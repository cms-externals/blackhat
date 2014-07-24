//============================================================================
// Name        : SplittingTest.cpp
// Author      : Carola F. Berger
// Version     : 8/25/08
// Copyright   : Resistance is futile!
// Description : tests ALL pure glue amplitudes up to 6 legs
//============================================================================




#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include<iostream>
#include<iomanip>
#include<vector>

#include "factorization.h"
#include "OneLoopHelAmpl.h"
#include "BH_utilities.h"
#include "settings_reader.h"

#include "gmp_r.h"

using namespace BH;
using namespace std;


#define _CHECK_SILENT2(check,pro,err) if ( check != 0 ){  _MESSAGE4(pro," NOT OK: ",check," collinear limits failed");err ++;	}


template <class T> T convert(const RGMP& g);

template <class T> std::complex<T> convert_complex(const std::complex<RGMP>& cg){
	return std::complex<T>(convert<T>(cg.real()),convert<T>(cg.imag()));
};

template <> RVHP convert<RVHP>(const RGMP& g){
	return RVHP(g.to_string().c_str());
};



void check_all(process PRO, color_structure cs, double mass,double tol, int arithm, int part, int & err)
{
 	int test;
	switch(arithm){
		case 0: test = TestCollinear(PRO, cs, R(mass),tol,part); break;
		case 1: test = TestCollinear(PRO, cs, RHP(mass),tol,part); break;
		case 2: test = TestCollinear(PRO, cs, RVHP(mass),tol,part); break;
		case 3: test = TestCollinear(PRO, cs, RGMP(mass),tol,part); break;
	}
	err+=test;
	return;
};
	 

int main(int argc,char *argv[]){


	string in_process;
	string in_cs;
	string in_mass;
	string in_seed;
	string in_arithm;
	string in_part;
	string in_tol;

	//settings::use_setting("SET_ALL_RAT_TO_ZERO no");
	//settings::use_setting("USE_KNOWN_FORMULAE no");
	settings::use_setting("USE_PARENT_FILES yes");
	settings::read_settings_from_file("BHsettings");


	if (argc < 15) { 
	    std::cout << endl; 
	    std::cout << "--------------------------------"<<endl;
	    std::cout << "* Usage:\n\n";
	    std::cout <<"\t./CollinearLimit\n\t\t -pro <process>\n\t\t -cs <color_structure>\n\t\t -mass <scale for limit>\n\t\t -seed <random-seed>\n\t\t -arithm <R,RHP,RVHP,RGMP>\n\t\t -part <full,cut,rat> \n\n"; // Inform the user of how to use the program
    	    std::cout << "\tExample: \n\t" ;
	    // std::cout << "\t ./CollinearLimit -pro qp_m_p_qbm_lm_lbp -cs leading_color -mass 0.00001 -seed 1111"<<endl; 
	    std::cout << "\t ./CollinearLimit -pro qp_m_p_qbm_lm_lbp -cs leading_color -mass 1e-80 -tol 1e-40 -seed 1111 -arithm RGMP -part cut"<<endl; 
	    std::cout << "\tVerbose: \n\t"; 
	    std::cout << "\t add 'splittingtest.cpp' in BH_debug.dat"<<endl; 
	    std::cout << "\tfor GMP: \n\t"; 
	    std::cout << "\t set e.g.: export BH_GMP_PRECISION=300 " <<endl; 
	    std::cout << "--------------------------------"<<endl;
    	    exit(0);
	} else {
    		for (int i = 1; i < argc; i++) { 
            if (strcmp(argv[i],       "-pro" )==0 && (i + 1 != argc) ) {
                in_process=string(argv[++i]);
            } else if (strcmp(argv[i],"-cs"  )==0 && (i + 1 != argc)) {
                in_cs=string(argv[++i]);
            } else if (strcmp(argv[i],"-mass")==0 && (i + 1 != argc)) {
                in_mass=string(argv[++i]);
            } else if (strcmp(argv[i],"-seed")==0 && (i + 1 != argc)) {
                in_seed=string(argv[++i]);
            } else if (strcmp(argv[i],"-arithm")==0 && (i + 1 != argc)) {
                in_arithm=string(argv[++i]);
            } else if (strcmp(argv[i],"-part")==0 && (i + 1 != argc)) {
                in_part=string(argv[++i]);
            } else if (strcmp(argv[i],"-tol")==0 && (i + 1 != argc)) {
                in_tol=string(argv[++i]);
	    }
    		}
	}

	// color_structure
	color_structure cs;
	if(in_cs=="leading_color") cs=leading_color;
	else if (in_cs=="nf")      cs=nf;
	else {cout<<"unknown cs; exit"<<endl; return 0;}

	// part
	int part;
	if(in_part=="full") part=0;
	else if(in_part=="cut") part=1;
	else if(in_part=="rat") part=2;
	else {cout<<"unknown part; must be either of [full,cut,rat]; exit"<<endl; return 0;}

	// arithm
	int arithm;
	if(in_arithm=="R") arithm=0;
	else if(in_arithm=="RHP") arithm=1;
	else if(in_arithm=="RVHP") arithm=2;
	else if(in_arithm=="RGMP") arithm=3;
	else {cout<<"unknown arithm; must be either of [R,RHP,RVHP,RGMP]; exit"<<endl; return 0;}



	// mass
	R mass;
	stringstream s_mass(in_mass);
	s_mass>>mass;

	// tol
	R tol;
	stringstream s_tol(in_tol);
	s_tol>>tol;

	// seed
	unsigned seed;
	stringstream s_seed(in_seed);
	s_seed>>seed;

	//process	
	vector<particle_ID> v_pro;
	size_t pos(0),pos1(0);
	while(pos1!=string::npos){
		pos1=in_process.find("_",pos);
		string str(in_process.substr(pos,pos1-pos));
		pos=pos1+1;
		if(str=="m") v_pro.push_back(m);
		else if(str=="p") v_pro.push_back(p);	
		else if(str=="m") v_pro.push_back(m);	
		//
		else if(str=="qp") v_pro.push_back(qp);	
		else if(str=="qbp") v_pro.push_back(qbp);	
		else if(str=="qm") v_pro.push_back(qm);	
		else if(str=="qbm") v_pro.push_back(qbm);	
		//	
		else if(str=="Gp") v_pro.push_back(Gp);	
		else if(str=="Gbp") v_pro.push_back(Gbp);	
		else if(str=="Gm") v_pro.push_back(Gm);	
		else if(str=="Gbm") v_pro.push_back(Gbm);	
		//	
		else if(str=="G2p") v_pro.push_back(G2p);	
		else if(str=="Gb2p") v_pro.push_back(Gb2p);	
		else if(str=="G2m") v_pro.push_back(G2m);	
		else if(str=="Gb2m") v_pro.push_back(Gb2m);	
		//	
		else if(str=="G3p") v_pro.push_back(G3p);	
		else if(str=="Gb3p") v_pro.push_back(Gb3p);	
		else if(str=="G3m") v_pro.push_back(G3m);	
		else if(str=="Gb3m") v_pro.push_back(Gb3m);	
		//
		else if(str=="lp") v_pro.push_back(lp);	
		else if(str=="lbp") v_pro.push_back(lbp);	
		else if(str=="lm") v_pro.push_back(lm);	
		else if(str=="lbm") v_pro.push_back(lbm);	
		//
		else if(str=="yp") v_pro.push_back(yp);	
		else if(str=="ym") v_pro.push_back(ym);	
		else {cout<<"need to update spectrum; exit"; exit(0);}
	}

	cout<<"========================================"<<endl;
	
	//return 0;
	
        //srand(unsigned(time(0)));
        srand(seed);
        //Read in the settings file
        //settings::read_settings_from_file("BHsettings");

	int err=0;
	check_all(process(v_pro),cs,mass,tol,arithm,part,err);
	
	cout<<"========================================"<<endl;
    
	cout<< "* pro: "<<process(v_pro)<<" cs: "<<cs<<" mass: "<<mass<<" err: "<<err<<endl;
	cout<<"\t tolerance: "<<in_tol<<endl;
	cout<<"\t seed: "<<seed<<endl;
	cout<<"\t arithm: "<<in_arithm<<endl;
	cout<<"\t part: "<<in_part<<endl;
	cout<<"========================================"<<endl;

    	return err;
}


