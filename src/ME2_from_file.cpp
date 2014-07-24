/*
 * ME2_from_file.cpp
 *
 *  Created on: Oct 4, 2009
 *      Author: harald
 */

#define _VERBOSE 1
#define MAP_LOOP_PRO 0 // [0,1] [map,don't map] to min set of amplitudes 

#ifndef BH_PUBLIC
#include "parent_diagrams.h"
#endif
#include <bitset>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include "constants.h"
#include "particles.h"
#include "BH_utilities.h"
#include "partitions.h"
#include "iterators.h"
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include "BH_debug.h"
#include "BH_typedefs.h"
#include "settings.h"
#include "eval_param.h"
#include "BHpath.h"
#include <sys/stat.h>
#include "color_algebra.h"
#include "cached_OLHA.h"
#include "cached_THA.h"
#include "assembly.h"
#include "primitive_ampl_map.h"
//#include "ME2.h"

#include "ME2_from_file.h"
#include "from_file.h"

#include "timing.h"
using namespace std;

namespace BH {
/*
string NoSpaces(string s)
{
  while(s.find(" ") != string::npos)
  {
    s.replace(s.find(" "), 1, "");
  }
  return s;
}
*/

bool color_info_match(vector<std::string> s1, vector<std::string>s2,int offset){

	int n=s1.size();
	for(int i=0;i<n;i++) if(s1[(i+offset)%n]!=s2[i]) return false;

	return true;
}

std::pair<int,int> pa_label_from_string(const std::string& s){
	int pos0=s.find('[',0);
	int pos=s.find('(',0);
	int pos2=s.find(')',0);

	std::string label=s.substr(pos0+1,pos0-pos-2);
	std::string index=s.substr(pos+1,pos2-pos-1);
	int pa_label; std::stringstream(label) >> pa_label;
	int ind; std::stringstream(index) >> ind;

//	std::pair<int,int> test_pair=make_pair(pa_label,ind);
//	cout<<"pa_label_from_string: "<<test_pair.first<<":"<<test_pair.second<<"\n";

	return make_pair(pa_label,ind);

}


color_constant  color_constant_from_string(const std::string& s){
	int pos_Nc=s.find('N',0);
    if(pos_Nc==string::npos){
        color_constant cc(0,0);
        return cc;
    }
	int pos_start=s.find('(',0)+1;
    //check for rational expression
    int pos_start_rat=s.find('(',pos_start)+1;
    int pos_end; pos_end=s.find(')',0);
    
    if(pos_start==string::npos){pos_start=0;}
	if(pos_end==string::npos){pos_end=s.size();}
    
    int pos_m=s.find('-',pos_start);
	if(pos_m==string::npos) pos_m=pos_start;

    int pos_ot=s.find('*',pos_start);
    int pos_ot_rat(0);
    int pos_frac=s.find('/',pos_start);
	int sign_pow(1);
	if(pos_frac!=string::npos){
		sign_pow=-1;
        pos_ot_rat=s.find('*',pos_frac);
        pos_ot=pos_frac;
	}
	int pos_pow=s.find('^',pos_start);


	std::string coeff_str=s.substr(pos_m,pos_ot-pos_m);
	int coeff; std::stringstream(coeff_str) >> coeff;

    int denom(1);
    if(pos_ot_rat!=string::npos && pos_ot_rat!=0){
    	std::string denom_str=s.substr(pos_start_rat,pos_pow-2-pos_start_rat);
        std::stringstream(denom_str) >> denom;
    }

    std::string pow_str=s.substr(pos_pow+1,pos_pow-s.size()-1);
	int pow(0); std::stringstream(pow_str) >> pow;pow=pow*sign_pow;

	//cout<<"color_constant_from_string: ("<<coeff<<","<<pow<<")\n";
	//cout<<s<<endl;

   /*
   cout<<"-------"<<endl;
   cout<< s <<endl;
   cout<< "denom "<<denom<<endl;
   cout<< "pow "<<pow<<endl;
   cout<< "coeff_str "<<coeff_str<<endl;
   cout<< "coeff "<<coeff<<endl;
   */
   //cout<< "factor "<<factor<<endl;
   //cout<< "res "<<factor*std::pow(NC,double(pow))<<endl;
   //cout<<"-------"<<endl;

   if(denom==1){
       color_constant cc(coeff,pow);
       return cc;
   }

       color_constant cc(coeff,denom,pow);
       return cc;

}


void helicity_orbit(const std::vector<plabel>& process_l,
	const std::vector<plabel>& process_r, 
	std::vector<std::vector<plabel> >& process_l_new,
	std::vector<std::vector<plabel> >& process_r_new,
    bool is_virtual){

	int n(process_l.size());
	std::vector<ph_type> hl(n),hr(n);
	std::vector<int> ll_map(n),rr_map(n);

	std::vector<int> helicity_map_l,helicity_map_r;
	std::vector<pair<int,int> > fermion_lines_l,fermion_lines_r;

	//identify fermion lines
	//identify independent helicities
	for(int i=0;i<n;i++){
		ll_map[process_l[i].ind()-1]=i;
		hl[i]=process_l[i];
		if(hl[i].is_a(gluon)){helicity_map_l.push_back(i);}
		else if(hl[i].is_a(quark)&&(!hl[i].is_anti_particle())){helicity_map_l.push_back(i);}
		else if(hl[i].is_a(gluino)&&(!hl[i].is_anti_particle())){helicity_map_l.push_back(i);}
		else if(hl[i].is_a(lepton)&&(!hl[i].is_anti_particle())){helicity_map_l.push_back(i);}
		else if(hl[i].is_a(photon)){helicity_map_l.push_back(i);}
		//if anti-fermion find matching fermion
		else if(hl[i].is_a(quark)||hl[i].is_a(gluino)||hl[i].is_a(lepton)){
			for(int pos=0;pos<n;pos++){
				if((*process_l[pos].type()==*hl[i].type())&&(process_l[pos].flavor()==hl[i].flavor())&&(pos!=i)){
					fermion_lines_l.push_back(make_pair(i,pos));
					break;	
				}
			};
		};
		//right process mapping
		//and quark line identification
		rr_map[process_r[i].ind()-1]=i;
		hr[i]=process_r[i];
		if((hr[i].is_a(quark)||hr[i].is_a(gluino)||hr[i].is_a(lepton))&&((hr[i].is_anti_particle()))){
			for(int pos=0;pos<n;pos++){
				if(((*process_r[pos].type()==*hr[i].type())&&(process_r[pos].flavor()==hr[i].flavor())&&(pos!=i))){
				fermion_lines_r.push_back(make_pair(i,pos));	
				}
			}
		}
	}

	//count fermions
	int nbr_fermions(fermion_lines_l.size());
	

	//independent helicities
	int all_hels(n-nbr_fermions);
	//sum over helicities
	bitset<16> helicities;
	//CHECK LENGTH OF ORBIT
	int max_helicity(int(std::pow(2.,all_hels)));
	for(int i=0;i<max_helicity;i++){
		helicities=i;
		for(int j=0;j<all_hels;j++){
			if(helicities[j])
				hl[helicity_map_l[j]].set_helicity(1);
			else 
				hl[helicity_map_l[j]].set_helicity(-1);
		};
		for(int j=0;j<nbr_fermions;j++){
				hl[fermion_lines_l[j].first].set_helicity(-hl[fermion_lines_l[j].second].helicity());
		}

		BH_DEBUG_PRINT(ll_map);
		BH_DEBUG_PRINT(rr_map);
        

    	//drop all plus and single minus amplitudes here
        int hel_sum(0), hel_l(0);
		for(int j=0;j<n;j++){
			hr[rr_map[j]].set_helicity(hl[ll_map[j]].helicity());
			hel_sum=hel_sum+hl[ll_map[j]].helicity();
		};
	    if(std::abs(hel_sum)>n-3) continue;
      
        // pick subclasses of helicity configurations to allow easier 
        // full color high multiplicity runs
        int keep=BH::settings::BH_interface_settings::s_same_helicity_projection;
        if(is_virtual && keep>-1){
            if(std::abs(hel_sum)!=keep) continue;
        }


        bool consistent(true);
		for(int j=0;(j<nbr_fermions)&&consistent;j++){
			consistent=consistent&&(hr[fermion_lines_r[j].first].helicity()==-hr[fermion_lines_r[j].second].helicity());
		}

		if(consistent){
			std::vector<plabel> plabel_l,plabel_r;
			for(int i=0;i<hl.size();i++){
				plabel_l.push_back(plabel(hl[i],process_l[i].ind()));	
				plabel_r.push_back(plabel(hr[i],process_r[i].ind()));	
			}
	
			process_l_new.push_back(plabel_l);
			process_r_new.push_back(plabel_r);
		};
    BH_DEBUG_PRINT(consistent);
	}
	
	return;
}


std::map<std::string,color_structure> init_cs_string_map(){
	std::map<std::string,color_structure> m;
	m.insert(std::make_pair("glue",glue));
	m.insert(std::make_pair("nf",nf));
	m.insert(std::make_pair("leading_color",leading_color));
	m.insert(std::make_pair("L",LT));
	m.insert(std::make_pair("R",RT));
	m.insert(std::make_pair("LL",LLT));
	m.insert(std::make_pair("LR",LRT));
	m.insert(std::make_pair("RL",RLT));
	m.insert(std::make_pair("RR",RRT));
	m.insert(std::make_pair("LLL",LLLT));
	m.insert(std::make_pair("RLL",RLLT));
	m.insert(std::make_pair("LRL",LRLT));
	m.insert(std::make_pair("LLR",LLRT));
	m.insert(std::make_pair("RRL",RRLT));
	m.insert(std::make_pair("RLR",RLRT));
	m.insert(std::make_pair("LRR",LRRT));
	m.insert(std::make_pair("RRR",RRRT));
	m.insert(std::make_pair("nfL",nfLT));
	m.insert(std::make_pair("nfR",nfRT));
	m.insert(std::make_pair("nfLL",nfLLT));
	m.insert(std::make_pair("nfLR",nfLRT));
	m.insert(std::make_pair("nfRL",nfRLT));
	m.insert(std::make_pair("nfRR",nfRRT));
//
	m.insert(std::make_pair("nfLLL",nfLLLT));
	m.insert(std::make_pair("nfRLL",nfRLLT));
	m.insert(std::make_pair("nfLRL",nfLRLT));
	m.insert(std::make_pair("nfLLR",nfLLRT));
	m.insert(std::make_pair("nfRRL",nfRRLT));
	m.insert(std::make_pair("nfRLR",nfRLRT));
	m.insert(std::make_pair("nfLRR",nfLRRT));
	m.insert(std::make_pair("nfRRR",nfRRRT));
	return m;
}


bool has_nf_loop(color_structure cs){
	switch(cs){
		case nf:case nfLT:case nfRT:case nfLLT:case nfLRT:case nfRLT:case nfRRT:
		case nfLLLT:case nfRLLT:case nfLRLT:case nfLLRT:case nfRRLT:case nfRLRT:case nfLRRT:case nfRRRT:{
			return true;} break;
		default: return false;
	}
}




// flip process and indices to take adventage of flip symmetry and of LT analytics
// needed to call correct analytics
// needed to speed up computation

void flip_pro_ind(process & pro, vector<int>& ind, double & sign, string & cs){
    //flips pro and ind: {12345...n} to {1 n ... 3 2} with (-1)^n
    //flips color structures
    size_t n=ind.size();
    if(n%2==1) sign=sign*(-1);
    reverse(ind.begin()+1,ind.end());
    vector<particle_ID> pro_flip(n);
    pro_flip[0]=pro[0];
    for(int i=1;i<n;i++) pro_flip[i]=pro[n-i];
    pro=process(pro_flip);

    flip_cs(cs);
}



void rot_qm_pro_ind(process & pro, vector<int>& ind, string & cs, bool qbm_first=true){
    //rotates pro and ind: {... qm ...} to {qm .... }
    //rotates pro and ind: {... qbm ...} to {qbm ... }
    //for case of (qm .... qbp) return ( qbm ... qp ) with R/L...T -> L/R ... T
    if(pro[0].is_a(quark)&&pro[0].helicity()==-1&&(qbm_first==pro[0].is_anti_particle())) return;

    size_t offset(0), qm_offset(0), qp_offset(0);
    while(!pro[offset].is_a(quark)) offset++;
        if(pro[offset].helicity()==-1){qm_offset=offset;}
        else {qp_offset=offset;};
        offset++;
    while(!pro[offset].is_a(quark)) offset++;
        if(pro[offset].helicity()==-1){qm_offset=offset;}
        else {qp_offset=offset;};

    rotate(ind.begin(),ind.begin()+qm_offset,ind.end());
    size_t n=ind.size();
    vector<particle_ID> pro_rot(n);
    for(int i=0;i<n;i++) pro_rot[i]=pro[(i+qm_offset)%n]; 
   
    if(qbm_first!=pro_rot[0].is_anti_particle()){
        qp_offset=(n+qp_offset-qm_offset)%n;
        pro_rot[0].set_is_antiparticle(qbm_first);
        pro_rot[qp_offset].set_is_antiparticle(!qbm_first);
        flip_cs_at(0,cs); //flip first L/R-turner
    }

    pro=process(pro_rot);
    return;
}


void rot_qm_pro_ind_ll(process & pro, vector<int>& ind, double & sign, string & cs){
    BH_DEBUG_MESSAGE("rot_qm_pro_ind_ll ");

    //rotates pro and ind: {... qm ...} to {qm .... }
    //rotates pro and ind: {... qbm ...} to {qbm ... }
    //returns true if first particle is qm !!!
    //for case of (qbm .... qbp) return ( qm ... qp ) with L...T -> R ... T
    //for leptons reorders to have (lm lbp)
    size_t n(ind.size());
    if(pro[0].is_a(quark)&&pro[0].helicity()==-1
            &&!pro[0].is_anti_particle()){
        if(pro[n-2].helicity()==-1) return;
        else{
            int ind_l(ind[n-2]);
            ind[n-2]=ind[n-1];
            ind[n-1]=ind_l;
            vector<particle_ID> pro_sort(n);
            for(size_t i=0;i<n-2;i++) pro_sort[i]=pro[i];
            pro_sort[n-2]=pro[n-1];
            pro_sort[n-1]=pro[n-2];
            pro_sort[n-2].set_is_antiparticle(false);
            pro_sort[n-1].set_is_antiparticle(true);
            pro=process(pro_sort); 
            sign*=-1;
            return;
        }
    } 

    //assume leptons at positions n-1 and n.
    vector<particle_ID> pro_short;
    vector<int> ind_short(ind.begin(),ind.end()-2); 
    for(size_t i=0;i<n-2;i++) pro_short.push_back(pro[i]);

    //reverse qcd part of process
    reverse(pro_short.begin(),pro_short.end());
    reverse(ind_short.begin(),ind_short.end());
    flip_cs(cs);
    if(n%2==1) sign=sign*(-1);

    //rotate qm into first position
    process pro_rot(pro_short);
    
    BH_DEBUG_MESSAGE7(pro_rot," ",ind_short," ",sign," ",cs);
    rot_qm_pro_ind(pro_rot,ind_short,cs,false);
    BH_DEBUG_MESSAGE7(pro_rot," ",ind_short," ",sign," ",cs);
    
    
    for(size_t i=0;i<n-2;i++) pro_short[i]=pro_rot[i];
    if(pro[n-2].helicity()==-1){
        sign*=-1;
        ind_short.push_back(ind[n-2]);
        ind_short.push_back(ind[n-1]);
        pro_short.push_back(pro[n-2]);
        pro_short.push_back(pro[n-1]);
    }
    else if(pro[n-2].helicity()==1){
        sign*=1;
        ind_short.push_back(ind[n-1]);
        ind_short.push_back(ind[n-2]);
        pro_short.push_back(pro[n-1]);
        pro_short.push_back(pro[n-2]);
        pro_short[n-2].set_is_antiparticle(false);
        pro_short[n-1].set_is_antiparticle(true);
    }

    pro=process(pro_short);
    ind=ind_short;
    return;
}

void rot_qm_pro_ind_glue(process & pro, vector<int>& ind, double & sign){
    
    size_t offset(0),n(ind.size());
    bool flip(false), rot(false);
    while(ind[offset]!=1) offset++;
    
    if(offset!=0){rot=true;}
    if(ind[(n+offset-1)%n]<ind[(n+offset+1)%n]) flip=true;
    
    BH_DEBUG_MESSAGE6("flip: ",flip," rot ",rot, ind[(n+offset-1)%n], ind[(n+offset+1)%n]);
    if(!flip && !rot) return;

    vector<particle_ID> pro_new(n);
    for(int i=0;i<n;i++) pro_new[i]=pro[(i+offset)%n];
    if(rot) rotate(ind.begin(),ind.begin()+offset,ind.end());

    if(!flip){pro=process(pro_new); return; }

    reverse(pro_new.begin()+1,pro_new.end());
    reverse(ind.begin()+1,ind.end());
    pro=process(pro_new);
    BH_DEBUG_MESSAGE(pro);
    if(n%2==1) sign*=-1;
    return;
}



void arrange_quarks(process& pro_in, string& cs){
    if(cs=="glue"||cs=="nf") return; //nothing needs to be done
    
    size_t extq(cs.size()-1);//number of external quarks
    vector<particle_ID> pro;
    if(cs==""){
        extq=0;
        for(size_t i=0;i<pro_in.n();i++){ 
            if(pro_in[i].is_a(quark)||pro_in[i].is_a(gluino)) extq++;
        }
        extq=extq/2;
    }
    else if(cs[0]=='n') extq-=2;
    BH_DEBUG_PRINT(extq);

    vector<bool> q_is_anti(extq,true); //vector to specify which cs/flavors have to be flipped
    size_t flavor(0);
    for(size_t i=0;i<pro_in.n();i++){
        pro.push_back(pro_in[i]);
        flavor=pro_in[i].flavor();
        BH_DEBUG_PRINT(pro[i]);
        if(flavor>1){
            if(q_is_anti[flavor-2]){
                q_is_anti[flavor-2]=false;
                if(!pro[i].is_anti_particle()) {
                    pro[i].set_is_antiparticle(true);   
                    flip_cs_at(flavor-1,cs);
                }
            }
            else if(pro[i].is_anti_particle()) pro[i].set_is_antiparticle(false);
        }
    }
    pro_in=process(pro);
    return;
}

bool has_leptons(process & pro){
    int n=pro.n();
    for(int i=0;i<n;i++){
        if(pro[n-i-1].is_a(lepton)) return true;
    }
    return false;
}




// ordering rules
// *) no quarks: - rotate leg=1 to first place 
//      - (leg-1) > (leg+1): flip
//      - first leg hel=-: conjugate

// sets first flavored gluino to  antifermion
// sets gluinos flavors to flavors according to their apperance in the plabel vector
// max 4 flavors
// for 2q2Gng2l processes
void arrange_quarks_leptons(vector<plabel>& pro_in){
   
    size_t n(pro_in.size());//number of partons
    size_t extq[]={0,0,0,0};//maximal 4 quark flavors 
    size_t flavor(2), flavor_loc(0);//flavor counter starts from 2; the lowest flavor
    
    for(size_t i=0;i<n;i++){ 
        if(pro_in[i].is_a(gluino)){
            flavor_loc=pro_in[i].flavor();
            switch(extq[flavor_loc]){
                case 0: {
                    extq[flavor_loc]=flavor;
                    pro_in[i].set_flavor(flavor); flavor++;
                    //pro_in[i].set_is_antiparticle(true);
                } break;
                default: {
                    pro_in[i].set_flavor(extq[flavor_loc]);
                    //pro_in[i].set_is_antiparticle(false);
                } break; 
            }
        }
    }
    return;
}



void read_primitive_amplitude(const std::string& input,
		const vector<plabel> & pro,
		const std::vector<int>& perm, 
		//CachedOLHA::partial_amplitude_cached* PA,
		//partial_amplitude_cached* PA,
        std::vector< std::pair<int,double> >& partial_ampl, 
        Squared_ME* ME,
		std::vector<kinematic_function*> prop_hel_fn){
	BH_DEBUG_MESSAGE("---------------------------");
	BH_DEBUG_MESSAGE("reading primitive amplitude");
	BH_DEBUG_MESSAGE2("input: ",input);
	BH_DEBUG(cout<<"perm: ";for(size_t i=0;i<perm.size();i++){ cout<<perm[i]<<" ";};cout<<endl;);

	int expr_start = input.find('{',1)+1; if(expr_start==string::npos) expr_start=0;
	int expr_sep_1 = input.find('|',expr_start); if(expr_sep_1==string::npos){BH_DEBUG_MESSAGE("empty partial amplitude data"); return;};
	int expr_sep_2 = input.find('|',expr_sep_1+1);
	int expr_end = input.find('}',expr_sep_2+2);if(expr_end==string::npos) expr_end=input.size()-1;

	//BH_DEBUG_MESSAGE(expr_start);
	//BH_DEBUG_MESSAGE(expr_sep_1);
	//BH_DEBUG_MESSAGE(expr_sep_2);
	//BH_DEBUG_MESSAGE(expr_end);

	//read color_structure
	string color_str=NoSpaces(input.substr(expr_start+1,expr_sep_1-expr_start-2));
	

	//read process
	vector<plabel> plprocess;
	std::string corner_str=input.substr(expr_sep_1+1,expr_sep_2-expr_sep_1-2);
	std::stringstream ss_corner(corner_str);
	while (ss_corner.good()){
		std::string partic;
		ss_corner >> partic;
		if ( !partic.empty() ){
//			BH_DEBUG_MESSAGE3("partic: \'",partic,"\'");
			plprocess.push_back(plabel_from_string(partic));
		}
	}
	BH_DEBUG_MESSAGE2("reading pro: ",plprocess);
	vector<int> ind_loc;
	for(int i=0;i<plprocess.size();i++){
		ind_loc.push_back(pro[perm[plprocess[i].ind()-1]-1].ind());
		plprocess[i].set_helicity(pro[perm[plprocess[i].ind()-1]-1].helicity());
		//plprocess[i].set_helicity(pro[perm[i]-1].helicity());
	}
	BH_DEBUG_MESSAGE2("input pro: ",pro);
	BH_DEBUG_MESSAGE2("primitive: ",plprocess);
	BH_DEBUG_MESSAGE2("ind:       ",ind_loc);


	//read color constant
	color_constant cc(0,0);
	int bracket_open=expr_sep_2+1,bracket_end;
	while ( bracket_open  != expr_end ){
		bracket_end = input.find('+',bracket_open+1);
		if (bracket_end == string::npos){
			bracket_end = expr_end;
		}
		std::string cc_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
//		_PRINT(bracket_end);
//		_PRINT(bracket_open);
//		BH_DEBUG_MESSAGE3("corner_str: \'",corner_str,"\'");
		std::stringstream ss_cc(cc_str);
		while (ss_cc.good()){
			std::string partic;
			ss_cc >> partic;
			if ( !partic.empty() ){
//				BH_DEBUG_MESSAGE3("partic: \'",partic,"\'");
				cc+=color_constant_from_string(partic);
			}
		}
		bracket_open = bracket_end;
	}

    

    // flip indices and process for LT terms for caching and for analytics
    double sign=1;
   
    short conjQ(0);
    //_MESSAGE9(pro_loc," ",ind_loc," ",sign," ",color_str," ",conjQ);


	for(int i=0;i<plprocess.size();i++) plprocess[i]=plabel(plprocess[i],ind_loc[i]);
#if MAP_LOOP_PRO==0 
	canonical_pro(plprocess,sign,conjQ,color_str);
#endif
    process pro_loc(plprocess);
    for(size_t i=0;i<plprocess.size();i++) ind_loc[i]=plprocess[i].ind();

    
    //_MESSAGE("---");
    //convert color structure string to cs 
    static std::map<std::string,color_structure> cs_string_map=init_cs_string_map();
    std::map<std::string,color_structure>::iterator it=cs_string_map.find(color_str);
    if(it==cs_string_map.end()){BH_DEBUG_MESSAGE3("Color structure ",color_str ," not known");}
	color_structure cs=it->second;
    BH_DEBUG_MESSAGE2("color structure: ",cs);

    //color projection
    //assumes leading color pprimitive in partial comes with factor of Nc or nf. Both contributions are kept
    switch (settings::BH_interface_settings::s_BH_color_mode){
	    case settings::BH_interface_settings::full_color: break;
        case settings::BH_interface_settings::leading_color: { 
	        if(has_nf_loop(cs)){ cc.project_to_Nc_powers(0,0);}
            else{                cc.project_to_Nc_powers(1,1);}
        } break; 
		case settings::BH_interface_settings::full_minus_leading_color:{
            if(has_nf_loop(cs)){ cc.project_to_Nc_powers(-1,-10);}
            else{                cc.project_to_Nc_powers(0,-10);}
        } break; 
    }
    if(cc.is_zero()) return;


    //weight for interference term from color_constant
	/*************************
	 * here nf-terms get a weight "nf" with nf the global setting of numbers of flavors
	 * **********************/
	double cc_double(cc.eval());
	if(has_nf_loop(cs)) cc_double *=(BH::settings::BH_interface_settings::s_nf);
	BH_DEBUG_MESSAGE2("weight of primitive: ",cc_double);

	BH_DEBUG_MESSAGE("++++++++++++++++++++++++++++");
	BH_DEBUG_MESSAGE(process(plprocess));
	BH_DEBUG_MESSAGE(cs);
	BH_DEBUG_MESSAGE(ind_loc);
	BH_DEBUG_MESSAGE(cc_double);
	BH_DEBUG_MESSAGE("++++++++++++++++++++++++++++");

    if(cc_double==0.0) return;

	//insert primitive amplitude into partial
	size_t loop_ampl=ME->add(pro_loc,ind_loc,cs,conjQ,prop_hel_fn);
    partial_ampl.push_back(make_pair(loop_ampl,cc_double*sign));
	
	return;
}


double ordered_permutation_orbit(const vector<plabel>& plabels,
        const vector<int>& perm,
        vector<vector<int> >& perm_orbit,
        vector<std::string >& color_info
        );

void read_cross_term(const std::string& input,
	const std::vector<std::pair<int, int> >& pa_labels, 
	const vector<int>& perm,
	construction_cache * constr_cache, 
	Squared_ME* ME,
	bool born_or_virt){

	BH_DEBUG_MESSAGE("---------------------");
    BH_DEBUG_MESSAGE2("   * interference: ",input);
    BH_DEBUG(cout<<"      permutation: ";for(size_t i=0;i<perm.size();i++) cout<<perm[i]<<" "; cout<<endl;);

	int expr_start_l = input.find('[',0);
	int expr_middle_l = input.find('|',expr_start_l);
	int expr_end_l = input.find('*',expr_start_l)-2;
	int expr_start_r = input.find('[',expr_end_l);
	int expr_middle_r = input.find('|',expr_start_r);
	int expr_end_r = input.find('*',expr_start_r)-2;
	int expr_start_w = input.find('[',expr_end_r);
	int expr_end_w = input.find(']',expr_start_w);
	int bracket_open=0, bracket_end=-1;

	vector<plabel> process_l,process_r;
	vector<vector<pair<int,int> > > coupling_l,coupling_r;

	bracket_open = expr_start_l;bracket_end = expr_middle_l;
	std::string process_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
	bracket_open = expr_middle_l;bracket_end = expr_end_l;
	std::string coupling_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
	BH_DEBUG_MESSAGE2("      coupling:  ",coupling_str);
	
	//vector<size_t> color_info;
	std::vector<std::string > color_info_l, color_info_r;
	

	/*
	switch(born_or_virt){
		case true: process_coupling_from_string(process_str,coupling_str,process_l,coupling_l,color_info_l); break;
		case false : partial_process_coupling_from_string(process_str,coupling_str,process_l,coupling_l,color_info_l); break;
	}
	*/
	partial_process_coupling_from_string(process_str,coupling_str,process_l,coupling_l,color_info_l);
    

	bracket_open = expr_start_r;bracket_end = expr_middle_r;
	process_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
	bracket_open = expr_middle_r;bracket_end = expr_end_r;
	coupling_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
	partial_process_coupling_from_string(process_str,coupling_str,process_r,coupling_r,color_info_r); 

	//read color constant
	bracket_open = expr_start_w;bracket_end = expr_end_w;
	color_constant cc(0,0);
	while ( bracket_open  != expr_end_w ){
		bracket_end = input.find('+',bracket_open+1);
		if (bracket_end == string::npos){
			bracket_end = expr_end_w;
		}
		std::string cc_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
        //_PRINT(cc_str);
		std::stringstream ss_cc(cc_str);
		while (ss_cc.good()){
			std::string partic;
			ss_cc >> partic;
			if ( !partic.empty() ){
				cc+=color_constant_from_string(partic);
			}
		}
		bracket_open = bracket_end;
	}
    //_PRINT(cc);

	vector<vector<int> > perm_orbit;
    double weight(1.);
	if(BH::settings::BH_interface_settings::s_born_final_state_gluon_symmetrization){
        //use this for normal operation mode
        permutation_orbit(process_l,perm,perm_orbit,color_info_l);
    }
    else{
        weight=ordered_permutation_orbit(process_l,perm,perm_orbit,color_info_l);
    }
	
    size_t n(perm_orbit.size());
	if(BH::settings::BH_interface_settings::s_use_symmetrized_assembly_files){
		BH_DEBUG_MESSAGE("      WARNING: expecting symmetrized input data");
		n=0;
	}

#if 0
    cout<<"process_l: ";
    for (int i=0;i<process_l.size();i++) cout<<process_l[i]<<" ";
    cout<<endl;
#endif


	//gernerate possible helicities of primitive amplitudes
	std::vector<std::vector<plabel> > process_l_new, process_r_new;
	
    BH_START_TIMER(helicity_orbit);
    helicity_orbit(process_l,process_r,process_l_new,process_r_new,!born_or_virt);
    BH_STOP_TIMER(helicity_orbit);
	int ampl_l(-1),ampl_r(-1);

	//weight for interference term from color_constant
	double cc_double(cc.eval());
	BH_DEBUG_MESSAGE2("cc_double: ",cc_double);
    //for symmerization of born pure gluon amplitudes
    cc_double*=weight;

	multi_precision_fraction factor(1),factor_virt(1);
	//charge conventions
	for(int i=0;i<process_l_new.size();i++){
        double sign_l(1),sign_r(1); //sign for flip symmetry and charge conjugation switch
		switch(born_or_virt){
			case true: ampl_l=constr_cache->new_cross_term_entry(process_l_new[i],perm,coupling_l,color_info_l,sign_l); break; //born
			case false: ampl_l=constr_cache->new_loop_cross_term_entry(process_l_new[i],perm,coupling_l,color_info_l);break; //loop
		};
		if(ampl_l!=-1){
			ampl_r=constr_cache->new_cross_term_entry(process_r_new[i],perm,coupling_r,color_info_r,sign_r);	
			if(ampl_r!=-1){
				switch(born_or_virt){
					case true:{
						//born
                            ME->add_tree(ampl_l,ampl_r,R(factor)*cc_double*sign_l*sign_r);
                        for(size_t j=0;j<n;j++){ //gluon permutation sum
                            sign_l=1; sign_r=1; 
							ampl_l=constr_cache->new_cross_term_entry(process_l_new[i],perm_orbit[j],coupling_l,color_info_l,sign_l);
							ampl_r=constr_cache->new_cross_term_entry(process_r_new[i],perm_orbit[j],coupling_r,color_info_r,sign_r);	
                            ME->add_tree(ampl_l,ampl_r,R(factor)*cc_double*sign_l*sign_r);
						};
					}; break;
					case false:{ 
						//loop
						BH_DEBUG_MESSAGE2("add: ",cc_double*sign_r);
						ME->add_loop(ampl_l,ampl_r,R(factor)*cc_double*sign_r);
						for(size_t j=0;j<n;j++){ //gluon permutation sum
                            sign_l=1; sign_r=1;
							ampl_l=constr_cache->new_loop_cross_term_entry(process_l_new[i],perm_orbit[j],coupling_l,color_info_l);
							ampl_r=constr_cache->new_cross_term_entry(process_r_new[i],perm_orbit[j],coupling_r,color_info_r,sign_r);	
						    ME->add_loop(ampl_l,ampl_r,R(factor)*cc_double*sign_r);
						};
					}; break;
				}
			}
		}
	}
}




void read_pa_labels(const std::string& input, vector<vector<std::pair<int,int> > >& pa_labels){
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
		pa_labels.push_back(std::vector<std::pair<int,int> >());

		std::string corner_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);

//		_PRINT(bracket_end);
//		_PRINT(bracket_open);
//		BH_DEBUG_MESSAGE3("corner_str: \'",corner_str,"\'");
		std::stringstream ss_corner(corner_str);

		while (ss_corner.good()){
			std::string partic;
			ss_corner >> partic;
			if ( !partic.empty() ){
//				BH_DEBUG_MESSAGE3("partic: \'",partic,"\'");
				pa_labels.back().push_back(pa_label_from_string(partic));
			}
		}


		bracket_open = bracket_end;

	}

}

struct plabel_has_same_type_and_ap {
	bool operator()(const plabel& p1, const particle_ID& p2){
		return (p1.type() == p2.type()) && (p1.is_anti_particle() == p2.is_anti_particle() && p1.flavor() == p2.flavor() ); }
};



// check if color traces have lowest label in first position
bool sorted_permutation(const vector<int> perm, const vector<int> trace_begin, const vector<int> trace_end){
    vector<int > perm_loc=perm;
    
    BH_DEBUG(
        cout<<"trace begin ";
        for(int i=0;i<trace_begin.size();i++) cout<<trace_begin[i]<<" ";
        cout<<endl<<"trace end   ";
        for(int i=0;i<trace_end.size();i++) cout<<trace_end[i]<<" ";
        cout<<endl;
    )

    for(int i=0;i<trace_begin.size();i++){
        sort(perm_loc.begin()+trace_begin[i],perm_loc.begin()+trace_end[i]);
        BH_DEBUG_PRINT(perm_loc);

        // make sure that there is no overcounting if we have two traces have same length
        if(trace_begin.size()==2&&(
                trace_end[0]-trace_begin[0]==
                trace_end[1]-trace_begin[1])){
            if((perm_loc[trace_begin[0]]>perm_loc[trace_begin[1]])||
               (perm[trace_begin[i]]!=perm_loc[trace_begin[i]])) return false;
        }
        else {
            if(perm[trace_begin[i]]!=perm_loc[trace_begin[i]]) return false;
        }
   
        BH_DEBUG(
            cout<<"perm ";
            for(int i=0;i<perm.size();i++) cout<<perm[i]<<" ";
            cout<<endl<<"perm_loc   ";
            for(int i=0;i<perm_loc.size();i++) cout<<perm_loc[i]<<" ";
            cout<<endl;
        )
    }   
    return true;
}

/* born only
 *
 * switch off permutation sum when using this function
 * 
 * 1) need to make sure that canonical permutation 
 *    is not added in perm_orbit 
 */


double ordered_permutation_orbit(const vector<plabel>& plabels,
        const vector<int>& perm,
        vector<vector<int> >& perm_orbit,
        vector<std::string >& color_info
        ){
	vector<int> orbit,orbit_ref,trace_begin,trace_end,trace_info;
	size_t m=plabels.size();
    vector<int> perm_new(m);perm_new=perm;	

    //invert index and collect gluon positions in to orbit vector
    //std::vector<std::pair<int,int> > pa_labels_loc(m);
	size_t ind1,ind2;
    for(int i=0;i<m;++i){
        if(plabels[i].is_a(gluon)){
            orbit_ref.push_back(plabels[i].ind()-1);
            if(perm[plabels[i].ind()-1]==1){
                ind1=plabels[i].ind()-1;
            }
            else if(perm[plabels[i].ind()-1]==2){
                ind2=plabels[i].ind()-1;
            }
            else{
                orbit.push_back(plabels[i].ind()-1);
            }
        }
    }
	if(orbit.size()==0) return double(1);

    BH_DEBUG_PRINT(plabels);
    BH_DEBUG(
            for(int i=0;i<plabels.size();i++) cout<<plabels[i].ind()<<" ";
            cout<<endl;)
    
    //orbit_ref=orbit;

	size_t n(orbit_ref.size());
    size_t n1(n), n2_min(1), n2_max(n);
    if(plabels.size()==n){
        n1=1;
        n2_min=2;
    }

	std::vector<int> orbit_loc(n);
	//for non-trivial partials we have to drop some permutations

            for(size_t i=0;i<n1;++i){
                for(size_t j=n2_min;j<n2_max;++j){
                    if(i!=j){
                       int shift(0);
                       for(size_t k=0;k<n;++k){
                           if(k==i){
                                orbit_loc[k]=ind1; 
                                --shift;
                            }
                            else if(k==j){
                                orbit_loc[k]=ind2; 
                                --shift;
                            }
                            else{
                                //_PRINT(k+shift);
                                orbit_loc[k]=orbit[k+shift];
                            }
                        }     
           
            for(size_t i=0;i<n;++i) perm_new[orbit_ref[i]]=perm[orbit_loc[i]];
    	
            perm_orbit.push_back(perm_new);
                    
                    }
                }
            }

    double weight(1);
    for(int i=1;i<n-1;i++) weight*=double(i);

    //for(size_t i=0;i<perm_orbit.size();++i) _PRINT(perm_orbit[i]);

    BH_DEBUG_PRINT(perm_new);
	return weight;
}

void permutation_orbit(const vector<plabel>& plabels,
        const vector<int>& perm,
        vector<vector<int> >& perm_orbit,
        vector<std::string >& color_info
        ){
	vector<int> orbit,orbit_ref,trace_begin,trace_end,trace_info;
	size_t m=plabels.size();
    vector<int> perm_new(m);perm_new=perm;	

    //invert index and collect gluon positions in to orbit vector
    //std::vector<std::pair<int,int> > pa_labels_loc(m);
	for(int i=0;i<m;i++){
        if(plabels[i].is_a(gluon)) orbit.push_back(plabels[i].ind()-1);
    }
	if(orbit.size()==0) return;

    BH_DEBUG_PRINT(plabels);
    BH_DEBUG(
            for(int i=0;i<plabels.size();i++) cout<<plabels[i].ind()<<" ";
            cout<<endl;)

    //arrange for supressing redundant permutations of partials
    if(color_info.size()>plabels.size()){
        // trace_info should have first, in each trace break or "," position.
        int shift(0),glue(-1);
        trace_info.push_back(0);
        for(int i=0;i<color_info.size();i++){
            if(color_info[i]==","||color_info[i]==";"){
                trace_info.push_back(i-shift);
                shift++;
            }
        }
        trace_info.push_back(m);
        for(int i=0;i<plabels.size();i++){
            if(plabels[i].is_a(gluon)){
                glue++;
                for(int j=0;j<trace_info.size();j++){
                    if(trace_info[j]==i){
                    trace_begin.push_back(glue);
                    trace_end.push_back(glue+trace_info[j+1]-trace_info[j]);
                    }
                }
            }
        }
    if(trace_end.size()>0){
    BH_DEBUG(
            _PRINT(perm_new);
            _PRINT(color_info);
            _PRINT(color_info.size());
            _PRINT(trace_info);
            _PRINT(orbit);
            _PRINT(trace_begin);
            _PRINT(trace_end);
    );
    }
    }
    else if(orbit.size()==plabels.size()){
        //WARNING: in case of pure glue amplitudes we can fix first leg for born and leading-color only!!
        orbit.erase(orbit.begin());
    }
    orbit_ref=orbit;

	size_t n(orbit.size());
	std::vector<int>::iterator begin,end;
	begin=orbit.begin();end=orbit.end();
	//for non-trivial partials we have to drop some permutations
    if(trace_begin.size()==0){
        while(next_permutation(begin,end)){
                for(int i=0;i<n;i++) perm_new[orbit_ref[i]]=perm[orbit[i]];
    		    perm_orbit.push_back(perm_new);
        }
    }
    else {
        while(next_permutation(begin,end)){
            for(int i=0;i<n;i++){
    		    	perm_new[orbit_ref[i]]=perm[orbit[i]];
    		}
            //reduce permutations    
            if(sorted_permutation(orbit,trace_begin,trace_end)){
                perm_orbit.push_back(perm_new);
            BH_DEBUG(
                    cout<<"perm:     ";
                    for(int i=0;i<perm.size();i++) cout<< perm[i];
                    cout<<endl;
                    cout<<"perm_new: ";
                    for(int i=0;i<perm_new.size();i++) cout<< perm_new[i];
                    cout<<endl;
                    cout <<"keep permutation " <<sorted_permutation(perm_new,trace_begin,trace_end)<<endl;
                    )
            }
        }
    }
    BH_DEBUG_PRINT(perm_new);
	return;
}

bool find_pa_labels_match(vector<pair<int,int> > ref_pa_labels,vector<pair<int,int> > pa_labels,vector<int>& perm){
	assert( perm.size() == pa_labels.size());
	if (pa_labels.size() != ref_pa_labels.size()){
		return false;
	}
	sort(pa_labels.begin(),pa_labels.end());
	sort(ref_pa_labels.begin(),ref_pa_labels.end());

	for (int pos=0;pos<ref_pa_labels.size();pos++){
		if (ref_pa_labels[pos].first==pa_labels[pos].first){
			perm[ref_pa_labels[pos].second-1]=pa_labels[pos].second;
		}
		else{ return false;}
	}
	return true;
}

//cross_term_entry
bool operator<(const cross_term_entry& cte1, const cross_term_entry& cte2){
	if(cte1.m_color_info.size()==0){
		if(cte1.m_coupling.size()==0&&cte2.m_coupling.size()==0){
			return(cte1.m_plabel<cte2.m_plabel);
		}	
		else{	
			if(cte1.m_plabel<cte2.m_plabel) return true;
            else if(cte1.m_plabel>cte2.m_plabel) return false;
            else return (cte1.m_coupling<cte2.m_coupling);
		}
	}
	else{
		if(cte1.m_coupling.size()==0&&cte2.m_coupling.size()==0){
			if(cte1.m_plabel<cte2.m_plabel) return true;
            else if(cte1.m_plabel>cte2.m_plabel) return false;
            else return (cte1.m_color_info<cte2.m_color_info);
		}	
		else{	
			if(cte1.m_plabel<cte2.m_plabel) return true;
            else if(cte1.m_plabel>cte2.m_plabel) return false;
            else {
                if(cte1.m_coupling<cte2.m_coupling) return true;
                else if (cte1.m_coupling>cte2.m_coupling) return false;
                else return (cte1.m_color_info<cte2.m_color_info);
            }
        }
	}
};

bool operator==(const cross_term_entry& cte1, const cross_term_entry& cte2){
	if(cte1.m_color_info.size()!=cte2.m_color_info.size()){
		return false;
	}
	else if(cte1.m_color_info.size()==0){
		return ((cte1.m_plabel.size()==cte2.m_plabel.size())&&
			(cte1.m_plabel==cte2.m_plabel)&&
			(cte1.m_coupling==cte2.m_coupling));
	}
	else{
		return ((cte1.m_plabel.size()==cte2.m_plabel.size())&&
			(cte1.m_plabel==cte2.m_plabel)&&
			(cte1.m_coupling==cte2.m_coupling)&&
			(cte1.m_color_info==cte2.m_color_info)
			);
	}			
};

bool compr(const cross_term_entry* cte1, const cross_term_entry* cte2){
	return (*cte1<*cte2);
};
bool equal(const cross_term_entry* cte1, const cross_term_entry* cte2){
	if(cte1==0||cte2==0){return false;}
	return (*cte1==*cte2);
};




/*

void process_coupling_from_string(const std::string& process_str, 
		const std::string& coupling_str,
		vector<plabel>& process,
		std::vector<std::vector<pair<int,int> > >& coupling){
	{	
//	std::string process_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
	BH_DEBUG_MESSAGE3("process_str: \'",process_str,"\'");
	std::stringstream ss_process(process_str);
	while (ss_process.good()){
		std::string partic;
		ss_process >> partic;
		if ( !partic.empty() ){
			BH_DEBUG_MESSAGE3("partic: \'",partic,"\'");
			process.push_back(plabel_from_string(partic));
		}
	}
	}
	{
	//std::string coupling_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
	BH_DEBUG_MESSAGE3("coupling_str: \'",coupling_str,"\'");
	std::stringstream ss_coupling(coupling_str);
	vector<pair<int,int> > single_coupling;
	bool has_coupling(false);
	while (ss_coupling.good()){
		std::string partic;
		ss_coupling >> partic;
		if ( !partic.empty() ){
			has_coupling=true;			
			if(partic==";"){
				BH_DEBUG_MESSAGE3("partic: \'",partic,"\'");
				coupling.push_back(single_coupling);
				single_coupling.clear();
			}
			else{
				single_coupling.push_back(pa_label_from_string(partic));
			}
		}
	}
	if(has_coupling) coupling.push_back(single_coupling);
	}
	return;
}
*/

void partial_process_coupling_from_string(const std::string& process_str, 
		const std::string& coupling_str,
		vector<plabel>& process,
		vector<vector<pair<int,int> > >& coupling,
//		vector<size_t> & color_info
		vector<std::string >& color_info
		){
	{
	color_info.clear();	
	bool pos(true);	
//	std::string process_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
	//BH_DEBUG_MESSAGE3("process_str: \'",process_str,"\'");
	std::stringstream ss_process(process_str);
	while (ss_process.good()){
		std::string partic;
		ss_process >> partic;
		if ( !partic.empty() ){
//				BH_DEBUG_MESSAGE3("partial_process_coupling: partic: \'",partic,"\'");
			if(partic[0]==','||partic[0]==';'){
                //H pos=false;
                pos=true;
                color_info.push_back(partic);
            }
			else{
				if(pos) color_info.push_back("_");
				pos=true;
//				BH_DEBUG_MESSAGE3("partic: \'",partic,"\'");
				process.push_back(plabel_from_string(partic));}
			}
		}
//	BH_DEBUG_MESSAGE2("color info in process match:", color_info.size());
//	for(int i=0;i<color_info.size();i++) cout<<color_info[i]; cout<<endl;
	}
	{
	//std::string coupling_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);
	//BH_DEBUG_MESSAGE3("coupling_str: \'",coupling_str,"\'");
	std::stringstream ss_coupling(coupling_str);
	vector<pair<int,int> > single_coupling;
	bool has_coupling(false);
	while (ss_coupling.good()){
		std::string partic;
		ss_coupling >> partic;
		if ( !partic.empty() ){
//			BH_DEBUG_MESSAGE3("partic: \'",partic,"\'");
			has_coupling=true;			
			if(partic==";"){
					coupling.push_back(single_coupling);
					single_coupling.clear();
			}
			else{
				single_coupling.push_back(pa_label_from_string(partic));
			}
		}
	}
	if(has_coupling) coupling.push_back(single_coupling);
	
	}
	return;
}


cross_term_entry::cross_term_entry(const std::vector<plabel>& labels,
		const std::vector<int>& perm,
		const std::vector<std::vector<pair<int,int> > >& coupling){

	std::vector<pair<int,int> > single_coupling;
/*	for(int i=0;i<labels.size();i++){
		cout<<labels[i];
	}
	cout<<endl;
*/
	for(int i=0;i<labels.size();i++){
		m_plabel.push_back(plabel(labels[i],perm[labels[i].ind()-1]));
	}	
	for(size_t h=0;h<coupling.size();h++){
		single_coupling.clear();
		for(int i=0;i<coupling[h].size();i++){
			single_coupling.push_back(make_pair(coupling[h][i].first,perm[coupling[h][i].second-1]));
		}
		m_coupling.push_back(single_coupling);
	}
	//_PRINT(m_plabel);	
	BH_DEBUG(cout<<"cross term: ";
		for(int i=0;i<labels.size();i++) cout<<m_plabel[i]<<" "; cout<<"  ";
		for(size_t h=0;h<m_coupling.size();h++){
			for(int j=0;j<4;j++) cout<<"("<<m_coupling[h][j].first<<","<<m_coupling[h][j].second<<") "; cout<<endl;
		};)
};

cross_term_entry::cross_term_entry(const std::vector<plabel>& labels,
		const std::vector<int>& perm,
		const std::vector<std::vector<pair<int,int> > >& coupling,
		//const std::vector<size_t>& color_info){
		const vector<std::string >& color_info){
/*	for(int i=0;i<labels.size();i++){
		cout<<labels[i];
	}
	cout<<endl;
*/
	std::vector<pair<int,int> > single_coupling;
	for(int i=0;i<labels.size();i++){
		m_plabel.push_back(plabel(labels[i],perm[labels[i].ind()-1]));
	}	

	for(size_t h=0;h<coupling.size();h++){
		single_coupling.clear();
		for(int i=0;i<coupling[h].size();i++){
			single_coupling.push_back(make_pair(coupling[h][i].first,perm[coupling[h][i].second-1]));
		}
		m_coupling.push_back(single_coupling);
	}
	
	m_color_info=color_info;
	//_PRINT(m_plabel);	
	BH_DEBUG(cout<<"cross term (loop): ";
		for(int i=0;i<labels.size();i++) cout<<m_plabel[i]<<" "; cout<<"  ";
		for(size_t h=0;h<m_coupling.size();h++){
			for(int j=0;j<4;j++) cout<<"("<<m_coupling[h][j].first<<","<<m_coupling[h][j].second<<") "; cout<<endl;
		}; )
};


std::vector<kinematic_function*> cross_term_entry::m_coupling_function(){
	
//	_PRINT(coupling_function(m_coupling,m_plabel));
	
	return coupling_function(m_coupling,m_plabel);
	}


//cross_term
kinematic_function* coupling_function_4(std::vector<pair<int,int> >& coupling, const std::vector<plabel>& labels){
	
	int coupling_type(coupling.size());
//	_PRINT(coupling_type);

	//photonZW
	int photonZW;
	switch(coupling_type){
		case 0:	return 0;break;
		case 3: photonZW=0; break;
		case 4:	//check if neutrino or charged leptons
			int lepton1,lepton2;
			lepton1=coupling[2].first;
			lepton2=coupling[3].first;
			switch(lepton1%2){
				case 0: {
					switch((-lepton2)%2){
						case 0: photonZW=2;break; //Z only
						case 1: photonZW=3;break; //W case
					}
				} break;
				case 1:{
					switch((-lepton2)%2){
						case 0: photonZW=3;break; //W case
						case 1: photonZW=1;break; //Z+photon case
					};
				} break;
			}; break;
	};

	//leading_vect_ax needs to be implemented properly; dropping vect_ax here
	int leading_vect_ax(0);

	
	//index i & j
	int ind_i(1),ind_j(1),label_quark(coupling[0].second-1),label_lepton(coupling[2].second-1);

//	_PRINT(label_quark);
//	_PRINT(label_lepton);

	std::vector<ph_type> ql_ph_type;
	if(photonZW!=0) {
		ind_i=coupling[2].second;
		ind_j=coupling[3].second;

		for(int j=0;j<labels.size();j++){
			if(labels[j].ind()==label_quark+1){	
				ql_ph_type.push_back(labels[j]);
                if(coupling[0].first>0) ql_ph_type.back().set_is_antiparticle(false);
                else ql_ph_type.back().set_is_antiparticle(true);
			}
		}
		for(int j=0;j<labels.size();j++){
			if(labels[j].ind()==label_lepton+1){
				ql_ph_type.push_back(labels[j]);
                if(coupling[2].first>0) ql_ph_type.back().set_is_antiparticle(false);
                else ql_ph_type.back().set_is_antiparticle(true);
			}
		}
	}

	//up_down_quark
	bool up_down_quark(false);
	//W keeps default value of up_down_quark=0
	if(photonZW!=3){up_down_quark=(coupling[0].first%2==1);}

	BH_DEBUG_PRINT(labels);
	//BH_DEBUG(for(int i=0;i<coupling.size();i++){
	//	cout<<"("<<coupling[i].first<<","<<coupling[i].second<<")";
	//};
	//cout<<endl;);

	BH_DEBUG_PRINT(up_down_quark);
	BH_DEBUG_PRINT(label_quark);
	BH_DEBUG_PRINT(photonZW);
	BH_DEBUG_PRINT(leading_vect_ax);
	BH_DEBUG_PRINT(ind_i);
	BH_DEBUG_PRINT(ind_j);
	BH_DEBUG_PRINT(ql_ph_type);

	prop_hel_fn* local_prop_hel_fn=(new prop_hel_fn(up_down_quark,photonZW,leading_vect_ax,ind_i,ind_j,ql_ph_type));

	//check if prop_hel_fn will be zero; if yes return zero pointer
	if(local_prop_hel_fn->selection_rule_is_zero()){
//		std::cout<<"set to zero:"<<std::endl;
		delete local_prop_hel_fn;return 0;
	};
//		std::cout<<"not set to zero:"<<std::endl;
	return local_prop_hel_fn;
};


//coupling function for two on-shell photons coupled to single quark line
kinematic_function* coupling_function_diphoton(std::vector<pair<int,int> >& coupling, const std::vector<plabel>& labels){
	
	//returns the charge squared of quark line
	prop_hel_fn_diphoton* local_prop_hel_fn=(new prop_hel_fn_diphoton(coupling[0].first));

	return local_prop_hel_fn;
};

vector<kinematic_function*> coupling_function(std::vector<std::vector<pair<int,int> > >& coupling, const std::vector<plabel >& labels){

	vector<kinematic_function* > coupling_prod;

	switch(coupling.size()){
		case 1: {
			if((coupling)[0].size()<5 && !(coupling[0].size()==4 && coupling[0][2].first == 22)){
				//single V-boson cases
				coupling_prod.push_back(coupling_function_4(coupling[0],labels));
				return coupling_prod;
			}
			else if(coupling[0].size()==4 && coupling[0][2].first==22 && coupling[0][3].first==22){
				//di photon 4-point coupling
				coupling_prod.push_back(coupling_function_diphoton(coupling[0],labels));
				return coupling_prod;
			}
			else{
				_MESSAGE("can only do single V-boson interaction.");
			}
				} break;
		case 2: {

			if((coupling)[0].size()<5 &&
				(coupling)[1].size()<5 ){
				coupling_prod.push_back(coupling_function_4(coupling[0],labels));
				coupling_prod.push_back(coupling_function_4(coupling[1],labels));
				return coupling_prod;
			}
			else{
				_MESSAGE("can only do single V-boson interaction.");
			}
				} break;
				
		}
};









coupling_process::coupling_process(std::vector<std::vector<pair<int,int> > >& coupling, const std::vector<plabel>& labels): m_coupling(coupling) {
    //reduce labels to the once that match the particles in couplings
	std::vector<plabel> single_plabels;
    for(int h=0;h<coupling.size();h++){
		single_plabels.clear();
    	for(int i=0;i<coupling[h].size();i++){
        	for(int j=0;j<labels.size();j++){
           		if(coupling[h][i].second==labels[j].ind()) single_plabels.push_back(labels[j]);
        	}
    	}
		m_plabels.push_back(single_plabels);
	}
};

bool operator<(const coupling_process& cpro1, const coupling_process& cpro2){
    if(cpro1.m_coupling<cpro2.m_coupling) return true;
    else if(cpro1.m_coupling>cpro2.m_coupling) return false;
    else if(cpro1.m_coupling==cpro2.m_coupling) return (cpro1.m_plabels<cpro2.m_plabels);
    
    _MESSAGE("problem in comparison of coupling_process");
};




int construction_cache::new_cross_term_entry(const std::vector<plabel>& labels_in,
	const std::vector<int>& perm_in,
	const std::vector<std::vector<pair<int,int> > >& coupling_in,
	const std::vector<std::string > & color_info, 
    double & sign){

/* for now color info is not used and we do not read 
	the partial files in for this. */


    //include permutation in labels
    vector<plabel> labels;
    vector<int> perm;
	std::vector<std::vector<pair<int,int> > > coupling;
	for(size_t i=0;i<labels_in.size();i++){
		labels.push_back(plabel(labels_in[i],perm_in[labels_in[i].ind()-1]));
	}	
	std::vector<pair<int,int> > single_coupling;
	for(size_t h=0;h<coupling_in.size();h++){
		single_coupling.clear();
		for(size_t i=0;i<coupling_in[h].size();i++){
	    	single_coupling.push_back(make_pair(coupling_in[h][i].first,perm_in[coupling_in[h][i].second-1]));
		}
		coupling.push_back(single_coupling);
	}
	for(size_t i=0;i<perm_in.size();i++) perm.push_back(i+1);


    // convert plabel vector to canonical order using 
    // the sorting functions for the virtual part
    sign=1; short conjQ=0;
    /* labels used for tree ampls
     * does not use charge conj 
    */
  
    canonical_pro_tree(labels,sign,conjQ);

    vector<int> ind;
	for(int i=0;i<labels.size();i++){
		ind.push_back(perm[labels[i].ind()-1]);
	}

	cross_term_entry* cte=new cross_term_entry(labels,perm,coupling);

    std::map<cross_term_entry*,int,CrossTermLessThan>::const_iterator it = m_cross_terms.find(cte);       
	if(it==m_cross_terms.end()){
		if(coupling.size()==0){
            		vector<kinematic_function*> local_kin_fn;
            		int res=m_ME->add(process(labels),ind,conjQ,local_kin_fn);
			m_cross_terms[cte]=res;
            		return res;
		}
		else{
	        vector<kinematic_function*> local_kin_fn_test;
	        local_kin_fn_test=(cte->m_coupling_function());
			for(size_t i=0;i<local_kin_fn_test.size();i++){
				if(local_kin_fn_test[i]==0){
                    delete cte;
				    return (-1);
			    }
			}
			{
                		// This allows to resort the assembly and improve the 
	                // efficiency of the born ME2.
	                vector<kinematic_function*> local_kin_fn;
	                coupling_process coupling_pro(cte->m_coupling,cte->m_plabel);
	                std::map<coupling_process,vector<kinematic_function*> >::const_iterator it_kin_fn = m_kinematic_functions.find(coupling_pro);       
		        	if(it_kin_fn==m_kinematic_functions.end()){
	                    		local_kin_fn=cte->m_coupling_function();
	                    		m_kinematic_functions[coupling_pro]=local_kin_fn;
	                	}
	                	else local_kin_fn=(*it_kin_fn).second;
               
                	int res=m_ME->add(process(labels),ind,conjQ,local_kin_fn);
                	m_cross_terms[cte]=res;
				return res;
			}
		}
	}
	else{
		delete cte;
		return (*it).second;
	}	
}

//reads processes in partial amplitude format, where color traces are indicated by commas "," and leptons are separated by semi-colons ";".
void read_PA_processes(const std::string& input,vector<vector<plabel> >& labels, vector<vector<std::string > >& color_info){
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
		vector<std::string > color_info_loc;
		color_info.push_back(std::vector<std::string >());

		std::string corner_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);

		std::stringstream ss_corner(corner_str);

		bool pos=true;
		while (ss_corner.good()){
			std::string partic;
			ss_corner >> partic;
			if ( !partic.empty() ){
				if(partic[0]==','||partic[0]==';'){
					pos=true;
                    color_info.back().push_back(partic);
				}
				else{
					if(pos) color_info.back().push_back("_");
					pos=true;
					labels.back().push_back(plabel_from_string(partic));
				}
			}
		}

		bracket_open = bracket_end;

	}
}

//check cross term entry printout
int construction_cache::new_loop_cross_term_entry(const std::vector<plabel>& labels,
	const std::vector<int>& perm,
	const std::vector<std::vector<pair<int,int> > >& coupling,
	const vector<std::string>& color_info){

	cross_term_entry* cte=new cross_term_entry(labels,perm,coupling,color_info);
//	vector<int> ind;
	vector<plabel> labels_new;
	for(int i=0;i<labels.size();i++){
		//_PRINT(perm[i]);
		labels_new.push_back(plabel(labels[i],perm[labels[i].ind()-1]));
	}
	//BH_DEBUG_MESSAGE2("  -- cross term entry (loop): ",labels_new);	
	BH_DEBUG_MESSAGE2("  -- cross term entry (loop): ",labels_new);	
	
	//if coupling vanishes, do not continue with the partial amplitude	
	if(coupling.size()!=0){
	        vector<kinematic_function*> local_kin_fn_test;
	        local_kin_fn_test=(cte->m_coupling_function());
			for(size_t i=0;i<local_kin_fn_test.size();i++){
				if(local_kin_fn_test[i]==0){
                    BH_DEBUG_MESSAGE("coupling vanishes, do not continue with this partial");
                    delete cte;
				    return (-1);
			    }
			}
	}

    std::map<cross_term_entry*,int,CrossTermLessThan>::const_iterator it = m_cross_terms_loop.find(cte);       
	if(it==m_cross_terms_loop.end()){
		if(coupling.size()==0){
			//pure QCD amplitude
			std::vector<kinematic_function*> local_kin_fn;
            //int res=PA_from_file(m_ME,this->m_PA_filename,labels_new,coupling,color_info,0);
	    string part_type("loop");
            int res=PA_from_file(m_ME,this->m_PA_filename,part_type,labels_new,coupling,color_info,local_kin_fn);
			m_cross_terms_loop[cte]=res;
			return res;
		}
		else{
			//QCD amplitude including Z/W/gam emissions
			if((cte->m_coupling_function()).size()==0){//vanishing coupling
				delete cte;
				return (-1);
			}
			else{
			//non-vanishing coupling
            //resort the assembly and improve the efficiency of the born ME2.
				std::vector<kinematic_function*> local_kin_fn;
                coupling_process coupling_pro(cte->m_coupling,cte->m_plabel);
                std::map<coupling_process,std::vector<kinematic_function*> >::const_iterator it_kin_fn = m_kinematic_functions.find(coupling_pro);       
	            if(it_kin_fn==m_kinematic_functions.end()){
                    local_kin_fn=cte->m_coupling_function();
                    m_kinematic_functions[coupling_pro]=local_kin_fn;
                }
                else local_kin_fn=(*it_kin_fn).second;
            
	    			string part_type("loop");
				int res=PA_from_file(m_ME,this->m_PA_filename,part_type,labels_new,coupling,color_info,local_kin_fn);
			    m_cross_terms_loop[cte]=res;
				return res;
			}
		}
	}
	else{
		//amplitude found in cache
		delete cte;
		return (*it).second;
	}	
}

construction_cache::~construction_cache(){
	std::map<cross_term_entry*,int,CrossTermLessThan>::iterator it;
	it=m_cross_terms.begin();
	while(it!=m_cross_terms.end()){ delete it->first; ++it; }
	it=m_cross_terms_loop.begin();
	while(it!=m_cross_terms_loop.end()){ delete it->first; ++it; }

}

std::string GetAssemblyDataDirectory(){
	if (settings::general::s_assembly_data_path != string("not set")){
		BH_DEBUG_MESSAGE3("Using ",settings::general::s_assembly_data_path," as the assembly data path.");
		return settings::general::s_assembly_data_path;
	} else {
		struct stat st;
		string datapath=string(BH_INSTALL_PATH)+"/share/blackhat/datafiles/assembly/";
		if( stat(datapath.c_str(),&st) == 0){
			BH_DEBUG_MESSAGE3("Using ",datapath," as the assembly data path.");
			return datapath;
		} else {
			std::string datafiles="/datafiles/assembly/";
			BH_DEBUG_MESSAGE3("Using ",BH_SOURCE_PATH+datafiles," as the assembly data path.");
			return BH_SOURCE_PATH+datafiles;
		}
	}
}




std::string ME2_file_name(const std::vector<pair<int,int> >& particle_labels){

	
	//decide on filename
	//if file exists parents will be appended
	std::string sfilename="ME2_";
	std::stringstream fileext;
	size_t n_quarks(0);
	size_t n_gluons(0);
	size_t n_photons(0);
	size_t n_leptons(0);

	int n(particle_labels.size());

	for(int i=0;i<n;i++){
		int pdb_label=particle_labels[i].second;
		if(pdb_label==21){n_gluons++;}
		else if(0<pdb_label&&pdb_label<7){n_quarks++;}
		else if(0>pdb_label&&pdb_label>-7){n_quarks++;}
		else if(10<pdb_label&&pdb_label<17){n_leptons++;}
		else if(-10>pdb_label&&pdb_label>-17){n_leptons++;}
		else if(pdb_label==22){n_photons++;}
		else {BH_DEBUG_MESSAGE("Particle not implemented in ME2_from_file:" );
			BH_DEBUG_MESSAGE(pdb_label);
		};
	};

	if(n_gluons>0){fileext<<n_gluons<<"g";};
	if(n_quarks>0){fileext<<n_quarks<<"q";};
	if(n_photons>0){fileext<<n_photons<<"y";};
	if(n_leptons>0){fileext<<n_leptons<<"l";};
	fileext<<".dat";
	sfilename+=fileext.str();
	//_PRINT(sfilename.c_str());
//	const char * ME2_filename=sfilename.c_str();
//	return ME2_filename;
	return sfilename.c_str();
};


std::string PA_file_name(const std::vector<pair<int,int> >& particle_labels,bool verbose){
	
	//decide on filename
	//if file exists parents will be appended
//	const char * partials_filename;
	std::string sfilename="partials_";
	std::stringstream fileext;
	size_t n_quarks(0);
	size_t n_gluons(0);
	size_t n_photons(0);
	size_t n_leptons(0);

	int n(particle_labels.size());

	for(int i=0;i<n;i++){
		int pdb_label=particle_labels[i].second;
		if(pdb_label<0){}
		else if(pdb_label==21){n_gluons++;}
		else if(0<pdb_label&&pdb_label<7){n_quarks++;}
		else if(10<pdb_label&&pdb_label<17){n_leptons++;}
		else if(pdb_label==22){n_photons++;}
		else {BH_DEBUG_MESSAGE("Particle not implemented in PartialAmplitude_from_file:" );
			BH_DEBUG_MESSAGE(pdb_label);
		};
	};

	if(n_gluons>0){fileext<<n_gluons<<"g";};
	if(n_quarks>0){fileext<<2*n_quarks<<"q";};
	if(n_photons>0){fileext<<n_photons<<"y";};
	if(n_leptons>0){fileext<<2*n_leptons<<"l";};
	fileext<<".dat";
	sfilename+=fileext.str();


	ifstream ifile;
	string fullfilename= GetParentDataDirectory() + string("/") + sfilename;
	//BH_DEBUG_MESSAGE2("fullpath:",fullfilename);
	/*
	ifile.open(fullfilename.c_str());
	if (verbose && !ifile){
		_WARNING("trying to read partial amplitude data");
		_WARNING4("problems opening ",sfilename, " in dir: ", settings::general::s_parent_data_path);
		throw ;
	}
	ifile.close();
	*/
	//return sfilename.c_str();
	//return fullfilename.c_str();
	return fullfilename;
}


/*  int is a reference to a vector of integer that has to have the right size (the size of PRO)
 *  perm will be filled with a map relating the labels of the reference process to indices in PRO
 *  For example if
 *    ref_process = { qp(1) qbp(2) yp(4) g(3)  }
 *  and
 *    PRO = { p qp qbp yp }
 * then perm will be
 *
 *   perm = { 2 3 1 4 }
 * because the particle labelled with k=1 in ref_process is the (perm[k-1])th particle in the process
 *
 * */
bool find_PA_process_match(const vector<plabel>& ref_process,const process& PRO,const vector<std::string >& color_info_ref, const vector<std::string >& color_info ,vector<int>& perm){

	assert( perm.size() == PRO.n());
	if (PRO.n() != ref_process.size()){
		return false;
	}
//	BH_DEBUG_MESSAGE2("color info in process match:", color_info.size());
//	for(int i=0;i<color_info.size();i++) cout<<color_info[i]; cout<<endl;
	BH_DEBUG_PRINT(PRO);;
	BH_DEBUG(cout << "ref: "; copy(ref_process.begin(),ref_process.end(),ostream_iterator<plabel>(cout," ")); cout << "\n";)
	for (int pos=1;pos<=ref_process.size();pos++){
			cyclic_iterator<particle_ID,process > cit(PRO,1,pos);
			//_PRINT(equal(ref_process.begin(),ref_process.end(),cit,plabel_has_same_type_and_ap()));
			//_PRINT(color_info_match(color_info_ref,color_info,pos-1));
			if ( equal(ref_process.begin(),ref_process.end(),cit,plabel_has_same_type_and_ap() )&&color_info_match(color_info_ref,color_info,pos-1)) {
				cyclic_iterator<particle_ID,process > cit2(PRO,1,pos);
				for (int j=0;j<ref_process.size();j++){
					int pos_in_PRO = cit2.position();  // this is the current position in PRO
					int label_in_ref = ref_process[j].ind();  // this is the label of the current particle in the referencee

					perm[label_in_ref-1]=pos_in_PRO;
					++cit2;
				}
				return true;
			}
		}

	return false;
}



//inserts into PA
//CachedOLHA::partial_amplitude_cached* 
int PA_from_file(
	    Squared_ME* ME,
		const std::string& filename,
		std::string part_type,
		const std::vector<plabel> & pro,
		const std::vector<vector<pair<int,int> > >& coupling,
		const vector<std::string >& color_info,
		std::vector<kinematic_function*> prop_hel_fn){
	
	BH_DEBUG_MESSAGE("**************************************");
	BH_DEBUG_MESSAGE("PA_from_file called");
	BH_DEBUG_MESSAGE2("file: ",filename);

	
	//new assembly CachedOLHA::partial_amplitude_cached* PA=new CachedOLHA::partial_amplitude_cached();
	//CachedOLHA::partial_amplitude_cached* PA=0;
	//partial_amplitude_cached* PA=new partial_amplitude_cached();
    std::vector<std::pair<int,double> > partial_ampl;
	if(prop_hel_fn.size()!=0){
		//new assembly PA->set_prefactor(*prop_hel_fn);
	};


	ifstream ifile;
	ifile.open(filename.c_str());
	
	int line_nbr=0;
	int n=pro.size();

	while ( ifile ){
		char line[2500];
		ifile.getline(line,2500);
		line_nbr++;
		switch ( line[0] ){
			case '\0':  break;
			case '*' : {
			BH_DEBUG_MESSAGE(" process analyse");
			string line_str(line);
				vector<vector<plabel> > processes;
				vector<vector<std::string > > color_info_file;
				read_PA_processes(string(line),processes,color_info_file);
				
				vector<int> perm(pro.size());
				for (int i=0;i<processes.size();i++){
					BH_DEBUG_MESSAGE2("ref process:  ",pro);
					BH_DEBUG_MESSAGE2("read process: ",processes[i]);
					//for(int j=0;j<color_info.size();j++) cout<<color_info[j];cout<<endl;
					//for(int j=0;j<color_info.size();j++) cout<<color_info_file[i][j];cout<<endl;
					if ( find_PA_process_match(processes[i],process(pro),color_info_file[i],color_info,perm) ){
						//BH_DEBUG_MESSAGE2("ref process:  ",pro);
						//BH_DEBUG_MESSAGE2("read process: ",processes[i]);
						//cout<<"pos in ref pro: ";for(int j=0;j<perm.size();j++) cout<<perm[j];cout<<endl;
						//combine input index with permutiation form PA-file
						BH_DEBUG_MESSAGE2("PA_MATCH at line ",line_nbr);
						char line[250];
						ifile.getline(line,250);
						//check that coupling matches
						//if(part=="phZW") ...
						while (line[0] != '*' && ifile ){
							if (line[0] != '\0' && line[0]!= '#'  ){
								string line_str(line);
								int part_end=line_str.find('{',1);
								string part=NoSpaces(line_str.substr(1,part_end-2));
								BH_DEBUG_MESSAGE2("part: ",part);
								//if(part=="loop"){
								if(part==part_type){
									while (line[0] != '*' //|| 
										//string::npos==line_str.find_last_of('}')
										&& ifile ){
										if(line[0] != '\0' && line[0] != '#')
//										read_primitive_amplitude(string(line),pro,ind_perm,PA);
										read_primitive_amplitude(string(line),pro,perm/*,PA*/,partial_ampl,ME,prop_hel_fn);
										ifile.getline(line,250);
									}
									BH_DEBUG_MESSAGE("**************************************");
									ifile.close();
                                    
				                                    return ME->add_partial(partial_ampl);
								}	
							}
							ifile.getline(line,250);
						}
						BH_DEBUG_MESSAGE("WARNING: failed to find loop part in partials data file");
						ifile.close();
						return (-1); //return if a matching partial is found and used.
					}
				}
			}; break;
			default: break;
		}
	}
	BH_DEBUG_MESSAGE("**************************************");
	ifile.close();
	return 0;
}



//needed
Virtual_SME* get_ME2_from_file(/*const string& filename, const std::string& ME2_order_str*/ QCDorder lo_or_nlo, const std::vector<std::pair<int, int> >& pa_labels_old){
	BH_DEBUG_MESSAGE("***************************************");
	BH_DEBUG_MESSAGE(" start reading ME2-files");
	std::string filename=ME2_file_name(pa_labels_old);
	BH_DEBUG_MESSAGE2(" file:      ",filename);
	cout<<" process:   "; for(int i=0;i<2;i++) {
        if(abs(pa_labels_old[i].second)>20){    
            cout<<pa_labels_old[i].second<<" ";}
        else{ cout<<-pa_labels_old[i].second<<" ";}
    };
	cout<<"-> "; for(int i=2;i<pa_labels_old.size();i++) {cout<<pa_labels_old[i].second<<" ";}; cout<<endl;
       	vector<std::string> ME2_order_str;
    Squared_ME* ME=new Squared_ME_Struct(lo_or_nlo);
	switch(lo_or_nlo){
		case lo: {
			BH_DEBUG_MESSAGE(" QCD order: LO"); 
			ME2_order_str.push_back("BSM");
            //ME = new Squared_ME_Struct(lo);
		    //
        } break;
		case nlo: {
			BH_DEBUG_MESSAGE(" QCD order: NLO"); 
			ME2_order_str.push_back("BSM");
			ME2_order_str.push_back("VSM");
			//ME = new Squared_ME(nlo);
			//ME = new Squared_ME_Struct(nlo);
		} break;
	}

	//new conventions for pa_labels
	std::vector<std::pair<int, int> > pa_labels;
	for(int i=0;i<pa_labels_old.size();i++){
		pa_labels.push_back(make_pair(pa_labels_old[i].second,pa_labels_old[i].first));
	}
	construction_cache * constr_cache= new construction_cache(ME,pa_labels_old);
	Virtual_SME* VSM=new Virtual_SME() ;

	std::string fullfilename=GetAssemblyDataDirectory()+filename;
	BH_DEBUG_MESSAGE(" start reading ME2-files");
	BH_DEBUG_MESSAGE2(" abs file name:  ",fullfilename);
	
	//loop over QCD approximations: LO, NLO, ... (possibly color-approximations and the like )
	for(int j=0;j<ME2_order_str.size();j++){
		BH_DEBUG_MESSAGE(" -----");
		bool found_match(false);
	
	//open file
	ifstream ifile;
	ifile.open(fullfilename.c_str());
	if (!ifile){
		_WARNING4("Problems opening file ",filename, " with full path ", fullfilename.c_str());
		//return zero pointer in order to use automated assembly to construct assembly file
		return 0;
	}
	
	//reading file
	int line_nbr=0;
	while ( !found_match && ifile ){
		char line[2500];
		ifile.getline(line,2500);
		line_nbr++;
		switch ( line[0] ){
		case '\0' :{
//			BH_DEBUG_MESSAGE(" blank line");
			   } break;
		case '*' : {
//			BH_DEBUG_MESSAGE(" process");
			string line_str(line);
			int ME2_order_end=line_str.find('{',1);
			string ME2_order=NoSpaces(line_str.substr(1,ME2_order_end-2));
			bool born_or_virt=(ME2_order_str[j]=="BSM");
				
			vector<vector<pair<int,int> > > pa_labels_ref;
			if (ME2_order == ME2_order_str[j]){
				read_pa_labels(string(line),pa_labels_ref);
				
				vector<int> perm(pa_labels.size());
				for (int i=0;i<pa_labels_ref.size();i++){
					if ( find_pa_labels_match(pa_labels_ref[i],pa_labels,perm) ){
						BH_DEBUG_MESSAGE4("* ME2 MATCH at line: ",line_nbr," for ",ME2_order_str[j]);
						BH_DEBUG_MESSAGE2("  permutation: ",perm);
						found_match=true;
						char line[250];
						ifile.getline(line,250);
						while (line[0] != '*' && ifile ){
							if (line[0] != '\0' && line[0]!= '#'  ){
							read_cross_term(string(line),pa_labels,perm,constr_cache, ME,born_or_virt);
							}
							ifile.getline(line,250);
						}
						/*
						 * found match
						 * add Virtual_SME
						 */
						break;
					} 
				};
			};
		} break;
		default: break;
		}
	}
	ifile.close();
	//if no suitable assembly entry is found return 0, 
	//so that it gets created
		if ( !found_match && !ifile ){
			BH_DEBUG_MESSAGE("No matching entry found in ME2 file - return zero-pointer");
			_MESSAGE("No matching entry found in ME2 file - return zero-pointer");
			delete constr_cache;
			return 0;
		}
	}
	BH_DEBUG_MESSAGE("-----");
	delete constr_cache;
	//if(lo_or_nlo==lo){
    //}
    //else if(lo_or_nlo==nlo){
    //}
    
    	ME->complete_construction();        
   	// ME->complete_virt_construction();        
	VSM->add(ME);
	return VSM;
}



}
