/*
 * process_utils_massive.cpp
 *
 *  Created on: 28-Jul-2008
 *      Author: Darren
 */

#include <vector>
#include "particles.h"
#include "process.h"
#include "process_utils.h"
#include "BH_utilities.h"
#include "process_utils_massive.h"
#include "iterators.h"
#include "eval_param.h"
#include <algorithm>
using namespace std;

namespace BH {

#if 0
//merge OCT23
//From process_utils.cpp
vector<int> find_all_flavors(const process& pro,const particle& type1, const particle& type2, const particle& type3){
	size_t count=0;
  	vector<int> flavors;
  	for(int i=1;i<=pro.n();i++)
  	{
  		if( ((*pro.p(i).type()) == type1)||((*pro.p(i).type()) == type2 )||((*pro.p(i).type()) == type3 )) {
  			vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
  	    	if ( pos == flavors.end() ){
  	    		count++;
  	    		flavors.push_back(pro.p(i).flavor());
  	    	}
		}
  	}
  	return flavors;
}
#endif
process replace_massive_gluino_with_massive_quark(const process& pro){
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(gluino_massive)){
			pp.push_back(particle_ID(quark_massive,pro.p(i).helicity(),pro.p(i).flavor(),pro.p(i).is_anti_particle()));
		} else {
			pp.push_back(pro.p(i));
		}
	}
	return process(pp);
}
#if 0
//merge OCT23
vector<int> find_all_flavors(const process& pro,const particle& type1, const particle& type2){
	size_t count=0;
  	vector<int> flavors;
  	for(int i=1;i<=pro.n();i++)
  	{
  		if( ((*pro.p(i).type()) == type1)||((*pro.p(i).type()) == type2 )) {
  			vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
  	    	if ( pos == flavors.end() ){
  	    		count++;
  	    		flavors.push_back(pro.p(i).flavor());
  	    	}
		}
  	}
  	return flavors;
}
#endif

process arrange_flavors_qq_massive_LT(const process& pro,vector<particle_ID>& propagators){
	propagators.push_back(gsc);

	process::const_iterator pos_q=find_if(pro.begin(),pro.end(),is_of_type_and_pap(quark,true));
	propagators.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor()+100,true));
	propagators.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor()+100,true));

	return pro;
}

process arrange_flavors_qq_massive_RT(const process& pro,vector<particle_ID>& propagators){
	propagators.push_back(gsc);

	process::const_iterator pos_q=find_if(pro.begin(),pro.end(),is_of_type_and_pap(quark,false));
	propagators.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor()+100,false));
	propagators.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor()+100,false));

	return pro;
}

process arrange_flavors_massive(const process& pro,size_t qb_pos,size_t q_pos,short increment,vector<particle_ID>& propagators){
	vector<particle_ID> temp;
	for (int i=1;i<=pro.n();i++){
		temp.push_back(pro.p(i));
	}

	// Next we compute the separation of the quarks on the quark line, i.e. how many gluonic three point
	//  vertices we will "pass" going around the loop on the quark line
	size_t distance;
	switch (increment){
		case 1 :{
			distance=(pro.n()+q_pos-qb_pos)%pro.n(); break;
		}
		case -1:	{
			distance=(pro.n()+qb_pos-q_pos)%pro.n();break;
		}
	}

	// For each three point quark-quark-gluon vertex we change the "flavour" of the gluon in the external leg and
	//  insert quarks into the propagator carrying this flavour.

	//In contrast to the massless case we have a "flavour" changing scalar at the three point
	//  vertex of the incoming quark and so we give it a new flavour rather than starting at zero.
	for (size_t j=1;j<distance;j++){
		temp[(qb_pos+pro.n()+j*increment-1)%pro.n()]=particle_ID(gluon,pro.p((qb_pos+pro.n()+j*increment-1)%pro.n()+1).helicity(),-j,false);
		propagators.push_back(particle_ID(quark_massive,1,j,true));
		propagators.push_back(particle_ID(quark_massive,-1,j,true));
		propagators.push_back(particle_ID(quark_massive,1,j,false));
		propagators.push_back(particle_ID(quark_massive,-1,j,false));
	}
	// For the final leg we simply add in the extra flavour propagator
	propagators.push_back(particle_ID(quark_massive,1,distance,true));
	propagators.push_back(particle_ID(quark_massive,-1,distance,true));
	propagators.push_back(particle_ID(quark_massive,1,distance,false));
	propagators.push_back(particle_ID(quark_massive,-1,distance,false));

	// We finally fix up the flavours of the two external quarks to match the change of
	//  "flavours" for the gluons.
	temp[qb_pos-1]=particle_ID(quark,pro.p(qb_pos).helicity(),1,true);
	temp[q_pos-1]=particle_ID(quark,pro.p(q_pos).helicity(),distance,false);

	return process(temp);
}

process arrange_flavors_2q2e_massive(const process& pro,vector<particle_ID>& propagators){
	vector<particle_ID> temp;
	size_t pos_qb=0;
	size_t pos_q=0;

	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(quark)){
			if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
		}
		temp.push_back(pro.p(i));
	}

	particle_type_match is_gluon(gluon);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);
	particle_type_match is_lepton(lepton);


	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"qb");
	rule R3(&is_lepton,"e");

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);


	std::string type=string_gen(pro,rules);
	type+=type;
	std::string::size_type loc;

	loc = type.find("qbeeq",0);
	if (loc != std::string::npos){
		propagators.push_back(particle_ID(quark_massive,+1,-1,true));
		propagators.push_back(particle_ID(quark_massive,-1,-1,true));
		propagators.push_back(particle_ID(quark_massive,+1,-2,true));
		propagators.push_back(particle_ID(quark_massive,-1,-2,true));
		propagators.push_back(gsc);
		temp[pos_qb-1]=particle_ID(quark,pro.p(pos_qb).helicity(),-1,true);
		temp[pos_q -1]=particle_ID(quark,pro.p(pos_q ).helicity(),-2,false);
	}

	loc = type.find("qeeqb",0);
	if (loc != std::string::npos){
		propagators.push_back(particle_ID(quark_massive,+1,-1,false));
		propagators.push_back(particle_ID(quark_massive,-1,-1,false));
		propagators.push_back(particle_ID(quark_massive,+1,-2,false));
		propagators.push_back(particle_ID(quark_massive,-1,-2,false));
		propagators.push_back(gsc);
		temp[pos_qb-1]=particle_ID(quark,pro.p(pos_qb).helicity(),-1,true);
		temp[pos_q -1]=particle_ID(quark,pro.p(pos_q ).helicity(),-2,false);
	}

	return process(temp);
}

process arrange_flavors_4q2e_massive(const process& pro,vector<particle_ID>& propagators){
//	_PRINT(pro);
	vector<particle_ID> temp(pro.begin(),pro.end());
	size_t pos_qb=0;
	size_t pos_q=0;
	size_t pos_qb2=0;
	size_t pos_q2=0;

	bool first_quark_is_antiquark=false;
	bool first_quark2_is_antiquark=false;
	int first_quark_helicity;
	int second_flavor;

	process::const_iterator lepton_pos=find_if(pro.begin(),pro.end(),is_of_type(lepton));

	cyclic_iterator<particle_ID,process> cy(pro,1,lepton_pos);
	cyclic_iterator<particle_ID,process> cy_end(pro);

	cyclic_iterator<particle_ID,process> first_q=find_if(cy,cy_end,is_of_type(quark));

	first_quark_helicity=(*first_q).helicity();
	int first_flavor=(*first_q).flavor();
	if ((*first_q).is_anti_particle()){
		pos_qb=first_q.position();
		first_quark_is_antiquark=true;
	} else {
		pos_q=first_q.position();
		first_quark_is_antiquark=false;
	}


	cyclic_iterator<particle_ID,process> second_q=find_if(++first_q,cy_end,is_of_type(quark));
	second_flavor=(*second_q).flavor();
	bool sub_leading_configuration = ( first_flavor == second_flavor );
	if ( !sub_leading_configuration){
		if ((*second_q).is_anti_particle()){
			pos_qb2=second_q.position();
			first_quark2_is_antiquark=true;
		} else {
			pos_q2=second_q.position();
			first_quark2_is_antiquark=false;
		}
	} else {
		if ((*second_q).is_anti_particle()){
			pos_qb=second_q.position();
		} else {
			pos_q=second_q.position();
		}
	}
	cyclic_iterator<particle_ID,process> third_q=find_if(++second_q,cy_end,is_of_type(quark));

	second_flavor=(*third_q).flavor();
	if (!sub_leading_configuration){
		if ((*third_q).is_anti_particle()){
			pos_qb2=third_q.position();
		} else {
			pos_q2=third_q.position();
		}
	} else {
		if ((*third_q).is_anti_particle()){
			pos_qb2=third_q.position();
			first_quark2_is_antiquark=true;
		} else {
			pos_q2=third_q.position();
			first_quark2_is_antiquark=false;
		}

	}
	cyclic_iterator<particle_ID,process> fourth_q=find_if(++third_q,cy_end,is_of_type(quark));

	if ( !sub_leading_configuration){
		if ((*fourth_q).is_anti_particle()){
			pos_qb=fourth_q.position();
		} else {
			pos_q=fourth_q.position();
		}
	}
	else {
		if ((*fourth_q).is_anti_particle()){
			pos_qb2=fourth_q.position();
		} else {
			pos_q2=fourth_q.position();
		}
	}

#if _VERBOSE
	_PRINT(pos_q);
	_PRINT(pos_qb);
	_PRINT(pos_q2);
	_PRINT(pos_qb2);
#endif



	if (! first_quark_is_antiquark){
		propagators.push_back(particle_ID(quark_massive,-first_quark_helicity,-1,true));
		propagators.push_back(particle_ID(quark_massive,-first_quark_helicity,-2,true));
		propagators.push_back(gsc);
		temp[pos_qb-1]=particle_ID(quark,pro.p(pos_qb).helicity(),-1,true);
		temp[pos_q -1]=particle_ID(quark,pro.p(pos_q ).helicity(),-2,false);
	}

	if (first_quark_is_antiquark){
		propagators.push_back(particle_ID(quark_massive,-first_quark_helicity,-1,false));
		propagators.push_back(particle_ID(quark_massive,-first_quark_helicity,-2,false));
		propagators.push_back(gsc);
		temp[pos_qb-1]=particle_ID(quark,pro.p(pos_qb).helicity(),-1,true);
		temp[pos_q -1]=particle_ID(quark,pro.p(pos_q ).helicity(),-2,false);
	}

	if (!sub_leading_configuration){
		if (first_quark2_is_antiquark){
			propagators.push_back(particle_ID(quark_massive,+1,second_flavor,true));
			propagators.push_back(particle_ID(quark_massive,-1,second_flavor,true));
		} else {
			propagators.push_back(particle_ID(quark_massive,+1,second_flavor,false));
			propagators.push_back(particle_ID(quark_massive,-1,second_flavor,false));
		}
	}
	else {
		if (first_quark2_is_antiquark){
			propagators.push_back(particle_ID(quark_massive,+1,second_flavor+100,false));
			propagators.push_back(particle_ID(quark_massive,-1,second_flavor+100,false));
		} else {
			propagators.push_back(particle_ID(quark_massive,+1,second_flavor+100,true));
			propagators.push_back(particle_ID(quark_massive,-1,second_flavor+100,true));
		}
	}

//	copy(propagators.begin(),propagators.end(),ostream_iterator<particle_ID>(cout," "));
//	cout << endl;
	return process(temp);
}

process arrange_flavors_2q2G2e_massive(const process& pro,vector<particle_ID>& propagators){
//	_PRINT(pro);
	vector<particle_ID> temp;
		size_t pos_qb=0;
		size_t pos_q=0;

		for (int i=1;i<=pro.n();i++){
			if (pro.p(i).is_a(quark)){
				if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
			}
			temp.push_back(pro.p(i));
		}

		particle_type_match is_a_gluon(gluon);
		particle_type_match is_a_gluino(gluino);
		particle_type_match is_a_quark(quark);
		particle_type_match is_a_lepton(lepton);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);




//	_PRINT(type);

	std::string::size_type loc,loc1,loc2,loc3;

	propagators.push_back(gsc);

	process::const_iterator lepton_pos=find_if(pro.begin(),pro.end(),is_of_type(lepton));

	cyclic_iterator<particle_ID,process> cy(pro,1,lepton_pos);
	cyclic_iterator<particle_ID,process> cy_end(pro);

	cyclic_iterator<particle_ID,process> first_G=find_if(cy,cy_end,is_of_either_type(gluino,gluino_massive));

	bool sub_leading_configuration;

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"qb");
	rule R3(&is_a_lepton,"e");
	rule R4(&is_a_gluino,"G");
	rule R5(&is_a_quark,"q");

	vector<rule> rules_L_or_SL;
	rules_L_or_SL.push_back(R3);
	rules_L_or_SL.push_back(R4);
	rules_L_or_SL.push_back(R5);

	std::string type_L_or_SL=string_gen(pro,rules_L_or_SL);
	type_L_or_SL+=type_L_or_SL;
	loc1 = type_L_or_SL.find("qqGGee",0);
	loc2 = type_L_or_SL.find("GGqqee",0);
	loc3 = type_L_or_SL.find("qGGqee",0);



	if (loc1 != std::string::npos || loc2 != std::string::npos){
		sub_leading_configuration=true;
	} else 	if (loc3 != std::string::npos ){
		sub_leading_configuration=false;
	} else {
		_WARNING("no type found in arrange_flavors_2q2G2e_massive");
	}

	if ( ! sub_leading_configuration ){
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),true));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),false));
		}
	}
	else {
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,false));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,true));
		}
	}

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string type=string_gen(pro,rules);
	type+=type;
	loc = type.find("qbeeq",0);
	if (loc != std::string::npos){
		//Include both helicites here unlike in the massless case as we can have helicity flip vertices
		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
	}

	loc = type.find("qeeqb",0);
	if (loc != std::string::npos){
		//Include both helicites here unlike in the massless case as we can have helicity flip vertices
		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_qb).helicity(),pro.p(pos_qb).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,pro.p(pos_qb).helicity(),pro.p(pos_qb).flavor()+100,false));
	}
	return process(temp);
}

process arrange_flavors_2q2G2e_q_massive(const process& pro,vector<particle_ID>& propagators){
//	_PRINT(pro);
	vector<particle_ID> temp;
		size_t pos_qb=0;
		size_t pos_q=0;

		for (int i=1;i<=pro.n();i++){
			if (pro.p(i).is_a(quark)){
				if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
			}
			temp.push_back(pro.p(i));
		}

		particle_type_match is_a_gluon(gluon);
		particle_type_match is_a_gluino(gluino);
		particle_type_match is_a_quark(quark);
		particle_type_match is_a_lepton(lepton);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);




//	_PRINT(type);

	std::string::size_type loc,loc1,loc2,loc3;

	propagators.push_back(gsc);

	process::const_iterator lepton_pos=find_if(pro.begin(),pro.end(),is_of_type(lepton));

	cyclic_iterator<particle_ID,process> cy(pro,1,lepton_pos);
	cyclic_iterator<particle_ID,process> cy_end(pro);

	cyclic_iterator<particle_ID,process> first_G=find_if(cy,cy_end,is_of_either_type(gluino,gluino_massive));

	bool sub_leading_configuration;

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"qb");
	rule R3(&is_a_lepton,"e");
	rule R4(&is_a_gluino,"G");
	rule R5(&is_a_quark,"q");

	vector<rule> rules_L_or_SL;
	rules_L_or_SL.push_back(R3);
	rules_L_or_SL.push_back(R4);
	rules_L_or_SL.push_back(R5);

	std::string type_L_or_SL=string_gen(pro,rules_L_or_SL);
	type_L_or_SL+=type_L_or_SL;
	loc1 = type_L_or_SL.find("qqGGee",0);
	loc2 = type_L_or_SL.find("GGqqee",0);
	loc3 = type_L_or_SL.find("qGGqee",0);



	if (loc1 != std::string::npos || loc2 != std::string::npos){
		sub_leading_configuration=true;
	} else 	if (loc3 != std::string::npos ){
		sub_leading_configuration=false;
	} else {
		_WARNING("no type found in arrange_flavors_2q2G2e");
	}

	if ( ! sub_leading_configuration ){
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),true));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),false));
		}
	}

/*
 * want no gluinos in the loop!!
 * else {
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino,+1,(*first_G).flavor()+100,false));
			propagators.push_back(particle_ID(gluino,-1,(*first_G).flavor()+100,false));
//						propagators.push_back(Gbp);
//						propagators.push_back(Gbm);
		} else {
			propagators.push_back(particle_ID(gluino,+1,(*first_G).flavor()+100,true));
			propagators.push_back(particle_ID(gluino,-1,(*first_G).flavor()+100,true));
//			propagators.push_back(Gp);
//			propagators.push_back(Gm);
		}
	}
*/

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string type=string_gen(pro,rules);
	type+=type;
	loc = type.find("qbeeq",0);
	if (loc != std::string::npos){
		//Include both helicites here unlike in the massless case as we can have helicity flip vertices
		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
	}

	loc = type.find("qeeqb",0);
	if (loc != std::string::npos){
		//Include both helicites here unlike in the massless case as we can have helicity flip vertices
		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_qb).helicity(),pro.p(pos_qb).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,pro.p(pos_qb).helicity(),pro.p(pos_qb).flavor()+100,false));
	}
	return process(temp);
}


process arrange_flavors_2q2G1y_q_massive(const process& pro,vector<particle_ID>& propagators){
//	_PRINT(pro);
	vector<particle_ID> temp;
		size_t pos_qb=0;
		size_t pos_q=0;

		for (int i=1;i<=pro.n();i++){
			if (pro.p(i).is_a(quark)){
				if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
			}
			temp.push_back(pro.p(i));
		}

		particle_type_match is_a_gluon(gluon);
		particle_type_match is_a_gluino(gluino);
		particle_type_match is_a_quark(quark);
		particle_type_match is_a_photon(photon);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);




//	_PRINT(type);

	std::string::size_type loc,loc1,loc2,loc3;

	propagators.push_back(gsc);

	process::const_iterator photon_pos=find_if(pro.begin(),pro.end(),is_of_type(photon));

	cyclic_iterator<particle_ID,process> cy(pro,1,photon_pos);
	cyclic_iterator<particle_ID,process> cy_end(pro);

	cyclic_iterator<particle_ID,process> first_G=find_if(cy,cy_end,is_of_either_type(gluino,gluino_massive));

	bool sub_leading_configuration;

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"qb");
	rule R3(&is_a_photon,"y");
	rule R4(&is_a_gluino,"G");
	rule R5(&is_a_quark,"q");

	vector<rule> rules_L_or_SL;
	rules_L_or_SL.push_back(R3);
	rules_L_or_SL.push_back(R4);
	rules_L_or_SL.push_back(R5);

	std::string type_L_or_SL=string_gen(pro,rules_L_or_SL);
	type_L_or_SL+=type_L_or_SL;
	loc1 = type_L_or_SL.find("qqGGy",0);
	loc2 = type_L_or_SL.find("GGqqy",0);
	loc3 = type_L_or_SL.find("qGGqy",0);



	if (loc1 != std::string::npos || loc2 != std::string::npos){
		sub_leading_configuration=true;
	} else 	if (loc3 != std::string::npos ){
		sub_leading_configuration=false;
	} else {
		_WARNING("no type found in arrange_flavors_2q2G1y");
	}

	if ( ! sub_leading_configuration ){
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),true));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),false));
		}
	}

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string type=string_gen(pro,rules);
	type+=type;
	loc = type.find("qbyq",0);
	if (loc != std::string::npos){
		//Include both helicites here unlike in the massless case as we can have helicity flip vertices
		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
	}

	loc = type.find("qyqb",0);
	if (loc != std::string::npos){
		//Include both helicites here unlike in the massless case as we can have helicity flip vertices
		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_qb).helicity(),pro.p(pos_qb).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,pro.p(pos_qb).helicity(),pro.p(pos_qb).flavor()+100,false));
	}
	return process(temp);
}


process arrange_flavors_2q2G2e_G_massive(const process& pro,vector<particle_ID>& propagators){
//	_PRINT(pro);
	vector<particle_ID> temp;
		size_t pos_qb=0;
		size_t pos_q=0;

		for (int i=1;i<=pro.n();i++){
			if (pro.p(i).is_a(quark)){
				if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
			}
			temp.push_back(pro.p(i));
		}


		particle_type_match is_a_gluon(gluon);
		particle_type_match is_a_gluino(gluino);
		particle_type_match is_a_quark(quark);
		particle_type_match is_a_lepton(lepton);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);




//	_PRINT(type);

	std::string::size_type loc,loc1,loc2,loc3;

	propagators.push_back(gsc);

	process::const_iterator lepton_pos=find_if(pro.begin(),pro.end(),is_of_type(lepton));

	cyclic_iterator<particle_ID,process> cy(pro,1,lepton_pos);
	cyclic_iterator<particle_ID,process> cy_end(pro);

	cyclic_iterator<particle_ID,process> first_G=find_if(cy,cy_end,is_of_type(gluino));

	bool sub_leading_configuration;

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"qb");
	rule R3(&is_a_lepton,"e");
	rule R4(&is_a_gluino,"G");
	rule R5(&is_a_quark,"q");

	vector<rule> rules_L_or_SL;
	rules_L_or_SL.push_back(R3);
	rules_L_or_SL.push_back(R4);
	rules_L_or_SL.push_back(R5);

	std::string type_L_or_SL=string_gen(pro,rules_L_or_SL);
	type_L_or_SL+=type_L_or_SL;
	loc1 = type_L_or_SL.find("qqGGee",0);
	loc2 = type_L_or_SL.find("GGqqee",0);
	loc3 = type_L_or_SL.find("qGGqee",0);



	if (loc1 != std::string::npos || loc2 != std::string::npos){
		sub_leading_configuration=true;
	} else 	if (loc3 != std::string::npos ){
		sub_leading_configuration=false;
	} else {
		_WARNING("no type found in arrange_flavors_2q2G2e");
	}

	if ( ! sub_leading_configuration ){
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),true));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),false));
		}
	}
	else {
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,false));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,true));
		}
	}
/* want no quarks in the loop
 *
 *
	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string type=string_gen(pro,rules);
	type+=type;
	loc = type.find("qbeeq",0);
	if (loc != std::string::npos){
		propagators.push_back(particle_ID(quark,-pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
//		propagators.push_back(particle_ID(quark,-pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
//		temp[pos_qb-1]=particle_ID(quark,pro.p(pos_qb).helicity(),-1,true);
//		temp[pos_q -1]=particle_ID(quark,pro.p(pos_q ).helicity(),-2,false);
	}

	loc = type.find("qeeqb",0);
	if (loc != std::string::npos){
		propagators.push_back(particle_ID(quark,-pro.p(pos_qb).helicity(),pro.p(pos_qb).flavor()+100,false));
//		propagators.push_back(particle_ID(quark,-pro.p(pos_qb).helicity(),pro.p(pos_qb).flavor()+100,false));
//		temp[pos_qb-1]=particle_ID(quark,pro.p(pos_qb).helicity(),-1,true);
//		temp[pos_q -1]=particle_ID(quark,pro.p(pos_q ).helicity(),-2,false);
	}
*/
	return process(temp);
}


process arrange_flavors_2q2G1y_G_massive(const process& pro,vector<particle_ID>& propagators){
//	_PRINT(pro);
	vector<particle_ID> temp;
		size_t pos_qb=0;
		size_t pos_q=0;

		for (int i=1;i<=pro.n();i++){
			if (pro.p(i).is_a(quark)){
				if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
			}
			temp.push_back(pro.p(i));
		}


		particle_type_match is_a_gluon(gluon);
		particle_type_match is_a_gluino(gluino);
		particle_type_match is_a_quark(quark);
		particle_type_match is_a_photon(photon);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);




//	_PRINT(type);

	std::string::size_type loc,loc1,loc2,loc3;

	propagators.push_back(gsc);

	process::const_iterator photon_pos=find_if(pro.begin(),pro.end(),is_of_type(photon));

	cyclic_iterator<particle_ID,process> cy(pro,1,photon_pos);
	cyclic_iterator<particle_ID,process> cy_end(pro);

	cyclic_iterator<particle_ID,process> first_G=find_if(cy,cy_end,is_of_type(gluino));

	bool sub_leading_configuration;

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"qb");
	rule R3(&is_a_photon,"y");
	rule R4(&is_a_gluino,"G");
	rule R5(&is_a_quark,"q");

	vector<rule> rules_L_or_SL;
	rules_L_or_SL.push_back(R3);
	rules_L_or_SL.push_back(R4);
	rules_L_or_SL.push_back(R5);

	std::string type_L_or_SL=string_gen(pro,rules_L_or_SL);
	type_L_or_SL+=type_L_or_SL;
	loc1 = type_L_or_SL.find("qqGGy",0);
	loc2 = type_L_or_SL.find("GGqqy",0);
	loc3 = type_L_or_SL.find("qGGqy",0);



	if (loc1 != std::string::npos || loc2 != std::string::npos){
		sub_leading_configuration=true;
	} else 	if (loc3 != std::string::npos ){
		sub_leading_configuration=false;
	} else {
		_WARNING("no type found in arrange_flavors_2q2G1y");
	}

	if ( ! sub_leading_configuration ){
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),true));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),false));
		}
	}
	else {
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,false));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,true));
		}
	}
	return process(temp);
}

process arrange_flavors_2q4G2e_massive(const process& pro,vector<particle_ID>& propagators){
	vector<particle_ID> temp;
		for (int i=1;i<=pro.n();i++){
			temp.push_back(pro.p(i));
		}

	propagators.push_back(gsc);

	process::const_iterator quark_pos=find_if(pro.begin(),pro.end(),is_of_type(quark));

	cyclic_iterator<particle_ID,process > c1(pro,2,quark_pos);

	// Any of the quarks is our starting point
	// Next we find a gluino1
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > first_G1(c1);
	// Find the second gluino1
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > second_G1(c1);
	// then the gluino2
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > first_G2(c1);
	// then the gluino2
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > second_G2(c1);
	// then the first quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > first_q(c1);
	// then the second quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);
	
	propagators.push_back(particle_ID(gluino_massive,(*first_G1).helicity(),(*first_G1).flavor(),((*first_G1).is_anti_particle())));
	propagators.push_back(particle_ID(gluino_massive,-(*first_G1).helicity(),(*first_G1).flavor(),((*first_G1).is_anti_particle())));
	propagators.push_back(particle_ID(gluino_massive,(*first_G2).helicity(),(*first_G2).flavor(),((*first_G2).is_anti_particle())));
	propagators.push_back(particle_ID(gluino_massive,-(*first_G2).helicity(),(*first_G2).flavor(),((*first_G2).is_anti_particle())));
	propagators.push_back(particle_ID(quark_massive,(*first_q).helicity(),(*first_q).flavor()+100,((*first_q).is_anti_particle())));
	propagators.push_back(particle_ID(quark_massive,-(*first_q).helicity(),(*first_q).flavor()+100,((*first_q).is_anti_particle())));
		
//	_PRINT(pro);
//	_PRINT(propagators);

	return process(temp);
}





//process arrange_flavors_2q2G1y_massive(const process& pro,vector<particle_ID>& propagators){
////	_PRINT(pro);
//	vector<particle_ID> temp;
//		size_t pos_qb=0;
//		size_t pos_q=0;
//
//		for (int i=1;i<=pro.n();i++){
//			if (pro.p(i).is_a(quark)){
//				if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
//			}
//			temp.push_back(pro.p(i));
//		}
//
//		particle_type_match is_a_gluon(gluon);
//		particle_type_match is_a_gluino(gluino);
//		particle_type_match is_a_quark(quark);
//		particle_type_match is_a_photon(photon);
//		type_and_pap_match is_quark(quark,false);
//		type_and_pap_match is_antiquark(quark,true);
//
//
//
//
////	_PRINT(type);
//
//	std::string::size_type loc,loc1,loc2,loc3;
//
//	propagators.push_back(gsc);
//
//	process::const_iterator photon_pos=find_if(pro.begin(),pro.end(),is_of_type(photon));
//
//	cyclic_iterator<particle_ID,process> cy(pro,1,photon_pos);
//	cyclic_iterator<particle_ID,process> cy_end(pro);
//
//	cyclic_iterator<particle_ID,process> first_G=find_if(cy,cy_end,is_of_type(gluino));
//
//	bool sub_leading_configuration;
//
//	rule R1(&is_quark,"q");
//	rule R2(&is_antiquark,"qb");
//	rule R3(&is_a_photon,"y");
//	rule R4(&is_a_gluino,"G");
//	rule R5(&is_a_quark,"q");
//
//	vector<rule> rules_L_or_SL;
//	rules_L_or_SL.push_back(R3);
//	rules_L_or_SL.push_back(R4);
//	rules_L_or_SL.push_back(R5);
//
//	std::string type_L_or_SL=string_gen(pro,rules_L_or_SL);
//	type_L_or_SL+=type_L_or_SL;
//	loc1 = type_L_or_SL.find("qqGGy",0);
//	loc2 = type_L_or_SL.find("GGqqy",0);
//	loc3 = type_L_or_SL.find("qGGqy",0);
//
//
//
//	if (loc1 != std::string::npos || loc2 != std::string::npos){
//		sub_leading_configuration=true;
//	} else 	if (loc3 != std::string::npos ){
//		sub_leading_configuration=false;
//	} else {
//		_WARNING("no type found in arrange_flavors_2q2G1y_massive");
//	}
//
//	if ( ! sub_leading_configuration ){
//		if ((*first_G).is_anti_particle()){
//			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),true));
//			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),true));
////			propagators.push_back(Gbp);
////			propagators.push_back(Gbm);
//		} else {
//			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor(),false));
//			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor(),false));
////			propagators.push_back(Gp);
////			propagators.push_back(Gm);
//		}
//	}
//	else {
//		if ((*first_G).is_anti_particle()){
//			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,false));
//			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,false));
//		} else {
//			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,true));
//			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,true));
//		}
//	}
//
//	vector<rule> rules;
//	rules.push_back(R1);
//	rules.push_back(R2);
//	rules.push_back(R3);
//
//	std::string type=string_gen(pro,rules);
//	type+=type;
//	loc = type.find("qbyq",0);
//	if (loc != std::string::npos){
//		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
//		propagators.push_back(particle_ID(quark_massive,pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
//	}
//
//	loc = type.find("qyqb",0);
//	if (loc != std::string::npos){
//		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_qb).helicity(),pro.p(pos_q).flavor()+100,false));
//		propagators.push_back(particle_ID(quark_massive,pro.p(pos_qb).helicity(),pro.p(pos_q).flavor()+100,false));
//	}
//	return process(temp);
//}

void arrange_possible_props_2q1y_massive(const process& pro, vector<particle_ID>& possible_props){
	possible_props.push_back(gsc);

	particle_type_match is_photon(photon);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"Q");
	rule R3(&is_photon,"y");

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string process_name=string_gen(pro,rules);
	process_name+=process_name;
	string::size_type loc;
	loc=process_name.find("qyQ",0);
	process::const_iterator pos_q=find_if(pro.begin(),pro.end(),is_of_type_and_pap(quark,false));
	if (loc != string::npos){
		possible_props.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor()+100,false));
		possible_props.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor()+100,false));
	} else {
		possible_props.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor()+100,true));
		possible_props.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor()+100,true));
	}
}



process arrange_flavors_2q2G_LC_massive(const process& pro,vector<particle_ID>& propagators){
	vector<particle_ID> temp;
		for (int i=1;i<=pro.n();i++){
			temp.push_back(pro.p(i));
		}

	propagators.push_back(gsc);

	process::const_iterator quark_pos=find_if(pro.begin(),pro.end(),is_of_type(quark));

	cyclic_iterator<particle_ID,process > c1(pro,2,quark_pos);

	// Any of the quarks is our starting point
	// Next we find a gluino1
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > first_G1(c1);
	// Find the second gluino1
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > second_G1(c1);
	// then the first quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > first_q(c1);
	// then the second quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);

	
	propagators.push_back(particle_ID(gluino_massive,(*first_G1).helicity(),(*first_G1).flavor(),((*first_G1).is_anti_particle())));
	propagators.push_back(particle_ID(gluino_massive,-(*first_G1).helicity(),(*first_G1).flavor(),((*first_G1).is_anti_particle())));
	propagators.push_back(particle_ID(quark_massive,(*first_q).helicity(),(*first_q).flavor()+100,((*first_q).is_anti_particle())));
	propagators.push_back(particle_ID(quark_massive,-(*first_q).helicity(),(*first_q).flavor()+100,((*first_q).is_anti_particle())));
	
//	_PRINT(pro);
//	_PRINT(propagators);

	return process(temp);
}



process arrange_flavors_2q4G_LC_massive(const process& pro,vector<particle_ID>& propagators){
	vector<particle_ID> temp;
		for (int i=1;i<=pro.n();i++){
			temp.push_back(pro.p(i));
		}

	propagators.push_back(gsc);

	process::const_iterator quark_pos=find_if(pro.begin(),pro.end(),is_of_type(quark));

	cyclic_iterator<particle_ID,process > c1(pro,2,quark_pos);

	// Any of the quarks is our starting point
	// Next we find a gluino1
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > first_G1(c1);
	// Find the second gluino1
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > second_G1(c1);
	// then the gluino2
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > first_G2(c1);
	// then the gluino2
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > second_G2(c1);
	// then the first quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > first_q(c1);
	// then the second quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);

	
	propagators.push_back(particle_ID(gluino_massive,(*first_G1).helicity(),(*first_G1).flavor(),((*first_G1).is_anti_particle())));
	propagators.push_back(particle_ID(gluino_massive,-(*first_G1).helicity(),(*first_G1).flavor(),((*first_G1).is_anti_particle())));
	propagators.push_back(particle_ID(gluino_massive,(*first_G2).helicity(),(*first_G2).flavor(),((*first_G2).is_anti_particle())));
	propagators.push_back(particle_ID(gluino_massive,-(*first_G2).helicity(),(*first_G2).flavor(),((*first_G2).is_anti_particle())));
	propagators.push_back(particle_ID(quark_massive,(*first_q).helicity(),(*first_q).flavor()+100,((*first_q).is_anti_particle())));
	propagators.push_back(particle_ID(quark_massive,-(*first_q).helicity(),(*first_q).flavor()+100,((*first_q).is_anti_particle())));
	
//	_PRINT(pro);
//	_PRINT(propagators);

	return process(temp);
}

process arrange_flavors_4q_massive_LC(const process& pro,vector<particle_ID>& propagators){
	// assumed are 2 distinct flavors
	vector<particle_ID> temp;
		for (int i=1;i<=pro.n();i++){
			temp.push_back(pro.p(i));
		}

	propagators.push_back(gsc);

	process::const_iterator quark_pos=find_if(pro.begin(),pro.end(),is_of_type(quark));

	cyclic_iterator<particle_ID,process > c1(pro,2,quark_pos);


	// Any of the quarks is our starting point
	// then the first quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > first_q(c1);
	// Find the second quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);
	// then the third quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > third_q(c1);
	// then the fourth quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > fourth_q(c1);

	if((*first_q).flavor()==(*second_q).flavor()){
	propagators.push_back(particle_ID(quark_massive,(*first_q).helicity(),(*first_q).flavor(),((*first_q).is_anti_particle())));
	propagators.push_back(particle_ID(quark_massive,(*third_q).helicity(),(*third_q).flavor(),((*third_q).is_anti_particle())));
	}
	else {
	propagators.push_back(particle_ID(quark_massive,(*second_q).helicity(),(*second_q).flavor(),((*second_q).is_anti_particle())));
	propagators.push_back(particle_ID(quark_massive,(*fourth_q).helicity(),(*fourth_q).flavor(),((*fourth_q).is_anti_particle())));
	}

//	_PRINT(pro);
//	_PRINT(propagators);

	return process(temp);
}

process arrange_flavors_6q_massive_LC(const process& pro,vector<particle_ID>& propagators){
	// assumed are 3 distinct flavors
	vector<particle_ID> temp;
		for (int i=1;i<=pro.n();i++){
			temp.push_back(pro.p(i));
		}

	propagators.push_back(gsc);

	process::const_iterator quark_pos=find_if(pro.begin(),pro.end(),is_of_type(quark));

	cyclic_iterator<particle_ID,process > c1(pro,2,quark_pos);

	// Any of the quarks is our starting point
	// Next we first quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > first_q(c1);
	// Find the second quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);
	// then the third quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > third_q(c1);
	// then the fourth quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > fourth_q(c1);
	// then the fifth quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > fifth_q(c1);
	// then the sixth quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > sixth_q(c1);

	//check if flavors match in pairs like: 12,34,56 or like: 23,45,61
	if((*first_q).flavor()==(*second_q).flavor()){
		propagators.push_back(particle_ID(quark_massive,(*first_q).helicity(),(*first_q).flavor(),((*first_q).is_anti_particle())));
		propagators.push_back(particle_ID(quark_massive,(*third_q).helicity(),(*third_q).flavor(),((*third_q).is_anti_particle())));
		propagators.push_back(particle_ID(quark_massive,(*fifth_q).helicity(),(*fifth_q).flavor(),((*fifth_q).is_anti_particle())));
	}
	else {
		propagators.push_back(particle_ID(quark_massive,(*second_q).helicity(),(*second_q).flavor(),((*second_q).is_anti_particle())));
		propagators.push_back(particle_ID(quark_massive,(*fourth_q).helicity(),(*fourth_q).flavor(),((*fourth_q).is_anti_particle())));
		propagators.push_back(particle_ID(quark_massive,(*sixth_q).helicity(),(*sixth_q).flavor(),((*sixth_q).is_anti_particle())));
	}
//	_PRINT(pro);
//	_PRINT(propagators);

	return process(temp);
}

process arrange_flavors_2q2G1y_massive(const process& pro,vector<particle_ID>& propagators){
//	_PRINT(pro);
	vector<particle_ID> temp;
		size_t pos_qb=0;
		size_t pos_q=0;

		for (int i=1;i<=pro.n();i++){
			if (pro.p(i).is_a(quark)){
				if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
			}
			temp.push_back(pro.p(i));
		}

		particle_type_match is_a_gluon(gluon);
		particle_type_match is_a_gluino(gluino);
		particle_type_match is_a_quark(quark);
		particle_type_match is_a_photon(photon);
		type_and_pap_match is_quark(quark,false);
		type_and_pap_match is_antiquark(quark,true);




//	_PRINT(type);

	std::string::size_type loc,loc1,loc2,loc3;

	propagators.push_back(gsc);

	process::const_iterator photon_pos=find_if(pro.begin(),pro.end(),is_of_type(photon));

	cyclic_iterator<particle_ID,process> cy(pro,1,photon_pos);
	cyclic_iterator<particle_ID,process> cy_end(pro);

	cyclic_iterator<particle_ID,process> first_G=find_if(cy,cy_end,is_of_type(gluino));

	bool sub_leading_configuration;

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"qb");
	rule R3(&is_a_photon,"y");
	rule R4(&is_a_gluino,"G");
	rule R5(&is_a_quark,"q");

	vector<rule> rules_L_or_SL;
	rules_L_or_SL.push_back(R3);
	rules_L_or_SL.push_back(R4);
	rules_L_or_SL.push_back(R5);

	std::string type_L_or_SL=string_gen(pro,rules_L_or_SL);
	type_L_or_SL+=type_L_or_SL;
	loc1 = type_L_or_SL.find("qqGGy",0);
	loc2 = type_L_or_SL.find("GGqqy",0);
	loc3 = type_L_or_SL.find("qGGqy",0);



	if (loc1 != std::string::npos || loc2 != std::string::npos){
		sub_leading_configuration=true;
	} else 	if (loc3 != std::string::npos ){
		sub_leading_configuration=false;
	} else {
		_WARNING("no type found in arrange_flavors_2q2G1y");
	}

	if ( ! sub_leading_configuration ){
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,true));
//			propagators.push_back(Gbp);
//			propagators.push_back(Gbm);
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,false));
//			propagators.push_back(Gp);
//			propagators.push_back(Gm);
		}
	}
	else {
		if ((*first_G).is_anti_particle()){
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,false));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,false));
		} else {
			propagators.push_back(particle_ID(gluino_massive,+1,(*first_G).flavor()+100,true));
			propagators.push_back(particle_ID(gluino_massive,-1,(*first_G).flavor()+100,true));
		}
	}

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string type=string_gen(pro,rules);
	type+=type;
	loc = type.find("qbyq",0);
	if (loc != std::string::npos){
		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_q).helicity(),pro.p(pos_q).flavor()+100,true));
	}

	loc = type.find("qyqb",0);
	if (loc != std::string::npos){
		propagators.push_back(particle_ID(quark_massive,-pro.p(pos_qb).helicity(),pro.p(pos_q).flavor()+100,false));
	}
	return process(temp);
}

process arrange_flavors_qqX_massive_LC(const process& pro,vector<particle_ID>& propagators){
	process::const_iterator pos_q=find_if(pro.begin(),pro.end(),is_of_type_and_pap(quark,false));
	propagators.push_back(gsc);

	particle_type_any_unordered is_photon;
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"Q");
	rule R3(&is_photon,"y");

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string process_name=string_gen(pro,rules);
	string::size_type loc;
	loc=process_name.find("qyQ",0);
	if (loc != string::npos){
		propagators.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor(),false));
		propagators.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor(),false));
	} else {
		propagators.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor(),true));
		propagators.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor(),true));
	}

	return pro;
}

process arrange_flavors_qqQQX_massive_LC(const process& pro,vector<particle_ID>& propagators){
	propagators.push_back(gsc);

	// Find the anti-quark as this is our starting point
	process::const_iterator ph=find_if(pro.begin(),pro.end(),is_of_type_and_pap(quark,true));
	cyclic_iterator<particle_ID,process > c1(pro,4,ph);
	cyclic_iterator<particle_ID,process > first_q(c1);
	// Find the second quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);
	// then the 3rd quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > third_q(c1);
	// Find the forth quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > forth_q(c1);

	// Then we need to find which quarks pair up as quark lines, we can have either (1st-2nd)(3rd-4th)
	//  or (2nd-3rd)(4th-1st)
	cyclic_iterator<particle_ID,process > ql1t(c1),ql1b(c1),ql2t(c1),ql2b(c1); // NOTE: We only set these to c1 as there is no () constructor? Is this true Daniel?
	if((*first_q).flavor()==(*second_q).flavor()){ // (1st-2nd)(3rd-4th)
		ql1t=first_q;
		ql1b=second_q;
		ql2t=third_q;
		ql2b=forth_q;
	}
	else{ // (2nd-3rd)(4th-1st)
		ql1t=second_q;
		ql1b=third_q;
		ql2t=forth_q;
		ql2b=first_q;
	}

	if ((*ql1t).is_anti_particle()){
//		propagators.push_back(particle_ID(quark,+1,(*ql1t).flavor()+100,true));
//		propagators.push_back(particle_ID(quark,-1,(*ql1t).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,+1,(*ql1t).flavor(),true));
		propagators.push_back(particle_ID(quark_massive,-1,(*ql1t).flavor(),true));
	} else {
//		propagators.push_back(particle_ID(quark,+1,(*ql1t).flavor()+100,false));
//		propagators.push_back(particle_ID(quark,-1,(*ql1t).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,+1,(*ql1t).flavor(),false));
		propagators.push_back(particle_ID(quark_massive,-1,(*ql1t).flavor(),false));
	}

	if ((*ql2t).is_anti_particle()){
//		propagators.push_back(particle_ID(quark,+1,(*ql2t).flavor()+100,true));
//		propagators.push_back(particle_ID(quark,-1,(*ql2t).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,+1,(*ql2t).flavor(),true));
		propagators.push_back(particle_ID(quark_massive,-1,(*ql2t).flavor(),true));
	} else {
//		propagators.push_back(particle_ID(quark,+1,(*ql2t).flavor()+100,false));
//		propagators.push_back(particle_ID(quark,-1,(*ql2t).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,+1,(*ql2t).flavor(),false));
		propagators.push_back(particle_ID(quark_massive,-1,(*ql2t).flavor(),false));
	}

	return pro;
}


process arrange_flavors_qqX_massive_SLC(const process& pro,vector<particle_ID>& propagators){
	process::const_iterator pos_q=find_if(pro.begin(),pro.end(),is_of_type_and_pap(quark,false));
	propagators.push_back(gsc);

	particle_type_any_unordered is_photon;
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"Q");
	rule R3(&is_photon,"y");

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string process_name=string_gen(pro,rules);
	string::size_type loc;
	loc=process_name.find("qyQ",0);
	process_name+=process_name;
	if (loc != string::npos){
		propagators.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor()+100,false));
	} else {
		propagators.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor()+100,true));
	}

	return pro;
}

process arrange_flavors_qqee_massive_SLC(const process& pro,vector<particle_ID>& propagators){
	process::const_iterator pos_q=find_if(pro.begin(),pro.end(),is_of_type_and_pap(quark,false));
	process::const_iterator pos_G=find_if(pro.begin(),pro.end(),is_of_type_and_pap(gluino,false));
	propagators.push_back(gsc);

	particle_type_match is_photon(photon);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);
	type_and_pap_match is_gluino(gluino,false);
	type_and_pap_match is_antigluino(gluino,true);

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"Q");
	rule R3(&is_gluino,"G");
	rule R4(&is_antigluino,"Gb");

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);
	rules.push_back(R4);

	std::string process_name=string_gen(pro,rules);
	string::size_type loc1,loc2,loc3,loc4;
	loc1=process_name.find("qQGGb",0);
	loc2=process_name.find("QqGGb",0);
	loc3=process_name.find("qQGbG",0);
	loc4=process_name.find("QqGbG",0);
	if (loc2 != string::npos || loc4 != string::npos ){
		propagators.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor()+100,false));
	} else {
		propagators.push_back(particle_ID(quark_massive,+1,(*pos_q).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,-1,(*pos_q).flavor()+100,true));
	}
	if (loc3 != string::npos || loc4 != string::npos ){
		propagators.push_back(particle_ID(gluino_massive,+1,(*pos_G).flavor()+100,false));
		propagators.push_back(particle_ID(gluino_massive,-1,(*pos_G).flavor()+100,false));
	} else {
		propagators.push_back(particle_ID(gluino_massive,+1,(*pos_G).flavor()+100,true));
		propagators.push_back(particle_ID(gluino_massive,-1,(*pos_G).flavor()+100,true));
	}

	return pro;
}


process arrange_flavors_qqQQX_massive_SLC(const process& pro,vector<particle_ID>& propagators){
	propagators.push_back(gsc);

	// Find the anti-quark as this is our starting point
	process::const_iterator ph=find_if(pro.begin(),pro.end(),is_of_type_and_pap(quark,true));
	cyclic_iterator<particle_ID,process > c1(pro,4,ph);
	cyclic_iterator<particle_ID,process > first_q(c1);
	// Find the second quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);
	// then the 3rd quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > third_q(c1);
	// Find the forth quark
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > forth_q(c1);

	// Then we need to find which quarks pair up as quark lines, we can have either (1st-2nd)(3rd-4th)
	//  or (2nd-3rd)(4th-1st)
	cyclic_iterator<particle_ID,process > ql1t(c1),ql1b(c1),ql2t(c1),ql2b(c1); // NOTE: We only set these to c1 as there is no () constructor? Is this true Daniel?
	if((*first_q).flavor()==(*second_q).flavor()){ // (1st-2nd)(3rd-4th)
		ql1t=first_q;
		ql1b=second_q;
		ql2t=third_q;
		ql2b=forth_q;
	}
	else{ // (2nd-3rd)(4th-1st)
		ql1t=second_q;
		ql1b=third_q;
		ql2t=forth_q;
		ql2b=first_q;
	}

	if (!(*ql1t).is_anti_particle()){
		propagators.push_back(particle_ID(quark_massive,+1,(*ql1t).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,-1,(*ql1t).flavor()+100,true));
//		propagators.push_back(particle_ID(quark_massive,+1,(*ql1t).flavor(),true));
//		propagators.push_back(particle_ID(quark_massive,-1,(*ql1t).flavor(),true));
	} else {
		propagators.push_back(particle_ID(quark_massive,+1,(*ql1t).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,-1,(*ql1t).flavor()+100,false));
//		propagators.push_back(particle_ID(quark_massive,+1,(*ql1t).flavor(),false));
//		propagators.push_back(particle_ID(quark_massive,-1,(*ql1t).flavor(),false));
	}

	if (!(*ql2t).is_anti_particle()){
		propagators.push_back(particle_ID(quark_massive,+1,(*ql2t).flavor()+100,true));
		propagators.push_back(particle_ID(quark_massive,-1,(*ql2t).flavor()+100,true));
//		propagators.push_back(particle_ID(quark_massive,+1,(*ql2t).flavor(),true));
//		propagators.push_back(particle_ID(quark_massive,-1,(*ql2t).flavor(),true));
	} else {
		propagators.push_back(particle_ID(quark_massive,+1,(*ql2t).flavor()+100,false));
		propagators.push_back(particle_ID(quark_massive,-1,(*ql2t).flavor()+100,false));
//		propagators.push_back(particle_ID(quark_massive,+1,(*ql2t).flavor(),false));
//		propagators.push_back(particle_ID(quark_massive,-1,(*ql2t).flavor(),false  ));
	}

	return pro;
}

}

