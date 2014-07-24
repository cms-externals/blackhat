/*
 * process_utils_basic.cpp
 *
 *  Created on: 02.08.2012
 *      Author: daniel
 */

#include "process_utils.h"
#include <vector>
#include "particles.h"
#include "process.h"
#include "partitions.h"
#include "BH_utilities.h"
#include "process_utils_massive.h"
#ifndef BH_PUBLIC
#include "options.h"
#include "parent_diagrams.h"
#endif
#include "iterators.h"
#include "arrange_flavors.h"
#include <algorithm>
#ifndef BH_PUBLIC
#include "parent_diagrams_from_file.h"
#endif
#include "settings.h"
#include "BH_debug.h"

using namespace std;

namespace BH {

size_t count_ph(const process& p, const particle& t); //In amplitudes.cpp


vector<int> find_all_flavors(const process& pro,const particle& type){
	size_t count=0;
  	vector<int> flavors;
  	for(int i=1;i<=pro.n();i++)
  	{
  		if( (*pro.p(i).type()) == type ) {
  			vector<int>::iterator pos= find ( flavors.begin(),flavors.end(),pro.p(i).flavor() ) ;
  	    	if ( pos == flavors.end() ){
  	    		count++;
  	    		flavors.push_back(pro.p(i).flavor());
  	    	}
		}
  	}
  	return flavors;
}

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

process replace_one_photon_with_gluon(const process& pro){
	bool done=false;
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(photon)){
			if (!done){
				pp.push_back(particle_ID(gluon,pro.p(i).helicity(),1,false));
				done=true;
			} else {
				pp.push_back(pro.p(i));
			}
		} else {
			pp.push_back(pro.p(i));
		}
	}
	return process(pp);
}


process replace_lepton_with_quark(const process& pro){
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(lepton)){
			pp.push_back(particle_ID(quark,pro.p(i).helicity(),pro.p(i).flavor(),pro.p(i).is_anti_particle()));
		} else {
			pp.push_back(pro.p(i));
		}
	}
	return process(pp);
}


process replace_photon_with_gluon(const process& pro){
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(photon)){
			pp.push_back(particle_ID(gluon,pro.p(i).helicity(),1,false));
		} else {
			pp.push_back(pro.p(i));
		}
	}
	return process(pp);
}



process replace_gluino_with_quark(const process& pro){
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(gluino)){
			pp.push_back(particle_ID(quark,pro.p(i).helicity(),pro.p(i).flavor(),pro.p(i).is_anti_particle()));
		} else if (pro.p(i).is_a(gluino_massive)){
			pp.push_back(particle_ID(quark_massive,pro.p(i).helicity(),pro.p(i).flavor(),pro.p(i).is_anti_particle()));
		} else {
			pp.push_back(pro.p(i));
		}
	}
	return process(pp);
}
process replace_gluino_with_quark(const process& pro,int flavor){
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(gluino)){
			pp.push_back(particle_ID(quark,pro.p(i).helicity(),flavor,pro.p(i).is_anti_particle()));
		} else if (pro.p(i).is_a(gluino_massive)){
			pp.push_back(particle_ID(quark_massive,pro.p(i).helicity(),flavor,pro.p(i).is_anti_particle()));
		} else {
			pp.push_back(pro.p(i));
		}
	}
#if _VERBOSE
	_MESSAGE3(pro ," changed into ",process(pp));
#endif
	return process(pp);
}

process replace_scalar_massive_with_gluon(const process& pro){
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(gluon_massive_scalar)){
			pp.push_back(particle_ID(scalar_massive,pro.p(i).helicity(),pro.p(i).flavor(),pro.p(i).is_anti_particle()));
		} else {
			pp.push_back(pro.p(i));
		}
	}
	return process(pp);
}

process replace_massive_scalars_with_gluons(const process& pro){
    vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(gluon_massive_scalar)||pro.p(i).is_a(gluon_massive)){
			pp.push_back(particle_ID(gluon,pro.p(i).helicity(),pro.p(i).flavor(),pro.p(i).is_anti_particle()));
		} else {
			pp.push_back(pro.p(i));
		}
	}
	return process(pp);
}

process replace_gluon_massive_with_scalar(const process& pro){
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(gluon_massive_scalar)||pro.p(i).is_a(gluon_massive)){
			pp.push_back(particle_ID(scalar_massive,pro.p(i).helicity(),pro.p(i).flavor(),pro.p(i).is_anti_particle()));
		} else {
			pp.push_back(pro.p(i));
		}
	}
	return process(pp);
}



process fix_flavors_gluons(const process& pro){
//	cout << "trying: " << pro<< ":" ;
	size_t nbr_gluons=0;
	size_t nbr_quarks=0;
	short q_flavor,qb_flavor;
	vector<particle_ID> temp;
	bool need_fixing=false;

	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(gluon) && pro.p(i).flavor() < 0 ) { need_fixing = true; }
		if(pro.p(i).is_a(quark)||pro.p(i).is_a(quark_massive)||pro.p(i).is_a(gluino_massive)){ nbr_quarks++; }
		temp.push_back(pro.p(i));
	}
	if ( need_fixing == false ) return pro;
	switch (nbr_quarks){
	case 2 : {
		size_t pos_qb=0;
		size_t pos_q=0;
		size_t pos_qb2=0;
		size_t pos_q2=0;
		nbr_quarks=2; nbr_gluons=pro.n()-2;
		for (int i=1;i<=pro.n();i++){
			if ((pro.p(i).is_a(quark)||pro.p(i).is_a(quark_massive)||pro.p(i).is_a(gluino_massive)) && pro.p(i).is_anti_particle() ) { pos_qb = i ; qb_flavor=pro.p(i).flavor();  }
			if ((pro.p(i).is_a(quark)||pro.p(i).is_a(quark_massive)||pro.p(i).is_a(gluino_massive)) && (!(pro.p(i).is_anti_particle()))) { pos_q = i ; q_flavor=pro.p(i).flavor(); }

		}
		if (pos_qb == 0 || pos_q == 0) return pro;

		size_t previous=(pro.n()+pos_qb-2)%pro.n()+1;
		size_t next    =(pro.n()+pos_qb)%pro.n()+1;

		if (pro.p(next).is_a(gluon) && pro.p(next).flavor()==-qb_flavor ){  //
			size_t distance = (pos_q-pos_qb+pro.n())%pro.n();
			if ( q_flavor-qb_flavor != distance-1 ) return pro;

			for (int j=0;j<distance-1;j++){
				size_t pos=(pro.n()+pos_qb+j)%pro.n()+1;
				if (pro.p(pos).is_a(gluon) && pro.p(pos).flavor()==-qb_flavor - j ){
					temp[pos-1]=particle_ID(gluon,pro.p(pos).helicity(),1,false);
				}
				else {return pro;}
			}
		}

		else { if (pro.p(previous).is_a(gluon) && pro.p(previous).flavor()==-qb_flavor ){  //
			size_t distance = (pos_qb-pos_q+pro.n())%pro.n();
			if ( q_flavor-qb_flavor != distance-1 ) return pro;

			for (int j=0;j<distance-1;j++){
				size_t pos=(pro.n()+pos_qb-j-2)%pro.n()+1;
//				_PRINT(pos);
				if (pro.p(pos).is_a(gluon) && pro.p(pos).flavor()==-qb_flavor - j ){
					temp[pos-1]=particle_ID(gluon,pro.p(pos).helicity(),1,false);
				}
				else {return pro;}
			}
		}

		else { return pro; }  // neither the particle before or after is a gluon with negative flavor
		}
		temp[pos_q-1]=particle_ID(pro.p(pos_q).type(),pro.p(pos_q).helicity(),1,false);
		temp[pos_qb-1]=particle_ID(pro.p(pos_qb).type(),pro.p(pos_qb).helicity(),1,true);
		return process(temp);
		break;
		}
	case 4 : {
		nbr_quarks=4; nbr_gluons=pro.n()-4;
		vector<size_t> quarks_pos;
		vector<size_t> antiquarks_pos;
		for (int i=1;i<=pro.n();i++){
			if (pro.p(i).is_a(quark)||pro.p(i).is_a(quark_massive)||pro.p(i).is_a(gluino_massive)){
				if ( pro.p(i).is_anti_particle() ) { antiquarks_pos.push_back(i) ;  } else {  quarks_pos.push_back(i) ;  }
			}
		}
		if (antiquarks_pos.size() !=2 || quarks_pos.size() !=2 ) return pro;
		size_t aq1_pos, aq2_pos, q1_pos, q2_pos, q_flavor1,q_flavor2,qb_flavor1,qb_flavor2;

		if (pro.p(antiquarks_pos[0]).flavor() < pro.p(antiquarks_pos[1]).flavor()){
			temp[antiquarks_pos[0]-1]=particle_ID(temp[antiquarks_pos[0]-1].type(),temp[antiquarks_pos[0]-1].helicity(),1,true);
			temp[antiquarks_pos[1]-1]=particle_ID(temp[antiquarks_pos[1]-1].type(),temp[antiquarks_pos[1]-1].helicity(),2,true);
			aq1_pos=antiquarks_pos[0];
			aq2_pos=antiquarks_pos[1];
			qb_flavor1=pro.p(antiquarks_pos[0]).flavor();
			qb_flavor2=pro.p(antiquarks_pos[1]).flavor();
		} else {
			temp[antiquarks_pos[0]-1]=particle_ID(temp[antiquarks_pos[0]-1].type(),temp[antiquarks_pos[0]-1].helicity(),2,true);
			temp[antiquarks_pos[1]-1]=particle_ID(temp[antiquarks_pos[1]-1].type(),temp[antiquarks_pos[1]-1].helicity(),1,true);
			aq1_pos=antiquarks_pos[1];
			aq2_pos=antiquarks_pos[0];
			qb_flavor1=pro.p(antiquarks_pos[1]).flavor();
			qb_flavor2=pro.p(antiquarks_pos[0]).flavor();
		}
		if (pro.p(quarks_pos[0]).flavor() < pro.p(quarks_pos[1]).flavor()){
			temp[quarks_pos[0]-1]=particle_ID(temp[quarks_pos[0]-1].type(),temp[quarks_pos[0]-1].helicity(),1,false);
			temp[quarks_pos[1]-1]=particle_ID(temp[quarks_pos[1]-1].type(),temp[quarks_pos[1]-1].helicity(),2,false);
			q1_pos=quarks_pos[0];
			q2_pos=quarks_pos[1];
			q_flavor1=pro.p(quarks_pos[0]).flavor();
			q_flavor2=pro.p(quarks_pos[1]).flavor();
		} else {
			temp[quarks_pos[0]-1]=particle_ID(temp[quarks_pos[0]-1].type(),temp[quarks_pos[0]-1].helicity(),2,false);
			temp[quarks_pos[1]-1]=particle_ID(temp[quarks_pos[1]-1].type(),temp[quarks_pos[1]-1].helicity(),1,false);
			q1_pos=quarks_pos[1];
			q2_pos=quarks_pos[0];
			q_flavor1=pro.p(quarks_pos[1]).flavor();
			q_flavor2=pro.p(quarks_pos[0]).flavor();
		}

//		_PRINT(q1_pos);_PRINT(q2_pos);
//		_PRINT(aq1_pos);_PRINT(aq2_pos);
		size_t previous1=(pro.n()+aq1_pos-2)%pro.n()+1;
		size_t next1    =(pro.n()+aq1_pos)%pro.n()+1;

		if (pro.p(next1).is_a(gluon) && pro.p(next1).flavor()==-qb_flavor1 ){  //
			size_t distance = (q1_pos-aq1_pos+pro.n())%pro.n();
			if (  (q_flavor2 != qb_flavor2) && (q_flavor1-qb_flavor1 != distance-1) ) return pro;
			size_t next_gluon_index=-qb_flavor1;
			for (int j=0;j<distance-1;j++){
				size_t pos=(pro.n()+aq1_pos+j)%pro.n()+1;
				if (pro.p(pos).is_a(gluon) && pro.p(pos).flavor()==next_gluon_index ){
					temp[pos-1]=particle_ID(gluon,pro.p(pos).helicity(),1,false);
					next_gluon_index--;
				}
				else if ((pro.p(pos).is_a(quark)||pro.p(pos).is_a(quark_massive)) && (q_flavor2 == qb_flavor2) ){
					// do nothing
				}
				else {return pro;}
			}
		}

	 if (pro.p(previous1).is_a(gluon) && pro.p(previous1).flavor()==-qb_flavor1 ){  //
			size_t distance = (aq1_pos-q1_pos+pro.n())%pro.n();
			if ( (q_flavor2 != qb_flavor2) &&( q_flavor1-qb_flavor1 != distance-1 )) return pro;
			size_t next_gluon_index=-qb_flavor1;

			for (int j=0;j<distance-1;j++){
				size_t pos=(pro.n()+aq1_pos-j-2)%pro.n()+1;
				if (pro.p(pos).is_a(gluon) && pro.p(pos).flavor()==next_gluon_index ){
					temp[pos-1]=particle_ID(gluon,pro.p(pos).helicity(),1,false);
					next_gluon_index--;
				}
				else {return pro;}
			}
		}

	 size_t previous2=(pro.n()+aq2_pos-2)%pro.n()+1;
		size_t next2    =(pro.n()+aq2_pos)%pro.n()+1;
//		_PRINT(previous2);
//		_PRINT(next2);
		if (pro.p(next2).is_a(gluon) && pro.p(next2).flavor()==-qb_flavor2 ){  //
			size_t distance = (q2_pos-aq2_pos+pro.n())%pro.n();
			if ( (q_flavor1 != qb_flavor1) && (q_flavor2-qb_flavor2 != distance-1) ) return pro;
			size_t next_gluon_index=-qb_flavor2;

			for (int j=0;j<distance-1;j++){
				size_t pos=(pro.n()+aq2_pos+j)%pro.n()+1;
				if (pro.p(pos).is_a(gluon) && pro.p(pos).flavor()==next_gluon_index ){
					temp[pos-1]=particle_ID(gluon,pro.p(pos).helicity(),1,false);
					next_gluon_index--;
				}
				else if ((pro.p(pos).is_a(quark)||pro.p(pos).is_a(quark_massive)) && (q_flavor1 == qb_flavor1)) {
					// do nothing
				}
				else {return pro;}
			}
		}

	 if (pro.p(previous2).is_a(gluon) && pro.p(previous2).flavor()==-qb_flavor2 ){  //
			size_t distance = (aq2_pos-q2_pos+pro.n())%pro.n();
//			_PRINT(distance);
			if ( (q_flavor1 != qb_flavor1) && (q_flavor2-qb_flavor2 != distance-1) ) return pro;
			size_t next_gluon_index=-qb_flavor2;
			for (int j=0;j<distance-1;j++){
				size_t pos=(pro.n()+aq2_pos-j-2)%pro.n()+1;
//				_PRINT(pos);
				if (pro.p(pos).is_a(gluon) && pro.p(pos).flavor()==next_gluon_index ){
					temp[pos-1]=particle_ID(gluon,pro.p(pos).helicity(),1,false);
				next_gluon_index--;
				}
				else if ((pro.p(pos).is_a(quark)||pro.p(pos).is_a(quark_massive)) && (q_flavor1 == qb_flavor1)) {
					// do nothing
				}
				else {return pro;}
			}
		}

//		_MESSAGE3(pro," --> ",process(temp));
		return process(temp);
	}

	default: return pro;
	}

	//We will never get here but to avoid compiler warnings
	return pro;
}

process fix_flavors_2q2e(const process& pro){
	vector<particle_ID> temp;

//	cout << "Fix flavours 2q2l: " << pro<< ":" ;

	size_t pos_qb=0;
	size_t pos_q=0;
	size_t pos_lb=0;
	size_t pos_l=0;
	for (int i=1;i<=pro.n();i++){
		if (pro.p(i).is_a(quark)||pro.p(i).is_a(quark_massive)){
			if ( pro.p(i).is_anti_particle() ) { pos_qb = i ;  } else {  pos_q = i ;  }
		}
		if (pro.p(i).is_a(lepton)){
			if ( pro.p(i).is_anti_particle() ) { pos_lb = i ;  } else {  pos_l = i ;  }
		}
		temp.push_back(pro.p(i));
	}
//	_MESSAGE8("Found all legs",pos_q,", ",pos_qb,", ",pos_l,", ",pos_lb);
	// if we have a process with qb qb or q q we don't need to fix anything
	if ( pos_q == 0 || pos_qb == 0 ) return pro;

	//Otherwise
	if ( (pro.p(pos_qb).flavor() < 0) && (pro.p(pos_q).flavor() < 0)  ) {

		if ( pro.p(pos_qb).flavor() == -1 && pro.p(pos_q).flavor() == -2) {
			temp[pos_qb-1]=particle_ID(pro.p(pos_qb).type(),pro.p(pos_qb).helicity(),1,true);
			temp[pos_q-1]=particle_ID(pro.p(pos_q).type(),pro.p(pos_q).helicity(),1,false);
			return process(temp);
		}


	}
//	_MESSAGE2("Fixed Process=",process(temp));
	return process(temp);


}


process fix_flavors_4q2e(const process& pro){
	vector<particle_ID> temp;


	size_t pos_lb=0;
	size_t pos_l=0;

	for (int i=1;i<=pro.n();i++){

		if (pro.p(i).is_a(lepton)){
			if ( pro.p(i).is_anti_particle() ) { pos_lb = i ;  } else {  pos_l = i ;  }
		}
		temp.push_back(pro.p(i));
	}

	size_t nbr_neg_flavors = count_if (pro.begin(),pro.end(),has_negative_flavor());

	if (nbr_neg_flavors == 4 ){
	int offset =0;
	while ( ! (pro.p((pos_lb+offset)%pro.n()+1).is_a(quark)||pro.p((pos_lb+offset)%pro.n()+1).is_a(quark_massive)  )){++offset;}
	size_t pos_q1=(pos_lb+offset)%pro.n()+1;
	offset =-2;
	while ( ! (pro.p((pro.n()+pos_l+offset)%pro.n()+1).is_a(quark)||pro.p((pro.n()+pos_l+offset)%pro.n()+1).is_a(quark_massive) )){--offset;}
	size_t pos_q2=(pro.n()+pos_l+offset)%pro.n()+1;

		if (	(
			(pro.p(pos_q1).flavor() == -1 && pro.p(pos_q2).flavor() == -2)
			||
			(pro.p(pos_q1).flavor() == -2 && pro.p(pos_q2).flavor() == -1)
			)
		&&
			(pro.p(pos_q1).is_anti_particle()!= pro.p(pos_q2).is_anti_particle())
		&&   // the lepton pair has to change the flavor of the a quark going into the loop and not of a lepton that "bounces" off the loop, since this factors out a closed fermion loop
			(pos_q2 == 1 || pos_q1 == pro.n() )

	){
				temp[pos_q1-1]=particle_ID(pro.p(pos_q1).type(),pro.p(pos_q1).helicity(),1,pro.p(pos_q1).is_anti_particle());
				temp[pos_q2-1]=particle_ID(pro.p(pos_q2).type(),pro.p(pos_q2).helicity(),1,pro.p(pos_q2).is_anti_particle());
				return process(temp);
	}
	}

	if (nbr_neg_flavors == 2 ){
		size_t pos_q1=0;
		size_t pos_q2=0;
		size_t pos_qb1=0;
		size_t pos_qb2=0;

		process::const_iterator lepton_pos=find_if(pro.begin(),pro.end(),is_of_type(lepton));

		cyclic_iterator<particle_ID,process> cy(pro,1,lepton_pos);
		cyclic_iterator<particle_ID,process> cy_end(pro);

		cyclic_iterator<particle_ID,process> first_q=find_if(cy,cy_end,has_negative_flavor());

		if ((*first_q).is_anti_particle()){
			pos_qb1=first_q.position();
		} else {
			pos_q1=first_q.position();
		}


		cyclic_iterator<particle_ID,process> second_q=find_if(++first_q,cy_end,has_negative_flavor());
		if ((*second_q).is_anti_particle()){
			pos_qb1=second_q.position();
		} else {
			pos_q1=second_q.position();
		}





		if ( pos_q1 == 0 || pos_qb1 == 0 ) return process(pro);
			if (
				(pro.p(pos_qb1).flavor() == -1 && pro.p(pos_q1).flavor() == -2)
//			&&   // the lepton pair has to change the flavor of the a quark going into the loop and not of a lepton that "bounces" off the loop, since this factors out a closed fermion loop
//				(pos_q1 == 1 || pos_q1 == pro.n() || pos_qb1 == 1 || pos_qb1 == pro.n() || ( (!pro.p(1).is_a(quark)) || !(pro.p(pro.n()).is_a(quark))) )

				){
				temp[pos_q1-1]=particle_ID(pro.p(pos_q1).type(),pro.p(pos_q1).helicity(),1,pro.p(pos_q1).is_anti_particle());
				temp[pos_qb1-1]=particle_ID(pro.p(pos_qb1).type(),pro.p(pos_qb1).helicity(),1,pro.p(pos_qb1).is_anti_particle());
//				copy(temp.begin(),temp.end(),ostream_iterator<particle_ID>(cout," "));
				return process(temp);
			}
	}


	return process(temp);


}

process fix_flavors_4q(const process& pro){
	vector<particle_ID> temp;
	copy(pro.begin(),pro.end(),back_inserter(temp));

	vector<int> all_flavors=find_all_flavors(pro,quark,quark_massive);

	// The cases we can have that we want to change are
	// q10x qx qx q10x->q1x q1x q2x q2x
	// qx qy qy q10x->q1x qy qy q1x
	// q10x qx qy q10y->q1x q1x q2y q2y
	if (all_flavors.size() < 2 ) return fix_flavors_gluons(pro);
	//if (all_flavors.size() < 2 || all_flavors.size() > 3  ) return fix_flavors_gluons(pro);

	if(all_flavors.size()==2){
		int large_flavor,small_flavor;
		if (all_flavors[0] > 100){
			large_flavor=all_flavors[0];
			small_flavor=all_flavors[1];
		} else if (all_flavors[1] > 100){
			large_flavor=all_flavors[1];
			small_flavor=all_flavors[0];
		} else return fix_flavors_gluons(pro);

		if ( (large_flavor % 100) != small_flavor) return fix_flavors_gluons(pro);

		process::const_iterator first_quark=find_if(pro.begin(),pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator second_quark=find_if(first_quark+1,pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator third_quark=find_if(second_quark+1,pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator fourth_quark=find_if(third_quark+1,pro.end(),is_of_either_type(quark,quark_massive));

		// 0 based:
		size_t q1_pos=first_quark-pro.begin();
		size_t q2_pos=second_quark-pro.begin();
		size_t q3_pos=third_quark-pro.begin();
		size_t q4_pos=fourth_quark-pro.begin();


		if ( (*first_quark).flavor() != (*second_quark).flavor() ){
			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark).helicity(),10+small_flavor,(*second_quark).is_anti_particle());
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark ).helicity(),20+small_flavor,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark).helicity(),20+small_flavor,(*fourth_quark).is_anti_particle());
		} else {
			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark).helicity(),20+small_flavor,(*second_quark).is_anti_particle());
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark ).helicity(),20+small_flavor,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark).helicity(),10+small_flavor,(*fourth_quark).is_anti_particle());
		}
	}
	else if(all_flavors.size()==3){
		sort(all_flavors.begin(),all_flavors.end());
		 //Only proceed if there is a possible 10x<->x match
		if(all_flavors[2]<100){return fix_flavors_gluons(pro);}
		if(all_flavors[1]>100){return fix_flavors_gluons(pro);}
		// Find if either of the flavours match the 100 flavour
		int small_flavor=0;
		if ( (all_flavors[2] % 100) == all_flavors[0]){
			small_flavor=all_flavors[0];
		}else if ( (all_flavors[2] % 100) == all_flavors[1]){
			small_flavor=all_flavors[1];
		}else{
			return fix_flavors_gluons(pro);
		}

		process::const_iterator first_quark=find_if(pro.begin(),pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator second_quark=find_if(first_quark+1,pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator third_quark=find_if(second_quark+1,pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator fourth_quark=find_if(third_quark+1,pro.end(),is_of_either_type(quark,quark_massive));

		// 0 based:
		size_t q1_pos=first_quark-pro.begin();
		size_t q2_pos=second_quark-pro.begin();
		size_t q3_pos=third_quark-pro.begin();
		size_t q4_pos=fourth_quark-pro.begin();

		if( (*first_quark).flavor() == (*second_quark).flavor() && ((*first_quark).flavor()!=small_flavor)){
			//Then it is the 3rd and 4th quarks
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark ).helicity(),10+small_flavor,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark).helicity(),10+small_flavor,(*fourth_quark).is_anti_particle());
		} else if( (*second_quark).flavor() == (*third_quark).flavor() && ((*second_quark).flavor()!=small_flavor) ){
			//Then it is the 1st and 4th quarks
			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark).helicity(),10+small_flavor,(*fourth_quark).is_anti_particle());
		} else if( (*third_quark).flavor() == (*fourth_quark).flavor()  && ((*third_quark).flavor()!=small_flavor)){
			//Then it is the 1st and 2nd quarks
			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark).helicity(),10+small_flavor,(*second_quark).is_anti_particle());
		}

		//If it's none of these cases something messed (i.e. we had 10x,x,x,y) up but we dont care and just return at this point
	}
	else if(all_flavors.size()==4){
		process::const_iterator first_quark=find_if(pro.begin(),pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator second_quark=find_if(first_quark+1,pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator third_quark=find_if(second_quark+1,pro.end(),is_of_either_type(quark,quark_massive));
		process::const_iterator fourth_quark=find_if(third_quark+1,pro.end(),is_of_either_type(quark,quark_massive));

		// 0 based:
		size_t q1_pos=first_quark-pro.begin();
		size_t q2_pos=second_quark-pro.begin();
		size_t q3_pos=third_quark-pro.begin();
		size_t q4_pos=fourth_quark-pro.begin();

		if( ((*first_quark).flavor()%100 == (*second_quark).flavor()%100)&&
			((*third_quark).flavor()%100 == (*fourth_quark).flavor()%100)){

			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+(*first_quark).flavor()%100,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark ).helicity(),10+(*second_quark).flavor()%100,(*second_quark ).is_anti_particle());
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark).helicity(),20+(*third_quark).flavor()%100,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark ).helicity(),20+(*fourth_quark).flavor()%100,(*fourth_quark ).is_anti_particle());

		}
		else if( ((*first_quark).flavor()%100 == (*fourth_quark).flavor()%100)&&
			((*third_quark).flavor()%100 == (*second_quark).flavor()%100)){

			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+(*first_quark).flavor()%100,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark ).helicity(),20+(*second_quark).flavor()%100,(*second_quark ).is_anti_particle());
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark).helicity(),20+(*third_quark).flavor()%100,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark ).helicity(),10+(*fourth_quark).flavor()%100,(*fourth_quark ).is_anti_particle());

		}
	}
	return fix_flavors_gluons(process(temp));

}

process fix_flavors_4G(const process& pro){
	vector<particle_ID> temp;
	copy(pro.begin(),pro.end(),back_inserter(temp));

	vector<int> all_flavors=find_all_flavors(pro,gluino,gluino_massive);

	// The cases we can have that we want to change are
	// q10x qx qx q10x->q1x q1x q2x q2x
	// qx qy qy q10x->q1x qy qy q1x
	// q10x qx qy q10y->q1x q1x q2y q2y
	//if (all_flavors.size() < 2 || all_flavors.size() > 3  ) return pro;
	if (all_flavors.size() < 2 ) return pro;

	if(all_flavors.size()==2){
		int large_flavor,small_flavor;
		if (all_flavors[0] > 100){
			large_flavor=all_flavors[0];
			small_flavor=all_flavors[1];
		} else if (all_flavors[1] > 100){
			large_flavor=all_flavors[1];
			small_flavor=all_flavors[0];
		} else return pro;

		if ( (large_flavor % 100) != small_flavor) return pro;

		process::const_iterator first_quark=find_if(pro.begin(),pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator second_quark=find_if(first_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator third_quark=find_if(second_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator fourth_quark=find_if(third_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));

		// 0 based:
		size_t q1_pos=first_quark-pro.begin();
		size_t q2_pos=second_quark-pro.begin();
		size_t q3_pos=third_quark-pro.begin();
		size_t q4_pos=fourth_quark-pro.begin();


		if ( (*first_quark).flavor() != (*second_quark).flavor() ){
			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark).helicity(),10+small_flavor,(*second_quark).is_anti_particle());
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark ).helicity(),20+small_flavor,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark).helicity(),20+small_flavor,(*fourth_quark).is_anti_particle());
		} else {
			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark).helicity(),20+small_flavor,(*second_quark).is_anti_particle());
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark ).helicity(),20+small_flavor,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark).helicity(),10+small_flavor,(*fourth_quark).is_anti_particle());
		}
	}
	else if(all_flavors.size()==3){
		sort(all_flavors.begin(),all_flavors.end());
		 //Only proceed if there is a possible 10x<->x match
		if(all_flavors[2]<100){return pro;}
		if(all_flavors[1]>100){return pro;}
		// Find if either of the flavours match the 100 flavour
		int small_flavor=0;
		if ( (all_flavors[2] % 100) == all_flavors[0]){
			small_flavor=all_flavors[0];
		}else if ( (all_flavors[2] % 100) == all_flavors[1]){
			small_flavor=all_flavors[1];
		}else{
			return pro;
		}

		process::const_iterator first_quark=find_if(pro.begin(),pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator second_quark=find_if(first_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator third_quark=find_if(second_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator fourth_quark=find_if(third_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));

		// 0 based:
		size_t q1_pos=first_quark-pro.begin();
		size_t q2_pos=second_quark-pro.begin();
		size_t q3_pos=third_quark-pro.begin();
		size_t q4_pos=fourth_quark-pro.begin();

		if( (*first_quark).flavor() == (*second_quark).flavor() && ((*first_quark).flavor()!=small_flavor)){
			//Then it is the 3rd and 4th quarks
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark ).helicity(),10+small_flavor,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark).helicity(),10+small_flavor,(*fourth_quark).is_anti_particle());
		} else if( (*second_quark).flavor() == (*third_quark).flavor() && ((*second_quark).flavor()!=small_flavor) ){
			//Then it is the 1st and 4th quarks
			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark).helicity(),10+small_flavor,(*fourth_quark).is_anti_particle());
		} else if( (*third_quark).flavor() == (*fourth_quark).flavor()  && ((*third_quark).flavor()!=small_flavor)){
			//Then it is the 1st and 2nd quarks
			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark).helicity(),10+small_flavor,(*second_quark).is_anti_particle());
		}

		//If it's none of these cases something messed (i.e. we had 10x,x,x,y) up but we dont care and just return at this point
	}
	else if(all_flavors.size()==4){
		process::const_iterator first_quark=find_if(pro.begin(),pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator second_quark=find_if(first_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator third_quark=find_if(second_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator fourth_quark=find_if(third_quark+1,pro.end(),is_of_either_type(gluino,gluino_massive));

		// 0 based:
		size_t q1_pos=first_quark-pro.begin();
		size_t q2_pos=second_quark-pro.begin();
		size_t q3_pos=third_quark-pro.begin();
		size_t q4_pos=fourth_quark-pro.begin();

		if( ((*first_quark).flavor()%100 == (*second_quark).flavor()%100)&&
			((*third_quark).flavor()%100 == (*fourth_quark).flavor()%100)){

			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+(*first_quark).flavor()%100,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark ).helicity(),10+(*second_quark).flavor()%100,(*second_quark ).is_anti_particle());
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark).helicity(),20+(*third_quark).flavor()%100,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark ).helicity(),20+(*fourth_quark).flavor()%100,(*fourth_quark ).is_anti_particle());

		}
		else if( ((*first_quark).flavor()%100 == (*fourth_quark).flavor()%100)&&
			((*third_quark).flavor()%100 == (*second_quark).flavor()%100)){

			temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+(*first_quark).flavor()%100,(*first_quark ).is_anti_particle());
			temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark ).helicity(),20+(*second_quark).flavor()%100,(*second_quark ).is_anti_particle());
			temp[q3_pos]=particle_ID((*third_quark).type(), (*third_quark).helicity(),20+(*third_quark).flavor()%100,(*third_quark ).is_anti_particle());
			temp[q4_pos]=particle_ID((*fourth_quark).type(), (*fourth_quark ).helicity(),10+(*fourth_quark).flavor()%100,(*fourth_quark ).is_anti_particle());

		}
	}
	return process(temp);
}


process fix_flavors_6G(const process& pro){
	vector<particle_ID> temp;
	copy(pro.begin(),pro.end(),back_inserter(temp));

	vector<int> all_flavors=find_all_flavors(pro,gluino,gluino_massive);



	// The cases we can have that we want to change are
	// q10x qx qz qz qx q10x->q1x q1x qz qz q2x q2x
	// qx qy qy qz qz q10x->q1x qy qy qz qz q1x
	// q10x qx qz qz qy q10y->q1x q1x qz qz q2y q2y
	//if (all_flavors.size() < 2 || all_flavors.size() > 3  ) return pro;
	if (all_flavors.size() < 2 ) return pro;


		process::const_iterator q1=find_if(pro.begin(),pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator q2=find_if(q1+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator q3=find_if(q2+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator q4=find_if(q3+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator q5=find_if(q4+1,pro.end(),is_of_either_type(gluino,gluino_massive));
		process::const_iterator q6=find_if(q5+1,pro.end(),is_of_either_type(gluino,gluino_massive));

		// 0 based:
		size_t q1_pos=q1-pro.begin();
		size_t q2_pos=q2-pro.begin();
		size_t q3_pos=q3-pro.begin();
		size_t q4_pos=q4-pro.begin();
		size_t q5_pos=q5-pro.begin();
		size_t q6_pos=q6-pro.begin();



	//1) q10x qx qz qz qx q10x -> q1x q1x qz qz q2x q2x
	//2) q10x qx qx qz qz q10x -> q1x q1x q2y qz qz q2y
	//3) q10x  qz qz qx qx q10x-> q1x qz qz q1x q2y q2y
	//4) q10x qy qx qx qy q10x -> q1x qy q1x q2x qy q2x    (vanishing, but we keep this here)
	if(all_flavors.size()==3){
		//1) q10x qx qz qz qx q10x->q1x q1x qz qz q2x q2x
		if((*q1).flavor()%100==(*q2).flavor()&&
			(*q3).flavor()==(*q4).flavor()&&
			(*q5).flavor()==(*q6).flavor()%100){

			temp[q1_pos]=particle_ID((*q1).type(), (*q1).helicity(),10+(*q1).flavor()%100,(*q1).is_anti_particle());
			temp[q2_pos]=particle_ID((*q2).type(), (*q2).helicity(),10+(*q2).flavor(),    (*q2).is_anti_particle());
			temp[q5_pos]=particle_ID((*q5).type(), (*q5).helicity(),20+(*q5).flavor(),    (*q5).is_anti_particle());
			temp[q6_pos]=particle_ID((*q6).type(), (*q6).helicity(),20+(*q6).flavor()%100,(*q6).is_anti_particle());
		}
		//2) q10x qx qx qz qz q10x->q1x q1x q2y qz qz q2y
		if((*q1).flavor()%100==(*q2).flavor()&&
			(*q4).flavor()==(*q5).flavor()&&
			(*q3).flavor()==(*q6).flavor()%100){

			temp[q1_pos]=particle_ID((*q1).type(), (*q1).helicity(),10+(*q1).flavor()%100,(*q1).is_anti_particle());
			temp[q2_pos]=particle_ID((*q2).type(), (*q2).helicity(),10+(*q2).flavor(),    (*q2).is_anti_particle());
			temp[q3_pos]=particle_ID((*q3).type(), (*q3).helicity(),20+(*q3).flavor(),    (*q3).is_anti_particle());
			temp[q6_pos]=particle_ID((*q6).type(), (*q6).helicity(),20+(*q6).flavor()%100,(*q6).is_anti_particle());
		}
		//3) q10x  qz qz qx qx q10x->q1x qz qz q1x q2y q2y
		if((*q1).flavor()%100==(*q4).flavor()&&
			(*q2).flavor()==(*q3).flavor()&&
			(*q5).flavor()==(*q6).flavor()%100){

			temp[q1_pos]=particle_ID((*q1).type(), (*q1).helicity(),10+(*q1).flavor()%100,(*q1).is_anti_particle());
			temp[q4_pos]=particle_ID((*q4).type(), (*q4).helicity(),10+(*q4).flavor(),    (*q4).is_anti_particle());
			temp[q5_pos]=particle_ID((*q5).type(), (*q5).helicity(),20+(*q5).flavor(),    (*q5).is_anti_particle());
			temp[q6_pos]=particle_ID((*q6).type(), (*q6).helicity(),20+(*q6).flavor()%100,(*q6).is_anti_particle());
		}
		//4) q10x qy qx qx qy q10x -> q1x qy q1x q2x qy q2x    (vanishing, but we keep this here)
		//unfinished
		if((*q1).flavor()%100==(*q3).flavor()&&
			(*q2).flavor()==(*q5).flavor()&&
			(*q4).flavor()==(*q6).flavor()%100){

			temp[q1_pos]=particle_ID((*q1).type(), (*q1).helicity(),10+(*q1).flavor()%100,(*q1).is_anti_particle());
			temp[q3_pos]=particle_ID((*q3).type(), (*q3).helicity(),10+(*q3).flavor(),    (*q3).is_anti_particle());
			temp[q4_pos]=particle_ID((*q4).type(), (*q4).helicity(),20+(*q4).flavor(),    (*q4).is_anti_particle());
			temp[q6_pos]=particle_ID((*q6).type(), (*q6).helicity(),20+(*q6).flavor()%100,(*q6).is_anti_particle());
		}



	}
	//1) qx qy qy qz qz q10x->q1x qy qy qz qz q1x
	else if(all_flavors.size()==4){
		//1) qx qy qy qz qz q10x->q1x qy qy qz qz q1x
		if((*q1).flavor()==(*q6).flavor()%100&&
			(((*q2).flavor()==(*q3).flavor()&&(*q4).flavor()==(*q5).flavor())
			 ||((*q2).flavor()==(*q3).flavor()&&(*q4).flavor()==(*q5).flavor()))){

			temp[q1_pos]=particle_ID((*q1).type(), (*q1).helicity(),10+(*q1).flavor(),(*q1).is_anti_particle());
			temp[q6_pos]=particle_ID((*q6).type(), (*q6).helicity(),10+(*q6).flavor()%100,(*q6).is_anti_particle());
		}
	}
	//1) q10x qx qz qz qy q10y->q1x q1x qz qz q2y q2y
	else if(all_flavors.size()==5){
		if((*q1).flavor()%100==(*q2).flavor()&&
			(*q3).flavor()==(*q4).flavor()&&
			(*q5).flavor()==(*q6).flavor()%100){


			temp[q1_pos]=particle_ID((*q1).type(), (*q1).helicity(),10+(*q1).flavor()%100,(*q1).is_anti_particle());
			temp[q2_pos]=particle_ID((*q2).type(), (*q2).helicity(),10+(*q2).flavor(),    (*q2).is_anti_particle());
			temp[q5_pos]=particle_ID((*q5).type(), (*q5).helicity(),10+(*q5).flavor(),    (*q5).is_anti_particle());
			temp[q6_pos]=particle_ID((*q6).type(), (*q6).helicity(),10+(*q6).flavor()%100,(*q6).is_anti_particle());
		}
	}
	return process(temp);
}





process fix_flavors_2q(const process& pro){
	vector<particle_ID> temp(pro.begin(),pro.end());

	vector<int> all_flavors=find_all_flavors(pro,quark,quark_massive);

	if (all_flavors.size() != 2 ) return fix_flavors_gluons(pro);


	int large_flavor,small_flavor;

	if (all_flavors[0] > 100){
		large_flavor=all_flavors[0];
		small_flavor=all_flavors[1];
	} else if (all_flavors[1] > 100){
		large_flavor=all_flavors[1];
		small_flavor=all_flavors[0];
	} else return fix_flavors_gluons(pro);

	if ( (large_flavor % 100) != small_flavor) return fix_flavors_gluons(pro);

	process::const_iterator first_quark=find_if(pro.begin(),pro.end(),is_of_either_type(quark,quark_massive));
	process::const_iterator second_quark=find_if(first_quark+1,pro.end(),is_of_either_type(quark,quark_massive));

// 0 based:
	size_t q1_pos=first_quark-pro.begin();
	size_t q2_pos=second_quark-pro.begin();


	temp[q1_pos]=particle_ID((*first_quark).type(), (*first_quark ).helicity(),10+small_flavor,(*first_quark ).is_anti_particle());
	temp[q2_pos]=particle_ID((*second_quark).type(), (*second_quark).helicity(),10+small_flavor,(*second_quark).is_anti_particle());


return fix_flavors_gluons(process(temp));

}

process fix_flavors_2G(const process& pro){
	vector<particle_ID> temp;
	copy(pro.begin(),pro.end(),back_inserter(temp));

	vector<int> all_flavors=find_all_flavors(pro,gluino,gluino_massive);

	if (all_flavors.size() != 2 ) return pro;

	int large_flavor,small_flavor;

	if (all_flavors[0] > 100){
		large_flavor=all_flavors[0];
		small_flavor=all_flavors[1];
	} else if (all_flavors[1] > 100){
		large_flavor=all_flavors[1];
		small_flavor=all_flavors[0];
	} else return pro;

	if ( (large_flavor % 100) != small_flavor) return pro;

	process::const_iterator first_gluino=find_if(pro.begin(),pro.end(),is_of_either_type(gluino,gluino_massive));
	process::const_iterator second_gluino=find_if(first_gluino+1,pro.end(),is_of_either_type(gluino,gluino_massive));

// 0 based:
	size_t q1_pos=first_gluino-pro.begin();
	size_t q2_pos=second_gluino-pro.begin();


	temp[q1_pos]=particle_ID((*first_gluino).type(), (*first_gluino ).helicity(),10+small_flavor,(*first_gluino ).is_anti_particle());
	temp[q2_pos]=particle_ID((*second_gluino).type(), (*second_gluino).helicity(),10+small_flavor,(*second_gluino).is_anti_particle());

return process(temp);

}

// converts a process with gluons with negative flavor ("flavor changing" gluons) into the right kind of process, with the right fermion lines connected
process fix_flavors(const process& pro){

	// Use a reduced pcode form to find out the particle content
	int ordering_code=count_ph(pro,lepton)*100
					+(count_ph(pro,gluino)+count_ph(pro,gluino_massive))*10
					+(count_ph(pro,quark)+count_ph(pro,quark_massive));

	switch(ordering_code){
	case 2:
		return fix_flavors_2q(pro);
	case 4:
		return fix_flavors_4q(pro);

	case 20:
		return fix_flavors_2G(pro);
	case 40:
		return fix_flavors_4G(pro);
	case 60:
		return fix_flavors_6G(pro);
	case 22:
		return fix_flavors_2G(fix_flavors_2q(pro));
	case 24:
		return fix_flavors_2G(fix_flavors_4q(pro));
	case 42:
		return fix_flavors_4G(fix_flavors_2q(pro));
	case 44:
		return fix_flavors_4G(fix_flavors_4q(pro));
	case 62:
		return fix_flavors_6G(fix_flavors_2q(pro));
	case 202:
		return fix_flavors_2q2e(fix_flavors_2q(pro));
	case 204:
		return fix_flavors_4q2e(fix_flavors_4q(pro));
	case 220:
		return fix_flavors_2q2e(fix_flavors_2G(pro));
	case 240:
		return fix_flavors_4q2e(fix_flavors_4G(pro));
	case 222:
		return fix_flavors_2q2e(fix_flavors_2q(fix_flavors_2G(pro)));
	case 224:
		return fix_flavors_4q2e(fix_flavors_4q(fix_flavors_2G(pro)));
	case 242:
		return fix_flavors_2q2e(fix_flavors_2q(fix_flavors_4G(pro)));
	case 244:
		return fix_flavors_4q2e(fix_flavors_4q(fix_flavors_4G(pro)));
	}

	return fix_flavors_gluons(pro);
}


std::string string_gen(const std::vector<particle_ID>& ps,rule_list& rules){
	std::string res;
	for (int i=0;i<ps.size();i++){
		rule_list::iterator pos;
		for (pos=rules.begin();pos!=rules.end();pos++){
			if ( (pos->first)->match(ps[i]) ) { res += pos->second;}
		};
	}
	return res;
}


std::string string_gen(const process& pro,rule_list& rules){
//	std::string res;
//	for (int i=1;i<=pro.n();i++){
//		rule_list::iterator pos;
//		for (pos=rules.begin();pos!=rules.end();pos++){
//			if ( (pos->first)->match(pro.p(i)) ) { res += pos->second;}
//		};
//	}
	return string_gen(pro.particle_IDs(),rules);
}


color_structure right_direction(const process& pro,color_structure cs){

	particle_type_match is_gluon(gluon);
	particle_type_match is_photon(photon);
	type_and_pap_match is_quark(quark,false);
	type_and_pap_match is_antiquark(quark,true);

	// leptons don't count for the ordering
	// particle_type_match is_lepton(lepton);

	rule R1(&is_quark,"q");
	rule R2(&is_antiquark,"qb");
	rule R3(&is_photon,"y");

	vector<rule> rules;
	rules.push_back(R1);
	rules.push_back(R2);
	rules.push_back(R3);

	std::string type=string_gen(pro,rules);
	type+=type;
	std::string::size_type loc;

	switch (cs) {
	case leading_color: case sub_leading_color : {
		loc = type.find("qbyq",0);
		if (loc != std::string::npos){
			return LT;
		}
		loc = type.find("qyqb",0);
		if (loc != std::string::npos){
			return RT;
		}

		loc = type.find("qbq",0);
		if (loc != std::string::npos){
			return LT;
		}
		loc = type.find("qqb",0);
		if (loc != std::string::npos){
			return RT;
		}
		// neither qbq nor qqb
		throw  BHerror("Process inconsistency");
		}
	case nf: return nf;
	default : throw BHerror("Unhandled color structure");
	}
}



} /* BH */


