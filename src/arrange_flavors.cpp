/*
 * arrange_flavors.cpp
 *
 *  Created on: 22-Apr-2009
 *      Author: daniel
 */

#include "arrange_flavors.h"
#include "particles.h"
#include "amplitudes.h"
#include <vector>
#include <memory>
#include <algorithm>
#include "process_utils.h"
#include "iterators.h"

#define _VERBOSE 0

using namespace std;

namespace BH {
process fix_flavors(const process& pro);



process arrange_flavors(const process& pro,size_t qb_pos,size_t q_pos,short increment,vector<particle_ID>& propagators){
//	_PRINT(pro);
	vector<particle_ID> temp;
	for (int i=1;i<=pro.n();i++){
		temp.push_back(pro.p(i));
	}
	size_t distance;
//	_PRINT(q_pos);_PRINT(qb_pos);
	switch (increment){
	case 1 :{
		distance=(pro.n()+q_pos-qb_pos)%pro.n(); break;
	}
	case -1:	{
		distance=(pro.n()+qb_pos-q_pos)%pro.n();break;
	}
	}
//	_PRINT(distance);
	for (size_t j=1;j<distance;j++){
		temp[(qb_pos+pro.n()+j*increment-1)%pro.n()]=particle_ID(*pro.p((qb_pos+pro.n()+j*increment-1)%pro.n()+1).type(),pro.p((qb_pos+pro.n()+j*increment-1)%pro.n()+1).helicity(),-j,pro.p((qb_pos+pro.n()+j*increment-1)%pro.n()+1).is_anti_particle()?true:false);
//		temp[(qb_pos+pro.n()+j*increment-1)%pro.n()]=particle_ID(gluon,pro.p((qb_pos+pro.n()+j*increment-1)%pro.n()+1).helicity(),-j,false);
		propagators.push_back(particle_ID(quark,1,j,true));
		propagators.push_back(particle_ID(quark,-1,j,true));
		propagators.push_back(particle_ID(quark,1,j,false));
		propagators.push_back(particle_ID(quark,-1,j,false));
	}
	propagators.push_back(particle_ID(quark,1,distance,true));
	propagators.push_back(particle_ID(quark,-1,distance,true));
	propagators.push_back(particle_ID(quark,1,distance,false));
	propagators.push_back(particle_ID(quark,-1,distance,false));
	temp[qb_pos-1]=particle_ID(quark,pro.p(qb_pos).helicity(),1,true);
	temp[q_pos-1]=particle_ID(quark,pro.p(q_pos).helicity(),distance,false);
//	_PRINT(process(temp));
	return process(temp);
}

vector<int> get_unordered_gluons_2q1y(const process& PRO){

	process::const_iterator ph=find_if(PRO.begin(),PRO.end(),is_of_type(photon));
	cyclic_iterator<particle_ID,process > c1(PRO,4,ph);

	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > first_q(c1);

	cyclic_iterator<particle_ID,process > c_pos=first_q;

	vector<int> gluons;

	while ( c_pos != second_q){
		if ( (*c_pos).is_a(gluon) ) {
			gluons.push_back(c_pos.position());
		}
		++c_pos;
	};


	return gluons;

}
vector<int> get_unordered_particles_2q2G1y(const process& PRO){

	process::const_iterator ph=find_if(PRO.begin(),PRO.end(),is_of_type(photon));
	cyclic_iterator<particle_ID,process > c1(PRO,4,ph);

	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > first_q(c1);

	cyclic_iterator<particle_ID,process > c_pos=first_q;

	vector<int> gluons;

	while ( ++c_pos != second_q){
			gluons.push_back(c_pos.position());
	};


	return gluons;

}



vector<int> get_unordered_gluons_2q2e(const process& PRO){

	process::const_iterator lep=find_if(PRO.begin(),PRO.end(),is_of_type(lepton));
	cyclic_iterator<particle_ID,process > c1(PRO,4,lep);

	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > second_q(c1);
	while ( !(*(++c1)).is_a(quark) ){};
	cyclic_iterator<particle_ID,process > first_q(c1);

	cyclic_iterator<particle_ID,process > c_pos=first_q;

	vector<int> gluons;

	while ( c_pos != second_q){
		if ( (*c_pos).is_a(gluon) ) {
			gluons.push_back(c_pos.position());
		}
		++c_pos;
	};


	return gluons;

}


vector<int> get_gluons_between_gluinos_y(const process& PRO){

	process::const_iterator ph=find_if(PRO.begin(),PRO.end(),is_of_type(photon));
	cyclic_iterator<particle_ID,process > c1(PRO,4,ph);

	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > first_G(c1);
	while ( !(*(++c1)).is_a(gluino) ){};
	cyclic_iterator<particle_ID,process > second_G(c1);

	cyclic_iterator<particle_ID,process > c_pos=first_G;

	vector<int> gluons;

	while ( c_pos != second_G){
		if ( (*c_pos).is_a(gluon) ) {
			gluons.push_back(c_pos.position());
		}
		++c_pos;
	};


	return gluons;

}


}
