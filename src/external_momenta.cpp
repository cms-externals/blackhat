/*
 * external_momenta.cpp
 *
 *  Created on: Jan 8, 2009
 *      Author: daniel
 */

#include "external_momenta.h"
#include <vector>
#include "cached_integral.h"
#include "mom_conf.h"
#include <iostream>

using namespace std;

namespace BH {


int External_Momenta_factory::new_external_momentum(const vector<int>& v){
	long s=index_combination_code(v);

	vector<long>::iterator it;
	it=find(d_indices_codes.begin(),d_indices_codes.end(),s);

	if ( it != d_indices_codes.end() ){
		int loc=it-d_indices_codes.begin();
		++d_load[loc];
		return loc;
	} else {
		d_indices_codes.push_back(s);
		d_indices.push_back(0);
		d_load.push_back(1);
		d_index_vectors.push_back(v);
		++d_length;
		return d_indices.size()-1;
	}
}
/** checks whether the mom_conf is still the same and that no new momenta has been added since the last preparation.
 *  only recompute the new momenta if the number of momenta has changed, but not the mc.
 */
void External_Momenta_factory::prepare(mom_conf& mc,int nbr_mom){
	int start;
	if ( d_last_mc_ID == mc.get_ID() ){
			start=d_last_length;
	} else {
		d_last_mc_ID= mc.get_ID();
		start=0;
	}
	for (int i=start;i<d_index_vectors.size();i++){
		int n=d_index_vectors[i].size();
		if ( n > 1 ){
			momentum<C> K;
			for (int k=0;k<d_index_vectors[i].size();k++){
				K+=mc.mom(d_index_vectors[i][k]);
			}
			if (n < nbr_mom-1 ) {
				d_indices[i]=mc.insert(K,_mt_massive);
			} else {
				d_indices[i]=mc.insert(K,_mt_massless);
			}
		} else {
			d_indices[i]=d_index_vectors[i][0];
		}
	}
	d_last_length=d_length;

}

void External_Momenta_factory::print_state(){
	for (int i=0;i<d_length;i++){
		momentum<C> K;
		for (int k=0;k<d_index_vectors[i].size();k++){
			cout << d_index_vectors[i][k];
		}
		cout << " used " << d_load[i] << " times." << endl;
	}
}

External_Momenta_factory External_Momenta_factory::s_global_external_factory;


}
