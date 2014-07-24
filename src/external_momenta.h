/*
 * external_momenta.h
 *
 *  Created on: Jan 8, 2009
 *      Author: daniel
 */

#ifndef EXTERNAL_MOMENTA_H_
#define EXTERNAL_MOMENTA_H_

#include <vector>

namespace BH {

template <class T> class momentum_configuration;

class External_Momenta_factory {
	long d_last_mc_ID;
	int d_last_length;
	int d_length;
	std::vector<std::vector<int> > d_index_vectors;
	std::vector<long> d_indices_codes;
	std::vector<int> d_indices;
	std::vector<int> d_load;
public:
	External_Momenta_factory() : d_last_mc_ID(1), d_last_length(0), d_length(0) {};
	int new_external_momentum(const std::vector<int>& ind);
	int get_momentum_index(int index){return d_indices[index];};
	// prepares the momenta in the storage of the factory. nbr_mom specifies how many momenta satisfy the momentum consrvation p_1+...+p_n=0
	void prepare(momentum_configuration<double>& mc,int nbr_mom);
	void print_state();
	static External_Momenta_factory s_global_external_factory;
};


}

#endif /* EXTERNAL_MOMENTA_H_ */
