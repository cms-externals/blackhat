/*
 * GMP_random.h
 *
 *  Created on: 3 Mar 2011
 *      Author: daniel
 */

#ifndef GMP_RANDOM_H_
#define GMP_RANDOM_H_

#include "gmp_r.h"

class GMP_random {
	gmp_randstate_t d_state;
public:
	GMP_random();
	void seed(unsigned long int seed);
	RGMP random();

};


#endif /* GMP_RANDOM_H_ */
