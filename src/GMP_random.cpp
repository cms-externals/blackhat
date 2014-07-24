/*
 * GMP_random.cpp
 *
 *  Created on: 3 Mar 2011
 *      Author: daniel
 */

#include "GMP_random.h"
#include "gmp_r.h"

GMP_random::GMP_random(){
	gmp_randinit_default(d_state );
}

RGMP GMP_random::random(){
	int old_prec=RGMP::get_current_precision();
	RGMP::set_precision(20480);
	mpfr::mpreal x=mpfr::urandom(d_state);
	RGMP::set_precision(old_prec);
	return RGMP(x);

}
void GMP_random::seed(unsigned long int seed ){
	gmp_randseed_ui(d_state,seed);
}
