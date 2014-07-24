/*
 * recrational_factory.h
 *
 *  Created on: 22-Apr-2009
 *      Author: daniel
 */

#ifndef RECRATIONAL_FACTORY_H_
#define RECRATIONAL_FACTORY_H_

#include "BH_typedefs.h"
#include "rational_factory.h"


namespace BH {

class process;
class option;



//! factory class for recursive rational terms
class Rec_Rational_factory : public Rational_factory<Rational_base> {
public :
	//! new rational with an option
//	virtual Rational_base* new_rational(const process& pro, option* opt);
	//! new rational with option and propagator content specified
	virtual Rational_base* new_rational(const process& pro,const std::vector<particle_ID>& pp,option* opt);
	//! new rational for process pro and color structure cs
	virtual Rational_base* new_rational(const process& pro,color_structure cs);
	static Rational_factory<Rational_base>* s_default_RRF;

};



}


#endif /* RECRATIONAL_FACTORY_H_ */
