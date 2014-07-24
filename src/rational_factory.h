/*
 * rational_factory.h
 *
 *  Created on: 8 Jul 2009
 *      Author: daniel
 */

#ifndef RATIONAL_FACTORY_H_
#define RATIONAL_FACTORY_H_


#include "BH_typedefs.h"
#include "rational.h"

namespace BH {

class process;
class option;

//! abstract class for rational part  factories
template <typename rational_type> class Rational_factory {
public:
	typedef rational_type Rational_Type;
	virtual rational_type* new_rational(const process& pro,color_structure cs)=0;
	virtual rational_type* new_rational(const process& pro,const std::vector<particle_ID>& possible_props,option* opt)=0;
	static Rational_factory<Rational_base>* default_rational_factory();
};

#ifdef SWIG
%template(RationalFactory) Rational_factory<Rational_base>;
#endif




//! factory class for recursive rational terms
class Known_Rational_factory : public Rational_factory<Rational_base> {
public :
	//! new rational with option and propagator content specified
	virtual Rational_base* new_rational(const process& pro,const std::vector<particle_ID>& pp,option* opt){throw BHerror("Not implemented");};
	//! new rational for process pro and color structure cs
	virtual Rational_base* new_rational(const process& pro,color_structure cs);
	static Rational_factory<Rational_base>* s_default_KRF;

};

}

#endif /* RATIONAL_FACTORY_H_ */
