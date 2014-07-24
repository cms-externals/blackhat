/*
 * rat_options.h
 *
 *  Created on: 14-Jul-2008
 *      Author: Darren
 */

#ifndef RAT_OPTIONS_H_
#define RAT_OPTIONS_H_

#include "options.h"

namespace BH {

//! option for cutD with closed quark loops
/** returns true if the cutD object has a closed quark loop */
class cutD_has_massive_quark_loop : public option {
private:
	int _quark_flavor;
	public:
		cutD_has_massive_quark_loop(int fl =0 ): _quark_flavor(fl ){};
		virtual bool operator()(cutD*);
		virtual option* clone() const ;
		virtual ~cutD_has_massive_quark_loop(){};
};

//! option for cutD with an external fermion turning right into the loop
/** a cut diagram is called "right turner" if the fermion line starting at the antifermion enters the loopand turns right */
class  is_massive_right_turner : public option {
private:
	size_t _quark_flavor;
public:
	is_massive_right_turner(size_t f=0) : _quark_flavor(f) {};
	virtual bool operator()(cutD*);
	virtual option* clone() const;
	virtual ~is_massive_right_turner(){};

};

//! option for cutD with an external fermion turning left into the loop
/** a cut diagram is called "left turner" if the fermion line starting at the antifermion enters the loop and turns left */
class  is_massive_left_turner : public option {
private:
	size_t _quark_flavor;
public:
	is_massive_left_turner(size_t f=0) : _quark_flavor(f) {};
	virtual bool operator()(cutD*);
	virtual option* clone() const;
	virtual ~is_massive_left_turner(){};

};

//! option for cutD with closed gluon loops
/** returns true if the cutD object has a closed gluon loop */
class cutD_has_massive_scalar_loop : public option {
	public:
		virtual bool operator()(cutD*);
		virtual option* clone() const;
		virtual ~cutD_has_massive_scalar_loop(){};
};

//! option for cutD with closed gluino loops
/** returns true if the cutD object has a closed gluino loop */
class cutD_has_massive_gluino_loop : public option {
private:
	int _gluino_flavor;
	public:
		cutD_has_massive_gluino_loop(int fl =0 ): _gluino_flavor(fl ){};
		virtual bool operator()(cutD*);
		virtual option* clone() const ;
		virtual ~cutD_has_massive_gluino_loop(){};
};

//! option for cutD with a closed massive quark loop that is not faked by a pinch
/** returns true if the cutD object has a closed massive quark loop */
class cutD_has_factorized_massive_quark_loop : public option {
private:
	int _quark_flavor;
	public:
		cutD_has_factorized_massive_quark_loop(int fl =0 ): _quark_flavor(fl ){};
		virtual bool operator()(cutD*);
		virtual option* clone() const ;
		virtual ~cutD_has_factorized_massive_quark_loop(){};
};

//! option for cutD with closed scalar loops
/** returns true if the cutD object has a closed massive loop with on external vertex containing only gluons from a specified set of indices*/
class cutD_has_massive_scalar_loop_with_external_gluons : public option {
	std::vector<int> _ind;
public:
		cutD_has_massive_scalar_loop_with_external_gluons(const std::vector<int>& ind): _ind(ind){ /*needs a sorted vector*/ sort(_ind.begin(),_ind.end());};
		virtual bool operator()(cutD*);
		virtual option* clone() const;
		virtual ~cutD_has_massive_scalar_loop_with_external_gluons(){};
};

//! option for cutD with a closed massive gluino loop that is not faked by a pinch
/** returns true if the cutD object has a closed massive gluino loop */
class cutD_has_factorized_massive_gluino_loop : public option {
private:
	int _gluino_flavor;
	public:
		cutD_has_factorized_massive_gluino_loop(int fl =0 ): _gluino_flavor(fl ){};
		virtual bool operator()(cutD*);
		virtual option* clone() const ;
		virtual ~cutD_has_factorized_massive_gluino_loop(){};
};

}


#endif /* RAT_OPTIONS_H_ */
