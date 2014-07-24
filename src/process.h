/*
 * process.h
 *
 *  Created on: 11-Nov-2008
 *      Author: daniel
 */

#ifndef PROCESS_H_
#define PROCESS_H_

#include <vector>
#include "particles.h"
#include "counted.h"
#include "BH_error.h"

namespace BH {

class process  //: private counted<process>
{
	size_t nbr;					/*!< number of particles in the process  */
	std::vector<particle_ID> particles;	/*!< vector containing the ptype of the particles  */
	long d_pcode;
	friend std::ostream& operator<<(std::ostream&, const process&);
public:
	/*! returns number of particles in the process  */
	size_t n() const {return nbr;};
	/*! \param n label of the particle
			    \return ptype of the nth particle
	*/

	typedef std::vector<particle_ID>::iterator iterator;
	typedef std::vector<particle_ID>::const_iterator const_iterator;
	const particle_ID& p(size_t n) const  {if (n<=nbr&&n>0) return particles[n-1]; else {std::cerr<< "Too large particle index in process::p with n=" << n << " for process="<< *this << std::endl;throw BHerror("Overflow in class process"); }};
	const std::vector<particle_ID>& particle_IDs() const {return particles;};
	process() : nbr(0) , particles(0), d_pcode(0) {};
	process(const std::vector<plabel>& pl);
	process(const std::vector<particle_ID>& pl);
	process(std::vector<particle_ID>& pl,long pcode); // this version takes over the vector of particle_IDs
	process(const std::vector<particle*>&,std::vector<short>);
	explicit process(particle_ID);
	std::vector<particle_ID>::const_iterator begin() const {return particles.begin();};
	std::vector<particle_ID>::const_iterator end() const {return particles.end();};
	size_t size() const {return nbr;};
	long pcode() const {return d_pcode;}
	const particle_ID& operator[](size_t n) const {return particles[n];}; // zero-based
	const particle_ID& front() const {return particles[0];};
	const particle_ID& back() const {return particles[nbr-1];};
	std::string print() const;
	process(particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID);
	process(particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID, particle_ID);
};

bool operator==(const process& p1, const process& p2);
bool operator<(const process& p1, const process& p2);

bool Tree_is_zero(const process& pro);

}

#endif /* PROCESS_H_ */
