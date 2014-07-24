/*
 * process_utils.h
 *
 *  Created on: Jul 17, 2008
 *      Author: daniel
 */

#ifndef PROCESS_UTILS_H_
#define PROCESS_UTILS_H_

#include "BH_typedefs.h"
#include "particles.h"
#include "process.h"
#include "iterators.h"
#include <vector>

#include "process_utils_basic.h"

namespace BH {

class process;
class option;

std::vector<size_t> find_indices(const process& pro,const particle& t);
bool are_quarks_adjacent(const process& pro);

class particle_type_match : public particle_pattern {
	particle& _type;
public:
	particle_type_match(particle& type): _type(type) {};
	virtual bool match(const particle_ID& p){ return (*(p.type()) == _type) ;};
};

class particle_type_ordered : public particle_pattern {
public:
	particle_type_ordered() {};
	virtual bool match(const particle_ID& p){ return p.ordered(particle::ord_type);};
};

class particle_type_unordered : public particle_pattern {
	int _level;
public:
	particle_type_unordered(int level) : _level(level ){};
	virtual bool match(const particle_ID& p){ return p.ordered(_level);};
};

class particle_type_any_unordered : public particle_pattern {
public:
	particle_type_any_unordered() {};
	virtual bool match(const particle_ID& p){ return !p.ordered(particle::ord_type);};
};

class particle_or_antiparticle_match : public particle_pattern {
	bool _anti_particle;
public:
	particle_or_antiparticle_match(bool ap): _anti_particle(ap) {};
	virtual bool match(const particle_ID& p){ return (p.is_anti_particle() == _anti_particle) ;};
};

class helicity_match : public particle_pattern {
	int d_helicity;
public:
	helicity_match(int hel): d_helicity(hel) {};
	virtual bool match(const particle_ID& p){ return (p.helicity() == d_helicity) ;};
};


class type_and_pap_match : public particle_pattern {
	bool _anti_particle;
	particle& _type;
public:
	type_and_pap_match(particle& type,bool ap): _type(type) ,_anti_particle(ap) {};
	virtual bool match(const particle_ID& p){ return  (*(p.type()) == _type) && (p.is_anti_particle() == _anti_particle) ;};
};

class type_and_flavor_match : public particle_pattern {
	int m_flavor;
	particle& m_type;
public:
	type_and_flavor_match(particle& type,int fl): m_type(type), m_flavor(fl) {};
	virtual bool match(const particle_ID& p){ return  (*(p.type()) == m_type) && (p.flavor() == m_flavor) ;};
};

class perfect_match : public particle_pattern {
	particle_ID _type;
public:
	perfect_match(particle_ID type): _type(type) {};
	virtual bool match(const particle_ID& p){ return ( p == _type)  ;};
	virtual ~perfect_match(){}
};




process arrange_flavors_2q2e(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G2e(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G2e_q(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G1y_q(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G2e_G(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G1y_G(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2e_SLC(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_4q2e(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G1y(const process& pro,std::vector<particle_ID>& propagators);
process order_2q2e(const process& pro,std::vector<int>& indices);

size_t is_a_rotation_of(const process& pro1,const process& pro2);
size_t helicity_is_a_rotation_of(const process& pro1,const process& pro2);



struct is_of_type : public std::unary_function<particle_ID,bool> {
	particle& m_particle;
	is_of_type(particle& pa): m_particle(pa){};
	bool operator()(particle_ID pa) const {return pa.is_a(m_particle);};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return (*ci).is_a(m_particle);};
};

struct is_of_either_type : public std::unary_function<particle_ID,bool> {
	particle& m_particle1;
	particle& m_particle2;
	is_of_either_type(particle& pa1, particle& pa2): m_particle1(pa1), m_particle2(pa2){};
	bool operator()(particle_ID pa) const {return (pa.is_a(m_particle1)||pa.is_a(m_particle2));};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return ((*ci).is_a(m_particle1)||(*ci).is_a(m_particle2));};
};

struct is_of_any_type : public std::unary_function<particle_ID,bool> {
	particle& m_particle1;
	particle& m_particle2;
	particle& m_particle3;
	is_of_any_type(particle& pa1, particle& pa2, particle& pa3): m_particle1(pa1), m_particle2(pa2), m_particle3(pa3){};
	bool operator()(particle_ID pa) const {return (pa.is_a(m_particle1)||pa.is_a(m_particle2)||pa.is_a(m_particle3));};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return ((*ci).is_a(m_particle1)||(*ci).is_a(m_particle2)||(*ci).is_a(m_particle3));};
};

struct is_of_type_and_pap : public std::unary_function<particle_ID,bool> {
	particle& m_particle;
	bool m_pap;
	is_of_type_and_pap(particle& pa,bool pap): m_particle(pa), m_pap(pap){};
	bool operator()(particle_ID pa) const {return ( pa.is_a(m_particle) && (pa.is_anti_particle() == m_pap ));};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return ( (*ci).is_a(m_particle) && ((*ci).is_anti_particle() == m_pap ));};
};

struct is_of_type_and_flavor : public std::unary_function<particle_ID,bool> {
	particle& m_particle;
	int m_flavor;
	is_of_type_and_flavor(particle& pa,int fl): m_particle(pa), m_flavor(fl){};
	bool operator()(particle_ID pa) const {return ( pa.is_a(m_particle) && (pa.flavor() == m_flavor ));};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return ( (*ci).is_a(m_particle) && ((*ci).flavor() == m_flavor ));};
};


struct is_of_type_and_flavor_mod : public std::unary_function<particle_ID,bool> {
	particle& m_particle;
	int m_flavor;
	int m_mod;
	is_of_type_and_flavor_mod(particle& pa,int fl, int mod): m_particle(pa), m_flavor(fl), m_mod(mod){};
	bool operator()(particle_ID pa) const {return ( pa.is_a(m_particle) && (pa.flavor()%m_mod == m_flavor%m_mod ));};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return ( (*ci).is_a(m_particle) && ((*ci).flavor() == m_flavor ));};
};


struct is_unordered : public std::unary_function<particle_ID,bool> {
	int _level;
	is_unordered(int level): _level(level){};
	bool operator()(particle_ID pa) const {return pa.ordered(_level);};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return (*ci).ordered(_level);};
};

struct is_any_unordered : public std::unary_function<particle_ID,bool> {
	bool operator()(particle_ID pa) const {return !pa.ordered(particle::ord_type);};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return !(*ci).ordered(particle::ord_type);};
};

struct is_ordered : public std::unary_function<particle_ID,bool> {
	bool operator()(particle_ID pa) const {return pa.ordered(particle::ord_type);};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return (*ci).ordered(particle::ord_type);};
};

struct has_positive_helicity : public std::unary_function<particle_ID,bool> {
	bool operator()(particle_ID pa) const {return  pa.helicity() > 0 ;};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return  (*ci).helicity() >0;};
};

struct has_negative_helicity : public std::unary_function<particle_ID,bool> {
	bool operator()(particle_ID pa) const {return pa.helicity() < 0 ;};
	bool operator()(cyclic_iterator<particle_ID,process > ci) const {return  (*ci).helicity() < 0;};
};

struct has_negative_flavor : public std::unary_function<particle_ID,bool> {
	bool operator()(particle_ID pa){return pa.flavor() < 0 ;};
};

struct binary_is_of_type : public std::binary_function<particle_ID,int,bool> {
	particle& m_particle;
	binary_is_of_type(particle& pa): m_particle(pa){};
	bool operator()(particle_ID pa,int dummy){return pa.is_a(m_particle);};
};

// sets possible_props to the right set of possible propagators and returns a suitable option to construct the tree of cutDs for the specified process and the specisied color structure and changes the process if necessary
option* get_possible_props_and_options(const process& pro,color_structure cs,std::vector<particle_ID>& possible_props,process& new_pro);



}
#endif /* PROCESS_UTILS_H_ */
