/*!\file particles.h
\brief Header for particles: ptype, plabel, ...
*/
#ifndef PARTICLES_H_
#define PARTICLES_H_

#include<iostream>
#include<vector>
#include "BH_typedefs.h"

namespace BH {




class plabel;
class process;
class particle_ID;
class mass_param;


std::ostream& operator<<(std::ostream&, const plabel& );
std::ostream& operator<<(std::ostream&, const process& );
//! class for particles
/** A particle object contains information about its statistics type (fermion, boson), itss mass and its name  */
class particle {
	friend class particle_ID;
	friend bool operator<(const particle& p1, const particle& p2);
	friend bool operator>(const particle& p1, const particle& p2);
	friend long scode(const process& p);
public:
	enum statistics_type {fermion, boson};
	// The different ways a unordered particle can couple to the ordered particles,
	//  does it couple to quarks only, qc_type, gluons only, gc_type, or both qgc_type
	//  a value of 0 corresponds to ordered.
	enum unordered_type {ord_type=0, qc_type=1, gc_type=2, qgc_type=2};
//! constuctor
	particle(statistics_type stat,std::string name,double mass,int pdg,int ordered=particle::ord_type,int mass_index=0);
	//! mass of the particle
	/** the return type is double, use T(m()) for other precision. \return mass of the particle*/
	double m() const {return _mass;};
	//! name of the particle
	const std::string& name() const {return _name;};
//! statistcs type
	statistics_type stat() const {return _stat;};
	//! returns the PDG code of the particle
	int pdg_code() const {return _pdg_code;};
	int get_ordered() const {return _ordered;};
	int get_mass_index() const {return d_mass_label;};
	void set_mass_index(int new_mass) {d_mass_label=new_mass;};
	bool ordered(int level) const {return _ordered==level;};
private:
	statistics_type _stat;
	double _mass;
	std::string _name;
	size_t _particle_nbr;
	static size_t _next_particle_nbr;
	size_t particle_nbr() const {return _particle_nbr;};
	int _pdg_code;
	int _ordered;
	int d_mass_label;
};

bool operator==(const particle& p1, const particle& p2);
bool operator!=(const particle& p1, const particle& p2);
bool operator<(const particle& p1, const particle& p2);
bool operator>(const particle& p1, const particle& p2);
//! Class for the quantum numbers of particles.
/** A particle_ID object contains the type, the helicity, the flavor and the particle/anti-particle nature of a particle. */
class particle_ID {
public:
	//! default constructor
	particle_ID(){};
	//! constructor
	particle_ID(particle& type,short helicity=0,short flavor=1,bool is_anti_particle=false);
	//! constructor
	particle_ID(particle* type,short helicity=0,short flavor=1,bool is_anti_particle=false);
	//! name of the particle
	const std::string& name() const {return _type->name();};
	//! type of the particle
	/** \return pointer to the particle objects that represents the type of the particle. */
	particle* type() const  {return _type;};
	//! helicity of the particle
	short helicity() const {return _helicity;};
	//! flavor of the particle
	short flavor() const {return _flavor;};
	//! particle/anti-particle nature
	/** \return true if the particle is an antiparticle, false if it is a particle */
	bool is_anti_particle() const {return _is_anti_particle;};
	//! comparaison with a particle
	/** \param p reference to a particle object \return true if the particle_ID has the type p provided as argument */
	bool is_a(const particle& p) const ;
	//! comparaison with a particle_ID
	/** \param p reference to a particle_ID object \return true if the particle_ID is the same as the one provided as argument */
	bool is_a(const particle_ID& p) const;
	//! comparaison with a particle
	/** \param p reference to a particle object \return true if the particle_ID has the type p provided as argument */
	bool is_not_a(const particle& p) const;
	//! comparaison with a particle_ID
	/** \param p reference to a particle_ID object \return true if the particle_ID is not the same as the one provided as argument */
	bool is_not_a(const particle_ID& p) const;
	/**\return true if the particle is a fermion */
	bool is_fermion() const {return _type->stat()==particle::fermion;};
	/**\return a label that is 0 for massless particles and a different non vanishing integer for each different massive particle type. */
	int get_ordered() const {return _type->get_ordered();};
	bool ordered(int level) const {return _type->ordered(level);};
	int mass_label() const {return _type->get_mass_index();};
//	void set_mass_label(size_t mass) {_mass=mass;}
	//! particle mass
	/** \retun the mass of the particle */
	double mass() const {return _type->m();};
	//! particle mass squared
	/** \retun the mass squared of the particle */
	double mass2() const {return _type->m()*_type->m();};
	//! particle conjugation
	/** \return particle_ID of the same type, with helicity flipped, and particle/anti-particle flipped if the particle is a fermion. \sa helicity_conjugate()  */
	particle_ID conjugate() const ;
	//! particle helicity flip
	/** \return particle_ID of the same type, same flavor with helicity flipped, but particle/anti-particle flipped if the particle is a fermion. \sa conjugate()  */
	particle_ID helicity_conjugate() const ;

	void set_helicity(int h){_helicity=h;};
	void set_flavor(int fl){_flavor=fl;};
	void set_is_antiparticle(bool iap){ _is_anti_particle=iap;};

private:
	particle*  _type;
	short _helicity;
	short _flavor;
	bool _is_anti_particle;
};
bool operator==(const particle_ID& p1, const particle_ID& p2);
bool operator!=(const particle_ID& p1, const particle_ID& p2);
bool operator<(const particle_ID& p1, const particle_ID& p2);
bool operator>(const particle_ID& p1, const particle_ID& p2);

std::ostream& operator<<(std::ostream& s, const particle_ID& r);

typedef particle_ID ph_type;

extern particle gluon;
extern particle quark;
extern particle lepton;
extern particle scalar;
extern particle photon;
extern particle gluino;
extern particle higgs;

// Particles needed for the rational extraction code
extern particle gluon_massive_scalar;
extern particle gluon_massive;
extern particle quark_massive;
extern particle gluino_massive;
extern particle scalar_massive;

extern particle_ID m;
extern particle_ID p;
extern particle_ID qp;
extern particle_ID qm;
extern particle_ID q2p;
extern particle_ID q2m;
extern particle_ID q3p;
extern particle_ID q3m;
extern particle_ID q4p;
extern particle_ID q4m;
extern particle_ID qbp;
extern particle_ID qbm;
extern particle_ID qb2p;
extern particle_ID qb2m;
extern particle_ID qb3p;
extern particle_ID qb3m;
extern particle_ID qb4p;
extern particle_ID qb4m;
extern particle_ID lp;
extern particle_ID lm;
extern particle_ID lbp;
extern particle_ID lbm;
extern particle_ID l2p;
extern particle_ID l2m;
extern particle_ID lb2p;
extern particle_ID lb2m;
extern particle_ID s0;
extern particle_ID yp;
extern particle_ID ym;
extern particle_ID y2p;
extern particle_ID y2m;
extern particle_ID Gp;
extern particle_ID Gm;
extern particle_ID Gbp;
extern particle_ID Gbm;
extern particle_ID G2p;
extern particle_ID G2m;
extern particle_ID Gb2p;
extern particle_ID Gb2m;
extern particle_ID G3p;
extern particle_ID G3m;
extern particle_ID Gb3p;
extern particle_ID Gb3m;

// Particles needed for the rational extraction code
extern particle_ID gsc;
extern particle_ID gpM;
extern particle_ID gmM;
extern particle_ID gD;
extern particle_ID scM;
extern particle_ID Qp;
extern particle_ID Qm;
extern particle_ID Qbp;
extern particle_ID Qbm;
extern particle_ID Q2p;
extern particle_ID Q2m;
extern particle_ID Qb2p;
extern particle_ID Qb2m;
extern particle_ID Lp;
extern particle_ID Lm;
extern particle_ID Lbp;
extern particle_ID Lbm;
extern particle_ID L2p;
extern particle_ID L2m;
extern particle_ID Lb2p;
extern particle_ID Lb2m;



// Higgs particles
extern particle_ID ph; //complex scalar field component of H
extern particle_ID phd; //complex scalar field component of H
extern particle_ID H; //complex scalar field component of H

//The mass_params to define the massive particles
extern mass_param rat_ext_mass;
extern mass_param top_ext_mass;
extern mass_param top_mass;
extern mass_param higgs_mass;

//! class for particle_ID that carry and index
class plabel : public particle_ID
{
	int label;
	friend std::ostream& operator<<(std::ostream&, const plabel&);

public:
	//! index
	int ind() const {return label;};
	//! default constructor
//	plabel(): particle_ID(no_particle,0), label(-1)  {};
	//! constructor
	plabel(const particle_ID&,int);
	//! constructor
//	plabel(particle_type,short,int);
	//! constructor
//	plabel(ph_type,int);
};

bool operator==(const plabel& p1, const plabel& p2);
bool operator!=(const plabel& p1, const plabel& p2);
bool operator<(const plabel& p1, const plabel& p2);
bool operator>(const plabel& p1, const plabel& p2);
//! process represents a collection of particles with helicities


std::vector<int> Indices(const std::vector<plabel>& p);

}
#endif /*PARTICLES_H_*/
