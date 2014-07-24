/*!\file particles.cpp
\brief Implementation of the classes that deal with the particle types and the process class


*/

#include "particles.h"
#include "BH_typedefs.h"
#include "BH_error.h"
#include "eval_param.h"
#include <algorithm>

using namespace std;

namespace BH {


plabel::plabel(const particle_ID& pt, int i): particle_ID(pt), label(i) {}


bool operator==(const plabel& p1, const plabel& p2){
	if (p1.helicity() != p2.helicity()) return false;
	if (p1.is_anti_particle() != p2.is_anti_particle()) return false;
	if (p1.flavor() != p2.flavor()) return false;
	if (p1.ind() != p2.ind()) return false;
	if (p1.type() != p2.type()) return false;
	return true;
}

bool operator!=(const plabel& p1, const plabel& p2){
	if (p1 == p2) return false; else return true;
}

bool operator<(const plabel& p1, const plabel& p2){
	const particle_ID& rp1=p1;
	const particle_ID& rp2=p2;
	if (p1.ind() < p2.ind()) return true;
	if (p1.ind() > p2.ind()) return false;
	if ( rp1 < rp2 ) return true;
	return false;
}

bool operator>(const plabel& p1, const plabel& p2){
	const particle_ID& rp1=p1;
	const particle_ID& rp2=p2;
	if (p1.ind() > p2.ind()) return true;
	if (p1.ind() < p2.ind()) return false;
	if ( rp1 > rp2 ) return true;
	return false;
}

ostream& operator<<(ostream& s, const plabel& r){
	const particle_ID& rpID=r;
	return s << rpID <<"(" << r.ind() <<")";
}

// makes all flavors 1 and all particles particle and not anti-particles



vector<int> Indices(const vector<plabel>& p)
{vector<int> indices;
  for (size_t j = 0;  j < p.size();  j += 1)
    indices.push_back(p[j].ind());
 return(indices);
}


bool operator==(const particle& p1, const particle& p2){
	if (p1.stat() != p2.stat()) return false;
	if (p1.m() != p2.m()) return false;
	if (p1.pdg_code() != p2.pdg_code()) return false;
// we don't compare the name of the particles, the pdg code should e enough
//	if (p1.name() != p2.name()) return false;
	return true;
}

bool operator!=(const particle& p1, const particle& p2){
	if (p1 == p2) return false; else return true;
}

bool operator<(const particle& p1, const particle& p2){
	if (p1.particle_nbr() < p2.particle_nbr()) return true;
	if (p1.particle_nbr() > p2.particle_nbr()) return false;
	return false;
}

bool operator>(const particle& p1, const particle& p2){
	if (p1.particle_nbr() > p2.particle_nbr()) return true;
	if (p1.particle_nbr() < p2.particle_nbr()) return false;
	return false;
}

size_t particle::_next_particle_nbr=1;

particle::particle(enum statistics_type stat,std::string name, double mass,int pdg,int ordered,int mass_label): _stat(stat),_mass(mass),_name(name) ,_particle_nbr(particle::_next_particle_nbr), _pdg_code(pdg), _ordered(ordered), d_mass_label(mass_label) {particle::_next_particle_nbr++;}

particle_ID::particle_ID(particle& type,short helicity,short flavor,bool is_anti_particle): _type(&type), _helicity(helicity), _flavor(flavor), _is_anti_particle(is_anti_particle) {}
particle_ID::particle_ID(particle* type,short helicity,short flavor,bool is_anti_particle): _type(type), _helicity(helicity), _flavor(flavor), _is_anti_particle(is_anti_particle) {}


bool particle_ID::is_a(const particle& p) const {
	return ( (*_type) == p );
}

bool particle_ID::is_not_a(const particle& p) const {
	return ( (*_type) != p );
}

bool particle_ID::is_a(const particle_ID& p) const {
	return ( (*this) == p );
}

bool particle_ID::is_not_a(const particle_ID& p) const {
	return ( (*this) != p );
}

//size_t particle_ID::mass_label() const {
//	if (mass() == 0.) {
//		return 0;
//	}
//	else {
//		return _type->particle_nbr();
//		}
//}

bool operator==(const particle_ID& p1, const particle_ID& p2){
	if (p1.type() != p2.type()) return false;
	if (p1.helicity() != p2.helicity()) return false;
	if (p1.flavor() != p2.flavor()) return false;
	if (p1.is_anti_particle() != p2.is_anti_particle()) return false;
	if (p1.mass_label() != p2.mass_label()) return false;
	return true;
}

bool operator!=(const particle_ID& p1, const particle_ID& p2){
	if (p1 == p2) return false; else return true;
}

bool operator<(const particle_ID& p1, const particle_ID& p2){
	if ((*p1.type()) < (*p2.type())) return true;
	if ((*p1.type()) > (*p2.type())) return false;
	if (p2.is_anti_particle() && !( p1.is_anti_particle())) return true;
	if (p1.is_anti_particle() && !( p2.is_anti_particle())) return false;
	if (p1.helicity() < p2.helicity()) return true;
	if (p1.helicity() > p2.helicity()) return false;
	if (p1.flavor() < p2.flavor()) return true;
	if (p1.flavor() > p2.flavor()) return false;
	return false;
}
// not exactly equal to !operator< since we want it to return false for equal (or indifferenciable) elements
bool operator>(const particle_ID& p1, const particle_ID& p2){
	if ((*p1.type()) > (*p2.type())) return true;
	if ((*p1.type()) < (*p2.type())) return false;
	if (p1.is_anti_particle() && !( p2.is_anti_particle())) return true;
	if (p2.is_anti_particle() && !( p1.is_anti_particle())) return false;
	if (p1.helicity() > p2.helicity()) return true;
	if (p1.helicity() < p2.helicity()) return false;
	if (p1.flavor() > p2.flavor()) return true;
	if (p1.flavor() < p2.flavor()) return false;
	return false;
}


particle_ID particle_ID::conjugate() const {
	switch (_type->stat()) {
	case particle::boson: {
		return particle_ID(_type,-_helicity,_flavor,false);
	}
	case particle::fermion: {
		return particle_ID(_type,-_helicity,_flavor,!(_is_anti_particle));
	}
	default: return *this;
	}
}

particle_ID particle_ID::helicity_conjugate() const {
	return particle_ID(_type,-_helicity,_flavor,_is_anti_particle);
}

ostream& operator<<(ostream& s, const particle_ID& r){
	s<<r.name();
	if ( r.is_anti_particle() ){
			s << "b";
		}
	if ( r.flavor() != 1 ){
			s << r.flavor();
		}
		switch (r.helicity()) {
		case 1 : s << "+"; break;
		case -1 :s << "-";break;
		case 0 :s << "0";break;
		default :s << "";
	}

	if(r.mass_label()!=0){
		C mass;
		eval_param<R>::get_mass_param(r.mass_label()).mass(mass);
		s<< "[" << mass << "]";
}

	return s;

}

	
mass_param_library setup_masses_original(){
	mass_param_library returnres;
	#ifdef BH_USE_GMP
	GMP_INIT
	#endif
	//The mass_params to define the massive particles
	mass_param rat_ext_mass_cre(C(0,0),CHP(0,0),CVHP(0,0));
	mass_param top_ext_mass_cre(C(172,0),CHP(172,0),CVHP(172,0));
	mass_param top_mass_cre(C(172,0),CHP(172,0),CVHP(172,0));
	mass_param higgs_mass_cre(C(1,0),CHP(1,0),CVHP(1,0)); // Chosen for debugging purposes, unless the higgs really does turn out to have a 1GeV mass	
	
	// Add all mass types to this list so that BH knows about them
	returnres.add(rat_ext_mass_cre);
	returnres.add(top_ext_mass_cre);
	returnres.add(top_mass_cre);
	returnres.add(higgs_mass_cre);
		
	return returnres;
}

// The mass_param_labels have already been created so we need to just add the correct ones in
mass_param_library setup_masses(){
	
	// If this is the first time we call this function then the R mass_param_library will be empty otherwise we have already created the mass_params
	if(eval_param<R>::nbr_masses()==0){
		return setup_masses_original();
	}
	else{ // We have created the mass_params and so we simply fill up the library with them
		mass_param_library returnres;
		
		// Add all mass types to this list so that BH knows about them
		returnres.add(eval_param<R>::modify_mass_param(1));
		returnres.add(eval_param<R>::modify_mass_param(2));
		returnres.add(eval_param<R>::modify_mass_param(3));
		returnres.add(eval_param<R>::modify_mass_param(4));
		
		return returnres;
	}
}

template <class T> mass_param_library eval_param<T>::_masses=setup_masses();
template mass_param_library eval_param<R>::_masses;
template mass_param_library eval_param<RHP>::_masses;
template mass_param_library eval_param<RVHP>::_masses;
	
#if BH_USE_GMP
template mass_param_library eval_param<RGMP>::_masses;
#endif
	// Setup the masses
	
// The mass_params to define the massive particles for use globaly, we assume that these were added to the mass_param_library in this order starting from zero
mass_param rat_ext_mass(eval_param<R>::get_mass_param(1));
mass_param top_ext_mass(eval_param<R>::get_mass_param(2));
mass_param top_mass(eval_param<R>::get_mass_param(3));
mass_param higgs_mass(eval_param<R>::get_mass_param(4));
	
	//Set up the particles
	
// please note that different particles must have different pdg_code!

particle quark(particle::fermion,"q",0.,1);
particle lepton(particle::fermion,"l",0.,11,particle::qc_type);
particle gluon(particle::boson,"g",0.,21);
particle photon(particle::boson,"y",0.,8,particle::qc_type);
particle gluino(particle::fermion,"G",0.,1000);
particle scalar(particle::boson,"s",0.,-5);
particle higgs(particle::boson,"h",1.,25,particle::gc_type,higgs_mass.label());

// Particles needed for the rational extraction code, we pass a dummy mass
particle gluon_massive_scalar(particle::boson,"Rsc",1.,-1,particle::ord_type,rat_ext_mass.label());
particle gluon_massive(particle::boson,"R",1.,-6,particle::ord_type,rat_ext_mass.label());
particle quark_massive(particle::fermion,"Q",1.,-2,particle::ord_type,rat_ext_mass.label());
particle gluino_massive(particle::fermion,"L",1.,-3,particle::ord_type,rat_ext_mass.label());
particle scalar_massive(particle::boson,"S",1.,-4,particle::ord_type,rat_ext_mass.label());
	
particle_ID m(gluon,-1);
particle_ID p(gluon,1);
particle_ID ym(photon,-1);
particle_ID yp(photon,1);
particle_ID y2m(photon,-1,2);
particle_ID y2p(photon,1,2);
particle_ID qp(quark,1);
particle_ID qm(quark,-1);
particle_ID q2p(quark,1,2);
particle_ID q2m(quark,-1,2);
particle_ID q3p(quark,1,3);
particle_ID q3m(quark,-1,3);
particle_ID q4p(quark,1,4);
particle_ID q4m(quark,-1,4);
particle_ID qbp(quark,1,1,true);
particle_ID qbm(quark,-1,1,true);
particle_ID qb2p(quark,1,2,true);
particle_ID qb2m(quark,-1,2,true);
particle_ID qb3p(quark,1,3,true);
particle_ID qb3m(quark,-1,3,true);
particle_ID qb4p(quark,1,4,true);
particle_ID qb4m(quark,-1,4,true);
particle_ID lp(lepton,1);
particle_ID lm(lepton,-1);
particle_ID lbp(lepton,1,1,true);
particle_ID lbm(lepton,-1,1,true);
particle_ID l2p(lepton,1,2,false);
particle_ID l2m(lepton,-1,2,false);
particle_ID lb2p(lepton,1,2,true);
particle_ID lb2m(lepton,-1,2,true);
particle_ID s0(scalar,0);
particle_ID Gp(gluino,1);
particle_ID Gm(gluino,-1);
particle_ID Gbp(gluino,1,1,true);
particle_ID Gbm(gluino,-1,1,true);
particle_ID G2p(gluino,1,2,false);
particle_ID G2m(gluino,-1,2,false);
particle_ID Gb2p(gluino,1,2,true);
particle_ID Gb2m(gluino,-1,2,true);
particle_ID G3p(gluino,1,3,false);
particle_ID G3m(gluino,-1,3,false);
particle_ID Gb3p(gluino,1,3,true);
particle_ID Gb3m(gluino,-1,3,true);

// Particles needed for the rational extraction code
particle_ID gsc(gluon_massive,0,1,false);
particle_ID gpM(gluon_massive,1,1,false);
particle_ID gmM(gluon_massive,-1,1,false);
particle_ID gD(gluon_massive_scalar,0,1,false);
particle_ID scM(scalar_massive,0,1,false);
particle_ID Qp(quark_massive,1,1,false);
particle_ID Qm(quark_massive,-1,1,false);
particle_ID Qbp(quark_massive,1,1,true);
particle_ID Qbm(quark_massive,-1,1,true);
particle_ID Q2p(quark_massive,1,2,false);
particle_ID Q2m(quark_massive,-1,2,false);
particle_ID Qb2p(quark_massive,1,2,true);
particle_ID Qb2m(quark_massive,-1,2,true);
particle_ID Lp(gluino_massive,1,1,false);
particle_ID Lm(gluino_massive,-1,1,false);
particle_ID Lbp(gluino_massive,1,1,true);
particle_ID Lbm(gluino_massive,-1,1,true);
particle_ID Lp2(gluino_massive,1,2,false);
particle_ID Lm2(gluino_massive,-1,2,false);
particle_ID Lb2p(gluino_massive,1,2,true);
particle_ID Lb2m(gluino_massive,-1,2,true);



// Higgs particles
particle_ID ph(higgs,0,1,false); //complex scalar field component of H
particle_ID phd(higgs,0,1,true); //complex scalar field component of H
particle_ID H(higgs,0,0,false); //The full higgs H
	

}
