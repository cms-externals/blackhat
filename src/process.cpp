/*
 * process.cpp
 *
 *  Created on: 11-Nov-2008
 *      Author: daniel
 */

#include <vector>
#include "process.h"
#include "BH_error.h"
#include <sstream>

using namespace std;

namespace BH {

long compute_pcode(const vector<particle_ID>& p){
  	int nbr_q=0;
  	int nbr_l=0;
  	int nbr_g=0;
  	int nbr_higgs=0;
  	int nbr_gsc_massive=0;
	int nbr_gsc_D_massive=0;
  	int nbr_q_massive=0;
  	int nbr_photon=0;
  	int nbr_gluino=0;
  	int nbr_gluino_massive=0;
	int nbr_sc_massless=0;
	int nbr_sc_massive=0;

  	for (size_t i=0;i<p.size();i++) {
		switch (p[i].type()->pdg_code()){
			case 1: nbr_q++; break; // quark
			case 11: nbr_l++; break; // lepton
			case 21: nbr_g++; break; // gluon
			case 25: nbr_higgs++; break; // higgs
			case -1: nbr_gsc_D_massive++; break; // massive gluon
            case -6: nbr_gsc_massive++; break; // massive D dim gluon helicities
			case -2: nbr_q_massive++; break; // massive quark
			case 8: nbr_photon++; break; // photon
			case 1000: nbr_gluino++; break; // gluinos
			case -3: nbr_gluino_massive++; break; // massive gluinos
			case -4: nbr_sc_massive++; break; // massive scalars
			case -5: nbr_sc_massless++; break; // massless scalars
  	  	}
  	}

  	return
  		nbr_g
  		+10*nbr_q
  		+100*nbr_l
  		+1000*nbr_sc_massive
  		+10000*nbr_q_massive
		+100000*nbr_photon
		+1000000*nbr_gluino
		+10000000*nbr_gluino_massive
		+100000000*nbr_higgs
		+1000000000*nbr_sc_massless
        +10000000000*nbr_gsc_massive
		+100000000000*nbr_gsc_D_massive;
}

process::process(const vector<particle*>& p,vector<short> hel){
if (p.size()!=hel.size()) throw BHerror("vector<particle_type> p and vector<short> hel in process constructor are not of the same size.");
	for (size_t j=0;j<p.size();j++){
		particles.push_back(particle_ID(p[j],hel[j],1,false));
	}
	nbr=p.size();
	d_pcode=compute_pcode(particles);
}

process::process(const vector<plabel>& pl) : particles(pl.begin(),pl.end()) {
//	for (size_t j=0;j<pl.size();j++){
//		particles.push_back(particle_ID(pl[j].type(),pl[j].helicity(),pl[j].flavor(),pl[j].is_anti_particle()));
//	}
	nbr=pl.size();
	d_pcode=compute_pcode(particles);
}

process::process(vector<particle_ID>& pl,long pcode) : d_pcode(pcode) {
	nbr=pl.size();
	particles.swap(pl);
}

process::process(const vector<particle_ID>& pl){
	copy(pl.begin(),pl.end(),back_inserter(particles));
//	for (size_t j=0;j<pl.size();j++){
//		particles.push_back(pl[j]);
//	}
	nbr=pl.size();
	d_pcode=compute_pcode(particles);
}

process::process(particle_ID p){
	particles.push_back(p);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2){
	particles.push_back(p1);
	particles.push_back(p2);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4,particle_ID p5){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4,particle_ID p5,particle_ID p6){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	particles.push_back(p6);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4,particle_ID p5,particle_ID p6,particle_ID p7){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	particles.push_back(p6);
	particles.push_back(p7);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4,particle_ID p5,particle_ID p6,particle_ID p7,particle_ID p8){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	particles.push_back(p6);
	particles.push_back(p7);
	particles.push_back(p8);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4,particle_ID p5,particle_ID p6,particle_ID p7,particle_ID p8,particle_ID p9){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	particles.push_back(p6);
	particles.push_back(p7);
	particles.push_back(p8);
	particles.push_back(p9);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4,particle_ID p5,particle_ID p6,particle_ID p7,particle_ID p8,particle_ID p9,particle_ID p10){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	particles.push_back(p6);
	particles.push_back(p7);
	particles.push_back(p8);
	particles.push_back(p9);
	particles.push_back(p10);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4,particle_ID p5,particle_ID p6,particle_ID p7,particle_ID p8,particle_ID p9,particle_ID p10,particle_ID p11){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	particles.push_back(p6);
	particles.push_back(p7);
	particles.push_back(p8);
	particles.push_back(p9);
	particles.push_back(p10);
	particles.push_back(p11);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}
process::process(particle_ID p1,particle_ID p2,particle_ID p3,particle_ID p4,particle_ID p5,particle_ID p6,particle_ID p7,particle_ID p8,particle_ID p9,particle_ID p10,particle_ID p11,particle_ID p12){
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);
	particles.push_back(p4);
	particles.push_back(p5);
	particles.push_back(p6);
	particles.push_back(p7);
	particles.push_back(p8);
	particles.push_back(p9);
	particles.push_back(p10);
	particles.push_back(p11);
	particles.push_back(p12);
	nbr=particles.size();
	d_pcode=compute_pcode(particles);
}

std::string process::print() const {
	stringstream s;
	s << *this;
	return s.str();
}

ostream& operator<<(ostream& s, const process& p){
	s<<"(";
	for (size_t j=0;j<p.particles.size()-1;j++){
		s << p.particles[j] << ",";
	}
	return s << p.particles.back() <<")";
}

bool operator==(const process& p1, const process& p2){
	if (p1.n() != p2.n() ) return false;
	for (int i=1;i<=p1.n();i++){
		if (p1.p(i) != p2.p(i)) return false;
	}
	return true;
}

bool operator<(const process& p1, const process& p2){
	if ( p1.n() != p2.n() ) return ( p1.n() < p2.n() ) ;
	for (int i=1;i<=p1.n();i++){
		if (p1.p(i) < p2.p(i)) return true;
		if (p2.p(i) < p1.p(i)) return false;
	}
	return false;
}


}
