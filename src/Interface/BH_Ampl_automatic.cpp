/*
 * matrix elements for lm lbp -> 2q2g
 *
 *  Created on: Oct 21, 2008
 *
 */



#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "cached_OLHA.h"
#include "BH_interface_impl.h"
//#include "ME2.h"

#define _VERBOSE 0


using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;
namespace BH {

#ifndef BH_PUBLIC

BH_Ampl_automated::BH_Ampl_automated(
		std::vector<std::pair<int,int> > particle_labels,
		ME2_factory* me2_factory,
		approx born_or_virt,
		BH_interface_impl* bhi,
		int NbrExtParticles,
		int NbrPowersOfAlphaS,
		int NbrPowersOfAlphaQED,
		int GeVdim,
		double scheme_shift,
		double renormalization_shift):
	BH_Ampl_impl(me2_factory->new_ME2(born_or_virt,particle_labels),	
			bhi, 
			NbrExtParticles,
			NbrPowersOfAlphaS,
			NbrPowersOfAlphaQED,
			GeVdim,
			1.,
			scheme_shift,
			renormalization_shift){
		for(int i=0;i<particle_labels.size();i++){
			momenta_assignment.push_back(particle_labels[i].first);
			particles.push_back(particle_labels[i].second);
		};
};
#endif

BH_Ampl_data::BH_Ampl_data(
		std::vector<std::pair<int,int> > particle_labels,
		QCDorder lo_or_nlo,
		BH_interface_impl* bhi,
		int NbrExtParticles,
		int NbrPowersOfAlphaS,
		int NbrPowersOfAlphaQED,
		int GeVdim,
		double scheme_shift,
		double renormalization_shift_l,
        double d_lc_l,
        double d_fc_l
        ):
	BH_Ampl_impl(
			get_ME2_from_file(/*const string& filename,*/ lo_or_nlo,particle_labels),
//			me2_factory->new_ME2(born_or_virt,particle_labels),	
			bhi, 
			NbrExtParticles,
			NbrPowersOfAlphaS,
			NbrPowersOfAlphaQED,
			GeVdim,
			1.,
			scheme_shift,
			renormalization_shift_l
			),
    renormalization_shift(renormalization_shift_l),
    d_lc(d_lc_l), d_fc(d_fc_l){
		for(int i=0;i<particle_labels.size();i++){
			momenta_assignment.push_back(particle_labels[i].first);
			particles.push_back(particle_labels[i].second);
		};
        
};

double BH_Ampl_data::getScaleVariationCoefficient(int logMuOrder){
    // returns bare matrix element squared

	switch (logMuOrder) {
	case 0: return get_finite();
	case 1:
		switch(settings::BH_interface_settings::s_BH_color_mode){
		    //Full color and Leading color
            case 0: case 1: return get_single_pole()+renormalization_shift;
		    //Full minus Leading color
		    case 2: return get_single_pole();
        }
	case 2:	return get_double_pole();
	}
};

}
