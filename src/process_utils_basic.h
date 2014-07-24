/*
 * process_utils_basic.h
 *
 *  Created on: 02.08.2012
 *      Author: daniel
 */

#ifndef PROCESS_UTILS_BASIC_H_
#define PROCESS_UTILS_BASIC_H_

#include "BH_typedefs.h"
#include "particles.h"
#include "process.h"
#include <vector>


namespace BH {


std::vector<int> find_all_flavors(const process& pro,const particle& type);
std::vector<int> find_all_flavors(const process& pro,const particle& type1, const particle& type2);
std::vector<int> find_all_flavors(const process& pro,const particle& type1, const particle& type2, const particle& type3);

process replace_one_photon_with_gluon(const process& pro);
process replace_lepton_with_quark(const process& pro);
process replace_photon_with_gluon(const process& pro);
process replace_gluino_with_quark(const process& pro);
process replace_gluino_with_quark(const process& pro,int flavor);
process replace_scalar_massive_with_gluon(const process& pro);
process replace_massive_scalars_with_gluons(const process& pro);
process replace_gluon_massive_with_scalar(const process& pro);

process fix_flavors(const process& pro);
process fix_flavors_2q2e(const process& pro);

class particle_pattern {
public:
	virtual bool match(const particle_ID&)=0;
	virtual ~particle_pattern(){};
};

typedef std::pair<particle_pattern*,std::string> rule;
typedef std::vector<rule> rule_list;


std::string string_gen(const std::vector<particle_ID>& ps,rule_list& rules);
std::string string_gen(const process& pro,rule_list& rules);

color_structure right_direction(const process& pro,color_structure cs);


} /* BH */
#endif /* PROCESS_UTILS_BASIC_H_ */
