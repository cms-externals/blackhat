/*
 * process_utils_massive.h
 *
 *  Created on: 28-Jul-2008
 *      Author: Darren
 */

#ifndef PROCESS_UTILS_MASSIVE_H_
#define PROCESS_UTILS_MASSIVE_H_

namespace BH {

// Contained in cut_part_factory.cpp
bool are_quarks_adjacent(const process& pro);
std::vector<size_t> find_indices(const process& pro,const particle& t);

// Replaces massive gluinos with massive quarks so that we can reuse massive quark code
process replace_massive_gluino_with_massive_quark(const process& pro);

// Arranges the initial flavours of the process
process arrange_flavors_massive(const process& pro,size_t qb_pos,size_t q_pos,short increment,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2e_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_4q2e_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G2e_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q4G2e_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G2e_q_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G1y_q_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G2e_G_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G1y_G_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q2G1y_massive(const process& pro,std::vector<particle_ID>& propagators);
void arrange_possible_props_2q1y_massive(const process& pro,std::vector<particle_ID>& possible_props);
process arrange_flavors_2q2G_LC_massive(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_2q4G_LC_massive(const process& pro,std::vector<particle_ID>& propagators);

process arrange_flavors_qq_massive_LT(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_qq_massive_RT(const process& pro,std::vector<particle_ID>& propagators);

process arrange_flavors_4q_massive_LC(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_6q_massive_LC(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_qqX_massive_LC(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_qqQQX_massive_LC(const process& pro,std::vector<particle_ID>& propagators);

process arrange_flavors_qqX_massive_SLC(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_qqee_massive_SLC(const process& pro,std::vector<particle_ID>& propagators);
process arrange_flavors_qqQQX_massive_SLC(const process& pro,std::vector<particle_ID>& propagators);

}

#endif /* PROCESS_UTILS_MASSIVE_H_ */
