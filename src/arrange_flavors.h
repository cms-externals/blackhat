/*
 * arrange_flavors.h
 *
 *  Created on: 22-Apr-2009
 *      Author: daniel
 */

#ifndef ARRANGE_FLAVORS_H_
#define ARRANGE_FLAVORS_H_

#include <stddef.h>
#include <vector>

namespace BH {

class process;
class particle_ID;

process arrange_flavors(const process& pro,size_t qb_pos,size_t q_pos,short increment,std::vector<particle_ID>& propagators);
std::vector<int> get_unordered_gluons_2q1y(const process& PRO);
std::vector<int> get_unordered_gluons_2q2e(const process& PRO);

}
#endif /* ARRANGE_FLAVORS_H_ */
