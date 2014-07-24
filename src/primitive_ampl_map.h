/*
 * primitive_ampl_map.h
 *
 *  Created on: 15 Jan 2011
 *      Author: 
 */

#ifndef PRIMITIVE_AMPL_MAP_H_
#define PRIMITIVE_AMPL_MAP_H_

#include <vector>
#include <string>
#include "process.h"
#include "particles.h"
#include "process_utils.h"

using namespace std;

namespace BH {

// Useful helper functions
void flip_cs_at(size_t pos, string& cs);
void flip_cs(string& cs);
bool is_Ltype_cs(string cs);


// Mapping functions
void conjugateQ(vector<plabel>& plabels , double & sign,short & conjQ, string & cs);

void canonical_pro(vector<plabel>& plabels, double & sign, short & conjQ, string & cs);
void canonical_pro_tree(vector<plabel>& plabels, double & sign, short & conjQ);

}

#endif /* PRIMITIVE_AMPL_MAP_H_ */
