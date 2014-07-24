/*
 * Known_Cut_wCI.h
 *
 *  Created on: 23 Sep 2009
 *      Author: daniel
 */

#ifndef KNOWN_CUT_WCI_H_
#define KNOWN_CUT_WCI_H_

#include "BH_typedefs.h"
#include <vector>

namespace BH {

class process;

namespace CachedIntegral {

class Cut_Part_wCI;

Cut_Part_wCI* CwCI_Ptr(const process& pro,color_structure cs,const std::vector<int>& ind);

} /* CachedIntegral */

} /* BH */


#endif /* KNOWN_CUT_WCI_H_ */
