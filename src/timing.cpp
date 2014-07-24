/*
 * timing.cpp
 *
 *  Created on: Oct 30, 2009
 *      Author: daniel
 */

#include "BH_utilities.h"

#include "timing.h"

namespace BH {

namespace timing {

bool counter::s_is_open=false;
std::ofstream counter::ofile;

counter::~counter(){
	if ( ! is_open()){ _MESSAGE("opening file");};
	ofile << d_name << ": " << d_nbr_calls << " calls, total: " << double(d_total_time)/double(CLOCKS_PER_SEC) << " s (" << double(d_total_time)/double(CLOCKS_PER_SEC)/double(d_nbr_calls) << " s per call)\n";
}

} /* timing */
} /* BH */
