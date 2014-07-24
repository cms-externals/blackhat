/*
 * filename.h
 *
 *  Created on: 7 Jul 2009
 *      Author: daniel
 */

#ifndef FILENAME_H_
#define FILENAME_H_

#include <iosfwd>
#include <string>
#include "BH_typedefs.h"

namespace BH {
	class process;
}

namespace BH {

namespace ratext {

std::string rat_filename(const process& pro,color_structure cs);

}
}



#endif /* FILENAME_H_ */
