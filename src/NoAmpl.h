/*
 * NoAmpl.h
 *
 *  Created on: 5 Jul 2009
 *      Author: daniel
 */

#ifndef NOAMPL_H_
#define NOAMPL_H_

#include "BH_error.h"
#include "process.h"
#include <sstream>

namespace BH {


class NoAmpl : public BHerror {
	process d_process;
	color_structure d_cs;
public:
	NoAmpl(process PRO,color_structure cs) : d_process(), d_cs(),  BHerror()	{
	std::string msg;
	std::stringstream ss;
	ss << "Amplitude for process " << PRO << " and color structure " << cs << " is not supported.";
	set_descr(ss.str());

	};
};

class NoTreeAmpl : public BHerror {
	process d_process;
public:
	NoTreeAmpl(process PRO) : d_process(),   BHerror()	{;
	std::string msg;
	std::stringstream ss;
	ss << "Amplitude for process " << PRO << " is not supported." ;
	set_descr(ss.str());
	};
};


}


#endif /* NOAMPL_H_ */
