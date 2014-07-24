/*
 * filename.cpp
 *
 *  Created on: 7 Jul 2009
 *      Author: daniel
 */

#ifndef BH_PUBLIC
#include "ratext/data_files.h"
#endif
#include "cut_worker.h"
#include "OneLoopHelAmpl.h"
#include <iostream>
#include <cassert>
#include "ratext/ratext_part.h"
#include "settings.h"
#include "from_file.h"
#include <cstdlib>
#include <unistd.h>
using namespace std;

namespace BH {
	string string_name(const process& pro);
}


namespace BH {

namespace ratext {


string rat_filename(const process& pro,color_structure cs){
	stringstream ss;
	ss << get_worker_dir("rat");
	ss << "/" << pro.n() ;
	if ( access( ss.str().c_str(), 0 ) != 0 ){
			_WARNING3("Data path ",ss.str(),"not present. Please create it. ");
			throw BHerror("Missing path");
		}
	ss << "/rat_";
	ss << pro;
	ss << "_";
	ss << cs;
	ss << ".dat";
	return ss.str();
}




}

}
