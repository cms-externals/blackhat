/*
 * path.cpp
 *
 *  Created on: 3 Mar 2010
 *      Author: daniel
 */

#include <string>
#include "BHpath.h"

using namespace std;

std::string GetDataPath(){
	static string path = string(BH_INSTALL_PATH) + string ("/share/blackhat/") ;
	return path;
}
std::string GetSrcPath(){
	static string path = string(BH_SOURCE_PATH)  ;
	return path;
}
