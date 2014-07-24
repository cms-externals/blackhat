/*
 * settings_reader.h
 *
 *  Created on: 25-Feb-2009
 *      Author: daniel
 */

#ifndef SETTINGS_READER_H_
#define SETTINGS_READER_H_

#include <string>

namespace BH {

namespace settings {


void read_settings_from_file(const std::string& path,bool printWarning=false);
void use_setting(const std::string& option);
bool read_from_stream(std::istream& in);

}



}

#endif /* SETTINGS_READER_H_ */
