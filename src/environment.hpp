#ifndef _H_ENVIRONMENT_HPP
#define _H_ENVIRONMENT_HPP

#include <cstdlib>
#include <string>

namespace BH {

namespace Tools {

template <class Target> Target stringConvert(const std::string&);

template <> int stringConvert<int>(const std::string& s) {
	return atoi(s.c_str());
};


template <class Type> Type readFromEnvironment(const char* name,Type dflt) {
	std::string value("");
	char* ptr = std::getenv (name);
	if ( ptr ) {
		value=ptr;
		return stringConvert<Type>(value);
	} else {
		return dflt;
	};
}



} /* Tools */



} /* BH */


#endif  /* _H-ENVIRONMENT_HPP*/