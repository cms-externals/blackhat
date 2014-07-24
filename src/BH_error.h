/*!\file BH_error.h
  \brief Header for error handling in BH
*/
#ifndef BH_ERROR_H_
#define BH_ERROR_H_

#include <iostream>
#include <string>

namespace BH {

class BHerror {
	std::string d_descr;
public:
	BHerror() {};
	BHerror(const char* d) : d_descr(d) {};
	BHerror(const std::string& d) : d_descr(d) {};
#ifdef SWIG
%rename(display) print;
#endif
	void print() { std::cerr << d_descr << std::endl;}
	void set_descr(const std::string& str){ d_descr=str; };
};

}
#endif /*BH_ERROR_H_*/
