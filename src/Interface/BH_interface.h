/*
 * BH_interface.h
 *
 *  Created on: 06-Mar-2009
 *      Author: daniel
 */

#ifndef BH_INTERFACE_H_
#define BH_INTERFACE_H_


#include <vector>
#include <string>

#include "setable.h"
#define BH_INTERFACE_BHSETTINGS


namespace BH {



class BH_interface_impl;
class BH_Ampl;
class BHinput;

class BH_interface : public setable {
public:
	BH_interface();
	BH_interface(const std::string&);
	virtual ~BH_interface();
	//! Return a pointer to a new BH_Ampl object.
	BH_Ampl* new_ampl(const std::vector<int>&);
	BH_Ampl* new_tree_ampl(const std::vector<int>&);
	//! Sets the named setting to the value provided.
	virtual bool set(const std::string& name,double value);
	virtual bool set(const std::string& name,int value);
	virtual bool set(const std::string& name,std::string value);
	virtual bool set(const std::string& name,bool value);
	template <class T> bool apply_setting(const std::string& name, T& param);
	//! Prints a table of the settings and their current values.
	void print_settings(std::ostream& os);
	void operator()(BHinput& in);
	void print_banner();
private:
	BH_interface_impl* d_impl;
};



}
#endif /* BH_INTERFACE_H_ */
