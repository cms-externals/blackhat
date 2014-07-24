/*
 * BH_interface_impl.h
 *
 *  Created on: 09-Mar-2009
 *      Author: daniel
 */

#ifndef BH_INTERFACE_IMPL_H_
#define BH_INTERFACE_IMPL_H_

#include "constants.h"
#include "BH_typedefs.h"
#include "settings.h"

namespace BH {

template <class T>  class momentum_configuration;

class BH_interface_impl {
protected:
	settings_table* d_settings_p;
	momentum_configuration<double>* d_mc_p;
	double d_mu;
public:
	//! Return a pointer to a new BH_Ampl object.
	virtual BH_Ampl* new_ampl(const std::vector<int>&, QCDorder lo_or_nlo = nlo )=0;
	virtual void operator()(BHinput& in)=0;
	template <class T> void set(const std::string& name,T value);
	double get_mu(){return d_mu;}
	momentum_configuration<double>* get_mc(){return d_mc_p;};
	// updates the values of the constant
	void update_constants(){constants::update_constants(d_settings_p);};
	settings_table* getSettingsTable(){ return d_settings_p;};
	void useSettingsTable(settings_table* settings_table){d_settings_p=settings_table;};
	template <class T> bool apply_setting(const std::string& name, T& param){return d_settings_p->apply_setting(name,param);};
	BH_interface_impl();
	virtual void remove_last(){};

	virtual ~BH_interface_impl();
};

}

#endif /* BH_INTERFACE_IMPL_H_ */
