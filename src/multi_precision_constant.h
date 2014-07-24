/*
 * multi_precision_constant.h
 *
 *  Created on: 11 Mar 2010
 *      Author: daniel
 */

#ifndef MULTI_PRECISION_CONSTANT_H_
#define MULTI_PRECISION_CONSTANT_H_

#if BH_USE_GMP
#include "gmp_r.h"
#endif
#include "qd_suppl.h"

namespace BH {

class multi_precision_constant {
	R m_value;
	RHP m_value_HP;
	RVHP m_value_VHP;
#if BH_USE_GMP
	RGMP m_value_RGMP;
#endif
public:
#if BH_USE_GMP
	explicit multi_precision_constant(int i): m_value(i), m_value_HP(i) ,m_value_VHP(i), m_value_RGMP(i) {};
	explicit multi_precision_constant(double d): m_value(d), m_value_HP(d) ,m_value_VHP(d), m_value_RGMP(d) {};
	explicit multi_precision_constant(const RHP& value): m_value(to_double(value)), m_value_HP(value) ,m_value_VHP(value), m_value_RGMP(value.to_string().c_str()) {};
	explicit multi_precision_constant(const RVHP& value): m_value(to_double(value)), m_value_HP(to_HP(value)) ,m_value_VHP(value) , m_value_RGMP(value.to_string().c_str()){};
	multi_precision_constant(R val,const RHP& val_HP,const RVHP& val_VHP,const RGMP& val_GMP): m_value(val), m_value_HP(val_HP) ,m_value_VHP(val_VHP), m_value_RGMP(val_GMP) {};
#else
	explicit multi_precision_constant(int i): m_value(i), m_value_HP(i) ,m_value_VHP(i) {};
	explicit multi_precision_constant(double d): m_value(d), m_value_HP(d) ,m_value_VHP(d) {};
	explicit multi_precision_constant(RHP& value): m_value(to_double(value)), m_value_HP(value) ,m_value_VHP(value) {};
	explicit multi_precision_constant(RVHP& value): m_value(to_double(value)), m_value_HP(to_HP(value)) ,m_value_VHP(value) {};
	multi_precision_constant(R val,RHP val_HP,RVHP val_VHP): m_value(val), m_value_HP(val_HP) ,m_value_VHP(val_VHP){};
#endif
	operator R() const {return m_value;};
	operator RHP() const {return m_value_HP;};
	operator RVHP() const {return m_value_VHP;};
#if BH_USE_GMP
	operator RGMP() const {return m_value_RGMP;};
	void set(R val,const RHP& val_HP,const RVHP& val_VHP,const RGMP& val_GMP){ m_value=val, m_value_HP=val_HP ,m_value_VHP=val_VHP; m_value_RGMP=val_GMP;};
	void set(R val){ m_value=val, m_value_HP=val ,m_value_VHP=val; m_value_RGMP=val;};
#else
	void set(R val,RHP val_HP,RVHP val_VHP){ m_value=val, m_value_HP=val_HP ,m_value_VHP=val_VHP;};
	void set(R val){ m_value=val, m_value_HP=val ,m_value_VHP=val;};
#endif
};

}
#endif /* MULTI_PRECISION_CONSTANT_H_ */
