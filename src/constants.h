/*
 * constants.h
 *
 *  Created on: May 5, 2009
 *      Author: daniel
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <iosfwd>
#include "BH_typedefs.h"
#include "BH_utilities.h"
#include "multi_precision_constant.h"

namespace BH {

class settings_table;

int gcd(int a, int b);

inline RHP toRHP(long i){ return RHP(static_cast<double>(i),0.0);}
inline RVHP toRVHP(long i){ return RVHP(static_cast<double>(i),0.0,0.0,0.0);}


class multi_precision_fraction {
	friend std::ostream& operator<<(std::ostream& s,const multi_precision_fraction& mpf);
	friend const multi_precision_fraction operator*(const multi_precision_fraction& mpf1,const multi_precision_fraction& mpf2);
	friend const multi_precision_fraction operator*(int,const multi_precision_fraction& mpf2);
	friend const multi_precision_fraction operator*(const multi_precision_fraction& mpf1,int);
	friend const multi_precision_fraction operator+(const multi_precision_fraction& mpf1,const multi_precision_fraction& mpf2);
	friend const multi_precision_fraction operator-(const multi_precision_fraction& mpf1,const multi_precision_fraction& mpf2);
public:
	long m_num;
	long m_den;
	multi_precision_fraction(long num,long den=1): m_num(num), m_den(den){ normalize();};
	operator R() const {return R(m_num)/R(m_den);};
	operator RHP() const {return toRHP(m_num)/toRHP(m_den);};
	operator RVHP() const {return toRVHP(m_num)/toRVHP(m_den);};
//	void set(R val,RHP val_HP,RVHP val_VHP){ m_value=val, m_value_HP=val_HP ,m_value_VHP=val_VHP;};
	const multi_precision_fraction operator-(){return multi_precision_fraction(-m_num,m_den);};
template <class T> T get(){return R(m_num)/R(m_den);}
private:
	void normalize();
};

struct constants {
	static multi_precision_constant MZ;
	static multi_precision_constant GZ;
	static multi_precision_constant MW;
	static multi_precision_constant GW;
	static multi_precision_constant Mtop;
	static multi_precision_constant G3_Lambda2;
	static R sin_th_2;
	static R sin_2th;
	static R vel;
	static R ver;
	static R vnuel;
	static R vnuer;
	static R vupl;
	static R vupr;
	static R vdownl;
	static R vdownr;
	static R alpha_S;
	static R alpha_QED;
	static void update_constants(settings_table*);
	static void set_GeV(R GeV){s_GeV=GeV;};
	// this static variable transforms between GeV units and O(1) units
	static R s_GeV;

//	static R get_GeV(){return s_GeV;};
//	//static R get_GeV();
//	static multi_precision_constant get_Mtop(){return m_Mtop;};
//	//static multi_precision_constant get_Mtop();

};

}

#endif /* CONSTANTS_H_ */
