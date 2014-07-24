/*
 * constants.cpp
 *
 *  Created on: May 5, 2009
 *      Author: daniel
 */

#include "constants.h"
#include "settings.h"
#include <iostream>
#include <cstdlib>
#include "multi_precision_constant.h"
using namespace std;

namespace BH {

const multi_precision_fraction operator*(const multi_precision_fraction& mpf1,const multi_precision_fraction& mpf2){
	return multi_precision_fraction(mpf1.m_num*mpf2.m_num,mpf1.m_den*mpf2.m_den);
}
const multi_precision_fraction operator*(const multi_precision_fraction& mpf1,int i){
	return multi_precision_fraction(mpf1.m_num*i,mpf1.m_den);
}
const multi_precision_fraction operator*(int i,const multi_precision_fraction& mpf2){
	return multi_precision_fraction(i*mpf2.m_num,mpf2.m_den);
}
const multi_precision_fraction operator+(const multi_precision_fraction& mpf1,const multi_precision_fraction& mpf2){
	return multi_precision_fraction(mpf1.m_num*mpf2.m_den+mpf2.m_num*mpf1.m_den,mpf1.m_den*mpf2.m_den);
}
const multi_precision_fraction operator-(const multi_precision_fraction& mpf1,const multi_precision_fraction& mpf2){
	return multi_precision_fraction(mpf1.m_num*mpf2.m_den-mpf2.m_num*mpf1.m_den,mpf1.m_den*mpf2.m_den);
}

ostream& operator<<(ostream& s,const multi_precision_fraction& mpf){
	return s << "(" << mpf.m_num << " " << mpf.m_den << ")";
;
}



//R constants::sin_th_2=0.231;
R constants::sin_th_2=0.23; //matches Sherpa
//R constants::sin_2th=0.8429448380528823;
R constants::sin_2th=0.8416650165000326; //matches Sherpa
//multi_precision_constant constants::m_MZ(91.1876,RHP(91.1876),RVHP(91.1876)); //Z-boson

/* older values
#if BH_USE_GMP
multi_precision_constant constants::MZ(91.188,RHP(91.188),RVHP(91.188),RGMP("91.188")); //Z-boson like Sherpa
multi_precision_constant constants::GZ(2.49,RHP(2.49),RVHP(2.49),RGMP("2.49")); //Z-boson like Sherpa
multi_precision_constant constants::MW(80.419,RHP(80.419),RVHP(80.419),RGMP("80.419")); //W-boson
multi_precision_constant constants::GW(2.06,RHP(2.06),RVHP(2.06),RGMP("2.06")); //W-boson
multi_precision_constant constants::Mtop(131.1,RHP(131.1),RVHP(131.1),RGMP("131.1")); //top-mass
multi_precision_constant constants::G3_Lambda2(1e5,RHP(1e5),RVHP(1e5),RGMP("1e	5")); //tr(G^3) term coupling
#else
multi_precision_constant constants::MZ(91.188,RHP(91.188),RVHP(91.188)); //Z-boson like Sherpa
multi_precision_constant constants::GZ(2.49,RHP(2.49),RVHP(2.49)); //Z-boson like Sherpa
multi_precision_constant constants::MW(80.419,RHP(80.419),RVHP(80.419)); //W-boson
multi_precision_constant constants::GW(2.06,RHP(2.06),RVHP(2.06)); //W-boson
multi_precision_constant constants::Mtop(131.1,RHP(131.1),RVHP(131.1)); //top-mass
multi_precision_constant constants::G3_Lambda2(1e10,RHP(1e10),RVHP(1e10)); //tr(G^3) term coupling
//multi_precision_constant constants::G3_Lambda2(9,RHP(9),RVHP(9)); //tr(G^3) term coupling
#endif
*/

// PDG 2012 values
#if BH_USE_GMP
multi_precision_constant constants::MZ(91.1876,RHP(91.1876),RVHP(91.1876),RGMP("91.1876")); //Z-boson like Sherpa
multi_precision_constant constants::GZ(2.4952,RHP(2.4952),RVHP(2.4952),RGMP("2.4952")); //Z-boson like Sherpa
multi_precision_constant constants::MW(80.385,RHP(80.385),RVHP(80.385),RGMP("80.385")); //W-boson
multi_precision_constant constants::GW(2.085,RHP(2.085),RVHP(2.085),RGMP("2.085")); //W-boson
multi_precision_constant constants::Mtop(173.15,RHP(173.15),RVHP(173.15),RGMP("173.15")); //top-mass
multi_precision_constant constants::G3_Lambda2(1e5,RHP(1e5),RVHP(1e5),RGMP("1e	5")); //tr(G^3) term coupling
#else
multi_precision_constant constants::MZ(91.1876,RHP(91.1876),RVHP(91.1876)); //Z-boson like Sherpa
multi_precision_constant constants::GZ(2.4952,RHP(2.4952),RVHP(2.4952)); //Z-boson like Sherpa
multi_precision_constant constants::MW(80.385,RHP(80.385),RVHP(80.385)); //W-boson
multi_precision_constant constants::GW(2.085,RHP(2.085),RVHP(2.085)); //W-boson
multi_precision_constant constants::Mtop(173.15,RHP(173.15),RVHP(173.15)); //top-mass
multi_precision_constant constants::G3_Lambda2(1e10,RHP(1e10),RVHP(1e10)); //tr(G^3) term coupling
//multi_precision_constant constants::G3_Lambda2(9,RHP(9),RVHP(9)); //tr(G^3) term coupling
#endif
R constants::vel=(-1.+2.*sin_th_2)/(sin_2th);
R constants::ver=(2.*sin_th_2)/(sin_2th);
R constants::vnuel=1./(sin_2th);
R constants::vnuer=0.;
R constants::vupl=(1.-2.*2./3.*sin_th_2)/(sin_2th);
R constants::vupr=-(2.*2./3.*sin_th_2)/(sin_2th);
R constants::vdownl=(-1.+2.*1./3.*sin_th_2)/(sin_2th);
R constants::vdownr=(2.*1./3.*sin_th_2)/(sin_2th);
// default value 1. --- no change in units
R constants::s_GeV=1.;
R constants::alpha_S=0.118;//NLO alpha_S(MZ)
R constants::alpha_QED=1./128.802;




void constants::update_constants(settings_table* st){
	st->apply_setting("sin_th_2",constants::sin_th_2);
	//st->apply_setting("sin_2th",constants::sin_2th);
	constants::sin_2th= sin(2*asin(sqrt(constants::sin_th_2)));
	R MZ,GZ,MW,GW,Mtop,G3_Lambda2;
	st->apply_setting("alpha_S",alpha_S);
 	st->apply_setting("alpha_QED",alpha_QED);	
	st->apply_setting("Z_mass",MZ);
	st->apply_setting("W_mass",MW);
	st->apply_setting("Z_width",GZ);
	st->apply_setting("W_width",GW);
	st->apply_setting("Top_mass",Mtop);
    //
	st->apply_setting("G3_Lambda2",G3_Lambda2);
#if BH_USE_GMP
	constants::MZ= multi_precision_constant(MZ,RHP(MZ),RVHP(MZ),RGMP(MZ));
	constants::MW= multi_precision_constant(MW,RHP(MW),RVHP(MW),RGMP(MW));
	constants::GZ= multi_precision_constant(GZ,RHP(GZ),RVHP(GZ),RGMP(GZ));
	constants::GW= multi_precision_constant(GW,RHP(GW),RVHP(GW),RGMP(GW));
	constants::Mtop= multi_precision_constant(Mtop,RHP(Mtop),RVHP(Mtop),RGMP(Mtop));
	//
    constants::G3_Lambda2= multi_precision_constant(G3_Lambda2,RHP(G3_Lambda2),RVHP(G3_Lambda2),RGMP(G3_Lambda2));
#else
	constants::MZ= multi_precision_constant(MZ,RHP(MZ),RVHP(MZ));
	constants::MW= multi_precision_constant(MW,RHP(MW),RVHP(MW));
	constants::GZ= multi_precision_constant(GZ,RHP(GZ),RVHP(GZ));
	constants::GW= multi_precision_constant(GW,RHP(GW),RVHP(GW));
	constants::Mtop= multi_precision_constant(Mtop,RHP(Mtop),RVHP(Mtop));
    //
    constants::G3_Lambda2= multi_precision_constant(G3_Lambda2,RHP(G3_Lambda2),RVHP(G3_Lambda2));
#endif
	constants::vel=(-1.+2.*sin_th_2)/(sin_2th);
	constants::ver=(2.*sin_th_2)/(sin_2th);
	constants::vnuel=1./(sin_2th);
	constants::vnuer=0.;
	constants::vupl=(1.-2.*2./3.*sin_th_2)/(sin_2th);
	constants::vupr=-(2.*2./3.*sin_th_2)/(sin_2th);
	constants::vdownl=(-1.+2.*1./3.*sin_th_2)/(sin_2th);
	constants::vdownr=(2.*1./3.*sin_th_2)/(sin_2th);
	
	constants::alpha_S=alpha_S;
	constants::alpha_QED=alpha_QED;


}

long gcd(long a, long b){
	if (b == 0) return a;
	else return gcd(b, a % b);
}


void multi_precision_fraction::normalize(){
	long g=gcd(abs(m_num), abs(m_den));
	m_num/=g;
	m_den/=g;

};


}
