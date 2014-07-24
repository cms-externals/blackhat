/*!\file rat_ampl.h
\brief Header for the rational vertices of the recursion for the rational terms
*/
#ifndef RAT_AMP_H_
#define RAT_AMP_H_

#define X 2.
//#define _VERBOSE 0  // 0 is silent

#include "tree_amp.h"
namespace BH {


// We set all 3-pt one loop terms to zero
/*template <int i1, int i2, class T> std::complex<T> R3gm(momentum_configuration<T>& mc, const vector<int>& ind)
{
	return (std::complex<T>(0.,-1.)*pow(mc.spa(ind.at(i1),ind.at(i2)),4))/(mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(2),ind.at(0)));
}

template <int i1, int i2, class T> std::complex<T> R3gp(momentum_configuration<T>& mc, const vector<int>& ind)
{
	return (std::complex<T>(0.,-1.)*pow(mc.spb(ind.at(i1),ind.at(i2)),4))/(mc.spb(ind.at(0),ind.at(1))*mc.spb(ind.at(1),ind.at(2))*mc.spb(ind.at(2),ind.at(0)));
}*/

#define _ONLY_X_PART 1

template <class T> std::complex<T> R4g0(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R4 : ----");
#endif
	return std::complex<T>(0,0)
#if _ONLY_X_PART
	-std::complex<T>(0.,1.)/std::complex<T>(3.,0.)*mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(2),ind.at(3))/(mc.spb(ind.at(0),ind.at(1))*mc.spb(ind.at(2),ind.at(3)))
#endif
;
}

template <class T> std::complex<T> R4g15(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R4 : ++++");
#endif
	return std::complex<T>(0,0)
	#if _ONLY_X_PART
	-std::complex<T>(0.,1.)/std::complex<T>(3.,0.)*mc.spb(ind.at(0),ind.at(1))*mc.spb(ind.at(2),ind.at(3))/(mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(2),ind.at(3)))
#endif
;
}

template <int i1, int i2, class T> std::complex<T> R4g2m2p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R4 : --++");
#endif
	return std::complex<T>(X,0.)/std::complex<T>(9.,0.)*A04g<i1,i2>(mc,ind);
}

template <int i1, int i2, int i3, int i4, class T> std::complex<T> R4g2m2p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R4 : -+-+");
#endif
	return std::complex<T>(X,0.)/std::complex<T>(9.,0.)*A04g<i1,i3>(mc,ind)
#if _ONLY_X_PART
	-std::complex<T>(0.,1.)*pow(mc.spa(ind[i1],ind[i3]),2)*mc.spb(ind[i1],ind[i2])*mc.spb(ind[i2],ind[i3])/(mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i1])*pow(mc.spb(ind[i1],ind[i3]),2))
#endif
;
}

template <int i1, int i2, int i3, int i4, class T> std::complex<T> R4g1p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R4 : +---");
#endif
	return std::complex<T>(0,0)
#if _ONLY_X_PART
	+std::complex<T>(0.,1.)/std::complex<T>(3.,0.)*mc.spb(ind.at(i2),ind.at(i4))*pow(mc.spa(ind.at(i2),ind.at(i4)),3)/(mc.spa(ind.at(i1),ind.at(i2))*mc.spb(ind.at(i2),ind.at(i3))*mc.spb(ind.at(i3),ind.at(i4))*mc.spa(ind.at(i4),ind.at(i1)))
#endif
	;
}

template <int i1, int i2, int i3, int i4, class T> std::complex<T> R4g1m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R4 : -+++");
#endif

	return std::complex<T>(0,0)
	#if _ONLY_X_PART
	+ std::complex<T>(0.,1.)/std::complex<T>(3.,0.)*mc.spa(ind.at(i2),ind.at(i4))*pow(mc.spb(ind.at(i2),ind.at(i4)),3)/(mc.spb(ind.at(i1),ind.at(i2))*mc.spa(ind.at(i2),ind.at(i3))*mc.spa(ind.at(i3),ind.at(i4))*mc.spb(ind.at(i4),ind.at(i1)))
#endif
;
}



template <class T> std::complex<T> R5g0(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R5 : -----");
#endif
	return std::complex<T>(0,0)
	#if _ONLY_X_PART
	+ std::complex<T>(1.,0.)/std::complex<T>(0.,-3.)*(
			mc.spa(ind[0],ind[1])*mc.spa(ind[1],ind[2])/(mc.spb(ind[2],ind[3])*mc.spb(ind[3],ind[4])*mc.spb(ind[4],ind[0]))
			+mc.spa(ind[3],ind[4])*mc.spa(ind[4],ind[0])/(mc.spb(ind[0],ind[1])*mc.spb(ind[1],ind[2])*mc.spb(ind[2],ind[3]))
			+mc.spa(ind[1],ind[4])*mc.spa(ind[2],ind[3])/(mc.spb(ind[0],ind[1])*mc.spb(ind[2],ind[3])*mc.spb(ind[4],ind[0])))
#endif
;
}

template <class T> std::complex<T> R5g31(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R5 : +++++");
#endif
	return std::complex<T>(0,0)
	#if _ONLY_X_PART
	+ std::complex<T>(1.,0.)/std::complex<T>(0.,-3.)*(
			mc.spb(ind[0],ind[1])*mc.spb(ind[1],ind[2])/(mc.spa(ind[2],ind[3])*mc.spa(ind[3],ind[4])*mc.spa(ind[4],ind[0]))
			+mc.spb(ind[3],ind[4])*mc.spb(ind[4],ind[0])/(mc.spa(ind[0],ind[1])*mc.spa(ind[1],ind[2])*mc.spa(ind[2],ind[3]))
			+mc.spb(ind[1],ind[4])*mc.spb(ind[2],ind[3])/(mc.spa(ind[0],ind[1])*mc.spa(ind[2],ind[3])*mc.spa(ind[4],ind[0])))
#endif
	;
}

template <int i1, int i2, int i3, int i4, int i5, class T> std::complex<T> R5g1p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R5 : +----");
#endif
	return std::complex<T>(0,0)
	#if _ONLY_X_PART
	+ std::complex<T>(0.,1.)/(std::complex<T>(3.,0.)*pow(mc.spb(ind[i3],ind[i4]),2))*(
					-pow(mc.spa(ind[i2],ind[i5]),3)/(mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i1]))
					+pow(mc.spb(ind[i1],ind[i4]),3)*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i3],ind[i5])/(mc.spb(ind[i1],ind[i2])*mc.spb(ind[i2],ind[i3])*pow(mc.spb(ind[i4],ind[i5]),2))
					-pow(mc.spb(ind[i1],ind[i3]),3)*mc.spa(ind[i3],ind[i2])*mc.spb(ind[i4],ind[i2])/(mc.spb(ind[i1],ind[i5])*mc.spb(ind[i5],ind[i4])*pow(mc.spb(ind[i3],ind[i2]),2)))
#endif
	;
}

template <int i1, int i2, int i3, int i4, int i5, class T> std::complex<T> R5g1m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE6("R5 : -++++",i1,i2,i3,i4,i5);
#endif
	return std::complex<T>(0,0)
	#if _ONLY_X_PART
	+ std::complex<T>(0.,1.)/(std::complex<T>(3.,0.)*pow(mc.spa(ind[i3],ind[i4]),2))*(
					-pow(mc.spb(ind[i2],ind[i5]),3)/(mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i1]))
					+pow(mc.spa(ind[i1],ind[i4]),3)*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i3],ind[i5])/(mc.spa(ind[i1],ind[i2])*mc.spa(ind[i2],ind[i3])*pow(mc.spa(ind[i4],ind[i5]),2))
					-pow(mc.spa(ind[i1],ind[i3]),3)*mc.spb(ind[i3],ind[i2])*mc.spa(ind[i4],ind[i2])/(mc.spa(ind[i1],ind[i5])*mc.spa(ind[i5],ind[i4])*pow(mc.spa(ind[i3],ind[i2]),2)))

#endif
	;
}

template <int i1, int i2, int i3, int i4, int i5, class T> std::complex<T> R5g2m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
std::complex<T> result(0,0);
std::complex<T> Xpart=std::complex<T>(X,0.)/std::complex<T>(9.,0.)*A05g2m<i1,i2>(mc,ind);
	result+=Xpart;
#if _ONLY_X_PART

	std::complex<T> Rpart=-(std::complex<T>(0.,1.)/std::complex<T>(3.,0.)*(
				mc.spa(ind[i3],ind[i5])*pow(mc.spb(ind[i3],ind[i5]),3)/(mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1])*mc.spb(ind[i1],ind[i2])*mc.spb(ind[i2],ind[i3]))
				-pow(mc.spb(ind[i3],ind[i5]),2)*mc.spa(ind[i1],ind[i2])/(mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1])*mc.spb(ind[i2],ind[i3]))
				+std::complex<T>(1.,0.)/std::complex<T>(2.,0.)*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i1])*mc.spa(ind[i4],ind[i2])*mc.spa(ind[i1],ind[i2])/(mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.s(ind[i5],ind[i1])*mc.s(ind[i2],ind[i3]))));
				/*CRhat*/
std::complex<T> Cpart=-(-std::complex<T>(0.,1.)/std::complex<T>(6.,0.)*(mc.spb(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spb(ind[i4],ind[i5])*(mc.spa(ind[i2],ind[i3])*mc.spb(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i1])+mc.spa(ind[i2],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1]))
				/(mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5]))*(mc.s(ind[i2],ind[i3])/mc.s(ind[i5],ind[i1])-mc.s(ind[i5],ind[i1])/mc.s(ind[i2],ind[i3]))/pow(mc.s(ind[i5],ind[i1])-mc.s(ind[i2],ind[i3]),3)));
#if _VERBOSE
	_MESSAGE6("R5 : --+++",i1,i2,i3,i4,i5);
	_PRINT(Xpart);
	_PRINT(Rpart);
	_PRINT(Cpart);
#endif

	result+=Rpart;
	result+=Cpart;

	#endif
	;

	return result;

}

template <int i1, int i2, int i3, int i4, int i5, class T> std::complex<T> R5g2p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE6("R5 : ++---",i1,i2,i3,i4,i5);
#endif
	return 			+std::complex<T>(X,0.)/std::complex<T>(9.,0.)*A05g2p<i1,i2>(mc,ind)
#if _ONLY_X_PART


	+(std::complex<T>(0.,1.)/std::complex<T>(3.,0.)*(
				mc.spb(ind[i3],ind[i5])*pow(mc.spa(ind[i3],ind[i5]),3)/(mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1])*mc.spa(ind[i1],ind[i2])*mc.spa(ind[i2],ind[i3]))
				-pow(mc.spa(ind[i3],ind[i5]),2)*mc.spb(ind[i1],ind[i2])/(mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1])*mc.spa(ind[i2],ind[i3]))
				+std::complex<T>(1.,0.)/std::complex<T>(2.,0.)*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i1])*mc.spb(ind[i4],ind[i2])*mc.spb(ind[i1],ind[i2])/(mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.s(ind[i5],ind[i1])*mc.s(ind[i2],ind[i3])))
					/*CRhat*/
				-std::complex<T>(0.,1.)/std::complex<T>(6.,0.)*(mc.spa(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spa(ind[i4],ind[i5])*(mc.spb(ind[i2],ind[i3])*mc.spa(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i1])+mc.spb(ind[i2],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1]))
				/(mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5]))*(mc.s(ind[i2],ind[i3])/mc.s(ind[i5],ind[i1])-mc.s(ind[i5],ind[i1])/mc.s(ind[i2],ind[i3]))/pow(mc.s(ind[i5],ind[i1])-mc.s(ind[i2],ind[i3]),3)));
#endif
	;

}

template <int i1, int i2, int i3, int i4, int i5, class T> std::complex<T> R5gmpm(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE6("R5 : -+-++",i1,i2,i3,i4,i5);
#endif
	return (std::complex<T>(0.,1.)/std::complex<T>(3.,0.)*(
			std::complex<T>(X,0.)/std::complex<T>(3.,0.)*pow(mc.spa(ind[i1],ind[i3]),4)/(mc.spa(ind[i1],ind[i2])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1]))
#if _ONLY_X_PART
			+pow(mc.spb(ind[i2],ind[i4]),2)*pow(mc.spb(ind[i2],ind[i5]),2)/(mc.spb(ind[i1],ind[i2])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1]))
			-mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i1]),2)*pow(mc.spb(ind[i2],ind[i4]),3)/(mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i3],ind[i4])*mc.s(ind[i5],ind[i1]))
			+mc.spa(ind[i3],ind[i2])*pow(mc.spa(ind[i5],ind[i3]),2)*pow(mc.spb(ind[i2],ind[i5]),3)/(mc.spa(ind[i5],ind[i4])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i2],ind[i5])*mc.spb(ind[i2],ind[i1])*mc.spb(ind[i1],ind[i5])*mc.s(ind[i3],ind[i4]))
			+std::complex<T>(1.,0.)/std::complex<T>(2.,0.)*pow(mc.spa(ind[i1],ind[i3]),2)*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i2],ind[i5])/(mc.s(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.s(ind[i5],ind[i1]))
#endif

	)
#if _ONLY_X_PART
	/*CRhat*/
			+std::complex<T>(0.,1.)*(-(((T(1.)/(T(1.) - mc.s(ind[i2],ind[i3])/mc.s(ind[i5],ind[i1])) + T(1.)/(T(1.) - mc.s(ind[i3],ind[i4])/mc.s(ind[i5],ind[i1])))*mc.spa(ind[i1],ind[i2])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i3],ind[i4])*
	        pow(mc.spa(ind[i4],ind[i1]),2)*pow(mc.spb(ind[i2],ind[i4]),2))/
	      (pow(mc.s(ind[i5],ind[i1]),2)*pow(mc.spa(ind[i2],ind[i4]),2)*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1]))) -
	   ((mc.s(ind[i2],ind[i3])/mc.s(ind[i5],ind[i1]) - mc.s(ind[i5],ind[i1])/mc.s(ind[i2],ind[i3]))*pow(mc.spa(ind[i2],ind[i3]),2)*pow(mc.spa(ind[i4],ind[i1]),3)*
	      pow(mc.spb(ind[i2],ind[i4]),3))/(T(3.)*pow(mc.s(ind[i5],ind[i1]) - mc.s(ind[i2],ind[i3]),3)*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1])) +
	   ((T(1.)/(T(1.) - mc.s(ind[i1],ind[i2])/mc.s(ind[i3],ind[i4])) + T(1.)/(T(1.) - mc.s(ind[i5],ind[i1])/mc.s(ind[i3],ind[i4])))*mc.spa(ind[i1],ind[i5])*mc.spa(ind[i2],ind[i1])*mc.spa(ind[i3],ind[i2])*
	      pow(mc.spa(ind[i5],ind[i3]),2)*pow(mc.spb(ind[i2],ind[i5]),2))/(pow(mc.s(ind[i3],ind[i4]),2)*pow(mc.spa(ind[i2],ind[i5]),2)*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i5],ind[i4]))
	     + ((mc.s(ind[i1],ind[i2])/mc.s(ind[i3],ind[i4]) - mc.s(ind[i3],ind[i4])/mc.s(ind[i1],ind[i2]))*pow(mc.spa(ind[i2],ind[i1]),2)*pow(mc.spa(ind[i5],ind[i3]),3)*
	      pow(mc.spb(ind[i2],ind[i5]),3))/(T(3.)*pow(mc.s(ind[i3],ind[i4]) - mc.s(ind[i1],ind[i2]),3)*mc.spa(ind[i2],ind[i5])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i5],ind[i4])) -
	   ((mc.s(ind[i3],ind[i4])/mc.s(ind[i5],ind[i1]) - mc.s(ind[i5],ind[i1])/mc.s(ind[i3],ind[i4]))*
	      ((T(2.)*pow(mc.spa(ind[i1],ind[i2]),2)*pow(mc.spa(ind[i3],ind[i4]),2)*mc.spa(ind[i4],ind[i1])*pow(mc.spb(ind[i2],ind[i4]),3))/
	         (T(3.)*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1])) -
	        (T(2.)*pow(mc.spa(ind[i1],ind[i5]),2)*pow(mc.spa(ind[i3],ind[i2]),2)*mc.spa(ind[i5],ind[i3])*pow(mc.spb(ind[i2],ind[i5]),3))/
	         (T(3.)*mc.spa(ind[i2],ind[i5])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i5],ind[i4])) +
	        (mc.spa(ind[i1],ind[i3])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i2],ind[i5])*
	           (-(mc.spa(ind[i2],ind[i1])*mc.spa(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i2])) + mc.spa(ind[i1],ind[i5])*mc.spa(ind[i2],ind[i3])*mc.spb(ind[i5],ind[i2])))/(T(3.)*mc.spa(ind[i4],ind[i5]))))/
	    (T(2.)*pow(mc.s(ind[i5],ind[i1]) - mc.s(ind[i3],ind[i4]),3)))
#endif
	   ) ;
}

template <int i1, int i2, int i3, int i4, int i5, class T> std::complex<T> R5gpmp(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE6("R5 : +-+--",i1,i2,i3,i4,i5);
#endif
	return -(std::complex<T>(0.,1.)/std::complex<T>(3.,0.)*(
			std::complex<T>(X,0.)/std::complex<T>(3.,0.)*pow(mc.spb(ind[i1],ind[i3]),4)/(mc.spb(ind[i1],ind[i2])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1]))
#if _ONLY_X_PART
			+pow(mc.spa(ind[i2],ind[i4]),2)*pow(mc.spa(ind[i2],ind[i5]),2)/(mc.spa(ind[i1],ind[i2])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i1]))
			-mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i1]),2)*pow(mc.spa(ind[i2],ind[i4]),3)/(mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i3],ind[i4])*mc.s(ind[i5],ind[i1]))
			+mc.spb(ind[i3],ind[i2])*pow(mc.spb(ind[i5],ind[i3]),2)*pow(mc.spa(ind[i2],ind[i5]),3)/(mc.spb(ind[i5],ind[i4])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i2],ind[i5])*mc.spa(ind[i2],ind[i1])*mc.spa(ind[i1],ind[i5])*mc.s(ind[i3],ind[i4]))
			+std::complex<T>(1.,0.)/std::complex<T>(2.,0.)*pow(mc.spb(ind[i1],ind[i3]),2)*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i2],ind[i5])/(mc.s(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.s(ind[i5],ind[i1]))
#endif
			)
#if _ONLY_X_PART
		/*CRhat*/
			+std::complex<T>(0.,1.)*(-(((T(1.)/(T(1.) - mc.s(ind[i2],ind[i3])/mc.s(ind[i5],ind[i1])) + T(1.)/(T(1.) - mc.s(ind[i3],ind[i4])/mc.s(ind[i5],ind[i1])))*mc.spb(ind[i1],ind[i2])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i3],ind[i4])*
		        pow(mc.spb(ind[i4],ind[i1]),2)*pow(mc.spa(ind[i2],ind[i4]),2))/
		      (pow(mc.s(ind[i5],ind[i1]),2)*pow(mc.spb(ind[i2],ind[i4]),2)*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1]))) -
		   ((mc.s(ind[i2],ind[i3])/mc.s(ind[i5],ind[i1]) - mc.s(ind[i5],ind[i1])/mc.s(ind[i2],ind[i3]))*pow(mc.spb(ind[i2],ind[i3]),2)*pow(mc.spb(ind[i4],ind[i1]),3)*
		      pow(mc.spa(ind[i2],ind[i4]),3))/(T(3.)*pow(mc.s(ind[i5],ind[i1]) - mc.s(ind[i2],ind[i3]),3)*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1])) +
		   ((T(1.)/(T(1.) - mc.s(ind[i1],ind[i2])/mc.s(ind[i3],ind[i4])) + T(1.)/(T(1.) - mc.s(ind[i5],ind[i1])/mc.s(ind[i3],ind[i4])))*mc.spb(ind[i1],ind[i5])*mc.spb(ind[i2],ind[i1])*mc.spb(ind[i3],ind[i2])*
		      pow(mc.spb(ind[i5],ind[i3]),2)*pow(mc.spa(ind[i2],ind[i5]),2))/(pow(mc.s(ind[i3],ind[i4]),2)*pow(mc.spb(ind[i2],ind[i5]),2)*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i5],ind[i4]))
		     + ((mc.s(ind[i1],ind[i2])/mc.s(ind[i3],ind[i4]) - mc.s(ind[i3],ind[i4])/mc.s(ind[i1],ind[i2]))*pow(mc.spb(ind[i2],ind[i1]),2)*pow(mc.spb(ind[i5],ind[i3]),3)*
		      pow(mc.spa(ind[i2],ind[i5]),3))/(T(3.)*pow(mc.s(ind[i3],ind[i4]) - mc.s(ind[i1],ind[i2]),3)*mc.spb(ind[i2],ind[i5])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i5],ind[i4])) -
		   ((mc.s(ind[i3],ind[i4])/mc.s(ind[i5],ind[i1]) - mc.s(ind[i5],ind[i1])/mc.s(ind[i3],ind[i4]))*
		      ((T(2.)*pow(mc.spb(ind[i1],ind[i2]),2)*pow(mc.spb(ind[i3],ind[i4]),2)*mc.spb(ind[i4],ind[i1])*pow(mc.spa(ind[i2],ind[i4]),3))/
		         (T(3.)*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i1])) -
		        (T(2.)*pow(mc.spb(ind[i1],ind[i5]),2)*pow(mc.spb(ind[i3],ind[i2]),2)*mc.spb(ind[i5],ind[i3])*pow(mc.spa(ind[i2],ind[i5]),3))/
		         (T(3.)*mc.spb(ind[i2],ind[i5])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i5],ind[i4])) +
		        (mc.spb(ind[i1],ind[i3])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i2],ind[i5])*
		           (-(mc.spb(ind[i2],ind[i1])*mc.spb(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i2])) + mc.spb(ind[i1],ind[i5])*mc.spb(ind[i2],ind[i3])*mc.spa(ind[i5],ind[i2])))/(T(3.)*mc.spb(ind[i4],ind[i5]))))/
		    (T(2.)*pow(mc.s(ind[i5],ind[i1]) - mc.s(ind[i3],ind[i4]),3)))
#endif
    );
}


// added by Daniel
template <int i1, int i2, int i3, int i4, int i5, int i6, class T> std::complex<T> R6g1m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R6 : -+++++");
#endif
	return std::complex<T>(0,0)
	#if _ONLY_X_PART
+std::complex<T>(0,1)/std::complex<T>(3,0)*
		((pow(mc.spa(ind.at(i1),ind.at(i4)),3)
		*mc.spa(ind.at(i3),ind.at(i5))*(mc.spab(ind.at(i1),ind.at(i2),ind.at(i4)) + mc.spab(ind.at(i1),ind.at(i3),ind.at(i4))))/(mc.spa(ind.at(i1),ind.at(i2))*mc.spa(ind.at(i1),ind.at(i6))*mc.spa(ind.at(i2),ind.at(i3))*pow(mc.spa(ind.at(i3),ind.at(i4)),2)*pow(mc.spa(ind.at(i4),ind.at(i5)),2)*mc.spa(ind.at(i5),ind.at(i6))) + pow(mc.spab(ind.at(i1),ind.at(i2),ind.at(i6)) + 	mc.spab(ind.at(i1),ind.at(i3),ind.at(i6)),3)/(mc.s(ind.at(i1),ind.at(i2),ind.at(i3))*mc.spa(ind.at(i1),ind.at(i2))*mc.spa(ind.at(i2),ind.at(i3))*pow(mc.spa(ind.at(i4),ind.at(i5)),2)*(mc.spab(ind.at(i3),ind.at(i1),ind.at(i6)) +mc.spab(ind.at(i3),ind.at(i2),ind.at(i6)))) - pow(mc.spab(ind.at(i1),ind.at(i3),ind.at(i2)) + mc.spab(ind.at(i1),ind.at(i4),ind.at(i2)),3)/(mc.s(ind.at(i2),ind.at(i3),ind.at(i4))*mc.spa(ind.at(i1),ind.at(i6))*pow(mc.spa(ind.at(i3),ind.at(i4)),2)*mc.spa(ind.at(i5),ind.at(i6))*(mc.spab(ind.at(i5),ind.at(i3),ind.at(i2)) +mc.spab(ind.at(i5),ind.at(i4),ind.at(i2)))) - (pow(mc.spa(ind.at(i1),ind.at(i3)),3)*mc.spa(ind.at(i2),ind.at(i4))*mc.spb(ind.at(i3),ind.at(i2)))/(mc.spa(ind.at(i1),ind.at(i6))*pow(mc.spa(ind.at(i2),ind.at(i3)),2)*pow(mc.spa(ind.at(i3),ind.at(i4)),2)*mc.spa(ind.at(i4),ind.at(i5))*mc.spa(ind.at(i5),ind.at(i6))) - (pow(mc.spa(ind.at(i1),ind.at(i5)),3)*mc.spa(ind.at(i4),ind.at(i6))*mc.spb(ind.at(i6),ind.at(i5)))/(mc.spa(ind.at(i1),ind.at(i2))*mc.spa(ind.at(i2),ind.at(i3))*mc.spa(ind.at(i3),ind.at(i4))*pow(mc.spa(ind.at(i4),ind.at(i5)),2)*pow(mc.spa(ind.at(i5),ind.at(i6)),2)) + (pow(mc.spb(ind.at(i6),ind.at(i2)),3)*((mc.spb(ind.at(i3),ind.at(i2))*mc.spb(ind.at(i4),ind.at(i3)))/(mc.spa(ind.at(i4),ind.at(i5))*(mc.spab(ind.at(i5),ind.at(i3),ind.at(i2)) + mc.spab(ind.at(i5),ind.at(i4),ind.at(i2)))) - mc.spb(ind.at(i5),ind.at(i3))/(mc.spa(ind.at(i3),ind.at(i4))*mc.spa(ind.at(i4),ind.at(i5))) - (mc.spb(ind.at(i5),ind.at(i4))*mc.spb(ind.at(i6),ind.at(i5)))/(mc.spa(ind.at(i3),ind.at(i4))*(mc.spab(ind.at(i3),ind.at(i1),ind.at(i6)) + mc.spab(ind.at(i3),ind.at(i2),ind.at(i6))))))/(mc.s(ind.at(i3),ind.at(i4),ind.at(i5))*mc.spb(ind.at(i2),ind.at(i1))*mc.spb(ind.at(i6),ind.at(i1)))
		);

		#endif
;
}

// added by Daniel
template <class T> std::complex<T> R6g6p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
#if _VERBOSE
	_MESSAGE("R6 : ++++++");
#endif
	return std::complex<T>(0,0)
	#if _ONLY_X_PART
+(T(1)/std::complex<T>(0,-3)*(
		- mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(2),ind.at(3))*mc.spb(ind.at(2),ind.at(1))*mc.spb(ind.at(3),ind.at(0))
		- mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(2),ind.at(4))*mc.spb(ind.at(2),ind.at(1))*mc.spb(ind.at(4),ind.at(0))
		- mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(3),ind.at(4))*mc.spb(ind.at(3),ind.at(1))*mc.spb(ind.at(4),ind.at(0))
		- mc.spa(ind.at(0),ind.at(2))*mc.spa(ind.at(3),ind.at(4))*mc.spb(ind.at(3),ind.at(2))*mc.spb(ind.at(4),ind.at(0))
		- mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(3),ind.at(4))*mc.spb(ind.at(3),ind.at(2))*mc.spb(ind.at(4),ind.at(1))
		- mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(2),ind.at(5))*mc.spb(ind.at(2),ind.at(1))*mc.spb(ind.at(5),ind.at(0))
		- mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(3),ind.at(5))*mc.spb(ind.at(3),ind.at(1))*mc.spb(ind.at(5),ind.at(0))
		- mc.spa(ind.at(0),ind.at(2))*mc.spa(ind.at(3),ind.at(5))*mc.spb(ind.at(3),ind.at(2))*mc.spb(ind.at(5),ind.at(0))
		- mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(4),ind.at(5))*mc.spb(ind.at(4),ind.at(1))*mc.spb(ind.at(5),ind.at(0))
		- mc.spa(ind.at(0),ind.at(2))*mc.spa(ind.at(4),ind.at(5))*mc.spb(ind.at(4),ind.at(2))*mc.spb(ind.at(5),ind.at(0))
		- mc.spa(ind.at(0),ind.at(3))*mc.spa(ind.at(4),ind.at(5))*mc.spb(ind.at(4),ind.at(3))*mc.spb(ind.at(5),ind.at(0))
		- mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(3),ind.at(5))*mc.spb(ind.at(3),ind.at(2))*mc.spb(ind.at(5),ind.at(1))
		- mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(4),ind.at(5))*mc.spb(ind.at(4),ind.at(2))*mc.spb(ind.at(5),ind.at(1))
		- mc.spa(ind.at(1),ind.at(3))*mc.spa(ind.at(4),ind.at(5))*mc.spb(ind.at(4),ind.at(3))*mc.spb(ind.at(5),ind.at(1))
		- mc.spa(ind.at(2),ind.at(3))*mc.spa(ind.at(4),ind.at(5))*mc.spb(ind.at(4),ind.at(3))*mc.spb(ind.at(5),ind.at(2))
		))/(mc.spa(ind.at(0),ind.at(1))*mc.spa(ind.at(0),ind.at(5))*mc.spa(ind.at(1),ind.at(2))*mc.spa(ind.at(2),ind.at(3))
				*mc.spa(ind.at(3),ind.at(4))*mc.spa(ind.at(4),ind.at(5)))
#endif
	;

}















template <int i1, int i2, int i3, int i4, int i5, int i6, class T> std::complex<T> R6g2m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
  return  std::complex<T>(0,1.)*((std::complex<T>(-X,0.)/std::complex<T>(9.,0.)*pow(mc.spa(ind[i1],ind[i2]),3))/(mc.spa(ind[i1],ind[i6])*mc.spa(ind[i2],ind[i3])*
mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6]))
#if _ONLY_X_PART
+
(std::complex<T>(1.,0.)/std::complex<T>(3.,0.)*pow(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i2]))+mc.spa(ind[i1],ind[i4])*mc.spb(ind[i4],ind[i3]),3)*
mc.spa(ind[i3],ind[i5]))/(mc.s(ind[i2],ind[i3],ind[i4])*mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])*
mc.spb(ind[i3],ind[i2])*(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2]))-mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(6.,0.)*mc.spa(ind[i1],ind[i2])*(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i2]))+mc.spa(ind[i1],ind[i4])*mc.spb(ind[i4],ind[i3]))*
(mc.spa(ind[i1],ind[i4])*mc.spb(ind[i4],ind[i3])+std::complex<T>(2.,0)*(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i2])) +
   mc.spa(ind[i1],ind[i4])*mc.spb(ind[i4],ind[i3]))))/(mc.s(ind[i2],ind[i3],ind[i4])*mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*
mc.spa(ind[i5],ind[i6])*mc.spb(ind[i3],ind[i2]))+(std::complex<T>(-1./3.,0)*pow(mc.spb(ind[i6],ind[i3]),3))/
 (pow(mc.spa(ind[i4],ind[i5]),2)*mc.spb(ind[i2],ind[i1])*mc.spb(ind[i3],ind[i2])*mc.spb(ind[i6],ind[i1])) +
(std::complex<T>(1.,0.)/std::complex<T>(3.,0.)*pow(mc.spa(ind[i1],ind[i2]),3)*pow(mc.spb(ind[i6],ind[i4]),2)*(mc.s(ind[i4],ind[i5])+mc.s(ind[i5],ind[i6])))/
 (mc.s(ind[i1],ind[i2],ind[i3])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])*(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2])) -
 mc.spa(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i3]))*(mc.spa(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i3])*mc.spb(ind[i6],ind[i2]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(3.,0.)*mc.s(ind[i3],ind[i5])*(mc.spa(ind[i1],ind[i4])*mc.spb(ind[i3],ind[i1])+mc.spa(ind[i2],ind[i4])*mc.spb(ind[i3],ind[i2]))*
(mc.spa(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i4])*mc.spb(ind[i6],ind[i2]))*(mc.spa(ind[i1],ind[i5])*mc.spb(ind[i6],ind[i1]) +
 mc.spa(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i2])))/(pow(mc.spa(ind[i3],ind[i4]),2)*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spb(ind[i2],ind[i1])*
(mc.spa(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i1])+mc.spa(ind[i2],ind[i6])*mc.spb(ind[i3],ind[i2]))*(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2])) -
 mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))*mc.spb(ind[i6],ind[i1])) +
(std::complex<T>(-1.,0.)/std::complex<T>(3.,0.)*pow(mc.spa(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i4])*mc.spb(ind[i6],ind[i2]),2)*mc.spa(ind[i3],ind[i5])*
mc.spb(ind[i6],ind[i3]))/(pow(mc.spa(ind[i3],ind[i4]),2)*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spb(ind[i2],ind[i1])*
(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2]))-mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))*mc.spb(ind[i6],ind[i1])) +
(std::complex<T>(1.,0.)/std::complex<T>(3.,0.)*pow(mc.spa(ind[i1],ind[i5]),2)*pow(mc.spb(ind[i4],ind[i3]),2)*
(-(mc.spa(ind[i1],ind[i6])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i3]))-mc.spa(ind[i5],ind[i6])*
  (-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i2]))+mc.spa(ind[i1],ind[i4])*mc.spb(ind[i4],ind[i3])))*mc.spb(ind[i6],ind[i5]))/
 (pow(mc.spa(ind[i5],ind[i6]),2)*mc.s(ind[i2],ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i3],ind[i2])*
(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2]))-mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))*
(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2]))-mc.spa(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i3]))) +
(std::complex<T>(1.,0.)/std::complex<T>(6.,0.)*pow(mc.spa(ind[i1],ind[i5])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i2]),2)*
(mc.s(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])+std::complex<T>(-2.,0)*(mc.s(ind[i1],ind[i5])*mc.spa(ind[i4],ind[i5]) +
   mc.s(ind[i2],ind[i5])*mc.spa(ind[i4],ind[i5])-mc.spa(ind[i1],ind[i5])*mc.spa(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i1]) -
   mc.spa(ind[i2],ind[i5])*mc.spa(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i2])))*mc.spb(ind[i6],ind[i5]))/
 (pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i2],ind[i1])*
(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2]))-mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))*mc.spb(ind[i6],ind[i1])*
(mc.spa(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i3])*mc.spb(ind[i6],ind[i2]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(6.,0.)*pow(mc.spa(ind[i1],ind[i2]),3)*mc.spa(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i4])*mc.spb(ind[i6],ind[i5]))/
 (mc.spa(ind[i2],ind[i3])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])*(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2])) -
 mc.spa(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i3]))*(mc.spa(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i3])*mc.spb(ind[i6],ind[i2]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(6.,0.)*mc.spa(ind[i1],ind[i2])*mc.spa(ind[i1],ind[i5])*mc.spb(ind[i4],ind[i3])*(mc.s(ind[i3],ind[i5])*mc.spa(ind[i1],ind[i5]) +
 mc.s(ind[i4],ind[i5])*mc.spa(ind[i1],ind[i5])+mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i3]) +
 mc.spa(ind[i1],ind[i6])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i6],ind[i4]))*mc.spb(ind[i6],ind[i5]))/(mc.s(ind[i2],ind[i3],ind[i4])*mc.spa(ind[i3],ind[i4])*
mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])*(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2]))-mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))*
(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2]))-mc.spa(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i3]))) +
(std::complex<T>(1.,0.)/std::complex<T>(3.,0.)*pow(mc.spa(ind[i1],ind[i2]),3)*pow(mc.spb(ind[i5],ind[i3]),2)*(-mc.s(ind[i3],ind[i4])-mc.s(ind[i4],ind[i5])))/
 (mc.s(ind[i3],ind[i4],ind[i5])*mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*(mc.spa(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i1]) +
 mc.spa(ind[i2],ind[i6])*mc.spb(ind[i3],ind[i2]))*(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i1])+mc.spa(ind[i2],ind[i6])*mc.spb(ind[i6],ind[i5]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(3.,0.)*pow(mc.spa(ind[i1],ind[i2]),2)*pow(mc.spb(ind[i5],ind[i3]),2)*(-(mc.s(ind[i3],ind[i5])*mc.spa(ind[i1],ind[i5])) +
 mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2])+mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2])))/
 (mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*(mc.spa(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i1]) +
 mc.spa(ind[i2],ind[i6])*mc.spb(ind[i3],ind[i2]))*(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2]))-mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))*
(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i1])+mc.spa(ind[i2],ind[i6])*mc.spb(ind[i6],ind[i5]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(3.,0.)*pow(mc.spb(ind[i5],ind[i3]),2)*mc.spa(ind[i1],ind[i2])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i3],ind[i5])*
(mc.spa(ind[i1],ind[i5])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i2]))*mc.spb(ind[i6],ind[i5]))/
 (pow(mc.spa(ind[i3],ind[i4]),2)*mc.spa(ind[i4],ind[i5])*(mc.spa(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i1]) +
 mc.spa(ind[i2],ind[i6])*mc.spb(ind[i3],ind[i2]))*(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2]))-mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))*
mc.spb(ind[i6],ind[i1])*(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i1])+mc.spa(ind[i2],ind[i6])*mc.spb(ind[i6],ind[i5]))) +
(std::complex<T>(1.,0.)/std::complex<T>(6.,0.)*pow(mc.spa(ind[i1],ind[i2]),2)*(-((mc.spa(ind[i1],ind[i4])*mc.spb(ind[i4],ind[i3]))/mc.spb(ind[i3],ind[i2])) -
 (mc.spa(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i5]))/mc.spb(ind[i6],ind[i1])))/(mc.spa(ind[i1],ind[i6])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i3],ind[i4])*
mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])) +
(std::complex<T>(-1.,0.)/std::complex<T>(6.,0.)*(-(((-pow(mc.s(ind[i2],ind[i3],ind[i4]),2)+pow(mc.spa(ind[i2],ind[i3]),2)*
pow(mc.spb(ind[i3],ind[i2]),2))*mc.spa(ind[i1],ind[i4])*mc.spa(ind[i2],ind[i4])*mc.spb(ind[i4],ind[i3])*
(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2]))-mc.spa(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i3]))*
(-(mc.spa(ind[i1],ind[i4])*mc.spa(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3]))+mc.spa(ind[i2],ind[i4])*
(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2]))-mc.spa(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i3]))))/
   (pow(-mc.s(ind[i2],ind[i4])-mc.s(ind[i3],ind[i4]),3)*mc.s(ind[i2],ind[i3],ind[i4])*mc.spb(ind[i3],ind[i2]))) -
 (((pow(mc.spa(ind[i1],ind[i6]),2)*pow(mc.spb(ind[i6],ind[i1]),2))/mc.s(ind[i2],ind[i3],ind[i4])-mc.s(ind[i2],ind[i3],ind[i4]))*
   mc.spa(ind[i1],ind[i5])*mc.spa(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i5])*(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i1]) +
mc.spa(ind[i2],ind[i6])*mc.spb(ind[i6],ind[i5]))*(mc.spa(ind[i1],ind[i6])*mc.spa(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i5]) +
mc.spa(ind[i1],ind[i5])*(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i1])+mc.spa(ind[i2],ind[i6])*mc.spb(ind[i6],ind[i5]))))/
  (pow(mc.s(ind[i1],ind[i6])-mc.s(ind[i2],ind[i3],ind[i4]),3)*mc.spb(ind[i6],ind[i1]))))/(mc.spa(ind[i1],ind[i6])*mc.spa(ind[i2],ind[i3])*
mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6]))
#endif
);


}

template <int i1, int i2, int i3, int i4, int i5, int i6, class T> std::complex<T> R6g2p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
  return  std::complex<T>(0,1.)*((std::complex<T>(-X,0.)/std::complex<T>(9.,0.)*pow(mc.spb(ind[i1],ind[i2]),3))/(mc.spb(ind[i1],ind[i6])*mc.spb(ind[i2],ind[i3])*
mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6]))
#if _ONLY_X_PART
+
(std::complex<T>(1.,0.)/std::complex<T>(3.,0.)*pow(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i2]))+mc.spb(ind[i1],ind[i4])*mc.spa(ind[i4],ind[i3]),3)*
mc.spb(ind[i3],ind[i5]))/(mc.s(ind[i2],ind[i3],ind[i4])*mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])*
mc.spa(ind[i3],ind[i2])*(-(mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2]))-mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(6.,0.)*mc.spb(ind[i1],ind[i2])*(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i2]))+mc.spb(ind[i1],ind[i4])*mc.spa(ind[i4],ind[i3]))*
(mc.spb(ind[i1],ind[i4])*mc.spa(ind[i4],ind[i3])+std::complex<T>(2.,0)*(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i2])) +
   mc.spb(ind[i1],ind[i4])*mc.spa(ind[i4],ind[i3]))))/(mc.s(ind[i2],ind[i3],ind[i4])*mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*
mc.spb(ind[i5],ind[i6])*mc.spa(ind[i3],ind[i2]))+(std::complex<T>(-1./3.,0)*pow(mc.spa(ind[i6],ind[i3]),3))/
 (pow(mc.spb(ind[i4],ind[i5]),2)*mc.spa(ind[i2],ind[i1])*mc.spa(ind[i3],ind[i2])*mc.spa(ind[i6],ind[i1])) +
(std::complex<T>(1.,0.)/std::complex<T>(3.,0.)*pow(mc.spb(ind[i1],ind[i2]),3)*pow(mc.spa(ind[i6],ind[i4]),2)*(mc.s(ind[i4],ind[i5])+mc.s(ind[i5],ind[i6])))/
 (mc.s(ind[i1],ind[i2],ind[i3])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])*(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i2])) -
 mc.spb(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i3]))*(mc.spb(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i1])+mc.spb(ind[i2],ind[i3])*mc.spa(ind[i6],ind[i2]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(3.,0.)*mc.s(ind[i3],ind[i5])*(mc.spb(ind[i1],ind[i4])*mc.spa(ind[i3],ind[i1])+mc.spb(ind[i2],ind[i4])*mc.spa(ind[i3],ind[i2]))*
(mc.spb(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i1])+mc.spb(ind[i2],ind[i4])*mc.spa(ind[i6],ind[i2]))*(mc.spb(ind[i1],ind[i5])*mc.spa(ind[i6],ind[i1]) +
 mc.spb(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i2])))/(pow(mc.spb(ind[i3],ind[i4]),2)*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spa(ind[i2],ind[i1])*
(mc.spb(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i1])+mc.spb(ind[i2],ind[i6])*mc.spa(ind[i3],ind[i2]))*(-(mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2])) -
 mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2]))*mc.spa(ind[i6],ind[i1])) +
(std::complex<T>(-1.,0.)/std::complex<T>(3.,0.)*pow(mc.spb(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i1])+mc.spb(ind[i2],ind[i4])*mc.spa(ind[i6],ind[i2]),2)*mc.spb(ind[i3],ind[i5])*
mc.spa(ind[i6],ind[i3]))/(pow(mc.spb(ind[i3],ind[i4]),2)*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spa(ind[i2],ind[i1])*
(-(mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2]))-mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2]))*mc.spa(ind[i6],ind[i1])) +
(std::complex<T>(1.,0.)/std::complex<T>(3.,0.)*pow(mc.spb(ind[i1],ind[i5]),2)*pow(mc.spa(ind[i4],ind[i3]),2)*
(-(mc.spb(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i3]))-mc.spb(ind[i5],ind[i6])*
  (-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i2]))+mc.spb(ind[i1],ind[i4])*mc.spa(ind[i4],ind[i3])))*mc.spa(ind[i6],ind[i5]))/
 (pow(mc.spb(ind[i5],ind[i6]),2)*mc.s(ind[i2],ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i3],ind[i2])*
(-(mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2]))-mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2]))*
(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i2]))-mc.spb(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i3]))) +
(std::complex<T>(1.,0.)/std::complex<T>(6.,0.)*pow(mc.spb(ind[i1],ind[i5])*mc.spa(ind[i6],ind[i1])+mc.spb(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i2]),2)*
(mc.s(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])+std::complex<T>(-2.,0)*(mc.s(ind[i1],ind[i5])*mc.spb(ind[i4],ind[i5]) +
   mc.s(ind[i2],ind[i5])*mc.spb(ind[i4],ind[i5])-mc.spb(ind[i1],ind[i5])*mc.spb(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i1]) -
   mc.spb(ind[i2],ind[i5])*mc.spb(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i2])))*mc.spa(ind[i6],ind[i5]))/
 (pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i2],ind[i1])*
(-(mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2]))-mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2]))*mc.spa(ind[i6],ind[i1])*
(mc.spb(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i1])+mc.spb(ind[i2],ind[i3])*mc.spa(ind[i6],ind[i2]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(6.,0.)*pow(mc.spb(ind[i1],ind[i2]),3)*mc.spb(ind[i3],ind[i5])*mc.spa(ind[i6],ind[i4])*mc.spa(ind[i6],ind[i5]))/
 (mc.spb(ind[i2],ind[i3])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])*(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i2])) -
 mc.spb(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i3]))*(mc.spb(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i1])+mc.spb(ind[i2],ind[i3])*mc.spa(ind[i6],ind[i2]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(6.,0.)*mc.spb(ind[i1],ind[i2])*mc.spb(ind[i1],ind[i5])*mc.spa(ind[i4],ind[i3])*(mc.s(ind[i3],ind[i5])*mc.spb(ind[i1],ind[i5]) +
 mc.s(ind[i4],ind[i5])*mc.spb(ind[i1],ind[i5])+mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i5])*mc.spa(ind[i6],ind[i3]) +
 mc.spb(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i6],ind[i4]))*mc.spa(ind[i6],ind[i5]))/(mc.s(ind[i2],ind[i3],ind[i4])*mc.spb(ind[i3],ind[i4])*
mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])*(-(mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2]))-mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2]))*
(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i2]))-mc.spb(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i3]))) +
(std::complex<T>(1.,0.)/std::complex<T>(3.,0.)*pow(mc.spb(ind[i1],ind[i2]),3)*pow(mc.spa(ind[i5],ind[i3]),2)*(-mc.s(ind[i3],ind[i4])-mc.s(ind[i4],ind[i5])))/
 (mc.s(ind[i3],ind[i4],ind[i5])*mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*(mc.spb(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i1]) +
 mc.spb(ind[i2],ind[i6])*mc.spa(ind[i3],ind[i2]))*(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i1])+mc.spb(ind[i2],ind[i6])*mc.spa(ind[i6],ind[i5]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(3.,0.)*pow(mc.spb(ind[i1],ind[i2]),2)*pow(mc.spa(ind[i5],ind[i3]),2)*(-(mc.s(ind[i3],ind[i5])*mc.spb(ind[i1],ind[i5])) +
 mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2])+mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2])))/
 (mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*(mc.spb(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i1]) +
 mc.spb(ind[i2],ind[i6])*mc.spa(ind[i3],ind[i2]))*(-(mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2]))-mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2]))*
(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i1])+mc.spb(ind[i2],ind[i6])*mc.spa(ind[i6],ind[i5]))) +
(std::complex<T>(-1.,0.)/std::complex<T>(3.,0.)*pow(mc.spa(ind[i5],ind[i3]),2)*mc.spb(ind[i1],ind[i2])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i3],ind[i5])*
(mc.spb(ind[i1],ind[i5])*mc.spa(ind[i6],ind[i1])+mc.spb(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i2]))*mc.spa(ind[i6],ind[i5]))/
 (pow(mc.spb(ind[i3],ind[i4]),2)*mc.spb(ind[i4],ind[i5])*(mc.spb(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i1]) +
 mc.spb(ind[i2],ind[i6])*mc.spa(ind[i3],ind[i2]))*(-(mc.spb(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i2]))-mc.spb(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i2]))*
mc.spa(ind[i6],ind[i1])*(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i1])+mc.spb(ind[i2],ind[i6])*mc.spa(ind[i6],ind[i5]))) +
(std::complex<T>(1.,0.)/std::complex<T>(6.,0.)*pow(mc.spb(ind[i1],ind[i2]),2)*(-((mc.spb(ind[i1],ind[i4])*mc.spa(ind[i4],ind[i3]))/mc.spa(ind[i3],ind[i2])) -
 (mc.spb(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i5]))/mc.spa(ind[i6],ind[i1])))/(mc.spb(ind[i1],ind[i6])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i3],ind[i4])*
mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])) +
(-std::complex<T>(1.,0.)/std::complex<T>(6.,0.)*(-(((-pow(mc.s(ind[i2],ind[i3],ind[i4]),2)+pow(mc.spb(ind[i2],ind[i3]),2)*
pow(mc.spa(ind[i3],ind[i2]),2))*mc.spb(ind[i1],ind[i4])*mc.spb(ind[i2],ind[i4])*mc.spa(ind[i4],ind[i3])*
(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i2]))-mc.spb(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i3]))*
(-(mc.spb(ind[i1],ind[i4])*mc.spb(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3]))+mc.spb(ind[i2],ind[i4])*
(-(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i2]))-mc.spb(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i3]))))/
   (pow(-mc.s(ind[i2],ind[i4])-mc.s(ind[i3],ind[i4]),3)*mc.s(ind[i2],ind[i3],ind[i4])*mc.spa(ind[i3],ind[i2]))) -
 (((pow(mc.spb(ind[i1],ind[i6]),2)*pow(mc.spa(ind[i6],ind[i1]),2))/mc.s(ind[i2],ind[i3],ind[i4])-mc.s(ind[i2],ind[i3],ind[i4]))*
   mc.spb(ind[i1],ind[i5])*mc.spb(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i5])*(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i1]) +
mc.spb(ind[i2],ind[i6])*mc.spa(ind[i6],ind[i5]))*(mc.spb(ind[i1],ind[i6])*mc.spb(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i5]) +
mc.spb(ind[i1],ind[i5])*(mc.spb(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i1])+mc.spb(ind[i2],ind[i6])*mc.spa(ind[i6],ind[i5]))))/
  (pow(mc.s(ind[i1],ind[i6])-mc.s(ind[i2],ind[i3],ind[i4]),3)*mc.spa(ind[i6],ind[i1]))))/(mc.spb(ind[i1],ind[i6])*mc.spb(ind[i2],ind[i3])*
mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6]))
#endif
);
}

//This result is not flip symetric and so we will not use it for now
/*template <int i1, int i2, int i3, int i4, int i5, int i6, class T> std::complex<T> R6g3ma(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
 return(
(-std::complex<T>(0.,1.)/std::complex<T>(6.,0.))*((-((pow(mc.spa(ind[i1],ind[i3]),3)*mc.spa(ind[i2],ind[i5])*mc.spb(ind[i3],ind[i2]))/
  (mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5]))) +
(pow(mc.spa(ind[i1],ind[i3]),2)*(std::complex<T>(-3.,0)*mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2]) -
   mc.spa(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i3])))/(mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])) -
(pow(mc.spb(ind[i6],ind[i4]),3)*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i5],ind[i2]))/(mc.spb(ind[i2],ind[i1])*mc.spb(ind[i4],ind[i3])*
  mc.spb(ind[i6],ind[i1]))-(pow(mc.spb(ind[i6],ind[i4]),ind[i2])*(std::complex<T>(3.,0)*mc.spa(ind[i1],ind[i5])*mc.spb(ind[i5],ind[i4]) +
   mc.spa(ind[i1],ind[i6])*mc.spb(ind[i6],ind[i4])))/(mc.spb(ind[i4],ind[i3])*mc.spb(ind[i6],ind[i1])) +
(pow(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2]))-mc.spa(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i3]),ind[i2])*
  (mc.spa(ind[i1],ind[i3])/mc.spa(ind[i3],ind[i4])+(-(mc.spa(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i2])) -
 mc.spa(ind[i1],ind[i5])*mc.spb(ind[i5],ind[i4]))/mc.s(ind[i2],ind[i3],ind[i4])+mc.spb(ind[i6],ind[i4])/mc.spb(ind[i6],ind[i1])))/
 (mc.spa(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i3])))/(mc.spa(ind[i5],ind[i6])*mc.spb(ind[i3],ind[i2])*
(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2]))-mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2]))) -
(((mc.s(ind[i2],ind[i3])+mc.s(ind[i2],ind[i4])+std::complex<T>(2.,0)*mc.s(ind[i3],ind[i4]))*mc.spa(ind[i1],ind[i2])*mc.spa(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i2])*
  (mc.spa(ind[i1],ind[i3])*mc.spb(ind[i3],ind[i2])+mc.spa(ind[i1],ind[i4])*mc.spb(ind[i4],ind[i2]))*
  (mc.s(ind[i2],ind[i3],ind[i4])*mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i2]) +
   mc.s(ind[i3],ind[i4])*(-(mc.s(ind[i2],ind[i3])*mc.spa(ind[i1],ind[i3]))-mc.spa(ind[i1],ind[i4])*mc.spa(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i2]))))/
 (pow(-mc.s(ind[i2],ind[i3])-mc.s(ind[i2],ind[i4]),2)*mc.s(ind[i2],ind[i3],ind[i4])*mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i6])*
  mc.spb(ind[i3],ind[i2]))+((std::complex<T>(2.,0)*mc.s(ind[i3],ind[i4])+mc.s(ind[i3],ind[i5])+mc.s(ind[i4],ind[i5]))*mc.spa(ind[i3],ind[i5])*
  mc.spb(ind[i5],ind[i4])*(mc.spa(ind[i1],ind[i5])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i2]))*mc.spb(ind[i6],ind[i5])*
  (mc.s(ind[i3],ind[i4])*(-(mc.spa(ind[i3],ind[i4])*mc.spb(ind[i6],ind[i4]))-mc.spa(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i5])) +
   mc.s(ind[i3],ind[i4],ind[i5])*(mc.spa(ind[i3],ind[i4])*mc.spb(ind[i6],ind[i4])-mc.spa(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i5]))))/
 (pow(-mc.s(ind[i3],ind[i5])-mc.s(ind[i4],ind[i5]),2)*mc.s(ind[i3],ind[i4],ind[i5])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i2],ind[i1])*
  mc.spb(ind[i6],ind[i1])))/(mc.s(ind[i3],ind[i4])*(-(mc.spa(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i2])) -
 mc.spa(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i2])))+((pow(mc.spa(ind[i1],ind[i3]),3)*mc.spa(ind[i2],ind[i5])*mc.spb(ind[i2],ind[i1]))/
 (mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i6]))+(pow(mc.spb(ind[i6],ind[i4]),3)*mc.spa(ind[i4],ind[i5])*
  mc.spb(ind[i5],ind[i2]))/(mc.spb(ind[i3],ind[i2])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i6],ind[i1])) +
(pow(mc.spa(ind[i1],ind[i3]),2)*(mc.spa(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])+std::complex<T>(3.,0)*mc.spa(ind[i2],ind[i3])*
mc.spb(ind[i6],ind[i2])))/(mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i4])) -
(pow(mc.spb(ind[i6],ind[i4]),2)*(-(mc.spa(ind[i3],ind[i4])*mc.spb(ind[i6],ind[i4]))+std::complex<T>(-3.,0)*mc.spa(ind[i3],ind[i5])*
mc.spb(ind[i6],ind[i5])))/(mc.spb(ind[i4],ind[i3])*mc.spb(ind[i6],ind[i1])) +
(pow(mc.spa(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])+mc.spa(ind[i2],ind[i3])*mc.spb(ind[i6],ind[i2]),2)*
  (-(mc.spa(ind[i1],ind[i3])/mc.spa(ind[i1],ind[i6]))-mc.spb(ind[i6],ind[i4])/mc.spb(ind[i4],ind[i3]) +
   (mc.spa(ind[i2],ind[i3])*mc.spb(ind[i6],ind[i2])+mc.spa(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i5]))/mc.s(ind[i1],ind[i2],ind[i6])))/
 (mc.spa(ind[i3],ind[i4])*mc.spb(ind[i6],ind[i1])))/(mc.spa(ind[i4],ind[i5])*mc.spb(ind[i2],ind[i1])*
(mc.spa(ind[i1],ind[i5])*mc.spb(ind[i2],ind[i1]) +
      mc.spa(ind[i5],ind[i6])*mc.spb(ind[i6],ind[i2]))))  );
}

template <int i1, int i2, int i3, int i4, int i5, int i6, class T> std::complex<T> R6g3m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
 return R6g3ma<0,1,2,3,4,5>(mc,ind)+R6g3ma<2,1,0,5,4,3>(mc,ind);
}*/

template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> std::complex<T> R7gap(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
  return std::complex<T>(0.,0.);

}

template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> std::complex<T> R7g1m(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
  return std::complex<T>(0.,1.)/T(3.)*(((mc.spa(ind[i1],ind[i2])*mc.spa(ind[i1],ind[i3])*(mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i1])*mc.spb(ind[i2],ind[i3]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i1])*mc.spb(ind[i2],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i1])*mc.spb(ind[i2],ind[i5]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i6]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i7]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i1])*mc.spb(ind[i3],ind[i4]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i1])*mc.spb(ind[i3],ind[i5]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i3],ind[i6]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i7])))/
      mc.spa(ind[i2],ind[i3]) + (mc.spa(ind[i1],ind[i3])*mc.spa(ind[i1],ind[i4])*
        (mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i1])*mc.spb(ind[i3],ind[i4]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i1])*mc.spb(ind[i3],ind[i5]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i3],ind[i6]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i7]) +
          mc.spa(ind[i1],ind[i4])*mc.spa(ind[i5],ind[i1])*mc.spb(ind[i4],ind[i5]) + mc.spa(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i4],ind[i6]) +
          mc.spa(ind[i1],ind[i4])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i4],ind[i7])))/mc.spa(ind[i3],ind[i4]) +
     (mc.spa(ind[i1],ind[i4])*mc.spa(ind[i1],ind[i5])*(mc.spa(ind[i1],ind[i4])*mc.spa(ind[i5],ind[i1])*mc.spb(ind[i4],ind[i5]) +
          mc.spa(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i4],ind[i6]) + mc.spa(ind[i1],ind[i4])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i4],ind[i7]) +
          mc.spa(ind[i1],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i5],ind[i6]) + mc.spa(ind[i1],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i5],ind[i7])))/
      mc.spa(ind[i4],ind[i5]) + (pow(mc.spa(ind[i1],ind[i6]),2)*mc.spa(ind[i1],ind[i7])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i6],ind[i7]))/mc.spa(ind[i6],ind[i7]) +
     (mc.spa(ind[i1],ind[i5])*mc.spa(ind[i1],ind[i6])*(mc.spa(ind[i1],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i5],ind[i6]) +
          mc.spa(ind[i1],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i5],ind[i7]) + mc.spa(ind[i1],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i6],ind[i7])))/
      mc.spa(ind[i5],ind[i6]) + (mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i5])*
        pow(mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i1])*mc.spb(ind[i3],ind[i5]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i3],ind[i6]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i7]) + mc.spa(ind[i1],ind[i4])*mc.spa(ind[i5],ind[i1])*mc.spb(ind[i4],ind[i5]) +
          mc.spa(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i4],ind[i6]) + mc.spa(ind[i1],ind[i4])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i4],ind[i7]),3)*
        (mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i4]),2)*mc.spa(ind[i5],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i4],ind[i5]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i4]),2)*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i4],ind[i6]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i4]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i4],ind[i7])))/
      ((mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spb(ind[i2],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i4])*mc.spb(ind[i2],ind[i4]))*
        (mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spb(ind[i2],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i2],ind[i4]))*
        (mc.spa(ind[i1],ind[i5])*mc.spa(ind[i3],ind[i2])*mc.spb(ind[i5],ind[i3]) + mc.spa(ind[i1],ind[i5])*mc.spa(ind[i4],ind[i2])*mc.spb(ind[i5],ind[i4]) +
          mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i2])*mc.spb(ind[i6],ind[i3]) + mc.spa(ind[i1],ind[i6])*mc.spa(ind[i4],ind[i2])*mc.spb(ind[i6],ind[i4]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i3],ind[i2])*mc.spb(ind[i7],ind[i3]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i2])*mc.spb(ind[i7],ind[i4]))*
        (mc.spa(ind[i1],ind[i5])*mc.spa(ind[i3],ind[i3])*mc.spb(ind[i5],ind[i3]) + mc.spa(ind[i1],ind[i5])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i5],ind[i4]) +
          mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i3])*mc.spb(ind[i6],ind[i3]) + mc.spa(ind[i1],ind[i6])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i6],ind[i4]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i3],ind[i3])*mc.spb(ind[i7],ind[i3]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i7],ind[i4]))*
        mc.s(ind[i3],ind[i4])) + (mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i6])*
        pow(mc.spa(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i4],ind[i6]) + mc.spa(ind[i1],ind[i4])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i4],ind[i7]) +
          mc.spa(ind[i1],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i5],ind[i6]) + mc.spa(ind[i1],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i5],ind[i7]),3)*
        (mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i6]) +
          mc.spa(ind[i1],ind[i3])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i6]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i3])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i7])))/
      ((mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i2],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i5])*mc.spb(ind[i2],ind[i5]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i3],ind[i4]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i5])*mc.spb(ind[i3],ind[i5]))*
        (mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i6])*mc.spb(ind[i2],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i2],ind[i5]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i6])*mc.spb(ind[i3],ind[i4]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i3],ind[i5]))*
        (mc.spa(ind[i1],ind[i6])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i6],ind[i4]) + mc.spa(ind[i1],ind[i6])*mc.spa(ind[i5],ind[i3])*mc.spb(ind[i6],ind[i5]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i7],ind[i4]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i3])*mc.spb(ind[i7],ind[i5]))*
        (mc.spa(ind[i1],ind[i6])*mc.spa(ind[i4],ind[i4])*mc.spb(ind[i6],ind[i4]) + mc.spa(ind[i1],ind[i6])*mc.spa(ind[i5],ind[i4])*mc.spb(ind[i6],ind[i5]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i4])*mc.spb(ind[i7],ind[i4]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i4])*mc.spb(ind[i7],ind[i5]))*
        mc.s(ind[i4],ind[i5])) + (mc.spa(ind[i4],ind[i5])*mc.spa(ind[i6],ind[i7])*
        pow(mc.spa(ind[i1],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i5],ind[i7]) + mc.spa(ind[i1],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i6],ind[i7]),3)*
        (mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i5],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i3])*pow(mc.spa(ind[i5],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i4])*pow(mc.spa(ind[i5],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7])))/
      ((mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i2],ind[i5]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i6],ind[i6])*mc.spb(ind[i2],ind[i6]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i3],ind[i5]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i6])*mc.spb(ind[i3],ind[i6]) +
          mc.spa(ind[i1],ind[i4])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i4],ind[i5]) + mc.spa(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i6])*mc.spb(ind[i4],ind[i6]))*
        (mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i7])*mc.spb(ind[i2],ind[i5]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i6],ind[i7])*mc.spb(ind[i2],ind[i6]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i7])*mc.spb(ind[i3],ind[i5]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i7])*mc.spb(ind[i3],ind[i6]) +
          mc.spa(ind[i1],ind[i4])*mc.spa(ind[i5],ind[i7])*mc.spb(ind[i4],ind[i5]) + mc.spa(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i7])*mc.spb(ind[i4],ind[i6]))*
        (mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i4])*mc.spb(ind[i7],ind[i5]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i6],ind[i4])*mc.spb(ind[i7],ind[i6]))*
        (mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i5])*mc.spb(ind[i7],ind[i5]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i6],ind[i5])*mc.spb(ind[i7],ind[i6]))*
        mc.s(ind[i5],ind[i6])) + (mc.spa(ind[i2],ind[i3])*mc.spa(ind[i5],ind[i6])*
        pow(mc.spa(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i3],ind[i6]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i7]) +
          mc.spa(ind[i1],ind[i4])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i4],ind[i6]) + mc.spa(ind[i1],ind[i4])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i4],ind[i7]) +
          mc.spa(ind[i1],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i5],ind[i6]) + mc.spa(ind[i1],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i5],ind[i7]),3)*
        (mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i4]),2)*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i4],ind[i6]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i4]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i4],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i6])*
           mc.spb(ind[i5],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i4])*
           mc.spb(ind[i4],ind[i6])*mc.spb(ind[i5],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i7])*mc.spb(ind[i5],ind[i3]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i4],ind[i7])*
           mc.spb(ind[i5],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i3])*
           mc.spb(ind[i4],ind[i3])*mc.spb(ind[i5],ind[i6]) + mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i5]),2)*mc.spa(ind[i6],ind[i1])*
           mc.spb(ind[i2],ind[i3])*mc.spb(ind[i5],ind[i3])*mc.spb(ind[i5],ind[i6]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i3])*
           mc.spb(ind[i5],ind[i6]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i6],ind[i1])*mc.spb(ind[i2],ind[i3])*
           mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i6]) + mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i6],ind[i1])*
           mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i6]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*
           mc.spb(ind[i5],ind[i7]) + mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i5]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*
           mc.spb(ind[i5],ind[i3])*mc.spb(ind[i5],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i3])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i5],ind[i4])*
           mc.spb(ind[i5],ind[i7]) + mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*
           mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i7])))/
      ((mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spb(ind[i2],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spb(ind[i2],ind[i4]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i5])*mc.spb(ind[i2],ind[i5]))*
        (mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i6])*mc.spb(ind[i2],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i6])*mc.spb(ind[i2],ind[i4]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i2],ind[i5]))*
        (mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i2])*mc.spb(ind[i6],ind[i3]) + mc.spa(ind[i1],ind[i6])*mc.spa(ind[i4],ind[i2])*mc.spb(ind[i6],ind[i4]) +
          mc.spa(ind[i1],ind[i6])*mc.spa(ind[i5],ind[i2])*mc.spb(ind[i6],ind[i5]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i3],ind[i2])*mc.spb(ind[i7],ind[i3]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i2])*mc.spb(ind[i7],ind[i4]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i2])*mc.spb(ind[i7],ind[i5]))*
        (mc.spa(ind[i1],ind[i6])*mc.spa(ind[i3],ind[i3])*mc.spb(ind[i6],ind[i3]) + mc.spa(ind[i1],ind[i6])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i6],ind[i4]) +
          mc.spa(ind[i1],ind[i6])*mc.spa(ind[i5],ind[i3])*mc.spb(ind[i6],ind[i5]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i3],ind[i3])*mc.spb(ind[i7],ind[i3]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i7],ind[i4]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i3])*mc.spb(ind[i7],ind[i5]))*
        mc.s(ind[i3],ind[i4],ind[i5])) + (mc.spa(ind[i3],ind[i4])*mc.spa(ind[i6],ind[i7])*
        pow(mc.spa(ind[i1],ind[i4])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i4],ind[i7]) + mc.spa(ind[i1],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i6],ind[i7]),3)*
        (mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i3])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i7])*
           mc.spb(ind[i6],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i5])*
           mc.spb(ind[i5],ind[i7])*mc.spb(ind[i6],ind[i4]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i7])*mc.spb(ind[i6],ind[i4]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i5],ind[i7])*
           mc.spb(ind[i6],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*
           mc.spb(ind[i5],ind[i4])*mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i6],ind[i4])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i4])*
           mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i3])*pow(mc.spa(ind[i4],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i4])*
           mc.spb(ind[i6],ind[i4])*mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i4])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i6],ind[i5])*
           mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i5],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i5])*
           mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i3],ind[i4])*mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i3])*pow(mc.spa(ind[i5],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7])))/
      ((mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i6])*mc.spb(ind[i2],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i2],ind[i5]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i6],ind[i6])*mc.spb(ind[i2],ind[i6]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i6])*mc.spb(ind[i3],ind[i4]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i3],ind[i5]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i6])*mc.spb(ind[i3],ind[i6]))*
        (mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i7])*mc.spb(ind[i2],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i7])*mc.spb(ind[i2],ind[i5]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i6],ind[i7])*mc.spb(ind[i2],ind[i6]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i4],ind[i7])*mc.spb(ind[i3],ind[i4]) +
          mc.spa(ind[i1],ind[i3])*mc.spa(ind[i5],ind[i7])*mc.spb(ind[i3],ind[i5]) + mc.spa(ind[i1],ind[i3])*mc.spa(ind[i6],ind[i7])*mc.spb(ind[i3],ind[i6]))*
        (mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i7],ind[i4]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i3])*mc.spb(ind[i7],ind[i5]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i6],ind[i3])*mc.spb(ind[i7],ind[i6]))*
        (mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i4])*mc.spb(ind[i7],ind[i4]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i4])*mc.spb(ind[i7],ind[i5]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i6],ind[i4])*mc.spb(ind[i7],ind[i6]))*mc.s(ind[i4],ind[i5],ind[i6])) +
     (mc.spa(ind[i2],ind[i3])*mc.spa(ind[i6],ind[i7])*pow(mc.spa(ind[i1],ind[i3])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i3],ind[i7]) +
          mc.spa(ind[i1],ind[i4])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i4],ind[i7]) + mc.spa(ind[i1],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i6],ind[i7]),3)*
        (mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i4]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i4],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i7])*
           mc.spb(ind[i5],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*
           mc.spb(ind[i4],ind[i7])*mc.spb(ind[i5],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i5]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i5],ind[i3])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i3])*
           mc.spb(ind[i5],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*
           mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i7]) + mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i5]),2)*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i5],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i7])*
           mc.spb(ind[i6],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*
           mc.spb(ind[i4],ind[i7])*mc.spb(ind[i6],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i5])*mc.spb(ind[i4],ind[i7])*mc.spb(ind[i6],ind[i3]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i5],ind[i7])*
           mc.spb(ind[i6],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*
           mc.spb(ind[i5],ind[i7])*mc.spb(ind[i6],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i5])*mc.spb(ind[i5],ind[i7])*mc.spb(ind[i6],ind[i3]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i5],ind[i7])*
           mc.spb(ind[i6],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*
           mc.spb(ind[i5],ind[i7])*mc.spb(ind[i6],ind[i4]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i5])*mc.spb(ind[i5],ind[i7])*mc.spb(ind[i6],ind[i4]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i3])*
           mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*
           mc.spb(ind[i4],ind[i5])*mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i3])*mc.spb(ind[i5],ind[i3])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i3])*
           mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*
           mc.spb(ind[i5],ind[i4])*mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i4])*mc.spb(ind[i5],ind[i4])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i3],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i6],ind[i3])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i6],ind[i3])*
           mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i5])*
           mc.spb(ind[i6],ind[i3])*mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i3])*mc.spb(ind[i6],ind[i4])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i4],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i4])*mc.spb(ind[i6],ind[i4])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i4])*
           mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i6])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i3])*
           mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i6])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i7],ind[i1])*
           mc.spb(ind[i2],ind[i4])*mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7]) +
          mc.spa(ind[i1],ind[i2])*pow(mc.spa(ind[i5],ind[i6]),2)*mc.spa(ind[i7],ind[i1])*mc.spb(ind[i2],ind[i5])*mc.spb(ind[i6],ind[i5])*mc.spb(ind[i6],ind[i7])))/
      ((mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i6])*mc.spb(ind[i2],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i6])*mc.spb(ind[i2],ind[i4]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i6])*mc.spb(ind[i2],ind[i5]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i6],ind[i6])*mc.spb(ind[i2],ind[i6]))*
        (mc.spa(ind[i1],ind[i2])*mc.spa(ind[i3],ind[i7])*mc.spb(ind[i2],ind[i3]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i4],ind[i7])*mc.spb(ind[i2],ind[i4]) +
          mc.spa(ind[i1],ind[i2])*mc.spa(ind[i5],ind[i7])*mc.spb(ind[i2],ind[i5]) + mc.spa(ind[i1],ind[i2])*mc.spa(ind[i6],ind[i7])*mc.spb(ind[i2],ind[i6]))*
        (mc.spa(ind[i1],ind[i7])*mc.spa(ind[i3],ind[i2])*mc.spb(ind[i7],ind[i3]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i2])*mc.spb(ind[i7],ind[i4]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i2])*mc.spb(ind[i7],ind[i5]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i6],ind[i2])*mc.spb(ind[i7],ind[i6]))*
        (mc.spa(ind[i1],ind[i7])*mc.spa(ind[i3],ind[i3])*mc.spb(ind[i7],ind[i3]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i4],ind[i3])*mc.spb(ind[i7],ind[i4]) +
          mc.spa(ind[i1],ind[i7])*mc.spa(ind[i5],ind[i3])*mc.spb(ind[i7],ind[i5]) + mc.spa(ind[i1],ind[i7])*mc.spa(ind[i6],ind[i3])*mc.spb(ind[i7],ind[i6]))*
        mc.s(ind[i3],ind[i4],ind[i5],ind[i6])))/
   (mc.spa(ind[i1],ind[i2])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i5],ind[i6])*mc.spa(ind[i6],ind[i7])*
     mc.spa(ind[i7],ind[i1])));
}

template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> std::complex<T> R7g1p(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
  return std::complex<T>(0.,1.)/T(3.)*(((mc.spb(ind[i1],ind[i2])*mc.spb(ind[i1],ind[i3])*(mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i1])*mc.spa(ind[i2],ind[i3]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i1])*mc.spa(ind[i2],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i1])*mc.spa(ind[i2],ind[i5]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i6]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i7]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i1])*mc.spa(ind[i3],ind[i4]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i1])*mc.spa(ind[i3],ind[i5]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i3],ind[i6]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i7])))/
      mc.spb(ind[i2],ind[i3]) + (mc.spb(ind[i1],ind[i3])*mc.spb(ind[i1],ind[i4])*
        (mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i1])*mc.spa(ind[i3],ind[i4]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i1])*mc.spa(ind[i3],ind[i5]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i3],ind[i6]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i7]) +
          mc.spb(ind[i1],ind[i4])*mc.spb(ind[i5],ind[i1])*mc.spa(ind[i4],ind[i5]) + mc.spb(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i4],ind[i6]) +
          mc.spb(ind[i1],ind[i4])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i4],ind[i7])))/mc.spb(ind[i3],ind[i4]) +
     (mc.spb(ind[i1],ind[i4])*mc.spb(ind[i1],ind[i5])*(mc.spb(ind[i1],ind[i4])*mc.spb(ind[i5],ind[i1])*mc.spa(ind[i4],ind[i5]) +
          mc.spb(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i4],ind[i6]) + mc.spb(ind[i1],ind[i4])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i4],ind[i7]) +
          mc.spb(ind[i1],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i5],ind[i6]) + mc.spb(ind[i1],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i5],ind[i7])))/
      mc.spb(ind[i4],ind[i5]) + (pow(mc.spb(ind[i1],ind[i6]),2)*mc.spb(ind[i1],ind[i7])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i6],ind[i7]))/mc.spb(ind[i6],ind[i7]) +
     (mc.spb(ind[i1],ind[i5])*mc.spb(ind[i1],ind[i6])*(mc.spb(ind[i1],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i5],ind[i6]) +
          mc.spb(ind[i1],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i5],ind[i7]) + mc.spb(ind[i1],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i6],ind[i7])))/
      mc.spb(ind[i5],ind[i6]) + (mc.spb(ind[i2],ind[i3])*mc.spb(ind[i4],ind[i5])*
        pow(mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i1])*mc.spa(ind[i3],ind[i5]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i3],ind[i6]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i7]) + mc.spb(ind[i1],ind[i4])*mc.spb(ind[i5],ind[i1])*mc.spa(ind[i4],ind[i5]) +
          mc.spb(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i4],ind[i6]) + mc.spb(ind[i1],ind[i4])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i4],ind[i7]),3)*
        (mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i4]),2)*mc.spb(ind[i5],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i4],ind[i5]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i4]),2)*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i4],ind[i6]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i4]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i4],ind[i7])))/
      ((mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spa(ind[i2],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i4])*mc.spa(ind[i2],ind[i4]))*
        (mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spa(ind[i2],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i2],ind[i4]))*
        (mc.spb(ind[i1],ind[i5])*mc.spb(ind[i3],ind[i2])*mc.spa(ind[i5],ind[i3]) + mc.spb(ind[i1],ind[i5])*mc.spb(ind[i4],ind[i2])*mc.spa(ind[i5],ind[i4]) +
          mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i2])*mc.spa(ind[i6],ind[i3]) + mc.spb(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i2])*mc.spa(ind[i6],ind[i4]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i3],ind[i2])*mc.spa(ind[i7],ind[i3]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i2])*mc.spa(ind[i7],ind[i4]))*
        (mc.spb(ind[i1],ind[i5])*mc.spb(ind[i3],ind[i3])*mc.spa(ind[i5],ind[i3]) + mc.spb(ind[i1],ind[i5])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i5],ind[i4]) +
          mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i3])*mc.spa(ind[i6],ind[i3]) + mc.spb(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i6],ind[i4]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i3],ind[i3])*mc.spa(ind[i7],ind[i3]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i7],ind[i4]))*
        mc.s(ind[i3],ind[i4])) + (mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i6])*
        pow(mc.spb(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i4],ind[i6]) + mc.spb(ind[i1],ind[i4])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i4],ind[i7]) +
          mc.spb(ind[i1],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i5],ind[i6]) + mc.spb(ind[i1],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i5],ind[i7]),3)*
        (mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i6]) +
          mc.spb(ind[i1],ind[i3])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i6]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i3])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i7])))/
      ((mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i2],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i5])*mc.spa(ind[i2],ind[i5]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i3],ind[i4]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i5])*mc.spa(ind[i3],ind[i5]))*
        (mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i6])*mc.spa(ind[i2],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i2],ind[i5]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i6])*mc.spa(ind[i3],ind[i4]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i3],ind[i5]))*
        (mc.spb(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i6],ind[i4]) + mc.spb(ind[i1],ind[i6])*mc.spb(ind[i5],ind[i3])*mc.spa(ind[i6],ind[i5]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i7],ind[i4]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i3])*mc.spa(ind[i7],ind[i5]))*
        (mc.spb(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i4])*mc.spa(ind[i6],ind[i4]) + mc.spb(ind[i1],ind[i6])*mc.spb(ind[i5],ind[i4])*mc.spa(ind[i6],ind[i5]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i4])*mc.spa(ind[i7],ind[i4]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i4])*mc.spa(ind[i7],ind[i5]))*
        mc.s(ind[i4],ind[i5])) + (mc.spb(ind[i4],ind[i5])*mc.spb(ind[i6],ind[i7])*
        pow(mc.spb(ind[i1],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i5],ind[i7]) + mc.spb(ind[i1],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i6],ind[i7]),3)*
        (mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i5],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i3])*pow(mc.spb(ind[i5],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i4])*pow(mc.spb(ind[i5],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i4],ind[i5])*mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7])))/
      ((mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i2],ind[i5]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i6],ind[i6])*mc.spa(ind[i2],ind[i6]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i3],ind[i5]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i6])*mc.spa(ind[i3],ind[i6]) +
          mc.spb(ind[i1],ind[i4])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i4],ind[i5]) + mc.spb(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i6])*mc.spa(ind[i4],ind[i6]))*
        (mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i7])*mc.spa(ind[i2],ind[i5]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i6],ind[i7])*mc.spa(ind[i2],ind[i6]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i7])*mc.spa(ind[i3],ind[i5]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i7])*mc.spa(ind[i3],ind[i6]) +
          mc.spb(ind[i1],ind[i4])*mc.spb(ind[i5],ind[i7])*mc.spa(ind[i4],ind[i5]) + mc.spb(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i7])*mc.spa(ind[i4],ind[i6]))*
        (mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i4])*mc.spa(ind[i7],ind[i5]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i6],ind[i4])*mc.spa(ind[i7],ind[i6]))*
        (mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i5])*mc.spa(ind[i7],ind[i5]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i6],ind[i5])*mc.spa(ind[i7],ind[i6]))*
        mc.s(ind[i5],ind[i6])) + (mc.spb(ind[i2],ind[i3])*mc.spb(ind[i5],ind[i6])*
        pow(mc.spb(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i3],ind[i6]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i7]) +
          mc.spb(ind[i1],ind[i4])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i4],ind[i6]) + mc.spb(ind[i1],ind[i4])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i4],ind[i7]) +
          mc.spb(ind[i1],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i5],ind[i6]) + mc.spb(ind[i1],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i5],ind[i7]),3)*
        (mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i4]),2)*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i4],ind[i6]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i4]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i4],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i6])*
           mc.spa(ind[i5],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i4])*
           mc.spa(ind[i4],ind[i6])*mc.spa(ind[i5],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i7])*mc.spa(ind[i5],ind[i3]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i4],ind[i7])*
           mc.spa(ind[i5],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i3])*
           mc.spa(ind[i4],ind[i3])*mc.spa(ind[i5],ind[i6]) + mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i5]),2)*mc.spb(ind[i6],ind[i1])*
           mc.spa(ind[i2],ind[i3])*mc.spa(ind[i5],ind[i3])*mc.spa(ind[i5],ind[i6]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i3])*
           mc.spa(ind[i5],ind[i6]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i6],ind[i1])*mc.spa(ind[i2],ind[i3])*
           mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i6]) + mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i6],ind[i1])*
           mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i6]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*
           mc.spa(ind[i5],ind[i7]) + mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i5]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*
           mc.spa(ind[i5],ind[i3])*mc.spa(ind[i5],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i3])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i5],ind[i4])*
           mc.spa(ind[i5],ind[i7]) + mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*
           mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i7])))/
      ((mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spa(ind[i2],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spa(ind[i2],ind[i4]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i5])*mc.spa(ind[i2],ind[i5]))*
        (mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i6])*mc.spa(ind[i2],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i6])*mc.spa(ind[i2],ind[i4]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i2],ind[i5]))*
        (mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i2])*mc.spa(ind[i6],ind[i3]) + mc.spb(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i2])*mc.spa(ind[i6],ind[i4]) +
          mc.spb(ind[i1],ind[i6])*mc.spb(ind[i5],ind[i2])*mc.spa(ind[i6],ind[i5]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i3],ind[i2])*mc.spa(ind[i7],ind[i3]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i2])*mc.spa(ind[i7],ind[i4]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i2])*mc.spa(ind[i7],ind[i5]))*
        (mc.spb(ind[i1],ind[i6])*mc.spb(ind[i3],ind[i3])*mc.spa(ind[i6],ind[i3]) + mc.spb(ind[i1],ind[i6])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i6],ind[i4]) +
          mc.spb(ind[i1],ind[i6])*mc.spb(ind[i5],ind[i3])*mc.spa(ind[i6],ind[i5]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i3],ind[i3])*mc.spa(ind[i7],ind[i3]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i7],ind[i4]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i3])*mc.spa(ind[i7],ind[i5]))*
        mc.s(ind[i3],ind[i4],ind[i5])) + (mc.spb(ind[i3],ind[i4])*mc.spb(ind[i6],ind[i7])*
        pow(mc.spb(ind[i1],ind[i4])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i4],ind[i7]) + mc.spb(ind[i1],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i6],ind[i7]),3)*
        (mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i3])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i7])*
           mc.spa(ind[i6],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i5])*
           mc.spa(ind[i5],ind[i7])*mc.spa(ind[i6],ind[i4]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i7])*mc.spa(ind[i6],ind[i4]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i5],ind[i7])*
           mc.spa(ind[i6],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*
           mc.spa(ind[i5],ind[i4])*mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i3],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i6],ind[i4])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i4])*
           mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i3])*pow(mc.spb(ind[i4],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i4])*
           mc.spa(ind[i6],ind[i4])*mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i3],ind[i5])*mc.spa(ind[i6],ind[i4])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i6],ind[i5])*
           mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i5],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i5])*
           mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i3],ind[i4])*mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i3])*pow(mc.spb(ind[i5],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i5])*mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7])))/
      ((mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i6])*mc.spa(ind[i2],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i2],ind[i5]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i6],ind[i6])*mc.spa(ind[i2],ind[i6]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i6])*mc.spa(ind[i3],ind[i4]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i3],ind[i5]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i6])*mc.spa(ind[i3],ind[i6]))*
        (mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i7])*mc.spa(ind[i2],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i7])*mc.spa(ind[i2],ind[i5]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i6],ind[i7])*mc.spa(ind[i2],ind[i6]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i4],ind[i7])*mc.spa(ind[i3],ind[i4]) +
          mc.spb(ind[i1],ind[i3])*mc.spb(ind[i5],ind[i7])*mc.spa(ind[i3],ind[i5]) + mc.spb(ind[i1],ind[i3])*mc.spb(ind[i6],ind[i7])*mc.spa(ind[i3],ind[i6]))*
        (mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i7],ind[i4]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i3])*mc.spa(ind[i7],ind[i5]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i6],ind[i3])*mc.spa(ind[i7],ind[i6]))*
        (mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i4])*mc.spa(ind[i7],ind[i4]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i4])*mc.spa(ind[i7],ind[i5]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i6],ind[i4])*mc.spa(ind[i7],ind[i6]))*mc.s(ind[i4],ind[i5],ind[i6])) +
     (mc.spb(ind[i2],ind[i3])*mc.spb(ind[i6],ind[i7])*pow(mc.spb(ind[i1],ind[i3])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i3],ind[i7]) +
          mc.spb(ind[i1],ind[i4])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i4],ind[i7]) + mc.spb(ind[i1],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i6],ind[i7]),3)*
        (mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i4]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i4],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i7])*
           mc.spa(ind[i5],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*
           mc.spa(ind[i4],ind[i7])*mc.spa(ind[i5],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i5]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i5],ind[i3])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i3])*
           mc.spa(ind[i5],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*
           mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i7]) + mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i5]),2)*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i5],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i7])*
           mc.spa(ind[i6],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*
           mc.spa(ind[i4],ind[i7])*mc.spa(ind[i6],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i5])*mc.spa(ind[i4],ind[i7])*mc.spa(ind[i6],ind[i3]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i5],ind[i7])*
           mc.spa(ind[i6],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*
           mc.spa(ind[i5],ind[i7])*mc.spa(ind[i6],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i5])*mc.spa(ind[i5],ind[i7])*mc.spa(ind[i6],ind[i3]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i5],ind[i7])*
           mc.spa(ind[i6],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*
           mc.spa(ind[i5],ind[i7])*mc.spa(ind[i6],ind[i4]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i5])*mc.spa(ind[i5],ind[i7])*mc.spa(ind[i6],ind[i4]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i4],ind[i3])*
           mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*
           mc.spa(ind[i4],ind[i5])*mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i3])*mc.spa(ind[i5],ind[i3])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i3])*
           mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i5])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*
           mc.spa(ind[i5],ind[i4])*mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i4])*mc.spa(ind[i5],ind[i4])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i3],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*mc.spa(ind[i6],ind[i3])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i6],ind[i3])*
           mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i5])*
           mc.spa(ind[i6],ind[i3])*mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i3])*mc.spa(ind[i6],ind[i4])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i4],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i4])*mc.spa(ind[i6],ind[i4])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i4])*
           mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i6])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i3])*
           mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i6])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i7],ind[i1])*
           mc.spa(ind[i2],ind[i4])*mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7]) +
          mc.spb(ind[i1],ind[i2])*pow(mc.spb(ind[i5],ind[i6]),2)*mc.spb(ind[i7],ind[i1])*mc.spa(ind[i2],ind[i5])*mc.spa(ind[i6],ind[i5])*mc.spa(ind[i6],ind[i7])))/
      ((mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i6])*mc.spa(ind[i2],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i6])*mc.spa(ind[i2],ind[i4]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i6])*mc.spa(ind[i2],ind[i5]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i6],ind[i6])*mc.spa(ind[i2],ind[i6]))*
        (mc.spb(ind[i1],ind[i2])*mc.spb(ind[i3],ind[i7])*mc.spa(ind[i2],ind[i3]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i4],ind[i7])*mc.spa(ind[i2],ind[i4]) +
          mc.spb(ind[i1],ind[i2])*mc.spb(ind[i5],ind[i7])*mc.spa(ind[i2],ind[i5]) + mc.spb(ind[i1],ind[i2])*mc.spb(ind[i6],ind[i7])*mc.spa(ind[i2],ind[i6]))*
        (mc.spb(ind[i1],ind[i7])*mc.spb(ind[i3],ind[i2])*mc.spa(ind[i7],ind[i3]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i2])*mc.spa(ind[i7],ind[i4]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i2])*mc.spa(ind[i7],ind[i5]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i6],ind[i2])*mc.spa(ind[i7],ind[i6]))*
        (mc.spb(ind[i1],ind[i7])*mc.spb(ind[i3],ind[i3])*mc.spa(ind[i7],ind[i3]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i4],ind[i3])*mc.spa(ind[i7],ind[i4]) +
          mc.spb(ind[i1],ind[i7])*mc.spb(ind[i5],ind[i3])*mc.spa(ind[i7],ind[i5]) + mc.spb(ind[i1],ind[i7])*mc.spb(ind[i6],ind[i3])*mc.spa(ind[i7],ind[i6]))*
        mc.s(ind[i3],ind[i4],ind[i5],ind[i6])))/
   (mc.spb(ind[i1],ind[i2])*mc.spb(ind[i2],ind[i3])*mc.spb(ind[i3],ind[i4])*mc.spb(ind[i4],ind[i5])*mc.spb(ind[i5],ind[i6])*mc.spb(ind[i6],ind[i7])*
     mc.spb(ind[i7],ind[i1])));
}

}

#endif /*RAT_AMP_H_*/
