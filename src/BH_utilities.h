#ifndef BH_UTILITIES_H_
#define BH_UTILITIES_H_

#include <iterator>
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include "BH_typedefs.h"
#if BH_USE_GMP
#include "gmp_r.h"
#endif

#define _DO_DEBUG 0

#if _DO_DEBUG==1
	#define _DEBUG(X) cout << "DEBUG: " << (X) << endl
#else
	#define _DEBUG(X) // cout << "DEBUG: " << (X) << endl
#endif
#if _DO_DEBUG==1
	#define _DO_IF_DEBUG(X) (X)
#else
	#define _DO_IF_DEBUG(X) // (X)
#endif



//! Macro that prints its argument and its value (useful to remember what we intended to print)
#define _PRINT(X) std::cout  <<(#X) << ": "  <<  (X) <<  std::endl
//! Macro that prints a message. The difference with _PRINT() is that it will only print the message, _PRit() will print it twice.
#define _MESSAGE(X) std::cout <<  (X) << std::endl
#define _MESSAGE2(X,Y) std::cout << (X)  << (Y) << std::endl
#define _MESSAGE3(X,Y,Z) std::cout << (X) << (Y) <<  (Z) << std::endl
#define _MESSAGE4(W,X,Y,Z) std::cout << (W) << (X) << (Y) <<  (Z) << std::endl
#define _MESSAGE5(V,W,X,Y,Z) std::cout << (V) << (W) << (X) << (Y) <<  (Z) << std::endl
//! Macro that a warning (with cerr instead of cout.)
#define _WARNING(X) std::cerr << (X) << std::endl
#define _WARNING2(X,Y) std::cerr << (X)  << (Y) << std::endl
#define _WARNING3(X,Y,Z) std::cerr << (X) << (Y) <<  (Z) << std::endl
#define _WARNING4(X,Y,Z,A) std::cerr << (X) << (Y) <<  (Z) << (A) << std::endl
#define _WARNING5(X,Y,Z,A,B) std::cerr << (X) << (Y) <<  (Z) << (A) << (B) << std::endl
#define _MESSAGE6(V,W,X,Y,Z,Z1) std::cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << std::endl
#define _MESSAGE7(V,W,X,Y,Z,Z1,Z2) std::cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << std::endl
#define _MESSAGE8(V,W,X,Y,Z,Z1,Z2,Z3) std::cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << std::endl
#define _MESSAGE9(V,W,X,Y,Z,Z1,Z2,Z3,Z4) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << std::endl
#define _MESSAGE10(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << std::endl
#define _MESSAGE11(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << std::endl
#define _MESSAGE12(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << std::endl
#define _MESSAGE13(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << std::endl
#define _MESSAGE14(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << (Z9) << std::endl
#define _MESSAGE15(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,U1) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << (Z9) << (Z10) << std::endl
#define _MESSAGE16(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,U1,U2) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << (Z9) << (U1) << (U2) << std::endl
#define _MESSAGE17(V,W,X,Y,Z,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,U1,U2,U3) cout << (V) << (W) << (X) << (Y) <<  (Z) << (Z1) << (Z2) << (Z3) << (Z4) << (Z5) << (Z6) << (Z7) << (Z8) << (Z9) << (U1) << (U2) << (U3) << std::endl


namespace BH {



// A function for printing a vector
template <class T> void _vector_cout(const std::vector<T>& out){
	std::cout << "{";
	copy(out.begin(), out.end()-1, std::ostream_iterator<T>(std::cout, ","));
	copy(out.end()-1, out.end(), std::ostream_iterator<T>(std::cout));
	std::cout << "}\n";
}

// A function that returns zero if we have nan, as this will be a case where we have 0/0 and so the amplitude should be zero
template <class T> inline std::complex<T> Amp_safe(const std::complex<T>& res)
{
	if (isnan(res.real()) || isinf(res.real())) {
		return std::complex<T>(0.,0.);
	}
	else return res;
}

template <> inline std::complex<R> Amp_safe<R>(const std::complex<R>& res)
{
	if (std::isnan(res.real()) || std::isinf(res.real())) {
		return std::complex<R>(0.,0.);
	}
	else return res;
}


// utility to print vectors

template <class T> std::ostream& operator<<(std::ostream& s, std::vector<T>& v){
std::cout << "{" ;
for (int i=0;i<v.size()-1;i++){
	std::cout << v[i] <<"," ;
}
return std::cout << v[v.size()-1] <<"}" ;
}



/*
 *
 *
 * Functions for setting the distance to zero in a type safe way
 *
 *
 */

template <class T> inline T DeltaZero();

template<> inline R DeltaZero(){return R(1e-13);}
template<> inline RHP DeltaZero(){return RHP(1e-29);}
template<> inline RVHP DeltaZero(){return RVHP(1e-61);}


template <class T> inline T DeltaZeroSqr();

template<> inline R DeltaZeroSqr(){return R(1e-7);}
template<> inline RHP DeltaZeroSqr(){return RHP(1e-15);}
template<> inline RVHP DeltaZeroSqr(){return RVHP(1e-31);}

/*
 *
 *
 * A function to return the maximum number of digits a particular type has to work with
 *
 *
 */

template <class T> inline T MaxDigits();

template<> inline R MaxDigits(){return R(16);}
template<> inline RHP MaxDigits(){return RHP(32);}
template<> inline RVHP MaxDigits(){return RVHP(64);}

#if BH_USE_GMP
template<> inline RGMP MaxDigits(){return RGMP(exp10(RGMP::get_current_nbr_digits()));}
#endif

}
#endif /*BH_UTILITIES_H_*/
