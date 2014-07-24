
#ifndef _H_RGMP_FP
#define _H_RGMP_FP


#include <iostream>
#include <iomanip>

#include <stdio.h>
#include <complex>

#include <gmpxx.h>
#include "mpreal.h"
#include "BH_typedefs.h"

int get_default_precision();

#if BH_USE_GMP
#define GMP_INIT static int prec = (RGMP::set_precision(get_default_precision()),RGMP::get_current_precision());
#endif

inline long minprec(const mpfr::mpreal& x,const mpfr::mpreal& y){
	long px=x.get_prec();
	long py=y.get_prec();
	return px<py?px:py;
}

class RGMP  {
//	mpf_class d_x;
	mpfr::mpreal d_x;
	static mp_prec_t s_default_precision;
	inline friend RGMP operator*(const RGMP& g1,const RGMP& g2 );
	inline friend RGMP operator*(int g1,const RGMP& g2 );
	inline friend RGMP operator*(const RGMP& g2, int g1);
	inline friend RGMP operator/(const RGMP& g1,const RGMP& g2 );
	inline friend RGMP operator/(const RGMP& g1,int g2 );
	inline friend RGMP operator-(const RGMP& g1,const RGMP& g2 );
	inline friend RGMP operator+(const RGMP& g1,const RGMP& g2 );

	inline friend bool operator<=(const RGMP& g1,const RGMP& g2 );
	inline friend bool operator>=(const RGMP& g1,const RGMP& g2 );
	inline friend bool operator<(const RGMP& g1,const RGMP& g2 );
	inline friend bool operator>(const RGMP& g1,const RGMP& g2 );
	inline friend bool operator==(const RGMP& g1,const RGMP& g2 );

	inline friend std::ostream& operator<<(std::ostream& os,const RGMP& g );
	inline friend std::istream& operator>>(std::istream& os, RGMP& g );

	friend bool  isinf(const RGMP& g1);
	friend bool  isnan(const RGMP& g1);
	friend RGMP exp(const RGMP& g1);
	friend RGMP exp10(const RGMP& g1);
	friend RGMP sin(const RGMP& g1);
	friend RGMP cos(const RGMP& g1);
	friend RGMP atan2(const RGMP& g1,const RGMP& g2);
	friend RGMP sqrt(const RGMP& g1);
	friend RGMP log(const RGMP& g1);
	friend RGMP log10(const RGMP& g1);
	friend RGMP abs(const RGMP& g1);
	friend RGMP fabs(const RGMP& g1);
public:
	RGMP(): d_x(0, s_default_precision) {}
	RGMP(int i): d_x(i, s_default_precision) {}
	RGMP(long int i): d_x(i, s_default_precision) {}
	RGMP(double x): d_x(x, s_default_precision) {}
	RGMP(long double x): d_x(x, s_default_precision) {}
	RGMP(const RGMP& g) {d_x.set_prec(g.d_x.get_prec());d_x=g.d_x; }
	RGMP(const RGMP& g,mp_prec_t prec): d_x(g.d_x) {d_x.set_prec(prec);}
	RGMP(const mpfr::mpreal& g): d_x(g) {}
	RGMP(const char* str): d_x(str) {}

	RGMP& operator=(const RGMP& g1){d_x.set_prec(g1.get_precision()); d_x=g1.d_x; return *this; };

	static void set_precision(int prec){s_default_precision=prec;mpfr::mpreal::set_default_prec(prec);};
	static int get_current_precision(){ return s_default_precision;};
	int get_precision()const { return d_x.get_prec();};
	int get_nbr_digits(){ return 39*((d_x.get_prec()+63)/64)/2;};
	RGMP& operator*=(const RGMP& g1){ d_x*=g1.d_x; d_x.set_prec(minprec(d_x,g1.d_x));  return *this; };
	RGMP& operator/=(const RGMP& g1){ d_x/=g1.d_x; d_x.set_prec(minprec(d_x,g1.d_x));return *this; };
	RGMP& operator+=(const RGMP& g1){ d_x+=g1.d_x; d_x.set_prec(minprec(d_x,g1.d_x));return *this; };
	RGMP& operator-=(const RGMP& g1){ d_x-=g1.d_x; d_x.set_prec(minprec(d_x,g1.d_x));return *this; };
	RGMP operator-() const { return -d_x;} ;
	double to_double() const { return double(d_x);}
	std::string to_string() const { return d_x.to_string();}
	static int get_current_nbr_digits();
	void set_prec(mp_prec_t prec){ d_x.set_prec(prec,mpfr::mpreal::default_rnd); }
};


inline RGMP operator*(const RGMP& g1,const RGMP& g2 ){return RGMP(g1.d_x*g2.d_x,minprec(g1.d_x,g2.d_x));};
inline RGMP operator*(int g1,const RGMP& g2 ){return  g2*RGMP(g1);};
inline RGMP operator*(const RGMP& g2, int g1){return  g2*RGMP(g1);};
inline RGMP operator/(const RGMP& g1,const RGMP& g2 ){return  RGMP(g1.d_x/g2.d_x,minprec(g1.d_x,g2.d_x));};
inline RGMP operator/(const RGMP& g1,int g2 ){return  g1.d_x/RGMP(g2);};
inline RGMP operator-(const RGMP& g1,const RGMP& g2 ){return RGMP(g1.d_x-g2.d_x,minprec(g1.d_x,g2.d_x));};
inline RGMP operator+(const RGMP& g1,const RGMP& g2 ){return RGMP(g1.d_x+g2.d_x,minprec(g1.d_x,g2.d_x));};

inline bool operator<(const RGMP& g1,const RGMP& g2 ){ return g1.d_x < g2.d_x; };
inline bool operator>(const RGMP& g1,const RGMP& g2 ){ return g1.d_x > g2.d_x; };
inline bool operator<=(const RGMP& g1,const RGMP& g2 ){ return g1.d_x <= g2.d_x; };
inline bool operator>=(const RGMP& g1,const RGMP& g2 ){ return g1.d_x >= g2.d_x; };
inline bool operator==(const RGMP& g1,const RGMP& g2 ){ return g1.d_x == g2.d_x; };

inline bool isinf(const RGMP& g1){return mpfr::_isinf(g1.d_x);};
inline bool isnan(const RGMP& g1){return mpfr::_isnan(g1.d_x);};
inline RGMP exp(const RGMP& g1){return RGMP(exp(g1.d_x));};
inline RGMP exp10(const RGMP& g1){return RGMP(exp10(g1.d_x));};
inline RGMP sin(const RGMP& g1){return RGMP(sin(g1.d_x));};
inline RGMP cos(const RGMP& g1){return RGMP(cos(g1.d_x));};
inline RGMP atan2(const RGMP& g1,const RGMP& g2){return RGMP(atan2(g1.d_x,g2.d_x));};
inline RGMP log(const RGMP& g1){return RGMP(log(g1.d_x));};
inline RGMP log10(const RGMP& g1){return RGMP(log10(g1.d_x));};
inline RGMP sqrt(const RGMP& g1){return RGMP(sqrt(g1.d_x));};
inline RGMP abs(const RGMP& g1){return RGMP(abs(g1.d_x));};
inline RGMP fabs(const RGMP& g1){return RGMP(fabs(g1.d_x));};

inline std::ostream& operator<<(std::ostream& os,const RGMP& g ){
  return os << g.d_x;
}

inline std::istream& operator>>(std::istream& is,RGMP& g ){
  return is >> g.d_x;
}

std::complex<RGMP> pow(std::complex<RGMP>  __x,int __n) ;



template <int N> class RGMP_FP;

template <int N> RGMP_FP<N> operator*(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
// bad because it creates a lot of ambigous resolution
//template <int N,class T> RGMP_FP<N> operator*(const T& g1,const RGMP_FP<N>& g2 );
//template <int N,class T> RGMP_FP<N> operator*(const RGMP_FP<N>& g2, const T& g1);
template <int N> RGMP_FP<N> operator*(int g1,const RGMP_FP<N>& g2 ){return  g2*RGMP_FP<N>(g1);};
template <int N> RGMP_FP<N> operator*(const RGMP_FP<N>& g2, int g1){return  g2*RGMP_FP<N>(g1);};
template <int N> RGMP_FP<N> operator/(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
template <int N> RGMP_FP<N> operator/(const RGMP_FP<N>& g1,int g2 ){return  g1/RGMP_FP<N>(g2);};
template <int N> RGMP_FP<N> operator-(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
template <int N> RGMP_FP<N> operator+(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );

template <int N> bool operator<(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
template <int N> bool operator>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
template <int N> bool operator<=(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
template <int N> bool operator>=(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
template <int N> bool operator==(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
template <int N> RGMP_FP<N> sqrt(const RGMP_FP<N>& g1){return RGMP_FP<N>(sqrt(g1.d_x));};
template <int N> RGMP_FP<N> abs(const RGMP_FP<N>& g1){return RGMP_FP<N>(abs(g1.d_x));};
template <int N> RGMP_FP<N> fabs(const RGMP_FP<N>& g1){return RGMP_FP<N>(fabs(g1.d_x));};

template <int N> std::ostream& operator<<(std::ostream& os,const RGMP_FP<N>& g );
template <int N> std::istream& operator>>(std::istream& os,RGMP_FP<N>& g );

template <int N> class RGMP_FP {
 mpf_class d_x;
 friend std::ostream& operator<< <> (std::ostream& os,const RGMP_FP<N>& g );
 friend std::istream& operator>> <> (std::istream& is,RGMP_FP<N>& g );
  friend  RGMP_FP<N> operator*<>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend  RGMP_FP<N> operator/<>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend  RGMP_FP<N> operator+<>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend  RGMP_FP<N> operator-<>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend  bool operator< <>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend  bool operator><>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend  bool operator<= <>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend  bool operator>=<>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend  bool operator==<>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 );
  friend RGMP_FP<N> sqrt<>(const RGMP_FP<N>& g1);
  friend RGMP_FP<N> abs<>(const RGMP_FP<N>& g1);
  friend RGMP_FP<N> fabs<>(const RGMP_FP<N>& g1);

public:
  RGMP_FP(): d_x(0,N){};
  template <class Type> RGMP_FP(Type x): d_x(x,N){};
  RGMP_FP<N>& operator*=(const RGMP_FP<N>& g1);
  RGMP_FP<N>& operator/=(const RGMP_FP<N>& g1);
  RGMP_FP<N>& operator+=(const RGMP_FP<N>& g1);
  RGMP_FP<N>& operator-=(const RGMP_FP<N>& g1);
  RGMP_FP<N> operator-() const;
};


template <int N> RGMP_FP<N> operator*(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return  RGMP_FP<N>(g1.d_x*g2.d_x);
} 

template <int N> RGMP_FP<N> operator/(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return  RGMP_FP<N>(g1.d_x/g2.d_x);
} 
template <int N> RGMP_FP<N> operator+(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return  RGMP_FP<N>(g1.d_x+g2.d_x);
} 
template <int N> RGMP_FP<N> operator-(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return  RGMP_FP<N>(g1.d_x-g2.d_x);
} 


template <int N> RGMP_FP<N>& RGMP_FP<N>::operator*=(const RGMP_FP<N>& g1){
  d_x*=g1.d_x;
  return *this;
} 
template <int N> RGMP_FP<N>& RGMP_FP<N>::operator/=(const RGMP_FP<N>& g1){
  d_x/=g1.d_x;
  return *this;
} 
template <int N> RGMP_FP<N>& RGMP_FP<N>::operator-=(const RGMP_FP<N>& g1){
  d_x-=g1.d_x;
  return *this;
} 
template <int N> RGMP_FP<N>& RGMP_FP<N>::operator+=(const RGMP_FP<N>& g1){
  d_x+=g1.d_x;
  return *this;
} 

template <int N> RGMP_FP<N> RGMP_FP<N>::operator-() const {
  return RGMP_FP<N>(-d_x);
}


template <int N> std::ostream& operator<<(std::ostream& os,const RGMP_FP<N>& g ){
  return os << g.d_x;
}

template <int N> std::istream& operator>>(std::istream& is,RGMP_FP<N>& g ){
  return is >> g.d_x;
}

template <int N> bool operator<(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return g1.d_x < g2.d_x;
}

template <int N> bool operator>(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return g1.d_x > g2.d_x;
}

template <int N> bool operator<=(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return g1.d_x <= g2.d_x;
}

template <int N> bool operator>=(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return g1.d_x >= g2.d_x;
}

template <int N> bool operator==(const RGMP_FP<N>& g1,const RGMP_FP<N>& g2 ){
  return g1.d_x == g2.d_x;
}

namespace BH {

template <class T> inline T DeltaZero();
template <class T> inline T DeltaZeroSqr();

template <> inline RGMP_FP<200> DeltaZero(){return RGMP_FP<200>(1e-61);};
template <> inline RGMP_FP<200> DeltaZeroSqr(){return RGMP_FP<200>(1e-31);}

template <> inline RGMP DeltaZero(){ return RGMP(exp10(RGMP(-20*((RGMP::get_current_precision()+63)/64))));};
template <> inline RGMP DeltaZeroSqr(){ return RGMP(exp10(RGMP(-10*((RGMP::get_current_precision()+63)/64))));}

RGMP get_pi();

inline double to_double(const RGMP& g){ return g.to_double(); };
inline RHP to_HP(const RGMP& g){ return RHP(g.to_string().c_str()); };
inline RVHP to_VHP(const RGMP& g){ return RVHP(g.to_string().c_str()); };
typedef std::complex<RGMP> CGMP;

inline std::complex<double> to_double(const std::complex<RGMP>& c){return std::complex<double>(to_double(c.real()),to_double(c.imag()));}
inline std::complex<dd_real> to_HP(const std::complex<RGMP>& c){return std::complex<dd_real>(to_HP(c.real()),to_HP(c.imag()));}
inline std::complex<qd_real> to_VHP(const std::complex<RGMP>& c){return std::complex<qd_real>(to_VHP(c.real()),to_VHP(c.imag()));}

}


#endif /* _H_RGMP_FP */
