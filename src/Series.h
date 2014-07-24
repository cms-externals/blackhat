/* Series.h */

/*  David A. Kosower, November 15, 2007  */

/* Definition of an object that holds several orders of a Laurent series in
   a parameter, taken around 0.  The series is of the form
      a_{-m}/x^m + a_{-m+1}/x^(m-1) +...+a_n x^n + O(x^(n+1))
*/

#ifndef SeriesDefined
#define SeriesDefined 1

#if 0
namespace std {
static C quiet_NaN() throw()
      { return C(__builtin_nan (""),__builtin_nan ("")); }
  //template<class T> numeric_limits<T>::quiet_NaN() throw();
};
#endif

#include <vector>
#include <complex>


#ifdef BH_USE_GMP
#include "gmp_r.h"
#endif

namespace BH {
template <class T> class Series;

template <class T> Series<T> operator+(const Series<T>& s1,const Series<T>& s2);
template <class T> Series<T> operator+(const Series<T>& s1,const T& v2);
template <class T> Series<T> operator+(const T& v1,const Series<T>& s2);

template <class T> Series<T> operator-(const Series<T>& s1,const Series<T>& s2);
template <class T> Series<T> operator-(const Series<T>& s1,const T& v2);
template <class T> Series<T> operator-(const T& v1,const Series<T>& s2);

template <class T> Series<T> operator-(const Series<T>& s1);

template <class T> Series<T> operator*(const Series<T>& s1,const Series<T>& s2);
template <class T> Series<T> operator*(const Series<T>& s1,const T& v2);
template <class T> Series<T> operator*(const T& v1,const Series<T>& s2);

template <class T> Series<T> operator/(const Series<T>& s,const T& v);

template <class T> Series<T> operator^(const Series<T>& s,unsigned int p);

template <class T> std::ostream& operator<<(std::ostream&, const Series<T>& );



template <class T> class Series {
protected:
  short int startOrder, endOrder;
  std::vector<T> term;
  std::string parameter; // optional, just for printing

#if defined(DefiningSeries)|defined(ImplementingSeries)
  T& operator[](int order) {
#if 0 // No range checking!
    if (order < startOrder) return(0); // term is truly absent
    if (order > endOrder) return(std::numeric_limits<T>::quiet_NaN); // undefined
#endif
    return term[order-startOrder];}
#endif

public:
#ifdef SWIG
%ignore operator=;
#endif
  inline Series<T> operator=(const Series<T>& s)
  {this->term = s.term;
  this->startOrder = s.startOrder;
  this->endOrder = s.endOrder;
  this->parameter = s.parameter;

  return *this;
  }
#ifndef SWIG
  friend Series<T> operator+<>(const Series<T>&,const Series<T>&);
  friend Series<T> operator+<>(const Series<T>&,const T&);
  friend Series<T> operator+<>(const T& v1,const Series<T>&);

  friend Series<T> operator-<>(const Series<T>&,const Series<T>&);
  friend Series<T> operator-<>(const Series<T>&,const T&);
  friend Series<T> operator-<>(const T&,const Series<T>&);

  Series<T>& operator+=(const Series<T>&);
  Series<T>& operator+=(const T&);

  Series<T> operator-=(const Series<T>&);
  Series<T> operator-=(const T&);

  Series<T> operator*=(const Series<T>&);
  Series<T> operator*=(const T&);

  friend Series<T> operator-<>(const Series<T>&);

  friend Series<T> operator*<>(const Series<T>&,const Series<T>&);
  friend Series<T> operator*<>(const Series<T>&,const T&);
  friend Series<T> operator*<>(const T&,const Series<T>&);

  friend Series<T> operator/<>(const Series<T>&,const T&);
  Series<T> operator/=(const T&);

  friend Series<T> operator^<>(const Series<T>&,unsigned int);
  Series<T> operator^=(unsigned int p);

  friend std::ostream& operator<<<>(std::ostream&, const Series<T>&);
#endif
  Series<T>() {}

  // Define a zero series object with given range of orders
  Series<T>(int s, int e) : startOrder(s), endOrder(e), term(e-s+1) {}

  // Define a series object starting with the given order
  // We would like to use the following,
#if 0 // but can't, because T may be non-POD and g++ can't handle that
  Series<T>(int s, int e, ...) : startOrder(s), endOrder(e) {
    va_list coeff; short int order = s;

    va_start(coeff, e);
    while (order++ <= e) term.push_back(va_arg(coeff,T));
    va_end(coeff);
  }
#endif
  // So instead we have a whole long list of functions
  Series<T>(int s, int e, T v1) : startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
  }
  Series<T>(int s, int e, T v1, T v2) :
    startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
    if (order++ <= e) term.push_back(v2);
  }
  Series<T>(int s, int e, T v1, T v2, T v3) :
    startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
    if (order++ <= e) term.push_back(v2);
    if (order++ <= e) term.push_back(v3);
  }
  Series<T>(int s, int e, T v1, T v2, T v3,
            T v4) :
    startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
    if (order++ <= e) term.push_back(v2);
    if (order++ <= e) term.push_back(v3);
    if (order++ <= e) term.push_back(v4);
  }
  Series<T>(int s, int e, T v1, T v2, T v3,
            T v4, T v5) :
    startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
    if (order++ <= e) term.push_back(v2);
    if (order++ <= e) term.push_back(v3);
    if (order++ <= e) term.push_back(v4);
    if (order++ <= e) term.push_back(v5);
  }
  Series<T>(int s, int e, T v1, T v2, T v3,
            T v4, T v5, T v6) :
    startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
    if (order++ <= e) term.push_back(v2);
    if (order++ <= e) term.push_back(v3);
    if (order++ <= e) term.push_back(v4);
    if (order++ <= e) term.push_back(v5);
    if (order++ <= e) term.push_back(v6);
  }
  Series<T>(int s, int e, T v1, T v2, T v3,
            T v4, T v5, T v6, T v7) :
    startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
    if (order++ <= e) term.push_back(v2);
    if (order++ <= e) term.push_back(v3);
    if (order++ <= e) term.push_back(v4);
    if (order++ <= e) term.push_back(v5);
    if (order++ <= e) term.push_back(v6);
    if (order++ <= e) term.push_back(v7);
  }
  Series<T>(int s, int e, T v1, T v2, T v3,
            T v4, T v5, T v6, T v7,
            T v8) :
    startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
    if (order++ <= e) term.push_back(v2);
    if (order++ <= e) term.push_back(v3);
    if (order++ <= e) term.push_back(v4);
    if (order++ <= e) term.push_back(v5);
    if (order++ <= e) term.push_back(v6);
    if (order++ <= e) term.push_back(v7);
    if (order++ <= e) term.push_back(v8);
  }
  Series<T>(int s, int e, T v1, T v2, T v3,
            T v4, T v5, T v6, T v7,
            T v8, T v9) :
    startOrder(s), endOrder(e) {
    short int order = s;

    if (order++ <= e) term.push_back(v1);
    if (order++ <= e) term.push_back(v2);
    if (order++ <= e) term.push_back(v3);
    if (order++ <= e) term.push_back(v4);
    if (order++ <= e) term.push_back(v5);
    if (order++ <= e) term.push_back(v6);
    if (order++ <= e) term.push_back(v7);
    if (order++ <= e) term.push_back(v8);
    if (order++ <= e) term.push_back(v9);
  }
  // and a catch-all for more complicated cases
  Series<T>(int s, int e, const std::vector<T>& v) : startOrder(s), endOrder(e) {
    short int order = s, i = 0;

    while (order++ <= e) term.push_back(v[i++]);
  }

  // Extract term of given order, can't return const T& because
  // of zero/NaN return values
#ifdef SWIG
%ignore operator[];
  #endif
  T const& operator[](int order) const {
    if (order < startOrder) return(zero); // term is truly absent
    if (order > endOrder)
       return(infinity); // undefined
    return term[order-startOrder];}

  // Same thing
  T const & operator()(int order) const {
    if (order < startOrder) return zero; // term is truly absent
    if (order > endOrder)
       return(infinity); // undefined
    return term[order-startOrder];}

  int leading() const {return startOrder;}
  int last() const {return endOrder;}
  std::vector<T> get_term() const {return term;};
#ifdef SWIG
	%extend {
		T __getitem__(int order) const {
			return (*$self)[order];
		}
		std::string __repr__() {
			//make sure you include sstream in the SWIG interface file
			std::ostringstream oss(std::ostringstream::out);
			oss << (*$self);
			return oss.str();
		}
		std::string __str__() {
			//make sure you include sstream in the SWIG interface file
			std::ostringstream oss(std::ostringstream::out);
			oss << (*$self);
			return oss.str();
		}
	}
#endif

private:
	static T zero;
	static T infinity;
};

typedef Series<C> CSeries;
typedef Series<CHP> CHPSeries;
typedef Series<CVHP> CVHPSeries;

template <class T> class SeriesC;
template <class T> SeriesC<T> operator*(const T& v1,const SeriesC<T>& s2);
template <class T> SeriesC<T> operator*(const SeriesC<T>& s2,const T& v1);

template <class T> class SeriesC : public Series<std::complex<T> > {
	friend SeriesC<R> to_double(SeriesC<R>);
	friend SeriesC<R> to_double(SeriesC<RHP>);
	friend SeriesC<R> to_double(SeriesC<RVHP>);
	friend SeriesC<RHP> to_HP(SeriesC<RVHP>);
	friend SeriesC<RHP> to_HP(SeriesC<RHP>);
	friend SeriesC<RHP> to_HP(SeriesC<R>);
	friend SeriesC<RVHP> to_VHP(SeriesC<RVHP>);
	friend SeriesC<RVHP> to_VHP(SeriesC<RHP>);
	friend SeriesC<RVHP> to_VHP(SeriesC<R>);
#ifdef BH_USE_GMP
	friend SeriesC<R> to_double(SeriesC<RGMP>);
	friend SeriesC<RHP> to_HP(SeriesC<RGMP>);
	friend SeriesC<RVHP> to_VHP(SeriesC<RGMP>);
#endif
#ifndef SWIG
	friend SeriesC<T> operator*<>(const T& v1,const SeriesC<T>& s2);
	friend SeriesC<T> operator*<>(const SeriesC<T>& s2,const T& v1);
#endif
public:
	SeriesC(){};
	 // Define a zero series object with given range of orders
	  SeriesC(int s, int e) :  Series<std::complex<T> >(s,e) {};
	  SeriesC(int s, int e, T v1) : Series<std::complex<T> >(s,e,v1) {};
	  SeriesC(int s, int e, T v1, T v2) : Series<std::complex<T> >(s,e,v1,v2) {};
	  SeriesC(int s, int e, T v1, T v2, T v3) : Series<std::complex<T> >(s,e,v1,v2,v3) {};
	  SeriesC(int s, int e, T v1, T v2, T v3, T v4) : Series<std::complex<T> >(s,e,v1,v2,v3,v4) {};
	  SeriesC(int s, int e, T v1, T v2, T v3, T v4, T v5) : Series<std::complex<T> >(s,e,v1,v2,v3,v4,v5) {};
	  SeriesC(int s, int e, T v1, T v2, T v3, T v4, T v5, T v6) : Series<std::complex<T> >(s,e,v1,v2,v3,v4,v5,v6) {};
	  SeriesC(int s, int e, T v1, T v2, T v3, T v4, T v5, T v6, T v7) : Series<std::complex<T> >(s,e,v1,v2,v3,v4,v5,v6,v7) {};
	  SeriesC(int s, int e, T v1, T v2, T v3, T v4, T v5, T v6, T v7, T v8) : Series<std::complex<T> >(s,e,v1,v2,v3,v4,v5,v6,v7,v8) {};
	  SeriesC(int s, int e, T v1, T v2, T v3, T v4, T v5, T v6, T v7, T v8, T v9) : Series<std::complex<T> >(s,e,v1,v2,v3,v4,v5,v6,v7,v8,v9) {};

	  SeriesC(int s, int e, const std::vector<std::complex<T> >& v) : Series<std::complex<T> >(s,e,v){};
	  SeriesC(Series<std::complex<T> > s) : Series<std::complex<T> >(s.leading(),s.last(),s.get_term()){};
};




#ifdef BH_USE_GMP
SeriesC<R> to_double(SeriesC<RGMP>);
SeriesC<RHP> to_HP(SeriesC<RGMP>);
SeriesC<RVHP> to_VHP(SeriesC<RGMP>);
#endif


}

#endif /* SeriesDefined */
