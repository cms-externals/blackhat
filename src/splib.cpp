/*!\file BlackHat/src/splib.cpp
\brief Implemetation for SpLib
This file implements the four-momenta and the spinor products. The conventions are the same as in the Mathematica package S\@M (aka SP\@M) so thatn a direct numerical comparaison is straight forward.
*/
#include <iostream>
#include <iomanip>
#include <complex>
#include "spinor.h"
#include "BH_typedefs.h"
#include "BH_utilities.h"
#include <cmath>
#include "BH_error.h"
#define _SAFER_LA_NORMALIZATION 0

using std::complex;

namespace BH {


//template <class T> momentum<T>::operator momentum<complex<T> > (){
//return momentum<complex<T> > (e,x,y,z);
//};

#if _DO_INLINE

#else
#include "spinor_inline.cpp"
#endif


template <class T> std::ostream& operator<<(std::ostream& s, const momentum<T>& p){
	return s << '(' << p.e << ',' << p.x <<',' << p.y <<',' << p.z << ')';
}

template <class T> Cmom<T> Cmom<T>::operator+=(const Cmom<T>& pp){
	p+=pp.P();
	La=lambda<T>(p);
	Lat=lambdat<T>(p);
	_type=_mt_unknown;
	return *this;
}
template <class T> Cmom<T> Cmom<T>::operator-=(const Cmom<T>& pp){
	p-=pp.P();
	La=lambda<T>(p);
	Lat=lambdat<T>(p);
	_type=_mt_unknown;
	return *this;
 }
template <class T> Cmom<T> operator*(const T& c,const Cmom<T>& p){
  if (c==T(0.)) return Cmom<T>(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),p.type());
  if (c>T(0.)) {
	  if (p.type()==_mt_massive) {
		  return Cmom<T>(complex<T>(c)*p.P(),_mt_massive);
	  }
	  else {
		  return Cmom<T>(complex<T>(c)*p.P(),sqrt(c)*p.L(),sqrt(c)*p.Lt(),p.type());
	  }
  }
  if (c<T(0.))
	  if (p.type()==_mt_massive) {
		  return Cmom<T>(complex<T>(c)*p.P(),_mt_massive);
	  }
	  else {
		  return Cmom<T>(complex<T>(c)*p.P(),sqrt(-c)*p.L(),(-sqrt(-c))*p.Lt(),p.type());
	  }

  _WARNING("no Cmom returned in  Cmom<T> operator*(const T& c,momentum<T> p), returned 0.");
  return Cmom<T>(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
}

template <class T> Cmom<T> Cmom<T>::operator*=(const T&  c){
	p*=complex<T>(c,0);
	La=lambda<T>(p);
	Lat=lambdat<T>(p);
	return *this;
}

template <class T> Cmom<T> operator/(const Cmom<T>& p,const T& c){
  if (c==T(0.)) {_WARNING("Division of a vector by zero"); throw BHerror("Momentum error");}
  if (c>T(0.))
	  if (p.type()==_mt_massive) {
		  return Cmom<T>(p.P()*(complex<T>(1,0)/c),_mt_massive);
	  }
	  else {
		  return Cmom<T>(p.P()*(complex<T>(1.0)/c),sqrt(T(1.)/c)*p.L(),sqrt(T(1.)/c)*p.Lt(),p.type());
	  }



  if (c<T(0.))
	  if (p.type()==_mt_massive) {
		  return Cmom<T>(p.P()*(complex<T>(1,0)/c),_mt_massive);
	  }
	  else {
		  return Cmom<T>((complex<T>(T(1.0)/c))*p.P(),sqrt(T(-1.)/c)*p.L(),-sqrt(T(-1.)/c)*p.Lt(),p.type());
	  }



  _WARNING("no Cmom returned in  Cmom<T> operator*(const T&c,momentum<T> p), returned 0.");
  return Cmom<T>(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
}

template <class T> Cmom<T> Cmom<T>::operator/=(const T& c){
	  if (c==T(0.)) {_WARNING("Division of a vector by zero"); throw BHerror("Momentum error"); }
	  if (c>T(0.)) {
		  p=p*(complex<T>(1.0)/c);
		  La=lambda<T>(sqrt(T(1.)/c)*La);
		  Lat=lambdat<T>(sqrt(T(1.)/c)*Lat);
	  }
	  if (c<T(0.)) {
		  p=p*(complex<T>(1.0)/c);
		  La=sqrt(T(-1.)/c)*La;
		  Lat=-sqrt(T(-1.)/c)*Lat;
	  }
  return *this;
}


template <class T> Cmom<T> operator/(const Cmom<T>& p,const complex<T>& c){
  if (c==complex<T>(0.,0.)) {_WARNING("Division of a vector by zero. Returned zero-momentum");throw BHerror("Momentum error");}
  if (c.imag()==T(0.)) return (T(1.0)/c.real())*p;

  if (p.type()==_mt_massive) {
	  return Cmom<T>((complex<T>(1.0,0.)/c)*p.P(),_mt_massive);
  }
  else {
	  return Cmom<T>((complex<T>(1.0,0.)/c)*p.P(),sqrt(complex<T>(1.0,0.)/c)*p.L(),sqrt(complex<T>(1.0,0.)/c)*p.Lt(),p.type());
  }
}

template <class T> Cmom<T> Cmom<T>::operator/=(const complex<T>& c){
	  if (c==complex<T>(0.,0.)) {_WARNING("Division of a vector by zero.");throw BHerror("Momentum error"); }
	  if (c.imag()==T(0.)) {
		  if (c.real()>T(0.)) {
			  p=p*(complex<T>(1.0)/c);
			  La=lambda<T>(sqrt(T(1.)/c.real())*La);
			  Lat=lambdat<T>(sqrt(T(1.)/c.real())*Lat);
		  }
		  if (c.real()<T(0.)) {
			  p=p*(complex<T>(1.0)/c.real());
			  La=sqrt(T(-1.)/c.real())*La;
			  Lat=-sqrt(T(-1.)/c.real())*Lat;
		  }
	  }
	  else {
		  p=(complex<T>(1.0,0.)/c)*p;
		  La=sqrt(complex<T>(1.0,0.)/c)*La;
		  Lat=sqrt(complex<T>(1.0,0.)/c)*Lat;
	  }
	  return *this;
}

template <class T> Cmom<T> Cmom<T>::operator*=(const complex<T>& c){
	  if (c==T(0.)) {
		  p=momentum<complex<T> >(complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.),complex<T>(0.,0.));
		  La=lambda<T>(0,0);
		  Lat=lambdat<T>(0,0);
	  }
	  p*=c;
	  if (type()!=_mt_massive) {
		  if (c.imag()==0. && c.real() < T(0.)) {
			  La=sqrt(-c)*La;
			  Lat=(-sqrt(-c))*Lat;
		  }
		  else //	  if (c>T(0.))
			    La=sqrt(c)*La;
			    Lat=sqrt(c)*Lat;
		  }
	return (*this);
}

template <class T> std::ostream& operator<<(std::ostream& s, const Cmom<T>& p){
	return s << '(' << p.p.E() << ',' << p.p.X() <<',' << p.p.Y() <<',' << p.p.Z() << ')';
}


template <class T> Cmom<T>::Cmom(const T& E,const T& X,const T& Y,const T& Z,momentum_type type) : p(complex<T>(E,0.),complex<T>(X,0.),complex<T>(Y,0.),complex<T>(Z,0.)) ,La(complex<T>(0.,0.),complex<T>(0.,0.)) ,Lat(complex<T>(0.,0.),complex<T>(0.,0.)), _type(type) {
	if ((type!=_mt_massive)) { LaLat(p,La,Lat);}
}
//! constructor
template <class T> Cmom<T>::Cmom(const complex<T>& E,const complex<T>& X,const complex<T>& Y,const complex<T>& Z,momentum_type type) : p(E,X,Y,Z) ,La(complex<T>(0.,0.),complex<T>(0.,0.)) ,Lat(complex<T>(0.,0.),complex<T>(0.,0.)) , _type(type) {
	if (type!=_mt_massive) {LaLat(p,La,Lat);}
}



//! constructor
template <class T> Cmom<T>::Cmom(const momentum<T>& pp,momentum_type type) :p(complex<T>(pp.E(),0.),complex<T>(pp.X(),0.),complex<T>(pp.Y(),0.),complex<T>(pp.Z(),0.)) , La(0,0), Lat(0,0), _type(type) {
	if (type!=_mt_massive) { LaLat(pp,La,Lat);};
}
//! constructor
template <class T> Cmom<T>::Cmom(const momentum<complex<T> >& pp,momentum_type type) : p(pp), _type(type)  , La(0,0), Lat(0,0) {  if (type!=_mt_massive) {
	LaLat(pp,La,Lat);}
}






template <class T> smatrix<T>::smatrix(const momentum<T>& p): a11(p.minus()) ,a22(p.plus()), a12(-p.X(),p.Y()), a21(-p.X(),-p.Y())  {
//	a11=p.minus();a22=p.plus();a12=-(p.X()-complex<T>(0.,1.)*p.Y());a21=-(p.X()+complex<T>(0.,1.)*p.Y());
}

template <class T> smatrix<T>::smatrix(const momentum<complex<T> >& p) : a11(p.minus()) ,a22(p.plus()), a12(-(p.X()-complex<T>(0.,1.)*p.Y())), a21(-(p.X()+complex<T>(0.,1.)*p.Y())) {
//	a11=p.minus();a22=p.plus();a12=-(p.X()-complex<T>(0.,1.)*p.Y());a21=-(p.X()+complex<T>(0.,1.)*p.Y());
}

//template <class T> smatrix<T> operator*(smatrix<T> m1,smatrix<T> m2){
//	return smatrix<T>(
//			m1.a11*m2.a11+m1.a12*m2.a21,
//			m1.a11*m2.a21+m1.a12*m2.a22,
//			m1.a21*m2.a11+m1.a12*m2.a21,
//			m1.a21*m2.a21+m1.a22*m2.a22);
//};



#if _SAFER_LA_NORMALIZATION
template <class T> spinor<T> la(momentum<complex<T> > p){
	complex<T> c1(-1.1,0.2);
	complex<T> c2(-0.7,0.41);
	complex<T> c3(-0.3,-0.87);
	complex<T> c4(0.6,-1.4);
	complex<T> pp=p.plus();
	complex<T> ptminus=p.X()-complex<T>(0.,1.)*p.Y();
	complex<T> ptplus=p.X()+complex<T>(0.,1.)*p.Y();
	complex<T> pm=p.minus();
	complex<T> AA=(pp + complex<T>(0.,1.2)*ptminus);
	complex<T> BB=(pp + complex<T>(0.,1.2)*ptplus);
	complex<T> CC=(complex<T>(0.,1.2)*pm + ptplus);
	complex<T> DD=( complex<T>(0.,1.2)*pm +ptminus);
	complex<T> norm=(c1*(AA*BB) +c2* (CC*DD) + c3*(AA*DD) + c4*(CC*BB))/(c1*pp + c2*pm +c3* ptminus +c4* ptplus);

	return spinor<T>(AA,CC)*(complex<T>(1.,0.)/sqrt(norm));
	}
#else

template <class T> spinor<T> la(const momentum<complex<T> >& p){

			if ((p.plus()*conj(p.plus())).real()<DeltaZero<T>())
				if ((p.minus()*conj(p.minus())).real()<DeltaZero<T>()) return spinor<T>((p.X()-complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()+complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
				else {complex<T> sqm=sqrt(p.minus()); return spinor<T>(complex<T>(1.,0.)/sqm*(p.X()-complex<T>(0.,1.)*p.Y()),sqm);}
			else {


#if _OLD_PHASE_CONVENTION
				complex<T> sqk1p=sqrt(p.plus());
				return spinor<T>(sqk1p,complex<T>(1.,0.)/sqk1p*(p.X()+complex<T>(0.,1.)*p.Y()));

#else
		T sqk1p=sqrt(abs(p.plus())); //T inv=T(1.)/sqk1p;
		return lambda<T>(sqk1p,sqk1p*(p.X()+complex<T>(0.,1.)*p.Y())/p.plus());
#endif

			}
		}
#endif /*(_SAFER_LA_NORMALIZATION)*/


//template <class T> spinor<T> la(momentum<T> *p){
//	return la(*p);
//}
//
//template <class T> spinor<T> la(momentum<complex<T> > *p){
//	return la(*p);
//}

#if _SAFER_LA_NORMALIZATION

template <class T> spinor<T> lat(momentum<complex<T> > p){
	complex<T> c1(-1.1,0.2);
	complex<T> c2(-0.7,0.41);
	complex<T> c3(-0.3,-0.87);
	complex<T> c4(0.6,-1.4);
	complex<T> pp=p.plus();
	complex<T> ptminus=p.X()-complex<T>(0.,1.)*p.Y();
	complex<T> ptplus=p.X()+complex<T>(0.,1.)*p.Y();
	complex<T> pm=p.minus();
	complex<T> AA=(pp + complex<T>(0.,1.2)*ptminus);
	complex<T> BB=(pp + complex<T>(0.,1.2)*ptplus);
	complex<T> CC=(complex<T>(0.,1.2)*pm + ptplus);
	complex<T> DD=(complex<T>(0.,1.2)*pm + ptminus);
	complex<T> norm=(c1*(AA*BB) +c2*(CC*DD) + c3* (AA*DD) +c4* (CC*BB))/(c1*pp + c2*pm +c3* ptminus + c4*ptplus);
	return spinor<T>(BB,DD)*(complex<T>(1.,0.)/(sqrt(norm)));
	}
#else

template <class T> spinor<T> lat(const momentum<complex<T> >& p){
	if ((p.plus()*conj(p.plus())).real()<DeltaZero<T>())
		if ((p.minus()*conj(p.minus())).real()<DeltaZero<T>()) return spinor<T>((p.X()+complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()-complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
		else {complex<T> sqm=sqrt(p.minus()); return spinor<T>(complex<T>(1.,0.)/sqm*(p.X()+complex<T>(0.,1.)*p.Y()),sqm);}
	else {
#if _OLD_PHASE_CONVENTION
		complex<T> sqk1p=sqrt(p.plus());
		return spinor<T>(sqk1p,complex<T>(1.,0.)/sqk1p*(p.X()-complex<T>(0.,1.)*p.Y()));
#else
		T sqk1p=sqrt(abs(p.plus())); T inv=T(1.)/sqk1p;
		return lambdat<T>(p.plus()*inv,inv*(p.X()-complex<T>(0.,1.)*p.Y()));
#endif

	}
}

#endif /*(_SAFER_LA_NORMALIZATION)*/



template <class T> std::ostream& operator<<(std::ostream& s, const spinor<T>& x){
	return s << '(' << x.a1 << ',' << x.a2  << ')';
}
template <class T> std::ostream& operator<<(std::ostream& s, const smatrix<T>& x){
	return s << '(' << x.a11 << ',' << x.a12  << ','<<std::endl  << x.a21  <<',' << x.a22  << ')';
}


//! spinor product < * | * | * |  * > for real momenta (overloaded)
/**
 \param p1 first (real) momentum
 \param p2 first slashed (real) momentum
 \param p3 second slashed (real) momentum
 \param p2 last (real) momentum
\return spinor product <p1|p2|p3|p4>
\sa spaa(momentum<T> p1,momentum<T> p2)
\sa spaa(momentum<complex<T> > p1,momentum<complex<T> > p2)
*/

template <class T> complex<T> spaa(const momentum<T>& p1,const momentum<T>& q1,const momentum<T>& q2,const momentum<T>& p2){
	spinor<T> b1=smatrix<T>(q1)*la(p1);
	spinor<T> b2=smatrix<T>(q2)*la(p2);
	return b1*b2;
}

//! spinor product < * | * | * |  * > for complex momenta (overloaded)
/**
 \param p1 first (complex) momentum
 \param p2 first slashed (complex) momentum
 \param p3 second slashed (complex) momentum
 \param p2 last (complex) momentum
\return spinor product <p1|p2|p3|p4>
\sa spaa(momentum<T> p1,momentum<T> p2)
\sa spaa(momentum<complex<T> > p1,momentum<complex<T> > p2)
*/

template <class T> complex<T> spaa(const momentum<complex<T> >& p1,const momentum<complex<T> >& q1,const momentum<complex<T> >& q2,const momentum<complex<T> >& p2){
	spinor<T> b1=smatrix<T>(q1)*la(p1);
	spinor<T> b2=smatrix<T>(q2)*la(p2);
	return b1*b2;
}

//template <class T>
//	complex<T> spaa(momentum<T> &p1,momentum<T> &p2){
//	return la(p1)*la(p2);
//}






//! spinor product [ * | * | * |  * ] for real momenta (overloaded)
/**
 \param p1 first (real) momentum
 \param p2 first slashed (real) momentum
 \param p3 second slashed (real) momentum
 \param p2 last (real) momentum
\return spinor product [p1|p2|p3|p4]
\sa spbb(momentum<T> p1,momentum<T> p2)
\sa spbb(momentum<complex<T> > p1,momentum<complex<T> > p2)
*/

template <class T> complex<T> spbb(const momentum<T>& p1,const momentum<T>& q1,const momentum<T>& q2,const momentum<T>& p2){
	spinor<T> b1=lat(p1)*smatrix<T>(q1);
	spinor<T> b2=lat(p2)*smatrix<T>(q2);
	return -(b1*b2);
		}
//! spinor product [ * | * | * |  * ] for complex momenta (overloaded)
/**
 \param p1 first (complex) momentum
 \param p2 first slashed (complex) momentum
 \param p3 second slashed (complex) momentum
 \param p2 last (complex) momentum
\return spinor product [p1|p2|p3|p4]
\sa spbb(momentum<T> p1,momentum<T> p2)
\sa spbb(momentum<complex<T> > p1,momentum<complex<T> > p2)
*/

template <class T> complex<T> spbb(const momentum<complex<T> >& p1,const momentum<complex<T> >& q1,const momentum<complex<T> >& q2,const momentum<complex<T> >& p2){
	spinor<T> b1=lat(p1)*smatrix<T>(q1);
	spinor<T> b2=lat(p2)*smatrix<T>(q2);
	return -(b1*b2);
}






//     EXPICIT INSTANCIATION

template class momentum<double>;
//template class momentum<double> operator+(const momentum<double>&,const momentum<double>&);
//template class momentum<double> operator-(const momentum<double>&,const momentum<double>&);
//template class momentum<double> operator-(const momentum<double>&);

template class momentum<complex<double> >;
//template class momentum<complex<double> > operator+(const momentum<complex<double> >&,const momentum<complex<double> >&);
//template class momentum<complex<double> > operator-(const momentum<complex<double> >&,const momentum<complex<double> >&);
//template class momentum<complex<double> > operator-(const momentum<complex<double> >&);

//template class momentum<double> operator*(const double&,const momentum<double>&);
//template class momentum<double> operator*(const momentum<double>&,const double&);
//template class momentum<double> operator/(const momentum<double>&,const double&);
//template double operator*(const momentum<double>&,const momentum<double>&);


//template class momentum<complex<double> > operator*(const complex<double>&,const momentum<complex<double> >&);
//template class momentum<complex<double> > operator*(const momentum<complex<double> >&,const complex<double>&);
//template class momentum<complex<double> > operator/(const momentum<complex<double> >&,const complex<double>&);
//template complex<double> operator*(const momentum<complex<double> >&,const momentum<complex<double> >&);




template class spinor<double>;
template class lambda<double>;

template class Cmom<double>;
//template class Cmom<double> operator+(const Cmom<double>&,const Cmom<double>&);
//template class Cmom<double> operator-(const Cmom<double>&,const Cmom<double>&);
//template class Cmom<double> operator-(const Cmom<double>&);
template Cmom<double> operator*(const double&,const Cmom<double>&);
//template class Cmom<double> operator*(const Cmom<double>&,const double&);
//template class Cmom<double> operator*(const Cmom<double>&,const complex<double>&);
//template class Cmom<double> operator*(const complex<double>&,const Cmom<double>&);
template Cmom<double> operator/(const Cmom<double>&,const double&);
template Cmom<double> operator/(const Cmom<double>&,const complex<double>&);

//template complex<double> operator*(const Cmom<double>&,const Cmom<double>&);
template std::ostream& operator<<(std::ostream& s, const Cmom<double>& p);

template class smatrix<double>;


template std::ostream& operator<<(std::ostream& s, const momentum<double>& p);
template std::ostream& operator<<(std::ostream& s, const momentum<complex<double> >& p);
template std::ostream& operator<<(std::ostream& s, const smatrix<double>& p);
template std::ostream& operator<<(std::ostream& s, const spinor<double>& p);


//template class spinor<double> la(const momentum<double>&);
//template class spinor<double> lat(const momentum<double>&);

template spinor<double> la(const momentum<complex<double> >&);
template spinor<double> lat<double>(const momentum<complex<double> >&);

//template class spinor<double> operator+(const spinor<double>&,const spinor<double>&);
//template class spinor<double> operator-(const spinor<double>&,const spinor<double>&);
//template class spinor<double> operator-(const spinor<double>&);
//template class complex<double> operator*(const spinor<double>&,const spinor<double>&);
//template class spinor<double> operator*(const smatrix<double>& m,const spinor<double>& b);
//template class spinor<double> operator*(const spinor<double>& b,const smatrix<double>& m);
//
//template class  spinor<double> operator*(const spinor<double>&,const complex<double>&);
//template class  spinor<double> operator*(const complex<double>&,const spinor<double>&);
//template class  spinor<double> operator*(const spinor<double>&,const double&);
//template class  spinor<double> operator*(const double&,const spinor<double>&);


//template class complex<double> spaa(const momentum<double>&,const momentum<double>&);
template complex<double> spaa(const momentum<double>&,const momentum<double>&,const momentum<double>&,const momentum<double>&);
//template class complex<double> spaa(const momentum<complex<double> >&,const momentum<complex<double> >&);
template complex<double> spaa(const momentum<complex<double> >&,const momentum<complex<double> >&,const momentum<complex<double> >&,const momentum<complex<double> >&);

//template class complex<double> spbb(const momentum<double>&,const momentum<double>&);

template complex<double> spbb(const momentum<double>&,const momentum<double>&,const momentum<double>&,const momentum<double>&);

//template class complex<double> spbb(const momentum<complex<double> >&,const momentum<complex<double> >&);
template complex<double> spbb(const momentum<complex<double> >&,const momentum<complex<double> >&,const momentum<complex<double> >&,const momentum<complex<double> >&);


//template class complex<double> spab(const momentum<double>&,const momentum<double>&,const momentum<double>&);


//template class complex<double> spab(const momentum<complex<double> >&,const momentum<complex<double> >&,const momentum<complex<double> >&);

//template double s<double>(const momentum<double>&,const momentum<double>&);
//template double s<double>(const momentum<double>&,const momentum<double>&,const momentum<double>&);
//template double s<double>(const momentum<double>&,const momentum<double>&,const momentum<double>&,const momentum<double>&);
//
//template complex<double> s(const momentum<complex<double> >&,const momentum<complex<double> >&);
//template complex<double> s(const momentum<complex<double> >&,const momentum<complex<double> >&,const momentum<complex<double> >&);
//template complex<double> s(const momentum<complex<double> >&,const momentum<complex<double> >&,const momentum<complex<double> >&,const momentum<complex<double> >&);
//
//template class lambda<double> operator+(const lambda<double>& b1,const lambda<double>& b2);
//template class lambda<double> operator-(const lambda<double>& b1,const lambda<double>& b2);
//template class lambda<double> operator-(const lambda<double>& b);
//template class lambdat<double> operator+(const lambdat<double>& b1,const lambdat<double>& b2);
//template class lambdat<double> operator-(const lambdat<double>& b1,const lambdat<double>& b2);
//template class lambdat<double> operator-(const lambdat<double>& b);
//template class  lambda<double> operator*(const complex<double>& c,const lambda<double>& b);
//template class  lambda<double> operator*(const lambda<double>& b,const complex<double>& c);
//template class  lambda<double> operator*(const double& c,const lambda<double>& b);
//template class lambda<double> operator*(const lambda<double>& b,const double& c);
//template class  lambdat<double> operator*(const complex<double>& c,const lambdat<double>& b);
//template class  lambdat<double> operator*(const lambdat<double>& b,const complex<double>& c);
//template class  lambdat<double> operator*(const double& c,const lambdat<double>& b);
//template class lambdat<double> operator*(const lambdat<double>& b,const double& c);
//template class lambdat<double> operator*(const smatrix<double>& m,const lambda<double>& b);
//template class lambdat<double> operator*(const lambda<double>& b,const smatrix<double>& m);
//template class lambda<double> operator*(const smatrix<double>& m,const lambdat<double>& b);
//template class lambda<double> operator*(const lambdat<double>& b,const smatrix<double>& m);
//template class complex<double> operator*(const lambdat<double>& b1,const lambdat<double>& b2);
//
//template class momentum<complex<double> > PfLLt(const lambdat<double>&,const lambda<double>&);



//     EXPICIT INSTANCIATION RHP

template class momentum<RHP>;
//template class momentum<RHP> operator+(const momentum<RHP>&,const momentum<RHP>&);
//template class momentum<RHP> operator-(const momentum<RHP>&,const momentum<RHP>&);
//template class momentum<RHP> operator-(const momentum<RHP>&);

template class momentum<complex<RHP> >;
//template class momentum<complex<RHP> > operator+(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);
//template class momentum<complex<RHP> > operator-(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);
//template class momentum<complex<RHP> > operator-(const momentum<complex<RHP> >&);

//template class momentum<RHP> operator*(const RHP&,const momentum<RHP>&);
//template class momentum<RHP> operator*(const momentum<RHP>&,const RHP&);
//template class momentum<RHP> operator/(const momentum<RHP>&,const RHP&);
//template RHP operator*(const momentum<RHP>&,const momentum<RHP>&);


//template class momentum<complex<RHP> > operator*(const complex<RHP>&,const momentum<complex<RHP> >&);
//template class momentum<complex<RHP> > operator*(const momentum<complex<RHP> >&,const complex<RHP>&);
//template class momentum<complex<RHP> > operator/(const momentum<complex<RHP> >&,const complex<RHP>&);
//template complex<RHP> operator*(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);

template class Cmom<RHP>;
//template class Cmom<RHP> operator+(const Cmom<RHP>&,const Cmom<RHP>&);
//template class Cmom<RHP> operator-(const Cmom<RHP>&,const Cmom<RHP>&);
//template class Cmom<RHP> operator-(const Cmom<RHP>&);
template Cmom<RHP> operator*(const RHP&,const Cmom<RHP>&);
//template class Cmom<RHP> operator*(const Cmom<RHP>&,const RHP&);
//template class Cmom<RHP> operator*(const Cmom<RHP>&,const complex<RHP>&);
//template class Cmom<RHP> operator*(const complex<RHP>&,const Cmom<RHP>&);
template Cmom<RHP> operator/(const Cmom<RHP>&,const RHP&);
template Cmom<RHP> operator/(const Cmom<RHP>&,const complex<RHP>&);

//template complex<RHP> operator*(const Cmom<RHP>&,const Cmom<RHP>&);
template std::ostream& operator<<(std::ostream& s, const Cmom<RHP>& p);

template class spinor<RHP>;
template class smatrix<RHP>;


template std::ostream& operator<<(std::ostream& s, const momentum<RHP>& p);
template std::ostream& operator<<(std::ostream& s, const momentum<complex<RHP> >& p);
template std::ostream& operator<<(std::ostream& s, const smatrix<RHP>& p);
template std::ostream& operator<<(std::ostream& s, const spinor<RHP>& p);


//template class spinor<RHP> la(const momentum<RHP>&);
//template class spinor<RHP> lat(const momentum<RHP>&);


template spinor<RHP> la(const momentum<complex<RHP> >&);

template spinor<RHP> lat<RHP>(const momentum<complex<RHP> >&);




//template class spinor<RHP> operator-(const spinor<RHP>&,const spinor<RHP>&);
//template class spinor<RHP> operator-(const spinor<RHP>&);
//template class spinor<RHP> operator+(const spinor<RHP>&,const spinor<RHP>&);
template complex<RHP> operator*(const spinor<RHP>&,const spinor<RHP>&);

//template class spinor<RHP> operator*(const smatrix<RHP>& m,const spinor<RHP>& b);
//template class spinor<RHP> operator*(const spinor<RHP>& b,const smatrix<RHP>& m);

//template class  spinor<RHP> operator*(const spinor<RHP>&,const complex<RHP>&);
//template class  spinor<RHP> operator*(const complex<RHP>&,const spinor<RHP>&);
//template class  spinor<RHP> operator*(const spinor<RHP>&,const RHP&);
//template class  spinor<RHP> operator*(const RHP&,const spinor<RHP>&);


//template class complex<RHP> spaa(const momentum<RHP>&,const momentum<RHP>&);
template complex<RHP> spaa(const momentum<RHP>&,const momentum<RHP>&,const momentum<RHP>&,const momentum<RHP>&);

//template class complex<RHP> spaa(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);
template  complex<RHP> spaa(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&,const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);

//template class complex<RHP> spbb(const momentum<RHP>&,const momentum<RHP>&);

template complex<RHP> spbb(const momentum<RHP>&,const momentum<RHP>&,const momentum<RHP>&,const momentum<RHP>&);

//template class complex<RHP> spbb(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);
template complex<RHP> spbb(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&,const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);


//template class complex<RHP> spab(const momentum<RHP>&,const momentum<RHP>&,const momentum<RHP>&);


//template class complex<RHP> spab(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);

//template RHP s<RHP>(const momentum<RHP>&,const momentum<RHP>&);
//template RHP s<RHP>(const momentum<RHP>&,const momentum<RHP>&,const momentum<RHP>&);
//template RHP s<RHP>(const momentum<RHP>&,const momentum<RHP>&,const momentum<RHP>&,const momentum<RHP>&);
//
//template complex<RHP> s(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);
//template complex<RHP> s(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);
//template complex<RHP> s(const momentum<complex<RHP> >&,const momentum<complex<RHP> >&,const momentum<complex<RHP> >&,const momentum<complex<RHP> >&);

//template class lambda<RHP> operator+(const lambda<RHP>& b1,const lambda<RHP>& b2);
//template class lambda<RHP> operator-(const lambda<RHP>& b1,const lambda<RHP>& b2);
//template class lambda<RHP> operator-(const lambda<RHP>& b);
//template class lambdat<RHP> operator+(const lambdat<RHP>& b1,const lambdat<RHP>& b2);
//template class lambdat<RHP> operator-(const lambdat<RHP>& b1,const lambdat<RHP>& b2);
//template class lambdat<RHP> operator-(const lambdat<RHP>& b);
//template class  lambda<RHP> operator*(const complex<RHP>& c,const lambda<RHP>& b);
//template class  lambda<RHP> operator*(const lambda<RHP>& b,const complex<RHP>& c);
//template class  lambda<RHP> operator*(const RHP& c,const lambda<RHP>& b);
//template class lambda<RHP> operator*(const lambda<RHP>& b,const RHP& c);
//template class  lambdat<RHP> operator*(const complex<RHP>& c,const lambdat<RHP>& b);
//template class  lambdat<RHP> operator*(const lambdat<RHP>& b,const complex<RHP>& c);
//template class  lambdat<RHP> operator*(const RHP& c,const lambdat<RHP>& b);
//template class lambdat<RHP> operator*(const lambdat<RHP>& b,const RHP& c);
//template class lambdat<RHP> operator*(const smatrix<RHP>& m,const lambda<RHP>& b);
//template class lambdat<RHP> operator*(const lambda<RHP>& b,const smatrix<RHP>& m);
//template class lambda<RHP> operator*(const smatrix<RHP>& m,const lambdat<RHP>& b);
//template class lambda<RHP> operator*(const lambdat<RHP>& b,const smatrix<RHP>& m);
//template class complex<RHP> operator*(const lambdat<RHP>& b1,const lambdat<RHP>& b2);
//
//template class momentum<complex<RHP> > PfLLt(const lambdat<RHP>&,const lambda<RHP>&);




//     EXPICIT INSTANCIATION RVHP

template class momentum<RVHP>;
//template class momentum<RVHP> operator+(const momentum<RVHP>&,const momentum<RVHP>&);
//template class momentum<RVHP> operator-(const momentum<RVHP>&,const momentum<RVHP>&);
//template class momentum<RVHP> operator-(const momentum<RVHP>&);

template class momentum<complex<RVHP> >;
//template class momentum<complex<RVHP> > operator+(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);
//template class momentum<complex<RVHP> > operator-(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);
//template class momentum<complex<RVHP> > operator-(const momentum<complex<RVHP> >&);

//template class momentum<RVHP> operator*(const RVHP&,const momentum<RVHP>&);
//template class momentum<RVHP> operator*(const momentum<RVHP>&,const RVHP&);
//template class momentum<RVHP> operator/(const momentum<RVHP>&,const RVHP&);
//template RVHP operator*(const momentum<RVHP>&,const momentum<RVHP>&);


//template class momentum<complex<RVHP> > operator*(const complex<RVHP>&,const momentum<complex<RVHP> >&);
//template class momentum<complex<RVHP> > operator*(const momentum<complex<RVHP> >&,const complex<RVHP>&);
//template class momentum<complex<RVHP> > operator/(const momentum<complex<RVHP> >&,const complex<RVHP>&);
//template complex<RVHP> operator*(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);

template class Cmom<RVHP>;
//template class Cmom<RVHP> operator+(const Cmom<RVHP>&,const Cmom<RVHP>&);
//template class Cmom<RVHP> operator-(const Cmom<RVHP>&,const Cmom<RVHP>&);
//template class Cmom<RVHP> operator-(const Cmom<RVHP>&);
template Cmom<RVHP> operator*(const RVHP&,const Cmom<RVHP>&);
//template class Cmom<RVHP> operator*(const Cmom<RVHP>&,const RVHP&);
//template class Cmom<RVHP> operator*(const Cmom<RVHP>&,const complex<RVHP>&);
//template class Cmom<RVHP> operator*(const complex<RVHP>&,const Cmom<RVHP>&);
template Cmom<RVHP> operator/(const Cmom<RVHP>&,const RVHP&);
template Cmom<RVHP> operator/(const Cmom<RVHP>&,const complex<RVHP>&);

//template complex<RVHP> operator*(const Cmom<RVHP>&,const Cmom<RVHP>&);
template std::ostream& operator<<(std::ostream& s, const Cmom<RVHP>& p);

template class spinor<RVHP>;
template class smatrix<RVHP>;


template std::ostream& operator<<(std::ostream& s, const momentum<RVHP>& p);
template std::ostream& operator<<(std::ostream& s, const momentum<complex<RVHP> >& p);
template std::ostream& operator<<(std::ostream& s, const smatrix<RVHP>& p);
template std::ostream& operator<<(std::ostream& s, const spinor<RVHP>& p);


//template class spinor<RVHP> la(const momentum<RVHP>&);

//template class spinor<RVHP> lat(const momentum<RVHP>&);


template spinor<RVHP> la(const momentum<complex<RVHP> >&);
//template class spinor<RVHP> la(momentum<complex<RVHP> >*);
template spinor<RVHP> lat<RVHP>(const momentum<complex<RVHP> >&);




//template class spinor<RVHP> operator-(const spinor<RVHP>&,const spinor<RVHP>&);
//template class spinor<RVHP> operator-(const spinor<RVHP>&);
//template class spinor<RVHP> operator+(const spinor<RVHP>&,const spinor<RVHP>&);
template complex<RVHP> operator*(const spinor<RVHP>&,const spinor<RVHP>&);

//template class spinor<RVHP> operator*(const smatrix<RVHP>& m,const spinor<RVHP>& b);
//template class spinor<RVHP> operator*(const spinor<RVHP>& b,const smatrix<RVHP>& m);
//
//template class  spinor<RVHP> operator*(const spinor<RVHP>&,const complex<RVHP>&);
//template class  spinor<RVHP> operator*(const complex<RVHP>&,const spinor<RVHP>&);
//template class  spinor<RVHP> operator*(const spinor<RVHP>&,const RVHP&);
//template class  spinor<RVHP> operator*(const RVHP&,const spinor<RVHP>&);
//
//
//template class complex<RVHP> spaa(const momentum<RVHP>&,const momentum<RVHP>&);
template complex<RVHP> spaa(const momentum<RVHP>&,const momentum<RVHP>&,const momentum<RVHP>&,const momentum<RVHP>&);
//template class complex<RVHP> spaa(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);
template complex<RVHP> spaa(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);

//template class complex<RVHP> spbb(const momentum<RVHP>&,const momentum<RVHP>&);

template complex<RVHP> spbb(const momentum<RVHP>&,const momentum<RVHP>&,const momentum<RVHP>&,const momentum<RVHP>&);

//template class complex<RVHP> spbb(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);
template complex<RVHP> spbb(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);


//template class complex<RVHP> spab(const momentum<RVHP>&,const momentum<RVHP>&,const momentum<RVHP>&);


//template class complex<RVHP> spab(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);

//template RVHP s<RVHP>(const momentum<RVHP>&,const momentum<RVHP>&);
//template RVHP s<RVHP>(const momentum<RVHP>&,const momentum<RVHP>&,const momentum<RVHP>&);
//template RVHP s<RVHP>(const momentum<RVHP>&,const momentum<RVHP>&,const momentum<RVHP>&,const momentum<RVHP>&);
//
//template complex<RVHP> s(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);
//template complex<RVHP> s(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);
//template complex<RVHP> s(const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&,const momentum<complex<RVHP> >&);
//
//template class lambda<RVHP> operator+(const lambda<RVHP>& b1,const lambda<RVHP>& b2);
//template class lambda<RVHP> operator-(const lambda<RVHP>& b1,const lambda<RVHP>& b2);
//template class lambda<RVHP> operator-(const lambda<RVHP>& b);
//template class lambdat<RVHP> operator+(const lambdat<RVHP>& b1,const lambdat<RVHP>& b2);
//template class lambdat<RVHP> operator-(const lambdat<RVHP>& b1,const lambdat<RVHP>& b2);
//template class lambdat<RVHP> operator-(const lambdat<RVHP>& b);
//template class  lambda<RVHP> operator*(const complex<RVHP>& c,const lambda<RVHP>& b);
//template class  lambda<RVHP> operator*(const lambda<RVHP>& b,const complex<RVHP>& c);
//template class  lambda<RVHP> operator*(const RVHP& c,const lambda<RVHP>& b);
//template class lambda<RVHP> operator*(const lambda<RVHP>& b,const RVHP& c);
//template class  lambdat<RVHP> operator*(const complex<RVHP>& c,const lambdat<RVHP>& b);
//template class  lambdat<RVHP> operator*(const lambdat<RVHP>& b,const complex<RVHP>& c);
//template class  lambdat<RVHP> operator*(const RVHP& c,const lambdat<RVHP>& b);
//template class lambdat<RVHP> operator*(const lambdat<RVHP>& b,const RVHP& c);
//template class lambdat<RVHP> operator*(const smatrix<RVHP>& m,const lambda<RVHP>& b);
//template class lambdat<RVHP> operator*(const lambda<RVHP>& b,const smatrix<RVHP>& m);
//template class lambda<RVHP> operator*(const smatrix<RVHP>& m,const lambdat<RVHP>& b);
//template class lambda<RVHP> operator*(const lambdat<RVHP>& b,const smatrix<RVHP>& m);
//template class complex<RVHP> operator*(const lambdat<RVHP>& b1,const lambdat<RVHP>& b2);
//
//template class momentum<complex<RVHP> > PfLLt(const lambdat<RVHP>&,const lambda<RVHP>&);



}
