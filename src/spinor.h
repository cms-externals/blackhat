/*!\file BlackHat/src/spinor.h
\brief Header file for SpLib
*/
#ifndef _H_spinor
#define _H_spinor

#include <iostream>
#include <iomanip>
#include <complex>
#include "qd_suppl.h"


#define _OLD_PHASE_CONVENTION 0   // with 1 the old phase convention is used
#define _DO_INLINE 1
namespace BH {


#ifdef SWIG
using namespace std;
#endif

enum momentum_type {_mt_massive,_mt_massless,_mt_unknown};


template <class T> class momentum ;
template <class T> class Cmom ;

template <class T> class spinor;
template <class T> class smatrix;
template <class T> class lambda;
template <class T> class lambdat;

template <class T> momentum<T> operator+(const momentum<T>& p1,const momentum<T>& p2);
template <class T> momentum<T> operator-(const momentum<T>& p1,const momentum<T>& p2);
template <class T> momentum<T> operator-(const momentum<T>&);
template <class T> T operator*(const momentum<T>& p1,const momentum<T>& p2);
template <class T> momentum<T> operator*(const T& c ,const momentum<T>& p);
template <class T> momentum<T> operator*(const momentum<T>& p,const T& c );
template <class T> momentum<T> operator/(const momentum<T>& p,const T& c );

template <class T> Cmom<T> operator+(const Cmom<T>&,const Cmom<T>&);
template <class T> Cmom<T> operator-(const Cmom<T>&,const Cmom<T>&);
template <class T> Cmom<T> operator-(const Cmom<T>&);
template <class T> std::complex<T> operator*(const Cmom<T>&,const Cmom<T>&);
template <class T> Cmom<T> operator*(const T&,const Cmom<T>&);
template <class T> Cmom<T> operator*(const Cmom<T>&,const T&);
template <class T> Cmom<T> operator*(const std::complex<T>&,const Cmom<T>&);
template <class T> Cmom<T> operator*(const Cmom<T>&,const std::complex<T>&);
template <class T> std::ostream& operator<<(std::ostream&, const Cmom<T>& );
template <class T> Cmom<T> operator/(const Cmom<T>& p,const T& c);
template <class T> Cmom<T> operator/(const Cmom<T>& p,const std::complex<T>& c);
template <class T> Cmom<T> operator/=(const Cmom<T>& p,const std::complex<T>& c);



template <class T> spinor<T> operator*(const spinor<T>&,const std::complex<T>&);
template <class T> spinor<T> operator*(const std::complex<T>&,const spinor<T>&);
template <class T> spinor<T> operator*(const spinor<T>&,const T&);
template <class T> spinor<T> operator*(const T&,const spinor<T>&);

template <class T> std::ostream& operator<<(std::ostream&, const spinor<T>& );

template <class T> std::complex<T> operator*(const spinor<T>&,const spinor<T>&);
template <class T> spinor<T> operator+(const spinor<T>&, const spinor<T>&);
template <class T> spinor<T> operator-(const spinor<T>&, const spinor<T>&);
template <class T> spinor<T> operator-(const spinor<T>&);
template <class T> spinor<T> operator*(const spinor<T>&,const smatrix<T>&);
template <class T> spinor<T> operator*(const smatrix<T>&,const spinor<T>&);
//template <class T> smatrix<T> operator*(smatrix<T>,smatrix<T>);

template <class T> std::ostream& operator<<(std::ostream&, const momentum<T>& );
#ifndef SWIG
template <class T> void Spinor_to_momentum(const lambdat<T>& l1,const lambda<T>& l2, momentum<std::complex<T> >& ret);
template <class T> void Momentum_to_spinor(const momentum<std::complex<T> >& p,lambda<T>& la,lambdat<T>& lat);
#endif

//! Template class for momenta
/**
 The momentum class has overloaded operators for (m1 and m2 are of type momentum<T>)
 \n - addition: m1+m2
 \n - subtraction: m1-m2, -m1+m2
 \n - scalar multiplication: m1*c, c*m1 where c is of the type T (or can be converted to it)
 \n - Minkovsky product: m1*m2 (with signature (1,-1,-1,-1))

 the momenta can be displayed using standard ostreams */

template <class T> class momentum {
protected:
	T e;
	T x;
	T y;
	T z;
    //friend complex<double> spaa(momentum p1,momentum p2);
    //friend ostream& operator<< (ostream& os, const momentum<type>& p);
public:
#ifndef SWIG
    friend momentum<T> operator+<>(const momentum<T>&,const momentum<T>&);
    friend momentum<T> operator-<>(const momentum<T>&,const momentum<T>&);
    friend momentum<T> operator-<>(const momentum<T>&);
    friend T operator*<>(const momentum<T>&,const momentum<T>&);
    friend momentum<T> operator*<>(const T&,const momentum<T>&);
    friend momentum<T> operator*<>(const momentum<T>&,const T&);
    friend momentum<T> operator/<>(const momentum<T>&,const T&);
    friend std::ostream& operator<<<>(std::ostream&, const momentum<T>& );
#endif

    momentum(){};
    momentum<T>(const T& E,const T& X,const T& Y,const T& Z): e(E), x(X), y(Y), z(Z) {};
    momentum<T> operator+=(const momentum<T>& p);
    momentum<T> operator-=(const momentum<T>& p);
    momentum<T> operator*=(const T& c);
    momentum<T> operator/=(const T& c);
    template <class U> momentum(const U& E, const U& X, const U& Y, const U& Z): e(T(E)), x(T(X)), y(T(Y)), z(T(Z)){}

    //operator momentum<complex<T> > ();
    //! \return the (Minkowski) square of the momentum p^2=e^2-x^2-y^2-z^2;
    T square() const { return e*e-x*x-y*y-z*z;};
    T plus() const { return e+z;};
    T minus() const { return e-z;};
    //! \return the x component of the vector
    const T& X() const { return x;};
    //! \return the y component of the vector
    const T& Y() const { return y;};
    //! \return the z component of the vector
    const T& Z() const { return z;};
    //! \return the E component of the vector
    const T& E() const { return e;};
    //template<typename _Tp, typename _CharT, class _Traits>
    //friend std::basic_ostream<_CharT, _Traits>&
    //    operator<<(std::basic_ostream<_CharT, _Traits>& , const momentum<_Tp>& );
#ifndef SWIG    
    void add_to(const momentum<T>& p) {e+=p.E();x+=p.X();y+=p.Y();z+=p.Z();};
    void set_to(const momentum<T>& p) {e=p.E();x=p.X();y=p.Y();z=p.Z();};
    void set_to(const T& E, const T& X, const T& Y, const T& Z) {e=E;x=X;y=Y;z=Z;};
    void mult_by(const T& p) {e*=p;x*=p;y*=p;z*=p;};
#endif
};
    
//! Template class for complex momenta Cmom

template <class T> class Cmom {
public:
protected:
	momentum<std::complex<T> > p;
	lambda<T> La;
	lambdat<T> Lat;
	momentum_type _type;
#ifndef SWIG
	friend Cmom<T> operator+<>(const Cmom<T>&,const Cmom<T>&);
	friend Cmom<T> operator-<>(const Cmom<T>&,const Cmom<T>&);
	friend Cmom<T> operator-<>(const Cmom<T>&);
	friend std::complex<T> operator*<>(const Cmom<T>&,const Cmom<T>&);
	friend Cmom<T> operator*<>(const T&,const Cmom<T>&);
	friend Cmom<T> operator*<>(const std::complex<T>&,const Cmom<T>&);
	friend Cmom<T> operator*<>(const Cmom<T>&,const T&);
	friend Cmom<T> operator*<>(const Cmom<T>&,const std::complex<T>&);
	friend std::ostream& operator<<<>(std::ostream&, const Cmom<T>& );
#endif
public:
	//! explicit four momentum
	/**
	\return the four momentum.
	*/
	const momentum<std::complex<T> >& P() const {return p;};
	//! explicit energy component of the  four vector
	/**
	\return the energy component of the four momentum.
	*/
	std::complex<T> E() const {return p.E();};
	//! explicit x component of the four momentum
	/**
	\return the x component of the four momentum.
	*/
	std::complex<T> X() const {return p.X();};
	//! explicit y component of the  four vector
	/**
	\return the y component of the four momentum.
	*/
	std::complex<T> Y() const {return p.Y();};
	//! explicit z component of the four momentum
	/**
	\return the x component of the four momentum.
	*/
	std::complex<T> Z() const {return p.Z();};

	//! lambda spinor
	/**
	\return the lambda spinor corresponding to the momentum.
	*/
	const lambda<T>& L() const { return La; };
	//! lambda tilde spinor
	/**
	\return the lambda spinor corresponding to the momentum.
	*/
	const lambdat<T>& Lt() const { return Lat;};
	//! slashed matrix
	/**
	\return the (2 dimensional) slashed matrix corresponding to the momentum.
	*/
	momentum_type type() const {return _type;};

	smatrix<T> Sm() const { return smatrix<T>(p);};
	//! Default constructeur
	Cmom() : p(std::complex<T>(0.,0.),std::complex<T>(0.,0.),std::complex<T>(0.,0.),std::complex<T>(0.,0.)), La(std::complex<T>(0.,0.),std::complex<T>(0.,0.)), Lat(std::complex<T>(0.,0.),std::complex<T>(0.,0.)), _type(_mt_unknown) {};
	//! constructeur
	Cmom(const T& E,const T& X,const T& Y,const T& Z,momentum_type type=_mt_unknown) ;
	//! constructeur
	Cmom(const std::complex<T>& E,const std::complex<T>& X,const std::complex<T>& Y,const std::complex<T>& Z,momentum_type type=_mt_unknown);
	//! constructeur
	Cmom(const momentum<T>& pp,momentum_type type=_mt_unknown) ;
	//! constructeur
	Cmom(const momentum<std::complex<T> >& pp,momentum_type type=_mt_unknown);
	//! constructeur
	Cmom(const lambdat<T>& lt,const lambda<T>& l) : p(PfLLt(lt,l)), La(l), Lat(lt), _type(_mt_massless) {};
	//! constructeur
	Cmom(const lambda<T>& l,const lambdat<T>& lt) : p(PfLLt(lt,l)), La(l), Lat(lt), _type(_mt_massless) {};
	//! constructeur
	Cmom(const momentum<std::complex< T> >& pp, const lambda<T>& l,const lambdat<T>& lt,momentum_type type=_mt_unknown) : p(pp), La(l), Lat(lt), _type(type) {};
	Cmom<T> operator+=(const Cmom<T>&);
	Cmom<T> operator-=(const Cmom<T>&);
	Cmom<T> operator*=(const T&);
	Cmom<T> operator*=(const std::complex<T>&);
	Cmom<T> operator/=(const T& c);
	Cmom<T> operator/=(const std::complex<T>& c);
	// type conversion
	template <class U> Cmom(const Cmom<std::complex<U> >& pp) : p(pp.P()) , La(pp.L()), Lat(pp.Lt()) {}
	template <class U> Cmom(const std::complex<U>& E, const std::complex<U>& X, const std::complex<U>& Y, const std::complex<U>& Z,momentum_type type=_mt_unknown) : p(E,X,Y,Z) ,La(std::complex<T>(0.,0.),std::complex<T>(0.,0.)) ,Lat(std::complex<T>(0.,0.),std::complex<T>(0.,0.)) , _type(type) {
		LaLat(p,La,Lat); }
    
    void add_to_M(const Cmom<T>& pin){p.add_to(pin.P());};
    void add_to_U(const Cmom<T>& pin){p.add_to(pin.P());Momentum_to_spinor(pin.P(),La,Lat);};
    void set_to(const Cmom<T>& pin){p.set_to(pin.P());La.set_to(pin.L());Lat.set_to(pin.Lt());};
    void set_to_M(const momentum<std::complex<T> >& pin) {p.set_to(pin);};
    void set_to_U(const momentum<std::complex<T> >& pin) {p.set_to(pin);Momentum_to_spinor(pin,La,Lat);};
    void set_to(const lambdat<T>& lt,const lambda<T>& l){Spinor_to_momentum(lt,l,p);La.set_to(l);Lat.set_to(lt);};
    void set_to(const lambda<T>& l,const lambdat<T>& lt){Spinor_to_momentum(lt,l,p);La.set_to(l);Lat.set_to(lt);};
    void mult_by_M(const std::complex<T>& c){p.mult_by(c);};
    void mult_by_U(const std::complex<T>& c){p.mult_by(c);La.mult_by(sqrt(c));Lat.mult_by(sqrt(c));};
#if SWIG
	%extend {
		std::string __str__() {
		  //make sure you include sstream in the SWIG interface file
		  std::ostringstream oss(std::ostringstream::out);
		  oss << (*$self);
		  return oss.str();
		}
	}
#endif
};


//! Template class for spinorial objects
template <class T> class spinor {
protected:
	std::complex<T> a1;
	std::complex<T> a2;
#ifndef SWIG
	friend spinor<T> operator+<>(const spinor<T>& b1,const spinor<T>& b2);
	friend spinor<T> operator-<>(const spinor<T>& b1,const spinor<T>& b2);
	friend spinor<T> operator-<>(const spinor<T>& b);
	friend spinor<T> operator*<>(const std::complex<T>& c,const spinor<T>& b);
	friend spinor<T> operator*<>(const spinor<T>& b,const std::complex<T>& c);
	friend spinor<T> operator*<>(const T& c,const spinor<T>& b);
	friend spinor<T> operator*<>(const spinor<T>& b,const T& c);
	friend std::complex<T> operator*<>(const spinor<T>& b1,const spinor<T>& b2);
	friend spinor<T> operator*<>(const smatrix<T>& m,const spinor<T>& b);
	friend std::ostream& operator<<<>(std::ostream&, const spinor<T>& );
	friend spinor<T> operator*<>(const spinor<T>& b,const smatrix<T>& m);
#endif
	template <class U> friend class lambda;
	template <class U> friend class lambdat;
public:
    //! Constructor (usually not needed, use the constructors for lambda and lambdat instead)
    /** \sa lambda  lambdat smatrix */
		spinor(const std::complex<T>& A1,const std::complex<T>& A2): a1(A1), a2(A2) {};
    spinor<T> conjugate(){return spinor(-a2,a1);};
		template <class U> spinor(const std::complex<U>& A1,const std::complex<U>& A2) : a1(std::complex<T>(A1)),a2(std::complex<T>(A2)) {}

    ~spinor(){};
    
    void add_to(const spinor<T>& p){a1+=p.a1;a2+=p.a2;};
    void set_to(const spinor<T>& p){a1=p.a1;a2=p.a2;};
    void set_to(const std::complex<T>& A1,const std::complex<T>& A2){a1=A1;a2=A2;};
    void mult_by(const std::complex<T>& p){a1*=p;a2*=p;};
};

template <class T> std::ostream& operator<<(std::ostream&, const smatrix<T>& );
//! Template class for slashed matrices
/**
 The slashed matrices are 2x2 matrices they correspond to the CLat[p].CLa[p] objects in S\@M.
Smatrix objects can be multipied with lambda and lambdat objects. They can't be multipied together since a smatrix object is supposed to have a dotted and an undotted index, whereas the product of two slashed matrices has two dotted or two undotted indices.

  \sa lambda lambdat */
template <class T> class smatrix {
protected:
	std::complex<T> a11;
	std::complex<T> a12;
	std::complex<T> a21;
	std::complex<T> a22;
#ifndef SWIG
	friend spinor<T> operator*<>(const spinor<T>& b,const smatrix<T>& m);
	friend spinor<T> operator*<>(const smatrix<T>& m,const spinor<T>& b);
//	friend smatrix<T> operator*<>(smatrix<T> m1,smatrix<T> m2);
	friend std::ostream& operator<<<>(std::ostream&, const smatrix<T>& );
	//template<typename _Tp, typename _CharT, class _Traits>
	//friend std::basic_ostream<_CharT, _Traits>&
	//operator<<(std::basic_ostream<_CharT, _Traits>& , const spinor<_Tp>& );

#endif
	public:
//! explicit constructor
		smatrix(std::complex<T> A11,std::complex<T> A12,std::complex<T> A21,std::complex<T> A22) : a11(A11), a12(A12), a21(A21), a22(A22){};
//! Constructor with momentum
/** \param p momentum */
    smatrix(const momentum<T>& p);
		smatrix(const momentum<std::complex<T> >& p);

};

#ifndef SWIG
template <class T> spinor<T> la(const momentum<T>& p);
#endif

//const complex<double> II=complex<double>(0,1);
//const complex<qd_real> II=complex<qd_real>(0,1);

template <class T> spinor<T> la(const momentum<T>& p);
template <class T> spinor<T> lat(const momentum<T>& p);
template <class T> spinor<T> la(const momentum<std::complex<T> >& p);
template <class T> spinor<T> lat(const momentum<std::complex<T> >& p);
//template <class T> spinor<T> la(Cmom<T> p);

//template <class T> std::complex<T> spaa(momentum<T> *p1,momentum<T> *p2);
template <class T> std::complex<T> spaa(const momentum<T>&,const momentum<T>& ,const momentum<T>& ,const momentum<T>& );
template <class T> std::complex<T> spaa(const momentum<T>&,const momentum<T>&);
template <class T> std::complex<T> spaa(const momentum<std::complex<T> >&,const momentum<std::complex<T> >);
//template <class T> std::complex<T> spaa(momentum<std::complex<T> > *p1,momentum<std::complex<T> > *p2);
template <class T> std::complex<T> spaa(const momentum<std::complex<T> >&,const momentum<std::complex<T> >& ,const momentum<std::complex<T> >& ,const momentum<std::complex<T> >& );

//template <class T> std::complex<T> spaa(momentum<T>&,momentum<T> &);

//template <class T> std::complex<T> spbb(momentum<T> *p1,momentum<T> *p2);
template <class T> std::complex<T> spbb(const momentum<T>& p1,const momentum<T>& p2);
template <class T> std::complex<T> spbb(const momentum<T>&,const momentum<T>& ,const momentum<T>& ,const momentum<T>& );

//template <class T> std::complex<T> spbb(momentum<std::complex<T> > *p1,momentum<std::complex<T> > *p2);
template <class T> std::complex<T> spbb(const momentum<std::complex<T> >& p1,const momentum<std::complex<T> >& p2);
template <class T> std::complex<T> spbb(const momentum<std::complex<T> >&,const momentum<std::complex<T> >& ,const momentum<std::complex<T> >& ,const momentum<std::complex<T> >& );

//template <class T> std::complex<T> spab(momentum<T> *,momentum<T> *,momentum<T> *);
template <class T> std::complex<T> spab(const momentum<T>&,const momentum<T>&,const momentum<T>&);

//template <class T> std::complex<T> spab(momentum<std::complex<T> > *,momentum<std::complex<T> > *,momentum<std::complex<T> > *);
template <class T> std::complex<T> spab(const momentum<std::complex<T> >&,const momentum<std::complex<T> >&,const momentum<std::complex<T> >&);

template <class T> T s(const momentum<T>&,const momentum<T>&);
template <class T> T s(const momentum<T>&,const momentum<T>&,const momentum<T>&);
template <class T> T s(const momentum<T>&,const momentum<T>&,const momentum<T>&,const momentum<T>&);




template <class T> lambda<T> operator+(const lambda<T>& b1,const lambda<T>& b2);
template <class T> lambda<T> operator-(const lambda<T>& b1,const lambda<T>& b2);
template <class T> lambda<T> operator-(const lambda<T>& b);
template <class T> lambdat<T> operator+(const lambdat<T>& b1,const lambdat<T>& b2);
template <class T> lambdat<T> operator-(const lambdat<T>& b1,const lambdat<T>& b2);
template <class T> lambdat<T> operator-(const lambdat<T>& b);

template <class T> lambda<T> operator*(const std::complex<T>& c,const lambda<T>& b);
template <class T> lambda<T> operator*(const lambda<T>& b,const std::complex<T>& c);
template <class T> lambda<T> operator*(const T& c,const lambda<T>& b);
template <class T> lambda<T> operator*(const lambda<T>& b,const T& c);

template <class T> lambdat<T> operator*(const std::complex<T>& c,const lambdat<T>& b);
template <class T> lambdat<T> operator*(const lambdat<T>& b,const std::complex<T>& c);
template <class T> lambdat<T> operator*(const T& c,const lambdat<T>& b);
template <class T> lambdat<T> operator*(const lambdat<T>& b,const T& c);


template <class T> lambdat<T> operator*(const smatrix<T>& m1,const lambda<T>& l);
template <class T> lambda<T> operator*(const smatrix<T>& m1,const lambdat<T>& l);
template <class T> lambdat<T> operator*(const lambda<T>&,const smatrix<T>&);
template <class T> lambda<T> operator*(const lambdat<T>&,const smatrix<T>&);
template <class T> std::complex<T> operator*(const lambdat<T>& b1,const lambdat<T>& b2);

template <class T> momentum<std::complex<T> > PfLLt(const lambdat<T>& l1,const lambda<T>& l2);


//! Template for spinorial objects of the type "lambda"
template <class T> class lambda : public spinor<T> {
#ifndef SWIG
	friend lambda<T> operator+<>(const lambda<T>& b1,const lambda<T>& b2);
	friend lambda<T> operator-<>(const lambda<T>& b1,const lambda<T>& b2);
	friend lambda<T> operator-<>(const lambda<T>& b);
	friend lambda<T> operator*<>(const std::complex<T>& c,const lambda<T>& b);
	friend lambda<T> operator*<>(const lambda<T>& b,const std::complex<T>& c);
	friend lambda<T> operator*<>(const T& c,const lambda<T>& b);
	friend lambda<T> operator*<>(const lambda<T>& b,const T& c);
    friend lambdat<T> operator*<>(const smatrix<T>& m1,const lambda<T>& l);
	friend lambda<T> operator*<>(const smatrix<T>& m1,const lambdat<T>& l);
	friend lambdat<T> operator*<>(const lambda<T>& l,const smatrix<T>& m1);
	friend lambda<T> operator*<>(const lambdat<T>& l,const smatrix<T>& m1);
	friend momentum<std::complex<T> > PfLLt<>(const lambdat<T>&,const lambda<T>&);
    friend void Spinor_to_momentum<>(const lambdat<T>& l1,const lambda<T>& l2, momentum<std::complex<T> >& ret);
#endif
public:
//	lambda(Cmom<T>& p) : spinor<T>(la(p)) {};
	lambda(const momentum<std::complex<T> >& p) : spinor<T>(la(p)) {};
	lambda(const momentum<T>& p) : spinor<T>(la(p)) {};
	lambda(const spinor<T>& s) : spinor<T>(s) {};
	lambda(const T& A1,const T& A2) : spinor<T>(A1,A2) {};
	lambda(const std::complex<T> A1, const std::complex<T> A2) : spinor<T>(A1,A2) {};
	template <class U> lambda(const U& A1,const U& A2) : spinor<T>(std::complex<T>(A1),std::complex<T>(A2)) {}
	template <class U> lambda(const lambda<U>& L) : spinor<T>(std::complex<T>(L.a1),std::complex<T>(L.a2)) {}

	//	lambda<T> operator-(){return lambda<T>(-this->a1,-this->a2);};

};

//! Template for spinorial objects of the type "lambda tilde"
/** \sa smatrix  */
template <class T> class lambdat : public spinor<T> {
#ifndef SWIG
	friend lambdat<T> operator+<>(const lambdat<T>& b1,const lambdat<T>& b2);
	friend lambdat<T> operator-<>(const lambdat<T>& b1,const lambdat<T>& b2);
	friend lambdat<T> operator-<>(const lambdat<T>& b1);
	friend lambdat<T> operator*<>(const std::complex<T>& c,const lambdat<T>& b);
	friend lambdat<T> operator*<>(const lambdat<T>& b,const std::complex<T>& c);
	friend lambdat<T> operator*<>(const T& c,const lambdat<T>& b);
	friend lambdat<T> operator*<>(const lambdat<T>& b,const T& c);
	friend lambda<T> operator*<>(const smatrix<T>& m1,const lambdat<T>& l);
	friend lambdat<T> operator*<>(const smatrix<T>& m1,const lambda<T>& l);
	friend lambda<T> operator*<>(const lambdat<T>& l,const smatrix<T>& m1);
	friend lambdat<T> operator*<>(const lambda<T>& l,const smatrix<T>& m1);
	friend std::complex<T> operator*<>(const lambdat<T>& b1,const lambdat<T>& b2);
	friend momentum<std::complex<T> > PfLLt<>(const lambdat<T>&,const lambda<T>&);
    friend void Spinor_to_momentum<>(const lambdat<T>& l1,const lambda<T>& l2, momentum<std::complex<T> >& ret);
#endif
public:
	lambdat(const momentum<std::complex<T> >& p) : spinor<T>(lat(p)) {};
	lambdat(const momentum<T>& p) : spinor<T>(lat(p)) {};
	lambdat(const spinor<T> s) : spinor<T>(s) {};
	lambdat(const T A1,const  T A2) : spinor<T>(A1,A2) {};
	lambdat(const std::complex<T> A1,const  std::complex<T> A2) : spinor<T>(A1,A2) {};
	template <class U> lambdat(const U& A1,const U& A2) : spinor<T>(std::complex<T>(A1),std::complex<T>(A2)) {}
	template <class U> lambdat(const lambdat<U>& L) : spinor<T>(L.a1,L.a2) {}

	//lambdat<T> operator-(){return lambdat<T>(-this->a1,-this->a2);};
};

#ifndef SWIG



template <class T> lambdat<T> operator*(const smatrix<T>& m,const lambda<T>& b);

template <class T> lambdat<T> operator*(const lambda<T>& b,const smatrix<T>& m);

template <class T> lambda<T> operator*(const smatrix<T>& m,const lambdat<T>& b);

template <class T> lambda<T> operator*(const lambdat<T>& b,const smatrix<T>& m);

template <class T> std::complex<T> operator*(const lambdat<T>& b1,const lambdat<T>& b2);


#endif

// typedefs

typedef Cmom<double> momC;
typedef momentum<double> mom;
//! std::complex momentum with higher precision (2 x double)
typedef momentum<dd_real> momCHP;
typedef momentum<dd_real> momHP;
//! std::complex momentum with high precision (4 x double)
typedef momentum<qd_real> momCVHP;


#if _DO_INLINE
#include "spinor_inline.h"
#endif
}



#endif
