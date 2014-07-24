#if _DO_INLINE
#define _INLINE inline
#else
#define _INLINE
#endif

#include <cmath>
using std::abs;

// Defined in BH_utilities.h
template <class T> inline T DeltaZero();
template <class T> inline T DeltaZeroSqr();



template <class T> _INLINE momentum<T> momentum<T>::operator+=(const momentum<T>& p1){
	e+=p1.E(); x+=p1.X(); y+=p1.Y(); z+=p1.Z();
	return *this;
}
template <class T> _INLINE momentum<T> operator+(const momentum<T>& p1,const momentum<T>& p2){
  return momentum<T>(p1.e+p2.e,p1.x+p2.x,p1.y+p2.y,p1.z+p2.z);
}
template <class T> _INLINE momentum<T> momentum<T>::operator-=(const momentum<T>& p1){
	e-=p1.E(); x-=p1.X(); y-=p1.Y(); z-=p1.Z();
	return *this;
}
template <class T> _INLINE momentum<T> operator-(const momentum<T>& p){
  return momentum<T>(-p.e,-p.x,-p.y,-p.z);
}
template <class T> _INLINE momentum<T> operator-(const momentum<T>& p1,const momentum<T>& p2){
  return momentum<T>(p1.e-p2.e,p1.x-p2.x,p1.y-p2.y,p1.z-p2.z);
}
template <class T> _INLINE momentum<T> momentum<T>::operator*=(const T& c){
	e*=c; x*=c; y*=c; z*=c;
	return *this;
}
template <class T> _INLINE T operator*(const momentum<T>& p1,const momentum<T>& p2){
  return (p1.e*p2.e - p1.x*p2.x - p1.y*p2.y - p1.z*p2.z);
}
template <class T> _INLINE momentum<T> operator*(const T& c,const momentum<T>& p){
  return momentum<T>(c*p.e, c*p.x ,c*p.y ,c*p.z);
}
template <class T> _INLINE momentum<T> operator*(const momentum<T>& p,const T& c){
  return momentum<T>(c*p.e, c*p.x ,c*p.y ,c*p.z);
}
template <class T> _INLINE momentum<T> operator/(const momentum<T>& p,const T& c){
  T f=T(1)/c;
	return momentum<T>(p.e*f, p.x*f ,p.y*f ,p.z*f);
}

template <class T> _INLINE momentum<T> momentum<T>::operator/=(const T& c){
	T factor=T(1)/c;
	e*=factor; x*=factor; y*=factor; z*=factor;
	return *this;
}

template <class T> _INLINE Cmom<T> operator+(const Cmom<T>& p1,const Cmom<T>& p2){
  return Cmom<T>(p1.P()+p2.P());
}

template <class T> _INLINE Cmom<T> operator-(const Cmom<T>& p){
  return Cmom<T>(-p.P(),p.L(),-p.Lt());
}
template <class T> _INLINE Cmom<T> operator-(const Cmom<T>& p1,const Cmom<T>& p2){
  return Cmom<T>(p1.P()-p2.P(),_mt_unknown);
}

template <class T> _INLINE std::complex<T> operator*(const Cmom<T>& p1,const Cmom<T>& p2){
  return (p1.P()*p2.P());
}

template <class T> _INLINE Cmom<T> operator*(const Cmom<T>& p,const T& c){
	return c*p;
}

template <class T> _INLINE Cmom<T> operator*(const std::complex<T>& c,const Cmom<T>& p){
  if (c==std::complex<T>(0.,0.)) return Cmom<T>(std::complex<T>(0.,0.),std::complex<T>(0.,0.),std::complex<T>(0.,0.),std::complex<T>(0.,0.));
  if (c.imag()==T(0.)) return c.real()*p;
  return Cmom<T>(c*p.P(),sqrt(c)*p.L(),sqrt(c)*p.Lt());
}

template <class T> _INLINE Cmom<T> operator*(const Cmom<T>& p,const std::complex<T>& c){
  return c*p;
}

template <class T> _INLINE Cmom<T> operator*(const momentum<std::complex<T> >& p,const T& c){
  return c*p;
}

template <class T> _INLINE spinor<T> operator+(const spinor<T>& b1,const spinor<T>& b2){
	return spinor<T>(b1.a1+b2.a1,b1.a2+b2.a2);
}
template <class T> _INLINE spinor<T> operator-(const spinor<T>& b1,const spinor<T>& b2){
	return spinor<T>(b1.a1-b2.a1,b1.a2-b2.a2);
}
template <class T> _INLINE spinor<T> operator-(const spinor<T>& b){
	return spinor<T>(-b.a1,-b.a2);
}
template <class T> _INLINE spinor<T> operator*(const std::complex<T>& c,const spinor<T>& b){
	return spinor<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE spinor<T> operator*(const spinor<T>& b,const std::complex<T>& c){
	return spinor<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE spinor<T> operator*(const T& c,const spinor<T>& b){
	return spinor<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE spinor<T> operator*(const spinor<T>& b,const T& c){
	return spinor<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE spinor<T> operator*(const spinor<T>& b,const smatrix<T>& m){
	return spinor<T>(b.a1*m.a11+b.a2*m.a21,b.a1*m.a12+b.a2*m.a22);
}

template <class T> _INLINE spinor<T> operator*(const smatrix<T>& m,const spinor<T>& b){
	return spinor<T>(m.a11*b.a1+m.a12*b.a2,m.a21*b.a1+m.a22*b.a2);
}

template <class T> _INLINE spinor<T> la(const momentum<T>& p){
	return la<T>(momentum<std::complex<T> >(std::complex<T>(p.E(),0.), std::complex<T>(p.X(),0.), std::complex<T>(p.Y(),0.), std::complex<T>(p.Z(),0.)));
}

template <class T> _INLINE spinor<T> la(const Cmom<T>& p){
	return la<T>(p.P());
}

template <class T> _INLINE spinor<T> lat(const momentum<T>& p){
	return lat<T>(momentum<std::complex<T> >(std::complex<T>(p.E(),0.), std::complex<T>(p.X(),0.), std::complex<T>(p.Y(),0.), std::complex<T>(p.Z(),0.)));
}

//! spinor product < * | * > for real momenta (overloaded)
/**
 \param p1 first (real) momentum
 \param p2 second (real) momentum
\return spinor product <p1|p2>
\sa spaa(momentum<std::complex<T> > p1,momentum<std::complex<T> > p2)
*/

template <class T> _INLINE std::complex<T> spaa(const momentum<T>& p1,const momentum<T>& p2){
	return la(p1)*la(p2);
}
//! spinor product < * | * > for std::complex momenta (overloaded)
/**
 \param p1 first (std::complex valued) momentum
 \param p2 second (std::complex valued)  momentum
\return spinor product <p1|p2>
\sa spaa(momentum<T> p1,momentum<T> p2)
*/

template <class T> _INLINE std::complex<T> spaa(const momentum<std::complex<T> >& p1,const momentum<std::complex<T> >& p2){
	return la(p1)*la(p2);
}

//! spinor product [ * | * ] for real momenta (overloaded)
/**
 \param p1 first (real) momentum
 \param p2 second (real) momentum
\return spinor product [p1|p2]
\sa spbb(momentum<std::complex<T> > p1,momentum<std::complex<T> > p2)
*/

template <class T> _INLINE std::complex<T> spbb(const momentum<T>& p1,const momentum<T>& p2){
return lat(p2)*lat(p1);
}
//! spinor product [ * | * ] for std::complex momenta (overloaded)
/**
 \param p1 first (std::complex valued) momentum
 \param p2 second (std::complex valued)  momentum
\return spinor product [p1|p2]
\sa spbb(momentum<T> p1,momentum<T> p2)
*/

template <class T> _INLINE std::complex<T> spbb(const momentum<std::complex<T> >& p1,const momentum<std::complex<T> >& p2){
return lat(p2)*lat(p1);
}



//! spinor product < * | * | * ] for complex momenta (overloaded)
/**
 \param p1 first (complex) momentum
 \param p2 slashed (complex) momentum
 \param p3 last (complex) momentum
\return spinor product <p1|p2|p3]
\sa spaa(momentum<T> p1,momentum<T> p2)
\sa spbb(momentum<T> p1,momentum<T> p2)
\sa spaa(momentum<complex<T> > p1,momentum<complex<T> > p2)
\sa spbb(momentum<complex<T> > p1,momentum<complex<T> > p2)
*/

template <class T> _INLINE std::complex<T> spab(const momentum<T>& p1,const momentum<T>& q,const momentum<T>& p2){
	spinor<T> b=lat(p2)*smatrix<T>(q);
	return b.conjugate()*la(p1);
}
//! spinor product < * | * | * ] for real momenta (overloaded)
/**
 \param p1 first (real) momentum
 \param p2 slashed (real) momentum
 \param p3 last (real) momentum
\return spinor product <p1|p2|p3]
\sa spaa(momentum<T> p1,momentum<T> p2)
\sa spbb(momentum<T> p1,momentum<T> p2)
\sa spaa(momentum<complex<T> > p1,momentum<complex<T> > p2)
\sa spbb(momentum<complex<T> > p1,momentum<complex<T> > p2)
*/
template <class T> _INLINE std::complex<T> spab(const momentum<std::complex<T> >& p1,const momentum<std::complex<T> >& q,const momentum<std::complex<T> >& p2){
	spinor<T> b=lat(p2)*smatrix<T>(q);
	return b.conjugate()*la(p1);
}

//! Invariant s(p1,p2)
/**
\param p1 momentum
\param p2 momentum
\return invariant (p1+p2)^2 which is of type T
*/
template <class T> _INLINE T s(const momentum<T>& p1,const momentum<T>& p2){
	momentum<T> sum=p1+p2;
	return sum*sum;
}
//! Invariant s(p1,p2,p3)
/**
\param p1 momentum
\param p2 momentum
\param p3 momentum
\return invariant (p1+p2+p3)^2 which is of type T
*/
template <class T> _INLINE T s(const momentum<T>& p1,const momentum<T>& p2,const momentum<T>& p3){
	momentum<T> sum=p1+p2+p3;
	return sum*sum;
}
//! Invariant s(p1,p2,p3,p4)
/**
\param p1 momentum
\param p2 momentum
\param p3 momentum
\param p4 momentum
\return invariant (p1+p2+p3+p4)^2 which is of type T
*/

template <class T> _INLINE T s(const momentum<T>& p1,const momentum<T>& p2,const momentum<T>& p3,const momentum<T>& p4){
	momentum<T> sum=p1+p2+p3+p4;
	return sum*sum;
}

template <class T> _INLINE lambda<T> operator+(const lambda<T>& b1,const lambda<T>& b2){
	return lambda<T>(b1.a1+b2.a1,b1.a2+b2.a2);
}
template <class T> _INLINE lambda<T> operator-(const lambda<T>& b1,const lambda<T>& b2){
	return lambda<T>(b1.a1-b2.a1,b1.a2-b2.a2);
}
template <class T> _INLINE lambda<T> operator-(const lambda<T>& b){
	return lambda<T>(-b.a1,-b.a2);
}
template <class T> _INLINE lambdat<T> operator+(const lambdat<T>& b1,const lambdat<T>& b2){
	return lambdat<T>(b1.a1+b2.a1,b1.a2+b2.a2);
}
template <class T> _INLINE lambdat<T> operator-(const lambdat<T>& b1,const lambdat<T>& b2){
	return lambdat<T>(b1.a1-b2.a1,b1.a2-b2.a2);
}
template <class T> _INLINE lambdat<T> operator-(const lambdat<T>& b){
	return lambdat<T>(-b.a1,-b.a2);
}

template <class T> _INLINE lambda<T> operator*(const std::complex<T>& c,const lambda<T>& b){
	return lambda<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE lambda<T> operator*(const lambda<T>& b,const std::complex<T>& c){
	return lambda<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE lambda<T> operator*(const T& c,const lambda<T>& b){
	return lambda<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE lambda<T> operator*(const lambda<T>& b,const T& c){
	return lambda<T>(c*b.a1,c*b.a2);
}

template <class T> _INLINE lambdat<T> operator*(const std::complex<T>& c,const lambdat<T>& b){
	return lambdat<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE lambdat<T> operator*(const lambdat<T>& b,const std::complex<T>& c){
	return lambdat<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE lambdat<T> operator*(const T& c,const lambdat<T>& b){
	return lambdat<T>(c*b.a1,c*b.a2);
}
template <class T> _INLINE lambdat<T> operator*(const lambdat<T>& b,const T& c){
	return lambdat<T>(c*b.a1,c*b.a2);
}

template <class T> _INLINE std::complex<T> operator*(const spinor<T>& b1,const spinor<T>& b2){
return b2.a1*b1.a2-b1.a1*b2.a2;
}

template <class T> _INLINE lambdat<T> operator*(const smatrix<T>& m,const lambda<T>& b){
	return (m*spinor<T>(-b.a1,-b.a2)).conjugate();
}

template <class T> _INLINE lambdat<T> operator*(const lambda<T>& b,const smatrix<T>& m){
	return lambda<T>((m*spinor<T>(-b.a1,-b.a2)).conjugate());
}

template <class T> _INLINE lambda<T> operator*(const smatrix<T>& m,const lambdat<T>& b){
	return ((spinor<T>(b.a1,b.a2))*m).conjugate();
}

template <class T> _INLINE lambda<T> operator*(const lambdat<T>& b,const smatrix<T>& m){
	return ((spinor<T>(b.a1,b.a2))*m).conjugate();
}

template <class T> _INLINE std::complex<T> operator*(const lambdat<T>& b1,const lambdat<T>& b2){
	return ((spinor<T>(b2.a1,b2.a2))*(spinor<T>(b1.a1,b1.a2)));
}

template <class T> _INLINE momentum<std::complex<T> > PfLLt(const lambdat<T>& l1,const lambda<T>& l2){
    T inv2=T(1)/T(2);
	std::complex<T> A11=inv2*l1.a1*l2.a1;
	std::complex<T> A12=inv2*l1.a1*l2.a2;
	std::complex<T> A21=inv2*l1.a2*l2.a1;
	std::complex<T> A22=inv2*l1.a2*l2.a2;
	return momentum<std::complex<T> >(A22+A11,A21+A12,std::complex<T>(0,-1)*(A12-A21),A11-A22);
}

template <class T> _INLINE void Spinor_to_momentum(const lambdat<T>& l1,const lambda<T>& l2, momentum<std::complex<T> >& ret){
    T inv2=T(1)/T(2);
	std::complex<T> A11=inv2*l1.a1*l2.a1;
	std::complex<T> A12=inv2*l1.a1*l2.a2;
	std::complex<T> A21=inv2*l1.a2*l2.a1;
	std::complex<T> A22=inv2*l1.a2*l2.a2;
	ret.set_to(A22+A11,A21+A12,std::complex<T>(0,-1)*(A12-A21),A11-A22);
}

template <class T> _INLINE void LaLat(const momentum<std::complex<T> >& p,lambda<T>& la,lambdat<T>& lat){

	if ((p.plus()*conj(p.plus())).real()<DeltaZero<T>()){
//		std::cout<< "USING DANGEROUS PHASE!"<< std::endl;

		#if _OLD_PHASE_CONVENTION
		if ((p.minus()*conj(p.minus())).real()<DeltaZero<T>()) {
			lat=lambdat<T>((p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
			la =lambda<T> ((p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
		}
		else {
			std::complex<T> sqm=sqrt(p.minus());
			lat=lambdat<T>(std::complex<T>(1.,0.)/sqm*(p.X()+std::complex<T>(0.,1.)*p.Y()),sqm);
			la=lambda<T>(std::complex<T>(1.,0.)/sqm*(p.X()-std::complex<T>(0.,1.)*p.Y()),sqm);
			}
#else
		if ((p.minus()*conj(p.minus())).real()<DeltaZero<T>()) {
			lat=lambdat<T>((p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
			la =lambda<T> ((p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
		}
		else {
			std::complex<T> sqm=sqrt(p.minus());
			lat=lambdat<T>(std::complex<T>(1.,0.)/sqm*(p.X()+std::complex<T>(0.,1.)*p.Y()),sqm);
			la=lambda<T>(std::complex<T>(1.,0.)/sqm*(p.X()-std::complex<T>(0.,1.)*p.Y()),sqm);
			}
#endif
	}
	else {
#if _OLD_PHASE_CONVENTION
		std::complex<T> sqk1p=sqrt(p.plus());std::complex<T> inv=std::complex<T>(1.,0.)/sqk1p;
		lat=lambdat<T>(sqk1p,inv*(p.X()-std::complex<T>(0.,1.)*p.Y()));
		la=lambda<T>(sqk1p,inv*(p.X()+std::complex<T>(0.,1.)*p.Y()));
#else
		T sqk1p=sqrt(abs(p.plus())); T inv=T(1.)/sqk1p;
		lat=lambdat<T>(p.plus()*inv,inv*(p.X()-std::complex<T>(0.,1.)*p.Y()));
		la=lambda<T>(std::complex<T>(sqk1p),sqk1p*(p.X()+std::complex<T>(T(0),T(1))*p.Y())/p.plus());
#endif
	}
}

template <class T> _INLINE void LaLat(const momentum<T>& p,lambda<T>& la,lambdat<T>& lat){
	if (abs(p.plus())<DeltaZeroSqr<T>()){
//		std::cout << "USING DANGEROUS PHASE!" << std::endl;
		if (abs(p.minus())<DeltaZeroSqr<T>()) {
			lat=lambdat<T>((p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
			la =lambda<T> ((p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
		}
		else {
			std::complex<T> sqm=sqrt(std::complex<T>(p.minus(),0));
			lat=lambdat<T>(std::complex<T>(1.,0.)/sqm*(p.X()+std::complex<T>(0.,1.)*p.Y()),sqm);
			la=lambda<T>(std::complex<T>(1.,0.)/sqm*(p.X()-std::complex<T>(0.,1.)*p.Y()),sqm);
			}
	}
	else {
#if _OLD_PHASE_CONVENTION
		std::complex<T> sqk1p=sqrt(std::complex<T>(p.plus(),0));std::complex<T> inv=std::complex<T>(1.,0.)/sqk1p;
		lat=lambdat<T>(sqk1p,inv*(p.X()-std::complex<T>(0.,1.)*p.Y()));
		la=lambda<T>(sqk1p,inv*(p.X()+std::complex<T>(0.,1.)*p.Y()));
#else
		T sqk1p=sqrt(abs(p.plus())); T inv=T(1)/sqk1p;
		lat=lambdat<T>(std::complex<T>(p.plus(),T(0))*inv,std::complex<T>(inv*p.X(),-inv*p.Y()));
		la=lambda<T>(std::complex<T>(sqk1p,T(0)),sqk1p*std::complex<T>(p.X(),p.Y())/p.plus());
#endif
    }
}

template <class T> _INLINE void Momentum_to_spinor(const momentum<std::complex<T> >& p,lambda<T>& la,lambdat<T>& lat){
    
	if ((p.plus()*conj(p.plus())).real()<DeltaZero<T>()){
        if ((p.minus()*conj(p.minus())).real()<DeltaZero<T>()) {
			lat.set_to((p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
			la.set_to((p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
		}
		else {
			std::complex<T> sqm=sqrt(p.minus());
			lat.set_to(std::complex<T>(1.,0.)/sqm*(p.X()+std::complex<T>(0.,1.)*p.Y()),sqm);
			la.set_to(std::complex<T>(1.,0.)/sqm*(p.X()-std::complex<T>(0.,1.)*p.Y()),sqm);
        }
	}
	else {
		T sqk1p=sqrt(abs(p.plus())); T inv=T(1.)/sqk1p;
		lat.set_to(p.plus()*inv,inv*(p.X()-std::complex<T>(0.,1.)*p.Y()));
		la.set_to(std::complex<T>(sqk1p),sqk1p*(p.X()+std::complex<T>(T(0),T(1))*p.Y())/p.plus());
	}
}

template <class T> _INLINE void Momentum_to_spinor(const momentum<T>& p,lambda<T>& la,lambdat<T>& lat){
	if (abs(p.plus())<DeltaZeroSqr<T>()){
		if (abs(p.minus())<DeltaZeroSqr<T>()) {
			lat.set_to((p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
			la.set_to((p.X()-std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()),(p.X()+std::complex<T>(0.,1.)*p.Y())/sqrt(T(2)*p.X()));
		}
		else {
			std::complex<T> sqm=sqrt(std::complex<T>(p.minus(),0));
			lat.set_to(std::complex<T>(1.,0.)/sqm*(p.X()+std::complex<T>(0.,1.)*p.Y()),sqm);
			la.set_to(std::complex<T>(1.,0.)/sqm*(p.X()-std::complex<T>(0.,1.)*p.Y()),sqm);
        }
	}
	else {
		T sqk1p=sqrt(abs(p.plus())); T inv=T(1)/sqk1p;
		lat.set_to(std::complex<T>(p.plus(),T(0))*inv,std::complex<T>(inv*p.X(),-inv*p.Y()));
		la.set_to(std::complex<T>(sqk1p,T(0)),sqk1p*std::complex<T>(p.X(),p.Y())/p.plus());
    }
}

