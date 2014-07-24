/* Series.cc */

/*  David A. Kosower, November 15, 2007  */

/* Implementation of an object that holds several orders of a Laurent series
   in a parameter, taken around 0.  The series is of the form
      a_{-m}/x^m + a_{-m+1}/x^(m-1) +...+a_n x^n + O(x^(n+1))
*/

//#include "header.h"
#define ImplementingSeries 1
#include "BH_typedefs.h"
#include "qd_suppl.h"
#include "math.h"
#include "Series.h"
#if BH_USE_GMP
#include "gmp_r.h"
#endif

#define IsOdd(number) ((number)&1)

using namespace std;
namespace BH {

#if 0
template <class T> Series<T> operator=(const Series<T>& s)
{this->term = s.term;
 this->startOrder = s.startOrder;
 this->endOrder = s.endOrder;
 this->parameter = s.parameter;

 return *this;
}
#endif

template<class T> Series<T> operator+(const Series<T>& s1,const Series<T>& s2)
{Series<T> r(min(s1.leading(),s2.leading()),min(s1.last(),s2.last()));
 int j;

 for (j = s1.startOrder; j < s2.startOrder; j += 1)  r[j] = s1[j];
 for (j = s2.startOrder; j < s1.startOrder; j += 1)  r[j] = s2[j];
 for (j = max(s1.startOrder,s2.startOrder); j <= min(s1.endOrder,s2.endOrder);
      j += 1)
   r[j] = s1[j] + s2[j];

 return r;
}

template <class T> Series<T> operator+(const Series<T>& s1,const T& v2)
{Series<T> r(s1);

 if (r.startOrder <= 0 and r.endOrder >= 0) r[0] += v2;
 return r;
}

template <class T> Series<T> operator+(const T& v1,const Series<T>& s2)
{Series<T> r(s2);

 if (r.startOrder <= 0 and r.endOrder >= 0) r[0] += v1;
 return r;
}

template<class T> Series<T>& Series<T>::operator+=(const Series<T>& s)
{Series<T> r(*this+s);

 this->term = r.term;
 this->startOrder = r.startOrder;
 this->endOrder = r.endOrder;
 return *this;
}

template<class T> Series<T>& Series<T>::operator+=(const T& v)
{if (this->startOrder <= 0 and this->endOrder >= 0) (*this)[0] += v;
 return *this;
}

template<class T> Series<T> operator-(const Series<T>& s1,const Series<T>& s2)
{Series<T> r(min(s1.startOrder,s2.startOrder),min(s1.endOrder,s2.endOrder));
 int j;

 for (j = s1.startOrder; j < s2.startOrder; j += 1)  r[j] = s1[j];
 for (j = s2.startOrder; j < s1.startOrder; j += 1)  r[j] = -s2[j];
 for (j = max(s1.startOrder,s2.startOrder); j <= min(s1.endOrder,s2.endOrder);
      j += 1)
   r[j] = s1[j] - s2[j];

 return r;
}


template <class T> Series<T> operator-(const Series<T>& s1,const T& v2)
{Series<T> r(s1);

 if (r.startOrder <= 0 and r.endOrder >= 0) r[0] -= v2;
 return r;
}

template <class T> Series<T> operator-(const T& v1,const Series<T>& s2)
{Series<T> r(-s2);

 if (r.startOrder <= 0 and r.endOrder >= 0) r[0] += v1;
 return r;
}

template<class T> Series<T> operator-(const Series<T>& s)
{Series<T> r(s);
 int j;

 for (j = r.startOrder; j <= r.endOrder; j += 1)  r[j] = -r[j];

 return r;
}

template<class T> Series<T> Series<T>::operator-=(const Series<T>& s)
{Series<T> r(*this-s);

 this->term = r.term;
 this->startOrder = r.startOrder;
 this->endOrder = r.endOrder;
 return *this;
}

template<class T> Series<T> Series<T>::operator-=(const T& v)
{if (this->startOrder <= 0 and this->endOrder >= 0) (*this)[0] -= v;
 return *this;
}

template<class T> Series<T> operator*(const Series<T>& s1,const Series<T>& s2)
{Series<T> r(s1.startOrder+s2.startOrder,
             min(s1.startOrder+s2.endOrder,s2.startOrder+s1.endOrder));

 for (int j1 = s1.startOrder;  j1 <= s1.endOrder;  j1 += 1)
 for (int j2 = s2.startOrder;  j2 <= s2.endOrder;  j2 += 1)
   {if (j1+j2 > r.endOrder) continue;
   r[j1+j2] += s1[j1]*s2[j2];}

 return r;
}

template <class T> Series<T> operator*(const Series<T>& s1,const T& v2)
{Series<T> r(s1);

 for (int j = r.startOrder;  j <= r.endOrder;  j += 1)
   r[j] *= v2;

 return r;
}

template <class T> Series<T> operator*(const T& v1,const Series<T>& s2)
{Series<T> r(s2);

 for (int j = r.startOrder;  j <= r.endOrder;  j += 1)
   r[j] *= v1;
 return r;
}

template<class T> Series<T> Series<T>::operator*=(const Series<T>& s)
{Series<T> r(*this * s);

 this->term = r.term;
 this->startOrder = r.startOrder;
 this->endOrder = r.endOrder;
 return *this;
}

template<class T> Series<T> Series<T>::operator*=(const T& v)
{ for (int j = this->startOrder;  j <= this->endOrder;  j += 1)
    (*this)[j] *= v;
 return *this;
}

template <class T> Series<T> operator/(const Series<T>& s1,const T& v2)
{Series<T> r(s1);

 for (int j = r.startOrder;  j <= r.endOrder;  j += 1)  r[j] /= v2;
 return r;
}

template<class T> Series<T> Series<T>::operator/=(const T& v)
{ for (int j = this->startOrder;  j <= this->endOrder;  j += 1)
    (*this)[j] /= v;
 return *this;
}

template <class T> Series<T> operator^(const Series<T>& s,unsigned int p)
{Series<T> r(p*s.startOrder,(p-1)*s.startOrder+s.endOrder);
#if 0
 double coeff = 1; // Really want the real base type of T
 int index;
#endif

 switch (p) {
 case 0:
   r.startOrder = r.endOrder = 0;
   r.term.push_back(T(1));
   r.term.push_back(T(1.0));


   break;
 case 1:
   r.term = s.term;
   break;
 case 2:
   // Quadratic terms
   for (int j = s.startOrder;  j <= s.endOrder;  j += 1)
     {if (2*j > r.endOrder) break;
     r[2*j] += s[j]*s[j];}
   // Cross terms
   for (int j1 = s.startOrder;  j1 <= s.endOrder;  j1 += 1)
     for (int j2 = j1+1;  j2 <= s.endOrder;  j2 += 1)
        {if (j1+j2 > r.endOrder) break;
        r[j1+j2] += T(2)*s[j1]*s[j2];}
   break;
 default:
   // Divide and conquer
   if (IsOdd(p))
      r = ((s^((p-1)/2))^2) * s;
   else r = (s^(p/2))^2;
   break;
 }

 return r;
}

template<class T> Series<T> Series<T>::operator^=(unsigned int p)
{Series<T> r(*this ^ p);

 this->term = r.term;
 this->startOrder = r.startOrder;
 this->endOrder = r.endOrder;
 return *this;
}

template <class T> std::ostream& operator<<(std::ostream& file,
                                            const Series<T>& s)
{for (int j = s.startOrder;  j < s.endOrder;  j += 1)
  file << j << ":" << s[j] << " ";
 file << s.endOrder << ":" << s[s.endOrder];
 return file;
}

#if 0
template <class T> std::ostream& operator>>(std::istream& file,
                                            Series<T>& s)
{for (int j = s.startOrder;  j < s.endOrder;  j += 1)
  file << j << ":" << s[j];
 file << s.endOrder << ":" << s[s.endOrder];
 return file;
}
#endif


SeriesC<R> to_double(SeriesC<R> s){
	vector<C> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_double(s[j]));
	}
return SeriesC<R>(s.startOrder,s.endOrder,values);
}

SeriesC<R> to_double(SeriesC<RHP> s){
	vector<C> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_double(s[j]));
	}
return SeriesC<R>(s.startOrder,s.endOrder,values);
}

SeriesC<R> to_double(SeriesC<RVHP> s){
	vector<C> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_double(s[j]));
	}
return SeriesC<R>(s.startOrder,s.endOrder,values);
}

#if BH_USE_GMP
SeriesC<R> to_double(SeriesC<RGMP> s){
	vector<C> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_double(s[j]));
	}
return SeriesC<R>(s.startOrder,s.endOrder,values);
}
#endif


SeriesC<RHP> to_HP(SeriesC<RVHP> s){
	vector<CHP> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_HP(s[j]));
	}
return SeriesC<RHP>(s.startOrder,s.endOrder,values);
}

SeriesC<RHP> to_HP(SeriesC<RHP> s){ return s;}

SeriesC<RHP> to_HP(SeriesC<R> s){
	vector<CHP> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_HP(s[j]));
	}
return SeriesC<RHP>(s.startOrder,s.endOrder,values);
}

SeriesC<RVHP> to_VHP(SeriesC<RVHP> s){ return s;}

SeriesC<RVHP> to_VHP(SeriesC<RHP> s){
	vector<CVHP> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_VHP(s[j]));
	}
return SeriesC<RVHP>(s.startOrder,s.endOrder,values);
}


SeriesC<RVHP> to_VHP(SeriesC<R> s){
	vector<CVHP> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_VHP(s[j]));
	}
return SeriesC<RVHP>(s.startOrder,s.endOrder,values);
}




#if BH_USE_GMP
SeriesC<RHP> to_HP(SeriesC<RGMP> s){
	vector<CHP> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_HP(s[j]));
	}
return SeriesC<RHP>(s.startOrder,s.endOrder,values);
}

SeriesC<RVHP> to_VHP(SeriesC<RGMP> s){
	vector<CVHP> values;
	for  (int j = s.startOrder;  j <= s.endOrder;  j++){
		values.push_back(BH::to_VHP(s[j]));
	}
return SeriesC<RVHP>(s.startOrder,s.endOrder,values);
}

#endif



template <class T> SeriesC<T> operator*(const SeriesC<T>& s1,const T& v2)
{SeriesC<T> r(s1);

 for (int j = r.startOrder;  j <= r.endOrder;  j += 1)
   r[j] *= v2;

 return r;
}

template <class T> SeriesC<T> operator*(const T& v2,const SeriesC<T>& s1)
{SeriesC<T> r(s1);

 for (int j = r.startOrder;  j <= r.endOrder;  j += 1)
   r[j] *= v2;

 return r;
}

template class SeriesC<R>;
template class SeriesC<RHP>;
template class SeriesC<RVHP>;

template SeriesC<R> operator*(const R& v1,const SeriesC<R>& s2);
template SeriesC<R> operator*(const SeriesC<R>& s2,const R& v1);
template SeriesC<RHP> operator*(const RHP& v1,const SeriesC<RHP>& s2);
template SeriesC<RHP> operator*(const SeriesC<RHP>& s2,const RHP& v1);
template SeriesC<RVHP> operator*(const RVHP& v1,const SeriesC<RVHP>& s2);
template SeriesC<RVHP> operator*(const SeriesC<RVHP>& s2,const RVHP& v1);



// g++ still requires explicit instantiation -- we do it for complex series

#define INSTANTIATE(TYPE)  \
template class Series<TYPE>;\
template Series<TYPE> operator+(const Series<TYPE>&,const Series<TYPE>&);\
template Series<TYPE> operator+(const Series<TYPE>&,const TYPE&);\
template Series<TYPE> operator+(const TYPE& v1,const Series<TYPE>& s2);\
template Series<TYPE> operator-(const Series<TYPE>& s1,const Series<TYPE>& s2);\
template Series<TYPE> operator-(const Series<TYPE>& s1,const TYPE& v2);\
template Series<TYPE> operator-(const TYPE& v1,const Series<TYPE>& s2);\
template Series<TYPE> operator-(const Series<TYPE>& s1);\
template Series<TYPE> operator*(const Series<TYPE>& s1,const Series<TYPE>& s2);\
template Series<TYPE> operator*(const Series<TYPE>& s1,const TYPE& v2);\
template Series<TYPE> operator*(const TYPE& v1,const Series<TYPE>& s2);\
template Series<TYPE> operator/(const Series<TYPE>& s,const TYPE& v);\
template Series<TYPE> operator^(const Series<TYPE>& s,unsigned int p);\
template std::ostream& operator<<(std::ostream&, const Series<TYPE>&);

#define INSTANTIATE_STATIC(TYPE)  \
template <> TYPE Series<TYPE>::zero=TYPE(0);\
template <> TYPE Series<TYPE>::infinity=std::numeric_limits<TYPE>::quiet_NaN();

INSTANTIATE(C)
INSTANTIATE(CHP)
INSTANTIATE(CVHP)
INSTANTIATE_STATIC(C)
INSTANTIATE_STATIC(CHP)
INSTANTIATE_STATIC(CVHP)

#if BH_USE_GMP
INSTANTIATE(std::complex<RGMP>)
template <> std::complex<RGMP> Series<std::complex<RGMP> >::zero= (RGMP::set_precision(get_default_precision()),std::complex<RGMP>(0));
template <> std::complex<RGMP> Series<std::complex<RGMP> >::infinity=std::numeric_limits<std::complex<RGMP> >::quiet_NaN();
#endif



}

#if 0
static void dummy() {
 static CSeries dummy1, dummy2;
 static C v;
 dummy1 = dummy1 + v;
 dummy1 = dummy1 + dummy2;
}


#endif

