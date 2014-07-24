/*
*C_2q2Q.cpp
*
* Created on 11/26, 2008
*      Author: Zvi's script
*/
 
#include "C_2q2Q_eval.h"
#include "integrals_ep.h"
 
using namespace std;
 
namespace BH  {
 
 
#define _VERBOSE 0
 
#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}


 
 
template <class T> SeriesC<T> C2q2Q_qmqpQmQp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qm, Qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQmQp nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:52 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c34);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qmqpQpQm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qp, Qm}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQpQm nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:52 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c34);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qmqpQmQp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qm, Qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQmQp L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:53 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> s13 = -(spa13*spb13);
complex<T> s23 = S(2,3);
complex<T> t4 = square(spa13); 
complex<T> t6 = square(s13) + square(s23); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> d2 = s13*spa34; d2 = T(1)/d2;
complex<T> d3 = spa34*square(spb13); d3 = T(1)/d3;
complex<T> d4 = spa12*spa34; d4 = T(1)/d4;
complex<T> d5 = spa12*spa34*square(spb13); d5 = T(1)/d5;
complex<T> d6 = spa34*square(spb13)*T(2); d6 = T(1)/d6;
complex<T> t2 = s23*t6; 
complex<T> t8 = d2*spb12; 
complex<T> t1 = d5*t2 - d4*s23*t4*T(2); 
complex<T> t3 = -(t4*(t8 + d1*T(2))); 
complex<T> co1 = t4*t8; 
complex<T> co2 = -(d3*spb12*t6); 
complex<T> co3 = -(d6*spb12*t2); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t3*Int(ep,mu,c12,c34) + co1*Int(ep,mu,c23,c14) + co2*Int(ep,mu,c1,c2,c34) + t1*Int(ep,mu,c2,c3,c14) + co3*Int(ep,mu,c1,c2,c3,c4));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qmqpQpQm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qp, Qm}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQpQm L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:53 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> s23 = S(2,3);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(2)); 
complex<T> co2 = -(d2*spb12*t1*T(2)); 
complex<T> co3 = -(d2*s23*spb12*t1); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(co1*Int(ep,mu,c12,c34) + co2*Int(ep,mu,c1,c2,c34) + co3*Int(ep,mu,c1,c2,c3,c4));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qmqpQmQp_SLC
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qm, Qp}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQmQp SLC");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:54 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34; d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb12*t1*T(2); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c12,c34) + co2*Int(ep,mu,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qmqpQpQm_SLC
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qp, Qm}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQpQm SLC");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:54 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa34; d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb12*t1*T(2); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c12,c34) + co2*Int(ep,mu,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qpqmQpQm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qp, Qm}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qpqmQpQm nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:54 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> t1 = square(spb13); 
complex<T> d1 = spb12*spb34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c34);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qpqmQmQp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qm, Qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qpqmQmQp nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:55 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> t1 = square(spb14); 
complex<T> d1 = spb12*spb34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c34);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qpqmQpQm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qp, Qm}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qpqmQpQm L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:55 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> s13 = -(spa13*spb13);
complex<T> s23 = S(2,3);
complex<T> t4 = square(spb13); 
complex<T> t6 = square(s13) + square(s23); 
complex<T> d1 = s13*spb34; d1 = T(1)/d1;
complex<T> d2 = spb12*spb34*T(3); d2 = T(1)/d2;
complex<T> d3 = spb34*square(spa13); d3 = T(1)/d3;
complex<T> d4 = spb12*spb34*square(spa13); d4 = T(1)/d4;
complex<T> d5 = spb12*spb34; d5 = T(1)/d5;
complex<T> d6 = spb34*square(spa13)*T(2); d6 = T(1)/d6;
complex<T> t2 = s23*t6; 
complex<T> t8 = d1*spa12; 
complex<T> t1 = d4*t2 - d5*s23*t4*T(2); 
complex<T> t3 = -(t4*(t8 + d2*T(2))); 
complex<T> co1 = t4*t8; 
complex<T> co2 = -(d3*spa12*t6); 
complex<T> co3 = -(d6*spa12*t2); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t3*Int(ep,mu,c12,c34) + co1*Int(ep,mu,c23,c14) + co2*Int(ep,mu,c1,c2,c34) + t1*Int(ep,mu,c2,c3,c14) + co3*Int(ep,mu,c1,c2,c3,c4));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qpqmQmQp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qm, Qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qpqmQmQp L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:56 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> s23 = S(2,3);
complex<T> t1 = square(spb14); 
complex<T> d1 = spb12*spb34*T(3); d1 = T(1)/d1;
complex<T> d2 = spb34; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(2)); 
complex<T> co2 = -(d2*spa12*t1*T(2)); 
complex<T> co3 = -(d2*s23*spa12*t1); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(co1*Int(ep,mu,c12,c34) + co2*Int(ep,mu,c1,c2,c34) + co3*Int(ep,mu,c1,c2,c3,c4));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qpqmQpQm_SLC
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qp, Qm}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qpqmQpQm SLC");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:56 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spb13); 
complex<T> d1 = spb12*spb34; d1 = T(1)/d1;
complex<T> d2 = spb34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spa12*t1*T(2); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c12,c34) + co2*Int(ep,mu,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q_qpqmQmQp_SLC
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qm, Qp}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qpqmQmQp SLC");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Wed 26 Nov 2008 10:06:56 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spb14); 
complex<T> d1 = spb12*spb34; d1 = T(1)/d1;
complex<T> d2 = spb34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spa12*t1*T(2); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c12,c34) + co2*Int(ep,mu,c1,c2,c34));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qmqpQmQp_nf C2q2Q_1237_nf
#define _C_qmqpQpQm_nf C2q2Q_1057_nf
#define _C_qmqpQmQp_L C2q2Q_1237_L
#define _C_qmqpQpQm_L C2q2Q_1057_L
#define _C_qmqpQmQp_SLC C2q2Q_1237_SLC
#define _C_qmqpQpQm_SLC C2q2Q_1057_SLC
#define _C_qpqmQpQm_nf C2q2Q_1052_nf
#define _C_qpqmQmQp_nf C2q2Q_1232_nf
#define _C_qpqmQpQm_L C2q2Q_1052_L
#define _C_qpqmQmQp_L C2q2Q_1232_L
#define _C_qpqmQpQm_SLC C2q2Q_1052_SLC
#define _C_qpqmQmQp_SLC C2q2Q_1232_SLC
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmqpQmQp_nf case 1237 : \
          return &C2q2Q_1237_nf
#define _CASE_qmqpQpQm_nf case 1057 : \
          return &C2q2Q_1057_nf
#define _CASE_qmqpQmQp_L case 1237 : \
          return &C2q2Q_1237_L
#define _CASE_qmqpQpQm_L case 1057 : \
          return &C2q2Q_1057_L
#define _CASE_qmqpQmQp_SLC case 1237 : \
          return &C2q2Q_1237_SLC
#define _CASE_qmqpQpQm_SLC case 1057 : \
          return &C2q2Q_1057_SLC
#define _CASE_qpqmQpQm_nf case 1052 : \
          return &C2q2Q_1052_nf
#define _CASE_qpqmQmQp_nf case 1232 : \
          return &C2q2Q_1232_nf
#define _CASE_qpqmQpQm_L case 1052 : \
          return &C2q2Q_1052_L
#define _CASE_qpqmQmQp_L case 1232 : \
          return &C2q2Q_1232_L
#define _CASE_qpqmQpQm_SLC case 1052 : \
          return &C2q2Q_1052_SLC
#define _CASE_qpqmQmQp_SLC case 1232 : \
          return &C2q2Q_1232_SLC
 
 
 // *************** function definitions using macros ************* 
 
template <class T> SeriesC<T> _C_qmqpQmQp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qmqpQmQp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQpQm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qmqpQpQm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQmQp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qmqpQmQp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQpQm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qmqpQpQm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQmQp_SLC(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qmqpQmQp_SLC(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQpQm_SLC(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qmqpQpQm_SLC(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQpQm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qpqmQpQm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQmQp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qpqmQmQp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQpQm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qpqmQpQm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQmQp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qpqmQmQp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQpQm_SLC(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qpqmQpQm_SLC(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQmQp_SLC(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q_qpqmQmQp_SLC(ep,mu);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> SeriesC<T> ( *C2q2Q_L_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqpQmQp_L;
       _CASE_qmqpQpQm_L;
       _CASE_qpqmQpQm_L;
       _CASE_qpqmQmQp_L;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q2Q_nf_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqpQmQp_nf;
       _CASE_qmqpQpQm_nf;
       _CASE_qpqmQpQm_nf;
       _CASE_qpqmQmQp_nf;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q2Q_SLC_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqpQmQp_SLC;
       _CASE_qmqpQpQm_SLC;
       _CASE_qpqmQpQm_SLC;
       _CASE_qpqmQmQp_SLC;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template SeriesC<R> ( *C2q2Q_L_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q_L_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q2Q_nf_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q2Q_SLC_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q_SLC_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q_SLC_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);




}
