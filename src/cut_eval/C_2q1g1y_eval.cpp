/*
*C_2q1g1y.cpp
*
* Created on 12/13, 2008
*      Author: Zvi's script
*/
 
#include "amplitudes_cut_eval.h"
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



 
 
template <class T> SeriesC<T> C2q1g1y_qmpqpgap_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmpqpgap L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmpqpgam_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmpqpgam L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Sat 13 Dec 2008 18:25:59 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> s14 = S(1,4);
complex<T> s34 = S(3,4);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa23*T(2); d1 = T(1)/d1;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d1*s14*s34*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c34,c12) + co2*Int(ep,mu,c1,c4,c3,c2));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g1y_qmmqpgap_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmmqpgap L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Sat 13 Dec 2008 18:26:00 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> s14 = S(1,4);
complex<T> s34 = S(3,4);
complex<T> t1 = square(spb34); 
complex<T> d1 = spb12*spb23*T(2); d1 = T(1)/d1;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d1*s14*s34*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c14,c23) + co2*Int(ep,mu,c3,c4,c1,c2));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g1y_qmmqpgam_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmmqpgam L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmgapqpp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, gap, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmgapqpp L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmgapqpm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, gap, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmgapqpm L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Sat 13 Dec 2008 18:26:01 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> s12 = S(1,2);
complex<T> s23 = S(2,3);
complex<T> t1 = square(spa12); 
complex<T> d1 = spa14*spa34*T(2); d1 = T(1)/d1;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d1*s12*s23*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c23,c14) + co2*Int(ep,mu,c1,c2,c3,c4));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g1y_qmgamqpp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, gam, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmgamqpp L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Sat 13 Dec 2008 18:26:01 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> s12 = S(1,2);
complex<T> s23 = S(2,3);
complex<T> t1 = square(spa12); 
complex<T> d1 = spa14*spa34*T(2); d1 = T(1)/d1;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d1*s12*s23*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c23,c14) + co2*Int(ep,mu,c1,c2,c3,c4));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g1y_qmgamqpm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, gam, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmgamqpm L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmpqpgap_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmpqpgap nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmpqpgam_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmpqpgam nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmmqpgap_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmmqpgap nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmmqpgam_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmmqpgam nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmgapqpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, gap, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmgapqpp nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmgapqpm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, gap, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmgapqpm nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmgamqpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, gam, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmgamqpp nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmgamqpm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, gam, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmgamqpm nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmqppgap_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmqppgap L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmqppgam_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmqppgam L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Sat 13 Dec 2008 18:26:04 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> s14 = -(spa14*spb14);
complex<T> s12 = -(spa12*spb12);
complex<T> s24 = S(2,4);
complex<T> s13 = S(1,3);
complex<T> t2 = square(spa14); 
complex<T> t3 = s12*spb12; 
complex<T> t5 = cube(spa14); 
complex<T> t9 = spa34*T(2); 
complex<T> t11 = s12 - s24; 
complex<T> t12 = square(s24); 
complex<T> t14 = -(spa14*T(3)); 
complex<T> d3 = s14; d3 = T(1)/d3;
complex<T> d6 = spa13*spa34; d6 = T(1)/d6;
complex<T> d7 = spa23*spa34; d7 = T(1)/d7;
complex<T> d8 = s13*spa23*spa34; d8 = T(1)/d8;
complex<T> d9 = spa13*spa34*square(spb14); d9 = T(1)/d9;
complex<T> d11 = spa13*spa34*spb14; d11 = T(1)/d11;
complex<T> d12 = s13*spa34; d12 = T(1)/d12;
complex<T> t7 = -(t11*T(3)); 
complex<T> t13 = -(d3*s24*T(2)) + T(3); 
complex<T> t16 = spa24*t2; 
complex<T> t19 = t11*T(3); 
complex<T> t22 = spb12*t2; 
complex<T> d1 = s14*spa12*spa13*t9; d1 = T(1)/d1;
complex<T> d2 = spa12*spa23*t9; d2 = T(1)/d2;
complex<T> d4 = spa13*spb14*t9; d4 = T(1)/d4;
complex<T> d5 = spa12*spa13*t9; d5 = T(1)/d5;
complex<T> d10 = t2*square(spb14); d10 = T(1)/d10;
complex<T> d13 = s13*t9; d13 = T(1)/d13;
complex<T> d14 = spa13*t9*square(spb14); d14 = T(1)/d14;
complex<T> t6 = d10*square(s12) - d3*s24*T(3); 
complex<T> t20 = d4*t13; 
complex<T> t27 = spb12*t16; 
complex<T> t29 = t16*t3; 
complex<T> t4 = d11*t6; 
complex<T> t8 = -t20; 
complex<T> t23 = t3*(d14*t12*t14 + d4*s24*t2*t6); 
complex<T> t26 = t20*t22 + d1*t5*t7 + d2*t16*T(3); 
complex<T> t1 = d9*spb12*t12*t14 + s24*t22*t4; 
complex<T> t21 = d7*t27 + d8*t29 + d9*s24*t14*t3 + t2*t3*t4 - d6*spb12*t5; 
complex<T> t28 = d1*t19*t5 + t22*t8 - d5*t5*T(3); 
complex<T> co1 = -(d12*spb23*t27); 
complex<T> co2 = -(d13*spb23*t29); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t26*Int(ep,mu,c12,c34) + t28*Int(ep,mu,c24,c13) + t21*Int(ep,mu,c1,c2,c34) + co1*Int(ep,mu,c2,c3,c14) + t1*Int(ep,mu,c2,c4,c13) + co2*Int(ep,mu,c1,c2,c3,c4) + t23*Int(ep,mu,c1,c2,c4,c3));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g1y_qmqpmgap_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmqpmgap L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Sat 13 Dec 2008 18:26:06 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb24 = SPB(2,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = S(2,3);
complex<T> s14 = S(1,4);
complex<T> t2 = square(spa13); 
complex<T> t6 = cube(spa13); 
complex<T> t9 = spa34*T(2); 
complex<T> t11 = s12 - s23; 
complex<T> t12 = square(s23); 
complex<T> t14 = s12*spb12; 
complex<T> t15 = spa13*T(3); 
complex<T> d3 = s13; d3 = T(1)/d3;
complex<T> d6 = spa14*spa34; d6 = T(1)/d6;
complex<T> d7 = spa24*spa34; d7 = T(1)/d7;
complex<T> d8 = s14*spa24*spa34; d8 = T(1)/d8;
complex<T> d9 = spa14*spa34*square(spb13); d9 = T(1)/d9;
complex<T> d11 = spa14*spa34*spb13; d11 = T(1)/d11;
complex<T> d12 = s14*spa34; d12 = T(1)/d12;
complex<T> t3 = -t14; 
complex<T> t7 = -(t11*T(3)); 
complex<T> t13 = -(d3*s23*T(2)) + T(3); 
complex<T> t20 = t11*T(3); 
complex<T> t24 = spa23*t2; 
complex<T> t26 = spb12*t2; 
complex<T> d1 = s13*spa12*spa14*t9; d1 = T(1)/d1;
complex<T> d2 = spa12*spa24*t9; d2 = T(1)/d2;
complex<T> d4 = spa14*spb13*t9; d4 = T(1)/d4;
complex<T> d5 = spa12*spa14*t9; d5 = T(1)/d5;
complex<T> d10 = t2*square(spb13); d10 = T(1)/d10;
complex<T> d13 = spa14*t9*square(spb13); d13 = T(1)/d13;
complex<T> d14 = s14*t9; d14 = T(1)/d14;
complex<T> t5 = d10*square(s12) - d3*s23*T(3); 
complex<T> t16 = -t24; 
complex<T> t21 = d4*t13; 
complex<T> t8 = -t21; 
complex<T> t22 = d11*t5; 
complex<T> t23 = d13*t12*t14*t15 + d4*s23*t2*t3*t5; 
complex<T> t27 = t21*t26 + d1*t6*t7 + d5*t6*T(3); 
complex<T> t1 = d9*spb12*t12*t15 - s23*t22*t26; 
complex<T> t25 = d9*s23*t14*t15 + d7*spb12*t16 + t2*t22*t3 + d8*t24*t3 + d6*spb12*t6; 
complex<T> t29 = d1*t20*t6 + t26*t8 - d2*t24*T(3); 
complex<T> co1 = d12*spb12*spb24*t24; 
complex<T> co2 = d14*spb24*t14*t24; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t29*Int(ep,mu,c12,c34) + t27*Int(ep,mu,c23,c14) + t25*Int(ep,mu,c1,c2,c34) + t1*Int(ep,mu,c2,c3,c14) + co1*Int(ep,mu,c2,c4,c13) + t23*Int(ep,mu,c1,c2,c3,c4) + co2*Int(ep,mu,c1,c2,c4,c3));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g1y_qmqpmgam_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmqpmgam L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmqpgapp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, gap, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmqpgapp L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q1g1y_qmqpgapm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, gap, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmqpgapm L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Sat 13 Dec 2008 18:26:09 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb24 = SPB(2,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = S(2,3);
complex<T> s14 = S(1,4);
complex<T> t2 = square(spa13); 
complex<T> t6 = cube(spa13); 
complex<T> t9 = spa34*T(2); 
complex<T> t11 = s12 - s23; 
complex<T> t12 = square(s23); 
complex<T> t14 = s12*spb12; 
complex<T> t15 = spa13*T(3); 
complex<T> d3 = s13; d3 = T(1)/d3;
complex<T> d6 = spa14*spa34; d6 = T(1)/d6;
complex<T> d7 = spa24*spa34; d7 = T(1)/d7;
complex<T> d8 = s14*spa24*spa34; d8 = T(1)/d8;
complex<T> d9 = spa14*spa34*square(spb13); d9 = T(1)/d9;
complex<T> d11 = spa14*spa34*spb13; d11 = T(1)/d11;
complex<T> d12 = s14*spa34; d12 = T(1)/d12;
complex<T> t3 = -t14; 
complex<T> t7 = -(t11*T(3)); 
complex<T> t13 = -(d3*s23*T(2)) + T(3); 
complex<T> t20 = t11*T(3); 
complex<T> t24 = spa23*t2; 
complex<T> t26 = spb12*t2; 
complex<T> d1 = s13*spa12*spa14*t9; d1 = T(1)/d1;
complex<T> d2 = spa12*spa24*t9; d2 = T(1)/d2;
complex<T> d4 = spa14*spb13*t9; d4 = T(1)/d4;
complex<T> d5 = spa12*spa14*t9; d5 = T(1)/d5;
complex<T> d10 = t2*square(spb13); d10 = T(1)/d10;
complex<T> d13 = spa14*t9*square(spb13); d13 = T(1)/d13;
complex<T> d14 = s14*t9; d14 = T(1)/d14;
complex<T> t5 = d10*square(s12) - d3*s23*T(3); 
complex<T> t16 = -t24; 
complex<T> t21 = d4*t13; 
complex<T> t8 = -t21; 
complex<T> t22 = d11*t5; 
complex<T> t23 = d13*t12*t14*t15 + d4*s23*t2*t3*t5; 
complex<T> t27 = t21*t26 + d1*t6*t7 + d5*t6*T(3); 
complex<T> t1 = d9*spb12*t12*t15 - s23*t22*t26; 
complex<T> t25 = d9*s23*t14*t15 + d7*spb12*t16 + t2*t22*t3 + d8*t24*t3 + d6*spb12*t6; 
complex<T> t29 = d1*t20*t6 + t26*t8 - d2*t24*T(3); 
complex<T> co1 = d12*spb12*spb24*t24; 
complex<T> co2 = d14*spb24*t14*t24; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t29*Int(ep,mu,c12,c34) + t27*Int(ep,mu,c23,c14) + t25*Int(ep,mu,c1,c2,c34) + t1*Int(ep,mu,c2,c3,c14) + co1*Int(ep,mu,c2,c4,c13) + t23*Int(ep,mu,c1,c2,c3,c4) + co2*Int(ep,mu,c1,c2,c4,c3));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g1y_qmqpgamp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, gam, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmqpgamp L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Sat 13 Dec 2008 18:26:11 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb24 = SPB(2,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = S(2,3);
complex<T> s14 = S(1,4);
complex<T> t2 = square(spa13); 
complex<T> t6 = cube(spa13); 
complex<T> t9 = spa34*T(2); 
complex<T> t11 = s12 - s23; 
complex<T> t12 = square(s23); 
complex<T> t14 = s12*spb12; 
complex<T> t15 = spa13*T(3); 
complex<T> d3 = s13; d3 = T(1)/d3;
complex<T> d6 = spa14*spa34; d6 = T(1)/d6;
complex<T> d7 = spa24*spa34; d7 = T(1)/d7;
complex<T> d8 = s14*spa24*spa34; d8 = T(1)/d8;
complex<T> d9 = spa14*spa34*square(spb13); d9 = T(1)/d9;
complex<T> d11 = spa14*spa34*spb13; d11 = T(1)/d11;
complex<T> d12 = s14*spa34; d12 = T(1)/d12;
complex<T> t3 = -t14; 
complex<T> t7 = -(t11*T(3)); 
complex<T> t13 = -(d3*s23*T(2)) + T(3); 
complex<T> t20 = t11*T(3); 
complex<T> t24 = spa23*t2; 
complex<T> t26 = spb12*t2; 
complex<T> d1 = s13*spa12*spa14*t9; d1 = T(1)/d1;
complex<T> d2 = spa12*spa24*t9; d2 = T(1)/d2;
complex<T> d4 = spa14*spb13*t9; d4 = T(1)/d4;
complex<T> d5 = spa12*spa14*t9; d5 = T(1)/d5;
complex<T> d10 = t2*square(spb13); d10 = T(1)/d10;
complex<T> d13 = spa14*t9*square(spb13); d13 = T(1)/d13;
complex<T> d14 = s14*t9; d14 = T(1)/d14;
complex<T> t5 = d10*square(s12) - d3*s23*T(3); 
complex<T> t16 = -t24; 
complex<T> t21 = d4*t13; 
complex<T> t8 = -t21; 
complex<T> t22 = d11*t5; 
complex<T> t23 = d13*t12*t14*t15 + d4*s23*t2*t3*t5; 
complex<T> t27 = t21*t26 + d1*t6*t7 + d5*t6*T(3); 
complex<T> t1 = d9*spb12*t12*t15 - s23*t22*t26; 
complex<T> t25 = d9*s23*t14*t15 + d7*spb12*t16 + t2*t22*t3 + d8*t24*t3 + d6*spb12*t6; 
complex<T> t29 = d1*t20*t6 + t26*t8 - d2*t24*T(3); 
complex<T> co1 = d12*spb12*spb24*t24; 
complex<T> co2 = d14*spb24*t14*t24; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t29*Int(ep,mu,c12,c34) + t27*Int(ep,mu,c23,c14) + t25*Int(ep,mu,c1,c2,c34) + t1*Int(ep,mu,c2,c3,c14) + co1*Int(ep,mu,c2,c4,c13) + t23*Int(ep,mu,c1,c2,c3,c4) + co2*Int(ep,mu,c1,c2,c4,c3));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g1y_qmqpgamm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, gam, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qmqpgamm L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
 // *************** table of switch values ************* 
 
#define _C_qmpqpgap_L C2q1g1y_1171_L
#define _C_qmpqpgam_L C2q1g1y_955_L
#define _C_qmmqpgap_L C2q1g1y_1153_L
#define _C_qmmqpgam_L C2q1g1y_937_L
#define _C_qmgapqpp_L C2q1g1y_751_L
#define _C_qmgapqpm_L C2q1g1y_103_L
#define _C_qmgamqpp_L C2q1g1y_745_L
#define _C_qmgamqpm_L C2q1g1y_97_L
#define _C_qmpqpgap_nf C2q1g1y_1171_nf
#define _C_qmpqpgam_nf C2q1g1y_955_nf
#define _C_qmmqpgap_nf C2q1g1y_1153_nf
#define _C_qmmqpgam_nf C2q1g1y_937_nf
#define _C_qmgapqpp_nf C2q1g1y_751_nf
#define _C_qmgapqpm_nf C2q1g1y_103_nf
#define _C_qmgamqpp_nf C2q1g1y_745_nf
#define _C_qmgamqpm_nf C2q1g1y_97_nf
#define _C_qmqppgap_L C2q1g1y_1201_L
#define _C_qmqppgam_L C2q1g1y_985_L
#define _C_qmqpmgap_L C2q1g1y_1093_L
#define _C_qmqpmgam_L C2q1g1y_877_L
#define _C_qmqpgapp_L C2q1g1y_841_L
#define _C_qmqpgapm_L C2q1g1y_193_L
#define _C_qmqpgamp_L C2q1g1y_805_L
#define _C_qmqpgamm_L C2q1g1y_157_L
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmpqpgap_L case 1171 : \
          return &C2q1g1y_1171_L
#define _CASE_qmpqpgam_L case 955 : \
          return &C2q1g1y_955_L
#define _CASE_qmmqpgap_L case 1153 : \
          return &C2q1g1y_1153_L
#define _CASE_qmmqpgam_L case 937 : \
          return &C2q1g1y_937_L
#define _CASE_qmgapqpp_L case 751 : \
          return &C2q1g1y_751_L
#define _CASE_qmgapqpm_L case 103 : \
          return &C2q1g1y_103_L
#define _CASE_qmgamqpp_L case 745 : \
          return &C2q1g1y_745_L
#define _CASE_qmgamqpm_L case 97 : \
          return &C2q1g1y_97_L
#define _CASE_qmpqpgap_nf case 1171 : \
          return &C2q1g1y_1171_nf
#define _CASE_qmpqpgam_nf case 955 : \
          return &C2q1g1y_955_nf
#define _CASE_qmmqpgap_nf case 1153 : \
          return &C2q1g1y_1153_nf
#define _CASE_qmmqpgam_nf case 937 : \
          return &C2q1g1y_937_nf
#define _CASE_qmgapqpp_nf case 751 : \
          return &C2q1g1y_751_nf
#define _CASE_qmgapqpm_nf case 103 : \
          return &C2q1g1y_103_nf
#define _CASE_qmgamqpp_nf case 745 : \
          return &C2q1g1y_745_nf
#define _CASE_qmgamqpm_nf case 97 : \
          return &C2q1g1y_97_nf
#define _CASE_qmqppgap_L case 1201 : \
          return &C2q1g1y_1201_L
#define _CASE_qmqppgam_L case 985 : \
          return &C2q1g1y_985_L
#define _CASE_qmqpmgap_L case 1093 : \
          return &C2q1g1y_1093_L
#define _CASE_qmqpmgam_L case 877 : \
          return &C2q1g1y_877_L
#define _CASE_qmqpgapp_L case 841 : \
          return &C2q1g1y_841_L
#define _CASE_qmqpgapm_L case 193 : \
          return &C2q1g1y_193_L
#define _CASE_qmqpgamp_L case 805 : \
          return &C2q1g1y_805_L
#define _CASE_qmqpgamm_L case 157 : \
          return &C2q1g1y_157_L
 
 
 // *************** function definitions using macros ************* 
 
template <class T> SeriesC<T> _C_qmpqpgap_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmpqpgap_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpgam_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmpqpgam_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqpgap_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmmqpgap_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqpgam_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmmqpgam_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmgapqpp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmgapqpp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmgapqpm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmgapqpm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmgamqpp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmgamqpp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmgamqpm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmgamqpm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpgap_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmpqpgap_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpgam_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmpqpgam_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqpgap_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmmqpgap_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqpgam_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmmqpgam_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmgapqpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmgapqpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmgapqpm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmgapqpm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmgamqpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmgamqpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmgamqpm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmgamqpm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppgap_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmqppgap_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppgam_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmqppgam_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmgap_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmqpmgap_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmgam_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmqpmgam_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpgapp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmqpgapp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpgapm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmqpgapm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpgamp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmqpgamp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpgamm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g1y_qmqpgamm_L(ep,mu);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> SeriesC<T> ( *C2q1g1y_SLC_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmpqpgap_L;
       _CASE_qmpqpgam_L;
       _CASE_qmmqpgap_L;
       _CASE_qmmqpgam_L;
       _CASE_qmgapqpp_L;
       _CASE_qmgapqpm_L;
       _CASE_qmgamqpp_L;
       _CASE_qmgamqpm_L;
       _CASE_qmqppgap_L;
       _CASE_qmqppgam_L;
       _CASE_qmqpmgap_L;
       _CASE_qmqpmgam_L;
       _CASE_qmqpgapp_L;
       _CASE_qmqpgapm_L;
       _CASE_qmqpgamp_L;
       _CASE_qmqpgamm_L;
 
       default: return 0;
        }
 }

template <class T> SeriesC<T> ( *C2q1g1y_L_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmpqpgap_L;
       _CASE_qmpqpgam_L;
       _CASE_qmmqpgap_L;
       _CASE_qmmqpgam_L;
       _CASE_qmgapqpp_L;
       _CASE_qmgapqpm_L;
       _CASE_qmgamqpp_L;
       _CASE_qmgamqpm_L;
       _CASE_qmqppgap_L;
       _CASE_qmqppgam_L;
       _CASE_qmqpmgap_L;
       _CASE_qmqpmgam_L;
       _CASE_qmqpgapp_L;
       _CASE_qmqpgapm_L;
       _CASE_qmqpgamp_L;
       _CASE_qmqpgamm_L;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q1g1y_nf_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmpqpgap_nf;
       _CASE_qmpqpgam_nf;
       _CASE_qmmqpgap_nf;
       _CASE_qmmqpgam_nf;
       _CASE_qmgapqpp_nf;
       _CASE_qmgapqpm_nf;
       _CASE_qmgamqpp_nf;
       _CASE_qmgamqpm_nf;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template SeriesC<R> ( *C2q1g1y_SLC_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q1g1y_SLC_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q1g1y_SLC_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);

template SeriesC<R> ( *C2q1g1y_L_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q1g1y_L_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q1g1y_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q1g1y_nf_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q1g1y_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q1g1y_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);




}
