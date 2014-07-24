/*
 * C_2q1g2l.cpp
 *
 *  Created on: Aug 20, 2008
 *      Author: daniel
 */



#include "C_2q1g2l_eval.h"
#include "tree_amp.h"
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


 
 
template <class T> SeriesC<T> C2q1g2l_qppqmemep_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, p, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qppqmemep L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
	                      c451.push_back(1-1);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(5-1);
                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(4-1);
                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
                            c134.push_back(4-1);    
         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
                            c235.push_back(5-1);    
         vector<int> c245; c245.push_back(2-1); c245.push_back(4-1);
                            c245.push_back(5-1);    
 // #define TimeStamp "Sat 13 Dec 2008 16:03:33 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t4 = square(spa34); 
complex<T> t6 = square(spa13); 
complex<T> t7 = square(spb15); 
complex<T> t9 = spa12*T(2); 
complex<T> d2 = (s23 - s45)*spa12*spa23; d2 = T(1)/d2;
complex<T> d4 = spa12*spa23; d4 = T(1)/d4;
complex<T> d5 = spa45*T(2); d5 = T(1)/d5;
complex<T> d1 = spa23*spa45*t9; d1 = T(1)/d1;
complex<T> d3 = spa23*t9*square(s23 - s45); d3 = T(1)/d3;
complex<T> t2 = -(d2*spa13*spa34*spb15) - d3*spa45*t6*t7 + d1*t4*T(3); 
complex<T> t3 = d2*spa13*spa34*spb15 + d3*spa45*t6*t7; 
complex<T> co1 = d4*spb45*t4; 
complex<T> co2 = d5*spb12*spb23*t4; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*Int(ep,mu,c23,c145) + t3*Int(ep,mu,c45,c123) + co1*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g2l_qpmqmemep_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, m, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpmqmemep L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
	                      c451.push_back(1-1);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(5-1);
                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(4-1);
                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
                            c134.push_back(4-1);    
         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
                            c235.push_back(5-1);    
         vector<int> c245; c245.push_back(2-1); c245.push_back(4-1);
                            c245.push_back(5-1);    
 // #define TimeStamp "Sat 13 Dec 2008 16:03:34 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t4 = square(spb15); 
complex<T> t6 = square(spa34); 
complex<T> t7 = square(spb13); 
complex<T> t10 = spb12*T(2); 
complex<T> d1 = (s12 - s45)*spb12*spb23; d1 = T(1)/d1;
complex<T> d4 = spb12*spb23; d4 = T(1)/d4;
complex<T> d5 = spb45*T(2); d5 = T(1)/d5;
complex<T> d2 = spb23*spb45*t10; d2 = T(1)/d2;
complex<T> d3 = spb23*t10*square(s12 - s45); d3 = T(1)/d3;
complex<T> t2 = -(d1*spa34*spb13*spb15) - d3*spb45*t6*t7; 
complex<T> t3 = d1*spa34*spb13*spb15 + d3*spb45*t6*t7 - d2*t4*T(3); 
complex<T> co1 = -(d4*spa45*t4); 
complex<T> co2 = -(d5*spa12*spa23*t4); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*Int(ep,mu,c12,c345) + t2*Int(ep,mu,c45,c123) + co1*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c1,c2,c3,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g2l_qpqmmemep_SLC
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpqmmemep SLC");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
	                      c451.push_back(1-1);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(5-1);
                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(4-1);
                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
                            c134.push_back(4-1);    
         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
                            c235.push_back(5-1);    
         vector<int> c245; c245.push_back(2-1); c245.push_back(4-1);
                            c245.push_back(5-1);    
 // #define TimeStamp "Sat 13 Dec 2008 16:03:37 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa13 = SPA(1,3);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = S(1,2);
complex<T> t1 = spb23*spb45; 
complex<T> t11 = square(spb25); 
complex<T> t12 = square(spb15); 
complex<T> t14 = -(spb13*spb25) - spb12*spb35; 
complex<T> t16 = square(spa23); 
complex<T> t26 = spa23*spb25; 
complex<T> d1 = spb23*square(s12 - s45); d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spb45*square(spb23); d2 = T(1)/d2;
complex<T> d4 = (s13 - s45)*spb45*square(spb23); d4 = T(1)/d4;
complex<T> d7 = spb45*cube(spb23); d7 = T(1)/d7;
complex<T> d8 = spb13*spb45*square(spb23); d8 = T(1)/d8;
complex<T> d9 = spb45*square(spb23); d9 = T(1)/d9;
complex<T> d10 = spb13*spb45; d10 = T(1)/d10;
complex<T> d11 = spb13*spb23; d11 = T(1)/d11;
complex<T> d12 = cube(spb23); d12 = T(1)/d12;
complex<T> d13 = spb13*square(spb23); d13 = T(1)/d13;
complex<T> d14 = spb45*cube(spb23)*T(2); d14 = T(1)/d14;
complex<T> d15 = spb45*square(spb23)*T(2); d15 = T(1)/d15;
complex<T> d16 = spb13*spb45*T(2); d16 = T(1)/d16;
complex<T> t20 = -(spb13*t11); 
complex<T> t22 = spb15*t14; 
complex<T> t29 = d1*spa34; 
complex<T> d3 = (s13 - s45)*t1; d3 = T(1)/d3;
complex<T> d5 = t1*square(s13 - s45)*T(2); d5 = T(1)/d5;
complex<T> d6 = spb13*t1*T(2); d6 = T(1)/d6;
complex<T> t4 = s12*(d7*t20 - d8*t22); 
complex<T> t5 = d7*s13*t20 + d9*spa13*t22; 
complex<T> t6 = -(spa45*(d11*t12 - d12*t20 + d13*t22)); 
complex<T> t9 = d4*spa23*spb13*t11 + d5*spb13*t11*t16 - d3*spb15*t26*T(2); 
complex<T> t10 = d4*spa23*t20 + d5*t16*t20 - d2*spb12*spb35*t26 - spa23*spb12*spb35*t29 + d3*spb15*t26*T(2) - d6*t12*T(3); 
complex<T> t19 = s12*(d14*s13*t20 + d15*spa13*t22); 
complex<T> t25 = spb12*spb35*(d2*t26 + spa23*t29); 
complex<T> co1 = d10*spa23*t12; 
complex<T> co2 = d16*s12*spa23*t12; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t25*Int(ep,mu,c12,c345) + t9*Int(ep,mu,c13,c245) + t10*Int(ep,mu,c45,c123) + t4*Int(ep,mu,c1,c2,c345) + t5*Int(ep,mu,c1,c3,c245) + co1*Int(ep,mu,c2,c3,c145) + t6*Int(ep,mu,c4,c5,c123) + t19*Int(ep,mu,c3,c1,c2,c45) + co2*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g2l_qpqmpemep_SLC
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpqmpemep SLC");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
	                      c451.push_back(1-1);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(5-1);
                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(4-1);
                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
                            c134.push_back(4-1);    
         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
                            c235.push_back(5-1);    
         vector<int> c245; c245.push_back(2-1); c245.push_back(4-1);
                            c245.push_back(5-1);    
 // #define TimeStamp "Sat 13 Dec 2008 16:03:40 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = spa13*spa45; 
complex<T> t12 = square(spa14); 
complex<T> t13 = square(spa24); 
complex<T> t15 = -(spa14*spa23) + spa12*spa34; 
complex<T> t20 = -(spa12*spa34); 
complex<T> t27 = spa14*spb13; 
complex<T> d1 = (s12 - s45)*spa45*square(spa13); d1 = T(1)/d1;
complex<T> d2 = spa13*square(s12 - s45); d2 = T(1)/d2;
complex<T> d3 = (s23 - s45)*spa45*square(spa13); d3 = T(1)/d3;
complex<T> d7 = spa45*cube(spa13); d7 = T(1)/d7;
complex<T> d8 = spa23*spa45*square(spa13); d8 = T(1)/d8;
complex<T> d9 = spa23*spa45; d9 = T(1)/d9;
complex<T> d10 = spa45*square(spa13); d10 = T(1)/d10;
complex<T> d11 = cube(spa13); d11 = T(1)/d11;
complex<T> d12 = spa13*spa23; d12 = T(1)/d12;
complex<T> d13 = spa23*square(spa13); d13 = T(1)/d13;
complex<T> d14 = spa23*spa45*T(2); d14 = T(1)/d14;
complex<T> d15 = spa45*cube(spa13)*T(2); d15 = T(1)/d15;
complex<T> d16 = spa45*square(spa13)*T(2); d16 = T(1)/d16;
complex<T> t8 = -(spa24*t15); 
complex<T> t19 = spa23*t12; 
complex<T> t28 = d2*spb35; 
complex<T> d4 = (s23 - s45)*t1; d4 = T(1)/d4;
complex<T> d5 = t1*square(s23 - s45)*T(2); d5 = T(1)/d5;
complex<T> d6 = spa23*t1*T(2); d6 = T(1)/d6;
complex<T> t4 = s12*(d8*spa24*t15 + d7*t19); 
complex<T> t5 = d7*s23*t19 + d10*spb23*t8; 
complex<T> t6 = spb45*(d12*t13 + d13*spa24*t15 + d11*t19); 
complex<T> t9 = -(d3*spb13*t19) - d5*t19*square(spb13) + d4*spa24*t27*T(2); 
complex<T> t10 = d3*spb13*t19 + d1*t20*t27 + spa12*spa34*spb13*t28 + d5*t19*square(spb13) - d4*spa24*t27*T(2) + d6*t13*T(3); 
complex<T> t11 = s12*(d15*s23*t19 + d16*spb23*t8); 
complex<T> t25 = d1*spa12*spa34*t27 + spb13*t20*t28; 
complex<T> co1 = -(d9*spb13*t13); 
complex<T> co2 = -(d14*s12*spb13*t13); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t25*Int(ep,mu,c12,c345) + t9*Int(ep,mu,c23,c145) + t10*Int(ep,mu,c45,c123) + t4*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c3,c245) + t5*Int(ep,mu,c2,c3,c145) + t6*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c3,c1,c2,c45) + t11*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g2l_qmqppemep_AX
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qmqppemep AX");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
	                      c451.push_back(1-1);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(5-1);
                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(4-1);
                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
                            c134.push_back(4-1);    
         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
                            c235.push_back(5-1);    
         vector<int> c245; c245.push_back(2-1); c245.push_back(4-1);
                            c245.push_back(5-1);    
 // #define TimeStamp "Sat 13 Dec 2008 16:03:41 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> d1 = square(s12 - s45); d1 = T(1)/d1;
complex<T> co1 = d1*spa14*spb23*spb35; 
complex<T> co2 = -(d1*spa14*spb23*spb35); 
complex<T> co3 = -co1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345) + co3*Int(ep,mu,c45,c123);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g2l_qmqpmemep_AX
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qmqpmemep AX");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
	                      c451.push_back(1-1);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(5-1);
                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(4-1);
                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
                            c134.push_back(4-1);    
         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
                            c235.push_back(5-1);    
         vector<int> c245; c245.push_back(2-1); c245.push_back(4-1);
                            c245.push_back(5-1);    
 // #define TimeStamp "Sat 13 Dec 2008 16:03:41 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> d1 = square(s12 - s45); d1 = T(1)/d1;
complex<T> co1 = -(d1*spa13*spa34*spb25); 
complex<T> co2 = d1*spa13*spa34*spb25; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345) + co2*Int(ep,mu,c45,c123);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g2l_qpqmpemep_AX
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, p, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpqmpemep AX");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
	                      c451.push_back(1-1);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(5-1);
                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(4-1);
                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
                            c134.push_back(4-1);    
         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
                            c235.push_back(5-1);    
         vector<int> c245; c245.push_back(2-1); c245.push_back(4-1);
                            c245.push_back(5-1);    
 // #define TimeStamp "Sat 13 Dec 2008 16:03:42 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> d1 = square(s12 - s45); d1 = T(1)/d1;
complex<T> co1 = d1*spa24*spb13*spb35; 
complex<T> co2 = -(d1*spa24*spb13*spb35); 
complex<T> co3 = -co1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345) + co3*Int(ep,mu,c45,c123);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q1g2l_qpqmmemep_AX
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, m, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpqmmemep AX");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
	                      c451.push_back(1-1);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(5-1);
                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(4-1);
                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
                            c134.push_back(4-1);    
         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
                            c235.push_back(5-1);    
         vector<int> c245; c245.push_back(2-1); c245.push_back(4-1);
                            c245.push_back(5-1);    
 // #define TimeStamp "Sat 13 Dec 2008 16:03:42 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> d1 = square(s12 - s45); d1 = T(1)/d1;
complex<T> co1 = -(d1*spa23*spa34*spb15); 
complex<T> co2 = d1*spa23*spa34*spb15; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345) + co2*Int(ep,mu,c45,c123);  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qppqmemep_L C2q1g2l_7400_L
#define _C_qpmqmemep_L C2q1g2l_7382_L
#define _C_qpqmmemep_SLC C2q1g2l_7352_SLC
#define _C_qpqmpemep_SLC C2q1g2l_7460_SLC
#define _C_qmqppemep_AX C2q1g2l_7465_AX
#define _C_qmqpmemep_AX C2q1g2l_7357_AX
#define _C_qpqmpemep_AX C2q1g2l_7460_AX
#define _C_qpqmmemep_AX C2q1g2l_7352_AX
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qppqmemep_L case 7400 : \
          return &C2q1g2l_7400_L
#define _CASE_qpmqmemep_L case 7382 : \
          return &C2q1g2l_7382_L
#define _CASE_qpqmmemep_SLC case 7352 : \
          return &C2q1g2l_7352_SLC
#define _CASE_qpqmpemep_SLC case 7460 : \
          return &C2q1g2l_7460_SLC
#define _CASE_qmqppemep_AX case 7465 : \
          return &C2q1g2l_7465_AX
#define _CASE_qmqpmemep_AX case 7357 : \
          return &C2q1g2l_7357_AX
#define _CASE_qpqmpemep_AX case 7460 : \
          return &C2q1g2l_7460_AX
#define _CASE_qpqmmemep_AX case 7352 : \
          return &C2q1g2l_7352_AX
 
 
 // *************** function definitions using macros ************* 
 
template <class T> SeriesC<T> _C_qppqmemep_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g2l_qppqmemep_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qpmqmemep_L(
        const eval_param<T>& ep, const T& mu){
          return C2q1g2l_qpmqmemep_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmmemep_SLC(
        const eval_param<T>& ep, const T& mu){
          return C2q1g2l_qpqmmemep_SLC(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmpemep_SLC(
        const eval_param<T>& ep, const T& mu){
          return C2q1g2l_qpqmpemep_SLC(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppemep_AX(
        const eval_param<T>& ep, const T& mu){
          return C2q1g2l_qmqppemep_AX(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmemep_AX(
        const eval_param<T>& ep, const T& mu){
          return C2q1g2l_qmqpmemep_AX(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmpemep_AX(
        const eval_param<T>& ep, const T& mu){
          return C2q1g2l_qpqmpemep_AX(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmmemep_AX(
        const eval_param<T>& ep, const T& mu){
          return C2q1g2l_qpqmmemep_AX(ep,mu);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> SeriesC<T> ( *C2q1g2l_AX_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqppemep_AX;
       _CASE_qmqpmemep_AX;
       _CASE_qpqmpemep_AX;
       _CASE_qpqmmemep_AX;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q1g2l_L_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qppqmemep_L;
       _CASE_qpmqmemep_L;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q1g2l_SLC_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qpqmmemep_SLC;
       _CASE_qpqmpemep_SLC;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template SeriesC<R> ( *C2q1g2l_AX_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q1g2l_AX_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q1g2l_AX_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q1g2l_L_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q1g2l_L_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q1g2l_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q1g2l_SLC_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q1g2l_SLC_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q1g2l_SLC_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);




}
