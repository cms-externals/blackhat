/*
*C_2q2Q1y.cpp
*
* Created on 11/17, 2008
*      Author: Zvi's script
*/
 
#include "C_2q2Q1y_eval.h"
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

 
 
template <class T> SeriesC<T> C2q2Q1y_qpQmQpqmgap_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, Qm, Qp, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQmQpqmgap nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:49:50 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa24); 
complex<T> d1 = spa15*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c23,c145);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qpQpQmqmgap_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, Qp, Qm, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQpQmqmgap nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:49:50 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa34); 
complex<T> d1 = spa15*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c23,c145);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qmQpQmqpgap_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, Qp, Qm, qp, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQpQmqpgap nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:49:51 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa15*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c23,c145);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qmQmQpqpgap_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, Qm, Qp, qp, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQmQpqpgap nf");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:49:51 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa12); 
complex<T> d1 = spa15*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c23,c145);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qpQmQpqmgap_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, Qm, Qp, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQmQpqmgap L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:53:05 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = -(spa35*spb35);
complex<T> s13 = -(spa13*spb13);
complex<T> t3 = spa15*T(2); 
complex<T> t4 = spa12*spa34; 
complex<T> t5 = square(s12 - s45); 
complex<T> t13 = square(spa24); 
complex<T> t14 = -(spa14*spb13); 
complex<T> t15 = -(spa13*T(3)); 
complex<T> t16 = spa12*spa45; 
complex<T> t17 = (s12 - s45)*spa15; 
complex<T> t24 = square(spb35); 
complex<T> t25 = square(spa23); 
complex<T> t26 = square(spa25); 
complex<T> t27 = square(spa34); 
complex<T> t28 = square(spa45); 
complex<T> t30 = spa13*spa24; 
complex<T> t31 = -(spa14*spa23) - spa13*spa24; 
complex<T> t33 = s23 - s45; 
complex<T> t34 = -(spa14*spb45); 
complex<T> t35 = square(spb13); 
complex<T> t36 = -(s45*T(3)) + s23*T(4); 
complex<T> t37 = spa12*spb13 - spa24*spb34; 
complex<T> t41 = s35*spa24 + spa12*spa34*spb13*T(2); 
complex<T> t47 = cube(spb35); 
complex<T> t53 = spb15*spb45; 
complex<T> t62 = spa14*spa24; 
complex<T> t66 = cube(spa24); 
complex<T> t76 = -(spb13*T(3)); 
complex<T> t81 = s25*spb25; 
complex<T> t90 = spa25*T(3); 
complex<T> t91 = -(spa23*spa45); 
complex<T> t92 = -(spb15*spb35); 
complex<T> t97 = spa13*spb23; 
complex<T> t106 = s12*s23; 
complex<T> t113 = spa24*spb35; 
complex<T> t125 = spb12*spb15; 
complex<T> t163 = spa45*spb35; 
complex<T> d7 = (s12 - s34)*spa15*spa34; d7 = T(1)/d7;
complex<T> d11 = spa15*spa34*square(s12 - s34); d11 = T(1)/d11;
complex<T> d12 = (s12 - s34)*s35*spa12*spa15; d12 = T(1)/d12;
complex<T> d26 = s12 - s34; d26 = T(1)/d26;
complex<T> d28 = spa12*T(2); d28 = T(1)/d28;
complex<T> d30 = square(s12 - s34); d30 = T(1)/d30;
complex<T> d31 = spa15*spa23*spa45*T(4); d31 = T(1)/d31;
complex<T> d32 = s23*spa15*spa45*T(4); d32 = T(1)/d32;
complex<T> d33 = s23*spa15*T(4); d33 = T(1)/d33;
complex<T> d34 = spa15*spa23*spa45*T(6); d34 = T(1)/d34;
complex<T> d42 = spa13*spa15*spa45; d42 = T(1)/d42;
complex<T> d43 = spa15*spa45*square(spa13); d43 = T(1)/d43;
complex<T> d44 = spa15*square(spa35); d44 = T(1)/d44;
complex<T> d45 = spa15*spa34*spa45; d45 = T(1)/d45;
complex<T> d46 = spa45*cube(spa35); d46 = T(1)/d46;
complex<T> d47 = spa15*spa34*cube(spa35); d47 = T(1)/d47;
complex<T> d48 = spa23*spa45; d48 = T(1)/d48;
complex<T> d49 = s23*spa45; d49 = T(1)/d49;
complex<T> d50 = s23; d50 = T(1)/d50;
complex<T> d51 = spa15*spa45; d51 = T(1)/d51;
complex<T> d52 = spa15; d52 = T(1)/d52;
complex<T> d53 = spa12*spa15*square(spa35); d53 = T(1)/d53;
complex<T> d56 = spa12*spa15*cube(spa35); d56 = T(1)/d56;
complex<T> d58 = s23*spa15; d58 = T(1)/d58;
complex<T> d59 = spa13*spa15; d59 = T(1)/d59;
complex<T> d60 = spa15*spa23; d60 = T(1)/d60;
complex<T> d61 = spa15*square(spa13); d61 = T(1)/d61;
complex<T> d63 = spa12*cube(spa35); d63 = T(1)/d63;
complex<T> d66 = spa34*square(s12 + s15 - s34)*T(2); d66 = T(1)/d66;
complex<T> d71 = s23*T(2); d71 = T(1)/d71;
complex<T> d74 = spa34*spa45*T(2); d74 = T(1)/d74;
complex<T> d75 = spa23*spa34*spa45*T(2); d75 = T(1)/d75;
complex<T> d78 = spa12*cube(spa35)*T(2); d78 = T(1)/d78;
complex<T> t43 = -(spa12*spa34*spb13) + spb35*t91; 
complex<T> t54 = d66*s12; 
complex<T> t69 = -t125; 
complex<T> t73 = t24*t28; 
complex<T> t74 = spa13*t25; 
complex<T> t75 = t26*t27; 
complex<T> t82 = d43*t31; 
complex<T> t88 = d48*spb15*t13 + d49*spa24*spb15*t14 + d50*spa24*t92; 
complex<T> t89 = spa14*t13; 
complex<T> t105 = d31*T(3); 
complex<T> t121 = t30*t90; 
complex<T> t123 = d28*spa23; 
complex<T> t130 = spa25*t24; 
complex<T> t132 = d42*T(3); 
complex<T> t137 = spb13*t62; 
complex<T> t154 = spa14*t25; 
complex<T> d1 = spa15*spa45*t4; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa15*spa45*t4; d2 = T(1)/d2;
complex<T> d3 = spa45*t17*t4; d3 = T(1)/d3;
complex<T> d4 = spa23*spa45*t3*t4; d4 = T(1)/d4;
complex<T> d5 = spa45*t17; d5 = T(1)/d5;
complex<T> d6 = s13*spa45*t17; d6 = T(1)/d6;
complex<T> d8 = t16*t17; d8 = T(1)/d8;
complex<T> d9 = (s12 - s34)*t16; d9 = T(1)/d9;
complex<T> d10 = t16*t5*T(2); d10 = T(1)/d10;
complex<T> d13 = s35*spa12*t17; d13 = T(1)/d13;
complex<T> d14 = t3*square(s12 - s34); d14 = T(1)/d14;
complex<T> d15 = t3*t4*t5; d15 = T(1)/d15;
complex<T> d16 = spa35*t3*t4*square(s12 - s34); d16 = T(1)/d16;
complex<T> d17 = (s12 - s34)*s35*spa15*spa35*t4; d17 = T(1)/d17;
complex<T> d18 = spa35*t3*t4*t5; d18 = T(1)/d18;
complex<T> d19 = s35*spa35*t17*t4; d19 = T(1)/d19;
complex<T> d20 = (s12 - s34)*t16*square(s35); d20 = T(1)/d20;
complex<T> d21 = s35*t16*square(s12 - s34)*T(2); d21 = T(1)/d21;
complex<T> d22 = s35*t16*t5*T(2); d22 = T(1)/d22;
complex<T> d23 = (s12 - s45)*t16*square(s35); d23 = T(1)/d23;
complex<T> d24 = spa34*t3; d24 = T(1)/d24;
complex<T> d25 = t16*T(2); d25 = T(1)/d25;
complex<T> d27 = (s12 - s45)*spa23*spa45*t3*t4; d27 = T(1)/d27;
complex<T> d29 = t3*t4; d29 = T(1)/d29;
complex<T> d35 = s23*spa45*t3; d35 = T(1)/d35;
complex<T> d36 = spa15*spa45*t33; d36 = T(1)/d36;
complex<T> d37 = s13*spa15*spa45*t33; d37 = T(1)/d37;
complex<T> d38 = s23*t3; d38 = T(1)/d38;
complex<T> d39 = s23*t3*t33; d39 = T(1)/d39;
complex<T> d40 = t3*square(t33); d40 = T(1)/d40;
complex<T> d41 = s23*t3*square(t33); d41 = T(1)/d41;
complex<T> d54 = t16*cube(spa35); d54 = T(1)/d54;
complex<T> d55 = spa15*t16; d55 = T(1)/d55;
complex<T> d57 = spa15*t4*cube(spa35); d57 = T(1)/d57;
complex<T> d62 = spa15*t4; d62 = T(1)/d62;
complex<T> d64 = spa13*spa45*t3; d64 = T(1)/d64;
complex<T> d65 = spa34*spa45*t3; d65 = T(1)/d65;
complex<T> d67 = t16*t3; d67 = T(1)/d67;
complex<T> d68 = spa45*t3*square(spa13); d68 = T(1)/d68;
complex<T> d69 = spa12*t3; d69 = T(1)/d69;
complex<T> d70 = spa12*spa23*t3; d70 = T(1)/d70;
complex<T> d72 = t4*T(2); d72 = T(1)/d72;
complex<T> d73 = spa23*t4*T(2); d73 = T(1)/d73;
complex<T> d76 = spa12*t3*square(spa35); d76 = T(1)/d76;
complex<T> d77 = spa12*t3*cube(spa35); d77 = T(1)/d77;
complex<T> t8 = d44*spb12*t121 + s12*t132*t62 + d47*spb12*t28*t74 + d46*spb12*t75 + s12*spa14*t82 - d45*spb12*t89; 
complex<T> t9 = -(d52*t113) + d51*spb23*t13 + d51*spa24*t14 + s23*t132*t62 + s23*spa14*t82; 
complex<T> t40 = -(d65*spb12); 
complex<T> t42 = -(d67*spb34); 
complex<T> t48 = -t74; 
complex<T> t49 = -t75; 
complex<T> t56 = d24*spb25; 
complex<T> t57 = d69*spb34; 
complex<T> t64 = d41*(s45*T(3) - s23*T(4)); 
complex<T> t79 = d29*s45; 
complex<T> t100 = t54*t81; 
complex<T> t116 = -t123; 
complex<T> t120 = t105*t13 + d32*t62*t76 - d33*t113*T(3); 
complex<T> t151 = d71*t137*t53 - d73*spa13*t53*t66 + d72*t53*t89 + d71*s45*spa24*t92; 
complex<T> t158 = d4*t15; 
complex<T> t169 = t130*t15; 
complex<T> t1 = d25*d26*spb35*t13 + d9*spb35*t13 + d14*spa34*t130 + d12*spa24*t169 + d30*spa24*t116*t24 + d11*spb15*t163*t30 + d7*spb15*t62 + d16*t48*t73 + d17*t48*t73 + d20*t47*t75 + d21*t47*t75 + d30*t13*t34*t79 + d2*s35*t89 + d1*t89*T(2) - d26*t13*t56*T(3); 
complex<T> t2 = -(d25*d26*spb35*t13) - d9*spb35*t13 - d14*spa34*t130 + d12*t121*t24 + d13*t121*t24 + d30*spa24*t123*t24 - d15*spa45*t154*t24 + d6*spa14*t31*t35 + d10*t113*t41 + d27*t13*t15*t43 + d20*t47*t49 + d21*t47*t49 + d22*t47*t49 + d23*t47*t49 - d7*spb15*t62 - d8*t37*t62 + t158*t66 + d16*t73*t74 + d17*t73*t74 + d18*t73*t74 + d19*t73*t74 + d5*t62*t76 + d1*t89 - d2*s35*t89 - d3*s35*t89 + d30*spb45*t79*t89 + d11*spa45*t30*t92 + d26*t13*t56*T(3); 
complex<T> t10 = d58*s45*t113 + d53*s45*t121 + d60*spb45*t13 + d58*spa24*spb45*t14 + d61*spa14*spb45*t31 + d63*spb45*t49 + d57*s45*t28*t74 + d62*spb45*t89 + d59*spb45*t62*T(3); 
complex<T> t11 = t120 + d13*spa24*t169 + d15*spa45*t154*t24 - d37*spa14*t31*t35 - d6*spa14*t31*t35 + d39*t113*t36 - d10*t113*t41 + d8*t37*t62 + spa24*spb45*t14*t64 + t158*t66 + d18*t48*t73 + d19*t48*t73 + d22*t47*t75 + d23*t47*t75 + d3*s35*t89 + d40*t16*t92 + d5*t137*T(3) + d27*spa13*t13*t43*T(3) + d36*t137*T(4); 
complex<T> t12 = d40*spb15*spb35*t16 + d37*spa14*t31*t35 - d39*t113*t36 + spb45*t137*t64 + d38*t113*T(3) + d35*t137*T(3) - d36*t137*T(4) + d34*t13*T(13); 
complex<T> t22 = d53*s34*t121 + d56*spb34*t28*t48 + d54*s34*t75 + d55*spa14*spb34*square(spa24); 
complex<T> t23 = d76*s34*s45*t121 + d77*s45*spb34*t28*t48 + d78*s34*spb45*t49 + t34*t57*square(spa24); 
complex<T> t134 = spb15*t100; 
complex<T> t142 = -(d70*spa13*spb34*spb45*t66) + spb45*t57*t89; 
complex<T> t156 = s23*t40*t89 + t40*t66*t97 + d64*t106*t62*T(3); 
complex<T> t160 = t42*(s23*t89 + t66*t97); 
complex<T> t175 = t13*t134; 
complex<T> t135 = t175 + d75*spa13*t66*t69 + d74*t125*t89; 
complex<T> co1 = -t175; 
complex<T> co2 = d68*spa14*t106*t31; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*Int(ep,mu,c12,c345) + t120*Int(ep,mu,c15,c234) + t12*Int(ep,mu,c23,c145) + t1*Int(ep,mu,c34,c125) + t11*Int(ep,mu,c45,c123) + t8*Int(ep,mu,c1,c2,c345) + t88*Int(ep,mu,c1,c5,c234) + t9*Int(ep,mu,c2,c3,c145) + t22*Int(ep,mu,c3,c4,c125) + t10*Int(ep,mu,c4,c5,c123) + t156*Int(ep,mu,c1,c2,c3,c45) + co1*Int(ep,mu,c2,c1,c5,c34) + t160*Int(ep,mu,c2,c3,c4,c15) + co2*Int(ep,mu,c3,c2,c1,c45) + t142*Int(ep,mu,c3,c4,c5,c12) + t151*Int(ep,mu,c4,c5,c1,c23) + t135*Int(ep,mu,c5,c1,c2,c34) + t23*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qpQpQmqmgap_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, Qp, Qm, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQpQmqmgap L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:53:10 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = spa15*T(2); 
complex<T> t10 = square(spa34); 
complex<T> t11 = -(s45*T(3)) + s23*T(4); 
complex<T> t15 = -(spb15*spb45); 
complex<T> t16 = spa14*spb12; 
complex<T> t18 = -(s23*T(4)); 
complex<T> t27 = spa34*spb25; 
complex<T> d1 = spa15*spa23*spa45*T(4); d1 = T(1)/d1;
complex<T> d2 = s23*spa15*spa45*T(4); d2 = T(1)/d2;
complex<T> d3 = s23*spa15*T(4); d3 = T(1)/d3;
complex<T> d4 = spa15*spa23*spa45*T(3); d4 = T(1)/d4;
complex<T> d11 = spa23*spa45; d11 = T(1)/d11;
complex<T> d12 = s23*spa45; d12 = T(1)/d12;
complex<T> d13 = s23; d13 = T(1)/d13;
complex<T> d14 = spa15*spa45; d14 = T(1)/d14;
complex<T> d15 = spa15; d15 = T(1)/d15;
complex<T> d16 = s23*spa15; d16 = T(1)/d16;
complex<T> d17 = spa15*spa23; d17 = T(1)/d17;
complex<T> d20 = s23*T(2); d20 = T(1)/d20;
complex<T> d21 = spa23*T(2); d21 = T(1)/d21;
complex<T> d22 = spa23*spa45*T(2); d22 = T(1)/d22;
complex<T> t5 = -t16; 
complex<T> t17 = -t27; 
complex<T> t21 = d20*spb15; 
complex<T> t23 = spa34*t16; 
complex<T> t30 = d1*T(3); 
complex<T> t34 = spb23*t10; 
complex<T> t35 = d2*T(3); 
complex<T> t39 = d3*T(3); 
complex<T> t41 = spb15*t10; 
complex<T> t44 = spb45*t10; 
complex<T> d5 = s23*spa45*t1; d5 = T(1)/d5;
complex<T> d6 = (s23 - s45)*spa23*spa45*t1; d6 = T(1)/d6;
complex<T> d7 = s23*t1; d7 = T(1)/d7;
complex<T> d8 = s23*(s23 - s45)*t1; d8 = T(1)/d8;
complex<T> d9 = t1*square(s23 - s45); d9 = T(1)/d9;
complex<T> d10 = s23*t1*square(s23 - s45); d10 = T(1)/d10;
complex<T> d18 = spa45*t1; d18 = T(1)/d18;
complex<T> d19 = spa23*t1; d19 = T(1)/d19;
complex<T> t24 = d10*(t18 + s45*T(3)); 
complex<T> t26 = d15*t27 + d14*(t23 + t34); 
complex<T> t31 = d9*spb15; 
complex<T> t33 = d16*s45*t17 + d16*spb45*t23 + d17*t44; 
complex<T> t38 = d12*spb15*t23 + d13*spb15*t27 + d11*t41; 
complex<T> t42 = t10*t30 + t23*t35 + t27*t39; 
complex<T> t45 = d21*t10*t15 + s45*t21*t27 + spa34*spb45*t21*t5; 
complex<T> t3 = d8*t11*t27 - spa13*spa45*spb25*t31 + spa34*spb45*t24*t5 + d4*t10*T(2) + d6*spa34*(spa13*spa45*spb15 + spa23*t16)*T(3) - d5*t23*T(3) - d7*t27*T(3); 
complex<T> t4 = d8*t11*t17 + spb45*t23*t24 + spa13*spa45*spb25*t31 + t23*t35 + t27*t39 - d1*t10*T(3) - d6*spa34*(spa13*spa45*spb15 + spa23*t16)*T(3); 
complex<T> co1 = d18*s12*t34; 
complex<T> co2 = d18*s34*t34; 
complex<T> co3 = d19*s34*t44; 
complex<T> co4 = d22*s12*t41; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t42*Int(ep,mu,c15,c234) + t3*Int(ep,mu,c23,c145) + t4*Int(ep,mu,c45,c123) + t38*Int(ep,mu,c1,c5,c234) + t26*Int(ep,mu,c2,c3,c145) + t33*Int(ep,mu,c4,c5,c123) + co1*Int(ep,mu,c1,c2,c3,c45) + co2*Int(ep,mu,c2,c3,c4,c15) + co3*Int(ep,mu,c3,c4,c5,c12) + t45*Int(ep,mu,c4,c5,c1,c23) + co4*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qmQpQmqpgap_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, Qp, Qm, qp, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQpQmqpgap L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:56:50 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> t4 = spa45*T(2); 
complex<T> t6 = spa12*spa15; 
complex<T> t14 = square(spa13); 
complex<T> t15 = -(spa14*spb24); 
complex<T> t16 = -(spa24*T(3)); 
complex<T> t17 = spa12*spa34; 
complex<T> t20 = (s12 - s34)*spa45; 
complex<T> t29 = square(spb25); 
complex<T> t30 = square(spa23); 
complex<T> t31 = square(spa12); 
complex<T> t32 = square(spa15); 
complex<T> t33 = square(spa35); 
complex<T> t34 = spa13*spa24; 
complex<T> t36 = square(s25); 
complex<T> t37 = -(spa14*spa23) - spa13*spa24; 
complex<T> t38 = s15 - s34; 
complex<T> t40 = s15*T(3); 
complex<T> t41 = s15*T(3) - s23*T(4); 
complex<T> t42 = square(spb24); 
complex<T> t44 = -(spa13*spb12) + spa34*spb24; 
complex<T> t46 = -(s25*spa13) - spa12*spa34*spb24*T(2); 
complex<T> t49 = cube(spb25); 
complex<T> t54 = -(spb24*T(3)); 
complex<T> t55 = spa25*spa45; 
complex<T> t57 = spb15*spb45; 
complex<T> t69 = cube(spa13); 
complex<T> t70 = spa14*spb15; 
complex<T> t71 = -(s15*T(3)) + s23*T(4); 
complex<T> t78 = spa13*spa14; 
complex<T> t99 = s23*s34; 
complex<T> t100 = spa24*spb23; 
complex<T> t115 = spa13*spb25; 
complex<T> t128 = spb34*spb45; 
complex<T> t155 = spa15*spb25; 
complex<T> d3 = (s12 - s34)*spa15*spa34; d3 = T(1)/d3;
complex<T> d11 = spa34*T(2); d11 = T(1)/d11;
complex<T> d12 = square(s12 - s34); d12 = T(1)/d12;
complex<T> d13 = spa15*spa34*T(2); d13 = T(1)/d13;
complex<T> d15 = s12 - s34; d15 = T(1)/d15;
complex<T> d17 = spa12*spa45*square(s12 - s34); d17 = T(1)/d17;
complex<T> d18 = spa15*spa23*spa45*T(4); d18 = T(1)/d18;
complex<T> d21 = (s15 - s23)*spa15*spa45; d21 = T(1)/d21;
complex<T> d22 = s23*spa15*spa45*T(4); d22 = T(1)/d22;
complex<T> d25 = (s15 - s23)*s24*spa15*spa45; d25 = T(1)/d25;
complex<T> d28 = s23*spa45*T(4); d28 = T(1)/d28;
complex<T> d39 = spa15*spa23*spa45*T(6); d39 = T(1)/d39;
complex<T> d46 = spa23*spa45; d46 = T(1)/d46;
complex<T> d47 = spa24*spa45; d47 = T(1)/d47;
complex<T> d48 = spa45*square(spa24); d48 = T(1)/d48;
complex<T> d50 = s23*spa45; d50 = T(1)/d50;
complex<T> d53 = spa15*spa24*spa45; d53 = T(1)/d53;
complex<T> d54 = spa15*spa45*square(spa24); d54 = T(1)/d54;
complex<T> d55 = spa15*spa45; d55 = T(1)/d55;
complex<T> d56 = spa45; d56 = T(1)/d56;
complex<T> d61 = spa15*spa23; d61 = T(1)/d61;
complex<T> d62 = s23*spa15; d62 = T(1)/d62;
complex<T> d63 = s23; d63 = T(1)/d63;
complex<T> d65 = s23*T(2); d65 = T(1)/d65;
complex<T> d68 = spa12*spa35*T(2); d68 = T(1)/d68;
complex<T> t5 = spa34*t38; 
complex<T> t18 = square(t38); 
complex<T> t35 = s12 + t38; 
complex<T> t39 = -t70; 
complex<T> t60 = d68*s34; 
complex<T> t67 = -t115; 
complex<T> t75 = -t128; 
complex<T> t77 = spb12*t36; 
complex<T> t79 = t29*t32; 
complex<T> t81 = spa24*t30; 
complex<T> t87 = d21*spb24; 
complex<T> t94 = spa35*t29; 
complex<T> t105 = spa14*t14; 
complex<T> t106 = t31*t33; 
complex<T> t110 = d18*T(3); 
complex<T> t113 = d54*t37; 
complex<T> t118 = -(d28*T(3)); 
complex<T> t124 = spa13*t15; 
complex<T> t125 = d53*T(3); 
complex<T> t133 = t34*T(3); 
complex<T> t134 = d11*spa23; 
complex<T> t151 = spa14*t30; 
complex<T> t157 = d22*t54; 
complex<T> d1 = spa34*spa45*t6; d1 = T(1)/d1;
complex<T> d4 = spa25*t17*t4*square(s12 - s34); d4 = T(1)/d4;
complex<T> d6 = t4*square(s12 - s34); d6 = T(1)/d6;
complex<T> d10 = t17*t4; d10 = T(1)/d10;
complex<T> d14 = spa12*t4; d14 = T(1)/d14;
complex<T> d16 = spa12*t20; d16 = T(1)/d16;
complex<T> d20 = spa23*spa34*t4*t6; d20 = T(1)/d20;
complex<T> d23 = spa15*spa45*t38; d23 = T(1)/d23;
complex<T> d24 = s23*t4*square(s15 - s23); d24 = T(1)/d24;
complex<T> d26 = s24*spa15*spa45*t38; d26 = T(1)/d26;
complex<T> d29 = s23*(-s15 + s23)*t4; d29 = T(1)/d29;
complex<T> d38 = t4*square(s15 - s23); d38 = T(1)/d38;
complex<T> d40 = s23*spa15*t4; d40 = T(1)/d40;
complex<T> d41 = s23*t4; d41 = T(1)/d41;
complex<T> d64 = spa15*spa34*t4; d64 = T(1)/d64;
complex<T> d66 = spa15*t4*square(spa24); d66 = T(1)/d66;
complex<T> d67 = t4*t6; d67 = T(1)/d67;
complex<T> d69 = t6*T(2); d69 = T(1)/d69;
complex<T> d70 = spa23*t6*T(2); d70 = T(1)/d70;
complex<T> d71 = spa15*spa24*t4; d71 = T(1)/d71;
complex<T> d72 = t17*T(2); d72 = T(1)/d72;
complex<T> d73 = spa23*t17*T(2); d73 = T(1)/d73;
complex<T> d74 = spa34*t4; d74 = T(1)/d74;
complex<T> d76 = spa23*spa34*t4; d76 = T(1)/d76;
complex<T> t3 = square(t35); 
complex<T> t8 = s23*spa14*t113 + d55*t124 + d55*spb23*t14 + d56*t67 + s23*t125*t78; 
complex<T> t10 = -(d67*spb34*(s23*t105 + t100*t69)) + d66*spa14*t37*t99; 
complex<T> t45 = -(d64*spb12); 
complex<T> t51 = -t81; 
complex<T> t59 = d10*s15; 
complex<T> t61 = d14*spb35; 
complex<T> t80 = -t106; 
complex<T> t121 = d38*spb45; 
complex<T> t122 = spb45*(d62*t124 + d61*t14 + d63*t67); 
complex<T> t138 = spb45*t60; 
complex<T> t139 = t57*(d72*t105 - d73*spa24*t69); 
complex<T> t158 = t115*t118 + t110*t14 + t157*t78; 
complex<T> t160 = d20*t16; 
complex<T> t162 = d65*(s15*spb45*t67 + spb24*t57*t78); 
complex<T> t172 = t16*t94; 
complex<T> t174 = t40*t94; 
complex<T> d2 = spa34*t20*t35*t6; d2 = T(1)/d2;
complex<T> d5 = spa25*t17*t20*t35; d5 = T(1)/d5;
complex<T> d7 = spa34*t20*t35; d7 = T(1)/d7;
complex<T> d9 = spa15*spa34*t35*square(s12 - s34)*T(2); d9 = T(1)/d9;
complex<T> d19 = spa45*t35*t5*t6; d19 = T(1)/d19;
complex<T> d27 = spa15*spa45*t5; d27 = T(1)/d27;
complex<T> d30 = spa15*spa34*t18*T(2); d30 = T(1)/d30;
complex<T> d31 = t17*t18*t4; d31 = T(1)/d31;
complex<T> d32 = spa25*t17*t18*t4; d32 = T(1)/d32;
complex<T> d33 = spa12*t35*t5*t55; d33 = T(1)/d33;
complex<T> d34 = spa45*t35*t5; d34 = T(1)/d34;
complex<T> d36 = spa15*spa34*t18*t35*T(2); d36 = T(1)/d36;
complex<T> d37 = spa23*t4*t5*t6; d37 = T(1)/d37;
complex<T> d45 = spa15*spa34*cube(t35); d45 = T(1)/d45;
complex<T> d52 = spa34*cube(t35); d52 = T(1)/d52;
complex<T> d60 = spa15*cube(t35); d60 = T(1)/d60;
complex<T> d78 = spa34*cube(t35)*T(2); d78 = T(1)/d78;
complex<T> t13 = spa34*t121*t155 - d25*spa14*t37*t42 + d24*spa13*spb24*t41*t70 + d29*t67*t71 + d41*t115*T(3) + d40*spb24*t78*T(3) + t78*t87*T(4) + d39*t14*T(13); 
complex<T> t26 = d37*(spa12*spa34*spb24 + spa23*t155); 
complex<T> t107 = t51*t79; 
complex<T> t131 = t45*(s23*t105 + t100*t69); 
complex<T> t132 = t49*t80; 
complex<T> t178 = t138*t14; 
complex<T> d8 = (s12 - s34)*spa15*spa34*t3; d8 = T(1)/d8;
complex<T> d35 = spa15*t3*t5; d35 = T(1)/d35;
complex<T> d42 = spa15*spa34*spa45*t3; d42 = T(1)/d42;
complex<T> d43 = spa34*spa45*t3; d43 = T(1)/d43;
complex<T> d44 = spa34*t3*t55; d44 = T(1)/d44;
complex<T> d49 = spa45*t17*t3; d49 = T(1)/d49;
complex<T> d51 = t17*t3*t55; d51 = T(1)/d51;
complex<T> d57 = spa45*t3*t6; d57 = T(1)/d57;
complex<T> d58 = spa12*t3*t55; d58 = T(1)/d58;
complex<T> d59 = spa45*t3; d59 = T(1)/d59;
complex<T> d75 = spa34*t3*t4; d75 = T(1)/d75;
complex<T> d77 = spa25*spa34*t3*t4; d77 = T(1)/d77;
complex<T> t1 = d1*t105 + d13*d15*spb25*t14 + d3*spb25*t14 + t14*t16*t26 + d12*spa13*t134*t29 - d31*spa15*t151*t29 - d17*spb45*t155*t34 - d19*t105*t36 - d2*t105*t36 - d26*spa14*t37*t42 + d35*t106*t49 + d36*t106*t49 + d8*t106*t49 + d9*t106*t49 + d30*t46*t67 + t160*t69 + d12*t14*t59*t70 + d16*spb45*t78 + d27*t44*t78 + d32*t79*t81 + d33*t79*t81 + d4*t79*t81 + d5*t79*t81 - d6*spa12*t94 + d34*t133*t94 + d7*t133*t94 - d15*t14*t61*T(3) + d23*spb24*t78*T(3); 
complex<T> t2 = d4*t107 + d5*t107 + d8*t132 + d9*t132 - d13*d15*spb25*t14 - d3*spb25*t14 + d7*spa13*t172 - d12*spa13*t134*t29 + d17*spb45*t155*t34 + d2*t105*t36 + d12*t14*t39*t59 - d16*spb45*t78 + d6*spa12*t94 + d1*t105*T(2) + d15*t14*t61*T(3); 
complex<T> t9 = d50*s15*t115 + d50*spb15*t124 + d52*spb15*t132 + d46*spb15*t14 + d43*t174*t34 + d49*t14*t36*t70 + d48*t37*t70 + d51*s15*t79*t81 + d47*spa13*t70*T(3); 
complex<T> t11 = s34*(spa14*t113 + t125*t78) + spb34*(-(d57*t105*t36) + d60*t106*t49 + d58*t79*t81 + d59*t133*t94); 
complex<T> t12 = d32*t107 + d33*t107 + d35*t132 + d36*t132 - spa34*t121*t155 + t158 + d34*spa13*t172 + d31*spa15*t151*t29 + d19*t105*t36 + d24*spb15*t124*t41 + d25*spa14*t37*t42 + d26*spa14*t37*t42 + d30*t115*t46 + t160*t69 + d29*t115*t71 - d27*t44*t78 + d23*t54*t78 + spa24*t14*t26*T(3) - t78*t87*T(4); 
complex<T> t27 = d44*spb12*t107 + d45*s12*t106*t49 + d43*s12*t133*t94 + d42*spa14*t77*square(spa13); 
complex<T> t28 = d77*s15*spb12*t107 + d78*s12*spb15*t132 + d75*s12*t174*t34 - d76*spa24*spb12*spb15*t69 + d74*spb12*t70*square(spa13) + d75*t39*t77*square(spa13); 
complex<T> t149 = d69*t105*t128 + t178 + d70*spa24*t69*t75; 
complex<T> co1 = d71*t78*t99*T(3); 
complex<T> co2 = -t178; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*Int(ep,mu,c12,c345) + t12*Int(ep,mu,c15,c234) + t13*Int(ep,mu,c23,c145) + t1*Int(ep,mu,c34,c125) + t158*Int(ep,mu,c45,c123) + t27*Int(ep,mu,c1,c2,c345) + t9*Int(ep,mu,c1,c5,c234) + t8*Int(ep,mu,c2,c3,c145) + t11*Int(ep,mu,c3,c4,c125) + t122*Int(ep,mu,c4,c5,c123) + t131*Int(ep,mu,c1,c2,c3,c45) + t162*Int(ep,mu,c1,c5,c4,c23) + t10*Int(ep,mu,c2,c3,c4,c15) + t149*Int(ep,mu,c3,c4,c5,c12) + co1*Int(ep,mu,c4,c3,c2,c15) + t139*Int(ep,mu,c4,c5,c1,c23) + t28*Int(ep,mu,c5,c1,c2,c34) + co2*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qmQmQpqpgap_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, Qm, Qp, qp, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQmQpqpgap L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:56:55 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa14 = SPA(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s15 = -(spa15*spb15);
complex<T> t9 = square(spa12); 
complex<T> t10 = s15*T(3) - s23*T(4); 
complex<T> t14 = spa14*spb34; 
complex<T> t15 = spa12*T(3); 
complex<T> t22 = -(s15*T(3)) + s23*T(4); 
complex<T> t25 = spa12*spb35; 
complex<T> t27 = spa15*spa24; 
complex<T> t33 = spa12*spb15; 
complex<T> d1 = spa15*spa23*spa45*T(4); d1 = T(1)/d1;
complex<T> d2 = s23*spa15*spa45*T(4); d2 = T(1)/d2;
complex<T> d3 = s23*spa45*square(s15 - s23)*T(2); d3 = T(1)/d3;
complex<T> d4 = s23*spa45*T(4); d4 = T(1)/d4;
complex<T> d5 = s23*(-s15 + s23)*spa45*T(2); d5 = T(1)/d5;
complex<T> d6 = spa45*square(s15 - s23)*T(2); d6 = T(1)/d6;
complex<T> d7 = (s15 - s23)*spa15*spa23*spa45*T(2); d7 = T(1)/d7;
complex<T> d8 = spa15*spa23*spa45*T(3); d8 = T(1)/d8;
complex<T> d9 = s23*spa15*spa45*T(2); d9 = T(1)/d9;
complex<T> d10 = s23*spa45*T(2); d10 = T(1)/d10;
complex<T> d11 = spa23*spa45; d11 = T(1)/d11;
complex<T> d12 = s23*spa45; d12 = T(1)/d12;
complex<T> d13 = spa15*spa45; d13 = T(1)/d13;
complex<T> d14 = spa45; d14 = T(1)/d14;
complex<T> d15 = spa15*spa23; d15 = T(1)/d15;
complex<T> d16 = s23*spa15; d16 = T(1)/d16;
complex<T> d17 = s23; d17 = T(1)/d17;
complex<T> d18 = spa15*spa45*T(2); d18 = T(1)/d18;
complex<T> d19 = s23*T(2); d19 = T(1)/d19;
complex<T> d20 = spa15*spa23*T(2); d20 = T(1)/d20;
complex<T> d21 = spa23*T(2); d21 = T(1)/d21;
complex<T> d22 = spa23*spa45*T(2); d22 = T(1)/d22;
complex<T> t5 = -t14; 
complex<T> t20 = -t25; 
complex<T> t29 = d1*T(3); 
complex<T> t31 = d2*t14; 
complex<T> t36 = d4*spb35; 
complex<T> t38 = spb15*t9; 
complex<T> t43 = spb23*t9; 
complex<T> t51 = spb45*t9; 
complex<T> t3 = d5*t20*t22 + d6*spb35*spb45*t27 + t15*t31 + d3*t10*t14*t33 + t15*t36 + d7*spa12*(spb45*t27 - spa23*t5)*T(3) - d1*t9*T(3); 
complex<T> t4 = d5*t22*t25 - d6*spb35*spb45*t27 - d7*spb45*t15*t27 + d7*spa23*t15*t5 + d3*t10*t33*t5 + d8*t9*T(2) - d9*spa12*t14*T(3) - d10*t25*T(3); 
complex<T> t24 = d19*spb45*(s15*t25 + t33*t5); 
complex<T> t32 = d12*s15*t20 + d12*t14*t33 + d11*t38; 
complex<T> t37 = d14*t25 + d13*(spa12*t14 + t43); 
complex<T> t42 = d16*spa12*spb45*t14 + d17*spb45*t25 + d15*t51; 
complex<T> t45 = t15*(t31 + t36) + t29*t9; 
complex<T> co1 = d18*s12*t43; 
complex<T> co2 = d18*s34*t43; 
complex<T> co3 = d20*s34*t51; 
complex<T> co4 = -(d21*spb45*t38); 
complex<T> co5 = d22*s12*t38; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t3*Int(ep,mu,c15,c234) + t4*Int(ep,mu,c23,c145) + t45*Int(ep,mu,c45,c123) + t32*Int(ep,mu,c1,c5,c234) + t37*Int(ep,mu,c2,c3,c145) + t42*Int(ep,mu,c4,c5,c123) + co1*Int(ep,mu,c1,c2,c3,c45) + t24*Int(ep,mu,c1,c5,c4,c23) + co2*Int(ep,mu,c2,c3,c4,c15) + co3*Int(ep,mu,c3,c4,c5,c12) + co4*Int(ep,mu,c4,c5,c1,c23) + co5*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qpqmQpQmgap_sl
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qp, Qm, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpqmQpQmgap sl");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:57:03 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> s12 = -(spa12*spb12);
complex<T> t3 = spa34*T(2); 
complex<T> t16 = square(spa24); 
complex<T> t17 = square(spa12); 
complex<T> t18 = square(spa45); 
complex<T> t19 = cube(spb15); 
complex<T> t20 = s12 + s25 - s34; 
complex<T> t21 = -(s25*spb25); 
complex<T> t23 = -(s15*spa24) - spa12*spa34*spb13*T(2); 
complex<T> t27 = square(spb35); 
complex<T> d3 = spa12*spa35*spa45*T(2); d3 = T(1)/d3;
complex<T> d4 = (s12 - s34)*spa15*spa34; d4 = T(1)/d4;
complex<T> d7 = (s12 - s34)*(s12 + s15 - s34)*spa15*spa34; d7 = T(1)/d7;
complex<T> d8 = spa15*square(s12 - s34)*T(2); d8 = T(1)/d8;
complex<T> d9 = (s15 - s34)*spa15*spa34; d9 = T(1)/d9;
complex<T> d10 = (s15 - s34)*(s12 + s15 - s34)*spa15*spa34; d10 = T(1)/d10;
complex<T> d15 = spa25*spa34*spa35; d15 = T(1)/d15;
complex<T> d16 = spa15*spa34*spa45; d16 = T(1)/d16;
complex<T> d17 = spa35*spa45; d17 = T(1)/d17;
complex<T> d19 = spa15*spa34*square(s12 + s15 - s34); d19 = T(1)/d19;
complex<T> d20 = spa34*square(s12 + s15 - s34); d20 = T(1)/d20;
complex<T> d22 = spa15*spa25; d22 = T(1)/d22;
complex<T> d24 = spa15*square(s12 + s15 - s34); d24 = T(1)/d24;
complex<T> t37 = t17*t18; 
complex<T> t38 = spb25*t16; 
complex<T> t39 = spa25*t3; 
complex<T> d2 = spa12*spa15*spa45*t3; d2 = T(1)/d2;
complex<T> d5 = (s12 - s34)*spa25*spa34*square(t20); d5 = T(1)/d5;
complex<T> d12 = (s25 - s34)*spa25*spa34*square(t20); d12 = T(1)/d12;
complex<T> d18 = spa25*spa34*cube(t20); d18 = T(1)/d18;
complex<T> d21 = spa34*cube(t20); d21 = T(1)/d21;
complex<T> d23 = spa25*cube(t20); d23 = T(1)/d23;
complex<T> d25 = t3*square(s12 + s15 - s34); d25 = T(1)/d25;
complex<T> d26 = t3*cube(t20); d26 = T(1)/d26;
complex<T> t13 = -(d17*spb12*t16) + d16*spa14*spb12*t16 - d15*spa23*spb12*t16 + d19*s12*t16*t21 + d18*s12*t19*t37; 
complex<T> t14 = spb34*(d22*t16 + d24*t16*t21 + d23*t19*t37); 
complex<T> t29 = -t37; 
complex<T> t51 = s25*t38; 
complex<T> d1 = spa12*spa35*t39; d1 = T(1)/d1;
complex<T> d6 = t20*t39*square(s12 - s34); d6 = T(1)/d6;
complex<T> d11 = t39*square(s25 - s34); d11 = T(1)/d11;
complex<T> d13 = t20*t39*square(s25 - s34); d13 = T(1)/d13;
complex<T> d14 = spa15*t39; d14 = T(1)/d14;
complex<T> t7 = d4*spa14*spa24*spb15 + d10*t16*t21 + d7*t16*t21 - d11*spa24*spb15*t23 + d8*spa25*spa34*t27 + d12*t19*t37 + d13*t19*t37 + d5*t19*t37 + d6*t19*t37 - d9*t38 + d14*t16*T(3); 
complex<T> t36 = d9*t38 + d10*t51; 
complex<T> t42 = t19*t29; 
complex<T> t6 = d11*spa24*spb15*t23 + (d12 + d13)*t42; 
complex<T> t15 = -(d4*spa14*spa24*spb15) - d8*spa25*spa34*t27 + d5*t42 + d6*t42 + d7*t51 - d3*t16*T(3) + d2*spa14*t16*T(3) - d1*spa23*t16*T(3); 
complex<T> co1 = d20*spb15*t51; 
complex<T> co2 = d21*spb25*t42; 
complex<T> co3 = d25*s12*spb15*t51; 
complex<T> co4 = d26*s12*spb25*t42; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t15*Int(ep,mu,c12,c345) + t36*Int(ep,mu,c15,c234) + t6*Int(ep,mu,c25,c134) + t7*Int(ep,mu,c34,c125) + t13*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + co2*Int(ep,mu,c2,c5,c134) + t14*Int(ep,mu,c3,c4,c125) + co3*Int(ep,mu,c5,c1,c2,c34) + co4*Int(ep,mu,c5,c2,c1,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qpqmQmQpgap_sl
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qm, Qp, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpqmQmQpgap sl");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:57:10 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> s12 = -(spa12*spb12);
complex<T> t3 = spa34*T(2); 
complex<T> t16 = square(spa23); 
complex<T> t17 = square(spa12); 
complex<T> t18 = square(spa35); 
complex<T> t19 = cube(spb15); 
complex<T> t20 = s12 + s25 - s34; 
complex<T> t21 = -(s25*spb25); 
complex<T> t23 = -(s15*spa23) + spa12*spa34*spb14*T(2); 
complex<T> t27 = square(spb45); 
complex<T> d3 = spa12*spa35*spa45*T(2); d3 = T(1)/d3;
complex<T> d4 = (s12 - s34)*spa15*spa34; d4 = T(1)/d4;
complex<T> d7 = (s12 - s34)*(s12 + s15 - s34)*spa15*spa34; d7 = T(1)/d7;
complex<T> d8 = spa15*square(s12 - s34)*T(2); d8 = T(1)/d8;
complex<T> d9 = (s15 - s34)*spa15*spa34; d9 = T(1)/d9;
complex<T> d10 = (s15 - s34)*(s12 + s15 - s34)*spa15*spa34; d10 = T(1)/d10;
complex<T> d15 = spa15*spa34*spa35; d15 = T(1)/d15;
complex<T> d16 = spa25*spa34*spa45; d16 = T(1)/d16;
complex<T> d17 = spa35*spa45; d17 = T(1)/d17;
complex<T> d19 = spa15*spa34*square(s12 + s15 - s34); d19 = T(1)/d19;
complex<T> d20 = spa34*square(s12 + s15 - s34); d20 = T(1)/d20;
complex<T> d22 = spa15*spa25; d22 = T(1)/d22;
complex<T> d24 = spa15*square(s12 + s15 - s34); d24 = T(1)/d24;
complex<T> t37 = t17*t18; 
complex<T> t38 = spb25*t16; 
complex<T> t39 = spa25*t3; 
complex<T> d1 = spa12*spa15*spa35*t3; d1 = T(1)/d1;
complex<T> d5 = (s12 - s34)*spa25*spa34*square(t20); d5 = T(1)/d5;
complex<T> d12 = (s25 - s34)*spa25*spa34*square(t20); d12 = T(1)/d12;
complex<T> d18 = spa25*spa34*cube(t20); d18 = T(1)/d18;
complex<T> d21 = spa34*cube(t20); d21 = T(1)/d21;
complex<T> d23 = spa25*cube(t20); d23 = T(1)/d23;
complex<T> d25 = t3*square(s12 + s15 - s34); d25 = T(1)/d25;
complex<T> d26 = t3*cube(t20); d26 = T(1)/d26;
complex<T> t13 = d17*spb12*t16 + d15*spa13*spb12*t16 - d16*spa24*spb12*t16 + d19*s12*t16*t21 + d18*s12*t19*t37; 
complex<T> t14 = spb34*(d22*t16 + d24*t16*t21 + d23*t19*t37); 
complex<T> t29 = -t37; 
complex<T> t51 = s25*t38; 
complex<T> d2 = spa12*spa45*t39; d2 = T(1)/d2;
complex<T> d6 = t20*t39*square(s12 - s34); d6 = T(1)/d6;
complex<T> d11 = t39*square(s25 - s34); d11 = T(1)/d11;
complex<T> d13 = t20*t39*square(s25 - s34); d13 = T(1)/d13;
complex<T> d14 = spa15*t39; d14 = T(1)/d14;
complex<T> t7 = d4*spa13*spa23*spb15 + d10*t16*t21 + d7*t16*t21 - d11*spa23*spb15*t23 + d8*spa25*spa34*t27 + d12*t19*t37 + d13*t19*t37 + d5*t19*t37 + d6*t19*t37 - d9*t38 + d14*t16*T(3); 
complex<T> t36 = d9*t38 + d10*t51; 
complex<T> t42 = t19*t29; 
complex<T> t6 = d11*spa23*spb15*t23 + (d12 + d13)*t42; 
complex<T> t15 = -(d4*spa13*spa23*spb15) - d8*spa25*spa34*t27 + d5*t42 + d6*t42 + d7*t51 + d3*t16*T(3) + d1*spa13*t16*T(3) - d2*spa24*t16*T(3); 
complex<T> co1 = d20*spb15*t51; 
complex<T> co2 = d21*spb25*t42; 
complex<T> co3 = d25*s12*spb15*t51; 
complex<T> co4 = d26*s12*spb25*t42; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t15*Int(ep,mu,c12,c345) + t36*Int(ep,mu,c15,c234) + t6*Int(ep,mu,c25,c134) + t7*Int(ep,mu,c34,c125) + t13*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + co2*Int(ep,mu,c2,c5,c134) + t14*Int(ep,mu,c3,c4,c125) + co3*Int(ep,mu,c5,c1,c2,c34) + co4*Int(ep,mu,c5,c2,c1,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qmqpQpQmgap_sl
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qp, Qm, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmqpQpQmgap sl");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:57:14 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa45 = SPA(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> t2 = spa34*T(2); 
complex<T> t8 = spa15*spa34; 
complex<T> t14 = square(spa12); 
complex<T> t15 = square(spa14); 
complex<T> t16 = square(spa45); 
complex<T> t17 = cube(spb25); 
complex<T> t18 = s12 + s15 - s34; 
complex<T> t22 = -(s25*spa14) + spa12*spa34*spb23*T(2); 
complex<T> t23 = square(spb35); 
complex<T> t24 = s12*spb15; 
complex<T> d3 = spa12*spa35*spa45*T(2); d3 = T(1)/d3;
complex<T> d6 = (s12 - s34)*spa25; d6 = T(1)/d6;
complex<T> d7 = spa25*square(s12 - s34)*T(2); d7 = T(1)/d7;
complex<T> d12 = spa25*spa34*spa35; d12 = T(1)/d12;
complex<T> d14 = spa35*spa45; d14 = T(1)/d14;
complex<T> d15 = (s12 + s25 - s34)*spa25*spa34; d15 = T(1)/d15;
complex<T> d18 = (s12 + s25 - s34)*spa34; d18 = T(1)/d18;
complex<T> d19 = spa15*spa25; d19 = T(1)/d19;
complex<T> d20 = (s12 + s25 - s34)*spa25; d20 = T(1)/d20;
complex<T> t1 = square(t18); 
complex<T> t27 = spa15*t2; 
complex<T> t33 = t14*t16; 
complex<T> t45 = spb15*t15; 
complex<T> t46 = t15*t24; 
complex<T> d1 = spa12*spa25*spa35*t2; d1 = T(1)/d1;
complex<T> d13 = spa45*t8; d13 = T(1)/d13;
complex<T> d16 = t8*cube(t18); d16 = T(1)/d16;
complex<T> d17 = spa34*cube(t18); d17 = T(1)/d17;
complex<T> d21 = spa15*cube(t18); d21 = T(1)/d21;
complex<T> d22 = t2*cube(t18); d22 = T(1)/d22;
complex<T> d23 = (s12 + s25 - s34)*t2; d23 = T(1)/d23;
complex<T> t11 = -(d14*spb12*t15) - d12*spa23*spb12*t15 + d16*s12*t17*t33 + d15*t46 + d13*spb12*cube(spa14); 
complex<T> t12 = spb34*(d19*t15 + d21*t17*t33 + d20*t45); 
complex<T> t26 = -t33; 
complex<T> d2 = spa12*spa45*t27; d2 = T(1)/d2;
complex<T> d4 = (s12 - s34)*t1*t8; d4 = T(1)/d4;
complex<T> d5 = t18*t27*square(s12 - s34); d5 = T(1)/d5;
complex<T> d8 = t27*square(s15 - s34); d8 = T(1)/d8;
complex<T> d9 = (s15 - s34)*t1*t8; d9 = T(1)/d9;
complex<T> d10 = t18*t27*square(s15 - s34); d10 = T(1)/d10;
complex<T> d11 = spa25*t27; d11 = T(1)/d11;
complex<T> t6 = -(d6*spa14*spb35) - d8*spa14*spb25*t22 + d10*t17*t33 + d4*t17*t33 + d5*t17*t33 + d9*t17*t33 + d7*t23*t8 + d11*t15*T(3); 
complex<T> t37 = t17*t26; 
complex<T> t5 = d8*spa14*spb25*t22 + (d10 + d9)*t37; 
complex<T> t13 = d6*spa14*spb35 - d7*spa15*spa34*t23 + d4*t37 + d5*t37 - d3*t15*T(3) - d1*spa23*t15*T(3) + d2*cube(spa14)*T(3); 
complex<T> co1 = d17*spb15*t37; 
complex<T> co2 = -(d18*spb25*t45); 
complex<T> co3 = d22*t24*t37; 
complex<T> co4 = -(d23*spb25*t46); 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t13*Int(ep,mu,c12,c345) + t5*Int(ep,mu,c15,c234) + t6*Int(ep,mu,c34,c125) + t11*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + co2*Int(ep,mu,c2,c5,c134) + t12*Int(ep,mu,c3,c4,c125) + co3*Int(ep,mu,c5,c1,c2,c34) + co4*Int(ep,mu,c5,c2,c1,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1y_qmqpQmQpgap_sl
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qm, Qp, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmqpQmQpgap sl");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);
	 vector<int> c5;  c5.push_back(5-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
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
 // #define TimeStamp "Mon 17 Nov 2008 23:57:19 on n303"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> t2 = spa34*T(2); 
complex<T> t8 = spa15*spa34; 
complex<T> t14 = square(spa12); 
complex<T> t15 = square(spa13); 
complex<T> t16 = square(spa35); 
complex<T> t17 = cube(spb25); 
complex<T> t18 = s12 + s15 - s34; 
complex<T> t22 = -(s25*spa13) - spa12*spa34*spb24*T(2); 
complex<T> t23 = square(spb45); 
complex<T> t24 = s12*spb15; 
complex<T> d3 = spa12*spa35*spa45*T(2); d3 = T(1)/d3;
complex<T> d6 = (s12 - s34)*spa25; d6 = T(1)/d6;
complex<T> d7 = spa25*square(s12 - s34)*T(2); d7 = T(1)/d7;
complex<T> d13 = spa25*spa34*spa45; d13 = T(1)/d13;
complex<T> d14 = spa35*spa45; d14 = T(1)/d14;
complex<T> d15 = (s12 + s25 - s34)*spa25*spa34; d15 = T(1)/d15;
complex<T> d18 = (s12 + s25 - s34)*spa34; d18 = T(1)/d18;
complex<T> d19 = spa15*spa25; d19 = T(1)/d19;
complex<T> d20 = (s12 + s25 - s34)*spa25; d20 = T(1)/d20;
complex<T> t1 = square(t18); 
complex<T> t27 = spa15*t2; 
complex<T> t33 = t14*t16; 
complex<T> t47 = spb15*t15; 
complex<T> t48 = t15*t24; 
complex<T> d2 = spa12*spa25*spa45*t2; d2 = T(1)/d2;
complex<T> d12 = spa35*t8; d12 = T(1)/d12;
complex<T> d16 = t8*cube(t18); d16 = T(1)/d16;
complex<T> d17 = spa34*cube(t18); d17 = T(1)/d17;
complex<T> d21 = spa15*cube(t18); d21 = T(1)/d21;
complex<T> d22 = t2*cube(t18); d22 = T(1)/d22;
complex<T> d23 = (s12 + s25 - s34)*t2; d23 = T(1)/d23;
complex<T> t11 = d14*spb12*t15 - d13*spa24*spb12*t15 + d16*s12*t17*t33 + d15*t48 + d12*spb12*cube(spa13); 
complex<T> t12 = spb34*(d19*t15 + d21*t17*t33 + d20*t47); 
complex<T> t26 = -t33; 
complex<T> d1 = spa12*spa35*t27; d1 = T(1)/d1;
complex<T> d4 = (s12 - s34)*t1*t8; d4 = T(1)/d4;
complex<T> d5 = t18*t27*square(s12 - s34); d5 = T(1)/d5;
complex<T> d8 = t27*square(s15 - s34); d8 = T(1)/d8;
complex<T> d9 = (s15 - s34)*t1*t8; d9 = T(1)/d9;
complex<T> d10 = t18*t27*square(s15 - s34); d10 = T(1)/d10;
complex<T> d11 = spa25*t27; d11 = T(1)/d11;
complex<T> t6 = d6*spa13*spb45 - d8*spa13*spb25*t22 + d10*t17*t33 + d4*t17*t33 + d5*t17*t33 + d9*t17*t33 + d7*t23*t8 + d11*t15*T(3); 
complex<T> t37 = t17*t26; 
complex<T> t5 = d8*spa13*spb25*t22 + (d10 + d9)*t37; 
complex<T> t13 = -(d6*spa13*spb45) - d7*spa15*spa34*t23 + d4*t37 + d5*t37 + d3*t15*T(3) - d2*spa24*t15*T(3) + d1*cube(spa13)*T(3); 
complex<T> co1 = d17*spb15*t37; 
complex<T> co2 = -(d18*spb25*t47); 
complex<T> co3 = d22*t24*t37; 
complex<T> co4 = -(d23*spb25*t48); 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t13*Int(ep,mu,c12,c345) + t5*Int(ep,mu,c15,c234) + t6*Int(ep,mu,c34,c125) + t11*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + co2*Int(ep,mu,c2,c5,c134) + t12*Int(ep,mu,c3,c4,c125) + co3*Int(ep,mu,c5,c1,c2,c34) + co4*Int(ep,mu,c5,c2,c1,c34));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qpQmQpqmgap_nf C2q2Q1y_4310_nf
#define _C_qpQpQmqmgap_nf C2q2Q1y_4280_nf
#define _C_qmQpQmqpgap_nf C2q2Q1y_4495_nf
#define _C_qmQmQpqpgap_nf C2q2Q1y_4525_nf
#define _C_qpQmQpqmgap_L C2q2Q1y_4310_L
#define _C_qpQpQmqmgap_L C2q2Q1y_4280_L
#define _C_qmQpQmqpgap_L C2q2Q1y_4495_L
#define _C_qmQmQpqpgap_L C2q2Q1y_4525_L
#define _C_qpqmQpQmgap_sl C2q2Q1y_4940_sl
#define _C_qpqmQmQpgap_sl C2q2Q1y_5120_sl
#define _C_qmqpQpQmgap_sl C2q2Q1y_4945_sl
#define _C_qmqpQmQpgap_sl C2q2Q1y_5125_sl
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qpQmQpqmgap_nf case 4310 : \
          return &C2q2Q1y_4310_nf
#define _CASE_qpQpQmqmgap_nf case 4280 : \
          return &C2q2Q1y_4280_nf
#define _CASE_qmQpQmqpgap_nf case 4495 : \
          return &C2q2Q1y_4495_nf
#define _CASE_qmQmQpqpgap_nf case 4525 : \
          return &C2q2Q1y_4525_nf
#define _CASE_qpQmQpqmgap_L case 4310 : \
          return &C2q2Q1y_4310_L
#define _CASE_qpQpQmqmgap_L case 4280 : \
          return &C2q2Q1y_4280_L
#define _CASE_qmQpQmqpgap_L case 4495 : \
          return &C2q2Q1y_4495_L
#define _CASE_qmQmQpqpgap_L case 4525 : \
          return &C2q2Q1y_4525_L
#define _CASE_qpqmQpQmgap_sl case 4940 : \
          return &C2q2Q1y_4940_sl
#define _CASE_qpqmQmQpgap_sl case 5120 : \
          return &C2q2Q1y_5120_sl
#define _CASE_qmqpQpQmgap_sl case 4945 : \
          return &C2q2Q1y_4945_sl
#define _CASE_qmqpQmQpgap_sl case 5125 : \
          return &C2q2Q1y_5125_sl
 
 
 // *************** function definitions using macros ************* 
 
template <class T> SeriesC<T> _C_qpQmQpqmgap_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qpQmQpqmgap_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qpQpQmqmgap_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qpQpQmqmgap_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmQpQmqpgap_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qmQpQmqpgap_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmQmQpqpgap_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qmQmQpqpgap_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qpQmQpqmgap_L(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qpQmQpqmgap_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qpQpQmqmgap_L(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qpQpQmqmgap_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmQpQmqpgap_L(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qmQpQmqpgap_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmQmQpqpgap_L(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qmQmQpqpgap_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQpQmgap_sl(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qpqmQpQmgap_sl(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQmQpgap_sl(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qpqmQmQpgap_sl(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQpQmgap_sl(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qmqpQpQmgap_sl(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQmQpgap_sl(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1y_qmqpQmQpgap_sl(ep,mu);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> SeriesC<T> ( *C2q2Q1y_L_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qpQmQpqmgap_L;
       _CASE_qpQpQmqmgap_L;
       _CASE_qmQpQmqpgap_L;
       _CASE_qmQmQpqpgap_L;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q2Q1y_nf_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qpQmQpqmgap_nf;
       _CASE_qpQpQmqmgap_nf;
       _CASE_qmQpQmqpgap_nf;
       _CASE_qmQmQpqpgap_nf;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q2Q1y_sl_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qpqmQpQmgap_sl;
       _CASE_qpqmQmQpgap_sl;
       _CASE_qmqpQpQmgap_sl;
       _CASE_qmqpQmQpgap_sl;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template SeriesC<R> ( *C2q2Q1y_L_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q1y_L_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q1y_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q2Q1y_nf_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q1y_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q1y_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q2Q1y_sl_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q1y_sl_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q1y_sl_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);




}
