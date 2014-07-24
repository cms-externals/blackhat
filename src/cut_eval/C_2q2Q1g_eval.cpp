/*
*C_2q2Q1g.cpp
*
* Created on 11/26, 2008
*      Author: Zvi's script
*/
 
#include "C_2q2Q1g_eval.h"
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


 
 
template <class T> SeriesC<T> C2q2Q1g_qmqpQpQmp_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qp, Qm, p}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpQmp PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:22:25 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = S(2,3);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> t4 = spa15*spa34; 
complex<T> t5 = spa12*spa45; 
complex<T> t6 = square(s15 - s34); 
complex<T> t7 = square(s12 - s45); 
complex<T> t9 = square(spa14); 
complex<T> t20 = cube(spa14); 
complex<T> t21 = square(spa13); 
complex<T> t22 = square(spa24); 
complex<T> t23 = square(spb25); 
complex<T> t24 = square(spa15); 
complex<T> t25 = square(spa45); 
complex<T> t26 = square(spb35); 
complex<T> t27 = s12 - s34; 
complex<T> t28 = square(s25); 
complex<T> t29 = s12 + s15 - s34; 
complex<T> t31 = cube(spa35); 
complex<T> t34 = spa12*T(2); 
complex<T> t39 = cube(spa13); 
complex<T> t40 = cube(spa24); 
complex<T> t44 = s15 + s45; 
complex<T> t47 = spb12*spb15; 
complex<T> t61 = spa12*spa35; 
complex<T> t79 = -(spa14*spa45); 
complex<T> t94 = spa14*spa15; 
complex<T> d43 = spa34*spa45*T(2); d43 = T(1)/d43;
complex<T> t3 = square(s15 + t27); 
complex<T> t10 = spa34*t5; 
complex<T> t32 = -(spb35*t4) - spb25*t5; 
complex<T> t37 = -t47; 
complex<T> t69 = t23*t24; 
complex<T> t70 = t25*t39; 
complex<T> t71 = spb45*t20; 
complex<T> t85 = t21*t26; 
complex<T> t93 = s35*t61; 
complex<T> t105 = spb12*t20; 
complex<T> d1 = t4*t5*T(3); d1 = T(1)/d1;
complex<T> d2 = t27*(s15 + t27)*t4*t5; d2 = T(1)/d2;
complex<T> d3 = t27*t4*t5; d3 = T(1)/d3;
complex<T> d4 = (s12 - s45)*t4*t5; d4 = T(1)/d4;
complex<T> d5 = t27*t4*t5*T(2); d5 = T(1)/d5;
complex<T> d8 = t4*T(2); d8 = T(1)/d8;
complex<T> d10 = square(t27); d10 = T(1)/d10;
complex<T> d11 = (s12 - s45)*t5; d11 = T(1)/d11;
complex<T> d12 = t34*t4*t7; d12 = T(1)/d12;
complex<T> d13 = spa35*t34*t4*square(t27); d13 = T(1)/d13;
complex<T> d15 = spa35*t34*t4*t7; d15 = T(1)/d15;
complex<T> d17 = t27*t4*t5*T(6); d17 = T(1)/d17;
complex<T> d18 = cube(t27)*T(3); d18 = T(1)/d18;
complex<T> d19 = t5*T(2); d19 = T(1)/d19;
complex<T> d20 = t34*t4; d20 = T(1)/d20;
complex<T> d21 = (s15 - s34)*(s15 + t27)*t4*t5; d21 = T(1)/d21;
complex<T> d22 = (s15 - s34)*t4; d22 = T(1)/d22;
complex<T> d26 = spa45*t4; d26 = T(1)/d26;
complex<T> d28 = t31*t4; d28 = T(1)/d28;
complex<T> d33 = spa12*spa15*t31; d33 = T(1)/d33;
complex<T> d35 = spa12*t31*t4; d35 = T(1)/d35;
complex<T> d36 = spa12*t4; d36 = T(1)/d36;
complex<T> d37 = spa45*t4*T(2); d37 = T(1)/d37;
complex<T> d40 = spa15*t5*T(2); d40 = T(1)/d40;
complex<T> d41 = spa15*t34; d41 = T(1)/d41;
complex<T> d42 = spa34*t34; d42 = T(1)/d42;
complex<T> d44 = spa15*t31*t34; d44 = T(1)/d44;
complex<T> t18 = d35*s45*t70 + d36*t71; 
complex<T> t49 = d41*spb34; 
complex<T> t58 = -t69; 
complex<T> t59 = -t70; 
complex<T> t73 = d19*spb35; 
complex<T> d6 = spa25*t10*square(t27)*T(2); d6 = T(1)/d6;
complex<T> d7 = spa25*t10*t27*(s15 + t27); d7 = T(1)/d7;
complex<T> d9 = t10*T(2); d9 = T(1)/d9;
complex<T> d14 = t27*t4*t93; d14 = T(1)/d14;
complex<T> d16 = (s12 - s45)*t4*t93; d16 = T(1)/d16;
complex<T> d23 = t10*t6*T(2); d23 = T(1)/d23;
complex<T> d24 = spa25*t10*t6*T(2); d24 = T(1)/d24;
complex<T> d25 = (s15 - s34)*spa25*t10*(s15 + t27); d25 = T(1)/d25;
complex<T> d27 = spa45*t3*t4; d27 = T(1)/d27;
complex<T> d29 = spa25*spa34*spa45*t3; d29 = T(1)/d29;
complex<T> d30 = t10*t3; d30 = T(1)/d30;
complex<T> d31 = spa25*t10*t3; d31 = T(1)/d31;
complex<T> d32 = spa15*t3*t5; d32 = T(1)/d32;
complex<T> d34 = spa25*t3*t5; d34 = T(1)/d34;
complex<T> d38 = spa34*spa45*t3*T(2); d38 = T(1)/d38;
complex<T> d39 = spa25*spa34*spa45*t3*T(2); d39 = T(1)/d39;
complex<T> t14 = d30*spb15*t20*t28 + d31*s15*t40*t69; 
complex<T> t16 = spb34*(-(d32*t20*t28) + d33*t59 + d34*t40*t69); 
complex<T> t17 = d4*s35*t20 + d15*t26*t59 + d16*t26*t59 + d12*spa14*spa45*t85 + d11*spb35*square(spa14); 
complex<T> t38 = -t49; 
complex<T> t65 = -t73; 
complex<T> t78 = t40*t58; 
complex<T> t1 = -(d1*t20) - d3*s35*t20 - d4*s35*t20 + d2*t20*t28 - d18*spb25*spb35*t32 - d5*t20*t44 + d13*t26*t70 + d14*t26*t70 + d15*t26*t70 + d16*t26*t70 + d6*t78 + d7*t78 + d12*t79*t85 + d10*d20*t79*t85 - d10*d8*s35*spb25*t9 - d11*spb35*t9 + d17*t32*t9 + d10*s25*t73*t9 + d10*d9*t22*t23*t94; 
complex<T> t2 = -(d1*t20) + d3*s35*t20 - d2*t20*t28 - d21*t20*t28 + d18*spb25*spb35*t32 + d5*t20*t44 + d13*t26*t59 + d14*t26*t59 + d24*t40*t69 + d25*t40*t69 + d6*t40*t69 + d7*t40*t69 + d10*d20*spa14*spa45*t85 + d22*spb25*t9 + d10*d8*s35*spb25*t9 - d17*t32*t9 + d10*s25*t65*t9 - d23*t22*t23*t94 - d10*d9*t22*t23*t94; 
complex<T> t13 = d21*t20*t28 + d24*t78 + d25*t78 + d23*t22*t23*t94 - d22*spb25*square(spa14); 
complex<T> t15 = d27*t105*t28 + d28*spb12*t70 + d29*spb12*t78 - d26*t105*T(2); 
complex<T> t56 = d44*s45*spb34*t59 + t38*t71; 
complex<T> t57 = d38*t20*t28*t37 + d39*s15*spb12*t78; 
complex<T> co1 = -(d37*s23*t105); 
complex<T> co2 = -(d40*s23*spb34*t20); 
complex<T> co3 = t49*t71; 
complex<T> co4 = d42*spb15*t71; 
complex<T> co5 = d43*t20*t47; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t1*Int(ep,mu,c12,c345) + t13*Int(ep,mu,c15,c234) + t2*Int(ep,mu,c34,c125) + t17*Int(ep,mu,c45,c123) + t15*Int(ep,mu,c1,c2,c345) + t14*Int(ep,mu,c1,c5,c234) + t16*Int(ep,mu,c3,c4,c125) + t18*Int(ep,mu,c4,c5,c123) + co1*Int(ep,mu,c1,c2,c3,c45) + t57*Int(ep,mu,c2,c1,c5,c34) + co2*Int(ep,mu,c2,c3,c4,c15) + co3*Int(ep,mu,c3,c4,c5,c12) + co4*Int(ep,mu,c4,c5,c1,c23) + co5*Int(ep,mu,c5,c1,c2,c34) + t56*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmqpQpQmp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qp, Qm, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpQmp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:22:26 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb15 = SPB(1,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t5 = spa15*spa34; 
complex<T> t6 = spa12*spa45; 
complex<T> t10 = cube(spa14); 
complex<T> d3 = cube(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> t8 = -(spb35*t5) - spb25*t6; 
complex<T> d1 = t5*t6*T(3); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*t5*t6*T(3); d2 = T(1)/d2;
complex<T> t2 = d1*t10 - t8*(d3*spb15*spb35 + d2*square(spa14)); 
complex<T> t3 = d1*t10 + d3*spb15*spb35*t8 + d2*t8*square(spa14); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*Int(ep,mu,c12,c345) + t2*Int(ep,mu,c34,c125));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qpqmQpQmp_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qp, Qm, p}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qpqmQpQmp PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:23:04 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> s23 = S(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s25 = -(spa25*spb25);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t3 = spa15*T(2); 
complex<T> t5 = spa12*spa34; 
complex<T> t6 = square(s12 - s45); 
complex<T> t20 = square(spa24); 
complex<T> t21 = square(spb35); 
complex<T> t22 = square(spa23); 
complex<T> t23 = square(spa45); 
complex<T> t30 = spa12*spb13 - spa24*spb34; 
complex<T> t31 = -(spa14*spb45); 
complex<T> t45 = s25*spb25; 
complex<T> t100 = spa13*spa24; 
complex<T> t101 = spa45*spb35; 
complex<T> d4 = s13*(s12 - s45)*spa15*spa45; d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*(s12 + s15 - s34)*spa15*spa34; d5 = T(1)/d5;
complex<T> d6 = (s12 - s45)*spa12*spa15*spa45; d6 = T(1)/d6;
complex<T> d7 = spa15*spa34*square(s12 - s34); d7 = T(1)/d7;
complex<T> d8 = spa12*spa45*square(s12 - s34)*T(3); d8 = T(1)/d8;
complex<T> d14 = cube(s12 - s34)*T(3); d14 = T(1)/d14;
complex<T> d16 = spa12*spa45*T(2); d16 = T(1)/d16;
complex<T> d17 = s12 - s34; d17 = T(1)/d17;
complex<T> d18 = spa12*T(2); d18 = T(1)/d18;
complex<T> d20 = square(s12 - s34); d20 = T(1)/d20;
complex<T> d21 = (s15 - s34)*spa15*spa34; d21 = T(1)/d21;
complex<T> d22 = (s15 - s34)*(s12 + s15 - s34)*spa15*spa34; d22 = T(1)/d22;
complex<T> d23 = (s23 - s45)*spa15*spa45; d23 = T(1)/d23;
complex<T> d24 = s13*(s23 - s45)*spa15*spa45; d24 = T(1)/d24;
complex<T> d26 = spa15*spa45*square(spa13); d26 = T(1)/d26;
complex<T> d27 = spa15*spa34*spa45; d27 = T(1)/d27;
complex<T> d28 = spa15*spa34*cube(spa35); d28 = T(1)/d28;
complex<T> d29 = spa15*spa34*square(s12 + s15 - s34); d29 = T(1)/d29;
complex<T> d30 = spa34*square(s12 + s15 - s34); d30 = T(1)/d30;
complex<T> d31 = spa12*spa15*cube(spa35); d31 = T(1)/d31;
complex<T> d32 = spa15*square(s12 + s15 - s34); d32 = T(1)/d32;
complex<T> d34 = spa15*square(spa13); d34 = T(1)/d34;
complex<T> d37 = spa34*square(s12 + s15 - s34)*T(2); d37 = T(1)/d37;
complex<T> d42 = spa34*spa45*T(2); d42 = T(1)/d42;
complex<T> t24 = -(spa14*spa23) - t100; 
complex<T> t26 = -t45; 
complex<T> t40 = d14*spa25; 
complex<T> t61 = spa14*t20; 
complex<T> t62 = t21*t23; 
complex<T> t63 = spa13*t22; 
complex<T> t64 = spa35*t5; 
complex<T> t80 = spa24*t21; 
complex<T> t81 = d18*spa23; 
complex<T> t98 = spa45*t22; 
complex<T> d1 = spa45*t3*t5; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa15*spa45*t5; d2 = T(1)/d2;
complex<T> d3 = (s12 - s45)*spa15*spa45*t5; d3 = T(1)/d3;
complex<T> d9 = t3*t5*t6; d9 = T(1)/d9;
complex<T> d15 = spa34*t3; d15 = T(1)/d15;
complex<T> d19 = t3*t5; d19 = T(1)/d19;
complex<T> d25 = spa15*spa45*t5*T(6); d25 = T(1)/d25;
complex<T> d33 = spa15*t5*cube(spa35); d33 = T(1)/d33;
complex<T> d35 = spa15*t5; d35 = T(1)/d35;
complex<T> d36 = spa34*spa45*t3; d36 = T(1)/d36;
complex<T> d38 = spa12*spa45*t3; d38 = T(1)/d38;
complex<T> d39 = spa45*t3*square(spa13); d39 = T(1)/d39;
complex<T> d40 = spa12*t3; d40 = T(1)/d40;
complex<T> d41 = t5*T(2); d41 = T(1)/d41;
complex<T> d43 = spa12*t3*cube(spa35); d43 = T(1)/d43;
complex<T> t35 = -(d15*T(3)); 
complex<T> t42 = -t61; 
complex<T> t43 = -t63; 
complex<T> t47 = d19*s45; 
complex<T> t48 = d40*spb34; 
complex<T> t49 = d24*t24; 
complex<T> t51 = d15*T(3); 
complex<T> t69 = d26*t24; 
complex<T> t72 = -t80; 
complex<T> t73 = spa14*t24; 
complex<T> t77 = spb15*t40; 
complex<T> t85 = t20*t26; 
complex<T> t92 = spb45*t61; 
complex<T> t110 = spb12*t61; 
complex<T> d10 = t3*t64*square(s12 - s34); d10 = T(1)/d10;
complex<T> d11 = (s12 - s34)*s35*spa15*t64; d11 = T(1)/d11;
complex<T> d12 = t3*t6*t64; d12 = T(1)/d12;
complex<T> d13 = s35*(s12 - s45)*spa15*t64; d13 = T(1)/d13;
complex<T> t1 = d7*spb15*t100*t101 + d21*spb25*t20 + d16*d17*spb35*t20 + d17*spb25*t20*t35 + d25*t42 + d22*t20*t45 + d5*t20*t45 + d20*t20*t31*t47 + d2*s35*t61 + d10*t43*t62 + d11*t43*t62 + d8*spa25*spa34*t72 + d20*t72*t81 + spa34*t21*t77*T(2); 
complex<T> t2 = -(d7*spb15*t100*t101) - d16*d17*spb35*t20 - d6*spa14*spa24*t30 + d1*t42 + d2*s35*t42 + d3*s35*t42 + d17*spb25*t20*t51 + d10*t62*t63 + d11*t62*t63 + d12*t62*t63 + d13*t62*t63 + d8*spa25*spa34*t80 + d20*t80*t81 + d5*t85 + d20*t47*t92 - d9*spa14*t21*t98 + d4*t73*square(spb13) - spa34*t21*t77*T(2); 
complex<T> t8 = spa14*spb13*(-(d23*spa24) + spb13*t49); 
complex<T> t10 = d23*spa14*spa24*spb13 + d6*spa14*spa24*t30 + d3*s35*t61 + d12*t43*t62 + d13*t43*t62 + d9*spa14*t21*t98 - spa14*t49*square(spb13) - d4*t73*square(spb13); 
complex<T> t11 = d33*s45*t23*t63 + d34*spb45*t73 + d35*t92; 
complex<T> t19 = spb34*(d31*t23*t43 + d32*t20*t45); 
complex<T> t59 = -(d21*spb25*t20) + d22*t85; 
complex<T> t60 = d43*s45*spb34*t23*t43 + t20*t31*t48; 
complex<T> t99 = spa14*t69; 
complex<T> t9 = d29*s12*t20*t45 + d28*spb12*t23*t63 + s12*t99 - d27*t110*T(2); 
complex<T> co1 = d30*spb15*t85; 
complex<T> co2 = s23*t99; 
complex<T> co3 = d36*s23*spb12*t42; 
complex<T> co4 = d37*s12*spb15*t85; 
complex<T> co5 = d38*s23*spb34*t42; 
complex<T> co6 = d39*s12*s23*t73; 
complex<T> co7 = t48*t92; 
complex<T> co8 = d41*spb15*t92; 
complex<T> co9 = d42*spb15*t110; 
complex<T> co10 = Complex(0,1); 
SeriesC<T> result = co10*(t2*Int(ep,mu,c12,c345) + t59*Int(ep,mu,c15,c234) + t8*Int(ep,mu,c23,c145) + t1*Int(ep,mu,c34,c125) + t10*Int(ep,mu,c45,c123) + t9*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + co2*Int(ep,mu,c2,c3,c145) + t19*Int(ep,mu,c3,c4,c125) + t11*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c1,c2,c3,c45) + co4*Int(ep,mu,c2,c1,c5,c34) + co5*Int(ep,mu,c2,c3,c4,c15) + co6*Int(ep,mu,c3,c2,c1,c45) + co7*Int(ep,mu,c3,c4,c5,c12) + co8*Int(ep,mu,c4,c5,c1,c23) + co9*Int(ep,mu,c5,c1,c2,c34) + t60*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qpqmQpQmp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qp, Qm, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qpqmQpQmp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:23:05 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t3 = spa12*spa45; 
complex<T> t6 = square(spb35); 
complex<T> t11 = spa25*spa34; 
complex<T> d3 = cube(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> t16 = t11*t6; 
complex<T> d1 = (s12 - s34)*t3; d1 = T(1)/d1;
complex<T> d2 = t3*square(s12 - s34)*T(3); d2 = T(1)/d2;
complex<T> d4 = spa15*spa34*t3*T(3); d4 = T(1)/d4;
complex<T> t4 = d2*spa24*t16 + d1*spb35*square(spa24) - d3*spb15*t16*T(2) + d4*spa14*square(spa24)*T(2); 
complex<T> t5 = -(d2*spa24*t16) - d1*spb35*square(spa24) + d3*spb15*t16*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t5*Int(ep,mu,c12,c345) + t4*Int(ep,mu,c34,c125));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qpqmQmQpp_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qm, Qp, p}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qpqmQmQpp PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:23:16 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> t4 = square(s12 - s34); 
complex<T> t5 = spa34*T(2); 
complex<T> t13 = square(spa23); 
complex<T> t14 = -s15 + s45; 
complex<T> t15 = -(s25*spb25); 
complex<T> t16 = spa25*spb45; 
complex<T> t17 = spa35*spb15; 
complex<T> t24 = -(spa12*T(2)); 
complex<T> t37 = spa13*spb15; 
complex<T> t38 = spa24*spb45; 
complex<T> t49 = spa14*spb12; 
complex<T> t55 = spa14*spb34; 
complex<T> d1 = spa12*spa15*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> d3 = (s12 - s34)*(s12 + s15 - s34)*spa15*spa34; d3 = T(1)/d3;
complex<T> d4 = (s12 - s34)*spa12*spa45; d4 = T(1)/d4;
complex<T> d7 = spa45*T(2); d7 = T(1)/d7;
complex<T> d8 = spa15*spa34; d8 = T(1)/d8;
complex<T> d9 = spa12*spa45; d9 = T(1)/d9;
complex<T> d10 = spa15*T(2); d10 = T(1)/d10;
complex<T> d11 = s12 - s34; d11 = T(1)/d11;
complex<T> d12 = spa15*spa34*T(6); d12 = T(1)/d12;
complex<T> d13 = spa12*spa45*T(6); d13 = T(1)/d13;
complex<T> d14 = spa15*spa45*T(3); d14 = T(1)/d14;
complex<T> d15 = cube(s12 - s34); d15 = T(1)/d15;
complex<T> d16 = (s15 - s34)*spa15*spa34; d16 = T(1)/d16;
complex<T> d17 = (s15 - s34)*(s12 + s15 - s34)*spa15*spa34; d17 = T(1)/d17;
complex<T> d18 = spa15*spa34*spa45; d18 = T(1)/d18;
complex<T> d19 = spa35*spa45; d19 = T(1)/d19;
complex<T> d20 = spa15*spa34*square(s12 + s15 - s34); d20 = T(1)/d20;
complex<T> d21 = spa34*square(s12 + s15 - s34); d21 = T(1)/d21;
complex<T> d22 = spa12*spa35*spa45; d22 = T(1)/d22;
complex<T> d23 = spa12*spa15*spa45; d23 = T(1)/d23;
complex<T> d24 = spa15*square(s12 + s15 - s34); d24 = T(1)/d24;
complex<T> d25 = spa12*spa35; d25 = T(1)/d25;
complex<T> d28 = spa12*spa15*spa45*T(2); d28 = T(1)/d28;
complex<T> d29 = spa12*spa15*T(2); d29 = T(1)/d29;
complex<T> d32 = spa12*spa35*T(2); d32 = T(1)/d32;
complex<T> t7 = -t17; 
complex<T> t19 = -t49; 
complex<T> t21 = d14*(s12 + s34); 
complex<T> t26 = -t37; 
complex<T> t27 = -t38; 
complex<T> t28 = -t55; 
complex<T> t29 = -(d10*T(3)); 
complex<T> t31 = d14*spb23; 
complex<T> t35 = spb25*t13; 
complex<T> t39 = d10*T(3); 
complex<T> t40 = -(d1*spa14); 
complex<T> t45 = t13*t15; 
complex<T> t46 = t16*t17; 
complex<T> t48 = -(d7*T(3)); 
complex<T> t50 = d9*spb35; 
complex<T> t54 = d7*T(3); 
complex<T> t69 = s34*t13; 
complex<T> d2 = spa15*spa45*t4*T(6); d2 = T(1)/d2;
complex<T> d5 = spa15*t4*T(2); d5 = T(1)/d5;
complex<T> d6 = spa45*t4*T(2); d6 = T(1)/d6;
complex<T> d26 = spa15*spa45*t5; d26 = T(1)/d26;
complex<T> d27 = t5*square(s12 + s15 - s34); d27 = T(1)/d27;
complex<T> d30 = spa12*t5; d30 = T(1)/d30;
complex<T> d31 = spa45*t5; d31 = T(1)/d31;
complex<T> t23 = -t35; 
complex<T> t25 = d2*t14; 
complex<T> t41 = -t50; 
complex<T> t53 = s25*t35; 
complex<T> t70 = t13*t19; 
complex<T> t71 = t13*t28; 
complex<T> t1 = -(d4*spb35*t13) - spa34*t16*t25 + d5*t16*t26 + d16*t35 + d6*t17*t38 + t13*t40 + d15*spa34*t24*t31*t46 + d17*t53 + d3*t53 + d11*(d8*t23 + d12*spa23*t26 + d13*spa23*t27 + spa23*spb45*t39 + t13*t41 + spa23*spb15*t54) + d15*spa14*t16*t21*t7 + spa12*t25*t7; 
complex<T> t2 = d4*spb35*t13 + spa34*t16*t25 + spa12*t17*t25 + d11*spa23*spb45*t29 + d11*d8*t35 + d11*d12*spa23*t37 + d5*t16*t37 + d11*d13*spa23*t38 + t13*t40 + d3*t45 + d15*spa14*t21*t46 + d11*spa23*spb15*t48 + d15*spa12*t31*t46*t5 + d11*t13*t50 + d6*t38*t7; 
complex<T> t12 = d16*t23 + d17*t45; 
complex<T> t52 = -(d19*spb12*t13) + d20*s12*t53 + d18*t70; 
complex<T> t62 = d24*spb34*t53 - d22*t69 + d23*t71; 
complex<T> co1 = d21*spb15*t45; 
complex<T> co2 = d25*spb45*t13; 
complex<T> co3 = d26*s23*t70; 
complex<T> co4 = d27*s12*spb15*t45; 
complex<T> co5 = d28*s23*t71; 
complex<T> co6 = d29*spb45*t13*t55; 
complex<T> co7 = d30*spa14*spb15*spb45*t13; 
complex<T> co8 = d31*spb15*t13*t49; 
complex<T> co9 = d32*spb45*t69; 
complex<T> co10 = Complex(0,1); 
SeriesC<T> result = co10*(t2*Int(ep,mu,c12,c345) + t12*Int(ep,mu,c15,c234) + t1*Int(ep,mu,c34,c125) + t52*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + t62*Int(ep,mu,c3,c4,c125) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c1,c2,c3,c45) + co4*Int(ep,mu,c2,c1,c5,c34) + co5*Int(ep,mu,c2,c3,c4,c15) + co6*Int(ep,mu,c3,c4,c5,c12) + co7*Int(ep,mu,c4,c5,c1,c23) + co8*Int(ep,mu,c5,c1,c2,c34) + co9*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qpqmQmQpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qm, Qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qpqmQmQpp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:23:19 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> s15 = -(spa15*spb15);
complex<T> s45 = -(spa45*spb45);
complex<T> t5 = square(spa23); 
complex<T> t9 = spa12*T(2); 
complex<T> t11 = -s15 + s45; 
complex<T> t12 = spa25*spb45; 
complex<T> t13 = spa35*spb15; 
complex<T> t14 = -(spa12*spa45*spb25) - spa15*spa34*spb35; 
complex<T> t18 = -(spa12*T(2)); 
complex<T> d1 = spa12*spa15*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> d2 = spa15*spa45*square(s12 - s34)*T(6); d2 = T(1)/d2;
complex<T> d4 = spa15*spa34*T(6); d4 = T(1)/d4;
complex<T> d5 = spa12*spa45*T(6); d5 = T(1)/d5;
complex<T> d6 = s12 - s34; d6 = T(1)/d6;
complex<T> d7 = spa15*spa45*T(3); d7 = T(1)/d7;
complex<T> d8 = cube(s12 - s34); d8 = T(1)/d8;
complex<T> t7 = -t13; 
complex<T> t15 = d7*(s12 + s34); 
complex<T> t17 = spa34*t12; 
complex<T> t22 = d1*spa14; 
complex<T> t24 = d7*spb23; 
complex<T> d3 = (s12 - s34)*spa15*spa34*spa45*t9; d3 = T(1)/d3;
complex<T> t33 = spa14*t15; 
complex<T> t40 = t13*t17; 
complex<T> t1 = -(d4*d6*spa13*spa23*spb15) - d5*d6*spa23*spa24*spb45 - d2*t11*t17 + d8*t18*t24*t40 + d3*t14*t5 + t22*t5 + d2*spa12*t11*t7 + d8*t12*t33*t7; 
complex<T> t2 = d4*d6*spa13*spa23*spb15 + d5*d6*spa23*spa24*spb45 + d2*spa12*t11*t13 + d2*t11*t17 + d8*t12*t13*t33 - d3*t14*t5 + t22*t5 + d8*t24*t40*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*Int(ep,mu,c12,c345) + t2*Int(ep,mu,c34,c125));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmpqpQpQm_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, Qp, Qm}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpqpQpQm PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:23:24 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> t6 = -(spa15*T(3)); 
complex<T> t9 = square(spa15); 
complex<T> t10 = s12*T(3) - s45*T(4); 
complex<T> t13 = spb12*spb23; 
complex<T> t15 = spa13*spb34; 
complex<T> t18 = -(s12*T(3)); 
complex<T> t22 = spa15*spb24; 
complex<T> d1 = spa12*spa23*spa45*T(4); d1 = T(1)/d1;
complex<T> d2 = s45*spa23*T(4); d2 = T(1)/d2;
complex<T> d3 = s45*(-s12 + s45)*spa23*T(2); d3 = T(1)/d3;
complex<T> d4 = spa23*square(s12 - s45)*T(2); d4 = T(1)/d4;
complex<T> d5 = s45*spa12*spa23*T(4); d5 = T(1)/d5;
complex<T> d6 = s45*spa23*square(s12 - s45)*T(2); d6 = T(1)/d6;
complex<T> d7 = (s12 - s45)*spa12*spa23*spa45*T(2); d7 = T(1)/d7;
complex<T> d8 = spa12*spa23*spa45*T(3); d8 = T(1)/d8;
complex<T> d9 = s45*spa23*T(2); d9 = T(1)/d9;
complex<T> d10 = s45*spa12*spa23*T(2); d10 = T(1)/d10;
complex<T> d11 = spa23*spa45; d11 = T(1)/d11;
complex<T> d12 = s45*spa23; d12 = T(1)/d12;
complex<T> d13 = spa12*spa45; d13 = T(1)/d13;
complex<T> d14 = s45; d14 = T(1)/d14;
complex<T> d15 = s45*spa12; d15 = T(1)/d15;
complex<T> d16 = spa23; d16 = T(1)/d16;
complex<T> d17 = spa12*spa23; d17 = T(1)/d17;
complex<T> d18 = spa45*T(2); d18 = T(1)/d18;
complex<T> d19 = s45*T(2); d19 = T(1)/d19;
complex<T> d20 = spa12*spa45*T(2); d20 = T(1)/d20;
complex<T> d21 = spa12*spa23*T(2); d21 = T(1)/d21;
complex<T> d22 = spa23*spa45*T(2); d22 = T(1)/d22;
complex<T> t5 = -t15; 
complex<T> t12 = -(d1*T(3)); 
complex<T> t17 = -t22; 
complex<T> t20 = d6*t10; 
complex<T> t29 = d4*spb23; 
complex<T> t32 = d2*spb24; 
complex<T> t37 = d5*t15; 
complex<T> t40 = -(spb45*t9); 
complex<T> t44 = spa15*t15; 
complex<T> t49 = -(spb12*t9); 
complex<T> t50 = -(spb23*t9); 
complex<T> t3 = spa12*spa35*spb24*t29 + spb12*t20*t44 - d7*spa12*spa35*spb23*t6 + d7*spa45*t5*t6 - d8*t9*T(2) + d9*t22*T(3) + d10*t44*T(3) + d3*t17*(t18 + s45*T(4)); 
complex<T> t28 = spa15*t5; 
complex<T> t43 = t32*t6 + t37*t6 + t12*t9; 
complex<T> t45 = d19*s12*spb23*t17 + d19*t13*t44 + d18*t13*t9; 
complex<T> t4 = spb12*t20*t28 - spa12*spa35*spb24*t29 + t32*t6 + t37*t6 + d7*spa15*(-(spa12*spa35*spb23) + spa45*t5)*T(3) + d1*t9*T(3) + d3*t22*(t18 + s45*T(4)); 
complex<T> t26 = d12*s12*t22 + d12*spb12*t28 + d11*t49; 
complex<T> t34 = d14*spb23*t17 + d15*spb23*t28 + d13*t50; 
complex<T> t39 = d16*t17 + d17*(t28 + t40); 
complex<T> co1 = d20*s34*t50; 
complex<T> co2 = d21*s34*t40; 
complex<T> co3 = d21*s15*t40; 
complex<T> co4 = d22*s15*t49; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t4*Int(ep,mu,c12,c345) + t43*Int(ep,mu,c23,c145) + t3*Int(ep,mu,c45,c123) + t26*Int(ep,mu,c1,c2,c345) + t34*Int(ep,mu,c2,c3,c145) + t39*Int(ep,mu,c4,c5,c123) + t45*Int(ep,mu,c1,c2,c3,c45) + co1*Int(ep,mu,c2,c3,c4,c15) + co2*Int(ep,mu,c3,c4,c5,c12) + co3*Int(ep,mu,c4,c5,c1,c23) + co4*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmpqpQmQp_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, Qm, Qp}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpqpQmQp PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:00 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s35 = -(spa35*spb35);
complex<T> t4 = spa23*T(2); 
complex<T> t6 = spa12*spa15; 
complex<T> t13 = square(spa14); 
complex<T> t14 = -(spa13*spb35); 
complex<T> t15 = -(spa35*T(3)); 
complex<T> t16 = spa15*spa34; 
complex<T> t19 = (s15 - s34)*spa23; 
complex<T> t28 = square(spb25); 
complex<T> t29 = square(spa45); 
complex<T> t30 = square(spa12); 
complex<T> t31 = square(spa15); 
complex<T> t32 = square(spa24); 
complex<T> t33 = spa14*spa35; 
complex<T> t35 = square(s25); 
complex<T> t36 = spa14*spa35 + spa13*spa45; 
complex<T> t37 = s12 - s34; 
complex<T> t39 = -(s12*T(3)); 
complex<T> t40 = s12*T(3) - s45*T(4); 
complex<T> t41 = square(spb35); 
complex<T> t48 = cube(spb25); 
complex<T> t54 = spa23*spa25; 
complex<T> t68 = cube(spa14); 
complex<T> t69 = spa13*spb12; 
complex<T> t75 = spb23*spb34; 
complex<T> t78 = spa13*spa14; 
complex<T> t83 = spb35*T(3); 
complex<T> t93 = spa14*spb25; 
complex<T> t95 = s45*spa13; 
complex<T> t97 = spa34*spb35; 
complex<T> t101 = spa35*spb12; 
complex<T> t114 = spa35*spb45; 
complex<T> d2 = spa12*spa23*spa45*T(4); d2 = T(1)/d2;
complex<T> d4 = s45*spa23*T(4); d4 = T(1)/d4;
complex<T> d14 = (s12 - s45)*spa12*spa23; d14 = T(1)/d14;
complex<T> d15 = s45*spa12*spa23*T(4); d15 = T(1)/d15;
complex<T> d18 = s35*(s12 - s45)*spa12*spa23; d18 = T(1)/d18;
complex<T> d25 = (s15 - s34)*spa12*spa34; d25 = T(1)/d25;
complex<T> d26 = spa15*spa23*square(s15 - s34); d26 = T(1)/d26;
complex<T> d34 = spa12*spa34*T(2); d34 = T(1)/d34;
complex<T> d35 = s15 - s34; d35 = T(1)/d35;
complex<T> d37 = spa34*T(2); d37 = T(1)/d37;
complex<T> d38 = square(s15 - s34); d38 = T(1)/d38;
complex<T> d39 = spa12*spa23*spa45*T(6); d39 = T(1)/d39;
complex<T> d43 = spa23*spa35; d43 = T(1)/d43;
complex<T> d44 = spa23*spa45; d44 = T(1)/d44;
complex<T> d45 = spa23*square(spa35); d45 = T(1)/d45;
complex<T> d46 = s45*spa23; d46 = T(1)/d46;
complex<T> d53 = spa12*spa45; d53 = T(1)/d53;
complex<T> d54 = s45; d54 = T(1)/d54;
complex<T> d55 = s45*spa12; d55 = T(1)/d55;
complex<T> d56 = spa12*spa23*spa35; d56 = T(1)/d56;
complex<T> d57 = spa12*spa23*square(spa35); d57 = T(1)/d57;
complex<T> d62 = spa23; d62 = T(1)/d62;
complex<T> d63 = spa12*spa23; d63 = T(1)/d63;
complex<T> d66 = s45*T(2); d66 = T(1)/d66;
complex<T> d70 = spa15*spa24*T(2); d70 = T(1)/d70;
complex<T> t5 = spa34*t37; 
complex<T> t12 = -(d63*spb45*t13) + d63*spb35*t78 + d62*t93 + d57*t36*t95 - d56*s45*t78*T(3); 
complex<T> t17 = square(t37); 
complex<T> t34 = s15 + t37; 
complex<T> t38 = -t69; 
complex<T> t43 = -(spa14*spb15) + t97; 
complex<T> t46 = -(s25*spa14) - spa15*t97*T(2); 
complex<T> t58 = d70*s34; 
complex<T> t71 = t39 + s45*T(4); 
complex<T> t79 = t28*t30; 
complex<T> t81 = spa35*t29; 
complex<T> t94 = t31*t32; 
complex<T> t99 = -(d2*T(3)); 
complex<T> t106 = spa24*t28; 
complex<T> t107 = -(spa13*t13); 
complex<T> t110 = d4*T(3); 
complex<T> t113 = d66*spb23; 
complex<T> t127 = spa14*t15; 
complex<T> t144 = spb23*(-(d53*t13) + d55*spb35*t78 + d54*t93); 
complex<T> t147 = spb12*t14; 
complex<T> d3 = spa34*spa45*t4*t6; d3 = T(1)/d3;
complex<T> d5 = s45*(-s12 + s45)*t4; d5 = T(1)/d5;
complex<T> d6 = t4*square(s12 - s45); d6 = T(1)/d6;
complex<T> d13 = spa12*spa23*t37; d13 = T(1)/d13;
complex<T> d16 = s45*t4*square(s12 - s45); d16 = T(1)/d16;
complex<T> d17 = s35*spa12*spa23*t37; d17 = T(1)/d17;
complex<T> d22 = spa23*spa34*t6; d22 = T(1)/d22;
complex<T> d24 = spa15*t19; d24 = T(1)/d24;
complex<T> d27 = t4*square(s15 - s34); d27 = T(1)/d27;
complex<T> d29 = spa25*t16*t4*square(s15 - s34); d29 = T(1)/d29;
complex<T> d33 = spa15*t4; d33 = T(1)/d33;
complex<T> d36 = t16*t4; d36 = T(1)/d36;
complex<T> d40 = s45*t4; d40 = T(1)/d40;
complex<T> d41 = s45*spa12*t4; d41 = T(1)/d41;
complex<T> d64 = t16*T(2); d64 = T(1)/d64;
complex<T> d65 = spa45*t16*T(2); d65 = T(1)/d65;
complex<T> d71 = t6*T(2); d71 = T(1)/d71;
complex<T> d72 = spa45*t6*T(2); d72 = T(1)/d72;
complex<T> d73 = spa12*spa35*t4; d73 = T(1)/d73;
complex<T> d74 = t4*t6; d74 = T(1)/d74;
complex<T> d75 = spa12*spa34*t4; d75 = T(1)/d75;
complex<T> d76 = spa34*t4; d76 = T(1)/d76;
complex<T> d77 = spa34*spa45*t4; d77 = T(1)/d77;
complex<T> d78 = spa12*t4*square(spa35); d78 = T(1)/d78;
complex<T> t3 = square(t34); 
complex<T> t10 = -(d6*spa12*spa34*spb23*spb25) + d16*spa14*t147*t40 - d18*spa13*t36*t41 + d5*t71*t93 - d41*spb35*t78*T(3) - d40*t93*T(3) - d14*spb35*t78*T(4) - d39*t13*T(13); 
complex<T> t44 = d75*spb15; 
complex<T> t45 = d74*spb34; 
complex<T> t50 = -t81; 
complex<T> t57 = d36*s12; 
complex<T> t60 = d33*spb24; 
complex<T> t73 = d3*T(3); 
complex<T> t80 = -t94; 
complex<T> t124 = spb23*t58; 
complex<T> t133 = t13*t38; 
complex<T> t162 = d15*t78*t83 + t110*t93 + t13*t99; 
complex<T> d1 = spa23*t34*t5*t6; d1 = T(1)/d1;
complex<T> d7 = spa23*t34*t5; d7 = T(1)/d7;
complex<T> d8 = t16*t17*t4; d8 = T(1)/d8;
complex<T> d9 = spa25*t16*t17*t4; d9 = T(1)/d9;
complex<T> d10 = spa15*t34*t5*t54; d10 = T(1)/d10;
complex<T> d12 = spa12*spa34*t17*t34*T(2); d12 = T(1)/d12;
complex<T> d19 = spa12*spa23*t5; d19 = T(1)/d19;
complex<T> d20 = spa12*spa34*t17*T(2); d20 = T(1)/d20;
complex<T> d21 = spa45*t4*t5*t6; d21 = T(1)/d21;
complex<T> d23 = spa34*t19*t34*t6; d23 = T(1)/d23;
complex<T> d28 = spa34*t19*t34; d28 = T(1)/d28;
complex<T> d30 = spa25*t16*t19*t34; d30 = T(1)/d30;
complex<T> d32 = spa12*spa34*t34*square(s15 - s34)*T(2); d32 = T(1)/d32;
complex<T> d49 = spa34*cube(t34); d49 = T(1)/d49;
complex<T> d52 = spa12*spa34*cube(t34); d52 = T(1)/d52;
complex<T> d61 = spa12*cube(t34); d61 = T(1)/d61;
complex<T> d69 = spa34*cube(t34)*T(2); d69 = T(1)/d69;
complex<T> t125 = spb15*(d76*t133 + d77*t101*t68); 
complex<T> t137 = t124*t13 + d71*t107*t75 + d72*spa35*t68*t75; 
complex<T> t148 = d64*spb23*t133 + spa14*t113*t147 + d65*spb23*t101*t68 + s12*t113*t93; 
complex<T> t153 = t44*(t114*t68 + t13*t95); 
complex<T> t158 = t114*t45*t68 + t13*t45*t95 - d73*s34*s45*t78*T(3); 
complex<T> d11 = spa12*t3*t5; d11 = T(1)/d11;
complex<T> d31 = (s15 - s34)*spa12*spa34*t3; d31 = T(1)/d31;
complex<T> d42 = spa23*t16*t3; d42 = T(1)/d42;
complex<T> d47 = spa23*spa34*t3; d47 = T(1)/d47;
complex<T> d48 = t16*t3*t54; d48 = T(1)/d48;
complex<T> d50 = spa12*spa23*spa34*t3; d50 = T(1)/d50;
complex<T> d51 = spa34*t3*t54; d51 = T(1)/d51;
complex<T> d58 = spa23*t3*t6; d58 = T(1)/d58;
complex<T> d59 = spa23*t3; d59 = T(1)/d59;
complex<T> d60 = spa15*t3*t54; d60 = T(1)/d60;
complex<T> d67 = spa34*t3*t4; d67 = T(1)/d67;
complex<T> d68 = spa25*spa34*t3*t4; d68 = T(1)/d68;
complex<T> t1 = -(d27*spa15*t106) + d25*spb25*t13 + d34*d35*spb25*t13 + d37*d38*spa14*spa45*t28 - d26*spa12*spb23*spb25*t33 + d23*t107*t35 + d38*t13*t57*t69 + d24*spb23*t78 + d29*t79*t81 + d30*t79*t81 + d31*t48*t94 + d32*t48*t94 - d22*spa13*t13*T(2) + d28*t106*t33*T(3) - d35*t13*t60*T(3); 
complex<T> t2 = d27*spa15*t106 + d22*t107 + d28*t106*t127 + d7*t106*t127 - d25*spb25*t13 - d34*d35*spb25*t13 - d37*d38*spa14*spa45*t28 + d8*spa12*spa13*t28*t29 + d26*spa12*spb23*spb25*t33 + d1*spa13*t13*t35 + d23*spa13*t13*t35 - d17*spa13*t36*t41 + d38*t133*t57 + spa35*t68*t73 - d24*spb23*t78 - d19*t43*t78 + d10*t50*t79 + d29*t50*t79 + d30*t50*t79 + d9*t50*t79 + d11*t48*t80 + d12*t48*t80 + d31*t48*t80 + d32*t48*t80 + d20*t46*t93 + d21*spa12*spa35*spa45*spb25*t13*T(3) + d21*spa35*spb35*t13*t16*T(3) + d35*t13*t60*T(3) - d13*spb35*t78*T(3); 
complex<T> t8 = d59*spb34*t106*t127 + d58*spa13*spb34*t13*t35 + d57*s34*spa13*t36 + d60*spb34*t50*t79 + d61*spb34*t48*t80 - d56*s34*t78*T(3); 
complex<T> t9 = d47*s12*t106*t127 - d44*spb12*t13 + d42*t133*t35 + d46*spa14*spb35*t69 + d45*t36*t69 + d48*s12*t50*t79 - d46*s12*t93 + d49*spb12*t48*t94 - d43*spa14*t69*T(3); 
complex<T> t11 = d6*spa12*spa34*spb23*spb25 + d21*t13*t15*(spa12*spa45*spb25 + spb35*t16) + t162 - d8*spa12*spa13*t28*t29 + d1*t107*t35 + d17*spa13*t36*t41 + d18*spa13*t36*t41 + d16*spa14*spb35*t40*t69 + spa35*t68*t73 + d19*t43*t78 + d10*t79*t81 + d9*t79*t81 + d13*t78*t83 - d20*t46*t93 - d5*t71*t93 + d11*t48*t94 + d12*t48*t94 + d7*t106*t33*T(3) + d14*spb35*t78*T(4); 
complex<T> t23 = d52*s15*t48*t80 + d51*spb15*t79*t81 - d50*spa13*spb15*t35*square(spa14) - d47*s15*t106*t33*T(3); 
complex<T> t24 = d67*s15*t106*t33*t39 + d68*s12*spb15*t79*t81 + d69*s15*spb12*t48*t94 + d67*spb15*t35*t69*square(spa14); 
complex<T> co1 = -(t124*t13); 
complex<T> co2 = d78*s34*t36*t95; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t11*Int(ep,mu,c12,c345) + t1*Int(ep,mu,c15,c234) + t162*Int(ep,mu,c23,c145) + t2*Int(ep,mu,c34,c125) + t10*Int(ep,mu,c45,c123) + t9*Int(ep,mu,c1,c2,c345) + t23*Int(ep,mu,c1,c5,c234) + t144*Int(ep,mu,c2,c3,c145) + t8*Int(ep,mu,c3,c4,c125) + t12*Int(ep,mu,c4,c5,c123) + t148*Int(ep,mu,c1,c2,c3,c45) + t24*Int(ep,mu,c2,c1,c5,c34) + t137*Int(ep,mu,c2,c3,c4,c15) + t158*Int(ep,mu,c3,c4,c5,c12) + co1*Int(ep,mu,c4,c3,c2,c15) + t153*Int(ep,mu,c4,c5,c1,c23) + t125*Int(ep,mu,c5,c1,c2,c34) + co2*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qppqmQmQp_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, p, qm, Qm, Qp}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qppqmQmQp BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:01 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t4 = square(spa34); 
complex<T> t7 = -(s23*T(3)) + s45*T(4); 
complex<T> t9 = spa12*T(2); 
complex<T> t12 = spa14*spa23; 
complex<T> d4 = spa12*spa23; d4 = T(1)/d4;
complex<T> d5 = spa45*T(2); d5 = T(1)/d5;
complex<T> t11 = spb23*t4; 
complex<T> d1 = spa45*t9*square(s23 - s45); d1 = T(1)/d1;
complex<T> d2 = t9*square(s23 - s45); d2 = T(1)/d2;
complex<T> d3 = spa23*spa45*t9; d3 = T(1)/d3;
complex<T> t2 = -(d2*spb12*spb25*t12) + d2*spa34*spb25*t7 - d1*t11*t7 + d3*t4*T(3); 
complex<T> t3 = d1*t11*t7 + d2*spb25*(spb12*t12 - spa34*t7); 
complex<T> co1 = d4*spb45*t4; 
complex<T> co2 = d5*spb12*t11; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*Int(ep,mu,c23,c145) + t2*Int(ep,mu,c45,c123) + co1*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmqppQmQp_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, Qm, Qp}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqppQmQp BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:02 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t6 = square(spa14); 
complex<T> t9 = square(spb35); 
complex<T> d1 = spa12*spa23*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa23*spa45; d2 = T(1)/d2;
complex<T> d3 = spa23*square(s12 - s45)*T(2); d3 = T(1)/d3;
complex<T> d4 = spa13*spa23*spa45; d4 = T(1)/d4;
complex<T> d5 = spa23*spa34*spa45; d5 = T(1)/d5;
complex<T> d6 = spa13*spa45; d6 = T(1)/d6;
complex<T> d7 = spa13*spa23; d7 = T(1)/d7;
complex<T> d8 = spa13*spa45*T(2); d8 = T(1)/d8;
complex<T> t4 = -(d2*spa14*spa24*spb23) - d2*spb13*t6 - d3*spa13*spa45*t9 + d1*spa24*t6*T(3); 
complex<T> t5 = d2*spa14*spa24*spb23 + d2*spb13*t6 + d3*spa13*spa45*t9; 
complex<T> t17 = s12*t6; 
complex<T> t3 = d4*t17 + d5*spa24*spb12*t6; 
complex<T> co1 = -(d6*spb23*t6); 
complex<T> co2 = d7*spb45*t6; 
complex<T> co3 = -(d8*spb23*t17); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t4*Int(ep,mu,c12,c345) + t5*Int(ep,mu,c45,c123) + t3*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c2,c3,c145) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmqppQpQm_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, Qp, Qm}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqppQpQm BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:03 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t5 = square(spa15); 
complex<T> t8 = square(spb34); 
complex<T> d1 = spa12*spa23*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa23; d2 = T(1)/d2;
complex<T> d3 = spa23*square(s12 - s45)*T(2); d3 = T(1)/d3;
complex<T> d4 = spa13*spa23*spa45; d4 = T(1)/d4;
complex<T> d5 = spa23*spa34*spa45; d5 = T(1)/d5;
complex<T> d6 = spa13*spa45; d6 = T(1)/d6;
complex<T> d7 = spa13*spa23; d7 = T(1)/d7;
complex<T> d8 = spa13*spa45*T(2); d8 = T(1)/d8;
complex<T> t3 = -(d2*spa15*spb34) - d3*spa13*spa45*t8 + d1*spa24*t5*T(3); 
complex<T> t4 = d2*spa15*spb34 + d3*spa13*spa45*t8; 
complex<T> t13 = s12*t5; 
complex<T> t2 = d4*t13 + d5*spa24*spb12*t5; 
complex<T> co1 = -(d6*spb23*t5); 
complex<T> co2 = d7*spb45*t5; 
complex<T> co3 = -(d8*spb23*t13); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t3*Int(ep,mu,c12,c345) + t4*Int(ep,mu,c45,c123) + t2*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c2,c3,c145) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qpqmpQmQp_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, p, Qm, Qp}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qpqmpQmQp BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:05 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> t2 = spa45*T(2); 
complex<T> t5 = square(s13); 
complex<T> t6 = spa23*spa45; 
complex<T> t9 = square(spa12); 
complex<T> t10 = square(spa34); 
complex<T> t11 = cube(spb13); 
complex<T> t19 = -(s13*spa24); 
complex<T> d9 = spa45*cube(spa13); d9 = T(1)/d9;
complex<T> d10 = spa23*cube(spa13); d10 = T(1)/d10;
complex<T> t16 = t19 + spa12*spa45*spb15*T(2); 
complex<T> t18 = spa23*t2; 
complex<T> t23 = t10*t9; 
complex<T> d3 = (s12 - s45)*t5*t6; d3 = T(1)/d3;
complex<T> d5 = (s23 - s45)*t5*t6; d5 = T(1)/d5;
complex<T> d7 = t6*cube(spa13); d7 = T(1)/d7;
complex<T> d8 = spa34*t6; d8 = T(1)/d8;
complex<T> d11 = t2*cube(spa13); d11 = T(1)/d11;
complex<T> t7 = d7*s12*t23 + d8*spb12*cube(spa24); 
complex<T> t17 = -t23; 
complex<T> d1 = spa12*spa34*t18; d1 = T(1)/d1;
complex<T> d2 = s13*t18*square(s12 - s45); d2 = T(1)/d2;
complex<T> d4 = s13*t18*square(s23 - s45); d4 = T(1)/d4;
complex<T> d6 = t18*square(s23 - s45); d6 = T(1)/d6;
complex<T> t4 = -(d6*spa24*spb13*t16) + t11*(d2*t17 + d4*t17 + (d3 + d5)*t23); 
complex<T> t8 = d3*t11*t17 + d2*t11*t23 + d1*cube(spa24)*T(3); 
complex<T> t21 = d5*t11*t17 + d6*spa24*spb13*(t19 + spa12*spb15*t2) + d4*t11*t23; 
complex<T> co1 = d9*spb23*t17; 
complex<T> co2 = d10*spb45*t23; 
complex<T> co3 = d11*s12*spb23*t17; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t8*Int(ep,mu,c12,c345) + t21*Int(ep,mu,c23,c145) + t4*Int(ep,mu,c45,c123) + t7*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c2,c3,c145) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qpqmpQpQm_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, p, Qp, Qm}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qpqmpQpQm BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:06 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> t2 = spa45*T(2); 
complex<T> t5 = square(s13); 
complex<T> t6 = spa23*spa45; 
complex<T> t10 = square(spa12); 
complex<T> t11 = square(spa35); 
complex<T> t12 = cube(spb13); 
complex<T> t16 = square(spa25); 
complex<T> t17 = -(s13*spa25) - spa12*spa45*spb14*T(2); 
complex<T> t20 = spa25*spb13; 
complex<T> d9 = spa45*cube(spa13); d9 = T(1)/d9;
complex<T> d10 = spa23*cube(spa13); d10 = T(1)/d10;
complex<T> t19 = spa23*t2; 
complex<T> t23 = t10*t11; 
complex<T> d3 = (s12 - s45)*t5*t6; d3 = T(1)/d3;
complex<T> d5 = (s23 - s45)*t5*t6; d5 = T(1)/d5;
complex<T> d7 = t6*cube(spa13); d7 = T(1)/d7;
complex<T> d8 = spa34*t6; d8 = T(1)/d8;
complex<T> d11 = t2*cube(spa13); d11 = T(1)/d11;
complex<T> t7 = d8*spa24*spb12*t16 + d7*s12*t23; 
complex<T> t18 = -t23; 
complex<T> d1 = spa12*spa34*t19; d1 = T(1)/d1;
complex<T> d2 = s13*t19*square(s12 - s45); d2 = T(1)/d2;
complex<T> d4 = s13*t19*square(s23 - s45); d4 = T(1)/d4;
complex<T> d6 = t19*square(s23 - s45); d6 = T(1)/d6;
complex<T> t8 = d3*t12*t18 + d2*t12*t23 + d1*spa24*t16*T(3); 
complex<T> t21 = d6*t17; 
complex<T> t4 = d2*t12*t18 + d4*t12*t18 - t20*t21 + d3*t12*t23 + d5*t12*t23; 
complex<T> t22 = d5*t12*t18 + t20*t21 + d4*t12*t23; 
complex<T> co1 = d9*spb23*t18; 
complex<T> co2 = d10*spb45*t23; 
complex<T> co3 = d11*s12*spb23*t18; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t8*Int(ep,mu,c12,c345) + t22*Int(ep,mu,c23,c145) + t4*Int(ep,mu,c45,c123) + t7*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c2,c3,c145) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmqpQmpQp_TriP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qm, p, Qp}, TriP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmpQp TriP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:07 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = spa34*spa45; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb12*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c12,c345) + co2*Int(ep,mu,c1,c2,c345));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmqpQppQm_TriP
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qp, p, Qm}, TriP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQppQm TriP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:07 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spa15); 
complex<T> d1 = spa12*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = spa34*spa45; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb12*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c12,c345) + co2*Int(ep,mu,c1,c2,c345));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmqpQmpQp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qm, p, Qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmpQp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:08 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qmqpQppQm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, Qp, p, Qm}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQppQm nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:08 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa15); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qpqmQmpQp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qm, p, Qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qpqmQmpQp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:09 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa23); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_qpqmQppQm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, Qp, p, Qm}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qpqmQppQm nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:09 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa25); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmQpqpqmp_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, Qp, qp, qm, p}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmQpqpqmp PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:48 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = S(2,3);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> t4 = spa15*spa34; 
complex<T> t5 = spa12*spa45; 
complex<T> t6 = square(s15 - s34); 
complex<T> t7 = square(s12 - s45); 
complex<T> t9 = square(spa14); 
complex<T> t20 = cube(spa14); 
complex<T> t21 = square(spa13); 
complex<T> t22 = square(spa24); 
complex<T> t23 = square(spb25); 
complex<T> t24 = square(spa15); 
complex<T> t25 = square(spa45); 
complex<T> t26 = square(spb35); 
complex<T> t27 = s12 - s34; 
complex<T> t28 = square(s25); 
complex<T> t29 = s12 + s15 - s34; 
complex<T> t31 = cube(spa35); 
complex<T> t34 = spa12*T(2); 
complex<T> t39 = cube(spa13); 
complex<T> t40 = cube(spa24); 
complex<T> t44 = s15 + s45; 
complex<T> t47 = spb12*spb15; 
complex<T> t61 = spa12*spa35; 
complex<T> t79 = -(spa14*spa45); 
complex<T> t94 = spa14*spa15; 
complex<T> d43 = spa34*spa45*T(2); d43 = T(1)/d43;
complex<T> t3 = square(s15 + t27); 
complex<T> t10 = spa34*t5; 
complex<T> t32 = -(spb35*t4) - spb25*t5; 
complex<T> t37 = -t47; 
complex<T> t69 = t23*t24; 
complex<T> t70 = t25*t39; 
complex<T> t71 = spb45*t20; 
complex<T> t85 = t21*t26; 
complex<T> t93 = s35*t61; 
complex<T> t105 = spb12*t20; 
complex<T> d1 = t4*t5*T(3); d1 = T(1)/d1;
complex<T> d2 = t27*(s15 + t27)*t4*t5; d2 = T(1)/d2;
complex<T> d3 = t27*t4*t5; d3 = T(1)/d3;
complex<T> d4 = (s12 - s45)*t4*t5; d4 = T(1)/d4;
complex<T> d5 = t27*t4*t5*T(2); d5 = T(1)/d5;
complex<T> d8 = t4*T(2); d8 = T(1)/d8;
complex<T> d10 = square(t27); d10 = T(1)/d10;
complex<T> d11 = (s12 - s45)*t5; d11 = T(1)/d11;
complex<T> d12 = t34*t4*t7; d12 = T(1)/d12;
complex<T> d13 = spa35*t34*t4*square(t27); d13 = T(1)/d13;
complex<T> d15 = spa35*t34*t4*t7; d15 = T(1)/d15;
complex<T> d17 = t27*t4*t5*T(6); d17 = T(1)/d17;
complex<T> d18 = cube(t27)*T(3); d18 = T(1)/d18;
complex<T> d19 = t5*T(2); d19 = T(1)/d19;
complex<T> d20 = t34*t4; d20 = T(1)/d20;
complex<T> d21 = (s15 - s34)*(s15 + t27)*t4*t5; d21 = T(1)/d21;
complex<T> d22 = (s15 - s34)*t4; d22 = T(1)/d22;
complex<T> d26 = spa45*t4; d26 = T(1)/d26;
complex<T> d28 = t31*t4; d28 = T(1)/d28;
complex<T> d33 = spa12*spa15*t31; d33 = T(1)/d33;
complex<T> d35 = spa12*t31*t4; d35 = T(1)/d35;
complex<T> d36 = spa12*t4; d36 = T(1)/d36;
complex<T> d37 = spa45*t4*T(2); d37 = T(1)/d37;
complex<T> d40 = spa15*t5*T(2); d40 = T(1)/d40;
complex<T> d41 = spa15*t34; d41 = T(1)/d41;
complex<T> d42 = spa34*t34; d42 = T(1)/d42;
complex<T> d44 = spa15*t31*t34; d44 = T(1)/d44;
complex<T> t18 = d35*s45*t70 + d36*t71; 
complex<T> t49 = d41*spb34; 
complex<T> t58 = -t69; 
complex<T> t59 = -t70; 
complex<T> t73 = d19*spb35; 
complex<T> d6 = spa25*t10*square(t27)*T(2); d6 = T(1)/d6;
complex<T> d7 = spa25*t10*t27*(s15 + t27); d7 = T(1)/d7;
complex<T> d9 = t10*T(2); d9 = T(1)/d9;
complex<T> d14 = t27*t4*t93; d14 = T(1)/d14;
complex<T> d16 = (s12 - s45)*t4*t93; d16 = T(1)/d16;
complex<T> d23 = t10*t6*T(2); d23 = T(1)/d23;
complex<T> d24 = spa25*t10*t6*T(2); d24 = T(1)/d24;
complex<T> d25 = (s15 - s34)*spa25*t10*(s15 + t27); d25 = T(1)/d25;
complex<T> d27 = spa45*t3*t4; d27 = T(1)/d27;
complex<T> d29 = spa25*spa34*spa45*t3; d29 = T(1)/d29;
complex<T> d30 = t10*t3; d30 = T(1)/d30;
complex<T> d31 = spa25*t10*t3; d31 = T(1)/d31;
complex<T> d32 = spa15*t3*t5; d32 = T(1)/d32;
complex<T> d34 = spa25*t3*t5; d34 = T(1)/d34;
complex<T> d38 = spa34*spa45*t3*T(2); d38 = T(1)/d38;
complex<T> d39 = spa25*spa34*spa45*t3*T(2); d39 = T(1)/d39;
complex<T> t14 = d30*spb15*t20*t28 + d31*s15*t40*t69; 
complex<T> t16 = spb34*(-(d32*t20*t28) + d33*t59 + d34*t40*t69); 
complex<T> t17 = d4*s35*t20 + d15*t26*t59 + d16*t26*t59 + d12*spa14*spa45*t85 + d11*spb35*square(spa14); 
complex<T> t38 = -t49; 
complex<T> t65 = -t73; 
complex<T> t78 = t40*t58; 
complex<T> t1 = -(d1*t20) - d3*s35*t20 - d4*s35*t20 + d2*t20*t28 - d18*spb25*spb35*t32 - d5*t20*t44 + d13*t26*t70 + d14*t26*t70 + d15*t26*t70 + d16*t26*t70 + d6*t78 + d7*t78 + d12*t79*t85 + d10*d20*t79*t85 - d10*d8*s35*spb25*t9 - d11*spb35*t9 + d17*t32*t9 + d10*s25*t73*t9 + d10*d9*t22*t23*t94; 
complex<T> t2 = -(d1*t20) + d3*s35*t20 - d2*t20*t28 - d21*t20*t28 + d18*spb25*spb35*t32 + d5*t20*t44 + d13*t26*t59 + d14*t26*t59 + d24*t40*t69 + d25*t40*t69 + d6*t40*t69 + d7*t40*t69 + d10*d20*spa14*spa45*t85 + d22*spb25*t9 + d10*d8*s35*spb25*t9 - d17*t32*t9 + d10*s25*t65*t9 - d23*t22*t23*t94 - d10*d9*t22*t23*t94; 
complex<T> t13 = d21*t20*t28 + d24*t78 + d25*t78 + d23*t22*t23*t94 - d22*spb25*square(spa14); 
complex<T> t15 = d27*t105*t28 + d28*spb12*t70 + d29*spb12*t78 - d26*t105*T(2); 
complex<T> t56 = d44*s45*spb34*t59 + t38*t71; 
complex<T> t57 = d38*t20*t28*t37 + d39*s15*spb12*t78; 
complex<T> co1 = -(d37*s23*t105); 
complex<T> co2 = -(d40*s23*spb34*t20); 
complex<T> co3 = t49*t71; 
complex<T> co4 = d42*spb15*t71; 
complex<T> co5 = d43*t20*t47; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t1*Int(ep,mu,c12,c345) + t13*Int(ep,mu,c15,c234) + t2*Int(ep,mu,c34,c125) + t17*Int(ep,mu,c45,c123) + t15*Int(ep,mu,c1,c2,c345) + t14*Int(ep,mu,c1,c5,c234) + t16*Int(ep,mu,c3,c4,c125) + t18*Int(ep,mu,c4,c5,c123) + co1*Int(ep,mu,c1,c2,c3,c45) + t57*Int(ep,mu,c2,c1,c5,c34) + co2*Int(ep,mu,c2,c3,c4,c15) + co3*Int(ep,mu,c3,c4,c5,c12) + co4*Int(ep,mu,c4,c5,c1,c23) + co5*Int(ep,mu,c5,c1,c2,c34) + t56*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmQpqpqmp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, Qp, qp, qm, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmQpqpqmp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:26:48 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb15 = SPB(1,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t5 = spa15*spa34; 
complex<T> t6 = spa12*spa45; 
complex<T> t10 = cube(spa14); 
complex<T> d3 = cube(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> t8 = -(spb35*t5) - spb25*t6; 
complex<T> d1 = t5*t6*T(3); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*t5*t6*T(3); d2 = T(1)/d2;
complex<T> t2 = d1*t10 - t8*(d3*spb15*spb35 + d2*square(spa14)); 
complex<T> t3 = d1*t10 + d3*spb15*spb35*t8 + d2*t8*square(spa14); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*Int(ep,mu,c12,c345) + t2*Int(ep,mu,c34,c125));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QpQmqpqmp_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, Qm, qp, qm, p}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QpQmqpqmp PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:27:26 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> s23 = S(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s25 = -(spa25*spb25);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t3 = spa15*T(2); 
complex<T> t5 = spa12*spa34; 
complex<T> t6 = square(s12 - s45); 
complex<T> t20 = square(spa24); 
complex<T> t21 = square(spb35); 
complex<T> t22 = square(spa23); 
complex<T> t23 = square(spa45); 
complex<T> t30 = spa12*spb13 - spa24*spb34; 
complex<T> t31 = -(spa14*spb45); 
complex<T> t45 = s25*spb25; 
complex<T> t100 = spa13*spa24; 
complex<T> t101 = spa45*spb35; 
complex<T> d4 = s13*(s12 - s45)*spa15*spa45; d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*(s12 + s15 - s34)*spa15*spa34; d5 = T(1)/d5;
complex<T> d6 = (s12 - s45)*spa12*spa15*spa45; d6 = T(1)/d6;
complex<T> d7 = spa15*spa34*square(s12 - s34); d7 = T(1)/d7;
complex<T> d8 = spa12*spa45*square(s12 - s34)*T(3); d8 = T(1)/d8;
complex<T> d14 = cube(s12 - s34)*T(3); d14 = T(1)/d14;
complex<T> d16 = spa12*spa45*T(2); d16 = T(1)/d16;
complex<T> d17 = s12 - s34; d17 = T(1)/d17;
complex<T> d18 = spa12*T(2); d18 = T(1)/d18;
complex<T> d20 = square(s12 - s34); d20 = T(1)/d20;
complex<T> d21 = (s15 - s34)*spa15*spa34; d21 = T(1)/d21;
complex<T> d22 = (s15 - s34)*(s12 + s15 - s34)*spa15*spa34; d22 = T(1)/d22;
complex<T> d23 = (s23 - s45)*spa15*spa45; d23 = T(1)/d23;
complex<T> d24 = s13*(s23 - s45)*spa15*spa45; d24 = T(1)/d24;
complex<T> d26 = spa15*spa45*square(spa13); d26 = T(1)/d26;
complex<T> d27 = spa15*spa34*spa45; d27 = T(1)/d27;
complex<T> d28 = spa15*spa34*cube(spa35); d28 = T(1)/d28;
complex<T> d29 = spa15*spa34*square(s12 + s15 - s34); d29 = T(1)/d29;
complex<T> d30 = spa34*square(s12 + s15 - s34); d30 = T(1)/d30;
complex<T> d31 = spa12*spa15*cube(spa35); d31 = T(1)/d31;
complex<T> d32 = spa15*square(s12 + s15 - s34); d32 = T(1)/d32;
complex<T> d34 = spa15*square(spa13); d34 = T(1)/d34;
complex<T> d37 = spa34*square(s12 + s15 - s34)*T(2); d37 = T(1)/d37;
complex<T> d42 = spa34*spa45*T(2); d42 = T(1)/d42;
complex<T> t24 = -(spa14*spa23) - t100; 
complex<T> t26 = -t45; 
complex<T> t40 = d14*spa25; 
complex<T> t61 = spa14*t20; 
complex<T> t62 = t21*t23; 
complex<T> t63 = spa13*t22; 
complex<T> t64 = spa35*t5; 
complex<T> t80 = spa24*t21; 
complex<T> t81 = d18*spa23; 
complex<T> t98 = spa45*t22; 
complex<T> d1 = spa45*t3*t5; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa15*spa45*t5; d2 = T(1)/d2;
complex<T> d3 = (s12 - s45)*spa15*spa45*t5; d3 = T(1)/d3;
complex<T> d9 = t3*t5*t6; d9 = T(1)/d9;
complex<T> d15 = spa34*t3; d15 = T(1)/d15;
complex<T> d19 = t3*t5; d19 = T(1)/d19;
complex<T> d25 = spa15*spa45*t5*T(6); d25 = T(1)/d25;
complex<T> d33 = spa15*t5*cube(spa35); d33 = T(1)/d33;
complex<T> d35 = spa15*t5; d35 = T(1)/d35;
complex<T> d36 = spa34*spa45*t3; d36 = T(1)/d36;
complex<T> d38 = spa12*spa45*t3; d38 = T(1)/d38;
complex<T> d39 = spa45*t3*square(spa13); d39 = T(1)/d39;
complex<T> d40 = spa12*t3; d40 = T(1)/d40;
complex<T> d41 = t5*T(2); d41 = T(1)/d41;
complex<T> d43 = spa12*t3*cube(spa35); d43 = T(1)/d43;
complex<T> t35 = -(d15*T(3)); 
complex<T> t42 = -t61; 
complex<T> t43 = -t63; 
complex<T> t47 = d19*s45; 
complex<T> t48 = d40*spb34; 
complex<T> t49 = d24*t24; 
complex<T> t51 = d15*T(3); 
complex<T> t69 = d26*t24; 
complex<T> t72 = -t80; 
complex<T> t73 = spa14*t24; 
complex<T> t77 = spb15*t40; 
complex<T> t85 = t20*t26; 
complex<T> t92 = spb45*t61; 
complex<T> t110 = spb12*t61; 
complex<T> d10 = t3*t64*square(s12 - s34); d10 = T(1)/d10;
complex<T> d11 = (s12 - s34)*s35*spa15*t64; d11 = T(1)/d11;
complex<T> d12 = t3*t6*t64; d12 = T(1)/d12;
complex<T> d13 = s35*(s12 - s45)*spa15*t64; d13 = T(1)/d13;
complex<T> t1 = d7*spb15*t100*t101 + d21*spb25*t20 + d16*d17*spb35*t20 + d17*spb25*t20*t35 + d25*t42 + d22*t20*t45 + d5*t20*t45 + d20*t20*t31*t47 + d2*s35*t61 + d10*t43*t62 + d11*t43*t62 + d8*spa25*spa34*t72 + d20*t72*t81 + spa34*t21*t77*T(2); 
complex<T> t2 = -(d7*spb15*t100*t101) - d16*d17*spb35*t20 - d6*spa14*spa24*t30 + d1*t42 + d2*s35*t42 + d3*s35*t42 + d17*spb25*t20*t51 + d10*t62*t63 + d11*t62*t63 + d12*t62*t63 + d13*t62*t63 + d8*spa25*spa34*t80 + d20*t80*t81 + d5*t85 + d20*t47*t92 - d9*spa14*t21*t98 + d4*t73*square(spb13) - spa34*t21*t77*T(2); 
complex<T> t8 = spa14*spb13*(-(d23*spa24) + spb13*t49); 
complex<T> t10 = d23*spa14*spa24*spb13 + d6*spa14*spa24*t30 + d3*s35*t61 + d12*t43*t62 + d13*t43*t62 + d9*spa14*t21*t98 - spa14*t49*square(spb13) - d4*t73*square(spb13); 
complex<T> t11 = d33*s45*t23*t63 + d34*spb45*t73 + d35*t92; 
complex<T> t19 = spb34*(d31*t23*t43 + d32*t20*t45); 
complex<T> t59 = -(d21*spb25*t20) + d22*t85; 
complex<T> t60 = d43*s45*spb34*t23*t43 + t20*t31*t48; 
complex<T> t99 = spa14*t69; 
complex<T> t9 = d29*s12*t20*t45 + d28*spb12*t23*t63 + s12*t99 - d27*t110*T(2); 
complex<T> co1 = d30*spb15*t85; 
complex<T> co2 = s23*t99; 
complex<T> co3 = d36*s23*spb12*t42; 
complex<T> co4 = d37*s12*spb15*t85; 
complex<T> co5 = d38*s23*spb34*t42; 
complex<T> co6 = d39*s12*s23*t73; 
complex<T> co7 = t48*t92; 
complex<T> co8 = d41*spb15*t92; 
complex<T> co9 = d42*spb15*t110; 
complex<T> co10 = Complex(0,1); 
SeriesC<T> result = co10*(t2*Int(ep,mu,c12,c345) + t59*Int(ep,mu,c15,c234) + t8*Int(ep,mu,c23,c145) + t1*Int(ep,mu,c34,c125) + t10*Int(ep,mu,c45,c123) + t9*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + co2*Int(ep,mu,c2,c3,c145) + t19*Int(ep,mu,c3,c4,c125) + t11*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c1,c2,c3,c45) + co4*Int(ep,mu,c2,c1,c5,c34) + co5*Int(ep,mu,c2,c3,c4,c15) + co6*Int(ep,mu,c3,c2,c1,c45) + co7*Int(ep,mu,c3,c4,c5,c12) + co8*Int(ep,mu,c4,c5,c1,c23) + co9*Int(ep,mu,c5,c1,c2,c34) + t60*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QpQmqpqmp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, Qm, qp, qm, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QpQmqpqmp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:27:27 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t3 = spa12*spa45; 
complex<T> t6 = square(spb35); 
complex<T> t11 = spa25*spa34; 
complex<T> d3 = cube(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> t16 = t11*t6; 
complex<T> d1 = (s12 - s34)*t3; d1 = T(1)/d1;
complex<T> d2 = t3*square(s12 - s34)*T(3); d2 = T(1)/d2;
complex<T> d4 = spa15*spa34*t3*T(3); d4 = T(1)/d4;
complex<T> t4 = d2*spa24*t16 + d1*spb35*square(spa24) - d3*spb15*t16*T(2) + d4*spa14*square(spa24)*T(2); 
complex<T> t5 = -(d2*spa24*t16) - d1*spb35*square(spa24) + d3*spb15*t16*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t5*Int(ep,mu,c12,c345) + t4*Int(ep,mu,c34,c125));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QpQmqmqpp_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, Qm, qm, qp, p}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QpQmqmqpp PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:27:38 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> t4 = square(s12 - s34); 
complex<T> t5 = spa34*T(2); 
complex<T> t13 = square(spa23); 
complex<T> t14 = -s15 + s45; 
complex<T> t15 = -(s25*spb25); 
complex<T> t16 = spa25*spb45; 
complex<T> t17 = spa35*spb15; 
complex<T> t24 = -(spa12*T(2)); 
complex<T> t37 = spa13*spb15; 
complex<T> t38 = spa24*spb45; 
complex<T> t49 = spa14*spb12; 
complex<T> t55 = spa14*spb34; 
complex<T> d1 = spa12*spa15*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> d3 = (s12 - s34)*(s12 + s15 - s34)*spa15*spa34; d3 = T(1)/d3;
complex<T> d4 = (s12 - s34)*spa12*spa45; d4 = T(1)/d4;
complex<T> d7 = spa45*T(2); d7 = T(1)/d7;
complex<T> d8 = spa15*spa34; d8 = T(1)/d8;
complex<T> d9 = spa12*spa45; d9 = T(1)/d9;
complex<T> d10 = spa15*T(2); d10 = T(1)/d10;
complex<T> d11 = s12 - s34; d11 = T(1)/d11;
complex<T> d12 = spa15*spa34*T(6); d12 = T(1)/d12;
complex<T> d13 = spa12*spa45*T(6); d13 = T(1)/d13;
complex<T> d14 = spa15*spa45*T(3); d14 = T(1)/d14;
complex<T> d15 = cube(s12 - s34); d15 = T(1)/d15;
complex<T> d16 = (s15 - s34)*spa15*spa34; d16 = T(1)/d16;
complex<T> d17 = (s15 - s34)*(s12 + s15 - s34)*spa15*spa34; d17 = T(1)/d17;
complex<T> d18 = spa15*spa34*spa45; d18 = T(1)/d18;
complex<T> d19 = spa35*spa45; d19 = T(1)/d19;
complex<T> d20 = spa15*spa34*square(s12 + s15 - s34); d20 = T(1)/d20;
complex<T> d21 = spa34*square(s12 + s15 - s34); d21 = T(1)/d21;
complex<T> d22 = spa12*spa35*spa45; d22 = T(1)/d22;
complex<T> d23 = spa12*spa15*spa45; d23 = T(1)/d23;
complex<T> d24 = spa15*square(s12 + s15 - s34); d24 = T(1)/d24;
complex<T> d25 = spa12*spa35; d25 = T(1)/d25;
complex<T> d28 = spa12*spa15*spa45*T(2); d28 = T(1)/d28;
complex<T> d29 = spa12*spa15*T(2); d29 = T(1)/d29;
complex<T> d32 = spa12*spa35*T(2); d32 = T(1)/d32;
complex<T> t7 = -t17; 
complex<T> t19 = -t49; 
complex<T> t21 = d14*(s12 + s34); 
complex<T> t26 = -t37; 
complex<T> t27 = -t38; 
complex<T> t28 = -t55; 
complex<T> t29 = -(d10*T(3)); 
complex<T> t31 = d14*spb23; 
complex<T> t35 = spb25*t13; 
complex<T> t39 = d10*T(3); 
complex<T> t40 = -(d1*spa14); 
complex<T> t45 = t13*t15; 
complex<T> t46 = t16*t17; 
complex<T> t48 = -(d7*T(3)); 
complex<T> t50 = d9*spb35; 
complex<T> t54 = d7*T(3); 
complex<T> t69 = s34*t13; 
complex<T> d2 = spa15*spa45*t4*T(6); d2 = T(1)/d2;
complex<T> d5 = spa15*t4*T(2); d5 = T(1)/d5;
complex<T> d6 = spa45*t4*T(2); d6 = T(1)/d6;
complex<T> d26 = spa15*spa45*t5; d26 = T(1)/d26;
complex<T> d27 = t5*square(s12 + s15 - s34); d27 = T(1)/d27;
complex<T> d30 = spa12*t5; d30 = T(1)/d30;
complex<T> d31 = spa45*t5; d31 = T(1)/d31;
complex<T> t23 = -t35; 
complex<T> t25 = d2*t14; 
complex<T> t41 = -t50; 
complex<T> t53 = s25*t35; 
complex<T> t70 = t13*t19; 
complex<T> t71 = t13*t28; 
complex<T> t1 = -(d4*spb35*t13) - spa34*t16*t25 + d5*t16*t26 + d16*t35 + d6*t17*t38 + t13*t40 + d15*spa34*t24*t31*t46 + d17*t53 + d3*t53 + d11*(d8*t23 + d12*spa23*t26 + d13*spa23*t27 + spa23*spb45*t39 + t13*t41 + spa23*spb15*t54) + d15*spa14*t16*t21*t7 + spa12*t25*t7; 
complex<T> t2 = d4*spb35*t13 + spa34*t16*t25 + spa12*t17*t25 + d11*spa23*spb45*t29 + d11*d8*t35 + d11*d12*spa23*t37 + d5*t16*t37 + d11*d13*spa23*t38 + t13*t40 + d3*t45 + d15*spa14*t21*t46 + d11*spa23*spb15*t48 + d15*spa12*t31*t46*t5 + d11*t13*t50 + d6*t38*t7; 
complex<T> t12 = d16*t23 + d17*t45; 
complex<T> t52 = -(d19*spb12*t13) + d20*s12*t53 + d18*t70; 
complex<T> t62 = d24*spb34*t53 - d22*t69 + d23*t71; 
complex<T> co1 = d21*spb15*t45; 
complex<T> co2 = d25*spb45*t13; 
complex<T> co3 = d26*s23*t70; 
complex<T> co4 = d27*s12*spb15*t45; 
complex<T> co5 = d28*s23*t71; 
complex<T> co6 = d29*spb45*t13*t55; 
complex<T> co7 = d30*spa14*spb15*spb45*t13; 
complex<T> co8 = d31*spb15*t13*t49; 
complex<T> co9 = d32*spb45*t69; 
complex<T> co10 = Complex(0,1); 
SeriesC<T> result = co10*(t2*Int(ep,mu,c12,c345) + t12*Int(ep,mu,c15,c234) + t1*Int(ep,mu,c34,c125) + t52*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + t62*Int(ep,mu,c3,c4,c125) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c1,c2,c3,c45) + co4*Int(ep,mu,c2,c1,c5,c34) + co5*Int(ep,mu,c2,c3,c4,c15) + co6*Int(ep,mu,c3,c4,c5,c12) + co7*Int(ep,mu,c4,c5,c1,c23) + co8*Int(ep,mu,c5,c1,c2,c34) + co9*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QpQmqmqpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, Qm, qm, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QpQmqmqpp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:27:40 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> s15 = -(spa15*spb15);
complex<T> s45 = -(spa45*spb45);
complex<T> t5 = square(spa23); 
complex<T> t9 = spa12*T(2); 
complex<T> t11 = -s15 + s45; 
complex<T> t12 = spa25*spb45; 
complex<T> t13 = spa35*spb15; 
complex<T> t14 = -(spa12*spa45*spb25) - spa15*spa34*spb35; 
complex<T> t18 = -(spa12*T(2)); 
complex<T> d1 = spa12*spa15*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> d2 = spa15*spa45*square(s12 - s34)*T(6); d2 = T(1)/d2;
complex<T> d4 = spa15*spa34*T(6); d4 = T(1)/d4;
complex<T> d5 = spa12*spa45*T(6); d5 = T(1)/d5;
complex<T> d6 = s12 - s34; d6 = T(1)/d6;
complex<T> d7 = spa15*spa45*T(3); d7 = T(1)/d7;
complex<T> d8 = cube(s12 - s34); d8 = T(1)/d8;
complex<T> t7 = -t13; 
complex<T> t15 = d7*(s12 + s34); 
complex<T> t17 = spa34*t12; 
complex<T> t22 = d1*spa14; 
complex<T> t24 = d7*spb23; 
complex<T> d3 = (s12 - s34)*spa15*spa34*spa45*t9; d3 = T(1)/d3;
complex<T> t33 = spa14*t15; 
complex<T> t40 = t13*t17; 
complex<T> t1 = -(d4*d6*spa13*spa23*spb15) - d5*d6*spa23*spa24*spb45 - d2*t11*t17 + d8*t18*t24*t40 + d3*t14*t5 + t22*t5 + d2*spa12*t11*t7 + d8*t12*t33*t7; 
complex<T> t2 = d4*d6*spa13*spa23*spb15 + d5*d6*spa23*spa24*spb45 + d2*spa12*t11*t13 + d2*t11*t17 + d8*t12*t13*t33 - d3*t14*t5 + t22*t5 + d8*t24*t40*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*Int(ep,mu,c12,c345) + t2*Int(ep,mu,c34,c125));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmpQpqpqm_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, p, Qp, qp, qm}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmpQpqpqm PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:27:46 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> t6 = -(spa15*T(3)); 
complex<T> t9 = square(spa15); 
complex<T> t10 = s12*T(3) - s45*T(4); 
complex<T> t13 = spb12*spb23; 
complex<T> t15 = spa13*spb34; 
complex<T> t18 = -(s12*T(3)); 
complex<T> t22 = spa15*spb24; 
complex<T> d1 = spa12*spa23*spa45*T(4); d1 = T(1)/d1;
complex<T> d2 = s45*spa23*T(4); d2 = T(1)/d2;
complex<T> d3 = s45*(-s12 + s45)*spa23*T(2); d3 = T(1)/d3;
complex<T> d4 = spa23*square(s12 - s45)*T(2); d4 = T(1)/d4;
complex<T> d5 = s45*spa12*spa23*T(4); d5 = T(1)/d5;
complex<T> d6 = s45*spa23*square(s12 - s45)*T(2); d6 = T(1)/d6;
complex<T> d7 = (s12 - s45)*spa12*spa23*spa45*T(2); d7 = T(1)/d7;
complex<T> d8 = spa12*spa23*spa45*T(3); d8 = T(1)/d8;
complex<T> d9 = s45*spa23*T(2); d9 = T(1)/d9;
complex<T> d10 = s45*spa12*spa23*T(2); d10 = T(1)/d10;
complex<T> d11 = spa23*spa45; d11 = T(1)/d11;
complex<T> d12 = s45*spa23; d12 = T(1)/d12;
complex<T> d13 = spa12*spa45; d13 = T(1)/d13;
complex<T> d14 = s45; d14 = T(1)/d14;
complex<T> d15 = s45*spa12; d15 = T(1)/d15;
complex<T> d16 = spa23; d16 = T(1)/d16;
complex<T> d17 = spa12*spa23; d17 = T(1)/d17;
complex<T> d18 = spa45*T(2); d18 = T(1)/d18;
complex<T> d19 = s45*T(2); d19 = T(1)/d19;
complex<T> d20 = spa12*spa45*T(2); d20 = T(1)/d20;
complex<T> d21 = spa12*spa23*T(2); d21 = T(1)/d21;
complex<T> d22 = spa23*spa45*T(2); d22 = T(1)/d22;
complex<T> t5 = -t15; 
complex<T> t12 = -(d1*T(3)); 
complex<T> t17 = -t22; 
complex<T> t20 = d6*t10; 
complex<T> t29 = d4*spb23; 
complex<T> t32 = d2*spb24; 
complex<T> t37 = d5*t15; 
complex<T> t40 = -(spb45*t9); 
complex<T> t44 = spa15*t15; 
complex<T> t49 = -(spb12*t9); 
complex<T> t50 = -(spb23*t9); 
complex<T> t3 = spa12*spa35*spb24*t29 + spb12*t20*t44 - d7*spa12*spa35*spb23*t6 + d7*spa45*t5*t6 - d8*t9*T(2) + d9*t22*T(3) + d10*t44*T(3) + d3*t17*(t18 + s45*T(4)); 
complex<T> t28 = spa15*t5; 
complex<T> t43 = t32*t6 + t37*t6 + t12*t9; 
complex<T> t45 = d19*s12*spb23*t17 + d19*t13*t44 + d18*t13*t9; 
complex<T> t4 = spb12*t20*t28 - spa12*spa35*spb24*t29 + t32*t6 + t37*t6 + d7*spa15*(-(spa12*spa35*spb23) + spa45*t5)*T(3) + d1*t9*T(3) + d3*t22*(t18 + s45*T(4)); 
complex<T> t26 = d12*s12*t22 + d12*spb12*t28 + d11*t49; 
complex<T> t34 = d14*spb23*t17 + d15*spb23*t28 + d13*t50; 
complex<T> t39 = d16*t17 + d17*(t28 + t40); 
complex<T> co1 = d20*s34*t50; 
complex<T> co2 = d21*s34*t40; 
complex<T> co3 = d21*s15*t40; 
complex<T> co4 = d22*s15*t49; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t4*Int(ep,mu,c12,c345) + t43*Int(ep,mu,c23,c145) + t3*Int(ep,mu,c45,c123) + t26*Int(ep,mu,c1,c2,c345) + t34*Int(ep,mu,c2,c3,c145) + t39*Int(ep,mu,c4,c5,c123) + t45*Int(ep,mu,c1,c2,c3,c45) + co1*Int(ep,mu,c2,c3,c4,c15) + co2*Int(ep,mu,c3,c4,c5,c12) + co3*Int(ep,mu,c4,c5,c1,c23) + co4*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmpQpqmqp_PentP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, p, Qp, qm, qp}, PentP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmpQpqmqp PentP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:18 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s35 = -(spa35*spb35);
complex<T> t4 = spa23*T(2); 
complex<T> t6 = spa12*spa15; 
complex<T> t13 = square(spa14); 
complex<T> t14 = -(spa13*spb35); 
complex<T> t15 = -(spa35*T(3)); 
complex<T> t16 = spa15*spa34; 
complex<T> t19 = (s15 - s34)*spa23; 
complex<T> t28 = square(spb25); 
complex<T> t29 = square(spa45); 
complex<T> t30 = square(spa12); 
complex<T> t31 = square(spa15); 
complex<T> t32 = square(spa24); 
complex<T> t33 = spa14*spa35; 
complex<T> t35 = square(s25); 
complex<T> t36 = spa14*spa35 + spa13*spa45; 
complex<T> t37 = s12 - s34; 
complex<T> t39 = -(s12*T(3)); 
complex<T> t40 = s12*T(3) - s45*T(4); 
complex<T> t41 = square(spb35); 
complex<T> t48 = cube(spb25); 
complex<T> t54 = spa23*spa25; 
complex<T> t68 = cube(spa14); 
complex<T> t69 = spa13*spb12; 
complex<T> t75 = spb23*spb34; 
complex<T> t78 = spa13*spa14; 
complex<T> t83 = spb35*T(3); 
complex<T> t93 = spa14*spb25; 
complex<T> t95 = s45*spa13; 
complex<T> t97 = spa34*spb35; 
complex<T> t101 = spa35*spb12; 
complex<T> t114 = spa35*spb45; 
complex<T> d2 = spa12*spa23*spa45*T(4); d2 = T(1)/d2;
complex<T> d4 = s45*spa23*T(4); d4 = T(1)/d4;
complex<T> d14 = (s12 - s45)*spa12*spa23; d14 = T(1)/d14;
complex<T> d15 = s45*spa12*spa23*T(4); d15 = T(1)/d15;
complex<T> d18 = s35*(s12 - s45)*spa12*spa23; d18 = T(1)/d18;
complex<T> d25 = (s15 - s34)*spa12*spa34; d25 = T(1)/d25;
complex<T> d26 = spa15*spa23*square(s15 - s34); d26 = T(1)/d26;
complex<T> d34 = spa12*spa34*T(2); d34 = T(1)/d34;
complex<T> d35 = s15 - s34; d35 = T(1)/d35;
complex<T> d37 = spa34*T(2); d37 = T(1)/d37;
complex<T> d38 = square(s15 - s34); d38 = T(1)/d38;
complex<T> d39 = spa12*spa23*spa45*T(6); d39 = T(1)/d39;
complex<T> d43 = spa23*spa35; d43 = T(1)/d43;
complex<T> d44 = spa23*spa45; d44 = T(1)/d44;
complex<T> d45 = spa23*square(spa35); d45 = T(1)/d45;
complex<T> d46 = s45*spa23; d46 = T(1)/d46;
complex<T> d53 = spa12*spa45; d53 = T(1)/d53;
complex<T> d54 = s45; d54 = T(1)/d54;
complex<T> d55 = s45*spa12; d55 = T(1)/d55;
complex<T> d56 = spa12*spa23*spa35; d56 = T(1)/d56;
complex<T> d57 = spa12*spa23*square(spa35); d57 = T(1)/d57;
complex<T> d62 = spa23; d62 = T(1)/d62;
complex<T> d63 = spa12*spa23; d63 = T(1)/d63;
complex<T> d66 = s45*T(2); d66 = T(1)/d66;
complex<T> d70 = spa15*spa24*T(2); d70 = T(1)/d70;
complex<T> t5 = spa34*t37; 
complex<T> t12 = -(d63*spb45*t13) + d63*spb35*t78 + d62*t93 + d57*t36*t95 - d56*s45*t78*T(3); 
complex<T> t17 = square(t37); 
complex<T> t34 = s15 + t37; 
complex<T> t38 = -t69; 
complex<T> t43 = -(spa14*spb15) + t97; 
complex<T> t46 = -(s25*spa14) - spa15*t97*T(2); 
complex<T> t58 = d70*s34; 
complex<T> t71 = t39 + s45*T(4); 
complex<T> t79 = t28*t30; 
complex<T> t81 = spa35*t29; 
complex<T> t94 = t31*t32; 
complex<T> t99 = -(d2*T(3)); 
complex<T> t106 = spa24*t28; 
complex<T> t107 = -(spa13*t13); 
complex<T> t110 = d4*T(3); 
complex<T> t113 = d66*spb23; 
complex<T> t127 = spa14*t15; 
complex<T> t144 = spb23*(-(d53*t13) + d55*spb35*t78 + d54*t93); 
complex<T> t147 = spb12*t14; 
complex<T> d3 = spa34*spa45*t4*t6; d3 = T(1)/d3;
complex<T> d5 = s45*(-s12 + s45)*t4; d5 = T(1)/d5;
complex<T> d6 = t4*square(s12 - s45); d6 = T(1)/d6;
complex<T> d13 = spa12*spa23*t37; d13 = T(1)/d13;
complex<T> d16 = s45*t4*square(s12 - s45); d16 = T(1)/d16;
complex<T> d17 = s35*spa12*spa23*t37; d17 = T(1)/d17;
complex<T> d22 = spa23*spa34*t6; d22 = T(1)/d22;
complex<T> d24 = spa15*t19; d24 = T(1)/d24;
complex<T> d27 = t4*square(s15 - s34); d27 = T(1)/d27;
complex<T> d29 = spa25*t16*t4*square(s15 - s34); d29 = T(1)/d29;
complex<T> d33 = spa15*t4; d33 = T(1)/d33;
complex<T> d36 = t16*t4; d36 = T(1)/d36;
complex<T> d40 = s45*t4; d40 = T(1)/d40;
complex<T> d41 = s45*spa12*t4; d41 = T(1)/d41;
complex<T> d64 = t16*T(2); d64 = T(1)/d64;
complex<T> d65 = spa45*t16*T(2); d65 = T(1)/d65;
complex<T> d71 = t6*T(2); d71 = T(1)/d71;
complex<T> d72 = spa45*t6*T(2); d72 = T(1)/d72;
complex<T> d73 = spa12*spa35*t4; d73 = T(1)/d73;
complex<T> d74 = t4*t6; d74 = T(1)/d74;
complex<T> d75 = spa12*spa34*t4; d75 = T(1)/d75;
complex<T> d76 = spa34*t4; d76 = T(1)/d76;
complex<T> d77 = spa34*spa45*t4; d77 = T(1)/d77;
complex<T> d78 = spa12*t4*square(spa35); d78 = T(1)/d78;
complex<T> t3 = square(t34); 
complex<T> t10 = -(d6*spa12*spa34*spb23*spb25) + d16*spa14*t147*t40 - d18*spa13*t36*t41 + d5*t71*t93 - d41*spb35*t78*T(3) - d40*t93*T(3) - d14*spb35*t78*T(4) - d39*t13*T(13); 
complex<T> t44 = d75*spb15; 
complex<T> t45 = d74*spb34; 
complex<T> t50 = -t81; 
complex<T> t57 = d36*s12; 
complex<T> t60 = d33*spb24; 
complex<T> t73 = d3*T(3); 
complex<T> t80 = -t94; 
complex<T> t124 = spb23*t58; 
complex<T> t133 = t13*t38; 
complex<T> t162 = d15*t78*t83 + t110*t93 + t13*t99; 
complex<T> d1 = spa23*t34*t5*t6; d1 = T(1)/d1;
complex<T> d7 = spa23*t34*t5; d7 = T(1)/d7;
complex<T> d8 = t16*t17*t4; d8 = T(1)/d8;
complex<T> d9 = spa25*t16*t17*t4; d9 = T(1)/d9;
complex<T> d10 = spa15*t34*t5*t54; d10 = T(1)/d10;
complex<T> d12 = spa12*spa34*t17*t34*T(2); d12 = T(1)/d12;
complex<T> d19 = spa12*spa23*t5; d19 = T(1)/d19;
complex<T> d20 = spa12*spa34*t17*T(2); d20 = T(1)/d20;
complex<T> d21 = spa45*t4*t5*t6; d21 = T(1)/d21;
complex<T> d23 = spa34*t19*t34*t6; d23 = T(1)/d23;
complex<T> d28 = spa34*t19*t34; d28 = T(1)/d28;
complex<T> d30 = spa25*t16*t19*t34; d30 = T(1)/d30;
complex<T> d32 = spa12*spa34*t34*square(s15 - s34)*T(2); d32 = T(1)/d32;
complex<T> d49 = spa34*cube(t34); d49 = T(1)/d49;
complex<T> d52 = spa12*spa34*cube(t34); d52 = T(1)/d52;
complex<T> d61 = spa12*cube(t34); d61 = T(1)/d61;
complex<T> d69 = spa34*cube(t34)*T(2); d69 = T(1)/d69;
complex<T> t125 = spb15*(d76*t133 + d77*t101*t68); 
complex<T> t137 = t124*t13 + d71*t107*t75 + d72*spa35*t68*t75; 
complex<T> t148 = d64*spb23*t133 + spa14*t113*t147 + d65*spb23*t101*t68 + s12*t113*t93; 
complex<T> t153 = t44*(t114*t68 + t13*t95); 
complex<T> t158 = t114*t45*t68 + t13*t45*t95 - d73*s34*s45*t78*T(3); 
complex<T> d11 = spa12*t3*t5; d11 = T(1)/d11;
complex<T> d31 = (s15 - s34)*spa12*spa34*t3; d31 = T(1)/d31;
complex<T> d42 = spa23*t16*t3; d42 = T(1)/d42;
complex<T> d47 = spa23*spa34*t3; d47 = T(1)/d47;
complex<T> d48 = t16*t3*t54; d48 = T(1)/d48;
complex<T> d50 = spa12*spa23*spa34*t3; d50 = T(1)/d50;
complex<T> d51 = spa34*t3*t54; d51 = T(1)/d51;
complex<T> d58 = spa23*t3*t6; d58 = T(1)/d58;
complex<T> d59 = spa23*t3; d59 = T(1)/d59;
complex<T> d60 = spa15*t3*t54; d60 = T(1)/d60;
complex<T> d67 = spa34*t3*t4; d67 = T(1)/d67;
complex<T> d68 = spa25*spa34*t3*t4; d68 = T(1)/d68;
complex<T> t1 = -(d27*spa15*t106) + d25*spb25*t13 + d34*d35*spb25*t13 + d37*d38*spa14*spa45*t28 - d26*spa12*spb23*spb25*t33 + d23*t107*t35 + d38*t13*t57*t69 + d24*spb23*t78 + d29*t79*t81 + d30*t79*t81 + d31*t48*t94 + d32*t48*t94 - d22*spa13*t13*T(2) + d28*t106*t33*T(3) - d35*t13*t60*T(3); 
complex<T> t2 = d27*spa15*t106 + d22*t107 + d28*t106*t127 + d7*t106*t127 - d25*spb25*t13 - d34*d35*spb25*t13 - d37*d38*spa14*spa45*t28 + d8*spa12*spa13*t28*t29 + d26*spa12*spb23*spb25*t33 + d1*spa13*t13*t35 + d23*spa13*t13*t35 - d17*spa13*t36*t41 + d38*t133*t57 + spa35*t68*t73 - d24*spb23*t78 - d19*t43*t78 + d10*t50*t79 + d29*t50*t79 + d30*t50*t79 + d9*t50*t79 + d11*t48*t80 + d12*t48*t80 + d31*t48*t80 + d32*t48*t80 + d20*t46*t93 + d21*spa12*spa35*spa45*spb25*t13*T(3) + d21*spa35*spb35*t13*t16*T(3) + d35*t13*t60*T(3) - d13*spb35*t78*T(3); 
complex<T> t8 = d59*spb34*t106*t127 + d58*spa13*spb34*t13*t35 + d57*s34*spa13*t36 + d60*spb34*t50*t79 + d61*spb34*t48*t80 - d56*s34*t78*T(3); 
complex<T> t9 = d47*s12*t106*t127 - d44*spb12*t13 + d42*t133*t35 + d46*spa14*spb35*t69 + d45*t36*t69 + d48*s12*t50*t79 - d46*s12*t93 + d49*spb12*t48*t94 - d43*spa14*t69*T(3); 
complex<T> t11 = d6*spa12*spa34*spb23*spb25 + d21*t13*t15*(spa12*spa45*spb25 + spb35*t16) + t162 - d8*spa12*spa13*t28*t29 + d1*t107*t35 + d17*spa13*t36*t41 + d18*spa13*t36*t41 + d16*spa14*spb35*t40*t69 + spa35*t68*t73 + d19*t43*t78 + d10*t79*t81 + d9*t79*t81 + d13*t78*t83 - d20*t46*t93 - d5*t71*t93 + d11*t48*t94 + d12*t48*t94 + d7*t106*t33*T(3) + d14*spb35*t78*T(4); 
complex<T> t23 = d52*s15*t48*t80 + d51*spb15*t79*t81 - d50*spa13*spb15*t35*square(spa14) - d47*s15*t106*t33*T(3); 
complex<T> t24 = d67*s15*t106*t33*t39 + d68*s12*spb15*t79*t81 + d69*s15*spb12*t48*t94 + d67*spb15*t35*t69*square(spa14); 
complex<T> co1 = -(t124*t13); 
complex<T> co2 = d78*s34*t36*t95; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t11*Int(ep,mu,c12,c345) + t1*Int(ep,mu,c15,c234) + t162*Int(ep,mu,c23,c145) + t2*Int(ep,mu,c34,c125) + t10*Int(ep,mu,c45,c123) + t9*Int(ep,mu,c1,c2,c345) + t23*Int(ep,mu,c1,c5,c234) + t144*Int(ep,mu,c2,c3,c145) + t8*Int(ep,mu,c3,c4,c125) + t12*Int(ep,mu,c4,c5,c123) + t148*Int(ep,mu,c1,c2,c3,c45) + t24*Int(ep,mu,c2,c1,c5,c34) + t137*Int(ep,mu,c2,c3,c4,c15) + t158*Int(ep,mu,c3,c4,c5,c12) + co1*Int(ep,mu,c4,c3,c2,c15) + t153*Int(ep,mu,c4,c5,c1,c23) + t125*Int(ep,mu,c5,c1,c2,c34) + co2*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QppQmqmqp_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, p, Qm, qm, qp}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QppQmqmqp BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:19 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t4 = square(spa34); 
complex<T> t7 = -(s23*T(3)) + s45*T(4); 
complex<T> t9 = spa12*T(2); 
complex<T> t12 = spa14*spa23; 
complex<T> d4 = spa12*spa23; d4 = T(1)/d4;
complex<T> d5 = spa45*T(2); d5 = T(1)/d5;
complex<T> t11 = spb23*t4; 
complex<T> d1 = spa45*t9*square(s23 - s45); d1 = T(1)/d1;
complex<T> d2 = t9*square(s23 - s45); d2 = T(1)/d2;
complex<T> d3 = spa23*spa45*t9; d3 = T(1)/d3;
complex<T> t2 = -(d2*spb12*spb25*t12) + d2*spa34*spb25*t7 - d1*t11*t7 + d3*t4*T(3); 
complex<T> t3 = d1*t11*t7 + d2*spb25*(spb12*t12 - spa34*t7); 
complex<T> co1 = d4*spb45*t4; 
complex<T> co2 = d5*spb12*t11; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*Int(ep,mu,c23,c145) + t2*Int(ep,mu,c45,c123) + co1*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmQppqmqp_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, Qp, p, qm, qp}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmQppqmqp BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:20 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t6 = square(spa14); 
complex<T> t9 = square(spb35); 
complex<T> d1 = spa12*spa23*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa23*spa45; d2 = T(1)/d2;
complex<T> d3 = spa23*square(s12 - s45)*T(2); d3 = T(1)/d3;
complex<T> d4 = spa13*spa23*spa45; d4 = T(1)/d4;
complex<T> d5 = spa23*spa34*spa45; d5 = T(1)/d5;
complex<T> d6 = spa13*spa45; d6 = T(1)/d6;
complex<T> d7 = spa13*spa23; d7 = T(1)/d7;
complex<T> d8 = spa13*spa45*T(2); d8 = T(1)/d8;
complex<T> t4 = -(d2*spa14*spa24*spb23) - d2*spb13*t6 - d3*spa13*spa45*t9 + d1*spa24*t6*T(3); 
complex<T> t5 = d2*spa14*spa24*spb23 + d2*spb13*t6 + d3*spa13*spa45*t9; 
complex<T> t17 = s12*t6; 
complex<T> t3 = d4*t17 + d5*spa24*spb12*t6; 
complex<T> co1 = -(d6*spb23*t6); 
complex<T> co2 = d7*spb45*t6; 
complex<T> co3 = -(d8*spb23*t17); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t4*Int(ep,mu,c12,c345) + t5*Int(ep,mu,c45,c123) + t3*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c2,c3,c145) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmQppqpqm_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, Qp, p, qp, qm}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmQppqpqm BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:21 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t5 = square(spa15); 
complex<T> t8 = square(spb34); 
complex<T> d1 = spa12*spa23*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa23; d2 = T(1)/d2;
complex<T> d3 = spa23*square(s12 - s45)*T(2); d3 = T(1)/d3;
complex<T> d4 = spa13*spa23*spa45; d4 = T(1)/d4;
complex<T> d5 = spa23*spa34*spa45; d5 = T(1)/d5;
complex<T> d6 = spa13*spa45; d6 = T(1)/d6;
complex<T> d7 = spa13*spa23; d7 = T(1)/d7;
complex<T> d8 = spa13*spa45*T(2); d8 = T(1)/d8;
complex<T> t3 = -(d2*spa15*spb34) - d3*spa13*spa45*t8 + d1*spa24*t5*T(3); 
complex<T> t4 = d2*spa15*spb34 + d3*spa13*spa45*t8; 
complex<T> t13 = s12*t5; 
complex<T> t2 = d4*t13 + d5*spa24*spb12*t5; 
complex<T> co1 = -(d6*spb23*t5); 
complex<T> co2 = d7*spb45*t5; 
complex<T> co3 = -(d8*spb23*t13); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t3*Int(ep,mu,c12,c345) + t4*Int(ep,mu,c45,c123) + t2*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c2,c3,c145) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QpQmpqmqp_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, Qm, p, qm, qp}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QpQmpqmqp BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:22 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> t2 = spa45*T(2); 
complex<T> t5 = square(s13); 
complex<T> t6 = spa23*spa45; 
complex<T> t9 = square(spa12); 
complex<T> t10 = square(spa34); 
complex<T> t11 = cube(spb13); 
complex<T> t19 = -(s13*spa24); 
complex<T> d9 = spa45*cube(spa13); d9 = T(1)/d9;
complex<T> d10 = spa23*cube(spa13); d10 = T(1)/d10;
complex<T> t16 = t19 + spa12*spa45*spb15*T(2); 
complex<T> t18 = spa23*t2; 
complex<T> t23 = t10*t9; 
complex<T> d3 = (s12 - s45)*t5*t6; d3 = T(1)/d3;
complex<T> d5 = (s23 - s45)*t5*t6; d5 = T(1)/d5;
complex<T> d7 = t6*cube(spa13); d7 = T(1)/d7;
complex<T> d8 = spa34*t6; d8 = T(1)/d8;
complex<T> d11 = t2*cube(spa13); d11 = T(1)/d11;
complex<T> t7 = d7*s12*t23 + d8*spb12*cube(spa24); 
complex<T> t17 = -t23; 
complex<T> d1 = spa12*spa34*t18; d1 = T(1)/d1;
complex<T> d2 = s13*t18*square(s12 - s45); d2 = T(1)/d2;
complex<T> d4 = s13*t18*square(s23 - s45); d4 = T(1)/d4;
complex<T> d6 = t18*square(s23 - s45); d6 = T(1)/d6;
complex<T> t4 = -(d6*spa24*spb13*t16) + t11*(d2*t17 + d4*t17 + (d3 + d5)*t23); 
complex<T> t8 = d3*t11*t17 + d2*t11*t23 + d1*cube(spa24)*T(3); 
complex<T> t21 = d5*t11*t17 + d6*spa24*spb13*(t19 + spa12*spb15*t2) + d4*t11*t23; 
complex<T> co1 = d9*spb23*t17; 
complex<T> co2 = d10*spb45*t23; 
complex<T> co3 = d11*s12*spb23*t17; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t8*Int(ep,mu,c12,c345) + t21*Int(ep,mu,c23,c145) + t4*Int(ep,mu,c45,c123) + t7*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c2,c3,c145) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QpQmpqpqm_BoxP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, Qm, p, qp, qm}, BoxP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QpQmpqpqm BoxP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:24 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> t2 = spa45*T(2); 
complex<T> t5 = square(s13); 
complex<T> t6 = spa23*spa45; 
complex<T> t10 = square(spa12); 
complex<T> t11 = square(spa35); 
complex<T> t12 = cube(spb13); 
complex<T> t16 = square(spa25); 
complex<T> t17 = -(s13*spa25) - spa12*spa45*spb14*T(2); 
complex<T> t20 = spa25*spb13; 
complex<T> d9 = spa45*cube(spa13); d9 = T(1)/d9;
complex<T> d10 = spa23*cube(spa13); d10 = T(1)/d10;
complex<T> t19 = spa23*t2; 
complex<T> t23 = t10*t11; 
complex<T> d3 = (s12 - s45)*t5*t6; d3 = T(1)/d3;
complex<T> d5 = (s23 - s45)*t5*t6; d5 = T(1)/d5;
complex<T> d7 = t6*cube(spa13); d7 = T(1)/d7;
complex<T> d8 = spa34*t6; d8 = T(1)/d8;
complex<T> d11 = t2*cube(spa13); d11 = T(1)/d11;
complex<T> t7 = d8*spa24*spb12*t16 + d7*s12*t23; 
complex<T> t18 = -t23; 
complex<T> d1 = spa12*spa34*t19; d1 = T(1)/d1;
complex<T> d2 = s13*t19*square(s12 - s45); d2 = T(1)/d2;
complex<T> d4 = s13*t19*square(s23 - s45); d4 = T(1)/d4;
complex<T> d6 = t19*square(s23 - s45); d6 = T(1)/d6;
complex<T> t8 = d3*t12*t18 + d2*t12*t23 + d1*spa24*t16*T(3); 
complex<T> t21 = d6*t17; 
complex<T> t4 = d2*t12*t18 + d4*t12*t18 - t20*t21 + d3*t12*t23 + d5*t12*t23; 
complex<T> t22 = d5*t12*t18 + t20*t21 + d4*t12*t23; 
complex<T> co1 = d9*spb23*t18; 
complex<T> co2 = d10*spb45*t23; 
complex<T> co3 = d11*s12*spb23*t18; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t8*Int(ep,mu,c12,c345) + t22*Int(ep,mu,c23,c145) + t4*Int(ep,mu,c45,c123) + t7*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c2,c3,c145) + co2*Int(ep,mu,c4,c5,c123) + co3*Int(ep,mu,c3,c2,c1,c45));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmQpqmpqp_TriP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, Qp, qm, p, qp}, TriP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmQpqmpqp TriP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:25 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = spa34*spa45; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb12*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c12,c345) + co2*Int(ep,mu,c1,c2,c345));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmQpqppqm_TriP
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, Qp, qp, p, qm}, TriP}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmQpqppqm TriP");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:25 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spa15); 
complex<T> d1 = spa12*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = spa34*spa45; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb12*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c12,c345) + co2*Int(ep,mu,c1,c2,c345));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmQpqmpqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, Qp, qm, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmQpqmpqp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:25 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QmQpqppqm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qm, Qp, qp, p, qm}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QmQpqppqm nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:26 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa15); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QpQmqmpqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, Qm, qm, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QpQmqmpqp nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:26 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa23); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345);  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q2Q1g_QpQmqppqm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{Qp, Qm, qp, p, qm}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  QpQmqppqm nf");
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
 // #define TimeStamp "Wed 26 Nov 2008 10:30:27 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa25); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*Int(ep,mu,c12,c345);  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qmqpQpQmp_PentP C2q2Q1g_4945_PentP
#define _C_qmqpQpQmp_nf C2q2Q1g_4945_nf
#define _C_qpqmQpQmp_PentP C2q2Q1g_4940_PentP
#define _C_qpqmQpQmp_nf C2q2Q1g_4940_nf
#define _C_qpqmQmQpp_PentP C2q2Q1g_5120_PentP
#define _C_qpqmQmQpp_nf C2q2Q1g_5120_nf
#define _C_qmpqpQpQm_PentP C2q2Q1g_6355_PentP
#define _C_qmpqpQmQp_PentP C2q2Q1g_7435_PentP
#define _C_qppqmQmQp_BoxP C2q2Q1g_7400_BoxP
#define _C_qmqppQmQp_BoxP C2q2Q1g_7465_BoxP
#define _C_qmqppQpQm_BoxP C2q2Q1g_6385_BoxP
#define _C_qpqmpQmQp_BoxP C2q2Q1g_7460_BoxP
#define _C_qpqmpQpQm_BoxP C2q2Q1g_6380_BoxP
#define _C_qmqpQmpQp_TriP C2q2Q1g_7285_TriP
#define _C_qmqpQppQm_TriP C2q2Q1g_6025_TriP
#define _C_qmqpQmpQp_nf C2q2Q1g_7285_nf
#define _C_qmqpQppQm_nf C2q2Q1g_6025_nf
#define _C_qpqmQmpQp_nf C2q2Q1g_7280_nf
#define _C_qpqmQppQm_nf C2q2Q1g_6020_nf
#define _C_QmQpqpqmp_PentP C2q2Q1g_4210_PentP
#define _C_QmQpqpqmp_nf C2q2Q1g_4210_nf
#define _C_QpQmqpqmp_PentP C2q2Q1g_4205_PentP
#define _C_QpQmqpqmp_nf C2q2Q1g_4205_nf
#define _C_QpQmqmqpp_PentP C2q2Q1g_4385_PentP
#define _C_QpQmqmqpp_nf C2q2Q1g_4385_nf
#define _C_QmpQpqpqm_PentP C2q2Q1g_1930_PentP
#define _C_QmpQpqmqp_PentP C2q2Q1g_3010_PentP
#define _C_QppQmqmqp_BoxP C2q2Q1g_2975_BoxP
#define _C_QmQppqmqp_BoxP C2q2Q1g_2950_BoxP
#define _C_QmQppqpqm_BoxP C2q2Q1g_1870_BoxP
#define _C_QpQmpqmqp_BoxP C2q2Q1g_2945_BoxP
#define _C_QpQmpqpqm_BoxP C2q2Q1g_1865_BoxP
#define _C_QmQpqmpqp_TriP C2q2Q1g_3310_TriP
#define _C_QmQpqppqm_TriP C2q2Q1g_2050_TriP
#define _C_QmQpqmpqp_nf C2q2Q1g_3310_nf
#define _C_QmQpqppqm_nf C2q2Q1g_2050_nf
#define _C_QpQmqmpqp_nf C2q2Q1g_3305_nf
#define _C_QpQmqppqm_nf C2q2Q1g_2045_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmqpQpQmp_PentP case 4945 : \
          return &C2q2Q1g_4945_PentP
#define _CASE_qmqpQpQmp_nf case 4945 : \
          return &C2q2Q1g_4945_nf
#define _CASE_qpqmQpQmp_PentP case 4940 : \
          return &C2q2Q1g_4940_PentP
#define _CASE_qpqmQpQmp_nf case 4940 : \
          return &C2q2Q1g_4940_nf
#define _CASE_qpqmQmQpp_PentP case 5120 : \
          return &C2q2Q1g_5120_PentP
#define _CASE_qpqmQmQpp_nf case 5120 : \
          return &C2q2Q1g_5120_nf
#define _CASE_qmpqpQpQm_PentP case 6355 : \
          return &C2q2Q1g_6355_PentP
#define _CASE_qmpqpQmQp_PentP case 7435 : \
          return &C2q2Q1g_7435_PentP
#define _CASE_qppqmQmQp_BoxP case 7400 : \
          return &C2q2Q1g_7400_BoxP
#define _CASE_qmqppQmQp_BoxP case 7465 : \
          return &C2q2Q1g_7465_BoxP
#define _CASE_qmqppQpQm_BoxP case 6385 : \
          return &C2q2Q1g_6385_BoxP
#define _CASE_qpqmpQmQp_BoxP case 7460 : \
          return &C2q2Q1g_7460_BoxP
#define _CASE_qpqmpQpQm_BoxP case 6380 : \
          return &C2q2Q1g_6380_BoxP
#define _CASE_qmqpQmpQp_TriP case 7285 : \
          return &C2q2Q1g_7285_TriP
#define _CASE_qmqpQppQm_TriP case 6025 : \
          return &C2q2Q1g_6025_TriP
#define _CASE_qmqpQmpQp_nf case 7285 : \
          return &C2q2Q1g_7285_nf
#define _CASE_qmqpQppQm_nf case 6025 : \
          return &C2q2Q1g_6025_nf
#define _CASE_qpqmQmpQp_nf case 7280 : \
          return &C2q2Q1g_7280_nf
#define _CASE_qpqmQppQm_nf case 6020 : \
          return &C2q2Q1g_6020_nf
#define _CASE_QmQpqpqmp_PentP case 4210 : \
          return &C2q2Q1g_4210_PentP
#define _CASE_QmQpqpqmp_nf case 4210 : \
          return &C2q2Q1g_4210_nf
#define _CASE_QpQmqpqmp_PentP case 4205 : \
          return &C2q2Q1g_4205_PentP
#define _CASE_QpQmqpqmp_nf case 4205 : \
          return &C2q2Q1g_4205_nf
#define _CASE_QpQmqmqpp_PentP case 4385 : \
          return &C2q2Q1g_4385_PentP
#define _CASE_QpQmqmqpp_nf case 4385 : \
          return &C2q2Q1g_4385_nf
#define _CASE_QmpQpqpqm_PentP case 1930 : \
          return &C2q2Q1g_1930_PentP
#define _CASE_QmpQpqmqp_PentP case 3010 : \
          return &C2q2Q1g_3010_PentP
#define _CASE_QppQmqmqp_BoxP case 2975 : \
          return &C2q2Q1g_2975_BoxP
#define _CASE_QmQppqmqp_BoxP case 2950 : \
          return &C2q2Q1g_2950_BoxP
#define _CASE_QmQppqpqm_BoxP case 1870 : \
          return &C2q2Q1g_1870_BoxP
#define _CASE_QpQmpqmqp_BoxP case 2945 : \
          return &C2q2Q1g_2945_BoxP
#define _CASE_QpQmpqpqm_BoxP case 1865 : \
          return &C2q2Q1g_1865_BoxP
#define _CASE_QmQpqmpqp_TriP case 3310 : \
          return &C2q2Q1g_3310_TriP
#define _CASE_QmQpqppqm_TriP case 2050 : \
          return &C2q2Q1g_2050_TriP
#define _CASE_QmQpqmpqp_nf case 3310 : \
          return &C2q2Q1g_3310_nf
#define _CASE_QmQpqppqm_nf case 2050 : \
          return &C2q2Q1g_2050_nf
#define _CASE_QpQmqmpqp_nf case 3305 : \
          return &C2q2Q1g_3305_nf
#define _CASE_QpQmqppqm_nf case 2045 : \
          return &C2q2Q1g_2045_nf
 
 
 // *************** function definitions using macros ************* 
 
template <class T> SeriesC<T> _C_qmqpQpQmp_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmqpQpQmp_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQpQmp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmqpQpQmp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQpQmp_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qpqmQpQmp_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQpQmp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qpqmQpQmp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQmQpp_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qpqmQmQpp_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQmQpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qpqmQmQpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpQpQm_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmpqpQpQm_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpQmQp_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmpqpQmQp_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_qppqmQmQp_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qppqmQmQp_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppQmQp_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmqppQmQp_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppQpQm_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmqppQpQm_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmpQmQp_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qpqmpQmQp_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmpQpQm_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qpqmpQpQm_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQmpQp_TriP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmqpQmpQp_TriP(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQppQm_TriP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmqpQppQm_TriP(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQmpQp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmqpQmpQp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpQppQm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qmqpQppQm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQmpQp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qpqmQmpQp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qpqmQppQm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_qpqmQppQm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_QmQpqpqmp_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmQpqpqmp_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_QmQpqpqmp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmQpqpqmp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_QpQmqpqmp_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QpQmqpqmp_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_QpQmqpqmp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QpQmqpqmp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_QpQmqmqpp_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QpQmqmqpp_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_QpQmqmqpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QpQmqmqpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_QmpQpqpqm_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmpQpqpqm_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_QmpQpqmqp_PentP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmpQpqmqp_PentP(ep,mu);}
 
template <class T> SeriesC<T> _C_QppQmqmqp_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QppQmqmqp_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_QmQppqmqp_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmQppqmqp_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_QmQppqpqm_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmQppqpqm_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_QpQmpqmqp_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QpQmpqmqp_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_QpQmpqpqm_BoxP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QpQmpqpqm_BoxP(ep,mu);}
 
template <class T> SeriesC<T> _C_QmQpqmpqp_TriP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmQpqmpqp_TriP(ep,mu);}
 
template <class T> SeriesC<T> _C_QmQpqppqm_TriP(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmQpqppqm_TriP(ep,mu);}
 
template <class T> SeriesC<T> _C_QmQpqmpqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmQpqmpqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_QmQpqppqm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QmQpqppqm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_QpQmqmpqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QpQmqmpqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_QpQmqppqm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q2Q1g_QpQmqppqm_nf(ep,mu);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> SeriesC<T> ( *C2q2Q1g_BoxP_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qppqmQmQp_BoxP;
       _CASE_qmqppQmQp_BoxP;
       _CASE_qmqppQpQm_BoxP;
       _CASE_qpqmpQmQp_BoxP;
       _CASE_qpqmpQpQm_BoxP;
       _CASE_QppQmqmqp_BoxP;
       _CASE_QmQppqmqp_BoxP;
       _CASE_QmQppqpqm_BoxP;
       _CASE_QpQmpqmqp_BoxP;
       _CASE_QpQmpqpqm_BoxP;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q2Q1g_nf_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqpQpQmp_nf;
       _CASE_qpqmQpQmp_nf;
       _CASE_qpqmQmQpp_nf;
       _CASE_qmqpQmpQp_nf;
       _CASE_qmqpQppQm_nf;
       _CASE_qpqmQmpQp_nf;
       _CASE_qpqmQppQm_nf;
       _CASE_QmQpqpqmp_nf;
       _CASE_QpQmqpqmp_nf;
       _CASE_QpQmqmqpp_nf;
       _CASE_QmQpqmpqp_nf;
       _CASE_QmQpqppqm_nf;
       _CASE_QpQmqmpqp_nf;
       _CASE_QpQmqppqm_nf;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q2Q1g_PentP_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqpQpQmp_PentP;
       _CASE_qpqmQpQmp_PentP;
       _CASE_qpqmQmQpp_PentP;
       _CASE_qmpqpQpQm_PentP;
       _CASE_qmpqpQmQp_PentP;
       _CASE_QmQpqpqmp_PentP;
       _CASE_QpQmqpqmp_PentP;
       _CASE_QpQmqmqpp_PentP;
       _CASE_QmpQpqpqm_PentP;
       _CASE_QmpQpqmqp_PentP;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q2Q1g_TriP_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqpQmpQp_TriP;
       _CASE_qmqpQppQm_TriP;
       _CASE_QmQpqmpqp_TriP;
       _CASE_QmQpqppqm_TriP;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template SeriesC<R> ( *C2q2Q1g_BoxP_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q1g_BoxP_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q1g_BoxP_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q2Q1g_nf_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q1g_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q1g_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q2Q1g_PentP_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q1g_PentP_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q1g_PentP_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q2Q1g_TriP_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q2Q1g_TriP_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q2Q1g_TriP_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);




}
