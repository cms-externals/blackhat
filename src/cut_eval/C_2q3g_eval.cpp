/*
*C_2q3g.cpp
*
* Created on 10/28, 2008
*      Author: Zvi's script
*/
 
#include "C_2q3g_eval.h"
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


 
 
template <class T> SeriesC<T> C2q3g_qmqpppp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpppp L");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmqpmpp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpmpp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:24:14 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = S(2,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = -(spa12*spb12);
complex<T> t5 = spa15*T(2); 
complex<T> t20 = square(spa13); 
complex<T> t21 = square(spb25); 
complex<T> t22 = square(spb24); 
complex<T> t24 = square(spa35); 
complex<T> t27 = s12 + s15 - s34; 
complex<T> t29 = square(spb45); 
complex<T> t30 = s15 - s23; 
complex<T> t31 = spa15*spa45; 
complex<T> t38 = cube(spb25); 
complex<T> t40 = cube(spb24); 
complex<T> t41 = spa14*spa34; 
complex<T> t43 = -(spa35*T(2)); 
complex<T> t49 = s23*s34; 
complex<T> t54 = cube(spa13); 
complex<T> t56 = square(spa23); 
complex<T> t60 = spa13*spa23; 
complex<T> t69 = spa12*spa23; 
complex<T> t72 = spa35*T(2); 
complex<T> t76 = -(spa14*T(2)); 
complex<T> d5 = spa12*spa45*square(s12 - s34)*T(6); d5 = T(1)/d5;
complex<T> d6 = spa45*cube(s12 - s34)*T(3); d6 = T(1)/d6;
complex<T> d8 = spa12*spa45; d8 = T(1)/d8;
complex<T> d9 = s12 - s34; d9 = T(1)/d9;
complex<T> d21 = spa34*spa45*T(2); d21 = T(1)/d21;
complex<T> d22 = spa12*spa34*spa45; d22 = T(1)/d22;
complex<T> d23 = square(s15 - s34); d23 = T(1)/d23;
complex<T> d28 = spa45*square(spa24); d28 = T(1)/d28;
complex<T> d29 = spa24*spa45; d29 = T(1)/d29;
complex<T> d30 = spa45*cube(spa24); d30 = T(1)/d30;
complex<T> d43 = spa12*spa34*T(2); d43 = T(1)/d43;
complex<T> t10 = spa45*t5; 
complex<T> t23 = -t69; 
complex<T> t32 = -t49; 
complex<T> t45 = d21*spb12; 
complex<T> t46 = spa15*t29; 
complex<T> t61 = t24*t38; 
complex<T> t62 = t40*t41; 
complex<T> t65 = d22*spa15; 
complex<T> t66 = d8*spb45; 
complex<T> t71 = spa14*t22; 
complex<T> t75 = t21*t60; 
complex<T> t81 = spb25*t20; 
complex<T> t98 = spa13*t21; 
complex<T> t103 = -(spb12*t54); 
complex<T> t105 = spb15*t54; 
complex<T> d1 = spa12*spa34*t31; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa34*spa45*t27; d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*spa34*spa45*square(t27); d3 = T(1)/d3;
complex<T> d4 = spa34*spa45*t27*square(s12 - s34)*T(2); d4 = T(1)/d4;
complex<T> d7 = spa34*t31; d7 = T(1)/d7;
complex<T> d12 = s24*t30*t31; d12 = T(1)/d12;
complex<T> d13 = s24*(s15 - s34)*t31; d13 = T(1)/d13;
complex<T> d14 = t30*t31*square(s24); d14 = T(1)/d14;
complex<T> d17 = (s15 - s34)*t31*square(s24); d17 = T(1)/d17;
complex<T> d18 = (s15 - s34)*spa34*spa45*t27; d18 = T(1)/d18;
complex<T> d19 = (s15 - s34)*spa34*spa45*square(t27); d19 = T(1)/d19;
complex<T> d20 = spa34*spa45*t27*square(s15 - s34)*T(2); d20 = T(1)/d20;
complex<T> d25 = spa34*t27*t31; d25 = T(1)/d25;
complex<T> d26 = spa34*spa45*square(t27); d26 = T(1)/d26;
complex<T> d27 = spa34*spa45*cube(t27); d27 = T(1)/d27;
complex<T> d31 = spa34*spa45*t27; d31 = T(1)/d31;
complex<T> d32 = t31*square(spa24); d32 = T(1)/d32;
complex<T> d33 = spa24*t31; d33 = T(1)/d33;
complex<T> d34 = t31*cube(spa24); d34 = T(1)/d34;
complex<T> d35 = t27*t31; d35 = T(1)/d35;
complex<T> d36 = spa45*square(t27); d36 = T(1)/d36;
complex<T> d37 = spa45*cube(t27); d37 = T(1)/d37;
complex<T> d40 = spa12*t5; d40 = T(1)/d40;
complex<T> d44 = spa34*spa45*t27*T(2); d44 = T(1)/d44;
complex<T> d45 = spa34*spa45*cube(t27)*T(2); d45 = T(1)/d45;
complex<T> t17 = t105*t45 + s12*(d45*s15*t61*t69 - d26*s15*spa35*t75 + d44*spa35*spb15*t81); 
complex<T> t33 = -t45; 
complex<T> t50 = d34*spa12; 
complex<T> t64 = d27*spa12; 
complex<T> t70 = -t81; 
complex<T> t73 = spa34*t46; 
complex<T> t99 = d32*t60; 
complex<T> d10 = (s15 - s34)*t10; d10 = T(1)/d10;
complex<T> d11 = t10*square(t30); d11 = T(1)/d11;
complex<T> d15 = s24*t10*square(t30); d15 = T(1)/d15;
complex<T> d16 = s24*t10*square(s15 - s34); d16 = T(1)/d16;
complex<T> d24 = spa12*spa34*t10; d24 = T(1)/d24;
complex<T> d38 = spa34*t10; d38 = T(1)/d38;
complex<T> d39 = spa12*t10; d39 = T(1)/d39;
complex<T> d41 = spa24*t10; d41 = T(1)/d41;
complex<T> d42 = t10*cube(spa24); d42 = T(1)/d42;
complex<T> t1 = -(d10*spb24*t20) + d1*t54 + d19*t23*t61 + d20*t23*t61 + d14*t23*t62 + d15*t23*t62 + d16*t23*t62 + d17*t23*t62 - d11*t60*t71 + d18*t72*t75 + d12*t22*t60*t76 + d13*t22*t60*t76 + d23*spa23*t33*t81 - d23*t56*t65*t98; 
complex<T> t2 = d6*spa34*spa35*spb25*t29*t5 - d1*t54 + d3*t23*t61 + d4*t23*t61 + d2*t72*t75 + d7*d9*spa35*t81 - d9*t20*t66*T(2) + d5*spa13*t73*T(5); 
complex<T> t3 = d10*spb24*t20 + d19*t61*t69 + d20*t61*t69 + d3*t61*t69 + d4*t61*t69 + d16*t62*t69 + d17*t62*t69 + d7*d9*spa35*t70 + d6*spb25*t43*t73 + d18*t43*t75 + d2*t43*t75 + d23*spa23*t45*t81 + d23*t56*t65*t98 + d9*t20*t66*T(2) + d13*t60*t71*T(2) + d24*t54*T(3) - d5*spa13*t73*T(5); 
complex<T> t14 = d14*t62*t69 + d15*t62*t69 + t60*t71*(d11 + d12*T(2)); 
complex<T> t59 = d41*t20*t32 + d42*t41*t49*t69 + spa14*t32*t99; 
complex<T> t93 = t41*t50; 
complex<T> t100 = t61*t64; 
complex<T> t13 = s23*(-(d33*t20) + spa23*t93 + t76*t99); 
complex<T> t15 = d7*t103 + s12*(spa23*t100 + d25*spa35*t70 + d26*t43*t75); 
complex<T> t16 = s15*(spa23*t100 + d26*t43*t75) + spb15*(-(d29*t20) + d30*t41*t69 + d28*t60*t76 + d31*spa35*t81); 
complex<T> t18 = -(d33*s34*t20) + d37*spb34*t61*t69 + d35*spa35*spb34*t70 + d36*spb34*t43*t75 + s34*spa23*t93 + s34*t76*t99; 
complex<T> co1 = d38*s23*t103; 
complex<T> co2 = -(d39*s23*spb34*t54); 
complex<T> co3 = d40*spb34*spb45*t54; 
complex<T> co4 = d43*spb45*t105; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t2*Int(ep,mu,c12,c345) + t1*Int(ep,mu,c15,c234) + t14*Int(ep,mu,c23,c145) + t3*Int(ep,mu,c34,c125) + t15*Int(ep,mu,c1,c2,c345) + t16*Int(ep,mu,c1,c5,c234) + t13*Int(ep,mu,c2,c3,c145) + t18*Int(ep,mu,c3,c4,c125) + co1*Int(ep,mu,c1,c2,c3,c45) + co2*Int(ep,mu,c2,c3,c4,c15) + co3*Int(ep,mu,c3,c4,c5,c12) + t59*Int(ep,mu,c4,c3,c2,c15) + co4*Int(ep,mu,c4,c5,c1,c23) + t17*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqppmp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqppmp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:25:50 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s35 = -(spa35*spb35);
complex<T> t1 = spa15*T(2); 
complex<T> t5 = square(s35); 
complex<T> t15 = square(spb35); 
complex<T> t16 = square(spa14); 
complex<T> t17 = square(spb25); 
complex<T> t18 = spa34*spa45; 
complex<T> t19 = spa15*spa23; 
complex<T> t29 = cube(spb35); 
complex<T> t30 = square(spa24); 
complex<T> t32 = cube(spa14); 
complex<T> t33 = -(spa15*T(2)); 
complex<T> t34 = cube(spb25); 
complex<T> t35 = s12 - s34; 
complex<T> t37 = s12 + s15 - s34; 
complex<T> t38 = s12 - s45; 
complex<T> t45 = s35*T(3); 
complex<T> t46 = square(square(spb35)); 
complex<T> t50 = spa23*T(2); 
complex<T> t54 = s34*s45; 
complex<T> t61 = spa24*spb15; 
complex<T> t73 = spa12*spa45; 
complex<T> t75 = spa13*T(2); 
complex<T> t84 = spa13*spa25; 
complex<T> t95 = spa15*spa24; 
complex<T> d30 = spa23*cube(spa35); d30 = T(1)/d30;
complex<T> d33 = square(square(spa35)); d33 = T(1)/d33;
complex<T> d37 = spa13*spa15*spa45; d37 = T(1)/d37;
complex<T> d39 = spa12*square(square(spa35)); d39 = T(1)/d39;
complex<T> t2 = square(s15 + t35); 
complex<T> t39 = spa14*(t19 + t84) + t75*t95; 
complex<T> t40 = spa14*spa25 + t95*T(2); 
complex<T> t43 = -(spa15*spa34*spb35) - spb25*t73; 
complex<T> t47 = spa13*t18; 
complex<T> t48 = -(spa12*t30); 
complex<T> t55 = d39*spa15; 
complex<T> t60 = spa14*t17; 
complex<T> t66 = t30*t34; 
complex<T> t67 = spa15*t46; 
complex<T> t87 = spa24*t16; 
complex<T> t93 = spa45*t34; 
complex<T> t94 = spa34*t50; 
complex<T> t101 = spa24*t32; 
complex<T> t102 = t33*t46; 
complex<T> t106 = spa14*t29; 
complex<T> t129 = s12*t32; 
complex<T> t133 = spb45*t32; 
complex<T> d1 = spa45*t19*t38; d1 = T(1)/d1;
complex<T> d2 = spa34*t19*t35*(s15 + t35); d2 = T(1)/d2;
complex<T> d3 = spa23*spa34*t1*square(t35); d3 = T(1)/d3;
complex<T> d6 = t19*t35; d6 = T(1)/d6;
complex<T> d7 = s35*spa12*t19*t35; d7 = T(1)/d7;
complex<T> d8 = s35*spa12*t19*t38; d8 = T(1)/d8;
complex<T> d9 = t50*square(t35); d9 = T(1)/d9;
complex<T> d10 = spa12*cube(t35)*T(3); d10 = T(1)/d10;
complex<T> d11 = spa12*cube(t38)*T(3); d11 = T(1)/d11;
complex<T> d12 = spa23*t35*t5; d12 = T(1)/d12;
complex<T> d13 = s35*t50*square(t35); d13 = T(1)/d13;
complex<T> d14 = s35*t50*square(t38); d14 = T(1)/d14;
complex<T> d15 = spa23*t38*t5; d15 = T(1)/d15;
complex<T> d16 = spa12*t35*cube(s35); d16 = T(1)/d16;
complex<T> d17 = spa12*t5*square(t35); d17 = T(1)/d17;
complex<T> d18 = spa12*t45*cube(t35); d18 = T(1)/d18;
complex<T> d19 = spa12*t45*cube(t38); d19 = T(1)/d19;
complex<T> d20 = spa12*t5*square(t38); d20 = T(1)/d20;
complex<T> d21 = spa12*t38*cube(s35); d21 = T(1)/d21;
complex<T> d22 = spa12*spa23*t1*t18*t35; d22 = T(1)/d22;
complex<T> d23 = (s15 - s34)*spa34*t19; d23 = T(1)/d23;
complex<T> d25 = (s15 - s34)*spa34*t19*(s15 + t35); d25 = T(1)/d25;
complex<T> d28 = spa12*spa23*t1*t18; d28 = T(1)/d28;
complex<T> d29 = spa13*spa45*t19; d29 = T(1)/d29;
complex<T> d31 = t19*square(spa35); d31 = T(1)/d31;
complex<T> d32 = t18*t19; d32 = T(1)/d32;
complex<T> d35 = spa23*spa34*cube(s15 + t35); d35 = T(1)/d35;
complex<T> d38 = spa12*t19*square(spa35); d38 = T(1)/d38;
complex<T> d41 = spa23*cube(s15 + t35); d41 = T(1)/d41;
complex<T> d42 = spa13*t19; d42 = T(1)/d42;
complex<T> d43 = t1*t18; d43 = T(1)/d43;
complex<T> d46 = t1*t73; d46 = T(1)/d46;
complex<T> d47 = spa13*spa45*t1; d47 = T(1)/d47;
complex<T> d48 = spa12*spa23*t1; d48 = T(1)/d48;
complex<T> d50 = t18*t50; d50 = T(1)/d50;
complex<T> d51 = spa12*spa23*t1*square(spa35); d51 = T(1)/d51;
complex<T> d52 = t50*cube(spa35); d52 = T(1)/d52;
complex<T> t7 = t54*(t47*(d52 + t55) - d51*spa14*(spa14*(t19 + t84) + t75*t95)); 
complex<T> t21 = -(spa24*t40); 
complex<T> t27 = -(d38*(spa14*(spa15*spa23 + t84) + t75*t95)); 
complex<T> t31 = -t47; 
complex<T> t81 = t66*t73; 
complex<T> t86 = t47*T(2); 
complex<T> t112 = t40*t60; 
complex<T> t113 = d3*spa24; 
complex<T> t117 = d11*spa45; 
complex<T> d4 = spa23*spa34*t2*t35; d4 = T(1)/d4;
complex<T> d5 = (s15 + t35)*t94*square(t35); d5 = T(1)/d5;
complex<T> d24 = t94*square(s15 - s34); d24 = T(1)/d24;
complex<T> d26 = (s15 - s34)*spa23*spa34*t2; d26 = T(1)/d26;
complex<T> d27 = (s15 + t35)*t94*square(s15 - s34); d27 = T(1)/d27;
complex<T> d34 = spa34*t19*t2; d34 = T(1)/d34;
complex<T> d36 = spa23*spa34*t2; d36 = T(1)/d36;
complex<T> d40 = t19*t2; d40 = T(1)/d40;
complex<T> d44 = t2*t94; d44 = T(1)/d44;
complex<T> d45 = t94*cube(s15 + t35); d45 = T(1)/d45;
complex<T> d49 = spa12*t94; d49 = T(1)/d49;
complex<T> t8 = d25*spa24*t112 + d24*t30*t60 - d26*t81 - d27*t81 - d23*spb25*t87; 
complex<T> t9 = -(d32*spb12*t101) - d29*t129 + d30*s12*t31 + d34*s12*t21*t60 + d35*s12*t81 + d33*spa15*spb12*t86 - d31*spa14*spb12*(spa14*(t19 + t84) + t75*t95); 
complex<T> t10 = s34*spa14*t27 + d30*s34*t47 + d40*spb34*t21*t60 + d41*spb34*t81 + s34*t55*t86; 
complex<T> t11 = -(d1*spb13*t32) + d19*t102*t47 + d21*t102*t47 + d14*t29*t47 + d15*t29*t47 + d20*t31*t67 - d1*spb23*t87 + d8*spa14*t15*(spa14*(t19 + t84) + t75*t95) - spa13*t106*t117*T(2); 
complex<T> t12 = d9*spa14*spa34*t15 + d6*spb35*t16 + d10*spa34*t106*t33 + d16*t102*t47 + d18*t102*t47 + d12*t29*t47 + d13*t29*t47 + d2*t21*t60 + d25*t21*t60 - d24*t30*t60 + d17*t31*t67 - t113*t60*t73 + d26*t81 + d27*t81 + d4*t81 + d5*t81 + d23*spb25*t87 + d7*spa14*t15*(spa14*(t19 + t84) + t75*t95) + d28*t101*T(3) - d22*t43*t87*T(3); 
complex<T> t13 = d2*spa24*t112 - d9*spa14*spa34*t15 - d6*spb35*t16 + d12*t29*t31 + d13*t29*t31 + d14*t29*t31 + d15*t29*t31 + d1*spb13*t32 + d17*t47*t67 + d20*t47*t67 + t113*t60*t73 + t106*t117*t75 - d4*t81 - d5*t81 + d16*t67*t86 + d18*t67*t86 + d19*t67*t86 + d21*t67*t86 + d1*spb23*t87 - d7*spa14*t15*(spa14*(t19 + t84) + t75*t95) - d8*spa14*t15*(spa14*(t19 + t84) + t75*t95) + d10*spa15*spa34*t106*T(2) + d22*t43*t87*T(3); 
complex<T> t14 = -(d42*t133) + s45*(spa14*t27 + d30*t47 + t55*t86); 
complex<T> t59 = d36*(spa14*spa25 + spa24*t1)*t60*t61 + d35*s15*t81; 
complex<T> t64 = s12*(d44*t112*t61 + d45*s15*t81); 
complex<T> co1 = d37*spb23*t32; 
complex<T> co2 = d43*spb12*spb23*t101; 
complex<T> co3 = d46*spb23*spb34*t101; 
complex<T> co4 = d47*spb23*t129; 
complex<T> co5 = d48*spb34*spb45*t101; 
complex<T> co6 = d49*t133*t61; 
complex<T> co7 = d50*spb12*t32*t61; 
complex<T> co8 = Complex(0,1); 
SeriesC<T> result = co8*(t13*Int(ep,mu,c12,c345) + t8*Int(ep,mu,c15,c234) + t12*Int(ep,mu,c34,c125) + t11*Int(ep,mu,c45,c123) + t9*Int(ep,mu,c1,c2,c345) + t59*Int(ep,mu,c1,c5,c234) + co1*Int(ep,mu,c2,c3,c145) + t10*Int(ep,mu,c3,c4,c125) + t14*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c1,c2,c3,c45) + t64*Int(ep,mu,c2,c1,c5,c34) + co3*Int(ep,mu,c2,c3,c4,c15) + co4*Int(ep,mu,c3,c2,c1,c45) + co5*Int(ep,mu,c3,c4,c5,c12) + co6*Int(ep,mu,c4,c5,c1,c23) + co7*Int(ep,mu,c5,c1,c2,c34) + t7*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqpppm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpppm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:25:59 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s14 = S(1,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t5 = spa23*spa34; 
complex<T> t11 = square(spa15); 
complex<T> t13 = square(spb34); 
complex<T> t14 = -(s14*spb14); 
complex<T> t18 = spa35*spb23; 
complex<T> t26 = cube(spb34); 
complex<T> t27 = spa25*spb12; 
complex<T> t31 = square(spa13); 
complex<T> t32 = square(spa45); 
complex<T> t36 = spa35*spb13; 
complex<T> d2 = spa12*spa34*square(s12 - s45); d2 = T(1)/d2;
complex<T> d4 = spa12*spa34*cube(s12 - s45)*T(3); d4 = T(1)/d4;
complex<T> d5 = spa12*spa34; d5 = T(1)/d5;
complex<T> d7 = s12 - s45; d7 = T(1)/d7;
complex<T> d15 = spa13*spa34*spa45; d15 = T(1)/d15;
complex<T> d16 = spa34*square(s15 - s23 + s45); d16 = T(1)/d16;
complex<T> d18 = spa34*spa45*T(2); d18 = T(1)/d18;
complex<T> d20 = spa12*spa45*T(2); d20 = T(1)/d20;
complex<T> d21 = spa13*spa34*spa45*T(2); d21 = T(1)/d21;
complex<T> d22 = spa12*spa23*T(2); d22 = T(1)/d22;
complex<T> t16 = -(spa12*t18) - spa13*spa45*spb34*T(3); 
complex<T> t17 = -t27; 
complex<T> t35 = spb14*t11; 
complex<T> t38 = -(d5*spb34); 
complex<T> t39 = t31*t32; 
complex<T> t46 = spa25*t11; 
complex<T> t48 = d2*t13; 
complex<T> t60 = t11*t18; 
complex<T> d1 = (s12 - s45)*spa45*t5; d1 = T(1)/d1;
complex<T> d3 = t5*square(s12 - s45)*T(2); d3 = T(1)/d3;
complex<T> d6 = spa12*spa45*t5*T(2); d6 = T(1)/d6;
complex<T> d8 = (s15 - s23)*t5; d8 = T(1)/d8;
complex<T> d9 = (s15 - s23)*(s15 - s23 + s45)*t5; d9 = T(1)/d9;
complex<T> d10 = (s23 - s45)*spa45*t5; d10 = T(1)/d10;
complex<T> d11 = (s23 - s45)*(s15 - s23 + s45)*t5; d11 = T(1)/d11;
complex<T> d12 = spa13*spa45*t5; d12 = T(1)/d12;
complex<T> d13 = spa45*t5; d13 = T(1)/d13;
complex<T> d14 = t5*square(s15 - s23 + s45); d14 = T(1)/d14;
complex<T> d17 = spa13*t5; d17 = T(1)/d17;
complex<T> d19 = t5*square(s15 - s23 + s45)*T(2); d19 = T(1)/d19;
complex<T> d23 = spa12*t5*T(2); d23 = T(1)/d23;
complex<T> d24 = spa45*t5*T(2); d24 = T(1)/d24;
complex<T> t24 = -t35; 
complex<T> t29 = d6*t16; 
complex<T> t30 = d14*s14; 
complex<T> t41 = d10*t11; 
complex<T> t43 = d3*spa35; 
complex<T> t50 = s14*t35; 
complex<T> t59 = t11*t17; 
complex<T> t1 = d5*d7*spb34*t11 - d7*spa15*spa25*t29 - d1*t11*t36 + t17*t41 - t36*t41 - spa13*spa45*t13*t43 - spa13*spa15*spa45*t48 + d11*t50 + d4*t26*t39*T(2) + d6*t46*T(3); 
complex<T> t2 = d7*spa15*spa25*t29 + d1*t11*t36 + d7*t11*t38 + spa13*spa45*t13*t43 + spa13*spa15*spa45*t48 - d4*t26*t39*T(2); 
complex<T> t34 = -(d12*s12*spa35*t11) + d13*t59; 
complex<T> t40 = d9*t11*t14 + d8*t24; 
complex<T> t45 = d16*spb23*t50 + d15*t60; 
complex<T> t54 = d11*t11*t14 + d8*t35 + t27*t41 + t36*t41 + d9*t50; 
complex<T> t61 = t30*t35; 
complex<T> t49 = -(d17*spa35*spb45*t11) + s45*t61; 
complex<T> co1 = s15*t61; 
complex<T> co2 = d18*spb23*t11*t27; 
complex<T> co3 = d19*s15*s45*t50; 
complex<T> co4 = d20*spb23*spb34*t46; 
complex<T> co5 = d21*s12*t60; 
complex<T> co6 = d22*spb34*spb45*t46; 
complex<T> co7 = -(d23*s15*spb45*t46); 
complex<T> co8 = d24*s15*t59; 
complex<T> co9 = Complex(0,1); 
SeriesC<T> result = co9*(t2*Int(ep,mu,c12,c345) + t40*Int(ep,mu,c15,c234) + t54*Int(ep,mu,c23,c145) + t1*Int(ep,mu,c45,c123) + t34*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + t45*Int(ep,mu,c2,c3,c145) + t49*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c1,c2,c3,c45) + co3*Int(ep,mu,c1,c5,c4,c23) + co4*Int(ep,mu,c2,c3,c4,c15) + co5*Int(ep,mu,c3,c2,c1,c45) + co6*Int(ep,mu,c3,c4,c5,c12) + co7*Int(ep,mu,c4,c5,c1,c23) + co8*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmpqppp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpqppp L");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmqppp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmqppp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:26:00 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> s12 = S(1,2);
complex<T> s23 = S(2,3);
complex<T> s15 = -(spa15*spb15);
complex<T> t2 = spa34*spa45; 
complex<T> t5 = square(spa12); 
complex<T> t7 = square(spa24); 
complex<T> t8 = square(spb45); 
complex<T> t12 = spa12*spa24; 
complex<T> d4 = spa15*T(2); d4 = T(1)/d4;
complex<T> d5 = spa34*T(2); d5 = T(1)/d5;
complex<T> d1 = (s15 - s23)*t2; d1 = T(1)/d1;
complex<T> d2 = t2*square(s15 - s23)*T(2); d2 = T(1)/d2;
complex<T> d3 = spa15*t2*T(2); d3 = T(1)/d3;
complex<T> t4 = d2*spa15*t7*t8 + d1*spb45*t12*T(2); 
complex<T> t16 = d3*t5; 
complex<T> t3 = -(d2*spa15*t7*t8) - d1*spb45*t12*T(2) - t16*T(3); 
complex<T> co1 = -(s12*s23*t16); 
complex<T> co2 = -(d4*spb34*spb45*t5); 
complex<T> co3 = -(d5*spb15*spb45*t5); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t4*Int(ep,mu,c15,c234) + t3*Int(ep,mu,c23,c145) + co1*Int(ep,mu,c1,c2,c3,c45) + co2*Int(ep,mu,c3,c4,c5,c12) + co3*Int(ep,mu,c4,c5,c1,c23));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmpqpmp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpqpmp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:26:55 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb12 = SPB(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s35 = -(spa35*spb35);
complex<T> t3 = spa12*spa23; 
complex<T> t5 = square(s12 - s34); 
complex<T> t17 = square(spb25); 
complex<T> t18 = square(spa14); 
complex<T> t19 = square(spb35); 
complex<T> t20 = square(spa13); 
complex<T> t22 = s12 + s15 - s34; 
complex<T> t24 = -(spa34*T(2)); 
complex<T> t31 = cube(spb25); 
complex<T> t32 = cube(spb35); 
complex<T> t33 = cube(spa14); 
complex<T> t34 = spa34*spa45; 
complex<T> t35 = -(spa24*T(2)); 
complex<T> t55 = spa12*spa45; 
complex<T> t62 = spa24*T(2); 
complex<T> t70 = spa34*T(2); 
complex<T> t74 = spa14*spb25; 
complex<T> d13 = spa15*spa23*spa25*T(2); d13 = T(1)/d13;
complex<T> d15 = s12 - s34; d15 = T(1)/d15;
complex<T> d16 = spa23*square(s15 - s34)*T(2); d16 = T(1)/d16;
complex<T> d21 = spa23*square(spa35); d21 = T(1)/d21;
complex<T> d22 = spa15*spa23*spa35; d22 = T(1)/d22;
complex<T> d23 = spa23*cube(spa35); d23 = T(1)/d23;
complex<T> d35 = spa23*spa45*T(2); d35 = T(1)/d35;
complex<T> t21 = -t55; 
complex<T> t49 = spa24*t31; 
complex<T> t56 = spa13*t19; 
complex<T> t60 = spa14*t17; 
complex<T> t61 = t20*t34; 
complex<T> t65 = spa14*t24; 
complex<T> t71 = -(spb45*t33); 
complex<T> t84 = spb25*t18; 
complex<T> t88 = spb15*t33; 
complex<T> d1 = spa15*spa45*t3*T(2); d1 = T(1)/d1;
complex<T> d2 = spa23*spa25*t5*T(2); d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*spa23*t22; d3 = T(1)/d3;
complex<T> d4 = (s12 - s34)*spa23*square(t22); d4 = T(1)/d4;
complex<T> d5 = spa23*t22*t5*T(2); d5 = T(1)/d5;
complex<T> d6 = t3*t5*T(2); d6 = T(1)/d6;
complex<T> d7 = (s12 - s34)*s35*t3; d7 = T(1)/d7;
complex<T> d8 = s35*(s12 - s45)*t3; d8 = T(1)/d8;
complex<T> d9 = (s12 - s34)*t3*square(s35); d9 = T(1)/d9;
complex<T> d10 = s35*t3*t5*T(2); d10 = T(1)/d10;
complex<T> d11 = s35*t3*square(s12 - s45)*T(2); d11 = T(1)/d11;
complex<T> d12 = (s12 - s45)*t3*square(s35); d12 = T(1)/d12;
complex<T> d14 = spa45*t3; d14 = T(1)/d14;
complex<T> d17 = (s15 - s34)*spa23*t22; d17 = T(1)/d17;
complex<T> d18 = (s15 - s34)*spa23*square(t22); d18 = T(1)/d18;
complex<T> d19 = spa23*t22*square(s15 - s34)*T(2); d19 = T(1)/d19;
complex<T> d20 = spa15*spa45*t3; d20 = T(1)/d20;
complex<T> d24 = spa15*spa23*t22; d24 = T(1)/d24;
complex<T> d25 = spa23*square(t22); d25 = T(1)/d25;
complex<T> d26 = spa23*cube(t22); d26 = T(1)/d26;
complex<T> d27 = spa23*t22; d27 = T(1)/d27;
complex<T> d28 = t3*square(spa35); d28 = T(1)/d28;
complex<T> d29 = spa15*spa35*t3; d29 = T(1)/d29;
complex<T> d30 = t3*cube(spa35); d30 = T(1)/d30;
complex<T> d31 = spa15*t3; d31 = T(1)/d31;
complex<T> d32 = spa15*t55*T(2); d32 = T(1)/d32;
complex<T> d33 = spa15*t3*T(2); d33 = T(1)/d33;
complex<T> d34 = t3*T(2); d34 = T(1)/d34;
complex<T> d36 = spa23*t22*T(2); d36 = T(1)/d36;
complex<T> d37 = spa23*cube(t22)*T(2); d37 = T(1)/d37;
complex<T> d38 = spa15*spa35*t3*T(2); d38 = T(1)/d38;
complex<T> d39 = t3*cube(spa35)*T(2); d39 = T(1)/d39;
complex<T> t12 = d18*t21*t49 + d19*t21*t49 + d16*spa24*t60 + d17*t60*t62; 
complex<T> t14 = d37*s12*s15*t49*t55 - d25*s12*s15*spa24*t60 + d36*s12*spb15*t84 + d35*spb12*t88; 
complex<T> t39 = d28*spa13; 
complex<T> t42 = d2*spb15; 
complex<T> t43 = d26*spa12; 
complex<T> t46 = d11*t32*t61 + d12*t32*t61 + d8*spa14*t56*t70; 
complex<T> t63 = d14*spb35; 
complex<T> t67 = -(d29*spa13); 
complex<T> t82 = d25*t35; 
complex<T> t1 = -(d1*t33) + d4*t21*t49 + d5*t21*t49 - d6*spa14*spa34*t56 - d10*t32*t61 - d11*t32*t61 - d12*t32*t61 - d9*t32*t61 + d3*t60*t62 + d15*t18*t24*t63 + d7*t56*t65 + d8*t56*t65 + t21*t42*t74 + d13*d15*t55*t74; 
complex<T> t2 = d18*t49*t55 + d19*t49*t55 + d4*t49*t55 + d5*t49*t55 + d6*spa14*spa34*t56 - d16*spa24*t60 + d17*t35*t60 + d3*t35*t60 + d10*t32*t61 + d9*t32*t61 + d7*spa14*t56*t70 + d15*t18*t63*t70 + d13*d15*t21*t74 + t42*t55*t74 + d20*t33*T(2); 
complex<T> t11 = -(s34*s45*(d38*spa13*t18 + spa14*spa34*t39 - d39*t61)); 
complex<T> t16 = d30*s45*t61 + s45*t39*t65 + s45*t18*t67 + d31*t71; 
complex<T> t47 = s34*(d26*t21*t49 + d30*t61 + d25*t60*t62 + t39*t65 + t18*t67 + d24*t84); 
complex<T> t78 = spa45*t43; 
complex<T> t13 = -(d22*spa13*spb12*t18) + d23*spb12*t61 + d21*spa13*spb12*t65 + s12*t49*t78 + s12*t60*t82 - d24*s12*t84; 
complex<T> t48 = s15*t49*t78 + s15*t60*t82 + d27*spb15*t84; 
complex<T> co1 = -(d32*s34*spb23*t33); 
complex<T> co2 = d33*s34*t71; 
complex<T> co3 = d34*spb45*t88; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t1*Int(ep,mu,c12,c345) + t12*Int(ep,mu,c15,c234) + t2*Int(ep,mu,c34,c125) + t46*Int(ep,mu,c45,c123) + t13*Int(ep,mu,c1,c2,c345) + t48*Int(ep,mu,c1,c5,c234) + t47*Int(ep,mu,c3,c4,c125) + t16*Int(ep,mu,c4,c5,c123) + co1*Int(ep,mu,c2,c3,c4,c15) + co2*Int(ep,mu,c3,c4,c5,c12) + co3*Int(ep,mu,c4,c5,c1,c23) + t14*Int(ep,mu,c5,c1,c2,c34) + t11*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmpqppm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpqppm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:27:00 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> t2 = spa23*spa34; 
complex<T> t7 = square(spa15); 
complex<T> t8 = spa35*spb34; 
complex<T> t11 = square(spa13); 
complex<T> t15 = square(spb34); 
complex<T> t19 = -(spa35*spb45); 
complex<T> t22 = s45*spa13; 
complex<T> t28 = spa13*spa15; 
complex<T> d4 = spa12*spa24*spa34; d4 = T(1)/d4;
complex<T> d6 = (s15 - s23 + s45)*spa12*spa34; d6 = T(1)/d6;
complex<T> d7 = spa12*spa24; d7 = T(1)/d7;
complex<T> d9 = spa12*spa24*T(2); d9 = T(1)/d9;
complex<T> d10 = spa12*spa45*T(2); d10 = T(1)/d10;
complex<T> d11 = spa12*spa23*T(2); d11 = T(1)/d11;
complex<T> t14 = -(s15*t7); 
complex<T> t27 = -((d4*s23 + d6*spa13*spb14*spb23)*t7); 
complex<T> t29 = spa35*t15; 
complex<T> t33 = t7*(d9*s23*spb34 + d10*spb23*t8); 
complex<T> d1 = spa12*spa45*t2*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa12*t2; d2 = T(1)/d2;
complex<T> d3 = spa12*t2*square(s12 - s45)*T(2); d3 = T(1)/d3;
complex<T> d5 = (s15 - s23 + s45)*spa12*t2; d5 = T(1)/d5;
complex<T> d8 = spa12*t2; d8 = T(1)/d8;
complex<T> d12 = (s15 - s23 + s45)*spa12*t2*T(2); d12 = T(1)/d12;
complex<T> d13 = spa12*t2*T(2); d13 = T(1)/d13;
complex<T> d14 = spa45*t2*T(2); d14 = T(1)/d14;
complex<T> t13 = -(d5*spb14); 
complex<T> t16 = -(d2*T(2)); 
complex<T> t18 = d3*spa45; 
complex<T> t20 = t14*(d13*spa35*spb45 + d12*spb14*t22); 
complex<T> t5 = s15*(d4 + spa13*t13)*t7; 
complex<T> t6 = -(t11*t18*t29) + d2*t28*t8*T(2); 
complex<T> t21 = t11*t18*t29 + t16*t28*t8 + d1*spa35*t7*T(3); 
complex<T> t30 = (d8*t19 + t13*t22)*t7; 
complex<T> co1 = d7*spb34*t7; 
complex<T> co2 = d11*spb45*t7*t8; 
complex<T> co3 = d14*spa35*spb12*t14; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t21*Int(ep,mu,c12,c345) + t6*Int(ep,mu,c45,c123) + t5*Int(ep,mu,c1,c5,c234) + t27*Int(ep,mu,c2,c3,c145) + co1*Int(ep,mu,c3,c4,c125) + t30*Int(ep,mu,c4,c5,c123) + t33*Int(ep,mu,c2,c3,c4,c15) + co2*Int(ep,mu,c3,c4,c5,c12) + t20*Int(ep,mu,c4,c5,c1,c23) + co3*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmppqpp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, p, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmppqpp L");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmpqpp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, p, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmpqpp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:27:02 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = -(spa15*spb15);
complex<T> t1 = spa15*spa34; 
complex<T> t3 = spa45*T(2); 
complex<T> t5 = spa12*spa24; 
complex<T> t8 = square(spa12); 
complex<T> t10 = square(spa14); 
complex<T> t11 = square(spb34); 
complex<T> t18 = -(spa24*spb45); 
complex<T> t25 = spa14*spb34; 
complex<T> d8 = spa15*spa35; d8 = T(1)/d8;
complex<T> d10 = spa23*spa34*T(2); d10 = T(1)/d10;
complex<T> d12 = spa15*spa35*T(2); d12 = T(1)/d12;
complex<T> t12 = spa15*t18 - spa23*t25; 
complex<T> t15 = -(spa14*t8); 
complex<T> t28 = spa24*t10; 
complex<T> d1 = (s15 - s23)*spa45*t1; d1 = T(1)/d1;
complex<T> d2 = t1*t3*square(s15 - s23); d2 = T(1)/d2;
complex<T> d3 = (s15 - s23)*spa23*t1*t3; d3 = T(1)/d3;
complex<T> d4 = spa23*t1*t3; d4 = T(1)/d4;
complex<T> d5 = spa35*t1; d5 = T(1)/d5;
complex<T> d6 = spa13*spa45*t1; d6 = T(1)/d6;
complex<T> d7 = spa45*t1; d7 = T(1)/d7;
complex<T> d9 = spa13*t1; d9 = T(1)/d9;
complex<T> d11 = spa13*t1*t3; d11 = T(1)/d11;
complex<T> t4 = d2*spa23*t11*t28 - d1*t25*t5*T(2) - d3*t12*t5*T(3); 
complex<T> t13 = -(d2*spa23); 
complex<T> t20 = s12*(d6*t15 + d5*t8); 
complex<T> t24 = d9*spb45*t15 - d5*s45*t8; 
complex<T> t27 = d6*s23*t15 - d7*spa24*spb23*t8; 
complex<T> t14 = t11*t13*t28 + d1*t25*t5*T(2) + d3*t12*t5*T(3) - d4*spa24*t8*T(3); 
complex<T> co1 = d8*spb34*t8; 
complex<T> co2 = d10*spb15*t18*t8; 
complex<T> co3 = d11*s12*s23*t15; 
complex<T> co4 = d12*s45*spb34*t8; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t4*Int(ep,mu,c15,c234) + t14*Int(ep,mu,c23,c145) + t20*Int(ep,mu,c1,c2,c345) + t27*Int(ep,mu,c2,c3,c145) + co1*Int(ep,mu,c3,c4,c125) + t24*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c1,c5,c4,c23) + co3*Int(ep,mu,c3,c2,c1,c45) + co4*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmpmqpp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, m, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpmqpp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:27:46 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = S(3,4);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = -(spa12*spb12);
complex<T> t4 = spa15*spa45; 
complex<T> t8 = square(spa13); 
complex<T> t18 = square(spb25); 
complex<T> t19 = square(spb24); 
complex<T> t20 = square(spa14); 
complex<T> t22 = -(spa15*spa23); 
complex<T> t24 = -(spa13*spa34); 
complex<T> t26 = s12 - s34; 
complex<T> t32 = cube(spb25); 
complex<T> t33 = cube(spb24); 
complex<T> t34 = spa23*spa34; 
complex<T> t35 = spa13*spa14; 
complex<T> t36 = spa45*T(2); 
complex<T> t38 = s12*spa35; 
complex<T> t39 = s23*s34; 
complex<T> t51 = spa15*spa35; 
complex<T> t57 = cube(spa13); 
complex<T> t74 = spa13*spb25; 
complex<T> t81 = spa34*spb24; 
complex<T> t88 = spb12*spb25; 
complex<T> d20 = s15 - s34; d20 = T(1)/d20;
complex<T> d25 = spa12*spa24*spa45; d25 = T(1)/d25;
complex<T> d26 = spa45*cube(spa24); d26 = T(1)/d26;
complex<T> d27 = spa45*square(spa24); d27 = T(1)/d27;
complex<T> d33 = spa12*spa23*T(2); d33 = T(1)/d33;
complex<T> t9 = t4*T(2); 
complex<T> t27 = -t38; 
complex<T> t29 = -(spa15*spa23*spb25) - spa12*t81; 
complex<T> t60 = t20*t34; 
complex<T> t61 = spa13*t18; 
complex<T> t67 = spa23*t32; 
complex<T> t68 = spa34*t35; 
complex<T> t73 = spa35*t22; 
complex<T> t90 = s15*t38; 
complex<T> d1 = t36*square(t26); d1 = T(1)/d1;
complex<T> d2 = spa45*t26*(s15 + t26); d2 = T(1)/d2;
complex<T> d3 = spa45*t26*square(s15 + t26); d3 = T(1)/d3;
complex<T> d4 = (s15 + t26)*t36*square(t26); d4 = T(1)/d4;
complex<T> d5 = spa12*spa23*t4; d5 = T(1)/d5;
complex<T> d6 = (s15 - s23)*s24*t4; d6 = T(1)/d6;
complex<T> d8 = s24*(s15 - s34)*t4; d8 = T(1)/d8;
complex<T> d9 = (s15 - s23)*t4*square(s24); d9 = T(1)/d9;
complex<T> d12 = (s15 - s34)*t4*square(s24); d12 = T(1)/d12;
complex<T> d13 = spa25*t36*square(s15 - s34); d13 = T(1)/d13;
complex<T> d14 = (s15 - s34)*spa45*(s15 + t26); d14 = T(1)/d14;
complex<T> d15 = (s15 - s34)*spa45*square(s15 + t26); d15 = T(1)/d15;
complex<T> d16 = (s15 + t26)*t36*square(s15 - s34); d16 = T(1)/d16;
complex<T> d18 = spa23*t4; d18 = T(1)/d18;
complex<T> d19 = spa12*spa25*t36; d19 = T(1)/d19;
complex<T> d22 = spa45*(s15 + t26); d22 = T(1)/d22;
complex<T> d23 = spa45*square(s15 + t26); d23 = T(1)/d23;
complex<T> d24 = spa45*cube(s15 + t26); d24 = T(1)/d24;
complex<T> d28 = spa12*spa45*(s15 + t26); d28 = T(1)/d28;
complex<T> d29 = spa12*spa24*t4; d29 = T(1)/d29;
complex<T> d30 = t4*cube(spa24); d30 = T(1)/d30;
complex<T> d31 = t4*square(spa24); d31 = T(1)/d31;
complex<T> d32 = spa12*t4; d32 = T(1)/d32;
complex<T> d34 = (s15 + t26)*t36; d34 = T(1)/d34;
complex<T> d35 = t36*cube(s15 + t26); d35 = T(1)/d35;
complex<T> d39 = t36*square(s15 + t26); d39 = T(1)/d39;
complex<T> t14 = -(d32*spb23*t57) + d30*s23*t60 + d31*s23*t68 - d29*s23*spa14*square(spa13); 
complex<T> t16 = d1*spa35*t61 - d2*spa35*t61 + (d3 + d4)*t32*t73; 
complex<T> t40 = d23*s15; 
complex<T> t41 = d13*spb12; 
complex<T> t43 = d24*spa15; 
complex<T> t48 = s34*(d30*t60 - d23*spa35*t61 + d31*t68 + d24*t32*t73 - d29*spa14*t8 + d28*spb25*t8); 
complex<T> t53 = -(d18*T(2)); 
complex<T> t54 = d6*t19; 
complex<T> t62 = d18*T(2); 
complex<T> d7 = t9*square(s15 - s34); d7 = T(1)/d7;
complex<T> d10 = s24*t9*square(s15 - s23); d10 = T(1)/d10;
complex<T> d11 = s24*t9*square(s15 - s34); d11 = T(1)/d11;
complex<T> d17 = (s15 - s34)*spa12*spa23*t9; d17 = T(1)/d17;
complex<T> d21 = spa12*spa23*t9; d21 = T(1)/d21;
complex<T> d36 = spa12*spa24*t9; d36 = T(1)/d36;
complex<T> d37 = t9*cube(spa24); d37 = T(1)/d37;
complex<T> d38 = t9*square(spa24); d38 = T(1)/d38;
complex<T> t1 = d8*spa14*t19*t24 + d21*t57 + d11*t33*t60 + d12*t33*t60 - d1*spa35*t61 + d14*spa35*t61 + d2*spa35*t61 + d15*t51*t67 + d16*t51*t67 + d3*t51*t67 + d4*t51*t67 + d7*t19*t68 + d19*d20*t22*t74 + spa15*spa23*t41*t74 + d20*t62*t8*t81 + d17*t29*t8*T(3); 
complex<T> t2 = d7*spa14*t19*t24 - d10*t33*t60 - d11*t33*t60 - d12*t33*t60 - d9*t33*t60 - d14*spa35*t61 + d8*t19*t68 + t54*t68 + d15*t32*t73 + d16*t32*t73 + d19*d20*spa15*spa23*t74 + t22*t41*t74 + d20*t53*t8*t81 - d5*t57*T(2) - d17*t29*t8*T(3); 
complex<T> t13 = t39*(d31*spa14*t24 + d37*t60 - d36*spa14*square(spa13)); 
complex<T> t17 = d26*spb15*t60 + spa35*t40*t61 + s15*spa35*t43*t67 + d27*spb15*t68 - d25*spa14*spb15*square(spa13) - d28*s15*spb25*square(spa13); 
complex<T> t58 = spa14*t24*t54 + (d10 + d9)*t33*t60; 
complex<T> t59 = d23*t38*t61 + t38*t43*t67 + d22*t8*t88; 
complex<T> t65 = t27*t40*t61 + d34*s15*t8*t88 + d35*spa15*t67*t90; 
complex<T> co1 = -(d33*spb15*spb45*t57); 
complex<T> co2 = d38*t39*t68*T(3); 
complex<T> co3 = d39*t61*t90*T(3); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t16*Int(ep,mu,c12,c345) + t2*Int(ep,mu,c15,c234) + t58*Int(ep,mu,c23,c145) + t1*Int(ep,mu,c34,c125) + t59*Int(ep,mu,c1,c2,c345) + t17*Int(ep,mu,c1,c5,c234) + t14*Int(ep,mu,c2,c3,c145) + t48*Int(ep,mu,c3,c4,c125) + co1*Int(ep,mu,c1,c5,c4,c23) + t65*Int(ep,mu,c2,c1,c5,c34) + t13*Int(ep,mu,c2,c3,c4,c15) + co2*Int(ep,mu,c4,c3,c2,c15) + co3*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmppqpm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, p, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmppqpm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:27:47 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa45 = SPA(4,5);
complex<T> s15 = S(1,5);
complex<T> s45 = S(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> t7 = square(spa15); 
complex<T> t9 = square(spa35); 
complex<T> t10 = square(spb23); 
complex<T> t18 = spa35*spb23; 
complex<T> d1 = spa12*spa23*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa23*spa34; d2 = T(1)/d2;
complex<T> d3 = spa23*spa34*square(s12 - s45)*T(2); d3 = T(1)/d3;
complex<T> d4 = (s12 - s45)*spa12*spa23*spa34*T(2); d4 = T(1)/d4;
complex<T> d5 = spa23*spa34*T(2); d5 = T(1)/d5;
complex<T> d6 = spa12*spa23*T(2); d6 = T(1)/d6;
complex<T> t11 = -(spa13*spa45*spb34) - spa12*t18; 
complex<T> t3 = -(d3*spa12*t10*t9) - d2*spa15*t18*T(2) - d4*spa15*t11*T(3); 
complex<T> t4 = d3*spa12*t10*t9 + d2*spa15*t18*T(2) + d4*spa15*t11*T(3) + d1*t7*T(3); 
complex<T> co1 = -(d5*s15*spb12*t7); 
complex<T> co2 = -(d6*s45*spb34*t7); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t4*Int(ep,mu,c12,c345) + t3*Int(ep,mu,c45,c123) + co1*Int(ep,mu,c2,c1,c5,c34) + co2*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmpppqp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, p, p, qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpppqp L");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmppqp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, p, p, qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmppqp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:27:52 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa14 = SPA(1,4);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s14 = -(spa14*spb14);
complex<T> s15 = -(spa15*spb15);
complex<T> s45 = -(spa45*spb45);
complex<T> t3 = square(s15 - s23 + s45); 
complex<T> t4 = spa34*spa45; 
complex<T> t13 = square(spa12); 
complex<T> t19 = square(spb34); 
complex<T> t21 = -(spa25*spb15); 
complex<T> t25 = s14*spa24; 
complex<T> t26 = spa14*spa23; 
complex<T> t31 = -(spa15*spb45); 
complex<T> t34 = spa12*spa25; 
complex<T> d3 = (s15 - s23)*spa15*spa34; d3 = T(1)/d3;
complex<T> d6 = spa15*spa34; d6 = T(1)/d6;
complex<T> d7 = s15 - s23; d7 = T(1)/d7;
complex<T> d15 = spa13*spa34; d15 = T(1)/d15;
complex<T> t8 = t4*T(2); 
complex<T> t14 = -t25; 
complex<T> t17 = spa24*t31 - spb34*t26*T(3); 
complex<T> t20 = -(spb34*t26) + spa24*t31; 
complex<T> t24 = spb14*t13; 
complex<T> t42 = -(d6*spb34); 
complex<T> d2 = (s15 - s23)*(s15 - s23 + s45)*spa23*t4; d2 = T(1)/d2;
complex<T> d8 = (s23 - s45)*t4; d8 = T(1)/d8;
complex<T> d9 = (s23 - s45)*(s15 - s23 + s45)*spa23*t4; d9 = T(1)/d9;
complex<T> d10 = (s23 - s45)*spa23*t4; d10 = T(1)/d10;
complex<T> d11 = spa13*t4; d11 = T(1)/d11;
complex<T> d12 = spa23*t3*t4; d12 = T(1)/d12;
complex<T> d13 = spa23*t4; d13 = T(1)/d13;
complex<T> d14 = t3*t4; d14 = T(1)/d14;
complex<T> d16 = spa23*spa34*t3; d16 = T(1)/d16;
complex<T> d18 = spa23*spa34*t3*T(2); d18 = T(1)/d18;
complex<T> t22 = -(d8*spb13); 
complex<T> t40 = t14*t24; 
complex<T> t44 = t24*t25; 
complex<T> t51 = -(d11*t13); 
complex<T> d1 = spa15*spa23*t8; d1 = T(1)/d1;
complex<T> d4 = t8*square(s15 - s23); d4 = T(1)/d4;
complex<T> d5 = (s15 - s23)*spa15*spa23*t8; d5 = T(1)/d5;
complex<T> d17 = spa13*t8; d17 = T(1)/d17;
complex<T> t27 = d1*t17; 
complex<T> t35 = d4*spa24; 
complex<T> t39 = d10*t13*t21 + t13*t22 + d9*t44; 
complex<T> t43 = d13*t13*t21 + d12*s15*t44; 
complex<T> t50 = d14*spb23*t44 + s23*t51; 
complex<T> t55 = spb45*t40; 
complex<T> t1 = d8*spb13*t13 + d10*spa25*spb15*t13 - d3*spb34*t13 + d6*d7*spb34*t13 - d7*t27*t34 - t19*t26*t35 + d9*t40 + d2*t44 + d5*t20*t34*T(3); 
complex<T> t2 = d3*spb34*t13 + d7*t27*t34 + t19*t26*t35 + d2*t40 + d7*t13*t42 - d1*spa25*t13*T(3) - d5*t20*t34*T(3); 
complex<T> t33 = -(d15*spb45*t13) + d16*t55; 
complex<T> co1 = s12*t51; 
complex<T> co2 = -(d17*s12*s23*t13); 
complex<T> co3 = d18*s15*t55; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t2*Int(ep,mu,c15,c234) + t1*Int(ep,mu,c23,c145) + t39*Int(ep,mu,c45,c123) + co1*Int(ep,mu,c1,c2,c345) + t43*Int(ep,mu,c1,c5,c234) + t50*Int(ep,mu,c2,c3,c145) + t33*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c1,c2,c3,c45) + co3*Int(ep,mu,c4,c5,c1,c23));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmpmpqp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, m, p, qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpmpqp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:29:30 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = -(spa14*spb14);
complex<T> s15 = -(spa15*spb15);
complex<T> s12 = -(spa12*spb12);
complex<T> s24 = -(spa24*spb24);
complex<T> t2 = spa12*T(2); 
complex<T> t3 = (s15 - s34)*spa45; 
complex<T> t4 = square(s15 - s23 + s45); 
complex<T> t5 = (s12 - s34)*spa34; 
complex<T> t15 = square(spb25); 
complex<T> t16 = square(spb24); 
complex<T> t17 = -(spa13*spa35); 
complex<T> t18 = spa12*spa45; 
complex<T> t27 = spa34*T(2); 
complex<T> t28 = square(spa13); 
complex<T> t29 = square(spa35); 
complex<T> t30 = cube(spb25); 
complex<T> t32 = s12 + s15 - s34; 
complex<T> t33 = cube(spb24); 
complex<T> t41 = spa45*T(2); 
complex<T> t46 = s12*s15; 
complex<T> t49 = cube(spa13); 
complex<T> t50 = spa15*spa23; 
complex<T> t52 = -(spa14*T(3)); 
complex<T> t54 = s14*spb14; 
complex<T> t64 = spa14*spa34; 
complex<T> t67 = -(spa12*T(2)); 
complex<T> t76 = spa14*spa35; 
complex<T> t77 = spa34*spa45; 
complex<T> t87 = -(spa14*spa25); 
complex<T> d10 = (s15 - s23)*s24*spa15; d10 = T(1)/d10;
complex<T> d11 = s24*(s15 - s34)*spa15; d11 = T(1)/d11;
complex<T> d13 = (s15 - s23)*s24*spa15*spa45; d13 = T(1)/d13;
complex<T> d17 = (s15 - s23)*spa45*square(s24); d17 = T(1)/d17;
complex<T> d31 = spa45*cube(spa24); d31 = T(1)/d31;
complex<T> d33 = square(spa24); d33 = T(1)/d33;
complex<T> d34 = spa45*square(spa24); d34 = T(1)/d34;
complex<T> d38 = spa15*square(spa24); d38 = T(1)/d38;
complex<T> d39 = spa15*spa45*square(spa24); d39 = T(1)/d39;
complex<T> d46 = spa15*square(spa24)*T(2); d46 = T(1)/d46;
complex<T> t1 = square(t32); 
complex<T> t31 = -(spa13*spa25) + spa35*t67; 
complex<T> t34 = -t50; 
complex<T> t35 = t67*t76 + spa13*(-t18 + t87); 
complex<T> t36 = -t64; 
complex<T> t38 = -t54; 
complex<T> t51 = t29*t30; 
complex<T> t55 = spa45*t27; 
complex<T> t65 = spa23*t33; 
complex<T> t75 = spa13*t15; 
complex<T> t78 = d31*spa23; 
complex<T> t90 = spa13*t16; 
complex<T> t91 = t76*T(3); 
complex<T> t99 = t15*t17; 
complex<T> d1 = t18*t5; d1 = T(1)/d1;
complex<T> d3 = spa45*t32*t5; d3 = T(1)/d3;
complex<T> d4 = t18*t32*t5; d4 = T(1)/d4;
complex<T> d7 = t2*t50*t77; d7 = T(1)/d7;
complex<T> d8 = (s15 - s23)*(s15 - s23 + s45)*spa23*t18; d8 = T(1)/d8;
complex<T> d9 = spa12*t3; d9 = T(1)/d9;
complex<T> d12 = t41*square(s15 - s34); d12 = T(1)/d12;
complex<T> d14 = s24*spa15*t3; d14 = T(1)/d14;
complex<T> d15 = (s15 - s23)*s24*spa15*t18; d15 = T(1)/d15;
complex<T> d16 = s24*spa12*spa15*t3; d16 = T(1)/d16;
complex<T> d18 = s24*t41*square(s15 - s23); d18 = T(1)/d18;
complex<T> d19 = s24*t41*square(s15 - s34); d19 = T(1)/d19;
complex<T> d20 = t3*square(s24); d20 = T(1)/d20;
complex<T> d21 = t2*t77*square(s15 - s34); d21 = T(1)/d21;
complex<T> d22 = spa34*t3*t32; d22 = T(1)/d22;
complex<T> d23 = spa12*spa34*t3*t32; d23 = T(1)/d23;
complex<T> d26 = (s15 - s23)*spa23*t18; d26 = T(1)/d26;
complex<T> d27 = (s23 - s45)*spa23*t18; d27 = T(1)/d27;
complex<T> d28 = (s23 - s45)*(s15 - s23 + s45)*spa23*t18; d28 = T(1)/d28;
complex<T> d30 = t77*cube(t32); d30 = T(1)/d30;
complex<T> d32 = spa23*t18*t4; d32 = T(1)/d32;
complex<T> d35 = spa23*spa34*t18; d35 = T(1)/d35;
complex<T> d36 = t18*square(spa24); d36 = T(1)/d36;
complex<T> d40 = spa15*t18*square(spa24); d40 = T(1)/d40;
complex<T> d41 = t18*t4; d41 = T(1)/d41;
complex<T> d44 = spa45*cube(t32); d44 = T(1)/d44;
complex<T> d45 = spa12*spa23*t4; d45 = T(1)/d45;
complex<T> d47 = t41*cube(spa24); d47 = T(1)/d47;
complex<T> d48 = spa15*t41*square(spa24); d48 = T(1)/d48;
complex<T> d49 = spa15*spa45*t2*square(spa24); d49 = T(1)/d49;
complex<T> d50 = spa23*t2*t4; d50 = T(1)/d50;
complex<T> t62 = t49*(d27*spb14 + d28*t54); 
complex<T> t83 = d40*t35; 
complex<T> t85 = t50*t51; 
complex<T> t95 = t29*t75; 
complex<T> t111 = t38*t49; 
complex<T> t119 = spa13*t91; 
complex<T> d2 = t55*square(s12 - s34); d2 = T(1)/d2;
complex<T> d5 = spa45*t1*t5; d5 = T(1)/d5;
complex<T> d6 = t32*t55*square(s12 - s34); d6 = T(1)/d6;
complex<T> d24 = spa34*t1*t3; d24 = T(1)/d24;
complex<T> d25 = t32*t55*square(s15 - s34); d25 = T(1)/d25;
complex<T> d29 = t1*t77; d29 = T(1)/d29;
complex<T> d37 = spa34*t1*t18; d37 = T(1)/d37;
complex<T> d42 = spa45*t1; d42 = T(1)/d42;
complex<T> d43 = t1*t18; d43 = T(1)/d43;
complex<T> d51 = t1*t55; d51 = T(1)/d51;
complex<T> d52 = t55*cube(t32); d52 = T(1)/d52;
complex<T> t7 = s23*s34*(d48*t119 + d46*t28 + d47*spa23*t64 + d49*spa13*(-(spa13*t18) + t67*t76 + spa13*t87)); 
complex<T> t8 = d39*s23*t119 + d38*s23*t28 + d41*spb23*t49*t54 + s23*t64*t78 + d40*s23*spa13*(-(spa13*t18) + t67*t76 + spa13*t87); 
complex<T> t9 = -(d1*spa35*spb25*t28) + d5*t34*t51 + d6*t34*t51 + d2*t95 + d4*t31*t99 - d3*t95*T(3); 
complex<T> t10 = d9*spb24*t28 + d1*spa35*spb25*t28 - d11*t16*t28 + d19*t64*t65 + d20*t64*t65 + d23*spa35*t31*t75 + d4*spa35*t31*t75 + d24*t85 + d25*t85 + d5*t85 + d6*t85 + d12*spa34*t90 + d14*spa35*t52*t90 - d16*(t67*t76 + spa13*(-t18 + t87))*t90 - d2*t95 + d21*t50*t99 + d22*t95*T(3) + d3*t95*T(3); 
complex<T> t12 = d39*s34*t119 + d38*s34*t28 + d43*spa35*spb34*t31*t75 + s34*t64*t78 + d44*spb34*t85 + d40*s34*spa13*(-(spa13*t18) + t67*t76 + spa13*t87) + d42*spb34*t95*T(3); 
complex<T> t13 = d28*t111 - d26*spa35*spb45*t28 - d10*t16*t28 - d27*spb14*t49 + d8*t49*t54 + d17*t64*t65 + d18*t64*t65 + d13*spa35*t52*t90 - d15*(t67*t76 + spa13*(-t18 + t87))*t90; 
complex<T> t14 = d8*t111 + d13*t119*t16 + d14*t119*t16 - d9*spb24*t28 + d26*spa35*spb45*t28 + d10*t16*t28 + d11*t16*t28 + d24*t34*t51 + d25*t34*t51 + d17*t36*t65 + d18*t36*t65 + d19*t36*t65 + d20*t36*t65 + d21*spa35*t50*t75 - d12*spa34*t90 + d15*(t67*t76 + spa13*(-t18 + t87))*t90 + d16*(t67*t76 + spa13*(-t18 + t87))*t90 + d23*t31*t99 - d7*spa35*t49*T(3) - d22*t95*T(3); 
complex<T> t74 = d52*t46*t85 + d51*s15*spb12*t31*t99 + d51*t46*t95*T(3); 
complex<T> t86 = d29*T(3); 
complex<T> t11 = d34*spb15*t119 + d33*spb15*t28 - d35*spa35*spb15*t49 + d32*s15*t49*t54 + d37*s15*spa35*t31*t75 + s15*t36*t78 + d30*s15*t85 + d36*spa13*spb15*(-(spa13*t18) + t67*t76 + spa13*t87) + s15*t86*t95; 
complex<T> t61 = d30*s12*t85 + s12*t86*t95 + d29*spb12*t31*t99; 
complex<T> co1 = d45*spb45*t111; 
complex<T> co2 = d50*s15*spb45*t111; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t9*Int(ep,mu,c12,c345) + t14*Int(ep,mu,c15,c234) + t13*Int(ep,mu,c23,c145) + t10*Int(ep,mu,c34,c125) + t62*Int(ep,mu,c45,c123) + t61*Int(ep,mu,c1,c2,c345) + t11*Int(ep,mu,c1,c5,c234) + t8*Int(ep,mu,c2,c3,c145) + t12*Int(ep,mu,c3,c4,c125) + co1*Int(ep,mu,c4,c5,c123) + t7*Int(ep,mu,c2,c3,c4,c15) + co2*Int(ep,mu,c4,c5,c1,c23) + t74*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmppmqp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, p, m, qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmppmqp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:30:23 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> s45 = S(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s35 = -(spa35*spb35);
complex<T> t5 = spa12*spa23; 
complex<T> t9 = square(spa14); 
complex<T> t12 = square(s35); 
complex<T> t21 = square(spb25); 
complex<T> t22 = square(spb35); 
complex<T> t23 = square(spa24); 
complex<T> t25 = s12 - s34; 
complex<T> t26 = s12 + s15 - s34; 
complex<T> t29 = -square(spa14); 
complex<T> t32 = square(spb23); 
complex<T> t36 = cube(spb25); 
complex<T> t38 = cube(spb35); 
complex<T> t40 = spa13*spa34; 
complex<T> t44 = s15*spa24; 
complex<T> t54 = square(spa45); 
complex<T> t57 = spa15*spa45; 
complex<T> t66 = spa24*spb25; 
complex<T> t69 = spa14*spa45; 
complex<T> t88 = spa14*spa34; 
complex<T> t98 = spa45*spb25; 
complex<T> d5 = spa23*spa34*T(2); d5 = T(1)/d5;
complex<T> d6 = spa15*spa23*spa34; d6 = T(1)/d6;
complex<T> d18 = (s15 - s34)*spa15*spa23; d18 = T(1)/d18;
complex<T> d19 = spa15*spa23*square(s15 - s34)*T(2); d19 = T(1)/d19;
complex<T> d23 = spa15*spa23; d23 = T(1)/d23;
complex<T> d25 = s15 - s34; d25 = T(1)/d25;
complex<T> d26 = spa23*spa35; d26 = T(1)/d26;
complex<T> d27 = spa23*cube(spa35); d27 = T(1)/d27;
complex<T> d28 = spa23*square(spa35); d28 = T(1)/d28;
complex<T> t24 = -t57; 
complex<T> t30 = t5*T(2); 
complex<T> t33 = -(spa15*spa34*spb35) - spa12*t98; 
complex<T> t37 = -t69; 
complex<T> t47 = -(d5*spb15); 
complex<T> t52 = -(spa14*t21); 
complex<T> t58 = t23*t36; 
complex<T> t59 = t38*t40; 
complex<T> t62 = d6*spa12; 
complex<T> t63 = -(d23*T(2)); 
complex<T> t70 = spa24*t21; 
complex<T> t71 = spa13*t22; 
complex<T> t73 = d23*T(2); 
complex<T> d2 = spa23*spa34*t25*(s15 + t25); d2 = T(1)/d2;
complex<T> d3 = spa23*spa34*t25*square(s15 + t25); d3 = T(1)/d3;
complex<T> d4 = spa23*spa34*(s15 + t25)*square(t25)*T(2); d4 = T(1)/d4;
complex<T> d7 = square(t25); d7 = T(1)/d7;
complex<T> d9 = s35*t25*t5; d9 = T(1)/d9;
complex<T> d11 = s35*(s12 - s45)*t5; d11 = T(1)/d11;
complex<T> d12 = t12*t25*t5; d12 = T(1)/d12;
complex<T> d15 = (s12 - s45)*t12*t5; d15 = T(1)/d15;
complex<T> d17 = spa15*spa34*t5; d17 = T(1)/d17;
complex<T> d20 = (s15 - s34)*spa23*spa34*(s15 + t25); d20 = T(1)/d20;
complex<T> d21 = (s15 - s34)*spa23*spa34*square(s15 + t25); d21 = T(1)/d21;
complex<T> d22 = spa23*spa34*(s15 + t25)*square(s15 - s34)*T(2); d22 = T(1)/d22;
complex<T> d24 = spa34*t5; d24 = T(1)/d24;
complex<T> d29 = spa23*spa34*(s15 + t25); d29 = T(1)/d29;
complex<T> d30 = spa23*spa34*square(s15 + t25); d30 = T(1)/d30;
complex<T> d31 = spa23*spa34*cube(s15 + t25); d31 = T(1)/d31;
complex<T> d32 = spa34*(s15 + t25)*t5; d32 = T(1)/d32;
complex<T> d33 = spa35*t5; d33 = T(1)/d33;
complex<T> d34 = t5*cube(spa35); d34 = T(1)/d34;
complex<T> d35 = t5*square(spa35); d35 = T(1)/d35;
complex<T> d36 = (s15 + t25)*t5; d36 = T(1)/d36;
complex<T> d37 = spa23*square(s15 + t25); d37 = T(1)/d37;
complex<T> d38 = spa23*cube(s15 + t25); d38 = T(1)/d38;
complex<T> d39 = spa23*spa34*(s15 + t25)*T(2); d39 = T(1)/d39;
complex<T> d40 = spa23*spa34*square(s15 + t25)*T(2); d40 = T(1)/d40;
complex<T> d41 = spa23*spa34*cube(s15 + t25)*T(2); d41 = T(1)/d41;
complex<T> t1 = d21*t24*t58 + d22*t24*t58 + d20*t37*t70 + d19*spa12*t32*t88 + d18*spb23*t9 + d25*spb23*t63*t9 + d24*d25*t66*t9 - d17*cube(spa14); 
complex<T> t46 = d34*spa15; 
complex<T> t76 = t57*t58; 
complex<T> t78 = d35*spa13; 
complex<T> d1 = spa15*spa34*t30; d1 = T(1)/d1;
complex<T> d8 = t25*t30; d8 = T(1)/d8;
complex<T> d10 = t30*square(s12 - s45); d10 = T(1)/d10;
complex<T> d13 = s35*t30*square(t25); d13 = T(1)/d13;
complex<T> d14 = s35*t30*square(s12 - s45); d14 = T(1)/d14;
complex<T> d16 = spa15*spa34*t25*t30; d16 = T(1)/d16;
complex<T> d42 = spa35*t30; d42 = T(1)/d42;
complex<T> d43 = t30*cube(spa35); d43 = T(1)/d43;
complex<T> d44 = t30*square(spa35); d44 = T(1)/d44;
complex<T> t2 = d8*spb35*t29 + d3*t24*t58 + d4*t24*t58 + d12*t24*t59 + d13*t24*t59 + d14*t24*t59 + d15*t24*t59 + d7*t52*t54*t62 + d2*t37*t70 + d10*t37*t71 + d11*t69*t71 + d9*t69*t71 + d5*d7*spb15*t29*t98 - d1*cube(spa14) - d16*t33*t9*T(3); 
complex<T> t3 = d18*spb23*t29 + d12*t57*t59 + d13*t57*t59 + d7*spa14*t21*t54*t62 + d24*d25*t29*t66 + d2*t69*t70 + d20*t69*t70 + d9*t37*t71 + d21*t76 + d22*t76 + d3*t76 + d4*t76 - d19*spa12*t32*t88 + d8*spb35*t9 + d25*spb23*t73*t9 + d5*d7*spb15*t9*t98 + d16*t33*t9*T(3); 
complex<T> t15 = d26*spb12*t29 + d27*spb12*t40*t57 + d28*spa13*spb12*t69 + d30*s12*t69*t70 + d31*s12*t76 + d29*spb12*t66*square(spa14); 
complex<T> t16 = d32*spb25*t29*t44 + d30*t21*t44*t69 + d31*s15*t76 - d24*spb15*cube(spa14); 
complex<T> t17 = d40*s12*t21*t44*t69 + d41*s12*s15*t76 + d39*spb12*spb25*t44*square(spa14); 
complex<T> t20 = d14*t57*t59 + d15*t57*t59 + d11*t37*t71 + d10*t69*t71; 
complex<T> t35 = s34*s45*(d43*t40*t57 + d44*spa13*t69 - d42*t9); 
complex<T> t96 = t40*t46; 
complex<T> t14 = s45*(d33*t29 + t69*t78 + spa45*t96); 
complex<T> t18 = d33*s34*t29 + d36*spb34*t29*t66 + d37*spb34*t69*t70 + d38*spb34*t76 + s34*t69*t78 + s34*spa45*t96; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*Int(ep,mu,c12,c345) + t1*Int(ep,mu,c15,c234) + t3*Int(ep,mu,c34,c125) + t20*Int(ep,mu,c45,c123) + t15*Int(ep,mu,c1,c2,c345) + t16*Int(ep,mu,c1,c5,c234) + t18*Int(ep,mu,c3,c4,c125) + t14*Int(ep,mu,c4,c5,c123) + t17*Int(ep,mu,c2,c1,c5,c34) + t35*Int(ep,mu,c3,c4,c5,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqpmmm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpmmm L");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmqppmm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqppmm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:30:33 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> s23 = S(2,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> t3 = square(s12 + s15 - s34); 
complex<T> t6 = spb15*spb45; 
complex<T> t13 = square(spb23); 
complex<T> t14 = square(spa45); 
complex<T> t17 = spa15*spb12*spb35 + spa45*spb25*spb34*T(3); 
complex<T> t26 = s25*spa25; 
complex<T> t27 = cube(spa45); 
complex<T> t29 = spa12*spb13; 
complex<T> t33 = square(spb25); 
complex<T> t34 = square(spb34); 
complex<T> t45 = spa34*spb13; 
complex<T> d1 = spb12*spb45; d1 = T(1)/d1;
complex<T> d3 = s12 - s34; d3 = T(1)/d3;
complex<T> d4 = spb12*spb45*square(s12 - s34); d4 = T(1)/d4;
complex<T> d5 = spb12*spb45*cube(s12 - s34)*T(3); d5 = T(1)/d5;
complex<T> d13 = spb24*spb45; d13 = T(1)/d13;
complex<T> d20 = spb12*spb15*T(2); d20 = T(1)/d20;
complex<T> d21 = spb12*spb34*T(2); d21 = T(1)/d21;
complex<T> d22 = spb34*spb45*T(2); d22 = T(1)/d22;
complex<T> t16 = -t26; 
complex<T> t18 = -t29; 
complex<T> t24 = spb35*t13; 
complex<T> t37 = -(s23*t13); 
complex<T> t46 = d1*spa45; 
complex<T> t47 = t33*t34; 
complex<T> t50 = d4*t14; 
complex<T> t52 = spa15*t13; 
complex<T> d2 = spb12*spb34*t6*T(2); d2 = T(1)/d2;
complex<T> d6 = (s12 - s34)*(s12 + s15 - s34)*spb34*t6; d6 = T(1)/d6;
complex<T> d7 = t6*square(s12 - s34)*T(2); d7 = T(1)/d7;
complex<T> d8 = (s15 - s34)*t6; d8 = T(1)/d8;
complex<T> d9 = (s15 - s34)*spb34*t6; d9 = T(1)/d9;
complex<T> d10 = (s15 - s34)*(s12 + s15 - s34)*spb34*t6; d10 = T(1)/d10;
complex<T> d11 = spb34*t6; d11 = T(1)/d11;
complex<T> d12 = spb34*t3*t6; d12 = T(1)/d12;
complex<T> d14 = spb34*spb45*t3; d14 = T(1)/d14;
complex<T> d15 = spb24*t6; d15 = T(1)/d15;
complex<T> d16 = t3*t6; d16 = T(1)/d16;
complex<T> d17 = spb34*t6*T(2); d17 = T(1)/d17;
complex<T> d18 = spb12*t6*T(2); d18 = T(1)/d18;
complex<T> d19 = spb24*t6*T(2); d19 = T(1)/d19;
complex<T> d23 = spb34*spb45*t3*T(2); d23 = T(1)/d23;
complex<T> t30 = d2*t17; 
complex<T> t32 = d8*spa24; 
complex<T> t36 = t37*(d19*s34 + d18*t45); 
complex<T> t40 = d7*spb35; 
complex<T> t41 = spa15*t16; 
complex<T> t44 = t24*t26; 
complex<T> t49 = -t52; 
complex<T> t55 = t13*t18; 
complex<T> t1 = d6*t16*t24 + spb25*spb34*t14*t40 - d3*(spb13*spb23*t30 + t13*t46) + spb23*spb25*spb34*t50 - d5*t27*t47*T(2); 
complex<T> t2 = d3*spb13*spb23*t30 - t13*t32 - spb25*spb34*t14*t40 + d10*t44 + d6*t44 + d3*t13*t46 - spb23*spb25*spb34*t50 + d9*t55 + d5*t27*t47*T(2) + d2*spb13*t13*T(3); 
complex<T> t43 = d12*s12*t44 + d11*t55; 
complex<T> t48 = d10*t16*t24 + t13*(d9*t29 + t32); 
complex<T> t51 = d23*s12*t24*t41 + d22*t29*t52; 
complex<T> t54 = -(d15*s34*t13) + d16*spa34*t44; 
complex<T> t60 = d14*t24*t41 + d13*t49; 
complex<T> co1 = d15*t37; 
complex<T> co2 = d17*s23*t55; 
complex<T> co3 = d20*spa45*t13*t45; 
complex<T> co4 = d21*spa45*spb13*t52; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t1*Int(ep,mu,c12,c345) + t48*Int(ep,mu,c15,c234) + t2*Int(ep,mu,c34,c125) + t43*Int(ep,mu,c1,c2,c345) + t60*Int(ep,mu,c1,c5,c234) + co1*Int(ep,mu,c2,c3,c145) + t54*Int(ep,mu,c3,c4,c125) + co2*Int(ep,mu,c1,c2,c3,c45) + t36*Int(ep,mu,c2,c3,c4,c15) + co3*Int(ep,mu,c3,c4,c5,c12) + co4*Int(ep,mu,c4,c5,c1,c23) + t51*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqpmpm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpmpm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:32:18 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb13 = SPB(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s35 = -(spa35*spb35);
complex<T> t1 = spb23*T(2); 
complex<T> t3 = spb15*spb45; 
complex<T> t6 = square(s35); 
complex<T> t16 = square(spa35); 
complex<T> t17 = square(spb24); 
complex<T> t18 = square(spa13); 
complex<T> t19 = spb23*spb34; 
complex<T> t29 = cube(spa35); 
complex<T> t30 = cube(spb24); 
complex<T> t31 = square(spb14); 
complex<T> t34 = s12 - s45; 
complex<T> t35 = cube(spa13); 
complex<T> t36 = -(spb15*spb23*spb24) - spb13*spb24*spb25 - spb14*spb23*spb25*T(2); 
complex<T> t38 = s12 - s34; 
complex<T> t40 = square(spb35); 
complex<T> t42 = spb13*spb24 + spb14*spb23*T(2); 
complex<T> t43 = s12 + s15 - s34; 
complex<T> t47 = s35*T(3); 
complex<T> t48 = square(square(spa35)); 
complex<T> t49 = spb25*spb45; 
complex<T> t50 = spb12*spb34; 
complex<T> t57 = s25*spa25; 
complex<T> t64 = spa23*spb14; 
complex<T> t71 = s35*spb15; 
complex<T> d37 = square(square(spb35)); d37 = T(1)/d37;
complex<T> d38 = spb15*cube(spb35); d38 = T(1)/d38;
complex<T> d43 = spb12*square(square(spb35)); d43 = T(1)/d43;
complex<T> d44 = spb15*spb23*square(spb13); d44 = T(1)/d44;
complex<T> d45 = spb15*cube(spb13); d45 = T(1)/d45;
complex<T> d52 = spb15*cube(spb35)*T(2); d52 = T(1)/d52;
complex<T> t32 = -t49; 
complex<T> t33 = -t50; 
complex<T> t41 = -t57; 
complex<T> t46 = spa35*spb23*spb45 + spa13*t50; 
complex<T> t52 = t31*t35; 
complex<T> t59 = d43*T(2); 
complex<T> t68 = t19*t49; 
complex<T> t69 = -(t48*T(2)); 
complex<T> t70 = spb24*t18; 
complex<T> t72 = spb14*t42; 
complex<T> t78 = spb34*t29; 
complex<T> t79 = t3*T(2); 
complex<T> t88 = spb14*t17; 
complex<T> t89 = t48*T(2); 
complex<T> t91 = d38*spb34; 
complex<T> t96 = spb15*t19; 
complex<T> t103 = spb14*t30; 
complex<T> t106 = t31*t50; 
complex<T> t140 = t30*t64; 
complex<T> d1 = spb15*spb23*t34; d1 = T(1)/d1;
complex<T> d2 = spb12*spb23*t38*t71; d2 = T(1)/d2;
complex<T> d3 = spb12*spb23*t34*t71; d3 = T(1)/d3;
complex<T> d6 = spb12*cube(t38)*T(3); d6 = T(1)/d6;
complex<T> d7 = s13*spb23*t3*t34; d7 = T(1)/d7;
complex<T> d9 = t3*t34*square(s13); d9 = T(1)/d9;
complex<T> d10 = t1*t3*square(t34); d10 = T(1)/d10;
complex<T> d11 = spb15*square(t34)*T(2); d11 = T(1)/d11;
complex<T> d12 = spb12*cube(t34)*T(3); d12 = T(1)/d12;
complex<T> d13 = spb15*t38*t6; d13 = T(1)/d13;
complex<T> d14 = t71*square(t38)*T(2); d14 = T(1)/d14;
complex<T> d15 = t71*square(t34)*T(2); d15 = T(1)/d15;
complex<T> d16 = spb15*t34*t6; d16 = T(1)/d16;
complex<T> d17 = spb12*t38*cube(s35); d17 = T(1)/d17;
complex<T> d18 = spb12*t6*square(t38); d18 = T(1)/d18;
complex<T> d19 = spb12*t47*cube(t38); d19 = T(1)/d19;
complex<T> d20 = spb12*t47*cube(t34); d20 = T(1)/d20;
complex<T> d21 = spb12*t6*square(t34); d21 = T(1)/d21;
complex<T> d22 = spb12*t34*cube(s35); d22 = T(1)/d22;
complex<T> d23 = t1*t3*t34*t50; d23 = T(1)/d23;
complex<T> d27 = (s23 - s45)*spb23*t3; d27 = T(1)/d27;
complex<T> d28 = s13*(s23 - s45)*spb23*t3; d28 = T(1)/d28;
complex<T> d30 = (s23 - s45)*t3*square(s13); d30 = T(1)/d30;
complex<T> d31 = t1*t3*t50; d31 = T(1)/d31;
complex<T> d33 = spb15*spb23*t40; d33 = T(1)/d33;
complex<T> d34 = spb23*t3*square(spb13); d34 = T(1)/d34;
complex<T> d35 = t19*t3; d35 = T(1)/d35;
complex<T> d36 = t3*cube(spb13); d36 = T(1)/d36;
complex<T> d39 = t19*square(s15 + t38); d39 = T(1)/d39;
complex<T> d40 = t3*square(spb13); d40 = T(1)/d40;
complex<T> d41 = spb15*spb23*square(s15 + t38); d41 = T(1)/d41;
complex<T> d42 = spb12*spb15*spb23*t40; d42 = T(1)/d42;
complex<T> d50 = spb12*spb15*t1; d50 = T(1)/d50;
complex<T> d51 = spb12*spb15*t1*t40; d51 = T(1)/d51;
complex<T> d53 = t1*t50; d53 = T(1)/d53;
complex<T> d54 = spb34*t1*square(s15 + t38); d54 = T(1)/d54;
complex<T> d55 = spb34*spb45*t1; d55 = T(1)/d55;
complex<T> t8 = d36*s23*t106 + d40*spb24*t42*t64; 
complex<T> t11 = d50*spa34*spa45*t103 + s34*s45*(d51*spb24*t36 + d52*spb34*t49 + d43*t68); 
complex<T> t63 = -t70; 
complex<T> t84 = d6*spb25; 
complex<T> t85 = d42*t36; 
complex<T> t94 = spa15*t41; 
complex<T> t113 = d12*spb45; 
complex<T> t124 = d10*spb14; 
complex<T> t134 = spa15*t103; 
complex<T> d4 = t38*t96; d4 = T(1)/d4;
complex<T> d5 = t38*(s15 + t38)*t96; d5 = T(1)/d5;
complex<T> d8 = s13*t79*square(t34); d8 = T(1)/d8;
complex<T> d24 = (s15 - s34)*t96; d24 = T(1)/d24;
complex<T> d25 = (s15 - s34)*(s15 + t38)*t96; d25 = T(1)/d25;
complex<T> d26 = t79*square(s23 - s45); d26 = T(1)/d26;
complex<T> d29 = s13*t79*square(s23 - s45); d29 = T(1)/d29;
complex<T> d32 = t96*square(s15 + t38); d32 = T(1)/d32;
complex<T> d46 = t79*square(spb13); d46 = T(1)/d46;
complex<T> d47 = spb34*t79; d47 = T(1)/d47;
complex<T> d48 = t79*cube(spb13); d48 = T(1)/d48;
complex<T> d49 = spb12*t79; d49 = T(1)/d49;
complex<T> t9 = d41*spa34*t30*t57 + s34*(t59*t68 + spb24*t85 + t49*t91); 
complex<T> t10 = d45*spa45*t106 - d44*spa45*spb24*t72 + s45*(t59*t68 + spb24*t85 + t49*t91); 
complex<T> t12 = -(d35*spa12*t103) + d36*s12*t106 + d33*spa12*spb24*t36 + d32*s12*t30*t57 - d34*s12*spb24*t72 + s12*t32*t91 + d37*spa12*t68*T(2); 
complex<T> t13 = d24*spa25*t30 - d2*spb24*t16*t36 + d18*t19*t32*t48 + d25*t30*t57 + d5*t30*t57 + d17*t68*t69 + d19*t68*t69 + d13*t49*t78 + d14*t49*t78 - d4*spa15*t88 - spb24*t78*t84*T(2); 
complex<T> t14 = -(d11*spb24*spb45*t16) - d1*spa35*t17 + d2*spb24*t16*t36 + d3*spb24*t16*t36 + d5*t30*t41 + d9*t33*t52 + d8*t50*t52 + d18*t48*t68 + d21*t48*t68 + t124*t50*t70 + d7*t63*t72 + d13*t32*t78 + d14*t32*t78 + d15*t32*t78 + d16*t32*t78 + d4*spa15*t88 + d17*t68*t89 + d19*t68*t89 + d20*t68*t89 + d22*t68*t89 + spb23*spb24*t113*t29*T(2) + spb24*t78*t84*T(2) - d23*t46*t88*T(3); 
complex<T> t15 = d11*spb24*spb45*t16 + d1*spa35*t17 - d3*spb24*t16*t36 + d21*t19*t32*t48 + d29*t33*t52 + d8*t33*t52 + d30*t50*t52 + d9*t50*t52 + d26*t31*t63 + d20*t68*t69 + d22*t68*t69 + t124*t33*t70 + d28*t70*t72 + d7*t70*t72 + d15*t49*t78 + d16*t49*t78 + d27*spa13*t88 - spb23*spb24*t113*t29*T(2) + d31*t103*T(3) + d23*t46*t88*T(3); 
complex<T> t62 = d48*s12*s23*t106 + d47*spa12*t140 + d46*s12*spb24*(spb13*spb24 + spb14*t1)*t64; 
complex<T> t66 = -(d24*spa25*t30) + d25*t30*t41; 
complex<T> t67 = d30*t33*t52 + d29*t50*t52 + d26*t31*t70 + d28*t63*t72 - d27*spa13*t88; 
complex<T> t150 = t30*t94; 
complex<T> t77 = d55*spa12*t134 + d54*s12*t150; 
complex<T> co1 = d39*t150; 
complex<T> co2 = d49*spa34*t140; 
complex<T> co3 = d53*spa45*t134; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t14*Int(ep,mu,c12,c345) + t66*Int(ep,mu,c15,c234) + t67*Int(ep,mu,c23,c145) + t13*Int(ep,mu,c34,c125) + t15*Int(ep,mu,c45,c123) + t12*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + t8*Int(ep,mu,c2,c3,c145) + t9*Int(ep,mu,c3,c4,c125) + t10*Int(ep,mu,c4,c5,c123) + t62*Int(ep,mu,c1,c2,c3,c45) + co2*Int(ep,mu,c2,c3,c4,c15) + t11*Int(ep,mu,c3,c4,c5,c12) + co3*Int(ep,mu,c4,c5,c1,c23) + t77*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqpmmp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpmmp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:33:41 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa14 = SPA(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s15 = S(1,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> t5 = spb34*spb45; 
complex<T> t7 = spb23*T(2); 
complex<T> t13 = (s15 - s23)*spb23; 
complex<T> t20 = square(spb25); 
complex<T> t21 = square(spa14); 
complex<T> t22 = square(spa13); 
complex<T> t24 = square(spb35); 
complex<T> t26 = s23 - s45; 
complex<T> t29 = square(spa34); 
complex<T> t37 = cube(spa14); 
complex<T> t39 = cube(spa13); 
complex<T> t40 = spb24*spb45; 
complex<T> t41 = -(spb35*T(2)); 
complex<T> t46 = s15*s45; 
complex<T> t53 = cube(spb25); 
complex<T> t55 = square(spb15); 
complex<T> t56 = spb15*spb25; 
complex<T> t66 = spb12*spb15; 
complex<T> t69 = spb23*spb45; 
complex<T> t75 = spb35*T(2); 
complex<T> d1 = spb12*spb34; d1 = T(1)/d1;
complex<T> d3 = s12 - s45; d3 = T(1)/d3;
complex<T> d8 = spb12*spb34*square(s12 - s45)*T(6); d8 = T(1)/d8;
complex<T> d9 = spb34*cube(s12 - s45)*T(3); d9 = T(1)/d9;
complex<T> d28 = spb23*spb34*square(s15 - s23 + s45); d28 = T(1)/d28;
complex<T> d29 = (s15 - s23 + s45)*spb23*spb34; d29 = T(1)/d29;
complex<T> d30 = spb23*spb34*cube(s15 - s23 + s45); d30 = T(1)/d30;
complex<T> d31 = spb34*square(s15 - s23 + s45); d31 = T(1)/d31;
complex<T> d32 = (s15 - s23 + s45)*spb34; d32 = T(1)/d32;
complex<T> d34 = spb34*cube(s15 - s23 + s45); d34 = T(1)/d34;
complex<T> d35 = spb34*square(spb13); d35 = T(1)/d35;
complex<T> d36 = spb13*spb23*spb34; d36 = T(1)/d36;
complex<T> d37 = spb34*cube(spb13); d37 = T(1)/d37;
complex<T> d38 = spb12*spb45*T(2); d38 = T(1)/d38;
complex<T> t10 = spb34*t26; 
complex<T> t23 = -t66; 
complex<T> t31 = -t46; 
complex<T> t33 = -(d1*T(2)); 
complex<T> t49 = d9*t29; 
complex<T> t52 = s13*t5; 
complex<T> t57 = t37*t40; 
complex<T> t58 = t24*t39; 
complex<T> t59 = d1*T(2); 
complex<T> t60 = d30*spb12; 
complex<T> t67 = spb24*t21; 
complex<T> t77 = -(spa14*t20); 
complex<T> t82 = spb35*t20; 
complex<T> t85 = t41*t56; 
complex<T> t90 = spa13*t20; 
complex<T> t93 = -(d28*T(2)); 
complex<T> t94 = t29*t69; 
complex<T> t97 = spb15*t24; 
complex<T> t98 = spb25*t22; 
complex<T> t105 = -(s15*t53); 
complex<T> d2 = spb23*t5; d2 = T(1)/d2;
complex<T> d4 = spb12*spb23*t5; d4 = T(1)/d4;
complex<T> d7 = (s12 - s45)*t5*square(s13); d7 = T(1)/d7;
complex<T> d10 = spb34*t7*square(s15 - s23); d10 = T(1)/d10;
complex<T> d11 = (s15 - s23 + s45)*spb34*t13; d11 = T(1)/d11;
complex<T> d12 = spb34*t13*square(s15 - s23 + s45); d12 = T(1)/d12;
complex<T> d13 = (s15 - s23 + s45)*spb34*t7*square(s15 - s23); d13 = T(1)/d13;
complex<T> d16 = spb12*t5; d16 = T(1)/d16;
complex<T> d17 = t5*T(2); d17 = T(1)/d17;
complex<T> d18 = square(t26); d18 = T(1)/d18;
complex<T> d21 = t26*t5*square(s13); d21 = T(1)/d21;
complex<T> d23 = (s15 - s23 + s45)*spb34*t7*square(t26); d23 = T(1)/d23;
complex<T> d24 = spb12*t5*t7; d24 = T(1)/d24;
complex<T> d25 = t5*square(spb13); d25 = T(1)/d25;
complex<T> d26 = spb13*spb23*t5; d26 = T(1)/d26;
complex<T> d27 = t5*cube(spb13); d27 = T(1)/d27;
complex<T> d33 = spb13*t5; d33 = T(1)/d33;
complex<T> d39 = spb13*t5*T(2); d39 = T(1)/d39;
complex<T> d40 = t5*cube(spb13)*T(2); d40 = T(1)/d40;
complex<T> d41 = spb12*t7; d41 = T(1)/d41;
complex<T> d42 = (s15 - s23 + s45)*spb34*t7; d42 = T(1)/d42;
complex<T> d43 = spb12*spb34*t7; d43 = T(1)/d43;
complex<T> d44 = spb34*t7*cube(s15 - s23 + s45); d44 = T(1)/d44;
complex<T> d45 = t5*t7; d45 = T(1)/d45;
complex<T> t16 = s12*(-(d25*s23*spb35*t56) + d40*s23*t24*t66 + d39*spa23*t82); 
complex<T> t44 = d17*spa12; 
complex<T> t47 = d27*spb12; 
complex<T> t68 = -t82; 
complex<T> t73 = t56*t67; 
complex<T> t87 = t57*t60; 
complex<T> d5 = (s12 - s45)*t52; d5 = T(1)/d5;
complex<T> d6 = t52*square(s12 - s45)*T(2); d6 = T(1)/d6;
complex<T> d14 = (s15 - s23 + s45)*spb23*t10; d14 = T(1)/d14;
complex<T> d15 = t10*t7; d15 = T(1)/d15;
complex<T> d19 = t26*t52; d19 = T(1)/d19;
complex<T> d20 = t52*square(t26)*T(2); d20 = T(1)/d20;
complex<T> d22 = spb23*t10*square(s15 - s23 + s45); d22 = T(1)/d22;
complex<T> t1 = d4*t53 + d22*t23*t57 + d21*t23*t58 + d12*t57*t66 + d13*t57*t66 + d23*t57*t66 + d20*t58*t66 - d10*t73 + d15*t77 + d19*t22*t85 - d18*spb15*t44*t90 - d16*d18*spb23*t55*t98 - d11*t73*T(2) + d14*t73*T(2); 
complex<T> t2 = d15*spa14*t20 + d23*t23*t57 + d20*t23*t58 + d6*t23*t58 + d3*spa34*t20*t59 + d22*t57*t66 + d21*t58*t66 + d7*t58*t66 + d2*d3*spa13*t68 + spa13*t41*t49*t69 + d19*t22*t56*t75 + d5*t22*t56*t75 + d18*spb15*t44*t90 + d16*d18*spb23*t55*t98 - d14*t73*T(2) + d24*t53*T(3) - d8*spb25*t94*T(5); 
complex<T> t3 = d3*spa34*t20*t33 - d4*t53 + d7*t23*t58 + d6*t58*t66 + spa13*spb35*spb45*t49*t7 + d2*d3*spa13*t82 + d5*t22*t85 + d8*spb25*t94*T(5); 
complex<T> t15 = -(d2*spa12*t53) + s12*(d26*t68 + d25*t85 + t47*t97); 
complex<T> t17 = s15*(d29*t77 + spb15*t87 + t73*t93); 
complex<T> t18 = d37*spa45*t24*t66 + d36*spa45*t68 + d29*s45*t77 + d35*spa45*t85 + s45*spb15*t87 + s45*t73*t93; 
complex<T> t19 = d34*spa23*t57*t66 + d32*spa23*t77 + d33*spa23*t82 + d25*s23*t85 + s23*t47*t97 - d31*spa23*t73*T(2); 
complex<T> t51 = d43*spa45*t105 + d42*spa14*t20*t31 + d44*t46*t57*t66 + d28*t31*t73; 
complex<T> t65 = d12*t23*t57 + d13*t23*t57 + t73*(d10 + d11*T(2)); 
complex<T> co1 = spa23*t44*t53; 
complex<T> co2 = d38*spa23*spa34*t53; 
complex<T> co3 = d41*spa34*spa45*t53; 
complex<T> co4 = d45*spa12*t105; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t3*Int(ep,mu,c12,c345) + t65*Int(ep,mu,c15,c234) + t1*Int(ep,mu,c23,c145) + t2*Int(ep,mu,c45,c123) + t15*Int(ep,mu,c1,c2,c345) + t17*Int(ep,mu,c1,c5,c234) + t19*Int(ep,mu,c2,c3,c145) + t18*Int(ep,mu,c4,c5,c123) + co1*Int(ep,mu,c1,c2,c3,c45) + co2*Int(ep,mu,c2,c3,c4,c15) + t16*Int(ep,mu,c3,c2,c1,c45) + co3*Int(ep,mu,c3,c4,c5,c12) + t51*Int(ep,mu,c4,c5,c1,c23) + co4*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmqpmm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmqpmm L");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpqpmm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpqpmm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:33:42 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = S(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s23 = S(2,3);
complex<T> t5 = square(spb23); 
complex<T> t7 = square(spa45); 
complex<T> t8 = square(spb25); 
complex<T> t11 = spa45*spb23; 
complex<T> d1 = (s12 - s34)*spb15*spb45; d1 = T(1)/d1;
complex<T> d2 = spb15*spb34*spb45*T(2); d2 = T(1)/d2;
complex<T> d3 = spb15*spb45*square(s12 - s34)*T(2); d3 = T(1)/d3;
complex<T> d4 = spb34*spb45*T(2); d4 = T(1)/d4;
complex<T> d5 = spb15*spb45*T(2); d5 = T(1)/d5;
complex<T> d6 = spb15*T(2); d6 = T(1)/d6;
complex<T> d7 = spb34*T(2); d7 = T(1)/d7;
complex<T> t4 = d3*spb34*t7*t8 - d1*spb25*t11*T(2); 
complex<T> t15 = d2*t5; 
complex<T> t3 = -(d3*spb34*t7*t8) + d1*spb25*t11*T(2) - t15*T(3); 
complex<T> co1 = -(s12*s23*t15); 
complex<T> co2 = -(d4*s12*spa15*t5); 
complex<T> co3 = d5*s23*spa34*t5; 
complex<T> co4 = -(d6*spa34*spa45*t5); 
complex<T> co5 = -(d5*s23*spa34*t5); 
complex<T> co6 = -(d7*spa15*spa45*t5); 
complex<T> co7 = d4*s12*spa15*t5; 
complex<T> co8 = -co3; 
complex<T> co9 = Complex(0,1); 
SeriesC<T> result = co9*(t3*Int(ep,mu,c12,c345) + t4*Int(ep,mu,c34,c125) + co1*Int(ep,mu,c1,c2,c3,c45) + co2*Int(ep,mu,c2,c1,c5,c34) + co3*Int(ep,mu,c2,c3,c4,c15) + co4*Int(ep,mu,c3,c4,c5,c12) + co8*Int(ep,mu,c4,c3,c2,c15) + co6*Int(ep,mu,c4,c5,c1,c23) + co7*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmqppm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmqppm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:33:44 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> t1 = spb12*spb15; 
complex<T> t6 = square(spb34); 
complex<T> t9 = square(spa15); 
complex<T> t10 = square(spb13); 
complex<T> t11 = -(spa45*spb14); 
complex<T> t24 = spa15*spb13; 
complex<T> t27 = spb14*spb34; 
complex<T> d4 = (s12 + s15 - s34)*spb15*spb23; d4 = T(1)/d4;
complex<T> d5 = spb15*spb23*spb35; d5 = T(1)/d5;
complex<T> d6 = (s12 + s15 - s34)*spb23; d6 = T(1)/d6;
complex<T> d9 = spb15*spb45*T(2); d9 = T(1)/d9;
complex<T> d10 = (s12 + s15 - s34)*spb23*T(2); d10 = T(1)/d10;
complex<T> d13 = spb12*spb23*T(2); d13 = T(1)/d13;
complex<T> d14 = spb23*spb45*T(2); d14 = T(1)/d14;
complex<T> t12 = -(d4*spa25); 
complex<T> t14 = spb14*t6; 
complex<T> t22 = s34*t6; 
complex<T> t29 = t10*t9; 
complex<T> t32 = s12*t6; 
complex<T> d1 = (s23 - s45)*spb23*t1; d1 = T(1)/d1;
complex<T> d2 = spb23*spb45*t1*T(2); d2 = T(1)/d2;
complex<T> d3 = spb23*t1*square(s23 - s45)*T(2); d3 = T(1)/d3;
complex<T> d7 = spb23*spb35*t1; d7 = T(1)/d7;
complex<T> d8 = spb23*t1; d8 = T(1)/d8;
complex<T> d11 = spb45*t1*T(2); d11 = T(1)/d11;
complex<T> d12 = spb23*t1*T(2); d12 = T(1)/d12;
complex<T> d15 = spb23*spb35*t1*T(2); d15 = T(1)/d15;
complex<T> t15 = -(d1*T(2)); 
complex<T> t17 = -(d7*spb13); 
complex<T> t19 = d3*spb45; 
complex<T> t20 = t12*t32 - d5*spa12*spb13*t6; 
complex<T> t4 = (d4*spa25 + t17)*t22; 
complex<T> t5 = -(spb14*t19*t29) + d1*t24*t27*T(2); 
complex<T> t21 = t15*t24*t27 + spb14*t19*t29 + d2*t14*T(3); 
complex<T> t26 = (d8*t11 + s45*t17)*t6; 
complex<T> co1 = d6*spa15*spa25*t6; 
complex<T> co2 = d9*spa12*spa23*t14; 
complex<T> co3 = d10*spa15*spa25*t32; 
complex<T> co4 = -(d11*s34*spa23*t14); 
complex<T> co5 = -(d9*spa12*spa23*t14); 
complex<T> co6 = d12*t11*t22; 
complex<T> co7 = d13*spa15*spa45*t14; 
complex<T> co8 = d14*spa12*spa15*t14; 
complex<T> co9 = -(d15*s45*spb13*t22); 
complex<T> co10 = -co2; 
complex<T> co11 = Complex(0,1); 
SeriesC<T> result = co11*(t21*Int(ep,mu,c23,c145) + t5*Int(ep,mu,c45,c123) + t20*Int(ep,mu,c1,c2,c345) + co1*Int(ep,mu,c1,c5,c234) + t4*Int(ep,mu,c3,c4,c125) + t26*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c1,c2,c3,c45) + co3*Int(ep,mu,c2,c1,c5,c34) + co4*Int(ep,mu,c2,c3,c4,c15) + co10*Int(ep,mu,c3,c2,c1,c45) + co6*Int(ep,mu,c3,c4,c5,c12) + co7*Int(ep,mu,c4,c5,c1,c23) + co8*Int(ep,mu,c5,c1,c2,c34) + co9*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmqpmp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmqpmp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:34:46 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb24 = SPB(2,4);
complex<T> s15 = S(1,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s24 = -(spa24*spb24);
complex<T> t4 = square(s15 - s23); 
complex<T> t5 = spb12*T(2); 
complex<T> t7 = (s23 - s45)*spb23; 
complex<T> t20 = square(spa14); 
complex<T> t21 = square(spb35); 
complex<T> t22 = square(spa24); 
complex<T> t23 = square(spb13); 
complex<T> t25 = s15 - s23 + s45; 
complex<T> t26 = -(spb15*T(2)); 
complex<T> t27 = s15 - s34; 
complex<T> t32 = cube(spa14); 
complex<T> t33 = cube(spb35); 
complex<T> t34 = cube(spa24); 
complex<T> t35 = spb23*spb45; 
complex<T> t36 = -(spb25*T(2)); 
complex<T> t49 = spb35*T(2); 
complex<T> t53 = spb15*spb45; 
complex<T> t54 = spa24*spb35; 
complex<T> d3 = (s15 - s23)*s24*spb12; d3 = T(1)/d3;
complex<T> d9 = (s15 - s23)*spb12*square(s24); d9 = T(1)/d9;
complex<T> d16 = s15 - s23; d16 = T(1)/d16;
complex<T> d22 = spb12*square(spb24); d22 = T(1)/d22;
complex<T> d24 = spb12*spb24*spb34; d24 = T(1)/d24;
complex<T> d26 = spb12*cube(spb24); d26 = T(1)/d26;
complex<T> d30 = spb12*spb24; d30 = T(1)/d30;
complex<T> d31 = spb12*spb23*spb34; d31 = T(1)/d31;
complex<T> d32 = spb34*spb45*T(2); d32 = T(1)/d32;
complex<T> t24 = -t35; 
complex<T> t46 = -t53; 
complex<T> t48 = spb25*t34; 
complex<T> t52 = t23*t32; 
complex<T> t57 = spb13*t20; 
complex<T> t58 = spb35*t26; 
complex<T> t63 = spb25*t22; 
complex<T> t72 = spa14*t21; 
complex<T> t77 = spb35*t36; 
complex<T> t86 = spa45*t33; 
complex<T> d1 = spb23*t4*t5; d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*spb12*spb23*t25; d2 = T(1)/d2;
complex<T> d4 = t5*square(t27); d4 = T(1)/d4;
complex<T> d5 = s24*spb12*t27; d5 = T(1)/d5;
complex<T> d6 = spb12*spb34*t35; d6 = T(1)/d6;
complex<T> d7 = (s15 - s23)*spb12*spb23*square(t25); d7 = T(1)/d7;
complex<T> d8 = spb23*t25*t4*t5; d8 = T(1)/d8;
complex<T> d10 = s24*t4*t5; d10 = T(1)/d10;
complex<T> d11 = s24*t5*square(t27); d11 = T(1)/d11;
complex<T> d12 = spb12*t27*square(s24); d12 = T(1)/d12;
complex<T> d13 = spb24*t4*t5; d13 = T(1)/d13;
complex<T> d14 = spb12*t35; d14 = T(1)/d14;
complex<T> d15 = spb24*spb34*t5; d15 = T(1)/d15;
complex<T> d17 = spb12*t25*t7; d17 = T(1)/d17;
complex<T> d18 = spb34*t35*t5; d18 = T(1)/d18;
complex<T> d19 = spb12*t7*square(t25); d19 = T(1)/d19;
complex<T> d20 = spb23*t25*t5*square(s23 - s45); d20 = T(1)/d20;
complex<T> d21 = spb12*spb23*square(t25); d21 = T(1)/d21;
complex<T> d23 = spb12*spb23*spb34*t25; d23 = T(1)/d23;
complex<T> d25 = spb12*spb23*cube(t25); d25 = T(1)/d25;
complex<T> d27 = spb12*square(t25); d27 = T(1)/d27;
complex<T> d28 = spb12*spb34*t25; d28 = T(1)/d28;
complex<T> d29 = spb12*cube(t25); d29 = T(1)/d29;
complex<T> d33 = spb45*t5; d33 = T(1)/d33;
complex<T> d34 = spb23*t5; d34 = T(1)/d34;
complex<T> d35 = spb24*t5; d35 = T(1)/d35;
complex<T> d36 = t5*cube(spb24); d36 = T(1)/d36;
complex<T> d37 = spb23*spb34*t25*t5; d37 = T(1)/d37;
complex<T> d38 = spb23*spb34*t5; d38 = T(1)/d38;
complex<T> d39 = spb23*t5*cube(t25); d39 = T(1)/d39;
complex<T> d40 = spb34*t35*T(2); d40 = T(1)/d40;
complex<T> t16 = d30*spa34*t21 + d26*s34*spb25*t35 + d22*s34*t77; 
complex<T> t17 = s23*(-(d22*s34*spb25*spb35) + d35*spa34*t21 + d36*s34*spb25*t35); 
complex<T> t18 = d11*t35*t48 + d12*t35*t48 + d4*spb35*t63 + d5*t49*t63; 
complex<T> t43 = d13*spa34; 
complex<T> t62 = t52*t53; 
complex<T> t66 = -t72; 
complex<T> t70 = t57*t58; 
complex<T> t73 = d23*spb13; 
complex<T> t79 = spb15*t57; 
complex<T> t1 = -(d18*t33) + d10*t35*t48 + d9*t35*t48 + d19*t46*t52 + d15*d16*t24*t54 + t24*t43*t54 + d20*t62 + d7*t62 + d8*t62 + d3*t49*t63 + d2*t70 - d1*spb35*t79 + d17*t49*t79 + d14*d16*spb15*t72*T(2); 
complex<T> t2 = d10*t24*t48 + d11*t24*t48 + d12*t24*t48 + d9*t24*t48 + d7*t46*t52 + d8*t46*t52 + d15*d16*t35*t54 + t35*t43*t54 - d4*spb35*t63 + d14*d16*t26*t72 + d3*t22*t77 + d5*t22*t77 + d1*spb35*t79 + d2*t49*t79 + d6*t33*T(2); 
complex<T> t12 = d20*t46*t52 + d19*t62 + d17*t70; 
complex<T> t13 = d25*s45*t62 + d21*s45*t70 + s45*t66*t73 - d31*t86; 
complex<T> t14 = s15*(d39*s45*t62 + d37*s45*spb13*t66 - d21*s45*spb35*t79 - d38*t86); 
complex<T> t15 = -(d24*s23*t21) + d26*s23*spb25*t35 + d29*spa23*t62 + d28*spa23*spb13*t66 + d27*spa23*t70 + d22*s23*t77; 
complex<T> t44 = s15*(d24*t21 + d26*spb25*t24 + d22*spb25*t49 + d25*t62 + d21*t70 + t66*t73); 
complex<T> co1 = d32*spa12*spa23*t33; 
complex<T> co2 = d33*spa23*spa34*t33; 
complex<T> co3 = -(d32*spa12*spa23*t33); 
complex<T> co4 = d34*spa34*t86; 
complex<T> co5 = -(d40*s15*spa12*t33); 
complex<T> co6 = -co1; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t2*Int(ep,mu,c15,c234) + t1*Int(ep,mu,c23,c145) + t18*Int(ep,mu,c34,c125) + t12*Int(ep,mu,c45,c123) + t44*Int(ep,mu,c1,c5,c234) + t15*Int(ep,mu,c2,c3,c145) + t16*Int(ep,mu,c3,c4,c125) + t13*Int(ep,mu,c4,c5,c123) + co1*Int(ep,mu,c1,c2,c3,c45) + co2*Int(ep,mu,c2,c3,c4,c15) + co6*Int(ep,mu,c3,c2,c1,c45) + co4*Int(ep,mu,c3,c4,c5,c12) + t17*Int(ep,mu,c4,c3,c2,c15) + t14*Int(ep,mu,c4,c5,c1,c23) + co5*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmmqpm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, m, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmmqpm L");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpmqpm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, m, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpmqpm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:35:16 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = S(1,2);
complex<T> s35 = -(spa35*spb35);
complex<T> t3 = spb15*spb45; 
complex<T> t8 = square(spb24); 
complex<T> t21 = square(spa13); 
complex<T> t22 = square(spa35); 
complex<T> t23 = square(spb14); 
complex<T> t24 = -(spb23*spb45); 
complex<T> t26 = -(spb12*spb24); 
complex<T> t29 = s23 - s45; 
complex<T> t30 = spb15*T(2); 
complex<T> t31 = s12*s23; 
complex<T> t35 = cube(spa13); 
complex<T> t36 = cube(spa35); 
complex<T> t37 = spb12*spb23; 
complex<T> t38 = spb14*spb24; 
complex<T> t39 = spb25*spb45; 
complex<T> t48 = -square(spb24); 
complex<T> t52 = s34*s45; 
complex<T> t58 = cube(spb24); 
complex<T> t61 = spb24*spb25; 
complex<T> t62 = spa35*spb45; 
complex<T> t81 = spa13*spb12; 
complex<T> d2 = (s12 - s34)*s35*spb15; d2 = T(1)/d2;
complex<T> d3 = s35*(s12 - s45)*spb15; d3 = T(1)/d3;
complex<T> d9 = (s12 - s34)*spb15*square(s35); d9 = T(1)/d9;
complex<T> d12 = (s12 - s45)*spb15*square(s35); d12 = T(1)/d12;
complex<T> d17 = s12 - s45; d17 = T(1)/d17;
complex<T> d22 = spb15*square(spb35); d22 = T(1)/d22;
complex<T> d23 = spb15*spb34*spb35; d23 = T(1)/d23;
complex<T> d27 = spb15*cube(spb35); d27 = T(1)/d27;
complex<T> d29 = spb15*spb35; d29 = T(1)/d29;
complex<T> d30 = spb15*cube(spb13); d30 = T(1)/d30;
complex<T> d31 = spb15*square(spb13); d31 = T(1)/d31;
complex<T> d32 = spb13*spb15*spb34; d32 = T(1)/d32;
complex<T> d38 = spb23*spb34*T(2); d38 = T(1)/d38;
complex<T> t9 = t3*T(2); 
complex<T> t32 = -(spb23*t62) - spb34*t81; 
complex<T> t60 = t23*t37; 
complex<T> t66 = spb12*t38; 
complex<T> t67 = spb23*t39; 
complex<T> t71 = spb25*t24; 
complex<T> d1 = t30*square(s12 - s34); d1 = T(1)/d1;
complex<T> d5 = (s12 - s45)*t3*square(s13); d5 = T(1)/d5;
complex<T> d7 = s13*(s12 - s45)*t3; d7 = T(1)/d7;
complex<T> d10 = s35*t30*square(s12 - s34); d10 = T(1)/d10;
complex<T> d11 = s35*t30*square(s12 - s45); d11 = T(1)/d11;
complex<T> d13 = spb35*t30*square(s12 - s45); d13 = T(1)/d13;
complex<T> d15 = spb23*t3; d15 = T(1)/d15;
complex<T> d16 = spb34*spb35*t30; d16 = T(1)/d16;
complex<T> d19 = t29*t3*square(s13); d19 = T(1)/d19;
complex<T> d20 = s13*t29*t3; d20 = T(1)/d20;
complex<T> d21 = spb23*spb34*t3; d21 = T(1)/d21;
complex<T> d24 = t3*cube(spb13); d24 = T(1)/d24;
complex<T> d25 = t3*square(spb13); d25 = T(1)/d25;
complex<T> d26 = spb13*spb34*t3; d26 = T(1)/d26;
complex<T> d28 = spb34*t3; d28 = T(1)/d28;
complex<T> d36 = spb35*t30; d36 = T(1)/d36;
complex<T> d37 = t30*cube(spb35); d37 = T(1)/d37;
complex<T> d39 = t30*square(spb35); d39 = T(1)/d39;
complex<T> t17 = d22*s34*t61 + d27*s34*t67 + d29*spa34*square(spb24); 
complex<T> t18 = d23*s45*t48 + d32*spa45*spb14*t48 + d30*spa45*t60 + d22*s45*t61 + d31*spa45*t66 + d27*s45*t67; 
complex<T> t19 = -(d22*t52*t61) + d37*t52*t67 + d36*s45*spa34*square(spb24); 
complex<T> t20 = d1*t22*t61 - d2*t22*t61 + (d10 + d9)*t36*t67; 
complex<T> t44 = -(d15*T(2)); 
complex<T> t45 = d13*spa34; 
complex<T> t50 = -t60; 
complex<T> t54 = d15*T(2); 
complex<T> t78 = d26*spb14; 
complex<T> d4 = s13*t9*square(s12 - s45); d4 = T(1)/d4;
complex<T> d6 = t9*square(s12 - s45); d6 = T(1)/d6;
complex<T> d8 = spb23*spb34*t9; d8 = T(1)/d8;
complex<T> d14 = (s12 - s45)*spb23*spb34*t9; d14 = T(1)/d14;
complex<T> d18 = s13*t9*square(t29); d18 = T(1)/d18;
complex<T> d33 = t9*square(spb13); d33 = T(1)/d33;
complex<T> d34 = t9*cube(spb13); d34 = T(1)/d34;
complex<T> d35 = spb13*spb34*t9; d35 = T(1)/d35;
complex<T> t1 = d16*d17*spa35*spb24*t24 + d20*spb14*t21*t26 + d6*spb14*t21*t26 + d7*spb14*t21*t26 + spa35*spb24*t24*t45 + d18*t35*t50 + d4*t35*t50 + d19*t35*t60 + d5*t35*t60 - d3*t22*t61 + d11*t36*t67 + d12*t36*t67 + d17*t54*t8*t81 - d21*t58*T(2) + d14*t32*t8*T(3); 
complex<T> t2 = d5*t35*t50 + d8*t58 + d4*t35*t60 - d1*t22*t61 + d2*t22*t61 + d3*t22*t61 + d16*d17*spb23*spb24*t62 + spb23*spb24*t45*t62 + d6*t21*t66 + d7*t21*t66 + d10*t36*t71 + d11*t36*t71 + d12*t36*t71 + d9*t36*t71 + d17*t44*t8*t81 - d14*t32*t8*T(3); 
complex<T> t13 = d19*t35*t50 + d18*t35*t60 + d20*t21*t66; 
complex<T> t14 = -(d28*spa23*t58) + s23*(d24*t60 + d25*t66 + t48*t78); 
complex<T> t16 = s12*(d24*t60 - d22*t61 + d25*t66 + d27*t71 + t48*t78 + d23*square(spb24)); 
complex<T> t59 = t31*(d25*spb14*t26 + d34*t60 - d35*spb14*t8); 
complex<T> co1 = d33*t31*t66*T(3); 
complex<T> co2 = -(d38*spa15*spa45*t58); 
complex<T> co3 = d39*t52*t61*T(3); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t2*Int(ep,mu,c12,c345) + t13*Int(ep,mu,c23,c145) + t20*Int(ep,mu,c34,c125) + t1*Int(ep,mu,c45,c123) + t16*Int(ep,mu,c1,c2,c345) + t14*Int(ep,mu,c2,c3,c145) + t17*Int(ep,mu,c3,c4,c125) + t18*Int(ep,mu,c4,c5,c123) + co1*Int(ep,mu,c1,c2,c3,c45) + t59*Int(ep,mu,c3,c2,c1,c45) + t19*Int(ep,mu,c3,c4,c5,c12) + co2*Int(ep,mu,c4,c5,c1,c23) + co3*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmpqpm_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, p, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmpqpm L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:35:19 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> s15 = -(spa15*spb15);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = spb12*spb45; 
complex<T> t3 = spb15*T(2); 
complex<T> t5 = spb13*spb34; 
complex<T> t8 = square(spb34); 
complex<T> t11 = square(spa12); 
complex<T> t12 = square(spb14); 
complex<T> t28 = spa12*spb14; 
complex<T> d5 = (s12 + s15 - s34)*spb45; d5 = T(1)/d5;
complex<T> d11 = spb12*spb23*T(2); d11 = T(1)/d11;
complex<T> d12 = (s12 + s15 - s34)*spb45*T(2); d12 = T(1)/d12;
complex<T> t13 = -(spa15*spb13*spb45) - spb23*t28; 
complex<T> t16 = -(spb14*t8); 
complex<T> t20 = spb23*t11; 
complex<T> t25 = -(spb13*t8); 
complex<T> t33 = s15*t8; 
complex<T> d1 = t1*t3*square(s23 - s45); d1 = T(1)/d1;
complex<T> d2 = (s23 - s45)*spb15*t1; d2 = T(1)/d2;
complex<T> d3 = spb23*t1*t3; d3 = T(1)/d3;
complex<T> d4 = (s23 - s45)*spb23*t1*t3; d4 = T(1)/d4;
complex<T> d6 = (s12 + s15 - s34)*t1; d6 = T(1)/d6;
complex<T> d7 = spb24*t1; d7 = T(1)/d7;
complex<T> d8 = spb15*t1; d8 = T(1)/d8;
complex<T> d9 = spb15*spb24*t1; d9 = T(1)/d9;
complex<T> d10 = spb24*t1*t3; d10 = T(1)/d10;
complex<T> t18 = d6*spa25; 
complex<T> t19 = d1*spb13; 
complex<T> t21 = t13*T(3); 
complex<T> t23 = d9*s23*t16 + d8*spa23*t25; 
complex<T> t4 = -(t12*t19*t20) - d2*t28*t5*T(2) - d4*t13*t5*T(3) - d3*spb13*t8*T(3); 
complex<T> t14 = -t18; 
complex<T> t22 = t12*t19*t20 + d4*t21*t5 + d2*t28*t5*T(2); 
complex<T> t30 = s34*(d9*t16 + t18*t8); 
complex<T> t27 = d7*spa15*t16 + t14*t33; 
complex<T> co1 = d5*spa12*spa25*t8; 
complex<T> co2 = d10*s23*s34*t16; 
complex<T> co3 = d11*spa15*spa45*t25; 
complex<T> co4 = d12*spa12*spa25*t33; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t4*Int(ep,mu,c23,c145) + t22*Int(ep,mu,c45,c123) + co1*Int(ep,mu,c1,c2,c345) + t27*Int(ep,mu,c1,c5,c234) + t23*Int(ep,mu,c2,c3,c145) + t30*Int(ep,mu,c3,c4,c125) + co2*Int(ep,mu,c2,c3,c4,c15) + co3*Int(ep,mu,c4,c5,c1,c23) + co4*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmmqpp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, m, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmmqpp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:35:20 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> s45 = S(4,5);
complex<T> s15 = S(1,5);
complex<T> s34 = -(spa34*spb34);
complex<T> t7 = square(spb45); 
complex<T> t9 = square(spa23); 
complex<T> t10 = square(spb25); 
complex<T> t13 = spa23*spb25; 
complex<T> d1 = spb12*spb23*square(s15 - s34)*T(2); d1 = T(1)/d1;
complex<T> d2 = (s15 - s34)*spb12*spb23; d2 = T(1)/d2;
complex<T> d3 = (s15 - s34)*spb12*spb23*spb34*T(2); d3 = T(1)/d3;
complex<T> d4 = spb12*spb23*spb34*T(2); d4 = T(1)/d4;
complex<T> d5 = spb12*spb23*T(2); d5 = T(1)/d5;
complex<T> d6 = spb23*spb34*T(2); d6 = T(1)/d6;
complex<T> t11 = spa12*spb15*spb24 + spb34*t13; 
complex<T> t3 = -(d1*spb34*t10*t9) + d2*spb45*t13*T(2) - d3*spb45*t11*T(3); 
complex<T> t4 = d1*spb34*t10*t9 - d2*spb45*t13*T(2) + d3*spb45*t11*T(3) + d4*t7*T(3); 
complex<T> co1 = -(d5*s45*spa34*t7); 
complex<T> co2 = -(d6*s15*spa12*t7); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*Int(ep,mu,c15,c234) + t4*Int(ep,mu,c34,c125) + co1*Int(ep,mu,c3,c4,c5,c12) + co2*Int(ep,mu,c5,c1,c2,c34));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmmmqp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, m, m, qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmmmqp L");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpmmqp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, m, m, qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpmmqp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:36:12 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> s15 = -(spa15*spb15);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s13 = -(spa13*spb13);
complex<T> s12 = S(1,2);
complex<T> t5 = spb34*spb45; 
complex<T> t9 = square(spb25); 
complex<T> t12 = square(s13); 
complex<T> t21 = square(spa14); 
complex<T> t22 = square(spa13); 
complex<T> t24 = square(spb24); 
complex<T> t25 = s23 - s45; 
complex<T> t26 = s15 - s23 + s45; 
complex<T> t29 = -square(spb25); 
complex<T> t32 = square(spa34); 
complex<T> t36 = cube(spa14); 
complex<T> t38 = cube(spa13); 
complex<T> t39 = spb23*spb35; 
complex<T> t44 = s15*spb24; 
complex<T> t53 = square(spb12); 
complex<T> t56 = spb12*spb15; 
complex<T> t65 = spa14*spb24; 
complex<T> t68 = spb12*spb25; 
complex<T> t87 = spb23*spb25; 
complex<T> t97 = spa14*spb12; 
complex<T> d8 = (s15 - s23)*spb15*spb34; d8 = T(1)/d8;
complex<T> d9 = spb15*spb34; d9 = T(1)/d9;
complex<T> d11 = s15 - s23; d11 = T(1)/d11;
complex<T> d13 = spb15*spb34*square(s15 - s23)*T(2); d13 = T(1)/d13;
complex<T> d22 = spb23*spb34*T(2); d22 = T(1)/d22;
complex<T> d23 = spb15*spb23*spb34; d23 = T(1)/d23;
complex<T> d35 = spb13*spb34; d35 = T(1)/d35;
complex<T> d37 = spb34*cube(spb13); d37 = T(1)/d37;
complex<T> d38 = spb34*square(spb13); d38 = T(1)/d38;
complex<T> t23 = -t56; 
complex<T> t30 = t5*T(2); 
complex<T> t33 = -(spa13*spb15*spb23) - spb45*t97; 
complex<T> t37 = -t68; 
complex<T> t45 = -(d22*spa15); 
complex<T> t52 = -(spb25*t21); 
complex<T> t57 = t24*t36; 
complex<T> t58 = t38*t39; 
complex<T> t61 = d23*spb45; 
complex<T> t62 = -(d9*T(2)); 
complex<T> t69 = spb24*t21; 
complex<T> t70 = spb35*t22; 
complex<T> t72 = d9*T(2); 
complex<T> d2 = (s12 - s45)*t12*t5; d2 = T(1)/d2;
complex<T> d4 = s13*(s12 - s45)*t5; d4 = T(1)/d4;
complex<T> d5 = (s15 - s23)*spb23*spb34*square(t26); d5 = T(1)/d5;
complex<T> d6 = spb23*spb34*t26*square(s15 - s23)*T(2); d6 = T(1)/d6;
complex<T> d7 = (s15 - s23)*spb23*spb34*t26; d7 = T(1)/d7;
complex<T> d10 = spb23*t5; d10 = T(1)/d10;
complex<T> d12 = spb15*spb23*t5; d12 = T(1)/d12;
complex<T> d14 = spb23*spb34*t25*square(t26); d14 = T(1)/d14;
complex<T> d15 = spb23*spb34*t26*square(t25)*T(2); d15 = T(1)/d15;
complex<T> d16 = spb23*spb34*t25*t26; d16 = T(1)/d16;
complex<T> d19 = t12*t25*t5; d19 = T(1)/d19;
complex<T> d20 = s13*t25*t5; d20 = T(1)/d20;
complex<T> d24 = square(t25); d24 = T(1)/d24;
complex<T> d26 = spb13*t5; d26 = T(1)/d26;
complex<T> d27 = t5*cube(spb13); d27 = T(1)/d27;
complex<T> d28 = t5*square(spb13); d28 = T(1)/d28;
complex<T> d29 = spb23*spb34*cube(t26); d29 = T(1)/d29;
complex<T> d30 = spb23*spb34*square(t26); d30 = T(1)/d30;
complex<T> d31 = spb23*t26*t5; d31 = T(1)/d31;
complex<T> d32 = spb34*cube(t26); d32 = T(1)/d32;
complex<T> d33 = spb34*square(t26); d33 = T(1)/d33;
complex<T> d34 = t26*t5; d34 = T(1)/d34;
complex<T> d36 = spb23*spb34*t26; d36 = T(1)/d36;
complex<T> d42 = spb23*spb34*cube(t26)*T(2); d42 = T(1)/d42;
complex<T> d43 = spb23*spb34*square(t26)*T(2); d43 = T(1)/d43;
complex<T> d44 = spb23*spb34*t26*T(2); d44 = T(1)/d44;
complex<T> t1 = d5*t23*t57 + d6*t23*t57 + d7*t37*t69 + d13*spb45*t32*t87 + d8*spa34*t9 + d11*spa34*t62*t9 + d10*d11*t65*t9 - d12*cube(spb25); 
complex<T> t47 = d27*spb15; 
complex<T> t75 = t56*t57; 
complex<T> t78 = d28*spb35; 
complex<T> d1 = s13*t30*square(s12 - s45); d1 = T(1)/d1;
complex<T> d3 = t30*square(s12 - s45); d3 = T(1)/d3;
complex<T> d17 = t25*t30; d17 = T(1)/d17;
complex<T> d18 = s13*t30*square(t25); d18 = T(1)/d18;
complex<T> d21 = spb15*spb23*t25*t30; d21 = T(1)/d21;
complex<T> d25 = spb15*spb23*t30; d25 = T(1)/d25;
complex<T> d39 = spb13*t30; d39 = T(1)/d39;
complex<T> d40 = t30*cube(spb13); d40 = T(1)/d40;
complex<T> d41 = t30*square(spb13); d41 = T(1)/d41;
complex<T> t2 = d15*t23*t57 + d1*t23*t58 + d18*t23*t58 + d19*t56*t58 + d2*t56*t58 + d24*t52*t53*t61 + d16*t68*t69 + d20*t37*t70 + d3*t37*t70 + d4*t37*t70 + d14*t75 + d17*spa13*t9 + d22*d24*spa15*t29*t97 - d25*cube(spb25) + d21*t33*t9*T(3); 
complex<T> t3 = d17*spa13*t29 + d8*spa34*t29 + d14*t23*t57 + d19*t23*t58 + d18*t56*t58 + d24*spb25*t21*t53*t61 + d10*d11*t29*t65 + d16*t37*t69 + d7*t68*t69 + d20*t68*t70 + d15*t75 + d5*t75 + d6*t75 - d13*spb45*t32*t87 + d11*spa34*t72*t9 + d22*d24*spa15*t9*t97 - d21*t33*t9*T(3); 
complex<T> t14 = d43*s45*t21*t44*t68 + d42*s15*s45*t75 + d44*spa14*spa45*t44*square(spb25); 
complex<T> t15 = d31*spa14*t29*t44 + d30*t21*t44*t68 + d29*s15*t75 - d10*spa15*cube(spb25); 
complex<T> t18 = d2*t23*t58 + d1*t56*t58 + (d3 + d4)*t68*t70; 
complex<T> t19 = d35*spa45*t29 + d37*spa45*t39*t56 + d38*spa45*spb35*t68 + d30*s45*t68*t69 + d29*s45*t75 + d36*spa45*t65*square(spb25); 
complex<T> t35 = s12*s23*(d40*t39*t56 + d41*spb35*t68 - d39*t9); 
complex<T> t95 = t39*t47; 
complex<T> t16 = s12*(d26*t29 + t68*t78 + spb12*t95); 
complex<T> t17 = d26*s23*t29 + d34*spa23*t29*t65 + d33*spa23*t68*t69 + d32*spa23*t75 + s23*t68*t78 + s23*spb12*t95; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t18*Int(ep,mu,c12,c345) + t1*Int(ep,mu,c15,c234) + t3*Int(ep,mu,c23,c145) + t2*Int(ep,mu,c45,c123) + t16*Int(ep,mu,c1,c2,c345) + t15*Int(ep,mu,c1,c5,c234) + t17*Int(ep,mu,c2,c3,c145) + t19*Int(ep,mu,c4,c5,c123) + t35*Int(ep,mu,c3,c2,c1,c45) + t14*Int(ep,mu,c4,c5,c1,c23));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmpmqp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, p, m, qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmpmqp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:37:59 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s15 = -(spa15*spb15);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> t2 = spb45*T(2); 
complex<T> t3 = (s15 - s23)*spb12; 
complex<T> t5 = (s23 - s45)*spb23; 
complex<T> t15 = square(spa14); 
complex<T> t16 = square(spa24); 
complex<T> t17 = -(spb13*spb35); 
complex<T> t18 = spb12*spb45; 
complex<T> t23 = s23*s34; 
complex<T> t27 = spb23*T(2); 
complex<T> t28 = square(spb35); 
complex<T> t29 = square(spb13); 
complex<T> t30 = cube(spa14); 
complex<T> t32 = s15 - s23 + s45; 
complex<T> t33 = cube(spa24); 
complex<T> t36 = -(spb23*spb25); 
complex<T> t37 = s15 - s34; 
complex<T> t41 = spb12*T(2); 
complex<T> t47 = s15*s45; 
complex<T> t49 = cube(spb35); 
complex<T> t50 = spb15*spb34; 
complex<T> t52 = -(spb25*T(3)); 
complex<T> t54 = s25*spa25; 
complex<T> t64 = spb13*spb25; 
complex<T> t65 = spb23*spb34; 
complex<T> t66 = -(spb45*T(2)); 
complex<T> t85 = spb12*spb23; 
complex<T> t87 = -(spb14*spb25); 
complex<T> d13 = (s15 - s23)*s24*spb15; d13 = T(1)/d13;
complex<T> d31 = spb12*cube(spb24); d31 = T(1)/d31;
complex<T> d33 = spb12*square(spb24); d33 = T(1)/d33;
complex<T> d34 = square(spb24); d34 = T(1)/d34;
complex<T> d41 = spb12*spb15*square(spb24); d41 = T(1)/d41;
complex<T> d42 = spb15*square(spb24); d42 = T(1)/d42;
complex<T> d51 = spb15*square(spb24)*T(2); d51 = T(1)/d51;
complex<T> t1 = square(t32); 
complex<T> t4 = square(s12 + t37); 
complex<T> t31 = -(spb14*spb35) + spb13*t66; 
complex<T> t34 = -t50; 
complex<T> t35 = t64*t66 + spb35*(-t18 + t87); 
complex<T> t38 = -t54; 
complex<T> t51 = t29*t30; 
complex<T> t55 = spb12*t27; 
complex<T> t75 = spb35*t15; 
complex<T> t76 = spb25*t65; 
complex<T> t89 = spb34*t36; 
complex<T> t90 = d41*T(3); 
complex<T> t93 = spb35*t16; 
complex<T> t97 = t15*t17; 
complex<T> d1 = (s12 - s34)*spb34*t18; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spb34*t18*(s12 + t37); d2 = T(1)/d2;
complex<T> d5 = t3*square(s24); d5 = T(1)/d5;
complex<T> d6 = s24*t41*square(s15 - s23); d6 = T(1)/d6;
complex<T> d7 = s24*t41*square(t37); d7 = T(1)/d7;
complex<T> d8 = spb12*t37*square(s24); d8 = T(1)/d8;
complex<T> d9 = spb23*t3*t32; d9 = T(1)/d9;
complex<T> d10 = t41*square(s15 - s23); d10 = T(1)/d10;
complex<T> d11 = s24*spb15*t3; d11 = T(1)/d11;
complex<T> d12 = s24*spb12*spb15*t37; d12 = T(1)/d12;
complex<T> d14 = s24*spb15*t37; d14 = T(1)/d14;
complex<T> d15 = t2*t85*square(s15 - s23); d15 = T(1)/d15;
complex<T> d16 = spb45*t3; d16 = T(1)/d16;
complex<T> d17 = spb34*t18*t37; d17 = T(1)/d17;
complex<T> d18 = spb34*t18*t37*(s12 + t37); d18 = T(1)/d18;
complex<T> d19 = t2*t50*t85; d19 = T(1)/d19;
complex<T> d20 = spb23*spb45*t3*t32; d20 = T(1)/d20;
complex<T> d21 = s24*spb15*spb45*t3; d21 = T(1)/d21;
complex<T> d22 = s24*spb15*t18*t37; d22 = T(1)/d22;
complex<T> d26 = spb12*t32*t5; d26 = T(1)/d26;
complex<T> d27 = t18*t5; d27 = T(1)/d27;
complex<T> d28 = t18*t32*t5; d28 = T(1)/d28;
complex<T> d30 = t85*cube(t32); d30 = T(1)/d30;
complex<T> d36 = t18*t65; d36 = T(1)/d36;
complex<T> d38 = t18*square(spb24); d38 = T(1)/d38;
complex<T> d39 = spb12*cube(t32); d39 = T(1)/d39;
complex<T> d44 = spb15*t18*square(spb24); d44 = T(1)/d44;
complex<T> d49 = t41*cube(spb24); d49 = T(1)/d49;
complex<T> d50 = spb15*t41*square(spb24); d50 = T(1)/d50;
complex<T> d52 = spb12*spb15*t2*square(spb24); d52 = T(1)/d52;
complex<T> t10 = t23*(d51*t28 + d49*t76 + d52*spb35*(t64*t66 + spb35*(-t18 + t87)) + d50*spb35*t64*T(3)); 
complex<T> t13 = -(d17*spa12*spb13*t28) - d14*t16*t28 + d1*spa25*t49 + d18*t49*t54 + d2*t49*t54 + d7*t33*t76 + d8*t33*t76 + d12*spb13*t52*t93 - d22*(t64*t66 + spb35*(-t18 + t87))*t93; 
complex<T> t83 = t50*t51; 
complex<T> t88 = t29*t75; 
complex<T> t111 = t38*t49; 
complex<T> d3 = spb23*t1*t3; d3 = T(1)/d3;
complex<T> d4 = t32*t55*square(s15 - s23); d4 = T(1)/d4;
complex<T> d23 = spb12*t1*t5; d23 = T(1)/d23;
complex<T> d24 = t32*t55*square(s23 - s45); d24 = T(1)/d24;
complex<T> d25 = t55*square(s23 - s45); d25 = T(1)/d25;
complex<T> d29 = spb34*spb45*t4; d29 = T(1)/d29;
complex<T> d32 = t1*t85; d32 = T(1)/d32;
complex<T> d35 = spb34*t18*t4; d35 = T(1)/d35;
complex<T> d37 = spb23*t1*t18; d37 = T(1)/d37;
complex<T> d40 = spb12*t1; d40 = T(1)/d40;
complex<T> d43 = t1*t18; d43 = T(1)/d43;
complex<T> d45 = t18*t4; d45 = T(1)/d45;
complex<T> d46 = t55*cube(t32); d46 = T(1)/d46;
complex<T> d47 = t1*t55; d47 = T(1)/d47;
complex<T> d48 = spb34*t2*t4; d48 = T(1)/d48;
complex<T> t7 = d27*spa14*spb13*t28 + d24*t34*t51 + d28*spb13*t31*t75 + d23*t83 + d25*t88 + d26*t88*T(3); 
complex<T> t8 = d42*s23*t28 + d43*spa23*spb13*t31*t75 + d31*s23*t76 + d39*spa23*t83 + d44*s23*spb35*(t64*t66 + spb35*(-t18 + t87)) + s23*spb35*t64*t90 + d40*spa23*t88*T(3); 
complex<T> t9 = d42*s34*t28 + d45*spa34*t49*t54 + d31*s34*t76 + d44*s34*spb35*(t64*t66 + spb35*(-t18 + t87)) + s34*spb35*t64*t90; 
complex<T> t12 = d16*spa24*t28 - d27*spa14*spb13*t28 - d13*t16*t28 + d23*t34*t51 + d20*spb13*t31*t75 + d5*t33*t76 + d6*t33*t76 + d24*t83 + d3*t83 + d4*t83 - d25*t88 + d10*spb23*t93 + d11*spb13*t52*t93 - d21*(t64*t66 + spb35*(-t18 + t87))*t93 + d28*t31*t97 + d15*t50*t97 - d26*t88*T(3) + d9*t88*T(3); 
complex<T> t14 = d18*t111 - d16*spa24*t28 + d17*spa12*spb13*t28 + d13*t16*t28 + d14*t16*t28 + d3*t34*t51 + d4*t34*t51 + d15*spb13*t50*t75 + d5*t33*t89 + d6*t33*t89 + d7*t33*t89 + d8*t33*t89 - d10*spb23*t93 + d21*(t64*t66 + spb35*(-t18 + t87))*t93 + d22*(t64*t66 + spb35*(-t18 + t87))*t93 + d20*t31*t97 - d19*spb13*t49*T(3) - d9*t88*T(3) + d11*t64*t93*T(3) + d12*t64*t93*T(3); 
complex<T> t61 = d2*t111 - d1*spa25*t49; 
complex<T> t74 = d46*t47*t83 + d47*s15*spa45*t31*t97 + d47*t47*t88*T(3); 
complex<T> t86 = d32*T(3); 
complex<T> t11 = d34*spa15*t28 - d36*spa15*spb13*t49 + d35*s15*t49*t54 + d37*s15*spb13*t31*t75 + d30*s15*t83 + d38*spa15*spb35*(t64*t66 + spb35*(-t18 + t87)) + s15*t86*t88 + d31*s15*t89 + d33*spa15*spb35*t64*T(3); 
complex<T> t62 = d30*s45*t83 + s45*t86*t88 + d32*spa45*t31*t97; 
complex<T> co1 = d29*spa12*t111; 
complex<T> co2 = d48*s15*spa12*t111; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t61*Int(ep,mu,c12,c345) + t14*Int(ep,mu,c15,c234) + t12*Int(ep,mu,c23,c145) + t13*Int(ep,mu,c34,c125) + t7*Int(ep,mu,c45,c123) + co1*Int(ep,mu,c1,c2,c345) + t11*Int(ep,mu,c1,c5,c234) + t8*Int(ep,mu,c2,c3,c145) + t9*Int(ep,mu,c3,c4,c125) + t62*Int(ep,mu,c4,c5,c123) + t74*Int(ep,mu,c1,c5,c4,c23) + co2*Int(ep,mu,c2,c1,c5,c34) + t10*Int(ep,mu,c4,c3,c2,c15));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmmpqp_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, m, p, qp}, L}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmmpqp L");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:38:04 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = S(4,5);
complex<T> t3 = square(s12 + s15 - s34); 
complex<T> t4 = spb12*spb23; 
complex<T> t13 = square(spb45); 
complex<T> t20 = square(spa23); 
complex<T> t25 = s25*spa25; 
complex<T> t26 = spb25*spb34; 
complex<T> t31 = -(spa12*spb15); 
complex<T> t34 = spb14*spb45; 
complex<T> t43 = spa15*spb14; 
complex<T> d6 = (s15 - s34)*spb15*spb23; d6 = T(1)/d6;
complex<T> d9 = spb15*spb23; d9 = T(1)/d9;
complex<T> d10 = s15 - s34; d10 = T(1)/d10;
complex<T> d12 = spb23*spb35; d12 = T(1)/d12;
complex<T> t7 = -t34; 
complex<T> t8 = t4*T(2); 
complex<T> t14 = -t25; 
complex<T> t18 = spb24*t31 - spa23*t26*T(3); 
complex<T> t21 = -(spa23*t26) + spb24*t31; 
complex<T> t22 = -t43; 
complex<T> t24 = spb24*t13; 
complex<T> t45 = d9*spa23; 
complex<T> d1 = (s12 - s34)*t4; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spb34*t4; d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*(s12 + s15 - s34)*spb34*t4; d3 = T(1)/d3;
complex<T> d8 = (s15 - s34)*(s12 + s15 - s34)*spb34*t4; d8 = T(1)/d8;
complex<T> d11 = spb23*spb34*t3; d11 = T(1)/d11;
complex<T> d13 = spb34*t4; d13 = T(1)/d13;
complex<T> d14 = spb34*t3*t4; d14 = T(1)/d14;
complex<T> d15 = t3*t4; d15 = T(1)/d15;
complex<T> d16 = spb35*t4; d16 = T(1)/d16;
complex<T> d17 = spb23*spb34*t3*T(2); d17 = T(1)/d17;
complex<T> t29 = d1*spa35; 
complex<T> t37 = spa12*t14; 
complex<T> t41 = t24*t25; 
complex<T> t46 = -(d16*t13); 
complex<T> d4 = t8*square(s15 - s34); d4 = T(1)/d4;
complex<T> d5 = (s15 - s34)*spb15*spb34*t8; d5 = T(1)/d5;
complex<T> d7 = spb15*spb34*t8; d7 = T(1)/d7;
complex<T> d18 = spb35*t8; d18 = T(1)/d18;
complex<T> t27 = d7*t18; 
complex<T> t36 = d4*spb24; 
complex<T> t40 = d13*t13*t22 + d14*s15*t41; 
complex<T> t44 = d15*spa34*t41 + s34*t46; 
complex<T> t48 = d3*t14*t24 + t13*(t29 + d2*t43); 
complex<T> t56 = t24*t37; 
complex<T> t1 = d6*spa23*t13 + d8*t14*t24 + d10*t27*t34 + t20*t26*t36 - d10*t13*t45 - d7*spb14*t13*T(3) - d5*t21*t34*T(3); 
complex<T> t2 = -(d6*spa23*t13) + d2*t13*t22 - t13*t29 - t20*t26*t36 + d3*t41 + d8*t41 + d10*t13*t45 + d10*t27*t7 + d5*t21*t34*T(3); 
complex<T> t33 = -(d12*spa12*t13) + d11*t56; 
complex<T> co1 = s45*t46; 
complex<T> co2 = d17*s15*t56; 
complex<T> co3 = -(d18*s34*s45*t13); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t48*Int(ep,mu,c12,c345) + t1*Int(ep,mu,c15,c234) + t2*Int(ep,mu,c34,c125) + t33*Int(ep,mu,c1,c2,c345) + t40*Int(ep,mu,c1,c5,c234) + t44*Int(ep,mu,c3,c4,c125) + co1*Int(ep,mu,c4,c5,c123) + co2*Int(ep,mu,c2,c1,c5,c34) + co3*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqpppp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpppp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmqpmpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpmpp nf");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:38:06 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t4 = square(spb45); 
complex<T> t6 = square(spa13); 
complex<T> t9 = spa35*spb25; 
complex<T> t10 = spa15*spa34; 
complex<T> d1 = (s12 - s34)*spa12*spa45; d1 = T(1)/d1;
complex<T> d2 = spa12*spa45*square(s12 - s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = spa45*cube(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> t11 = d3*T(2); 
complex<T> t12 = d2*spa13; 
complex<T> t14 = -(d1*spb45); 
complex<T> t16 = t10*t4; 
complex<T> t3 = -(t12*t16) + d1*spb45*t6 - d3*t16*t9*T(2); 
complex<T> t15 = t12*t16 + t14*t6 + t11*t16*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*Int(ep,mu,c12,c345) + t15*Int(ep,mu,c34,c125));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqppmp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqppmp nf");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:38:17 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> s34 = S(3,4);
complex<T> s45 = S(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s35 = -(spa35*spb35);
complex<T> t3 = square(s35); 
complex<T> t6 = cube(s35); 
complex<T> t10 = square(spb35); 
complex<T> t13 = square(spa14); 
complex<T> t17 = s34*s45; 
complex<T> t18 = cube(spb35); 
complex<T> t19 = spa13*spa45; 
complex<T> t20 = spa15*spa34; 
complex<T> t21 = spa12*T(3); 
complex<T> t25 = square(square(spb35)); 
complex<T> t30 = -(spa14*T(2)); 
complex<T> t40 = spa14*T(2); 
complex<T> d1 = (s12 - s34)*s35*spa12; d1 = T(1)/d1;
complex<T> d2 = s35*(s12 - s45)*spa12; d2 = T(1)/d2;
complex<T> d11 = square(spa35); d11 = T(1)/d11;
complex<T> d12 = square(square(spa35)); d12 = T(1)/d12;
complex<T> d13 = spa12*square(spa35); d13 = T(1)/d13;
complex<T> d14 = spa12*square(square(spa35)); d14 = T(1)/d14;
complex<T> d15 = spa12*square(spa35)*T(2); d15 = T(1)/d15;
complex<T> t11 = -t20; 
complex<T> t22 = -(d14*T(2)); 
complex<T> t26 = t19*t20; 
complex<T> d3 = t21*cube(s12 - s34); d3 = T(1)/d3;
complex<T> d4 = t21*cube(s12 - s45); d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*spa12*t6; d5 = T(1)/d5;
complex<T> d6 = spa12*t3*square(s12 - s34); d6 = T(1)/d6;
complex<T> d7 = s35*t21*cube(s12 - s34); d7 = T(1)/d7;
complex<T> d8 = s35*t21*cube(s12 - s45); d8 = T(1)/d8;
complex<T> d9 = spa12*t3*square(s12 - s45); d9 = T(1)/d9;
complex<T> d10 = (s12 - s45)*spa12*t6; d10 = T(1)/d10;
complex<T> t7 = s34*(d13*t13 + t22*t26); 
complex<T> t8 = spb12*(d11*t13 - d12*t26*T(2)); 
complex<T> t16 = -(d10*T(2)); 
complex<T> t29 = t25*t26; 
complex<T> t32 = s45*(d13*t13 + t22*t26); 
complex<T> t41 = t11*t19; 
complex<T> t9 = -(d1*t10*t13) + d6*t29 + d3*t18*t20*t40 + d5*t29*T(2) + d7*t29*T(2); 
complex<T> t24 = -(d2*t10*t13) + d9*t29 + d4*t18*t19*t40 + d10*t29*T(2) + d8*t29*T(2); 
complex<T> t28 = t17*(d15*t13 + d14*t41); 
complex<T> t33 = d1*t10*t13 + d2*t10*t13 + t16*t29 + d4*t18*t19*t30 + d3*t18*t20*t30 + d6*t25*t41 + d9*t25*t41 - d5*t29*T(2) - d7*t29*T(2) - d8*t29*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t33*Int(ep,mu,c12,c345) + t9*Int(ep,mu,c34,c125) + t24*Int(ep,mu,c45,c123) + t8*Int(ep,mu,c1,c2,c345) + t7*Int(ep,mu,c3,c4,c125) + t32*Int(ep,mu,c4,c5,c123) + t28*Int(ep,mu,c5,c4,c3,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqpppm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpppm nf");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:38:18 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa13 = SPA(1,3);
complex<T> spa45 = SPA(4,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> t5 = square(spb34); 
complex<T> t7 = square(spa13); 
complex<T> t8 = square(spa15); 
complex<T> t9 = square(spa45); 
complex<T> t11 = cube(spb34); 
complex<T> t16 = spa15*spa45; 
complex<T> d1 = (s12 - s45)*spa12*spa34; d1 = T(1)/d1;
complex<T> d2 = spa12*spa34*square(s12 - s45); d2 = T(1)/d2;
complex<T> d3 = spa12*spa34*cube(s12 - s45)*T(3); d3 = T(1)/d3;
complex<T> t12 = -(d2*spa13); 
complex<T> t14 = d3*T(2); 
complex<T> t20 = t7*t9; 
complex<T> t3 = d2*spa13*t16*t5 - d1*spb34*t8 - d3*t11*t20*T(2); 
complex<T> t15 = t11*t14*t20 + t12*t16*t5 + d1*spb34*t8; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t15*Int(ep,mu,c12,c345) + t3*Int(ep,mu,c45,c123));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmpqppp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpqppp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmqppp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmqppp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpqpmp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpqpmp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpqppm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpqppm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmppqpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, p, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmppqpp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmpqpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, p, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmpqpp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpmqpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, m, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpmqpp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmppqpm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, p, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmppqpm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpppqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, p, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpppqp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmppqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, p, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmppqp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpmpqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, m, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpmpqp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmppmqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, p, m, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmppmqp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmqpmmm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpmmm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmqppmm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqppmm nf");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:38:19 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t5 = square(spa45); 
complex<T> t7 = square(spb23); 
complex<T> t8 = square(spb25); 
complex<T> t9 = square(spb34); 
complex<T> t11 = cube(spa45); 
complex<T> t16 = spb25*spb34; 
complex<T> d1 = (s12 - s34)*spb12*spb45; d1 = T(1)/d1;
complex<T> d2 = spb12*spb45*square(s12 - s34); d2 = T(1)/d2;
complex<T> d3 = spb12*spb45*cube(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> t12 = -(d2*spb23); 
complex<T> t14 = d3*T(2); 
complex<T> t20 = t8*t9; 
complex<T> t3 = d2*spb23*t16*t5 - d1*spa45*t7 - d3*t11*t20*T(2); 
complex<T> t15 = t11*t14*t20 + t12*t16*t5 + d1*spa45*t7; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t15*Int(ep,mu,c12,c345) + t3*Int(ep,mu,c34,c125));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqpmpm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpmpm nf");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:38:31 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = S(3,4);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = S(4,5);
complex<T> t3 = square(s35); 
complex<T> t6 = cube(s35); 
complex<T> t10 = square(spa35); 
complex<T> t13 = square(spb24); 
complex<T> t17 = s34*s45; 
complex<T> t18 = cube(spa35); 
complex<T> t19 = spb23*spb45; 
complex<T> t20 = spb25*spb34; 
complex<T> t21 = spb12*T(3); 
complex<T> t25 = square(square(spa35)); 
complex<T> t30 = -(spb24*T(2)); 
complex<T> t40 = spb24*T(2); 
complex<T> d1 = (s12 - s34)*s35*spb12; d1 = T(1)/d1;
complex<T> d2 = s35*(s12 - s45)*spb12; d2 = T(1)/d2;
complex<T> d11 = square(spb35); d11 = T(1)/d11;
complex<T> d12 = square(square(spb35)); d12 = T(1)/d12;
complex<T> d13 = spb12*square(spb35); d13 = T(1)/d13;
complex<T> d14 = spb12*square(square(spb35)); d14 = T(1)/d14;
complex<T> d15 = spb12*square(spb35)*T(2); d15 = T(1)/d15;
complex<T> t11 = -t20; 
complex<T> t22 = -(d14*T(2)); 
complex<T> t26 = t19*t20; 
complex<T> d3 = t21*cube(s12 - s34); d3 = T(1)/d3;
complex<T> d4 = t21*cube(s12 - s45); d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*spb12*t6; d5 = T(1)/d5;
complex<T> d6 = spb12*t3*square(s12 - s34); d6 = T(1)/d6;
complex<T> d7 = s35*t21*cube(s12 - s34); d7 = T(1)/d7;
complex<T> d8 = s35*t21*cube(s12 - s45); d8 = T(1)/d8;
complex<T> d9 = spb12*t3*square(s12 - s45); d9 = T(1)/d9;
complex<T> d10 = (s12 - s45)*spb12*t6; d10 = T(1)/d10;
complex<T> t7 = s34*(d13*t13 + t22*t26); 
complex<T> t8 = spa12*(d11*t13 - d12*t26*T(2)); 
complex<T> t16 = -(d10*T(2)); 
complex<T> t29 = t25*t26; 
complex<T> t32 = s45*(d13*t13 + t22*t26); 
complex<T> t41 = t11*t19; 
complex<T> t9 = -(d1*t10*t13) + d6*t29 + d3*t18*t20*t40 + d5*t29*T(2) + d7*t29*T(2); 
complex<T> t24 = -(d2*t10*t13) + d9*t29 + d4*t18*t19*t40 + d10*t29*T(2) + d8*t29*T(2); 
complex<T> t28 = t17*(d15*t13 + d14*t41); 
complex<T> t33 = d1*t10*t13 + d2*t10*t13 + t16*t29 + d4*t18*t19*t30 + d3*t18*t20*t30 + d6*t25*t41 + d9*t25*t41 - d5*t29*T(2) - d7*t29*T(2) - d8*t29*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t33*Int(ep,mu,c12,c345) + t9*Int(ep,mu,c34,c125) + t24*Int(ep,mu,c45,c123) + t8*Int(ep,mu,c1,c2,c345) + t7*Int(ep,mu,c3,c4,c125) + t32*Int(ep,mu,c4,c5,c123) + t28*Int(ep,mu,c3,c4,c5,c12));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmqpmmp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, qp, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmqpmmp nf");
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
 // #define TimeStamp "Thu 11 Dec 2008 21:38:32 on n301"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> t4 = square(spa34); 
complex<T> t6 = square(spb25); 
complex<T> t9 = spa13*spb35; 
complex<T> t10 = spb23*spb45; 
complex<T> d1 = (s12 - s45)*spb12*spb34; d1 = T(1)/d1;
complex<T> d2 = spb12*spb34*square(s12 - s45)*T(3); d2 = T(1)/d2;
complex<T> d3 = spb34*cube(s12 - s45)*T(3); d3 = T(1)/d3;
complex<T> t11 = d3*T(2); 
complex<T> t12 = -(d1*spa34); 
complex<T> t13 = d2*spb25; 
complex<T> t16 = t10*t4; 
complex<T> t3 = -(t13*t16) + d1*spa34*t6 - d3*t16*t9*T(2); 
complex<T> t15 = t13*t16 + t12*t6 + t11*t16*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*Int(ep,mu,c12,c345) + t15*Int(ep,mu,c45,c123));  
 return(result);
} 
  
  
 
 
template <class T> SeriesC<T> C2q3g_qmmqpmm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmqpmm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpqpmm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, qp, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpqpmm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmqppm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmqppm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmqpmp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, qp, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmqpmp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmmqpm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, m, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmmqpm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpmqpm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, m, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpmqpm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmpqpm_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, p, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmpqpm nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmmqpp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, m, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmmqpp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmmmqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, m, m, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmmmqp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmpmmqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, p, m, m, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmpmmqp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmpmqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, p, m, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmpmqp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
template <class T> SeriesC<T> C2q3g_qmmmpqp_nf
      (const eval_param<T>& ep,
                 const T& mu){
//{{qm, m, m, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q3g :  qmmmpqp nf");
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
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
 // *************** table of switch values ************* 
 
#define _C_qmqpppp_L C2q3g_1017_L
#define _C_qmqpmpp_L C2q3g_969_L
#define _C_qmqppmp_L C2q3g_825_L
#define _C_qmqpppm_L C2q3g_249_L
#define _C_qmpqppp_L C2q3g_1005_L
#define _C_qmmqppp_L C2q3g_993_L
#define _C_qmpqpmp_L C2q3g_813_L
#define _C_qmpqppm_L C2q3g_237_L
#define _C_qmppqpp_L C2q3g_957_L
#define _C_qmmpqpp_L C2q3g_945_L
#define _C_qmpmqpp_L C2q3g_909_L
#define _C_qmppqpm_L C2q3g_189_L
#define _C_qmpppqp_L C2q3g_765_L
#define _C_qmmppqp_L C2q3g_753_L
#define _C_qmpmpqp_L C2q3g_717_L
#define _C_qmppmqp_L C2q3g_573_L
#define _C_qmqpmmm_L C2q3g_9_L
#define _C_qmqppmm_L C2q3g_57_L
#define _C_qmqpmpm_L C2q3g_201_L
#define _C_qmqpmmp_L C2q3g_777_L
#define _C_qmmqpmm_L C2q3g_33_L
#define _C_qmpqpmm_L C2q3g_45_L
#define _C_qmmqppm_L C2q3g_225_L
#define _C_qmmqpmp_L C2q3g_801_L
#define _C_qmmmqpm_L C2q3g_129_L
#define _C_qmpmqpm_L C2q3g_141_L
#define _C_qmmpqpm_L C2q3g_177_L
#define _C_qmmmqpp_L C2q3g_897_L
#define _C_qmmmmqp_L C2q3g_513_L
#define _C_qmpmmqp_L C2q3g_525_L
#define _C_qmmpmqp_L C2q3g_561_L
#define _C_qmmmpqp_L C2q3g_705_L
#define _C_qmqpppp_nf C2q3g_1017_nf
#define _C_qmqpmpp_nf C2q3g_969_nf
#define _C_qmqppmp_nf C2q3g_825_nf
#define _C_qmqpppm_nf C2q3g_249_nf
#define _C_qmpqppp_nf C2q3g_1005_nf
#define _C_qmmqppp_nf C2q3g_993_nf
#define _C_qmpqpmp_nf C2q3g_813_nf
#define _C_qmpqppm_nf C2q3g_237_nf
#define _C_qmppqpp_nf C2q3g_957_nf
#define _C_qmmpqpp_nf C2q3g_945_nf
#define _C_qmpmqpp_nf C2q3g_909_nf
#define _C_qmppqpm_nf C2q3g_189_nf
#define _C_qmpppqp_nf C2q3g_765_nf
#define _C_qmmppqp_nf C2q3g_753_nf
#define _C_qmpmpqp_nf C2q3g_717_nf
#define _C_qmppmqp_nf C2q3g_573_nf
#define _C_qmqpmmm_nf C2q3g_9_nf
#define _C_qmqppmm_nf C2q3g_57_nf
#define _C_qmqpmpm_nf C2q3g_201_nf
#define _C_qmqpmmp_nf C2q3g_777_nf
#define _C_qmmqpmm_nf C2q3g_33_nf
#define _C_qmpqpmm_nf C2q3g_45_nf
#define _C_qmmqppm_nf C2q3g_225_nf
#define _C_qmmqpmp_nf C2q3g_801_nf
#define _C_qmmmqpm_nf C2q3g_129_nf
#define _C_qmpmqpm_nf C2q3g_141_nf
#define _C_qmmpqpm_nf C2q3g_177_nf
#define _C_qmmmqpp_nf C2q3g_897_nf
#define _C_qmmmmqp_nf C2q3g_513_nf
#define _C_qmpmmqp_nf C2q3g_525_nf
#define _C_qmmpmqp_nf C2q3g_561_nf
#define _C_qmmmpqp_nf C2q3g_705_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmqpppp_L case 1017 : \
          return &C2q3g_1017_L
#define _CASE_qmqpmpp_L case 969 : \
          return &C2q3g_969_L
#define _CASE_qmqppmp_L case 825 : \
          return &C2q3g_825_L
#define _CASE_qmqpppm_L case 249 : \
          return &C2q3g_249_L
#define _CASE_qmpqppp_L case 1005 : \
          return &C2q3g_1005_L
#define _CASE_qmmqppp_L case 993 : \
          return &C2q3g_993_L
#define _CASE_qmpqpmp_L case 813 : \
          return &C2q3g_813_L
#define _CASE_qmpqppm_L case 237 : \
          return &C2q3g_237_L
#define _CASE_qmppqpp_L case 957 : \
          return &C2q3g_957_L
#define _CASE_qmmpqpp_L case 945 : \
          return &C2q3g_945_L
#define _CASE_qmpmqpp_L case 909 : \
          return &C2q3g_909_L
#define _CASE_qmppqpm_L case 189 : \
          return &C2q3g_189_L
#define _CASE_qmpppqp_L case 765 : \
          return &C2q3g_765_L
#define _CASE_qmmppqp_L case 753 : \
          return &C2q3g_753_L
#define _CASE_qmpmpqp_L case 717 : \
          return &C2q3g_717_L
#define _CASE_qmppmqp_L case 573 : \
          return &C2q3g_573_L
#define _CASE_qmqpmmm_L case 9 : \
          return &C2q3g_9_L
#define _CASE_qmqppmm_L case 57 : \
          return &C2q3g_57_L
#define _CASE_qmqpmpm_L case 201 : \
          return &C2q3g_201_L
#define _CASE_qmqpmmp_L case 777 : \
          return &C2q3g_777_L
#define _CASE_qmmqpmm_L case 33 : \
          return &C2q3g_33_L
#define _CASE_qmpqpmm_L case 45 : \
          return &C2q3g_45_L
#define _CASE_qmmqppm_L case 225 : \
          return &C2q3g_225_L
#define _CASE_qmmqpmp_L case 801 : \
          return &C2q3g_801_L
#define _CASE_qmmmqpm_L case 129 : \
          return &C2q3g_129_L
#define _CASE_qmpmqpm_L case 141 : \
          return &C2q3g_141_L
#define _CASE_qmmpqpm_L case 177 : \
          return &C2q3g_177_L
#define _CASE_qmmmqpp_L case 897 : \
          return &C2q3g_897_L
#define _CASE_qmmmmqp_L case 513 : \
          return &C2q3g_513_L
#define _CASE_qmpmmqp_L case 525 : \
          return &C2q3g_525_L
#define _CASE_qmmpmqp_L case 561 : \
          return &C2q3g_561_L
#define _CASE_qmmmpqp_L case 705 : \
          return &C2q3g_705_L
#define _CASE_qmqpppp_nf case 1017 : \
          return &C2q3g_1017_nf
#define _CASE_qmqpmpp_nf case 969 : \
          return &C2q3g_969_nf
#define _CASE_qmqppmp_nf case 825 : \
          return &C2q3g_825_nf
#define _CASE_qmqpppm_nf case 249 : \
          return &C2q3g_249_nf
#define _CASE_qmpqppp_nf case 1005 : \
          return &C2q3g_1005_nf
#define _CASE_qmmqppp_nf case 993 : \
          return &C2q3g_993_nf
#define _CASE_qmpqpmp_nf case 813 : \
          return &C2q3g_813_nf
#define _CASE_qmpqppm_nf case 237 : \
          return &C2q3g_237_nf
#define _CASE_qmppqpp_nf case 957 : \
          return &C2q3g_957_nf
#define _CASE_qmmpqpp_nf case 945 : \
          return &C2q3g_945_nf
#define _CASE_qmpmqpp_nf case 909 : \
          return &C2q3g_909_nf
#define _CASE_qmppqpm_nf case 189 : \
          return &C2q3g_189_nf
#define _CASE_qmpppqp_nf case 765 : \
          return &C2q3g_765_nf
#define _CASE_qmmppqp_nf case 753 : \
          return &C2q3g_753_nf
#define _CASE_qmpmpqp_nf case 717 : \
          return &C2q3g_717_nf
#define _CASE_qmppmqp_nf case 573 : \
          return &C2q3g_573_nf
#define _CASE_qmqpmmm_nf case 9 : \
          return &C2q3g_9_nf
#define _CASE_qmqppmm_nf case 57 : \
          return &C2q3g_57_nf
#define _CASE_qmqpmpm_nf case 201 : \
          return &C2q3g_201_nf
#define _CASE_qmqpmmp_nf case 777 : \
          return &C2q3g_777_nf
#define _CASE_qmmqpmm_nf case 33 : \
          return &C2q3g_33_nf
#define _CASE_qmpqpmm_nf case 45 : \
          return &C2q3g_45_nf
#define _CASE_qmmqppm_nf case 225 : \
          return &C2q3g_225_nf
#define _CASE_qmmqpmp_nf case 801 : \
          return &C2q3g_801_nf
#define _CASE_qmmmqpm_nf case 129 : \
          return &C2q3g_129_nf
#define _CASE_qmpmqpm_nf case 141 : \
          return &C2q3g_141_nf
#define _CASE_qmmpqpm_nf case 177 : \
          return &C2q3g_177_nf
#define _CASE_qmmmqpp_nf case 897 : \
          return &C2q3g_897_nf
#define _CASE_qmmmmqp_nf case 513 : \
          return &C2q3g_513_nf
#define _CASE_qmpmmqp_nf case 525 : \
          return &C2q3g_525_nf
#define _CASE_qmmpmqp_nf case 561 : \
          return &C2q3g_561_nf
#define _CASE_qmmmpqp_nf case 705 : \
          return &C2q3g_705_nf
 
 
 // *************** function definitions using macros ************* 
 
template <class T> SeriesC<T> _C_qmqpppp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpppp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmpp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpmpp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppmp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqppmp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpppm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpppm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqppp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpqppp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqppp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmqppp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpmp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpqpmp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqppm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpqppm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmppqpp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmppqpp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmpqpp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmpqpp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpmqpp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpmqpp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmppqpm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmppqpm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpppqp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpppqp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmppqp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmppqp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpmpqp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpmpqp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmppmqp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmppmqp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmmm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpmmm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppmm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqppmm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmpm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpmpm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmmp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpmmp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqpmm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmqpmm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpmm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpqpmm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqppm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmqppm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqpmp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmqpmp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmmqpm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmmqpm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpmqpm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpmqpm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmpqpm_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmpqpm_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmmqpp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmmqpp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmmmqp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmmmqp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpmmqp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpmmqp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmpmqp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmpmqp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmmpqp_L(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmmpqp_L(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpppp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpppp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpmpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppmp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqppmp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpppm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpppm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqppp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpqppp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqppp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmqppp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpmp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpqpmp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqppm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpqppm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmppqpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmppqpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmpqpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmpqpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpmqpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpmqpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmppqpm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmppqpm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpppqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpppqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmppqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmppqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpmpqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpmpqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmppmqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmppmqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmmm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpmmm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqppmm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqppmm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmpm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpmpm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmqpmmp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmqpmmp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqpmm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmqpmm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpqpmm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpqpmm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqppm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmqppm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmqpmp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmqpmp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmmqpm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmmqpm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpmqpm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpmqpm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmpqpm_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmpqpm_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmmqpp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmmqpp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmmmqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmmmqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmpmmqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmpmmqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmpmqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmpmqp_nf(ep,mu);}
 
template <class T> SeriesC<T> _C_qmmmpqp_nf(
        const eval_param<T>& ep, const T& mu){
          return C2q3g_qmmmpqp_nf(ep,mu);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> SeriesC<T> ( *C2q3g_L_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqpppp_L;
       _CASE_qmqpmpp_L;
       _CASE_qmqppmp_L;
       _CASE_qmqpppm_L;
       _CASE_qmpqppp_L;
       _CASE_qmmqppp_L;
       _CASE_qmpqpmp_L;
       _CASE_qmpqppm_L;
       _CASE_qmppqpp_L;
       _CASE_qmmpqpp_L;
       _CASE_qmpmqpp_L;
       _CASE_qmppqpm_L;
       _CASE_qmpppqp_L;
       _CASE_qmmppqp_L;
       _CASE_qmpmpqp_L;
       _CASE_qmppmqp_L;
       _CASE_qmqpmmm_L;
       _CASE_qmqppmm_L;
       _CASE_qmqpmpm_L;
       _CASE_qmqpmmp_L;
       _CASE_qmmqpmm_L;
       _CASE_qmpqpmm_L;
       _CASE_qmmqppm_L;
       _CASE_qmmqpmp_L;
       _CASE_qmmmqpm_L;
       _CASE_qmpmqpm_L;
       _CASE_qmmpqpm_L;
       _CASE_qmmmqpp_L;
       _CASE_qmmmmqp_L;
       _CASE_qmpmmqp_L;
       _CASE_qmmpmqp_L;
       _CASE_qmmmpqp_L;
 
       default: return 0;
        }
 }
 
template <class T> SeriesC<T> ( *C2q3g_nf_Ptr_eval( int hc))
     (const eval_param<T>&, const T&) {
       switch (hc) {
       _CASE_qmqpppp_nf;
       _CASE_qmqpmpp_nf;
       _CASE_qmqppmp_nf;
       _CASE_qmqpppm_nf;
       _CASE_qmpqppp_nf;
       _CASE_qmmqppp_nf;
       _CASE_qmpqpmp_nf;
       _CASE_qmpqppm_nf;
       _CASE_qmppqpp_nf;
       _CASE_qmmpqpp_nf;
       _CASE_qmpmqpp_nf;
       _CASE_qmppqpm_nf;
       _CASE_qmpppqp_nf;
       _CASE_qmmppqp_nf;
       _CASE_qmpmpqp_nf;
       _CASE_qmppmqp_nf;
       _CASE_qmqpmmm_nf;
       _CASE_qmqppmm_nf;
       _CASE_qmqpmpm_nf;
       _CASE_qmqpmmp_nf;
       _CASE_qmmqpmm_nf;
       _CASE_qmpqpmm_nf;
       _CASE_qmmqppm_nf;
       _CASE_qmmqpmp_nf;
       _CASE_qmmmqpm_nf;
       _CASE_qmpmqpm_nf;
       _CASE_qmmpqpm_nf;
       _CASE_qmmmqpp_nf;
       _CASE_qmmmmqp_nf;
       _CASE_qmpmmqp_nf;
       _CASE_qmmpmqp_nf;
       _CASE_qmmmpqp_nf;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template SeriesC<R> ( *C2q3g_L_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q3g_L_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q3g_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);


template SeriesC<R> ( *C2q3g_nf_Ptr_eval(int hc))
             (const eval_param<R>&, const R&);
template SeriesC<RHP> ( *C2q3g_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&, const RHP&);
template SeriesC<RVHP> ( *C2q3g_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&, const RVHP&);




}
