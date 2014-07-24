/*
*C_2q2g1y_wCI.cpp
*
* Created on 2/5, 2011
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#include <vector>
#include "integrals.h"
#include "cached_integral.h"
using namespace std;

using BH::CachedIntegral::Cached_Bubble_Integral_User;
using BH::CachedIntegral::Cached_Triangle_Integral_User;
using BH::CachedIntegral::Cached_Box_Integral_User;
 
namespace BH  {
 
class Index_Vector;
namespace CachedIntegral {
 
#define _VERBOSE 0
 
#define SPA(i,j) mc.spa(ind[i-1],ind[j-1])
#define SPB(i,j) mc.spb(ind[i-1],ind[j-1])
#define S(i,j) mc.s(ind[i-1],ind[j-1])
#define SS(i,j,k) mc.s(ind[i-1],ind[j-1],ind[k-1])

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}
 



GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmgamqpmm_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmgamqpmm_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmgamqpmp_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmgamqpmp_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmgamqppm_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmgamqppm_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmgamqppp_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmgamqppp_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmmmqpgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmmmqpgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmmpqpgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmmpqpgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmmqpmgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmmqppgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmpmqpgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmpmqpgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmppqpgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmppqpgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmpqpmgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmpqppgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmqpmmgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmqpmpgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmqppmgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qmqpppgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpgapqmmm_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpgapqmmm_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpgapqmmp_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpgapqmmp_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpgapqmpm_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpgapqmpm_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpgapqmpp_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpgapqmpp_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpmmqmgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpmmqmgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpmpqmgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpmpqmgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpmqmmgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpmqmpgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qppmqmgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qppmqmgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpppqmgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpppqmgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qppqmmgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qppqmpgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpqmmmgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpqmmpgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpqmpmgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g1y_qpqmppgap_SLC_wCI)

 

C2q2g1y_qmmmqpgam_L_wCI::\
C2q2g1y_qmmmqpgam_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmmmqpgam_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmmmqpgam L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1y_qmmpqpgam_L_wCI::\
C2q2g1y_qmmpqpgam_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmmpqpgam_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, p, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmmpqpgam L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:06:23 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = -(spa15*spb15);
complex<T> s12 = -(spa12*spb12);
complex<T> t1 = spb12*spb15; 
complex<T> t6 = square(spb34); 
complex<T> t10 = square(spb14); 
complex<T> t13 = spa12*spb13; 
complex<T> t21 = square(spa12); 
complex<T> t25 = -(spa12*spa25); 
complex<T> t27 = spb14*spb34; 
complex<T> d4 = (s12 + s15 - s34)*spb45; d4 = T(1)/d4;
complex<T> d5 = (s12 + s15 - s34)*spb12*spb45; d5 = T(1)/d5;
complex<T> d6 = spb12*spb24*spb45; d6 = T(1)/d6;
complex<T> d9 = spb15*spb45*T(2); d9 = T(1)/d9;
complex<T> d13 = (s12 + s15 - s34)*spb45*T(2); d13 = T(1)/d13;
complex<T> d14 = spb23*spb45*T(2); d14 = T(1)/d14;
complex<T> t8 = -t13; 
complex<T> t14 = s34*t6; 
complex<T> t16 = d5*spa25; 
complex<T> t22 = spb13*t10; 
complex<T> t32 = s15*t6; 
complex<T> d1 = spb45*t1*square(s23 - s45)*T(2); d1 = T(1)/d1;
complex<T> d2 = (s23 - s45)*spb45*t1; d2 = T(1)/d2;
complex<T> d3 = spb23*spb45*t1*T(2); d3 = T(1)/d3;
complex<T> d7 = spb45*t1; d7 = T(1)/d7;
complex<T> d8 = spb24*spb45*t1; d8 = T(1)/d8;
complex<T> d10 = spb45*t1*T(2); d10 = T(1)/d10;
complex<T> d11 = spb24*spb45*t1*T(2); d11 = T(1)/d11;
complex<T> d12 = spb23*t1*T(2); d12 = T(1)/d12;
complex<T> t4 = d1*spb23*t21*t22 + d2*t13*t27*T(2); 
complex<T> t11 = -t16; 
complex<T> t12 = -(d1*spb23); 
complex<T> t15 = -(d2*T(2)); 
complex<T> t17 = d8*spb14; 
complex<T> t23 = (d10*spa23*spb13 + d11*s23*spb14)*t14; 
complex<T> t33 = t16*t32 + d6*spa15*spb14*t6; 
complex<T> t35 = d13*t25*t32 + d14*spa15*t6*t8; 
complex<T> t5 = (d7*spa23*spb13 + s23*t17)*t6; 
complex<T> t20 = t12*t21*t22 + t13*t15*t27 - d3*spb13*t6*T(3); 
complex<T> t29 = t14*(t11 + t17); 
complex<T> co1 = d4*t25*t6; 
complex<T> co2 = d9*spa23*t6*t8; 
complex<T> co3 = d12*spa45*spb13*t14; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t20*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + t33*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + t29*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t23*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t35*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qmpmqpgam_L_wCI::\
C2q2g1y_qmpmqpgam_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmpmqpgam_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, m, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmpmqpgam L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:07:00 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> s13 = -(spa13*spb13);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = S(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t3 = spb15*spb45; 
complex<T> t4 = square(s12 - s45); 
complex<T> t19 = square(spb24); 
complex<T> t20 = square(spa13); 
complex<T> t21 = square(spa35); 
complex<T> t22 = square(spb14); 
complex<T> t25 = -(spb12*T(2)); 
complex<T> t29 = -(s12*s23); 
complex<T> t32 = cube(spa13); 
complex<T> t33 = cube(spa35); 
complex<T> t34 = cube(spb24); 
complex<T> t36 = -(spb25*T(2)); 
complex<T> t38 = s12*spb14; 
complex<T> t47 = spb24*T(2); 
complex<T> t48 = spb23*spb45; 
complex<T> t54 = spb12*spb23; 
complex<T> t56 = spa35*spb24; 
complex<T> d1 = spb15*square(s12 - s34)*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*s35*spb15; d2 = T(1)/d2;
complex<T> d3 = s35*(s12 - s45)*spb15; d3 = T(1)/d3;
complex<T> d9 = (s12 - s34)*spb15*square(s35); d9 = T(1)/d9;
complex<T> d10 = s35*spb15*square(s12 - s34)*T(2); d10 = T(1)/d10;
complex<T> d12 = (s12 - s45)*spb15*square(s35); d12 = T(1)/d12;
complex<T> d15 = spb15*spb34*spb35*T(2); d15 = T(1)/d15;
complex<T> d16 = s12 - s45; d16 = T(1)/d16;
complex<T> d21 = spb15*square(spb35); d21 = T(1)/d21;
complex<T> d22 = spb15*spb34*spb35; d22 = T(1)/d22;
complex<T> d26 = spb15*cube(spb35); d26 = T(1)/d26;
complex<T> d28 = spb15*spb35; d28 = T(1)/d28;
complex<T> d29 = spb15*cube(spb13); d29 = T(1)/d29;
complex<T> d30 = spb15*square(spb13); d30 = T(1)/d30;
complex<T> d31 = spb13*spb15*spb34; d31 = T(1)/d31;
complex<T> d36 = spb15*spb23*T(2); d36 = T(1)/d36;
complex<T> d37 = spb15*spb35*T(2); d37 = T(1)/d37;
complex<T> d38 = spb15*cube(spb35)*T(2); d38 = T(1)/d38;
complex<T> t23 = -t48; 
complex<T> t35 = -t54; 
complex<T> t49 = t22*t32; 
complex<T> t50 = spb25*t33; 
complex<T> t55 = spb14*t20; 
complex<T> t59 = spb12*t47; 
complex<T> t63 = spb24*t21; 
complex<T> t72 = spa23*t34; 
complex<T> t75 = spb25*t47; 
complex<T> d4 = s13*t3*t4*T(2); d4 = T(1)/d4;
complex<T> d5 = (s12 - s45)*t3*square(s13); d5 = T(1)/d5;
complex<T> d6 = t3*t4*T(2); d6 = T(1)/d6;
complex<T> d7 = s13*(s12 - s45)*t3; d7 = T(1)/d7;
complex<T> d8 = spb23*spb34*t3; d8 = T(1)/d8;
complex<T> d11 = s35*spb15*t4*T(2); d11 = T(1)/d11;
complex<T> d13 = spb15*spb35*t4*T(2); d13 = T(1)/d13;
complex<T> d14 = spb23*t3; d14 = T(1)/d14;
complex<T> d17 = s13*t3*square(s23 - s45)*T(2); d17 = T(1)/d17;
complex<T> d18 = (s23 - s45)*t3*square(s13); d18 = T(1)/d18;
complex<T> d19 = s13*(s23 - s45)*t3; d19 = T(1)/d19;
complex<T> d20 = spb23*spb34*t3*T(2); d20 = T(1)/d20;
complex<T> d23 = t3*cube(spb13); d23 = T(1)/d23;
complex<T> d24 = t3*square(spb13); d24 = T(1)/d24;
complex<T> d25 = spb13*spb34*t3; d25 = T(1)/d25;
complex<T> d27 = spb34*t3; d27 = T(1)/d27;
complex<T> d32 = spb34*t3*T(2); d32 = T(1)/d32;
complex<T> d33 = t3*T(2); d33 = T(1)/d33;
complex<T> d34 = t3*cube(spb13)*T(2); d34 = T(1)/d34;
complex<T> d35 = spb13*spb34*t3*T(2); d35 = T(1)/d35;
complex<T> d39 = spb34*t48*T(2); d39 = T(1)/d39;
complex<T> t11 = d17*t35*t49 + d18*t49*t54 + d19*t55*t59; 
complex<T> t15 = -(d28*spa34*t19) + d26*s34*spb25*t23 + d21*s34*t75; 
complex<T> t17 = d21*s34*s45*spb24*spb25 - d37*s45*spa34*t19 + d38*s34*s45*spb25*t23 - d36*spa34*spa45*t34; 
complex<T> t18 = d10*t23*t50 + d9*t23*t50 - d1*spb25*t63 + d2*t36*t63; 
complex<T> t41 = d24*s23; 
complex<T> t43 = d13*spa34; 
complex<T> t61 = d14*spa13; 
complex<T> t68 = spb24*t55; 
complex<T> t74 = t22*t35; 
complex<T> t1 = d4*t35*t49 + d10*t48*t50 + d11*t48*t50 + d12*t48*t50 + d9*t48*t50 + d5*t49*t54 + d15*d16*t23*t56 + t23*t43*t56 + d7*t55*t59 + d1*spb25*t63 - d6*spb12*t68 + d2*t21*t75 + d3*t21*t75 - d8*t34*T(2) + d16*spb12*t19*t61*T(2); 
complex<T> t2 = d20*t34 + d18*t35*t49 + d5*t35*t49 + d11*t23*t50 + d12*t23*t50 + d17*t49*t54 + d4*t49*t54 + d15*d16*t48*t56 + t43*t48*t56 + d16*t19*t25*t61 + d3*t36*t63 + d6*spb12*t68 + d19*t25*t68 + d7*t25*t68; 
complex<T> t12 = d25*s23*spb14*t19 + spb14*t41*t59 + d27*t72 + d23*s23*t74; 
complex<T> t14 = -(d22*s12*t19) + d21*s12*spb24*t36 + d25*t19*t38 + d26*s12*spb25*t48 + d24*t38*t59 + d23*s12*t74; 
complex<T> t16 = d22*s45*t19 + d31*spa45*spb14*t19 + d26*s45*spb25*t23 + d30*spa45*spb14*t59 + d29*spa45*t74 + d21*s45*t75; 
complex<T> t46 = d35*s23*t19*t38 + spb12*spb24*t38*t41 + d34*t22*t29*t54; 
complex<T> co1 = d32*s12*t72; 
complex<T> co2 = -(d33*spa34*t72); 
complex<T> co3 = d39*s12*spa15*t34; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t11*(*CI_users[1]->get_value(mc,ind,mu)) + t18*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t14*(*CI_users[4]->get_value(mc,ind,mu)) + t12*(*CI_users[5]->get_value(mc,ind,mu)) + t15*(*CI_users[6]->get_value(mc,ind,mu)) + t16*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + t46*(*CI_users[10]->get_value(mc,ind,mu)) + t17*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qmppqpgam_L_wCI::\
C2q2g1y_qmppqpgam_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmppqpgam_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, p, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmppqpgam L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:07:01 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> s15 = S(1,5);
complex<T> s45 = S(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> t5 = square(spa15); 
complex<T> t7 = square(spa35); 
complex<T> t8 = square(spb23); 
complex<T> t12 = spa15*spa35; 
complex<T> d1 = (s12 - s45)*spa23*spa34; d1 = T(1)/d1;
complex<T> d2 = spa23*spa34*square(s12 - s45)*T(2); d2 = T(1)/d2;
complex<T> d3 = spa12*spa23*spa34*T(2); d3 = T(1)/d3;
complex<T> d4 = spa34*T(2); d4 = T(1)/d4;
complex<T> d5 = spa23*spa34*T(2); d5 = T(1)/d5;
complex<T> d6 = spa12*T(2); d6 = T(1)/d6;
complex<T> d7 = spa12*spa23*T(2); d7 = T(1)/d7;
complex<T> t3 = -(d2*spa12*t7*t8) - d1*spb23*t12*T(2); 
complex<T> t15 = d3*t5; 
complex<T> t4 = d2*spa12*t7*t8 + d1*spb23*t12*T(2) + t15*T(3); 
complex<T> co1 = d4*spb12*spb23*t5; 
complex<T> co2 = d5*s15*spb12*t5; 
complex<T> co3 = d6*spb23*spb34*t5; 
complex<T> co4 = -(d7*s45*spb34*t5); 
complex<T> co5 = s15*s45*t15; 
complex<T> co6 = -(d5*s15*spb12*t5); 
complex<T> co7 = d7*s45*spb34*t5; 
complex<T> co8 = -co2; 
complex<T> co9 = Complex(0,1); 
SeriesC<T> result = co9*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)) + co8*(*CI_users[7]->get_value(mc,ind,mu)) + co7*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g1y_qmmqpmgam_SLC_wCI::\
C2q2g1y_qmmqpmgam_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmmqpmgam_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmmqpmgam SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1y_qmmqppgam_SLC_wCI::\
C2q2g1y_qmmqppgam_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c5, c124));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c2, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmmqppgam_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, qp, p, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmmqppgam SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:08:07 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> s14 = S(1,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s35 = -(spa35*spb35);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = S(3,4);
complex<T> t5 = (s14 - s23)*spb12; 
complex<T> t6 = (s23 - s45)*spb23; 
complex<T> t8 = spb15*T(2); 
complex<T> t10 = square(spb34); 
complex<T> t11 = spb12*spb23; 
complex<T> t16 = -(spb14*T(3)); 
complex<T> t23 = square(spa15); 
complex<T> t24 = square(spa25); 
complex<T> t25 = square(spb13); 
complex<T> t26 = s14 - s23 + s45; 
complex<T> t29 = -(spb14*spb34); 
complex<T> t31 = s14 - s35; 
complex<T> t34 = s14*s45; 
complex<T> t39 = spb12*T(2); 
complex<T> t40 = cube(spa15); 
complex<T> t41 = cube(spa25); 
complex<T> t43 = spb13*spb34; 
complex<T> t45 = s23*spb24; 
complex<T> t49 = spa12*spa23; 
complex<T> t54 = spa15*spb13; 
complex<T> t59 = spb14*spb45; 
complex<T> t60 = spb23*spb24; 
complex<T> t66 = s45*spb13; 
complex<T> t76 = spa15*spb14; 
complex<T> t77 = spb23*spb45; 
complex<T> t99 = spa25*spb34; 
complex<T> d17 = s14 - s23; d17 = T(1)/d17;
complex<T> d26 = (s12 + s15 - s34)*spb15*spb23; d26 = T(1)/d26;
complex<T> d27 = spb15*spb23*spb35; d27 = T(1)/d27;
complex<T> d34 = (s12 + s15 - s34)*spb23; d34 = T(1)/d34;
complex<T> d43 = (s12 + s15 - s34)*spb23*T(2); d43 = T(1)/d43;
complex<T> d46 = spb35*spb45*T(2); d46 = T(1)/d46;
complex<T> t27 = -t77; 
complex<T> t28 = -s23 + t31; 
complex<T> t36 = -t49; 
complex<T> t42 = -t59; 
complex<T> t44 = t11*T(2); 
complex<T> t57 = t25*t40; 
complex<T> t68 = spb34*t24; 
complex<T> t83 = t23*t29; 
complex<T> t88 = spa25*t10; 
complex<T> t89 = spb45*t60; 
complex<T> t90 = spb14*t23; 
complex<T> t102 = s35*t45; 
complex<T> t106 = s34*t10; 
complex<T> d2 = spb23*t26*t5; d2 = T(1)/d2;
complex<T> d3 = t39*square(t31); d3 = T(1)/d3;
complex<T> d7 = spb23*t5*square(t26); d7 = T(1)/d7;
complex<T> d13 = spb25*t39*square(s14 - s23); d13 = T(1)/d13;
complex<T> d14 = spb35*t5*t77*T(2); d14 = T(1)/d14;
complex<T> d15 = spb45*t11; d15 = T(1)/d15;
complex<T> d16 = spb25*spb35*t39; d16 = T(1)/d16;
complex<T> d18 = spb12*t26*t6; d18 = T(1)/d18;
complex<T> d19 = spb12*spb15*t6; d19 = T(1)/d19;
complex<T> d20 = spb35*spb45*t11; d20 = T(1)/d20;
complex<T> d21 = spb12*t6*square(t26); d21 = T(1)/d21;
complex<T> d23 = t11*t8*square(s23 - s45); d23 = T(1)/d23;
complex<T> d24 = spb12*spb45*t6*t8; d24 = T(1)/d24;
complex<T> d25 = spb45*t11*t8; d25 = T(1)/d25;
complex<T> d28 = t11*square(t26); d28 = T(1)/d28;
complex<T> d31 = spb35*t11*t26; d31 = T(1)/d31;
complex<T> d32 = t11*cube(t26); d32 = T(1)/d32;
complex<T> d35 = spb12*square(t26); d35 = T(1)/d35;
complex<T> d36 = spb12*spb35*t26; d36 = T(1)/d36;
complex<T> d37 = spb12*cube(t26); d37 = T(1)/d37;
complex<T> d38 = spb15*spb35*t11; d38 = T(1)/d38;
complex<T> d40 = spb15*t11; d40 = T(1)/d40;
complex<T> d41 = spb35*t11; d41 = T(1)/d41;
complex<T> d45 = spb45*t8; d45 = T(1)/d45;
complex<T> d51 = spb35*t11*t8; d51 = T(1)/d51;
complex<T> t17 = d24*(spa12*spb14*spb23 + spb45*t54); 
complex<T> t20 = d14*(spa25*t27 + spb35*t76); 
complex<T> t50 = d13*spa35; 
complex<T> t52 = d45*spb14*t10*t49 + d46*t36*cube(spb34); 
complex<T> t61 = -(d15*T(2)); 
complex<T> t69 = t57*t59; 
complex<T> t71 = d15*T(2); 
complex<T> t75 = -t88; 
complex<T> t81 = d27*spa12*spb13*t10 + d26*s12*t88; 
complex<T> t82 = t27*t41; 
complex<T> t86 = d23*t25; 
complex<T> d1 = t44*square(s14 - s23); d1 = T(1)/d1;
complex<T> d4 = t28*t5; d4 = T(1)/d4;
complex<T> d5 = spb12*t28*t31; d5 = T(1)/d5;
complex<T> d6 = spb35*spb45*t44; d6 = T(1)/d6;
complex<T> d8 = t26*t44*square(s14 - s23); d8 = T(1)/d8;
complex<T> d9 = t5*square(t28); d9 = T(1)/d9;
complex<T> d10 = spb12*t31*square(t28); d10 = T(1)/d10;
complex<T> d11 = t28*t39*square(s14 - s23); d11 = T(1)/d11;
complex<T> d12 = t28*t39*square(t31); d12 = T(1)/d12;
complex<T> d22 = t26*t44*square(s23 - s45); d22 = T(1)/d22;
complex<T> d29 = spb12*square(t28); d29 = T(1)/d29;
complex<T> d30 = spb12*spb35*t28; d30 = T(1)/d30;
complex<T> d33 = spb12*cube(t28); d33 = T(1)/d33;
complex<T> d39 = spb12*t28; d39 = T(1)/d39;
complex<T> d42 = t44*square(t26); d42 = T(1)/d42;
complex<T> d44 = t39*square(t28); d44 = T(1)/d44;
complex<T> d47 = t28*t39; d47 = T(1)/d47;
complex<T> d48 = t39*cube(t28); d48 = T(1)/d48;
complex<T> d49 = spb35*t26*t44; d49 = T(1)/d49;
complex<T> d50 = t44*cube(t26); d50 = T(1)/d50;
complex<T> t1 = spb34*t16*t17 + d21*t42*t57 + d4*spb24*t68 + d22*t69 + d7*t69 + d8*t69 + d17*t10*t71*t76 + d2*spb13*t83 + t23*t42*t86 + d11*t41*t89 + d9*t41*t89 + d1*t43*t90 + d18*t43*t90 + d16*d17*t27*t99 + t27*t50*t99 + d19*t43*t76*T(2) - d20*cube(spb34)*T(2) - t10*t20*T(3); 
complex<T> t2 = d7*t42*t57 + d8*t42*t57 + d3*spb24*t68 - d4*spb24*t68 - d5*spb24*t68 + d17*t10*t61*t76 + d10*spb24*t82 + d11*spb24*t82 + d12*spb24*t82 + d9*spb24*t82 + d1*spb13*t83 + d2*t43*t90 + d16*d17*t77*t99 + t50*t77*t99 + d6*cube(spb34) + t10*t20*T(3); 
complex<T> t9 = d22*t42*t57 + d21*t69 + d18*spb13*t83 + t23*t59*t86 - d19*t43*t76*T(2) + d25*spb14*t10*T(3) + spb14*spb34*t17*T(3); 
complex<T> t18 = d32*s45*t69 + d28*t66*t83 - d41*spa45*cube(spb34) + d40*spa45*spb14*square(spb34) + d31*s45*t54*square(spb34) + d38*t66*square(spb34); 
complex<T> t19 = t34*(d50*t69 + d28*t43*t90 + d49*t54*square(spb34)); 
complex<T> t21 = -(d3*spb24*t68) + d5*spb24*t68 + (d10 + d12)*t41*t89; 
complex<T> t22 = -(d29*t45*t68) + d37*spa23*t69 + d33*t45*t82 + d35*spa23*spb13*t83 - d30*s23*spa25*square(spb34) + d36*spa23*t54*square(spb34); 
complex<T> t48 = d29*s35; 
complex<T> t56 = s14*(d31*spa15*spb13*t10 + d29*spb24*t68 + d32*t69 + d28*spb13*t83 + d30*t88 + d33*t41*t89); 
complex<T> t87 = d38*spb13*t106 + d26*s34*t75; 
complex<T> t35 = -t48; 
complex<T> t74 = t45*t48*t68 + d48*t102*t82 + d47*s23*spa35*t88; 
complex<T> t67 = spb24*t35*t68 + d33*s35*spb24*t82 + d39*spa35*t88; 
complex<T> co1 = d34*spa15*t75; 
complex<T> co2 = d42*t16*t23*t34*t43; 
complex<T> co3 = d43*s12*spa15*t75; 
complex<T> co4 = -(d44*t102*t68*T(3)); 
complex<T> co5 = d51*t106*t66; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t21*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t81*(*CI_users[4]->get_value(mc,ind,mu)) + t56*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t22*(*CI_users[7]->get_value(mc,ind,mu)) + t87*(*CI_users[8]->get_value(mc,ind,mu)) + t67*(*CI_users[9]->get_value(mc,ind,mu)) + t18*(*CI_users[10]->get_value(mc,ind,mu)) + co2*(*CI_users[12]->get_value(mc,ind,mu)) + co3*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[15]->get_value(mc,ind,mu)) + t52*(*CI_users[16]->get_value(mc,ind,mu)) + t74*(*CI_users[22]->get_value(mc,ind,mu)) + t19*(*CI_users[23]->get_value(mc,ind,mu)) + co5*(*CI_users[24]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qmpqpmgam_SLC_wCI::\
C2q2g1y_qmpqpmgam_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c35));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c2, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmpqpmgam_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, m, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmpqpmgam SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:08:09 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa14 = SPA(1,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> s12 = S(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> s23 = S(2,3);
complex<T> t12 = square(spb23); 
complex<T> t13 = square(spa45); 
complex<T> t16 = square(spb24); 
complex<T> t17 = square(spb25); 
complex<T> d1 = (s12 - s35)*spb14*spb45; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spb15*spb45; d2 = T(1)/d2;
complex<T> d3 = spb15*spb45*square(s12 - s34)*T(2); d3 = T(1)/d3;
complex<T> d4 = spb14*spb45*square(s12 - s35)*T(2); d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*spb15*spb34*spb45*T(2); d5 = T(1)/d5;
complex<T> d6 = (s12 - s35)*spb14*spb35*spb45*T(2); d6 = T(1)/d6;
complex<T> d7 = spb15*spb34*spb45*T(2); d7 = T(1)/d7;
complex<T> d8 = spb14*spb35*spb45*T(2); d8 = T(1)/d8;
complex<T> d9 = spb35*spb45*T(2); d9 = T(1)/d9;
complex<T> d10 = spb34*spb45*T(2); d10 = T(1)/d10;
complex<T> d11 = spb15*spb45*T(2); d11 = T(1)/d11;
complex<T> d12 = spb14*spb45*T(2); d12 = T(1)/d12;
complex<T> t24 = d1*spb24; 
complex<T> t25 = d2*spb25; 
complex<T> t4 = -(d3*spb34*t13*t17) + spa45*spb23*t25*T(2) - d5*spb23*(spa45*spb25*spb34 + spa15*spb12*spb35)*T(3) - d7*t12*T(3); 
complex<T> t5 = d4*spb35*t13*t16 + spa45*spb23*t24*T(2) + d6*spb23*(spa14*spb12*spb34 - spa45*spb24*spb35)*T(3) + d8*t12*T(3); 
complex<T> t6 = -(d4*spb35*t13*t16) + d3*spb34*t13*t17 - spa45*spb23*t24*T(2) - spa45*spb23*t25*T(2) + d5*spb23*(spa45*spb25*spb34 + spa15*spb12*spb35)*T(3) + d6*(-(spa14*spb12*spb23*spb34*T(3)) + spa45*spb23*spb24*spb35*T(3)); 
complex<T> co1 = -(d9*s12*spa14*t12); 
complex<T> co2 = d10*s12*spa15*t12; 
complex<T> co3 = d11*s23*spa34*t12; 
complex<T> co4 = -(d12*s23*spa35*t12); 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t6*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qmpqppgam_SLC_wCI::\
C2q2g1y_qmpqppgam_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c5, c124));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c35));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c3, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmpqppgam_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, p, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmpqppgam SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:09:13 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s35 = S(3,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s14 = -(spa14*spb14);
complex<T> s34 = -(spa34*spb34);
complex<T> t4 = spa12*spa23; 
complex<T> t6 = spa34*T(2); 
complex<T> t10 = square(spa15); 
complex<T> t16 = -(spa35*T(3)); 
complex<T> t23 = square(spb24); 
complex<T> t24 = square(spb34); 
complex<T> t25 = square(spa13); 
complex<T> t28 = -(spa15*spa35); 
complex<T> t30 = s12 - s45; 
complex<T> t31 = s14 - s35; 
complex<T> t33 = s35*s45; 
complex<T> t40 = cube(spb24); 
complex<T> t41 = cube(spb34); 
complex<T> t43 = spa13*spa15; 
complex<T> t44 = spa23*T(2); 
complex<T> t46 = s12*spa25; 
complex<T> t49 = spb12*spb23; 
complex<T> t55 = s45*spa13; 
complex<T> t59 = spa35*spa45; 
complex<T> t68 = spa12*spa45; 
complex<T> t70 = -(spb34*T(2)); 
complex<T> t77 = spb34*T(2); 
complex<T> t85 = spa15*spb24; 
complex<T> d18 = s12 - s35; d18 = T(1)/d18;
complex<T> d26 = spa14*spa23*spa34; d26 = T(1)/d26;
complex<T> d27 = spa23*square(spa34); d27 = T(1)/d27;
complex<T> d28 = spa23*cube(spa34); d28 = T(1)/d28;
complex<T> d33 = spa12*spa24*spa34; d33 = T(1)/d33;
complex<T> d35 = (s15 - s23 + s45)*spa12*spa34; d35 = T(1)/d35;
complex<T> d36 = spa12*spa24; d36 = T(1)/d36;
complex<T> d42 = spa14*spa45*T(2); d42 = T(1)/d42;
complex<T> d45 = spa12*spa24*T(2); d45 = T(1)/d45;
complex<T> t11 = t4*T(2); 
complex<T> t27 = -t68; 
complex<T> t37 = -t49; 
complex<T> t42 = -t59; 
complex<T> t60 = t25*t41; 
complex<T> t67 = spa25*t40; 
complex<T> t69 = spa35*t43; 
complex<T> t74 = spa15*t23; 
complex<T> t75 = t25*t59; 
complex<T> t90 = spa35*t10; 
complex<T> t95 = -(spb14*t10); 
complex<T> t101 = s14*t46; 
complex<T> t106 = s15*t10; 
complex<T> t107 = s23*t10; 
complex<T> d1 = spa14*spa45*t4; d1 = T(1)/d1;
complex<T> d2 = spa24*t44*square(s12 - s35); d2 = T(1)/d2;
complex<T> d3 = (s12 - s35)*spa23*(s12 + t31); d3 = T(1)/d3;
complex<T> d4 = (s12 - s35)*spa23*square(s12 + t31); d4 = T(1)/d4;
complex<T> d5 = (s12 + t31)*t44*square(s12 - s35); d5 = T(1)/d5;
complex<T> d6 = spa34*t30*t4; d6 = T(1)/d6;
complex<T> d8 = s34*(s12 - s35)*t4; d8 = T(1)/d8;
complex<T> d9 = s34*t30*t4; d9 = T(1)/d9;
complex<T> d10 = t4*t6*square(t30); d10 = T(1)/d10;
complex<T> d12 = (s12 - s35)*t4*square(s34); d12 = T(1)/d12;
complex<T> d14 = t30*t4*square(s34); d14 = T(1)/d14;
complex<T> d16 = spa14*spa24*t44; d16 = T(1)/d16;
complex<T> d17 = spa45*t4; d17 = T(1)/d17;
complex<T> d19 = spa45*t30*t4*t6; d19 = T(1)/d19;
complex<T> d20 = t44*square(t31); d20 = T(1)/d20;
complex<T> d21 = spa23*t31*(s12 + t31); d21 = T(1)/d21;
complex<T> d22 = spa23*t31*square(s12 + t31); d22 = T(1)/d22;
complex<T> d23 = (s12 + t31)*t44*square(t31); d23 = T(1)/d23;
complex<T> d25 = spa45*t4*t6; d25 = T(1)/d25;
complex<T> d29 = spa14*spa23*(s12 + t31); d29 = T(1)/d29;
complex<T> d30 = spa23*square(s12 + t31); d30 = T(1)/d30;
complex<T> d31 = spa23*cube(s12 + t31); d31 = T(1)/d31;
complex<T> d32 = spa23*(s12 + t31); d32 = T(1)/d32;
complex<T> d34 = (s15 - s23 + s45)*spa34*t4; d34 = T(1)/d34;
complex<T> d37 = spa14*spa34*t4; d37 = T(1)/d37;
complex<T> d38 = t4*square(spa34); d38 = T(1)/d38;
complex<T> d39 = t4*cube(spa34); d39 = T(1)/d39;
complex<T> d40 = spa14*t4; d40 = T(1)/d40;
complex<T> d41 = spa34*t4; d41 = T(1)/d41;
complex<T> d43 = spa45*t6; d43 = T(1)/d43;
complex<T> d44 = t44*square(s12 + t31); d44 = T(1)/d44;
complex<T> d47 = (s12 + t31)*t44; d47 = T(1)/d47;
complex<T> d48 = t44*cube(s12 + t31); d48 = T(1)/d48;
complex<T> d49 = (s15 - s23 + s45)*t4*t6; d49 = T(1)/d49;
complex<T> d50 = spa14*t4*t6; d50 = T(1)/d50;
complex<T> t18 = d22*t27*t67 + d23*t27*t67 + (-d20 + d21)*spa25*t74; 
complex<T> t36 = d34*spb14; 
complex<T> t47 = d30*s14; 
complex<T> t48 = d2*spb14; 
complex<T> t51 = d31*spa12; 
complex<T> t53 = d43*t49*t90 + d42*t37*cube(spa15); 
complex<T> t57 = s35*(d37*spa13*t10 - d29*spb24*t10 + d38*spa13*t28 + d31*t27*t67 + d30*spa25*t74 + d39*t75); 
complex<T> t72 = d35*spa13*spb14*spb23*t10 + d33*t107; 
complex<T> d7 = t11*square(s12 - s35); d7 = T(1)/d7;
complex<T> d11 = s34*t11*square(s12 - s35); d11 = T(1)/d11;
complex<T> d13 = s34*t11*square(t30); d13 = T(1)/d13;
complex<T> d15 = (s12 - s35)*spa14*spa45*t11; d15 = T(1)/d15;
complex<T> d24 = spa14*spa45*t11; d24 = T(1)/d24;
complex<T> d46 = t11*square(spa34); d46 = T(1)/d46;
complex<T> d51 = t11*cube(spa34); d51 = T(1)/d51;
complex<T> t9 = -(d19*spa15*(spa12*spa35*spb23 + spa13*spa45*spb34)*t16) + d13*t59*t60 + d14*t59*t60 + d9*t24*t69 + d6*t69*t70 + d10*t24*t75 + d25*t90*T(3); 
complex<T> t17 = t33*(d38*t69 + d51*t75 + d50*spa13*square(spa15)); 
complex<T> t19 = d27*spa13*spb12*t28 + spa45*t40*t46*t51 - d30*t46*t74 + d28*spb12*t75 + d26*spa13*spb12*square(spa15) + d29*s12*spb24*square(spa15); 
complex<T> t20 = d15*(-(spa14*spa35*spb34) + spb24*t68); 
complex<T> t22 = d38*t28*t55 + d39*s45*t75 - d40*spb45*cube(spa15) + d41*spa35*spb45*square(spa15) + d37*t55*square(spa15) + t36*t55*square(spa15); 
complex<T> t35 = -t47; 
complex<T> t66 = d48*t101*t40*t68 + t46*t47*t74 + d47*s12*spb24*t95; 
complex<T> t79 = -(d33*t106) + spa13*t106*t36; 
complex<T> t1 = d7*spa13*t24*t28 + d11*t59*t60 + d12*t59*t60 + d22*t67*t68 + d23*t67*t68 + d4*t67*t68 + d5*t67*t68 + d8*t24*t69 + d20*spa25*t74 - d21*spa25*t74 - d3*spa25*t74 + d16*d18*t27*t85 + t48*t68*t85 + d17*d18*t77*t90 + d24*cube(spa15) + t10*t20*T(3); 
complex<T> t2 = d8*spa13*t24*t28 + d9*spa13*t24*t28 + d10*t24*t25*t42 + d11*t42*t60 + d12*t42*t60 + d13*t42*t60 + d14*t42*t60 + d4*t27*t67 + d5*t27*t67 + d7*t24*t69 + d3*spa25*t74 + d6*t69*t77 + t27*t48*t85 + d16*d18*t68*t85 + d17*d18*t70*t90 - d1*cube(spa15)*T(2) - d19*spa15*spa35*(spa12*spa35*spb23 + spa13*spa45*spb34)*T(3) - t10*t20*T(3); 
complex<T> t73 = s14*spa45*t51*t67 + spa25*t35*t74 + d32*spb24*t95; 
complex<T> co1 = -(d36*spb34*t10); 
complex<T> co2 = -(d44*t101*t74*T(3)); 
complex<T> co3 = -(d45*spb34*t107); 
complex<T> co4 = d46*t16*t33*t43; 
complex<T> co5 = d49*spb14*t106*t55; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t18*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t19*(*CI_users[4]->get_value(mc,ind,mu)) + t73*(*CI_users[5]->get_value(mc,ind,mu)) + t79*(*CI_users[6]->get_value(mc,ind,mu)) + t72*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t57*(*CI_users[9]->get_value(mc,ind,mu)) + t22*(*CI_users[10]->get_value(mc,ind,mu)) + t53*(*CI_users[11]->get_value(mc,ind,mu)) + co2*(*CI_users[12]->get_value(mc,ind,mu)) + co3*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + t66*(*CI_users[17]->get_value(mc,ind,mu)) + co5*(*CI_users[18]->get_value(mc,ind,mu)) + t17*(*CI_users[19]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g1y_qmqpmmgam_SLC_wCI::\
C2q2g1y_qmqpmmgam_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmqpmmgam_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmqpmmgam SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1y_qmqpmpgam_SLC_wCI::\
C2q2g1y_qmqpmpgam_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmqpmpgam_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, p, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmqpmpgam SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:24:15 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s14 = S(1,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = -(spa35*spb35);
complex<T> s13 = -(spa13*spb13);
complex<T> t6 = spb23*T(2); 
complex<T> t10 = square(s12 - s45); 
complex<T> t12 = square(s12 + s15 - s34); 
complex<T> t18 = spb15*spb34; 
complex<T> t25 = square(spb24); 
complex<T> t26 = square(spa13); 
complex<T> t27 = square(spa35); 
complex<T> t28 = spb25*spb34; 
complex<T> t29 = spb12*spb23; 
complex<T> t30 = -(s12*spb14); 
complex<T> t31 = spb35*T(2); 
complex<T> t32 = spb15*spb45; 
complex<T> t45 = square(spa15); 
complex<T> t47 = cube(spa13); 
complex<T> t50 = s23 - s45; 
complex<T> t51 = square(spb45); 
complex<T> t55 = cube(spa35); 
complex<T> t56 = -(spb15*spb23*spb24) - spb13*spb24*spb25 - spb14*spb23*spb25*T(2); 
complex<T> t57 = s25 - s34; 
complex<T> t60 = spb13*spb24 + spb14*spb23*T(2); 
complex<T> t61 = s14 - s23; 
complex<T> t62 = s14 - s25; 
complex<T> t69 = s14*s34; 
complex<T> t72 = square(s13); 
complex<T> t73 = square(spb14); 
complex<T> t74 = cube(spa15); 
complex<T> t75 = cube(spb24); 
complex<T> t76 = spb12*spb34; 
complex<T> t79 = -(spb25*T(3)); 
complex<T> t81 = s12*s23; 
complex<T> t82 = s25*spa25; 
complex<T> t83 = spb15*T(2); 
complex<T> t85 = s14*s45; 
complex<T> t92 = spb14*spb24; 
complex<T> t93 = spb34*spb35; 
complex<T> t99 = square(spb34); 
complex<T> t100 = -(spa15*spb45); 
complex<T> t106 = spb12*spb14; 
complex<T> t107 = spb25*spb45; 
complex<T> t110 = spb23*spb35; 
complex<T> t166 = spb25*T(3); 
complex<T> d1 = (s12 - s34)*s35*spb12; d1 = T(1)/d1;
complex<T> d2 = s35*(s12 - s45)*spb12; d2 = T(1)/d2;
complex<T> d3 = (s12 - s45)*spb15*spb23; d3 = T(1)/d3;
complex<T> d4 = (s12 - s34)*s35*spb12*spb15; d4 = T(1)/d4;
complex<T> d5 = s35*(s12 - s45)*spb12*spb15; d5 = T(1)/d5;
complex<T> d10 = (s12 - s34)*spb12*spb35; d10 = T(1)/d10;
complex<T> d11 = (s12 - s45)*spb12*spb35; d11 = T(1)/d11;
complex<T> d14 = spb12*spb35; d14 = T(1)/d14;
complex<T> d16 = s12 - s45; d16 = T(1)/d16;
complex<T> d24 = s13*(s12 - s45)*spb35*spb45; d24 = T(1)/d24;
complex<T> d28 = (s12 - s34)*spb15*square(s35); d28 = T(1)/d28;
complex<T> d31 = (s12 - s45)*spb15*square(s35); d31 = T(1)/d31;
complex<T> d37 = s12 - s34; d37 = T(1)/d37;
complex<T> d50 = spb12*spb35*spb45; d50 = T(1)/d50;
complex<T> d79 = square(spb35); d79 = T(1)/d79;
complex<T> d80 = spb15*square(spb35); d80 = T(1)/d80;
complex<T> d81 = spb15*spb23*square(spb35); d81 = T(1)/d81;
complex<T> d86 = spb35*spb45*square(spb13); d86 = T(1)/d86;
complex<T> d88 = spb35*spb45*cube(spb13); d88 = T(1)/d88;
complex<T> d89 = spb15*cube(spb35); d89 = T(1)/d89;
complex<T> d102 = spb13*spb35*spb45; d102 = T(1)/d102;
complex<T> d109 = spb12*square(spb35); d109 = T(1)/d109;
complex<T> d110 = spb12*spb15*square(spb35); d110 = T(1)/d110;
complex<T> d115 = spb15*square(spb13); d115 = T(1)/d115;
complex<T> d116 = spb15*spb23*square(spb13); d116 = T(1)/d116;
complex<T> d117 = spb15*cube(spb13); d117 = T(1)/d117;
complex<T> d118 = spb35*square(spb13); d118 = T(1)/d118;
complex<T> d120 = spb35*cube(spb13); d120 = T(1)/d120;
complex<T> d129 = spb12*square(spb35)*T(2); d129 = T(1)/d129;
complex<T> t14 = spb25*t57; 
complex<T> t16 = square(t62); 
complex<T> t46 = -t92; 
complex<T> t48 = -t76; 
complex<T> t52 = s12 + t57; 
complex<T> t53 = s34 + t62; 
complex<T> t59 = -t82; 
complex<T> t67 = spa15*spb14*spb25 + spa13*t76; 
complex<T> t78 = spb23*t26; 
complex<T> t84 = spb35*t50; 
complex<T> t94 = t26*(spb13*spb24 + spb14*t6); 
complex<T> t98 = -t106; 
complex<T> t111 = -(d14*T(2)); 
complex<T> t127 = t106*t74; 
complex<T> t129 = -(spb45*t28); 
complex<T> t130 = d14*T(2); 
complex<T> t146 = t47*t99; 
complex<T> t147 = spb45*t55; 
complex<T> t154 = spb25*t45; 
complex<T> t155 = spb24*t73; 
complex<T> t156 = spb34*t47; 
complex<T> t162 = spb14*t29; 
complex<T> t163 = spa15*t25; 
complex<T> t164 = spb34*t26; 
complex<T> t165 = spb45*t27; 
complex<T> t176 = spb45*t45; 
complex<T> t177 = t47*t73; 
complex<T> t178 = spa35*t25; 
complex<T> t185 = spa13*t25; 
complex<T> t210 = t79*t92; 
complex<T> t221 = t75*t82; 
complex<T> t244 = t106*t99; 
complex<T> t262 = t47*t76; 
complex<T> d6 = (s12 - s34)*s35*spb15*t29; d6 = T(1)/d6;
complex<T> d7 = s35*(s12 - s45)*spb15*t29; d7 = T(1)/d7;
complex<T> d8 = (s12 - s34)*spb23*t18; d8 = T(1)/d8;
complex<T> d9 = (s12 - s34)*(s12 + s15 - s34)*spb23*t18; d9 = T(1)/d9;
complex<T> d12 = spb12*spb35*t28; d12 = T(1)/d12;
complex<T> d13 = spb12*t31*square(s12 - s34); d13 = T(1)/d13;
complex<T> d15 = spb45*t110; d15 = T(1)/d15;
complex<T> d17 = s13*(s12 - s45)*t32; d17 = T(1)/d17;
complex<T> d18 = s13*(s12 - s45)*spb23*t32; d18 = T(1)/d18;
complex<T> d19 = spb12*spb45*t18*t6; d19 = T(1)/d19;
complex<T> d20 = s13*t10*t32*T(2); d20 = T(1)/d20;
complex<T> d21 = (s12 - s45)*t32*t72; d21 = T(1)/d21;
complex<T> d22 = t10*t32*t6; d22 = T(1)/d22;
complex<T> d23 = spb35*spb45*t29; d23 = T(1)/d23;
complex<T> d25 = s13*spb45*t10*t31; d25 = T(1)/d25;
complex<T> d26 = (s12 - s45)*spb35*spb45*t72; d26 = T(1)/d26;
complex<T> d27 = t10*t83; d27 = T(1)/d27;
complex<T> d29 = s35*t83*square(s12 - s34); d29 = T(1)/d29;
complex<T> d30 = s35*t10*t83; d30 = T(1)/d30;
complex<T> d32 = spb12*t10*t31; d32 = T(1)/d32;
complex<T> d36 = spb35*t28; d36 = T(1)/d36;
complex<T> d40 = spb35*t6*square(t61); d40 = T(1)/d40;
complex<T> d41 = t110*t61*(s45 + t61); d41 = T(1)/d41;
complex<T> d44 = t110*t61*square(s45 + t61); d44 = T(1)/d44;
complex<T> d45 = spb35*t6*(s45 + t61)*square(t61); d45 = T(1)/d45;
complex<T> d46 = (s15 - s34)*spb23*t18; d46 = T(1)/d46;
complex<T> d47 = (s15 - s34)*(s12 + s15 - s34)*spb23*t18; d47 = T(1)/d47;
complex<T> d49 = t110*t50*(s45 + t61); d49 = T(1)/d49;
complex<T> d51 = spb45*t31; d51 = T(1)/d51;
complex<T> d52 = square(t50); d52 = T(1)/d52;
complex<T> d53 = t32*square(t50)*T(2); d53 = T(1)/d53;
complex<T> d54 = s13*t32*t50; d54 = T(1)/d54;
complex<T> d55 = spb23*t32*t50; d55 = T(1)/d55;
complex<T> d56 = s13*spb23*t32*t50; d56 = T(1)/d56;
complex<T> d57 = s13*t32*square(t50)*T(2); d57 = T(1)/d57;
complex<T> d58 = t32*t50*t72; d58 = T(1)/d58;
complex<T> d59 = spb12*spb35*spb45*t6; d59 = T(1)/d59;
complex<T> d61 = s13*spb45*t31*square(t50); d61 = T(1)/d61;
complex<T> d63 = t110*t50*square(s45 + t61); d63 = T(1)/d63;
complex<T> d64 = spb35*t6*(s45 + t61)*square(t50); d64 = T(1)/d64;
complex<T> d66 = spb34*t31; d66 = T(1)/d66;
complex<T> d67 = spb35*t76; d67 = T(1)/d67;
complex<T> d68 = square(t57); d68 = T(1)/d68;
complex<T> d71 = spb12*t28*t31; d71 = T(1)/d71;
complex<T> d78 = spb23*t12*t18; d78 = T(1)/d78;
complex<T> d82 = t32*square(spb13); d82 = T(1)/d82;
complex<T> d83 = spb23*t32*square(spb13); d83 = T(1)/d83;
complex<T> d84 = spb23*spb45*t18; d84 = T(1)/d84;
complex<T> d85 = t32*cube(spb13); d85 = T(1)/d85;
complex<T> d87 = spb13*spb45*t110; d87 = T(1)/d87;
complex<T> d93 = t110*(s45 + t61); d93 = T(1)/d93;
complex<T> d96 = t110*square(s45 + t61); d96 = T(1)/d96;
complex<T> d98 = t110*cube(s45 + t61); d98 = T(1)/d98;
complex<T> d99 = spb23*spb34*t12; d99 = T(1)/d99;
complex<T> d100 = spb35*(s45 + t61); d100 = T(1)/d100;
complex<T> d101 = spb35*square(s45 + t61); d101 = T(1)/d101;
complex<T> d103 = spb35*cube(s45 + t61); d103 = T(1)/d103;
complex<T> d108 = spb15*spb23*t12; d108 = T(1)/d108;
complex<T> d111 = spb15*t29*square(spb35); d111 = T(1)/d111;
complex<T> d119 = spb13*t110; d119 = T(1)/d119;
complex<T> d121 = t32*square(spb13)*T(2); d121 = T(1)/d121;
complex<T> d122 = t32*cube(spb13)*T(2); d122 = T(1)/d122;
complex<T> d123 = spb45*t31*square(spb13); d123 = T(1)/d123;
complex<T> d124 = spb13*spb45*t31; d124 = T(1)/d124;
complex<T> d125 = spb45*t31*cube(spb13); d125 = T(1)/d125;
complex<T> d130 = spb12*t83*square(spb35); d130 = T(1)/d130;
complex<T> d131 = spb12*spb15*t6*square(spb35); d131 = T(1)/d131;
complex<T> d132 = t83*cube(spb35); d132 = T(1)/d132;
complex<T> d133 = spb34*t12*t6; d133 = T(1)/d133;
complex<T> d137 = spb35*t6*(s45 + t61); d137 = T(1)/d137;
complex<T> d138 = spb35*t6*square(s45 + t61); d138 = T(1)/d138;
complex<T> d139 = spb35*t6*cube(s45 + t61); d139 = T(1)/d139;
complex<T> t8 = square(t53); 
complex<T> t22 = s34*s45*(d132*t129 + d130*t210 - d129*t25 - d131*spb24*t56); 
complex<T> t34 = cube(t53); 
complex<T> t40 = -(d111*t56); 
complex<T> t80 = -t147; 
complex<T> t108 = -t163; 
complex<T> t109 = -t156; 
complex<T> t115 = d51*spa12; 
complex<T> t132 = d66*spa12; 
complex<T> t133 = -(d15*spa13); 
complex<T> t134 = d46*spa25; 
complex<T> t144 = -t155; 
complex<T> t153 = d124*s12*spa23*spb34*t25 + d125*t244*t81 + d123*spb34*t81*t92; 
complex<T> t159 = -(d82*T(3)); 
complex<T> t169 = t107*t127; 
complex<T> t172 = spa35*t111; 
complex<T> t184 = t127*t51; 
complex<T> t192 = t74*t98; 
complex<T> t199 = t154*t46; 
complex<T> t202 = d22*t26; 
complex<T> t235 = t48*t73; 
complex<T> d33 = (s12 - s34)*t52*t93; d33 = T(1)/d33;
complex<T> d34 = (s12 - s34)*t93*square(t52); d34 = T(1)/d34;
complex<T> d35 = spb34*t31*t52*square(s12 - s34); d35 = T(1)/d35;
complex<T> d38 = spb25*t16*t31; d38 = T(1)/d38;
complex<T> d39 = spb25*spb35*t53*t62; d39 = T(1)/d39;
complex<T> d43 = spb25*t16*t31*t53; d43 = T(1)/d43;
complex<T> d48 = t6*t84; d48 = T(1)/d48;
complex<T> d60 = s13*spb45*t84; d60 = T(1)/d60;
complex<T> d62 = spb45*t72*t84; d62 = T(1)/d62;
complex<T> d65 = spb12*spb45*t6*t84; d65 = T(1)/d65;
complex<T> d69 = spb35*t14*t53; d69 = T(1)/d69;
complex<T> d70 = t14*t31; d70 = T(1)/d70;
complex<T> d73 = spb25*t31*t53*square(t57); d73 = T(1)/d73;
complex<T> d74 = t14*t31*t76; d74 = T(1)/d74;
complex<T> d75 = t52*t57*t93; d75 = T(1)/d75;
complex<T> d76 = t57*t93*square(t52); d76 = T(1)/d76;
complex<T> d77 = spb34*t31*t52*square(t57); d77 = T(1)/d77;
complex<T> d90 = t93*square(t52); d90 = T(1)/d90;
complex<T> d91 = spb35*t28*t52; d91 = T(1)/d91;
complex<T> d92 = t93*cube(t52); d92 = T(1)/d92;
complex<T> d95 = spb25*spb35*t53; d95 = T(1)/d95;
complex<T> d105 = spb35*t53; d105 = T(1)/d105;
complex<T> d107 = t52*t93; d107 = T(1)/d107;
complex<T> d112 = spb35*square(t52); d112 = T(1)/d112;
complex<T> d113 = spb25*spb35*t52; d113 = T(1)/d113;
complex<T> d114 = spb35*cube(t52); d114 = T(1)/d114;
complex<T> d127 = spb25*t31*t53; d127 = T(1)/d127;
complex<T> d134 = spb34*t31*square(t52); d134 = T(1)/d134;
complex<T> d135 = spb34*t31*t52; d135 = T(1)/d135;
complex<T> d136 = spb34*t31*cube(t52); d136 = T(1)/d136;
complex<T> t4 = d8*spb14*t108 + d25*t106*t146 + d27*spb24*t165 + d32*spb23*spb24*t165 - d10*t178 - d11*t178 + d3*t178 + d16*t130*t178 + d37*t130*t178 + d15*d16*spb34*t185 + d9*t221 + d36*d37*t100*t25 + d4*t210*t27 + d5*t210*t27 - d1*t25*t27 - d2*t25*t27 + d28*t147*t28 + d29*t147*t28 + d30*t147*t28 + d31*t147*t28 + d13*spb24*t27*t28 + d20*t177*t48 + d34*t192*t51 + d35*t192*t51 - d6*spb24*t27*t56 - d7*spb24*t27*t56 - d12*t75 - d23*t75 + d21*t177*t76 + t202*t46*t76 + d24*t164*t92 + d33*t176*t92 + d18*t26*t60*t92 + d26*t146*t98 - d17*t155*t26*T(3) + d19*spb14*t75*T(3); 
complex<T> t19 = d121*spa23*spb24*t30*t60 + d122*t235*t81 - d121*t155*t81*T(3); 
complex<T> t20 = spa23*(d100*t163 + d103*t169 + d101*t199 + d102*spb34*t25 + d82*t46*t60) + s23*(t155*t159 + d85*t235 + d88*t244 + d86*spb34*t92); 
complex<T> t21 = d89*s45*t129 + d93*s45*t163 + d98*s45*t169 + d96*s45*t199 + d110*s45*t210 + d117*spa45*t235 + d120*spa45*t244 - d109*s45*t25 - d119*spa45*spb34*t25 + s45*spb24*t40 + d118*spa45*spb34*t92 + d116*spa45*t60*t92 - d115*spa45*t155*T(3); 
complex<T> t23 = d80*spa12*t210 - d79*spa12*t25 + d90*spb24*t176*t30 - d81*spa12*spb24*t56 - d15*spa12*t75 - d36*spa12*t75 + d84*spa12*spb14*t75 + s12*(t155*t159 + d91*spb45*t163 + d92*t184 + d85*t235 + d88*t244 - d87*spb34*t25 + d89*spb45*t28 + d78*t59*t75 + d86*spb34*t92 + d83*t60*t92); 
complex<T> t42 = d65*(spa13*spb14*spb23 + spb12*t100); 
complex<T> t87 = -t115; 
complex<T> t141 = d136*s12*s25*t184 + d135*s12*spa25*t100*t25 + d134*s25*spb24*t176*t30; 
complex<T> t175 = (d137*t163 + d139*t169 + d138*t199)*t85; 
complex<T> t191 = d47*t221 + t134*t75; 
complex<T> d42 = spb25*spb35*t62*t8; d42 = T(1)/d42;
complex<T> d72 = spb35*t14*t8; d72 = T(1)/d72;
complex<T> d94 = spb25*spb35*t8; d94 = T(1)/d94;
complex<T> d97 = spb25*spb35*t34; d97 = T(1)/d97;
complex<T> d104 = spb35*t8; d104 = T(1)/d104;
complex<T> d106 = spb35*t34; d106 = T(1)/d106;
complex<T> d126 = spb25*t31*t8; d126 = T(1)/d126;
complex<T> d128 = spb25*t31*t34; d128 = T(1)/d128;
complex<T> t1 = d68*spb14*t108*t132 + d67*d68*t144*t154 + d72*t109*t162 + d42*t156*t162 + d43*t156*t162 + d73*t156*t162 - d70*t185 + d76*t192*t51 + d77*t192*t51 - d71*t75 + d38*t46*t78 + d69*t46*t78 + d75*t176*t92 + d39*t78*t92 + d74*t25*t67*T(3); 
complex<T> t2 = d48*t108 + d26*t106*t146 + d62*t106*t146 - d27*spb24*t165 - d32*spb23*spb24*t165 + d63*t169 + d11*t178 - d3*t178 - d55*spb14*t185 + d52*spb14*t115*t185 + d64*t107*t192 + d49*t199 + d16*spb34*t133*t25 + d16*t172*t25 + d53*t155*t26 + d2*t25*t27 + d24*t164*t46 + d60*t164*t46 + d21*t177*t48 + d58*t177*t48 + d30*t129*t55 + d31*t129*t55 + d7*spb24*t27*t56 + d18*t26*t46*t60 + d56*t26*t46*t60 + d20*t177*t76 + d57*t177*t76 + d50*d52*t155*t78 + d5*t166*t27*t92 + t202*t76*t92 + d25*t146*t98 + d61*t146*t98 + d17*t155*t26*T(3) + d54*t155*t26*T(3) - t25*t42*T(3); 
complex<T> t3 = d61*t106*t146 + d48*t163 + d44*t169 + d45*t169 + d64*t169 + d55*spb14*t185 + d63*t107*t192 + d41*t199 + d53*t144*t26 + d57*t177*t48 - d59*t75 + d58*t177*t76 + d50*d52*t144*t78 + d52*spb14*t185*t87 + d40*t154*t92 + d49*t154*t92 + d60*t164*t92 + d56*t26*t60*t92 + d62*t146*t98 - d54*t155*t26*T(3) + t25*t42*T(3); 
complex<T> t5 = d67*d68*t154*t155 + d73*t109*t162 + d72*t156*t162 + d8*spb14*t163 + d36*d37*spb45*t163 + d68*spb14*t132*t163 + d10*t178 + d34*t184 + d35*t184 + d76*t184 + d77*t184 + d70*t185 + d37*t172*t25 + d1*t25*t27 - d13*spb24*t27*t28 + d33*t176*t46 + d75*t176*t46 + d28*t129*t55 + d29*t129*t55 + d6*spb24*t27*t56 - t134*t75 + d47*t59*t75 + d9*t59*t75 + d4*t166*t27*t92 + d69*t78*t92 - d74*t25*t67*T(3); 
complex<T> t24 = d89*s34*t129 + d97*s34*t156*t162 + d113*spa34*spb45*t163 + d114*spa34*t184 - d95*s34*t185 + d110*s34*t210 - d109*s34*t25 + s34*spb24*t40 + d112*spa34*t176*t46 + d108*spa34*t59*t75 + d94*s34*t78*t92; 
complex<T> t43 = s14*(d98*t169 + d96*t199 + d97*spb14*spb23*t262 + d94*spb23*t92*square(spa13) - d95*spa13*square(spb24) + d93*spa15*square(spb24)); 
complex<T> t44 = d92*s25*t184 + d106*spa25*spb14*spb23*t262 + d90*s25*t176*t46 + d104*spa25*spb23*t92*square(spa13) - d105*spa13*spa25*square(spb24) + d107*spa25*t100*square(spb24); 
complex<T> t126 = d42*t109*t162 + d43*t109*t162 + d44*t107*t192 + d45*t107*t192 + d40*t199 + d39*t46*t78 + d41*t154*t92 + d38*t78*t92; 
complex<T> t143 = t69*(d128*t156*t162 - d127*t185 + d126*t78*t92); 
complex<T> co1 = d99*spa15*t221; 
complex<T> co2 = d133*s12*spa15*t221; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t126*(*CI_users[1]->get_value(mc,ind,mu)) + t191*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t1*(*CI_users[4]->get_value(mc,ind,mu)) + t5*(*CI_users[5]->get_value(mc,ind,mu)) + t2*(*CI_users[6]->get_value(mc,ind,mu)) + t23*(*CI_users[7]->get_value(mc,ind,mu)) + t43*(*CI_users[8]->get_value(mc,ind,mu)) + co1*(*CI_users[9]->get_value(mc,ind,mu)) + t20*(*CI_users[10]->get_value(mc,ind,mu)) + t44*(*CI_users[11]->get_value(mc,ind,mu)) + t24*(*CI_users[12]->get_value(mc,ind,mu)) + t21*(*CI_users[13]->get_value(mc,ind,mu)) + t19*(*CI_users[14]->get_value(mc,ind,mu)) + t153*(*CI_users[19]->get_value(mc,ind,mu)) + t143*(*CI_users[20]->get_value(mc,ind,mu)) + t22*(*CI_users[21]->get_value(mc,ind,mu)) + co2*(*CI_users[25]->get_value(mc,ind,mu)) + t141*(*CI_users[26]->get_value(mc,ind,mu)) + t175*(*CI_users[28]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qmqppmgam_SLC_wCI::\
C2q2g1y_qmqppmgam_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c5, c124));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmqppmgam_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, m, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmqppmgam SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:29:28 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = S(2,3);
complex<T> s24 = -(spa24*spb24);
complex<T> s15 = -(spa15*spb15);
complex<T> t5 = spb25*T(2); 
complex<T> t13 = spb15*spb45; 
complex<T> t20 = square(spa45); 
complex<T> t21 = square(spa15); 
complex<T> t22 = -(spb13*spb23); 
complex<T> t23 = spb12*spb35; 
complex<T> t24 = spb14*spb25; 
complex<T> t27 = spb45*T(2); 
complex<T> t39 = square(spb23); 
complex<T> t40 = square(spb13); 
complex<T> t41 = s12 - s34 - s35; 
complex<T> t42 = cube(spa15); 
complex<T> t43 = cube(spa45); 
complex<T> t45 = s12 + s25 - s34; 
complex<T> t50 = s14 - s35; 
complex<T> t51 = s12 + s14 - s35; 
complex<T> t52 = s15 - s34; 
complex<T> t53 = s12 + s15 - s34; 
complex<T> t65 = spb14*T(2); 
complex<T> t66 = cube(spb23); 
complex<T> t67 = -(spb13*T(3)); 
complex<T> t70 = s25*spa25; 
complex<T> t71 = s24*spa24; 
complex<T> t72 = spb14*spb45; 
complex<T> t74 = -(spb25*T(2)); 
complex<T> t86 = spb13*T(3); 
complex<T> t90 = spb12*T(2); 
complex<T> t97 = spb24*spb35; 
complex<T> t114 = spb25*spb34; 
complex<T> t132 = spa14*spb12; 
complex<T> t149 = spa15*spb12; 
complex<T> t155 = -(spb15*spb24); 
complex<T> t161 = spb13*spb24; 
complex<T> d21 = spb12*spb45; d21 = T(1)/d21;
complex<T> d23 = s12 - s34; d23 = T(1)/d23;
complex<T> d25 = s12 - s35; d25 = T(1)/d25;
complex<T> d26 = (s12 - s34)*spb12*spb45; d26 = T(1)/d26;
complex<T> d27 = (s12 - s35)*spb12*spb45; d27 = T(1)/d27;
complex<T> d66 = spb24*spb45; d66 = T(1)/d66;
complex<T> t6 = square(s12 + t50); 
complex<T> t7 = spb14*t41; 
complex<T> t9 = square(s12 + t52); 
complex<T> t44 = -t97; 
complex<T> t46 = -t71; 
complex<T> t47 = -t70; 
complex<T> t48 = spb23*(t155 - t24) + t161*t74; 
complex<T> t49 = spb15*spb23 + spb13*t5; 
complex<T> t56 = spb34*t132 - spa45*t97*T(3); 
complex<T> t68 = -(t23*t42); 
complex<T> t69 = spb34*t43; 
complex<T> t84 = spb23*t21; 
complex<T> t98 = t40*t42; 
complex<T> t110 = spb14*t5; 
complex<T> t112 = spb23*t20; 
complex<T> t128 = spb34*t65; 
complex<T> t129 = -(t40*T(3)); 
complex<T> t146 = spb13*t39; 
complex<T> t152 = spb34*t39; 
complex<T> t159 = spb35*t39; 
complex<T> t160 = spa45*t114; 
complex<T> t166 = -(t20*t39); 
complex<T> t170 = spa15*t70; 
complex<T> t173 = spa45*t39; 
complex<T> t190 = s23*t39; 
complex<T> d1 = (s12 - s34)*spb12*t41; d1 = T(1)/d1;
complex<T> d2 = (s12 - s35)*spb12*t41; d2 = T(1)/d2;
complex<T> d5 = (s12 - s34)*t24; d5 = T(1)/d5;
complex<T> d8 = (s12 - s34)*spb14*spb34*t45; d8 = T(1)/d8;
complex<T> d9 = (s12 - s34)*spb34*t24*t45; d9 = T(1)/d9;
complex<T> d10 = t65*square(s12 - s34); d10 = T(1)/d10;
complex<T> d11 = (s12 - s35)*spb35*t24; d11 = T(1)/d11;
complex<T> d12 = (s12 - s35)*spb35*t24*(s12 + t50); d12 = T(1)/d12;
complex<T> d14 = (s12 - s34)*spb14*spb34*square(t45); d14 = T(1)/d14;
complex<T> d17 = (s12 - s34)*spb14*square(t41); d17 = T(1)/d17;
complex<T> d18 = (s12 - s35)*spb14*square(t41); d18 = T(1)/d18;
complex<T> d22 = spb34*t13*t90; d22 = T(1)/d22;
complex<T> d24 = spb14*t23*t27; d24 = T(1)/d24;
complex<T> d28 = (s12 - s35)*spb35*(s12 + t50)*t72; d28 = T(1)/d28;
complex<T> d29 = (s12 - s34)*spb34*t13*(s12 + t52); d29 = T(1)/d29;
complex<T> d30 = spb14*t27*square(s12 - s35); d30 = T(1)/d30;
complex<T> d31 = t13*square(s12 - s34)*T(2); d31 = T(1)/d31;
complex<T> d32 = (s12 - s34)*spb34*t13*t90; d32 = T(1)/d32;
complex<T> d33 = (s12 - s35)*spb14*t23*t27; d33 = T(1)/d33;
complex<T> d34 = spb35*t24*t50; d34 = T(1)/d34;
complex<T> d35 = spb35*t24*t50*(s12 + t50); d35 = T(1)/d35;
complex<T> d36 = (s14 - s23)*t72; d36 = T(1)/d36;
complex<T> d37 = (s14 - s23)*(-s23 + t50)*t72; d37 = T(1)/d37;
complex<T> d38 = t50*(-s23 + t50)*t72; d38 = T(1)/d38;
complex<T> d39 = spb35*t50*t72; d39 = T(1)/d39;
complex<T> d40 = spb35*t50*(s12 + t50)*t72; d40 = T(1)/d40;
complex<T> d41 = t13*t52; d41 = T(1)/d41;
complex<T> d42 = spb34*t13*t52; d42 = T(1)/d42;
complex<T> d43 = spb34*t13*t52*(s12 + t52); d43 = T(1)/d43;
complex<T> d45 = (s25 - s34)*spb14*spb34*t45; d45 = T(1)/d45;
complex<T> d46 = (s25 - s34)*spb34*t24; d46 = T(1)/d46;
complex<T> d47 = (s25 - s34)*spb34*t24*t45; d47 = T(1)/d47;
complex<T> d48 = (s25 - s34)*spb14*spb34*square(t45); d48 = T(1)/d48;
complex<T> d50 = square(t41); d50 = T(1)/d50;
complex<T> d51 = spb14*square(t41); d51 = T(1)/d51;
complex<T> d52 = t24*square(t41); d52 = T(1)/d52;
complex<T> d53 = spb14*spb34*square(t45); d53 = T(1)/d53;
complex<T> d54 = spb34*t24*square(t45); d54 = T(1)/d54;
complex<T> d56 = spb34*spb35*t24; d56 = T(1)/d56;
complex<T> d57 = spb14*spb34*cube(t45); d57 = T(1)/d57;
complex<T> d58 = spb14*cube(t41); d58 = T(1)/d58;
complex<T> d59 = spb34*t13; d59 = T(1)/d59;
complex<T> d60 = spb35*t72; d60 = T(1)/d60;
complex<T> d64 = spb45*square(s23 - t50); d64 = T(1)/d64;
complex<T> d68 = t72*square(s23 - t50); d68 = T(1)/d68;
complex<T> d69 = spb24*t13; d69 = T(1)/d69;
complex<T> d70 = spb14*square(t45); d70 = T(1)/d70;
complex<T> d71 = spb12*square(t41); d71 = T(1)/d71;
complex<T> d72 = spb12*spb14*square(t41); d72 = T(1)/d72;
complex<T> d73 = t24*square(t45); d73 = T(1)/d73;
complex<T> d74 = spb12*t24*square(t41); d74 = T(1)/d74;
complex<T> d75 = spb14*cube(t45); d75 = T(1)/d75;
complex<T> d81 = spb24*t13*T(2); d81 = T(1)/d81;
complex<T> d82 = spb14*t27*square(s23 - t50); d82 = T(1)/d82;
complex<T> d86 = t90*square(t41); d86 = T(1)/d86;
complex<T> d87 = spb12*t65*square(t41); d87 = T(1)/d87;
complex<T> d89 = t65*cube(t41); d89 = T(1)/d89;
complex<T> t37 = d33*(spb34*t132 + spa45*t44); 
complex<T> t59 = spa15*t23 + t160*T(3); 
complex<T> t77 = d24*t56; 
complex<T> t79 = d39*spa12; 
complex<T> t93 = t39*(d36*spa25 + d37*t47); 
complex<T> t96 = -t112; 
complex<T> t100 = spa14*t46; 
complex<T> t103 = -(d42*spa12); 
complex<T> t104 = -(d36*spa25); 
complex<T> t113 = t23*t98; 
complex<T> t115 = -(d53*T(3)); 
complex<T> t118 = -(d34*spa24); 
complex<T> t124 = t190*(d69 + d68*t70); 
complex<T> t126 = -t173; 
complex<T> t133 = d72*spb24; 
complex<T> t134 = d74*t48; 
complex<T> t135 = spa25*t49; 
complex<T> t139 = t69*t97; 
complex<T> t141 = -(d41*spa24); 
complex<T> t168 = d31*spb35; 
complex<T> t177 = t112*t67; 
complex<T> t178 = d30*spb34; 
complex<T> t183 = t40*t84; 
complex<T> t186 = spb13*t84; 
complex<T> t208 = t159*t170; 
complex<T> d3 = (s12 - s34)*spb12*t7; d3 = T(1)/d3;
complex<T> d4 = (s12 - s35)*spb12*t7; d4 = T(1)/d4;
complex<T> d6 = (s12 - s34)*spb12*spb25*t7; d6 = T(1)/d6;
complex<T> d7 = (s12 - s35)*spb12*spb25*t7; d7 = T(1)/d7;
complex<T> d13 = spb34*t110*t23; d13 = T(1)/d13;
complex<T> d15 = t128*t45*square(s12 - s34); d15 = T(1)/d15;
complex<T> d16 = spb34*t110*square(s12 - s34); d16 = T(1)/d16;
complex<T> d19 = t7*square(s12 - s34)*T(2); d19 = T(1)/d19;
complex<T> d20 = t7*square(s12 - s35)*T(2); d20 = T(1)/d20;
complex<T> d44 = t128*square(s25 - s34); d44 = T(1)/d44;
complex<T> d49 = t128*t45*square(s25 - s34); d49 = T(1)/d49;
complex<T> d55 = spb35*t24*t6; d55 = T(1)/d55;
complex<T> d61 = spb35*t6*t72; d61 = T(1)/d61;
complex<T> d62 = spb34*t13*t9; d62 = T(1)/d62;
complex<T> d63 = spb25*spb35*t6; d63 = T(1)/d63;
complex<T> d65 = spb35*spb45*t6; d65 = T(1)/d65;
complex<T> d67 = spb34*spb45*t9; d67 = T(1)/d67;
complex<T> d76 = t13*t9; d76 = T(1)/d76;
complex<T> d77 = t24*t6; d77 = T(1)/d77;
complex<T> d78 = t6*t72; d78 = T(1)/d78;
complex<T> d79 = t128*square(t45); d79 = T(1)/d79;
complex<T> d80 = t128*cube(t45); d80 = T(1)/d80;
complex<T> d83 = spb35*t5*t6; d83 = T(1)/d83;
complex<T> d84 = spb35*t27*t6; d84 = T(1)/d84;
complex<T> d85 = spb34*t27*t9; d85 = T(1)/d85;
complex<T> d88 = spb12*t110*square(t41); d88 = T(1)/d88;
complex<T> t1 = d14*t113 + d15*t113 + d48*t113 + d49*t113 + d42*spa12*t146 - d46*spa15*t146 + d26*t173 + d5*t173 + d44*t183 + t114*t168*t20 + d16*t21*t22*t23 + d41*spa24*t39 + d1*t20*t39 + d29*t159*t47 + d43*t159*t47 + d47*t186*t49 + d9*t186*t49 + d17*t44*t69 + d19*t44*t69 + d6*t112*(spb23*(t155 - t24) + t161*t74) + d45*t129*t84 + d8*t129*t84 + d3*spb24*t112*t86 + d32*spb23*(t160 + spa15*t23)*t86 + d10*spb34*t96 + d23*(d21*t126 + d22*t22*(spa15*t23 + t160*T(3))); 
complex<T> t2 = d10*spb34*t112 + d26*t126 + d27*t126 + d5*t126 + d17*t139 + d18*t139 + d19*t139 + d20*t139 + d11*spa14*t146 + d1*t166 + d2*t166 + d3*spb24*t177 + d4*spb24*t177 - t114*t168*t20 + d16*t186*t23 + d28*t152*t46 + d9*t21*t22*t49 + d12*t46*t66 + d32*spb23*(t160 + spa15*t23)*t67 + d24*t39*t67 + d13*t66*t67 + d14*t40*t68 + d15*t40*t68 + d29*t159*t70 + d25*(d21*t173 + t22*t77) + spb23*t37*t86 + d22*t39*t86 + d6*(spb23*(t155 - t24) + t161*t74)*t96 + d7*(spb23*(t155 - t24) + t161*t74)*t96 + t178*t20*t97 + d8*t183*T(3) + d23*(d21*t173 + d22*spb13*spb23*(spa15*t23 + t160*T(3))); 
complex<T> t3 = -(d11*spa14*t146) + d27*t173 + d2*t20*t39 + t178*t20*t44 + d38*t39*t47 + d34*spa24*t66 + spb23*t37*t67 + d18*t44*t69 + d20*t44*t69 + d28*t152*t71 + d40*t152*t71 + d12*t66*t71 + d35*t66*t71 + d7*t112*(spb23*(t155 - t24) + t161*t74) + d25*(d21*t126 + spb13*spb23*t77) - t146*t79 + d4*spb24*t112*t86; 
complex<T> t14 = d46*spa15*t146 - d44*t183 + d48*t40*t68 + d49*t40*t68 + d47*t21*t22*(spb15*spb23 + spb13*spb25*T(2)) + d45*t183*T(3); 
complex<T> t15 = d57*s12*t113 + d59*spa12*t146 - d60*spa12*t146 + d50*spa12*t166 + d51*spa12*spb24*t177 + s12*t115*t183 + d62*s12*t159*t47 - d56*spa12*spb13*t66 + d58*s12*t44*t69 + d61*s12*t152*t71 + d55*s12*t66*t71 + d52*spa12*(spb23*(t155 - t24) + t161*t74)*t96 + d54*s12*t186*(spb15*spb23 + spb13*spb25*T(2)); 
complex<T> t16 = d75*spa34*t113 + d58*s34*t139 + d71*s34*t166 + s34*t133*t177 + d69*s34*t39 + d76*spa34*t159*t47 + d70*spa34*t129*t84 + d74*s34*(spb23*(t155 - t24) + t161*t74)*t96 + d73*spa34*t186*(spb15*spb23 + spb13*spb25*T(2)); 
complex<T> t17 = d58*s35*t139 + d71*s35*t166 + s35*t133*t177 + d68*s35*t39*t70 + d78*spa35*t152*t71 + d77*spa35*t66*t71 + d74*s35*(spb23*(t155 - t24) + t161*t74)*t96; 
complex<T> t18 = s34*s35*(d89*t139 + d86*t166 + d87*spb24*t177 + d88*(spb23*(t155 - t24) + t161*t74)*t96); 
complex<T> t78 = d79*s12; 
complex<T> t83 = d22*t59; 
complex<T> t111 = s12*t100*(d84*t152 + d83*t66); 
complex<T> t137 = d65*t100*t152 + d63*t100*t66 + d64*spa14*t39*t70; 
complex<T> t145 = t103*t146 + t141*t39 + d43*t159*t70; 
complex<T> t151 = d67*t208 + d66*spa15*t39; 
complex<T> t162 = t135*t21; 
complex<T> t172 = t104*t39 + d40*t152*t46 + t118*t66 + d35*t46*t66 + d37*t39*t70 + d38*t39*t70 + t146*t79; 
complex<T> t95 = d57*s25*t113 + s25*t115*t183 + d53*t162*t22; 
complex<T> t125 = d80*s12*s25*t113 + t162*t22*t78 + s25*t129*t78*t84; 
complex<T> co1 = d81*s34*t190; 
complex<T> co2 = d82*s35*t190*t70; 
complex<T> co3 = d85*s12*t208; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t172*(*CI_users[1]->get_value(mc,ind,mu)) + t145*(*CI_users[2]->get_value(mc,ind,mu)) + t93*(*CI_users[3]->get_value(mc,ind,mu)) + t14*(*CI_users[4]->get_value(mc,ind,mu)) + t1*(*CI_users[5]->get_value(mc,ind,mu)) + t3*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + t137*(*CI_users[8]->get_value(mc,ind,mu)) + t151*(*CI_users[9]->get_value(mc,ind,mu)) + t124*(*CI_users[10]->get_value(mc,ind,mu)) + t95*(*CI_users[11]->get_value(mc,ind,mu)) + t16*(*CI_users[12]->get_value(mc,ind,mu)) + t17*(*CI_users[13]->get_value(mc,ind,mu)) + t125*(*CI_users[15]->get_value(mc,ind,mu)) + co1*(*CI_users[16]->get_value(mc,ind,mu)) + co2*(*CI_users[17]->get_value(mc,ind,mu)) + t111*(*CI_users[22]->get_value(mc,ind,mu)) + co3*(*CI_users[24]->get_value(mc,ind,mu)) + t18*(*CI_users[25]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qmqpppgam_SLC_wCI::\
C2q2g1y_qmqpppgam_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c5, c124));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c35));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c2, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c3, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmqpppgam_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, p, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmqpppgam SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:38:11 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = -(spa14*spb14);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s25 = S(2,5);
complex<T> s34 = -(spa34*spb34);
complex<T> t6 = spa14*T(2); 
complex<T> t8 = square(s12 - s35); 
complex<T> t9 = spa23*spa34; 
complex<T> t21 = square(spa15); 
complex<T> t22 = square(spb24); 
complex<T> t23 = square(spb34); 
complex<T> t24 = spa12*spa45; 
complex<T> t25 = -(spa25*T(3)); 
complex<T> t26 = spa14*spa23; 
complex<T> t28 = spa34*T(2); 
complex<T> t38 = square(spa25); 
complex<T> t39 = cube(spb24); 
complex<T> t40 = square(spb23); 
complex<T> t43 = square(spa45); 
complex<T> t44 = cube(spa15); 
complex<T> t45 = s14 - s35; 
complex<T> t47 = cube(spb34); 
complex<T> t50 = spa23*T(2); 
complex<T> t51 = s14 - s25; 
complex<T> t52 = s15 - s23 + s45; 
complex<T> t53 = s12*s14; 
complex<T> t65 = -(spa15*spa25); 
complex<T> t66 = cube(spb23); 
complex<T> t70 = spa14*spa34; 
complex<T> t76 = s25*s35; 
complex<T> t97 = spa12*spa25; 
complex<T> t98 = spa13*spa45; 
complex<T> t100 = -(spa35*spb23); 
complex<T> t101 = s14*spb14; 
complex<T> t115 = spa13*spa15; 
complex<T> t163 = spa13*spa35; 
complex<T> t169 = spa13*T(2); 
complex<T> t200 = spa14*spa25; 
complex<T> d15 = (s12 - s35)*spa12*spa34; d15 = T(1)/d15;
complex<T> d16 = (s12 - s45)*spa12*spa34; d16 = T(1)/d16;
complex<T> d17 = s34*(s12 - s35)*spa12; d17 = T(1)/d17;
complex<T> d18 = s34*(s12 - s45)*spa12; d18 = T(1)/d18;
complex<T> d19 = s34*(s12 - s35)*spa12*spa23; d19 = T(1)/d19;
complex<T> d20 = s34*(s12 - s45)*spa12*spa23; d20 = T(1)/d20;
complex<T> d27 = (s12 - s35)*spa23*square(s34); d27 = T(1)/d27;
complex<T> d29 = (s12 - s45)*spa23*square(s34); d29 = T(1)/d29;
complex<T> d31 = spa12*spa34; d31 = T(1)/d31;
complex<T> d32 = s12 - s35; d32 = T(1)/d32;
complex<T> d34 = s12 - s45; d34 = T(1)/d34;
complex<T> d55 = spa12*spa34*spa35; d55 = T(1)/d55;
complex<T> d63 = spa23*cube(spa34); d63 = T(1)/d63;
complex<T> d64 = square(spa34); d64 = T(1)/d64;
complex<T> d65 = spa23*square(spa34); d65 = T(1)/d65;
complex<T> d86 = spa12*square(spa34); d86 = T(1)/d86;
complex<T> d87 = spa12*spa23*square(spa34); d87 = T(1)/d87;
complex<T> d108 = spa12*square(spa34)*T(2); d108 = T(1)/d108;
complex<T> t10 = spa35*t45; 
complex<T> t41 = -t163; 
complex<T> t42 = s12 + t45; 
complex<T> t48 = spa15*spa24 + t200*T(2); 
complex<T> t49 = spa24*t115 + spa15*t26 + spa13*spa25*t6; 
complex<T> t57 = spa12*t100 - spb34*t98*T(3); 
complex<T> t58 = -t76; 
complex<T> t69 = -(spa12*t39); 
complex<T> t71 = -t101; 
complex<T> t73 = -(d31*T(2)); 
complex<T> t82 = spa15*t22; 
complex<T> t99 = -(spb24*t21); 
complex<T> t102 = -(d31*spb34); 
complex<T> t116 = -(t24*t38); 
complex<T> t118 = d31*T(2); 
complex<T> t120 = d55*spa14; 
complex<T> t143 = spa35*t98; 
complex<T> t148 = spa23*t6; 
complex<T> t153 = t39*t43; 
complex<T> t154 = spa13*t40; 
complex<T> t161 = spa25*t21; 
complex<T> t162 = t24*t39; 
complex<T> t184 = t163*t66; 
complex<T> t185 = spa25*t115; 
complex<T> t188 = t23*T(3); 
complex<T> t189 = t115*t25; 
complex<T> d1 = spa12*spa35*t70; d1 = T(1)/d1;
complex<T> d2 = t24*t9*T(2); d2 = T(1)/d2;
complex<T> d4 = (s12 - s45)*spa45*t26; d4 = T(1)/d4;
complex<T> d5 = (s12 - s45)*spa45*t9; d5 = T(1)/d5;
complex<T> d14 = (s12 - s35)*t26; d14 = T(1)/d14;
complex<T> d21 = s34*(s12 - s35)*spa12*t26; d21 = T(1)/d21;
complex<T> d22 = s34*(s12 - s45)*spa12*t26; d22 = T(1)/d22;
complex<T> d23 = t50*t8; d23 = T(1)/d23;
complex<T> d24 = spa12*t28*t8; d24 = T(1)/d24;
complex<T> d25 = t9*square(s12 - s45)*T(2); d25 = T(1)/d25;
complex<T> d26 = s34*t50*t8; d26 = T(1)/d26;
complex<T> d28 = s34*t50*square(s12 - s45); d28 = T(1)/d28;
complex<T> d30 = spa35*t70; d30 = T(1)/d30;
complex<T> d33 = (s12 - s45)*t24*t9*T(2); d33 = T(1)/d33;
complex<T> d35 = spa12*spa34*spa35*t6; d35 = T(1)/d35;
complex<T> d36 = spa34*t45*t6; d36 = T(1)/d36;
complex<T> d37 = spa34*t6*square(t51); d37 = T(1)/d37;
complex<T> d38 = (-s25 + t45)*t51*t70; d38 = T(1)/d38;
complex<T> d39 = t45*(-s25 + t45)*t70; d39 = T(1)/d39;
complex<T> d40 = t51*t70*square(s25 - t45); d40 = T(1)/d40;
complex<T> d41 = t45*t70*square(s25 - t45); d41 = T(1)/d41;
complex<T> d42 = spa34*(-s25 + t45)*t6*square(t51); d42 = T(1)/d42;
complex<T> d43 = spa34*(-s25 + t45)*t6*square(t45); d43 = T(1)/d43;
complex<T> d45 = spa35*t50*square(t45); d45 = T(1)/d45;
complex<T> d54 = spa35*t28; d54 = T(1)/d54;
complex<T> d56 = square(t45); d56 = T(1)/d56;
complex<T> d57 = (s15 - s23)*t9; d57 = T(1)/d57;
complex<T> d58 = (s15 - s23)*t52*t9; d58 = T(1)/d58;
complex<T> d59 = (s23 - s45)*spa45*t9; d59 = T(1)/d59;
complex<T> d60 = (s23 - s45)*t52*t9; d60 = T(1)/d60;
complex<T> d61 = t26*t98; d61 = T(1)/d61;
complex<T> d62 = t9*t98; d62 = T(1)/d62;
complex<T> d66 = t26*square(spa34); d66 = T(1)/d66;
complex<T> d67 = spa45*t9; d67 = T(1)/d67;
complex<T> d68 = spa35*spa45*t26; d68 = T(1)/d68;
complex<T> d75 = spa34*(-s25 + t45); d75 = T(1)/d75;
complex<T> d76 = spa34*square(s25 - t45); d76 = T(1)/d76;
complex<T> d77 = spa34*cube(-s25 + t45); d77 = T(1)/d77;
complex<T> d79 = t9*square(t52); d79 = T(1)/d79;
complex<T> d80 = spa14*t98; d80 = T(1)/d80;
complex<T> d81 = spa34*t98; d81 = T(1)/d81;
complex<T> d82 = spa34*square(t52); d82 = T(1)/d82;
complex<T> d83 = (-s25 + t45)*t70; d83 = T(1)/d83;
complex<T> d84 = t70*square(s25 - t45); d84 = T(1)/d84;
complex<T> d85 = t70*cube(-s25 + t45); d85 = T(1)/d85;
complex<T> d88 = spa12*t26*square(spa34); d88 = T(1)/d88;
complex<T> d95 = spa13*t26; d95 = T(1)/d95;
complex<T> d96 = spa13*t9; d96 = T(1)/d96;
complex<T> d97 = t9*square(t52)*T(2); d97 = T(1)/d97;
complex<T> d100 = t6*t98; d100 = T(1)/d100;
complex<T> d101 = t28*t98; d101 = T(1)/d101;
complex<T> d102 = spa34*(-s25 + t45)*t6; d102 = T(1)/d102;
complex<T> d103 = spa34*t6*square(s25 - t45); d103 = T(1)/d103;
complex<T> d104 = spa34*t6*cube(-s25 + t45); d104 = T(1)/d104;
complex<T> d109 = spa12*t50*square(spa34); d109 = T(1)/d109;
complex<T> d111 = t50*cube(spa34); d111 = T(1)/d111;
complex<T> t5 = square(t42); 
complex<T> t27 = spa35*t42; 
complex<T> t33 = d88*(spa14*spa15*spa23 + spa24*t115 + t169*t200); 
complex<T> t61 = -t120; 
complex<T> t68 = cube(t42); 
complex<T> t81 = d2*t57; 
complex<T> t94 = (d57*spb14 + d58*t101)*t21; 
complex<T> t105 = d24*spa14; 
complex<T> t106 = d54*spb12; 
complex<T> t114 = s12*(d101*t100*t21 + d100*spb23*t44); 
complex<T> t121 = -(d59*spb12); 
complex<T> t123 = -(d57*spb14); 
complex<T> t132 = t41*t66; 
complex<T> t142 = t38*t82; 
complex<T> t149 = d79*t71; 
complex<T> t167 = t162*t38; 
complex<T> t173 = spa45*t41; 
complex<T> t182 = t153*t97; 
complex<T> t187 = t43*t69; 
complex<T> t193 = t154*t65; 
complex<T> t199 = t21*t71; 
complex<T> d3 = spa35*t148*t24; d3 = T(1)/d3;
complex<T> d8 = spa35*t148*t8; d8 = T(1)/d8;
complex<T> d44 = t10*t26; d44 = T(1)/d44;
complex<T> d46 = spa23*t10*t42; d46 = T(1)/d46;
complex<T> d47 = t10*t26*t42; d47 = T(1)/d47;
complex<T> d48 = spa34*t10*t42; d48 = T(1)/d48;
complex<T> d53 = spa12*spa34*t10*t6; d53 = T(1)/d53;
complex<T> d89 = t42*t70; d89 = T(1)/d89;
complex<T> d110 = spa12*t148*square(spa34); d110 = T(1)/d110;
complex<T> t2 = d59*spb12*t161 - d4*spb23*t161 + d20*t185*t188 + d60*t199 + d5*spa35*spb13*t21 + d59*spa35*spb13*t21 + d16*spb34*t21 + d25*t143*t23 + d18*t21*t23 - d4*spb13*t44 + d28*t143*t47 + d29*t143*t47 - d22*spa15*t23*(spa24*t115 + spa15*t26 + spa13*spa25*t6) + d34*(t102*t21 + spa15*spa25*t81) + d33*spa15*t25*(spa12*t100 - spb34*t98); 
complex<T> t14 = s35*s45*(d111*t143 + d109*t189 - d108*t21 + d110*spa15*(spa24*t115 + t169*t200 + spa15*t26)); 
complex<T> t35 = d53*(spa12*spa35*spb23 - spb24*t200); 
complex<T> t77 = -t106; 
complex<T> t129 = t121*t161 + d58*t199 + (-(d59*spa35*spb13) + d60*t101 + t123)*t21; 
complex<T> t130 = d37*t193 + d38*t185*t40 + (d40 + d42)*t184*t97; 
complex<T> t141 = d82*spb23*t199 + d81*t100*t21 + d80*spb23*t44; 
complex<T> t178 = t132*t97; 
complex<T> t234 = t149*t21; 
complex<T> d6 = (s12 - s35)*spa23*t27; d6 = T(1)/d6;
complex<T> d7 = (s12 - s35)*t26*t27; d7 = T(1)/d7;
complex<T> d9 = (s12 - s35)*spa34*t27; d9 = T(1)/d9;
complex<T> d10 = (s12 - s35)*spa23*spa35*t5; d10 = T(1)/d10;
complex<T> d11 = t27*t50*t8; d11 = T(1)/d11;
complex<T> d12 = (s12 - s35)*spa34*spa35*t5; d12 = T(1)/d12;
complex<T> d13 = t27*t28*t8; d13 = T(1)/d13;
complex<T> d49 = spa23*t10*t5; d49 = T(1)/d49;
complex<T> d50 = t27*t50*square(t45); d50 = T(1)/d50;
complex<T> d51 = spa34*t10*t5; d51 = T(1)/d51;
complex<T> d52 = t27*t28*square(t45); d52 = T(1)/d52;
complex<T> d69 = t27*t70; d69 = T(1)/d69;
complex<T> d70 = spa23*spa35*t5; d70 = T(1)/d70;
complex<T> d71 = spa35*t26*t5; d71 = T(1)/d71;
complex<T> d72 = spa34*spa35*t5; d72 = T(1)/d72;
complex<T> d73 = spa23*spa35*t68; d73 = T(1)/d73;
complex<T> d74 = spa34*spa35*t68; d74 = T(1)/d74;
complex<T> d78 = spa34*t27; d78 = T(1)/d78;
complex<T> d90 = spa23*t5; d90 = T(1)/d90;
complex<T> d91 = t26*t5; d91 = T(1)/d91;
complex<T> d92 = spa34*t5; d92 = T(1)/d92;
complex<T> d93 = spa23*t68; d93 = T(1)/d93;
complex<T> d94 = spa34*t68; d94 = T(1)/d94;
complex<T> d98 = spa35*t5*t50; d98 = T(1)/d98;
complex<T> d99 = spa35*t50*t68; d99 = T(1)/d99;
complex<T> d105 = t27*t28; d105 = T(1)/d105;
complex<T> d106 = spa35*t28*t5; d106 = T(1)/d106;
complex<T> d107 = spa35*t28*t68; d107 = T(1)/d107;
complex<T> t1 = d4*spb23*t161 - d12*t182 - d13*t182 - d5*spa35*spb13*t21 + d14*spb34*t21 - d15*spb34*t21 - d16*spb34*t21 + d23*spa15*spa35*t23 + spa15*spa35*t105*t23 + d25*t173*t23 + d19*t189*t23 + d20*t189*t23 - d17*t21*t23 - d18*t21*t23 + d10*t116*t39 + d11*t116*t39 - d1*t44 + d4*spb13*t44 + d3*t25*t44 + d26*t173*t47 + d27*t173*t47 + d28*t173*t47 + d29*t173*t47 + d21*spa15*t23*(spa24*t115 + spa15*t26 + spa13*spa25*t6) + d22*spa15*t23*(spa24*t115 + spa15*t26 + spa13*spa25*t6) + d7*t22*t48*t65 + d34*(d31*spb34*t21 + t65*t81) + d9*spa25*spa45*t82 + d8*spa25*t24*t82 + d32*(spb34*t118*t21 + d30*spa45*t99) + d6*t142*T(3) + d2*t161*T(3) + d33*spa15*spa25*(spa12*t100 - spb34*t98)*T(3); 
complex<T> t3 = -(d45*t142) + d44*spb24*t161 + d40*t178 + d41*t178 + d42*t178 + d43*t178 - d51*t182 - d52*t182 + d38*t193 + d39*t193 + d36*spb23*t21 + d49*t116*t39 + d50*t116*t39 + d37*t185*t40 - d35*t44 + d56*t142*t61 + d47*t22*t48*t65 + d56*spb24*t161*t77 + d48*spa25*spa45*t82 + d46*t142*T(3) - t21*t35*T(3); 
complex<T> t4 = d45*t142 + d56*(t120*t142 + spb24*t106*t161) + d10*t167 + d11*t167 + d49*t167 + d50*t167 + d12*t182 + d13*t182 + d51*t182 + d52*t182 + d19*t185*t188 - d36*spb23*t21 - d14*spb34*t21 + d15*spb34*t21 - d23*spa15*spa35*t23 - spa15*spa35*t105*t23 + d17*t21*t23 + d39*t185*t40 + d26*t143*t47 + d27*t143*t47 - d21*spa15*t23*(spa24*t115 + spa15*t26 + spa13*spa25*t6) + d48*spa45*t22*t65 + d9*spa45*t22*t65 + d8*t22*t24*t65 + d32*t21*(d30*spa45*spb24 + spb34*t73) + d47*spa25*t48*t82 + d7*spa25*t48*t82 + d41*t184*t97 + d43*t184*t97 + d44*spa25*t99 - d46*t142*T(3) - d6*t142*T(3) + t21*t35*T(3); 
complex<T> t15 = d99*t167*t53 + d98*s12*spb14*t22*t48*t65 - d98*t142*t53*T(3); 
complex<T> t19 = d63*s35*t143 + d93*spb35*t167 + d85*s35*t178 + d94*spb35*t182 + d87*s35*t189 + d84*s35*t193 - d86*s35*t21 - d83*s35*spb23*t21 + d89*spa45*spb24*spb35*t21 + s35*spa15*t33 + d92*spa45*spb35*t22*t65 + d91*spa25*spb35*t48*t82 - d90*spb35*t142*T(3); 
complex<T> t20 = d63*s45*t143 + d87*s45*t189 - d86*s45*t21 + d96*spa35*spb45*t21 + s45*t234 + s45*spa15*t33 - d95*spb45*t44; 
complex<T> t34 = s25*(d85*t178 + d84*t193 - d83*spb23*square(spa15)); 
complex<T> t96 = d107*t182*t53 + d106*spa45*t22*t53*t65 + d105*s12*spa45*spb14*t99; 
complex<T> t125 = d72*t22; 
complex<T> t134 = -(d70*T(3)); 
complex<T> t151 = d102*spb23*t21*t58 + d103*t185*t40*t58 + d104*t178*t76; 
complex<T> t16 = d66*spa15*spb12*(spa24*t115 + t169*t200 + spa15*t26) + spb12*(d67*t161 + d65*t189 - d64*t21 - d30*t44 - d68*spa25*t44) + s12*(t134*t142 + d73*t167 + d63*t173 + d74*t182 + d62*spa35*t21 + d69*spa45*spb24*t21 - d61*t44 + spa45*t125*t65 + d71*spa25*t48*t82); 
complex<T> t17 = s14*(t134*t142 + d73*t167 + d74*t182 + spa45*t125*t65) + spb14*(d77*t178 + d76*t193 - d75*spb23*t21 + d70*t22*t48*t65 + d78*spa45*t99); 
complex<T> co1 = s15*t234; 
complex<T> co2 = d97*s15*s45*t199; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t94*(*CI_users[2]->get_value(mc,ind,mu)) + t129*(*CI_users[3]->get_value(mc,ind,mu)) + t130*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + t2*(*CI_users[6]->get_value(mc,ind,mu)) + t16*(*CI_users[7]->get_value(mc,ind,mu)) + t17*(*CI_users[8]->get_value(mc,ind,mu)) + co1*(*CI_users[9]->get_value(mc,ind,mu)) + t141*(*CI_users[10]->get_value(mc,ind,mu)) + t34*(*CI_users[11]->get_value(mc,ind,mu)) + t19*(*CI_users[12]->get_value(mc,ind,mu)) + t20*(*CI_users[13]->get_value(mc,ind,mu)) + co2*(*CI_users[16]->get_value(mc,ind,mu)) + t15*(*CI_users[17]->get_value(mc,ind,mu)) + t114*(*CI_users[21]->get_value(mc,ind,mu)) + t151*(*CI_users[24]->get_value(mc,ind,mu)) + t96*(*CI_users[26]->get_value(mc,ind,mu)) + t14*(*CI_users[28]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g1y_qmmmqpgam_nf_wCI::\
C2q2g1y_qmmmqpgam_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmmmqpgam_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmmmqpgam nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qmmpqpgam_nf_wCI::\
C2q2g1y_qmmpqpgam_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmmpqpgam_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmmpqpgam nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qmpmqpgam_nf_wCI::\
C2q2g1y_qmpmqpgam_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmpmqpgam_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmpmqpgam nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qmppqpgam_nf_wCI::\
C2q2g1y_qmppqpgam_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmppqpgam_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmppqpgam nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qmgamqpmm_L_wCI::\
C2q2g1y_qmgamqpmm_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmgamqpmm_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmgamqpmm L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1y_qmgamqpmp_L_wCI::\
C2q2g1y_qmgamqpmp_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmgamqpmp_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, gam, qp, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmgamqpmp L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:39:04 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co7*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t18*(*CI_users[2]->get_value(mc,ind,mu)) + t12*(*CI_users[3]->get_value(mc,ind,mu)) + t44*(*CI_users[4]->get_value(mc,ind,mu)) + t15*(*CI_users[5]->get_value(mc,ind,mu)) + t16*(*CI_users[6]->get_value(mc,ind,mu)) + t13*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + co6*(*CI_users[10]->get_value(mc,ind,mu)) + co4*(*CI_users[11]->get_value(mc,ind,mu)) + t17*(*CI_users[12]->get_value(mc,ind,mu)) + t14*(*CI_users[13]->get_value(mc,ind,mu)) + co5*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qmgamqppm_L_wCI::\
C2q2g1y_qmgamqppm_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmgamqppm_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, gam, qp, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmgamqppm L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:39:07 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co11*(t21*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t20*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + t4*(*CI_users[4]->get_value(mc,ind,mu)) + t26*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + co4*(*CI_users[8]->get_value(mc,ind,mu)) + co10*(*CI_users[9]->get_value(mc,ind,mu)) + co6*(*CI_users[10]->get_value(mc,ind,mu)) + co7*(*CI_users[11]->get_value(mc,ind,mu)) + co8*(*CI_users[12]->get_value(mc,ind,mu)) + co9*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qmgamqppp_L_wCI::\
C2q2g1y_qmgamqppp_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qmgamqppp_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, gam, qp, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmgamqppp L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:39:08 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co4*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g1y_qmgamqpmm_nf_wCI::\
C2q2g1y_qmgamqpmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmgamqpmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmgamqpmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qmgamqpmp_nf_wCI::\
C2q2g1y_qmgamqpmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmgamqpmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmgamqpmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qmgamqppm_nf_wCI::\
C2q2g1y_qmgamqppm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmgamqppm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmgamqppm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qmgamqppp_nf_wCI::\
C2q2g1y_qmgamqppp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qmgamqppp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qmgamqppp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qpppqmgap_L_wCI::\
C2q2g1y_qpppqmgap_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpppqmgap_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpppqmgap L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1y_qppmqmgap_L_wCI::\
C2q2g1y_qppmqmgap_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qppmqmgap_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, p, m, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qppmqmgap L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:39:11 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = S(3,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> t1 = spa12*spa15; 
complex<T> t6 = square(spa34); 
complex<T> t10 = square(spa14); 
complex<T> t13 = spa13*spb12; 
complex<T> t21 = square(spb12); 
complex<T> t25 = spa14*spa34; 
complex<T> t27 = -(spb12*spb25); 
complex<T> d4 = (s12 + s15 - s34)*spa45; d4 = T(1)/d4;
complex<T> d5 = spa12*spa24*spa45; d5 = T(1)/d5;
complex<T> d6 = (s12 + s15 - s34)*spa12*spa45; d6 = T(1)/d6;
complex<T> d9 = spa15*spa45*T(2); d9 = T(1)/d9;
complex<T> d13 = spa23*spa45*T(2); d13 = T(1)/d13;
complex<T> d14 = (s12 + s15 - s34)*spa45*T(2); d14 = T(1)/d14;
complex<T> t8 = -t13; 
complex<T> t14 = s34*t6; 
complex<T> t19 = d6*spb25; 
complex<T> t22 = spa13*t10; 
complex<T> t32 = s15*t6; 
complex<T> d1 = (s23 - s45)*spa45*t1; d1 = T(1)/d1;
complex<T> d2 = spa45*t1*square(s23 - s45)*T(2); d2 = T(1)/d2;
complex<T> d3 = spa23*spa45*t1*T(2); d3 = T(1)/d3;
complex<T> d7 = spa24*spa45*t1; d7 = T(1)/d7;
complex<T> d8 = spa45*t1; d8 = T(1)/d8;
complex<T> d10 = spa24*spa45*t1*T(2); d10 = T(1)/d10;
complex<T> d11 = spa45*t1*T(2); d11 = T(1)/d11;
complex<T> d12 = spa23*t1*T(2); d12 = T(1)/d12;
complex<T> t12 = -t19; 
complex<T> t16 = d7*spa14; 
complex<T> t17 = d2*spa23; 
complex<T> t23 = (d10*s23*spa14 + d11*spa13*spb23)*t14; 
complex<T> t33 = t19*t32 + d5*spa14*spb15*t6; 
complex<T> t35 = d14*t27*t32 + d13*spb15*t6*t8; 
complex<T> t4 = -(t17*t21*t22) - d1*t13*t25*T(2) - d3*spa13*t6*T(3); 
complex<T> t5 = (d8*spa13*spb23 + s23*t16)*t6; 
complex<T> t20 = t17*t21*t22 + d1*t13*t25*T(2); 
complex<T> t29 = t14*(t12 + t16); 
complex<T> co1 = d4*t27*t6; 
complex<T> co2 = d9*spb23*t6*t8; 
complex<T> co3 = d12*spa13*spb45*t14; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t20*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + t33*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + t29*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t23*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t35*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qpmpqmgap_L_wCI::\
C2q2g1y_qpmpqmgap_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpmpqmgap_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, p, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpmpqmgap L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:39:48 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t3 = spa15*spa45; 
complex<T> t4 = square(s12 - s45); 
complex<T> t18 = square(spa24); 
complex<T> t19 = square(spb13); 
complex<T> t20 = square(spb35); 
complex<T> t21 = square(spa14); 
complex<T> t24 = -(spa12*T(2)); 
complex<T> t31 = cube(spb13); 
complex<T> t32 = cube(spb35); 
complex<T> t33 = cube(spa24); 
complex<T> t35 = -(spa25*T(2)); 
complex<T> t37 = s12*spa14; 
complex<T> t47 = spa24*T(2); 
complex<T> t48 = spa23*spa45; 
complex<T> t54 = spa12*spa23; 
complex<T> t57 = spa24*spb35; 
complex<T> d7 = spa15*square(s12 - s34)*T(2); d7 = T(1)/d7;
complex<T> d8 = (s12 - s34)*s35*spa15; d8 = T(1)/d8;
complex<T> d9 = s35*(s12 - s45)*spa15; d9 = T(1)/d9;
complex<T> d10 = (s12 - s34)*spa15*square(s35); d10 = T(1)/d10;
complex<T> d11 = s35*spa15*square(s12 - s34)*T(2); d11 = T(1)/d11;
complex<T> d13 = (s12 - s45)*spa15*square(s35); d13 = T(1)/d13;
complex<T> d15 = spa15*spa34*spa35*T(2); d15 = T(1)/d15;
complex<T> d16 = s12 - s45; d16 = T(1)/d16;
complex<T> d21 = spa15*square(spa35); d21 = T(1)/d21;
complex<T> d22 = spa15*spa34*spa35; d22 = T(1)/d22;
complex<T> d26 = spa15*cube(spa35); d26 = T(1)/d26;
complex<T> d28 = spa15*spa35; d28 = T(1)/d28;
complex<T> d29 = spa15*cube(spa13); d29 = T(1)/d29;
complex<T> d30 = spa15*square(spa13); d30 = T(1)/d30;
complex<T> d31 = spa13*spa15*spa34; d31 = T(1)/d31;
complex<T> d36 = spa15*cube(spa35)*T(2); d36 = T(1)/d36;
complex<T> d37 = spa15*spa35*T(2); d37 = T(1)/d37;
complex<T> d38 = spa15*spa23*T(2); d38 = T(1)/d38;
complex<T> t22 = -t48; 
complex<T> t34 = -t54; 
complex<T> t49 = t21*t31; 
complex<T> t50 = spa25*t32; 
complex<T> t55 = spa14*t19; 
complex<T> t59 = spa12*t47; 
complex<T> t64 = spa25*t20; 
complex<T> t66 = spb23*t33; 
complex<T> t77 = spa24*t35; 
complex<T> d1 = spa23*spa34*t3; d1 = T(1)/d1;
complex<T> d2 = t3*t4*T(2); d2 = T(1)/d2;
complex<T> d3 = s13*(s12 - s45)*t3; d3 = T(1)/d3;
complex<T> d4 = s13*t3*t4*T(2); d4 = T(1)/d4;
complex<T> d5 = (s12 - s45)*t3*square(s13); d5 = T(1)/d5;
complex<T> d6 = spa15*spa35*t4*T(2); d6 = T(1)/d6;
complex<T> d12 = s35*spa15*t4*T(2); d12 = T(1)/d12;
complex<T> d14 = spa23*t3; d14 = T(1)/d14;
complex<T> d17 = s13*(s23 - s45)*t3; d17 = T(1)/d17;
complex<T> d18 = s13*t3*square(s23 - s45)*T(2); d18 = T(1)/d18;
complex<T> d19 = (s23 - s45)*t3*square(s13); d19 = T(1)/d19;
complex<T> d20 = spa23*spa34*t3*T(2); d20 = T(1)/d20;
complex<T> d23 = t3*cube(spa13); d23 = T(1)/d23;
complex<T> d24 = t3*square(spa13); d24 = T(1)/d24;
complex<T> d25 = spa13*spa34*t3; d25 = T(1)/d25;
complex<T> d27 = spa34*t3; d27 = T(1)/d27;
complex<T> d32 = spa34*t3*T(2); d32 = T(1)/d32;
complex<T> d33 = t3*T(2); d33 = T(1)/d33;
complex<T> d34 = t3*cube(spa13)*T(2); d34 = T(1)/d34;
complex<T> d35 = spa13*spa34*t3*T(2); d35 = T(1)/d35;
complex<T> d39 = spa34*t48*T(2); d39 = T(1)/d39;
complex<T> t11 = s23*(d34*s12*t21*t34 + d24*spa12*spa24*t37 + d35*t18*t37); 
complex<T> t13 = -(d28*spb34*t18) + s34*spa25*(d26*t22 + d21*t47); 
complex<T> t15 = d10*t22*t50 + d11*t22*t50 - d7*spa24*t64 + d8*t20*t77; 
complex<T> t16 = d22*s45*t18 + d31*spa14*spb45*t18 + d26*s45*spa25*t22 + d29*spb45*t21*t34 + d21*s45*spa25*t47 + d30*spa14*spb45*t59; 
complex<T> t17 = d21*s34*s45*spa24*spa25 - d37*s45*spb34*t18 + d36*s34*s45*spa25*t22 - d38*spb34*spb45*t33; 
complex<T> t42 = d6*spb34; 
complex<T> t45 = d18*t34*t49 + d19*t49*t54 + d17*t55*t59; 
complex<T> t62 = d14*spb13; 
complex<T> t69 = spa24*t55; 
complex<T> t75 = d23*t21; 
complex<T> t1 = d4*t34*t49 + d10*t48*t50 + d11*t48*t50 + d12*t48*t50 + d13*t48*t50 + d5*t49*t54 + d15*d16*t22*t57 + t22*t42*t57 + d3*t55*t59 + d7*spa24*t64 + d8*t47*t64 + d9*t47*t64 - d2*spa12*t69 - d1*t33*T(2) + d16*spa12*t18*t62*T(2); 
complex<T> t2 = d20*t33 + d19*t34*t49 + d5*t34*t49 + d12*t22*t50 + d13*t22*t50 + d18*t49*t54 + d4*t49*t54 + d15*d16*t48*t57 + t42*t48*t57 + d16*t18*t24*t62 + d2*spa12*t69 + d17*t24*t69 + d3*t24*t69 + d9*t20*t77; 
complex<T> t12 = d25*s23*spa14*t18 + d24*s23*spa14*t59 + d27*t66 + s23*t34*t75; 
complex<T> t46 = -(d22*s12*t18) + d25*t18*t37 + d26*s12*spa25*t48 + d24*t37*t59 + s12*t34*t75 + d21*s12*t77; 
complex<T> co1 = d32*s12*t66; 
complex<T> co2 = -(d33*spb34*t66); 
complex<T> co3 = d39*s12*spb15*t33; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t45*(*CI_users[1]->get_value(mc,ind,mu)) + t15*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t46*(*CI_users[4]->get_value(mc,ind,mu)) + t12*(*CI_users[5]->get_value(mc,ind,mu)) + t13*(*CI_users[6]->get_value(mc,ind,mu)) + t16*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + t11*(*CI_users[10]->get_value(mc,ind,mu)) + t17*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qpmmqmgap_L_wCI::\
C2q2g1y_qpmmqmgap_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpmmqmgap_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, m, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpmmqmgap L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:39:49 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> s45 = S(4,5);
complex<T> s15 = S(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> t5 = square(spb15); 
complex<T> t7 = square(spa23); 
complex<T> t8 = square(spb35); 
complex<T> t12 = spa23*spb15; 
complex<T> d1 = (s12 - s45)*spb23*spb34; d1 = T(1)/d1;
complex<T> d2 = spb23*spb34*square(s12 - s45)*T(2); d2 = T(1)/d2;
complex<T> d3 = spb12*spb23*spb34*T(2); d3 = T(1)/d3;
complex<T> d4 = spb34*T(2); d4 = T(1)/d4;
complex<T> d5 = spb23*spb34*T(2); d5 = T(1)/d5;
complex<T> d6 = spb12*T(2); d6 = T(1)/d6;
complex<T> d7 = spb12*spb23*T(2); d7 = T(1)/d7;
complex<T> t3 = -(d2*spb12*t7*t8) - d1*spb35*t12*T(2); 
complex<T> t15 = d3*t5; 
complex<T> t4 = d2*spb12*t7*t8 + d1*spb35*t12*T(2) + t15*T(3); 
complex<T> co1 = d4*spa12*spa23*t5; 
complex<T> co2 = d5*s15*spa12*t5; 
complex<T> co3 = d6*spa23*spa34*t5; 
complex<T> co4 = -(d7*s45*spa34*t5); 
complex<T> co5 = s15*s45*t15; 
complex<T> co6 = -(d5*s15*spa12*t5); 
complex<T> co7 = d7*s45*spa34*t5; 
complex<T> co8 = -co2; 
complex<T> co9 = Complex(0,1); 
SeriesC<T> result = co9*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)) + co8*(*CI_users[7]->get_value(mc,ind,mu)) + co7*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g1y_qppqmpgap_SLC_wCI::\
C2q2g1y_qppqmpgap_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qppqmpgap_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qppqmpgap SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1y_qppqmmgap_SLC_wCI::\
C2q2g1y_qppqmmgap_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c5, c124));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c2, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qppqmmgap_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, p, qm, m, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qppqmmgap SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:40:59 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = S(1,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s35 = -(spa35*spb35);
complex<T> t5 = (s14 - s23)*spa12; 
complex<T> t6 = (s23 - s45)*spa23; 
complex<T> t8 = spa15*T(2); 
complex<T> t10 = square(spa34); 
complex<T> t11 = spa12*spa23; 
complex<T> t16 = -(spa14*T(3)); 
complex<T> t23 = square(spb15); 
complex<T> t24 = square(spb25); 
complex<T> t25 = square(spa13); 
complex<T> t26 = s14 - s23 + s45; 
complex<T> t29 = -(spa14*spa34); 
complex<T> t31 = s14 - s35; 
complex<T> t34 = s14*s45; 
complex<T> t39 = spa12*T(2); 
complex<T> t40 = cube(spb15); 
complex<T> t41 = cube(spb25); 
complex<T> t43 = spa13*spa34; 
complex<T> t45 = s23*spa24; 
complex<T> t49 = spb12*spb23; 
complex<T> t54 = spa13*spb15; 
complex<T> t59 = spa14*spa45; 
complex<T> t60 = spa23*spa24; 
complex<T> t66 = s45*spa13; 
complex<T> t76 = spa23*spa45; 
complex<T> t105 = spa34*spb25; 
complex<T> d16 = s14 - s23; d16 = T(1)/d16;
complex<T> d26 = spa15*spa23*spa35; d26 = T(1)/d26;
complex<T> d27 = (s12 + s15 - s34)*spa15*spa23; d27 = T(1)/d27;
complex<T> d34 = (s12 + s15 - s34)*spa23; d34 = T(1)/d34;
complex<T> d43 = (s12 + s15 - s34)*spa23*T(2); d43 = T(1)/d43;
complex<T> d46 = spa35*spa45*T(2); d46 = T(1)/d46;
complex<T> t27 = -t76; 
complex<T> t28 = -s23 + t31; 
complex<T> t37 = -t49; 
complex<T> t42 = -t59; 
complex<T> t44 = t11*T(2); 
complex<T> t57 = t25*t40; 
complex<T> t68 = spa34*t24; 
complex<T> t83 = t23*t29; 
complex<T> t88 = spa14*t43; 
complex<T> t89 = spa45*t41; 
complex<T> t91 = t23*t34; 
complex<T> t92 = s35*t45; 
complex<T> t94 = spb25*t10; 
complex<T> t100 = spa14*t10; 
complex<T> t112 = s34*t10; 
complex<T> d3 = spa23*t26*t5; d3 = T(1)/d3;
complex<T> d4 = spa23*t5*square(t26); d4 = T(1)/d4;
complex<T> d6 = t39*square(t31); d6 = T(1)/d6;
complex<T> d13 = spa35*t5*t76*T(2); d13 = T(1)/d13;
complex<T> d14 = spa45*t11; d14 = T(1)/d14;
complex<T> d15 = spa25*spa35*t39; d15 = T(1)/d15;
complex<T> d17 = spa25*t39*square(s14 - s23); d17 = T(1)/d17;
complex<T> d18 = spa35*spa45*t11; d18 = T(1)/d18;
complex<T> d19 = spa12*spa15*t6; d19 = T(1)/d19;
complex<T> d20 = spa12*t26*t6; d20 = T(1)/d20;
complex<T> d21 = t11*t8*square(s23 - s45); d21 = T(1)/d21;
complex<T> d22 = spa12*t6*square(t26); d22 = T(1)/d22;
complex<T> d24 = spa12*spa45*t6*t8; d24 = T(1)/d24;
complex<T> d25 = spa45*t11*t8; d25 = T(1)/d25;
complex<T> d28 = spa35*t11*t26; d28 = T(1)/d28;
complex<T> d29 = t11*square(t26); d29 = T(1)/d29;
complex<T> d30 = t11*cube(t26); d30 = T(1)/d30;
complex<T> d35 = spa12*spa35*t26; d35 = T(1)/d35;
complex<T> d36 = spa12*square(t26); d36 = T(1)/d36;
complex<T> d37 = spa12*cube(t26); d37 = T(1)/d37;
complex<T> d38 = spa15*spa35*t11; d38 = T(1)/d38;
complex<T> d40 = spa15*t11; d40 = T(1)/d40;
complex<T> d41 = spa35*t11; d41 = T(1)/d41;
complex<T> d45 = spa45*t8; d45 = T(1)/d45;
complex<T> d51 = spa35*t11*t8; d51 = T(1)/d51;
complex<T> t17 = -(d24*(spa14*spa23*spb12 + spa45*t54)); 
complex<T> t19 = d13*(-(spa14*spa35*spb15) + spb25*t76); 
complex<T> t50 = d17*spb35; 
complex<T> t52 = d45*t100*t49 + d46*t37*cube(spa34); 
complex<T> t61 = -(d14*T(2)); 
complex<T> t69 = t57*t59; 
complex<T> t71 = d14*T(2); 
complex<T> t75 = -t94; 
complex<T> t81 = d26*spa13*spb12*t10 + d27*s12*t94; 
complex<T> t82 = t27*t41; 
complex<T> t86 = d21*t25; 
complex<T> d1 = spa35*spa45*t44; d1 = T(1)/d1;
complex<T> d2 = t44*square(s14 - s23); d2 = T(1)/d2;
complex<T> d5 = t26*t44*square(s14 - s23); d5 = T(1)/d5;
complex<T> d7 = t28*t5; d7 = T(1)/d7;
complex<T> d8 = spa12*t28*t31; d8 = T(1)/d8;
complex<T> d9 = t5*square(t28); d9 = T(1)/d9;
complex<T> d10 = spa12*t31*square(t28); d10 = T(1)/d10;
complex<T> d11 = t28*t39*square(s14 - s23); d11 = T(1)/d11;
complex<T> d12 = t28*t39*square(t31); d12 = T(1)/d12;
complex<T> d23 = t26*t44*square(s23 - s45); d23 = T(1)/d23;
complex<T> d31 = spa12*spa35*t28; d31 = T(1)/d31;
complex<T> d32 = spa12*square(t28); d32 = T(1)/d32;
complex<T> d33 = spa12*cube(t28); d33 = T(1)/d33;
complex<T> d39 = spa12*t28; d39 = T(1)/d39;
complex<T> d42 = t44*square(t26); d42 = T(1)/d42;
complex<T> d44 = t39*square(t28); d44 = T(1)/d44;
complex<T> d47 = t39*cube(t28); d47 = T(1)/d47;
complex<T> d48 = t28*t39; d48 = T(1)/d48;
complex<T> d49 = spa35*t26*t44; d49 = T(1)/d49;
complex<T> d50 = t44*cube(t26); d50 = T(1)/d50;
complex<T> t1 = d15*d16*t105*t27 + t105*t27*t50 + d22*t42*t57 + d7*spa24*t68 + d23*t69 + d4*t69 + d5*t69 + d16*spb15*t100*t71 + d3*spa13*t83 + t23*t42*t86 + d2*t23*t88 + d20*t23*t88 + d11*t60*t89 + d9*t60*t89 + d19*spb15*t88*T(2) - d18*cube(spa34)*T(2) + spa14*spa34*t17*T(3) + t10*t19*T(3); 
complex<T> t2 = d4*t42*t57 + d5*t42*t57 + d16*spb15*t100*t61 + d6*spa24*t68 - d7*spa24*t68 - d8*spa24*t68 + d15*d16*t105*t76 + t105*t50*t76 + d10*spa24*t82 + d11*spa24*t82 + d12*spa24*t82 + d9*spa24*t82 + d2*spa13*t83 + d3*t23*t88 + d1*cube(spa34) - t10*t19*T(3); 
complex<T> t9 = spa34*t16*t17 + d23*t42*t57 + d22*t69 + d20*spa13*t83 + t23*t59*t86 - d19*spb15*t88*T(2) + d25*t100*T(3); 
complex<T> t18 = d50*t34*t69 + d29*t88*t91 + d49*t34*t54*square(spa34); 
complex<T> t20 = -(d6*spa24*t68) + d8*spa24*t68 + (d10 + d12)*t60*t89; 
complex<T> t21 = -(d32*t45*t68) + d37*spb23*t69 + d33*t45*t82 + d36*spa13*spb23*t83 - d31*s23*spb25*square(spa34) + d35*spb23*t54*square(spa34); 
complex<T> t22 = d30*s45*t69 + d29*t66*t83 - d41*spb45*cube(spa34) + d40*spa14*spb45*square(spa34) + d28*s45*t54*square(spa34) + d38*t66*square(spa34); 
complex<T> t48 = d32*s35; 
complex<T> t56 = s14*(d28*spa13*spb15*t10 + d32*spa24*t68 + d30*t69 + d29*spa13*t83 + d33*t60*t89 + d31*t94); 
complex<T> t93 = d38*spa13*t112 + d27*s34*t75; 
complex<T> t35 = -t48; 
complex<T> t74 = t45*t48*t68 + d47*t82*t92 + d48*s23*spb35*t94; 
complex<T> t67 = spa24*t35*t68 + d33*s35*spa24*t82 + d39*spb35*t94; 
complex<T> co1 = d34*spb15*t75; 
complex<T> co2 = d42*t16*t43*t91; 
complex<T> co3 = d43*s12*spb15*t75; 
complex<T> co4 = -(d44*t68*t92*T(3)); 
complex<T> co5 = d51*t112*t66; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t20*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t81*(*CI_users[4]->get_value(mc,ind,mu)) + t56*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t21*(*CI_users[7]->get_value(mc,ind,mu)) + t93*(*CI_users[8]->get_value(mc,ind,mu)) + t67*(*CI_users[9]->get_value(mc,ind,mu)) + t22*(*CI_users[10]->get_value(mc,ind,mu)) + co2*(*CI_users[12]->get_value(mc,ind,mu)) + co3*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[15]->get_value(mc,ind,mu)) + t52*(*CI_users[16]->get_value(mc,ind,mu)) + t74*(*CI_users[22]->get_value(mc,ind,mu)) + t18*(*CI_users[23]->get_value(mc,ind,mu)) + co5*(*CI_users[24]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qpmqmpgap_SLC_wCI::\
C2q2g1y_qpmqmpgap_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c35));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c2, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpmqmpgap_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, qm, p, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpmqmpgap SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:41:02 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> s12 = S(1,2);
complex<T> s23 = S(2,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t12 = square(spa23); 
complex<T> t13 = square(spb45); 
complex<T> t16 = square(spa24); 
complex<T> t17 = square(spa25); 
complex<T> t28 = spa24*spb45; 
complex<T> d1 = (s12 - s35)*spa14*spa45; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa15*spa45; d2 = T(1)/d2;
complex<T> d3 = spa15*spa45*square(s12 - s34)*T(2); d3 = T(1)/d3;
complex<T> d4 = spa14*spa45*square(s12 - s35)*T(2); d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*spa15*spa34*spa45*T(2); d5 = T(1)/d5;
complex<T> d6 = (s12 - s35)*spa14*spa35*spa45*T(2); d6 = T(1)/d6;
complex<T> d7 = spa15*spa34*spa45*T(2); d7 = T(1)/d7;
complex<T> d8 = spa14*spa35*spa45*T(2); d8 = T(1)/d8;
complex<T> d9 = spa35*spa45*T(2); d9 = T(1)/d9;
complex<T> d10 = spa34*spa45*T(2); d10 = T(1)/d10;
complex<T> d11 = spa15*spa45*T(2); d11 = T(1)/d11;
complex<T> d12 = spa14*spa45*T(2); d12 = T(1)/d12;
complex<T> t6 = d4*spa35*t13*t16 + d1*spa23*t28*T(2) + d8*t12*T(3) + d6*spa23*(spa12*spa34*spb14 - spa35*t28)*T(3); 
complex<T> t24 = d2*spa25; 
complex<T> t4 = -(d3*spa34*t13*t17) + spa23*spb45*t24*T(2) - d5*spa23*(spa12*spa35*spb15 + spa25*spa34*spb45)*T(3) - d7*t12*T(3); 
complex<T> t5 = -(d4*spa35*t13*t16) + d3*spa34*t13*t17 - spa23*spb45*t24*T(2) - d1*spa23*t28*T(2) + d5*spa23*(spa12*spa35*spb15 + spa25*spa34*spb45)*T(3) + d6*(-(spa12*spa23*spa34*spb14*T(3)) + spa23*spa35*t28*T(3)); 
complex<T> co1 = -(d9*s12*spb14*t12); 
complex<T> co2 = d10*s12*spb15*t12; 
complex<T> co3 = d11*s23*spb34*t12; 
complex<T> co4 = -(d12*s23*spb35*t12); 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qpmqmmgap_SLC_wCI::\
C2q2g1y_qpmqmmgap_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c5, c124));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c35));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c3, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpmqmmgap_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, qm, m, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpmqmmgap SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:42:04 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa45 = SPA(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s14 = -(spa14*spb14);
complex<T> s35 = S(3,5);
complex<T> s15 = S(1,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = -(spa34*spb34);
complex<T> t4 = spb12*spb23; 
complex<T> t6 = spb34*T(2); 
complex<T> t10 = square(spb15); 
complex<T> t16 = -(spb35*T(3)); 
complex<T> t23 = square(spa24); 
complex<T> t24 = square(spa34); 
complex<T> t25 = square(spb13); 
complex<T> t28 = -(spb15*spb35); 
complex<T> t30 = s12 - s45; 
complex<T> t31 = s14 - s35; 
complex<T> t33 = s35*s45; 
complex<T> t40 = cube(spa24); 
complex<T> t41 = cube(spa34); 
complex<T> t43 = spb13*spb15; 
complex<T> t44 = spb23*T(2); 
complex<T> t46 = s12*spb25; 
complex<T> t49 = spa12*spa23; 
complex<T> t55 = s45*spb13; 
complex<T> t58 = spb35*spb45; 
complex<T> t68 = spb12*spb45; 
complex<T> t70 = -(spa34*T(2)); 
complex<T> t77 = spa34*T(2); 
complex<T> t89 = spa24*spb15; 
complex<T> d19 = s12 - s35; d19 = T(1)/d19;
complex<T> d28 = spb14*spb23*spb34; d28 = T(1)/d28;
complex<T> d29 = spb23*square(spb34); d29 = T(1)/d29;
complex<T> d31 = spb23*cube(spb34); d31 = T(1)/d31;
complex<T> d34 = spb12*spb24*spb34; d34 = T(1)/d34;
complex<T> d35 = (s15 - s23 + s45)*spb12*spb34; d35 = T(1)/d35;
complex<T> d36 = spb12*spb24; d36 = T(1)/d36;
complex<T> d42 = spb14*spb45*T(2); d42 = T(1)/d42;
complex<T> d45 = spb12*spb24*T(2); d45 = T(1)/d45;
complex<T> t11 = t4*T(2); 
complex<T> t27 = -t68; 
complex<T> t37 = -t49; 
complex<T> t42 = -t58; 
complex<T> t59 = t25*t41; 
complex<T> t67 = spb25*t40; 
complex<T> t69 = spb35*t43; 
complex<T> t75 = spb15*t23; 
complex<T> t76 = t25*t58; 
complex<T> t84 = -(spa14*t10); 
complex<T> t86 = spb13*t24; 
complex<T> t100 = s15*t10; 
complex<T> t104 = s14*t46; 
complex<T> t108 = s23*t10; 
complex<T> d1 = (s12 - s35)*spb23*(s12 + t31); d1 = T(1)/d1;
complex<T> d3 = s34*(s12 - s35)*t4; d3 = T(1)/d3;
complex<T> d4 = s34*t30*t4; d4 = T(1)/d4;
complex<T> d5 = spb34*t30*t4; d5 = T(1)/d5;
complex<T> d6 = spb14*spb45*t4; d6 = T(1)/d6;
complex<T> d7 = spb24*t44*square(s12 - s35); d7 = T(1)/d7;
complex<T> d8 = (s12 - s35)*spb23*square(s12 + t31); d8 = T(1)/d8;
complex<T> d9 = (s12 + t31)*t44*square(s12 - s35); d9 = T(1)/d9;
complex<T> d11 = (s12 - s35)*t4*square(s34); d11 = T(1)/d11;
complex<T> d13 = t30*t4*square(s34); d13 = T(1)/d13;
complex<T> d14 = t4*t6*square(t30); d14 = T(1)/d14;
complex<T> d16 = spb45*t30*t4*t6; d16 = T(1)/d16;
complex<T> d17 = spb45*t4; d17 = T(1)/d17;
complex<T> d18 = spb14*spb24*t44; d18 = T(1)/d18;
complex<T> d20 = t44*square(t31); d20 = T(1)/d20;
complex<T> d21 = spb23*t31*(s12 + t31); d21 = T(1)/d21;
complex<T> d22 = spb23*t31*square(s12 + t31); d22 = T(1)/d22;
complex<T> d23 = (s12 + t31)*t44*square(t31); d23 = T(1)/d23;
complex<T> d25 = spb45*t4*t6; d25 = T(1)/d25;
complex<T> d26 = spb14*spb23*(s12 + t31); d26 = T(1)/d26;
complex<T> d27 = spb23*square(s12 + t31); d27 = T(1)/d27;
complex<T> d30 = spb23*cube(s12 + t31); d30 = T(1)/d30;
complex<T> d32 = spb23*(s12 + t31); d32 = T(1)/d32;
complex<T> d33 = (s15 - s23 + s45)*spb34*t4; d33 = T(1)/d33;
complex<T> d37 = spb14*spb34*t4; d37 = T(1)/d37;
complex<T> d38 = t4*square(spb34); d38 = T(1)/d38;
complex<T> d39 = t4*cube(spb34); d39 = T(1)/d39;
complex<T> d40 = spb14*t4; d40 = T(1)/d40;
complex<T> d41 = spb34*t4; d41 = T(1)/d41;
complex<T> d43 = spb45*t6; d43 = T(1)/d43;
complex<T> d44 = t44*square(s12 + t31); d44 = T(1)/d44;
complex<T> d47 = (s12 + t31)*t44; d47 = T(1)/d47;
complex<T> d48 = t44*cube(s12 + t31); d48 = T(1)/d48;
complex<T> d49 = (s15 - s23 + s45)*t4*t6; d49 = T(1)/d49;
complex<T> d50 = spb14*t4*t6; d50 = T(1)/d50;
complex<T> t18 = d16*(spa23*spb12*spb35 + spa34*spb13*spb45); 
complex<T> t19 = d22*t27*t67 + d23*t27*t67 + (-d20 + d21)*spb25*t75; 
complex<T> t20 = s35*(d38*spb13*t28 + d30*t27*t67 + d27*spb25*t75 + d39*t76 - d26*spa24*square(spb15) + d37*spb13*square(spb15)); 
complex<T> t36 = d33*spa14; 
complex<T> t47 = d27*s14; 
complex<T> t48 = d7*spa14; 
complex<T> t51 = d30*spb12; 
complex<T> t53 = d43*spb35*t10*t49 + d42*t37*cube(spb15); 
complex<T> t80 = d35*spa14*spa23*spb13*t10 + d34*t108; 
complex<T> t101 = d17*spb35; 
complex<T> d2 = t11*square(s12 - s35); d2 = T(1)/d2;
complex<T> d10 = s34*t11*square(s12 - s35); d10 = T(1)/d10;
complex<T> d12 = s34*t11*square(t30); d12 = T(1)/d12;
complex<T> d15 = (s12 - s35)*spb14*spb45*t11; d15 = T(1)/d15;
complex<T> d24 = spb14*spb45*t11; d24 = T(1)/d24;
complex<T> d46 = t11*square(spb34); d46 = T(1)/d46;
complex<T> d51 = t11*cube(spb34); d51 = T(1)/d51;
complex<T> t9 = d12*t58*t59 + d13*t58*t59 + d4*t24*t69 + d5*t69*t70 + d14*t24*t76 + d25*spb35*t10*T(3) + spb15*spb35*t18*T(3); 
complex<T> t17 = d15*(spa34*spb14*spb35 + spa24*t27); 
complex<T> t21 = d38*t28*t55 + d39*s45*t76 - d40*spa45*cube(spb15) + d41*spa45*spb35*square(spb15) + d37*t55*square(spb15) + t36*t55*square(spb15); 
complex<T> t22 = d29*spa12*spb13*t28 + spb45*t40*t46*t51 - d27*t46*t75 + d31*spa12*t76 + d26*s12*spa24*square(spb15) + d28*spa12*spb13*square(spb15); 
complex<T> t35 = -t47; 
complex<T> t81 = d48*t104*t40*t68 + t46*t47*t75 + d47*s12*spa24*t84; 
complex<T> t82 = t33*(d50*spb13*t10 + d38*t69 + d51*t76); 
complex<T> t88 = -(d34*t100) + spb13*t100*t36; 
complex<T> t1 = d10*t58*t59 + d11*t58*t59 + d22*t67*t68 + d23*t67*t68 + d8*t67*t68 + d9*t67*t68 + d3*t24*t69 - d1*spb25*t75 + d20*spb25*t75 - d21*spb25*t75 + d19*t10*t101*t77 + d2*t28*t86 + d18*d19*t27*t89 + t48*t68*t89 + d24*cube(spb15) - t10*t17*T(3); 
complex<T> t2 = spb15*t16*t18 + d14*t24*t25*t42 + d10*t42*t59 + d11*t42*t59 + d12*t42*t59 + d13*t42*t59 + d8*t27*t67 + d9*t27*t67 + d2*t24*t69 + d19*t10*t101*t70 + d1*spb25*t75 + d5*t69*t77 + d3*t28*t86 + d4*t28*t86 + t27*t48*t89 + d18*d19*t68*t89 - d6*cube(spb15)*T(2) + t10*t17*T(3); 
complex<T> t74 = s14*spb45*t51*t67 + spb25*t35*t75 + d32*spa24*t84; 
complex<T> co1 = -(d36*spa34*t10); 
complex<T> co2 = -(d44*t104*t75*T(3)); 
complex<T> co3 = -(d45*spa34*t108); 
complex<T> co4 = d46*t16*t33*t43; 
complex<T> co5 = d49*spa14*t100*t55; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t19*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t22*(*CI_users[4]->get_value(mc,ind,mu)) + t74*(*CI_users[5]->get_value(mc,ind,mu)) + t88*(*CI_users[6]->get_value(mc,ind,mu)) + t80*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t20*(*CI_users[9]->get_value(mc,ind,mu)) + t21*(*CI_users[10]->get_value(mc,ind,mu)) + t53*(*CI_users[11]->get_value(mc,ind,mu)) + co2*(*CI_users[12]->get_value(mc,ind,mu)) + co3*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + t81*(*CI_users[17]->get_value(mc,ind,mu)) + co5*(*CI_users[18]->get_value(mc,ind,mu)) + t82*(*CI_users[19]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g1y_qpqmppgap_SLC_wCI::\
C2q2g1y_qpqmppgap_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpqmppgap_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpqmppgap SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1y_qpqmpmgap_SLC_wCI::\
C2q2g1y_qpqmpmgap_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpqmpmgap_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, m, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpqmpmgap SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 21:55:43 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = S(1,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s23 = -(spa23*spb23);
complex<T> s13 = -(spa13*spb13);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s35 = -(spa35*spb35);
complex<T> t6 = spa23*T(2); 
complex<T> t10 = square(s12 - s45); 
complex<T> t12 = square(s12 + s15 - s34); 
complex<T> t18 = spa15*spa34; 
complex<T> t25 = square(spa24); 
complex<T> t26 = square(spb13); 
complex<T> t27 = square(spb35); 
complex<T> t28 = spa25*spa34; 
complex<T> t29 = spa12*spa23; 
complex<T> t30 = -(s12*spa14); 
complex<T> t31 = spa35*T(2); 
complex<T> t32 = spa15*spa45; 
complex<T> t45 = square(spb15); 
complex<T> t47 = cube(spb13); 
complex<T> t50 = s23 - s45; 
complex<T> t51 = square(spa45); 
complex<T> t55 = cube(spb35); 
complex<T> t56 = s25 - s34; 
complex<T> t61 = s14 - s23; 
complex<T> t62 = s14 - s25; 
complex<T> t69 = s14*s34; 
complex<T> t72 = square(s13); 
complex<T> t73 = square(spa14); 
complex<T> t74 = cube(spb15); 
complex<T> t75 = cube(spa24); 
complex<T> t76 = spa12*spa34; 
complex<T> t79 = -(spa25*T(3)); 
complex<T> t81 = s12*s23; 
complex<T> t82 = spa15*T(2); 
complex<T> t83 = s25*spb25; 
complex<T> t85 = s14*s45; 
complex<T> t92 = spa14*spa24; 
complex<T> t93 = spa34*spa35; 
complex<T> t99 = square(spa34); 
complex<T> t100 = -(spa45*spb15); 
complex<T> t107 = spa12*spa14; 
complex<T> t108 = spa25*spa45; 
complex<T> t111 = spa23*spa35; 
complex<T> t166 = spa25*T(3); 
complex<T> t186 = spa13*spa25; 
complex<T> t246 = spa14*spa23; 
complex<T> d7 = s13*(s12 - s45)*spa35*spa45; d7 = T(1)/d7;
complex<T> d17 = (s12 - s45)*spa15*spa23; d17 = T(1)/d17;
complex<T> d18 = (s12 - s34)*spa12*spa35; d18 = T(1)/d18;
complex<T> d19 = (s12 - s45)*spa12*spa35; d19 = T(1)/d19;
complex<T> d20 = (s12 - s34)*s35*spa12; d20 = T(1)/d20;
complex<T> d21 = s35*(s12 - s45)*spa12; d21 = T(1)/d21;
complex<T> d22 = (s12 - s34)*s35*spa12*spa15; d22 = T(1)/d22;
complex<T> d23 = s35*(s12 - s45)*spa12*spa15; d23 = T(1)/d23;
complex<T> d29 = (s12 - s34)*spa15*square(s35); d29 = T(1)/d29;
complex<T> d32 = (s12 - s45)*spa15*square(s35); d32 = T(1)/d32;
complex<T> d34 = spa12*spa35; d34 = T(1)/d34;
complex<T> d35 = s12 - s45; d35 = T(1)/d35;
complex<T> d37 = s12 - s34; d37 = T(1)/d37;
complex<T> d59 = spa12*spa35*spa45; d59 = T(1)/d59;
complex<T> d81 = spa35*spa45*square(spa13); d81 = T(1)/d81;
complex<T> d83 = spa35*spa45*cube(spa13); d83 = T(1)/d83;
complex<T> d84 = spa15*cube(spa35); d84 = T(1)/d84;
complex<T> d85 = square(spa35); d85 = T(1)/d85;
complex<T> d86 = spa15*square(spa35); d86 = T(1)/d86;
complex<T> d87 = spa15*spa23*square(spa35); d87 = T(1)/d87;
complex<T> d100 = spa13*spa35*spa45; d100 = T(1)/d100;
complex<T> d108 = spa12*square(spa35); d108 = T(1)/d108;
complex<T> d109 = spa12*spa15*square(spa35); d109 = T(1)/d109;
complex<T> d115 = spa15*square(spa13); d115 = T(1)/d115;
complex<T> d116 = spa15*spa23*square(spa13); d116 = T(1)/d116;
complex<T> d117 = spa15*cube(spa13); d117 = T(1)/d117;
complex<T> d118 = spa35*square(spa13); d118 = T(1)/d118;
complex<T> d120 = spa35*cube(spa13); d120 = T(1)/d120;
complex<T> d129 = spa12*square(spa35)*T(2); d129 = T(1)/d129;
complex<T> t14 = spa25*t56; 
complex<T> t16 = square(t62); 
complex<T> t46 = -t92; 
complex<T> t48 = -t76; 
complex<T> t52 = s12 + t56; 
complex<T> t53 = s34 + t62; 
complex<T> t58 = spa13*spa24 + t246*T(2); 
complex<T> t59 = -t83; 
complex<T> t60 = spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6; 
complex<T> t68 = -(spa14*spa25*spb15) - spb13*t76; 
complex<T> t78 = spa23*t26; 
complex<T> t84 = spa35*t50; 
complex<T> t94 = t26*(spa13*spa24 + spa14*t6); 
complex<T> t98 = -t107; 
complex<T> t112 = -(d34*T(2)); 
complex<T> t127 = t107*t74; 
complex<T> t129 = -(spa45*t28); 
complex<T> t130 = d34*T(2); 
complex<T> t146 = t47*t99; 
complex<T> t147 = spa45*t55; 
complex<T> t154 = spa25*t45; 
complex<T> t155 = spa24*t73; 
complex<T> t156 = spa34*t47; 
complex<T> t162 = spa14*t29; 
complex<T> t163 = spb15*t25; 
complex<T> t164 = spa34*t26; 
complex<T> t165 = spa45*t27; 
complex<T> t176 = spa45*t45; 
complex<T> t177 = t47*t73; 
complex<T> t185 = spb13*t25; 
complex<T> t210 = t79*t92; 
complex<T> t224 = t75*t83; 
complex<T> t239 = t107*t99; 
complex<T> t261 = t47*t76; 
complex<T> d1 = spa12*spa35*t28; d1 = T(1)/d1;
complex<T> d2 = spa12*spa45*t18*t6; d2 = T(1)/d2;
complex<T> d3 = spa35*spa45*t29; d3 = T(1)/d3;
complex<T> d4 = s13*(s12 - s45)*t32; d4 = T(1)/d4;
complex<T> d5 = s13*(s12 - s45)*spa23*t32; d5 = T(1)/d5;
complex<T> d6 = t10*t32*t6; d6 = T(1)/d6;
complex<T> d8 = s13*t10*t32*T(2); d8 = T(1)/d8;
complex<T> d9 = (s12 - s45)*t32*t72; d9 = T(1)/d9;
complex<T> d10 = s13*spa45*t10*t31; d10 = T(1)/d10;
complex<T> d11 = (s12 - s45)*spa35*spa45*t72; d11 = T(1)/d11;
complex<T> d12 = (s12 - s34)*spa23*t18; d12 = T(1)/d12;
complex<T> d16 = (s12 - s34)*(s12 + s15 - s34)*spa23*t18; d16 = T(1)/d16;
complex<T> d24 = (s12 - s34)*s35*spa15*t29; d24 = T(1)/d24;
complex<T> d25 = s35*(s12 - s45)*spa15*t29; d25 = T(1)/d25;
complex<T> d26 = spa12*t31*square(s12 - s34); d26 = T(1)/d26;
complex<T> d27 = t10*t82; d27 = T(1)/d27;
complex<T> d28 = spa12*t10*t31; d28 = T(1)/d28;
complex<T> d30 = s35*t82*square(s12 - s34); d30 = T(1)/d30;
complex<T> d31 = s35*t10*t82; d31 = T(1)/d31;
complex<T> d33 = spa45*t111; d33 = T(1)/d33;
complex<T> d36 = spa35*t28; d36 = T(1)/d36;
complex<T> d42 = spa35*t6*square(t61); d42 = T(1)/d42;
complex<T> d43 = t111*t61*(s45 + t61); d43 = T(1)/d43;
complex<T> d44 = t111*t61*square(s45 + t61); d44 = T(1)/d44;
complex<T> d45 = spa35*t6*(s45 + t61)*square(t61); d45 = T(1)/d45;
complex<T> d46 = (s15 - s34)*spa23*t18; d46 = T(1)/d46;
complex<T> d47 = (s15 - s34)*(s12 + s15 - s34)*spa23*t18; d47 = T(1)/d47;
complex<T> d48 = spa12*spa35*spa45*t6; d48 = T(1)/d48;
complex<T> d49 = spa23*t32*t50; d49 = T(1)/d49;
complex<T> d50 = t32*square(t50)*T(2); d50 = T(1)/d50;
complex<T> d51 = s13*t32*t50; d51 = T(1)/d51;
complex<T> d52 = s13*spa23*t32*t50; d52 = T(1)/d52;
complex<T> d54 = s13*t32*square(t50)*T(2); d54 = T(1)/d54;
complex<T> d55 = t32*t50*t72; d55 = T(1)/d55;
complex<T> d56 = s13*spa45*t31*square(t50); d56 = T(1)/d56;
complex<T> d58 = spa45*t31; d58 = T(1)/d58;
complex<T> d60 = square(t50); d60 = T(1)/d60;
complex<T> d62 = t111*t50*(s45 + t61); d62 = T(1)/d62;
complex<T> d63 = t111*t50*square(s45 + t61); d63 = T(1)/d63;
complex<T> d64 = spa35*t6*(s45 + t61)*square(t50); d64 = T(1)/d64;
complex<T> d66 = spa12*t28*t31; d66 = T(1)/d66;
complex<T> d75 = spa34*t31; d75 = T(1)/d75;
complex<T> d76 = spa35*t76; d76 = T(1)/d76;
complex<T> d77 = square(t56); d77 = T(1)/d77;
complex<T> d78 = t32*square(spa13); d78 = T(1)/d78;
complex<T> d79 = spa23*t32*square(spa13); d79 = T(1)/d79;
complex<T> d80 = t32*cube(spa13); d80 = T(1)/d80;
complex<T> d82 = spa13*spa45*t111; d82 = T(1)/d82;
complex<T> d88 = spa23*spa45*t18; d88 = T(1)/d88;
complex<T> d92 = spa23*t12*t18; d92 = T(1)/d92;
complex<T> d96 = t111*(s45 + t61); d96 = T(1)/d96;
complex<T> d97 = t111*square(s45 + t61); d97 = T(1)/d97;
complex<T> d98 = t111*cube(s45 + t61); d98 = T(1)/d98;
complex<T> d99 = spa23*spa34*t12; d99 = T(1)/d99;
complex<T> d101 = spa35*(s45 + t61); d101 = T(1)/d101;
complex<T> d102 = spa35*square(s45 + t61); d102 = T(1)/d102;
complex<T> d103 = spa35*cube(s45 + t61); d103 = T(1)/d103;
complex<T> d110 = spa15*t29*square(spa35); d110 = T(1)/d110;
complex<T> d114 = spa15*spa23*t12; d114 = T(1)/d114;
complex<T> d119 = spa13*t111; d119 = T(1)/d119;
complex<T> d121 = t32*square(spa13)*T(2); d121 = T(1)/d121;
complex<T> d122 = t32*cube(spa13)*T(2); d122 = T(1)/d122;
complex<T> d123 = spa45*t31*square(spa13); d123 = T(1)/d123;
complex<T> d124 = spa45*t31*cube(spa13); d124 = T(1)/d124;
complex<T> d125 = spa13*spa45*t31; d125 = T(1)/d125;
complex<T> d130 = spa12*t82*square(spa35); d130 = T(1)/d130;
complex<T> d131 = spa12*spa15*t6*square(spa35); d131 = T(1)/d131;
complex<T> d132 = t82*cube(spa35); d132 = T(1)/d132;
complex<T> d133 = spa34*t12*t6; d133 = T(1)/d133;
complex<T> d137 = spa35*t6*(s45 + t61); d137 = T(1)/d137;
complex<T> d138 = spa35*t6*square(s45 + t61); d138 = T(1)/d138;
complex<T> d139 = spa35*t6*cube(s45 + t61); d139 = T(1)/d139;
complex<T> t8 = square(t53); 
complex<T> t19 = s34*s45*(d132*t129 + d130*t210 - d129*t25 + d131*spa24*(spa15*spa23*spa24 + spa24*t186 + spa25*t246*T(2))); 
complex<T> t34 = cube(t53); 
complex<T> t39 = d110*(spa15*spa23*spa24 + spa24*t186 + spa25*t246*T(2)); 
complex<T> t80 = -t147; 
complex<T> t109 = -t163; 
complex<T> t110 = -t156; 
complex<T> t115 = d58*spb12; 
complex<T> t133 = d75*spb12; 
complex<T> t134 = -(d33*spb13); 
complex<T> t135 = d46*spb25; 
complex<T> t144 = -t155; 
complex<T> t153 = d125*s12*spa34*spb23*t25 + d124*t239*t81 + d123*spa34*t81*t92; 
complex<T> t159 = -(d78*T(3)); 
complex<T> t169 = t108*t127; 
complex<T> t172 = spb35*t112; 
complex<T> t180 = spb35*t130; 
complex<T> t184 = t127*t51; 
complex<T> t193 = t74*t98; 
complex<T> t200 = t154*t46; 
complex<T> t231 = t48*t73; 
complex<T> t233 = d6*t26; 
complex<T> d13 = (s12 - s34)*t52*t93; d13 = T(1)/d13;
complex<T> d14 = (s12 - s34)*t93*square(t52); d14 = T(1)/d14;
complex<T> d15 = spa34*t31*t52*square(s12 - s34); d15 = T(1)/d15;
complex<T> d38 = spa25*t16*t31; d38 = T(1)/d38;
complex<T> d39 = spa25*spa35*t53*t62; d39 = T(1)/d39;
complex<T> d41 = spa25*t16*t31*t53; d41 = T(1)/d41;
complex<T> d53 = s13*spa45*t84; d53 = T(1)/d53;
complex<T> d57 = spa45*t72*t84; d57 = T(1)/d57;
complex<T> d61 = t6*t84; d61 = T(1)/d61;
complex<T> d65 = spa12*spa45*t6*t84; d65 = T(1)/d65;
complex<T> d67 = t14*t31; d67 = T(1)/d67;
complex<T> d68 = spa35*t14*t53; d68 = T(1)/d68;
complex<T> d70 = spa25*t31*t53*square(t56); d70 = T(1)/d70;
complex<T> d71 = t52*t56*t93; d71 = T(1)/d71;
complex<T> d72 = t56*t93*square(t52); d72 = T(1)/d72;
complex<T> d73 = spa34*t31*t52*square(t56); d73 = T(1)/d73;
complex<T> d74 = t14*t31*t76; d74 = T(1)/d74;
complex<T> d89 = spa35*t28*t52; d89 = T(1)/d89;
complex<T> d90 = t93*square(t52); d90 = T(1)/d90;
complex<T> d91 = t93*cube(t52); d91 = T(1)/d91;
complex<T> d93 = spa25*spa35*t53; d93 = T(1)/d93;
complex<T> d104 = spa35*t53; d104 = T(1)/d104;
complex<T> d107 = t52*t93; d107 = T(1)/d107;
complex<T> d111 = spa25*spa35*t52; d111 = T(1)/d111;
complex<T> d112 = spa35*square(t52); d112 = T(1)/d112;
complex<T> d113 = spa35*cube(t52); d113 = T(1)/d113;
complex<T> d126 = spa25*t31*t53; d126 = T(1)/d126;
complex<T> d134 = spa34*t31*square(t52); d134 = T(1)/d134;
complex<T> d135 = spa34*t31*cube(t52); d135 = T(1)/d135;
complex<T> d136 = spa34*t31*t52; d136 = T(1)/d136;
complex<T> t5 = d12*spa14*t109 + d10*t107*t146 + d27*spa24*t165 + d28*spa23*spa24*t165 + d16*t224 + d17*spb35*t25 - d18*spb35*t25 - d19*spb35*t25 + d37*(d36*t100 + t180)*t25 + d35*(d33*spa34*t185 + t180*t25) + d22*t210*t27 + d23*t210*t27 - d20*t25*t27 - d21*t25*t27 + d29*t147*t28 + d30*t147*t28 + d31*t147*t28 + d32*t147*t28 + d26*spa24*t27*t28 + d8*t177*t48 + d14*t193*t51 + d15*t193*t51 + d24*spa24*t27*(spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6) + d25*spa24*t27*(spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6) - d1*t75 - d3*t75 + d9*t177*t76 + t233*t46*t76 + d7*t164*t92 + d13*t176*t92 + d5*t26*t58*t92 + d11*t146*t98 - d4*t155*t26*T(3) + d2*spa14*t75*T(3); 
complex<T> t20 = d121*spa24*spb23*t30*t58 + d122*t231*t81 - d121*t155*t81*T(3); 
complex<T> t21 = spb23*(d101*t163 + d103*t169 + d102*t200 + d100*spa34*t25 + d78*t46*t58) + s23*(t155*t159 + d80*t231 + d83*t239 + d81*spa34*t92); 
complex<T> t22 = d86*spb12*t210 - d85*spb12*t25 + d90*spa24*t176*t30 - d33*spb12*t75 - d36*spb12*t75 + d88*spa14*spb12*t75 + s12*(t155*t159 + d89*spa45*t163 + d91*t184 + d80*t231 + d83*t239 - d82*spa34*t25 + d84*spa45*t28 + d92*t59*t75 + d81*spa34*t92 + d79*t58*t92) + d87*spa24*spb12*(spa15*spa23*spa24 + spa24*t186 + spa25*t246*T(2)); 
complex<T> t24 = d84*s45*t129 + d96*s45*t163 + d98*s45*t169 + d97*s45*t200 + d109*s45*t210 + d117*spb45*t231 + d120*spb45*t239 - d108*s45*t25 - d119*spa34*spb45*t25 + s45*spa24*t39 + d118*spa34*spb45*t92 + d116*spb45*t58*t92 - d115*spb45*t155*T(3); 
complex<T> t41 = d65*(spa12*spa45*spb15 - spb13*t246); 
complex<T> t87 = -t115; 
complex<T> t141 = d135*s12*s25*t184 + d136*s12*spb25*t100*t25 + d134*s25*spa24*t176*t30; 
complex<T> t175 = (d137*t163 + d139*t169 + d138*t200)*t85; 
complex<T> t192 = d47*t224 + t135*t75; 
complex<T> d40 = spa25*spa35*t62*t8; d40 = T(1)/d40;
complex<T> d69 = spa35*t14*t8; d69 = T(1)/d69;
complex<T> d94 = spa25*spa35*t8; d94 = T(1)/d94;
complex<T> d95 = spa25*spa35*t34; d95 = T(1)/d95;
complex<T> d105 = spa35*t8; d105 = T(1)/d105;
complex<T> d106 = spa35*t34; d106 = T(1)/d106;
complex<T> d127 = spa25*t31*t8; d127 = T(1)/d127;
complex<T> d128 = spa25*t31*t34; d128 = T(1)/d128;
complex<T> t1 = d56*t107*t146 + d61*t163 + d44*t169 + d45*t169 + d64*t169 + d49*spa14*t185 + d63*t108*t193 + d43*t200 + d50*t144*t26 + d54*t177*t48 - d48*t75 + d55*t177*t76 + d59*d60*t144*t78 + d60*spa14*t185*t87 + d42*t154*t92 + d62*t154*t92 + d53*t164*t92 + d52*t26*t58*t92 + d57*t146*t98 - d51*t155*t26*T(3) - t25*t41*T(3); 
complex<T> t2 = d77*spa14*t109*t133 + d76*d77*t144*t154 + d69*t110*t162 + d40*t156*t162 + d41*t156*t162 + d70*t156*t162 - d67*t185 + d72*t193*t51 + d73*t193*t51 - d66*t75 + d38*t46*t78 + d68*t46*t78 + d71*t176*t92 + d39*t78*t92 - d74*t25*t68*T(3); 
complex<T> t3 = d61*t109 + d11*t107*t146 + d57*t107*t146 - d27*spa24*t165 - d28*spa23*spa24*t165 + d63*t169 - d49*spa14*t185 + d64*t108*t193 + d62*t200 - d17*spb35*t25 + d19*spb35*t25 + d35*(spa34*t134 + t172)*t25 + d50*t155*t26 + d21*t25*t27 + d53*t164*t46 + d7*t164*t46 + d55*t177*t48 + d9*t177*t48 + d31*t129*t55 + d32*t129*t55 + d5*t26*t46*t58 + d52*t26*t46*t58 - d25*spa24*t27*(spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6) + d54*t177*t76 + d8*t177*t76 + d60*(spa14*t115*t185 + d59*t155*t78) + d23*t166*t27*t92 + t233*t76*t92 + d10*t146*t98 + d56*t146*t98 + d4*t155*t26*T(3) + d51*t155*t26*T(3) + t25*t41*T(3); 
complex<T> t4 = d70*t110*t162 + d69*t156*t162 + d12*spa14*t163 + d77*(d76*t154*t155 + spa14*t133*t163) + d14*t184 + d15*t184 + d72*t184 + d73*t184 + d67*t185 + d18*spb35*t25 + d37*(d36*spa45*t163 + t172*t25) + d20*t25*t27 - d26*spa24*t27*t28 + d13*t176*t46 + d71*t176*t46 + d29*t129*t55 + d30*t129*t55 - d24*spa24*t27*(spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6) - t135*t75 + d16*t59*t75 + d47*t59*t75 + d22*t166*t27*t92 + d68*t78*t92 + d74*t25*t68*T(3); 
complex<T> t23 = d84*s34*t129 + d95*s34*t156*t162 + d111*spa45*spb34*t163 + d113*spb34*t184 - d93*s34*t185 + d109*s34*t210 - d108*s34*t25 + s34*spa24*t39 + d112*spb34*t176*t46 + d114*spb34*t59*t75 + d94*s34*t78*t92; 
complex<T> t42 = s14*(d98*t169 + d97*t200 + d95*t246*t261 - d93*spb13*square(spa24) + d96*spb15*square(spa24) + d94*spa23*t92*square(spb13)); 
complex<T> t44 = d91*s25*t184 + d106*spb25*t246*t261 + d90*s25*t176*t46 - d104*spb13*spb25*square(spa24) + d107*spb25*t100*square(spa24) + d105*spa23*spb25*t92*square(spb13); 
complex<T> t126 = d40*t110*t162 + d41*t110*t162 + d44*t108*t193 + d45*t108*t193 + d42*t200 + d39*t46*t78 + d43*t154*t92 + d38*t78*t92; 
complex<T> t143 = t69*(d128*t156*t162 - d126*t185 + d127*t78*t92); 
complex<T> co1 = d99*spb15*t224; 
complex<T> co2 = d133*s12*spb15*t224; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t126*(*CI_users[1]->get_value(mc,ind,mu)) + t192*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + t3*(*CI_users[6]->get_value(mc,ind,mu)) + t22*(*CI_users[7]->get_value(mc,ind,mu)) + t42*(*CI_users[8]->get_value(mc,ind,mu)) + co1*(*CI_users[9]->get_value(mc,ind,mu)) + t21*(*CI_users[10]->get_value(mc,ind,mu)) + t44*(*CI_users[11]->get_value(mc,ind,mu)) + t23*(*CI_users[12]->get_value(mc,ind,mu)) + t24*(*CI_users[13]->get_value(mc,ind,mu)) + t20*(*CI_users[14]->get_value(mc,ind,mu)) + t153*(*CI_users[19]->get_value(mc,ind,mu)) + t143*(*CI_users[20]->get_value(mc,ind,mu)) + t19*(*CI_users[21]->get_value(mc,ind,mu)) + co2*(*CI_users[25]->get_value(mc,ind,mu)) + t141*(*CI_users[26]->get_value(mc,ind,mu)) + t175*(*CI_users[28]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qpqmmpgap_SLC_wCI::\
C2q2g1y_qpqmmpgap_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c5, c124));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpqmmpgap_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, p, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpqmmpgap SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 22:00:57 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb45 = SPB(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> s23 = S(2,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s24 = -(spa24*spb24);
complex<T> s14 = -(spa14*spb14);
complex<T> s35 = -(spa35*spb35);
complex<T> s15 = -(spa15*spb15);
complex<T> t5 = spa25*T(2); 
complex<T> t8 = square(s12 - s34); 
complex<T> t13 = spa15*spa45; 
complex<T> t20 = square(spb45); 
complex<T> t21 = square(spb15); 
complex<T> t23 = spa12*spa35; 
complex<T> t24 = spa14*spa25; 
complex<T> t27 = spa45*T(2); 
complex<T> t30 = s34*s35; 
complex<T> t39 = square(spa23); 
complex<T> t40 = square(spa13); 
complex<T> t41 = s12 - s34 - s35; 
complex<T> t42 = cube(spb15); 
complex<T> t43 = cube(spb45); 
complex<T> t45 = s12 + s25 - s34; 
complex<T> t49 = s14 - s35; 
complex<T> t50 = s12 + s14 - s35; 
complex<T> t52 = s15 - s34; 
complex<T> t53 = s12 + s15 - s34; 
complex<T> t64 = spa14*T(2); 
complex<T> t65 = cube(spa23); 
complex<T> t66 = -(spa13*T(3)); 
complex<T> t69 = s25*spb25; 
complex<T> t70 = s24*spb24; 
complex<T> t71 = spa14*spa45; 
complex<T> t74 = spa24*T(2); 
complex<T> t90 = spa12*T(2); 
complex<T> t97 = spa24*spa35; 
complex<T> t115 = spa15*spa24; 
complex<T> t116 = -(spb45*T(3)); 
complex<T> t128 = -(spa12*spb14); 
complex<T> t155 = spa25*spa34; 
complex<T> t159 = spa13*spa23; 
complex<T> t170 = spa35*spb15; 
complex<T> d14 = (s12 - s34)*spa12*spa45; d14 = T(1)/d14;
complex<T> d15 = (s12 - s35)*spa12*spa45; d15 = T(1)/d15;
complex<T> d31 = spa12*spa45; d31 = T(1)/d31;
complex<T> d32 = s12 - s34; d32 = T(1)/d32;
complex<T> d33 = s12 - s35; d33 = T(1)/d33;
complex<T> d66 = spa24*spa45; d66 = T(1)/d66;
complex<T> t6 = square(s12 + t49); 
complex<T> t7 = spa14*t41; 
complex<T> t9 = square(s12 + t52); 
complex<T> t22 = -t159; 
complex<T> t44 = -t97; 
complex<T> t46 = -t70; 
complex<T> t47 = -t69; 
complex<T> t48 = spa15*spa23 + spa13*t5; 
complex<T> t51 = spa23*(t115 + t24) + spa13*spa24*t5; 
complex<T> t56 = spa34*t128 + spb45*t97*T(3); 
complex<T> t59 = t116*t155 - spb15*t23; 
complex<T> t67 = -(t23*t42); 
complex<T> t68 = spa34*t43; 
complex<T> t83 = spa23*t20; 
complex<T> t87 = t21*(spa15*spa23 + spa13*spa25*T(2)); 
complex<T> t98 = -(t40*T(3)); 
complex<T> t99 = -t155; 
complex<T> t109 = spa14*t5; 
complex<T> t112 = spa23*t21; 
complex<T> t113 = t23*t40; 
complex<T> t114 = spa34*t64; 
complex<T> t137 = t39*t69; 
complex<T> t139 = spa34*t20; 
complex<T> t144 = spa13*t39; 
complex<T> t146 = spa34*t70; 
complex<T> t148 = d31*spb45; 
complex<T> t152 = -(t20*t39); 
complex<T> t153 = -t170; 
complex<T> t154 = spa24*t66; 
complex<T> t167 = spa25*t74; 
complex<T> t178 = spa34*t39; 
complex<T> d2 = spa34*t13*t90; d2 = T(1)/d2;
complex<T> d3 = spa14*t23*t27; d3 = T(1)/d3;
complex<T> d4 = (s12 - s35)*spa35*t24; d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*spa14*spa34*t45; d5 = T(1)/d5;
complex<T> d6 = (s12 - s34)*spa34*t24*t45; d6 = T(1)/d6;
complex<T> d8 = (s12 - s34)*spa14*spa34*square(t45); d8 = T(1)/d8;
complex<T> d10 = (s12 - s35)*spa35*t24*(s12 + t49); d10 = T(1)/d10;
complex<T> d11 = (s12 - s35)*spa35*(s12 + t49)*t71; d11 = T(1)/d11;
complex<T> d12 = (s12 - s34)*spa34*t13*(s12 + t52); d12 = T(1)/d12;
complex<T> d13 = (s12 - s34)*t24; d13 = T(1)/d13;
complex<T> d16 = (s12 - s34)*spa12*t41; d16 = T(1)/d16;
complex<T> d17 = (s12 - s35)*spa12*t41; d17 = T(1)/d17;
complex<T> d22 = t64*t8; d22 = T(1)/d22;
complex<T> d23 = spa14*t27*square(s12 - s35); d23 = T(1)/d23;
complex<T> d24 = t13*t8*T(2); d24 = T(1)/d24;
complex<T> d25 = (s12 - s34)*spa14*square(t41); d25 = T(1)/d25;
complex<T> d26 = (s12 - s35)*spa14*square(t41); d26 = T(1)/d26;
complex<T> d29 = (s12 - s34)*spa34*t13*t90; d29 = T(1)/d29;
complex<T> d30 = (s12 - s35)*spa14*t23*t27; d30 = T(1)/d30;
complex<T> d34 = spa35*t49*t71; d34 = T(1)/d34;
complex<T> d35 = spa35*t24*t49; d35 = T(1)/d35;
complex<T> d36 = spa35*t24*t49*(s12 + t49); d36 = T(1)/d36;
complex<T> d37 = spa35*t49*(s12 + t49)*t71; d37 = T(1)/d37;
complex<T> d38 = (s14 - s23)*t71; d38 = T(1)/d38;
complex<T> d39 = (s14 - s23)*(-s23 + t49)*t71; d39 = T(1)/d39;
complex<T> d40 = t49*(-s23 + t49)*t71; d40 = T(1)/d40;
complex<T> d41 = spa34*t13*t52; d41 = T(1)/d41;
complex<T> d42 = t13*t52; d42 = T(1)/d42;
complex<T> d43 = spa34*t13*t52*(s12 + t52); d43 = T(1)/d43;
complex<T> d44 = (s25 - s34)*spa34*t24; d44 = T(1)/d44;
complex<T> d46 = (s25 - s34)*spa14*spa34*t45; d46 = T(1)/d46;
complex<T> d47 = (s25 - s34)*spa34*t24*t45; d47 = T(1)/d47;
complex<T> d48 = (s25 - s34)*spa14*spa34*square(t45); d48 = T(1)/d48;
complex<T> d50 = spa34*spa35*t24; d50 = T(1)/d50;
complex<T> d51 = spa34*t13; d51 = T(1)/d51;
complex<T> d52 = spa35*t71; d52 = T(1)/d52;
complex<T> d53 = spa14*spa34*square(t45); d53 = T(1)/d53;
complex<T> d54 = spa34*t24*square(t45); d54 = T(1)/d54;
complex<T> d55 = spa14*spa34*cube(t45); d55 = T(1)/d55;
complex<T> d59 = square(t41); d59 = T(1)/d59;
complex<T> d60 = spa14*square(t41); d60 = T(1)/d60;
complex<T> d61 = t24*square(t41); d61 = T(1)/d61;
complex<T> d62 = spa14*cube(t41); d62 = T(1)/d62;
complex<T> d65 = spa45*square(s23 - t49); d65 = T(1)/d65;
complex<T> d68 = spa24*t13; d68 = T(1)/d68;
complex<T> d69 = t71*square(s23 - t49); d69 = T(1)/d69;
complex<T> d70 = spa14*square(t45); d70 = T(1)/d70;
complex<T> d71 = t24*square(t45); d71 = T(1)/d71;
complex<T> d72 = spa14*cube(t45); d72 = T(1)/d72;
complex<T> d74 = spa12*square(t41); d74 = T(1)/d74;
complex<T> d75 = spa12*spa14*square(t41); d75 = T(1)/d75;
complex<T> d76 = spa12*t24*square(t41); d76 = T(1)/d76;
complex<T> d81 = t13*t74; d81 = T(1)/d81;
complex<T> d82 = spa14*t27*square(s23 - t49); d82 = T(1)/d82;
complex<T> d86 = t90*square(t41); d86 = T(1)/d86;
complex<T> d87 = spa12*t64*square(t41); d87 = T(1)/d87;
complex<T> d89 = t64*cube(t41); d89 = T(1)/d89;
complex<T> t36 = d30*(spa34*t128 + spb45*t97); 
complex<T> t76 = d3*t56; 
complex<T> t79 = d34*spb12; 
complex<T> t82 = d2*t59; 
complex<T> t94 = t39*(d38*spb25 + d39*t47); 
complex<T> t100 = -t148; 
complex<T> t101 = spb14*t46; 
complex<T> t103 = -(d41*spb12); 
complex<T> t104 = -(d38*spb25); 
complex<T> t120 = -(d35*spb24); 
complex<T> t123 = s23*(d69*t137 + d68*t39); 
complex<T> t126 = t68*t97; 
complex<T> t132 = spb25*t48; 
complex<T> t138 = t113*t42; 
complex<T> t175 = t112*t98; 
complex<T> t206 = t137*t170; 
complex<T> d1 = spa34*t109*t23; d1 = T(1)/d1;
complex<T> d7 = spa34*t109*t8; d7 = T(1)/d7;
complex<T> d9 = t114*t45*t8; d9 = T(1)/d9;
complex<T> d18 = (s12 - s34)*spa12*t7; d18 = T(1)/d18;
complex<T> d19 = (s12 - s35)*spa12*t7; d19 = T(1)/d19;
complex<T> d20 = (s12 - s34)*spa12*spa25*t7; d20 = T(1)/d20;
complex<T> d21 = (s12 - s35)*spa12*spa25*t7; d21 = T(1)/d21;
complex<T> d27 = t7*t8*T(2); d27 = T(1)/d27;
complex<T> d28 = t7*square(s12 - s35)*T(2); d28 = T(1)/d28;
complex<T> d45 = t114*square(s25 - s34); d45 = T(1)/d45;
complex<T> d49 = t114*t45*square(s25 - s34); d49 = T(1)/d49;
complex<T> d56 = spa35*t24*t6; d56 = T(1)/d56;
complex<T> d57 = spa35*t6*t71; d57 = T(1)/d57;
complex<T> d58 = spa34*t13*t9; d58 = T(1)/d58;
complex<T> d63 = spa25*spa35*t6; d63 = T(1)/d63;
complex<T> d64 = spa35*spa45*t6; d64 = T(1)/d64;
complex<T> d67 = spa34*spa45*t9; d67 = T(1)/d67;
complex<T> d73 = t13*t9; d73 = T(1)/d73;
complex<T> d77 = t24*t6; d77 = T(1)/d77;
complex<T> d78 = t6*t71; d78 = T(1)/d78;
complex<T> d79 = t114*square(t45); d79 = T(1)/d79;
complex<T> d80 = t114*cube(t45); d80 = T(1)/d80;
complex<T> d83 = spa35*t5*t6; d83 = T(1)/d83;
complex<T> d84 = spa35*t27*t6; d84 = T(1)/d84;
complex<T> d85 = spa34*t27*t9; d85 = T(1)/d85;
complex<T> d88 = spa12*t109*square(t41); d88 = T(1)/d88;
complex<T> t14 = d44*spb15*t144 - d45*t112*t40 + d48*t40*t67 + d49*t40*t67 + d47*t22*t87 + d46*t112*t40*T(3); 
complex<T> t15 = d55*s12*t138 + d51*spb12*t144 - d52*spb12*t144 + d59*spb12*t152 + d53*s12*t175 + d57*s12*t146*t39 + d58*s12*spa35*t39*t47 - d50*spa13*spb12*t65 + d62*s12*t44*t68 + d56*s12*t65*t70 + d60*spb12*t154*t83 + d61*spb12*(spa13*t167 + spa23*(t115 + t24))*t83 + d54*s12*t159*t87; 
complex<T> t16 = d62*s34*t126 + d72*spb34*t138 + d74*s34*t152 + d70*spb34*t175 + d68*s34*t39 + d73*spa35*spb34*t39*t47 + d75*s34*t154*t83 + d76*s34*(spa13*t167 + spa23*(t115 + t24))*t83 + d71*spb34*t159*t87; 
complex<T> t17 = d62*s35*t126 + d69*s35*t137 + d74*s35*t152 + d78*spb35*t146*t39 + d77*spb35*t65*t70 + d75*s35*t154*t83 + d76*s35*(spa13*t167 + spa23*(t115 + t24))*t83; 
complex<T> t18 = t30*(d89*t126 + d86*t152 + d87*t154*t83 + d88*(spa13*t167 + spa23*(t115 + t24))*t83); 
complex<T> t77 = d79*s12; 
complex<T> t110 = s12*t101*(d84*t178 + d83*t65); 
complex<T> t136 = d43*spa35*t137 + t103*t144 - d42*spb24*t39; 
complex<T> t140 = d18*spa24; 
complex<T> t143 = d67*t206 + d66*spb15*t39; 
complex<T> t149 = d19*spa24; 
complex<T> t151 = d65*spb14*t137 + d64*t101*t178 + d63*t101*t65; 
complex<T> t162 = t132*t21; 
complex<T> t168 = d39*t137 + d40*t137 + t104*t39 + d37*t178*t46 + t120*t65 + d36*t46*t65 + t144*t79; 
complex<T> t1 = d48*t138 + d49*t138 + d8*t138 + d9*t138 + d41*spb12*t144 - d44*spb15*t144 + d46*t175 + d5*t175 + d24*spa35*t155*t20 + d7*t21*t22*t23 + d42*spb24*t39 + d13*spb45*t39 + d14*spb45*t39 + d16*t20*t39 + d32*(d2*t159*(t116*t155 - spb15*t23) + t100*t39) + d45*t112*t40 + d12*spa35*t39*t47 + d43*spa35*t39*t47 + d47*t159*t21*t48 + d6*t159*t21*t48 + d25*t44*t68 + d27*t44*t68 - d22*spa34*t83 - d20*(spa23*(t115 + t24) + spa13*spa24*t5)*t83 + d29*spa23*t66*(-(spb15*t23) + spb45*t99) + t140*t159*t20*T(3); 
complex<T> t2 = -(d4*spb14*t144) + d15*spb45*t39 + d11*t146*t39 + d37*t146*t39 + d17*t20*t39 + d23*t139*t44 + d40*t39*t47 + d35*spb24*t65 + d26*t44*t68 + d28*t44*t68 + d10*t65*t70 + d36*t65*t70 + d33*(t100*t39 + t22*t76) - t144*t79 - d21*(spa23*(t115 + t24) + spa13*spa24*t5)*t83 + t149*t159*t20*T(3) + t159*t36*T(3); 
complex<T> t3 = d25*t126 + d26*t126 + d27*t126 + d28*t126 + d12*spa35*t137 + d4*spb14*t144 + d16*t152 + d17*t152 + d7*t159*t21*t23 - d13*spb45*t39 - d14*spb45*t39 - d15*spb45*t39 + d32*(d2*t22*(t116*t155 - spb15*t23) + t148*t39) + d11*t178*t46 + d6*t21*t22*t48 + d10*t46*t65 + spa23*t36*t66 + d3*t39*t66 + d1*t65*t66 + d8*t40*t67 + d9*t40*t67 + d33*(t148*t39 + t159*t76) + d22*spa34*t83 + d20*(spa23*(t115 + t24) + spa13*spa24*t5)*t83 + d21*(spa23*(t115 + t24) + spa13*spa24*t5)*t83 + t140*t66*t83 + t149*t66*t83 + d23*t139*t97 + d24*spa35*t20*t99 + d2*t144*T(3) + d5*t112*t40*T(3) + d29*t159*(-(spb15*t23) + spb45*t99)*T(3); 
complex<T> t95 = d55*s25*t138 + d53*s25*t175 + d53*t162*t22; 
complex<T> t124 = d80*s12*s25*t138 + s25*t175*t77 + t162*t22*t77; 
complex<T> co1 = d81*s23*s34*t39; 
complex<T> co2 = d82*s23*s35*t137; 
complex<T> co3 = d85*s12*t206; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t168*(*CI_users[1]->get_value(mc,ind,mu)) + t136*(*CI_users[2]->get_value(mc,ind,mu)) + t94*(*CI_users[3]->get_value(mc,ind,mu)) + t14*(*CI_users[4]->get_value(mc,ind,mu)) + t1*(*CI_users[5]->get_value(mc,ind,mu)) + t2*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + t151*(*CI_users[8]->get_value(mc,ind,mu)) + t143*(*CI_users[9]->get_value(mc,ind,mu)) + t123*(*CI_users[10]->get_value(mc,ind,mu)) + t95*(*CI_users[11]->get_value(mc,ind,mu)) + t16*(*CI_users[12]->get_value(mc,ind,mu)) + t17*(*CI_users[13]->get_value(mc,ind,mu)) + t124*(*CI_users[15]->get_value(mc,ind,mu)) + co1*(*CI_users[16]->get_value(mc,ind,mu)) + co2*(*CI_users[17]->get_value(mc,ind,mu)) + t110*(*CI_users[22]->get_value(mc,ind,mu)) + co3*(*CI_users[24]->get_value(mc,ind,mu)) + t18*(*CI_users[25]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qpqmmmgap_SLC_wCI::\
C2q2g1y_qpqmmmgap_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c235));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c35, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c235));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c5, c124));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c35));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c2, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c5, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c35));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c3, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpqmmmgap_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, m, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpqmmmgap SLC");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 22:10:09 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s14 = -(spa14*spb14);
complex<T> s45 = -(spa45*spb45);
complex<T> s25 = S(2,5);
complex<T> s35 = -(spa35*spb35);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> t6 = spb14*T(2); 
complex<T> t8 = square(s12 - s35); 
complex<T> t9 = spb23*spb34; 
complex<T> t21 = square(spb15); 
complex<T> t22 = square(spa24); 
complex<T> t23 = square(spa34); 
complex<T> t24 = spb12*spb45; 
complex<T> t25 = -(spb25*T(3)); 
complex<T> t26 = spb14*spb23; 
complex<T> t28 = spb34*T(2); 
complex<T> t38 = cube(spa24); 
complex<T> t39 = square(spb25); 
complex<T> t40 = square(spa23); 
complex<T> t41 = -(spb13*spb15); 
complex<T> t43 = square(spb45); 
complex<T> t45 = s14 - s35; 
complex<T> t47 = cube(spa34); 
complex<T> t49 = spb15*spb24 + spb14*spb25*T(2); 
complex<T> t50 = spb23*T(2); 
complex<T> t51 = s14 - s25; 
complex<T> t52 = s15 - s23 + s45; 
complex<T> t53 = s12*s14; 
complex<T> t65 = cube(spa23); 
complex<T> t67 = -(spb12*spb25); 
complex<T> t70 = spb14*spb34; 
complex<T> t72 = s25*s35; 
complex<T> t88 = cube(spb15); 
complex<T> t90 = -(spa12*spb25); 
complex<T> t97 = spb13*spb25; 
complex<T> t100 = -(spa23*spb35); 
complex<T> t101 = -(spb14*T(2)); 
complex<T> t103 = s14*spa14; 
complex<T> t114 = spb13*spb45; 
complex<T> t116 = spb12*spb35; 
complex<T> t152 = spb23*spb35; 
complex<T> t157 = spb25*spb45; 
complex<T> d1 = s34*(s12 - s35)*spb12; d1 = T(1)/d1;
complex<T> d2 = s34*(s12 - s45)*spb12; d2 = T(1)/d2;
complex<T> d4 = s34*(s12 - s35)*spb12*spb23; d4 = T(1)/d4;
complex<T> d5 = s34*(s12 - s45)*spb12*spb23; d5 = T(1)/d5;
complex<T> d8 = (s12 - s35)*spb12*spb34; d8 = T(1)/d8;
complex<T> d9 = (s12 - s45)*spb12*spb34; d9 = T(1)/d9;
complex<T> d24 = (s12 - s35)*spb23*square(s34); d24 = T(1)/d24;
complex<T> d26 = (s12 - s45)*spb23*square(s34); d26 = T(1)/d26;
complex<T> d31 = spb12*spb34; d31 = T(1)/d31;
complex<T> d33 = s12 - s35; d33 = T(1)/d33;
complex<T> d34 = s12 - s45; d34 = T(1)/d34;
complex<T> d61 = square(spb34); d61 = T(1)/d61;
complex<T> d62 = spb23*square(spb34); d62 = T(1)/d62;
complex<T> d73 = spb23*cube(spb34); d73 = T(1)/d73;
complex<T> d88 = spb12*square(spb34); d88 = T(1)/d88;
complex<T> d89 = spb12*spb23*square(spb34); d89 = T(1)/d89;
complex<T> d108 = spb12*square(spb34)*T(2); d108 = T(1)/d108;
complex<T> t10 = spb35*t45; 
complex<T> t42 = s12 + t45; 
complex<T> t44 = -t116; 
complex<T> t48 = -(spb13*spb15*spb24) - spb15*t26 + t101*t97; 
complex<T> t57 = spa23*t116 + spa34*t114*T(3); 
complex<T> t58 = -t72; 
complex<T> t66 = -t114; 
complex<T> t69 = -t103; 
complex<T> t74 = -(d31*T(2)); 
complex<T> t81 = spb15*t22; 
complex<T> t85 = spb25*(spb15*spb24 + spb25*t6); 
complex<T> t115 = spb15*t23; 
complex<T> t117 = spb25*t38; 
complex<T> t128 = t65*t97; 
complex<T> t130 = spb35*t47; 
complex<T> t131 = spb12*t43; 
complex<T> t140 = spb35*t114; 
complex<T> t150 = t24*t38; 
complex<T> t151 = spb25*t40; 
complex<T> t153 = d89*spb13; 
complex<T> t158 = spa24*t21; 
complex<T> t170 = spb35*t21; 
complex<T> t174 = t38*t43; 
complex<T> t185 = spb25*t21; 
complex<T> t187 = spb15*t25; 
complex<T> t200 = spb15*t40; 
complex<T> d3 = (s12 - s35)*t26; d3 = T(1)/d3;
complex<T> d6 = s34*(s12 - s35)*spb12*t26; d6 = T(1)/d6;
complex<T> d7 = s34*(s12 - s45)*spb12*t26; d7 = T(1)/d7;
complex<T> d12 = t116*t70; d12 = T(1)/d12;
complex<T> d13 = t50*t8; d13 = T(1)/d13;
complex<T> d14 = spb12*t28*t8; d14 = T(1)/d14;
complex<T> d15 = (s12 - s45)*spb45*t26; d15 = T(1)/d15;
complex<T> d16 = t24*t9*T(2); d16 = T(1)/d16;
complex<T> d17 = t152*t24*t6; d17 = T(1)/d17;
complex<T> d18 = (s12 - s45)*spb45*t9; d18 = T(1)/d18;
complex<T> d19 = t152*t6*t8; d19 = T(1)/d19;
complex<T> d23 = s34*t50*t8; d23 = T(1)/d23;
complex<T> d25 = s34*t50*square(s12 - s45); d25 = T(1)/d25;
complex<T> d27 = t9*square(s12 - s45)*T(2); d27 = T(1)/d27;
complex<T> d30 = (s12 - s45)*t24*t9*T(2); d30 = T(1)/d30;
complex<T> d32 = spb35*t70; d32 = T(1)/d32;
complex<T> d35 = spb34*t45*t6; d35 = T(1)/d35;
complex<T> d36 = spb34*t6*square(t51); d36 = T(1)/d36;
complex<T> d37 = (-s25 + t45)*t51*t70; d37 = T(1)/d37;
complex<T> d38 = t45*(-s25 + t45)*t70; d38 = T(1)/d38;
complex<T> d39 = spb35*t28; d39 = T(1)/d39;
complex<T> d40 = spb34*t116; d40 = T(1)/d40;
complex<T> d41 = square(t45); d41 = T(1)/d41;
complex<T> d43 = spb35*t50*square(t45); d43 = T(1)/d43;
complex<T> d46 = spb34*t116*t6; d46 = T(1)/d46;
complex<T> d47 = t51*t70*square(s25 - t45); d47 = T(1)/d47;
complex<T> d48 = t45*t70*square(s25 - t45); d48 = T(1)/d48;
complex<T> d49 = spb34*(-s25 + t45)*t6*square(t51); d49 = T(1)/d49;
complex<T> d50 = spb34*(-s25 + t45)*t6*square(t45); d50 = T(1)/d50;
complex<T> d57 = (s15 - s23)*t9; d57 = T(1)/d57;
complex<T> d58 = (s15 - s23)*t52*t9; d58 = T(1)/d58;
complex<T> d59 = (s23 - s45)*t52*t9; d59 = T(1)/d59;
complex<T> d60 = (s23 - s45)*spb45*t9; d60 = T(1)/d60;
complex<T> d63 = t26*square(spb34); d63 = T(1)/d63;
complex<T> d66 = t114*t26; d66 = T(1)/d66;
complex<T> d67 = spb45*t9; d67 = T(1)/d67;
complex<T> d68 = spb35*spb45*t26; d68 = T(1)/d68;
complex<T> d69 = t114*t9; d69 = T(1)/d69;
complex<T> d75 = spb34*(-s25 + t45); d75 = T(1)/d75;
complex<T> d76 = spb34*square(s25 - t45); d76 = T(1)/d76;
complex<T> d77 = spb34*cube(-s25 + t45); d77 = T(1)/d77;
complex<T> d79 = t9*square(t52); d79 = T(1)/d79;
complex<T> d80 = spb34*square(t52); d80 = T(1)/d80;
complex<T> d81 = spb14*t114; d81 = T(1)/d81;
complex<T> d82 = spb34*t114; d82 = T(1)/d82;
complex<T> d83 = (-s25 + t45)*t70; d83 = T(1)/d83;
complex<T> d84 = t70*square(s25 - t45); d84 = T(1)/d84;
complex<T> d85 = t70*cube(-s25 + t45); d85 = T(1)/d85;
complex<T> d90 = spb12*t26*square(spb34); d90 = T(1)/d90;
complex<T> d95 = spb13*t26; d95 = T(1)/d95;
complex<T> d96 = spb13*t9; d96 = T(1)/d96;
complex<T> d97 = t9*square(t52)*T(2); d97 = T(1)/d97;
complex<T> d100 = t114*t6; d100 = T(1)/d100;
complex<T> d101 = t114*t28; d101 = T(1)/d101;
complex<T> d102 = spb34*(-s25 + t45)*t6; d102 = T(1)/d102;
complex<T> d103 = spb34*t6*square(s25 - t45); d103 = T(1)/d103;
complex<T> d104 = spb34*t6*cube(-s25 + t45); d104 = T(1)/d104;
complex<T> d109 = spb12*t50*square(spb34); d109 = T(1)/d109;
complex<T> d110 = spb12*spb23*t6*square(spb34); d110 = T(1)/d110;
complex<T> d111 = t50*cube(spb34); d111 = T(1)/d111;
complex<T> t5 = square(t42); 
complex<T> t16 = s35*s45*(d111*t140 + d109*spb13*t187 - d108*t21 + d110*spb15*(spb15*t26 - spb24*t41 - t101*t97)); 
complex<T> t27 = spb35*t42; 
complex<T> t64 = -t81; 
complex<T> t68 = cube(t42); 
complex<T> t80 = d16*t57; 
complex<T> t99 = -t150; 
complex<T> t105 = d39*spa12; 
complex<T> t121 = -(d57*spa14); 
complex<T> t123 = d14*spb14; 
complex<T> t125 = (d57*spa14 + d58*t103)*t21; 
complex<T> t129 = -t158; 
complex<T> t134 = -(d60*spa13); 
complex<T> t136 = d40*spb14; 
complex<T> t137 = -(d90*t48); 
complex<T> t139 = t39*t81; 
complex<T> t146 = d79*t69; 
complex<T> t164 = t117*t131; 
complex<T> t169 = t150*t39; 
complex<T> t172 = s12*(d101*t100*t21 + d100*spa23*t88); 
complex<T> t173 = t128*t44; 
complex<T> t175 = t21*t69; 
complex<T> t178 = t151*t41; 
complex<T> t208 = t115*T(3); 
complex<T> t210 = spb13*t115; 
complex<T> d42 = t10*t26; d42 = T(1)/d42;
complex<T> d44 = spb23*t10*t42; d44 = T(1)/d44;
complex<T> d45 = t10*t26*t42; d45 = T(1)/d45;
complex<T> d51 = spb12*spb34*t10*t6; d51 = T(1)/d51;
complex<T> d54 = spb34*t10*t42; d54 = T(1)/d54;
complex<T> d92 = t42*t70; d92 = T(1)/d92;
complex<T> t1 = d18*spa13*t170 + d60*spa13*t170 + d59*t175 + d60*spa12*t185 - d15*spa23*t185 + d9*spa34*t21 + d27*t140*t23 + d2*t21*t23 + d25*t140*t47 + d26*t140*t47 - d34*(d31*spa34*t21 + spb15*spb25*t80) - d15*spa13*t88 + d5*t208*t97 + d7*t115*(-(spb15*t26) + spb24*t41 + t101*t97) + d30*spb15*spb25*(spa34*t114 + spa23*t116)*T(3); 
complex<T> t34 = d51*(spa24*spb14*spb25 + spa23*t44); 
complex<T> t35 = s25*(d85*t173 + d84*t178 - d83*spa23*square(spb15)); 
complex<T> t61 = -t136; 
complex<T> t76 = -t105; 
complex<T> t113 = d47*t116*t128 + d49*t116*t128 + d36*t178 + d37*t200*t97; 
complex<T> t127 = d102*spa23*t21*t58 + d104*t173*t72 + d103*t178*t72; 
complex<T> t168 = d80*spa23*t175 + d82*t100*t21 + d81*spa23*t88; 
complex<T> t177 = t134*t170 + d58*t175 + t21*(d59*t103 + t121 + d60*t90); 
complex<T> t182 = t157*t64; 
complex<T> t231 = t146*t21; 
complex<T> d10 = (s12 - s35)*spb23*t27; d10 = T(1)/d10;
complex<T> d11 = (s12 - s35)*t26*t27; d11 = T(1)/d11;
complex<T> d20 = (s12 - s35)*t152*t5; d20 = T(1)/d20;
complex<T> d21 = t27*t50*t8; d21 = T(1)/d21;
complex<T> d22 = (s12 - s35)*spb34*t27; d22 = T(1)/d22;
complex<T> d28 = (s12 - s35)*spb34*spb35*t5; d28 = T(1)/d28;
complex<T> d29 = t27*t28*t8; d29 = T(1)/d29;
complex<T> d52 = spb23*t10*t5; d52 = T(1)/d52;
complex<T> d53 = t27*t50*square(t45); d53 = T(1)/d53;
complex<T> d55 = spb34*t10*t5; d55 = T(1)/d55;
complex<T> d56 = t27*t28*square(t45); d56 = T(1)/d56;
complex<T> d64 = t152*t5; d64 = T(1)/d64;
complex<T> d65 = spb35*t26*t5; d65 = T(1)/d65;
complex<T> d70 = t152*t68; d70 = T(1)/d70;
complex<T> d71 = t27*t70; d71 = T(1)/d71;
complex<T> d72 = spb34*spb35*t5; d72 = T(1)/d72;
complex<T> d74 = spb34*spb35*t68; d74 = T(1)/d74;
complex<T> d78 = spb34*t27; d78 = T(1)/d78;
complex<T> d86 = spb23*t5; d86 = T(1)/d86;
complex<T> d87 = t26*t5; d87 = T(1)/d87;
complex<T> d91 = spb23*t68; d91 = T(1)/d91;
complex<T> d93 = spb34*t5; d93 = T(1)/d93;
complex<T> d94 = spb34*t68; d94 = T(1)/d94;
complex<T> d98 = spb35*t5*t50; d98 = T(1)/d98;
complex<T> d99 = spb35*t50*t68; d99 = T(1)/d99;
complex<T> d105 = t27*t28; d105 = T(1)/d105;
complex<T> d106 = spb35*t28*t5; d106 = T(1)/d106;
complex<T> d107 = spb35*t28*t68; d107 = T(1)/d107;
complex<T> t2 = d13*spb35*t115 + spb35*t115*t123 - d18*spa13*t170 + d15*spa23*t185 + d30*(spa34*t114 + spa23*t116)*t187 + d3*spa34*t21 - d8*spa34*t21 - d9*spa34*t21 + d4*spb13*t187*t23 + d5*spb13*t187*t23 - d1*t21*t23 - d2*t21*t23 + d11*spb25*t49*t64 + d23*t130*t66 + d24*t130*t66 + d25*t130*t66 + d26*t130*t66 + d27*spb35*t23*t66 + d28*t174*t67 + d29*t174*t67 + d34*(d31*spa34*t21 + spb15*spb25*t80) + d22*t157*t81 + d19*spb25*t24*t81 - d12*t88 + d15*spa13*t88 + d17*t25*t88 + d6*t115*(spb15*t26 - spb24*t41 - t101*t97) + d7*t115*(spb15*t26 - spb24*t41 - t101*t97) + d20*t39*t99 + d21*t39*t99 + d33*(d32*spb45*t129 + d31*spa34*t21*T(2)) + d10*t139*T(3) + d16*t185*T(3); 
complex<T> t3 = d42*spb25*t158 + d47*t173 + d48*t173 + d49*t173 + d50*t173 + d37*t178 + d38*t178 + d35*spa23*t21 + d41*t139*t61 + d43*t39*t64 + d45*spb25*t49*t64 + d55*t174*t67 + d56*t174*t67 + d41*spb25*t158*t76 + d54*t157*t81 - d46*t88 + d36*t200*t97 + d52*t39*t99 + d53*t39*t99 + d44*t139*T(3) + t21*t34*T(3); 
complex<T> t4 = -(d13*spb35*t115) - spb35*t115*t123 + d48*t116*t128 + d50*t116*t128 + d42*spb25*t129 + d43*t139 + d41*t136*t139 + d32*d33*spb45*t158 + d41*spb25*t105*t158 + d28*t164 + d29*t164 + d55*t164 + d56*t164 + d20*t169 + d21*t169 + d52*t169 + d53*t169 + d22*t182 + d54*t182 - d35*spa23*t21 - d3*spa34*t21 + d8*spa34*t21 + d1*t21*t23 - d6*spb15*t115*t26 + d6*spb24*t115*t41 + d23*t140*t47 + d24*t140*t47 + d19*spb25*t24*t64 + d33*spa34*t21*t74 + d11*spb25*t49*t81 + d45*spb25*t49*t81 + d6*t101*t115*t97 + d38*t200*t97 + d4*t208*t97 - d10*t139*T(3) - d44*t139*T(3) - t21*t34*T(3); 
complex<T> t14 = d99*t169*t53 + d98*s12*spa14*spb25*t49*t64 - d98*t139*t53*T(3); 
complex<T> t15 = d73*s45*t140 + d96*spa45*t170 + s45*t153*t187 - d88*s45*t21 + s45*t231 - d95*spa45*t88 + d90*s45*spb15*(spb15*t26 - spb24*t41 - t101*t97); 
complex<T> t19 = d73*s35*t140 + d92*spa35*spb45*t158 + d94*spa35*t164 + d91*spa35*t169 + d85*s35*t173 + d84*s35*t178 + d93*spa35*t182 + s35*t153*t187 - d88*s35*t21 - d83*s35*spa23*t21 + d87*spa35*spb25*t49*t81 + d90*s35*spb15*(spb15*t26 - spb24*t41 - t101*t97) - d86*spa35*t139*T(3); 
complex<T> t96 = d105*s12*spa14*spb45*t129 + d107*t164*t53 + d106*t182*t53; 
complex<T> t133 = -(d64*T(3)); 
complex<T> t17 = d67*spa12*t185 + d62*spa12*spb13*t187 - d61*spa12*t21 - d63*spa12*spb15*spb24*t41 - d32*spa12*t88 + s12*(t133*t139 + d71*spb45*t158 + d74*t164 + d70*t169 + d69*t170 + d72*t182 + d73*spb35*t66 + d65*spb25*t49*t81 - d66*t88) + d68*t88*t90 - d63*spa12*spb15*t101*t97 + d63*spa12*t26*square(spb15); 
complex<T> t18 = d78*spa14*spb45*t129 + s14*(t133*t139 + d74*t164 + d70*t169 + d72*t182) + spa14*(d77*t173 + d76*t178 - d75*spa23*t21 + d64*spb25*t49*t64); 
complex<T> co1 = s15*t231; 
complex<T> co2 = d97*s15*s45*t175; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t125*(*CI_users[2]->get_value(mc,ind,mu)) + t177*(*CI_users[3]->get_value(mc,ind,mu)) + t113*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + t1*(*CI_users[6]->get_value(mc,ind,mu)) + t17*(*CI_users[7]->get_value(mc,ind,mu)) + t18*(*CI_users[8]->get_value(mc,ind,mu)) + co1*(*CI_users[9]->get_value(mc,ind,mu)) + t168*(*CI_users[10]->get_value(mc,ind,mu)) + t35*(*CI_users[11]->get_value(mc,ind,mu)) + t19*(*CI_users[12]->get_value(mc,ind,mu)) + t15*(*CI_users[13]->get_value(mc,ind,mu)) + co2*(*CI_users[16]->get_value(mc,ind,mu)) + t14*(*CI_users[17]->get_value(mc,ind,mu)) + t172*(*CI_users[21]->get_value(mc,ind,mu)) + t127*(*CI_users[24]->get_value(mc,ind,mu)) + t96*(*CI_users[26]->get_value(mc,ind,mu)) + t16*(*CI_users[28]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g1y_qpppqmgap_nf_wCI::\
C2q2g1y_qpppqmgap_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpppqmgap_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpppqmgap nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qppmqmgap_nf_wCI::\
C2q2g1y_qppmqmgap_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qppmqmgap_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qppmqmgap nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qpmpqmgap_nf_wCI::\
C2q2g1y_qpmpqmgap_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpmpqmgap_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpmpqmgap nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qpmmqmgap_nf_wCI::\
C2q2g1y_qpmmqmgap_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpmmqmgap_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpmmqmgap nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qpgapqmpp_L_wCI::\
C2q2g1y_qpgapqmpp_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpgapqmpp_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpgapqmpp L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1y_qpgapqmpm_L_wCI::\
C2q2g1y_qpgapqmpm_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpgapqmpm_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, gap, qm, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpgapqmpm L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 22:11:03 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = -(spa34*spb34);
complex<T> s24 = -(spa24*spb24);
complex<T> t4 = square(s15 - s23); 
complex<T> t5 = spa12*T(2); 
complex<T> t7 = (s23 - s45)*spa23; 
complex<T> t20 = square(spb14); 
complex<T> t21 = square(spa35); 
complex<T> t22 = square(spb24); 
complex<T> t23 = square(spa13); 
complex<T> t25 = s15 - s23 + s45; 
complex<T> t26 = -(spa15*T(2)); 
complex<T> t27 = s15 - s34; 
complex<T> t32 = cube(spb14); 
complex<T> t33 = cube(spa35); 
complex<T> t34 = cube(spb24); 
complex<T> t35 = spa23*spa45; 
complex<T> t36 = -(spa25*T(2)); 
complex<T> t42 = s15*s45; 
complex<T> t49 = spa35*T(2); 
complex<T> t53 = spa15*spa45; 
complex<T> t54 = spa35*spb24; 
complex<T> d6 = (s15 - s23)*s24*spa12; d6 = T(1)/d6;
complex<T> d9 = (s15 - s23)*spa12*square(s24); d9 = T(1)/d9;
complex<T> d15 = s15 - s23; d15 = T(1)/d15;
complex<T> d21 = spa12*square(spa24); d21 = T(1)/d21;
complex<T> d22 = spa12*spa24*spa34; d22 = T(1)/d22;
complex<T> d23 = spa12*cube(spa24); d23 = T(1)/d23;
complex<T> d30 = spa12*spa24; d30 = T(1)/d30;
complex<T> d31 = spa12*spa23*spa34; d31 = T(1)/d31;
complex<T> d32 = spa34*spa45*T(2); d32 = T(1)/d32;
complex<T> t24 = -t35; 
complex<T> t30 = -t42; 
complex<T> t46 = -t53; 
complex<T> t48 = spa25*t34; 
complex<T> t52 = t23*t32; 
complex<T> t57 = spa13*t20; 
complex<T> t58 = spa35*t26; 
complex<T> t63 = spa25*t22; 
complex<T> t72 = spb14*t21; 
complex<T> t77 = spa35*t36; 
complex<T> t86 = spb45*t33; 
complex<T> d1 = spa12*spa34*t35; d1 = T(1)/d1;
complex<T> d2 = spa23*t4*t5; d2 = T(1)/d2;
complex<T> d3 = (s15 - s23)*spa12*spa23*t25; d3 = T(1)/d3;
complex<T> d4 = (s15 - s23)*spa12*spa23*square(t25); d4 = T(1)/d4;
complex<T> d5 = spa23*t25*t4*t5; d5 = T(1)/d5;
complex<T> d7 = t5*square(t27); d7 = T(1)/d7;
complex<T> d8 = s24*spa12*t27; d8 = T(1)/d8;
complex<T> d10 = s24*t4*t5; d10 = T(1)/d10;
complex<T> d11 = s24*t5*square(t27); d11 = T(1)/d11;
complex<T> d12 = spa12*t27*square(s24); d12 = T(1)/d12;
complex<T> d13 = spa12*t35; d13 = T(1)/d13;
complex<T> d14 = spa24*spa34*t5; d14 = T(1)/d14;
complex<T> d16 = spa24*t4*t5; d16 = T(1)/d16;
complex<T> d17 = spa34*t35*t5; d17 = T(1)/d17;
complex<T> d18 = spa12*t25*t7; d18 = T(1)/d18;
complex<T> d19 = spa12*t7*square(t25); d19 = T(1)/d19;
complex<T> d20 = spa23*t25*t5*square(s23 - s45); d20 = T(1)/d20;
complex<T> d24 = spa12*spa23*spa34*t25; d24 = T(1)/d24;
complex<T> d25 = spa12*spa23*square(t25); d25 = T(1)/d25;
complex<T> d26 = spa12*spa23*cube(t25); d26 = T(1)/d26;
complex<T> d27 = spa12*spa34*t25; d27 = T(1)/d27;
complex<T> d28 = spa12*square(t25); d28 = T(1)/d28;
complex<T> d29 = spa12*cube(t25); d29 = T(1)/d29;
complex<T> d33 = spa45*t5; d33 = T(1)/d33;
complex<T> d34 = spa23*t5; d34 = T(1)/d34;
complex<T> d35 = t5*cube(spa24); d35 = T(1)/d35;
complex<T> d36 = spa24*t5; d36 = T(1)/d36;
complex<T> d37 = spa23*spa34*t25*t5; d37 = T(1)/d37;
complex<T> d38 = spa23*t5*cube(t25); d38 = T(1)/d38;
complex<T> d39 = spa23*spa34*t5; d39 = T(1)/d39;
complex<T> d40 = spa34*t35*T(2); d40 = T(1)/d40;
complex<T> t16 = d11*t35*t48 + d12*t35*t48 + d7*spa35*t63 + d8*t49*t63; 
complex<T> t17 = d30*spb34*t21 + d23*s34*spa25*t35 + d21*s34*t77; 
complex<T> t18 = s23*(-(d21*s34*spa25*spa35) + d36*spb34*t21 + d35*s34*spa25*t35); 
complex<T> t43 = d16*spb34; 
complex<T> t62 = t52*t53; 
complex<T> t66 = -t72; 
complex<T> t70 = t57*t58; 
complex<T> t73 = d24*spa13; 
complex<T> t79 = spa15*t57; 
complex<T> t1 = -(d17*t33) + d10*t35*t48 + d9*t35*t48 + d19*t46*t52 + d14*d15*t24*t54 + t24*t43*t54 + d20*t62 + d4*t62 + d5*t62 + d6*t49*t63 + d3*t70 - d2*spa35*t79 + d18*t49*t79 + d13*d15*spa15*t72*T(2); 
complex<T> t2 = d10*t24*t48 + d11*t24*t48 + d12*t24*t48 + d9*t24*t48 + d4*t46*t52 + d5*t46*t52 + d14*d15*t35*t54 + t35*t43*t54 - d7*spa35*t63 + d13*d15*t26*t72 + d6*t22*t77 + d8*t22*t77 + d2*spa35*t79 + d3*t49*t79 + d1*t33*T(2); 
complex<T> t12 = d20*t46*t52 + d19*t62 + d18*t70; 
complex<T> t13 = s15*(d22*t21 + d23*spa25*t24 + d21*spa25*t49 + d26*t62 + d25*t70 + t66*t73); 
complex<T> t14 = -(d22*s23*t21) + d23*s23*spa25*t35 + d29*spb23*t62 + d27*spa13*spb23*t66 + d28*spb23*t70 + d21*s23*t77; 
complex<T> t19 = d26*s45*t62 + d25*s45*t70 + s45*t66*t73 - d31*t86; 
complex<T> t45 = d38*t42*t62 + d37*spa13*t30*t72 + d25*spa35*t30*t79 - d39*s15*t86; 
complex<T> co1 = d32*spb12*spb23*t33; 
complex<T> co2 = d33*spb23*spb34*t33; 
complex<T> co3 = -(d32*spb12*spb23*t33); 
complex<T> co4 = d34*spb34*t86; 
complex<T> co5 = -(d40*s15*spb12*t33); 
complex<T> co6 = -co1; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t16*(*CI_users[2]->get_value(mc,ind,mu)) + t12*(*CI_users[3]->get_value(mc,ind,mu)) + t13*(*CI_users[4]->get_value(mc,ind,mu)) + t14*(*CI_users[5]->get_value(mc,ind,mu)) + t17*(*CI_users[6]->get_value(mc,ind,mu)) + t19*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + co6*(*CI_users[10]->get_value(mc,ind,mu)) + co4*(*CI_users[11]->get_value(mc,ind,mu)) + t18*(*CI_users[12]->get_value(mc,ind,mu)) + t45*(*CI_users[13]->get_value(mc,ind,mu)) + co5*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qpgapqmmp_L_wCI::\
C2q2g1y_qpgapqmmp_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpgapqmmp_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, gap, qm, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpgapqmmp L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 22:11:05 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> t1 = spa12*spa15; 
complex<T> t6 = square(spa34); 
complex<T> t9 = square(spa13); 
complex<T> t10 = square(spb15); 
complex<T> t11 = -(spa14*spb45); 
complex<T> t24 = spa13*spa14; 
complex<T> t28 = spa34*spb15; 
complex<T> d4 = spa15*spa23*spa35; d4 = T(1)/d4;
complex<T> d5 = (s12 + s15 - s34)*spa15*spa23; d5 = T(1)/d5;
complex<T> d6 = (s12 + s15 - s34)*spa23; d6 = T(1)/d6;
complex<T> d9 = spa15*spa45*T(2); d9 = T(1)/d9;
complex<T> d10 = (s12 + s15 - s34)*spa23*T(2); d10 = T(1)/d10;
complex<T> d13 = spa12*spa23*T(2); d13 = T(1)/d13;
complex<T> d14 = spa23*spa45*T(2); d14 = T(1)/d14;
complex<T> t14 = spa14*t6; 
complex<T> t19 = d5*spb25; 
complex<T> t22 = s34*t6; 
complex<T> t29 = spa14*t9; 
complex<T> t32 = s12*t6; 
complex<T> d1 = spa23*spa45*t1*T(2); d1 = T(1)/d1;
complex<T> d2 = (s23 - s45)*spa23*t1; d2 = T(1)/d2;
complex<T> d3 = spa23*t1*square(s23 - s45)*T(2); d3 = T(1)/d3;
complex<T> d7 = spa23*spa35*t1; d7 = T(1)/d7;
complex<T> d8 = spa23*t1; d8 = T(1)/d8;
complex<T> d11 = spa45*t1*T(2); d11 = T(1)/d11;
complex<T> d12 = spa23*t1*T(2); d12 = T(1)/d12;
complex<T> d15 = spa23*spa35*t1*T(2); d15 = T(1)/d15;
complex<T> t5 = -(t19*t32) - d4*spa13*spb12*t6; 
complex<T> t15 = -(d2*T(2)); 
complex<T> t16 = -(d7*spa13); 
complex<T> t18 = d3*spa45; 
complex<T> t4 = -(t10*t18*t29) + d2*t24*t28*T(2); 
complex<T> t20 = (d8*t11 + s45*t16)*t6; 
complex<T> t21 = t15*t24*t28 + t10*t18*t29 + d1*t14*T(3); 
complex<T> t26 = (t16 + t19)*t22; 
complex<T> co1 = d6*spb15*spb25*t6; 
complex<T> co2 = d9*spb12*spb23*t14; 
complex<T> co3 = d10*spb15*spb25*t32; 
complex<T> co4 = -(d11*s34*spb23*t14); 
complex<T> co5 = -(d9*spb12*spb23*t14); 
complex<T> co6 = d12*t11*t22; 
complex<T> co7 = d13*spb15*spb45*t14; 
complex<T> co8 = d14*spb12*spb15*t14; 
complex<T> co9 = -(d15*s45*spa13*t22); 
complex<T> co10 = -co2; 
complex<T> co11 = Complex(0,1); 
SeriesC<T> result = co11*(t21*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + t26*(*CI_users[4]->get_value(mc,ind,mu)) + t20*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + co4*(*CI_users[8]->get_value(mc,ind,mu)) + co10*(*CI_users[9]->get_value(mc,ind,mu)) + co6*(*CI_users[10]->get_value(mc,ind,mu)) + co7*(*CI_users[11]->get_value(mc,ind,mu)) + co8*(*CI_users[12]->get_value(mc,ind,mu)) + co9*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g1y_qpgapqmmm_L_wCI::\
C2q2g1y_qpgapqmmm_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);
	 vector<int> c5;  c5.push_back(ind[5-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
	 vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
	 vector<int> c45;  c45.push_back(ind[4-1]); c45.push_back(ind[5-1]);
	 vector<int> c51;  c51.push_back(ind[5-1]); c51.push_back(ind[1-1]);
         vector<int> c15 = c51;
         vector<int> c14;  c14.push_back(ind[1-1]); c14.push_back(ind[4-1]);
	 vector<int> c25;  c25.push_back(ind[2-1]); c25.push_back(ind[5-1]);
	 vector<int> c35;  c35.push_back(ind[3-1]); c35.push_back(ind[5-1]);

	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(ind[i-1]);}
	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(ind[i-1]);}
	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(ind[i-1]);}
	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(ind[i-1]);}
	                      c451.push_back(ind[1-1]);
	 vector<int> c145 = c451;
       	 vector<int> c512;    c512.push_back(ind[5-1]);
                          for(int i = 1; i<=2; i++) {c512.push_back(ind[i-1]);}
         vector<int> c125 = c512;
       	 vector<int> c412;    c412.push_back(ind[4-1]);
                          for(int i = 1; i<=2; i++) {c412.push_back(ind[i-1]);}
         vector<int> c124 = c412;
         vector<int> c134; c134.push_back(ind[1-1]); c134.push_back(ind[3-1]);
                            c134.push_back(ind[4-1]);    
         vector<int> c235; c235.push_back(ind[2-1]); c235.push_back(ind[3-1]);
                            c235.push_back(ind[5-1]);    
         vector<int> c245; c245.push_back(ind[2-1]); c245.push_back(ind[4-1]);
                            c245.push_back(ind[5-1]);    

CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1y_qpgapqmmm_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, gap, qm, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpgapqmmm L");
#endif
 
//#define TimeStamp "Sat 5 Feb 2011 22:11:06 on n4024"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = S(2,3);
complex<T> s12 = S(1,2);
complex<T> t4 = spb34*spb45; 
complex<T> t5 = square(spb12); 
complex<T> t7 = square(spa45); 
complex<T> t8 = square(spb24); 
complex<T> t11 = spa45*spb12; 
complex<T> d4 = spb15*T(2); d4 = T(1)/d4;
complex<T> d5 = spb34*T(2); d5 = T(1)/d5;
complex<T> d1 = (s15 - s23)*t4; d1 = T(1)/d1;
complex<T> d2 = t4*square(s15 - s23)*T(2); d2 = T(1)/d2;
complex<T> d3 = spb15*t4*T(2); d3 = T(1)/d3;
complex<T> t3 = d2*spb15*t7*t8 + d1*spb24*t11*T(2); 
complex<T> t16 = d3*t5; 
complex<T> t2 = -(d2*spb15*t7*t8) - d1*spb24*t11*T(2) - t16*T(3); 
complex<T> co1 = -(s12*s23*t16); 
complex<T> co2 = -(d4*spa34*spa45*t5); 
complex<T> co3 = -(d5*spa15*spa45*t5); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g1y_qpgapqmpp_nf_wCI::\
C2q2g1y_qpgapqmpp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpgapqmpp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpgapqmpp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qpgapqmpm_nf_wCI::\
C2q2g1y_qpgapqmpm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpgapqmpm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpgapqmpm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qpgapqmmp_nf_wCI::\
C2q2g1y_qpgapqmmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpgapqmmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpgapqmmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1y_qpgapqmmm_nf_wCI::\
C2q2g1y_qpgapqmmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1y_qpgapqmmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1y :  qpgapqmmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
 // *************** table of switch values ************* 
 
#define _C_qmmmqpgam_L C2q2g1y_5617_L
#define _C_qmmpqpgam_L C2q2g1y_5725_L
#define _C_qmpmqpgam_L C2q2g1y_5635_L
#define _C_qmppqpgam_L C2q2g1y_5743_L
#define _C_qmmqpmgam_SLC C2q2g1y_5257_SLC
#define _C_qmmqppgam_SLC C2q2g1y_5905_SLC
#define _C_qmpqpmgam_SLC C2q2g1y_5275_SLC
#define _C_qmpqppgam_SLC C2q2g1y_5923_SLC
#define _C_qmqpmmgam_SLC C2q2g1y_5197_SLC
#define _C_qmqpmpgam_SLC C2q2g1y_5845_SLC
#define _C_qmqppmgam_SLC C2q2g1y_5305_SLC
#define _C_qmqpppgam_SLC C2q2g1y_5953_SLC
#define _C_qmmmqpgam_nf C2q2g1y_5617_nf
#define _C_qmmpqpgam_nf C2q2g1y_5725_nf
#define _C_qmpmqpgam_nf C2q2g1y_5635_nf
#define _C_qmppqpgam_nf C2q2g1y_5743_nf
#define _C_qmgamqpmm_L C2q2g1y_97_L
#define _C_qmgamqpmp_L C2q2g1y_3985_L
#define _C_qmgamqppm_L C2q2g1y_745_L
#define _C_qmgamqppp_L C2q2g1y_4633_L
#define _C_qmgamqpmm_nf C2q2g1y_97_nf
#define _C_qmgamqpmp_nf C2q2g1y_3985_nf
#define _C_qmgamqppm_nf C2q2g1y_745_nf
#define _C_qmgamqppp_nf C2q2g1y_4633_nf
#define _C_qpppqmgap_L C2q2g1y_6824_L
#define _C_qppmqmgap_L C2q2g1y_6716_L
#define _C_qpmpqmgap_L C2q2g1y_6806_L
#define _C_qpmmqmgap_L C2q2g1y_6698_L
#define _C_qppqmpgap_SLC C2q2g1y_7184_SLC
#define _C_qppqmmgap_SLC C2q2g1y_6536_SLC
#define _C_qpmqmpgap_SLC C2q2g1y_7166_SLC
#define _C_qpmqmmgap_SLC C2q2g1y_6518_SLC
#define _C_qpqmppgap_SLC C2q2g1y_7244_SLC
#define _C_qpqmpmgap_SLC C2q2g1y_6596_SLC
#define _C_qpqmmpgap_SLC C2q2g1y_7136_SLC
#define _C_qpqmmmgap_SLC C2q2g1y_6488_SLC
#define _C_qpppqmgap_nf C2q2g1y_6824_nf
#define _C_qppmqmgap_nf C2q2g1y_6716_nf
#define _C_qpmpqmgap_nf C2q2g1y_6806_nf
#define _C_qpmmqmgap_nf C2q2g1y_6698_nf
#define _C_qpgapqmpp_L C2q2g1y_4604_L
#define _C_qpgapqmpm_L C2q2g1y_716_L
#define _C_qpgapqmmp_L C2q2g1y_3956_L
#define _C_qpgapqmmm_L C2q2g1y_68_L
#define _C_qpgapqmpp_nf C2q2g1y_4604_nf
#define _C_qpgapqmpm_nf C2q2g1y_716_nf
#define _C_qpgapqmmp_nf C2q2g1y_3956_nf
#define _C_qpgapqmmm_nf C2q2g1y_68_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmmmqpgam_L case 5617
 
#define _CASE_qmmpqpgam_L case 5725
 
#define _CASE_qmpmqpgam_L case 5635
 
#define _CASE_qmppqpgam_L case 5743
 
#define _CASE_qmmqpmgam_SLC case 5257
 
#define _CASE_qmmqppgam_SLC case 5905
 
#define _CASE_qmpqpmgam_SLC case 5275
 
#define _CASE_qmpqppgam_SLC case 5923
 
#define _CASE_qmqpmmgam_SLC case 5197
 
#define _CASE_qmqpmpgam_SLC case 5845
 
#define _CASE_qmqppmgam_SLC case 5305
 
#define _CASE_qmqpppgam_SLC case 5953
 
#define _CASE_qmmmqpgam_nf case 5617
 
#define _CASE_qmmpqpgam_nf case 5725
 
#define _CASE_qmpmqpgam_nf case 5635
 
#define _CASE_qmppqpgam_nf case 5743
 
#define _CASE_qmgamqpmm_L case 97
 
#define _CASE_qmgamqpmp_L case 3985
 
#define _CASE_qmgamqppm_L case 745
 
#define _CASE_qmgamqppp_L case 4633
 
#define _CASE_qmgamqpmm_nf case 97
 
#define _CASE_qmgamqpmp_nf case 3985
 
#define _CASE_qmgamqppm_nf case 745
 
#define _CASE_qmgamqppp_nf case 4633
 
#define _CASE_qpppqmgap_L case 6824
 
#define _CASE_qppmqmgap_L case 6716
 
#define _CASE_qpmpqmgap_L case 6806
 
#define _CASE_qpmmqmgap_L case 6698
 
#define _CASE_qppqmpgap_SLC case 7184
 
#define _CASE_qppqmmgap_SLC case 6536
 
#define _CASE_qpmqmpgap_SLC case 7166
 
#define _CASE_qpmqmmgap_SLC case 6518
 
#define _CASE_qpqmppgap_SLC case 7244
 
#define _CASE_qpqmpmgap_SLC case 6596
 
#define _CASE_qpqmmpgap_SLC case 7136
 
#define _CASE_qpqmmmgap_SLC case 6488
 
#define _CASE_qpppqmgap_nf case 6824
 
#define _CASE_qppmqmgap_nf case 6716
 
#define _CASE_qpmpqmgap_nf case 6806
 
#define _CASE_qpmmqmgap_nf case 6698
 
#define _CASE_qpgapqmpp_L case 4604
 
#define _CASE_qpgapqmpm_L case 716
 
#define _CASE_qpgapqmmp_L case 3956
 
#define _CASE_qpgapqmmm_L case 68
 
#define _CASE_qpgapqmpp_nf case 4604
 
#define _CASE_qpgapqmpm_nf case 716
 
#define _CASE_qpgapqmmp_nf case 3956
 
#define _CASE_qpgapqmmm_nf case 68
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2g1y_L( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmmmqpgam_L: return new 
                       C2q2g1y_qmmmqpgam_L_wCI(ind);
    _CASE_qmmpqpgam_L: return new 
                       C2q2g1y_qmmpqpgam_L_wCI(ind);
    _CASE_qmpmqpgam_L: return new 
                       C2q2g1y_qmpmqpgam_L_wCI(ind);
    _CASE_qmppqpgam_L: return new 
                       C2q2g1y_qmppqpgam_L_wCI(ind);
    _CASE_qmgamqpmm_L: return new 
                       C2q2g1y_qmgamqpmm_L_wCI(ind);
    _CASE_qmgamqpmp_L: return new 
                       C2q2g1y_qmgamqpmp_L_wCI(ind);
    _CASE_qmgamqppm_L: return new 
                       C2q2g1y_qmgamqppm_L_wCI(ind);
    _CASE_qmgamqppp_L: return new 
                       C2q2g1y_qmgamqppp_L_wCI(ind);
    _CASE_qpppqmgap_L: return new 
                       C2q2g1y_qpppqmgap_L_wCI(ind);
    _CASE_qppmqmgap_L: return new 
                       C2q2g1y_qppmqmgap_L_wCI(ind);
    _CASE_qpmpqmgap_L: return new 
                       C2q2g1y_qpmpqmgap_L_wCI(ind);
    _CASE_qpmmqmgap_L: return new 
                       C2q2g1y_qpmmqmgap_L_wCI(ind);
    _CASE_qpgapqmpp_L: return new 
                       C2q2g1y_qpgapqmpp_L_wCI(ind);
    _CASE_qpgapqmpm_L: return new 
                       C2q2g1y_qpgapqmpm_L_wCI(ind);
    _CASE_qpgapqmmp_L: return new 
                       C2q2g1y_qpgapqmmp_L_wCI(ind);
    _CASE_qpgapqmmm_L: return new 
                       C2q2g1y_qpgapqmmm_L_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g1y_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmmmqpgam_nf: return new 
                       C2q2g1y_qmmmqpgam_nf_wCI(ind);
    _CASE_qmmpqpgam_nf: return new 
                       C2q2g1y_qmmpqpgam_nf_wCI(ind);
    _CASE_qmpmqpgam_nf: return new 
                       C2q2g1y_qmpmqpgam_nf_wCI(ind);
    _CASE_qmppqpgam_nf: return new 
                       C2q2g1y_qmppqpgam_nf_wCI(ind);
    _CASE_qmgamqpmm_nf: return new 
                       C2q2g1y_qmgamqpmm_nf_wCI(ind);
    _CASE_qmgamqpmp_nf: return new 
                       C2q2g1y_qmgamqpmp_nf_wCI(ind);
    _CASE_qmgamqppm_nf: return new 
                       C2q2g1y_qmgamqppm_nf_wCI(ind);
    _CASE_qmgamqppp_nf: return new 
                       C2q2g1y_qmgamqppp_nf_wCI(ind);
    _CASE_qpppqmgap_nf: return new 
                       C2q2g1y_qpppqmgap_nf_wCI(ind);
    _CASE_qppmqmgap_nf: return new 
                       C2q2g1y_qppmqmgap_nf_wCI(ind);
    _CASE_qpmpqmgap_nf: return new 
                       C2q2g1y_qpmpqmgap_nf_wCI(ind);
    _CASE_qpmmqmgap_nf: return new 
                       C2q2g1y_qpmmqmgap_nf_wCI(ind);
    _CASE_qpgapqmpp_nf: return new 
                       C2q2g1y_qpgapqmpp_nf_wCI(ind);
    _CASE_qpgapqmpm_nf: return new 
                       C2q2g1y_qpgapqmpm_nf_wCI(ind);
    _CASE_qpgapqmmp_nf: return new 
                       C2q2g1y_qpgapqmmp_nf_wCI(ind);
    _CASE_qpgapqmmm_nf: return new 
                       C2q2g1y_qpgapqmmm_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g1y_SLC( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_qmmqpmgam_SLC: return new 
                       C2q2g1y_qmmqpmgam_SLC_wCI(ind);
    _CASE_qmmqppgam_SLC: return new 
                       C2q2g1y_qmmqppgam_SLC_wCI(ind);
    _CASE_qmpqpmgam_SLC: return new 
                       C2q2g1y_qmpqpmgam_SLC_wCI(ind);
    _CASE_qmpqppgam_SLC: return new 
                       C2q2g1y_qmpqppgam_SLC_wCI(ind);
    _CASE_qmqpmmgam_SLC: return new 
                       C2q2g1y_qmqpmmgam_SLC_wCI(ind);
    _CASE_qmqpmpgam_SLC: return new 
                       C2q2g1y_qmqpmpgam_SLC_wCI(ind);
    _CASE_qmqppmgam_SLC: return new 
                       C2q2g1y_qmqppmgam_SLC_wCI(ind);
    _CASE_qmqpppgam_SLC: return new 
                       C2q2g1y_qmqpppgam_SLC_wCI(ind);
    _CASE_qppqmpgap_SLC: return new 
                       C2q2g1y_qppqmpgap_SLC_wCI(ind);
    _CASE_qppqmmgap_SLC: return new 
                       C2q2g1y_qppqmmgap_SLC_wCI(ind);
    _CASE_qpmqmpgap_SLC: return new 
                       C2q2g1y_qpmqmpgap_SLC_wCI(ind);
    _CASE_qpmqmmgap_SLC: return new 
                       C2q2g1y_qpmqmmgap_SLC_wCI(ind);
    _CASE_qpqmppgap_SLC: return new 
                       C2q2g1y_qpqmppgap_SLC_wCI(ind);
    _CASE_qpqmpmgap_SLC: return new 
                       C2q2g1y_qpqmpmgap_SLC_wCI(ind);
    _CASE_qpqmmpgap_SLC: return new 
                       C2q2g1y_qpqmmpgap_SLC_wCI(ind);
    _CASE_qpqmmmgap_SLC: return new 
                       C2q2g1y_qpqmmmgap_SLC_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
