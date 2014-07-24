/*
*C_2q1g2l_wCI.cpp
*
* Created on 9/24, 2009
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
 


GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qppqmemep_L_wCI)

GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qpmqmemep_L_wCI)
 
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qppqmemep_nf_wCI)
 
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qpmqmemep_nf_wCI)

GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qppqmemep_nf_top_wCI)

GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qpmqmemep_nf_top_wCI)
 
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qpqmmemep_SLC_wCI)

GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qpqmpemep_SLC_wCI)
 
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qmqppemep_AX_wCI)

GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qmqpmemep_AX_wCI)
 
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qpqmpemep_AX_wCI)

GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g2l_qpqmmemep_AX_wCI)



C2q1g2l_qppqmemep_L_wCI::C2q1g2l_qppqmemep_L_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
} 
  
  
template <class T> SeriesC<T> 
     C2q1g2l_qppqmemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, p, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qppqmemep L");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:20:11 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q1g2l_qpmqmemep_L_wCI::C2q1g2l_qpmqmemep_L_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
} 
  
  
template <class T> SeriesC<T> 
     C2q1g2l_qpmqmemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpmqmemep L");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:20:12 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co3*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q1g2l_qppqmemep_nf_wCI::C2q1g2l_qppqmemep_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g2l_qppqmemep_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qppqmemep nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g2l_qpmqmemep_nf_wCI::C2q1g2l_qpmqmemep_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g2l_qpmqmemep_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpmqmemep nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g2l_qppqmemep_nf_top_wCI::C2q1g2l_qppqmemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g2l_qppqmemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qppqmemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g2l_qpmqmemep_nf_top_wCI::C2q1g2l_qpmqmemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g2l_qpmqmemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpmqmemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q1g2l_qpqmmemep_SLC_wCI::C2q1g2l_qpqmmemep_SLC_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c13, c245));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c3, c245));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c1, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
} 
  
  
template <class T> SeriesC<T> 
     C2q1g2l_qpqmmemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpqmmemep SLC");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:20:15 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co3*(t25*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t10*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t6*(*CI_users[6]->get_value(mc,ind,mu)) + t19*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q1g2l_qpqmpemep_SLC_wCI::C2q1g2l_qpqmpemep_SLC_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c3, c245));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c1, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
} 
  
  
template <class T> SeriesC<T> 
     C2q1g2l_qpqmpemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpqmpemep SLC");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:20:18 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co3*(t25*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t10*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + t5*(*CI_users[5]->get_value(mc,ind,mu)) + t6*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t11*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q1g2l_qmqppemep_AX_wCI::C2q1g2l_qmqppemep_AX_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q1g2l_qmqppemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qmqppemep AX");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:20:18 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co3*(*CI_users[1]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q1g2l_qmqpmemep_AX_wCI::C2q1g2l_qmqpmemep_AX_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q1g2l_qmqpmemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qmqpmemep AX");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:20:18 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> d1 = square(s12 - s45); d1 = T(1)/d1;
complex<T> co1 = -(d1*spa13*spa34*spb25); 
complex<T> co2 = d1*spa13*spa34*spb25; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q1g2l_qpqmpemep_AX_wCI::C2q1g2l_qpqmpemep_AX_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q1g2l_qpqmpemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpqmpemep AX");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:20:19 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co3*(*CI_users[1]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q1g2l_qpqmmemep_AX_wCI::C2q1g2l_qpqmmemep_AX_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q1g2l_qpqmmemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q1g2l :  qpqmmemep AX");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:20:19 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> d1 = square(s12 - s45); d1 = T(1)/d1;
complex<T> co1 = -(d1*spa23*spa34*spb15); 
complex<T> co2 = d1*spa23*spa34*spb15; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qppqmemep_L C2q1g2l_7400_L
#define _C_qpmqmemep_L C2q1g2l_7382_L
#define _C_qppqmemep_nf C2q1g2l_7400_nf
#define _C_qpmqmemep_nf C2q1g2l_7382_nf
#define _C_qppqmemep_nf_top C2q1g2l_7400_nf_top
#define _C_qpmqmemep_nf_top C2q1g2l_7382_nf_top
#define _C_qpqmmemep_SLC C2q1g2l_7352_SLC
#define _C_qpqmpemep_SLC C2q1g2l_7460_SLC
#define _C_qmqppemep_AX C2q1g2l_7465_AX
#define _C_qmqpmemep_AX C2q1g2l_7357_AX
#define _C_qpqmpemep_AX C2q1g2l_7460_AX
#define _C_qpqmmemep_AX C2q1g2l_7352_AX
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qppqmemep_L case 7400
 
#define _CASE_qpmqmemep_L case 7382
 
#define _CASE_qppqmemep_nf case 7400
 
#define _CASE_qpmqmemep_nf case 7382
 
#define _CASE_qppqmemep_nf_top case 7400
 
#define _CASE_qpmqmemep_nf_top case 7382
 
#define _CASE_qpqmmemep_SLC case 7352
 
#define _CASE_qpqmpemep_SLC case 7460
 
#define _CASE_qmqppemep_AX case 7465
 
#define _CASE_qmqpmemep_AX case 7357
 
#define _CASE_qpqmpemep_AX case 7460
 
#define _CASE_qpqmmemep_AX case 7352
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q1g2l_AX( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmqppemep_AX: return new 
                       C2q1g2l_qmqppemep_AX_wCI(ind);
    _CASE_qmqpmemep_AX: return new 
                       C2q1g2l_qmqpmemep_AX_wCI(ind);
    _CASE_qpqmpemep_AX: return new 
                       C2q1g2l_qpqmpemep_AX_wCI(ind);
    _CASE_qpqmmemep_AX: return new 
                       C2q1g2l_qpqmmemep_AX_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q1g2l_L( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qppqmemep_L: return new 
                       C2q1g2l_qppqmemep_L_wCI(ind);
    _CASE_qpmqmemep_L: return new 
                       C2q1g2l_qpmqmemep_L_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q1g2l_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qppqmemep_nf: return new 
                       C2q1g2l_qppqmemep_nf_wCI(ind);
    _CASE_qpmqmemep_nf: return new 
                       C2q1g2l_qpmqmemep_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q1g2l_nf_top( int hc,const std::vector<int>& 
ind) { 
    switch (hc) {
    _CASE_qppqmemep_nf_top: return new 
                       C2q1g2l_qppqmemep_nf_top_wCI(ind);
    _CASE_qpmqmemep_nf_top: return new 
                       C2q1g2l_qpmqmemep_nf_top_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q1g2l_SLC( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_qpqmmemep_SLC: return new 
                       C2q1g2l_qpqmmemep_SLC_wCI(ind);
    _CASE_qpqmpemep_SLC: return new 
                       C2q1g2l_qpqmpemep_SLC_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
