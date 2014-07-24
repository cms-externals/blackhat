/*
*C_2q2Q1y_wCI.cpp
*
* Created on 2/6, 2011
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
 


GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpQmQpqmgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpQpQmqmgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmQpQmqpgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmQmQpqpgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpQmQpqmgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpQpQmqmgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmQpQmqpgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmQmQpqpgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpqmQpQmgap_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpqmQmQpgap_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmqpQpQmgap_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmqpQmQpgap_sl_wCI)

GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmQpQmqpgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmQmQpqpgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpQmQpqmgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpQpQmqmgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmQpQmqpgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmQmQpqpgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpQmQpqmgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpQpQmqmgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmqpQmQpgam_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qmqpQpQmgam_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpqmQmQpgam_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q1y_qpqmQpQmgam_sl_wCI)


C2q2Q1y_qmQpQmqpgam_nf_wCI::\
C2q2Q1y_qmQpQmqpgam_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmQpQmqpgam_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQpQmqpgam nf");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:01:26 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb24); 
complex<T> d1 = spb15*spb23*spb45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1y_qmQmQpqpgam_nf_wCI::\
C2q2Q1y_qmQmQpqpgam_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmQmQpqpgam_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQmQpqpgam nf");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:01:27 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb34); 
complex<T> d1 = spb15*spb23*spb45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1y_qpQmQpqmgam_nf_wCI::\
C2q2Q1y_qpQmQpqmgam_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpQmQpqmgam_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qm, Qp, qm, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQmQpqmgam nf");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:01:28 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb13); 
complex<T> d1 = spb15*spb23*spb45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1y_qpQpQmqmgam_nf_wCI::\
C2q2Q1y_qpQpQmqmgam_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpQpQmqmgam_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qp, Qm, qm, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQpQmqmgam nf");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:01:28 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb12); 
complex<T> d1 = spb15*spb23*spb45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1y_qmQpQmqpgam_L_wCI::\
C2q2Q1y_qmQpQmqpgam_L_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
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
     C2q2Q1y_qmQpQmqpgam_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQpQmqpgam L");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:04:19 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s25 = -(spa25*spb25);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s13 = -(spa13*spb13);
complex<T> s35 = -(spa35*spb35);
complex<T> t3 = spb15*T(2); 
complex<T> t4 = spb12*spb34; 
complex<T> t5 = square(s12 - s45); 
complex<T> t13 = square(spb24); 
complex<T> t14 = -(spa13*spb14); 
complex<T> t15 = -(spb13*T(3)); 
complex<T> t16 = spb12*spb45; 
complex<T> t17 = (s12 - s45)*spb15; 
complex<T> t24 = square(spa35); 
complex<T> t25 = square(spb23); 
complex<T> t26 = square(spb25); 
complex<T> t27 = square(spb34); 
complex<T> t28 = square(spb45); 
complex<T> t30 = spb13*spb24; 
complex<T> t31 = -(spb14*spb23) - spb13*spb24; 
complex<T> t33 = s23 - s45; 
complex<T> t34 = spa45*spb14; 
complex<T> t35 = square(spa13); 
complex<T> t36 = -(s45*T(3)) + s23*T(4); 
complex<T> t37 = spa13*spb12 - spa34*spb24; 
complex<T> t43 = -(s35*spb24) - spa13*spb12*spb34*T(2); 
complex<T> t47 = cube(spa35); 
complex<T> t63 = spa35*spb24; 
complex<T> t67 = cube(spb24); 
complex<T> t70 = spa12*spa15; 
complex<T> t74 = spb14*spb24; 
complex<T> t79 = spa13*T(3); 
complex<T> t84 = s25*spa25; 
complex<T> t96 = s12*s23; 
complex<T> t99 = spa23*spb13; 
complex<T> t113 = spa45*spb13; 
complex<T> t119 = spa35*spb45; 
complex<T> d1 = (s12 - s34)*s35*spb12*spb15; d1 = T(1)/d1;
complex<T> d3 = spb12*T(2); d3 = T(1)/d3;
complex<T> d5 = square(s12 - s34); d5 = T(1)/d5;
complex<T> d6 = (s12 - s34)*spb15*spb34; d6 = T(1)/d6;
complex<T> d10 = s12 - s34; d10 = T(1)/d10;
complex<T> d25 = spb15*spb34*square(s12 - s34); d25 = T(1)/d25;
complex<T> d31 = s23*spb15*T(4); d31 = T(1)/d31;
complex<T> d32 = s23*spb15*spb45*T(4); d32 = T(1)/d32;
complex<T> d33 = spb15*spb23*spb45*T(4); d33 = T(1)/d33;
complex<T> d39 = spb15*spb23*spb45*T(6); d39 = T(1)/d39;
complex<T> d42 = spb15*square(spb35); d42 = T(1)/d42;
complex<T> d43 = spb13*spb15*spb45; d43 = T(1)/d43;
complex<T> d44 = spb15*spb45*square(spb13); d44 = T(1)/d44;
complex<T> d45 = spb15*spb34*spb45; d45 = T(1)/d45;
complex<T> d46 = spb45*cube(spb35); d46 = T(1)/d46;
complex<T> d47 = spb15*spb34*cube(spb35); d47 = T(1)/d47;
complex<T> d48 = s23; d48 = T(1)/d48;
complex<T> d49 = s23*spb45; d49 = T(1)/d49;
complex<T> d50 = spb23*spb45; d50 = T(1)/d50;
complex<T> d51 = spb15; d51 = T(1)/d51;
complex<T> d52 = spb15*spb45; d52 = T(1)/d52;
complex<T> d53 = spb12*spb15*square(spb35); d53 = T(1)/d53;
complex<T> d56 = spb12*spb15*cube(spb35); d56 = T(1)/d56;
complex<T> d57 = s23*spb15; d57 = T(1)/d57;
complex<T> d58 = spb13*spb15; d58 = T(1)/d58;
complex<T> d59 = spb15*spb23; d59 = T(1)/d59;
complex<T> d60 = spb15*square(spb13); d60 = T(1)/d60;
complex<T> d62 = spb12*cube(spb35); d62 = T(1)/d62;
complex<T> d66 = spb34*square(s12 + s15 - s34)*T(2); d66 = T(1)/d66;
complex<T> d71 = s23*T(2); d71 = T(1)/d71;
complex<T> d74 = spb34*spb45*T(2); d74 = T(1)/d74;
complex<T> d75 = spb23*spb34*spb45*T(2); d75 = T(1)/d75;
complex<T> d76 = spb12*cube(spb35)*T(2); d76 = T(1)/d76;
complex<T> t44 = spa13*spb12*spb34 + spb23*t119; 
complex<T> t75 = t24*t28; 
complex<T> t76 = spb13*t25; 
complex<T> t77 = t26*t27; 
complex<T> t80 = d66*s12; 
complex<T> t105 = spb14*t13; 
complex<T> t106 = spb25*t24; 
complex<T> t110 = d31*T(3); 
complex<T> t114 = d3*spb23; 
complex<T> t116 = -(d44*t31); 
complex<T> t118 = -(t13*T(3)); 
complex<T> t123 = d53*spb25; 
complex<T> t128 = -(d43*T(3)); 
complex<T> t129 = d71*spa15; 
complex<T> t139 = spa15*t84; 
complex<T> t141 = spb24*t15; 
complex<T> t145 = spa15*(-(d50*t13) + d48*t63 + d49*spa13*t74); 
complex<T> t146 = spb24*t34; 
complex<T> t149 = d32*t79; 
complex<T> t151 = t13*t34; 
complex<T> t157 = spa45*t14; 
complex<T> t163 = spb14*t25; 
complex<T> d2 = s35*spb12*t17; d2 = T(1)/d2;
complex<T> d4 = t3*t4; d4 = T(1)/d4;
complex<T> d7 = t3*square(s12 - s34); d7 = T(1)/d7;
complex<T> d8 = spb34*t3; d8 = T(1)/d8;
complex<T> d9 = t16*T(2); d9 = T(1)/d9;
complex<T> d11 = spb45*t17; d11 = T(1)/d11;
complex<T> d12 = (s12 - s34)*t16; d12 = T(1)/d12;
complex<T> d13 = t16*t17; d13 = T(1)/d13;
complex<T> d14 = s13*spb45*t17; d14 = T(1)/d14;
complex<T> d15 = spb15*spb45*t4; d15 = T(1)/d15;
complex<T> d16 = (s12 - s34)*spb15*spb45*t4; d16 = T(1)/d16;
complex<T> d17 = spb45*t17*t4; d17 = T(1)/d17;
complex<T> d18 = spb23*spb45*t3*t4; d18 = T(1)/d18;
complex<T> d19 = (s12 - s34)*t16*square(s35); d19 = T(1)/d19;
complex<T> d20 = s35*t16*square(s12 - s34)*T(2); d20 = T(1)/d20;
complex<T> d21 = s35*t16*t5*T(2); d21 = T(1)/d21;
complex<T> d22 = (s12 - s45)*t16*square(s35); d22 = T(1)/d22;
complex<T> d23 = t16*t5*T(2); d23 = T(1)/d23;
complex<T> d24 = t3*t4*t5; d24 = T(1)/d24;
complex<T> d26 = spb35*t3*t4*square(s12 - s34); d26 = T(1)/d26;
complex<T> d27 = (s12 - s34)*s35*spb15*spb35*t4; d27 = T(1)/d27;
complex<T> d28 = spb35*t3*t4*t5; d28 = T(1)/d28;
complex<T> d29 = s35*spb35*t17*t4; d29 = T(1)/d29;
complex<T> d30 = (s12 - s45)*spb23*spb45*t3*t4; d30 = T(1)/d30;
complex<T> d34 = s23*t3; d34 = T(1)/d34;
complex<T> d35 = s23*t3*t33; d35 = T(1)/d35;
complex<T> d36 = s23*t3*square(t33); d36 = T(1)/d36;
complex<T> d37 = s23*spb45*t3; d37 = T(1)/d37;
complex<T> d38 = spb15*spb45*t33; d38 = T(1)/d38;
complex<T> d40 = s13*spb15*spb45*t33; d40 = T(1)/d40;
complex<T> d41 = t3*square(t33); d41 = T(1)/d41;
complex<T> d54 = spb15*t16; d54 = T(1)/d54;
complex<T> d55 = t16*cube(spb35); d55 = T(1)/d55;
complex<T> d61 = spb15*t4; d61 = T(1)/d61;
complex<T> d63 = spb15*t4*cube(spb35); d63 = T(1)/d63;
complex<T> d64 = spb13*spb45*t3; d64 = T(1)/d64;
complex<T> d65 = spb34*spb45*t3; d65 = T(1)/d65;
complex<T> d67 = t16*t3; d67 = T(1)/d67;
complex<T> d68 = spb45*t3*square(spb13); d68 = T(1)/d68;
complex<T> d69 = spb12*t3; d69 = T(1)/d69;
complex<T> d70 = spb12*spb23*t3; d70 = T(1)/d70;
complex<T> d72 = t4*T(2); d72 = T(1)/d72;
complex<T> d73 = spb23*t4*T(2); d73 = T(1)/d73;
complex<T> d77 = spb12*t3*square(spb35); d77 = T(1)/d77;
complex<T> d78 = spb12*t3*cube(spb35); d78 = T(1)/d78;
complex<T> t8 = s23*spb14*t116 - d52*spa23*t13 + d51*t63 + d52*spa13*t74 + s23*t128*t74; 
complex<T> t23 = d78*s45*spa34*t28*t76 + d76*s34*spa45*t77 + d69*spa34*t34*square(spb24) - d77*s34*s45*spb25*t30*T(3); 
complex<T> t39 = d18*T(3); 
complex<T> t41 = d65*spa12; 
complex<T> t42 = d67*spa34; 
complex<T> t49 = -t76; 
complex<T> t50 = -t77; 
complex<T> t53 = -t80; 
complex<T> t57 = d8*spa25; 
complex<T> t58 = -(d69*spa34); 
complex<T> t82 = d4*s45; 
complex<T> t91 = -t105; 
complex<T> t100 = -t114; 
complex<T> t111 = d38*spa13; 
complex<T> t126 = -t151; 
complex<T> t150 = d33*t118 + t110*t63 + t149*t74; 
complex<T> t178 = t13*t139; 
complex<T> t1 = -(d7*spb34*t106) - d12*spa35*t13 - d10*d9*spa35*t13 + d5*spb24*t114*t24 - d25*spa15*t119*t30 + d19*t47*t50 + d20*t47*t50 - d6*spa15*t74 + d26*t75*t76 + d27*t75*t76 + d5*t151*t82 + d16*s35*t91 - d15*t105*T(2) + d1*t106*t30*T(3) + d10*t13*t57*T(3); 
complex<T> t9 = -(d41*spa15*spa35*t16) - d40*spb14*t31*t35 + d35*t36*t63 + d36*s45*spb24*t157*T(3) - d34*t63*T(3) - d37*spa13*t74*T(3) - d36*s23*spb24*t157*T(4) + t111*t74*T(4) - d39*t13*T(13); 
complex<T> t10 = d61*t126 - d59*spa45*t13 + s45*t123*t141 + d57*spa13*t146 - d60*t31*t34 + d63*s45*t28*t49 - d57*s45*t63 + d62*spa45*t77 - d58*t146*T(3); 
complex<T> t11 = d45*spa12*t105 + s12*spb14*t116 + d42*spa12*spb25*t141 + d47*spa12*t28*t49 + d46*spa12*t50 + s12*t128*t74; 
complex<T> t22 = d55*s34*t50 + d56*spa34*t28*t76 - d54*spa34*spb14*square(spb24) - s34*t123*t30*T(3); 
complex<T> t55 = -t82; 
complex<T> t133 = d72*spa15*t126 + spb24*t129*t157 + s45*t129*t63 + d73*spa15*t113*t67; 
complex<T> t140 = t151*t58 + d70*spa34*t113*t67; 
complex<T> t161 = t178*t53 + d75*spb13*t67*t70 + d74*t70*t91; 
complex<T> t164 = spb13*t39; 
complex<T> t167 = s23*t105*t41 + t41*t67*t99 - d64*t74*t96*T(3); 
complex<T> t170 = t42*(s23*t105 + t67*t99); 
complex<T> t2 = d16*s35*t105 + d17*s35*t105 + d7*spb34*t106 + d12*spa35*t13 + d10*d9*spa35*t13 + d1*t106*t141 + d2*t106*t141 + d5*spb24*t100*t24 + d24*spb45*t163*t24 + d25*spa15*t119*t30 - d14*spb14*t31*t35 + d30*t13*t15*t44 + d5*t151*t55 + d10*t118*t57 + d23*t43*t63 + t164*t67 + d6*spa15*t74 + d13*t37*t74 + d26*t49*t75 + d27*t49*t75 + d28*t49*t75 + d29*t49*t75 + d19*t47*t77 + d20*t47*t77 + d21*t47*t77 + d22*t47*t77 + d11*t74*t79 + d15*t91; 
complex<T> t12 = t150 + d41*spa15*spa35*t16 - d24*spb45*t163*t24 + d14*spb14*t31*t35 + d40*spb14*t31*t35 + d21*t47*t50 + d22*t47*t50 - d35*t36*t63 - d23*t43*t63 + t164*t67 - d13*t37*t74 + d28*t75*t76 + d29*t75*t76 + d17*s35*t91 + d2*t106*t30*T(3) + d30*spb13*t13*t44*T(3) - d11*spa13*t74*T(3) - t111*t74*T(4) + d36*spa13*t146*(s45*T(3) - s23*T(4)); 
complex<T> co1 = t178*t80; 
complex<T> co2 = -(d68*spb14*t31*t96); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t150*(*CI_users[1]->get_value(mc,ind,mu)) + t9*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t12*(*CI_users[4]->get_value(mc,ind,mu)) + t11*(*CI_users[5]->get_value(mc,ind,mu)) + t145*(*CI_users[6]->get_value(mc,ind,mu)) + t8*(*CI_users[7]->get_value(mc,ind,mu)) + t22*(*CI_users[8]->get_value(mc,ind,mu)) + t10*(*CI_users[9]->get_value(mc,ind,mu)) + t167*(*CI_users[10]->get_value(mc,ind,mu)) + co1*(*CI_users[11]->get_value(mc,ind,mu)) + t170*(*CI_users[12]->get_value(mc,ind,mu)) + co2*(*CI_users[13]->get_value(mc,ind,mu)) + t140*(*CI_users[14]->get_value(mc,ind,mu)) + t133*(*CI_users[15]->get_value(mc,ind,mu)) + t161*(*CI_users[16]->get_value(mc,ind,mu)) + t23*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qmQmQpqpgam_L_wCI::\
C2q2Q1y_qmQmQpqpgam_L_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmQmQpqpgam_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQmQpqpgam L");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:04:24 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb13 = SPB(1,3);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = S(1,2);
complex<T> t1 = spb15*T(2); 
complex<T> t6 = -(spb34*T(3)); 
complex<T> t10 = square(spb34); 
complex<T> t11 = -(s45*T(3)) + s23*T(4); 
complex<T> t14 = spa15*spa45; 
complex<T> t17 = spa12*spb14; 
complex<T> t18 = -(s23*T(4)); 
complex<T> t27 = spa25*spb34; 
complex<T> d1 = s23*spb15*T(4); d1 = T(1)/d1;
complex<T> d2 = s23*spb15*spb45*T(4); d2 = T(1)/d2;
complex<T> d3 = spb15*spb23*spb45*T(4); d3 = T(1)/d3;
complex<T> d8 = spb15*spb23*spb45*T(3); d8 = T(1)/d8;
complex<T> d11 = s23; d11 = T(1)/d11;
complex<T> d12 = s23*spb45; d12 = T(1)/d12;
complex<T> d13 = spb23*spb45; d13 = T(1)/d13;
complex<T> d14 = spb15; d14 = T(1)/d14;
complex<T> d15 = spb15*spb45; d15 = T(1)/d15;
complex<T> d16 = s23*spb15; d16 = T(1)/d16;
complex<T> d17 = spb15*spb23; d17 = T(1)/d17;
complex<T> d20 = s23*T(2); d20 = T(1)/d20;
complex<T> d21 = spb23*T(2); d21 = T(1)/d21;
complex<T> d22 = spb23*spb45*T(2); d22 = T(1)/d22;
complex<T> t5 = -t17; 
complex<T> t13 = -(d3*T(3)); 
complex<T> t16 = -t27; 
complex<T> t32 = d2*t17; 
complex<T> t36 = d1*spa25; 
complex<T> t39 = -(spa23*t10); 
complex<T> t43 = -(spa15*t10); 
complex<T> t46 = -(spa45*t10); 
complex<T> d4 = s23*t1; d4 = T(1)/d4;
complex<T> d5 = s23*(s23 - s45)*t1; d5 = T(1)/d5;
complex<T> d6 = s23*t1*square(s23 - s45); d6 = T(1)/d6;
complex<T> d7 = s23*spb45*t1; d7 = T(1)/d7;
complex<T> d9 = t1*square(s23 - s45); d9 = T(1)/d9;
complex<T> d10 = (s23 - s45)*spb23*spb45*t1; d10 = T(1)/d10;
complex<T> d18 = spb45*t1; d18 = T(1)/d18;
complex<T> d19 = spb23*t1; d19 = T(1)/d19;
complex<T> t23 = spb34*t5; 
complex<T> t24 = d6*(t18 + s45*T(3)); 
complex<T> t25 = d21*t10*t14 + d20*s45*spa15*t16 + d20*spb34*t14*t17; 
complex<T> t29 = d9*spb45; 
complex<T> t45 = t10*t13 + (t32 + t36)*t6; 
complex<T> t3 = spa45*t23*t24 + d5*t11*t27 - spa15*spa25*spb13*t29 - d10*spa15*spb13*spb45*t6 + t32*t6 + t36*t6 + d10*spb23*t5*t6 + d3*t10*T(3); 
complex<T> t4 = d5*t11*t16 + spa45*spb34*t17*t24 + spa15*spa25*spb13*t29 - d8*t10*T(2) - d10*spa15*spb13*spb34*spb45*T(3) + d7*spb34*t17*T(3) + d4*t27*T(3) + d10*spb23*spb34*t5*T(3); 
complex<T> t34 = d14*t16 + d15*(t23 + t39); 
complex<T> t38 = d11*spa15*t16 + d12*spa15*t23 + d13*t43; 
complex<T> t42 = d16*spa45*t23 + d16*s45*t27 + d17*t46; 
complex<T> co1 = d18*s12*t39; 
complex<T> co2 = d18*s34*t39; 
complex<T> co3 = d19*s34*t46; 
complex<T> co4 = d22*s12*t43; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t45*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)) + t38*(*CI_users[3]->get_value(mc,ind,mu)) + t34*(*CI_users[4]->get_value(mc,ind,mu)) + t42*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t25*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qpQmQpqmgam_L_wCI::\
C2q2Q1y_qpQmQpqmgam_L_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpQmQpqmgam_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qm, Qp, qm, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQmQpqmgam L");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:07:41 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb12 = SPB(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = -(spa12*spb12);
complex<T> s25 = -(spa25*spb25);
complex<T> t4 = spb45*T(2); 
complex<T> t6 = spb12*spb15; 
complex<T> t14 = square(spb13); 
complex<T> t15 = -(spa24*spb14); 
complex<T> t16 = -(spb24*T(3)); 
complex<T> t17 = spb12*spb34; 
complex<T> t20 = (s12 - s34)*spb45; 
complex<T> t29 = square(spa25); 
complex<T> t30 = square(spb23); 
complex<T> t31 = square(spb12); 
complex<T> t32 = square(spb15); 
complex<T> t33 = square(spb35); 
complex<T> t34 = spb13*spb24; 
complex<T> t36 = square(s25); 
complex<T> t37 = -(spb14*spb23) - spb13*spb24; 
complex<T> t38 = s15 - s34; 
complex<T> t40 = -(s15*T(3)); 
complex<T> t41 = s15*T(3) - s23*T(4); 
complex<T> t42 = square(spa24); 
complex<T> t44 = -(spa12*spb13) + spa24*spb34; 
complex<T> t48 = s25*spb13 + spa24*spb12*spb34*T(2); 
complex<T> t49 = cube(spa25); 
complex<T> t55 = spb25*spb45; 
complex<T> t69 = cube(spb13); 
complex<T> t70 = spa15*spb14; 
complex<T> t76 = spa34*spa45; 
complex<T> t79 = spb13*spb14; 
complex<T> t84 = spa24*T(3); 
complex<T> t94 = spa25*spb13; 
complex<T> t102 = spa15*spb24; 
complex<T> t114 = spa23*spb24; 
complex<T> t170 = spa25*spb15; 
complex<T> d1 = (s12 - s34)*spb15*spb34; d1 = T(1)/d1;
complex<T> d4 = spb15*spb34*T(2); d4 = T(1)/d4;
complex<T> d6 = s12 - s34; d6 = T(1)/d6;
complex<T> d7 = spb34*T(2); d7 = T(1)/d7;
complex<T> d9 = square(s12 - s34); d9 = T(1)/d9;
complex<T> d11 = spb12*spb45*square(s12 - s34); d11 = T(1)/d11;
complex<T> d21 = s23*spb45*T(4); d21 = T(1)/d21;
complex<T> d24 = (s15 - s23)*spb15*spb45; d24 = T(1)/d24;
complex<T> d25 = s23*spb15*spb45*T(4); d25 = T(1)/d25;
complex<T> d27 = spb15*spb23*spb45*T(4); d27 = T(1)/d27;
complex<T> d28 = (s15 - s23)*s24*spb15*spb45; d28 = T(1)/d28;
complex<T> d41 = spb15*spb23*spb45*T(6); d41 = T(1)/d41;
complex<T> d47 = s23*spb45; d47 = T(1)/d47;
complex<T> d48 = spb23*spb45; d48 = T(1)/d48;
complex<T> d49 = spb24*spb45; d49 = T(1)/d49;
complex<T> d50 = spb45*square(spb24); d50 = T(1)/d50;
complex<T> d53 = spb45; d53 = T(1)/d53;
complex<T> d54 = spb15*spb45; d54 = T(1)/d54;
complex<T> d55 = spb15*spb24*spb45; d55 = T(1)/d55;
complex<T> d56 = spb15*spb45*square(spb24); d56 = T(1)/d56;
complex<T> d61 = s23; d61 = T(1)/d61;
complex<T> d62 = s23*spb15; d62 = T(1)/d62;
complex<T> d63 = spb15*spb23; d63 = T(1)/d63;
complex<T> d65 = s23*T(2); d65 = T(1)/d65;
complex<T> d70 = spb12*spb35*T(2); d70 = T(1)/d70;
complex<T> t5 = spb34*t38; 
complex<T> t18 = square(t38); 
complex<T> t35 = s12 + t38; 
complex<T> t39 = -t70; 
complex<T> t72 = t40 + s23*T(4); 
complex<T> t80 = t29*t32; 
complex<T> t82 = spb24*t30; 
complex<T> t88 = d70*s34; 
complex<T> t95 = t31*t33; 
complex<T> t100 = d21*T(3); 
complex<T> t107 = spb35*t29; 
complex<T> t110 = -(d27*T(3)); 
complex<T> t111 = d24*spa24; 
complex<T> t119 = spb14*t14; 
complex<T> t136 = spa45*t14; 
complex<T> t142 = spb13*t16; 
complex<T> t143 = d65*spa45; 
complex<T> t145 = spa15*t15; 
complex<T> t152 = spa24*t70; 
complex<T> t158 = spb14*t30; 
complex<T> t190 = s23*t79; 
complex<T> d5 = spb12*t4; d5 = T(1)/d5;
complex<T> d8 = t17*t4; d8 = T(1)/d8;
complex<T> d10 = spb12*t20; d10 = T(1)/d10;
complex<T> d12 = spb34*spb45*t6; d12 = T(1)/d12;
complex<T> d14 = spb25*t17*t4*square(s12 - s34); d14 = T(1)/d14;
complex<T> d16 = t4*square(s12 - s34); d16 = T(1)/d16;
complex<T> d22 = s23*(-s15 + s23)*t4; d22 = T(1)/d22;
complex<T> d23 = s23*t4*square(s15 - s23); d23 = T(1)/d23;
complex<T> d26 = spb15*spb45*t38; d26 = T(1)/d26;
complex<T> d29 = s24*spb15*spb45*t38; d29 = T(1)/d29;
complex<T> d32 = spb23*spb34*t4*t6; d32 = T(1)/d32;
complex<T> d35 = t4*square(s15 - s23); d35 = T(1)/d35;
complex<T> d39 = s23*t4; d39 = T(1)/d39;
complex<T> d40 = s23*spb15*t4; d40 = T(1)/d40;
complex<T> d64 = spb15*spb34*t4; d64 = T(1)/d64;
complex<T> d66 = t4*t6; d66 = T(1)/d66;
complex<T> d67 = spb15*t4*square(spb24); d67 = T(1)/d67;
complex<T> d68 = t6*T(2); d68 = T(1)/d68;
complex<T> d69 = spb23*t6*T(2); d69 = T(1)/d69;
complex<T> d71 = spb15*spb24*t4; d71 = T(1)/d71;
complex<T> d72 = t17*T(2); d72 = T(1)/d72;
complex<T> d73 = spb23*t17*T(2); d73 = T(1)/d73;
complex<T> d75 = spb34*t4; d75 = T(1)/d75;
complex<T> d77 = spb23*spb34*t4; d77 = T(1)/d77;
complex<T> t3 = square(t35); 
complex<T> t8 = -(d54*spa23*t14) - d56*s23*spb14*t37 + d54*spa24*t79 + d53*t94 - d55*t190*T(3); 
complex<T> t9 = -(d67*s23*s34*spb14*t37) + d66*spa34*(s23*t119 + t114*t69); 
complex<T> t45 = d64*spa12; 
complex<T> t51 = -t82; 
complex<T> t58 = d8*s15; 
complex<T> t59 = -t88; 
complex<T> t60 = d5*spa35; 
complex<T> t74 = d32*T(3); 
complex<T> t81 = -t95; 
complex<T> t106 = -t119; 
complex<T> t122 = d35*spa45; 
complex<T> t129 = -t136; 
complex<T> t134 = d72*t136*t39 + d73*spa45*t102*t69; 
complex<T> t135 = t80*t82; 
complex<T> t141 = t49*t95; 
complex<T> t154 = t107*t142; 
complex<T> t157 = t107*t34; 
complex<T> t161 = t110*t14 + d25*t79*t84 + t100*t94; 
complex<T> t165 = t143*(spb13*t145 + s15*t94); 
complex<T> d3 = spb15*spb34*t35*square(s12 - s34)*T(2); d3 = T(1)/d3;
complex<T> d13 = spb34*t20*t35*t6; d13 = T(1)/d13;
complex<T> d15 = spb25*t17*t20*t35; d15 = T(1)/d15;
complex<T> d17 = spb34*t20*t35; d17 = T(1)/d17;
complex<T> d18 = spb15*spb34*t18*T(2); d18 = T(1)/d18;
complex<T> d20 = spb15*spb34*t18*t35*T(2); d20 = T(1)/d20;
complex<T> d30 = spb45*t35*t5*t6; d30 = T(1)/d30;
complex<T> d31 = t17*t18*t4; d31 = T(1)/d31;
complex<T> d33 = spb25*t17*t18*t4; d33 = T(1)/d33;
complex<T> d34 = spb12*t35*t5*t55; d34 = T(1)/d34;
complex<T> d36 = spb15*spb45*t5; d36 = T(1)/d36;
complex<T> d37 = spb23*t4*t5*t6; d37 = T(1)/d37;
complex<T> d38 = spb45*t35*t5; d38 = T(1)/d38;
complex<T> d42 = spb15*spb34*cube(t35); d42 = T(1)/d42;
complex<T> d46 = spb34*cube(t35); d46 = T(1)/d46;
complex<T> d57 = spb15*cube(t35); d57 = T(1)/d57;
complex<T> d74 = spb34*cube(t35)*T(2); d74 = T(1)/d74;
complex<T> t10 = -(spb34*t122*t170) + d23*spb13*t145*t41 + d28*spb14*t37*t42 + d22*t72*t94 - d40*spa24*t79*T(3) - d39*t94*T(3) - t111*t79*T(4) - d41*t14*T(13); 
complex<T> t25 = -(d37*(spa24*spb12*spb34 + spb23*t170)); 
complex<T> t140 = t45*(s23*t119 + t114*t69); 
complex<T> t153 = t136*t59 + d68*t106*t76 + d69*spb24*t69*t76; 
complex<T> t156 = d63*t129 + d62*spa24*spa45*t79 + d61*spa45*t94; 
complex<T> t164 = spb24*t74; 
complex<T> d2 = (s12 - s34)*spb15*spb34*t3; d2 = T(1)/d2;
complex<T> d19 = spb15*t3*t5; d19 = T(1)/d19;
complex<T> d43 = spb15*spb34*spb45*t3; d43 = T(1)/d43;
complex<T> d44 = spb34*t3*t55; d44 = T(1)/d44;
complex<T> d45 = spb34*spb45*t3; d45 = T(1)/d45;
complex<T> d51 = spb45*t17*t3; d51 = T(1)/d51;
complex<T> d52 = t17*t3*t55; d52 = T(1)/d52;
complex<T> d58 = spb45*t3*t6; d58 = T(1)/d58;
complex<T> d59 = spb12*t3*t55; d59 = T(1)/d59;
complex<T> d60 = spb45*t3; d60 = T(1)/d60;
complex<T> d76 = spb34*t3*t4; d76 = T(1)/d76;
complex<T> d78 = spb25*spb34*t3*t4; d78 = T(1)/d78;
complex<T> t1 = d12*t106 + d16*spb12*t107 - d1*spa25*t14 - d4*d6*spa25*t14 + d17*t154 + d38*t154 + t14*t16*t25 - d7*d9*spb13*spb23*t29 + d31*spb15*t158*t29 + d11*spa45*t170*t34 + d13*t119*t36 + d30*t119*t36 + d29*spb14*t37*t42 + d9*t14*t39*t58 + t164*t69 - d10*spa45*t79 - d36*t44*t79 + d14*t51*t80 + d15*t51*t80 + d33*t51*t80 + d34*t51*t80 + d19*t49*t81 + d2*t49*t81 + d20*t49*t81 + d3*t49*t81 - d18*t48*t94 + d6*t14*t60*T(3) - d26*spa24*t79*T(3); 
complex<T> t2 = -(d16*spb12*t107) + d14*t135 + d15*t135 + d1*spa25*t14 + d4*d6*spa25*t14 + d2*t141 + d3*t141 + d7*d9*spb13*spb23*t29 - d11*spa45*t170*t34 + d13*t106*t36 + d9*t14*t58*t70 + d10*spa45*t79 - d12*t119*T(2) + d17*t157*T(3) - d6*t14*t60*T(3); 
complex<T> t11 = d33*t135 + d34*t135 + d19*t141 + d20*t141 + t161 + spb34*t122*t170 - d31*spb15*t158*t29 + d30*t106*t36 + d23*spb13*t152*t41 - d28*spb14*t37*t42 - d29*spb14*t37*t42 + t164*t69 + d36*t44*t79 + d26*t79*t84 + d18*t48*t94 - d22*t72*t94 + d38*t157*T(3) + spb24*t14*t25*T(3) + t111*t79*T(4); 
complex<T> t12 = -(d48*spa15*t14) + d46*spa15*t141 + d47*spb13*t152 + d45*s15*t154 + d51*t14*t36*t39 + d50*t37*t39 + d52*s15*t51*t80 - d47*s15*t94 - d49*spb13*t70*T(3); 
complex<T> t13 = d60*spa34*t154 + d58*spa34*t119*t36 - d56*s34*spb14*t37 + d59*spa34*t51*t80 + d57*spa34*t49*t81 - d55*s34*t79*T(3); 
complex<T> t27 = d44*spa12*t135 + d42*s12*t49*t81 - d43*spa12*spb14*t36*square(spb13) - d45*s12*t157*T(3); 
complex<T> t28 = d78*s15*spa12*t135 + d74*s12*spa15*t141 + d76*s12*t157*t40 + d77*spa12*t102*t69 + d75*spa12*t39*square(spb13) + d76*spa12*t36*t70*square(spb13); 
complex<T> co1 = -(d71*s34*t190*T(3)); 
complex<T> co2 = t136*t88; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t11*(*CI_users[1]->get_value(mc,ind,mu)) + t10*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t161*(*CI_users[4]->get_value(mc,ind,mu)) + t27*(*CI_users[5]->get_value(mc,ind,mu)) + t12*(*CI_users[6]->get_value(mc,ind,mu)) + t8*(*CI_users[7]->get_value(mc,ind,mu)) + t13*(*CI_users[8]->get_value(mc,ind,mu)) + t156*(*CI_users[9]->get_value(mc,ind,mu)) + t140*(*CI_users[10]->get_value(mc,ind,mu)) + t165*(*CI_users[11]->get_value(mc,ind,mu)) + t9*(*CI_users[12]->get_value(mc,ind,mu)) + t153*(*CI_users[13]->get_value(mc,ind,mu)) + co1*(*CI_users[14]->get_value(mc,ind,mu)) + t134*(*CI_users[15]->get_value(mc,ind,mu)) + t28*(*CI_users[16]->get_value(mc,ind,mu)) + co2*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qpQpQmqmgam_L_wCI::\
C2q2Q1y_qpQpQmqmgam_L_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpQpQmqmgam_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qp, Qm, qm, gam}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQpQmqmgam L");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:07:46 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = S(3,4);
complex<T> s12 = S(1,2);
complex<T> t6 = -(spb12*T(3)); 
complex<T> t9 = square(spb12); 
complex<T> t10 = s15*T(3) - s23*T(4); 
complex<T> t14 = spa34*spb14; 
complex<T> t17 = -(s15*T(3)); 
complex<T> t21 = spa35*spb12; 
complex<T> t27 = spa15*spb12; 
complex<T> d1 = s23*spb45*T(4); d1 = T(1)/d1;
complex<T> d2 = s23*(-s15 + s23)*spb45*T(2); d2 = T(1)/d2;
complex<T> d3 = s23*spb45*square(s15 - s23)*T(2); d3 = T(1)/d3;
complex<T> d4 = s23*spb15*spb45*T(4); d4 = T(1)/d4;
complex<T> d5 = spb15*spb23*spb45*T(4); d5 = T(1)/d5;
complex<T> d6 = spb45*square(s15 - s23)*T(2); d6 = T(1)/d6;
complex<T> d7 = (s15 - s23)*spb15*spb23*spb45*T(2); d7 = T(1)/d7;
complex<T> d8 = s23*spb45*T(2); d8 = T(1)/d8;
complex<T> d9 = s23*spb15*spb45*T(2); d9 = T(1)/d9;
complex<T> d10 = spb15*spb23*spb45*T(3); d10 = T(1)/d10;
complex<T> d11 = s23*spb45; d11 = T(1)/d11;
complex<T> d12 = spb23*spb45; d12 = T(1)/d12;
complex<T> d13 = spb45; d13 = T(1)/d13;
complex<T> d14 = spb15*spb45; d14 = T(1)/d14;
complex<T> d15 = s23; d15 = T(1)/d15;
complex<T> d16 = s23*spb15; d16 = T(1)/d16;
complex<T> d17 = spb15*spb23; d17 = T(1)/d17;
complex<T> d18 = spb15*spb45*T(2); d18 = T(1)/d18;
complex<T> d19 = s23*T(2); d19 = T(1)/d19;
complex<T> d20 = spb15*spb23*T(2); d20 = T(1)/d20;
complex<T> d21 = spb23*T(2); d21 = T(1)/d21;
complex<T> d22 = spb23*spb45*T(2); d22 = T(1)/d22;
complex<T> t5 = -t14; 
complex<T> t12 = -(d5*T(3)); 
complex<T> t16 = -t21; 
complex<T> t29 = d6*spb24; 
complex<T> t31 = d1*spa35; 
complex<T> t32 = d4*t14; 
complex<T> t43 = -(spa23*t9); 
complex<T> t45 = -(spa45*t9); 
complex<T> t51 = spa15*t9; 
complex<T> t3 = -(spa35*spa45*spb15*t29) + d3*t10*t27*t5 + d7*spa45*spb15*spb24*t6 + d7*spb23*t14*t6 + t31*t6 + t32*t6 + d5*t9*T(3) + d2*t21*(t17 + s23*T(4)); 
complex<T> t4 = d3*t10*t14*t27 + spa35*spa45*spb15*t29 - d10*t9*T(2) + d9*spb12*t14*T(3) + d7*spb12*(spa45*spb15*spb24 + spb23*t14)*T(3) + d8*t21*T(3) + d2*t16*(t17 + s23*T(4)); 
complex<T> t24 = d19*spa45*(s15*t16 + t14*t27); 
complex<T> t38 = d13*t16 + d14*(t43 + spb12*t5); 
complex<T> t39 = -t51; 
complex<T> t42 = d15*spa45*t16 + d17*t45 + d16*spa45*spb12*t5; 
complex<T> t44 = t31*t6 + t32*t6 + t12*t9; 
complex<T> t33 = d11*s15*t21 + d12*t39 + d11*t27*t5; 
complex<T> co1 = d18*s12*t43; 
complex<T> co2 = d18*s34*t43; 
complex<T> co3 = d20*s34*t45; 
complex<T> co4 = d21*spa45*t51; 
complex<T> co5 = d22*s12*t39; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t44*(*CI_users[2]->get_value(mc,ind,mu)) + t33*(*CI_users[3]->get_value(mc,ind,mu)) + t38*(*CI_users[4]->get_value(mc,ind,mu)) + t42*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t24*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu)) + co3*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)) + co5*(*CI_users[11]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qmqpQmQpgam_sl_wCI::\
C2q2Q1y_qmqpQmQpgam_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmqpQmQpgam_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, gam}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmqpQmQpgam sl");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:07:53 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb45 = SPB(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> t3 = spb34*T(2); 
complex<T> t16 = square(spb24); 
complex<T> t17 = cube(spa15); 
complex<T> t18 = square(spb12); 
complex<T> t19 = square(spb45); 
complex<T> t20 = s12 + s25 - s34; 
complex<T> t21 = -(s25*spa25); 
complex<T> t23 = s15*spb24 + spa13*spb12*spb34*T(2); 
complex<T> t26 = square(spa35); 
complex<T> d1 = (s12 - s34)*spb15*spb34; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*(s12 + s15 - s34)*spb15*spb34; d2 = T(1)/d2;
complex<T> d3 = spb15*square(s12 - s34)*T(2); d3 = T(1)/d3;
complex<T> d6 = spb12*spb35*spb45*T(2); d6 = T(1)/d6;
complex<T> d9 = (s15 - s34)*spb15*spb34; d9 = T(1)/d9;
complex<T> d10 = (s15 - s34)*(s12 + s15 - s34)*spb15*spb34; d10 = T(1)/d10;
complex<T> d15 = spb15*spb34*square(s12 + s15 - s34); d15 = T(1)/d15;
complex<T> d16 = spb25*spb34*spb35; d16 = T(1)/d16;
complex<T> d17 = spb15*spb34*spb45; d17 = T(1)/d17;
complex<T> d18 = spb35*spb45; d18 = T(1)/d18;
complex<T> d20 = spb34*square(s12 + s15 - s34); d20 = T(1)/d20;
complex<T> d22 = spb15*square(s12 + s15 - s34); d22 = T(1)/d22;
complex<T> d23 = spb15*spb25; d23 = T(1)/d23;
complex<T> t37 = t17*t18; 
complex<T> t38 = spa25*t16; 
complex<T> t39 = spb25*t3; 
complex<T> t45 = t16*t21; 
complex<T> d5 = spb12*spb15*spb45*t3; d5 = T(1)/d5;
complex<T> d7 = (s12 - s34)*spb25*spb34*square(t20); d7 = T(1)/d7;
complex<T> d12 = (s25 - s34)*spb25*spb34*square(t20); d12 = T(1)/d12;
complex<T> d19 = spb25*spb34*cube(t20); d19 = T(1)/d19;
complex<T> d21 = spb34*cube(t20); d21 = T(1)/d21;
complex<T> d24 = spb25*cube(t20); d24 = T(1)/d24;
complex<T> d25 = t3*square(s12 + s15 - s34); d25 = T(1)/d25;
complex<T> d26 = t3*cube(t20); d26 = T(1)/d26;
complex<T> t29 = -t37; 
complex<T> t30 = -t38; 
complex<T> t44 = t19*t37; 
complex<T> d4 = spb12*spb35*t39; d4 = T(1)/d4;
complex<T> d8 = t20*t39*square(s12 - s34); d8 = T(1)/d8;
complex<T> d11 = t39*square(s25 - s34); d11 = T(1)/d11;
complex<T> d13 = t20*t39*square(s25 - s34); d13 = T(1)/d13;
complex<T> d14 = spb15*t39; d14 = T(1)/d14;
complex<T> t6 = d11*spa15*spb24*t23 + (d12 + d13)*t44; 
complex<T> t7 = -(d1*spa15*spb14*spb24) - d11*spa15*spb24*t23 - d3*spb25*spb34*t26 + d12*t19*t29 + d13*t19*t29 + d7*t19*t29 + d8*t19*t29 + d9*t38 + d10*s25*t38 + d2*s25*t38 - d14*t16*T(3); 
complex<T> t13 = d1*spa15*spb14*spb24 + d3*spb25*spb34*t26 + d7*t44 + d8*t44 + d2*t45 + d6*t16*T(3) - d5*spb14*t16*T(3) + d4*spb23*t16*T(3); 
complex<T> t14 = d18*spa12*t16 - d17*spa12*spb14*t16 + d16*spa12*spb23*t16 + d19*s12*t19*t29 + d15*s12*s25*t38; 
complex<T> t15 = spa34*(-(d23*t16) + d24*t19*t29 + d22*s25*t38); 
complex<T> t36 = d9*t30 + d10*t45; 
complex<T> co1 = d20*spa15*t45; 
complex<T> co2 = d21*spa25*t44; 
complex<T> co3 = d25*s12*spa15*t45; 
complex<T> co4 = d26*s12*spa25*t44; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t13*(*CI_users[0]->get_value(mc,ind,mu)) + t36*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t7*(*CI_users[3]->get_value(mc,ind,mu)) + t14*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qmqpQpQmgam_sl_wCI::\
C2q2Q1y_qmqpQpQmgam_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmqpQpQmgam_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, gam}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmqpQpQmgam sl");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:07:59 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> t3 = spb34*T(2); 
complex<T> t16 = square(spb23); 
complex<T> t17 = cube(spa15); 
complex<T> t18 = square(spb12); 
complex<T> t19 = square(spb35); 
complex<T> t20 = s12 + s25 - s34; 
complex<T> t21 = -(s25*spa25); 
complex<T> t23 = s15*spb23 - spa14*spb12*spb34*T(2); 
complex<T> t26 = square(spa45); 
complex<T> d1 = (s12 - s34)*spb15*spb34; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*(s12 + s15 - s34)*spb15*spb34; d2 = T(1)/d2;
complex<T> d3 = spb15*square(s12 - s34)*T(2); d3 = T(1)/d3;
complex<T> d8 = spb12*spb35*spb45*T(2); d8 = T(1)/d8;
complex<T> d9 = (s15 - s34)*spb15*spb34; d9 = T(1)/d9;
complex<T> d10 = (s15 - s34)*(s12 + s15 - s34)*spb15*spb34; d10 = T(1)/d10;
complex<T> d15 = spb15*spb34*square(s12 + s15 - s34); d15 = T(1)/d15;
complex<T> d16 = spb15*spb34*spb35; d16 = T(1)/d16;
complex<T> d18 = spb25*spb34*spb45; d18 = T(1)/d18;
complex<T> d19 = spb35*spb45; d19 = T(1)/d19;
complex<T> d20 = spb34*square(s12 + s15 - s34); d20 = T(1)/d20;
complex<T> d22 = spb15*square(s12 + s15 - s34); d22 = T(1)/d22;
complex<T> d23 = spb15*spb25; d23 = T(1)/d23;
complex<T> t37 = t17*t18; 
complex<T> t38 = spa25*t16; 
complex<T> t39 = spb25*t3; 
complex<T> t45 = t16*t21; 
complex<T> d4 = spb12*spb15*spb35*t3; d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*spb25*spb34*square(t20); d5 = T(1)/d5;
complex<T> d12 = (s25 - s34)*spb25*spb34*square(t20); d12 = T(1)/d12;
complex<T> d17 = spb25*spb34*cube(t20); d17 = T(1)/d17;
complex<T> d21 = spb34*cube(t20); d21 = T(1)/d21;
complex<T> d24 = spb25*cube(t20); d24 = T(1)/d24;
complex<T> d25 = t3*square(s12 + s15 - s34); d25 = T(1)/d25;
complex<T> d26 = t3*cube(t20); d26 = T(1)/d26;
complex<T> t29 = -t37; 
complex<T> t30 = -t38; 
complex<T> t44 = t19*t37; 
complex<T> d6 = t20*t39*square(s12 - s34); d6 = T(1)/d6;
complex<T> d7 = spb12*spb45*t39; d7 = T(1)/d7;
complex<T> d11 = t39*square(s25 - s34); d11 = T(1)/d11;
complex<T> d13 = t20*t39*square(s25 - s34); d13 = T(1)/d13;
complex<T> d14 = spb15*t39; d14 = T(1)/d14;
complex<T> t6 = d11*spa15*spb23*t23 + (d12 + d13)*t44; 
complex<T> t7 = -(d1*spa15*spb13*spb23) - d11*spa15*spb23*t23 - d3*spb25*spb34*t26 + d12*t19*t29 + d13*t19*t29 + d5*t19*t29 + d6*t19*t29 + d9*t38 + d10*s25*t38 + d2*s25*t38 - d14*t16*T(3); 
complex<T> t13 = d1*spa15*spb13*spb23 + d3*spb25*spb34*t26 + d5*t44 + d6*t44 + d2*t45 - d8*t16*T(3) - d4*spb13*t16*T(3) + d7*spb24*t16*T(3); 
complex<T> t14 = -(d19*spa12*t16) - d16*spa12*spb13*t16 + d18*spa12*spb24*t16 + d17*s12*t19*t29 + d15*s12*s25*t38; 
complex<T> t15 = spa34*(-(d23*t16) + d24*t19*t29 + d22*s25*t38); 
complex<T> t36 = d9*t30 + d10*t45; 
complex<T> co1 = d20*spa15*t45; 
complex<T> co2 = d21*spa25*t44; 
complex<T> co3 = d25*s12*spa15*t45; 
complex<T> co4 = d26*s12*spa25*t44; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t13*(*CI_users[0]->get_value(mc,ind,mu)) + t36*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t7*(*CI_users[3]->get_value(mc,ind,mu)) + t14*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qpqmQmQpgam_sl_wCI::\
C2q2Q1y_qpqmQmQpgam_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpqmQmQpgam_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, Qm, Qp, gam}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpqmQmQpgam sl");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:08:04 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb45 = SPB(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> t2 = spb34*T(2); 
complex<T> t8 = spb15*spb34; 
complex<T> t14 = cube(spa25); 
complex<T> t15 = square(spb12); 
complex<T> t16 = square(spb14); 
complex<T> t17 = square(spb45); 
complex<T> t18 = s12 + s15 - s34; 
complex<T> t21 = s25*spb14 - spa23*spb12*spb34*T(2); 
complex<T> t23 = square(spa35); 
complex<T> t28 = s12*spa15; 
complex<T> d1 = (s12 - s34)*spb25; d1 = T(1)/d1;
complex<T> d2 = spb25*square(s12 - s34)*T(2); d2 = T(1)/d2;
complex<T> d5 = spb12*spb35*spb45*T(2); d5 = T(1)/d5;
complex<T> d12 = (s12 + s25 - s34)*spb25*spb34; d12 = T(1)/d12;
complex<T> d13 = spb25*spb34*spb35; d13 = T(1)/d13;
complex<T> d15 = spb35*spb45; d15 = T(1)/d15;
complex<T> d18 = (s12 + s25 - s34)*spb34; d18 = T(1)/d18;
complex<T> d19 = (s12 + s25 - s34)*spb25; d19 = T(1)/d19;
complex<T> d20 = spb15*spb25; d20 = T(1)/d20;
complex<T> t1 = square(t18); 
complex<T> t26 = spb15*t2; 
complex<T> t33 = t14*t15; 
complex<T> t49 = spa15*t16; 
complex<T> d3 = spb12*spb25*spb35*t2; d3 = T(1)/d3;
complex<T> d14 = spb45*t8; d14 = T(1)/d14;
complex<T> d16 = t8*cube(t18); d16 = T(1)/d16;
complex<T> d17 = spb34*cube(t18); d17 = T(1)/d17;
complex<T> d21 = spb15*cube(t18); d21 = T(1)/d21;
complex<T> d22 = t2*cube(t18); d22 = T(1)/d22;
complex<T> d23 = (s12 + s25 - s34)*t2; d23 = T(1)/d23;
complex<T> t25 = -t33; 
complex<T> t42 = t17*t33; 
complex<T> d4 = spb12*spb45*t26; d4 = T(1)/d4;
complex<T> d6 = (s12 - s34)*t1*t8; d6 = T(1)/d6;
complex<T> d7 = t18*t26*square(s12 - s34); d7 = T(1)/d7;
complex<T> d8 = t26*square(s15 - s34); d8 = T(1)/d8;
complex<T> d9 = (s15 - s34)*t1*t8; d9 = T(1)/d9;
complex<T> d10 = t18*t26*square(s15 - s34); d10 = T(1)/d10;
complex<T> d11 = spb25*t26; d11 = T(1)/d11;
complex<T> t5 = d1*spa35*spb14 - d8*spa25*spb14*t21 + d10*t17*t25 + d6*t17*t25 + d7*t17*t25 + d9*t17*t25 - d2*t23*t8 - d11*t16*T(3); 
complex<T> t6 = d8*spa25*spb14*t21 + (d10 + d9)*t42; 
complex<T> t11 = -(d1*spa35*spb14) + d2*spb15*spb34*t23 + d6*t42 + d7*t42 + d5*t16*T(3) + d3*spb23*t16*T(3) - d4*cube(spb14)*T(3); 
complex<T> t12 = d15*spa12*t16 + d13*spa12*spb23*t16 + d16*s12*t17*t25 - d12*t16*t28 - d14*spa12*cube(spb14); 
complex<T> t13 = -(spa34*(d20*t16 - d21*t17*t25 + d19*t49)); 
complex<T> co1 = d17*spa15*t42; 
complex<T> co2 = d18*spa25*t49; 
complex<T> co3 = d22*t28*t42; 
complex<T> co4 = d23*spa25*t16*t28; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t11*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + t12*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + t13*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + co4*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qpqmQpQmgam_sl_wCI::\
C2q2Q1y_qpqmQpQmgam_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpqmQpQmgam_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, Qp, Qm, gam}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpqmQpQmgam sl");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:08:09 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> t2 = spb34*T(2); 
complex<T> t8 = spb15*spb34; 
complex<T> t14 = cube(spa25); 
complex<T> t15 = square(spb12); 
complex<T> t16 = square(spb13); 
complex<T> t17 = square(spb35); 
complex<T> t18 = s12 + s15 - s34; 
complex<T> t21 = s25*spb13 + spa24*spb12*spb34*T(2); 
complex<T> t23 = square(spa45); 
complex<T> t28 = s12*spa15; 
complex<T> d1 = (s12 - s34)*spb25; d1 = T(1)/d1;
complex<T> d2 = spb25*square(s12 - s34)*T(2); d2 = T(1)/d2;
complex<T> d7 = spb12*spb35*spb45*T(2); d7 = T(1)/d7;
complex<T> d12 = (s12 + s25 - s34)*spb25*spb34; d12 = T(1)/d12;
complex<T> d15 = spb25*spb34*spb45; d15 = T(1)/d15;
complex<T> d16 = spb35*spb45; d16 = T(1)/d16;
complex<T> d18 = (s12 + s25 - s34)*spb34; d18 = T(1)/d18;
complex<T> d19 = (s12 + s25 - s34)*spb25; d19 = T(1)/d19;
complex<T> d20 = spb15*spb25; d20 = T(1)/d20;
complex<T> t1 = square(t18); 
complex<T> t26 = spb15*t2; 
complex<T> t33 = t14*t15; 
complex<T> t49 = spa15*t16; 
complex<T> d6 = spb12*spb25*spb45*t2; d6 = T(1)/d6;
complex<T> d13 = spb35*t8; d13 = T(1)/d13;
complex<T> d14 = t8*cube(t18); d14 = T(1)/d14;
complex<T> d17 = spb34*cube(t18); d17 = T(1)/d17;
complex<T> d21 = spb15*cube(t18); d21 = T(1)/d21;
complex<T> d22 = t2*cube(t18); d22 = T(1)/d22;
complex<T> d23 = (s12 + s25 - s34)*t2; d23 = T(1)/d23;
complex<T> t25 = -t33; 
complex<T> t42 = t17*t33; 
complex<T> d3 = spb12*spb35*t26; d3 = T(1)/d3;
complex<T> d4 = (s12 - s34)*t1*t8; d4 = T(1)/d4;
complex<T> d5 = t18*t26*square(s12 - s34); d5 = T(1)/d5;
complex<T> d8 = t26*square(s15 - s34); d8 = T(1)/d8;
complex<T> d9 = (s15 - s34)*t1*t8; d9 = T(1)/d9;
complex<T> d10 = t18*t26*square(s15 - s34); d10 = T(1)/d10;
complex<T> d11 = spb25*t26; d11 = T(1)/d11;
complex<T> t5 = -(d1*spa45*spb13) - d8*spa25*spb13*t21 + d10*t17*t25 + d4*t17*t25 + d5*t17*t25 + d9*t17*t25 - d2*t23*t8 - d11*t16*T(3); 
complex<T> t6 = d8*spa25*spb13*t21 + (d10 + d9)*t42; 
complex<T> t11 = d1*spa45*spb13 + d2*spb15*spb34*t23 + d4*t42 + d5*t42 - d7*t16*T(3) + d6*spb24*t16*T(3) - d3*cube(spb13)*T(3); 
complex<T> t12 = -(d16*spa12*t16) + d15*spa12*spb24*t16 + d14*s12*t17*t25 - d12*t16*t28 - d13*spa12*cube(spb13); 
complex<T> t13 = -(spa34*(d20*t16 - d21*t17*t25 + d19*t49)); 
complex<T> co1 = d17*spa15*t42; 
complex<T> co2 = d18*spa25*t49; 
complex<T> co3 = d22*t28*t42; 
complex<T> co4 = d23*spa25*t16*t28; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t11*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + t12*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + t13*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + co4*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qpQmQpqmgap_nf_wCI::\
C2q2Q1y_qpQmQpqmgap_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpQmQpqmgap_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qm, Qp, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQmQpqmgap nf");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:08:10 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa24); 
complex<T> d1 = spa15*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1y_qpQpQmqmgap_nf_wCI::\
C2q2Q1y_qpQpQmqmgap_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpQpQmqmgap_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qp, Qm, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQpQmqmgap nf");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:08:10 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa34); 
complex<T> d1 = spa15*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1y_qmQpQmqpgap_nf_wCI::\
C2q2Q1y_qmQpQmqpgap_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmQpQmqpgap_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQpQmqpgap nf");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:08:11 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa15*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1y_qmQmQpqpgap_nf_wCI::\
C2q2Q1y_qmQmQpqpgap_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmQmQpqpgap_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQmQpqpgap nf");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:08:11 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa12); 
complex<T> d1 = spa15*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1y_qpQmQpqmgap_L_wCI::\
C2q2Q1y_qpQmQpqmgap_L_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
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
     C2q2Q1y_qpQmQpqmgap_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qm, Qp, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQmQpqmgap L");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:11:06 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t120*(*CI_users[1]->get_value(mc,ind,mu)) + t12*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t11*(*CI_users[4]->get_value(mc,ind,mu)) + t8*(*CI_users[5]->get_value(mc,ind,mu)) + t88*(*CI_users[6]->get_value(mc,ind,mu)) + t9*(*CI_users[7]->get_value(mc,ind,mu)) + t22*(*CI_users[8]->get_value(mc,ind,mu)) + t10*(*CI_users[9]->get_value(mc,ind,mu)) + t156*(*CI_users[10]->get_value(mc,ind,mu)) + co1*(*CI_users[11]->get_value(mc,ind,mu)) + t160*(*CI_users[12]->get_value(mc,ind,mu)) + co2*(*CI_users[13]->get_value(mc,ind,mu)) + t142*(*CI_users[14]->get_value(mc,ind,mu)) + t151*(*CI_users[15]->get_value(mc,ind,mu)) + t135*(*CI_users[16]->get_value(mc,ind,mu)) + t23*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qpQpQmqmgap_L_wCI::\
C2q2Q1y_qpQpQmqmgap_L_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpQpQmqmgap_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qp, Qm, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpQpQmqmgap L");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:11:11 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co5*(t42*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t4*(*CI_users[2]->get_value(mc,ind,mu)) + t38*(*CI_users[3]->get_value(mc,ind,mu)) + t26*(*CI_users[4]->get_value(mc,ind,mu)) + t33*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t45*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qmQpQmqpgap_L_wCI::\
C2q2Q1y_qmQpQmqpgap_L_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmQpQmqpgap_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQpQmqpgap L");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:14:30 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t12*(*CI_users[1]->get_value(mc,ind,mu)) + t13*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t158*(*CI_users[4]->get_value(mc,ind,mu)) + t27*(*CI_users[5]->get_value(mc,ind,mu)) + t9*(*CI_users[6]->get_value(mc,ind,mu)) + t8*(*CI_users[7]->get_value(mc,ind,mu)) + t11*(*CI_users[8]->get_value(mc,ind,mu)) + t122*(*CI_users[9]->get_value(mc,ind,mu)) + t131*(*CI_users[10]->get_value(mc,ind,mu)) + t162*(*CI_users[11]->get_value(mc,ind,mu)) + t10*(*CI_users[12]->get_value(mc,ind,mu)) + t149*(*CI_users[13]->get_value(mc,ind,mu)) + co1*(*CI_users[14]->get_value(mc,ind,mu)) + t139*(*CI_users[15]->get_value(mc,ind,mu)) + t28*(*CI_users[16]->get_value(mc,ind,mu)) + co2*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qmQmQpqpgap_L_wCI::\
C2q2Q1y_qmQmQpqpgap_L_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmQmQpqpgap_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp, gap}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmQmQpqpgap L");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:14:35 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co6*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t45*(*CI_users[2]->get_value(mc,ind,mu)) + t32*(*CI_users[3]->get_value(mc,ind,mu)) + t37*(*CI_users[4]->get_value(mc,ind,mu)) + t42*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t24*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu)) + co3*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)) + co5*(*CI_users[11]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qpqmQpQmgap_sl_wCI::\
C2q2Q1y_qpqmQpQmgap_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpqmQpQmgap_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, Qp, Qm, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpqmQpQmgap sl");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:14:41 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co5*(t15*(*CI_users[0]->get_value(mc,ind,mu)) + t36*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t7*(*CI_users[3]->get_value(mc,ind,mu)) + t13*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t14*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qpqmQmQpgap_sl_wCI::\
C2q2Q1y_qpqmQmQpgap_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qpqmQmQpgap_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, Qm, Qp, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qpqmQmQpgap sl");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:14:48 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co5*(t15*(*CI_users[0]->get_value(mc,ind,mu)) + t36*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t7*(*CI_users[3]->get_value(mc,ind,mu)) + t13*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t14*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qmqpQpQmgap_sl_wCI::\
C2q2Q1y_qmqpQpQmgap_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmqpQpQmgap_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmqpQpQmgap sl");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:14:52 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co5*(t13*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t11*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + t12*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + co4*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1y_qmqpQmQpgap_sl_wCI::\
C2q2Q1y_qmqpQmQpgap_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c1, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1y_qmqpQmQpgap_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1y :  qmqpQmQpgap sl");
#endif
 
//#define TimeStamp "Sun 6 Feb 2011 19:14:57 on n2148"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co5*(t13*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t11*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + t12*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + co4*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qmQpQmqpgam_nf C2q2Q1y_607_nf
#define _C_qmQmQpqpgam_nf C2q2Q1y_637_nf
#define _C_qpQmQpqmgam_nf C2q2Q1y_422_nf
#define _C_qpQpQmqmgam_nf C2q2Q1y_392_nf
#define _C_qmQpQmqpgam_L C2q2Q1y_607_L
#define _C_qmQmQpqpgam_L C2q2Q1y_637_L
#define _C_qpQmQpqmgam_L C2q2Q1y_422_L
#define _C_qpQpQmqmgam_L C2q2Q1y_392_L
#define _C_qmqpQmQpgam_sl C2q2Q1y_1237_sl
#define _C_qmqpQpQmgam_sl C2q2Q1y_1057_sl
#define _C_qpqmQmQpgam_sl C2q2Q1y_1232_sl
#define _C_qpqmQpQmgam_sl C2q2Q1y_1052_sl
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
 
#define _CASE_qmQpQmqpgam_nf case 607
 
#define _CASE_qmQmQpqpgam_nf case 637
 
#define _CASE_qpQmQpqmgam_nf case 422
 
#define _CASE_qpQpQmqmgam_nf case 392
 
#define _CASE_qmQpQmqpgam_L case 607
 
#define _CASE_qmQmQpqpgam_L case 637
 
#define _CASE_qpQmQpqmgam_L case 422
 
#define _CASE_qpQpQmqmgam_L case 392
 
#define _CASE_qmqpQmQpgam_sl case 1237
 
#define _CASE_qmqpQpQmgam_sl case 1057
 
#define _CASE_qpqmQmQpgam_sl case 1232
 
#define _CASE_qpqmQpQmgam_sl case 1052
 
#define _CASE_qpQmQpqmgap_nf case 4310
 
#define _CASE_qpQpQmqmgap_nf case 4280
 
#define _CASE_qmQpQmqpgap_nf case 4495
 
#define _CASE_qmQmQpqpgap_nf case 4525
 
#define _CASE_qpQmQpqmgap_L case 4310
 
#define _CASE_qpQpQmqmgap_L case 4280
 
#define _CASE_qmQpQmqpgap_L case 4495
 
#define _CASE_qmQmQpqpgap_L case 4525
 
#define _CASE_qpqmQpQmgap_sl case 4940
 
#define _CASE_qpqmQmQpgap_sl case 5120
 
#define _CASE_qmqpQpQmgap_sl case 4945
 
#define _CASE_qmqpQmQpgap_sl case 5125
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2Q1y_L( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmQpQmqpgam_L: return new 
                       C2q2Q1y_qmQpQmqpgam_L_wCI(ind);
    _CASE_qmQmQpqpgam_L: return new 
                       C2q2Q1y_qmQmQpqpgam_L_wCI(ind);
    _CASE_qpQmQpqmgam_L: return new 
                       C2q2Q1y_qpQmQpqmgam_L_wCI(ind);
    _CASE_qpQpQmqmgam_L: return new 
                       C2q2Q1y_qpQpQmqmgam_L_wCI(ind);
    _CASE_qpQmQpqmgap_L: return new 
                       C2q2Q1y_qpQmQpqmgap_L_wCI(ind);
    _CASE_qpQpQmqmgap_L: return new 
                       C2q2Q1y_qpQpQmqmgap_L_wCI(ind);
    _CASE_qmQpQmqpgap_L: return new 
                       C2q2Q1y_qmQpQmqpgap_L_wCI(ind);
    _CASE_qmQmQpqpgap_L: return new 
                       C2q2Q1y_qmQmQpqpgap_L_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q1y_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmQpQmqpgam_nf: return new 
                       C2q2Q1y_qmQpQmqpgam_nf_wCI(ind);
    _CASE_qmQmQpqpgam_nf: return new 
                       C2q2Q1y_qmQmQpqpgam_nf_wCI(ind);
    _CASE_qpQmQpqmgam_nf: return new 
                       C2q2Q1y_qpQmQpqmgam_nf_wCI(ind);
    _CASE_qpQpQmqmgam_nf: return new 
                       C2q2Q1y_qpQpQmqmgam_nf_wCI(ind);
    _CASE_qpQmQpqmgap_nf: return new 
                       C2q2Q1y_qpQmQpqmgap_nf_wCI(ind);
    _CASE_qpQpQmqmgap_nf: return new 
                       C2q2Q1y_qpQpQmqmgap_nf_wCI(ind);
    _CASE_qmQpQmqpgap_nf: return new 
                       C2q2Q1y_qmQpQmqpgap_nf_wCI(ind);
    _CASE_qmQmQpqpgap_nf: return new 
                       C2q2Q1y_qmQmQpqpgap_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q1y_sl( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmqpQmQpgam_sl: return new 
                       C2q2Q1y_qmqpQmQpgam_sl_wCI(ind);
    _CASE_qmqpQpQmgam_sl: return new 
                       C2q2Q1y_qmqpQpQmgam_sl_wCI(ind);
    _CASE_qpqmQmQpgam_sl: return new 
                       C2q2Q1y_qpqmQmQpgam_sl_wCI(ind);
    _CASE_qpqmQpQmgam_sl: return new 
                       C2q2Q1y_qpqmQpQmgam_sl_wCI(ind);
    _CASE_qpqmQpQmgap_sl: return new 
                       C2q2Q1y_qpqmQpQmgap_sl_wCI(ind);
    _CASE_qpqmQmQpgap_sl: return new 
                       C2q2Q1y_qpqmQmQpgap_sl_wCI(ind);
    _CASE_qmqpQpQmgap_sl: return new 
                       C2q2Q1y_qmqpQpQmgap_sl_wCI(ind);
    _CASE_qmqpQmQpgap_sl: return new 
                       C2q2Q1y_qmqpQmQpgap_sl_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
