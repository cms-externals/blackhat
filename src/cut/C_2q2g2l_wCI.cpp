/*
*C_2q2g2l_wCI.cpp
*
* Created on 9/25, 2009
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
 



 
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpppqmemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qppmqmemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpmpqmemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpmmqmemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qmppqpemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qmpmqpemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qmmpqpemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qmmmqpemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qppqmpemep_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qppqmmemep_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpmqmpemep_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpmqmmemep_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmppemep_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmpmemep_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmmpemep_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmmmemep_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpppqmemep_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qppmqmemep_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpmpqmemep_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpmmqmemep_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpppqmemep_nf_top_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qppmqmemep_nf_top_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpmpqmemep_nf_top_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpmmqmemep_nf_top_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmppemep_VECT_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmpmemep_VECT_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmmpemep_VECT_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmmmemep_VECT_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmppemep_AX_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmpmemep_AX_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmmpemep_AX_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmmmemep_AX_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmppemep_AXSL_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmpmemep_AXSL_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmmpemep_AXSL_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2g2l_qpqmmmemep_AXSL_wCI)


C2q2g2l_qpppqmemep_L_wCI::\
C2q2g2l_qpppqmemep_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c4, c23, c1, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpppqmemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, p, p, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpppqmemep L");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:21:06 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb56 = SPB(5,6);
complex<T> spa25 = SPA(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb13 = SPB(1,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s56 = -(spa56*spb56);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> t4 = spa56*T(2); 
complex<T> t7 = spa12*spa23; 
complex<T> t11 = square(spa45); 
complex<T> t14 = square(spa25); 
complex<T> t15 = square(spb23); 
complex<T> t16 = s123*s234 - s23*s56; 
complex<T> t22 = spa24*spb12; 
complex<T> t3 = square(spa15*spa34*spb13 + spa15*t22); 
complex<T> d1 = (-s234 + s56)*spa34*spa56*t7; d1 = T(1)/d1;
complex<T> d2 = spa34*t4*t7*square(s234 - s56); d2 = T(1)/d2;
complex<T> d3 = (-s234 + s34)*spa56*t7; d3 = T(1)/d3;
complex<T> d4 = t4*t7*square(s234 - s34); d4 = T(1)/d4;
complex<T> d5 = spa34*t4*t7; d5 = T(1)/d5;
complex<T> d6 = spa34*t7; d6 = T(1)/d6;
complex<T> d7 = spa34*t4; d7 = T(1)/d7;
complex<T> d8 = spa12*t4; d8 = T(1)/d8;
complex<T> t9 = d1*(spa34*spb13 + t22); 
complex<T> t20 = d3*spa25; 
complex<T> t28 = d5*t11; 
complex<T> t1 = d2*t3 - spa15*spa45*t9*T(2) + t28*T(3); 
complex<T> t2 = d4*spa34*t14*t15 - d2*t3 + spa45*(spb23*t20 + spa15*t9)*T(2); 
complex<T> t10 = -(d4*spa34*t14*t15) - spa45*spb23*t20*T(2); 
complex<T> co1 = d6*spb56*t11; 
complex<T> co2 = d7*spb12*spb23*t11; 
complex<T> co3 = t16*t28; 
complex<T> co4 = d8*spb23*spb34*t11; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t10*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qppmqmemep_L_wCI::\
C2q2g2l_qppmqmemep_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c56, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c56, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qppmqmemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, p, m, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qppmqmemep L");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:21:34 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb16 = SPB(1,6);
complex<T> spa13 = SPA(1,3);
complex<T> spb26 = SPB(2,6);
complex<T> spb56 = SPB(5,6);
complex<T> spb46 = SPB(4,6);
complex<T> spb36 = SPB(3,6);
complex<T> spa25 = SPA(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s123 = SS(1,2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s56 = -(spa56*spb56);
complex<T> s234 = SS(2,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> t9 = -(s12*s34*square(spa45)) + square(spa14*spa35*spb13 + spa14*spa45*spb14 + spa24*spa35*spb23 + spa24*spa45*spb24); 
complex<T> t21 = square(spa34); 
complex<T> t22 = square(spb12); 
complex<T> t24 = spa12*(spa12*spb24 + spa13*spb34); 
complex<T> t25 = square(spa35*spb23 + spa45*spb24); 
complex<T> t26 = square(spa13*spb16 + spa23*spb26); 
complex<T> t27 = square(spa25*spb24 + spa35*spb34); 
complex<T> t28 = square(spa12*spb26 + spa13*spb36); 
complex<T> t37 = spb24*(spa35*spb23 + spa45*spb24); 
complex<T> t39 = spa56*(spa12*spb24 + spa13*spb34); 
complex<T> t43 = s234*spb34; 
complex<T> t44 = -(spa13*spb16) - spa23*spb26; 
complex<T> t46 = spa24*spb12 + spa34*spb13; 
complex<T> t52 = square(spa15); 
complex<T> t54 = spa14*spb13 + spa24*spb23; 
complex<T> t57 = square(spb46); 
complex<T> t72 = s23*s56; 
complex<T> t73 = spb34*T(2); 
complex<T> d12 = T(2); d12 = T(1)/d12;
complex<T> t12 = -(s12*s34*square(spb16)) + square(spa14*spb14*spb16 + spa24*spb14*spb26 - spb13*t44); 
complex<T> t40 = spb56*t46; 
complex<T> t69 = spb56*t24; 
complex<T> t70 = spa13*t44; 
complex<T> t74 = spa56*t46; 
complex<T> t86 = t22*t28; 
complex<T> t91 = t21*t27; 
complex<T> t100 = t22*t52; 
complex<T> t101 = t21*t57; 
complex<T> d5 = (s23 - s234)*t39*t43; d5 = T(1)/d5;
complex<T> d6 = t39*t43*square(s23 - s234)*T(2); d6 = T(1)/d6;
complex<T> d7 = spb23*t39*t73*square(s234 - s56); d7 = T(1)/d7;
complex<T> d8 = (-s234 + s56)*spb23*spb34*t39; d8 = T(1)/d8;
complex<T> d9 = s123*spa12*t72; d9 = T(1)/d9;
complex<T> d10 = t43*t72; d10 = T(1)/d10;
complex<T> d11 = spa12*spb34*t72; d11 = T(1)/d11;
complex<T> d13 = spb23*t39*t43; d13 = T(1)/d13;
complex<T> t78 = d8*spa15; 
complex<T> t82 = d13*spb24; 
complex<T> t92 = s234*t40; 
complex<T> d1 = s123*t69*square(s123 - s23)*T(2); d1 = T(1)/d1;
complex<T> d2 = s123*(-s123 + s23)*t69; d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spa23*t69; d3 = T(1)/d3;
complex<T> d4 = spa23*t69*square(s123 - s56)*T(2); d4 = T(1)/d4;
complex<T> d15 = s123*spb13*spb23*t74; d15 = T(1)/d15;
complex<T> d16 = s123*spa23*t69; d16 = T(1)/d16;
complex<T> d17 = s123*spb23*t54*t74; d17 = T(1)/d17;
complex<T> d18 = spa12*t54*t73*t74; d18 = T(1)/d18;
complex<T> d20 = spa12*(spa23*spb13 + spa24*spb14)*t40*t73; d20 = T(1)/d20;
complex<T> t18 = -(d7*s234*spb24*t100) + d6*spb23*spb24*t91 + d5*(s34*spa35 - spa25*spa34*spb24)*t37*T(2) - spb12*t37*t78*T(2); 
complex<T> t19 = -(d1*spa13*spa23*t86) - d6*spb23*spb24*t91 - d5*s34*spa35*t37*T(2) + d5*spa25*spa34*spb24*t37*T(2) + d2*s12*spb26*t70*T(2) - d2*spa13*spb12*spb36*t70*T(2); 
complex<T> t80 = d3*spa34; 
complex<T> t96 = d16*spa13; 
complex<T> d14 = spa23*(spa23*spb13 + spa24*spb14)*t92; d14 = T(1)/d14;
complex<T> d19 = spa23*spa24*t92; d19 = T(1)/d19;
complex<T> t8 = s23*(t25*t82 + d19*cube(spa34)*square(spb16)); 
complex<T> t10 = s12*(t26*t96 + d15*cube(spb12)*square(spa45)); 
complex<T> t13 = d7*s234*spb24*t100 + d4*s123*spa13*t101 + spb12*t37*t78*T(2) - spb46*t70*t80*T(2) - d10*d12*spa34*spb16*t37*T(3) + d11*d12*spa35*spb23*t44*T(3) + d11*d12*spa45*spb24*t44*T(3) + d12*d9*spa45*spb12*t70*T(3); 
complex<T> t14 = t25*t82 + d14*(-(spa23*spb12) + spa34*spb14)*t21*square(spb16); 
complex<T> t15 = t26*t96 + d17*(spa14*spb12 - spa34*spb23)*t22*square(spa45); 
complex<T> t20 = -(d4*s123*spa13*t101) + d1*spa13*spa23*t86 - d2*s12*spb26*t70*T(2) + d2*spa13*spb12*spb36*t70*T(2) + spb46*t70*t80*T(2); 
complex<T> t1 = (-s156 + s34)*(-t14 + t25*t82 + d19*cube(spa34)*square(spb16)); 
complex<T> t3 = s23*(t15 + d10*spa34*spb16*t37 - d11*spa35*spb23*t44 - d11*spa45*spb24*t44 - d9*spa45*spb12*t70) + t8; 
complex<T> t4 = (-s123 + s23)*(-t15 + t26*t96 + d15*cube(spb12)*square(spa45)); 
complex<T> t11 = s34*t15; 
complex<T> t71 = d12*t14; 
complex<T> t6 = -(s56*(d12*t15 + t71)); 
complex<T> t90 = d12*t11; 
complex<T> t108 = s12*t71; 
complex<T> t2 = t10 + t108 - s12*(d12*t15 - d10*spa34*spb16*t37 + d11*spa35*spb23*t44 + d11*spa45*spb24*t44 + d9*spa45*spb12*t70); 
complex<T> t5 = d10*s34*spa34*spb16*t37 - d11*s34*(spa35*spb23 + spa45*spb24)*t44 - d9*s34*spa45*spb12*t70 + s34*t71 + t90; 
complex<T> co1 = -(d18*spb12*t9); 
complex<T> co2 = -(d20*spa34*t12); 
complex<T> co3 = s234*t108; 
complex<T> co4 = d12*s34*t8; 
complex<T> co5 = d12*s23*t10; 
complex<T> co6 = s123*t90; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t20*(*CI_users[0]->get_value(mc,ind,mu)) + t19*(*CI_users[1]->get_value(mc,ind,mu)) + t18*(*CI_users[2]->get_value(mc,ind,mu)) + t13*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t3*(*CI_users[7]->get_value(mc,ind,mu)) + t1*(*CI_users[8]->get_value(mc,ind,mu)) + t5*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t6*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)) + co5*(*CI_users[14]->get_value(mc,ind,mu)) + co6*(*CI_users[15]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpmpqmemep_L_wCI::\
C2q2g2l_qpmpqmemep_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpmpqmemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, p, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpmpqmemep L");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:34:16 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb16 = SPB(1,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb36 = SPB(3,6);
complex<T> spb14 = SPB(1,4);
complex<T> spb46 = SPB(4,6);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa14 = SPA(1,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s24 = -(spa24*spb24);
complex<T> s134 = SS(1,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s124 = SS(1,2,4);
complex<T> t8 = -(spa15*spb34); 
complex<T> t21 = spa45*spb13; 
complex<T> t23 = spa12*spb46; 
complex<T> t46 = spa23*spb34; 
complex<T> t47 = square(spb46); 
complex<T> t48 = square(spa15); 
complex<T> t49 = spa12*spb26; 
complex<T> t56 = spa35*spb36; 
complex<T> t57 = spa56*T(2); 
complex<T> t68 = square(spa24); 
complex<T> t69 = square(spb13); 
complex<T> t70 = square(spa25*spb24 + spa35*spb34); 
complex<T> t75 = square(spa35*spb23 + spa45*spb24); 
complex<T> t76 = square(spa13*spb16 + spa23*spb26); 
complex<T> t85 = square(spa12); 
complex<T> t86 = (-s234 + s34)*spa56; 
complex<T> t87 = square(spb34); 
complex<T> t89 = s123*T(2); 
complex<T> t99 = s234*spb23; 
complex<T> t123 = -(spa12*spb16); 
complex<T> t124 = spa45*spb34; 
complex<T> t129 = spa25*spb23*T(2) - spa45*spb34*T(2); 
complex<T> t131 = spb13*(-(spa12*spb16) + spa23*spb36); 
complex<T> t133 = s123*(-(spa12*spb14) + spa23*spb34); 
complex<T> t135 = square(s234); 
complex<T> t140 = (s12 - s34 + s56)*spa34; 
complex<T> t141 = (s12 - s34 - s56)*spb12; 
complex<T> t146 = -s12 - s34 + s56; 
complex<T> t148 = spa12*spb23 - spa14*spb34; 
complex<T> t151 = -(spa23*spb13) - spa24*spb14; 
complex<T> t152 = -(spa14*spb13) - spa24*spb23; 
complex<T> t156 = spa24*spb12 + spa34*spb13; 
complex<T> t166 = cube(spa24); 
complex<T> t167 = square(spa45); 
complex<T> t168 = square(spb16); 
complex<T> t171 = spa12*spb26 + spa13*spb36; 
complex<T> t182 = -(spa25*spb24) - spa35*spb34; 
complex<T> t212 = spa12*spb24; 
complex<T> t216 = cube(spb13); 
complex<T> t224 = s12*s56; 
complex<T> t230 = s23*spa34; 
complex<T> t259 = spa13*spb34; 
complex<T> t268 = s34*s56; 
complex<T> t286 = s14*spa45; 
complex<T> t291 = spa14*spb24; 
complex<T> t316 = spa24*spb34; 
complex<T> t342 = s234*spb24; 
complex<T> t350 = spa12*spb34; 
complex<T> t355 = spa15*spb46; 
complex<T> d3 = (s12 - s123)*s123*spa13*(spa13*spb14 + spa23*spb24)*spb56; d3 = T(1)/d3;
complex<T> d10 = spa56; d10 = T(1)/d10;
complex<T> d11 = spa34*spb56; d11 = T(1)/d11;
complex<T> d14 = spb56*square(spa13*spb14 + spa23*spb24)*T(2); d14 = T(1)/d14;
complex<T> d15 = spa13*spb14 + spa23*spb24; d15 = T(1)/d15;
complex<T> d37 = spa56*spb12; d37 = T(1)/d37;
complex<T> d38 = spb56; d38 = T(1)/d38;
complex<T> d39 = spb12*spb56*T(2); d39 = T(1)/d39;
complex<T> d44 = s123*s23*s56*spb12; d44 = T(1)/d44;
complex<T> d45 = T(2); d45 = T(1)/d45;
complex<T> d57 = s234; d57 = T(1)/d57;
complex<T> d61 = square(s123)*T(2); d61 = T(1)/d61;
complex<T> d62 = spa34*spb12*T(2); d62 = T(1)/d62;
complex<T> d63 = spa34; d63 = T(1)/d63;
complex<T> d65 = square(spa13*spb14 + spa23*spb24); d65 = T(1)/d65;
complex<T> d75 = s123; d75 = T(1)/d75;
complex<T> d76 = spb12; d76 = T(1)/d76;
complex<T> t14 = d38*spa24*(-(spa14*spb16*spb34) - spa24*spb26*spb34 + spa13*spb13*spb36 + spa23*spb23*spb36 + spa14*spb13*spb46 + spa24*spb23*spb46)*(s134*spb46 - spa35*spb34*spb56) + d37*((spa13*spb14 + spa23*spb24)*(s12*t124 - s34*t124 + s56*t124 - spa56*spb36*t146) + (s14 - s23)*(s12*spa35*spb34 - s34*spa35*spb34 + s56*spa35*spb34 + spa56*spb46*t146))*t21 + (s14 - s23)*t316*(spa45*spb46 + t56)*T(2); 
complex<T> t50 = -t133; 
complex<T> t51 = -(s234*t148); 
complex<T> t71 = spb56*(t212 + t259); 
complex<T> t72 = square(t171); 
complex<T> t83 = -(spa35*spb13*spb56) - spb46*t156; 
complex<T> t114 = t146*t355 + spa35*spb26*t350*T(2); 
complex<T> t125 = spa56*(t212 + t259); 
complex<T> t126 = spb23*t70; 
complex<T> t127 = t151*t152; 
complex<T> t128 = spa15*t148; 
complex<T> t132 = -(spa23*t47); 
complex<T> t136 = d45*(s14 - s23); 
complex<T> t138 = -(d62*t146); 
complex<T> t145 = -(d3*(spa13*spb16 + spa23*spb26)); 
complex<T> t147 = square(s12) + square(s34) + square(s56) - s12*s34*T(2) - t224*T(2) - t268*T(2); 
complex<T> t154 = -(spa12*spb14) + t46; 
complex<T> t158 = -(spa25*spb23) + t124; 
complex<T> t160 = spa23*spb36 + t123; 
complex<T> t184 = -(spa25*spb26) - t56; 
complex<T> t215 = -(spb34*t68); 
complex<T> t231 = -(d44*spa45); 
complex<T> t244 = s234*(spa13*spb23 + t291); 
complex<T> t261 = spa12*t69; 
complex<T> t282 = (t212 + t259)*t57; 
complex<T> t287 = t133*square(spb46); 
complex<T> t301 = spa34*t156; 
complex<T> t327 = (t212 + t259)*t286; 
complex<T> t348 = spa56*t156; 
complex<T> t361 = spa13*t89; 
complex<T> d5 = spa34*t57; d5 = T(1)/d5;
complex<T> d6 = spa13*spb23 + t291; d6 = T(1)/d6;
complex<T> d9 = (t212 + t259)*(spa13*spb23 + t291); d9 = T(1)/d9;
complex<T> d17 = t212 + t259; d17 = T(1)/d17;
complex<T> d28 = (t212 + t259)*t86*square(spa13*spb23 + t291); d28 = T(1)/d28;
complex<T> d31 = (s23 - s234 + s34)*(t212 + t259)*t342*t86; d31 = T(1)/d31;
complex<T> d34 = t57*square(spa13*spb23 + t291); d34 = T(1)/d34;
complex<T> d41 = (spa13*spb14 + spa23*spb24)*(t212 + t259); d41 = T(1)/d41;
complex<T> d42 = s56*spb12*t230; d42 = T(1)/d42;
complex<T> d43 = s234*s56*t230; d43 = T(1)/d43;
complex<T> d58 = s123*(spa13*spb14 + spa23*spb24)*(t212 + t259); d58 = T(1)/d58;
complex<T> d64 = (spa13*spb14 + spa23*spb24)*(t212 + t259)*(spa13*spb23 + t291)*T(2); d64 = T(1)/d64;
complex<T> d74 = t135*T(2); d74 = T(1)/d74;
complex<T> d77 = square(spa13*spb23 + t291); d77 = T(1)/d77;
complex<T> d79 = (t212 + t259)*cube(spa13*spb14 + spa23*spb24)*T(2); d79 = T(1)/d79;
complex<T> d80 = (t212 + t259)*cube(spa13*spb23 + t291)*T(2); d80 = T(1)/d80;
complex<T> t185 = s123*t71; 
complex<T> t218 = -(t127*T(3)); 
complex<T> t221 = t47*t50; 
complex<T> t234 = d43*t158; 
complex<T> t262 = spa23*t72; 
complex<T> t263 = t126*t68; 
complex<T> t264 = t48*t51; 
complex<T> t267 = t158*T(2); 
complex<T> t278 = t136*t184; 
complex<T> t288 = spb16*t138; 
complex<T> t388 = d28*t99; 
complex<T> d1 = (spa13*spb14 + spa23*spb24)*spb56*t361*square(s12 - s123); d1 = T(1)/d1;
complex<T> d4 = (s12 - s123)*t71*square(spa13*spb14 + spa23*spb24); d4 = T(1)/d4;
complex<T> d7 = spa34*t147*t57; d7 = T(1)/d7;
complex<T> d8 = square(t147); d8 = T(1)/d8;
complex<T> d12 = (spa13*spb14 + spa23*spb24)*t147*(spa13*spb23 + t291); d12 = T(1)/d12;
complex<T> d13 = (spa13*spb14 + spa23*spb24)*square(t147); d13 = T(1)/d13;
complex<T> d16 = spa34*(spa13*spb14 + spa23*spb24)*spb56*t147*T(2); d16 = T(1)/d16;
complex<T> d18 = t361*t71*square(s123 - s23); d18 = T(1)/d18;
complex<T> d21 = (-s123 + s56)*spa23*(spa13*spb14 + spa23*spb24)*t71; d21 = T(1)/d21;
complex<T> d22 = (-s123 + s56)*t71*square(spa13*spb14 + spa23*spb24); d22 = T(1)/d22;
complex<T> d23 = spa23*(spa13*spb14 + spa23*spb24)*t71*square(s123 - s56)*T(2); d23 = T(1)/d23;
complex<T> d24 = t282*t342*square(s23 - s234); d24 = T(1)/d24;
complex<T> d25 = (s23 - s234)*(s23 - s234 + s34)*t125*t342; d25 = T(1)/d25;
complex<T> d26 = (s23 - s234)*t125*t342; d26 = T(1)/d26;
complex<T> d27 = spb24*t244*t57*square(s234 - s34); d27 = T(1)/d27;
complex<T> d29 = (-s234 + s56)*t125*square(spa13*spb23 + t291); d29 = T(1)/d29;
complex<T> d30 = spb23*t282*(spa13*spb23 + t291)*square(s234 - s56); d30 = T(1)/d30;
complex<T> d32 = spb24*t244*t86; d32 = T(1)/d32;
complex<T> d33 = (-s234 + s56)*spb23*t125*(spa13*spb23 + t291); d33 = T(1)/d33;
complex<T> d35 = spb12*t147*(spa13*spb23 + t291)*t57; d35 = T(1)/d35;
complex<T> d36 = (spa13*spb23 + t291)*square(t147); d36 = T(1)/d36;
complex<T> d40 = spb12*spb56*t147*T(2); d40 = T(1)/d40;
complex<T> d46 = t282*cube(spa13*spb23 + t291); d46 = T(1)/d46;
complex<T> d47 = t125*cube(spa13*spb23 + t291); d47 = T(1)/d47;
complex<T> d48 = t125*t99*cube(spa13*spb23 + t291); d48 = T(1)/d48;
complex<T> d49 = s234*spa23*spb56*t301; d49 = T(1)/d49;
complex<T> d50 = s123*spb12*spb23*t348; d50 = T(1)/d50;
complex<T> d54 = t71*cube(spa13*spb14 + spa23*spb24); d54 = T(1)/d54;
complex<T> d56 = t71*cube(spa13*spb14 + spa23*spb24)*T(2); d56 = T(1)/d56;
complex<T> d59 = spb12*spb23*t348; d59 = T(1)/d59;
complex<T> d60 = spa23*(spa13*spb14 + spa23*spb24)*t71; d60 = T(1)/d60;
complex<T> d66 = (spa13*spb14 + spa23*spb24)*t147; d66 = T(1)/d66;
complex<T> d67 = t147*(t212 + t259); d67 = T(1)/d67;
complex<T> d68 = t125*t342*square(s23 - s234 + s34); d68 = T(1)/d68;
complex<T> d69 = t125*t99*cube(spb24); d69 = T(1)/d69;
complex<T> d70 = s234*t125*cube(spb24); d70 = T(1)/d70;
complex<T> d71 = t244*(t212 + t259); d71 = T(1)/d71;
complex<T> d72 = spb23*t125*(spa13*spb23 + t291); d72 = T(1)/d72;
complex<T> d73 = spa23*spb56*t301; d73 = T(1)/d73;
complex<T> d78 = t147*(spa13*spb23 + t291); d78 = T(1)/d78;
complex<T> d81 = t282*t342*square(s23 - s234 + s34); d81 = T(1)/d81;
complex<T> d82 = t361*t71*square(s12 - s123 + s23); d82 = T(1)/d82;
complex<T> t12 = d12*(d10*spb13*(spa15*spa23*spb13 - spa12*spa35*spb13 + spa15*spa24*spb14 - spa12*spa45*spb14 + spa23*spa25*spb23 + spa24*spa25*spb24)*(s124*spa15 - spa56*t49) - d11*spa24*spb16*(-((s12*t123 - s34*t123 - s56*t123 - spa25*spb56*t146)*(spa13*spb23 + t291)) - (s14 - s23)*(spa15*spb56*t146 + (-s12 + s34 + s56)*t49)) + (s14 - s23)*spa12*spb13*(spa15*spb16 + spa25*spb26)*T(2)); 
complex<T> t13 = -(d34*spa15*t152) + d35*((spa15*spb12 - spa56*spb26)*t127 + d6*(-s13 - s14 + s23 + s24)*spa15*t141*t152 - (s12 - s34 + s56)*spb13*(spa24*spa56*spb26 + spa15*t156)) - d36*t218*(spa15*spb56*t146 + (-s12 + s34 + s56)*t49); 
complex<T> t15 = d14*spb46*t151 + d13*(s12*spa35*spb34 - s34*spa35*spb34 + s56*spa35*spb34 + spa56*spb46*t146)*t218 + d16*(-(spa34*spb46*t127) + spa35*spb56*t127 + d15*(s13 - s14 + s23 - s24)*spb46*t140*t151 + (s12 - s34 - s56)*spa24*t83); 
complex<T> t77 = d42*spa24*spb16*t21 + spb16*t234*square(spa24) + t160*t231*square(spb13); 
complex<T> t214 = -t262; 
complex<T> t228 = d56*s34; 
complex<T> t233 = d29*t148; 
complex<T> t248 = d54*spa23; 
complex<T> t249 = d4*spa12; 
complex<T> t255 = d47*t148; 
complex<T> t277 = d21*t154; 
complex<T> t302 = d50*t167; 
complex<T> t303 = d49*t168; 
complex<T> t305 = d8*t218; 
complex<T> t340 = t261*t262; 
complex<T> t341 = spb34*t263; 
complex<T> t367 = d33*spb13; 
complex<T> d2 = (s12 - s123)*(s12 - s123 + s23)*spa13*t185; d2 = T(1)/d2;
complex<T> d19 = (-s123 + s23)*(s12 - s123 + s23)*spa13*t185; d19 = T(1)/d19;
complex<T> d20 = (-s123 + s23)*spa13*t185; d20 = T(1)/d20;
complex<T> d51 = spa23*t185*cube(spa13); d51 = T(1)/d51;
complex<T> d52 = t185*cube(spa13); d52 = T(1)/d52;
complex<T> d53 = spa23*t185*cube(spa13*spb14 + spa23*spb24); d53 = T(1)/d53;
complex<T> d55 = spa13*t185*square(s12 - s123 + s23); d55 = T(1)/d55;
complex<T> t6 = -(d5*spa15*spa24*spb13) + t305*t350*(-((s12 - s34 + s56)*spb26) + spa15*spb12*spb56*T(2)) + d7*((s12 - s34 - s56)*spa12*spb13*(spa24*spa56*spb26 + spa15*t156) + s234*spa15*T(2)*((-s13 - s14 + s23 + s24)*(spa24*spb13 - d6*spa14*spb23*t151) + (spa14*spb13 - spa24*spb23)*t151*T(2)) - t127*(-(s12*spa15) + s56*spa15 + spa12*spa35*spb23 + spa12*spa45*spb24 - s34*spa15*T(3))); 
complex<T> t7 = -(d39*spa24*spb13*spb46) + t305*t350*((s12 - s34 - s56)*spa35 + spa34*spb46*t57) + d40*((s12 - s34 + s56)*t316*t83 + spb46*t89*((s13 - s14 + s23 - s24)*(spa24*spb13 - d15*spa23*spb14*t152) + (-(spa23*spb13) + spa24*spb14)*t152*T(2)) - t127*(spa13*spb16*spb34 + spa23*spb26*spb34 + spb46*(-s34 + s56 - s12*T(3)))); 
complex<T> t34 = s23*(-(d70*spb34*t126) + t166*t303 + d69*t75*cube(spb34)); 
complex<T> t35 = s12*(d52*spa12*t214 + t216*t302 + d51*t76*cube(spa12)); 
complex<T> t52 = t166*t303 - d48*t75*cube(t148) + t255*t99*square(spa15); 
complex<T> t53 = t248*t287 + t216*t302 + d53*t76*cube(spa12*spb14 - spa23*spb34); 
complex<T> t304 = d20*t171; 
complex<T> t306 = s234*t233; 
complex<T> t307 = t160*t277; 
complex<T> t1 = -(d17*t15*t23) + d2*t214*t261 + s123*spa23*spb13*t249*t47 + d9*spa15*t6 - d1*t69*t76*t85 - t12*T(2) - spa12*t131*t145*T(2); 
complex<T> t2 = d17*t15*t23 + d22*spa24*t287 + t128*t129*t367 + spb13*t306*t48 - d9*spa15*t6 + d23*t221*t68 + d30*t264*t69 - d41*spb46*t7 + d17*t13*t8 + t12*T(2) + d12*t14*T(2) + spa24*spb46*t307*T(2) - d45*(d42*spa24*spb16*t21 + spb16*t234*t68 + t160*t231*t69)*T(3); 
complex<T> t3 = d17*spa15*spb34*t13 + d31*t126*t215 + d32*(spa35*spb23 + spa45*spb24)*t129*t316 + t316*t388*t48 + d41*spb46*t7 - d27*t68*t75*t87 - d12*t14*T(2); 
complex<T> t19 = d47*spb23*t264*(d45*(s12 - s34 - s56) + d57*t268) + (d45*(s12 - s34 - s56) + d57*t268)*t52 + d64*(spa14*spa25*spb16*spb34*t212 - d38*spa24*spb36*spb46*t154*(spa13*spb23 + t291) - d11*(s14 + s24 + s34)*spa24*spb16*t23*(spa13*spb23 + t291) + t288*t327 + t278*t350 + spa15*(spa14*spb16*spb34 + (s234 + spa34*spb34)*spb36)*t46 + d63*spa56*spb16*spb46*t46*square(spa14)) + d58*t154*t21*t23*T(2) - d61*(-(d59*t167*t216) + d60*t154*square(t160))*(s12*s123 + s123*(-s34 + s56) - t224*T(2)) + d12*spb13*(s56*spa24*t114 - spa25*(s12*s123 - s123*(s34 + s56) + s56*t146)*t49*T(2) + spa12*spa15*spa45*spb56*T(2)*(spb14*t146 - spb12*t46*T(2))) + d67*t23*(-(spa24*spa56*spb13*spb36) + spa56*(d65*(-s13 + s14 - s23 + s24)*spa23*spb14*spb46*t152 + d66*spb34*(spb46*t140 + spa35*spb56*t146)*t218) - d15*t152*(spa35*spb34*t151 + spa23*spa56*spb13*spb46*T(3))); 
complex<T> t20 = t221*(-(d45*(s12 - s34 + s56)) + d75*t224)*t248 + (-(d45*(s12 - s34 + s56)) + d75*t224)*t53 + d64*(spa12*spa13*spb14*spb36*t124 - d10*spa25*spb13*(spa13*spb14 + spa23*spb24)*t128 - d37*(s12 + s13 + s14)*spa15*(spa13*spb14 + spa23*spb24)*spb34*t21 + (s123*spa25 + spa12*spa25*spb12 + spa12*spa45*spb14)*spb23*t23 + t288*t327 + t278*t350 + d76*spa12*spa15*spa45*spb23*spb56*square(spb14)) + d71*spb16*t128*t316*T(2) + d74*(-(d73*t166*t168) + d72*t148*square(t158))*(s12*s234 - s234*(s34 + s56) + t268*T(2)) - d12*spa24*(-(s56*spb13*t114) + spb34*(-(s234*(s12 - s34 + s56)) + s56*t146)*t56*T(2) + spb16*spb34*spb46*t57*(-(spa14*t146) + spa12*spa34*spb23*T(2))) + d67*t8*(spa24*spa25*spb13*spb56 + d77*(-s13 - s14 + s23 + s24)*spa14*spa15*spb23*spb56*t151 + d78*spa12*spb56*(spa15*t141 - spa56*spb26*t146)*t218 + d6*t151*(t152*t49 + spa15*spa24*spb23*spb56*T(3))); 
complex<T> t25 = s34*(d45*t34 + d81*s23*t341); 
complex<T> t26 = (-s156 + s34)*(-(d70*spb34*t126) + t166*t303 + d68*t341 - t52 + t255*t48*t99 + d69*t75*cube(spb34)); 
complex<T> t28 = (-s123 + s23)*(d52*spa12*t214 + t248*t287 + t216*t302 + d55*t340 - t53 + d51*t76*cube(spa12)); 
complex<T> t36 = s34*t53; 
complex<T> t60 = d26*t129*t182*t316 + d32*spa35*spb23*t267*t316 + d32*spa45*spb24*t267*t316 + d24*t341 + d25*t341 + d31*t341 + t128*t267*t367 + d27*t68*t75*t87 - (spb13*t306 + t316*t388 - d30*s234*t148*t69)*square(spa15); 
complex<T> t62 = d24*t126*t215 + d25*t126*t215 + d18*t214*t261 + d19*t214*t261 + d26*t182*t267*t316 + spa12*t131*t304*T(2); 
complex<T> t65 = -(d22*spa24*t287) + d18*t340 + d19*t340 + d2*t340 + d23*t287*t68 + d1*t69*t76*t85 - s123*spa23*spb13*t249*square(spb46) + spa12*t131*t145*T(2) + d21*spa24*(spa12*spb14 - spa23*spb34)*spb46*t160*T(2) - spa12*t131*t304*T(2); 
complex<T> t67 = s23*(d82*s12*t340 + d45*t35); 
complex<T> t207 = t34 + s23*(t221*t248 + d68*t341 + t53 + t77); 
complex<T> t223 = d45*t52; 
complex<T> t27 = d42*s12*spa24*spa45*spb13*spb16 + t35 + s12*(t223 + d46*spb23*t264 + d56*spa23*t287 + d55*t340 - d45*t53 + spb16*t234*t68 + t160*t231*t69); 
complex<T> t29 = d79*spa23*spa56*t221 + d80*spb23*spb56*t264 - s56*(t223 + d45*t53); 
complex<T> t208 = spa23*t221*t228 + d45*t36 + s34*(t223 + d46*spb23*t264 + t77); 
complex<T> t210 = s12*(s234*t223 - d46*spb23*t135*t148*t48); 
complex<T> t245 = s123*(s123*t132*t154*t228 + d45*t36); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t65*(*CI_users[1]->get_value(mc,ind,mu)) + t62*(*CI_users[2]->get_value(mc,ind,mu)) + t60*(*CI_users[3]->get_value(mc,ind,mu)) + t3*(*CI_users[4]->get_value(mc,ind,mu)) + t2*(*CI_users[5]->get_value(mc,ind,mu)) + t27*(*CI_users[6]->get_value(mc,ind,mu)) + t28*(*CI_users[7]->get_value(mc,ind,mu)) + t19*(*CI_users[8]->get_value(mc,ind,mu)) + t207*(*CI_users[9]->get_value(mc,ind,mu)) + t26*(*CI_users[10]->get_value(mc,ind,mu)) + t208*(*CI_users[11]->get_value(mc,ind,mu)) + t20*(*CI_users[12]->get_value(mc,ind,mu)) + t29*(*CI_users[13]->get_value(mc,ind,mu)) + t210*(*CI_users[14]->get_value(mc,ind,mu)) + t25*(*CI_users[15]->get_value(mc,ind,mu)) + t67*(*CI_users[16]->get_value(mc,ind,mu)) + t245*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpmmqmemep_L_wCI::\
C2q2g2l_qpmmqmemep_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpmmqmemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, m, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpmmqmemep L");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:34:21 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spa12 = SPA(1,2);
complex<T> spb36 = SPB(3,6);
complex<T> spa24 = SPA(2,4);
complex<T> spb46 = SPB(4,6);
complex<T> spb13 = SPB(1,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> t4 = spb56*T(2); 
complex<T> t7 = spb12*spb23; 
complex<T> t10 = square(spb16); 
complex<T> t13 = square(spa23); 
complex<T> t15 = square(spb36); 
complex<T> t16 = s123*s234 - s23*s56; 
complex<T> t21 = -(spa34*spb13); 
complex<T> d1 = (s12 - s123)*spb23*spb34*spb56; d1 = T(1)/d1;
complex<T> t3 = square(spa24*spb12*spb46 - spb46*t21); 
complex<T> t19 = d1*spa23; 
complex<T> d2 = spb23*spb34*t4*square(s12 - s123); d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spb34*spb56*t7; d3 = T(1)/d3;
complex<T> d4 = spb34*t4*t7*square(s123 - s56); d4 = T(1)/d4;
complex<T> d5 = spb34*t4*t7; d5 = T(1)/d5;
complex<T> d6 = spb34*t7; d6 = T(1)/d6;
complex<T> d7 = spb34*t4; d7 = T(1)/d7;
complex<T> d8 = spb12*t4; d8 = T(1)/d8;
complex<T> t1 = d2*spb12*t13*t15 - d4*t3 + d3*(spa24*spb12 + spa34*spb13)*spb16*spb46*T(2) + spb16*spb36*t19*T(2); 
complex<T> t9 = -(d2*spb12*t13*t15) - spb16*spb36*t19*T(2); 
complex<T> t31 = d5*t10; 
complex<T> t2 = d4*t3 - d3*(spa24*spb12 + spa34*spb13)*spb16*spb46*T(2) + t31*T(3); 
complex<T> co1 = d6*spa56*t10; 
complex<T> co2 = d7*spa12*spa23*t10; 
complex<T> co3 = t16*t31; 
complex<T> co4 = d8*spa23*spa34*t10; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t9*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qmppqpemep_L_wCI::\
C2q2g2l_qmppqpemep_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qmppqpemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, p, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qmppqpemep L");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:34:24 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb56 = SPB(5,6);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> t4 = spa56*T(2); 
complex<T> t7 = spa12*spa23; 
complex<T> t11 = square(spa15); 
complex<T> t14 = square(spa35); 
complex<T> t15 = square(spb23); 
complex<T> t16 = s123*s234 - s23*s56; 
complex<T> t23 = spa13*spb34; 
complex<T> d1 = (s12 - s123)*spa23*spa34*spa56; d1 = T(1)/d1;
complex<T> t3 = square(spa12*spa45*spb24 + spa45*t23); 
complex<T> t22 = d1*spa35; 
complex<T> d2 = spa23*spa34*t4*square(s12 - s123); d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spa34*spa56*t7; d3 = T(1)/d3;
complex<T> d4 = spa34*t4*t7*square(s123 - s56); d4 = T(1)/d4;
complex<T> d5 = spa34*t4*t7; d5 = T(1)/d5;
complex<T> d6 = spa34*t7; d6 = T(1)/d6;
complex<T> d7 = spa34*t4; d7 = T(1)/d7;
complex<T> d8 = spa12*t4; d8 = T(1)/d8;
complex<T> t9 = -(d2*spa12*t14*t15) - spa15*spb23*t22*T(2); 
complex<T> t10 = d3*(spa12*spb24 + t23); 
complex<T> t27 = d5*t11; 
complex<T> t1 = d2*spa12*t14*t15 - d4*t3 + spa15*(spa45*t10 + spb23*t22)*T(2); 
complex<T> t2 = d4*t3 - spa15*spa45*t10*T(2) + t27*T(3); 
complex<T> co1 = d6*spb56*t11; 
complex<T> co2 = d7*spb12*spb23*t11; 
complex<T> co3 = t16*t27; 
complex<T> co4 = d8*spb23*spb34*t11; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t9*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)) + co3*(*CI_users[5]->get_value(mc,ind,mu)) + co4*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qmpmqpemep_L_wCI::\
C2q2g2l_qmpmqpemep_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qmpmqpemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, m, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qmpmqpemep L");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:48:23 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb46 = SPB(4,6);
complex<T> spb56 = SPB(5,6);
complex<T> spa45 = SPA(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb26 = SPB(2,6);
complex<T> spb36 = SPB(3,6);
complex<T> spb16 = SPB(1,6);
complex<T> spb14 = SPB(1,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s24 = -(spa24*spb24);
complex<T> s134 = SS(1,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s124 = SS(1,2,4);
complex<T> t8 = -(spa45*spb12); 
complex<T> t22 = spa15*spb24; 
complex<T> t23 = spa34*spb16; 
complex<T> t45 = spa23*spb12; 
complex<T> t46 = square(spa45); 
complex<T> t47 = square(spb16); 
complex<T> t49 = spa34*spb36; 
complex<T> t50 = s234*(-(spa23*spb12) + spa34*spb14); 
complex<T> t55 = spa25*spb26; 
complex<T> t67 = square(spa13); 
complex<T> t68 = square(spb24); 
complex<T> t70 = spa56*(spa24*spb12 + spa34*spb13); 
complex<T> t71 = square(spa24*spb26 + spa34*spb36); 
complex<T> t74 = square(spa15*spb13 + spa25*spb23); 
complex<T> t75 = square(spa23*spb36 + spa24*spb46); 
complex<T> t84 = square(spa34); 
complex<T> t85 = square(spb12); 
complex<T> t89 = s234*T(2); 
complex<T> t100 = s123*spb23; 
complex<T> t124 = (spa24*spb12 + spa34*spb13)*spb56; 
complex<T> t125 = -(spa15*spb12); 
complex<T> t126 = spa34*spb46; 
complex<T> t130 = spa15*spb12*T(2) - spa35*spb23*T(2); 
complex<T> t132 = spb24*(-(spa23*spb26) + spa34*spb46); 
complex<T> t135 = s123*spb13; 
complex<T> t147 = -s12 - s34 + s56; 
complex<T> t149 = -(spa14*spb12) + spa34*spb23; 
complex<T> t150 = s12 - s34 - s56; 
complex<T> t151 = -s12 + s34 - s56; 
complex<T> t153 = -(spa13*spb23) - spa14*spb24; 
complex<T> t154 = -(spa13*spb14) - spa23*spb24; 
complex<T> t156 = spa12*spb24 + spa13*spb34; 
complex<T> t163 = -s234 + s34; 
complex<T> t164 = s12 - s123 + s23; 
complex<T> t168 = s23 - s234 + s34; 
complex<T> t169 = cube(spa13); 
complex<T> t170 = square(spa15); 
complex<T> t171 = spa25*spb12 + spa35*spb13; 
complex<T> t174 = square(spb46); 
complex<T> t181 = -square(spa45); 
complex<T> t182 = cube(spb13); 
complex<T> t185 = -(spa24*spb26) - spa34*spb36; 
complex<T> t187 = square(s123); 
complex<T> t216 = spa13*spb12; 
complex<T> t218 = spa12*T(2); 
complex<T> t224 = square(s234); 
complex<T> t227 = s12*s56; 
complex<T> t233 = s23*spb34; 
complex<T> t244 = spa13*spb16; 
complex<T> t247 = cube(spb24); 
complex<T> t259 = spa14*spb13; 
complex<T> t262 = spa23*spa34; 
complex<T> t264 = spa45*spb56; 
complex<T> t265 = spa13*spb46; 
complex<T> t267 = s34*s56; 
complex<T> t279 = s14*(spa24*spb12 + spa34*spb13); 
complex<T> t309 = spa45*spb16; 
complex<T> t347 = spa34*spb12; 
complex<T> d9 = spa24*spb12 + spa34*spb13; d9 = T(1)/d9;
complex<T> d10 = spa56*spb34; d10 = T(1)/d10;
complex<T> d11 = spb56; d11 = T(1)/d11;
complex<T> d14 = spb34*spb56*T(2); d14 = T(1)/d14;
complex<T> d15 = spa23*spb13 + spa24*spb14; d15 = T(1)/d15;
complex<T> d17 = (spa24*spb12 + spa34*spb13)*(spa23*spb13 + spa24*spb14); d17 = T(1)/d17;
complex<T> d37 = spa56; d37 = T(1)/d37;
complex<T> d38 = spa12*spb56; d38 = T(1)/d38;
complex<T> d40 = spb56*square(spa23*spb13 + spa24*spb14)*T(2); d40 = T(1)/d40;
complex<T> d42 = s123*s23*s56*spa12; d42 = T(1)/d42;
complex<T> d45 = T(2); d45 = T(1)/d45;
complex<T> d61 = s234; d61 = T(1)/d61;
complex<T> d63 = spb34; d63 = T(1)/d63;
complex<T> d71 = s123; d71 = T(1)/d71;
complex<T> d72 = s234*(spa24*spb12 + spa34*spb13)*(spa23*spb13 + spa24*spb14); d72 = T(1)/d72;
complex<T> d73 = spa12; d73 = T(1)/d73;
complex<T> d77 = square(spa23*spb13 + spa24*spb14); d77 = T(1)/d77;
complex<T> d79 = (spa24*spb12 + spa34*spb13)*cube(spa23*spb13 + spa24*spb14)*T(2); d79 = T(1)/d79;
complex<T> t14 = d10*((s14 - s23)*(spa56*spb16*t147 + spa25*spb12*t150) - (spa23*spb13 + spa24*spb14)*(spa56*spb26*t147 + t125*t150))*t22 + d11*spa13*(s124*spb16 + spa25*spb12*spb56)*(-(spa14*spb16*spb24) - spa23*spb23*spb26 - spa24*spb24*spb26 + spa13*spb12*spb36 + spa14*spb12*spb46 - spb23*t244) + (s14 - s23)*t216*(spa15*spb16 + t55)*T(2); 
complex<T> t57 = t153*t154; 
complex<T> t69 = square(t171); 
complex<T> t113 = t147*t309 + spa25*spb36*t347*T(2); 
complex<T> t129 = spa45*t149; 
complex<T> t133 = -(spa23*t47); 
complex<T> t137 = spa56*(spa24*spb23 + t259); 
complex<T> t138 = d45*(-s14 + s23); 
complex<T> t148 = square(s12) + square(s34) + square(s56) - s12*s34*T(2) - t227*T(2) - t267*T(2); 
complex<T> t157 = spa35*spb23 + t125; 
complex<T> t158 = -(spa34*spb14) + t45; 
complex<T> t161 = -(spa23*spb26) + t126; 
complex<T> t178 = s234*t124; 
complex<T> t186 = -(spa35*spb36) - t55; 
complex<T> t215 = -t262; 
complex<T> t225 = s123*t149; 
complex<T> t261 = t68*t71; 
complex<T> t290 = t47*t50; 
complex<T> t314 = spb34*t156; 
complex<T> t322 = t70*T(2); 
complex<T> t323 = -(t169*t174); 
complex<T> t331 = t135*t70; 
complex<T> t346 = spb56*t156; 
complex<T> t350 = spa24*t89; 
complex<T> d2 = (s12 - s123)*t70*square(spa24*spb23 + t259); d2 = T(1)/d2;
complex<T> d5 = spa56*square(spa24*spb23 + t259)*T(2); d5 = T(1)/d5;
complex<T> d6 = spa24*spb23 + t259; d6 = T(1)/d6;
complex<T> d21 = (-s123 + s56)*t70*square(spa24*spb23 + t259); d21 = T(1)/d21;
complex<T> d22 = (-s123 + s56)*spb23*(spa24*spb23 + t259)*t70; d22 = T(1)/d22;
complex<T> d27 = (-s234 + s56)*t124*square(spa23*spb13 + spa24*spb14); d27 = T(1)/d27;
complex<T> d28 = spa23*(spa23*spb13 + spa24*spb14)*t124*square(s234 - s56)*T(2); d28 = T(1)/d28;
complex<T> d29 = t124*t163*square(spa23*spb13 + spa24*spb14); d29 = T(1)/d29;
complex<T> d32 = (-s234 + s56)*spa23*(spa23*spb13 + spa24*spb14)*t124; d32 = T(1)/d32;
complex<T> d33 = s234*spa24*(spa23*spb13 + spa24*spb14)*spb56*t163; d33 = T(1)/d33;
complex<T> d34 = spa56*t218; d34 = T(1)/d34;
complex<T> d36 = (spa24*spb12 + spa34*spb13)*(spa24*spb23 + t259); d36 = T(1)/d36;
complex<T> d43 = s56*spa12*t233; d43 = T(1)/d43;
complex<T> d44 = s234*s56*t233; d44 = T(1)/d44;
complex<T> d48 = s123*t182*t70; d48 = T(1)/d48;
complex<T> d49 = t100*t182*t70; d49 = T(1)/d49;
complex<T> d51 = t70*cube(spa24*spb23 + t259); d51 = T(1)/d51;
complex<T> d52 = t100*t70*cube(spa24*spb23 + t259); d52 = T(1)/d52;
complex<T> d54 = t124*cube(spa23*spb13 + spa24*spb14); d54 = T(1)/d54;
complex<T> d56 = t124*cube(spa23*spb13 + spa24*spb14)*T(2); d56 = T(1)/d56;
complex<T> d57 = s123*(spa24*spb12 + spa34*spb13)*(spa24*spb23 + t259); d57 = T(1)/d57;
complex<T> d58 = spb23*(spa24*spb23 + t259)*t70; d58 = T(1)/d58;
complex<T> d60 = t187*T(2); d60 = T(1)/d60;
complex<T> d62 = spb34*t218; d62 = T(1)/d62;
complex<T> d64 = (spa24*spb12 + spa34*spb13)*(spa23*spb13 + spa24*spb14)*(spa24*spb23 + t259)*T(2); d64 = T(1)/d64;
complex<T> d65 = square(spa24*spb23 + t259); d65 = T(1)/d65;
complex<T> d75 = spa23*(spa23*spb13 + spa24*spb14)*t124; d75 = T(1)/d75;
complex<T> d76 = t224*T(2); d76 = T(1)/d76;
complex<T> d80 = (spa24*spb12 + spa34*spb13)*cube(spa24*spb23 + t259)*T(2); d80 = T(1)/d80;
complex<T> t48 = -t225; 
complex<T> t127 = spb23*t69; 
complex<T> t139 = d62*t147; 
complex<T> t146 = d33*(spa23*spb36 + spa24*spb46); 
complex<T> t210 = d56*s12; 
complex<T> t222 = -(t57*T(3)); 
complex<T> t232 = d44*spa15; 
complex<T> t266 = t157*T(2); 
complex<T> t274 = d32*t158; 
complex<T> t300 = t138*t186; 
complex<T> t315 = d21*t225; 
complex<T> t329 = t261*t262; 
complex<T> t339 = spa23*t290; 
complex<T> t343 = spb24*t129; 
complex<T> t356 = t158*t47; 
complex<T> d1 = (s12 - s123)*t164*t331; d1 = T(1)/d1;
complex<T> d3 = t135*t137*square(s12 - s123)*T(2); d3 = T(1)/d3;
complex<T> d4 = (s12 - s123)*t135*t137; d4 = T(1)/d4;
complex<T> d7 = spb34*t137*t148*T(2); d7 = T(1)/d7;
complex<T> d8 = (spa24*spb23 + t259)*square(t148); d8 = T(1)/d8;
complex<T> d12 = (spa23*spb13 + spa24*spb14)*t148*(spa24*spb23 + t259); d12 = T(1)/d12;
complex<T> d13 = square(t148); d13 = T(1)/d13;
complex<T> d16 = spb34*spb56*t148*T(2); d16 = T(1)/d16;
complex<T> d18 = t135*t322*square(s123 - s23); d18 = T(1)/d18;
complex<T> d19 = (-s123 + s23)*t164*t331; d19 = T(1)/d19;
complex<T> d20 = (-s123 + s23)*t331; d20 = T(1)/d20;
complex<T> d23 = spb23*(spa24*spb23 + t259)*t322*square(s123 - s56); d23 = T(1)/d23;
complex<T> d24 = t124*t350*square(s23 - s234); d24 = T(1)/d24;
complex<T> d25 = (s23 - s234)*spa24*(s23 + t163)*t178; d25 = T(1)/d25;
complex<T> d26 = (s23 - s234)*spa24*t178; d26 = T(1)/d26;
complex<T> d30 = spa24*t163*(s23 + t163)*t178; d30 = T(1)/d30;
complex<T> d31 = (spa23*spb13 + spa24*spb14)*spb56*t350*square(t163); d31 = T(1)/d31;
complex<T> d35 = spa56*t148*t218; d35 = T(1)/d35;
complex<T> d39 = (spa23*spb13 + spa24*spb14)*square(t148); d39 = T(1)/d39;
complex<T> d41 = (spa23*spb13 + spa24*spb14)*spb56*t148*t218; d41 = T(1)/d41;
complex<T> d46 = t331*square(t164); d46 = T(1)/d46;
complex<T> d47 = t322*cube(spa24*spb23 + t259); d47 = T(1)/d47;
complex<T> d50 = s123*spa12*spa23*t346; d50 = T(1)/d50;
complex<T> d53 = s234*spa56*spb23*t314; d53 = T(1)/d53;
complex<T> d55 = spa23*t178*cube(spa23*spb13 + spa24*spb14); d55 = T(1)/d55;
complex<T> d59 = spa12*spa23*t346; d59 = T(1)/d59;
complex<T> d66 = t148*(spa24*spb23 + t259); d66 = T(1)/d66;
complex<T> d67 = (spa24*spb12 + spa34*spb13)*t148; d67 = T(1)/d67;
complex<T> d68 = t178*cube(spa24); d68 = T(1)/d68;
complex<T> d69 = spa23*t178*cube(spa24); d69 = T(1)/d69;
complex<T> d70 = spa24*t178*square(s23 + t163); d70 = T(1)/d70;
complex<T> d74 = spa56*spb23*t314; d74 = T(1)/d74;
complex<T> d78 = (spa23*spb13 + spa24*spb14)*t148; d78 = T(1)/d78;
complex<T> d81 = t124*t350*square(s23 + t163); d81 = T(1)/d81;
complex<T> d82 = t135*t322*square(t164); d82 = T(1)/d82;
complex<T> t12 = d12*(-(d37*spb24*(spa13*spa35*spb13 - spa15*spa34*spb14 + spa13*spa45*spb14 + spa23*spa35*spb23 - spa25*spa34*spb24 + spa23*spa45*spb24)*(s134*spa45 + spa56*t49)) - d38*t265*((spa35*spb56*t147 - t126*t151)*(spa24*spb23 + t259) - (s14 - s23)*(t147*t264 + t151*t49)) + (s14 - s23)*spa34*spb24*(spa35*spb36 + spa45*spb46)*T(2)); 
complex<T> t13 = -(d5*spa45*t153) - d8*t222*(t147*t264 + t151*t49) + d7*(-(spa13*spa56*spb24*spb36*t150) + d6*(s13 - s14 + s23 - s24)*spa45*spb34*t151*t153 + spa45*spb24*t150*t156 + spa45*spb34*t57 + spa56*spb36*t57); 
complex<T> t15 = d40*spb16*t154 + d39*(spa56*spb16*t147 + spa25*spb12*t150)*t222 + d41*(d15*(s13 + s14 - s23 - s24)*spa12*spb16*t150*t154 + spa13*t151*(spa25*spb24*spb56 - spb16*t156) - (spa12*spb16 + spa25*spb56)*t57); 
complex<T> t51 = d51*t100*t149*t181 + d50*t323 + d52*t74*cube(t149); 
complex<T> t76 = d43*t22*t265 - d42*spb46*t157*square(spa13) + t161*t232*square(spb24); 
complex<T> t230 = -(d47*s34); 
complex<T> t255 = -(d53*t170); 
complex<T> t286 = spb12*t127; 
complex<T> t289 = t46*t48; 
complex<T> t296 = spa15*t139; 
complex<T> t299 = d26*t185; 
complex<T> t301 = d13*t222; 
complex<T> t6 = -(d34*spa13*spa45*spb24) + t301*t347*(-(spb36*t150) + spb34*t264*T(2)) + d35*(spa34*spb24*t151*(-(spa13*spa56*spb36) + spa45*t156) + s123*spa45*T(2)*((s13 - s14 + s23 - s24)*(spa13*spb24 - d6*spa14*spb23*t154) + (-(spa13*spb23) + spa14*spb24)*t154*T(2)) + t57*(s34*spa45 - s56*spa45 - spa15*spa34*spb13 - spa25*spa34*spb23 + s12*spa45*T(3))); 
complex<T> t7 = -(d14*spb24*t244) + (-(spa25*t151) + spa56*spb16*t218)*t301*t347 + d16*(t150*(-(spa25*spb24*spb56) + spb16*t156)*t216 + spb16*t89*((-s13 - s14 + s23 + s24)*(spa13*spb24 - d15*spa23*spb14*t153) + (spa13*spb14 - spa23*spb24)*t153*T(2)) + t57*(s12*spb16 - s56*spb16 - spa23*spb12*spb36 - spa24*spb12*spb46 + s34*spb16*T(3))); 
complex<T> t30 = -(d45*t51); 
complex<T> t36 = -(s12*(d48*t286 + d50*t323 - d49*t74*cube(spb12))); 
complex<T> t37 = -(s23*(t247*t255 + d68*t262*t71 - d69*t75*cube(spa34))); 
complex<T> t52 = t247*t255 + d55*t75*cube(spa23*spb12 - spa34*spb14) + d54*s234*spa23*(-(spa23*spb12) + spa34*spb14)*square(spb16); 
complex<T> t63 = d24*t329 + d25*t329 + d30*t329 + d31*t68*t75*t84 - d27*s234*spa13*(spa23*spb12 - spa34*spb14)*square(spb16) + d29*s234*spb24*t215*square(spb16) + d28*s234*(spa23*spb12 - spa34*spb14)*t67*square(spb16) + spa34*t132*t146*T(2) + d32*(spa23*spb12 - spa34*spb14)*t161*t244*T(2) - spa34*t132*t299*T(2); 
complex<T> t260 = -t286; 
complex<T> t321 = spb23*t289; 
complex<T> t328 = t286*t67; 
complex<T> t336 = spb46*t296; 
complex<T> t1 = -(d9*t15*t23) + d30*t215*t261 + d29*s234*spb24*t262*t47 + d36*spa45*t6 - d31*t68*t75*t84 + t12*T(2) - spa34*t132*t146*T(2); 
complex<T> t2 = d9*t15*t23 + d22*t266*t343 + d27*s234*spa13*t356 + spb24*t315*t46 - d36*spa45*t6 + d28*t290*t67 + d23*t289*t68 - d17*spb16*t7 + d9*t13*t8 - t12*T(2) - d12*t14*T(2) - t161*t244*t274*T(2) + d45*(d43*t22*t265 - d42*spb46*t157*t67 + t161*t232*t68)*T(3); 
complex<T> t3 = d9*spa45*spb12*t13 - d4*(spa15*spb13 + spa25*spb23)*t130*t216 + d2*t100*t216*t46 + d1*t260*t67 + d17*spb16*t7 - d3*t67*t74*t85 + d12*t14*T(2); 
complex<T> t18 = d51*(d45*t151 + d71*t227)*t321 - (d45*t151 + d71*t227)*t51 + d64*(-(spa35*spb12*t126*t259) - d11*spb26*t158*t244*(spa24*spb23 + t259) - d38*(s12 + s13 + s14)*t23*(spa24*spb23 + t259)*t265 + t279*t336 + t300*t347 - spa45*(s123*spb26 + spa12*spb12*spb26 + spa14*spb12*spb46)*t45 + d73*spa56*spb16*spb46*t45*square(spa14)) - d72*t158*t22*t23*T(2) + d76*(-(d74*t170*t247) + d75*t158*square(t161))*(s234*t150 + t267*T(2)) - d12*spb24*(s56*spa13*t113 - spa35*(s56*t147 + s234*t151)*t49*T(2) - spa15*spa34*t264*T(2)*(spb14*t147 - spb34*t45*T(2))) + d67*t23*(-(spa13*spa56*spb24*spb26) + d77*(s13 + s14 - s23 - s24)*spa23*spa56*spb14*spb16*t153 + d78*spa56*spb12*(spa25*spb56*t147 + spa12*spb16*t150)*t57*T(3) + d15*t153*(spa25*spb12*t154 - spa23*spa56*spb16*spb24*T(3))); 
complex<T> t20 = d54*(d45*t150 + d61*t267)*t339 - (d45*t150 + d61*t267)*t52 + d64*(spa24*spa34*spb14*spb26*t125 - d10*(s14 + s24 + s34)*spa45*spb12*(spa23*spb13 + spa24*spb14)*t22 - spb23*(s234*spa35 + spa15*spa34*spb14 + spa34*spa35*spb34)*t23 + t279*t336 - d37*spa35*(spa23*spb13 + spa24*spb14)*t343 + t300*t347 + d63*spa15*spa34*spb23*t264*square(spb14)) - d57*spb46*t129*t216*T(2) + d60*(d59*t323 + d58*t149*square(t157))*(s123*t151 + t227*T(2)) + d12*spa13*(-(s56*spb24*t113) + spa56*spb12*spb16*spb46*(spa14*t147 - spa34*spb23*t218)*T(2) + spb12*(s56*t147 + s123*t150)*t55*T(2)) + d67*t8*(spa13*spa35*spb24*spb56 + d65*(s13 - s14 + s23 - s24)*spa14*spb23*t154*t264 - d66*spa34*spb56*(spa56*spb36*t147 + spa45*spb34*t151)*t57*T(3) + d6*t154*(-(t153*t49) + spa13*spb23*t264*T(3))); 
complex<T> t25 = s23*(d82*s12*t328 + d45*t36); 
complex<T> t26 = (-s123 + s23)*(-(d48*t286) - d50*t323 + d46*t328 + d51*t100*t149*t46 + t51 + d49*t74*cube(spb12)); 
complex<T> t27 = (-s156 + s34)*(-(t247*t255) + d70*t329 + d54*s234*spa23*t356 + t52 - d68*t262*t71 + d69*t75*cube(spa34)); 
complex<T> t29 = d80*spb56*t321 + d79*spa56*t339 + d45*s56*(t51 + t52); 
complex<T> t58 = d20*t130*t171*t216 + d2*t100*t181*t216 - d4*(spa15*spb13 + spa25*spb23)*t216*t266 + spb24*t181*t315 + d1*t328 + d18*t328 + d19*t328 + d22*t130*t343 + d3*t67*t74*t85 + d23*t225*t68*square(spa45); 
complex<T> t62 = d24*t215*t261 + d25*t215*t261 + d20*t171*t216*t266 + d18*t260*t67 + d19*t260*t67 + spa34*t132*t299*T(2); 
complex<T> t66 = s34*(d81*s23*t329 + d45*t37); 
complex<T> t209 = d51*s23*t321 + d70*s23*t329 + t37 - s23*t51 - s23*t76; 
complex<T> t226 = -(d45*t52); 
complex<T> t256 = s123*s34*t30 + spb23*t149*t187*t230*t46; 
complex<T> t28 = s34*(t226 - d43*spa15*spb24*t265 + t30 + d47*t321 + d56*t339 + d42*spb46*t157*t67 - t161*t232*t68); 
complex<T> t123 = t210*t339 + t36 + s12*(t226 + d46*t328 + d47*t100*t149*t46 + d45*t51 - t76); 
complex<T> t281 = t133*t158*t210*t224 + s12*s234*t226; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t58*(*CI_users[1]->get_value(mc,ind,mu)) + t62*(*CI_users[2]->get_value(mc,ind,mu)) + t63*(*CI_users[3]->get_value(mc,ind,mu)) + t1*(*CI_users[4]->get_value(mc,ind,mu)) + t2*(*CI_users[5]->get_value(mc,ind,mu)) + t123*(*CI_users[6]->get_value(mc,ind,mu)) + t26*(*CI_users[7]->get_value(mc,ind,mu)) + t20*(*CI_users[8]->get_value(mc,ind,mu)) + t209*(*CI_users[9]->get_value(mc,ind,mu)) + t27*(*CI_users[10]->get_value(mc,ind,mu)) + t28*(*CI_users[11]->get_value(mc,ind,mu)) + t18*(*CI_users[12]->get_value(mc,ind,mu)) + t29*(*CI_users[13]->get_value(mc,ind,mu)) + t281*(*CI_users[14]->get_value(mc,ind,mu)) + t66*(*CI_users[15]->get_value(mc,ind,mu)) + t25*(*CI_users[16]->get_value(mc,ind,mu)) + t256*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qmmpqpemep_L_wCI::\
C2q2g2l_qmmpqpemep_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c56, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c56, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qmmpqpemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, p, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qmmpqpemep L");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:48:53 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spa14 = SPA(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spb46 = SPB(4,6);
complex<T> spb36 = SPB(3,6);
complex<T> spa24 = SPA(2,4);
complex<T> spb56 = SPB(5,6);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb26 = SPB(2,6);
complex<T> spb16 = SPB(1,6);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s123 = SS(1,2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s56 = -(spa56*spb56);
complex<T> s234 = SS(2,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> t11 = -(s12*s34*square(spb46)) + square(spa13*spb14*spb36 + spa23*spb24*spb36 + spa14*spb14*spb46 + spa24*spb24*spb46); 
complex<T> t20 = square(spa12); 
complex<T> t21 = square(spb34); 
complex<T> t23 = spa34*(spa24*spb12 + spa34*spb13); 
complex<T> t24 = square(spa25*spb12 + spa35*spb13); 
complex<T> t25 = square(spa15*spb13 + spa25*spb23); 
complex<T> t26 = square(spa24*spb26 + spa34*spb36); 
complex<T> t27 = square(spa23*spb36 + spa24*spb46); 
complex<T> t35 = spa56*(spa24*spb12 + spa34*spb13); 
complex<T> t38 = spa24*(spa23*spb36 + spa24*spb46); 
complex<T> t42 = s123*spb12; 
complex<T> t43 = -(spa15*spb13) - spa25*spb23; 
complex<T> t45 = spa12*spb24 + spa13*spb34; 
complex<T> t51 = square(spa45); 
complex<T> t54 = square(spb16); 
complex<T> t55 = spa13*spb23 + spa14*spb24; 
complex<T> t72 = s23*s56; 
complex<T> t83 = spb12*T(2); 
complex<T> d12 = T(2); d12 = T(1)/d12;
complex<T> t8 = -(s12*s34*square(spa15)) + square(spa14*spa15*spb14 + spa14*spa25*spb24 - spa13*t43); 
complex<T> t39 = spb56*t45; 
complex<T> t68 = spb56*t23; 
complex<T> t69 = spb13*t43; 
complex<T> t73 = spa56*t45; 
complex<T> t85 = t21*t26; 
complex<T> t90 = t20*t24; 
complex<T> t98 = t21*t51; 
complex<T> t99 = t20*t54; 
complex<T> d1 = t35*t42*square(s123 - s23)*T(2); d1 = T(1)/d1;
complex<T> d2 = (-s123 + s23)*t35*t42; d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spb12*spb23*t35; d3 = T(1)/d3;
complex<T> d4 = spb23*t35*t83*square(s123 - s56); d4 = T(1)/d4;
complex<T> d9 = t42*t72; d9 = T(1)/d9;
complex<T> d10 = spa34*spb12*t72; d10 = T(1)/d10;
complex<T> d11 = s234*spa34*t72; d11 = T(1)/d11;
complex<T> d13 = spb23*t35*t42; d13 = T(1)/d13;
complex<T> t79 = d3*spa45; 
complex<T> t89 = -(d13*spb13); 
complex<T> t91 = s123*t39; 
complex<T> d5 = s234*t68*square(s23 - s234)*T(2); d5 = T(1)/d5;
complex<T> d6 = (s23 - s234)*s234*t68; d6 = T(1)/d6;
complex<T> d7 = spa23*t68*square(s234 - s56)*T(2); d7 = T(1)/d7;
complex<T> d8 = (-s234 + s56)*spa23*t68; d8 = T(1)/d8;
complex<T> d16 = s234*spb23*t55*t73; d16 = T(1)/d16;
complex<T> d17 = s234*spa23*t68; d17 = T(1)/d17;
complex<T> d18 = spa34*(spa13*spb14 + spa23*spb24)*t39*t83; d18 = T(1)/d18;
complex<T> d19 = s234*spb23*spb24*t73; d19 = T(1)/d19;
complex<T> d20 = spa34*t55*t73*t83; d20 = T(1)/d20;
complex<T> t17 = d1*spb13*spb23*t90 - d4*s123*spb13*t98 - d2*s12*spa25*t69*T(2) + d2*spa12*spa35*spb13*t69*T(2) + spb34*t69*t79*T(2); 
complex<T> t18 = -(d5*spa23*spa24*t85) - d1*spb13*spb23*t90 + d6*(spa24*spb26*spb34 - s34*spb36)*t38*T(2) + d2*(s12*spa25 - spa12*spa35*spb13)*t69*T(2); 
complex<T> t77 = d8*spa12; 
complex<T> t88 = -(d17*spa24); 
complex<T> d14 = spa13*spa23*t91; d14 = T(1)/d14;
complex<T> d15 = spa23*(spa13*spb14 + spa23*spb24)*t91; d15 = T(1)/d15;
complex<T> t9 = s12*(-(t25*t89) + d14*cube(spa12)*square(spb46)); 
complex<T> t10 = s23*(-(t27*t88) + d19*cube(spb34)*square(spa15)); 
complex<T> t12 = d4*s123*spb13*t98 + d7*s234*spa24*t99 + spb16*t38*t77*T(2) - spb34*t69*t79*T(2) + d11*d12*spa15*spb34*t38*T(3) - d10*d12*(spa23*spb36 + spa24*spb46)*t43*T(3) - d12*d9*spa12*spb46*t69*T(3); 
complex<T> t13 = t25*t89 + d15*(-(spa12*spb14) + spa23*spb34)*t20*square(spb46); 
complex<T> t14 = t27*t88 + d16*(spa12*spb23 - spa14*spb34)*t21*square(spa15); 
complex<T> t19 = d5*spa23*spa24*t85 - d7*s234*spa24*t99 - d6*spa24*spb26*spb34*t38*T(2) + d6*s34*spb36*t38*T(2) - spb16*t38*t77*T(2); 
complex<T> t1 = (-s123 + s23)*(t13 - t25*t89 + d14*cube(spa12)*square(spb46)); 
complex<T> t2 = t10 + s23*(-t13 - d11*spa15*spb34*t38 + d10*spa23*spb36*t43 + d10*spa24*spb46*t43 + d9*spa12*spb46*t69); 
complex<T> t5 = (-s156 + s34)*(t14 - t27*t88 + d19*cube(spb34)*square(spa15)); 
complex<T> t6 = d12*s56*(t13 + t14); 
complex<T> t7 = -(d12*t13); 
complex<T> t70 = -(d12*t14); 
complex<T> t105 = s34*t7; 
complex<T> t112 = s12*t70; 
complex<T> t3 = t112 + d12*s12*t13 - d11*s12*spa15*spb34*t38 + d10*s12*spa23*spb36*t43 + d10*s12*spa24*spb46*t43 + d9*s12*spa12*spb46*t69 + t9; 
complex<T> t4 = t105 + s34*(-(d11*spa15*spb34*t38) + d10*spa23*spb36*t43 + d10*spa24*spb46*t43 + d9*spa12*spb46*t69 + t70); 
complex<T> co1 = -(d18*spa12*t11); 
complex<T> co2 = -(d20*spb34*t8); 
complex<T> co3 = s234*t112; 
complex<T> co4 = d12*s34*t10; 
complex<T> co5 = d12*s23*t9; 
complex<T> co6 = s123*t105; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t17*(*CI_users[0]->get_value(mc,ind,mu)) + t18*(*CI_users[1]->get_value(mc,ind,mu)) + t19*(*CI_users[2]->get_value(mc,ind,mu)) + t12*(*CI_users[3]->get_value(mc,ind,mu)) + t3*(*CI_users[4]->get_value(mc,ind,mu)) + t1*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t2*(*CI_users[7]->get_value(mc,ind,mu)) + t5*(*CI_users[8]->get_value(mc,ind,mu)) + t4*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t6*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)) + co5*(*CI_users[14]->get_value(mc,ind,mu)) + co6*(*CI_users[15]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qmmmqpemep_L_wCI::\
C2q2g2l_qmmmqpemep_L_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c4, c23, c1, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qmmmqpemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, m, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qmmmqpemep L");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:48:56 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb46 = SPB(4,6);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spa12 = SPA(1,2);
complex<T> spb26 = SPB(2,6);
complex<T> spb16 = SPB(1,6);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s56 = -(spa56*spb56);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> t4 = spb56*T(2); 
complex<T> t7 = spb12*spb23; 
complex<T> t10 = square(spb46); 
complex<T> t13 = square(spa23); 
complex<T> t14 = square(spb26); 
complex<T> t16 = s123*s234 - s23*s56; 
complex<T> t21 = -(spa12*spb24); 
complex<T> t3 = square(spa13*spb16*spb34 - spb16*t21); 
complex<T> d1 = t4*t7*square(s234 - s34); d1 = T(1)/d1;
complex<T> d2 = spb34*t4*t7*square(s234 - s56); d2 = T(1)/d2;
complex<T> d3 = (-s234 + s34)*spb56*t7; d3 = T(1)/d3;
complex<T> d4 = (-s234 + s56)*spb34*spb56*t7; d4 = T(1)/d4;
complex<T> d5 = spb34*t4*t7; d5 = T(1)/d5;
complex<T> d6 = spb34*t7; d6 = T(1)/d6;
complex<T> d7 = spb34*t4; d7 = T(1)/d7;
complex<T> d8 = spb12*t4; d8 = T(1)/d8;
complex<T> t19 = d3*spa23; 
complex<T> t32 = d5*t10; 
complex<T> t1 = d1*spb34*t13*t14 - d2*t3 + d4*spb16*(spa12*spb24 + spa13*spb34)*spb46*T(2) + spb26*spb46*t19*T(2); 
complex<T> t2 = d2*t3 - d4*spb16*(spa12*spb24 + spa13*spb34)*spb46*T(2) + t32*T(3); 
complex<T> t9 = -(d1*spb34*t13*t14) - spb26*spb46*t19*T(2); 
complex<T> co1 = d6*spa56*t10; 
complex<T> co2 = d7*spa12*spa23*t10; 
complex<T> co3 = t16*t32; 
complex<T> co4 = d8*spa23*spa34*t10; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)) + co3*(*CI_users[5]->get_value(mc,ind,mu)) + co4*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qppqmpemep_SLC_wCI::\
C2q2g2l_qppqmpemep_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c2356));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c356));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c4, c23, c1, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qppqmpemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, p, qm, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qppqmpemep SLC");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 18:49:35 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa56 = SPA(5,6);
complex<T> spb14 = SPB(1,4);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb46 = SPB(4,6);
complex<T> spb24 = SPB(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb16 = SPB(1,6);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s124 = SS(1,2,4);
complex<T> s234 = SS(2,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> t1 = spa14*spa23; 
complex<T> t2 = spa12*spa56; 
complex<T> t14 = square(spb24); 
complex<T> t15 = square(spa14); 
complex<T> t16 = square(spa15); 
complex<T> t17 = square(spa25); 
complex<T> t18 = square(spa35); 
complex<T> t19 = -(spa13*spa45); 
complex<T> t29 = -s234 + s56; 
complex<T> t30 = square(spa13); 
complex<T> t31 = -(spa15*spa34) + spa13*spa45; 
complex<T> t32 = spa23*spb12 - spa34*spb14; 
complex<T> t39 = -(spa13*spb14) - spa23*spb24; 
complex<T> t41 = s123*s234 - s23*s56; 
complex<T> t44 = -(spa25*spa34) + spa23*spa45; 
complex<T> t45 = square(spb16); 
complex<T> t46 = square(spb46); 
complex<T> t47 = spa34*spb24; 
complex<T> t49 = s123*s124 - s12*s56; 
complex<T> t59 = spa24*spa56; 
complex<T> t69 = cube(spa13); 
complex<T> t74 = spa13*spa35; 
complex<T> t85 = spa35*spb24; 
complex<T> t93 = spa13*spa34; 
complex<T> d25 = spa12*spa23*cube(spa14); d25 = T(1)/d25;
complex<T> d29 = spa14*spa34*spa56*T(2); d29 = T(1)/d29;
complex<T> t25 = -(t32*T(2)); 
complex<T> t48 = t2*T(2); 
complex<T> t55 = -t74; 
complex<T> t75 = spa34*t14; 
complex<T> t84 = -t93; 
complex<T> t107 = spa13*t18; 
complex<T> d1 = (-s123 + s23)*t15*t2; d1 = T(1)/d1;
complex<T> d2 = (-s123 + s56)*spa23*t15*t2; d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spa12*t1; d3 = T(1)/d3;
complex<T> d4 = spa12*t1*square(s123 - s56); d4 = T(1)/d4;
complex<T> d5 = (s23 - s234)*t15*t2; d5 = T(1)/d5;
complex<T> d6 = (s23 - s234)*spa14*t59; d6 = T(1)/d6;
complex<T> d7 = (s23 - s234)*(s23 - s234 + s34)*spa24*t2; d7 = T(1)/d7;
complex<T> d8 = spa14*t59*square(s23 - s234)*T(2); d8 = T(1)/d8;
complex<T> d10 = spa23*t15*t2*t29; d10 = T(1)/d10;
complex<T> d11 = spa34*t1*t2*t29; d11 = T(1)/d11;
complex<T> d12 = spa12*spa34*t1*t29; d12 = T(1)/d12;
complex<T> d13 = spa12*spa34*t1*square(t29)*T(2); d13 = T(1)/d13;
complex<T> d14 = (-s234 + s34)*spa24*t2; d14 = T(1)/d14;
complex<T> d16 = (-s234 + s34)*(s23 - s234 + s34)*spa24*t2; d16 = T(1)/d16;
complex<T> d17 = spa34*t1*t2; d17 = T(1)/d17;
complex<T> d18 = spa23*t2*cube(spa14); d18 = T(1)/d18;
complex<T> d19 = spa23*spa34*t15*t2; d19 = T(1)/d19;
complex<T> d20 = spa23*t59; d20 = T(1)/d20;
complex<T> d21 = spa34*t2*square(spa24); d21 = T(1)/d21;
complex<T> d22 = spa24*t2*square(s23 - s234 + s34); d22 = T(1)/d22;
complex<T> d23 = spa23*t15*t2; d23 = T(1)/d23;
complex<T> d24 = t1*t59; d24 = T(1)/d24;
complex<T> d26 = spa12*spa34*t1; d26 = T(1)/d26;
complex<T> d27 = spa12*spa23*spa34*t15; d27 = T(1)/d27;
complex<T> d28 = spa23*t59*T(2); d28 = T(1)/d28;
complex<T> t9 = (s12 - s124)*(-(d17*t107) + d24*t18); 
complex<T> t23 = -(d19*t31); 
complex<T> t56 = -t75; 
complex<T> t60 = d1*spb12; 
complex<T> t62 = d2*t39; 
complex<T> t64 = d6*spa45; 
complex<T> t76 = d3*spb46; 
complex<T> t77 = d18*t16; 
complex<T> t80 = d14*T(2); 
complex<T> t83 = spb56*(d26*t107 + d27*t31*t74 + d25*t16*t93); 
complex<T> t87 = d11*spa15; 
complex<T> t88 = d16*t17; 
complex<T> t94 = d4*t46; 
complex<T> t98 = d22*t17; 
complex<T> d9 = spa34*t1*t48; d9 = T(1)/d9;
complex<T> d15 = spa24*t48*square(s234 - s34); d15 = T(1)/d15;
complex<T> d30 = spa23*t48*cube(spa14); d30 = T(1)/d30;
complex<T> d31 = spa23*spa34*t15*t48; d31 = T(1)/d31;
complex<T> d32 = spa24*t48*square(s23 - s234 + s34); d32 = T(1)/d32;
complex<T> d33 = t48*square(spa24); d33 = T(1)/d33;
complex<T> t10 = (-s123 + s23)*(t23*t74 + t77*t84); 
complex<T> t11 = s23*(d21*spa35*t44 + t75*t98); 
complex<T> t12 = (-s156 + s34)*(d21*spa35*t44 + t23*t74 + t77*t84 + t75*t98); 
complex<T> t27 = d5*t47*square(spa15) - spa13*t60*square(spa15) - d7*spa34*square(spa25)*square(spb24) + d8*spa23*square(spa45)*square(spb24) + t64*t85*T(2); 
complex<T> t54 = -(d33*s23*spa35*spb34*t44) + d32*s23*s34*t17*t75; 
complex<T> t78 = d15*t17; 
complex<T> t91 = t41*(d31*t31*t74 + d30*t16*t93); 
complex<T> t99 = spa15*t62; 
complex<T> t101 = d23*spb34*t31*t55 + s34*t77*t93; 
complex<T> t13 = -(d9*t107) - d12*spa35*spb16*t30 - d10*spa13*t16*t32 - d5*t16*t47 - d13*spa56*t45*t69 + d7*t17*t75 + t75*t78 + t75*t88 - d8*spa23*t14*square(spa45) - d14*spa25*t85*T(2) - t64*t85*T(2) + t32*t74*t87*T(2); 
complex<T> t68 = d12*spa35*spb16*t30 + d10*spa13*t16*t32 + d13*spa56*t45*t69 + t55*t76 + t25*t74*t87 + spa56*t93*t94 + t19*t99 + d17*t107*T(2); 
complex<T> t92 = spa25*t80*t85 + t56*(t78 + t88); 
complex<T> t106 = spa13*t16*t60 + t74*t76 + spa56*t84*t94 + spa13*spa45*t99; 
complex<T> co1 = -(d20*spb14*t18); 
complex<T> co2 = -(d28*s12*spb14*t18); 
complex<T> co3 = d9*t107*t49; 
complex<T> co4 = d29*spb12*spb23*t107; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t106*(*CI_users[0]->get_value(mc,ind,mu)) + t27*(*CI_users[1]->get_value(mc,ind,mu)) + t13*(*CI_users[2]->get_value(mc,ind,mu)) + t92*(*CI_users[3]->get_value(mc,ind,mu)) + t68*(*CI_users[4]->get_value(mc,ind,mu)) + t10*(*CI_users[6]->get_value(mc,ind,mu)) + co1*(*CI_users[7]->get_value(mc,ind,mu)) + t11*(*CI_users[8]->get_value(mc,ind,mu)) + t12*(*CI_users[9]->get_value(mc,ind,mu)) + t101*(*CI_users[10]->get_value(mc,ind,mu)) + t9*(*CI_users[11]->get_value(mc,ind,mu)) + t83*(*CI_users[12]->get_value(mc,ind,mu)) + co2*(*CI_users[13]->get_value(mc,ind,mu)) + co3*(*CI_users[14]->get_value(mc,ind,mu)) + co4*(*CI_users[15]->get_value(mc,ind,mu)) + t91*(*CI_users[16]->get_value(mc,ind,mu)) + t54*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qppqmmemep_SLC_wCI::\
C2q2g2l_qppqmmemep_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c2356));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c2356));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c14, c23, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c23, c14, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c356));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c14, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c23, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qppqmmemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, p, qm, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qppqmmemep SLC");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 19:07:52 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb16 = SPB(1,6);
complex<T> spb26 = SPB(2,6);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb36 = SPB(3,6);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb46 = SPB(4,6);
complex<T> spa14 = SPA(1,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa56 = SPA(5,6);
complex<T> spa24 = SPA(2,4);
complex<T> spb13 = SPB(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s123 = SS(1,2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s14 = -(spa14*spb14);
complex<T> s124 = SS(1,2,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s24 = -(spa24*spb24);
complex<T> s134 = SS(1,3,4);
complex<T> s234 = SS(2,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> t18 = spa14*spb16; 
complex<T> t19 = spa35*spb12; 
complex<T> t35 = s124*s23; 
complex<T> t51 = square(spb36); 
complex<T> t52 = square(spa24); 
complex<T> t54 = spa25*spb26; 
complex<T> t56 = spa34*spb16; 
complex<T> t61 = -((spa24*spb12 + spa34*spb13)*(spa25*spb26 + spa35*spb36)); 
complex<T> t62 = spa45*spb46; 
complex<T> t74 = spa12*(spa12*spb13 + spa24*spb34); 
complex<T> t75 = square(spa13*spb16 + spa23*spb26); 
complex<T> t77 = spa13*spb23; 
complex<T> t97 = spa24*spb13; 
complex<T> t129 = spa35*spb23; 
complex<T> t131 = -(spa24*spb26); 
complex<T> t133 = spa13*spb14; 
complex<T> t135 = s14*s23; 
complex<T> t136 = square(spb12); 
complex<T> t137 = spa14*T(2); 
complex<T> t138 = square(spa34); 
complex<T> t150 = (spa25*spa34 + spa24*spa35)*spb16; 
complex<T> t153 = -((spa24*spb12 + spa34*spb13)*spb36); 
complex<T> t163 = -s14 + s23 - s56; 
complex<T> t164 = -s14 - s23 + s56; 
complex<T> t166 = s14 - s23 - s56; 
complex<T> t167 = -(spa25*spb24) - spa35*spb34; 
complex<T> t168 = spa12*spb26 + spa13*spb36; 
complex<T> t173 = square(s123); 
complex<T> t174 = s13 - s24; 
complex<T> t175 = cube(spa13); 
complex<T> t176 = square(spa35); 
complex<T> t178 = -(spa13*spb16) - spa23*spb26; 
complex<T> t179 = square(spb46); 
complex<T> t180 = -(spa12*spb16) + spa24*spb46; 
complex<T> t187 = -s234 + s56; 
complex<T> t189 = -(spa15*spb13) + spa45*spb34; 
complex<T> t191 = -s123 + s23; 
complex<T> t200 = -(spb16*spb23) + spb12*spb36; 
complex<T> t237 = -(s56*T(2)); 
complex<T> t238 = spa12*spb24; 
complex<T> t239 = spa13*spb34; 
complex<T> t241 = spa24*spb36; 
complex<T> t247 = spa35*spb13; 
complex<T> t248 = spa56*spb16; 
complex<T> t268 = cube(spa12*spb13 + spa24*spb34); 
complex<T> t269 = spa12*spb56; 
complex<T> t274 = spa12*spa23; 
complex<T> t301 = spa45*spb14; 
complex<T> t303 = spa34*spb23; 
complex<T> t304 = spa14*spa23; 
complex<T> t334 = spa45*spb24; 
complex<T> t335 = spa25*spb12; 
complex<T> t339 = spa56*spb46; 
complex<T> t357 = spa23*spb36; 
complex<T> t374 = s14*s56; 
complex<T> t383 = spa13*spa34; 
complex<T> t434 = spa24*spa34; 
complex<T> d9 = s124; d9 = T(1)/d9;
complex<T> d10 = spa12*spb13 + spa24*spb34; d10 = T(1)/d10;
complex<T> d13 = spa56; d13 = T(1)/d13;
complex<T> d14 = spb56; d14 = T(1)/d14;
complex<T> d17 = spa23*spb56; d17 = T(1)/d17;
complex<T> d19 = (spa12*spb13 + spa24*spb34)*spb56; d19 = T(1)/d19;
complex<T> d30 = spa56*spb14; d30 = T(1)/d30;
complex<T> d32 = spb14*spb56; d32 = T(1)/d32;
complex<T> d37 = s124*s56*spa12*spb14; d37 = T(1)/d37;
complex<T> d38 = s234*s56*spa23*spb34; d38 = T(1)/d38;
complex<T> d40 = T(2); d40 = T(1)/d40;
complex<T> d41 = s234*spa23*(spa23*spb13 + spa24*spb14)*spb56; d41 = T(1)/d41;
complex<T> d52 = spa23*square(s234); d52 = T(1)/d52;
complex<T> d53 = s123*s234*spa23; d53 = T(1)/d53;
complex<T> d54 = s234; d54 = T(1)/d54;
complex<T> d57 = spa56*spb14*T(2); d57 = T(1)/d57;
complex<T> d58 = spa12; d58 = T(1)/d58;
complex<T> d64 = square(s124)*T(2); d64 = T(1)/d64;
complex<T> d65 = s123; d65 = T(1)/d65;
complex<T> t20 = -(d13*spa45); 
complex<T> t21 = -(spb23*t163); 
complex<T> t53 = square(t168); 
complex<T> t58 = -(spa13*(spa23*spb24 + t133)); 
complex<T> t64 = d10*(spa13*spb12 + spa34*spb24); 
complex<T> t72 = spa23*t163; 
complex<T> t73 = (t238 + t239)*T(2); 
complex<T> t76 = square(t129 + t334); 
complex<T> t78 = spb56*(spa23*spb24 + t133); 
complex<T> t90 = d65*s12*s56 - d40*(s12 - s34 + s56); 
complex<T> t111 = spb26*t153 + t131*t200 - d10*spa24*spb13*(spa13*spb12 + spa34*spb24)*square(spb36); 
complex<T> t121 = d14*spb46*(spa24*spb26 + t18) + d13*spa15*(t247 + t301); 
complex<T> t132 = square(t238 + t239); 
complex<T> t140 = (-s12 + s34 + t174)*T(3); 
complex<T> t143 = -(d13*spb14*square(spa45)) - d14*spa14*square(spb16); 
complex<T> t156 = -(d19*(spa34*spb14*spb36 + spb56*t335 - spb12*t357)); 
complex<T> t159 = (spa24*spb12 - spa34*spb13)*(spa35*spb36 + t54); 
complex<T> t165 = t129 + t334; 
complex<T> t169 = t247 + t335; 
complex<T> t181 = t131 - t18; 
complex<T> t188 = -(spa15*spb12) + t334; 
complex<T> t190 = t51*t52; 
complex<T> t203 = spa23*spa45*spb13*spb26 + spa25*spb23*t56; 
complex<T> t204 = -(spa14*spb24) - t77; 
complex<T> t226 = -(spa14*t136); 
complex<T> t254 = d10*s23; 
complex<T> t257 = -(d14*spb26); 
complex<T> t273 = cube(t238 + t239); 
complex<T> t296 = spa15*(t247 + t301); 
complex<T> t311 = d10*s56; 
complex<T> t323 = t74*T(2); 
complex<T> t328 = -(spa13*spb16*spb23) - spb24*t18; 
complex<T> t360 = spa56*t164; 
complex<T> t388 = -(d14*spb36); 
complex<T> t393 = t241*T(2); 
complex<T> t468 = -(spa13*t75); 
complex<T> t477 = t75*t77; 
complex<T> d1 = (t238 + t239)*t274*square(s123 - s56); d1 = T(1)/d1;
complex<T> d4 = (-s124 + s56)*spb56*t74; d4 = T(1)/d4;
complex<T> d5 = (-s124 + s56)*t269*square(spa12*spb13 + spa24*spb34); d5 = T(1)/d5;
complex<T> d8 = s124*(-s124 + s14)*spb56*t74; d8 = T(1)/d8;
complex<T> d12 = s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2); d12 = T(1)/d12;
complex<T> d15 = t238 + t239; d15 = T(1)/d15;
complex<T> d18 = (spa12*spb13 + spa24*spb34)*(t238 + t239)*(s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2)); d18 = T(1)/d18;
complex<T> d21 = (s23 - s234)*spa56*spb34*(t238 + t239); d21 = T(1)/d21;
complex<T> d22 = s234*spa56*spb34*(t238 + t239); d22 = T(1)/d22;
complex<T> d25 = spa56*(s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2)); d25 = T(1)/d25;
complex<T> d26 = square(s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2)); d26 = T(1)/d26;
complex<T> d27 = spb56*(s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2)); d27 = T(1)/d27;
complex<T> d31 = spb14*spb56*(t238 + t239); d31 = T(1)/d31;
complex<T> d35 = s234*spb34*t187*(t238 + t239); d35 = T(1)/d35;
complex<T> d39 = s56*spb14*spb34*t274; d39 = T(1)/d39;
complex<T> d42 = s124*spa56*spb14*(spa23*spb24 + t133); d42 = T(1)/d42;
complex<T> d43 = t268*t269; d43 = T(1)/d43;
complex<T> d44 = s124*t268*t269; d44 = T(1)/d44;
complex<T> d46 = spb56*t268; d46 = T(1)/d46;
complex<T> d49 = s124*spa56*(spa23*spb24 + t133); d49 = T(1)/d49;
complex<T> d50 = t268*t269*T(2); d50 = T(1)/d50;
complex<T> d55 = s123*spa12*(spa23*spb24 + t133); d55 = T(1)/d55;
complex<T> d56 = s124*t74; d56 = T(1)/d56;
complex<T> d59 = spa12*(t238 + t239); d59 = T(1)/d59;
complex<T> d60 = t269*T(2); d60 = T(1)/d60;
complex<T> d61 = (spa12*spb13 + spa24*spb34)*(spa23*spb24 + t133); d61 = T(1)/d61;
complex<T> d62 = spa56*spb14*(spa23*spb24 + t133); d62 = T(1)/d62;
complex<T> d63 = spb56*t74; d63 = T(1)/d63;
complex<T> d68 = T(2)*(s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2)); d68 = T(1)/d68;
complex<T> d69 = t74*(s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2)); d69 = T(1)/d69;
complex<T> d73 = spa12*t268*T(2); d73 = T(1)/d73;
complex<T> d75 = s124*spa56*(spa23*spb24 + t133)*T(2); d75 = T(1)/d75;
complex<T> t4 = -(d31*(spb14*spb26 + spb12*spb46)*t168) + d12*(d32*spb16*spb24*(-(spa34*spb36) - spa45*spb56 + t131)*t163 - d14*spb23*t168*t56 + d15*t163*(-(spa15*spb16) + t62)*t77 + spa45*(spb56*t334 + spb16*t77) + d32*spb23*spb46*(d15*spa13*spa24*spb12*spb46*t163 + d15*spa13*spa34*spb13*spb46*t163 - spb12*spb16*t304 + spa34*spb16*t163*T(2)) + t20*(t204*t301 - spa15*spb12*t166*T(2) + spb26*t360*T(2)) + spa15*spb23*(d15*spa45*spb56*t133 + t56)*T(4)) + d26*spb23*t61*(spa13*t164 + spb24*t304*T(2))*T(6); 
complex<T> t9 = d18*(d13*spb12*(-(spa14*t169) + spa24*t188 + spa34*t189)*(s134*spa15 - spa14*t339) + d17*(-(spa15*spb56*t164*t174) + spa14*spb46*t166*t174 - (spa45*spb56*t164 + t166*t18)*(t238 + t239))*t56 + spa14*spb12*t174*(spa15*spb16 + t62)*T(2)); 
complex<T> t10 = -(d14*spa34*(s123*spb36 - spa25*spb23*spb56)*(-(spa13*spb16*spb23) + spa34*spb24*spb36 + spa34*spb23*spb46 - spb24*t131 + spb12*t168)) + d30*t19*(-((spa12*spb13 + spa24*spb34)*(t129*t163 + spb26*t360)) + t174*(spa25*t21 + spb36*t360)) - t174*t303*(spa35*spb36 + t54)*T(2); 
complex<T> t11 = d14*t18*(spa23*spb24*(spa34*spb36 - t131) + spa34*spb36*t133 - t131*t133 - spa35*spb56*t164*T(2) + spa34*spb46*t166*T(2)) + d15*spa13*(d14*spb46*t18 + d13*spa15*t301)*(-(s12*t166) + s34*t166 - t166*t174 + square(s14) + square(s23) + square(s56) - s23*s56*T(2) - t135*T(2) - t374*T(2)) + d13*spa45*(-(s14*spa13*t169) + spb24*t169*t304 + spa23*spb12*(spa15*t163 + t137*t339)*T(2) - spa45*t133*t166*T(3)) - d12*t304*t61*(spb24*t164 + spb14*t77*T(2))*T(6); 
complex<T> t33 = spa25*spa34*t169 - spa23*spa45*(t301 + t335) + t176*t97 + t64*t97*square(spa25); 
complex<T> t59 = d43*s124; 
complex<T> t79 = d9*t181 + d10*t393; 
complex<T> t125 = -(d13*spa45*t335) + t156*t393 + t388*t56; 
complex<T> t142 = d12*t164; 
complex<T> t145 = d22*t76; 
complex<T> t149 = d39*(t247 + t301); 
complex<T> t154 = d8*t180; 
complex<T> t157 = -(d4*(spa24*spb26 + t18)); 
complex<T> t201 = -(d13*spa25*t129) + t257*t357; 
complex<T> t262 = d21*t165; 
complex<T> t277 = d26*(-s12 + s34 + t174); 
complex<T> t309 = t140*t61; 
complex<T> t310 = -(d1*s123); 
complex<T> t318 = d42*t176; 
complex<T> t320 = d35*t204; 
complex<T> t345 = -(d38*t165); 
complex<T> t373 = t53*t58; 
complex<T> t430 = t248*t254; 
complex<T> t447 = t176*t226; 
complex<T> t465 = t311*t97; 
complex<T> d2 = t132*t191*t269; d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spa23*t132*t269; d3 = T(1)/d3;
complex<T> d6 = spb56*t323*square(s124 - s56); d6 = T(1)/d6;
complex<T> d7 = spb56*t323*square(s124 - s14); d7 = T(1)/d7;
complex<T> d11 = (-s124 + s14)*spb56*t323; d11 = T(1)/d11;
complex<T> d16 = t274*t73*(s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2)); d16 = T(1)/d16;
complex<T> d20 = t323*(s14*t237 + s23*t237 + square(s14) + square(s23) + square(s56) - t135*T(2)); d20 = T(1)/d20;
complex<T> d23 = spa56*spb34*t73*square(s23 - s234); d23 = T(1)/d23;
complex<T> d24 = (s23 - s234)*s234*spa56*spb34*t73; d24 = T(1)/d24;
complex<T> d28 = spb14*t323; d28 = T(1)/d28;
complex<T> d29 = spb14*(spa12*spb13 + spa24*spb34)*t73; d29 = T(1)/d29;
complex<T> d33 = spa12*t73; d33 = T(1)/d33;
complex<T> d34 = spa56*spb34*t187*t73; d34 = T(1)/d34;
complex<T> d36 = spb34*t73*square(t187); d36 = T(1)/d36;
complex<T> d45 = spa23*t73*t78; d45 = T(1)/d45;
complex<T> d47 = spa23*t269*t273; d47 = T(1)/d47;
complex<T> d48 = t273*t274*t78; d48 = T(1)/d48;
complex<T> d51 = spa23*t269*t273*T(2); d51 = T(1)/d51;
complex<T> d66 = spa12*(t238 + t239)*t78; d66 = T(1)/d66;
complex<T> d67 = t269*t273*T(2); d67 = T(1)/d67;
complex<T> d70 = t274*t73*t78; d70 = T(1)/d70;
complex<T> d71 = (t238 + t239)*t274*t78; d71 = T(1)/d71;
complex<T> d72 = (spa23*spb24 + t133)*t274*t73; d72 = T(1)/d72;
complex<T> d74 = t273*t274*T(2); d74 = T(1)/d74;
complex<T> d76 = t73*t78; d76 = T(1)/d76;
complex<T> d77 = s234*spa56*t73; d77 = T(1)/d77;
complex<T> t55 = t136*t318 - t59*square(spa24)*square(spb36) + d44*square(spa14*spb13 + spa24*spb23)*square(t180); 
complex<T> t57 = t145 + d41*t138*square(spb16); 
complex<T> t139 = -(d47*t53); 
complex<T> t161 = d11*t79; 
complex<T> t245 = -(d48*t173); 
complex<T> t250 = d23*t167; 
complex<T> t253 = d51*s14; 
complex<T> t255 = -(d70*s34); 
complex<T> t264 = d3*t168; 
complex<T> t266 = t154*t181; 
complex<T> t313 = d7*spb26; 
complex<T> t316 = d24*t165; 
complex<T> t368 = t277*t61; 
complex<T> t433 = t262*t334; 
complex<T> t467 = t241*t430; 
complex<T> t5 = (d40*t163 + d9*t374)*t55 + t190*(d40*t163 + d9*t374)*t59 + d69*spa14*(d68*s23*t163*t309 + t301*t303*t357 + d40*s23*t61 + spa12*spb36*t169*t64*t97 + spa23*(spa34*spb12*spb36*t189 + spa35*spb13*spb36*t303 + spb13*t303*t54 + spb36*t189*t64*t97) + spa24*(spa34*spb12*spb36*t167 - spa13*spb12*spb36*t169 + spa35*spb14*spb36*t303 + spb14*t303*t54 - spb36*t167*t64*t97 + spa23*spb12*spb36*t188*T(2))) + d18*spa34*(-(s56*spb12*(spa25*spb23*spb46*t137 + spa15*spb36*t164)) + spb23*(s234*t163 + s56*t164)*t54*T(2) - spb23*spb36*t248*T(2)*(spa13*t164 + spb24*t304*T(2))) + d59*spa13*(-(d65*spa45*spb12*t168) - d26*spa14*spb23*t61*(t133*t166 - spb24*t72)*T(3) + d12*(spa45*spb12*(-(spa15*spb56*t164) + spa14*spb46*t166) + spb23*(spa15*t163 - spa45*(t238 + t239) + t137*t339)*t56 + d15*(-s12 - s13 + s24 + s34)*spa15*(spa24*spb12 + spa34*spb13)*spb46*t77 - (spa24*spb12 + spa34*spb13)*(spa14*spb24*(spa35*spb36 + t54) + t62*t77*T(3)))); 
complex<T> t6 = -(d16*spa13*t11) + d7*spa14*spb12*t131*t180 + spb12*t137*t266 - d11*spa14*spb12*t180*(d9*t181 + d10*t393) - t9*T(2) - d20*spa14*(t142*t309 - t166*(t20*t335 + t156*t393 + t388*t56) - t159*T(2) + t203*T(4) + (d13*spa25*t129 + t257*t357)*(spa24*spb12 + spa34*spb13 + t465*T(4)) + t467*T(8)); 
complex<T> t14 = -(d55*spa13*spa45*spb12*t178) + d56*t137*t19*t241 + d47*(d65*s23*s56 + d40*t166)*t373 + (d65*s23*s56 + d40*t166)*(t175*t179*t245 + t139*t58) + d61*(-(d60*t181*(t163*t178 + spa13*spb16*t237)) - d57*t163*t19*(t247 + t301) + spb26*(spa24*t19 - spa34*t301) + d58*spb36*t19*t304 + (spa15*spb12 - t129)*t56 - d59*spa13*spb12*(spa12*spa13*spb13 + spa13*spa24*spb34 - spa12*t174)*t62) + d64*(-(d62*t136*t176) - d63*square(spa24*spb26 + t18))*(s124*t163 + t374*T(2)) + d18*spb12*(-(s56*spa34*(spa25*spb23*spb46*t137 + spa15*spb36*t164)) - spa14*(s56*t164 + s124*t166)*t62*T(2) + spa15*spa35*spb56*t137*(spb13*t164 + spa24*spb14*spb23*T(2))); 
complex<T> t24 = d53*spa13*(d15*(spa24*spb12 + spa34*spb13)*t178*(t129 + t334) - spa45*spb12*t56) + (d40*(s12 - s34 - s56) + d54*s34*s56)*t57 + d52*spb16*t138*(t129 + t334)*T(2); 
complex<T> t27 = (s12 - s124)*(t136*t318 - t55 - t190*t59); 
complex<T> t29 = d37*s23*spa35*spb12*(-(spa14*spb16) + t131) + d50*t190*t35 - d67*spb23*t373 + d66*t477 + s23*(t145 + t149*t178 + d40*t175*t179*t245 + d40*t55 - d38*t129*t56 - d38*t334*t56 + d40*t139*t58); 
complex<T> t30 = d40*t35*t55 + d50*s23*t190*square(s124); 
complex<T> t39 = s14*(t175*t179*t245 + t139*t58); 
complex<T> t46 = s124*(-(d6*t138) + d5*t434)*t51 - spa34*spb36*t157*T(2) + spa14*spb12*(t161*t180 + spa24*t180*t313 - t266*T(2)); 
complex<T> t120 = -(d36*spa15*spb12*t328) - t250*t303*t334 + d34*spa15*spb12*(t129 + t334) - d24*t167*t303*(t129 + t334) + d35*spb16*(spa13*spb23 + spa14*spb24)*(t129 + t334)*T(2) - d21*t334*(t129 + t334)*T(2); 
complex<T> t128 = -((s156 - s34)*(t145 - t57)); 
complex<T> t246 = d40*t57; 
complex<T> t271 = -(t175*t179*t191*t245); 
complex<T> t322 = t255*t75; 
complex<T> t329 = s123*spb46*t264*t383 + d2*spa13*spb12*t53 + t310*t383*t62; 
complex<T> t392 = t135*t368; 
complex<T> t1 = -(d29*spb12*(d14*spb46*(spa24*spb26 + t18) + d13*spa15*(t247 + t301))) + t250*t303*t334 + d24*t167*t303*(t129 + t334) - d33*spa13*t4 - d2*spa13*spb12*t53 - d18*t10*T(2) + t145*T(2) + d21*t334*(t129 + t334)*T(2) + d28*(d25*spb23*t163*t33 + t20*t335 - d12*t143*t35 + d27*t111*t72 + d10*d13*spa25*t129*t97 - d10*t257*t357*t97 + spb23*t142*t150*T(2) + t392*T(6)); 
complex<T> t2 = d16*spa13*t11 + d29*spb12*(d14*spb46*(spa24*spb26 + t18) + d13*spa15*(t247 + t301)) - d34*spa15*spb12*(t129 + t334) - s123*spb46*t264*t383 + d33*spa13*t4 + d6*s124*t138*t51 - d5*s124*t434*t51 + d1*s123*t383*t62 - d36*spa15*spb12*(spb24*t18 + spb16*t77) + d18*t10*T(2) - t145*T(2) + d4*spa34*spb36*t181*T(2) + spb16*t320*(t129 + t334)*T(2) + t9*T(2) - d40*(t149*t178 + d37*t181*t19 - d38*(t129 + t334)*t56)*T(3) + d28*(d25*t21*t33 + d13*spa45*t335 + d12*t143*t35 - d27*t111*t72 - d10*d13*spa25*t129*t97 + d10*t257*t357*t97 - spb23*t142*t150*T(2) - t392*T(6)) + d20*spa14*(t142*t309 - t166*(t20*t335 + t156*t393 + t388*t56) - t159*T(2) + t203*T(4) + (d13*spa25*t129 + t257*t357)*(spa24*spb12 + spa34*spb13 + t465*T(4)) + t467*T(8)); 
complex<T> t25 = s123*(t253*t373 + d40*t39); 
complex<T> t28 = -(d50*s124*s14*t190) + t253*t373 + d40*t39 + d49*t447 - d40*s14*t55; 
complex<T> t31 = d73*s124*spa56*t190 - s56*t246 + d74*spa56*t373 + d72*spa56*t468 - d40*s56*(t175*t179*t245 + t55 + t139*t58); 
complex<T> t347 = s12*t246; 
complex<T> t442 = spa13*t322; 
complex<T> t23 = s34*t246 + t442; 
complex<T> t26 = -(d46*s124*spb12*t190) + t347 + s12*t55 + s12*(d37*spa35*spb12*(-(spa14*spb16) + t131) + t149*t178 - d38*(t129 + t334)*t56) + d45*spa13*spb12*t75; 
complex<T> co1 = d71*t468*t90; 
complex<T> co2 = s234*t347; 
complex<T> co3 = d75*s12*t447; 
complex<T> co4 = -(d76*spb12*t477); 
complex<T> co5 = s123*t442; 
complex<T> co6 = -(d77*s23*spa34*t76); 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t329*(*CI_users[0]->get_value(mc,ind,mu)) + t46*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t120*(*CI_users[4]->get_value(mc,ind,mu)) + t2*(*CI_users[5]->get_value(mc,ind,mu)) + t26*(*CI_users[6]->get_value(mc,ind,mu)) + t271*(*CI_users[7]->get_value(mc,ind,mu)) + t28*(*CI_users[8]->get_value(mc,ind,mu)) + t24*(*CI_users[9]->get_value(mc,ind,mu)) + t14*(*CI_users[10]->get_value(mc,ind,mu)) + t29*(*CI_users[11]->get_value(mc,ind,mu)) + t128*(*CI_users[12]->get_value(mc,ind,mu)) + t5*(*CI_users[13]->get_value(mc,ind,mu)) + t23*(*CI_users[14]->get_value(mc,ind,mu)) + co1*(*CI_users[15]->get_value(mc,ind,mu)) + t27*(*CI_users[16]->get_value(mc,ind,mu)) + t31*(*CI_users[17]->get_value(mc,ind,mu)) + co2*(*CI_users[18]->get_value(mc,ind,mu)) + co3*(*CI_users[19]->get_value(mc,ind,mu)) + co4*(*CI_users[20]->get_value(mc,ind,mu)) + t30*(*CI_users[21]->get_value(mc,ind,mu)) + t25*(*CI_users[22]->get_value(mc,ind,mu)) + co5*(*CI_users[23]->get_value(mc,ind,mu)) + co6*(*CI_users[24]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpmqmpemep_SLC_wCI::\
C2q2g2l_qpmqmpemep_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c2356));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c14, c23, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c23, c14, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c14, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c356));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c23, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpmqmpemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, qm, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpmqmpemep SLC");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 19:25:27 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa56 = SPA(5,6);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb46 = SPB(4,6);
complex<T> spa13 = SPA(1,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb16 = SPB(1,6);
complex<T> spb56 = SPB(5,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb34 = SPB(3,4);
complex<T> spb36 = SPB(3,6);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s56 = -(spa56*spb56);
complex<T> s123 = SS(1,2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s124 = SS(1,2,4);
complex<T> s234 = SS(2,3,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s24 = -(spa24*spb24);
complex<T> s134 = SS(1,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> t18 = spa35*spb34; 
complex<T> t19 = spa23*spb16; 
complex<T> t35 = -(s12*s234); 
complex<T> t51 = square(spa15); 
complex<T> t52 = square(spa15*spb13 + spa25*spb23); 
complex<T> t53 = square(spb24); 
complex<T> t54 = spa25*spb26; 
complex<T> t55 = spa35*spb14; 
complex<T> t62 = spa45*spb46; 
complex<T> t74 = spb23*(spa13*spb23 + spa14*spb24); 
complex<T> t75 = square(spa25*spb12 + spa35*spb13); 
complex<T> t79 = spb56*T(2); 
complex<T> t131 = -(spa34*spb13); 
complex<T> t133 = s12*s34; 
complex<T> t134 = square(spb16); 
complex<T> t136 = square(spa23); 
complex<T> t137 = spa12*spb26; 
complex<T> t138 = spb34*T(2); 
complex<T> t139 = square(spb14); 
complex<T> t164 = s12 - s34 - s56; 
complex<T> t165 = -s12 - s34 + s56; 
complex<T> t167 = -s12 + s34 - s56; 
complex<T> t168 = spa24*spb12 + spa34*spb13; 
complex<T> t169 = -(spa14*spb16) - spa24*spb26; 
complex<T> t170 = -(spa15*spb13) - spa25*spb23; 
complex<T> t175 = s13 - s24; 
complex<T> t177 = -(spa13*spb14) - spa23*spb24; 
complex<T> t178 = spa35*spb23 + spa45*spb24; 
complex<T> t180 = s12 - s124; 
complex<T> t184 = -s124 + s56; 
complex<T> t185 = -s234 + s56; 
complex<T> t187 = spa13*spb36 + spa14*spb46; 
complex<T> t189 = -(spa13*spb16) + spa34*spb46; 
complex<T> t199 = spb16*spb24 + spb14*spb26; 
complex<T> t240 = -(s56*T(2)); 
complex<T> t241 = spa15*spb24; 
complex<T> t242 = -(spa34*spb46); 
complex<T> t245 = spb12*spb34; 
complex<T> t267 = cube(spa13*spb23 + spa14*spb24); 
complex<T> t268 = spb23*T(2); 
complex<T> t270 = spa12*spb13; 
complex<T> t281 = spb24*square(spa15); 
complex<T> t290 = spa14*spb13; 
complex<T> t291 = spa24*spb23; 
complex<T> t295 = spa12*spa35; 
complex<T> t297 = spa25*spb24; 
complex<T> t299 = spa45*spb56; 
complex<T> t329 = spa24*spb46; 
complex<T> t354 = spa12*spb14; 
complex<T> t355 = spa15*spb12; 
complex<T> t366 = s56*T(2); 
complex<T> t381 = spb13*spb14; 
complex<T> t382 = spa23*spb26; 
complex<T> t387 = spa35*spb56; 
complex<T> t410 = -(s56*T(4)); 
complex<T> d2 = spa56; d2 = T(1)/d2;
complex<T> d3 = spb56; d3 = T(1)/d3;
complex<T> d5 = spa13*spb23 + spa14*spb24; d5 = T(1)/d5;
complex<T> d15 = spa34*spb56; d15 = T(1)/d15;
complex<T> d19 = spa34*spa56; d19 = T(1)/d19;
complex<T> d31 = s234; d31 = T(1)/d31;
complex<T> d33 = spa56*(spa13*spb23 + spa14*spb24); d33 = T(1)/d33;
complex<T> d35 = spa56*spb12; d35 = T(1)/d35;
complex<T> d37 = s234*s56*spa34*spb23; d37 = T(1)/d37;
complex<T> d38 = s124*s56*spa14*spb12; d38 = T(1)/d38;
complex<T> d39 = s56*spa14*spa34*spb12*spb23; d39 = T(1)/d39;
complex<T> d40 = T(2); d40 = T(1)/d40;
complex<T> d49 = s124*spa56*spb12*(spa13*spb12 + spa34*spb24); d49 = T(1)/d49;
complex<T> d54 = s123; d54 = T(1)/d54;
complex<T> d59 = spb12*square(s124); d59 = T(1)/d59;
complex<T> d60 = s123*s124*spb12; d60 = T(1)/d60;
complex<T> d61 = s124; d61 = T(1)/d61;
complex<T> d67 = square(s234)*T(2); d67 = T(1)/d67;
complex<T> d68 = spb23; d68 = T(1)/d68;
complex<T> t20 = -(d3*spb46); 
complex<T> t21 = -(spa12*t164); 
complex<T> t56 = -(spb13*t168); 
complex<T> t61 = -((spa15*spb16 + spa25*spb26)*t177); 
complex<T> t64 = d5*(spa23*spb13 + spa24*spb14); 
complex<T> t72 = spb12*t164; 
complex<T> t73 = square(t290 + t291); 
complex<T> t78 = spa56*t168; 
complex<T> t93 = d40*(s14 - s23 - s56) + d54*s23*s56; 
complex<T> t130 = -t297; 
complex<T> t132 = cube(t290 + t291); 
complex<T> t140 = (s14 - s23 + t175)*T(3); 
complex<T> t144 = t168*t52; 
complex<T> t150 = d39*(spa25*spb12 + spa35*spb13); 
complex<T> t154 = d33*(spa15*spa34*spb14 - spa23*t355 + spa56*t382); 
complex<T> t161 = (spa13*spb14 - spa23*spb24)*(spa15*spb16 + t54); 
complex<T> t166 = -(spa12*spb16) + t329; 
complex<T> t171 = -(spa13*spb16) - t382; 
complex<T> t186 = -t290 - t291; 
complex<T> t188 = spa23*spb36 + t329; 
complex<T> t190 = t51*t53; 
complex<T> t196 = -(spa15*spa23) + t295; 
complex<T> t203 = spa13*spa25*spb12*spb46 - t137*t55; 
complex<T> t247 = d5*spa13; 
complex<T> t253 = spb36*t165; 
complex<T> t259 = -(d2*spa25); 
complex<T> t271 = spb12*(t290 + t291); 
complex<T> t273 = t74*T(2); 
complex<T> t278 = spa24*t18 + spa35*t270; 
complex<T> t287 = spb36*(spa13*spb16 + t242); 
complex<T> t332 = -(spb56*t165); 
complex<T> t336 = d40*s34; 
complex<T> t341 = spa45*t170; 
complex<T> t357 = d5*s12; 
complex<T> t363 = -(spa25*t177); 
complex<T> t373 = d2*spa15; 
complex<T> t380 = t241*T(2); 
complex<T> t433 = s23*t136; 
complex<T> t460 = -(t241*T(8)); 
complex<T> t462 = -(spb13*t75); 
complex<T> t468 = t270*t75; 
complex<T> d4 = spa34*(spa13*spb23 + spa14*spb24)*(t290 + t291)*T(2); d4 = T(1)/d4;
complex<T> d6 = spa56*(s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2)); d6 = T(1)/d6;
complex<T> d7 = square(s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2)); d7 = T(1)/d7;
complex<T> d8 = s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2); d8 = T(1)/d8;
complex<T> d9 = spb56*(s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2)); d9 = T(1)/d9;
complex<T> d11 = spa14*(t290 + t291)*t79*square(t180); d11 = T(1)/d11;
complex<T> d12 = s124*spa14*t180*(t290 + t291)*t79; d12 = T(1)/d12;
complex<T> d13 = spa14*spb56*t180*(t290 + t291); d13 = T(1)/d13;
complex<T> d14 = s124*spa14*spb56*(t290 + t291); d14 = T(1)/d14;
complex<T> d16 = (spa13*spb23 + spa14*spb24)*(t290 + t291)*(s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2)); d16 = T(1)/d16;
complex<T> d17 = spa34*spa56*(t290 + t291); d17 = T(1)/d17;
complex<T> d18 = t290 + t291; d18 = T(1)/d18;
complex<T> d20 = t268*(t290 + t291); d20 = T(1)/d20;
complex<T> d23 = spa14*(t290 + t291)*square(t184)*T(2); d23 = T(1)/d23;
complex<T> d24 = s124*spa14*t184*(t290 + t291); d24 = T(1)/d24;
complex<T> d25 = spa14*t184*(t290 + t291)*t79; d25 = T(1)/d25;
complex<T> d26 = spa56*spb23*t185*square(spa13*spb23 + spa14*spb24); d26 = T(1)/d26;
complex<T> d29 = spa56*t185*t74; d29 = T(1)/d29;
complex<T> d30 = s234*(-s234 + s34)*spa56*t74; d30 = T(1)/d30;
complex<T> d45 = spa56*t267*t268; d45 = T(1)/d45;
complex<T> d46 = spa56*spb23*t267; d46 = T(1)/d46;
complex<T> d47 = s234*spa56*spb23*t267; d47 = T(1)/d47;
complex<T> d48 = s234*spa34*spb56*t168; d48 = T(1)/d48;
complex<T> d51 = s124*spb56*(t290 + t291); d51 = T(1)/d51;
complex<T> d52 = T(2)*(s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2)); d52 = T(1)/d52;
complex<T> d53 = t74*(s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2)); d53 = T(1)/d53;
complex<T> d55 = spb23*(t290 + t291); d55 = T(1)/d55;
complex<T> d63 = s234*t74; d63 = T(1)/d63;
complex<T> d64 = s123*spb23*t168; d64 = T(1)/d64;
complex<T> d65 = spa56*t74; d65 = T(1)/d65;
complex<T> d66 = spa34*spb56*t168; d66 = T(1)/d66;
complex<T> d69 = spa56*t268; d69 = T(1)/d69;
complex<T> d70 = spa34*t79; d70 = T(1)/d70;
complex<T> d71 = (spa13*spb23 + spa14*spb24)*t168; d71 = T(1)/d71;
complex<T> d74 = t267*t268; d74 = T(1)/d74;
complex<T> d76 = s234*t168*t79; d76 = T(1)/d76;
complex<T> d77 = s124*(t290 + t291)*t79; d77 = T(1)/d77;
complex<T> t4 = d17*(-(spa25*spa34) + spa23*spa45)*t170 + d8*(-(d19*spa24*spa35*(-(spa15*spb14) + spa56*spb46 + t130)*t164) + spb46*(spa35*t270 - spa56*t329) - d2*spa12*t170*t55 + d18*spb13*t21*(-(spa35*spb36) + t62) + d19*spa12*spa45*(d18*spa45*spb13*t164*t177 - spa23*spa35*t245 + spa35*spb14*t164*T(2)) - d3*spb46*(spa24*spa34*spb34*spb46 + spa34*spb46*t270 - spa23*spb36*t167*T(2) - spa25*t332*T(2)) - spa12*spb36*(d18*spa56*spb46*t131 + t55)*T(4)) + d7*spa12*t61*(spb13*t165 + spa24*t245*T(2))*T(6); 
complex<T> t8 = d16*(d2*spb14*(-(s123*spa15) + spa56*t137)*(spa12*spa35*spb13 - spa15*spa24*spb14 + spa12*spa45*spb14 + spa24*t130 + spa23*t170) + d15*t19*(t175*(t137*t164 + spa15*t332) - (spa13*spb23 + spa14*spb24)*(spb16*t21 + spa25*t332)) + t175*t354*(spa15*spb16 + t54)*T(2)); 
complex<T> t10 = -(d3*spa23*(spb34*t171 - spb14*t187 - spb24*t188)*(s134*spb36 + spb34*t299)) + d35*(spa45*spb34*t167*t175 + spa56*spb46*t165*t186 - t167*t18*t186 + spa56*t175*t253)*t55 - spa23*t138*t175*(spa35*spb36 + t62); 
complex<T> t11 = d18*spb13*(d2*spa45*t18 - spa34*spb36*t20)*(s14*t167 - s23*t167 - t167*t175 - s12*t366 - s34*t366 + square(s12) + square(s34) + square(s56) - t133*T(2)) + d2*t18*((-(spa15*spb14) + t130)*t168 + spa56*spb16*t165*T(2) + spa45*spb14*t167*T(2)) + t20*(-(s34*spb13*t171) + spa24*t171*t245 + spa23*spb12*(-(spb36*t164*T(2)) + spa45*spb34*t79*T(2)) - spa34*spb13*spb46*t167*T(3)) - d8*t245*t61*(spa24*t165 + spa34*t270*T(2))*T(6); 
complex<T> t33 = spb14*spb26*t171 + spb12*spb46*(spa34*spb46 - t382) - spa13*spb24*(t134 + t64*square(spb26)); 
complex<T> t60 = d46*s234; 
complex<T> t76 = square(t166); 
complex<T> t104 = t130*t196 + spa23*spb13*t247*t281 + spa24*spb14*t247*t281 + spa15*t363; 
complex<T> t109 = spa15*t253 - spa45*spb34*t137*T(2); 
complex<T> t122 = d2*spa45*(t130 - t18) + d3*t287; 
complex<T> t125 = t154*t380 - d3*spb46*t382 + t373*t55; 
complex<T> t141 = d8*t165; 
complex<T> t143 = spa35*t166; 
complex<T> t152 = d30*t178; 
complex<T> t159 = -(d29*(t18 + t297)); 
complex<T> t179 = t130 - t18; 
complex<T> t274 = d7*(s14 - s23 + t175); 
complex<T> t301 = d11*t169; 
complex<T> t304 = t140*t61; 
complex<T> t339 = d13*t166; 
complex<T> t367 = t52*t56; 
complex<T> t376 = d38*t166; 
complex<T> t426 = spb24*t247; 
complex<T> t447 = t357*t387; 
complex<T> d1 = (s12 - s123)*spa56*spb23*t73; d1 = T(1)/d1;
complex<T> d10 = spa34*t273; d10 = T(1)/d10;
complex<T> d21 = (-s123 + s56)*spa56*spb12*spb23*t73; d21 = T(1)/d21;
complex<T> d22 = spb23*t271*square(s123 - s56); d22 = T(1)/d22;
complex<T> d27 = spa56*t273*square(t185); d27 = T(1)/d27;
complex<T> d28 = spa56*t273*square(s234 - s34); d28 = T(1)/d28;
complex<T> d32 = (-s234 + s34)*spa56*t273; d32 = T(1)/d32;
complex<T> d34 = t273*(s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2)); d34 = T(1)/d34;
complex<T> d36 = t268*t271*(s12*t240 + s34*t240 + square(s12) + square(s34) + square(s56) - t133*T(2)); d36 = T(1)/d36;
complex<T> d41 = spb23*(t290 + t291)*t78; d41 = T(1)/d41;
complex<T> d42 = spa56*t132*t268; d42 = T(1)/d42;
complex<T> d43 = spb12*spb23*t132*t78; d43 = T(1)/d43;
complex<T> d44 = spa56*spb12*spb23*t132; d44 = T(1)/d44;
complex<T> d50 = t268*t271*t78; d50 = T(1)/d50;
complex<T> d56 = spb23*t271*t78; d56 = T(1)/d56;
complex<T> d57 = t271*t78*T(2); d57 = T(1)/d57;
complex<T> d58 = spa56*spb12*t132; d58 = T(1)/d58;
complex<T> d62 = spa56*spb12*t132*t268; d62 = T(1)/d62;
complex<T> d72 = t168*t268*t271; d72 = T(1)/d72;
complex<T> d73 = spb12*t132*t268; d73 = T(1)/d73;
complex<T> d75 = (t290 + t291)*t78*T(2); d75 = T(1)/d75;
complex<T> t13 = -(d64*spa23*spb13*(spa25*spb12 + spa35*spb13)*spb46) + d63*t138*t19*t241 + d44*(d54*s12*s56 + d40*t167)*t367 + d71*(d70*t164*t189*t19 + spa25*(spa34*spb14*spb46 + spb24*t19) + d68*spa15*t19*t245 + d69*t179*(spa25*spb12*t164 + spa35*spb13*(t164 + t366)) + (-(spa12*spb16) + spa23*spb36)*t55 - d55*spa23*spb13*(spa13*spb13*spb23 + spa14*spb13*spb24 - spb23*t175)*t62) + (d54*s12*s56 + d40*t167)*(d44*spb13*t168*t52 - d43*cube(spb13)*square(s123)*square(spa45)) + d67*(s234*t164 + s34*t366)*(-(d66*t134*t136) - d65*square(t18 + t297)) - d16*spa23*(s56*spb14*t109 + t138*(s56*t165 + s234*t167)*t62 - spa56*spb16*spb34*spb36*T(2)*(spa13*t165 + spa12*spa34*spb24*T(2))); 
complex<T> t22 = d43*(-s123 + s23)*cube(spb13)*square(s123)*square(spa45); 
complex<T> t57 = -(d48*t134*t136) + t60*square(spa15)*square(spb24) - d47*square(spa12*spb24 + spa13*spb34)*square(t178); 
complex<T> t80 = d31*t179 + d5*t380; 
complex<T> t250 = d14*t76; 
complex<T> t256 = -(d50*s14); 
complex<T> t258 = d62*s34; 
complex<T> t262 = d44*t144; 
complex<T> t289 = d21*s123*t341*t381 + d1*spa23*spb13*t52 - d22*s123*t381*t62; 
complex<T> t309 = d28*spb34; 
complex<T> t313 = t141*t199; 
complex<T> t345 = t274*t61; 
complex<T> t424 = t329*t339; 
complex<T> t442 = -(spb14*t76); 
complex<T> t5 = d36*spb13*t11 + spa23*t178*t297*t309 + d32*spa23*spb34*t178*(d31*t179 + d5*t380) + d16*t10*T(2) - spa23*spb34*t152*t179*T(2) + d34*spb34*(t141*t304 + (d3*spb16*t137 + t259*t355)*(t177 + t410*t426) + t447*t460 - t167*(t154*t380 + t20*t382 + t373*t55) - t161*T(2) + t203*T(4)); 
complex<T> t6 = -((d31*s34*s56 + d40*t164)*t57) + (d31*s34*s56 + d40*t164)*t190*t60 - d53*spb34*(d52*s12*t164*t304 - t354*(spa15*spa34*spb16*spb24 + t242*t355 + spa34*spb24*t54 + spa13*spb12*(spa15*spb16 + t54)) + d40*s12*t61 + spa13*(spb24*t169 + spb23*t171 + spb12*t187)*t241*t64 - spa15*spa23*(spb14*spb24*t169 + spb13*spb24*t171 - spb12*spb14*t187 - spb12*spb24*t188*T(2))) + d16*spb14*(-(s56*spa23*t109) + spa12*(s124*t164 + s56*t165)*t54*T(2) - spa15*t295*t79*(spb13*t165 + spa24*t245*T(2))) - d55*spb13*(d54*spa23*spb46*t170 + d7*spa12*spb34*t61*(t131*t167 + spa24*t72)*T(3) + d8*(-(spa23*spb46*(spa45*spb34*t167 + spa56*t253)) + d18*(-s13 + s14 - s23 + s24)*spa45*(spa13*spb14 + spa23*spb24)*spb36*t270 + spa12*t55*(-(spb36*t164) + spb46*t186 + spa45*spb34*t79) - t177*(spa15*spa24*spb16*spb34 + spa24*spb34*t54 + t270*t62*T(3)))); 
complex<T> t23 = d40*t35*t57 + d45*s12*t190*square(s234); 
complex<T> t24 = d45*s234*s34*t190 + t258*t367 + d44*spb13*t168*t336*t52 - t336*t57 - d43*t336*cube(spb13)*square(s123)*square(spa45); 
complex<T> t25 = (-s156 + s34)*(d48*t134*t136 + t57 - t190*t60); 
complex<T> t58 = t250 + d49*t139*square(spa35); 
complex<T> t77 = spb13*t262 - d43*cube(spb13)*square(s123)*square(spa45); 
complex<T> t121 = -(d25*spa23*spb36*t166) + d23*spa23*spb36*t278 + d12*t166*t169*t354 + t301*t329*t354 - d24*t143*(spa24*spb34 + t270)*T(2) - t424*T(2); 
complex<T> t163 = d32*t80; 
complex<T> t315 = t256*t75; 
complex<T> t389 = t133*t345; 
complex<T> t441 = t295*t313; 
complex<T> t1 = -(d36*spb13*t11) + d25*spa23*spb36*t166 - d23*spa23*spb36*t278 - d26*s234*spb14*t281 + d4*spa23*(d2*spa45*(t130 - t18) + d3*t287) - d21*s123*t341*t381 - d20*spb13*t4 + d27*s234*t139*t51 + d22*s123*t381*t62 - d16*t10*T(2) + d29*spa15*spb14*t179*T(2) - t250*T(2) + d24*t143*(spa24*spb34 + t270)*T(2) - t8*T(2) - d40*(t150*t189 + d37*t179*t19 + t376*t55)*T(3) - d34*spb34*(t141*t304 + (d3*spb16*t137 + t259*t355)*(t177 + t410*t426) + t447*t460 - t167*(t154*t380 + t20*t382 + t373*t55) - t161*T(2) + t203*T(4)) + d10*(d9*spa12*t164*t33 + t20*t382 - d3*spb16*t137*t426 + t259*t355*t426 + d6*t104*t72 + d2*d8*spb34*t35*square(spa35) + d3*d8*spa34*t35*square(spb46) - t441*T(2) + t389*T(6)); 
complex<T> t2 = d2*d4*spa23*spa45*(-t130 + t18) - d3*d4*spa23*t287 - d12*t166*t169*t354 - t301*t329*t354 + d20*spb13*t4 - d1*spa23*spb13*t52 + t250*T(2) + t424*T(2) + t8*T(2) + d10*(d9*t21*t33 + d3*spb46*t382 + d3*spb16*t137*t426 - t259*t355*t426 - d6*t104*t72 + d2*d8*s12*s234*spb34*square(spa35) + d3*d8*s12*s234*spa34*square(spb46) + t441*T(2) - t389*T(6)); 
complex<T> t28 = d60*spb13*(-(d18*(spa25*spb12 + spa35*spb13)*t166*t177) + spa23*spb46*t55) + (d61*s14*s56 - d40*(s14 - s23 + s56))*t58 - d59*t139*t143*T(2); 
complex<T> t29 = d37*s12*spa23*spb16*(-(spa35*spb34) + t130) + s12*t150*t189 + d45*s12*s234*t190 + d42*spa12*t367 + d41*t468 + s12*t376*t55 + s12*t58 + d40*s12*(-(d44*spb13*t168*t52) - t57 + d43*cube(spb13)*square(s123)*square(spa45)); 
complex<T> t32 = spb56*(d74*s234*t190 + d73*t367 + d72*t462) + d40*s56*(-(d44*spb13*t168*t52) + t57 - t58 + d43*cube(spb13)*square(s123)*square(spa45)); 
complex<T> t39 = s23*t58; 
complex<T> t46 = -(spa23*spb34*t163*t178) + spa23*t138*t152*t179 + d26*s234*spb14*t281 + spa23*t130*t178*t309 - d27*s234*t139*t51 - spa15*spb14*t159*T(2); 
complex<T> t324 = s123*(t258*t367 + t336*t77); 
complex<T> t348 = t180*(t250 - t58); 
complex<T> t451 = spb13*t315; 
complex<T> t27 = d51*t442 + t451 - d40*s14*t58; 
complex<T> t364 = d40*t39; 
complex<T> t30 = d37*s23*spa23*spb16*(-(spa35*spb34) + t130) + s23*t150*t189 + t364 + d48*t134*t433 + d44*s23*spb13*t168*t52 + d58*spa23*spb13*t168*t52 + s23*t376*t55 + d57*spa23*spb13*t75 - d43*s23*cube(spb13)*square(s123)*square(spa45); 
complex<T> co1 = d56*t462*t93; 
complex<T> co2 = -(d75*spa23*t468); 
complex<T> co3 = -(d76*spb34*t134*t433); 
complex<T> co4 = s124*t364; 
complex<T> co5 = d77*s12*t442; 
complex<T> co6 = s123*t451; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t289*(*CI_users[1]->get_value(mc,ind,mu)) + t121*(*CI_users[2]->get_value(mc,ind,mu)) + t46*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + t1*(*CI_users[5]->get_value(mc,ind,mu)) + t29*(*CI_users[6]->get_value(mc,ind,mu)) + t22*(*CI_users[7]->get_value(mc,ind,mu)) + t27*(*CI_users[8]->get_value(mc,ind,mu)) + t6*(*CI_users[9]->get_value(mc,ind,mu)) + co1*(*CI_users[10]->get_value(mc,ind,mu)) + t30*(*CI_users[11]->get_value(mc,ind,mu)) + t25*(*CI_users[12]->get_value(mc,ind,mu)) + t28*(*CI_users[13]->get_value(mc,ind,mu)) + t24*(*CI_users[14]->get_value(mc,ind,mu)) + t13*(*CI_users[15]->get_value(mc,ind,mu)) + t348*(*CI_users[16]->get_value(mc,ind,mu)) + t32*(*CI_users[17]->get_value(mc,ind,mu)) + co2*(*CI_users[18]->get_value(mc,ind,mu)) + t23*(*CI_users[19]->get_value(mc,ind,mu)) + co3*(*CI_users[20]->get_value(mc,ind,mu)) + co4*(*CI_users[21]->get_value(mc,ind,mu)) + co5*(*CI_users[22]->get_value(mc,ind,mu)) + co6*(*CI_users[23]->get_value(mc,ind,mu)) + t324*(*CI_users[24]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpmqmmemep_SLC_wCI::\
C2q2g2l_qpmqmmemep_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c2356));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c2356));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c356));
CI_users.push_back(new Cached_Box_Integral_User(c4, c12, c3, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpmqmmemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, m, qm, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpmqmmemep SLC");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 19:26:07 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb16 = SPB(1,6);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb36 = SPB(3,6);
complex<T> spb46 = SPB(4,6);
complex<T> spa23 = SPA(2,3);
complex<T> spb56 = SPB(5,6);
complex<T> spb24 = SPB(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb26 = SPB(2,6);
complex<T> spa35 = SPA(3,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s56 = -(spa56*spb56);
complex<T> s14 = -(spa14*spb14);
complex<T> s124 = SS(1,2,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s156 = SS(1,5,6);
complex<T> s234 = SS(2,3,4);
complex<T> t1 = spb12*spb34; 
complex<T> t2 = spb23*spb56; 
complex<T> t15 = square(spa24); 
complex<T> t16 = square(spb34); 
complex<T> t17 = square(spb16); 
complex<T> t18 = square(spb26); 
complex<T> t19 = square(spb36); 
complex<T> t28 = s12 - s124; 
complex<T> t29 = -s124 + s56; 
complex<T> t30 = square(spb13); 
complex<T> t31 = spa23*spb12 - spa34*spb14; 
complex<T> t33 = -(spb14*spb36) - spb13*spb46; 
complex<T> t35 = s12 - s124 + s14; 
complex<T> t41 = s123*s124 - s12*s56; 
complex<T> t42 = square(spa35); 
complex<T> t43 = square(spa45); 
complex<T> t44 = spa24*spb12 + spa34*spb13; 
complex<T> t45 = square(spb46); 
complex<T> t46 = -(spb14*spb26) - spb12*spb46; 
complex<T> t47 = spa24*spb14; 
complex<T> t49 = s123*s234 - s23*s56; 
complex<T> t57 = spb24*spb56; 
complex<T> t67 = cube(spb13); 
complex<T> t73 = spb13*spb16; 
complex<T> t79 = spb36*spb46; 
complex<T> t94 = spb14*spb56; 
complex<T> t101 = spa24*spb16; 
complex<T> d24 = spb12*spb56*cube(spb34); d24 = T(1)/d24;
complex<T> d29 = spb12*spb23*cube(spb34); d29 = T(1)/d29;
complex<T> t20 = spb24*t2; 
complex<T> t23 = -(spb14*t28); 
complex<T> t24 = -(t31*T(2)); 
complex<T> t25 = -(spa14*t46); 
complex<T> t38 = spb12*t2; 
complex<T> t53 = -t73; 
complex<T> t74 = spb14*t15; 
complex<T> t93 = spb13*t17; 
complex<T> t100 = spb13*t19; 
complex<T> d2 = (s12 - s123)*t16*t2; d2 = T(1)/d2;
complex<T> d3 = t16*t2*t28; d3 = T(1)/d3;
complex<T> d4 = spb34*t28*t57; d4 = T(1)/d4;
complex<T> d5 = spb34*t57*square(t28)*T(2); d5 = T(1)/d5;
complex<T> d6 = (-s123 + s56)*spb23*t1; d6 = T(1)/d6;
complex<T> d8 = spb23*t1*square(s123 - s56); d8 = T(1)/d8;
complex<T> d9 = spb14*spb23*t1*t29; d9 = T(1)/d9;
complex<T> d13 = spb14*t1*t2*T(2); d13 = T(1)/d13;
complex<T> d14 = spb14*t1*t2*t29; d14 = T(1)/d14;
complex<T> d16 = spb14*spb23*t1*square(t29)*T(2); d16 = T(1)/d16;
complex<T> d17 = spb14*t1*t2; d17 = T(1)/d17;
complex<T> d21 = t2*square(spb24); d21 = T(1)/d21;
complex<T> d22 = t1*t94; d22 = T(1)/d22;
complex<T> d23 = t1*t57; d23 = T(1)/d23;
complex<T> d25 = spb12*t16*t94; d25 = T(1)/d25;
complex<T> d27 = spb14*t2*square(spb24); d27 = T(1)/d27;
complex<T> d28 = spb14*spb23*t1; d28 = T(1)/d28;
complex<T> d30 = spb12*spb14*spb23*t16; d30 = T(1)/d30;
complex<T> d31 = spb34*t94*T(2); d31 = T(1)/d31;
complex<T> d32 = spb12*t57*T(2); d32 = T(1)/d32;
complex<T> d34 = t2*square(spb24)*T(2); d34 = T(1)/d34;
complex<T> t54 = -t74; 
complex<T> t59 = d2*spa23; 
complex<T> t61 = d14*spb36; 
complex<T> t69 = d5*spb12; 
complex<T> t76 = d6*spa45; 
complex<T> t83 = spa56*(d29*spb14*t100 + d30*t33*t73 + d28*t93); 
complex<T> t85 = -t93; 
complex<T> t86 = t18*t74; 
complex<T> t90 = -(d8*t43); 
complex<T> t120 = spa23*t93; 
complex<T> d1 = t20*t28*(s14 + t28); d1 = T(1)/d1;
complex<T> d7 = (-s123 + s56)*t16*t38; d7 = T(1)/d7;
complex<T> d10 = (-s124 + s14)*t20; d10 = T(1)/d10;
complex<T> d11 = t20*square(s124 - s14)*T(2); d11 = T(1)/d11;
complex<T> d12 = (-s124 + s14)*t20*(s14 + t28); d12 = T(1)/d12;
complex<T> d15 = t16*t29*t38; d15 = T(1)/d15;
complex<T> d18 = t38*cube(spb34); d18 = T(1)/d18;
complex<T> d19 = spb14*t16*t38; d19 = T(1)/d19;
complex<T> d20 = t20*square(s14 + t28); d20 = T(1)/d20;
complex<T> d26 = spb14*t38; d26 = T(1)/d26;
complex<T> d33 = t20*square(s14 + t28)*T(2); d33 = T(1)/d33;
complex<T> d35 = t38*cube(spb34)*T(2); d35 = T(1)/d35;
complex<T> d36 = spb14*t16*t38*T(2); d36 = T(1)/d36;
complex<T> t9 = (-s156 + s34)*(d23*t17 + d17*t85); 
complex<T> t10 = d21*spb16*t25 + d20*s14*t86; 
complex<T> t11 = s12*(d34*spb16*t25 + d33*s14*t86); 
complex<T> t13 = -(d24*spa23*spb14*t100) + d22*t120 + d23*s23*t17 + d25*spa23*t33*t53; 
complex<T> t26 = -(d19*t33); 
complex<T> t27 = t45*t69*square(spa24) - d1*spb14*square(spa24)*square(spb26) + d3*t47*square(spb36) - spb13*t59*square(spb36) - d4*spb46*t101*T(2); 
complex<T> t65 = d7*t44; 
complex<T> t77 = d11*t18; 
complex<T> t78 = d10*T(2); 
complex<T> t92 = t41*(d35*spb14*t100 + d36*t33*t73); 
complex<T> t12 = (s123 - s23)*(d18*spb14*t100 - t26*t73); 
complex<T> t14 = -(d9*spa35*spb16*t30) - d15*t100*t31 - d3*t19*t47 - d16*spb56*t42*t67 - t15*t45*t69 + t74*t77 + d13*t85 + d1*t86 + d12*t86 - d10*spb26*t101*T(2) + d4*spb46*t101*T(2) + t31*t61*t73*T(2); 
complex<T> t22 = -t65; 
complex<T> t72 = d18*t100*t23 + t28*(d27*spb16*t46 + t26*t73 + d20*t86); 
complex<T> t84 = d12*t18*t54 + t54*t77 + spb26*t101*t78; 
complex<T> t99 = t100*t59 + t53*t76 + spb13*t65*t79 + spb13*t90*t94; 
complex<T> t66 = d9*spa35*spb16*t30 + d15*t100*t31 + d16*spb56*t42*t67 + t24*t61*t73 + t73*t76 + spb13*t22*t79 + d8*spb13*t43*t94 + d17*t93*T(2); 
complex<T> co1 = d26*spa34*t85; 
complex<T> co2 = d31*spa12*t120; 
complex<T> co3 = d13*t49*t93; 
complex<T> co4 = -(d32*s23*spa34*t17); 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t27*(*CI_users[0]->get_value(mc,ind,mu)) + t99*(*CI_users[1]->get_value(mc,ind,mu)) + t14*(*CI_users[2]->get_value(mc,ind,mu)) + t84*(*CI_users[3]->get_value(mc,ind,mu)) + t66*(*CI_users[4]->get_value(mc,ind,mu)) + t12*(*CI_users[6]->get_value(mc,ind,mu)) + t10*(*CI_users[7]->get_value(mc,ind,mu)) + t13*(*CI_users[8]->get_value(mc,ind,mu)) + t9*(*CI_users[9]->get_value(mc,ind,mu)) + co1*(*CI_users[10]->get_value(mc,ind,mu)) + t72*(*CI_users[11]->get_value(mc,ind,mu)) + t83*(*CI_users[12]->get_value(mc,ind,mu)) + co2*(*CI_users[13]->get_value(mc,ind,mu)) + co3*(*CI_users[14]->get_value(mc,ind,mu)) + co4*(*CI_users[15]->get_value(mc,ind,mu)) + t11*(*CI_users[16]->get_value(mc,ind,mu)) + t92*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmppemep_SLC_wCI::\
C2q2g2l_qpqmppemep_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c256));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c4, c356));
CI_users.push_back(new Cached_Box_Integral_User(c2, c14, c3, c56));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c4, c23, c1, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmppemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmppemep SLC");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 19:26:55 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa56 = SPA(5,6);
complex<T> spa25 = SPA(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb13 = SPB(1,3);
complex<T> spb12 = SPB(1,2);
complex<T> spb46 = SPB(4,6);
complex<T> spb24 = SPB(2,4);
complex<T> spb36 = SPB(3,6);
complex<T> spb16 = SPB(1,6);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s124 = SS(1,2,4);
complex<T> s14 = -(spa14*spb14);
complex<T> s134 = SS(1,3,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s156 = SS(1,5,6);
complex<T> s234 = SS(2,3,4);
complex<T> s256 = SS(2,5,6);
complex<T> t4 = spa14*spa34; 
complex<T> t5 = spa56*T(2); 
complex<T> t24 = square(spa14); 
complex<T> t26 = square(spa15); 
complex<T> t27 = square(spa24); 
complex<T> t28 = square(spa45); 
complex<T> t29 = spa23*spa56; 
complex<T> t30 = square(spa12); 
complex<T> t41 = square(spa25); 
complex<T> t42 = cube(spa14); 
complex<T> t43 = spa24*spa35 + spa23*spa45; 
complex<T> t44 = -(spa23*spb13) - spa24*spb14; 
complex<T> t48 = -(spa15*spa24) + spa12*spa45; 
complex<T> t55 = spa15*spa24 + spa14*spa25; 
complex<T> t56 = -(spa15*spa24) + spa14*spa25*T(2); 
complex<T> t57 = spa15*spa23 - spa12*spa35; 
complex<T> t59 = square(spb16); 
complex<T> t60 = -(spa35*(spa15*spb13 + spa25*spb23)); 
complex<T> t61 = -(spa15*spb14) - spa25*spb24; 
complex<T> t63 = -(spa12*spb14) + spa23*spb34; 
complex<T> t64 = square(spb36); 
complex<T> t65 = square(spb46); 
complex<T> t67 = -s256 + s34; 
complex<T> t68 = s124*s134 - s14*s56; 
complex<T> t90 = spa15*spa25; 
complex<T> t101 = spa25*spa45; 
complex<T> d2 = (s12 - s123)*spa13*spa56*square(spa34); d2 = T(1)/d2;
complex<T> d3 = (s12 - s124)*spa14*spa56*square(spa34); d3 = T(1)/d3;
complex<T> d8 = spa14*square(spa34); d8 = T(1)/d8;
complex<T> d9 = -s123 + s56; d9 = T(1)/d9;
complex<T> d10 = (-s124 + s56)*spa14*square(spa34); d10 = T(1)/d10;
complex<T> d21 = -s234 + s56; d21 = T(1)/d21;
complex<T> d25 = spa14*spa56*cube(spa34); d25 = T(1)/d25;
complex<T> d27 = spa14*spa56*square(spa34); d27 = T(1)/d27;
complex<T> d32 = spa14*cube(spa34); d32 = T(1)/d32;
complex<T> d33 = spa14*spa23*square(spa34); d33 = T(1)/d33;
complex<T> t17 = d25*s23*spa23*t28 - d27*spa25*spb23*t43; 
complex<T> t31 = -(spa25*t48); 
complex<T> t76 = -(spa23*t28); 
complex<T> t100 = -(t26*t27); 
complex<T> d1 = spa13*spa56*t4; d1 = T(1)/d1;
complex<T> d4 = t29*t4; d4 = T(1)/d4;
complex<T> d5 = (-s123 + s23)*spa13*t24*t29; d5 = T(1)/d5;
complex<T> d6 = spa23*t4*square(s123 - s56); d6 = T(1)/d6;
complex<T> d7 = spa23*spa34*t24; d7 = T(1)/d7;
complex<T> d11 = t4*square(s124 - s56); d11 = T(1)/d11;
complex<T> d12 = spa13*t24*t29; d12 = T(1)/d12;
complex<T> d13 = (s23 - s234)*spa56*t4; d13 = T(1)/d13;
complex<T> d14 = (s23 - s234)*spa34*spa56*t24; d14 = T(1)/d14;
complex<T> d15 = t4*t5*square(s23 - s234); d15 = T(1)/d15;
complex<T> d16 = spa23*t4*t5; d16 = T(1)/d16;
complex<T> d17 = (-s234 + s56)*t29*t4; d17 = T(1)/d17;
complex<T> d18 = spa23*t4*square(s234 - s56)*T(2); d18 = T(1)/d18;
complex<T> d19 = spa34*t24*t29; d19 = T(1)/d19;
complex<T> d20 = spa23*t4; d20 = T(1)/d20;
complex<T> d22 = spa56*t4*square(spa13); d22 = T(1)/d22;
complex<T> d23 = spa13*t29*t4; d23 = T(1)/d23;
complex<T> d24 = spa34*t29*t42; d24 = T(1)/d24;
complex<T> d26 = spa14*t29*square(spa34); d26 = T(1)/d26;
complex<T> d28 = t29*t42; d28 = T(1)/d28;
complex<T> d29 = spa14*t29; d29 = T(1)/d29;
complex<T> d30 = t24*t29; d30 = T(1)/d30;
complex<T> d31 = spa23*spa34*t42; d31 = T(1)/d31;
complex<T> d34 = spa23*spa34*t5; d34 = T(1)/d34;
complex<T> d35 = spa14*t5*cube(spa34); d35 = T(1)/d35;
complex<T> d36 = spa14*spa23*t5*square(spa34); d36 = T(1)/d36;
complex<T> d37 = t4*t5*square(spa13); d37 = T(1)/d37;
complex<T> d38 = spa13*t4*t5; d38 = T(1)/d38;
complex<T> d39 = spa23*spa34*t42*t5; d39 = T(1)/d39;
complex<T> d40 = spa23*spa34*t24*t5; d40 = T(1)/d40;
complex<T> t11 = (-s156 + s34)*(d24*t100 + d19*t31); 
complex<T> t12 = (s123*s234 - s23*s56)*(d39*t26*t27 + d40*spa25*t48); 
complex<T> t15 = (s123*s124 - s12*s56)*(d35*spa23*t28 + d36*spa25*t43); 
complex<T> t18 = d2*spa23*t60 - d3*spa24*spa45*t61 + d1*t90; 
complex<T> t19 = spb34*(d28*t100 + d30*t31 - d29*t41); 
complex<T> t20 = d14*spa45*spb34*t55 + d12*spa12*t90 - d5*t30*(spa15*spa35*spb13 + spb12*t90) + d15*spa23*t28*square(spb34) - d13*spb34*t101*T(2); 
complex<T> t21 = d10*spa24*spa35*spb36 + d3*spa24*spa45*t61 - d11*t29*t64; 
complex<T> t23 = spb56*(d31*t26*t27 + d32*spa23*t28 + d20*t41 + d33*spa25*t43 + d7*spa25*t48); 
complex<T> t33 = -(d7*t56); 
complex<T> t35 = -(d26*t43); 
complex<T> t39 = -(d4*t63); 
complex<T> t74 = s12*(d37*s23*spa23*t26 + d38*spa25*spb23*t57); 
complex<T> t110 = d22*spa23; 
complex<T> t115 = -(d4*t41); 
complex<T> t138 = d16*t41; 
complex<T> t1 = -(d20*d21*spa12*spa25*spb16) - t138 - d19*d21*spa24*t26*t44 - d14*spa45*spb34*t55 - d18*spa56*t30*t59 + d15*t76*square(spb34) + d13*spb34*t101*T(2) + d17*t44*t90*T(2); 
complex<T> t2 = -(d8*d9*spa24*spa45*spb46) + t115 + d5*spa15*spa35*spb13*t30 + d7*d9*spa24*spb46*t56 - d2*spa23*t60 + d4*d9*t101*t63 - d6*spa56*t27*t65 + d5*spb12*t30*t90; 
complex<T> t3 = d20*d21*spa12*spa25*spb16 - d10*spa24*spa35*spb36 + d8*d9*spa24*spa45*spb46 + d9*spa24*spb46*t33 + d9*t101*t39 + d19*d21*spa24*t26*t44 + d18*spa56*t30*t59 + d11*t29*t64 + d6*spa56*t27*t65 - d1*t90 - d12*spa12*t90 - d17*t44*t90*T(2) + d4*t41*T(3); 
complex<T> t10 = s12*(t110*t26 - d23*spa25*t57); 
complex<T> t13 = (s12 - s124)*(spa25*t35 + d25*t76); 
complex<T> t14 = (-s123 + s23)*(d24*t100 + t110*t26 + d19*t31 + spa25*t35 - d23*spa25*t57 + d25*t76); 
complex<T> co1 = t115*t67; 
complex<T> co2 = -(d34*s12*spb14*t41); 
complex<T> co3 = t138*t68; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t18*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t21*(*CI_users[2]->get_value(mc,ind,mu)) + t20*(*CI_users[3]->get_value(mc,ind,mu)) + t1*(*CI_users[4]->get_value(mc,ind,mu)) + t3*(*CI_users[5]->get_value(mc,ind,mu)) + t10*(*CI_users[6]->get_value(mc,ind,mu)) + t14*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t17*(*CI_users[9]->get_value(mc,ind,mu)) + t11*(*CI_users[10]->get_value(mc,ind,mu)) + t19*(*CI_users[11]->get_value(mc,ind,mu)) + t13*(*CI_users[12]->get_value(mc,ind,mu)) + t23*(*CI_users[13]->get_value(mc,ind,mu)) + co2*(*CI_users[14]->get_value(mc,ind,mu)) + co3*(*CI_users[15]->get_value(mc,ind,mu)) + t15*(*CI_users[16]->get_value(mc,ind,mu)) + t74*(*CI_users[17]->get_value(mc,ind,mu)) + t12*(*CI_users[18]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmpmemep_SLC_wCI::\
C2q2g2l_qpqmpmemep_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c134, c256));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c2356));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c256));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c2356));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c14, c23, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c23, c14, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c14, c56));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c256));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c356));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c23, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmpmemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmpmemep SLC");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 22:37:26 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb14 = SPB(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb12 = SPB(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb34 = SPB(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb36 = SPB(3,6);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb46 = SPB(4,6);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spa26 = SPA(2,6);
complex<T> s12 = -(spa12*spb12);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s34 = -(spa34*spb34);
complex<T> s124 = SS(1,2,4);
complex<T> s134 = SS(1,3,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s24 = -(spa24*spb24);
complex<T> s234 = SS(2,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s256 = SS(2,5,6);
complex<T> t37 = spa14*spb16; 
complex<T> t38 = -(spa25*spb23); 
complex<T> t41 = -(spa35*spb12); 
complex<T> t42 = -(spa12*spb46); 
complex<T> t44 = spa56*spb26; 
complex<T> t45 = spa15*spb56; 
complex<T> t67 = spa14*spb12; 
complex<T> t68 = spa12*spb23; 
complex<T> t76 = s123*s34; 
complex<T> t124 = square(s123); 
complex<T> t125 = square(spb46); 
complex<T> t126 = square(s124); 
complex<T> t127 = square(spa35); 
complex<T> t130 = square(spb26); 
complex<T> t131 = square(spa15); 
complex<T> t135 = spa45*spb14; 
complex<T> t138 = spa23*spb36; 
complex<T> t144 = spb12*spb34; 
complex<T> t169 = square(spb23); 
complex<T> t170 = square(spa24); 
complex<T> t171 = square(spa14); 
complex<T> t172 = square(spb13); 
complex<T> t173 = square(spa25*spb24 + spa35*spb34); 
complex<T> t176 = (spa13*spb12 + spa34*spb24)*spb56; 
complex<T> t177 = square(spa35*spb23 + spa45*spb24); 
complex<T> t178 = square(spa15*spb14 + spa25*spb24); 
complex<T> t179 = square(spa13*spb16 + spa23*spb26); 
complex<T> t208 = (-s134 + s56)*spa34; 
complex<T> t210 = (-s234 + s56)*spb34; 
complex<T> t231 = cube(spa12*spb23 - spa14*spb34); 
complex<T> t293 = spa13*spb14; 
complex<T> t294 = spa23*spb24; 
complex<T> t295 = s34*T(2); 
complex<T> t296 = spa45*spb34; 
complex<T> t298 = spa15*spb12; 
complex<T> t299 = -(spa24*spb26); 
complex<T> t300 = cube(spa13*spb23 + spa14*spb24); 
complex<T> t303 = spa34*spb23; 
complex<T> t304 = spa56*(spa13*spb12 + spa34*spb24); 
complex<T> t306 = spa14*spb46; 
complex<T> t311 = square(s234); 
complex<T> t318 = -((spa24*spb12 + spa34*spb13)*(spa25*spb26 + spa35*spb36)); 
complex<T> t338 = s14 - s23 - s56; 
complex<T> t339 = -s14 + s23 - s56; 
complex<T> t340 = s12 - s34 - s56; 
complex<T> t341 = -s14 - s23 + s56; 
complex<T> t342 = spa23*spb12 - spa34*spb14; 
complex<T> t344 = square(s56); 
complex<T> t345 = -s12 - s34 + s56; 
complex<T> t346 = spa15*spb13; 
complex<T> t347 = -(spa23*spb13) - spa24*spb14; 
complex<T> t348 = spa12*spb26 + spa13*spb36; 
complex<T> t350 = -s12 + s34 - s56; 
complex<T> t352 = square(s12); 
complex<T> t353 = square(s34); 
complex<T> t354 = -(spb23*T(2)); 
complex<T> t355 = -(spa12*spb26); 
complex<T> t357 = -s123 + s23; 
complex<T> t362 = spa25*spb12 + spa35*spb13; 
complex<T> t371 = square(spa13); 
complex<T> t372 = square(spb24); 
complex<T> t373 = -(spa13*spb12) - spa34*spb24; 
complex<T> t374 = spb12*spb36; 
complex<T> t384 = -s134 + s14 + s34; 
complex<T> t385 = s23 - s234 + s34; 
complex<T> t399 = square(spb16); 
complex<T> t401 = -square(s123); 
complex<T> t402 = -square(s124); 
complex<T> t413 = -s124 + s14; 
complex<T> t418 = s13 - s24; 
complex<T> t480 = -(s56*T(2)); 
complex<T> t484 = s12*s34; 
complex<T> t491 = s14*s23; 
complex<T> t495 = s13 + s23; 
complex<T> t496 = s14 + s24; 
complex<T> t543 = square(spa13*spb12 + spa34*spb24); 
complex<T> t549 = spa14*T(2); 
complex<T> t550 = spa24*spb23; 
complex<T> t551 = cube(spa13*spb12 + spa34*spb24); 
complex<T> t553 = spa25*spb13; 
complex<T> t554 = spa34*spb36; 
complex<T> t555 = spa15*spb16; 
complex<T> t556 = spa56*spb34; 
complex<T> t558 = spa13*(spa13*spb23 + spa14*spb24); 
complex<T> t560 = s134*spb14; 
complex<T> t568 = cube(spa24); 
complex<T> t569 = -square(spa15); 
complex<T> t621 = spa12*spb24; 
complex<T> t622 = spa13*spb34; 
complex<T> t624 = spa14*spb23; 
complex<T> t631 = spa34*spb56; 
complex<T> t640 = s23*s56; 
complex<T> t661 = spa56*T(2); 
complex<T> t662 = spb24*(spa13*spb23 + spa14*spb24); 
complex<T> t691 = -(spa45*spb56); 
complex<T> t693 = spa12*spb13; 
complex<T> t695 = spa25*spb16; 
complex<T> t731 = spa14*spb13; 
complex<T> t736 = spa35*spb56; 
complex<T> t737 = s234*spa23; 
complex<T> t761 = -(spa25*spb56); 
complex<T> t768 = spa34*spb46; 
complex<T> t789 = spa12*spb14; 
complex<T> t793 = spa56*spb16; 
complex<T> t820 = s56*T(2); 
complex<T> t821 = spa24*spb36; 
complex<T> t847 = spa13*spb24; 
complex<T> t855 = spa14*spb26; 
complex<T> t895 = spa14*spa34; 
complex<T> t897 = spb23*spb34; 
complex<T> t930 = spb13*spb36; 
complex<T> d3 = spa13*spb23 + spa14*spb24; d3 = T(1)/d3;
complex<T> d7 = spa56; d7 = T(1)/d7;
complex<T> d8 = spb56; d8 = T(1)/d8;
complex<T> d16 = spb34*spb56; d16 = T(1)/d16;
complex<T> d19 = spa34*spa56; d19 = T(1)/d19;
complex<T> d41 = spa23*spa56; d41 = T(1)/d41;
complex<T> d48 = spa23*spb56; d48 = T(1)/d48;
complex<T> d52 = spa13*spb12 + spa34*spb24; d52 = T(1)/d52;
complex<T> d57 = spb14*spb56; d57 = T(1)/d57;
complex<T> d58 = spa56*spb14; d58 = T(1)/d58;
complex<T> d69 = spa56*(spa13*spb23 + spa14*spb24); d69 = T(1)/d69;
complex<T> d72 = (spa13*spb23 + spa14*spb24)*spb56; d72 = T(1)/d72;
complex<T> d75 = s34*s56*spa23*spb14; d75 = T(1)/d75;
complex<T> d78 = T(2); d78 = T(1)/d78;
complex<T> d103 = s234; d103 = T(1)/d103;
complex<T> d107 = s124*spb14*(spa13*spb12 + spa34*spb24); d107 = T(1)/d107;
complex<T> d110 = spb14*square(s134); d110 = T(1)/d110;
complex<T> d114 = spa23; d114 = T(1)/d114;
complex<T> d116 = spb14; d116 = T(1)/d116;
complex<T> d118 = s134; d118 = T(1)/d118;
complex<T> d119 = spa34*spb14; d119 = T(1)/d119;
complex<T> d121 = spb14*spb34; d121 = T(1)/d121;
complex<T> d124 = (spa13*spb23 + spa14*spb24)*(spa13*spb12 + spa34*spb24)*T(2); d124 = T(1)/d124;
complex<T> d134 = spa23*spa34; d134 = T(1)/d134;
complex<T> d135 = spa23*spb34; d135 = T(1)/d135;
complex<T> d141 = s123; d141 = T(1)/d141;
complex<T> d143 = s124; d143 = T(1)/d143;
complex<T> t39 = -t296; 
complex<T> t40 = -t554; 
complex<T> t70 = d41*(-(s234*(spa15*spb14 + spa25*spb24)) + s124*(spa25*spb24 + spa35*spb34)) + spb24*spb36; 
complex<T> t71 = spa13*(spa45 + d57*(-s123 + s134)*spb16) + d57*(s134*spa23*spb26 + s123*t768); 
complex<T> t120 = -t550; 
complex<T> t121 = -t731; 
complex<T> t123 = square(t348); 
complex<T> t128 = square(t342); 
complex<T> t141 = spa12*t340; 
complex<T> t153 = spa34*(spb14*t345 - spa23*spb12*spb34*T(2)); 
complex<T> t154 = d41*(spa25*t338 + spa23*spb36*t661); 
complex<T> t162 = d57*(spb16*t339 - spa45*spb14*spb56*T(2)); 
complex<T> t168 = (t293 + t294)*T(2); 
complex<T> t174 = square(-(spa13*spb16) + t768); 
complex<T> t175 = spa56*(t621 + t622); 
complex<T> t180 = square(spa13*spb36 + t306); 
complex<T> t181 = spa35*spb13*t350 + t345*t793; 
complex<T> t182 = spa24*spb46*t350 + t345*t761; 
complex<T> t191 = s12*t295 + s56*t295 - t344 - t352 - t353 + t350*(-s14 + s23 + t418) + s12*t820; 
complex<T> t192 = s12*t295 + s56*t295 - t344 - t352 - t353 - t350*(s13 + s23 - t496) + s12*t820; 
complex<T> t195 = spb12*t341 - spb14*t303*T(2); 
complex<T> t196 = spa25*spb15*spb23 + spa26*spb16*spb23 + spb13*t339; 
complex<T> t203 = spa14*spa25*spb15 + spa24*t338 + spa26*t37; 
complex<T> t276 = t298*t340 - t345*t44; 
complex<T> t301 = -t693; 
complex<T> t305 = spb56*(t621 + t622); 
complex<T> t307 = s134*(spa34*spb23 - t67); 
complex<T> t308 = d7*spa25; 
complex<T> t312 = -(d78*s56); 
complex<T> t313 = t126*t127; 
complex<T> t319 = (spa45*spb46 + t555)*(spa24*spb34 + t693); 
complex<T> t325 = d3*(-(d7*spa15*spa25*spb12) + d8*spb16*t355); 
complex<T> t343 = spa23*spb34 - t789; 
complex<T> t351 = spa45*spb24 - t298; 
complex<T> t356 = cube(t293 + t294); 
complex<T> t359 = -(spa14*spb34) + t68; 
complex<T> t363 = -(spa25*spb26) - t555; 
complex<T> t365 = s14 - s23 + t418; 
complex<T> t379 = t299 - t37; 
complex<T> t380 = -(spa15*spb13) + t38; 
complex<T> t381 = t303 - t67; 
complex<T> t392 = -(spa23*t125); 
complex<T> t430 = -t550 - t731; 
complex<T> t458 = -(spa12*spa45*spb16*spb23) + t298*t821; 
complex<T> t460 = spa15*spa34*spb13 + spa24*t298 + spa56*t299; 
complex<T> t477 = spa25*spb12*t340 + t345*t793; 
complex<T> t485 = t296 + t38; 
complex<T> t486 = d8*spb36; 
complex<T> t494 = s134*t176; 
complex<T> t531 = d3*t347; 
complex<T> t544 = square(t621 + t622); 
complex<T> t552 = cube(t621 + t622); 
complex<T> t623 = -(t484*T(2)); 
complex<T> t625 = d7*spa45; 
complex<T> t690 = s12*t480; 
complex<T> t697 = d78*s34; 
complex<T> t701 = spa56*t341; 
complex<T> t720 = spb13*t171; 
complex<T> t729 = s34*t480; 
complex<T> t765 = d78*s12; 
complex<T> t785 = t169*t170; 
complex<T> t792 = d3*s12; 
complex<T> t810 = t171*t172; 
complex<T> t989 = t37*t553; 
complex<T> d1 = (s12 - s124)*spa56*spb24*square(t293 + t294); d1 = T(1)/d1;
complex<T> d10 = t295*t662; d10 = T(1)/d10;
complex<T> d11 = t295*t558; d11 = T(1)/d11;
complex<T> d12 = (s12 - s123)*spa13*spb56*square(t293 + t294); d12 = T(1)/d12;
complex<T> d14 = spb34*spb56*(t293 + t294); d14 = T(1)/d14;
complex<T> d15 = t293 + t294; d15 = T(1)/d15;
complex<T> d18 = spa34*spa56*(t293 + t294); d18 = T(1)/d18;
complex<T> d21 = spa23*(t293 + t294)*(t621 + t622)*square(s123 - s56); d21 = T(1)/d21;
complex<T> d26 = (-s124 + s56)*t304*square(t293 + t294); d26 = T(1)/d26;
complex<T> d27 = (-s124 + s56)*spa56*spb14*(t293 + t294)*t543; d27 = T(1)/d27;
complex<T> d28 = spa56*(t293 + t294)*t413*t543; d28 = T(1)/d28;
complex<T> d29 = spb24*(t293 + t294)*t304*t413; d29 = T(1)/d29;
complex<T> d30 = spb14*(spa13*spb12 + spa34*spb24)*(t293 + t294)*square(s124 - s56); d30 = T(1)/d30;
complex<T> d31 = (-s134 + s34)*t176*square(spa13*spb23 + spa14*spb24); d31 = T(1)/d31;
complex<T> d32 = t176*t208*square(spa13*spb23 + spa14*spb24); d32 = T(1)/d32;
complex<T> d33 = spa34*(spa13*spb23 + spa14*spb24)*t176*square(s134 - s56)*T(2); d33 = T(1)/d33;
complex<T> d34 = (spa13*spb23 + spa14*spb24)*t176*t208; d34 = T(1)/d34;
complex<T> d35 = s134*(-s134 + s34)*spb56*t558; d35 = T(1)/d35;
complex<T> d36 = s134*spb56*t558*square(s134 - s34)*T(2); d36 = T(1)/d36;
complex<T> d44 = (t293 + t294)*square(t344 + s14*t480 + s23*t480 + square(s14) + square(s23) - t491*T(2)); d44 = T(1)/d44;
complex<T> d45 = t344 + s14*t480 + s23*t480 + square(s14) + square(s23) - t491*T(2); d45 = T(1)/d45;
complex<T> d46 = t621 + t622; d46 = T(1)/d46;
complex<T> d49 = (spa13*spb12 + spa34*spb24)*(t621 + t622)*(t344 + s14*t480 + s23*t480 + square(s14) + square(s23) - t491*T(2)); d49 = T(1)/d49;
complex<T> d50 = spa23*t304; d50 = T(1)/d50;
complex<T> d51 = square(t344 + s14*t480 + s23*t480 + square(s14) + square(s23) - t491*T(2)); d51 = T(1)/d51;
complex<T> d62 = s234*t661*t662*square(s234 - s34); d62 = T(1)/d62;
complex<T> d67 = s234*(-s234 + s34)*spa56*t662; d67 = T(1)/d67;
complex<T> d76 = s34*s56*t737; d76 = T(1)/d76;
complex<T> d77 = s34*s56*t560; d77 = T(1)/d77;
complex<T> d80 = spa56*spb14*(t293 + t294)*t551; d80 = T(1)/d80;
complex<T> d86 = (spa23*spb13 + spa24*spb14)*t631*t737; d86 = T(1)/d86;
complex<T> d91 = (spa23*spb13 + spa24*spb14)*t556*t560; d91 = T(1)/d91;
complex<T> d92 = spa34*t176*t300; d92 = T(1)/d92;
complex<T> d94 = spa34*t176*t300*T(2); d94 = T(1)/d94;
complex<T> d100 = (t293 + t294)*t304*t372; d100 = T(1)/d100;
complex<T> d101 = spb14*(t293 + t294)*t304*t372; d101 = T(1)/d101;
complex<T> d104 = spa23*t311; d104 = T(1)/d104;
complex<T> d105 = t556*t737*T(2); d105 = T(1)/d105;
complex<T> d106 = t556*t560*T(2); d106 = T(1)/d106;
complex<T> d108 = (spa13*spb12 + spa34*spb24)*(t293 + t294); d108 = T(1)/d108;
complex<T> d109 = s123*spa23*(t621 + t622); d109 = T(1)/d109;
complex<T> d112 = (t293 + t294)*(t621 + t622); d112 = T(1)/d112;
complex<T> d120 = spb14*t556; d120 = T(1)/d120;
complex<T> d122 = spb14*t631; d122 = T(1)/d122;
complex<T> d123 = t631*square(s134); d123 = T(1)/d123;
complex<T> d125 = t631*t737*T(2); d125 = T(1)/d125;
complex<T> d126 = t560*t631*T(2); d126 = T(1)/d126;
complex<T> d132 = spa23*t556; d132 = T(1)/d132;
complex<T> d133 = t311*t556; d133 = T(1)/d133;
complex<T> d136 = spa23*t631; d136 = T(1)/d136;
complex<T> d137 = (spa13*spb23 + spa14*spb24)*(t621 + t622)*T(2); d137 = T(1)/d137;
complex<T> d139 = s123*(t293 + t294); d139 = T(1)/d139;
complex<T> d140 = (spa13*spb12 + spa34*spb24)*(t621 + t622); d140 = T(1)/d140;
complex<T> d142 = (t293 + t294)*(t344 + s14*t480 + s23*t480 + square(s14) + square(s23) - t491*T(2)); d142 = T(1)/d142;
complex<T> d148 = s124*(t293 + t294); d148 = T(1)/d148;
complex<T> d150 = t176*t300*T(2); d150 = T(1)/d150;
complex<T> d151 = spa34*(spa13*spb12 + spa34*spb24)*t300*T(2); d151 = T(1)/d151;
complex<T> d156 = spb34*t300*(t621 + t622)*T(2); d156 = T(1)/d156;
complex<T> t20 = d49*(d7*spb13*(s124*spa15 - spa56*t306)*(spa34*(t296 - t346) + spa24*t351 - spa14*t362) - d48*spa24*spb16*((s12 - s34)*(t306*t338 - t341*t45) + (t621 + t622)*(-(t338*t37) + t341*t691)) + (s12 - s34)*(spa45*spb46 + t555)*t731*T(2)); 
complex<T> t25 = -((s12 - s34)*spa24*(spa25*spb26 + spa35*spb36)*t354) + d58*t553*((s12 - s34)*(spa35*spb23*t339 + t341*t44) + t373*(t339*t38 + spb36*t701)) + d8*spa24*(-(spa12*spb16*spb23) + spa24*spb23*spb46 + spb34*t299 - spb13*t348 + spb34*t40)*(s123*spb26 + spb23*t736); 
complex<T> t26 = d46*spa12*(d7*spa15*t135 + d8*spb46*t37)*(-(s12*t338) + s34*t338 + t344 - t338*t418 - s14*t820 + square(s14) + square(s23) - t491*T(2) - t640*T(2)) + d8*t37*(t299*t343 + t343*t40 + spa24*spb46*t338*T(2) + t341*t761*T(2)) - t625*(s14*spa12*t362 + spa14*spa23*spb34*t362 + spa15*spa23*spb13*t339*T(2) + spa23*spb13*t306*t661*T(2) + spa12*t135*t338*T(3)) + d45*spa14*spa23*t318*(spb34*t341 - spb14*t68*T(2))*T(6); 
complex<T> t27 = -(d52*spb12*(d8*spb26*t138 - d7*spa35*t38)*(-(s12*t339) + s13*t339 - s24*t339 + s34*t339 + t344 - s14*t820 + square(s14) + square(s23) - t491*T(2) - t640*T(2))) + d7*t38*(t296*t342 - t342*t346 + spa35*spb13*t339*T(2) + t341*t793*T(2)) + t486*(s23*spb12*(spa12*spb16 - spa24*spb46) + spa12*spb14*spb16*t303 - spa24*spb14*spb46*t303 + spa24*spb14*spb26*t338*T(2) + spa24*spb14*t354*t736*T(2) - spb12*t138*t339*T(3)) - d45*spb14*spb23*t319*(spa34*t341 - spa23*t67*T(2))*T(6); 
complex<T> t43 = -(t430*T(3)); 
complex<T> t72 = s34*spb14*t379 + spa23*t144*t379 + spa24*spb12*T(2)*(spb46*t340 + spb34*t736*T(2)) - spb14*t350*t554*T(3); 
complex<T> t73 = t343*t379 + t182*T(2); 
complex<T> t74 = s34*spa23*t380 + spa34*t380*t789 - t693*(spa35*t340 + t661*t768)*T(2) + spa23*t296*t350*T(3); 
complex<T> t75 = t342*t380 + t181*T(2); 
complex<T> t122 = square(t351); 
complex<T> t129 = square(t343); 
complex<T> t132 = -t307; 
complex<T> t134 = s234*t359; 
complex<T> t137 = d80*(t402*square(spa35)*square(spb12) + square(t342)*square(t351)); 
complex<T> t143 = t363*t430; 
complex<T> t149 = -t531; 
complex<T> t197 = -square(t351); 
complex<T> t200 = -(d41*spa25*t362) + spb16*t486; 
complex<T> t201 = d57*spb16*(-(spa12*spb16) + spa24*spb46) - spa45*t308; 
complex<T> t239 = d69*spa15*(spa24*t298 + spa56*t299 + spa34*t346)*t354 + t299*t486 + t346*t625; 
complex<T> t279 = -(spb26*t141) + t345*t45; 
complex<T> t324 = d67*(spa35*spb23 + spa45*spb24); 
complex<T> t328 = d35*(spa13*spb36 + spa14*spb46); 
complex<T> t349 = -t37 + t40; 
complex<T> t390 = t120 + t731; 
complex<T> t391 = t120 + t121 - d3*s56*t624*T(4); 
complex<T> t394 = spa24*t485; 
complex<T> t400 = s234*t175; 
complex<T> t410 = d7*spa35*t296 + t486*t768; 
complex<T> t412 = d8*spa12*spb16*spb26 - t298*t308; 
complex<T> t478 = -(spb16*t141) + t345*t761; 
complex<T> t502 = -(d44*T(3)); 
complex<T> t503 = -(spa23*t123); 
complex<T> t510 = d34*spb26; 
complex<T> t527 = d51*t318; 
complex<T> t576 = d92*t171; 
complex<T> t582 = d31*s134; 
complex<T> t606 = -(t363*T(2)); 
complex<T> t612 = d86*t399; 
complex<T> t629 = -(t172*t174); 
complex<T> t632 = t121 + t550; 
complex<T> t633 = s123*t343; 
complex<T> t647 = d51*t319; 
complex<T> t652 = d27*t351; 
complex<T> t653 = d69*t460; 
complex<T> t671 = t307*square(spb26); 
complex<T> t712 = t162*t347; 
complex<T> t716 = (spa13*spb23 + spa14*spb24)*t168; 
complex<T> t749 = d30*t342; 
complex<T> t788 = spa13*t494; 
complex<T> t875 = t312*t341; 
complex<T> t880 = spb26*t301; 
complex<T> d4 = spa56*(t344 + t352 + t353 + t623 + t690 + t729); d4 = T(1)/d4;
complex<T> d5 = square(t344 + t352 + t353 + t623 + t690 + t729); d5 = T(1)/d5;
complex<T> d6 = t344 + t352 + t353 + t623 + t690 + t729; d6 = T(1)/d6;
complex<T> d9 = spb56*(t344 + t352 + t353 + t623 + t690 + t729); d9 = T(1)/d9;
complex<T> d17 = spa13*t168; d17 = T(1)/d17;
complex<T> d20 = spb24*t168; d20 = T(1)/d20;
complex<T> d22 = spa13*(t293 + t294)*t305*t357; d22 = T(1)/d22;
complex<T> d23 = spb56*(t293 + t294)*t357*t544; d23 = T(1)/d23;
complex<T> d24 = (-s123 + s56)*spa23*spb56*(t293 + t294)*t544; d24 = T(1)/d24;
complex<T> d25 = (-s123 + s56)*t305*square(t293 + t294); d25 = T(1)/d25;
complex<T> d42 = (spa13*spb12 + spa34*spb24)*t168*(t621 + t622); d42 = T(1)/d42;
complex<T> d43 = t168*(t344 + s14*t480 + s23*t480 + square(s14) + square(s23) - t491*T(2)); d43 = T(1)/d43;
complex<T> d47 = spa23*t168*(t621 + t622)*(t344 + s14*t480 + s23*t480 + square(s14) + square(s23) - t491*T(2)); d47 = T(1)/d47;
complex<T> d53 = (spa13*spb12 + spa34*spb24)*t168; d53 = T(1)/d53;
complex<T> d59 = spb14*(spa13*spb12 + spa34*spb24)*t168*(t344 + s14*t480 + s23*t480 + square(s14) + square(s23) - t491*T(2)); d59 = T(1)/d59;
complex<T> d60 = spb14*t305; d60 = T(1)/d60;
complex<T> d61 = t168*(t621 + t622); d61 = T(1)/d61;
complex<T> d63 = (-s234 + s34)*t175*square(spa13*spb23 + spa14*spb24); d63 = T(1)/d63;
complex<T> d64 = t175*t210*square(spa13*spb23 + spa14*spb24); d64 = T(1)/d64;
complex<T> d65 = (spa13*spb23 + spa14*spb24)*spb34*t175*square(s234 - s56)*T(2); d65 = T(1)/d65;
complex<T> d68 = (spa13*spb23 + spa14*spb24)*t175*t210; d68 = T(1)/d68;
complex<T> d70 = t662*(t344 + t352 + t353 + t623 + t690 + t729)*T(2); d70 = T(1)/d70;
complex<T> d71 = spa13*t168*(t344 + t352 + t353 + t623 + t690 + t729); d71 = T(1)/d71;
complex<T> d73 = t558*(t344 + t352 + t353 + t623 + t690 + t729)*T(2); d73 = T(1)/d73;
complex<T> d74 = spb24*t168*(t344 + t352 + t353 + t623 + t690 + t729); d74 = T(1)/d74;
complex<T> d79 = t304*t356*T(2); d79 = T(1)/d79;
complex<T> d81 = t304*t356; d81 = T(1)/d81;
complex<T> d82 = spb14*t304*t356; d82 = T(1)/d82;
complex<T> d83 = spb34*t175*t300*T(2); d83 = T(1)/d83;
complex<T> d84 = spb34*t175*t300; d84 = T(1)/d84;
complex<T> d87 = spa23*(t293 + t294)*t305*t371; d87 = T(1)/d87;
complex<T> d88 = (t293 + t294)*t305*t371; d88 = T(1)/d88;
complex<T> d89 = spa23*t305*t356; d89 = T(1)/d89;
complex<T> d90 = t305*t356; d90 = T(1)/d90;
complex<T> d93 = spa34*t300*t494; d93 = T(1)/d93;
complex<T> d95 = t305*t356*T(2); d95 = T(1)/d95;
complex<T> d96 = spa23*spb56*(t293 + t294)*t552; d96 = T(1)/d96;
complex<T> d97 = t494*cube(spa13); d97 = T(1)/d97;
complex<T> d99 = spa56*t168*t551; d99 = T(1)/d99;
complex<T> d102 = spa23*spb56*t168*t552; d102 = T(1)/d102;
complex<T> d111 = (spa13*spb12 + spa34*spb24)*(t344 + t352 + t353 + t623 + t690 + t729); d111 = T(1)/d111;
complex<T> d113 = (t293 + t294)*(t621 + t622)*(t344 + t352 + t353 + t623 + t690 + t729); d113 = T(1)/d113;
complex<T> d127 = (t621 + t622)*square(t344 + t352 + t353 + t623 + t690 + t729); d127 = T(1)/d127;
complex<T> d128 = (spa13*spb12 + spa34*spb24)*square(t344 + t352 + t353 + t623 + t690 + t729); d128 = T(1)/d128;
complex<T> d129 = (spa13*spb12 + spa34*spb24)*(t293 + t294)*(t344 + t352 + t353 + t623 + t690 + t729); d129 = T(1)/d129;
complex<T> d130 = (t293 + t294)*(t621 + t622)*square(t344 + t352 + t353 + t623 + t690 + t729); d130 = T(1)/d130;
complex<T> d131 = (spa13*spb12 + spa34*spb24)*(t293 + t294)*square(t344 + t352 + t353 + t623 + t690 + t729); d131 = T(1)/d131;
complex<T> d138 = (t621 + t622)*(t344 + t352 + t353 + t623 + t690 + t729); d138 = T(1)/d138;
complex<T> d144 = spa56*spb14*t168*t551; d144 = T(1)/d144;
complex<T> d147 = spb56*t168*t552; d147 = T(1)/d147;
complex<T> d149 = t175*t300*T(2); d149 = T(1)/d149;
complex<T> d152 = spa23*t168*t552; d152 = T(1)/d152;
complex<T> d153 = t356*(t621 + t622)*T(2); d153 = T(1)/d153;
complex<T> d154 = (spa13*spb12 + spa34*spb24)*t356*T(2); d154 = T(1)/d154;
complex<T> d155 = spb14*t168*t551; d155 = T(1)/d155;
complex<T> d157 = t168*t305*t371; d157 = T(1)/d157;
complex<T> d159 = t168*t304*t372; d159 = T(1)/d159;
complex<T> t9 = -(d60*(spb14*spb36 + spb13*spb46)*t348) + d45*(d8*spb16*t348*t550 - d46*t339*(spa45*spb46 - t555)*t68 + spa45*(spb56*t296 - spb16*t68) + d57*spb16*spb34*t339*(t299 + t40 + t691) + t625*(-(t135*t359) + t338*t346*T(2) - spb36*t701*T(2)) + d57*spb23*spb46*(d46*(spa24*spb12 + spa34*spb13)*t339*t42 - spb16*(spa23*t731 + spa24*t339*T(2))) - spa15*spb23*(spa24*spb16 + d46*spa12*spb56*t135)*T(4)) + spb23*t527*(-(spa12*t341) + spa23*spb34*t549)*T(6); 
complex<T> t12 = -(d50*(spa24*spa35 + spa23*spa45)*t351) + d45*(t121*t308*t351 + t138*t381*t486 + spa56*spb36*t554 + spa25*spb36*t67 + d52*spa25*spb26*t338*t67 - d52*spa35*spb36*t338*t67 + t299*t339*t486*T(2) + t341*t486*t691*T(2) + d41*(spa25*spa34*t338*(spa56*spb36 + t296 - t346) + spa14*spa35*(d52*spa24*spb34*t338*t41 - d52*t301*t338*t41 + spa25*spb14*t550 + spa25*spb13*t338*T(2))) + d52*spa56*spb12*t138*t855*T(4) - t553*t855*T(4)) + spa14*t195*t647*T(6); 
complex<T> t48 = s14*(d159*s12*spb14*t197 + d100*spb14*t122*t765 - d101*t178*t765*square(spb12)); 
complex<T> t78 = s23*t137; 
complex<T> t133 = d90*spa23*square(s123)*square(spb46) - d89*t179*square(t343); 
complex<T> t139 = d96*(t401*square(spa12)*square(spb46) + square(t343)*square(t348)); 
complex<T> t140 = t576*t671 - d93*t180*cube(t381) + d91*cube(spb13)*square(spa25); 
complex<T> t146 = -(d87*t179*square(spa12)) + d88*spa23*square(t348); 
complex<T> t147 = d97*(t171*t174 - t180*square(spa34)); 
complex<T> t152 = d6*(spa23*t345 - spa34*t789*T(2)); 
complex<T> t204 = spa24*spb16*(d76*spa24*t485 - d75*t553) + d77*spa25*t349*square(spb13); 
complex<T> t232 = d23*t343; 
complex<T> t240 = d6*(d7*spb34*square(spa45) + d8*spa34*square(spb36)); 
complex<T> t256 = t299*t486 + t346*t625 + d72*spb26*t549*(spb34*t299 + spb56*t346 + t880); 
complex<T> t309 = d6*t345; 
complex<T> t310 = t122*t128; 
complex<T> t314 = t123*t129; 
complex<T> t315 = t143*T(3); 
complex<T> t332 = d90*t124; 
complex<T> t337 = d68*t359; 
complex<T> t416 = -(d79*spb14); 
complex<T> t461 = spa15*spb13*spb56 + spb34*t299 + t880; 
complex<T> t462 = spa14*spa25*spb12*spb36 + spa45*t880; 
complex<T> t530 = d4*t340; 
complex<T> t536 = t363*t390; 
complex<T> t561 = d5*t143; 
complex<T> t639 = t130*t132; 
complex<T> t642 = d22*spa23; 
complex<T> t650 = d24*t348; 
complex<T> t670 = d81*spb14; 
complex<T> t675 = d64*spb23; 
complex<T> t696 = t131*t134; 
complex<T> t705 = -(d102*s14); 
complex<T> t706 = -(d144*s23); 
complex<T> t743 = d6*t143; 
complex<T> t777 = t349*t381; 
complex<T> t791 = t171*t629; 
complex<T> t800 = d21*t633; 
complex<T> t801 = t342*t652; 
complex<T> t808 = spb24*t400; 
complex<T> t816 = d83*s12; 
complex<T> t876 = t347*t502; 
complex<T> t888 = spa15*t653; 
complex<T> t978 = spb36*t712; 
complex<T> d2 = s34*spa56*t716; d2 = T(1)/d2;
complex<T> d13 = s34*spb56*t716; d13 = T(1)/d13;
complex<T> d37 = (-s134 + s14)*t788; d37 = T(1)/d37;
complex<T> d38 = t788*square(s134 - s14)*T(2); d38 = T(1)/d38;
complex<T> d39 = (-s134 + s14)*t384*t788; d39 = T(1)/d39;
complex<T> d40 = (-s134 + s34)*t384*t788; d40 = T(1)/d40;
complex<T> d85 = spb34*t300*t400; d85 = T(1)/d85;
complex<T> d98 = t788*square(t384); d98 = T(1)/d98;
complex<T> d115 = (t621 + t622)*t716; d115 = T(1)/d115;
complex<T> d117 = (spa13*spb12 + spa34*spb24)*t716; d117 = T(1)/d117;
complex<T> d146 = t400*cube(spb24); d146 = T(1)/d146;
complex<T> d158 = t788*square(t384)*T(2); d158 = T(1)/d158;
complex<T> t10 = -(d14*(spa13*spb16 + spa23*spb26)*(spb16*spb34 - spb13*spb46)) + d6*(spa45*(spb56*t135 - spb12*t138 + d15*spa23*spb12*spb46*t340) + d8*spa13*spa24*spb16*t374 + d8*spa23*spa24*spb26*t374 + d16*spb14*spb36*t340*t379 - t296*t342*t625 + d16*spb14*spb36*t340*t691 + d16*spa34*spb12*spb36*spb46*t693 + d16*spa24*spb12*spb36*spb46*t340*T(2) - t181*t625*T(2) + spa35*spb12*t821*T(4) + d15*spa23*spb12*(-(spa35*spb36*t340) + d16*t340*t430*square(spb46) + spa35*spb56*t296*T(4))) - spb12*t561*(spa23*t345 - spa34*t789*T(2))*T(6); 
complex<T> t11 = d18*(spa25*spa34 + spa24*spa35)*(spa15*spb14 + spa25*spb24) + d6*(spb36*(spa12*t135 - spa56*t138) - d15*spb14*(spa35*spb36 - spa45*spb46)*t141 + d19*spa23*spa45*t340*(spa56*spb36 - t346 + t38) - (spa15*spb14 + spa25*spb24)*t301*t625 + t486*(t343*t40 + t182*T(2)) - d19*spa12*spa35*(spa24*spa45*t144 + d15*spa35*spb14*t340*t430 + spa45*spb13*t340*T(2)) - spa12*spb46*(spa45*spb13 + d15*spa56*spb14*t554)*T(4)) + spa12*t561*(spb14*t345 - spa23*t144*T(2))*T(6); 
complex<T> t13 = t137*(d143*s14*s56 + d78*t339) + d80*t128*t197*(d143*s14*s56 + d78*t339) - d49*spa24*(spa35*spb36*(s234*t339 + s56*t341)*t354 - s56*spb13*(spa15*spb26*t341 + spa35*t306*t354) + spb16*t354*t44*(-(spa12*t341) + spa23*spb34*t549)) + d148*(spa25*spb12 + t135)*t821 + d140*(t299*t346 + d78*(spa25*spb26 + spa35*spb36)*(t120 + t731) + d15*t430*(spa45*spb46*t294 - t294*t555 + t695*t847)) - t624*(s56*t135*t138 + t138*t339*t793 - t695*t875)*t876 + d142*t347*(t37*t38 + d78*spa56*spb36*(-(spa45*spb56) + t37*T(2))) + d112*t343*(-(d141*spa45*spb13*t348) + t527*t624*(spa23*spb34*t339 + t338*t789)*T(3) + d45*(spa45*spb13*(t306*t338 - t341*t45) + spb16*t120*(spa15*t339 - spa45*(t621 + t622) + t306*t661) + d46*(s12 + s13 - s24 - s34)*spa15*t68*(spa24*spb12*spb46 + spb13*t768) - (spa24*spb12 + spa34*spb13)*(spa14*spb34*(spa25*spb26 + spa35*spb36) - spa45*spb46*t68*T(3)))); 
complex<T> t14 = d139*spa45*spb13*(-(spa12*spb16) + t138) + d142*t347*(d78*spa45*spb56*(-(spa56*spb36) + spa25*t354) + t37*t38) + t139*(d78*t338 + d141*t640) - d96*t314*(d78*t338 + d141*t640) + d140*(t299*t346 + d78*(spa45*spb46 + t555)*t632 + d15*t430*(-(spa25*spb26*t293) + spa35*spb36*t293 + t695*t847)) + t624*(-(s56*t135*t138) + spa25*spb56*t135*t338 + t695*t875)*t876 + d49*spb13*(-(spa25*t195*t45*t549) + s134*spa45*t306*t338*T(2) + s56*(spa15*spa24*spb26*t341 + spa24*spa35*t306*t354 + spa45*t306*t341*T(2))) - d108*t342*(-(d143*t351*t821) - (spa34*spb14*t338 + spa23*spb12*t339)*t624*t647*T(3) + d45*(-(d52*spa35*spb26*(-(spa24*spb34) + t301)*(-s12 + s34 + t418)*t67) + spa25*t121*(spb26*t338 - spb36*t373 + t354*t736) + (spa35*spb23*t339 + t341*t44)*t821 - (spa24*spb34 - t301)*(spa45*spb46*t303 + t303*t555 - spa35*spb36*t67*T(3)))); 
complex<T> t54 = t357*(-t133 - t139 + t146 + d96*t314 + spa23*t125*t332 + d88*t503); 
complex<T> t59 = (s256 - s34)*(-t147 - d98*t791 - d93*t180*cube(t381) + d91*cube(spb13)*square(spa25)); 
complex<T> t136 = t670*square(s124)*square(spa35) - d82*t178*square(t342); 
complex<T> t145 = d146*(t169*t173 - t177*square(spb34)); 
complex<T> t155 = -(d1*spa24*spb12*t178) + d29*spa24*spb14*t197 + d28*spa24*t197*t342 + s124*spa35*spb13*t801 - s124*spa35*t749*t930 + d26*spb13*t402*square(spa35); 
complex<T> t159 = d12*t179*t301 + spa24*spb46*(d25*spb46*t401 - t633*t650 + spa45*t800) + spb13*(t232 - t642)*square(t348); 
complex<T> t282 = s134*(t140*t765 - d94*s12*s134*t130*t381*square(spa14)); 
complex<T> t329 = d37*(-(spa13*spb16) + t768); 
complex<T> t506 = -(d80*t310); 
complex<T> t534 = t309*t365; 
complex<T> t563 = (-s13 - s14 + s23 + s24)*t309; 
complex<T> t564 = d78*t133*t76 + d95*s34*t392*cube(s123); 
complex<T> t654 = d72*t461; 
complex<T> t683 = s23*(d157*s12*t503 + t146*t765); 
complex<T> t709 = d85*t177; 
complex<T> t725 = s14*(t147*t697 + d158*s34*t791); 
complex<T> t751 = t337*t485; 
complex<T> t756 = s124*(t310*t706 + d78*t78); 
complex<T> t773 = t143*t152; 
complex<T> t794 = d78*t139; 
complex<T> t877 = t484*t561; 
complex<T> d54 = t808*square(s23 - s234)*T(2); d54 = T(1)/d54;
complex<T> d55 = (s23 - s234)*t385*t808; d55 = T(1)/d55;
complex<T> d56 = (s23 - s234)*t808; d56 = T(1)/d56;
complex<T> d66 = (-s234 + s34)*t385*t808; d66 = T(1)/d66;
complex<T> d145 = t808*square(t385); d145 = T(1)/d145;
complex<T> d160 = t808*square(t385)*T(2); d160 = T(1)/d160;
complex<T> t4 = d29*spa24*spb14*t122 + d53*t12*t342 + d28*spa24*t122*t342 + d47*t26*t343 - d42*spa24*spa35*t70 + d38*t174*t810 + d39*t174*t810 + t20*T(2) - t329*t349*t731*T(2) + d44*spa14*t347*(spb36*t135*t338 + spb16*t339*t38 + spb36*t341*t793 - spb56*t135*t38*T(2))*T(3) + d43*(spa14*spa25*spb15*t200 + spa24*t200*t338 + spa26*t200*t37 + spa45*t154*t347*T(2) + t550*t695*T(4)); 
complex<T> t22 = d126*spa24*spb13*spb16*t349 + d125*spb16*t170*t349 - d110*spa25*t172*t349 + d108*spa24*spa35*t374 + d107*spa24*(spa25*spb12 + t135)*t374 - d104*spb16*t170*t485 + d106*spa25*t172*t485 + d105*t394*t553 - s134*t130*(d118*s34*s56 + d78*t340)*t381*t576 + d112*spa45*spb46*t693 + d109*spa45*(-(spa12*spb16) + t138)*t693 + d84*t169*(d103*s34*s56 + d78*t340)*t696 + (d103*s34*s56 + d78*t340)*(d84*t134*t169*t569 + t568*t612 + t231*t709) + d131*s12*t43*t736*(spa56*spb12*(spb16*t141 + spa25*spb56*t345) - s123*(spa25*spb12*t340 + t345*t793)) + d130*s12*spa56*spb46*t43*(s124*(spb16*t141 + spa25*spb56*t345) - spa12*spb56*(spa25*spb12*t340 + t345*t793)) + (d118*s34*s56 + d78*t340)*(s134*t130*t381*t576 - d93*t180*cube(t381) + d91*cube(spb13)*square(spa25)) + d128*t430*t44*(t279*t531*t67 + spa24*spb12*(-(spb16*t141) + t345*t761))*T(3) + d127*t430*t45*(spa25*spb12*t301*t340 + t149*(t298*t340 - t345*t44)*t68 + t301*t345*t793)*T(3) + d129*t41*(s124*spb36*(spa23*t430 + spa24*t496 + spa34*t693) + t430*(s123*spa12*spb16 - s124*spa12*spb16 - s123*spa24*spb46 + s124*t138 - spb16*t141 + t345*t761) - d15*s124*spb46*(t120 + t121)*(spa24*(t293 + t294) - spa23*t495 + spa34*t789) - spa12*spb56*(s13*spa45*spb13 + s23*spa45*spb13 - spa45*spb13*t496 + spa25*spb12*t430*T(3) + t135*t430*T(3))) + d113*t42*(s123*spa45*(spa24*t144 + spb14*t430 + spb13*t495) - d15*s123*spa35*(t120 + t121)*(spa23*t144 + spb13*(t293 + t294) - spb14*t496) + t430*(s123*(spa25*spb12 + t135) + spa25*spb12*t340 - s124*t362 + t345*t793) + spa56*spb12*((-s14 + s23 + t418)*t821 + spa12*spb16*t430*T(3) - t138*t430*T(3))) + d138*spa15*(-(spb13*(spa45*spb56*t301 + spa24*spb56*t38 + spa12*spb16*t430 - t430*t761)) + spb23*t149*(t355*t430 + d3*(s13 + s14 - s23 - s24)*t45*t624 - t45*t550*T(3))) - d111*spb26*(spa24*(-(spa56*spb13*t37) + t430*(spa25*spb12 + t793)) + spa56*t374*square(spa24) + spa14*t149*(t298*t430 + d3*t365*t44*t624 + t44*t731*T(3))) + d137*(-(d135*(s134 + s23)*spa12*spb36*(spa25*spb23 - t346)) - d136*t130*t568*(t621 + t622) + d132*t172*t569*(s24*spa12 + t141 - spa24*t622) + d134*spa24*spa45*(s24*t355 + spa24*spb26*t622 - t37*(t621 + t622)) - d133*(s56*t295 + s234*t340)*t359*square(t485) + d114*spa12*spa45*(spb36*t350 + s24*spb36*T(2) + spb13*t45*T(2)) - d103*t359*t550*t555*T(4)) + d124*(-(d122*t130*t170*(s13*spb12 - spa34*spb13*spb24 + spb12*t340)) + d119*(s14 + s234)*spa45*spb12*(t299 + t37) + d121*(spa13*spa25*spb12*spb23 + s13*t298 + spa25*spb24*t303 - spa34*spb24*t346)*t930 + d120*t131*t373*cube(spb13) + d123*(s56*t295 + s134*t340)*t381*square(t37 + t554) + d116*t374*(spa45*t350 + s13*spa45*T(2) - spa24*t44*T(2)) + d118*t381*t553*t855*T(4)) - d117*(s14*spa35*spb23*t37 + s13*t351*t554 + s34*t555*t67 + s24*(spa35*spb23*t299 - t351*t554 + spa35*spb36*t67) + d116*spa23*spa25*spb26*t144*t67*T(2) + d116*spa23*spa45*spb46*t144*t67*T(2) + spa23*spb26*t39*t67*T(2) - spa25*spb26*t495*t67*T(2) - s13*spa45*spb46*t67*T(3) - spa25*spb24*t554*t67*T(3) - spa13*spa35*t374*t550*T(5)) + d115*(s23*spa25*spb23*t306 + t296*t348*t418 - s34*spa25*spb26*t68 - s13*(t306*t346 + spa45*spb46*t68) - spa15*spb14*t554*t68*T(2) + t496*t555*t68*T(2) - d114*spa34*spa35*spb36*t68*t789*T(2) - d114*spa34*t555*t68*t789*T(2) + s24*spa35*spb36*t68*T(3) + spa13*spb16*t39*t68*T(3) + spa45*spb13*t306*t621*T(5)); 
complex<T> t23 = s12*s234*t240 + t299*t486 + t325*t624 + spb12*t530*(spa12*spa45*t38 - spa15*(spa24*t38 + spa25*t430) + t131*t149*t624) + d9*t141*(spb13*spb26*t379 + t374*(t299 + t40) - (t130*t149 + t399)*t624) - spa12*spa45*(spb16*spb23 + spb13*spb26)*t309*T(2) + t365*t877*T(6); 
complex<T> t24 = s12*s134*t240 + t325*t624 + d9*t141*(spb13*spb26*t37 + t37*t374 - spb16*spb26*t430 + t130*t149*t624) + t346*t625 - spb12*t530*(spa12*spa45*(t296 - t346) + spa15*spa24*(t346 - t38) + t624*(t131*t149 + square(spa25))) + (spa15*spa24 + spa14*spa25)*t309*t374*T(2) + (-s13 - s14 + s23 + s24)*t877*T(6); 
complex<T> t47 = s124*t136*t697 + s34*t127*t416*cube(s124); 
complex<T> t49 = (s12 - s124)*(-t136 - t137 + d100*spb14*(t122 + t197) + d80*t310 + t313*t670 - d101*t178*square(spb12)); 
complex<T> t55 = d141*s12*s56*(t133 + d90*spa23*t125*t401) + d143*s12*s56*(t136 + t127*t402*t670) + d78*t350*(t133 + t136 + d90*spa23*t125*t401 + t127*t402*t670); 
complex<T> t56 = d150*s134*spb34*t130*t171*t381 + d95*s34*spa23*t125*t401 + d79*s34*spb14*t127*t402 + d149*spa34*t134*t169*t569 + t133*t697 + t136*t697 + d84*t134*t169*t569*t697 + s134*t130*t381*t576*t697 + t568*t612*t697 + t231*t697*t709 - d93*t180*t697*cube(t381) + d91*t697*cube(spb13)*square(spa25); 
complex<T> t58 = -(d78*s14*t137) + s14*t147 + d99*spa14*t128*t197 + d100*s14*spb14*(t122 + t197) + t314*t705 + d98*s14*t791 + s14*t794 - d101*s14*t178*square(spb12); 
complex<T> t62 = d155*spb56*t128*t197 + t133*t312 + t136*t312 + t137*t312 + t139*t312 - d152*spa56*t314 - d151*s134*spa56*t130*t171*t381 + d153*spa23*spa56*t125*t401 + d154*spb14*spb56*t127*t402 + d84*t134*t169*t312*t569 + s134*t130*t312*t381*t576 + t312*t568*t612 + d156*spb56*t169*t696 + t231*t312*t709 - d93*t180*t312*cube(t381) + d91*t312*cube(spb13)*square(spa25); 
complex<T> t80 = s12*(d84*s234*(spa12*spb23 - spa14*spb34)*t169*t569 + t568*t612 + t231*t709); 
complex<T> t90 = -(t239*t350) + t412*t430 + t315*t534 + t606*t632 + spa14*spa25*t374*T(4) - d3*s56*t412*t624*T(4) + spa45*t880*T(4) - spa15*spa45*spb23*spb56*t792*T(8); 
complex<T> t99 = -(t256*t350) + t315*t563 + t606*(t120 + t731) + t458*T(4) + (t298*t308 + d8*spb16*t355)*(t430 - d3*s56*t624*T(4)) + spa56*spb36*t792*t855*T(8); 
complex<T> t164 = spb13*(t328 + t329)*t349*t549 - d32*spa14*spa24*t671 + d33*t170*t671 + d38*t791 + d39*t791 + d40*t791 - d36*t172*t180*t895 + t582*t720*square(spb26) + spa24*t510*t777*T(2); 
complex<T> t166 = s23*(t145*t697 - d160*s34*t173*t785); 
complex<T> t326 = -(d56*(spa25*spb24 + spa35*spb34)); 
complex<T> t504 = d78*t136; 
complex<T> t514 = -(d145*t173); 
complex<T> t680 = s123*(t314*t705 + s14*t794); 
complex<T> t851 = spb26*t654; 
complex<T> t937 = t144*t773; 
complex<T> t1002 = t346*t751; 
complex<T> t1 = -(spb13*t123*t232) - d59*t27*t342 + t326*t354*t394 + spb13*t123*t642 - d42*spb13*spb46*t71 + d54*t173*t785 + d55*t173*t785 - d61*t343*t9 + d49*t25*T(2) + d44*spb23*t347*(-(t37*(spa25*t338 + t138*t661)) + spa45*(t138*t339 + t341*t761))*T(3) + d43*(t196*t201 - t978*T(2) + t989*T(4)); 
complex<T> t2 = d17*spa12*t10 - d20*spb12*t11 + d1*spa24*spb12*t178 - d10*spb23*t23 - d11*spa14*t24 - d13*spb46*t349*t430 + d2*spa35*t430*t485 + d12*t179*t693; 
complex<T> t3 = -(d17*spa12*t10) + d20*spb12*t11 + d25*spa24*t124*t125 + d10*spb23*t23 + d11*spa14*t24 + d26*spb13*t313 - d53*t12*t342 + d59*t27*t342 - d47*t26*t343 + d32*s134*spa14*spa24*t130*t381 - d33*s134*t130*t170*t381 + d13*spb46*t349*t430 - d2*spa35*t430*t485 + spa24*spb46*t633*t650 + spb13*t134*t569*t675 + d65*t172*t696 + d42*spa24*spa35*t70 + d42*spb13*spb46*t71 - spa24*spa45*spb46*t800 - s124*spa35*spb13*t801 + spb23*(-(t37*(spa25*t338 + t138*t661)) + spa45*(t138*t339 + t341*t761))*t876 + d61*t343*t9 + d70*spb23*t90 + s124*spa35*t749*t930 + d73*spa14*t99 + t1002*T(2) - t20*T(2) - d49*t25*T(2) - spa24*t510*t777*T(2) + spa14*t876*(spb36*t135*t338 + spb16*t339*t38 + spb36*t341*t793 - spb56*t135*t38*T(2)) + d78*(d77*spa25*t172*t349 + d76*spb16*t170*t485 - d75*spa24*spb16*t553)*T(3) - d43*(spa14*spa25*spb15*t200 + spa24*t200*t338 + spa26*t200*t37 + spa45*t154*t347*T(2) + t550*t695*T(4)) - d43*(t196*t201 - t978*T(2) + t989*T(4)) + d71*(d15*spa23*t191*(d7*spa35*t296 + t486*t768) + d8*t40*(t343*t379 + t182*T(2)) - t625*(s34*spa23*(-t346 + t38) + spa34*(-t346 + t38)*t789 - t693*(spa35*t340 + t661*t768)*T(2) + spa23*t296*t350*T(3)) + spa12*t153*t743*T(6)) + d74*(d15*spb14*t192*(d7*spa35*t296 + t486*t768) + d7*t39*(t342*(-t346 + t38) + t181*T(2)) + t486*(s34*spb14*t379 + spa23*t144*t379 + spa24*spb12*T(2)*(spb46*t340 + spb34*t736*T(2)) - spb14*t350*t554*T(3)) + t937*T(6)); 
complex<T> t15 = d63*s234*spa24*t169*t569 - t130*t582*t720 + d66*t173*t785 + d40*t174*t810 + d36*t172*t180*t895 + d62*t170*t177*t897 - d70*spb23*t90 - d73*spa14*t99 + t324*t485*t550*T(2) - t328*t349*t731*T(2) + d71*(-(d15*spa23*t191*(d7*spa35*t296 + t486*t768)) + spa34*t486*(t343*t379 + t182*T(2)) + t625*(s34*spa23*(-t346 + t38) + spa34*(-t346 + t38)*t789 - t693*(spa35*t340 + t661*t768)*T(2) + spa23*t296*t350*T(3)) - spa12*t153*t743*T(6)) + d74*(-(d15*spb14*t192*(d7*spa35*t296 + t486*t768)) + d7*t296*(t342*(-t346 + t38) + t181*T(2)) - t486*(s34*spb14*t379 + spa23*t144*t379 + spa24*spb12*T(2)*(spb46*t340 + spb34*t736*T(2)) - spb14*t350*t554*T(3)) - t937*T(6)); 
complex<T> t50 = -((s156 - s34)*(t145 - t568*t612 - t231*t709 + t514*t785)); 
complex<T> t53 = d147*spb23*t314 + d78*t78 + s23*(t133 + t145 + d144*t128*t197 + d90*spa23*t125*t401 + t514*t785 + t794); 
complex<T> t165 = t324*t354*t394 + d65*s234*t172*t359*t569 - d54*t173*t785 - d55*t173*t785 - d66*t173*t785 - d62*t170*t177*t897 + d63*s234*spa24*t169*square(spa15) + s234*spb13*t359*t675*square(spa15) - t1002*T(2) + t326*t485*t550*T(2); 
complex<T> t726 = d95*s12*spa23*t124*t125 + s12*(t137 + t146 - t204 + t313*t416 + d88*t503 + t504 + t506 + d94*t171*t639) - t133*t765 + t140*t765 + d78*t80 + t169*t696*t816; 
complex<T> t757 = d78*s234*t80 + t131*t169*t311*t359*t816; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t159*(*CI_users[1]->get_value(mc,ind,mu)) + t155*(*CI_users[2]->get_value(mc,ind,mu)) + t164*(*CI_users[3]->get_value(mc,ind,mu)) + t4*(*CI_users[4]->get_value(mc,ind,mu)) + t1*(*CI_users[5]->get_value(mc,ind,mu)) + t165*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + t3*(*CI_users[8]->get_value(mc,ind,mu)) + t726*(*CI_users[9]->get_value(mc,ind,mu)) + t54*(*CI_users[10]->get_value(mc,ind,mu)) + t59*(*CI_users[11]->get_value(mc,ind,mu)) + t58*(*CI_users[12]->get_value(mc,ind,mu)) + t22*(*CI_users[13]->get_value(mc,ind,mu)) + t14*(*CI_users[14]->get_value(mc,ind,mu)) + t53*(*CI_users[15]->get_value(mc,ind,mu)) + t50*(*CI_users[16]->get_value(mc,ind,mu)) + t13*(*CI_users[17]->get_value(mc,ind,mu)) + t56*(*CI_users[18]->get_value(mc,ind,mu)) + t55*(*CI_users[19]->get_value(mc,ind,mu)) + t49*(*CI_users[20]->get_value(mc,ind,mu)) + t62*(*CI_users[21]->get_value(mc,ind,mu)) + t757*(*CI_users[22]->get_value(mc,ind,mu)) + t282*(*CI_users[23]->get_value(mc,ind,mu)) + t683*(*CI_users[24]->get_value(mc,ind,mu)) + t756*(*CI_users[25]->get_value(mc,ind,mu)) + t725*(*CI_users[26]->get_value(mc,ind,mu)) + t47*(*CI_users[27]->get_value(mc,ind,mu)) + t48*(*CI_users[28]->get_value(mc,ind,mu)) + t680*(*CI_users[29]->get_value(mc,ind,mu)) + t564*(*CI_users[30]->get_value(mc,ind,mu)) + t166*(*CI_users[31]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmmpemep_SLC_wCI::\
C2q2g2l_qpqmmpemep_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c134, c256));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c256));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c2356));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c14, c23, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c23, c14, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c14, c56));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c256));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c356));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c23, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmmpemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmmpemep SLC");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:39:13 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa56 = SPA(5,6);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb16 = SPB(1,6);
complex<T> spb46 = SPB(4,6);
complex<T> spb56 = SPB(5,6);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb26 = SPB(2,6);
complex<T> spb36 = SPB(3,6);
complex<T> s12 = -(spa12*spb12);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s34 = -(spa34*spb34);
complex<T> s124 = SS(1,2,4);
complex<T> s134 = SS(1,3,4);
complex<T> s234 = SS(2,3,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s24 = -(spa24*spb24);
complex<T> s156 = SS(1,5,6);
complex<T> s256 = SS(2,5,6);
complex<T> t22 = -(spa34*spb46); 
complex<T> t23 = spa15*spa23; 
complex<T> t24 = spb14*spb26; 
complex<T> t44 = s12*s134; 
complex<T> t70 = square(spa15); 
complex<T> t71 = square(spb26); 
complex<T> t72 = square(spa15*spb13 + spa25*spb23); 
complex<T> t73 = square(spa14*spb16 + spa24*spb26); 
complex<T> t79 = -(spa12*spa35); 
complex<T> t80 = spb12*spb46; 
complex<T> t81 = -(spa34*T(2)); 
complex<T> t98 = (spa13*spb23 + spa14*spb24)*T(2); 
complex<T> t99 = square(spa25*spb12 + spa35*spb13); 
complex<T> t100 = square(spa12*spb16 - spa24*spb46); 
complex<T> t101 = spa12*spb14; 
complex<T> t102 = spa23*spb12; 
complex<T> t119 = (-s134 + s56)*spa34; 
complex<T> t120 = (-s234 + s56)*spa56; 
complex<T> t172 = -(spa13*spb16); 
complex<T> t173 = -(spa25*spb24); 
complex<T> t178 = -(spa23*spb26); 
complex<T> t180 = s12*s34; 
complex<T> t182 = -(spa35*spb36); 
complex<T> t186 = s134*(spa12*spb13 + spa24*spb34); 
complex<T> t187 = spa45*spb46; 
complex<T> t191 = spa56*spb23; 
complex<T> t201 = (-s134 + s34)*spb56; 
complex<T> t204 = s12 - s34 - s56; 
complex<T> t205 = spa24*spb12 + spa34*spb13; 
complex<T> t206 = -s12 + s34 - s56; 
complex<T> t208 = square(spa13); 
complex<T> t209 = square(spb24); 
complex<T> t211 = spa35*spb23 + spa45*spb24; 
complex<T> t212 = spa13*spb36 + spa14*spb46; 
complex<T> t213 = -(spa15*spb16) - spa25*spb26; 
complex<T> t214 = -(spa12*spb13) - spa24*spb34; 
complex<T> t215 = -s12 - s34 + s56; 
complex<T> t216 = -(spa15*spb13) - spa25*spb23; 
complex<T> t217 = -(spa14*spb16) - spa24*spb26; 
complex<T> t218 = square(s12); 
complex<T> t219 = square(s34); 
complex<T> t220 = square(s56); 
complex<T> t221 = spa15*spb14; 
complex<T> t222 = square(spa23); 
complex<T> t223 = square(spa25); 
complex<T> t224 = square(spb16); 
complex<T> t225 = s13 + s14 - s23 - s24; 
complex<T> t226 = spa23*spb13 + spa24*spb14; 
complex<T> t228 = -(spb56*T(2)); 
complex<T> t310 = -(s56*T(2)); 
complex<T> t312 = cube(spa13); 
complex<T> t313 = cube(spb24); 
complex<T> t316 = square(spb14); 
complex<T> t322 = -(spa12*spb26); 
complex<T> t324 = square(s134); 
complex<T> t325 = -(spa25*spb12); 
complex<T> t326 = s124*spb36; 
complex<T> t359 = cube(spa13*spb23 + spa14*spb24); 
complex<T> t360 = spa14*spb56; 
complex<T> t366 = spa13*spb26; 
complex<T> t367 = spa15*spb24; 
complex<T> t371 = spa12*spb36; 
complex<T> t400 = spa14*spb13; 
complex<T> t401 = spa24*spb23; 
complex<T> t406 = spa23*spb34; 
complex<T> t410 = spa25*spb16; 
complex<T> t445 = s34*s56; 
complex<T> t446 = spa35*spb34; 
complex<T> t447 = spa45*spb14; 
complex<T> t452 = spa56*spb36; 
complex<T> t470 = s12*T(2); 
complex<T> t491 = spa13*spb14; 
complex<T> t492 = spa23*spb24; 
complex<T> t499 = spa45*spb34; 
complex<T> t522 = spb16*spb24; 
complex<T> t539 = spb12*spb34; 
complex<T> t545 = spa13*spa25; 
complex<T> t556 = spa13*spb24; 
complex<T> t559 = s12*spa35; 
complex<T> t579 = spb14*spb16; 
complex<T> t587 = spa23*spa25; 
complex<T> d3 = spa13*spb23 + spa14*spb24; d3 = T(1)/d3;
complex<T> d7 = spa56; d7 = T(1)/d7;
complex<T> d8 = spb56; d8 = T(1)/d8;
complex<T> d16 = spb34*spb56; d16 = T(1)/d16;
complex<T> d19 = spa34*spa56; d19 = T(1)/d19;
complex<T> d30 = s134; d30 = T(1)/d30;
complex<T> d37 = s234; d37 = T(1)/d37;
complex<T> d39 = spa56*(spa13*spb23 + spa14*spb24); d39 = T(1)/d39;
complex<T> d42 = (spa13*spb23 + spa14*spb24)*spb56; d42 = T(1)/d42;
complex<T> d48 = T(2); d48 = T(1)/d48;
complex<T> d68 = s234*spb23*(spa13*spb23 + spa14*spb24); d68 = T(1)/d68;
complex<T> d69 = s134*spa14*(spa13*spb23 + spa14*spb24); d69 = T(1)/d69;
complex<T> d78 = s124*spa14*(spa12*spb13 + spa24*spb34); d78 = T(1)/d78;
complex<T> d80 = spb23; d80 = T(1)/d80;
complex<T> d82 = spa14; d82 = T(1)/d82;
complex<T> d85 = spa34; d85 = T(1)/d85;
complex<T> d86 = spa14*spa34; d86 = T(1)/d86;
complex<T> d87 = spa56*spb34; d87 = T(1)/d87;
complex<T> d92 = spb34; d92 = T(1)/d92;
complex<T> d94 = spb23*spb34; d94 = T(1)/d94;
complex<T> d95 = spa34*spb56; d95 = T(1)/d95;
complex<T> d97 = s124; d97 = T(1)/d97;
complex<T> d99 = s123; d99 = T(1)/d99;
complex<T> t21 = -t446; 
complex<T> t25 = spa12*t204; 
complex<T> t68 = -t491; 
complex<T> t69 = -t492; 
complex<T> t77 = -(t213*(t491 + t492)); 
complex<T> t86 = -(d3*(spa23*spb13 + spa24*spb14)); 
complex<T> t97 = square(t400 + t401); 
complex<T> t103 = spb23*t205; 
complex<T> t105 = spb13*t215 + spa24*t539*T(2); 
complex<T> t113 = (-s13 + s14 - s23 + s24)*t206 - t218 - t219 - t220 + s56*t470 + t180*T(2) + t445*T(2); 
complex<T> t115 = spa56*spb16*t215 + t206*t447; 
complex<T> t116 = spa23*spb36*t206 - spa25*spb56*t215; 
complex<T> t177 = -t221; 
complex<T> t179 = cube(t400 + t401); 
complex<T> t185 = s234*t205; 
complex<T> t188 = d7*spb34*square(spa35) + d8*spa34*square(spb46); 
complex<T> t207 = spa34*spb46 + t172; 
complex<T> t227 = -t400 - t401; 
complex<T> t243 = -(spa34*spb14) + t102; 
complex<T> t275 = spa15*spb14*spb56 + spb14*t322 + spb26*t406; 
complex<T> t300 = t522*t79 + spa15*spa23*t80; 
complex<T> t302 = spb14*spb26*t79 + t545*t80; 
complex<T> t315 = t172 + t178; 
complex<T> t321 = d48*s12; 
complex<T> t323 = -(t205*t72); 
complex<T> t332 = t214*t73; 
complex<T> t340 = d3*spa56; 
complex<T> t365 = (t400 + t401)*T(2); 
complex<T> t402 = -(t180*T(2)); 
complex<T> t403 = d8*spb46; 
complex<T> t404 = d7*spa35; 
complex<T> t408 = s234*t70; 
complex<T> t418 = d3*spb56; 
complex<T> t475 = -(t223*t316); 
complex<T> t477 = d7*spa15; 
complex<T> t478 = d8*spb16; 
complex<T> t533 = t366*T(2); 
complex<T> t567 = -(t367*T(2)); 
complex<T> t596 = t191*t98; 
complex<T> t620 = spb24*t221; 
complex<T> t623 = t360*t98; 
complex<T> t625 = s23*t222; 
complex<T> d2 = s34*spa56*(t400 + t401)*t98; d2 = T(1)/d2;
complex<T> d10 = s34*spb23*t98; d10 = T(1)/d10;
complex<T> d11 = s34*spa14*t98; d11 = T(1)/d11;
complex<T> d13 = s34*spb56*(t400 + t401)*t98; d13 = T(1)/d13;
complex<T> d14 = spb34*spb56*(t400 + t401); d14 = T(1)/d14;
complex<T> d15 = t400 + t401; d15 = T(1)/d15;
complex<T> d18 = spa34*spa56*(t400 + t401); d18 = T(1)/d18;
complex<T> d22 = spb23*(t400 + t401)*square(s123 - s56); d22 = T(1)/d22;
complex<T> d23 = spa14*(t400 + t401)*square(s124 - s56); d23 = T(1)/d23;
complex<T> d25 = t119*t360*square(spa13*spb23 + spa14*spb24); d25 = T(1)/d25;
complex<T> d28 = (spa13*spb23 + spa14*spb24)*t119*t360; d28 = T(1)/d28;
complex<T> d29 = s134*spa14*(spa13*spb23 + spa14*spb24)*t201; d29 = T(1)/d29;
complex<T> d31 = spa14*t201*t98; d31 = T(1)/d31;
complex<T> d33 = spb23*spb34*t120*square(spa13*spb23 + spa14*spb24); d33 = T(1)/d33;
complex<T> d35 = s234*(-s234 + s34)*(spa13*spb23 + spa14*spb24)*t191; d35 = T(1)/d35;
complex<T> d36 = spb23*(spa13*spb23 + spa14*spb24)*spb34*t120; d36 = T(1)/d36;
complex<T> d45 = s234*spb23*t445; d45 = T(1)/d45;
complex<T> d46 = s134*spa14*t445; d46 = T(1)/d46;
complex<T> d47 = spa14*spb23*t445; d47 = T(1)/d47;
complex<T> d53 = spb34*t191*t359*T(2); d53 = T(1)/d53;
complex<T> d54 = spb34*t191*t359; d54 = T(1)/d54;
complex<T> d55 = s234*spb34*t191*t359; d55 = T(1)/d55;
complex<T> d59 = spa56*spb34*t186*t226; d59 = T(1)/d59;
complex<T> d60 = spa34*t359*t360; d60 = T(1)/d60;
complex<T> d61 = s134*spa34*t359*t360; d61 = T(1)/d61;
complex<T> d62 = spa34*t359*t360*T(2); d62 = T(1)/d62;
complex<T> d64 = (spa12*spb13 + spa24*spb34)*t360*(t400 + t401); d64 = T(1)/d64;
complex<T> d65 = spa56*spb13*t186; d65 = T(1)/d65;
complex<T> d71 = spa34*t186*T(2); d71 = T(1)/d71;
complex<T> d72 = t205*square(s234); d72 = T(1)/d72;
complex<T> d73 = spa56*spb34*t186*T(2); d73 = T(1)/d73;
complex<T> d79 = (spa12*spb13 + spa24*spb34)*t324; d79 = T(1)/d79;
complex<T> d81 = t205*(t400 + t401)*t98; d81 = T(1)/d81;
complex<T> d83 = (spa12*spb13 + spa24*spb34)*(t400 + t401)*t98; d83 = T(1)/d83;
complex<T> d88 = spa34*t360; d88 = T(1)/d88;
complex<T> d89 = (spa12*spb13 + spa24*spb34)*t98; d89 = T(1)/d89;
complex<T> d93 = spb34*t191; d93 = T(1)/d93;
complex<T> d96 = t205*t98; d96 = T(1)/d96;
complex<T> d98 = spa14*(t400 + t401); d98 = T(1)/d98;
complex<T> d100 = spb23*(t400 + t401); d100 = T(1)/d100;
complex<T> d101 = spa14*(spa12*spb13 + spa24*spb34)*t191*t226*T(2); d101 = T(1)/d101;
complex<T> d107 = t191*t359*T(2); d107 = T(1)/d107;
complex<T> d108 = t359*t360*T(2); d108 = T(1)/d108;
complex<T> d109 = spa14*spa34*t359*T(2); d109 = T(1)/d109;
complex<T> d114 = spb23*spb34*t359*T(2); d114 = T(1)/d114;
complex<T> d115 = spa56*spb13*t186*T(2); d115 = T(1)/d115;
complex<T> t40 = t214*t315 + t116*T(2); 
complex<T> t43 = s34*spb13*t315 - spa24*t315*t539 + spb36*t102*t204*T(2) + t102*t228*t499*T(2) + spa34*spb13*spb46*t206*T(3); 
complex<T> t52 = -(s14*s23*t224) + square(-(spa12*spb12*spb16) + spa24*spb12*spb46 + spb13*t207); 
complex<T> t83 = d60*s134; 
complex<T> t84 = d54*s234; 
complex<T> t106 = d30*t207 + d3*t533; 
complex<T> t183 = t77*T(3); 
complex<T> t198 = d3*(t325*t477 + t322*t478); 
complex<T> t229 = d7*spa45*t21 + d8*spb36*t22; 
complex<T> t252 = t325*t477 + spa12*spb26*t478; 
complex<T> t253 = t492 + t68; 
complex<T> t290 = spa15*spa23*spb12 + spa34*t177 + spa56*t178; 
complex<T> t309 = t173 + t21; 
complex<T> t311 = t68 + t69; 
complex<T> t314 = t173 + t177; 
complex<T> t346 = d28*t207; 
complex<T> t378 = t491 + t69; 
complex<T> t422 = d35*t211; 
complex<T> t423 = d29*t212; 
complex<T> t425 = d42*t275; 
complex<T> t432 = d62*t71; 
complex<T> t540 = s12*t340; 
complex<T> t580 = d25*spa23; 
complex<T> t590 = d33*spb14; 
complex<T> t612 = t418*t559; 
complex<T> d1 = (s12 - s123)*t191*t97; d1 = T(1)/d1;
complex<T> d4 = spa56*(t218 + t219 + t220 + s12*t310 + s34*t310 + t402); d4 = T(1)/d4;
complex<T> d5 = square(t218 + t219 + t220 + s12*t310 + s34*t310 + t402); d5 = T(1)/d5;
complex<T> d6 = t218 + t219 + t220 + s12*t310 + s34*t310 + t402; d6 = T(1)/d6;
complex<T> d9 = spb56*(t218 + t219 + t220 + s12*t310 + s34*t310 + t402); d9 = T(1)/d9;
complex<T> d12 = (s12 - s124)*t360*t97; d12 = T(1)/d12;
complex<T> d17 = spa14*t365; d17 = T(1)/d17;
complex<T> d20 = spb23*t365; d20 = T(1)/d20;
complex<T> d21 = (-s123 + s56)*t191*t97; d21 = T(1)/d21;
complex<T> d24 = (-s124 + s56)*t360*t97; d24 = T(1)/d24;
complex<T> d26 = spa34*t623*square(s134 - s56); d26 = T(1)/d26;
complex<T> d27 = t623*square(s134 - s34); d27 = T(1)/d27;
complex<T> d32 = t596*square(s234 - s34); d32 = T(1)/d32;
complex<T> d34 = spb34*t596*square(s234 - s56); d34 = T(1)/d34;
complex<T> d38 = (-s234 + s34)*t596; d38 = T(1)/d38;
complex<T> d40 = spb23*(t218 + t219 + t220 + s12*t310 + s34*t310 + t402)*t98; d40 = T(1)/d40;
complex<T> d41 = spa14*t365*(t218 + t219 + t220 + s12*t310 + s34*t310 + t402); d41 = T(1)/d41;
complex<T> d43 = spa14*(t218 + t219 + t220 + s12*t310 + s34*t310 + t402)*t98; d43 = T(1)/d43;
complex<T> d44 = spb23*t365*(t218 + t219 + t220 + s12*t310 + s34*t310 + t402); d44 = T(1)/d44;
complex<T> d49 = spa56*t103*(t400 + t401); d49 = T(1)/d49;
complex<T> d50 = t179*t191*T(2); d50 = T(1)/d50;
complex<T> d51 = spa56*t103*t179; d51 = T(1)/d51;
complex<T> d52 = t179*t191; d52 = T(1)/d52;
complex<T> d56 = spa34*spb56*t185*t226; d56 = T(1)/d56;
complex<T> d57 = t179*t360; d57 = T(1)/d57;
complex<T> d58 = (spa12*spb13 + spa24*spb34)*t179*t360; d58 = T(1)/d58;
complex<T> d63 = t179*t360*T(2); d63 = T(1)/d63;
complex<T> d66 = spa56*t103*t365; d66 = T(1)/d66;
complex<T> d67 = (spa12*spb13 + spa24*spb34)*spb56*t365; d67 = T(1)/d67;
complex<T> d70 = spb34*t185*T(2); d70 = T(1)/d70;
complex<T> d74 = spb34*t596*square(s234); d74 = T(1)/d74;
complex<T> d75 = (t218 + t219 + t220 + s12*t310 + s34*t310 + t402)*T(2); d75 = T(1)/d75;
complex<T> d76 = spa14*(spa13*spb23 + spa14*spb24)*(t218 + t219 + t220 + s12*t310 + s34*t310 + t402); d76 = T(1)/d76;
complex<T> d77 = s123*t103; d77 = T(1)/d77;
complex<T> d84 = spb23*(spa13*spb23 + spa14*spb24)*(t218 + t219 + t220 + s12*t310 + s34*t310 + t402); d84 = T(1)/d84;
complex<T> d90 = spa34*spb56*t185*T(2); d90 = T(1)/d90;
complex<T> d91 = spa34*t324*t623; d91 = T(1)/d91;
complex<T> d102 = spa56*t205*t365; d102 = T(1)/d102;
complex<T> d103 = spa56*t179; d103 = T(1)/d103;
complex<T> d104 = spa24*spb56*t185; d104 = T(1)/d104;
complex<T> d105 = (spa12*spb13 + spa24*spb34)*t360*t365; d105 = T(1)/d105;
complex<T> d106 = t103*t226*t360*T(2); d106 = T(1)/d106;
complex<T> d110 = spa14*t179*T(2); d110 = T(1)/d110;
complex<T> d111 = spa14*(spa12*spb13 + spa24*spb34)*t365; d111 = T(1)/d111;
complex<T> d112 = t103*t365; d112 = T(1)/d112;
complex<T> d113 = spb23*t179*T(2); d113 = T(1)/d113;
complex<T> d116 = spa24*spb56*t185*T(2); d116 = T(1)/d116;
complex<T> t39 = t205*t314 + t115*T(2); 
complex<T> t42 = s34*spa24*t314 - spa12*spa34*spb13*t314 - t101*(spa45*t204 + t452*t81)*T(2) - spa24*t206*t446*T(3); 
complex<T> t61 = -(s14*s23*t223) + square(spa12*spa25*spb12 + spa12*spa35*spb13 - spa24*t309); 
complex<T> t74 = -(d51*square(s123)*square(spa45)*square(spb13)) + d52*t205*square(spa15*spb13 + spa25*spb23); 
complex<T> t75 = d56*t222*t224*t243 - t313*t84*square(spa15) + d55*spb24*square(spa12*spb24 + spa13*spb34)*square(t211); 
complex<T> t76 = d59*(-t101 + t406)*t475 - t312*t83*square(spb26) + d61*spa13*square(spa13*spb12 + spa34*spb24)*square(t212); 
complex<T> t85 = d57*t214; 
complex<T> t89 = d6*(spa24*t215 + spa12*spa34*spb13*T(2)); 
complex<T> t90 = d23*s124*spa23*t182 - d24*spa23*t217*t326 + d12*t101*square(spa14*spb16 + spa24*spb26); 
complex<T> t91 = -(d22*s123*spb14*t187) + d21*s123*t216*t447 + d1*t102*square(spa15*spb13 + spa25*spb23); 
complex<T> t107 = d37*t309 + d3*t367*T(2); 
complex<T> t117 = d47*t207*t309 + d46*spa25*t207*t68 + d45*spb16*t309*t69; 
complex<T> t157 = t178*t403 + t221*t404 + t425*t533; 
complex<T> t184 = d6*t215; 
complex<T> t244 = t311 - d3*s56*t556*T(4); 
complex<T> t327 = d31*t106; 
complex<T> t334 = -(d66*s14); 
complex<T> t335 = -(d105*s23); 
complex<T> t339 = d6*spa34; 
complex<T> t341 = d67*t100; 
complex<T> t353 = t213*t253; 
complex<T> t356 = d102*t99; 
complex<T> t369 = d5*t77; 
complex<T> t417 = d4*spb12; 
complex<T> t426 = d39*t290; 
complex<T> t429 = t213*t378; 
complex<T> t453 = t183*t225; 
complex<T> t461 = d36*t309; 
complex<T> t476 = d50*s34; 
complex<T> t484 = t309*t422; 
complex<T> t485 = t207*t423; 
complex<T> t503 = d75*t204; 
complex<T> t551 = d34*spb24; 
complex<T> t570 = d26*spa13; 
complex<T> t585 = spb46*t540; 
complex<T> t634 = t367*t612; 
complex<T> t5 = d14*(spb16*spb34 + spb14*spb36)*t217 + d6*(d15*spa24*spb12*(spa35*spb36 - t187)*t204 + d16*spb13*spb46*t204*(spa35*spb56 - t315) - d8*spa23*t217*t80 + spa35*(spa35*spb13*spb56 - spa24*t80) - t404*(t205*t21 + t115*T(2)) - d16*spb12*spb36*(-(spa34*spb46*t101) + d15*spa24*spb36*t204*t311 + spa23*spb46*t204*T(2)) + spa45*spb12*(spa23*spb46 + d15*spa24*spb56*t21)*T(4)) - spb12*t369*(spa24*t215 + spa12*spa34*spb13*T(2))*T(6); 
complex<T> t6 = d18*(-(spa25*spa34) + spa23*spa45)*t216 + d6*(spb46*(spa12*spa35*spb13 - spa24*spa56*spb46) - d15*spb13*(t182 + t187)*t25 - d19*spa24*spa35*t204*(spa56*spb46 + t314) + d7*spb14*t216*t79 + t403*(spa34*spb46*t214 + t116*T(2)) + d19*spa12*spa45*(-(spa35*spb34*t102) + d15*spa45*spb13*t204*t311 + spa35*spb14*t204*T(2)) - (spa35*spb14 + d15*spa56*spb13*t22)*t371*T(4)) + spa12*t105*t369*T(6); 
complex<T> t8 = d90*spa23*spb16*t207*t243 + d72*spa23*spb16*t243*t309 + d79*spa25*spb14*t207*(t101 - t406) + d73*spa25*spb14*t309*(t101 - t406) + d71*spa23*t410*t68 + d70*spb14*t410*t69 + (d48*t204 + d37*t445)*t75 + (d48*t204 + d30*t445)*t76 + d78*spb14*(-(spa12*spb16) + spa24*spb46)*t79 - d77*spa23*(spa25*spb12 + spa35*spb13)*t80 + t312*(d48*t204 + d30*t445)*t71*t83 + t313*(d48*t204 + d37*t445)*t70*t84 + d68*spb16*t209*t23*T(2) + d69*spa25*t208*t24*T(2) + d84*spb24*(-(spa34*spb24*t101*t213) + spa15*spb12*t101*t22 - spb12*spb14*t212*t23 + spb14*spb24*t217*t23 + spb13*spb24*t23*t315 - s12*t453*t503 - t321*t77 + spa13*(-(spb12*t101*t213) + spb12*t212*t367*t86 + spb24*t217*t367*t86 + spb23*t315*t367*t86) - spa23*spb12*spb24*spb36*t23*T(2) - spa24*spb12*spb24*spb46*t23*T(2)) - d91*spa13*square(t207)*(s134*t204 + t445*T(2)) - d74*spb24*square(spa25*spb24 + t446)*(s234*t204 + t445*T(2)) + d96*(-(d95*spa23*spb16*(spa34*spb24*t178 + spa13*spb12*(spa35*spb56 + t178))) + d80*spa35*spb12*t228*t309 - d93*spb24*t204*t309*t447 - d94*(-(s123*t173) + s234*t21 + s56*(t173 + t21))*t522 + d92*spb24*t410*t68 - spa35*spb14*(t178 + spa34*spb46*T(2))) + d89*(d87*spa25*spb14*(spa12*spb24*(spa56*spb46 - t177) - spa13*spb34*t177) + d88*spa13*spa23*spb36*t204*t207 - d86*((s134 + s56)*spa34*spb46 + (-s124 + s56)*t172)*t545 + d85*spa13*t410*t69 + d82*spa12*spa56*spb46*t207*T(2) + spa23*spb46*(t177 - t446*T(2))) - d76*spa13*(t102*(spa12*spb26*t21 + spa12*spb24*t213 + spa13*spb34*t213) + s12*(-s13 - s14 + s23 + s24)*t183*t503 + t321*t77 + spb24*(spa12*t211 - spa13*t216 - spa14*t314)*t366*t86 - t24*(spa13*(spa23*t216 + spa24*t314) + spa12*(spa23*t211 + spa13*(spa35*spb13 + t447)*T(2)))) - d100*spb12*(d5*spa12*spb34*t183*(spa24*spb12*t204 - spa34*spb13*t206) + d99*spa23*spb46*t216 + d6*(spa24*spb34*t213*t311 - spa23*spb46*(t215*t452 + t206*t499) + d15*(s13 - s14 + s23 - s24)*spa45*spb13*t371*(t68 + t69) + spb14*(spb36*t204 - spb46*t227 + t228*t499)*t79 - spa12*spb13*t187*t311*T(3))) + d98*spa12*(d97*spa35*spb14*t217 + d5*spa34*spb12*t183*(spa24*spb34*t206 - spb13*t25) + d6*(spa35*spb14*(spa34*spb36*t206 + spa45*spb56*t215) + d15*(s13 - s14 + s23 - s24)*spa24*spa45*spb12*spb36*(t68 + t69) + spa23*t80*(spa45*t204 - spa35*t227 + t452*t81) + t311*(-(spa34*spb13*t213) + spa24*spa35*spb12*spb36*T(3)))) + d83*(s34*spa35*spb36*t101 - s24*spa23*spb36*t314 - spa24*spb46*t211*(t101 + t406) + t187*(-(s24*t101) + s12*t406) + s14*(spa23*spb36 + spa24*spb46)*t446 + (spa25*spb23*spb46 - spb36*t177)*(s13*spa23 + spa24*t68) - spb36*(s24*spa25 + spa24*t221)*(t68 + t69) + spa25*t406*(s14*spb26 + s24*spb26 + spb26*t470 + spa14*t522)*T(2) + spb14*t214*(spa25*t366 + spa15*(t172 + t22)*T(2)) + d82*spa34*spb23*t101*(spa12*spa25*spb26 - spa13*spa25*spb36 - spa35*t371*T(2)) + s23*spa15*spb16*(t406 - t101*T(3)) - spa12*spa24*t24*t446*T(4)) + d81*((-(s12*spa34*spb14) + s13*t102)*t182 - spa35*spb13*(spa34*spb14 + t102)*t212 + s13*t315*t447 + s23*spa34*spb46*(spa35*spb13 + t447) + spa45*(s13*spb16 + spa23*spb13*spb26)*(t68 + t69) + (-(spa14*spa35*spb16) + spa45*t178)*(s24*spb14 + spb13*t69) + s34*spa23*spa45*t80 + (-(spa15*(s13 + s23 + t470)) - spb23*t545)*t579*t81 - spa23*t205*(spb16*t367 + spb26*(t173 + t446)*T(2)) + d80*spa14*spb34*t102*(spa15*spb12*spb16 + spa45*spb16*spb24 - spa45*t80*T(2)) + s14*spa25*spb26*(spa34*spb14 - t102*T(3)) + spa34*spb13*t23*t80*T(4)); 
complex<T> t26 = (-s123 + s23)*(d52*t205*t72 - t74); 
complex<T> t28 = s234*(d53*s12*s234*t313*t70 + t321*t75); 
complex<T> t29 = (s156 - s34)*(d104*t222*t224 + t75 + t313*t70*t84); 
complex<T> t34 = -((s256 - s34)*(d65*t475 - t76 - t312*t71*t83)); 
complex<T> t78 = -(t85*square(spa14*spb16 + spa24*spb26)) - d58*square(s124)*square(spa24)*square(spb36); 
complex<T> t147 = t178*t403 + t221*t404 + t426*t567; 
complex<T> t328 = d38*t107; 
complex<T> t372 = (-s13 - s14 + s23 + s24)*t184; 
complex<T> t428 = t100*t335; 
complex<T> t433 = t334*t99; 
complex<T> t441 = s12*t312*t324*t432 + d48*t44*t76; 
complex<T> t454 = d48*t74; 
complex<T> t464 = t105*t339; 
complex<T> t543 = t77*t89; 
complex<T> t552 = t180*t369; 
complex<T> t633 = spa23*t356; 
complex<T> t636 = spb14*t341; 
complex<T> t642 = t484*t492; 
complex<T> t643 = t366*t585; 
complex<T> t649 = t485*t491; 
complex<T> t662 = t461*t620; 
complex<T> t11 = -(d6*s12*s234*t188) + t178*t403 + t198*t556 + t204*t417*(spa12*spa35*t173 - t173*t23 - spa15*spa25*t311 - t556*t70*t86) + d9*t25*(t24*t315 - t224*t556 + spa34*spb46*t80 + t178*t80 + t556*t71*t86) - spa12*spa35*t184*(t24 + t522)*T(2) + t225*t552*T(6); 
complex<T> t12 = t221*t404 - d6*t188*t44 + t198*t556 + t204*t417*(t23*t314 - t223*t556 + t177*t79 + t21*t79 + t556*t70*t86) - d9*t25*(spb16*spb26*t311 + t172*(t24 + t80) + t556*t71*t86) + t184*(t23 + t545)*t80*T(2) + (-s13 - s14 + s23 + s24)*t552*T(6); 
complex<T> t30 = -((s12 - s124)*(t78 + t73*t85)); 
complex<T> t31 = d52*(d99*s12*s56 + d48*t206)*t323 + d99*s12*s56*t74 + d48*t206*t74 + d97*s12*s56*t78 + d48*t206*t78 + d97*s12*s56*t73*t85 + d48*t206*t73*t85; 
complex<T> t32 = t428 - d104*t224*t625 + t633 + d103*spa23*t205*t72 + s23*t74; 
complex<T> t41 = -t454; 
complex<T> t47 = s34*t78; 
complex<T> t51 = t212*t327*t491 + d27*t208*t212*t579 - s134*t222*t570*t71 + s134*t208*t580*t71 - spa23*t346*t366*T(2) - t649*T(2); 
complex<T> t53 = -(t147*t206) + t184*t453 - t213*(t492 + t68)*T(2) + t302*T(4) + t252*(t311 - d3*s56*t556*T(4)) - t634*T(8); 
complex<T> t57 = -(t157*t206) + t183*t372 - t213*(t491 + t69)*T(2) + t300*T(4) + (spa25*spb12*t477 + t322*t478)*(t311 - d3*s56*t556*T(4)) + t643*T(8); 
complex<T> t60 = -(t316*t408*t551) - d32*t209*t211*t587 + t209*t408*t590 + t211*t328*t69 + t642*T(2) - t662*T(2); 
complex<T> t92 = t433 + d65*s14*t475 + t636; 
complex<T> t164 = s123*(s34*t454 + t323*t476); 
complex<T> t531 = t464*t77; 
complex<T> t611 = t539*t543; 
complex<T> t2 = d10*spb24*t11 + d11*spa13*t12 + d13*spb36*t207*t311 - d2*spa45*t309*t311 - d17*spa12*t5 + d20*spb12*t6 - d1*t102*t72 - d12*t101*t73; 
complex<T> t33 = d63*s34*t332 - d107*spa34*t313*t408 + s34*t454 + d48*t47 + t323*t476 - d108*s134*spb34*t312*t71 + d48*s34*t75 + d48*s34*t76; 
complex<T> t36 = -(d111*spa56*t100) + d113*spb56*t323 + d110*spa56*t332 + d114*spb56*t313*t408 + s56*t41 + d109*s134*spa56*t312*t71 - d48*s56*t75 - d48*s56*t76 - d48*s56*t78 - d112*spb56*t99; 
complex<T> t95 = d48*s124*t47 + d63*s124*s34*t214*square(spa14*spb16 + spa24*spb26); 
complex<T> t396 = -(d64*s12*t100) + t312*t432*t44 + t321*t75 + t321*t76 + t321*t78 + s12*(t117 + d63*t332 + d53*t313*t408 + t41 + d50*t205*t72 - d49*t99); 
complex<T> t622 = spa12*t531; 
complex<T> t1 = -(d10*spb24*t11) - d11*spa13*t12 + d22*s123*spb14*t187 - d13*spb36*t207*t311 + d2*spa45*t309*t311 + d23*spa23*spa35*t326 + d24*spa23*t217*t326 - d21*s123*t216*t447 + d17*spa12*t5 - d40*spb24*t53 + spa23*t346*t533 + t316*t408*t551 - d43*spa13*t57 - t209*t408*t590 - d20*spb12*t6 + s134*t222*t570*t71 - s134*t208*t580*t71 + t662*T(2) - d48*(d47*t207*t309 + d46*spa25*t207*t68 + d45*spb16*t309*t69)*T(3) + d44*(d7*t21*t39 - t403*t43 + d15*spb13*t229*(-(s13*t206) + s14*t206 - s23*t206 + s24*t206 + t218 + t219 + t220 - s56*t470 - t180*T(2) - t445*T(2)) + t611*T(6)) + d41*(-(d15*spa24*t113*t229) + d8*t22*t40 + t404*t42 + t622*T(6)); 
complex<T> t7 = d40*spb24*t53 + d43*spa13*t57 - d27*t208*t212*t579 + d32*t209*t211*t587 + t212*t327*t68 - t642*T(2) + t649*T(2) + d38*t211*t492*(d37*t309 + d3*t367*T(2)) + d44*(spb34*t39*t404 + t403*t43 + d15*spb13*t229*((s13 - s14 + s23 - s24)*t206 - t218 - t219 - t220 + s56*t470 + t180*T(2) + t445*T(2)) - t611*T(6)) + d41*(d15*spa24*t113*t229 + spa34*t40*t403 - t404*t42 - t622*T(6)); 
complex<T> co1 = -(d101*spb14*t61); 
complex<T> co2 = -(d106*spa23*t52); 
complex<T> co3 = s12*t633; 
complex<T> co4 = s124*t428; 
complex<T> co5 = d115*s14*s34*t475; 
complex<T> co6 = s12*t636; 
complex<T> co7 = s123*t433; 
complex<T> co8 = -(d116*s34*t224*t625); 
complex<T> co9 = Complex(0,1); 
SeriesC<T> result = co9*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t91*(*CI_users[1]->get_value(mc,ind,mu)) + t90*(*CI_users[2]->get_value(mc,ind,mu)) + t51*(*CI_users[3]->get_value(mc,ind,mu)) + t60*(*CI_users[4]->get_value(mc,ind,mu)) + t7*(*CI_users[5]->get_value(mc,ind,mu)) + t1*(*CI_users[6]->get_value(mc,ind,mu)) + t396*(*CI_users[7]->get_value(mc,ind,mu)) + t26*(*CI_users[8]->get_value(mc,ind,mu)) + t34*(*CI_users[9]->get_value(mc,ind,mu)) + t92*(*CI_users[10]->get_value(mc,ind,mu)) + t8*(*CI_users[11]->get_value(mc,ind,mu)) + co1*(*CI_users[12]->get_value(mc,ind,mu)) + t32*(*CI_users[13]->get_value(mc,ind,mu)) + t29*(*CI_users[14]->get_value(mc,ind,mu)) + co2*(*CI_users[15]->get_value(mc,ind,mu)) + t33*(*CI_users[16]->get_value(mc,ind,mu)) + t31*(*CI_users[17]->get_value(mc,ind,mu)) + t30*(*CI_users[18]->get_value(mc,ind,mu)) + t36*(*CI_users[19]->get_value(mc,ind,mu)) + t28*(*CI_users[20]->get_value(mc,ind,mu)) + t441*(*CI_users[21]->get_value(mc,ind,mu)) + co3*(*CI_users[22]->get_value(mc,ind,mu)) + co4*(*CI_users[23]->get_value(mc,ind,mu)) + co5*(*CI_users[24]->get_value(mc,ind,mu)) + t95*(*CI_users[25]->get_value(mc,ind,mu)) + co6*(*CI_users[26]->get_value(mc,ind,mu)) + co7*(*CI_users[27]->get_value(mc,ind,mu)) + t164*(*CI_users[28]->get_value(mc,ind,mu)) + co8*(*CI_users[29]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmmmemep_SLC_wCI::\
C2q2g2l_qpqmmmemep_SLC_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c134, c256));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c2356));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c256));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c2356));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c3, c14, c2, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c356));
CI_users.push_back(new Cached_Box_Integral_User(c4, c12, c3, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmmmemep_SLC_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmmmemep SLC");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:40:09 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb14 = SPB(1,4);
complex<T> spb16 = SPB(1,6);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb26 = SPB(2,6);
complex<T> spb36 = SPB(3,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb46 = SPB(4,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spa14 = SPA(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s124 = SS(1,2,4);
complex<T> s134 = SS(1,3,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s156 = SS(1,5,6);
complex<T> s234 = SS(2,3,4);
complex<T> s256 = SS(2,5,6);
complex<T> t4 = spb23*spb34; 
complex<T> t5 = spb56*T(2); 
complex<T> t24 = square(spb23); 
complex<T> t26 = square(spb26); 
complex<T> t27 = square(spb13); 
complex<T> t28 = spb14*spb56; 
complex<T> t29 = square(spb36); 
complex<T> t31 = square(spb12); 
complex<T> t37 = spa14*(spb14*spb26 + spb12*spb46); 
complex<T> t41 = square(spb16); 
complex<T> t42 = cube(spb23); 
complex<T> t43 = -(spb13*spb26) - spb12*spb36; 
complex<T> t44 = -(spa23*spb13) - spa24*spb14; 
complex<T> t46 = spb14*spb36 + spb13*spb46; 
complex<T> t50 = -s256 + s34; 
complex<T> t54 = square(spa25); 
complex<T> t56 = square(spa35); 
complex<T> t57 = square(spa45); 
complex<T> t58 = spa23*spb12 - spa34*spb14; 
complex<T> t60 = -(spa13*spb16) - spa23*spb26; 
complex<T> t61 = -(spb13*spb26) + spb16*spb23*T(2); 
complex<T> t62 = spb16*spb23 + spb13*spb26; 
complex<T> t65 = (spa14*spb16 + spa24*spb26)*spb46; 
complex<T> t68 = -s156 + s34; 
complex<T> t69 = s123*s234 - s23*s56; 
complex<T> t94 = spa34*spb36; 
complex<T> t103 = spb16*spb26; 
complex<T> d2 = (s12 - s123)*spb23*spb56*square(spb34); d2 = T(1)/d2;
complex<T> d3 = (s12 - s124)*spb24*spb56*square(spb34); d3 = T(1)/d3;
complex<T> d4 = (-s123 + s56)*spb23*square(spb34); d4 = T(1)/d4;
complex<T> d7 = spb23*square(spb34); d7 = T(1)/d7;
complex<T> d9 = -s124 + s56; d9 = T(1)/d9;
complex<T> d14 = -s134 + s56; d14 = T(1)/d14;
complex<T> d23 = spb23*spb56*cube(spb34); d23 = T(1)/d23;
complex<T> d28 = spb56*cube(spb34); d28 = T(1)/d28;
complex<T> d35 = spb23*cube(spb34); d35 = T(1)/d35;
complex<T> d36 = spb14*spb23*square(spb34); d36 = T(1)/d36;
complex<T> t30 = -(spb16*t43); 
complex<T> t77 = -(spb14*t29); 
complex<T> t92 = -(t26*t27); 
complex<T> d1 = spb24*spb56*t4; d1 = T(1)/d1;
complex<T> d5 = t4*square(s123 - s56); d5 = T(1)/d5;
complex<T> d6 = spb14*spb34*t24; d6 = T(1)/d6;
complex<T> d8 = t28*t4; d8 = T(1)/d8;
complex<T> d10 = (-s124 + s14)*spb24*t24*t28; d10 = T(1)/d10;
complex<T> d11 = spb14*t4*square(s124 - s56); d11 = T(1)/d11;
complex<T> d12 = spb14*t4; d12 = T(1)/d12;
complex<T> d13 = spb34*t24*t28; d13 = T(1)/d13;
complex<T> d15 = spb14*t4*t5; d15 = T(1)/d15;
complex<T> d16 = (-s134 + s56)*t28*t4; d16 = T(1)/d16;
complex<T> d17 = (-s134 + s14)*spb56*t4; d17 = T(1)/d17;
complex<T> d18 = (-s134 + s14)*spb34*spb56*t24; d18 = T(1)/d18;
complex<T> d19 = t4*t5*square(s134 - s14); d19 = T(1)/d19;
complex<T> d20 = spb14*t4*square(s134 - s56)*T(2); d20 = T(1)/d20;
complex<T> d21 = spb24*t24*t28; d21 = T(1)/d21;
complex<T> d22 = spb34*t28*t42; d22 = T(1)/d22;
complex<T> d24 = spb23*t28*square(spb34); d24 = T(1)/d24;
complex<T> d25 = spb34*spb56*t42; d25 = T(1)/d25;
complex<T> d26 = spb56*t4*square(spb24); d26 = T(1)/d26;
complex<T> d27 = spb34*spb56*t24; d27 = T(1)/d27;
complex<T> d29 = t28*square(spb34); d29 = T(1)/d29;
complex<T> d30 = spb23*t28; d30 = T(1)/d30;
complex<T> d31 = t28*t42; d31 = T(1)/d31;
complex<T> d32 = t24*t28; d32 = T(1)/d32;
complex<T> d33 = spb24*t28*t4; d33 = T(1)/d33;
complex<T> d34 = spb14*spb34*t42; d34 = T(1)/d34;
complex<T> d37 = spb14*spb34*t5; d37 = T(1)/d37;
complex<T> d38 = spb14*spb34*t42*t5; d38 = T(1)/d38;
complex<T> d39 = spb14*spb34*t24*t5; d39 = T(1)/d39;
complex<T> d40 = t4*t5*square(spb24); d40 = T(1)/d40;
complex<T> d41 = spb24*t4*t5; d41 = T(1)/d41;
complex<T> d42 = spb23*t5*cube(spb34); d42 = T(1)/d42;
complex<T> d43 = spb14*spb23*t5*square(spb34); d43 = T(1)/d43;
complex<T> t3 = d12*d14*spa25*spb12*spb16 + d7*d9*spa35*spb13*spb36 - d4*spa45*spb13*spb46 - d1*t103 - d21*spb12*t103 + d13*d14*spb13*t26*t44 + d20*spb56*t31*t54 + d11*spb56*t27*t56 + d5*t28*t57 - d8*d9*spb16*spb36*t58 + d6*d9*spa35*spb13*t61 - d16*t103*t44*T(2) + d8*t41*T(3); 
complex<T> t12 = s12*(d22*t26*t27 + d13*spb16*t43); 
complex<T> t13 = (s124*s134 - s14*s56)*(d38*t26*t27 + d39*spb16*t43); 
complex<T> t14 = spa34*(d32*t30 - d30*t41 + d31*t92); 
complex<T> t15 = d26*s14*spb14*t26 + d25*spa14*t26*t27 + d1*spb16*t37 + d27*spa14*spb16*t43; 
complex<T> t18 = (s123*s124 - s12*s56)*(d42*spb14*t29 - d43*spb16*t46); 
complex<T> t19 = spa23*(d29*spb16*t46 + d28*t77); 
complex<T> t20 = spa56*(d34*t26*t27 + d35*spb14*t29 + d12*t41 + d6*spb16*t43 - d36*spb16*t46); 
complex<T> t21 = d21*spb12*t103 + d10*spa24*spb26*spb46*t31 - d10*spa12*t103*t31 - d18*t62*t94 + d19*spb14*t29*square(spa34) + d17*spb16*t94*T(2); 
complex<T> t22 = d1*t103 + d2*spb13*spb36*t60 + d3*spb14*t65; 
complex<T> t23 = d4*spa45*spb13*spb46 - d5*t28*t57 - d2*spb13*spb36*t60; 
complex<T> t35 = -(d6*t61); 
complex<T> t38 = d24*t46; 
complex<T> t73 = s12*(d40*s14*spb14*t26 + d41*spb16*t37); 
complex<T> t91 = t50*(d13*t30 + d22*t92); 
complex<T> t111 = -(d8*t41); 
complex<T> t140 = d15*t41; 
complex<T> t1 = -(d12*d14*spa25*spb12*spb16) - t140 - d13*d14*spb13*t26*t44 - d20*spb56*t31*t54 + d18*t62*t94 + d19*t77*square(spa34) + d16*t103*t44*T(2) - d17*spb16*t94*T(2); 
complex<T> t2 = -(d7*d9*spa35*spb13*spb36) + t111 - d10*spa24*spb26*spb46*t31 + d10*spa12*t103*t31 + d9*spa35*spb13*t35 - d11*spb56*t27*t56 + d8*d9*spb16*spb36*t58 - d3*spb14*t65; 
complex<T> t16 = (s12 - s124)*(-(d33*spb16*(spb14*spb26 + spb12*spb46)) + d26*spb14*t26 + d13*t30 + spb16*t38 + d23*t77 + d22*t92); 
complex<T> t17 = (-s123 + s23)*(spb16*t38 + d23*t77); 
complex<T> co1 = t111*t68; 
complex<T> co2 = -(d37*s12*spa23*t41); 
complex<T> co3 = t140*t69; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t22*(*CI_users[0]->get_value(mc,ind,mu)) + t23*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t21*(*CI_users[4]->get_value(mc,ind,mu)) + t3*(*CI_users[5]->get_value(mc,ind,mu)) + t12*(*CI_users[6]->get_value(mc,ind,mu)) + t17*(*CI_users[7]->get_value(mc,ind,mu)) + t91*(*CI_users[8]->get_value(mc,ind,mu)) + t15*(*CI_users[9]->get_value(mc,ind,mu)) + t19*(*CI_users[10]->get_value(mc,ind,mu)) + co1*(*CI_users[11]->get_value(mc,ind,mu)) + t14*(*CI_users[12]->get_value(mc,ind,mu)) + t16*(*CI_users[13]->get_value(mc,ind,mu)) + t20*(*CI_users[14]->get_value(mc,ind,mu)) + co2*(*CI_users[15]->get_value(mc,ind,mu)) + co3*(*CI_users[16]->get_value(mc,ind,mu)) + t13*(*CI_users[17]->get_value(mc,ind,mu)) + t73*(*CI_users[18]->get_value(mc,ind,mu)) + t18*(*CI_users[19]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g2l_qpppqmemep_nf_wCI::\
C2q2g2l_qpppqmemep_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g2l_qpppqmemep_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpppqmemep nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g2l_qppmqmemep_nf_wCI::\
C2q2g2l_qppmqmemep_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g2l_qppmqmemep_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qppmqmemep nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g2l_qpmpqmemep_nf_wCI::\
C2q2g2l_qpmpqmemep_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g2l_qpmpqmemep_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpmpqmemep nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g2l_qpmmqmemep_nf_wCI::\
C2q2g2l_qpmmqmemep_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g2l_qpmmqmemep_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpmmqmemep nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g2l_qpppqmemep_nf_top_wCI::\
C2q2g2l_qpppqmemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g2l_qpppqmemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpppqmemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g2l_qppmqmemep_nf_top_wCI::\
C2q2g2l_qppmqmemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g2l_qppmqmemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qppmqmemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g2l_qpmpqmemep_nf_top_wCI::\
C2q2g2l_qpmpqmemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g2l_qpmpqmemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpmpqmemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g2l_qpmmqmemep_nf_top_wCI::\
C2q2g2l_qpmmqmemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g2l_qpmmqmemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpmmqmemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g2l_qpqmppemep_VECT_wCI::\
C2q2g2l_qpqmppemep_VECT_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c356, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c456, c123));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c12, c3, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmppemep_VECT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, p, em, ep}, VECT}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmppemep VECT");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:40:21 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa56 = SPA(5,6);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb56 = SPB(5,6);
complex<T> spb46 = SPB(4,6);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb36 = SPB(3,6);
complex<T> s23 = S(2,3);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s12 = S(1,2);
complex<T> s124 = SS(1,2,4);
complex<T> t1 = square(spa34); 
complex<T> t13 = spa24*spa35; 
complex<T> t14 = spa23*spa45; 
complex<T> t15 = square(spa25); 
complex<T> t22 = s12 - s124; 
complex<T> t23 = -s124 + s56; 
complex<T> t29 = spa12*spa56; 
complex<T> t35 = spa24*spa35 + spa23*spa45; 
complex<T> t38 = spa23*spb36; 
complex<T> t39 = spa24*spb46; 
complex<T> t40 = square(spa35); 
complex<T> t41 = square(spb13); 
complex<T> t48 = square(spa23); 
complex<T> t49 = square(spa24); 
complex<T> t50 = square(spa45); 
complex<T> t51 = square(spb14); 
complex<T> t52 = square(spb36); 
complex<T> t53 = square(spb46); 
complex<T> d2 = (s12 - s123)*spa56*cube(spa34); d2 = T(1)/d2;
complex<T> d8 = (-s123 + s56)*spa12*cube(spa34); d8 = T(1)/d8;
complex<T> d17 = spa12*square(square(spa34)); d17 = T(1)/d17;
complex<T> t57 = t13*t14; 
complex<T> d1 = (s12 - s123)*spa56*t1; d1 = T(1)/d1;
complex<T> d3 = spa56*t1*square(s12 - s123); d3 = T(1)/d3;
complex<T> d4 = spa56*t1*t22; d4 = T(1)/d4;
complex<T> d5 = spa56*t22*cube(spa34); d5 = T(1)/d5;
complex<T> d6 = spa56*t1*square(t22); d6 = T(1)/d6;
complex<T> d7 = t29*cube(spa34); d7 = T(1)/d7;
complex<T> d9 = spa12*t23*cube(spa34); d9 = T(1)/d9;
complex<T> d10 = spa12*t1*t23; d10 = T(1)/d10;
complex<T> d11 = spa12*t1*square(t23); d11 = T(1)/d11;
complex<T> d12 = (-s123 + s56)*spa12*t1; d12 = T(1)/d12;
complex<T> d13 = spa12*t1*square(s123 - s56); d13 = T(1)/d13;
complex<T> d14 = t1*t29; d14 = T(1)/d14;
complex<T> d15 = t29*square(square(spa34)); d15 = T(1)/d15;
complex<T> d16 = spa12*t1; d16 = T(1)/d16;
complex<T> d18 = t29*square(square(spa34))*T(2); d18 = T(1)/d18;
complex<T> d19 = t1*t29*T(2); d19 = T(1)/d19;
complex<T> t8 = -(d1*spa25*spa35*spb13) - d4*spa25*spa45*spb14 + d2*spa35*spb13*t35 - d5*spa45*spb14*t35 + d3*spa12*t40*t41 + d6*spa12*t50*t51; 
complex<T> t9 = d4*spa25*spa45*spb14 - d7*spa25*t35 + d5*spa45*spb14*t35 + d9*t35*t38 - d6*spa12*t50*t51; 
complex<T> t10 = d1*spa25*spa35*spb13 + d7*spa25*t35 - d2*spa35*spb13*t35 - d8*t35*t39 - d3*spa12*t40*t41; 
complex<T> t11 = -(d10*spa25*t38) - d9*t35*t38 - d12*spa25*t39 + d8*t35*t39 + d11*spa56*t48*t52 + d13*spa56*t49*t53; 
complex<T> t16 = d18*(s123*s124 - s12*s56); 
complex<T> t19 = s23*(d14*square(spa25) + d15*spa23*spa24*spa35*spa45*T(2)); 
complex<T> t21 = spb56*(d16*square(spa25) + d17*spa23*spa24*spa35*spa45*T(2)); 
complex<T> t42 = -(d15*T(2)); 
complex<T> t46 = d10*spa25*t38 - d11*spa56*t48*t52; 
complex<T> t47 = d12*spa25*t39 - d13*spa56*t49*t53; 
complex<T> t6 = (s123 - s23)*(d14*t15 - t42*t57); 
complex<T> t56 = -(d14*t15*t22) + t22*t42*t57; 
complex<T> t80 = t16*t57; 
complex<T> t7 = d19*(s123*s124 - s12*s56)*t15 + t80; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t8*(*CI_users[0]->get_value(mc,ind,mu)) + t10*(*CI_users[1]->get_value(mc,ind,mu)) + t9*(*CI_users[2]->get_value(mc,ind,mu)) + t46*(*CI_users[3]->get_value(mc,ind,mu)) + t47*(*CI_users[4]->get_value(mc,ind,mu)) + t11*(*CI_users[5]->get_value(mc,ind,mu)) + t6*(*CI_users[6]->get_value(mc,ind,mu)) + t19*(*CI_users[7]->get_value(mc,ind,mu)) + t56*(*CI_users[8]->get_value(mc,ind,mu)) + t21*(*CI_users[9]->get_value(mc,ind,mu)) + t80*(*CI_users[10]->get_value(mc,ind,mu)) + t7*(*CI_users[11]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmpmemep_VECT_wCI::\
C2q2g2l_qpqmpmemep_VECT_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c356, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c456, c123));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c56, c34, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmpmemep_VECT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, m, em, ep}, VECT}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmpmemep VECT");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:45:43 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa56 = SPA(5,6);
complex<T> spb46 = SPB(4,6);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa36 = SPA(3,6);
complex<T> spa25 = SPA(2,5);
complex<T> spa26 = SPA(2,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spa13 = SPA(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb16 = SPB(1,6);
complex<T> spb36 = SPB(3,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa46 = SPA(4,6);
complex<T> s56 = -(spa56*spb56);
complex<T> s123 = SS(1,2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> s124 = SS(1,2,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s24 = -(spa24*spb24);
complex<T> s35 = -(spa35*spb35);
complex<T> s36 = -(spa36*spb36);
complex<T> s45 = -(spa45*spb45);
complex<T> s46 = -(spa46*spb46);
complex<T> t12 = -(spa23*(spa13*spb16 + spa23*spb26)); 
complex<T> t13 = -(spb46*(spa25*spb45 + spa26*spb46)); 
complex<T> t15 = -(spb14*(spa15*spb14 + spa25*spb24)); 
complex<T> t34 = spa23*T(2); 
complex<T> t35 = spa35*spb14; 
complex<T> t57 = spa45*spb14; 
complex<T> t58 = spa24*spb46; 
complex<T> t59 = spa35*spb13; 
complex<T> t60 = spa23*spb12; 
complex<T> t63 = (spa23*spb13 + spa24*spb14)*(spa35*spb36 + spa45*spb46); 
complex<T> t91 = spa13*spb14 + spa23*spb24; 
complex<T> t94 = -s12 - s34 + s56; 
complex<T> t95 = -(spa35*spb15) - spa36*spb16; 
complex<T> t97 = -s35 - s36 + s45 + s46; 
complex<T> t98 = s12 - s34 - s56; 
complex<T> t99 = square(spa23); 
complex<T> t100 = square(spa24); 
complex<T> t101 = square(spa35); 
complex<T> t102 = square(spa45); 
complex<T> t103 = square(spb13); 
complex<T> t104 = square(spb14); 
complex<T> t107 = square(spb36); 
complex<T> t109 = square(spb46); 
complex<T> t110 = s13 - s14 + s23 - s24; 
complex<T> t113 = -(spa45*spb35) - spa46*spb36; 
complex<T> t146 = -(s12*T(2)); 
complex<T> t150 = spa25*spb12; 
complex<T> t151 = spa12*spb16; 
complex<T> t154 = s123*spb46; 
complex<T> t168 = (s35 + s36 - s45 - s46)*T(2); 
complex<T> t181 = -(s34*T(2)); 
complex<T> t215 = spa24*spa56; 
complex<T> t218 = spa25*spb16; 
complex<T> d4 = spa56; d4 = T(1)/d4;
complex<T> d8 = spb56; d8 = T(1)/d8;
complex<T> d19 = spa12; d19 = T(1)/d19;
complex<T> d21 = spb12; d21 = T(1)/d21;
complex<T> d23 = spb12*square(s124 - s56)*square(spa35*spb45 + spa36*spb46); d23 = T(1)/d23;
complex<T> d24 = (-s124 + s56)*spb12*square(spa35*spb45 + spa36*spb46); d24 = T(1)/d24;
complex<T> d25 = (-s124 + s56)*spb12*cube(spa35*spb45 + spa36*spb46); d25 = T(1)/d25;
complex<T> d26 = spa12*square(s123 - s56)*square(spa35*spb45 + spa36*spb46); d26 = T(1)/d26;
complex<T> d27 = (-s123 + s56)*spa12*square(spa35*spb45 + spa36*spb46); d27 = T(1)/d27;
complex<T> d28 = (-s123 + s56)*spa12*cube(spa35*spb45 + spa36*spb46); d28 = T(1)/d28;
complex<T> d40 = T(2); d40 = T(1)/d40;
complex<T> d41 = s124; d41 = T(1)/d41;
complex<T> d42 = s123; d42 = T(1)/d42;
complex<T> d47 = s124*(spa35*spb45 + spa36*spb46)*T(2); d47 = T(1)/d47;
complex<T> d48 = s123*(spa35*spb45 + spa36*spb46)*T(2); d48 = T(1)/d48;
complex<T> t8 = -t151; 
complex<T> t9 = -t150; 
complex<T> t10 = square(-(spa23*spb36) + t151); 
complex<T> t11 = square(t150 + t57); 
complex<T> t14 = spa35*t95; 
complex<T> t16 = t110*T(2); 
complex<T> t29 = square(t91); 
complex<T> t37 = -(d8*(s124*spb16 + spb56*t150)*t98) + (-s13 + s14 - s23 + s24)*t57*T(2); 
complex<T> t51 = -(d48*spa45*spb13*(spa25*spb56 + t58)) + d47*spa24*spb36*(spa56*spb16 - t59); 
complex<T> t65 = d4*spb12*square(spa25) + d8*spa12*square(spb16) - t218*T(2); 
complex<T> t76 = -(d19*t94); 
complex<T> t84 = -(d25*(spa35*spb34 + spa56*spb46)); 
complex<T> t85 = d28*(spa34*spb46 + spa35*spb56); 
complex<T> t112 = s123*spa25 + spa56*t151; 
complex<T> t137 = d24*spb56; 
complex<T> t155 = s124*t35; 
complex<T> t156 = spa56*t100; 
complex<T> t160 = t113*t63; 
complex<T> t185 = t101*t103; 
complex<T> t186 = t102*t104; 
complex<T> t187 = t12*t13; 
complex<T> t198 = t107*t99; 
complex<T> t232 = spb46*t60; 
complex<T> t246 = t13*t215; 
complex<T> t261 = t58*t97; 
complex<T> d3 = (s12 - s124)*spa56*cube(t91); d3 = T(1)/d3;
complex<T> d6 = cube(t91)*(s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34)); d6 = T(1)/d6;
complex<T> d7 = t91*square(s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34)); d7 = T(1)/d7;
complex<T> d9 = t91*(s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34)); d9 = T(1)/d9;
complex<T> d10 = (s12 - s123)*spb56*cube(t91); d10 = T(1)/d10;
complex<T> d13 = spa12*spb56*cube(t91); d13 = T(1)/d13;
complex<T> d15 = spa56*spb12*cube(t91); d15 = T(1)/d15;
complex<T> d17 = (spa35*spb45 + spa36*spb46)*square(s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34)); d17 = T(1)/d17;
complex<T> d18 = cube(spa35*spb45 + spa36*spb46)*(s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34)); d18 = T(1)/d18;
complex<T> d20 = (s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34))*square(spa35*spb45 + spa36*spb46); d20 = T(1)/d20;
complex<T> d22 = (spa35*spb45 + spa36*spb46)*(s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34)); d22 = T(1)/d22;
complex<T> d30 = spa56*square(square(t91)); d30 = T(1)/d30;
complex<T> d32 = spb56*square(square(t91)); d32 = T(1)/d32;
complex<T> d33 = spa12*spb56*square(square(t91)); d33 = T(1)/d33;
complex<T> d34 = s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34); d34 = T(1)/d34;
complex<T> d35 = s124*t91*T(2); d35 = T(1)/d35;
complex<T> d36 = s123*t91*T(2); d36 = T(1)/d36;
complex<T> d38 = spa56*spb12*square(square(t91)); d38 = T(1)/d38;
complex<T> d44 = spa12*square(square(t91)); d44 = T(1)/d44;
complex<T> d46 = spb12*square(square(t91)); d46 = T(1)/d46;
complex<T> t19 = -(d9*(spa14*spb13 + spa24*spb23)); 
complex<T> t28 = d33*s123; 
complex<T> t45 = t168*t59 - (s124*spa25 + spa56*t151)*t76; 
complex<T> t52 = -(d21*(s123*spb16 + spb56*t150)*t94) + t261*T(2); 
complex<T> t67 = d6*(-s13 + s14 - s23 + s24); 
complex<T> t71 = d18*(s35 + s36 - s45 - s46); 
complex<T> t72 = d17*t98; 
complex<T> t78 = d3*(spa23*spb12 - spa34*spb14); 
complex<T> t79 = d4*t112; 
complex<T> t80 = -(d7*(spa14*spb13 + spa24*spb23)); 
complex<T> t81 = -(d38*(spa15*spb14 + spa25*spb24)); 
complex<T> t82 = d10*(-(spa12*spb14) + spa23*spb34); 
complex<T> t125 = -(s124*spa25) + spa56*t8; 
complex<T> t126 = -(s123*spb16) + spb56*t9; 
complex<T> t159 = d6*t112; 
complex<T> t163 = spb56*t84; 
complex<T> t165 = d18*t97; 
complex<T> t189 = d38*s124; 
complex<T> t200 = t109*t156; 
complex<T> t254 = t246*t85; 
complex<T> t256 = t137*t185; 
complex<T> d1 = spa56*t29*square(s12 - s124); d1 = T(1)/d1;
complex<T> d2 = (s12 - s124)*spa56*t29; d2 = T(1)/d2;
complex<T> d5 = t29*(s34*t146 + s56*(s56 + t146 + t181) + square(s12) + square(s34)); d5 = T(1)/d5;
complex<T> d11 = spb56*t29*square(s12 - s123); d11 = T(1)/d11;
complex<T> d12 = (s12 - s123)*spb56*t29; d12 = T(1)/d12;
complex<T> d14 = spa12*spb56*t29; d14 = T(1)/d14;
complex<T> d16 = spa56*spb12*t29; d16 = T(1)/d16;
complex<T> d29 = spa56*t29*T(2); d29 = T(1)/d29;
complex<T> d31 = spb56*t29*T(2); d31 = T(1)/d31;
complex<T> d37 = spa56*spb12*t29*T(2); d37 = T(1)/d37;
complex<T> d39 = spa12*spb56*t29*T(2); d39 = T(1)/d39;
complex<T> d43 = spa12*t29*T(2); d43 = T(1)/d43;
complex<T> d45 = spb12*t29*T(2); d45 = T(1)/d45;
complex<T> t20 = (s12 - s124)*(-(d16*square(t150 + t57)) - t155*t81*t95*T(2)); 
complex<T> t21 = s124*s34*(s124*t35*t81*t95 + d37*square(t150 + t57)); 
complex<T> t24 = s23*(spa13*spb16 + spa23*spb26)*spb46*(spa25*spb45 + spa26*spb46)*t28*t34 + d14*s23*square(-(spa23*spb36) + t151); 
complex<T> t25 = d32*(spa13*spb16 + spa23*spb26)*(spa25*spb45 + spa26*spb46)*t154*t60 + d30*spa12*(spa15*spb14 + spa25*spb24)*t155*t95 + d31*spb12*square(-(spa23*spb36) + t151) - d29*spa12*square(t150 + t57); 
complex<T> t26 = d44*spa23*spa56*(spa13*spb16 + spa23*spb26)*(spa25*spb45 + spa26*spb46)*t154 - d46*(spa15*spb14 + spa25*spb24)*spb56*t155*t95 + d43*spa56*square(-(spa23*spb36) + t151) + d45*spb56*square(t150 + t57); 
complex<T> t38 = spa23*spb36*t16 + t79*t98; 
complex<T> t89 = spb56*t125; 
complex<T> t90 = spa56*t126; 
complex<T> t192 = d21*t126; 
complex<T> t195 = t160*t72; 
complex<T> t201 = t35*t67; 
complex<T> t208 = t63*t80; 
complex<T> t214 = d39*s34; 
complex<T> t216 = d12*spb12; 
complex<T> t219 = t14*t163; 
complex<T> t222 = t187*t28; 
complex<T> t229 = t14*t189; 
complex<T> t230 = spa45*t78; 
complex<T> t242 = d27*t200; 
complex<T> t247 = spb12*t82; 
complex<T> t1 = d1*s12*spa12*t186 - d2*spa12*t186*T(2) - spa12*t15*t230*T(2) - d15*t155*(t150 + t59)*T(2) + d16*(-square(t150 + t59) - spa45*spb13*t35*T(2)); 
complex<T> t3 = d14*(d42*s12*s56 - d40*(s12 - s34 + s56))*t10 + d16*(d41*s12*s56 - d40*(s12 - s34 + s56))*t11 - (d40*(s12 - s34 + s56)*(t222 + t15*t229) - s12*s56*(d42*t222 + d41*t15*t229))*T(2); 
complex<T> t4 = -(d35*spa24*spb36*(t150 + t57)) - d6*spa23*(spa14*spb13 + spa24*spb23)*spb46*t16*(spa56*t232 + spa13*spb14*(t150 + t57) + spa23*spb24*(t150 + t57)) - d36*spa45*spb13*(spa23*spb36 + t8) + d5*spb46*t34*t35*square(spa14*spb13 + spa24*spb23) + d34*spa24*spa45*spb13*spb36*T(2) - (spa14*spb13 + spa24*spb23)*t201*(spa12*spb56*t35 + spa23*spb24*t8 + spa13*spb14*(spa23*spb36 + t8) + spb24*spb36*square(spa23))*T(2) + t19*(spa23*spa45*spb13*spb46 + spa24*spb36*t35 + d34*s34*(s12 - s34 + s56)*t63)*T(6); 
complex<T> t53 = (s123 - s23)*(d14*t10 + t222*T(2)); 
complex<T> t87 = s123*(d33*s123*s34*t187 + t10*t214); 
complex<T> t139 = d37*s34*t11 + t10*t214 + s34*(t222 + t15*t229); 
complex<T> t179 = d26*s56*t200 - (t242 + t254)*T(2); 
complex<T> t220 = t208*t94; 
complex<T> t260 = t198*t216; 
complex<T> t263 = spb13*t219; 
complex<T> t2 = d11*s12*spb12*t198 - spb36*t12*t247*T(2) - t260*T(2) - d13*spa23*t154*(t58 + t8)*T(2) + d14*(-square(t58 + t8) - spa23*spb36*t58*T(2)); 
complex<T> t5 = -(d23*s56*spb56*t185) - d26*s56*t200 + d20*spa25*spb56*t45 - d20*spa56*spb16*t52 + spb46*t165*t34*t90 + t242*T(2) + t254*T(2) + t256*T(2) + t263*T(2) - t35*t71*t89*T(2) + d22*t113*(d19*spb56*square(spa25) + d21*spa56*square(spb16) - t218*T(2)) + t195*T(6); 
complex<T> t6 = -(d1*s12*spa12*t186) - d11*s12*spb12*t198 + t159*t16*t232 + t19*t65 + d5*t37*t8 + d5*t150*(spa23*spb36*t16 + t79*t98) + d2*spa12*t186*T(2) - spa12*(s124*spb16 + spb56*t150)*t201*T(2) + spa12*t15*t230*T(2) + spb36*t12*t247*T(2) + t260*T(2) + t220*T(6); 
complex<T> t7 = d5*t151*t37 - d20*spa25*spb56*t45 + d20*spa56*spb16*t52 - t19*t65 + d13*t154*t34*(t58 + t8) + d5*t9*(spa23*spb36*t16 + t79*t98) + d14*(spb36*t34*t58 + square(t58 + t8)) + spa12*(s124*spb16 + spb56*t150)*t201*T(2) - t110*t159*t232*T(2) + d15*t155*(t150 + t59)*T(2) + t35*t71*t89*T(2) - spa23*spb46*t165*t90*T(2) - d22*t113*(d19*spb56*square(spa25) + d21*spa56*square(spb16) - t218*T(2)) + d16*(square(t150 + t59) + spa45*spb13*t35*T(2)) - t195*T(6) - t220*T(6); 
complex<T> t177 = d23*s56*spb56*t185 - (t256 + t263)*T(2); 
complex<T> co1 = Complex(0,1)*t6; 
complex<T> co2 = Complex(0,1)*t2; 
complex<T> co3 = Complex(0,1)*t1; 
complex<T> co4 = Complex(0,1)*t7; 
complex<T> co5 = Complex(0,1)*t177; 
complex<T> co6 = Complex(0,1)*t179; 
complex<T> co7 = Complex(0,1)*t5; 
complex<T> co8 = Complex(0,1)*t25; 
complex<T> co9 = Complex(0,1)*t53; 
complex<T> co10 = Complex(0,1)*t4; 
complex<T> co11 = Complex(0,1)*t24; 
complex<T> co12 = Complex(0,1)*t139; 
complex<T> co13 = Complex(0,1)*t3; 
complex<T> co14 = Complex(0,1)*t20; 
complex<T> co15 = Complex(0,1)*t26; 
complex<T> co16 = Complex(0,1)*t51; 
complex<T> co17 = Complex(0,1)*t21; 
complex<T> co18 = Complex(0,1)*t87; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + co3*(*CI_users[2]->get_value(mc,ind,mu)) + co4*(*CI_users[3]->get_value(mc,ind,mu)) + co5*(*CI_users[4]->get_value(mc,ind,mu)) + co6*(*CI_users[5]->get_value(mc,ind,mu)) + co7*(*CI_users[6]->get_value(mc,ind,mu)) + co8*(*CI_users[7]->get_value(mc,ind,mu)) + co9*(*CI_users[8]->get_value(mc,ind,mu)) + co10*(*CI_users[9]->get_value(mc,ind,mu)) + co11*(*CI_users[10]->get_value(mc,ind,mu)) + co12*(*CI_users[11]->get_value(mc,ind,mu)) + co13*(*CI_users[12]->get_value(mc,ind,mu)) + co14*(*CI_users[13]->get_value(mc,ind,mu)) + co15*(*CI_users[14]->get_value(mc,ind,mu)) + co16*(*CI_users[15]->get_value(mc,ind,mu)) + co17*(*CI_users[16]->get_value(mc,ind,mu)) + co18*(*CI_users[17]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g2l_qpqmmpemep_VECT_wCI::\
C2q2g2l_qpqmmpemep_VECT_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c356, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c456, c123));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c56, c34, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmmpemep_VECT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, p, em, ep}, VECT}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmmpemep VECT");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:50:48 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa46 = SPA(4,6);
complex<T> spb36 = SPB(3,6);
complex<T> spb56 = SPB(5,6);
complex<T> spb15 = SPB(1,5);
complex<T> spb16 = SPB(1,6);
complex<T> spb34 = SPB(3,4);
complex<T> spa56 = SPA(5,6);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa26 = SPA(2,6);
complex<T> spa13 = SPA(1,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb46 = SPB(4,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb45 = SPB(4,5);
complex<T> spa36 = SPA(3,6);
complex<T> s56 = -(spa56*spb56);
complex<T> s123 = SS(1,2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> s124 = SS(1,2,4);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s24 = -(spa24*spb24);
complex<T> s35 = -(spa35*spb35);
complex<T> s36 = -(spa36*spb36);
complex<T> s45 = -(spa45*spb45);
complex<T> s46 = -(spa46*spb46);
complex<T> t12 = -(spa45*(spa45*spb15 + spa46*spb16)); 
complex<T> t16 = (s13 - s14 + s23 - s24)*T(2); 
complex<T> t35 = spa24*T(2); 
complex<T> t36 = spa45*spb13; 
complex<T> t57 = spa23*spb36; 
complex<T> t58 = spa24*spb46; 
complex<T> t59 = spa45*spb14; 
complex<T> t62 = (spa23*spb13 + spa24*spb14)*(spa35*spb36 + spa45*spb46); 
complex<T> t63 = -(s123*(spa45*spb15 + spa46*spb16)); 
complex<T> t90 = spa14*spb13 + spa24*spb23; 
complex<T> t93 = -s12 - s34 + s56; 
complex<T> t95 = -s35 - s36 + s45 + s46; 
complex<T> t96 = square(spa23); 
complex<T> t97 = square(spa24); 
complex<T> t98 = square(spa35); 
complex<T> t99 = square(spa45); 
complex<T> t100 = square(spb13); 
complex<T> t101 = square(spb14); 
complex<T> t102 = -(spa15*spb13) - spa25*spb23; 
complex<T> t104 = -(spa14*spb16) - spa24*spb26; 
complex<T> t105 = square(spb36); 
complex<T> t106 = -(spa25*spb35) - spa26*spb36; 
complex<T> t107 = square(spb46); 
complex<T> t109 = -(spa56*T(2)); 
complex<T> t111 = -(spa35*spb45) - spa36*spb46; 
complex<T> t144 = -(s12*T(2)); 
complex<T> t149 = spa25*spb12; 
complex<T> t150 = spa12*spb16; 
complex<T> t154 = -(spa14*spb13); 
complex<T> t166 = (s35 + s36 - s45 - s46)*T(2); 
complex<T> t179 = -(s34*T(2)); 
complex<T> t180 = spa24*spb36; 
complex<T> t216 = spa25*spb16; 
complex<T> t236 = spb14*spb56; 
complex<T> d5 = spa56; d5 = T(1)/d5;
complex<T> d8 = spb56; d8 = T(1)/d8;
complex<T> d19 = spa12; d19 = T(1)/d19;
complex<T> d21 = spb12; d21 = T(1)/d21;
complex<T> d23 = spa12*square(s124 - s56)*square(spa45*spb35 + spa46*spb36); d23 = T(1)/d23;
complex<T> d24 = (-s124 + s56)*spa12*square(spa45*spb35 + spa46*spb36); d24 = T(1)/d24;
complex<T> d25 = (-s124 + s56)*spa12*cube(spa45*spb35 + spa46*spb36); d25 = T(1)/d25;
complex<T> d26 = spb12*square(s123 - s56)*square(spa45*spb35 + spa46*spb36); d26 = T(1)/d26;
complex<T> d27 = (-s123 + s56)*spb12*square(spa45*spb35 + spa46*spb36); d27 = T(1)/d27;
complex<T> d28 = (-s123 + s56)*spb12*cube(spa45*spb35 + spa46*spb36); d28 = T(1)/d28;
complex<T> d40 = T(2); d40 = T(1)/d40;
complex<T> d41 = s123; d41 = T(1)/d41;
complex<T> d42 = s124; d42 = T(1)/d42;
complex<T> d47 = s123*(spa45*spb35 + spa46*spb36)*T(2); d47 = T(1)/d47;
complex<T> d48 = s124*(spa45*spb35 + spa46*spb36)*T(2); d48 = T(1)/d48;
complex<T> t8 = -t149; 
complex<T> t9 = -t150; 
complex<T> t10 = square(spa35*spb13 + t149); 
complex<T> t11 = square(t150 - t58); 
complex<T> t13 = spb13*t102; 
complex<T> t14 = spa24*t104; 
complex<T> t15 = spb36*t106; 
complex<T> t29 = square(t90); 
complex<T> t31 = d5*(s12 - s34 - s56)*(s124*spa25 + spa56*t150) + (-s13 + s14 - s23 + s24)*spb46*t35; 
complex<T> t37 = -(d8*(s12 - s34 - s56)*(s123*spb16 + spb56*t149)) + spa35*spb13*t16; 
complex<T> t51 = -(d48*spa35*spb14*(spa25*spb56 + t57)) + d47*spa23*spb46*(spa56*spb16 - t59); 
complex<T> t52 = t166*t57 - d21*(s124*spb16 + spb56*t149)*t93; 
complex<T> t64 = d5*spb12*square(spa25) + d8*spa12*square(spb16) - t216*T(2); 
complex<T> t75 = -(d19*t93); 
complex<T> t78 = spb12*(s124*spa25 + spa56*t150); 
complex<T> t82 = d28*(spa45*spb34 - spa56*spb36); 
complex<T> t84 = d25*(-(spa34*spb36) + spa45*spb56); 
complex<T> t85 = spa12*(s123*spb16 + spb56*t149); 
complex<T> t158 = spb56*t99; 
complex<T> t161 = t111*t62; 
complex<T> t183 = t100*t98; 
complex<T> t200 = t105*t96; 
complex<T> t203 = spb12*t97; 
complex<T> t207 = t102*t63; 
complex<T> t261 = t59*t95; 
complex<T> d3 = (s12 - s123)*spa56*cube(t90); d3 = T(1)/d3;
complex<T> d4 = cube(t90)*(s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34)); d4 = T(1)/d4;
complex<T> d7 = t90*square(s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34)); d7 = T(1)/d7;
complex<T> d9 = t90*(s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34)); d9 = T(1)/d9;
complex<T> d10 = (s12 - s124)*spb56*cube(t90); d10 = T(1)/d10;
complex<T> d13 = spa56*spb12*cube(t90); d13 = T(1)/d13;
complex<T> d15 = spa12*spb56*cube(t90); d15 = T(1)/d15;
complex<T> d17 = (spa45*spb35 + spa46*spb36)*square(s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34)); d17 = T(1)/d17;
complex<T> d18 = cube(spa45*spb35 + spa46*spb36)*(s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34)); d18 = T(1)/d18;
complex<T> d20 = (s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34))*square(spa45*spb35 + spa46*spb36); d20 = T(1)/d20;
complex<T> d22 = (spa45*spb35 + spa46*spb36)*(s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34)); d22 = T(1)/d22;
complex<T> d30 = spa56*square(square(t90)); d30 = T(1)/d30;
complex<T> d31 = spb56*square(square(t90)); d31 = T(1)/d31;
complex<T> d33 = spa56*spb12*square(square(t90)); d33 = T(1)/d33;
complex<T> d34 = s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34); d34 = T(1)/d34;
complex<T> d35 = s123*t90*T(2); d35 = T(1)/d35;
complex<T> d36 = s124*t90*T(2); d36 = T(1)/d36;
complex<T> d38 = spa12*spb56*square(square(t90)); d38 = T(1)/d38;
complex<T> d43 = spa12*square(square(t90)); d43 = T(1)/d43;
complex<T> d46 = spb12*square(square(t90)); d46 = T(1)/d46;
complex<T> t21 = -(d9*(spa13*spb14 + spa23*spb24)); 
complex<T> t28 = d38*s124; 
complex<T> t45 = -((s123*spa25 + spa56*t150)*t75) + t261*T(2); 
complex<T> t66 = d4*(-s13 + s14 - s23 + s24); 
complex<T> t70 = d18*(s35 + s36 - s45 - s46); 
complex<T> t71 = d17*(s12 - s34 - s56); 
complex<T> t76 = d3*(spa24*spb12 + spa34*spb13); 
complex<T> t77 = d13*(spa45*spb14 + t149); 
complex<T> t80 = -(d7*(spa13*spb14 + spa23*spb24)); 
complex<T> t81 = -(d10*(spa12*spb13 + spa24*spb34)); 
complex<T> t88 = d15*(spa23*spb36 + t9); 
complex<T> t124 = -(s124*spb16) + spb56*t8; 
complex<T> t125 = -(s123*spa25) + spa56*t9; 
complex<T> t153 = d33*s123; 
complex<T> t198 = d18*t95; 
complex<T> t214 = t14*t15; 
complex<T> t225 = t101*t158; 
complex<T> t239 = t15*t84; 
complex<T> t248 = t12*t82; 
complex<T> d1 = spa56*t29*square(s12 - s123); d1 = T(1)/d1;
complex<T> d2 = (s12 - s123)*spa56*t29; d2 = T(1)/d2;
complex<T> d6 = t29*(s34*t144 + s56*(s56 + t144 + t179) + square(s12) + square(s34)); d6 = T(1)/d6;
complex<T> d11 = spb56*t29*square(s12 - s124); d11 = T(1)/d11;
complex<T> d12 = (s12 - s124)*spb56*t29; d12 = T(1)/d12;
complex<T> d14 = spa56*spb12*t29; d14 = T(1)/d14;
complex<T> d16 = spa12*spb56*t29; d16 = T(1)/d16;
complex<T> d29 = spa56*t29*T(2); d29 = T(1)/d29;
complex<T> d32 = spb56*t29*T(2); d32 = T(1)/d32;
complex<T> d37 = spa56*spb12*t29*T(2); d37 = T(1)/d37;
complex<T> d39 = spa12*spb56*t29*T(2); d39 = T(1)/d39;
complex<T> d44 = spa12*t29*T(2); d44 = T(1)/d44;
complex<T> d45 = spb12*t29*T(2); d45 = T(1)/d45;
complex<T> t4 = -(d35*spa23*spb46*(spa35*spb13 + t149)) + (-(spa13*spb14) - spa23*spb24)*spb36*((spa35*spb13 + t149)*(spa24*spb23 - t154) + spa56*spb12*t180)*t35*t66 - d36*spa35*spb14*(t58 + t9) + d4*(-(spa13*spb14) - spa23*spb24)*t16*t36*(spa12*spb56*t36 + (spa24*spb23 - t154)*(t58 + t9)) + d6*spb36*t35*t36*square(spa13*spb14 + spa23*spb24) + d34*spa23*spa35*spb14*spb46*T(2) + t21*(spa35*spb14*t180 + spa23*spb46*t36 + d34*s34*(s12 - s34 + s56)*t62)*T(6); 
complex<T> t19 = d14*s23*square(spa35*spb13 + t149) + d33*s23*t207*t36*T(2); 
complex<T> t20 = s123*s34*(-(d33*s123*(spa45*spb15 + spa46*spb16)*t102*t36) + d37*square(spa35*spb13 + t149)); 
complex<T> t24 = (s12 - s124)*(-(d16*square(t150 - t58)) - t104*t106*t180*t28*T(2)); 
complex<T> t25 = -(d31*s124*spb12*t104*t106*t180) + d30*spa12*t207*t36 + d29*spa12*square(spa35*spb13 + t149) - d32*spb12*square(t150 - t58); 
complex<T> t26 = d43*s124*spa56*t104*t106*t180 + d46*spb56*t207*t36 + d45*spb56*square(spa35*spb13 + t149) + d44*spa56*square(t150 - t58); 
complex<T> t87 = spa56*t124; 
complex<T> t89 = spb56*t125; 
complex<T> t129 = d2*spa12; 
complex<T> t177 = d23*s56*spa56*t200 + d24*t109*t200 + spa23*t109*t239; 
complex<T> t184 = t12*t153; 
complex<T> t190 = d12*t107; 
complex<T> t191 = d39*t11; 
complex<T> t192 = d21*t124; 
complex<T> t195 = t161*t71; 
complex<T> t208 = spa12*t76; 
complex<T> t209 = t62*t80; 
complex<T> t210 = spb46*t81; 
complex<T> t252 = d27*t225; 
complex<T> t259 = t236*t248; 
complex<T> t5 = -(d23*s56*spa56*t200) - d26*s56*t225 + d20*spa25*spb56*t45 - d20*spa56*spb16*t52 + spb36*t35*t70*t87 + d24*spa56*t200*T(2) + spa23*spa56*t239*T(2) + t252*T(2) + t259*T(2) - t198*t36*t89*T(2) + d22*t111*(d19*spb56*square(spa25) + d21*spa56*square(spb16) - t216*T(2)) + t195*T(6); 
complex<T> t86 = s124*s34*(t191 + d38*s124*t214); 
complex<T> t175 = d26*s56*t225 - (t252 + t259)*T(2); 
complex<T> t215 = t13*t184; 
complex<T> t219 = spb12*t210; 
complex<T> t223 = t209*t93; 
complex<T> t251 = t129*t183; 
complex<T> t254 = t190*t203; 
complex<T> t1 = d1*s12*spa12*t183 - (spa35*t13*t208 + t251 + s123*t36*t77)*T(2) - d14*(square(t149) + square(t59) + spa35*spb14*t36*T(2) + t149*t59*T(2)); 
complex<T> t3 = d14*(d41*s12*s56 - d40*(s12 - s34 + s56))*t10 + d16*(d42*s12*s56 - d40*(s12 - s34 + s56))*t11 - (d40*(s12 - s34 + s56)*(t215 + t214*t28) - s12*s56*(d41*t215 + d42*t214*t28))*T(2); 
complex<T> t7 = d6*t150*(-(d8*(s12 - s34 - s56)*(s123*spb16 + spb56*t149)) + spa35*spb13*t16) - d20*spa25*spb56*t45 + d20*spa56*spb16*t52 - t21*t64 + t109*t124*t180*t70 + d6*t31*t8 + d4*t16*t36*t85 + d15*s124*spb36*t35*(t57 + t9) + d16*(spb46*t35*t57 + square(t57 + t9)) + s123*t36*t77*T(2) - t180*t66*t78*T(2) + t198*t36*t89*T(2) - d22*t111*(d19*spb56*square(spa25) + d21*spa56*square(spb16) - t216*T(2)) + d14*(square(t149 + t59) + spa35*spb14*t36*T(2)) - t195*T(6) - t223*T(6); 
complex<T> t53 = (s123 - s23)*(d14*t10 + t215*T(2)); 
complex<T> t138 = s34*(d37*t10 + t191 + t215 + t214*t28); 
complex<T> t265 = t14*t219; 
complex<T> t2 = d11*s12*t107*t203 - t254*T(2) - t265*T(2) - d15*s124*t180*(t57 + t9)*T(2) + d16*(-square(t57 + t9) - t57*t58*T(2)); 
complex<T> t6 = -(d1*s12*spa12*t183) - d11*s12*t107*t203 + d6*t149*t31 + t21*t64 + spb36*t35*t66*t78 + d6*(-(d8*(s12 - s34 - s56)*(s123*spb16 + spb56*t149)) + spa35*spb13*t16)*t9 + spa35*t13*t208*T(2) + t251*T(2) + t254*T(2) + t265*T(2) - d4*(s13 - s14 + s23 - s24)*t36*t85*T(2) + t223*T(6); 
complex<T> co1 = Complex(0,1)*t6; 
complex<T> co2 = Complex(0,1)*t1; 
complex<T> co3 = Complex(0,1)*t2; 
complex<T> co4 = Complex(0,1)*t7; 
complex<T> co5 = Complex(0,1)*t177; 
complex<T> co6 = Complex(0,1)*t175; 
complex<T> co7 = Complex(0,1)*t5; 
complex<T> co8 = Complex(0,1)*t25; 
complex<T> co9 = Complex(0,1)*t53; 
complex<T> co10 = Complex(0,1)*t4; 
complex<T> co11 = Complex(0,1)*t19; 
complex<T> co12 = Complex(0,1)*t138; 
complex<T> co13 = Complex(0,1)*t3; 
complex<T> co14 = Complex(0,1)*t24; 
complex<T> co15 = Complex(0,1)*t26; 
complex<T> co16 = Complex(0,1)*t51; 
complex<T> co17 = Complex(0,1)*t86; 
complex<T> co18 = Complex(0,1)*t20; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + co3*(*CI_users[2]->get_value(mc,ind,mu)) + co4*(*CI_users[3]->get_value(mc,ind,mu)) + co5*(*CI_users[4]->get_value(mc,ind,mu)) + co6*(*CI_users[5]->get_value(mc,ind,mu)) + co7*(*CI_users[6]->get_value(mc,ind,mu)) + co8*(*CI_users[7]->get_value(mc,ind,mu)) + co9*(*CI_users[8]->get_value(mc,ind,mu)) + co10*(*CI_users[9]->get_value(mc,ind,mu)) + co11*(*CI_users[10]->get_value(mc,ind,mu)) + co12*(*CI_users[11]->get_value(mc,ind,mu)) + co13*(*CI_users[12]->get_value(mc,ind,mu)) + co14*(*CI_users[13]->get_value(mc,ind,mu)) + co15*(*CI_users[14]->get_value(mc,ind,mu)) + co16*(*CI_users[15]->get_value(mc,ind,mu)) + co17*(*CI_users[16]->get_value(mc,ind,mu)) + co18*(*CI_users[17]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g2l_qpqmmmemep_VECT_wCI::\
C2q2g2l_qpqmmmemep_VECT_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c356, c124));
CI_users.push_back(new Cached_Bubble_Integral_User(c456, c123));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c12, c3, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmmmemep_VECT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, m, em, ep}, VECT}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmmmemep VECT");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:51:00 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb36 = SPB(3,6);
complex<T> spb46 = SPB(4,6);
complex<T> spb56 = SPB(5,6);
complex<T> spa45 = SPA(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> s23 = S(2,3);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s12 = S(1,2);
complex<T> s124 = SS(1,2,4);
complex<T> t1 = square(spb34); 
complex<T> t13 = spb14*spb36; 
complex<T> t14 = spb13*spb46; 
complex<T> t15 = square(spb16); 
complex<T> t22 = s12 - s124; 
complex<T> t23 = -s124 + s56; 
complex<T> t24 = -(spa35*spb13); 
complex<T> t25 = -(spa45*spb14); 
complex<T> t29 = spb12*spb56; 
complex<T> t35 = spb14*spb36 + spb13*spb46; 
complex<T> t40 = square(spa23); 
complex<T> t41 = square(spb36); 
complex<T> t48 = square(spa24); 
complex<T> t49 = square(spa35); 
complex<T> t50 = square(spa45); 
complex<T> t51 = square(spb13); 
complex<T> t52 = square(spb14); 
complex<T> t53 = square(spb46); 
complex<T> d5 = (s12 - s123)*spb56*cube(spb34); d5 = T(1)/d5;
complex<T> d7 = (-s123 + s56)*spb12*cube(spb34); d7 = T(1)/d7;
complex<T> d17 = spb12*square(square(spb34)); d17 = T(1)/d17;
complex<T> t57 = t13*t14; 
complex<T> d1 = (s12 - s123)*spb56*t1; d1 = T(1)/d1;
complex<T> d2 = spb56*t1*square(s12 - s123); d2 = T(1)/d2;
complex<T> d3 = spb56*t1*t22; d3 = T(1)/d3;
complex<T> d4 = spb56*t1*square(t22); d4 = T(1)/d4;
complex<T> d6 = spb56*t22*cube(spb34); d6 = T(1)/d6;
complex<T> d8 = t29*cube(spb34); d8 = T(1)/d8;
complex<T> d9 = spb12*t23*cube(spb34); d9 = T(1)/d9;
complex<T> d10 = spb12*t1*t23; d10 = T(1)/d10;
complex<T> d11 = spb12*t1*square(t23); d11 = T(1)/d11;
complex<T> d12 = (-s123 + s56)*spb12*t1; d12 = T(1)/d12;
complex<T> d13 = spb12*t1*square(s123 - s56); d13 = T(1)/d13;
complex<T> d14 = t1*t29; d14 = T(1)/d14;
complex<T> d15 = t29*square(square(spb34)); d15 = T(1)/d15;
complex<T> d16 = spb12*t1; d16 = T(1)/d16;
complex<T> d18 = t1*t29*T(2); d18 = T(1)/d18;
complex<T> d19 = t29*square(square(spb34))*T(2); d19 = T(1)/d19;
complex<T> t8 = -(d1*spa23*spb16*spb36) + d7*spa45*spb14*t35 + d8*spb16*t35 + d5*spa23*spb36*t35 - d2*spb12*t40*t41; 
complex<T> t9 = -(d3*spa24*spb16*spb46) - d8*spb16*t35 - d6*spa24*spb46*t35 + d9*t24*t35 - d4*spb12*t48*t53; 
complex<T> t10 = d1*spa23*spb16*spb36 + d3*spa24*spb16*spb46 - d5*spa23*spb36*t35 + d6*spa24*spb46*t35 + d2*spb12*t40*t41 + d4*spb12*t48*t53; 
complex<T> t11 = d10*spa35*spb13*spb16 + d12*spa45*spb14*spb16 + d9*spa35*spb13*t35 + d7*t25*t35 + d11*spb56*t49*t51 + d13*spb56*t50*t52; 
complex<T> t16 = d19*(s123*s124 - s12*s56); 
complex<T> t19 = s23*(d14*square(spb16) + d15*spb13*spb14*spb36*spb46*T(2)); 
complex<T> t20 = spa56*(d16*square(spb16) + d17*spb13*spb14*spb36*spb46*T(2)); 
complex<T> t42 = -(d15*T(2)); 
complex<T> t46 = d10*spb16*t24 - d11*spb56*t49*t51; 
complex<T> t47 = d12*spb16*t25 - d13*spb56*t50*t52; 
complex<T> t6 = (s123 - s23)*(d14*t15 - t42*t57); 
complex<T> t56 = -(d14*t15*t22) + t22*t42*t57; 
complex<T> t80 = t16*t57; 
complex<T> t7 = d18*(s123*s124 - s12*s56)*t15 + t80; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t10*(*CI_users[0]->get_value(mc,ind,mu)) + t8*(*CI_users[1]->get_value(mc,ind,mu)) + t9*(*CI_users[2]->get_value(mc,ind,mu)) + t46*(*CI_users[3]->get_value(mc,ind,mu)) + t47*(*CI_users[4]->get_value(mc,ind,mu)) + t11*(*CI_users[5]->get_value(mc,ind,mu)) + t6*(*CI_users[6]->get_value(mc,ind,mu)) + t19*(*CI_users[7]->get_value(mc,ind,mu)) + t56*(*CI_users[8]->get_value(mc,ind,mu)) + t20*(*CI_users[9]->get_value(mc,ind,mu)) + t7*(*CI_users[10]->get_value(mc,ind,mu)) + t80*(*CI_users[11]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmppemep_AX_wCI::\
C2q2g2l_qpqmppemep_AX_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c4, c12, c3, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmppemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, p, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmppemep AX");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:51:05 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spb56 = SPB(5,6);
complex<T> spb13 = SPB(1,3);
complex<T> spa13 = SPA(1,3);
complex<T> spb46 = SPB(4,6);
complex<T> spb34 = SPB(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> spb36 = SPB(3,6);
complex<T> s23 = S(2,3);
complex<T> s12 = S(1,2);
complex<T> s123 = SS(1,2,3);
complex<T> s14 = S(1,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s124 = SS(1,2,4);
complex<T> t4 = square(spa34); 
complex<T> t6 = spa24*spa35 + spa23*spa45; 
complex<T> t7 = cube(spa34); 
complex<T> t9 = s12 - s124; 
complex<T> t14 = -s123 + s23; 
complex<T> t16 = s123*s124 - s12*s56; 
complex<T> t17 = s14 + s34; 
complex<T> t27 = square(spa25); 
complex<T> t28 = spa23*spa25; 
complex<T> t37 = spa25*spb46; 
complex<T> d5 = spa13*spa34*square(s123 - s56); d5 = T(1)/d5;
complex<T> d6 = spa12*spa34*square(s123 - s56); d6 = T(1)/d6;
complex<T> d8 = spa24*spa34*square(s124 - s56); d8 = T(1)/d8;
complex<T> d9 = spa12*spa34*square(s124 - s56); d9 = T(1)/d9;
complex<T> t8 = -t37; 
complex<T> t18 = -t28; 
complex<T> t24 = d6*spb34; 
complex<T> t29 = d8*spb13; 
complex<T> t31 = d5*t17; 
complex<T> t34 = d9*spb34; 
complex<T> d1 = (s12 - s123)*spa56*t4; d1 = T(1)/d1;
complex<T> d2 = spa56*t4*t9; d2 = T(1)/d2;
complex<T> d3 = spa12*spa56*t4; d3 = T(1)/d3;
complex<T> d4 = (-s123 + s56)*spa12*t4; d4 = T(1)/d4;
complex<T> d7 = (-s124 + s56)*spa12*t4; d7 = T(1)/d7;
complex<T> d10 = spa12*spa56*t7*T(2); d10 = T(1)/d10;
complex<T> d11 = spa12*t7*T(2); d11 = T(1)/d11;
complex<T> d12 = spa12*spa56*t7*T(4); d12 = T(1)/d12;
complex<T> t20 = d4*spa24; 
complex<T> t21 = -(d1*spa35); 
complex<T> t30 = d7*spb36; 
complex<T> t33 = d2*spa45; 
complex<T> t38 = spa24*t24; 
complex<T> t39 = spb36*t29; 
complex<T> t5 = d3*t27 + t18*t30 - spa25*spb14*t33 + spb36*t28*t34 + t28*t39; 
complex<T> t26 = d1*spa25*spa35*spb13 - d3*t27 + t20*t37 + t31*t8 + t38*t8; 
complex<T> t36 = spa25*(spb13*t21 + spb14*t33); 
complex<T> t41 = t28*t30 + spb36*t18*t34 + t31*t37 + t37*t38 + t18*t39 + t20*t8; 
complex<T> co1 = d10*spa25*t14*t6; 
complex<T> co2 = -(d10*s23*spa25*t6); 
complex<T> co3 = d10*spa25*t6*t9; 
complex<T> co4 = -(d11*spa25*spb56*t6); 
complex<T> co5 = -(d12*spa25*t16*t6); 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t36*(*CI_users[0]->get_value(mc,ind,mu)) + t26*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + t41*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)) + co4*(*CI_users[7]->get_value(mc,ind,mu)) + co5*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmpmemep_AX_wCI::\
C2q2g2l_qpqmpmemep_AX_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmpmemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, m, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmpmemep AX");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:52:44 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spa14 = SPA(1,4);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb36 = SPB(3,6);
complex<T> spb34 = SPB(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb46 = SPB(4,6);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa36 = SPA(3,6);
complex<T> spb35 = SPB(3,5);
complex<T> spa46 = SPA(4,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb26 = SPB(2,6);
complex<T> spa15 = SPA(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s35 = -(spa35*spb35);
complex<T> s36 = -(spa36*spb36);
complex<T> s45 = -(spa45*spb45);
complex<T> s46 = -(spa46*spb46);
complex<T> s123 = SS(1,2,3);
complex<T> s124 = SS(1,2,4);
complex<T> t11 = square(spa13*spb16 + spa23*spb26)*square(spa12*spb14 - spa23*spb34) - square(s123)*square(spa23)*square(spb46); 
complex<T> t12 = -(square(s124)*square(spa35)*square(spb14)) + square(spa23*spb12 - spa34*spb14)*square(spa15*spb14 + spa25*spb24); 
complex<T> t22 = square(spa35*spb45 + spa36*spb46); 
complex<T> t23 = spa12*T(2); 
complex<T> t44 = -(spa12*spb16); 
complex<T> t45 = spa25*spb12; 
complex<T> t54 = -(spa23*spa45); 
complex<T> t58 = -(spb14*spb36); 
complex<T> t61 = -s12 + s34 - s56; 
complex<T> t62 = spa13*spb14 + spa23*spb24; 
complex<T> t63 = -(spa14*spb13) - spa24*spb23; 
complex<T> t65 = -(spa45*spb35) - spa46*spb36; 
complex<T> t66 = -(spa15*spb14) - spa25*spb24; 
complex<T> t67 = -(spa13*spb12) - spa34*spb24; 
complex<T> t68 = -(spa13*spb16) - spa23*spb26; 
complex<T> t69 = spa12*spb24 + spa13*spb34; 
complex<T> t73 = s12 - s124; 
complex<T> t83 = -s123 + s23; 
complex<T> t94 = -(s12*T(2)); 
complex<T> t95 = square(s34); 
complex<T> t96 = square(s56); 
complex<T> t97 = spa56*spb12; 
complex<T> t98 = spa24*spb46; 
complex<T> t99 = square(spa45); 
complex<T> t100 = spa35*spb13; 
complex<T> t102 = spa12*spb56; 
complex<T> t111 = -(spa23*spb24); 
complex<T> t125 = square(spa24); 
complex<T> t126 = square(spb13); 
complex<T> t135 = s34*s56; 
complex<T> t138 = square(spb36); 
complex<T> t139 = -(spa45*spb14); 
complex<T> d33 = T(2); d33 = T(1)/d33;
complex<T> d34 = s124; d34 = T(1)/d34;
complex<T> d36 = s123; d36 = T(1)/d36;
complex<T> t21 = square(t62); 
complex<T> t26 = spa25*t61 - spa12*spa56*spb16*T(2); 
complex<T> t27 = s13*spa35 - s14*spa35 + s23*spa35 - s24*spa35 + spa13*t139*T(3) + spb24*t54*T(3); 
complex<T> t29 = -(s13*spb46) + s14*spb46 - s23*spb46 + s24*spb46 + spb36*t111*T(3) + spa13*t58*T(3); 
complex<T> t30 = (-s35 - s36 + s45 + s46)*spa23 + (spa24*spa35*spb45 + spa36*t98)*T(3); 
complex<T> t31 = (s35 + s36 - s45 - s46)*spb14 + (spa36*spb13*spb46 + spb45*t100)*T(3); 
complex<T> t32 = -(spb16*t61) + spb56*t45*T(2); 
complex<T> t43 = square(square(t62)); 
complex<T> t46 = -t100; 
complex<T> t49 = -((spa25*spb12 + spa45*spb14)*t66); 
complex<T> t59 = spa23*spb36 + t44; 
complex<T> t60 = spa45*spb14 + t45; 
complex<T> t64 = -(spa56*spb16) + t100; 
complex<T> t70 = spa25*spb56 + t98; 
complex<T> t170 = spb12*t63; 
complex<T> d1 = spa34*spa56*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d1 = T(1)/d1;
complex<T> d5 = t62*square(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d5 = T(1)/d5;
complex<T> d6 = spb34*spb56*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d6 = T(1)/d6;
complex<T> d11 = spa13*t102*t62*square(s123 - s56); d11 = T(1)/d11;
complex<T> d12 = spb24*t62*t97*square(s124 - s56); d12 = T(1)/d12;
complex<T> d15 = spa12*spa34*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d15 = T(1)/d15;
complex<T> d16 = spb12*spb34*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d16 = T(1)/d16;
complex<T> d17 = spb12*spb34*t22; d17 = T(1)/d17;
complex<T> d18 = spb12*spb34*t22*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d18 = T(1)/d18;
complex<T> d20 = (spa35*spb45 + spa36*spb46)*square(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d20 = T(1)/d20;
complex<T> d21 = spa12*spa34*t22; d21 = T(1)/d21;
complex<T> d22 = spa12*spa34*t22*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d22 = T(1)/d22;
complex<T> d26 = s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2); d26 = T(1)/d26;
complex<T> d27 = s124*t62*T(2); d27 = T(1)/d27;
complex<T> d29 = s123*t62*T(2); d29 = T(1)/d29;
complex<T> d30 = t62*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d30 = T(1)/d30;
complex<T> t18 = d5*(spb16*t61 - spb56*t45*T(2)); 
complex<T> t19 = t60*t63; 
complex<T> t20 = spb14*t59; 
complex<T> t25 = square(t60); 
complex<T> t28 = square(t59); 
complex<T> t51 = t65*t70; 
complex<T> t103 = t64*t65; 
complex<T> t105 = d5*t59; 
complex<T> t108 = d11*s123; 
complex<T> t109 = d12*s124; 
complex<T> t118 = d20*t32; 
complex<T> t140 = t60*t66; 
complex<T> t141 = t59*t68; 
complex<T> t145 = d20*spb56; 
complex<T> d2 = spa34*spa56*t21; d2 = T(1)/d2;
complex<T> d3 = spa56*t21*t73; d3 = T(1)/d3;
complex<T> d4 = spa34*spa56*t21*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d4 = T(1)/d4;
complex<T> d7 = (s12 - s123)*spb56*t21; d7 = T(1)/d7;
complex<T> d8 = spb34*spb56*t21; d8 = T(1)/d8;
complex<T> d9 = spb34*spb56*t21*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d9 = T(1)/d9;
complex<T> d10 = (-s123 + s56)*spa13*t102*t21; d10 = T(1)/d10;
complex<T> d13 = (-s124 + s56)*spb24*t21*t97; d13 = T(1)/d13;
complex<T> d14 = t21*t97; d14 = T(1)/d14;
complex<T> d19 = t102*t21; d19 = T(1)/d19;
complex<T> d23 = spa56*t43*T(4); d23 = T(1)/d23;
complex<T> d24 = spb56*t43*T(4); d24 = T(1)/d24;
complex<T> d25 = spb56*t23*t43; d25 = T(1)/d25;
complex<T> d28 = t21*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t135*T(2)); d28 = T(1)/d28;
complex<T> d31 = t43*t97*T(4); d31 = T(1)/d31;
complex<T> d32 = t102*t43*T(4); d32 = T(1)/d32;
complex<T> d35 = t43*t97*T(2); d35 = T(1)/d35;
complex<T> d37 = spa12*t43*T(4); d37 = T(1)/d37;
complex<T> d38 = spb12*t43*T(4); d38 = T(1)/d38;
complex<T> t3 = -(d24*spb12*t11) + d23*spa12*t12; 
complex<T> t4 = -(d37*spa56*t11) - d38*spb56*t12; 
complex<T> t6 = -(d31*t12); 
complex<T> t7 = d28*(-s13 + s14 - s23 + s24)*spa23*spb46*t19 + d29*spa45*spb13*t59 + d27*spa24*spb36*t60 + d28*(s13 - s14 + s23 - s24)*spa35*t20*t63 - d26*spa24*spa45*spb13*spb36*T(2) + d5*t19*(s34*(spa25*spb56 - t44) + s12*(spa25*spb56 + t44) - s56*(spa25*spb56 + t44))*t61*T(3) + t105*(-((s12 + s34 - s56)*spa56*spb16) + (s12 - s34 - s56)*t45)*t61*t63*T(3) + d30*t63*(s12*spa25*spb16 + s56*spa25*spb16 + spb14*spb36*t54 + t46*t98)*T(3); 
complex<T> t106 = d13*t67; 
complex<T> t107 = d10*t69; 
complex<T> t112 = d25*t11; 
complex<T> t146 = spa14*t108; 
complex<T> t147 = spb23*t109; 
complex<T> t148 = -(d32*t11); 
complex<T> t149 = d35*t12; 
complex<T> t150 = spa12*t18; 
complex<T> t169 = t118*t51; 
complex<T> t183 = t103*t145; 
complex<T> t5 = -t112; 
complex<T> t156 = spb13*t106*t140 + d3*spa24*t49 + t100*t147*t60; 
complex<T> t157 = d7*spb13*t141 + spa24*t107*t141 + t146*t59*t98; 
complex<T> t184 = s34*t148; 
complex<T> t205 = spa56*t169; 
complex<T> t207 = s34*t6; 
complex<T> t216 = t150*t19; 
complex<T> t1 = t184 + t207; 
complex<T> t2 = -(d34*s12*s56*t149) + d36*s12*s56*t5 + d33*(-t149 + t5)*t61; 
complex<T> t8 = d6*spa24*spb13*t138 + d3*spa24*t140 - d7*spb13*t141 - d4*spa23*t19*t27 + d8*spb13*spb46*t59 + d2*spa24*spa35*t60 - d9*t20*t29*t63 + d1*spa24*spb13*t99 + t216*T(6) - t105*t170*t26*T(6); 
complex<T> t9 = -(d15*spa45*spb36*t125) - d16*spa45*spb36*t126 - spa24*t107*t141 - d14*t25 - d19*t28 - d18*spb46*t103*t31 + spb13*t106*t49 - d22*spa35*t30*t51 + t147*t46*t60 + d17*t58*t64 + d21*t54*t70 - t146*t59*t98 - t205*T(6) + spa56*spb16*t183*t23*T(6) - spa25*t183*t61*T(6); 
complex<T> t10 = d15*spa45*spb36*t125 + d16*spa45*spb36*t126 - d6*spa24*spb13*t138 + d14*t25 + d4*spa23*t19*t27 + d19*t28 + d18*spb46*t103*t31 + d22*spa35*t30*t51 - d8*spb13*spb46*t59 - d2*spa24*spa35*t60 + d9*t20*t29*t63 + d17*spb14*spb36*t64 + d21*spa23*spa45*t70 - d1*spa24*spb13*t99 + t205*T(6) - t216*T(6) - spa56*spb16*t183*t23*T(6) + t105*t170*t26*T(6) + spa25*t183*t61*T(6); 
complex<T> co1 = t112*t83; 
complex<T> co2 = s23*t5; 
complex<T> co3 = t149*t73; 
complex<T> co4 = s124*t207; 
complex<T> co5 = s123*t184; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t8*(*CI_users[0]->get_value(mc,ind,mu)) + t157*(*CI_users[1]->get_value(mc,ind,mu)) + t156*(*CI_users[2]->get_value(mc,ind,mu)) + t10*(*CI_users[3]->get_value(mc,ind,mu)) + t9*(*CI_users[4]->get_value(mc,ind,mu)) + t3*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t7*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu)) + t1*(*CI_users[9]->get_value(mc,ind,mu)) + t2*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t4*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)) + co5*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmmpemep_AX_wCI::\
C2q2g2l_qpqmmpemep_AX_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c12, c56));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmmpemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, p, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmmpemep AX");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:54:25 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spa13 = SPA(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb46 = SPB(4,6);
complex<T> spb34 = SPB(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb36 = SPB(3,6);
complex<T> spa45 = SPA(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa46 = SPA(4,6);
complex<T> spb45 = SPB(4,5);
complex<T> spa36 = SPA(3,6);
complex<T> spa34 = SPA(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb26 = SPB(2,6);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s35 = -(spa35*spb35);
complex<T> s36 = -(spa36*spb36);
complex<T> s45 = -(spa45*spb45);
complex<T> s46 = -(spa46*spb46);
complex<T> s123 = SS(1,2,3);
complex<T> s124 = SS(1,2,4);
complex<T> t12 = square(spa14*spb16 + spa24*spb26)*square(spa12*spb13 + spa24*spb34) - square(s124)*square(spa24)*square(spb36); 
complex<T> t23 = spa12*T(2); 
complex<T> t44 = -(spa12*spb16); 
complex<T> t45 = spa25*spb12; 
complex<T> t46 = square(spa23); 
complex<T> t47 = square(spb14); 
complex<T> t49 = -(spa45*(spa25*spb12 + spa35*spb13)); 
complex<T> t54 = -(spa24*spa35); 
complex<T> t59 = -(spb13*spb46); 
complex<T> t62 = -s12 + s34 - s56; 
complex<T> t63 = spa14*spb13 + spa24*spb23; 
complex<T> t64 = -(spa13*spb14) - spa23*spb24; 
complex<T> t65 = spa45*spb14 - spa56*spb16; 
complex<T> t66 = -(spa15*spb13) - spa25*spb23; 
complex<T> t67 = -(spa35*spb45) - spa36*spb46; 
complex<T> t68 = spa24*spb12 + spa34*spb13; 
complex<T> t70 = -(spa14*spb16) - spa24*spb26; 
complex<T> t71 = -(spa12*spb13) - spa24*spb34; 
complex<T> t72 = spa23*spb36 + spa25*spb56; 
complex<T> t73 = s12 - s124; 
complex<T> t79 = square(spa35); 
complex<T> t81 = square(spb46); 
complex<T> t85 = -s123 + s23; 
complex<T> t94 = -(s12*T(2)); 
complex<T> t95 = square(s34); 
complex<T> t96 = square(s56); 
complex<T> t97 = spa56*spb12; 
complex<T> t100 = spa12*spb56; 
complex<T> t103 = spa45*spb35; 
complex<T> t108 = -(spa14*spb13); 
complex<T> t110 = -(spa24*spb23); 
complex<T> t136 = s34*s56; 
complex<T> t141 = spa46*spb36; 
complex<T> d33 = T(2); d33 = T(1)/d33;
complex<T> d34 = s123; d34 = T(1)/d34;
complex<T> d35 = s124; d35 = T(1)/d35;
complex<T> t11 = -(square(s123)*square(spa45)*square(spb13)) + square(spa15*spb13 + spa25*spb23)*square(t68); 
complex<T> t21 = square(t63); 
complex<T> t22 = square(t103 + t141); 
complex<T> t26 = spa25*t62 - spa12*spa56*spb16*T(2); 
complex<T> t27 = -(s13*spa45) + s14*spa45 - s23*spa45 + s24*spa45 + spa35*t108*T(3) + spb23*t54*T(3); 
complex<T> t28 = (s35 + s36 - s45 - s46)*spa24 + spa23*(t103 + t141)*T(3); 
complex<T> t29 = (-s35 - s36 + s45 + s46)*spb13 + spb14*(t103 + t141)*T(3); 
complex<T> t31 = s13*spb36 - s14*spb36 + s23*spb36 - s24*spb36 + spb46*t110*T(3) + spa14*t59*T(3); 
complex<T> t43 = square(square(t63)); 
complex<T> t50 = (spa12*spb16 - spa24*spb46)*t70; 
complex<T> t51 = t67*t72; 
complex<T> t60 = spa24*spb46 + t44; 
complex<T> t61 = spa35*spb13 + t45; 
complex<T> t92 = spa56*spb16*t23 - spa25*t62; 
complex<T> t102 = t65*t67; 
complex<T> t198 = spb14*t66; 
complex<T> t203 = spb56*t45; 
complex<T> d1 = spa34*spa56*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d1 = T(1)/d1;
complex<T> d5 = t63*square(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d5 = T(1)/d5;
complex<T> d6 = spb34*spb56*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d6 = T(1)/d6;
complex<T> d10 = spb13*t63*t97*square(s123 - s56); d10 = T(1)/d10;
complex<T> d13 = spa24*t100*t63*square(s124 - s56); d13 = T(1)/d13;
complex<T> d15 = spa12*spa34*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d15 = T(1)/d15;
complex<T> d16 = spb12*spb34*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d16 = T(1)/d16;
complex<T> d20 = (t103 + t141)*square(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d20 = T(1)/d20;
complex<T> d27 = s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2); d27 = T(1)/d27;
complex<T> d28 = s123*t63*T(2); d28 = T(1)/d28;
complex<T> d29 = s124*t63*T(2); d29 = T(1)/d29;
complex<T> d30 = t63*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d30 = T(1)/d30;
complex<T> t18 = d5*(spb16*t62 - t203*T(2)); 
complex<T> t19 = t61*t64; 
complex<T> t20 = spb13*t60; 
complex<T> t25 = square(t61); 
complex<T> t30 = square(t60); 
complex<T> t32 = -(spb16*t62) + t203*T(2); 
complex<T> t101 = -(spb36*t60); 
complex<T> t106 = d10*s123; 
complex<T> t140 = d13*s124; 
complex<T> t143 = d5*spb12; 
complex<T> t146 = d20*spb56; 
complex<T> t211 = t60*t64; 
complex<T> d2 = spa34*spa56*t21; d2 = T(1)/d2;
complex<T> d3 = (s12 - s123)*spa56*t21; d3 = T(1)/d3;
complex<T> d4 = spa34*spa56*t21*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d4 = T(1)/d4;
complex<T> d7 = spb56*t21*t73; d7 = T(1)/d7;
complex<T> d8 = spb34*spb56*t21; d8 = T(1)/d8;
complex<T> d9 = spb34*spb56*t21*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d9 = T(1)/d9;
complex<T> d11 = (-s123 + s56)*spb13*t21*t97; d11 = T(1)/d11;
complex<T> d12 = (-s124 + s56)*spa24*t100*t21; d12 = T(1)/d12;
complex<T> d14 = t21*t97; d14 = T(1)/d14;
complex<T> d17 = spb12*spb34*t22; d17 = T(1)/d17;
complex<T> d18 = spb12*spb34*t22*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d18 = T(1)/d18;
complex<T> d19 = t100*t21; d19 = T(1)/d19;
complex<T> d21 = spa12*spa34*t22; d21 = T(1)/d21;
complex<T> d22 = spa12*spa34*t22*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d22 = T(1)/d22;
complex<T> d23 = spa56*t43*T(4); d23 = T(1)/d23;
complex<T> d24 = spb56*t43*T(4); d24 = T(1)/d24;
complex<T> d25 = t43*t97*T(2); d25 = T(1)/d25;
complex<T> d26 = t21*(s34*t94 + s56*t94 + t95 + t96 + square(s12) - t136*T(2)); d26 = T(1)/d26;
complex<T> d31 = t43*t97*T(4); d31 = T(1)/d31;
complex<T> d32 = t100*t43*T(4); d32 = T(1)/d32;
complex<T> d36 = spb56*t23*t43; d36 = T(1)/d36;
complex<T> d37 = spa12*t43*T(4); d37 = T(1)/d37;
complex<T> d38 = spb12*t43*T(4); d38 = T(1)/d38;
complex<T> t3 = d23*spa12*t11 - d24*spb12*t12; 
complex<T> t4 = d38*spb56*t11 + d37*spa56*t12; 
complex<T> t5 = d25*t11; 
complex<T> t6 = d32*t12; 
complex<T> t7 = -(d26*(s13 - s14 + s23 - s24)*spa24*spb36*t19) - d29*spa35*spb14*t60 - d28*spa23*spb46*t61 + d26*(s13 - s14 + s23 - s24)*spa45*t20*t64 + d27*spa23*spa35*spb14*spb46*T(2) - d5*t19*(s34*(spa25*spb56 - t44) + s12*(spa25*spb56 + t44) - s56*(spa25*spb56 + t44))*t62*T(3) + d5*t211*(s12*spa56*spb16 + s34*spa56*spb16 - s56*spa56*spb16 - s12*t45 + s34*t45 + s56*t45)*t62*T(3) - d30*(s12*spa25*spb16 + s56*spa25*spb16 - spa23*spa45*spb14*spb36 + spb13*spb46*t54)*t64*T(3); 
complex<T> t55 = -(d11*t68); 
complex<T> t105 = d12*t71; 
complex<T> t111 = d31*t11; 
complex<T> t112 = d36*t12; 
complex<T> t117 = d20*t32; 
complex<T> t119 = d3*t66; 
complex<T> t162 = spa12*t18; 
complex<T> t172 = t102*t146; 
complex<T> t2 = d35*s12*s56*t112 + d34*s12*s56*t5 + d33*(t112 + t5)*t62; 
complex<T> t8 = -(d4*spa24*t19*t27) + d8*spb14*spb36*t60 + d2*spa23*spa45*t61 - spa23*t119*t61 - d9*t20*t31*t64 + d7*spb14*t60*t70 + d1*spa23*spb14*t79 + d6*spa23*spb14*t81 - t162*t19*T(6) + t143*t211*t26*T(6); 
complex<T> t151 = t101*t140*t46 + d7*spb14*t50 + spa23*t105*t50; 
complex<T> t153 = t106*t47*t49 + spa23*t119*t61 + t198*t55*t61; 
complex<T> t182 = s34*t111; 
complex<T> t183 = spa56*t117; 
complex<T> t197 = s34*t6; 
complex<T> t212 = t172*t92; 
complex<T> t1 = t182 + t197; 
complex<T> t196 = t183*t51; 
complex<T> t9 = d8*spb14*t101 - d14*t25 + d4*spa24*t19*t27 + d18*spb36*t102*t29 - d19*t30 + d15*spa35*spb46*t46 + d16*spa35*spb46*t47 + d2*spa23*t49 + d22*spa45*t28*t51 + d9*t20*t31*t64 + d17*spb13*spb46*t65 + d21*spa24*spa35*t72 - d1*spa23*spb14*t79 - d6*spa23*spb14*t81 + t162*t19*T(6) - t196*T(6) + t212*T(6) - t143*t211*t26*T(6); 
complex<T> t10 = d14*t25 - d18*spb36*t102*t29 + d19*t30 - d15*spa35*spb46*t46 - d16*spa35*spb46*t47 - d22*spa45*t28*t51 + spb36*t140*t46*t60 + spa45*t106*t47*t61 + d17*t59*t65 + d11*t198*t61*t68 + spa23*t105*t60*t70 + d21*t54*t72 + t196*T(6) - t212*T(6); 
complex<T> co1 = -(t5*t85); 
complex<T> co2 = s23*t5; 
complex<T> co3 = -(t112*t73); 
complex<T> co4 = s124*t197; 
complex<T> co5 = s123*t182; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t8*(*CI_users[0]->get_value(mc,ind,mu)) + t153*(*CI_users[1]->get_value(mc,ind,mu)) + t151*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t10*(*CI_users[4]->get_value(mc,ind,mu)) + t3*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t7*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu)) + t1*(*CI_users[9]->get_value(mc,ind,mu)) + t2*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t4*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)) + co5*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmmmemep_AX_wCI::\
C2q2g2l_qpqmmmemep_AX_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c3456));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c563, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c4, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmmmemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, m, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmmmemep AX");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:54:30 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spb34 = SPB(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> spb36 = SPB(3,6);
complex<T> spb13 = SPB(1,3);
complex<T> spb46 = SPB(4,6);
complex<T> spb56 = SPB(5,6);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb24 = SPB(2,4);
complex<T> s23 = S(2,3);
complex<T> s12 = S(1,2);
complex<T> s123 = SS(1,2,3);
complex<T> s56 = -(spa56*spb56);
complex<T> s124 = SS(1,2,4);
complex<T> s34 = -(spa34*spb34);
complex<T> t5 = square(spb34); 
complex<T> t8 = spb14*spb36 + spb13*spb46; 
complex<T> t9 = -s23 - s34; 
complex<T> t10 = cube(spb34); 
complex<T> t11 = s12 - s124; 
complex<T> t14 = square(spb16); 
complex<T> t17 = -s123 + s23; 
complex<T> t18 = s123*s124 - s12*s56; 
complex<T> t19 = -(spa45*spb14); 
complex<T> t29 = spa24*spb16; 
complex<T> t36 = spa35*spb16; 
complex<T> d4 = spb12*spb34*square(s123 - s56); d4 = T(1)/d4;
complex<T> d5 = spb13*spb34*square(s123 - s56); d5 = T(1)/d5;
complex<T> d8 = spb12*spb34*square(s124 - s56); d8 = T(1)/d8;
complex<T> d9 = spb24*spb34*square(s124 - s56); d9 = T(1)/d9;
complex<T> t6 = -t36; 
complex<T> t20 = -t29; 
complex<T> t22 = d9*t9; 
complex<T> t30 = d4*spa34; 
complex<T> t33 = d8*spb13; 
complex<T> d1 = (s12 - s123)*spb56*t5; d1 = T(1)/d1;
complex<T> d2 = spb56*t11*t5; d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spb12*t5; d3 = T(1)/d3;
complex<T> d6 = spb12*spb56*t5; d6 = T(1)/d6;
complex<T> d7 = (-s124 + s56)*spb12*t5; d7 = T(1)/d7;
complex<T> d10 = spb12*spb56*t10*T(2); d10 = T(1)/d10;
complex<T> d11 = spb12*t10*T(2); d11 = T(1)/d11;
complex<T> d12 = spb12*spb56*t10*T(4); d12 = T(1)/d12;
complex<T> t23 = -(d1*spa23); 
complex<T> t25 = d7*spb13; 
complex<T> t27 = d3*spa45*spb14*spb16 + d1*spa23*spb16*spb36 + d6*t14 + d5*t19*t29 + spb16*t19*t30; 
complex<T> t31 = d2*spb46; 
complex<T> t37 = spa34*t33; 
complex<T> t4 = d3*spb16*t19 + d5*spa45*spb14*t29 + spa45*spb14*spb16*t30 + t22*t36 + t25*t36 + t37*t6; 
complex<T> t35 = spb16*spb36*t23 + t29*t31; 
complex<T> t39 = -(d6*t14) + t20*t31 + t36*t37 + t22*t6 + t25*t6; 
complex<T> co1 = -(d10*spb16*t17*t8); 
complex<T> co2 = d10*s23*spb16*t8; 
complex<T> co3 = -(d10*spb16*t11*t8); 
complex<T> co4 = d11*spa56*spb16*t8; 
complex<T> co5 = d12*spb16*t18*t8; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t35*(*CI_users[0]->get_value(mc,ind,mu)) + t27*(*CI_users[1]->get_value(mc,ind,mu)) + t39*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)) + co4*(*CI_users[7]->get_value(mc,ind,mu)) + co5*(*CI_users[8]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmppemep_AXSL_wCI::\
C2q2g2l_qpqmppemep_AXSL_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmppemep_AXSL_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, p, em, ep}, AXSL}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmppemep AXSL");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:54:33 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb46 = SPB(4,6);
complex<T> spa14 = SPA(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb36 = SPB(3,6);
complex<T> s56 = S(5,6);
complex<T> s123 = SS(1,2,3);
complex<T> s124 = SS(1,2,4);
complex<T> t6 = -(spa12*spb14) + spa23*spb34; 
complex<T> t7 = -(spa12*spb13) - spa24*spb34; 
complex<T> t10 = spa25*spb36; 
complex<T> t13 = spa25*spb46; 
complex<T> d1 = spa13*spa23*square(s123 - s56); d1 = T(1)/d1;
complex<T> d2 = spa14*spa24*square(s124 - s56); d2 = T(1)/d2;
complex<T> t11 = d1*t6; 
complex<T> t12 = d2*t7; 
complex<T> t3 = -(t10*t12) - t11*t13; 
complex<T> co1 = t11*t13; 
complex<T> co2 = t10*t12; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmpmemep_AXSL_wCI::\
C2q2g2l_qpqmpmemep_AXSL_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmpmemep_AXSL_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, p, m, em, ep}, AXSL}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmpmemep AXSL");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:54:36 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spb36 = SPB(3,6);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> s56 = S(5,6);
complex<T> s123 = SS(1,2,3);
complex<T> s124 = SS(1,2,4);
complex<T> t6 = spa25*spb12 + spa45*spb14; 
complex<T> t7 = -(spa12*spb16) + spa23*spb36; 
complex<T> t10 = spa24*spa45; 
complex<T> t11 = spb13*spb36; 
complex<T> d1 = spa13*spa23*square(s123 - s56); d1 = T(1)/d1;
complex<T> d2 = spb14*spb24*square(s124 - s56); d2 = T(1)/d2;
complex<T> t12 = d2*t6; 
complex<T> t13 = d1*t7; 
complex<T> t3 = -(t11*t12) - t10*t13; 
complex<T> co1 = t10*t13; 
complex<T> co2 = t11*t12; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmmpemep_AXSL_wCI::\
C2q2g2l_qpqmmpemep_AXSL_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmmpemep_AXSL_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, p, em, ep}, AXSL}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmmpemep AXSL");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:54:39 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb13 = SPB(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb46 = SPB(4,6);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> s56 = S(5,6);
complex<T> s123 = SS(1,2,3);
complex<T> s124 = SS(1,2,4);
complex<T> t6 = spa25*spb12 + spa35*spb13; 
complex<T> t7 = -(spa12*spb16) + spa24*spb46; 
complex<T> t10 = spa23*spa35; 
complex<T> t11 = spb14*spb46; 
complex<T> d1 = spb13*spb23*square(s123 - s56); d1 = T(1)/d1;
complex<T> d2 = spa14*spa24*square(s124 - s56); d2 = T(1)/d2;
complex<T> t12 = d1*t6; 
complex<T> t13 = d2*t7; 
complex<T> t3 = -(t11*t12) - t10*t13; 
complex<T> co1 = t11*t12; 
complex<T> co2 = t10*t13; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g2l_qpqmmmemep_AXSL_wCI::\
C2q2g2l_qpqmmmemep_AXSL_wCI
      (const std::vector<int>& ind){
	//define corner vectors -- assume ordered
	 vector<int> c1;  c1.push_back(ind.at(0));
	 vector<int> c2;  c2.push_back(ind.at(1));
	 vector<int> c3;  c3.push_back(ind.at(2));
	 vector<int> c4;  c4.push_back(ind.at(3));
	 vector<int> c5;  c5.push_back(ind.at(4));
	 vector<int> c6;  c6.push_back(ind.at(5));

	 vector<int> c12;  c12.push_back(ind.at(0)); c12.push_back(ind.at(1));
	 vector<int> c23;  c23.push_back(ind.at(1)); c23.push_back(ind.at(2));
	 vector<int> c34;  c34.push_back(ind.at(2)); c34.push_back(ind.at(3));
	 vector<int> c45;  c45.push_back(ind.at(3)); c45.push_back(ind.at(4));
	 vector<int> c56;  c56.push_back(ind.at(4)); c56.push_back(ind.at(5));
	 vector<int> c16;  c16.push_back(ind.at(5)); c16.push_back(ind.at(0));
	 vector<int> c61;  c61.push_back(ind.at(5)); c61.push_back(ind.at(0));
	 vector<int> c41;  c41.push_back(ind.at(3)); c41.push_back(ind.at(0));
	 vector<int> c14 = c41; 
	 vector<int> c51;  c51.push_back(ind.at(4)); c51.push_back(ind.at(0));
	 vector<int> c15 = c51; 

	 vector<int> c123;  for(int i = 1; i<=3; i++) {c123.push_back(ind.at(i-1));}
	 vector<int> c234;  for(int i = 2; i<=4; i++) {c234.push_back(ind.at(i-1));}
	 vector<int> c345;  for(int i = 3; i<=5; i++) {c345.push_back(ind.at(i-1));}
	 vector<int> c456;  for(int i = 4; i<=6; i++) {c456.push_back(ind.at(i-1));}
	 vector<int> c356;  c356.push_back(ind.at(2));
	                    for(int i = 5; i<=6; i++) {c356.push_back(ind.at(i-1));}
	 vector<int> c563 = c356;
	 vector<int> c561;  for(int i = 5; i<=6; i++) {c561.push_back(ind.at(i-1));}
	                      c561.push_back(ind.at(0));
	 vector<int> c156;for(int i = 5; i<=6; i++) {c156.push_back(ind.at(i-1));}
	                      c156.push_back(ind.at(0));
	 vector<int> c256;for(int i = 5; i<=6; i++) {c256.push_back(ind.at(i-1));}
	                      c256.push_back(ind.at(1));

	 vector<int> c126;  c126.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c126.push_back(ind.at(i-1));}
	 vector<int> c612;  c612.push_back(ind.at(5)) ;
	                    for(int i = 1; i<=2; i++) {c612.push_back(ind.at(i-1));}
	 vector<int> c124;  c124.push_back(ind.at(3)) ;
	                    for(int i = 1; i<=2; i++) {c124.push_back(ind.at(i-1));}
	 vector<int> c134;  c134.push_back(ind.at(0)) ;
	                    for(int i = 3; i<=4; i++) {c134.push_back(ind.at(i-1));}

	 vector<int> c1234;  for(int i = 1; i<=4; i++) {c1234.push_back(ind.at(i-1));}
	 vector<int> c2345;  for(int i = 2; i<=5; i++) {c2345.push_back(ind.at(i-1));}
	 vector<int> c3456;  for(int i = 3; i<=6; i++) {c3456.push_back(ind.at(i-1));}
	 vector<int> c1456;  for(int i = 4; i<=6; i++) {c1456.push_back(ind.at(i-1));}
	                     c1456.push_back(ind.at(0));
	 vector<int> c4561 = c1456;
	 vector<int> c1256;  for(int i = 5; i<=6; i++) {c1256.push_back(ind.at(i-1));}
	                     c1256.push_back(ind.at(0)); c1256.push_back(ind.at(1));
	 vector<int> c5612 = c1256;
	 vector<int> c1236;  c1236.push_back(ind.at(5));
	                     for(int i = 1; i<=3; i++) {c1236.push_back(ind.at(i-1));}
	 vector<int> c2356;  for(int i = 5; i<=6; i++) {c2356.push_back(ind.at(i-1));}
	                     c2356.push_back(ind.at(1)); c2356.push_back(ind.at(2));
CI_users.push_back(new Cached_Bubble_Integral_User(c123, c456));
CI_users.push_back(new Cached_Bubble_Integral_User(c124, c356));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g2l_qpqmmmemep_AXSL_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, m, m, em, ep}, AXSL}
 
#if _VERBOSE
  _MESSAGE("C2q2g2l :  qpqmmmemep AXSL");
#endif
 
//#define TimeStamp "Fri 25 Sep 2009 23:54:42 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spb16 = SPB(1,6);
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> s56 = S(5,6);
complex<T> s123 = SS(1,2,3);
complex<T> s124 = SS(1,2,4);
complex<T> t6 = spa24*spb12 + spa34*spb13; 
complex<T> t7 = spa23*spb12 - spa34*spb14; 
complex<T> t10 = spa35*spb16; 
complex<T> t13 = spa45*spb16; 
complex<T> d1 = spb13*spb23*square(s123 - s56); d1 = T(1)/d1;
complex<T> d2 = spb14*spb24*square(s124 - s56); d2 = T(1)/d2;
complex<T> t11 = d1*t6; 
complex<T> t12 = d2*t7; 
complex<T> t3 = -(t10*t12) - t11*t13; 
complex<T> co1 = t11*t13; 
complex<T> co2 = t10*t12; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qpppqmemep_L C2q2g2l_44408_L
#define _C_qppmqmemep_L C2q2g2l_44300_L
#define _C_qpmpqmemep_L C2q2g2l_44390_L
#define _C_qpmmqmemep_L C2q2g2l_44282_L
#define _C_qmppqpemep_L C2q2g2l_44623_L
#define _C_qmpmqpemep_L C2q2g2l_44515_L
#define _C_qmmpqpemep_L C2q2g2l_44605_L
#define _C_qmmmqpemep_L C2q2g2l_44497_L
#define _C_qppqmpemep_SLC C2q2g2l_44768_SLC
#define _C_qppqmmemep_SLC C2q2g2l_44120_SLC
#define _C_qpmqmpemep_SLC C2q2g2l_44750_SLC
#define _C_qpmqmmemep_SLC C2q2g2l_44102_SLC
#define _C_qpqmppemep_SLC C2q2g2l_44828_SLC
#define _C_qpqmpmemep_SLC C2q2g2l_44180_SLC
#define _C_qpqmmpemep_SLC C2q2g2l_44720_SLC
#define _C_qpqmmmemep_SLC C2q2g2l_44072_SLC
#define _C_qpppqmemep_nf C2q2g2l_44408_nf
#define _C_qppmqmemep_nf C2q2g2l_44300_nf
#define _C_qpmpqmemep_nf C2q2g2l_44390_nf
#define _C_qpmmqmemep_nf C2q2g2l_44282_nf
#define _C_qpppqmemep_nf_top C2q2g2l_44408_nf_top
#define _C_qppmqmemep_nf_top C2q2g2l_44300_nf_top
#define _C_qpmpqmemep_nf_top C2q2g2l_44390_nf_top
#define _C_qpmmqmemep_nf_top C2q2g2l_44282_nf_top
#define _C_qpqmppemep_VECT C2q2g2l_44828_VECT
#define _C_qpqmpmemep_VECT C2q2g2l_44180_VECT
#define _C_qpqmmpemep_VECT C2q2g2l_44720_VECT
#define _C_qpqmmmemep_VECT C2q2g2l_44072_VECT
#define _C_qpqmppemep_AX C2q2g2l_44828_AX
#define _C_qpqmpmemep_AX C2q2g2l_44180_AX
#define _C_qpqmmpemep_AX C2q2g2l_44720_AX
#define _C_qpqmmmemep_AX C2q2g2l_44072_AX
#define _C_qpqmppemep_AXSL C2q2g2l_44828_AXSL
#define _C_qpqmpmemep_AXSL C2q2g2l_44180_AXSL
#define _C_qpqmmpemep_AXSL C2q2g2l_44720_AXSL
#define _C_qpqmmmemep_AXSL C2q2g2l_44072_AXSL
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qpppqmemep_L case 44408
 
#define _CASE_qppmqmemep_L case 44300
 
#define _CASE_qpmpqmemep_L case 44390
 
#define _CASE_qpmmqmemep_L case 44282
 
#define _CASE_qmppqpemep_L case 44623
 
#define _CASE_qmpmqpemep_L case 44515
 
#define _CASE_qmmpqpemep_L case 44605
 
#define _CASE_qmmmqpemep_L case 44497
 
#define _CASE_qppqmpemep_SLC case 44768
 
#define _CASE_qppqmmemep_SLC case 44120
 
#define _CASE_qpmqmpemep_SLC case 44750
 
#define _CASE_qpmqmmemep_SLC case 44102
 
#define _CASE_qpqmppemep_SLC case 44828
 
#define _CASE_qpqmpmemep_SLC case 44180
 
#define _CASE_qpqmmpemep_SLC case 44720
 
#define _CASE_qpqmmmemep_SLC case 44072
 
#define _CASE_qpppqmemep_nf case 44408
 
#define _CASE_qppmqmemep_nf case 44300
 
#define _CASE_qpmpqmemep_nf case 44390
 
#define _CASE_qpmmqmemep_nf case 44282
 
#define _CASE_qpppqmemep_nf_top case 44408
 
#define _CASE_qppmqmemep_nf_top case 44300
 
#define _CASE_qpmpqmemep_nf_top case 44390
 
#define _CASE_qpmmqmemep_nf_top case 44282
 
#define _CASE_qpqmppemep_VECT case 44828
 
#define _CASE_qpqmpmemep_VECT case 44180
 
#define _CASE_qpqmmpemep_VECT case 44720
 
#define _CASE_qpqmmmemep_VECT case 44072
 
#define _CASE_qpqmppemep_AX case 44828
 
#define _CASE_qpqmpmemep_AX case 44180
 
#define _CASE_qpqmmpemep_AX case 44720
 
#define _CASE_qpqmmmemep_AX case 44072
 
#define _CASE_qpqmppemep_AXSL case 44828
 
#define _CASE_qpqmpmemep_AXSL case 44180
 
#define _CASE_qpqmmpemep_AXSL case 44720
 
#define _CASE_qpqmmmemep_AXSL case 44072
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2g2l_AX( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qpqmppemep_AX: return new 
                       C2q2g2l_qpqmppemep_AX_wCI(ind);
    _CASE_qpqmpmemep_AX: return new 
                       C2q2g2l_qpqmpmemep_AX_wCI(ind);
    _CASE_qpqmmpemep_AX: return new 
                       C2q2g2l_qpqmmpemep_AX_wCI(ind);
    _CASE_qpqmmmemep_AX: return new 
                       C2q2g2l_qpqmmmemep_AX_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g2l_AXSL( int hc,const std::vector<int>& ind) 
{ 
    switch (hc) {
    _CASE_qpqmppemep_AXSL: return new 
                       C2q2g2l_qpqmppemep_AXSL_wCI(ind);
    _CASE_qpqmpmemep_AXSL: return new 
                       C2q2g2l_qpqmpmemep_AXSL_wCI(ind);
    _CASE_qpqmmpemep_AXSL: return new 
                       C2q2g2l_qpqmmpemep_AXSL_wCI(ind);
    _CASE_qpqmmmemep_AXSL: return new 
                       C2q2g2l_qpqmmmemep_AXSL_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g2l_L( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qpppqmemep_L: return new 
                       C2q2g2l_qpppqmemep_L_wCI(ind);
    _CASE_qppmqmemep_L: return new 
                       C2q2g2l_qppmqmemep_L_wCI(ind);
    _CASE_qpmpqmemep_L: return new 
                       C2q2g2l_qpmpqmemep_L_wCI(ind);
    _CASE_qpmmqmemep_L: return new 
                       C2q2g2l_qpmmqmemep_L_wCI(ind);
    _CASE_qmppqpemep_L: return new 
                       C2q2g2l_qmppqpemep_L_wCI(ind);
    _CASE_qmpmqpemep_L: return new 
                       C2q2g2l_qmpmqpemep_L_wCI(ind);
    _CASE_qmmpqpemep_L: return new 
                       C2q2g2l_qmmpqpemep_L_wCI(ind);
    _CASE_qmmmqpemep_L: return new 
                       C2q2g2l_qmmmqpemep_L_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g2l_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qpppqmemep_nf: return new 
                       C2q2g2l_qpppqmemep_nf_wCI(ind);
    _CASE_qppmqmemep_nf: return new 
                       C2q2g2l_qppmqmemep_nf_wCI(ind);
    _CASE_qpmpqmemep_nf: return new 
                       C2q2g2l_qpmpqmemep_nf_wCI(ind);
    _CASE_qpmmqmemep_nf: return new 
                       C2q2g2l_qpmmqmemep_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g2l_nf_top( int hc,const std::vector<int>& ind) 
{ 
    switch (hc) {
    _CASE_qpppqmemep_nf_top: return new 
                       C2q2g2l_qpppqmemep_nf_top_wCI(ind);
    _CASE_qppmqmemep_nf_top: return new 
                       C2q2g2l_qppmqmemep_nf_top_wCI(ind);
    _CASE_qpmpqmemep_nf_top: return new 
                       C2q2g2l_qpmpqmemep_nf_top_wCI(ind);
    _CASE_qpmmqmemep_nf_top: return new 
                       C2q2g2l_qpmmqmemep_nf_top_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g2l_SLC( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_qppqmpemep_SLC: return new 
                       C2q2g2l_qppqmpemep_SLC_wCI(ind);
    _CASE_qppqmmemep_SLC: return new 
                       C2q2g2l_qppqmmemep_SLC_wCI(ind);
    _CASE_qpmqmpemep_SLC: return new 
                       C2q2g2l_qpmqmpemep_SLC_wCI(ind);
    _CASE_qpmqmmemep_SLC: return new 
                       C2q2g2l_qpmqmmemep_SLC_wCI(ind);
    _CASE_qpqmppemep_SLC: return new 
                       C2q2g2l_qpqmppemep_SLC_wCI(ind);
    _CASE_qpqmpmemep_SLC: return new 
                       C2q2g2l_qpqmpmemep_SLC_wCI(ind);
    _CASE_qpqmmpemep_SLC: return new 
                       C2q2g2l_qpqmmpemep_SLC_wCI(ind);
    _CASE_qpqmmmemep_SLC: return new 
                       C2q2g2l_qpqmmmemep_SLC_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g2l_VECT( int hc,const std::vector<int>& ind) 
{ 
    switch (hc) {
    _CASE_qpqmppemep_VECT: return new 
                       C2q2g2l_qpqmppemep_VECT_wCI(ind);
    _CASE_qpqmpmemep_VECT: return new 
                       C2q2g2l_qpqmpmemep_VECT_wCI(ind);
    _CASE_qpqmmpemep_VECT: return new 
                       C2q2g2l_qpqmmpemep_VECT_wCI(ind);
    _CASE_qpqmmmemep_VECT: return new 
                       C2q2g2l_qpqmmmemep_VECT_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
