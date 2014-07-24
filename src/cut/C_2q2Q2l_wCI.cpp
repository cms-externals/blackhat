/*
*C_2q2Q2l_wCI.cpp
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
 


GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpQmQpqmemep_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpQpQmqmemep_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmQpQmqpemep_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmQmQpqpemep_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpQmQpqmemep_nf_top_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpQpQmqmemep_nf_top_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmQpQmqpemep_nf_top_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmQmQpqpemep_nf_top_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpQmQpqmemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpQpQmqmemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmQpQmqpemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmQmQpqpemep_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpqmQpQmemep_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpqmQmQpemep_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmqpQpQmemep_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmqpQmQpemep_sl_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpqmQpQmemep_AX_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmqpQpQmemep_AX_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qpqmQmQpemep_AX_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q2Q2l_qmqpQmQpemep_AX_wCI)

C2q2Q2l_qpQmQpqmemep_nf_wCI::C2q2Q2l_qpQmQpqmemep_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qpQmQpqmemep_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qm, Qp, qm, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpQmQpqmemep nf");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:21:55 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spa23 = SPA(2,3);
complex<T> spb36 = SPB(3,6);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s56 = S(5,6);
complex<T> s123 = SS(1,2,3);
complex<T> s234 = SS(2,3,4);
complex<T> d1 = s23*s234*s56; d1 = T(1)/d1;
complex<T> d2 = s123*s23*s56; d2 = T(1)/d2;
complex<T> d3 = T(3); d3 = T(1)/d3;
complex<T> t1 = d1*spa24*spb16*(spa25*spb23 - spa45*spb34) + d2*spa45*spb13*(spa12*spb16 - spa23*spb36); 
complex<T> co1 = Complex(0,2)*d3*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q2l_qpQpQmqmemep_nf_wCI::C2q2Q2l_qpQpQmqmemep_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qpQpQmqmemep_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qp, Qm, qm, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpQpQmqmemep nf");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:21:56 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spb16 = SPB(1,6);
complex<T> spa23 = SPA(2,3);
complex<T> spb26 = SPB(2,6);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s56 = S(5,6);
complex<T> s123 = SS(1,2,3);
complex<T> s234 = SS(2,3,4);
complex<T> d1 = s23*s234*s56; d1 = T(1)/d1;
complex<T> d2 = s123*s23*s56; d2 = T(1)/d2;
complex<T> d3 = T(3); d3 = T(1)/d3;
complex<T> t1 = -(d1*spa34*spb16*(spa35*spb23 + spa45*spb24)) + d2*spa45*spb12*(spa13*spb16 + spa23*spb26); 
complex<T> co1 = Complex(0,-2)*d3*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q2l_qmQpQmqpemep_nf_wCI::C2q2Q2l_qmQpQmqpemep_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qmQpQmqpemep_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmQpQmqpemep nf");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:21:57 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb46 = SPB(4,6);
complex<T> spb24 = SPB(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb26 = SPB(2,6);
complex<T> spa34 = SPA(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s56 = S(5,6);
complex<T> s123 = SS(1,2,3);
complex<T> s234 = SS(2,3,4);
complex<T> d1 = s123*s23*s56; d1 = T(1)/d1;
complex<T> d2 = s23*s234*s56; d2 = T(1)/d2;
complex<T> d3 = T(3); d3 = T(1)/d3;
complex<T> t1 = d1*spa13*(spa15*spb12 - spa35*spb23)*spb46 + d2*spa15*spb24*(spa23*spb26 - spa34*spb46); 
complex<T> co1 = Complex(0,2)*d3*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q2l_qmQmQpqpemep_nf_wCI::C2q2Q2l_qmQmQpqpemep_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qmQmQpqpemep_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmQmQpqpemep nf");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:21:57 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb46 = SPB(4,6);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb36 = SPB(3,6);
complex<T> spa24 = SPA(2,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s56 = S(5,6);
complex<T> s123 = SS(1,2,3);
complex<T> s234 = SS(2,3,4);
complex<T> d1 = s123*s23*s56; d1 = T(1)/d1;
complex<T> d2 = s23*s234*s56; d2 = T(1)/d2;
complex<T> d3 = T(3); d3 = T(1)/d3;
complex<T> t1 = d1*spa12*(spa15*spb13 + spa25*spb23)*spb46 - d2*spa15*spb34*(spa23*spb36 + spa24*spb46); 
complex<T> co1 = Complex(0,-2)*d3*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C2q2Q2l_qpQmQpqmemep_nf_top_wCI::C2q2Q2l_qpQmQpqmemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2Q2l_qpQmQpqmemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpQmQpqmemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2Q2l_qpQpQmqmemep_nf_top_wCI::C2q2Q2l_qpQpQmqmemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2Q2l_qpQpQmqmemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpQpQmqmemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2Q2l_qmQpQmqpemep_nf_top_wCI::C2q2Q2l_qmQpQmqpemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2Q2l_qmQpQmqpemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmQpQmqpemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2Q2l_qmQmQpqpemep_nf_top_wCI::C2q2Q2l_qmQmQpqpemep_nf_top_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2Q2l_qmQmQpqpemep_nf_top_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmQmQpqpemep nf_top");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2Q2l_qpQmQpqmemep_L_wCI::C2q2Q2l_qpQmQpqmemep_L_wCI
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
     C2q2Q2l_qpQmQpqmemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qm, Qp, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpQmQpqmemep L");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:23:05 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
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
complex<T> spb14 = SPB(1,4);
complex<T> spb36 = SPB(3,6);
complex<T> spb46 = SPB(4,6);
complex<T> spa25 = SPA(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa35 = SPA(3,5);
complex<T> s123 = SS(1,2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s124 = SS(1,2,4);
complex<T> s234 = SS(2,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s134 = SS(1,3,4);
complex<T> t21 = -(s12*T(2)); 
complex<T> t33 = square(spa24); 
complex<T> t34 = square(spb13); 
complex<T> t37 = square(spa13*spb16 + spa23*spb26); 
complex<T> t43 = square(spa25*spb24 + spa35*spb34); 
complex<T> t65 = square(spa15); 
complex<T> t66 = square(spb46); 
complex<T> t70 = -(spa12*spb14); 
complex<T> t73 = spa14*spb34; 
complex<T> t82 = spa12*spb24 + spa13*spb34; 
complex<T> t83 = -(spa25*spb23) + spa45*spb34; 
complex<T> t84 = -(spa12*spb16) + spa23*spb36; 
complex<T> t86 = spa35*spb23 + spa45*spb24; 
complex<T> t87 = spa12*spb26 + spa13*spb36; 
complex<T> t89 = square(spa45); 
complex<T> t91 = -s123 + s56; 
complex<T> t92 = -s234 + s56; 
complex<T> t93 = spa24*spb12 + spa34*spb13; 
complex<T> t94 = -(spa13*spb16) - spa23*spb26; 
complex<T> t95 = -(spa25*spb24) - spa35*spb34; 
complex<T> t123 = spa13*spb14; 
complex<T> t124 = spa23*spb24; 
complex<T> t131 = s34*s56; 
complex<T> t158 = s123*spa23; 
complex<T> d1 = spa13*spb23 + spa14*spb24; d1 = T(1)/d1;
complex<T> d2 = spa56; d2 = T(1)/d2;
complex<T> d3 = spb56; d3 = T(1)/d3;
complex<T> d15 = s23*s234*s56; d15 = T(1)/d15;
complex<T> d16 = s123*s23*s56; d16 = T(1)/d16;
complex<T> d17 = T(6); d17 = T(1)/d17;
complex<T> d18 = s234*(-s234 + s34)*spa56*spb24*(spa13*spb23 + spa14*spb24); d18 = T(1)/d18;
complex<T> d24 = T(2); d24 = T(1)/d24;
complex<T> d33 = s234; d33 = T(1)/d33;
complex<T> d36 = square(s123)*T(2); d36 = T(1)/d36;
complex<T> d43 = square(s234)*T(2); d43 = T(1)/d43;
complex<T> d44 = s123; d44 = T(1)/d44;
complex<T> t35 = spa56*t82; 
complex<T> t36 = square(t86); 
complex<T> t39 = d1*(spa23*spb13 + spa24*spb14)*t65 + square(spa25); 
complex<T> t41 = -(spa25*spb16) + d1*spa15*(spa23*spb13 + spa24*spb14)*spb26; 
complex<T> t42 = square(spb16) + d1*(spa23*spb13 + spa24*spb14)*square(spb26); 
complex<T> t44 = square(t87); 
complex<T> t67 = spb56*t82; 
complex<T> t88 = s34*t21 + s56*t21 + square(s12) + square(s34) + square(s56) - t131*T(2); 
complex<T> d4 = (t123 + t124)*(s34*t21 + s56*t21 + square(s12) + square(s34) + square(s56) - t131*T(2)); d4 = T(1)/d4;
complex<T> d5 = (s12 - s123)*s123*spa13*spb56*(t123 + t124); d5 = T(1)/d5;
complex<T> d22 = t123 + t124; d22 = T(1)/d22;
complex<T> d23 = (spa13*spb23 + spa14*spb24)*(s34*t21 + s56*t21 + square(s12) + square(s34) + square(s56) - t131*T(2)); d23 = T(1)/d23;
complex<T> d26 = s234*spa23*spb56*t93; d26 = T(1)/d26;
complex<T> d27 = s123*spa56*spb23*t93; d27 = T(1)/d27;
complex<T> d30 = (t123 + t124)*t82*T(2)*(s34*t21 + s56*t21 + square(s12) + square(s34) + square(s56) - t131*T(2)); d30 = T(1)/d30;
complex<T> d31 = s123*t82; d31 = T(1)/d31;
complex<T> d32 = (t123 + t124)*t82*T(2); d32 = T(1)/d32;
complex<T> d34 = spa56*spb23*t93; d34 = T(1)/d34;
complex<T> d38 = s234*t82; d38 = T(1)/d38;
complex<T> d39 = (spa13*spb23 + spa14*spb24)*t82*T(2); d39 = T(1)/d39;
complex<T> d40 = (spa13*spb23 + spa14*spb24)*t82*T(2)*(s34*t21 + s56*t21 + square(s12) + square(s34) + square(s56) - t131*T(2)); d40 = T(1)/d40;
complex<T> d42 = spa23*spb56*t93; d42 = T(1)/d42;
complex<T> t1 = -(d5*spa12*spb13*t37) + d4*(d2*s124*spb12*t39 + d3*s123*spa12*t42 + s12*t41*T(2)); 
complex<T> t30 = s34*(-(spa45*spb36) + d22*spa35*(spa14*spb13 + spa24*spb23)*spb46); 
complex<T> t40 = t89 + d22*(spa14*spb13 + spa24*spb23)*square(spa35); 
complex<T> t45 = d22*(spa14*spb13 + spa24*spb23)*t66 + square(spb36); 
complex<T> t120 = -(d27*t34); 
complex<T> t156 = s234*t35; 
complex<T> d6 = s123*(-s123 + s23)*spa13*t67; d6 = T(1)/d6;
complex<T> d7 = s123*t67*square(s123 - s23)*T(2); d7 = T(1)/d7;
complex<T> d8 = s123*(-s123 + s23)*t67; d8 = T(1)/d8;
complex<T> d9 = spa23*t67*t91; d9 = T(1)/d9;
complex<T> d10 = spa23*t67*square(t91)*T(2); d10 = T(1)/d10;
complex<T> d11 = (t123 + t124)*t67*t91; d11 = T(1)/d11;
complex<T> d19 = spb23*t35*square(t92)*T(2); d19 = T(1)/d19;
complex<T> d20 = (spa13*spb23 + spa14*spb24)*t35*t92; d20 = T(1)/d20;
complex<T> d21 = spb23*t35*t92; d21 = T(1)/d21;
complex<T> d28 = t158*t67*square(spa13); d28 = T(1)/d28;
complex<T> d29 = t158*t67*square(t123 + t124); d29 = T(1)/d29;
complex<T> d35 = spa23*t67; d35 = T(1)/d35;
complex<T> d41 = spb23*t35; d41 = T(1)/d41;
complex<T> t2 = d18*spa24*spb34*t36 - d23*(d2*s234*spb34*t40 + d3*s134*spa34*t45 + t30*T(2)); 
complex<T> t13 = s12*(t120*t89 + d28*t37*square(spa12)); 
complex<T> t20 = t120*t89 + d29*t37*square(spa23*spb34 + t70); 
complex<T> t164 = d21*t83; 
complex<T> t165 = d9*t84; 
complex<T> d12 = (s23 - s234)*spb24*t156; d12 = T(1)/d12;
complex<T> d13 = t156*square(s23 - s234)*T(2); d13 = T(1)/d13;
complex<T> d14 = (s23 - s234)*t156; d14 = T(1)/d14;
complex<T> d25 = spb23*t156*square(spa13*spb23 + spa14*spb24); d25 = T(1)/d25;
complex<T> d37 = spb23*t156*square(spb24); d37 = T(1)/d37;
complex<T> t3 = -(d2*d4*s124*spb12*t39) + d2*d23*s234*spb34*t40 + d4*t21*t41 - d3*d4*s123*spa12*t42 + d23*d3*s134*spa34*t45 + d20*s234*spb13*t65 - d19*s234*t34*t65 + d10*s123*t33*t66 - spa15*spb13*t164*T(2) + d23*t30*T(2) - d16*d24*spa45*spb13*t84*T(3) - spa24*(d11*s123*t66 + spb46*t165*T(2) + d15*d24*spb16*t83*T(3)); 
complex<T> t6 = (-s123 + s23)*(-t20 + t120*t89 + d28*t37*square(spa12)); 
complex<T> t10 = (d44*s12*s56 - d24*(s12 - s34 + s56))*t20 - d40*(s12 + s34 - s56)*(spa35*spb36 + spa45*spb46)*(s234*spa12*spb23 - s234*t73 + s56*t73) + d39*(spa12*spa15*spb16*spb23 + spa12*spa25*spb23*spb26 + spa15*spb16*t73 - spa25*spb26*t73) + d38*spa15*spa24*spb16*spb34*T(2) - d43*(d42*t33*square(spb16) - d41*square(t83))*(s12*s234 - s234*(s34 + s56) + t131*T(2)); 
complex<T> t15 = s23*(d26*t33*square(spb16) - d37*t36*square(spb34)); 
complex<T> t16 = s34*t20; 
complex<T> t17 = d13*spb23*t33*t43 - d7*spa23*t34*t44 + d8*spb13*t84*t87 + d6*spa12*spb13*t87*t94 - d14*spa24*t83*t95 - d12*spa24*spb34*t86*t95 + d15*d17*spa24*spb16*t83*T(13) + d16*d17*spa45*spb13*t84*T(13); 
complex<T> t19 = d26*t33*square(spb16) - d25*t36*square(-(spa12*spb23) + t73); 
complex<T> t27 = -(d18*spa24*spb34*t36) - d13*spb23*t33*t43 - d20*s234*spb13*t65 + d19*s234*t34*t65 + d14*spa24*t83*t95 + d12*spa24*spb34*t86*t95 + spa15*spb13*t164*T(2); 
complex<T> t31 = d5*spa12*spb13*t37 + d7*spa23*t34*t44 + d11*s123*spa24*t66 - d10*s123*t33*t66 - d8*spb13*t84*t87 - d6*spa12*spb13*t87*t94 + spa24*spb46*t165*T(2); 
complex<T> t4 = (s156 - s34)*(t19 - d26*t33*square(spb16) + d37*t36*square(spb34)); 
complex<T> t5 = t15 + s23*t20; 
complex<T> t11 = (d24*(s12 - s34 - s56) + d33*t131)*t19 - d32*(spa23*spa35*spb34*spb36 + spa23*spa45*spb34*spb46 + spa35*spb36*t70 - spa45*spb46*t70) - d30*(s12 + s34 - s56)*(spa15*spb16 + spa25*spb26)*(s56*t70 - s123*(spa23*spb34 + t70)) - d31*spa12*spa45*spb13*spb46*T(2) + d36*(-(d34*t34*t89) + d35*square(t84))*(s123*(-s34 + s56) + s12*(s123 - s56*T(2))); 
complex<T> t127 = d24*t19; 
complex<T> t177 = d24*t16; 
complex<T> t7 = -(s56*(t127 + d24*t20)); 
complex<T> t9 = t177 + s34*(t127 + d15*spa24*spb16*t83 + d16*spa45*spb13*t84); 
complex<T> t159 = s12*t127; 
complex<T> t8 = t13 + t159 + s12*(-(d24*t20) + d15*spa24*spb16*t83 + d16*spa45*spb13*t84); 
complex<T> co1 = s234*t159; 
complex<T> co2 = d24*s34*t15; 
complex<T> co3 = d24*s23*t13; 
complex<T> co4 = s123*t177; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t31*(*CI_users[1]->get_value(mc,ind,mu)) + t17*(*CI_users[2]->get_value(mc,ind,mu)) + t27*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t3*(*CI_users[5]->get_value(mc,ind,mu)) + t8*(*CI_users[6]->get_value(mc,ind,mu)) + t6*(*CI_users[7]->get_value(mc,ind,mu)) + t11*(*CI_users[8]->get_value(mc,ind,mu)) + t5*(*CI_users[9]->get_value(mc,ind,mu)) + t4*(*CI_users[10]->get_value(mc,ind,mu)) + t9*(*CI_users[11]->get_value(mc,ind,mu)) + t10*(*CI_users[12]->get_value(mc,ind,mu)) + t7*(*CI_users[13]->get_value(mc,ind,mu)) + co1*(*CI_users[14]->get_value(mc,ind,mu)) + co2*(*CI_users[15]->get_value(mc,ind,mu)) + co3*(*CI_users[16]->get_value(mc,ind,mu)) + co4*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q2l_qpQpQmqmemep_L_wCI::C2q2Q2l_qpQpQmqmemep_L_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qpQpQmqmemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, Qp, Qm, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpQpQmqmemep L");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:23:12 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spb16 = SPB(1,6);
complex<T> spb26 = SPB(2,6);
complex<T> spa12 = SPA(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb46 = SPB(4,6);
complex<T> spb36 = SPB(3,6);
complex<T> spa35 = SPA(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s234 = SS(2,3,4);
complex<T> t11 = square(spa34); 
complex<T> t12 = square(spb12); 
complex<T> t14 = spa56*(spa12*spb24 + spa13*spb34); 
complex<T> t15 = square(spa25*spb24 + spa35*spb34); 
complex<T> t16 = square(spa12*spb26 + spa13*spb36); 
complex<T> t23 = (spa12*spb24 + spa13*spb34)*spb56; 
complex<T> t28 = spa35*spb23 + spa45*spb24; 
complex<T> t29 = -(spa13*spb16) - spa23*spb26; 
complex<T> t34 = square(spa15); 
complex<T> t38 = square(spb46); 
complex<T> d5 = s23*s234*s56; d5 = T(1)/d5;
complex<T> d6 = s123*s23*s56; d6 = T(1)/d6;
complex<T> d7 = T(6); d7 = T(1)/d7;
complex<T> d12 = T(2); d12 = T(1)/d12;
complex<T> d14 = s234*spa23*(spa24*spb12 + spa34*spb13)*spb56; d14 = T(1)/d14;
complex<T> d15 = s123*spa56*(spa24*spb12 + spa34*spb13)*spb23; d15 = T(1)/d15;
complex<T> t13 = -(d5*spa34*spb16*t28) - d6*spa45*spb12*t29; 
complex<T> d1 = s123*(-s123 + s23)*t23; d1 = T(1)/d1;
complex<T> d2 = s123*t23*square(s123 - s23)*T(2); d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spa23*t23; d3 = T(1)/d3;
complex<T> d4 = spa23*t23*square(s123 - s56)*T(2); d4 = T(1)/d4;
complex<T> d8 = (s23 - s234)*s234*t14; d8 = T(1)/d8;
complex<T> d9 = s234*t14*square(s23 - s234)*T(2); d9 = T(1)/d9;
complex<T> d10 = spb23*t14*square(s234 - s56)*T(2); d10 = T(1)/d10;
complex<T> d11 = (-s234 + s56)*spb23*t14; d11 = T(1)/d11;
complex<T> d13 = s234*spb23*t14; d13 = T(1)/d13;
complex<T> d16 = s123*spa23*t23; d16 = T(1)/d16;
complex<T> t7 = d14*t11*square(spb16) - d13*square(t28); 
complex<T> t8 = -(d15*t12*square(spa45)) + d16*square(spa13*spb16 + spa23*spb26); 
complex<T> t26 = -(d8*(spa25*spb24 + spa35*spb34)); 
complex<T> t27 = d1*(spa12*spb26 + spa13*spb36); 
complex<T> t51 = d11*spa15; 
complex<T> t53 = d3*spa34; 
complex<T> t2 = s23*(t7 + t8); 
complex<T> t4 = s12*t8; 
complex<T> t5 = d9*spb23*t11*t15 - d2*spa23*t12*t16 + spa34*t26*t28*T(2) - spb12*t27*t29*T(2) + d7*t13*T(13); 
complex<T> t6 = -(d10*s234*t12*t34) + d4*s123*t11*t38 - spb12*t28*t51*T(2) - spb46*t29*t53*T(2) - d12*t13*T(3); 
complex<T> t9 = -(d9*spb23*t11*t15) + d10*s234*t12*t34 - spa34*t26*t28*T(2) + spb12*t28*t51*T(2); 
complex<T> t10 = d2*spa23*t12*t16 - d4*s123*t11*t38 + spb12*t27*t29*T(2) + spb46*t29*t53*T(2); 
complex<T> t46 = d12*t7; 
complex<T> t59 = s34*t8; 
complex<T> t3 = -t46; 
complex<T> t72 = d12*t4; 
complex<T> t78 = s12*t46; 
complex<T> t82 = s34*t46; 
complex<T> t86 = d12*t59; 
complex<T> t1 = s12*t13 + t72 + t78; 
complex<T> t44 = s34*t13 + t82 + t86; 
complex<T> t45 = s56*(t3 - d12*t8); 
complex<T> co1 = s234*t78; 
complex<T> co2 = s23*t82; 
complex<T> co3 = s23*t72; 
complex<T> co4 = s123*t86; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t10*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t9*(*CI_users[2]->get_value(mc,ind,mu)) + t6*(*CI_users[3]->get_value(mc,ind,mu)) + t1*(*CI_users[4]->get_value(mc,ind,mu)) + t2*(*CI_users[5]->get_value(mc,ind,mu)) + t44*(*CI_users[6]->get_value(mc,ind,mu)) + t45*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + co3*(*CI_users[10]->get_value(mc,ind,mu)) + co4*(*CI_users[11]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q2l_qmQpQmqpemep_L_wCI::C2q2Q2l_qmQpQmqpemep_L_wCI
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
     C2q2Q2l_qmQpQmqpemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmQpQmqpemep L");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:24:15 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb46 = SPB(4,6);
complex<T> spb56 = SPB(5,6);
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb16 = SPB(1,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb36 = SPB(3,6);
complex<T> s123 = SS(1,2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s124 = SS(1,2,4);
complex<T> s234 = SS(2,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s134 = SS(1,3,4);
complex<T> t20 = -(s12*T(2)); 
complex<T> t32 = square(spa13); 
complex<T> t33 = square(spb24); 
complex<T> t35 = square(spa15*spb13 + spa25*spb23); 
complex<T> t43 = square(spa24*spb26 + spa34*spb36); 
complex<T> t64 = square(spa45); 
complex<T> t65 = square(spb16); 
complex<T> t69 = -(spa14*spb12); 
complex<T> t70 = spa34*spb14; 
complex<T> t81 = spa24*spb12 + spa34*spb13; 
complex<T> t82 = -(spa23*spb26) + spa34*spb46; 
complex<T> t83 = -(spa15*spb12) + spa35*spb23; 
complex<T> t85 = spa25*spb12 + spa35*spb13; 
complex<T> t86 = spa23*spb36 + spa24*spb46; 
complex<T> t89 = square(spb46); 
complex<T> t90 = -s123 + s56; 
complex<T> t91 = -s234 + s56; 
complex<T> t92 = -(spa15*spb13) - spa25*spb23; 
complex<T> t93 = spa12*spb24 + spa13*spb34; 
complex<T> t94 = -(spa24*spb26) - spa34*spb36; 
complex<T> t122 = spa14*spb13; 
complex<T> t123 = spa24*spb23; 
complex<T> t153 = s34*s56; 
complex<T> d2 = spa23*spb13 + spa24*spb14; d2 = T(1)/d2;
complex<T> d3 = spa56; d3 = T(1)/d3;
complex<T> d4 = spb56; d4 = T(1)/d4;
complex<T> d12 = s123*s23*s56; d12 = T(1)/d12;
complex<T> d13 = s23*s234*s56; d13 = T(1)/d13;
complex<T> d14 = T(6); d14 = T(1)/d14;
complex<T> d20 = s234*(-s234 + s34)*spa24*(spa23*spb13 + spa24*spb14)*spb56; d20 = T(1)/d20;
complex<T> d24 = T(2); d24 = T(1)/d24;
complex<T> d35 = square(s123)*T(2); d35 = T(1)/d35;
complex<T> d36 = s234; d36 = T(1)/d36;
complex<T> d41 = s123; d41 = T(1)/d41;
complex<T> d44 = square(s234)*T(2); d44 = T(1)/d44;
complex<T> t34 = spa56*t81; 
complex<T> t36 = square(t86); 
complex<T> t38 = square(t85); 
complex<T> t39 = square(spa15) + d2*(spa13*spb23 + spa14*spb24)*square(spa25); 
complex<T> t41 = d2*spa25*spb16*(spa13*spb23 + spa14*spb24) - spa15*spb26; 
complex<T> t42 = d2*(spa13*spb23 + spa14*spb24)*t65 + square(spb26); 
complex<T> t66 = spb56*t81; 
complex<T> t87 = s34*t20 + s56*t20 + square(s12) + square(s34) + square(s56) - t153*T(2); 
complex<T> t124 = d24*s34; 
complex<T> d1 = (s12 - s123)*s123*spa56*spb13*(t122 + t123); d1 = T(1)/d1;
complex<T> d5 = (t122 + t123)*(s34*t20 + s56*t20 + square(s12) + square(s34) + square(s56) - t153*T(2)); d5 = T(1)/d5;
complex<T> d22 = t122 + t123; d22 = T(1)/d22;
complex<T> d23 = (spa23*spb13 + spa24*spb14)*(s34*t20 + s56*t20 + square(s12) + square(s34) + square(s56) - t153*T(2)); d23 = T(1)/d23;
complex<T> d26 = s123*spa23*spb56*t93; d26 = T(1)/d26;
complex<T> d28 = s234*spa56*spb23*t93; d28 = T(1)/d28;
complex<T> d30 = (t122 + t123)*t81*T(2)*(s34*t20 + s56*t20 + square(s12) + square(s34) + square(s56) - t153*T(2)); d30 = T(1)/d30;
complex<T> d31 = s123*t81; d31 = T(1)/d31;
complex<T> d32 = (t122 + t123)*t81*T(2); d32 = T(1)/d32;
complex<T> d34 = spa23*spb56*t93; d34 = T(1)/d34;
complex<T> d38 = s234*t81; d38 = T(1)/d38;
complex<T> d39 = (spa23*spb13 + spa24*spb14)*t81*T(2); d39 = T(1)/d39;
complex<T> d40 = (spa23*spb13 + spa24*spb14)*t81*T(2)*(s34*t20 + s56*t20 + square(s12) + square(s34) + square(s56) - t153*T(2)); d40 = T(1)/d40;
complex<T> d42 = spa56*spb23*t93; d42 = T(1)/d42;
complex<T> t1 = d1*spa13*spb12*t35 - d5*(d3*s123*spb12*t39 + d4*s124*spa12*t42 - s12*t41*T(2)); 
complex<T> t27 = s34*(d22*spa45*(spa13*spb14 + spa23*spb24)*spb36 - spa35*spb46); 
complex<T> t40 = d22*(spa13*spb14 + spa23*spb24)*t64 + square(spa35); 
complex<T> t44 = t89 + d22*(spa13*spb14 + spa23*spb24)*square(spb36); 
complex<T> t118 = d26*t32; 
complex<T> t157 = s123*t34; 
complex<T> t169 = s234*t66; 
complex<T> d9 = (t122 + t123)*t34*t90; d9 = T(1)/d9;
complex<T> d10 = spb23*t34*t90; d10 = T(1)/d10;
complex<T> d11 = spb23*t34*square(t90)*T(2); d11 = T(1)/d11;
complex<T> d18 = spa23*t66*square(t91)*T(2); d18 = T(1)/d18;
complex<T> d19 = (spa23*spb13 + spa24*spb14)*t66*t91; d19 = T(1)/d19;
complex<T> d21 = spa23*t66*t91; d21 = T(1)/d21;
complex<T> d33 = spb23*t34; d33 = T(1)/d33;
complex<T> d43 = spa23*t66; d43 = T(1)/d43;
complex<T> t3 = -(d20*spa34*spb24*t36) + d23*(d3*s134*spb34*t40 + d4*s234*spa34*t44 - t27*T(2)); 
complex<T> t77 = d10*t83; 
complex<T> t164 = d21*t82; 
complex<T> d6 = t157*square(s123 - s23)*T(2); d6 = T(1)/d6;
complex<T> d7 = (-s123 + s23)*spb13*t157; d7 = T(1)/d7;
complex<T> d8 = (-s123 + s23)*t157; d8 = T(1)/d8;
complex<T> d15 = t169*square(s23 - s234)*T(2); d15 = T(1)/d15;
complex<T> d16 = (s23 - s234)*spa24*t169; d16 = T(1)/d16;
complex<T> d17 = (s23 - s234)*t169; d17 = T(1)/d17;
complex<T> d25 = spb23*t157*square(spb13); d25 = T(1)/d25;
complex<T> d27 = spb23*t157*square(t122 + t123); d27 = T(1)/d27;
complex<T> d29 = spa23*t169*square(spa23*spb13 + spa24*spb14); d29 = T(1)/d29;
complex<T> d37 = spa23*t169*square(spa24); d37 = T(1)/d37;
complex<T> t2 = d3*d5*s123*spb12*t39 + d5*t20*t41 + d4*d5*s124*spa12*t42 + d9*s123*spb24*t64 - d11*s123*t33*t64 + d18*s234*t32*t65 + spa45*spb24*t77*T(2) + d23*(-(d3*s134*spb34*t40) - d4*s234*spa34*t44 + t27*T(2)) - d13*d24*spa15*spb24*t82*T(3) + spa13*(-(d19*s234*t65) + spb16*t164*T(2) - d12*d24*spb46*t83*T(3)); 
complex<T> t14 = s12*(t118*t89 - d25*t35*square(spb12)); 
complex<T> t15 = -(d28*s23*t33*square(spa15)) + d37*s23*t36*square(spa34); 
complex<T> t17 = d6*spb23*t32*t38 - d15*spa23*t33*t43 - d8*spa13*t83*t85 - d7*spa13*spb12*t85*t92 + d17*spb24*t82*t94 + d16*spa34*spb24*t86*t94 + d13*d14*spa15*spb24*t82*T(13) + d12*d14*spa13*spb46*t83*T(13); 
complex<T> t18 = t118*t89 - d27*t35*square(spa34*spb23 + t69); 
complex<T> t19 = -(d28*t33*square(spa15)) + d29*t36*square(-(spa23*spb12) + t70); 
complex<T> t23 = -(d1*spa13*spb12*t35) - d6*spb23*t32*t38 - d9*s123*spb24*t64 + d11*s123*t33*t64 + d8*spa13*t83*t85 + d7*spa13*spb12*t85*t92 - spa45*spb24*t77*T(2); 
complex<T> t30 = d20*spa34*spb24*t36 + d15*spa23*t33*t43 + d19*s234*spa13*t65 - d18*s234*t32*t65 - d17*spb24*t82*t94 - d16*spa34*spb24*t86*t94 - spa13*spb16*t164*T(2); 
complex<T> t4 = (s123 - s23)*(t18 - t118*t89 + d25*t35*square(spb12)); 
complex<T> t5 = t15 + s23*t18; 
complex<T> t6 = (s156 - s34)*(t19 + d28*t33*square(spa15) - d37*t36*square(spa34)); 
complex<T> t10 = (d41*s12*s56 - d24*(s12 - s34 + s56))*t18 - d40*(s12 + s34 - s56)*(spa35*spb36 + spa45*spb46)*(s234*spa23*spb12 - s234*t70 + s56*t70) + d39*(spa25*spb26*(spa23*spb12 - t70) + spa15*spb16*(spa23*spb12 + t70)) + d38*spa15*spa34*spb16*spb24*T(2) + d44*(d42*t33*square(spa15) - d43*square(t82))*(s12*s234 - s234*(s34 + s56) + t153*T(2)); 
complex<T> t11 = (d24*(s12 - s34 - s56) + d36*t153)*t19 - d32*(spa34*spa35*spb23*spb36 + spa34*spa45*spb23*spb46 + spa35*spb36*t69 - spa45*spb46*t69) - d30*(s12 + s34 - s56)*(spa15*spb16 + spa25*spb26)*(s56*t69 - s123*(spa34*spb23 + t69)) - d31*spa13*spa45*spb12*spb46*T(2) + d35*(-(d34*t32*t89) + d33*square(t83))*(-(s123*(s12 - s34 + s56)) + s12*s56*T(2)); 
complex<T> t167 = d24*t19; 
complex<T> t192 = t124*t18; 
complex<T> t7 = -(s56*(t167 + d24*t18)); 
complex<T> t9 = t124*t19 + t192 + d13*s34*spa15*spb24*t82 + d12*s34*spa13*spb46*t83; 
complex<T> t184 = s12*t167; 
complex<T> t8 = t14 - d24*s12*t18 + t184 + d13*s12*spa15*spb24*t82 + d12*s12*spa13*spb46*t83; 
complex<T> co1 = s234*t184; 
complex<T> co2 = t124*t15; 
complex<T> co3 = d24*s23*t14; 
complex<T> co4 = s123*t192; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t23*(*CI_users[1]->get_value(mc,ind,mu)) + t17*(*CI_users[2]->get_value(mc,ind,mu)) + t30*(*CI_users[3]->get_value(mc,ind,mu)) + t3*(*CI_users[4]->get_value(mc,ind,mu)) + t2*(*CI_users[5]->get_value(mc,ind,mu)) + t8*(*CI_users[6]->get_value(mc,ind,mu)) + t4*(*CI_users[7]->get_value(mc,ind,mu)) + t11*(*CI_users[8]->get_value(mc,ind,mu)) + t5*(*CI_users[9]->get_value(mc,ind,mu)) + t6*(*CI_users[10]->get_value(mc,ind,mu)) + t9*(*CI_users[11]->get_value(mc,ind,mu)) + t10*(*CI_users[12]->get_value(mc,ind,mu)) + t7*(*CI_users[13]->get_value(mc,ind,mu)) + co1*(*CI_users[14]->get_value(mc,ind,mu)) + co2*(*CI_users[15]->get_value(mc,ind,mu)) + co3*(*CI_users[16]->get_value(mc,ind,mu)) + co4*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q2l_qmQmQpqpemep_L_wCI::C2q2Q2l_qmQmQpqpemep_L_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c456));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qmQmQpqpemep_L_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmQmQpqpemep L");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:24:22 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa56 = SPA(5,6);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb46 = SPB(4,6);
complex<T> spb56 = SPB(5,6);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb36 = SPB(3,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb16 = SPB(1,6);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s123 = SS(1,2,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s234 = SS(2,3,4);
complex<T> t11 = square(spa12); 
complex<T> t12 = square(spb34); 
complex<T> t14 = spa56*(spa24*spb12 + spa34*spb13); 
complex<T> t15 = square(spa25*spb12 + spa35*spb13); 
complex<T> t16 = square(spa24*spb26 + spa34*spb36); 
complex<T> t23 = (spa24*spb12 + spa34*spb13)*spb56; 
complex<T> t28 = spa23*spb36 + spa24*spb46; 
complex<T> t29 = -(spa15*spb13) - spa25*spb23; 
complex<T> t34 = square(spa45); 
complex<T> t37 = square(spb16); 
complex<T> d5 = s123*s23*s56; d5 = T(1)/d5;
complex<T> d6 = s23*s234*s56; d6 = T(1)/d6;
complex<T> d7 = T(6); d7 = T(1)/d7;
complex<T> d12 = T(2); d12 = T(1)/d12;
complex<T> d14 = s123*spa23*(spa12*spb24 + spa13*spb34)*spb56; d14 = T(1)/d14;
complex<T> d15 = s234*spa56*spb23*(spa12*spb24 + spa13*spb34); d15 = T(1)/d15;
complex<T> t13 = -(d6*spa15*spb34*t28) - d5*spa12*spb46*t29; 
complex<T> d1 = s123*t14*square(s123 - s23)*T(2); d1 = T(1)/d1;
complex<T> d2 = s123*(-s123 + s23)*t14; d2 = T(1)/d2;
complex<T> d3 = (-s123 + s56)*spb23*t14; d3 = T(1)/d3;
complex<T> d4 = spb23*t14*square(s123 - s56)*T(2); d4 = T(1)/d4;
complex<T> d8 = s234*t23*square(s23 - s234)*T(2); d8 = T(1)/d8;
complex<T> d9 = (s23 - s234)*s234*t23; d9 = T(1)/d9;
complex<T> d10 = spa23*t23*square(s234 - s56)*T(2); d10 = T(1)/d10;
complex<T> d11 = (-s234 + s56)*spa23*t23; d11 = T(1)/d11;
complex<T> d13 = s123*spb23*t14; d13 = T(1)/d13;
complex<T> d16 = s234*spa23*t23; d16 = T(1)/d16;
complex<T> t7 = -(d13*square(spa15*spb13 + spa25*spb23)) + d14*t11*square(spb46); 
complex<T> t8 = -(d15*t12*square(spa15)) + d16*square(t28); 
complex<T> t26 = d2*(spa25*spb12 + spa35*spb13); 
complex<T> t27 = -(d9*(spa24*spb26 + spa34*spb36)); 
complex<T> t51 = d11*spa12; 
complex<T> t53 = d3*spa45; 
complex<T> t2 = s23*(t7 + t8); 
complex<T> t4 = s12*t8; 
complex<T> t5 = -(d4*s123*t12*t34) + d10*s234*t11*t37 + spb16*t28*t51*T(2) + spb34*t29*t53*T(2) - d12*t13*T(3); 
complex<T> t6 = d1*spb23*t11*t15 - d8*spa23*t12*t16 - spb34*t27*t28*T(2) + spa12*t26*t29*T(2) + d7*t13*T(13); 
complex<T> t9 = -(d1*spb23*t11*t15) + d4*s123*t12*t34 - spa12*t26*t29*T(2) - spb34*t29*t53*T(2); 
complex<T> t10 = d8*spa23*t12*t16 - d10*s234*t11*t37 + spb34*t27*t28*T(2) - spb16*t28*t51*T(2); 
complex<T> t46 = d12*t7; 
complex<T> t59 = s34*t8; 
complex<T> t3 = -t46; 
complex<T> t72 = d12*t4; 
complex<T> t78 = s12*t46; 
complex<T> t82 = s34*t46; 
complex<T> t86 = d12*t59; 
complex<T> t1 = s12*t13 + t72 + t78; 
complex<T> t44 = s34*t13 + t82 + t86; 
complex<T> t45 = s56*(t3 - d12*t8); 
complex<T> co1 = s234*t72; 
complex<T> co2 = s23*t86; 
complex<T> co3 = s23*t78; 
complex<T> co4 = s123*t82; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t9*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t10*(*CI_users[2]->get_value(mc,ind,mu)) + t5*(*CI_users[3]->get_value(mc,ind,mu)) + t1*(*CI_users[4]->get_value(mc,ind,mu)) + t2*(*CI_users[5]->get_value(mc,ind,mu)) + t44*(*CI_users[6]->get_value(mc,ind,mu)) + t45*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + co3*(*CI_users[10]->get_value(mc,ind,mu)) + co4*(*CI_users[11]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q2l_qpqmQpQmemep_sl_wCI::C2q2Q2l_qpqmQpQmemep_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c134, c256));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c256, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c256));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c56, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c56, c12, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c56, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qpqmQpQmemep_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, Qp, Qm, em, ep}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpqmQpQmemep sl");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:27:06 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb16 = SPB(1,6);
complex<T> spa25 = SPA(2,5);
complex<T> spb26 = SPB(2,6);
complex<T> spa56 = SPA(5,6);
complex<T> spb34 = SPB(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb36 = SPB(3,6);
complex<T> spb12 = SPB(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spb46 = SPB(4,6);
complex<T> spa23 = SPA(2,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb25 = SPB(2,5);
complex<T> spa16 = SPA(1,6);
complex<T> spb35 = SPB(3,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa26 = SPA(2,6);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = S(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s134 = SS(1,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s15 = -(spa15*spb15);
complex<T> s16 = -(spa16*spb16);
complex<T> s25 = -(spa25*spb25);
complex<T> s26 = -(spa26*spb26);
complex<T> s234 = SS(2,3,4);
complex<T> s256 = SS(2,5,6);
complex<T> t13 = -(spa56*spb26); 
complex<T> t61 = spa14*spb16; 
complex<T> t62 = spa15*spb13; 
complex<T> t63 = spa25*spb23; 
complex<T> t64 = spa24*spb26; 
complex<T> t66 = (spa14*spb13 + spa24*spb23)*(spa15*spb16 + spa25*spb26); 
complex<T> t84 = s12 - s34 - s56; 
complex<T> t88 = -s12 - s34 + s56; 
complex<T> t89 = -s12 + s34 - s56; 
complex<T> t91 = square(spa14); 
complex<T> t92 = square(spa15); 
complex<T> t93 = square(spb23); 
complex<T> t94 = square(spb26); 
complex<T> t99 = -s156 + s34; 
complex<T> t100 = -s256 + s34; 
complex<T> t101 = square(spa24); 
complex<T> t102 = square(spa25); 
complex<T> t103 = square(spb13); 
complex<T> t104 = square(spb16); 
complex<T> t124 = spa13*spb23; 
complex<T> t125 = spa14*spb24; 
complex<T> t126 = spa15*spb25; 
complex<T> t127 = spa16*spb26; 
complex<T> t128 = spa45*spb34; 
complex<T> t129 = spa34*spb36; 
complex<T> t130 = s34*s56; 
complex<T> t156 = spa34*spb56; 
complex<T> t169 = spa14*spb16 - spa24*spb26; 
complex<T> t181 = spa56*spb34; 
complex<T> t199 = spa56*spb36; 
complex<T> t203 = spa14*spb26; 
complex<T> t204 = spa15*spb23; 
complex<T> t223 = -(s12*spa45); 
complex<T> t228 = spa14*spa24; 
complex<T> t229 = spa15*spa25; 
complex<T> t236 = spb13*spb23; 
complex<T> t237 = spb16*spb26; 
complex<T> d17 = T(2); d17 = T(1)/d17;
complex<T> d48 = s134; d48 = T(1)/d48;
complex<T> d49 = s156; d49 = T(1)/d49;
complex<T> t12 = -t66; 
complex<T> t30 = (t124 + t125)*T(2); 
complex<T> t31 = spb34*(t126 + t127); 
complex<T> t32 = square(t128 - t63); 
complex<T> t33 = square(t129 + t61); 
complex<T> t34 = square(t199 + t62); 
complex<T> t38 = -square(s12) - square(s34) - square(s56) + s12*s34*T(2) + s12*s56*T(2) + t130*T(2) + s12*t84*T(3); 
complex<T> t44 = (s13 + s14 - s23 - s24)*spa45 - spa56*t129*T(2); 
complex<T> t45 = -((s15 + s16 - s25 - s26)*spa45) + spa56*t129*T(2); 
complex<T> t51 = d48*t130 + d17*t84; 
complex<T> t52 = d49*t130 + d17*t84; 
complex<T> t67 = -t129; 
complex<T> t68 = (t126 + t127)*T(2); 
complex<T> t69 = spa56*(t124 + t125); 
complex<T> t72 = -(t89*T(3)); 
complex<T> t90 = t128 - t63; 
complex<T> t95 = spa45*spb56 - t64; 
complex<T> t105 = -t199 - t62; 
complex<T> t110 = t62 - t63; 
complex<T> t111 = -t61 + t64; 
complex<T> t166 = -t62 + t63; 
complex<T> t178 = -t204; 
complex<T> t182 = t66*T(3); 
complex<T> t211 = spb56*t128; 
complex<T> d1 = t124 + t125; d1 = T(1)/d1;
complex<T> d5 = t126 + t127; d5 = T(1)/d5;
complex<T> d6 = spa34*(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d6 = T(1)/d6;
complex<T> d7 = (t124 + t125)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d7 = T(1)/d7;
complex<T> d8 = t181*square(t124 + t125)*T(2); d8 = T(1)/d8;
complex<T> d10 = (t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d10 = T(1)/d10;
complex<T> d12 = t181*square(t126 + t127)*T(2); d12 = T(1)/d12;
complex<T> d14 = square(t124 + t125)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d14 = T(1)/d14;
complex<T> d15 = s234*t130; d15 = T(1)/d15;
complex<T> d16 = s134*t130; d16 = T(1)/d16;
complex<T> d18 = spb56*(t124 + t125)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d18 = T(1)/d18;
complex<T> d19 = t156*square(t124 + t125)*T(2); d19 = T(1)/d19;
complex<T> d20 = s134*(t124 + t125)*t156*T(4); d20 = T(1)/d20;
complex<T> d21 = square(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d21 = T(1)/d21;
complex<T> d22 = t156*square(t126 + t127)*T(2); d22 = T(1)/d22;
complex<T> d23 = s256*(t126 + t127)*t156*T(4); d23 = T(1)/d23;
complex<T> d24 = (-s134 + s56)*t156*square(t124 + t125); d24 = T(1)/d24;
complex<T> d26 = (-s134 + s56)*(t124 + t125)*t156; d26 = T(1)/d26;
complex<T> d27 = t181*t99*square(t126 + t127); d27 = T(1)/d27;
complex<T> d30 = (-s234 + s56)*t181*square(t124 + t125); d30 = T(1)/d30;
complex<T> d33 = t100*t156*square(t126 + t127); d33 = T(1)/d33;
complex<T> d35 = t100*(t126 + t127)*t156; d35 = T(1)/d35;
complex<T> d36 = t181*square(t126 + t127); d36 = T(1)/d36;
complex<T> d38 = t156*square(t126 + t127); d38 = T(1)/d38;
complex<T> d40 = t181*square(t124 + t125); d40 = T(1)/d40;
complex<T> d42 = t156*square(t124 + t125); d42 = T(1)/d42;
complex<T> d44 = s156*t181*cube(t126 + t127); d44 = T(1)/d44;
complex<T> d45 = s156*(spa25*spb15 + spa26*spb16)*t156; d45 = T(1)/d45;
complex<T> d46 = s134*(spa23*spb13 + spa24*spb14)*t181; d46 = T(1)/d46;
complex<T> d47 = s134*t156*cube(t124 + t125); d47 = T(1)/d47;
complex<T> d50 = square(t124 + t125); d50 = T(1)/d50;
complex<T> d52 = (t124 + t125)*square(s134); d52 = T(1)/d52;
complex<T> d54 = (t126 + t127)*square(s156); d54 = T(1)/d54;
complex<T> d55 = square(t126 + t127); d55 = T(1)/d55;
complex<T> t10 = -(d45*t101*t104) + d44*square(spa15*spb12 + spa56*spb26)*square(spa15*spb35 + spa16*spb36); 
complex<T> t11 = -(d46*t102*t103) + d47*square(spa14*spb12 - spa34*spb23)*square(spa13*spb36 + spa14*spb46); 
complex<T> t36 = square(t95); 
complex<T> t37 = (-s13 - s14 + s23 + s24)*spb36 - t211*T(2); 
complex<T> t40 = s134*t62 - s234*t63 + d1*s234*t204*t84; 
complex<T> t41 = -(s156*t62) + s256*t63 + d5*s156*t204*t84; 
complex<T> t42 = -(s134*t61) + s234*t64 + d1*s134*t203*t84; 
complex<T> t43 = s156*t61 - s256*t64 + d5*s256*t203*t84; 
complex<T> t46 = (s15 + s16 - s25 - s26)*spb36 + t211*T(2); 
complex<T> t80 = -(d19*(spa14*spb16 + t129)); 
complex<T> t81 = -(d12*(spa15*spb13 + t199)); 
complex<T> t82 = d22*(spa24*spb26 - spa45*spb56); 
complex<T> t96 = -t61 + t67; 
complex<T> t143 = d20*t33; 
complex<T> t198 = d35*t95; 
complex<T> t215 = d8*t90; 
complex<T> d2 = t69*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d2 = T(1)/d2;
complex<T> d3 = t30*square(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d3 = T(1)/d3;
complex<T> d4 = t68*square(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d4 = T(1)/d4;
complex<T> d9 = s234*spb34*t69*T(4); d9 = T(1)/d9;
complex<T> d11 = t31*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d11 = T(1)/d11;
complex<T> d13 = s156*spa56*t31*T(4); d13 = T(1)/d13;
complex<T> d25 = t156*t30*square(s134 - s56); d25 = T(1)/d25;
complex<T> d28 = spa56*t31*square(t99)*T(2); d28 = T(1)/d28;
complex<T> d29 = spa56*t31*t99; d29 = T(1)/d29;
complex<T> d31 = t181*t30*square(s234 - s56); d31 = T(1)/d31;
complex<T> d32 = (-s234 + s56)*spb34*t69; d32 = T(1)/d32;
complex<T> d34 = t156*t68*square(t100); d34 = T(1)/d34;
complex<T> d37 = s156*spa56*t31*T(2); d37 = T(1)/d37;
complex<T> d39 = s256*t156*t68; d39 = T(1)/d39;
complex<T> d41 = s234*t181*t30; d41 = T(1)/d41;
complex<T> d43 = s134*t156*t30; d43 = T(1)/d43;
complex<T> d51 = t30*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d51 = T(1)/d51;
complex<T> d53 = t68*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t130*T(2)); d53 = T(1)/d53;
complex<T> t5 = s12*t11; 
complex<T> t6 = d53*(s15 + s16 - s25 - s26)*spb36*t223 + d4*(s15 + s16 - s25 - s26)*t12*t38 + d55*t178*(spa45*spb56 - t61) + d21*t178*t84*(spa45*spb56*t88 + t129*t89) + d10*t156*square(t62 + t63) - d54*spa24*spb16*(spa15*spb35 + spa16*spb36)*(spa15*spb12 - t13)*T(2); 
complex<T> t26 = t203*t80 - d25*s134*t101*t94 + d24*s134*t228*t94 - d26*t64*t96*T(2) - t143*T(3); 
complex<T> t35 = d15*spa24*spb16*(-t128 + t63) + d16*spa25*spb13*t96; 
complex<T> t70 = d3*(s13 + s14 - s23 - s24); 
complex<T> t79 = d32*(-(spa25*spb23) + t128); 
complex<T> t132 = d17*t10; 
complex<T> t141 = d29*t105; 
complex<T> t144 = d23*t36; 
complex<T> t160 = d3*(-s13 - s14 + s23 + s24); 
complex<T> t192 = d9*t32; 
complex<T> t193 = d13*t34; 
complex<T> t2 = -(s56*(d17*t11 + t132)); 
complex<T> t4 = d51*(s13 + s14 - s23 - s24)*spb36*t223 - d50*t203*(t128 - t62) + t12*t38*t70 - d14*t203*t84*(t199*t88 + t128*t89) + d7*t181*square(t61 + t64) + d52*spa25*spb13*(spa14*spb12 - spa34*spb23)*(spa13*spb36 + spa14*spb46)*T(2); 
complex<T> t25 = -(t35*T(3)); 
complex<T> t27 = t204*t81 - d28*s156*t102*t93 + d27*s156*t229*t93 + d29*(t199 + t62)*t63*T(2) - t193*T(3); 
complex<T> t29 = t203*t82 + d34*s256*t104*t91 - d33*s256*t237*t91 + t198*t61*T(2) + t144*T(3); 
complex<T> t177 = s34*(d17*t11 + t132 + t35); 
complex<T> t233 = s12*t132; 
complex<T> t248 = t141*t63; 
complex<T> t251 = d17*t5; 
complex<T> t266 = t62*t79; 
complex<T> t1 = t233 + t251 + s12*t35; 
complex<T> t7 = d7*(t111*t128 + t110*t129) + d17*t25 + d14*spa34*t178*t37 + d2*spa45*t40 + d18*spb36*t42 - d14*spb34*t203*t44 - d36*t178*(t199 + t62) + t182*t70*t88 - d34*s256*t104*t91 + d33*s256*t237*t91 + d28*s156*t102*t93 - d27*s156*t229*t93 + d38*t203*t95 - t198*t61*T(2) - d29*t199*t63*T(2) - d29*t62*t63*T(2) + d37*t34*T(3) - d39*t36*T(3) - t160*t66*t88*T(3); 
complex<T> t8 = -(d7*t111*t128) - d10*spa45*spb56*t166 - d10*t169*t199 + d17*t25 + d14*spa34*t204*t37 - d2*spa45*t40 - d11*spb36*t41 - d18*spb36*t42 - d6*spa45*t43 + d14*spb34*t203*t44 + d21*spb56*t178*t45 + d21*spa14*t13*t46 + d8*t178*(t128 - t63) + d7*t110*t67 + d4*(s15 + s16 - s25 - s26)*t66*t72 + t203*t80 + t204*t81 + t203*t82 + t160*t182*t88 + d4*(-s15 - s16 + s25 + s26)*t182*t89 - t143*T(3) + t144*T(3) + t192*T(3) - t193*T(3) - t66*t70*t88*T(3); 
complex<T> t9 = d10*spa45*spb56*t166 + d10*t169*t199 + d11*spb36*t41 + d6*spa45*t43 + d21*spb56*t204*t45 + d21*spa56*t203*t46 + d40*t204*(t128 - t63) + d4*(-s15 - s16 + s25 + s26)*t66*t72 + d4*(s15 + s16 - s25 - s26)*t182*t89 - d31*s234*t103*t92 + d30*s234*t236*t92 + d25*s134*t101*t94 - d24*s134*t228*t94 - d42*t203*t96 - t266*T(2) + d26*t64*t96*T(2) - d41*t32*T(3) + d43*t33*T(3); 
complex<T> t22 = d8*t178*(t128 - t63) + d31*s234*t103*t92 - d30*s234*t236*t92 + t266*T(2) + t192*T(3); 
complex<T> co1 = Complex(0,1)*t8; 
complex<T> co2 = Complex(0,1)*t26; 
complex<T> co3 = Complex(0,1)*t27; 
complex<T> co4 = Complex(0,1)*t22; 
complex<T> co5 = Complex(0,1)*t29; 
complex<T> co6 = Complex(0,1)*t7; 
complex<T> co7 = Complex(0,1)*t9; 
complex<T> co8 = Complex(0,1)*t1; 
complex<T> co9 = Complex(0,-1)*t100*t11; 
complex<T> co10 = Complex(0,1)*t11*t51; 
complex<T> co11 = Complex(0,1)*t10*t52; 
complex<T> co12 = Complex(0,-1)*t10*t99; 
complex<T> co13 = Complex(0,1)*t177; 
complex<T> co14 = Complex(0,1)*t4; 
complex<T> co15 = Complex(0,1)*t2; 
complex<T> co16 = Complex(0,1)*t6; 
complex<T> co17 = Complex(0,1)*s134*t251; 
complex<T> co18 = Complex(0,1)*s156*t233; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + co3*(*CI_users[2]->get_value(mc,ind,mu)) + co4*(*CI_users[3]->get_value(mc,ind,mu)) + co5*(*CI_users[4]->get_value(mc,ind,mu)) + co6*(*CI_users[5]->get_value(mc,ind,mu)) + co7*(*CI_users[6]->get_value(mc,ind,mu)) + co8*(*CI_users[7]->get_value(mc,ind,mu)) + co9*(*CI_users[8]->get_value(mc,ind,mu)) + co10*(*CI_users[9]->get_value(mc,ind,mu)) + co11*(*CI_users[10]->get_value(mc,ind,mu)) + co12*(*CI_users[11]->get_value(mc,ind,mu)) + co13*(*CI_users[12]->get_value(mc,ind,mu)) + co14*(*CI_users[13]->get_value(mc,ind,mu)) + co15*(*CI_users[14]->get_value(mc,ind,mu)) + co16*(*CI_users[15]->get_value(mc,ind,mu)) + co17*(*CI_users[16]->get_value(mc,ind,mu)) + co18*(*CI_users[17]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q2l_qpqmQmQpemep_sl_wCI::C2q2Q2l_qpqmQmQpemep_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c134, c256));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c256, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c256));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c56, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c56, c12, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c56, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qpqmQmQpemep_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, Qm, Qp, em, ep}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpqmQmQpemep sl");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:29:51 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spa14 = SPA(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb16 = SPB(1,6);
complex<T> spa25 = SPA(2,5);
complex<T> spb26 = SPB(2,6);
complex<T> spa56 = SPA(5,6);
complex<T> spb34 = SPB(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb46 = SPB(4,6);
complex<T> spb12 = SPB(1,2);
complex<T> spa34 = SPA(3,4);
complex<T> spb36 = SPB(3,6);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb56 = SPB(5,6);
complex<T> spb25 = SPB(2,5);
complex<T> spa16 = SPA(1,6);
complex<T> spb45 = SPB(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa26 = SPA(2,6);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = S(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s134 = SS(1,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s15 = -(spa15*spb15);
complex<T> s16 = -(spa16*spb16);
complex<T> s25 = -(spa25*spb25);
complex<T> s26 = -(spa26*spb26);
complex<T> s234 = SS(2,3,4);
complex<T> s256 = SS(2,5,6);
complex<T> t13 = -(spa56*spb26); 
complex<T> t14 = -(spa34*spb24); 
complex<T> t61 = spa13*spb16; 
complex<T> t62 = spa15*spb14; 
complex<T> t63 = spa25*spb24; 
complex<T> t64 = spa23*spb26; 
complex<T> t66 = (spa13*spb14 + spa23*spb24)*(spa15*spb16 + spa25*spb26); 
complex<T> t84 = s12 - s34 - s56; 
complex<T> t88 = -s12 - s34 + s56; 
complex<T> t89 = -s12 + s34 - s56; 
complex<T> t91 = square(spa13); 
complex<T> t92 = square(spa15); 
complex<T> t93 = square(spb24); 
complex<T> t94 = square(spb26); 
complex<T> t99 = -s156 + s34; 
complex<T> t100 = -s256 + s34; 
complex<T> t101 = square(spa23); 
complex<T> t102 = square(spa25); 
complex<T> t103 = square(spb14); 
complex<T> t104 = square(spb16); 
complex<T> t124 = spa13*spb23; 
complex<T> t125 = spa14*spb24; 
complex<T> t126 = spa15*spb25; 
complex<T> t127 = spa16*spb26; 
complex<T> t129 = s34*s56; 
complex<T> t130 = spa35*spb34; 
complex<T> t158 = spa34*spb56; 
complex<T> t159 = spa56*spb46; 
complex<T> t172 = spa13*spb16 - spa23*spb26; 
complex<T> t201 = spa56*spb34; 
complex<T> t208 = spa13*spb26; 
complex<T> t209 = spa15*spb24; 
complex<T> t212 = -(s12*spa35); 
complex<T> t233 = spa13*spa23; 
complex<T> t234 = spa15*spa25; 
complex<T> t236 = spb14*spb24; 
complex<T> t243 = spb16*spb26; 
complex<T> d17 = T(2); d17 = T(1)/d17;
complex<T> d48 = s134; d48 = T(1)/d48;
complex<T> d49 = s156; d49 = T(1)/d49;
complex<T> t12 = -t66; 
complex<T> t30 = (t124 + t125)*T(2); 
complex<T> t31 = spb34*(t126 + t127); 
complex<T> t32 = square(t130 + t63); 
complex<T> t34 = square(t159 + t62); 
complex<T> t38 = -square(s12) - square(s34) - square(s56) + s12*s34*T(2) + s12*s56*T(2) + t129*T(2) + s12*t84*T(3); 
complex<T> t44 = -((s15 + s16 - s25 - s26)*spa35) - spa34*t159*T(2); 
complex<T> t45 = (s13 + s14 - s23 - s24)*spa35 + spa34*t159*T(2); 
complex<T> t51 = d48*t129 + d17*t84; 
complex<T> t52 = d49*t129 + d17*t84; 
complex<T> t65 = -t130; 
complex<T> t68 = (t126 + t127)*T(2); 
complex<T> t69 = spa56*(t124 + t125); 
complex<T> t72 = -(t89*T(3)); 
complex<T> t90 = spa34*spb46 - t61; 
complex<T> t95 = spa35*spb56 - t64; 
complex<T> t105 = -t159 - t62; 
complex<T> t109 = t62 - t63; 
complex<T> t110 = -t61 + t64; 
complex<T> t168 = -t62 + t63; 
complex<T> t181 = -t209; 
complex<T> t185 = t66*T(3); 
complex<T> t229 = spb56*t130; 
complex<T> d1 = t124 + t125; d1 = T(1)/d1;
complex<T> d5 = t126 + t127; d5 = T(1)/d5;
complex<T> d6 = spa34*(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d6 = T(1)/d6;
complex<T> d7 = (t124 + t125)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d7 = T(1)/d7;
complex<T> d8 = t201*square(t124 + t125)*T(2); d8 = T(1)/d8;
complex<T> d10 = (t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d10 = T(1)/d10;
complex<T> d12 = t201*square(t126 + t127)*T(2); d12 = T(1)/d12;
complex<T> d14 = square(t124 + t125)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d14 = T(1)/d14;
complex<T> d15 = s234*t129; d15 = T(1)/d15;
complex<T> d16 = s134*t129; d16 = T(1)/d16;
complex<T> d18 = spb56*(t124 + t125)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d18 = T(1)/d18;
complex<T> d19 = t158*square(t124 + t125)*T(2); d19 = T(1)/d19;
complex<T> d20 = s134*(t124 + t125)*t158*T(4); d20 = T(1)/d20;
complex<T> d21 = square(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d21 = T(1)/d21;
complex<T> d22 = t158*square(t126 + t127)*T(2); d22 = T(1)/d22;
complex<T> d23 = s256*(t126 + t127)*t158*T(4); d23 = T(1)/d23;
complex<T> d24 = (-s134 + s56)*t158*square(t124 + t125); d24 = T(1)/d24;
complex<T> d26 = (-s134 + s56)*(t124 + t125)*t158; d26 = T(1)/d26;
complex<T> d27 = t201*t99*square(t126 + t127); d27 = T(1)/d27;
complex<T> d30 = (-s234 + s56)*t201*square(t124 + t125); d30 = T(1)/d30;
complex<T> d33 = t100*t158*square(t126 + t127); d33 = T(1)/d33;
complex<T> d35 = t100*(t126 + t127)*t158; d35 = T(1)/d35;
complex<T> d36 = t201*square(t126 + t127); d36 = T(1)/d36;
complex<T> d38 = t158*square(t126 + t127); d38 = T(1)/d38;
complex<T> d40 = t201*square(t124 + t125); d40 = T(1)/d40;
complex<T> d42 = t158*square(t124 + t125); d42 = T(1)/d42;
complex<T> d44 = s156*t201*cube(t126 + t127); d44 = T(1)/d44;
complex<T> d45 = s156*(spa25*spb15 + spa26*spb16)*t158; d45 = T(1)/d45;
complex<T> d46 = s134*(spa23*spb13 + spa24*spb14)*t201; d46 = T(1)/d46;
complex<T> d47 = s134*t158*cube(t124 + t125); d47 = T(1)/d47;
complex<T> d50 = square(t124 + t125); d50 = T(1)/d50;
complex<T> d52 = (t124 + t125)*square(s134); d52 = T(1)/d52;
complex<T> d54 = (t126 + t127)*square(s156); d54 = T(1)/d54;
complex<T> d55 = square(t126 + t127); d55 = T(1)/d55;
complex<T> t10 = d45*t101*t104 - d44*square(spa15*spb12 + spa56*spb26)*square(spa15*spb45 + spa16*spb46); 
complex<T> t11 = d46*t102*t103 - d47*square(spa13*spb12 + spa34*spb24)*square(spa13*spb36 + spa14*spb46); 
complex<T> t33 = square(t90); 
complex<T> t36 = square(t95); 
complex<T> t37 = (-s13 - s14 + s23 + s24)*spb46 + t229*T(2); 
complex<T> t40 = s134*t62 - s234*t63 + d1*s234*t209*t84; 
complex<T> t41 = -(s156*t62) + s256*t63 + d5*s156*t209*t84; 
complex<T> t42 = -(s134*t61) + s234*t64 + d1*s134*t208*t84; 
complex<T> t43 = s156*t61 - s256*t64 + d5*s256*t208*t84; 
complex<T> t46 = (s15 + s16 - s25 - s26)*spb46 - t229*T(2); 
complex<T> t80 = d19*(-(spa13*spb16) + spa34*spb46); 
complex<T> t81 = -(d12*(spa15*spb14 + t159)); 
complex<T> t82 = d22*(spa23*spb26 - spa35*spb56); 
complex<T> t96 = -t63 + t65; 
complex<T> t199 = d26*t90; 
complex<T> t200 = d35*t95; 
complex<T> d2 = t69*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d2 = T(1)/d2;
complex<T> d3 = t30*square(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d3 = T(1)/d3;
complex<T> d4 = t68*square(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d4 = T(1)/d4;
complex<T> d9 = s234*spb34*t69*T(4); d9 = T(1)/d9;
complex<T> d11 = t31*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d11 = T(1)/d11;
complex<T> d13 = s156*spa56*t31*T(4); d13 = T(1)/d13;
complex<T> d25 = t158*t30*square(s134 - s56); d25 = T(1)/d25;
complex<T> d28 = spa56*t31*square(t99)*T(2); d28 = T(1)/d28;
complex<T> d29 = spa56*t31*t99; d29 = T(1)/d29;
complex<T> d31 = t201*t30*square(s234 - s56); d31 = T(1)/d31;
complex<T> d32 = (-s234 + s56)*spb34*t69; d32 = T(1)/d32;
complex<T> d34 = t158*t68*square(t100); d34 = T(1)/d34;
complex<T> d37 = s156*spa56*t31*T(2); d37 = T(1)/d37;
complex<T> d39 = s256*t158*t68; d39 = T(1)/d39;
complex<T> d41 = s234*t201*t30; d41 = T(1)/d41;
complex<T> d43 = s134*t158*t30; d43 = T(1)/d43;
complex<T> d51 = t30*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d51 = T(1)/d51;
complex<T> d53 = t68*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t129*T(2)); d53 = T(1)/d53;
complex<T> t6 = d53*(s15 + s16 - s25 - s26)*spb46*t212 + d4*(s15 + s16 - s25 - s26)*t12*t38 + d55*t181*(spa35*spb56 - t61) + d21*t181*t84*(spa35*spb56*t88 - spa34*spb46*t89) - d10*t158*square(t62 + t63) - d54*spa23*spb16*(spa15*spb45 + spa16*spb46)*(spa15*spb12 - t13)*T(2); 
complex<T> t35 = d16*spa25*spb14*t90 - d15*spa23*spb16*t96; 
complex<T> t70 = d3*(s13 + s14 - s23 - s24); 
complex<T> t79 = -(d32*(spa25*spb24 + t130)); 
complex<T> t143 = d29*t105; 
complex<T> t145 = d20*t33; 
complex<T> t146 = d23*t36; 
complex<T> t162 = d3*(-s13 - s14 + s23 + s24); 
complex<T> t188 = d17*t10; 
complex<T> t194 = d9*t32; 
complex<T> t195 = d13*t34; 
complex<T> t211 = d17*t11; 
complex<T> t232 = d8*t96; 
complex<T> t254 = t199*t64; 
complex<T> t1 = s56*(t188 + t211); 
complex<T> t3 = -t211; 
complex<T> t5 = d51*(s13 + s14 - s23 - s24)*spb46*t212 + d50*t208*(t62 - t65) + t12*t38*t70 - d14*t208*t84*(t159*t88 + t65*t89) - d7*t201*square(t61 + t64) + d52*spa25*spb14*(spa13*spb36 + spa14*spb46)*(spa13*spb12 - t14)*T(2); 
complex<T> t25 = t35*T(3); 
complex<T> t26 = t208*t80 - d25*s134*t101*t94 + d24*s134*t233*t94 - t254*T(2) - t145*T(3); 
complex<T> t27 = t209*t81 - d28*s156*t102*t93 + d27*s156*t234*t93 + d29*(t159 + t62)*t63*T(2) - t195*T(3); 
complex<T> t29 = t208*t82 + d34*s256*t104*t91 - d33*s256*t243*t91 + t200*t61*T(2) + t146*T(3); 
complex<T> t133 = -t188; 
complex<T> t253 = t143*t63; 
complex<T> t271 = t62*t79; 
complex<T> t2 = s12*(t133 + t3 - t35); 
complex<T> t7 = -(d10*spa35*spb56*t168) - d10*t159*t172 + d11*spb46*t41 + d6*spa35*t43 + d21*spb56*t181*t44 + d21*spa13*t13*t46 + d4*(s15 + s16 - s25 - s26)*t66*t72 + d4*(-s15 - s16 + s25 + s26)*t185*t89 - d42*t208*t90 - d31*s234*t103*t92 + d30*s234*t236*t92 + d25*s134*t101*t94 - d24*s134*t233*t94 + d40*t209*t96 + t254*T(2) - t271*T(2) - d41*t32*T(3) + d43*t33*T(3); 
complex<T> t8 = d7*(spa34*spb46*t109 + t110*t130) + d17*t25 + d14*spa15*t14*t37 - d2*spa35*t40 - d18*spb46*t42 - d14*spb34*t208*t45 - d36*t181*(t159 + t62) + t162*t185*t88 - d34*s256*t104*t91 + d33*s256*t243*t91 + d28*s156*t102*t93 - d27*s156*t234*t93 + d38*t208*t95 - t200*t61*T(2) - d29*t159*t63*T(2) - d29*t62*t63*T(2) + d37*t34*T(3) - d39*t36*T(3) - t66*t70*t88*T(3); 
complex<T> t9 = -(d7*spa34*spb46*t109) + d10*spa35*spb56*t168 + d10*t159*t172 + t181*t232 + d17*t25 + d14*spa34*t209*t37 + d2*spa35*t40 - d11*spb46*t41 + d18*spb46*t42 - d6*spa35*t43 + d21*spb56*t209*t44 + d14*spb34*t208*t45 + d21*spa56*t208*t46 + d7*t110*t65 + d4*(-s15 - s16 + s25 + s26)*t66*t72 + t208*t80 + t209*t81 + t208*t82 + t185*t70*t88 + d4*(s15 + s16 - s25 - s26)*t185*t89 - t145*T(3) + t146*T(3) + t194*T(3) - t195*T(3) - t162*t66*t88*T(3); 
complex<T> t22 = t181*t232 + d31*s234*t103*t92 - d30*s234*t236*t92 + t271*T(2) + t194*T(3); 
complex<T> t180 = s34*(t133 + t3 - t35); 
complex<T> co1 = Complex(0,1)*t9; 
complex<T> co2 = Complex(0,1)*t26; 
complex<T> co3 = Complex(0,1)*t27; 
complex<T> co4 = Complex(0,1)*t22; 
complex<T> co5 = Complex(0,1)*t29; 
complex<T> co6 = Complex(0,1)*t8; 
complex<T> co7 = Complex(0,1)*t7; 
complex<T> co8 = Complex(0,1)*t2; 
complex<T> co9 = Complex(0,1)*t100*t11; 
complex<T> co10 = Complex(0,-1)*t11*t51; 
complex<T> co11 = Complex(0,-1)*t10*t52; 
complex<T> co12 = Complex(0,1)*t10*t99; 
complex<T> co13 = Complex(0,1)*t180; 
complex<T> co14 = Complex(0,-1)*t5; 
complex<T> co15 = Complex(0,1)*t1; 
complex<T> co16 = Complex(0,-1)*t6; 
complex<T> co17 = Complex(0,-1)*s12*s134*t211; 
complex<T> co18 = Complex(0,-1)*s12*s156*t188; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + co3*(*CI_users[2]->get_value(mc,ind,mu)) + co4*(*CI_users[3]->get_value(mc,ind,mu)) + co5*(*CI_users[4]->get_value(mc,ind,mu)) + co6*(*CI_users[5]->get_value(mc,ind,mu)) + co7*(*CI_users[6]->get_value(mc,ind,mu)) + co8*(*CI_users[7]->get_value(mc,ind,mu)) + co9*(*CI_users[8]->get_value(mc,ind,mu)) + co10*(*CI_users[9]->get_value(mc,ind,mu)) + co11*(*CI_users[10]->get_value(mc,ind,mu)) + co12*(*CI_users[11]->get_value(mc,ind,mu)) + co13*(*CI_users[12]->get_value(mc,ind,mu)) + co14*(*CI_users[13]->get_value(mc,ind,mu)) + co15*(*CI_users[14]->get_value(mc,ind,mu)) + co16*(*CI_users[15]->get_value(mc,ind,mu)) + co17*(*CI_users[16]->get_value(mc,ind,mu)) + co18*(*CI_users[17]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q2l_qmqpQpQmemep_sl_wCI::C2q2Q2l_qmqpQpQmemep_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c134, c256));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c256, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c256));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c56, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c56, c12, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c56, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qmqpQpQmemep_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, em, ep}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmqpQpQmemep sl");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:32:29 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb16 = SPB(1,6);
complex<T> spa25 = SPA(2,5);
complex<T> spb26 = SPB(2,6);
complex<T> spa45 = SPA(4,5);
complex<T> spb36 = SPB(3,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spa35 = SPA(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa56 = SPA(5,6);
complex<T> spa13 = SPA(1,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa16 = SPA(1,6);
complex<T> spb15 = SPB(1,5);
complex<T> spa26 = SPA(2,6);
complex<T> spa46 = SPA(4,6);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = S(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s134 = SS(1,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s15 = -(spa15*spb15);
complex<T> s16 = -(spa16*spb16);
complex<T> s25 = -(spa25*spb25);
complex<T> s26 = -(spa26*spb26);
complex<T> s234 = SS(2,3,4);
complex<T> s256 = SS(2,5,6);
complex<T> t14 = -(spa24*spb34); 
complex<T> t62 = spa15*spb13; 
complex<T> t63 = spa14*spb16; 
complex<T> t64 = spa24*spb26; 
complex<T> t65 = spa25*spb23; 
complex<T> t67 = (spa14*spb13 + spa24*spb23)*(spa15*spb16 + spa25*spb26); 
complex<T> t86 = s12 - s34 - s56; 
complex<T> t90 = -s12 - s34 + s56; 
complex<T> t91 = -s12 + s34 - s56; 
complex<T> t93 = square(spa24); 
complex<T> t94 = square(spa25); 
complex<T> t95 = square(spb13); 
complex<T> t96 = square(spb16); 
complex<T> t99 = s13 + s14 - s23 - s24; 
complex<T> t100 = s15 + s16 - s25 - s26; 
complex<T> t101 = -s156 + s34; 
complex<T> t102 = -s256 + s34; 
complex<T> t103 = square(spa14); 
complex<T> t104 = square(spa15); 
complex<T> t105 = square(spb23); 
complex<T> t106 = square(spb26); 
complex<T> t126 = spa23*spb13; 
complex<T> t127 = spa24*spb14; 
complex<T> t128 = spa25*spb15; 
complex<T> t129 = spa26*spb16; 
complex<T> t131 = s34*s56; 
complex<T> t158 = spa56*spb34; 
complex<T> t173 = spa15*spb13 - spa25*spb23; 
complex<T> t184 = spa34*spb36; 
complex<T> t207 = spa24*spb16; 
complex<T> t208 = spa25*spb13; 
complex<T> t209 = spa34*spb56; 
complex<T> t233 = spa14*spa24; 
complex<T> t234 = spa15*spa25; 
complex<T> t248 = spb13*spb23; 
complex<T> t249 = spb16*spb26; 
complex<T> d17 = T(2); d17 = T(1)/d17;
complex<T> d48 = s134; d48 = T(1)/d48;
complex<T> d49 = s156; d49 = T(1)/d49;
complex<T> t12 = -t67; 
complex<T> t31 = (t126 + t127)*T(2); 
complex<T> t32 = spa34*(t128 + t129); 
complex<T> t34 = square(t184 + t64); 
complex<T> t35 = square(spa56*spb36 + t65); 
complex<T> t36 = (-s15 - s16 + s25 + s26)*spa45 - spa56*t184*T(2); 
complex<T> t37 = (s13 + s14 - s23 - s24)*spa45 + spa56*t184*T(2); 
complex<T> t40 = -square(s12) - square(s34) - square(s56) + s12*s34*T(2) + s12*s56*T(2) + t131*T(2) + s12*t86*T(3); 
complex<T> t46 = spb36*t100 - spa45*spb34*spb56*T(2); 
complex<T> t47 = -(spb36*t99) + spa45*spb34*spb56*T(2); 
complex<T> t52 = d48*t131 + d17*t86; 
complex<T> t53 = d49*t131 + d17*t86; 
complex<T> t66 = -t184; 
complex<T> t70 = (t128 + t129)*T(2); 
complex<T> t73 = -(t91*T(3)); 
complex<T> t92 = spa45*spb34 - t62; 
complex<T> t97 = spa45*spb56 - t63; 
complex<T> t107 = -(spa56*spb36) - t65; 
complex<T> t112 = t63 - t64; 
complex<T> t183 = -t208; 
complex<T> t185 = t67*T(3); 
complex<T> d1 = t126 + t127; d1 = T(1)/d1;
complex<T> d2 = spa56*(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d2 = T(1)/d2;
complex<T> d3 = t128 + t129; d3 = T(1)/d3;
complex<T> d7 = (t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d7 = T(1)/d7;
complex<T> d8 = t158*square(t126 + t127)*T(2); d8 = T(1)/d8;
complex<T> d9 = s134*(t126 + t127)*t158*T(4); d9 = T(1)/d9;
complex<T> d10 = (t128 + t129)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d10 = T(1)/d10;
complex<T> d11 = spb34*(t128 + t129)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d11 = T(1)/d11;
complex<T> d12 = t158*square(t128 + t129)*T(2); d12 = T(1)/d12;
complex<T> d13 = s256*(t128 + t129)*t158*T(4); d13 = T(1)/d13;
complex<T> d14 = square(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d14 = T(1)/d14;
complex<T> d15 = s134*t131; d15 = T(1)/d15;
complex<T> d16 = s234*t131; d16 = T(1)/d16;
complex<T> d18 = spb56*(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d18 = T(1)/d18;
complex<T> d19 = t209*square(t126 + t127)*T(2); d19 = T(1)/d19;
complex<T> d20 = s234*(t126 + t127)*t209*T(4); d20 = T(1)/d20;
complex<T> d21 = square(t128 + t129)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d21 = T(1)/d21;
complex<T> d22 = t209*square(t128 + t129)*T(2); d22 = T(1)/d22;
complex<T> d24 = (-s134 + s56)*t158*square(t126 + t127); d24 = T(1)/d24;
complex<T> d26 = (-s134 + s56)*(t126 + t127)*t158; d26 = T(1)/d26;
complex<T> d27 = t101*t209*square(t128 + t129); d27 = T(1)/d27;
complex<T> d30 = (-s234 + s56)*t209*square(t126 + t127); d30 = T(1)/d30;
complex<T> d32 = (-s234 + s56)*(t126 + t127)*t209; d32 = T(1)/d32;
complex<T> d33 = t102*t158*square(t128 + t129); d33 = T(1)/d33;
complex<T> d35 = t102*(t128 + t129)*t158; d35 = T(1)/d35;
complex<T> d36 = t158*square(t128 + t129); d36 = T(1)/d36;
complex<T> d38 = t209*square(t128 + t129); d38 = T(1)/d38;
complex<T> d40 = t158*square(t126 + t127); d40 = T(1)/d40;
complex<T> d42 = t209*square(t126 + t127); d42 = T(1)/d42;
complex<T> d44 = s134*t158*cube(t126 + t127); d44 = T(1)/d44;
complex<T> d45 = s134*(spa13*spb23 + spa14*spb24)*t209; d45 = T(1)/d45;
complex<T> d46 = s156*(spa15*spb25 + spa16*spb26)*t158; d46 = T(1)/d46;
complex<T> d47 = s156*t209*cube(t128 + t129); d47 = T(1)/d47;
complex<T> d50 = (t126 + t127)*square(s134); d50 = T(1)/d50;
complex<T> d52 = square(t126 + t127); d52 = T(1)/d52;
complex<T> d54 = square(t128 + t129); d54 = T(1)/d54;
complex<T> d55 = (t128 + t129)*square(s156); d55 = T(1)/d55;
complex<T> t10 = -(d46*t104*t105) + d47*square(spa45*spb15 + spa46*spb16)*square(spa12*spb16 - spa25*spb56); 
complex<T> t11 = -(d45*t103*t106) + d44*square(spa35*spb13 + spa45*spb14)*square(spa12*spb13 + spa24*spb34); 
complex<T> t33 = square(t92); 
complex<T> t39 = square(t97); 
complex<T> t42 = -(s134*t62) + s234*t65 + d1*s134*t208*t86; 
complex<T> t43 = s156*t62 - s256*t65 + d3*s256*t208*t86; 
complex<T> t44 = s134*t63 - s234*t64 + d1*s234*t207*t86; 
complex<T> t45 = -(s156*t63) + s256*t64 + d3*s156*t207*t86; 
complex<T> t80 = d26*(-(spa15*spb13) + spa45*spb34); 
complex<T> t81 = -(d19*(spa24*spb26 + t184)); 
complex<T> t82 = -(d12*(spa25*spb23 + spa56*spb36)); 
complex<T> t84 = d22*(spa14*spb16 - spa45*spb56); 
complex<T> t98 = -t64 + t66; 
complex<T> t146 = d35*t107; 
complex<T> t147 = d20*t34; 
complex<T> t195 = d13*t35; 
complex<T> t220 = d8*t92; 
complex<T> d4 = t32*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d4 = T(1)/d4;
complex<T> d5 = t31*square(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d5 = T(1)/d5;
complex<T> d6 = t70*square(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d6 = T(1)/d6;
complex<T> d23 = s156*spb56*t32*T(4); d23 = T(1)/d23;
complex<T> d25 = t158*t31*square(s134 - s56); d25 = T(1)/d25;
complex<T> d28 = spb56*t32*square(t101)*T(2); d28 = T(1)/d28;
complex<T> d29 = spb56*t101*t32; d29 = T(1)/d29;
complex<T> d31 = t209*t31*square(s234 - s56); d31 = T(1)/d31;
complex<T> d34 = t158*t70*square(t102); d34 = T(1)/d34;
complex<T> d37 = s256*t158*t70; d37 = T(1)/d37;
complex<T> d39 = s156*spb56*t32*T(2); d39 = T(1)/d39;
complex<T> d41 = s134*t158*t31; d41 = T(1)/d41;
complex<T> d43 = s234*t209*t31; d43 = T(1)/d43;
complex<T> d51 = t31*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d51 = T(1)/d51;
complex<T> d53 = t70*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t131*T(2)); d53 = T(1)/d53;
complex<T> t4 = -(d53*s12*spa45*spb36*t100) + d6*t100*t12*t40 + d54*t207*(spa56*spb36 + t62) + d21*t207*t86*(spa56*spb36*t90 + spa45*spb34*t91) + d10*t158*square(t63 + t64) - d55*spa15*(spa45*spb15 + spa46*spb16)*spb23*(spa12*spb16 - spa25*spb56)*T(2); 
complex<T> t27 = t208*t82 - d34*s256*t104*t95 + d33*s256*t234*t95 - t146*t62*T(2) - t195*T(3); 
complex<T> t38 = d15*spa14*spb26*t92 - d16*spa15*spb23*t98; 
complex<T> t71 = d5*t99; 
complex<T> t148 = d23*t39; 
complex<T> t162 = d5*(-s13 - s14 + s23 + s24); 
complex<T> t188 = d17*t10; 
complex<T> t194 = d9*t33; 
complex<T> t198 = d29*t97; 
complex<T> t199 = d32*t98; 
complex<T> t212 = d17*t11; 
complex<T> t265 = t65*t80; 
complex<T> t2 = s56*(t188 + t212); 
complex<T> t3 = -t212; 
complex<T> t6 = d52*t183*(-t63 + t66) + t12*t40*t71 + d14*t183*t86*(-(spa45*spb56*t90) + t66*t91) - d51*s12*spa45*spb36*t99 + d7*t209*square(t62 + t65) + d50*spa14*(spa35*spb13 + spa45*spb14)*spb26*(spa12*spb13 - t14)*T(2); 
complex<T> t22 = t183*t220 + d25*s134*t105*t94 - d24*s134*t248*t94 + t265*T(2) + t194*T(3); 
complex<T> t25 = t38*T(3); 
complex<T> t134 = -t188; 
complex<T> t261 = t199*t63; 
complex<T> t262 = t198*t64; 
complex<T> t1 = s12*(t134 + t3 - t38); 
complex<T> t7 = d10*spa45*spb56*t173 + d21*spb56*t183*t36 + d11*spb36*t43 + d4*spa45*t45 - d21*spa56*t207*t46 + d10*spa56*spb36*(-t63 + t64) + d6*t100*t67*t73 + d6*(-s15 - s16 + s25 + s26)*t185*t91 + d40*t208*t92 - d25*s134*t105*t94 + d24*s134*t248*t94 + d31*s234*t103*t96 - d30*s234*t233*t96 - d42*t207*t98 + t261*T(2) - t265*T(2) - d41*t33*T(3) + d43*t34*T(3); 
complex<T> t8 = -(d7*spa45*spb34*t112) - d10*spa45*spb56*t173 + t183*t220 + d17*t25 + d21*spb56*t208*t36 + d14*spb16*t14*t37 - d2*spa45*t42 - d11*spb36*t43 - d18*spb36*t44 - d4*spa45*t45 + d21*spa56*t207*t46 + d14*spa34*t183*t47 + d10*spa56*spb36*(t63 - t64) + d7*(-t62 + t65)*t66 + d6*(-s15 - s16 + s25 + s26)*t67*t73 + t207*t81 + t208*t82 + t207*t84 + t185*t71*t90 + d6*t100*t185*t91 - t147*T(3) + t148*T(3) + t194*T(3) - t195*T(3) - t162*t67*t90*T(3); 
complex<T> t9 = d7*spa45*spb34*t112 + d36*t107*t183 + d17*t25 + d14*spb34*t207*t37 + d2*spa45*t42 + d18*spb36*t44 + d14*spa34*t208*t47 + d7*t184*(-t62 + t65) + t162*t185*t90 - d28*s156*t106*t93 + d27*s156*t249*t93 + d34*s256*t104*t95 - d33*s256*t234*t95 + d38*t207*t97 - t262*T(2) + t146*t62*T(2) + d37*t35*T(3) - d39*t39*T(3) - t67*t71*t90*T(3); 
complex<T> t26 = t207*t81 - d31*s234*t103*t96 + d30*s234*t233*t96 - t261*T(2) - t147*T(3); 
complex<T> t30 = t207*t84 + d28*s156*t106*t93 - d27*s156*t249*t93 + t262*T(2) + t148*T(3); 
complex<T> t181 = s34*(t134 + t3 - t38); 
complex<T> co1 = Complex(0,1)*t8; 
complex<T> co2 = Complex(0,1)*t22; 
complex<T> co3 = Complex(0,1)*t30; 
complex<T> co4 = Complex(0,1)*t26; 
complex<T> co5 = Complex(0,1)*t27; 
complex<T> co6 = Complex(0,1)*t9; 
complex<T> co7 = Complex(0,1)*t7; 
complex<T> co8 = Complex(0,1)*t1; 
complex<T> co9 = Complex(0,1)*t102*t11; 
complex<T> co10 = Complex(0,-1)*t11*t52; 
complex<T> co11 = Complex(0,-1)*t10*t53; 
complex<T> co12 = Complex(0,1)*t10*t101; 
complex<T> co13 = Complex(0,1)*t181; 
complex<T> co14 = Complex(0,-1)*t6; 
complex<T> co15 = Complex(0,1)*t2; 
complex<T> co16 = Complex(0,-1)*t4; 
complex<T> co17 = Complex(0,-1)*s12*s134*t212; 
complex<T> co18 = Complex(0,-1)*s12*s156*t188; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + co3*(*CI_users[2]->get_value(mc,ind,mu)) + co4*(*CI_users[3]->get_value(mc,ind,mu)) + co5*(*CI_users[4]->get_value(mc,ind,mu)) + co6*(*CI_users[5]->get_value(mc,ind,mu)) + co7*(*CI_users[6]->get_value(mc,ind,mu)) + co8*(*CI_users[7]->get_value(mc,ind,mu)) + co9*(*CI_users[8]->get_value(mc,ind,mu)) + co10*(*CI_users[9]->get_value(mc,ind,mu)) + co11*(*CI_users[10]->get_value(mc,ind,mu)) + co12*(*CI_users[11]->get_value(mc,ind,mu)) + co13*(*CI_users[12]->get_value(mc,ind,mu)) + co14*(*CI_users[13]->get_value(mc,ind,mu)) + co15*(*CI_users[14]->get_value(mc,ind,mu)) + co16*(*CI_users[15]->get_value(mc,ind,mu)) + co17*(*CI_users[16]->get_value(mc,ind,mu)) + co18*(*CI_users[17]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q2l_qmqpQmQpemep_sl_wCI::C2q2Q2l_qmqpQmQpemep_sl_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c134, c256));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c256, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c256));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c56, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c5, c6, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c56, c12, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c34, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c56, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qmqpQmQpemep_sl_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, em, ep}, sl}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmqpQmQpemep sl");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:35:08 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb16 = SPB(1,6);
complex<T> spa25 = SPA(2,5);
complex<T> spb26 = SPB(2,6);
complex<T> spa35 = SPA(3,5);
complex<T> spb46 = SPB(4,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb56 = SPB(5,6);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa56 = SPA(5,6);
complex<T> spb23 = SPB(2,3);
complex<T> spa14 = SPA(1,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa16 = SPA(1,6);
complex<T> spb15 = SPB(1,5);
complex<T> spa26 = SPA(2,6);
complex<T> spa36 = SPA(3,6);
complex<T> s13 = -(spa13*spb13);
complex<T> s14 = -(spa14*spb14);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = S(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = -(spa56*spb56);
complex<T> s134 = SS(1,3,4);
complex<T> s156 = SS(1,5,6);
complex<T> s15 = -(spa15*spb15);
complex<T> s16 = -(spa16*spb16);
complex<T> s25 = -(spa25*spb25);
complex<T> s26 = -(spa26*spb26);
complex<T> s234 = SS(2,3,4);
complex<T> s256 = SS(2,5,6);
complex<T> t36 = (s13 + s14 - s23 - s24)*spa35 - spa34*spa56*spb46*T(2); 
complex<T> t37 = (-s15 - s16 + s25 + s26)*spa35 + spa34*spa56*spb46*T(2); 
complex<T> t62 = spa15*spb14; 
complex<T> t63 = spa13*spb16; 
complex<T> t64 = spa23*spb26; 
complex<T> t65 = spa25*spb24; 
complex<T> t67 = (spa13*spb14 + spa23*spb24)*(spa15*spb16 + spa25*spb26); 
complex<T> t86 = s12 - s34 - s56; 
complex<T> t90 = -s12 - s34 + s56; 
complex<T> t91 = -s12 + s34 - s56; 
complex<T> t93 = square(spa23); 
complex<T> t94 = square(spa25); 
complex<T> t95 = square(spb14); 
complex<T> t96 = square(spb16); 
complex<T> t99 = s13 + s14 - s23 - s24; 
complex<T> t100 = s15 + s16 - s25 - s26; 
complex<T> t101 = -s156 + s34; 
complex<T> t102 = -s256 + s34; 
complex<T> t103 = square(spa13); 
complex<T> t104 = square(spa15); 
complex<T> t105 = square(spb24); 
complex<T> t106 = square(spb26); 
complex<T> t126 = spa23*spb13; 
complex<T> t127 = spa24*spb14; 
complex<T> t128 = spa25*spb15; 
complex<T> t129 = spa26*spb16; 
complex<T> t131 = spa35*spb34; 
complex<T> t132 = s34*s56; 
complex<T> t158 = spa56*spb34; 
complex<T> t166 = spa34*spb56; 
complex<T> t174 = spa15*spb14 - spa25*spb24; 
complex<T> t208 = spa23*spb16; 
complex<T> t209 = spa25*spb14; 
complex<T> t234 = spa13*spa23; 
complex<T> t235 = spa15*spa25; 
complex<T> t239 = spb14*spb24; 
complex<T> t250 = spb16*spb26; 
complex<T> d17 = T(2); d17 = T(1)/d17;
complex<T> d48 = s134; d48 = T(1)/d48;
complex<T> d49 = s156; d49 = T(1)/d49;
complex<T> t12 = -t67; 
complex<T> t32 = spa34*(t128 + t129); 
complex<T> t33 = square(t131 + t62); 
complex<T> t35 = square(spa56*spb46 + t65); 
complex<T> t40 = -square(s12) - square(s34) - square(s56) + s12*s34*T(2) + s12*s56*T(2) + t132*T(2) + s12*t86*T(3); 
complex<T> t46 = -(spb46*t99) - spb56*t131*T(2); 
complex<T> t47 = spb46*t100 + spb56*t131*T(2); 
complex<T> t52 = d48*t132 + d17*t86; 
complex<T> t53 = d49*t132 + d17*t86; 
complex<T> t68 = -t131; 
complex<T> t70 = (t128 + t129)*T(2); 
complex<T> t73 = -(t91*T(3)); 
complex<T> t92 = spa34*spb46 - t64; 
complex<T> t97 = spa35*spb56 - t63; 
complex<T> t107 = -(spa56*spb46) - t65; 
complex<T> t113 = t63 - t64; 
complex<T> t114 = -t62 + t65; 
complex<T> t169 = -t63 + t64; 
complex<T> t181 = -t209; 
complex<T> t183 = t67*T(3); 
complex<T> d1 = t126 + t127; d1 = T(1)/d1;
complex<T> d2 = spa56*(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d2 = T(1)/d2;
complex<T> d3 = t128 + t129; d3 = T(1)/d3;
complex<T> d5 = (t126 + t127)*square(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2))*T(2); d5 = T(1)/d5;
complex<T> d7 = (t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d7 = T(1)/d7;
complex<T> d8 = t158*square(t126 + t127)*T(2); d8 = T(1)/d8;
complex<T> d9 = s134*(t126 + t127)*t158*T(4); d9 = T(1)/d9;
complex<T> d10 = (t128 + t129)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d10 = T(1)/d10;
complex<T> d11 = spb34*(t128 + t129)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d11 = T(1)/d11;
complex<T> d12 = t158*square(t128 + t129)*T(2); d12 = T(1)/d12;
complex<T> d13 = s256*(t128 + t129)*t158*T(4); d13 = T(1)/d13;
complex<T> d14 = square(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d14 = T(1)/d14;
complex<T> d15 = s134*t132; d15 = T(1)/d15;
complex<T> d16 = s234*t132; d16 = T(1)/d16;
complex<T> d18 = spb56*(t126 + t127)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d18 = T(1)/d18;
complex<T> d19 = t166*square(t126 + t127)*T(2); d19 = T(1)/d19;
complex<T> d20 = s234*(t126 + t127)*t166*T(4); d20 = T(1)/d20;
complex<T> d21 = square(t128 + t129)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d21 = T(1)/d21;
complex<T> d22 = t166*square(t128 + t129)*T(2); d22 = T(1)/d22;
complex<T> d24 = (-s134 + s56)*t158*square(t126 + t127); d24 = T(1)/d24;
complex<T> d25 = (t126 + t127)*t158*square(s134 - s56)*T(2); d25 = T(1)/d25;
complex<T> d26 = (-s134 + s56)*(t126 + t127)*t158; d26 = T(1)/d26;
complex<T> d27 = t101*t166*square(t128 + t129); d27 = T(1)/d27;
complex<T> d30 = (-s234 + s56)*t166*square(t126 + t127); d30 = T(1)/d30;
complex<T> d31 = (t126 + t127)*t166*square(s234 - s56)*T(2); d31 = T(1)/d31;
complex<T> d32 = (-s234 + s56)*(t126 + t127)*t166; d32 = T(1)/d32;
complex<T> d33 = t102*t158*square(t128 + t129); d33 = T(1)/d33;
complex<T> d35 = t102*(t128 + t129)*t158; d35 = T(1)/d35;
complex<T> d36 = t158*square(t128 + t129); d36 = T(1)/d36;
complex<T> d38 = t166*square(t128 + t129); d38 = T(1)/d38;
complex<T> d40 = t158*square(t126 + t127); d40 = T(1)/d40;
complex<T> d41 = s134*(t126 + t127)*t158*T(2); d41 = T(1)/d41;
complex<T> d42 = t166*square(t126 + t127); d42 = T(1)/d42;
complex<T> d43 = s234*(t126 + t127)*t166*T(2); d43 = T(1)/d43;
complex<T> d44 = s134*t158*cube(t126 + t127); d44 = T(1)/d44;
complex<T> d45 = s134*(spa13*spb23 + spa14*spb24)*t166; d45 = T(1)/d45;
complex<T> d46 = s156*(spa15*spb25 + spa16*spb26)*t158; d46 = T(1)/d46;
complex<T> d47 = s156*t166*cube(t128 + t129); d47 = T(1)/d47;
complex<T> d50 = (t126 + t127)*square(s134); d50 = T(1)/d50;
complex<T> d51 = (t126 + t127)*T(2)*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d51 = T(1)/d51;
complex<T> d52 = square(t126 + t127); d52 = T(1)/d52;
complex<T> d54 = square(t128 + t129); d54 = T(1)/d54;
complex<T> d55 = (t128 + t129)*square(s156); d55 = T(1)/d55;
complex<T> t10 = d45*t103*t106 - d44*square(spa35*spb13 + spa45*spb14)*square(spa12*spb14 - spa23*spb34); 
complex<T> t11 = d46*t104*t105 - d47*square(spa35*spb15 + spa36*spb16)*square(spa12*spb16 - spa25*spb56); 
complex<T> t34 = square(t92); 
complex<T> t39 = square(t97); 
complex<T> t42 = -(s134*t62) + s234*t65 + d1*s134*t209*t86; 
complex<T> t43 = s156*t62 - s256*t65 + d3*s256*t209*t86; 
complex<T> t44 = s134*t63 - s234*t64 + d1*s234*t208*t86; 
complex<T> t45 = -(s156*t63) + s256*t64 + d3*s156*t208*t86; 
complex<T> t71 = d5*t99; 
complex<T> t80 = -(d26*(spa15*spb14 + t131)); 
complex<T> t81 = d19*(-(spa23*spb26) + spa34*spb46); 
complex<T> t82 = -(d12*(spa25*spb24 + spa56*spb46)); 
complex<T> t84 = d22*(spa13*spb16 - spa35*spb56); 
complex<T> t98 = -t62 + t68; 
complex<T> t144 = d35*t107; 
complex<T> t162 = d5*(-s13 - s14 + s23 + s24); 
complex<T> t194 = d9*t33; 
complex<T> t195 = d13*t35; 
complex<T> t198 = d32*t92; 
complex<T> d4 = t32*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d4 = T(1)/d4;
complex<T> d6 = t70*square(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d6 = T(1)/d6;
complex<T> d23 = s156*spb56*t32*T(4); d23 = T(1)/d23;
complex<T> d28 = spb56*t32*square(t101)*T(2); d28 = T(1)/d28;
complex<T> d29 = spb56*t101*t32; d29 = T(1)/d29;
complex<T> d34 = t158*t70*square(t102); d34 = T(1)/d34;
complex<T> d37 = s256*t158*t70; d37 = T(1)/d37;
complex<T> d39 = s156*spb56*t32*T(2); d39 = T(1)/d39;
complex<T> d53 = t70*(square(s12) + square(s34) + square(s56) - s12*(s34 + s56)*T(2) - t132*T(2)); d53 = T(1)/d53;
complex<T> t4 = -(d53*s12*spa35*spb46*t100) + d6*t100*t12*t40 + d54*t208*(spa56*spb46 + t62) + d21*t208*t86*(spa56*spb46*t90 - t131*t91) - d10*t158*square(t63 + t64) - d55*spa15*(spa35*spb15 + spa36*spb16)*spb24*(spa12*spb16 - spa25*spb56)*T(2); 
complex<T> t5 = s12*t11; 
complex<T> t6 = d52*t181*(spa34*spb46 - t63) + t12*t40*t71 + d14*t181*t86*(-(spa35*spb56*t90) + spa34*spb46*t91) - d51*s12*spa35*spb46*t99 - d7*t166*square(t62 + t65) + d50*spa13*(spa35*spb13 + spa45*spb14)*spb26*(spa12*spb14 - spa23*spb34)*T(2); 
complex<T> t27 = t209*t82 - d34*s256*t104*t95 + d33*s256*t235*t95 - t144*t62*T(2) - t195*T(3); 
complex<T> t38 = -(d16*spa15*spb24*t92) + d15*spa13*spb26*t98; 
complex<T> t134 = d17*t10; 
complex<T> t146 = d20*t34; 
complex<T> t147 = d23*t39; 
complex<T> t199 = d29*t97; 
complex<T> t232 = d8*t98; 
complex<T> t261 = t198*t63; 
complex<T> t263 = t65*t80; 
complex<T> t2 = -(s56*(d17*t11 + t134)); 
complex<T> t9 = -(d10*spa56*spb46*t169) - d10*spa35*spb56*t174 + d21*spb56*t209*t37 + d11*spb46*t43 + d4*spa35*t45 + d21*spa56*t208*t47 + d6*(-s15 - s16 + s25 + s26)*t67*t73 + d6*t100*t183*t91 - d42*t208*t92 - d25*s134*t105*t94 + d24*s134*t239*t94 + d31*s234*t103*t96 - d30*s234*t234*t96 + d40*t209*t98 + t261*T(2) - t263*T(2) - d41*t33*T(3) + d43*t34*T(3); 
complex<T> t22 = t181*t232 + d25*s134*t105*t94 - d24*s134*t239*t94 + t263*T(2) + t194*T(3); 
complex<T> t25 = -(t38*T(3)); 
complex<T> t26 = t208*t81 - d31*s234*t103*t96 + d30*s234*t234*t96 - t261*T(2) - t146*T(3); 
complex<T> t179 = s34*(d17*t11 + t134 + t38); 
complex<T> t240 = s12*t134; 
complex<T> t247 = d17*t5; 
complex<T> t262 = t199*t64; 
complex<T> t1 = t240 + t247 + s12*t38; 
complex<T> t7 = d7*spa34*spb46*t114 + d7*t113*t131 + d36*t107*t181 + d17*t25 + d14*spb34*t208*t36 - d2*spa35*t42 - d18*spb46*t44 + d14*spa34*t209*t46 + t183*t71*t90 - d28*s156*t106*t93 + d27*s156*t250*t93 + d34*s256*t104*t95 - d33*s256*t235*t95 + d38*t208*t97 - t262*T(2) + t144*t62*T(2) + d37*t35*T(3) - d39*t39*T(3) - t162*t67*t90*T(3); 
complex<T> t8 = -(d7*spa34*spb46*t114) + d10*spa56*spb46*t169 + d10*spa35*spb56*t174 + t181*t232 + d17*t25 - d14*spb34*t208*t36 + d21*spb56*t181*t37 + d2*spa35*t42 - d11*spb46*t43 + d18*spb46*t44 - d4*spa35*t45 + d14*spa34*t181*t46 - d21*spa56*t208*t47 + d7*t113*t68 + d6*t100*t67*t73 + t208*t81 + t209*t82 + t208*t84 + t162*t183*t90 + d6*(-s15 - s16 + s25 + s26)*t183*t91 - t146*T(3) + t147*T(3) + t194*T(3) - t195*T(3) - t67*t71*t90*T(3); 
complex<T> t30 = t208*t84 + d28*s156*t106*t93 - d27*s156*t250*t93 + t262*T(2) + t147*T(3); 
complex<T> co1 = Complex(0,1)*t8; 
complex<T> co2 = Complex(0,1)*t22; 
complex<T> co3 = Complex(0,1)*t30; 
complex<T> co4 = Complex(0,1)*t26; 
complex<T> co5 = Complex(0,1)*t27; 
complex<T> co6 = Complex(0,1)*t7; 
complex<T> co7 = Complex(0,1)*t9; 
complex<T> co8 = Complex(0,1)*t1; 
complex<T> co9 = Complex(0,-1)*t10*t102; 
complex<T> co10 = Complex(0,1)*t10*t52; 
complex<T> co11 = Complex(0,1)*t11*t53; 
complex<T> co12 = Complex(0,-1)*t101*t11; 
complex<T> co13 = Complex(0,1)*t179; 
complex<T> co14 = Complex(0,1)*t6; 
complex<T> co15 = Complex(0,1)*t2; 
complex<T> co16 = Complex(0,1)*t4; 
complex<T> co17 = Complex(0,1)*s134*t240; 
complex<T> co18 = Complex(0,1)*s156*t247; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + co3*(*CI_users[2]->get_value(mc,ind,mu)) + co4*(*CI_users[3]->get_value(mc,ind,mu)) + co5*(*CI_users[4]->get_value(mc,ind,mu)) + co6*(*CI_users[5]->get_value(mc,ind,mu)) + co7*(*CI_users[6]->get_value(mc,ind,mu)) + co8*(*CI_users[7]->get_value(mc,ind,mu)) + co9*(*CI_users[8]->get_value(mc,ind,mu)) + co10*(*CI_users[9]->get_value(mc,ind,mu)) + co11*(*CI_users[10]->get_value(mc,ind,mu)) + co12*(*CI_users[11]->get_value(mc,ind,mu)) + co13*(*CI_users[12]->get_value(mc,ind,mu)) + co14*(*CI_users[13]->get_value(mc,ind,mu)) + co15*(*CI_users[14]->get_value(mc,ind,mu)) + co16*(*CI_users[15]->get_value(mc,ind,mu)) + co17*(*CI_users[16]->get_value(mc,ind,mu)) + co18*(*CI_users[17]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q2l_qpqmQpQmemep_AX_wCI::C2q2Q2l_qpqmQpQmemep_AX_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qpqmQpQmemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, Qp, Qm, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpqmQpQmemep AX");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:35:10 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb16 = SPB(1,6);
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb36 = SPB(3,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = S(5,6);
complex<T> t2 = -(s12*T(2)); 
complex<T> t12 = -(s34*T(2)); 
complex<T> t13 = (s12 - s34 - s56)*T(3); 
complex<T> t24 = s34*T(2); 
complex<T> d3 = spb12; d3 = T(1)/d3;
complex<T> d4 = spa12; d4 = T(1)/d4;
complex<T> d5 = s56; d5 = T(1)/d5;
complex<T> d6 = spa34; d6 = T(1)/d6;
complex<T> d7 = spb34; d7 = T(1)/d7;
complex<T> t17 = d3*spa45*spb13*spb16 + d4*spa24*spa25*spb36; 
complex<T> t18 = -(d6*spa24*spa45*spb16) - d7*spa25*spb13*spb36; 
complex<T> d1 = square(s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34)); d1 = T(1)/d1;
complex<T> d2 = (s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34))*T(2); d2 = T(1)/d2;
complex<T> d8 = s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34); d8 = T(1)/d8;
complex<T> t3 = d5*(-d2 + d1*s56*(-s12 - s34 + s56)*T(3)); 
complex<T> t14 = -(d1*(s12 - s34 + s56)); 
complex<T> t15 = s12*t17; 
complex<T> t28 = d1*t13; 
complex<T> t4 = t17*(-d2 + s12*t28); 
complex<T> t7 = -d8 + s12*t28; 
complex<T> t16 = t18*t3; 
complex<T> t33 = t14*T(3); 
complex<T> t37 = t15*T(2); 
complex<T> t5 = t18*(-d2 + s34*t33); 
complex<T> t9 = -d8 + s34*t33; 
complex<T> t11 = t12*t16 + t3*t37 + t4*T(2); 
complex<T> t1 = -((t4 + t5)*T(2)); 
complex<T> t22 = t16*t24 + t17*t2*t3 + t5*T(2); 
complex<T> co1 = t37*t9; 
complex<T> co2 = t18*t24*t7; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t22*(*CI_users[0]->get_value(mc,ind,mu)) + t11*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q2l_qmqpQpQmemep_AX_wCI::C2q2Q2l_qmqpQpQmemep_AX_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qmqpQpQmemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmqpQpQmemep AX");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:35:12 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb26 = SPB(2,6);
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb36 = SPB(3,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = S(5,6);
complex<T> t2 = -(s12*T(2)); 
complex<T> t12 = -(s34*T(2)); 
complex<T> t13 = (s12 - s34 - s56)*T(3); 
complex<T> d3 = spb12; d3 = T(1)/d3;
complex<T> d4 = spa12; d4 = T(1)/d4;
complex<T> d5 = s56; d5 = T(1)/d5;
complex<T> d6 = spa34; d6 = T(1)/d6;
complex<T> d7 = spb34; d7 = T(1)/d7;
complex<T> t16 = -(d6*spa14*spa45*spb26) - d7*spa15*spb23*spb36; 
complex<T> t17 = -(d3*spa45*spb23*spb26) - d4*spa14*spa15*spb36; 
complex<T> d1 = square(s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34)); d1 = T(1)/d1;
complex<T> d2 = (s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34))*T(2); d2 = T(1)/d2;
complex<T> d8 = s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34); d8 = T(1)/d8;
complex<T> t3 = d5*(-d2 + d1*s56*(-s12 - s34 + s56)*T(3)); 
complex<T> t14 = -(d1*(s12 - s34 + s56)); 
complex<T> t27 = d1*t13; 
complex<T> t4 = t17*(-d2 + s12*t27); 
complex<T> t7 = -d8 + s12*t27; 
complex<T> t15 = t17*t3; 
complex<T> t31 = t14*T(3); 
complex<T> t5 = t16*(-d2 + s34*t31); 
complex<T> t9 = -d8 + s34*t31; 
complex<T> t11 = t15*t2 + s34*t16*t3*T(2) - t4*T(2); 
complex<T> t1 = (t4 + t5)*T(2); 
complex<T> t21 = t12*t16*t3 + s12*t15*T(2) - t5*T(2); 
complex<T> co1 = t17*t2*t9; 
complex<T> co2 = t12*t16*t7; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t21*(*CI_users[0]->get_value(mc,ind,mu)) + t11*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q2l_qpqmQmQpemep_AX_wCI::C2q2Q2l_qpqmQmQpemep_AX_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qpqmQmQpemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qp, qm, Qm, Qp, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qpqmQmQpemep AX");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:35:15 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb16 = SPB(1,6);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb46 = SPB(4,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = S(5,6);
complex<T> t2 = -(s12*T(2)); 
complex<T> t12 = -(s34*T(2)); 
complex<T> t13 = (s12 - s34 - s56)*T(3); 
complex<T> d3 = spb12; d3 = T(1)/d3;
complex<T> d4 = spa12; d4 = T(1)/d4;
complex<T> d5 = s56; d5 = T(1)/d5;
complex<T> d6 = spa34; d6 = T(1)/d6;
complex<T> d7 = spb34; d7 = T(1)/d7;
complex<T> t16 = d6*spa23*spa35*spb16 + d7*spa25*spb14*spb46; 
complex<T> t17 = d3*spa35*spb14*spb16 + d4*spa23*spa25*spb46; 
complex<T> d1 = square(s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34)); d1 = T(1)/d1;
complex<T> d2 = (s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34))*T(2); d2 = T(1)/d2;
complex<T> d8 = s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34); d8 = T(1)/d8;
complex<T> t3 = d5*(-d2 + d1*s56*(-s12 - s34 + s56)*T(3)); 
complex<T> t14 = -(d1*(s12 - s34 + s56)); 
complex<T> t27 = d1*t13; 
complex<T> t4 = t17*(-d2 + s12*t27); 
complex<T> t7 = -d8 + s12*t27; 
complex<T> t15 = t17*t3; 
complex<T> t31 = t14*T(3); 
complex<T> t5 = t16*(-d2 + s34*t31); 
complex<T> t9 = -d8 + s34*t31; 
complex<T> t11 = t15*t2 + s34*t16*t3*T(2) - t4*T(2); 
complex<T> t1 = (t4 + t5)*T(2); 
complex<T> t21 = t12*t16*t3 + s12*t15*T(2) - t5*T(2); 
complex<T> co1 = t17*t2*t9; 
complex<T> co2 = t12*t16*t7; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t21*(*CI_users[0]->get_value(mc,ind,mu)) + t11*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q2l_qmqpQmQpemep_AX_wCI::C2q2Q2l_qmqpQmQpemep_AX_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c56, c1234));
CI_users.push_back(new Cached_Triangle_Integral_User(c12, c34, c56));
CI_users.push_back(new Cached_Triangle_Integral_User(c34, c12, c56));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q2l_qmqpQmQpemep_AX_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("C2q2Q2l :  qmqpQmQpemep AX");
#endif
 
//#define TimeStamp "Thu 24 Sep 2009 22:35:17 on n2001"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb26 = SPB(2,6);
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb46 = SPB(4,6);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s56 = S(5,6);
complex<T> t2 = -(s12*T(2)); 
complex<T> t12 = -(s34*T(2)); 
complex<T> t13 = (s12 - s34 - s56)*T(3); 
complex<T> t24 = s34*T(2); 
complex<T> d3 = spb12; d3 = T(1)/d3;
complex<T> d4 = spa12; d4 = T(1)/d4;
complex<T> d5 = s56; d5 = T(1)/d5;
complex<T> d6 = spa34; d6 = T(1)/d6;
complex<T> d7 = spb34; d7 = T(1)/d7;
complex<T> t17 = -(d3*spa35*spb24*spb26) - d4*spa13*spa15*spb46; 
complex<T> t18 = d6*spa13*spa35*spb26 + d7*spa15*spb24*spb46; 
complex<T> d1 = square(s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34)); d1 = T(1)/d1;
complex<T> d2 = (s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34))*T(2); d2 = T(1)/d2;
complex<T> d8 = s34*t2 + s56*(s56 + t12 + t2) + square(s12) + square(s34); d8 = T(1)/d8;
complex<T> t3 = d5*(-d2 + d1*s56*(-s12 - s34 + s56)*T(3)); 
complex<T> t14 = -(d1*(s12 - s34 + s56)); 
complex<T> t15 = s12*t17; 
complex<T> t28 = d1*t13; 
complex<T> t4 = t17*(-d2 + s12*t28); 
complex<T> t7 = -d8 + s12*t28; 
complex<T> t16 = t18*t3; 
complex<T> t33 = t14*T(3); 
complex<T> t37 = t15*T(2); 
complex<T> t5 = t18*(-d2 + s34*t33); 
complex<T> t9 = -d8 + s34*t33; 
complex<T> t11 = t12*t16 + t3*t37 + t4*T(2); 
complex<T> t1 = -((t4 + t5)*T(2)); 
complex<T> t22 = t16*t24 + t17*t2*t3 + t5*T(2); 
complex<T> co1 = t37*t9; 
complex<T> co2 = t18*t24*t7; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t22*(*CI_users[0]->get_value(mc,ind,mu)) + t11*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qpQmQpqmemep_nf C2q2Q2l_254161_nf
#define _C_qpQpQmqmemep_nf C2q2Q2l_254105_nf
#define _C_qmQpQmqpemep_nf C2q2Q2l_254616_nf
#define _C_qmQmQpqpemep_nf C2q2Q2l_254672_nf
#define _C_qpQmQpqmemep_nf_top C2q2Q2l_254161_nf_top
#define _C_qpQpQmqmemep_nf_top C2q2Q2l_254105_nf_top
#define _C_qmQpQmqpemep_nf_top C2q2Q2l_254616_nf_top
#define _C_qmQmQpqpemep_nf_top C2q2Q2l_254672_nf_top
#define _C_qpQmQpqmemep_L C2q2Q2l_254161_L
#define _C_qpQpQmqmemep_L C2q2Q2l_254105_L
#define _C_qmQpQmqpemep_L C2q2Q2l_254616_L
#define _C_qmQmQpqpemep_L C2q2Q2l_254672_L
#define _C_qpqmQpQmemep_sl C2q2Q2l_255169_sl
#define _C_qpqmQmQpemep_sl C2q2Q2l_255617_sl
#define _C_qmqpQpQmemep_sl C2q2Q2l_255176_sl
#define _C_qmqpQmQpemep_sl C2q2Q2l_255624_sl
#define _C_qpqmQpQmemep_AX C2q2Q2l_255169_AX
#define _C_qmqpQpQmemep_AX C2q2Q2l_255176_AX
#define _C_qpqmQmQpemep_AX C2q2Q2l_255617_AX
#define _C_qmqpQmQpemep_AX C2q2Q2l_255624_AX
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qpQmQpqmemep_nf case 254161
 
#define _CASE_qpQpQmqmemep_nf case 254105
 
#define _CASE_qmQpQmqpemep_nf case 254616
 
#define _CASE_qmQmQpqpemep_nf case 254672
 
#define _CASE_qpQmQpqmemep_nf_top case 254161
 
#define _CASE_qpQpQmqmemep_nf_top case 254105
 
#define _CASE_qmQpQmqpemep_nf_top case 254616
 
#define _CASE_qmQmQpqpemep_nf_top case 254672
 
#define _CASE_qpQmQpqmemep_L case 254161
 
#define _CASE_qpQpQmqmemep_L case 254105
 
#define _CASE_qmQpQmqpemep_L case 254616
 
#define _CASE_qmQmQpqpemep_L case 254672
 
#define _CASE_qpqmQpQmemep_sl case 255169
 
#define _CASE_qpqmQmQpemep_sl case 255617
 
#define _CASE_qmqpQpQmemep_sl case 255176
 
#define _CASE_qmqpQmQpemep_sl case 255624
 
#define _CASE_qpqmQpQmemep_AX case 255169
 
#define _CASE_qmqpQpQmemep_AX case 255176
 
#define _CASE_qpqmQmQpemep_AX case 255617
 
#define _CASE_qmqpQmQpemep_AX case 255624
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2Q2l_AX( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qpqmQpQmemep_AX: return new 
                       C2q2Q2l_qpqmQpQmemep_AX_wCI(ind);
    _CASE_qmqpQpQmemep_AX: return new 
                       C2q2Q2l_qmqpQpQmemep_AX_wCI(ind);
    _CASE_qpqmQmQpemep_AX: return new 
                       C2q2Q2l_qpqmQmQpemep_AX_wCI(ind);
    _CASE_qmqpQmQpemep_AX: return new 
                       C2q2Q2l_qmqpQmQpemep_AX_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q2l_L( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qpQmQpqmemep_L: return new 
                       C2q2Q2l_qpQmQpqmemep_L_wCI(ind);
    _CASE_qpQpQmqmemep_L: return new 
                       C2q2Q2l_qpQpQmqmemep_L_wCI(ind);
    _CASE_qmQpQmqpemep_L: return new 
                       C2q2Q2l_qmQpQmqpemep_L_wCI(ind);
    _CASE_qmQmQpqpemep_L: return new 
                       C2q2Q2l_qmQmQpqpemep_L_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q2l_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qpQmQpqmemep_nf: return new 
                       C2q2Q2l_qpQmQpqmemep_nf_wCI(ind);
    _CASE_qpQpQmqmemep_nf: return new 
                       C2q2Q2l_qpQpQmqmemep_nf_wCI(ind);
    _CASE_qmQpQmqpemep_nf: return new 
                       C2q2Q2l_qmQpQmqpemep_nf_wCI(ind);
    _CASE_qmQmQpqpemep_nf: return new 
                       C2q2Q2l_qmQmQpqpemep_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q2l_nf_top( int hc,const std::vector<int>& 
ind) { 
    switch (hc) {
    _CASE_qpQmQpqmemep_nf_top: return new 
                       C2q2Q2l_qpQmQpqmemep_nf_top_wCI(ind);
    _CASE_qpQpQmqmemep_nf_top: return new 
                       C2q2Q2l_qpQpQmqmemep_nf_top_wCI(ind);
    _CASE_qmQpQmqpemep_nf_top: return new 
                       C2q2Q2l_qmQpQmqpemep_nf_top_wCI(ind);
    _CASE_qmQmQpqpemep_nf_top: return new 
                       C2q2Q2l_qmQmQpqpemep_nf_top_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q2l_sl( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qpqmQpQmemep_sl: return new 
                       C2q2Q2l_qpqmQpQmemep_sl_wCI(ind);
    _CASE_qpqmQmQpemep_sl: return new 
                       C2q2Q2l_qpqmQmQpemep_sl_wCI(ind);
    _CASE_qmqpQpQmemep_sl: return new 
                       C2q2Q2l_qmqpQpQmemep_sl_wCI(ind);
    _CASE_qmqpQmQpemep_sl: return new 
                       C2q2Q2l_qmqpQmQpemep_sl_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
