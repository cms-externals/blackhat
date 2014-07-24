/*
*C_2q2g_wCI.cpp
*
* Created on 11/12, 2010
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
 



 
class C2q2g_qmpmqp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmpmqp_LT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmmpqp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmmpqp_LT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmpqpm_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmpqpm_LT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmmqpp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmmqpp_LT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmqppm_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmqppm_LT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmqpmp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmqpmp_LT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmpqpm_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmpqpm_nfLT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmmqpp_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmmqpp_nfLT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmqppm_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmqppm_nfLT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};

 
class C2q2g_qmqpmp_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g_qmqpmp_nfLT_wCI
      (const std::vector<int>&);
       SeriesC<R> eval(momentum_configuration<R>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#if BH_USE_GMP
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,
               const Index_Vector& ind,int mu)
        {return eval_fn(mc,ind,mu);};
#endif
private:
       template <class T> SeriesC<T>
          eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,
               int mu);
};


C2q2g_qmpmqp_LT_wCI::\
C2q2g_qmpmqp_LT_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
         vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
         vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
	 vector<int> c41;  c41.push_back(ind[4-1]); c41.push_back(ind[1-1]);
         vector<int> c14 = c41;
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c34));
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c23));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c4, c1, c2, c3));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g_qmpmqp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, m, qp}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g :  qmpmqp LT");
#endif
 
//#define TimeStamp "Fri 12 Nov 2010 14:03:40 on n48"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> s14 = -(spa14*spb14);
complex<T> s24 = -(spa24*spb24);
complex<T> t2 = cube(spa13); 
complex<T> t5 = cube(s14); 
complex<T> t6 = -(spb12*spb14); 
complex<T> t7 = spa23*T(2); 
complex<T> t9 = cube(spa24); 
complex<T> t10 = cube(spb24); 
complex<T> t11 = square(s14); 
complex<T> d3 = spa23*square(spa24)*square(spb24); d3 = T(1)/d3;
complex<T> d5 = spa12*spa23; d5 = T(1)/d5;
complex<T> t12 = spb12*t2; 
complex<T> t16 = t2*t6; 
complex<T> d1 = spa12*spa14*t7; d1 = T(1)/d1;
complex<T> d2 = s24*spa14*t7; d2 = T(1)/d2;
complex<T> d4 = spa23*t10*t9; d4 = T(1)/d4;
complex<T> d6 = spa12*spa23*t10*t9; d6 = T(1)/d6;
complex<T> d7 = t10*t7*t9; d7 = T(1)/d7;
complex<T> t1 = spb14*t2*(-d5 + d6*t5); 
complex<T> t8 = -(d2*T(3)); 
complex<T> t13 = d2*T(3); 
complex<T> t15 = d3*spb14*t12 + t12*t13 - d1*t2*T(3); 
complex<T> t18 = d3*t16 + t12*t8; 
complex<T> co1 = d4*t11*t16; 
complex<T> co2 = d7*t16*t5; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t15*(*CI_users[0]->get_value(mc,ind,mu)) + t18*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g_qmmpqp_LT_wCI::\
C2q2g_qmmpqp_LT_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
         vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
         vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
	 vector<int> c41;  c41.push_back(ind[4-1]); c41.push_back(ind[1-1]);
         vector<int> c14 = c41;
CI_users.push_back(new Cached_Bubble_Integral_User(c14, c23));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c4, c23));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c2));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g_qmmpqp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, p, qp}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g :  qmmpqp LT");
#endif
 
//#define TimeStamp "Fri 12 Nov 2010 14:03:41 on n48"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> s13 = S(1,3);
complex<T> s14 = -(spa14*spb14);
complex<T> t2 = square(spa12); 
complex<T> t4 = spa24*spb14; 
complex<T> d1 = spa14*spa23*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = spa23*spa34; d2 = T(1)/d2;
complex<T> d3 = s13*spa23*spa34; d3 = T(1)/d3;
complex<T> d4 = s13*spa23; d4 = T(1)/d4;
complex<T> d5 = s13*spa23*T(2); d5 = T(1)/d5;
complex<T> t6 = s14*t2; 
complex<T> t1 = -(t4*(d2*t2 + d3*t6)); 
complex<T> co1 = -(d1*spa24*t2*T(3)); 
complex<T> co2 = d4*spb34*t2*t4; 
complex<T> co3 = d5*spb34*t4*t6; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g_qmpqpm_LT_wCI::\
C2q2g_qmpqpm_LT_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
         vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
         vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
	 vector<int> c41;  c41.push_back(ind[4-1]); c41.push_back(ind[1-1]);
         vector<int> c14 = c41;
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c34));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c4));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g_qmpqpm_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, m}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g :  qmpqpm LT");
#endif
 
//#define TimeStamp "Fri 12 Nov 2010 14:03:41 on n48"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = S(1,2);
complex<T> s23 = S(2,3);
complex<T> t1 = square(spb23); 
complex<T> d1 = spb14*spb34*T(2); d1 = T(1)/d1;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d1*s12*s23*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g_qmmqpp_LT_wCI::\
C2q2g_qmmqpp_LT_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
         vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
         vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
	 vector<int> c41;  c41.push_back(ind[4-1]); c41.push_back(ind[1-1]);
         vector<int> c14 = c41;
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c14));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g_qmmqpp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, qp, p}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g :  qmmqpp LT");
#endif
 
//#define TimeStamp "Fri 12 Nov 2010 14:03:42 on n48"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
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
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g_qmqppm_LT_wCI::\
C2q2g_qmqppm_LT_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
         vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
         vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
	 vector<int> c41;  c41.push_back(ind[4-1]); c41.push_back(ind[1-1]);
         vector<int> c14 = c41;
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g_qmqppm_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, m}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g :  qmqppm LT");
#endif
 
//#define TimeStamp "Fri 12 Nov 2010 14:03:43 on n48"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s13 = S(1,3);
complex<T> t2 = square(spa14); 
complex<T> t4 = spa24*spb12; 
complex<T> d1 = spa12*spa23*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = spa23*spa34; d2 = T(1)/d2;
complex<T> d3 = s13*spa23*spa34; d3 = T(1)/d3;
complex<T> d4 = s13*spa34; d4 = T(1)/d4;
complex<T> d5 = spa34; d5 = T(1)/d5;
complex<T> d6 = s13*spa34*T(2); d6 = T(1)/d6;
complex<T> t3 = -t4; 
complex<T> t5 = spb23*t2; 
complex<T> t1 = (d2 + d3*s12)*t2*t3; 
complex<T> t8 = t4*t5; 
complex<T> t7 = (d5 + d6*s12)*t8; 
complex<T> co1 = d1*spa24*t2*T(3); 
complex<T> co2 = d4*t8; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + t7*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2g_qmqpmp_LT_wCI::\
C2q2g_qmqpmp_LT_wCI
      (const std::vector<int>& ind){
	//define corner vectors
	 vector<int> c1;  c1.push_back(ind[1-1]);
	 vector<int> c2;  c2.push_back(ind[2-1]);
	 vector<int> c3;  c3.push_back(ind[3-1]);
	 vector<int> c4;  c4.push_back(ind[4-1]);

	 vector<int> c12;  c12.push_back(ind[1-1]); c12.push_back(ind[2-1]);
 	 vector<int> c13;  c13.push_back(ind[1-1]); c13.push_back(ind[3-1]);
         vector<int> c23;  c23.push_back(ind[2-1]); c23.push_back(ind[3-1]);
	 vector<int> c34;  c34.push_back(ind[3-1]); c34.push_back(ind[4-1]);
         vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
	 vector<int> c41;  c41.push_back(ind[4-1]); c41.push_back(ind[1-1]);
         vector<int> c14 = c41;
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c34));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g_qmqpmp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, p}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g :  qmqpmp LT");
#endif
 
//#define TimeStamp "Fri 12 Nov 2010 14:03:44 on n48"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spa12 = SPA(1,2);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = S(2,3);
complex<T> t3 = square(spa13); 
complex<T> t7 = spa14*T(2); 
complex<T> t9 = s23*spb12; 
complex<T> t10 = cube(spa13); 
complex<T> t11 = cube(s12); 
complex<T> t17 = s12*spa13; 
complex<T> d2 = spa14*spa34*square(spb13); d2 = T(1)/d2;
complex<T> d4 = spa14*spa34; d4 = T(1)/d4;
complex<T> d5 = spa14*spa34*cube(spb13); d5 = T(1)/d5;
complex<T> t4 = -t9; 
complex<T> t8 = d2*T(3); 
complex<T> t19 = d2*spa13; 
complex<T> t21 = -(d4*t10); 
complex<T> d1 = spa12*spa34*t7; d1 = T(1)/d1;
complex<T> d3 = spa34*spb13*t7; d3 = T(1)/d3;
complex<T> d6 = spa34*t7*cube(spb13); d6 = T(1)/d6;
complex<T> d7 = spa34*t7*square(spb13); d7 = T(1)/d7;
complex<T> t1 = d5*t9*square(s12) + spa13*spb12*t8*square(s23); 
complex<T> t2 = t19*t9 + d1*t10*T(3) + d3*spb12*t3*T(3); 
complex<T> t14 = -(d3*T(3)); 
complex<T> t15 = d4*t10*t4 + d6*t11*t9 + d7*spb12*t17*square(s23)*T(3); 
complex<T> t22 = d5*spb12*t11 + spb12*t21 + t17*t8*t9; 
complex<T> t20 = spb12*t14*t3 + t19*t4; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t20*(*CI_users[1]->get_value(mc,ind,mu)) + t22*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t15*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2g_qmpqpm_nfLT_wCI::\
C2q2g_qmpqpm_nfLT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g_qmpqpm_nfLT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g :  qmpqpm nfLT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g_qmmqpp_nfLT_wCI::\
C2q2g_qmmqpp_nfLT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g_qmmqpp_nfLT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g :  qmmqpp nfLT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g_qmqppm_nfLT_wCI::\
C2q2g_qmqppm_nfLT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g_qmqppm_nfLT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g :  qmqppm nfLT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g_qmqpmp_nfLT_wCI::\
C2q2g_qmqpmp_nfLT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g_qmqpmp_nfLT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g :  qmqpmp nfLT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
 // *************** table of switch values ************* 
 
#define _C_qmpmqp_LT C2q2g_141_LT
#define _C_qmmpqp_LT C2q2g_177_LT
#define _C_qmpqpm_LT C2q2g_45_LT
#define _C_qmmqpp_LT C2q2g_225_LT
#define _C_qmqppm_LT C2q2g_57_LT
#define _C_qmqpmp_LT C2q2g_201_LT
#define _C_qmpqpm_nfLT C2q2g_45_nfLT
#define _C_qmmqpp_nfLT C2q2g_225_nfLT
#define _C_qmqppm_nfLT C2q2g_57_nfLT
#define _C_qmqpmp_nfLT C2q2g_201_nfLT
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmpmqp_LT case 141
 
#define _CASE_qmmpqp_LT case 177
 
#define _CASE_qmpqpm_LT case 45
 
#define _CASE_qmmqpp_LT case 225
 
#define _CASE_qmqppm_LT case 57
 
#define _CASE_qmqpmp_LT case 201
 
#define _CASE_qmpqpm_nfLT case 45
 
#define _CASE_qmmqpp_nfLT case 225
 
#define _CASE_qmqppm_nfLT case 57
 
#define _CASE_qmqpmp_nfLT case 201
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2g_LT( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmpmqp_LT: return new 
                       C2q2g_qmpmqp_LT_wCI(ind);
    _CASE_qmmpqp_LT: return new 
                       C2q2g_qmmpqp_LT_wCI(ind);
    _CASE_qmpqpm_LT: return new 
                       C2q2g_qmpqpm_LT_wCI(ind);
    _CASE_qmmqpp_LT: return new 
                       C2q2g_qmmqpp_LT_wCI(ind);
    _CASE_qmqppm_LT: return new 
                       C2q2g_qmqppm_LT_wCI(ind);
    _CASE_qmqpmp_LT: return new 
                       C2q2g_qmqpmp_LT_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g_nfLT( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmpqpm_nfLT: return new 
                       C2q2g_qmpqpm_nfLT_wCI(ind);
    _CASE_qmmqpp_nfLT: return new 
                       C2q2g_qmmqpp_nfLT_wCI(ind);
    _CASE_qmqppm_nfLT: return new 
                       C2q2g_qmqppm_nfLT_wCI(ind);
    _CASE_qmqpmp_nfLT: return new 
                       C2q2g_qmqpmp_nfLT_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
