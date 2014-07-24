/*
*C_4g_wCI.cpp
*
* Created on 11/8, 2010
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
 



 
class C4g_mmpp_G_wCI : public Cut_Part_wCI {
public:
       C4g_mmpp_G_wCI
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

 
class C4g_ppmm_G_wCI : public Cut_Part_wCI {
public:
       C4g_ppmm_G_wCI
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

 
class C4g_mpmp_G_wCI : public Cut_Part_wCI {
public:
       C4g_mpmp_G_wCI
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

 
class C4g_pmpm_G_wCI : public Cut_Part_wCI {
public:
       C4g_pmpm_G_wCI
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

 
class C4g_mppm_G_wCI : public Cut_Part_wCI {
public:
       C4g_mppm_G_wCI
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

 
class C4g_pmmp_G_wCI : public Cut_Part_wCI {
public:
       C4g_pmmp_G_wCI
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

 
class C4g_mmpp_nf_wCI : public Cut_Part_wCI {
public:
       C4g_mmpp_nf_wCI
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

 
class C4g_ppmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g_ppmm_nf_wCI
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

 
class C4g_mpmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g_mpmp_nf_wCI
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

 
class C4g_pmpm_nf_wCI : public Cut_Part_wCI {
public:
       C4g_pmpm_nf_wCI
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

 
class C4g_mppm_nf_wCI : public Cut_Part_wCI {
public:
       C4g_mppm_nf_wCI
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

 
class C4g_pmmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g_pmmp_nf_wCI
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


C4g_mmpp_G_wCI::\
C4g_mmpp_G_wCI
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
     C4g_mmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g :  mmpp G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:04 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> s12 = S(1,2);
complex<T> t1 = cube(spa12); 
complex<T> d1 = spa14*spa23*spa34*T(3); d1 = T(1)/d1;
complex<T> d2 = spa14*spa34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(11); 
complex<T> co2 = -(d2*s12*spb23*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C4g_ppmm_G_wCI::\
C4g_ppmm_G_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c1, c2));
} 
  
  
template <class T> SeriesC<T> 
     C4g_ppmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g :  ppmm G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:04 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> s34 = S(3,4);
complex<T> t1 = cube(spa34); 
complex<T> d1 = spa12*spa14*spa23*T(3); d1 = T(1)/d1;
complex<T> d2 = spa12*spa23; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(11); 
complex<T> co2 = -(d2*s34*spb14*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C4g_mpmp_G_wCI::\
C4g_mpmp_G_wCI
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
     C4g_mpmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g :  mpmp G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:06 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> t3 = square(spa13); 
complex<T> t4 = spb12*spb23; 
complex<T> t6 = -(s12*T(2)); 
complex<T> t7 = square(s23); 
complex<T> t9 = cube(spa13); 
complex<T> t12 = spa34*T(6); 
complex<T> t17 = square(square(spa13)); 
complex<T> d3 = spa14*spa34*cube(spb13); d3 = T(1)/d3;
complex<T> d5 = spa14*spa34*square(square(spb13)); d5 = T(1)/d5;
complex<T> d6 = spa14*spa34*square(spb13); d6 = T(1)/d6;
complex<T> d7 = spa14*spa34; d7 = T(1)/d7;
complex<T> t10 = -(s12*t4); 
complex<T> t20 = t4*T(2); 
complex<T> t21 = d5*t7; 
complex<T> t23 = d6*t4; 
complex<T> d1 = spa12*spa14*spa23*t12; d1 = T(1)/d1;
complex<T> d2 = spa14*spa23*spb13*t12; d2 = T(1)/d2;
complex<T> d4 = spa12*spa14*spb13*t12; d4 = T(1)/d4;
complex<T> t8 = d1*T(11); 
complex<T> t14 = -(d2*T(11)); 
complex<T> t28 = t23*t3; 
complex<T> t1 = d3*(s12 - s23)*spa13*t4 + t17*t8 + (d2*spb12 - d4*spb23)*t9*T(11); 
complex<T> t2 = s12*(d5*s12*s23*t20 - t28*T(4)); 
complex<T> t15 = d7*t17*t4 + s23*t28*t6 + t21*t4*square(s12); 
complex<T> t16 = d3*spa13*(t10 + s23*t4) + t17*t8 + spb12*t14*t9 + d4*spb23*t9*T(11); 
complex<T> t27 = s12*t20*t21 - s23*t28*T(4); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t16*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t27*(*CI_users[3]->get_value(mc,ind,mu)) + t15*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C4g_pmpm_G_wCI::\
C4g_pmpm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C4g_pmpm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g :  pmpm G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:08 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> t3 = square(spa24); 
complex<T> t4 = spb23*spb34; 
complex<T> t6 = -(s23*T(2)); 
complex<T> t7 = square(s34); 
complex<T> t9 = cube(spa24); 
complex<T> t12 = spa14*T(6); 
complex<T> t17 = square(square(spa24)); 
complex<T> d3 = spa12*spa14*cube(spb24); d3 = T(1)/d3;
complex<T> d5 = spa12*spa14*square(square(spb24)); d5 = T(1)/d5;
complex<T> d6 = spa12*spa14*square(spb24); d6 = T(1)/d6;
complex<T> d7 = spa12*spa14; d7 = T(1)/d7;
complex<T> t10 = -(s23*t4); 
complex<T> t20 = t4*T(2); 
complex<T> t21 = d5*t7; 
complex<T> t23 = d6*t4; 
complex<T> d1 = spa12*spa23*spa34*t12; d1 = T(1)/d1;
complex<T> d2 = spa12*spa34*spb24*t12; d2 = T(1)/d2;
complex<T> d4 = spa12*spa23*spb24*t12; d4 = T(1)/d4;
complex<T> t8 = d1*T(11); 
complex<T> t14 = -(d2*T(11)); 
complex<T> t28 = t23*t3; 
complex<T> t1 = d3*(s23 - s34)*spa24*t4 + t17*t8 + (d2*spb23 - d4*spb34)*t9*T(11); 
complex<T> t2 = s23*(d5*s23*s34*t20 - t28*T(4)); 
complex<T> t15 = d7*t17*t4 + s34*t28*t6 + t21*t4*square(s23); 
complex<T> t16 = d3*spa24*(t10 + s34*t4) + t17*t8 + spb23*t14*t9 + d4*spb34*t9*T(11); 
complex<T> t27 = s23*t20*t21 - s34*t28*T(4); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t16*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t27*(*CI_users[3]->get_value(mc,ind,mu)) + t15*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C4g_mppm_G_wCI::\
C4g_mppm_G_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
} 
  
  
template <class T> SeriesC<T> 
     C4g_mppm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g :  mppm G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:09 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> t1 = cube(spa14); 
complex<T> d1 = spa12*spa23*spa34*T(3); d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(11); 
complex<T> co2 = d2*spb12*spb23*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C4g_pmmp_G_wCI::\
C4g_pmmp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c12));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C4g_pmmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g :  pmmp G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:09 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> s23 = S(2,3);
complex<T> t1 = cube(spa23); 
complex<T> d1 = spa12*spa14*spa34*T(3); d1 = T(1)/d1;
complex<T> d2 = spa12*spa14; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(11); 
complex<T> co2 = -(d2*s23*spb34*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C4g_mmpp_nf_wCI::\
C4g_mmpp_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C4g_mmpp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g :  mmpp nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:10 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = cube(spa12); 
complex<T> d1 = spa14*spa23*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g_ppmm_nf_wCI::\
C4g_ppmm_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C4g_ppmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g :  ppmm nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:10 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = cube(spa34); 
complex<T> d1 = spa12*spa14*spa23*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g_mpmp_nf_wCI::\
C4g_mpmp_nf_wCI
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
     C4g_mpmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g :  mpmp nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:12 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> t3 = square(spa13); 
complex<T> t6 = spa14*spa34; 
complex<T> t7 = square(s12); 
complex<T> t8 = square(s23); 
complex<T> t10 = cube(spa13); 
complex<T> t11 = spb12*spb23; 
complex<T> t16 = square(square(spa13)); 
complex<T> t4 = -t11; 
complex<T> t13 = t6*T(3); 
complex<T> t19 = s23*t11; 
complex<T> d3 = t6*cube(spb13); d3 = T(1)/d3;
complex<T> d5 = t6*square(square(spb13)); d5 = T(1)/d5;
complex<T> d6 = t6*square(spb13); d6 = T(1)/d6;
complex<T> d7 = t6*square(spb13)*T(2); d7 = T(1)/d7;
complex<T> t9 = -(d5*T(2)); 
complex<T> t14 = d3*s12; 
complex<T> t18 = d7*s12*t19*t3 + d5*t4*t7*t8; 
complex<T> d1 = spa12*spa23*t13; d1 = T(1)/d1;
complex<T> d2 = spa23*spb13*t13; d2 = T(1)/d2;
complex<T> d4 = spa12*spb13*t13; d4 = T(1)/d4;
complex<T> t1 = d2*spb12*t10 - d4*spb23*t10 + spa13*t11*t14 - d1*t16 + d3*s23*spa13*t4; 
complex<T> t2 = d6*s12*t11*t3 + t19*t7*t9; 
complex<T> t15 = -(d2*spb12*t10) + d4*spb23*t10 - d1*t16 + d3*spa13*t19 + spa13*t14*t4; 
complex<T> t22 = d6*t19*t3 + s12*t11*t8*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t15*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t22*(*CI_users[3]->get_value(mc,ind,mu)) + t18*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C4g_pmpm_nf_wCI::\
C4g_pmpm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C4g_pmpm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g :  pmpm nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:13 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> t3 = square(spa24); 
complex<T> t6 = spa12*spa14; 
complex<T> t7 = square(s23); 
complex<T> t8 = square(s34); 
complex<T> t10 = cube(spa24); 
complex<T> t11 = spb23*spb34; 
complex<T> t16 = square(square(spa24)); 
complex<T> t4 = -t11; 
complex<T> t13 = t6*T(3); 
complex<T> t19 = s34*t11; 
complex<T> d3 = t6*cube(spb24); d3 = T(1)/d3;
complex<T> d5 = t6*square(square(spb24)); d5 = T(1)/d5;
complex<T> d6 = t6*square(spb24); d6 = T(1)/d6;
complex<T> d7 = t6*square(spb24)*T(2); d7 = T(1)/d7;
complex<T> t9 = -(d5*T(2)); 
complex<T> t14 = d3*s23; 
complex<T> t18 = d7*s23*t19*t3 + d5*t4*t7*t8; 
complex<T> d1 = spa23*spa34*t13; d1 = T(1)/d1;
complex<T> d2 = spa34*spb24*t13; d2 = T(1)/d2;
complex<T> d4 = spa23*spb24*t13; d4 = T(1)/d4;
complex<T> t1 = d2*spb23*t10 - d4*spb34*t10 + spa24*t11*t14 - d1*t16 + d3*s34*spa24*t4; 
complex<T> t2 = d6*s23*t11*t3 + t19*t7*t9; 
complex<T> t15 = -(d2*spb23*t10) + d4*spb34*t10 - d1*t16 + d3*spa24*t19 + spa24*t14*t4; 
complex<T> t22 = d6*t19*t3 + s23*t11*t8*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t15*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t22*(*CI_users[3]->get_value(mc,ind,mu)) + t18*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C4g_mppm_nf_wCI::\
C4g_mppm_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C4g_mppm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g :  mppm nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:13 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = cube(spa14); 
complex<T> d1 = spa12*spa23*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g_pmmp_nf_wCI::\
C4g_pmmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g_pmmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g :  pmmp nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 16:00:14 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = cube(spa23); 
complex<T> d1 = spa12*spa14*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_mmpp_G C4g_12_G
#define _C_ppmm_G C4g_3_G
#define _C_mpmp_G C4g_10_G
#define _C_pmpm_G C4g_5_G
#define _C_mppm_G C4g_6_G
#define _C_pmmp_G C4g_9_G
#define _C_mmpp_nf C4g_12_nf
#define _C_ppmm_nf C4g_3_nf
#define _C_mpmp_nf C4g_10_nf
#define _C_pmpm_nf C4g_5_nf
#define _C_mppm_nf C4g_6_nf
#define _C_pmmp_nf C4g_9_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_mmpp_G case 12
 
#define _CASE_ppmm_G case 3
 
#define _CASE_mpmp_G case 10
 
#define _CASE_pmpm_G case 5
 
#define _CASE_mppm_G case 6
 
#define _CASE_pmmp_G case 9
 
#define _CASE_mmpp_nf case 12
 
#define _CASE_ppmm_nf case 3
 
#define _CASE_mpmp_nf case 10
 
#define _CASE_pmpm_nf case 5
 
#define _CASE_mppm_nf case 6
 
#define _CASE_pmmp_nf case 9
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_4g_G( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_mmpp_G: return new 
                       C4g_mmpp_G_wCI(ind);
    _CASE_ppmm_G: return new 
                       C4g_ppmm_G_wCI(ind);
    _CASE_mpmp_G: return new 
                       C4g_mpmp_G_wCI(ind);
    _CASE_pmpm_G: return new 
                       C4g_pmpm_G_wCI(ind);
    _CASE_mppm_G: return new 
                       C4g_mppm_G_wCI(ind);
    _CASE_pmmp_G: return new 
                       C4g_pmmp_G_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_4g_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_mmpp_nf: return new 
                       C4g_mmpp_nf_wCI(ind);
    _CASE_ppmm_nf: return new 
                       C4g_ppmm_nf_wCI(ind);
    _CASE_mpmp_nf: return new 
                       C4g_mpmp_nf_wCI(ind);
    _CASE_pmpm_nf: return new 
                       C4g_pmpm_nf_wCI(ind);
    _CASE_mppm_nf: return new 
                       C4g_mppm_nf_wCI(ind);
    _CASE_pmmp_nf: return new 
                       C4g_pmmp_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
