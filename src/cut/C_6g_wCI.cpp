/*
*C_6g_wCI.cpp
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
 



 
class C6g_mmpppp_G_wCI : public Cut_Part_wCI {
public:
       C6g_mmpppp_G_wCI
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

 
class C6g_mmpppp_nf_wCI : public Cut_Part_wCI {
public:
       C6g_mmpppp_nf_wCI
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

 
class C6g_mpmppp_G_wCI : public Cut_Part_wCI {
public:
       C6g_mpmppp_G_wCI
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

 
class C6g_mpmppp_nf_wCI : public Cut_Part_wCI {
public:
       C6g_mpmppp_nf_wCI
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

 
class C6g_mppmpp_G_wCI : public Cut_Part_wCI {
public:
       C6g_mppmpp_G_wCI
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

 
class C6g_mppmpp_nf_wCI : public Cut_Part_wCI {
public:
       C6g_mppmpp_nf_wCI
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


// cc
 
class C6g_ppmmmm_G_wCI : public Cut_Part_wCI {
public:
       C6g_ppmmmm_G_wCI
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

 
class C6g_ppmmmm_nf_wCI : public Cut_Part_wCI {
public:
       C6g_ppmmmm_nf_wCI
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

 
class C6g_pmpmmm_G_wCI : public Cut_Part_wCI {
public:
       C6g_pmpmmm_G_wCI
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

 
class C6g_pmpmmm_nf_wCI : public Cut_Part_wCI {
public:
       C6g_pmpmmm_nf_wCI
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

 
class C6g_pmmpmm_G_wCI : public Cut_Part_wCI {
public:
       C6g_pmmpmm_G_wCI
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

 
class C6g_pmmpmm_nf_wCI : public Cut_Part_wCI {
public:
       C6g_pmmpmm_nf_wCI
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

//cc-end


C6g_mmpppp_G_wCI::\
C6g_mmpppp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c2, c34, c5, c16));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c126));
CI_users.push_back(new Cached_Box_Integral_User(c3, c45, c6, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c6, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c6, c1, c234));
CI_users.push_back(new Cached_Box_Integral_User(c6, c1, c2, c345));
} 
  
  
template <class T> SeriesC<T> 
     C6g_mmpppp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C6g :  mmpppp G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 17:38:58 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spb16 = SPB(1,6);
complex<T> spa16 = SPA(1,6);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb56 = SPB(5,6);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa26 = SPA(2,6);
complex<T> spa14 = SPA(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb46 = SPB(4,6);
complex<T> spb14 = SPB(1,4);
complex<T> s12 = S(1,2);
complex<T> s16 = -(spa16*spb16);
complex<T> s156 = SS(1,5,6);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s123 = SS(1,2,3);
complex<T> s345 = SS(3,4,5);
complex<T> t5 = spa16*spa34; 
complex<T> t6 = spa23*spa56; 
complex<T> t9 = square(spa12*spb15 + spa26*spb56); 
complex<T> t12 = square(spa15); 
complex<T> t14 = square(spa12); 
complex<T> t16 = -(spa24*(spa15*spb45 + spa16*spb46)); 
complex<T> t22 = cube(spa12); 
complex<T> t23 = -(spa12*spb15) - spa26*spb56; 
complex<T> t24 = square(spa25); 
complex<T> t25 = square(spb56); 
complex<T> t28 = spa15*spb45 + spa16*spb46; 
complex<T> t29 = spa12*spb14 + spa25*spb45 + spa26*spb46; 
complex<T> t32 = square(spa14); 
complex<T> t33 = square(spa24); 
complex<T> t36 = spa25*spb56; 
complex<T> t41 = -(spa26*spb46); 
complex<T> d10 = spa16*spa45*spa56*T(2); d10 = T(1)/d10;
complex<T> d15 = spa23*spa34*spa45*T(2); d15 = T(1)/d15;
complex<T> t15 = spa45*t6; 
complex<T> t17 = spa14*(-(spa12*spb14) - spa25*spb45 + t41); 
complex<T> t38 = spa15*t23; 
complex<T> t39 = square(t28); 
complex<T> t40 = square(t29); 
complex<T> t55 = s12*t22; 
complex<T> t71 = spb34*t22; 
complex<T> t77 = spb56*t22; 
complex<T> d7 = spa45*spa56*t5*T(2); d7 = T(1)/d7;
complex<T> d9 = spa45*t5*T(2); d9 = T(1)/d9;
complex<T> d12 = spa16*t6*T(2); d12 = T(1)/d12;
complex<T> d13 = t5*t6*T(2); d13 = T(1)/d13;
complex<T> d14 = spa23*t5*T(2); d14 = T(1)/d14;
complex<T> d1 = (-s156 + s16)*spa34*t15*T(6); d1 = T(1)/d1;
complex<T> d2 = (-s156 + s16)*t15*t5*T(6); d2 = T(1)/d2;
complex<T> d3 = spa34*t15*cube(-s156 + s16)*T(3); d3 = T(1)/d3;
complex<T> d4 = t15*t5*T(6); d4 = T(1)/d4;
complex<T> d5 = (s23 - s234)*t15*t5*T(6); d5 = T(1)/d5;
complex<T> d6 = t15*t5*cube(s23 - s234)*T(3); d6 = T(1)/d6;
complex<T> d8 = t15*t5*T(2); d8 = T(1)/d8;
complex<T> d11 = t15*T(2); d11 = T(1)/d11;
complex<T> d16 = spa34*t15*T(2); d16 = T(1)/d16;
complex<T> t30 = d4*T(11); 
complex<T> t31 = d8*s123; 
complex<T> t42 = d6*t32; 
complex<T> t45 = d8*s234*s345*t22 - d11*spb16*t71; 
complex<T> t52 = d3*t24; 
complex<T> t61 = d6*t33; 
complex<T> t62 = d3*t9; 
complex<T> t1 = -(t16*t40*t42) - t17*t39*t61 - d5*t14*(t16 + t17)*T(11); 
complex<T> t2 = t22*t30 + t16*t40*t42 + t17*t39*t61 + d5*t14*(t16 + t17)*T(11); 
complex<T> t54 = s234*t22*t31 - d9*spb23*t77; 
complex<T> t63 = s345*t22*t31 + d13*spb45*t55; 
complex<T> t76 = spa16*t52; 
complex<T> t3 = t22*t30 - t12*t36*t62 + t25*t38*t76 - d1*t14*t36*T(11) + d2*t14*t38*T(11); 
complex<T> t4 = t12*t36*t62 - t25*t38*t76 + d1*t14*t36*T(11) - d2*t14*t38*T(11); 
complex<T> co1 = -(d7*spb23*t55); 
complex<T> co2 = d10*spb23*t71; 
complex<T> co3 = d12*spb45*t71; 
complex<T> co4 = d14*spb45*t77; 
complex<T> co5 = d15*spb16*t77; 
complex<T> co6 = -(d16*spb16*t55); 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + t54*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t45*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t63*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)) + co5*(*CI_users[11]->get_value(mc,ind,mu)) + co6*(*CI_users[12]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_mmpppp_nf_wCI::\
C6g_mmpppp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
} 
  
  
template <class T> SeriesC<T> 
     C6g_mmpppp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C6g :  mmpppp nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 17:39:02 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa16 = SPA(1,6);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb56 = SPB(5,6);
complex<T> spb15 = SPB(1,5);
complex<T> spa26 = SPA(2,6);
complex<T> spa14 = SPA(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb46 = SPB(4,6);
complex<T> spb14 = SPB(1,4);
complex<T> s16 = S(1,6);
complex<T> s156 = SS(1,5,6);
complex<T> s23 = S(2,3);
complex<T> s234 = SS(2,3,4);
complex<T> t7 = square(spa12*spb15 + spa26*spb56); 
complex<T> t10 = square(spa15); 
complex<T> t12 = square(spa12); 
complex<T> t13 = square(spa15*spb45 + spa16*spb46); 
complex<T> t14 = square(spa12*spb14 + spa25*spb45 + spa26*spb46); 
complex<T> t16 = spa23*T(3); 
complex<T> t17 = spa34*spa45; 
complex<T> t21 = square(spa14); 
complex<T> t22 = square(spa24); 
complex<T> t23 = -(spa15*spb45) - spa16*spb46; 
complex<T> t24 = -(spa12*spb14) - spa25*spb45 - spa26*spb46; 
complex<T> t25 = -(spa12*spb15) - spa26*spb56; 
complex<T> t26 = square(spa25); 
complex<T> t27 = square(spb56); 
complex<T> t30 = cube(spa12); 
complex<T> t47 = t14*t21; 
complex<T> t48 = t26*t27; 
complex<T> d1 = (-s156 + s16)*spa56*t16*t17; d1 = T(1)/d1;
complex<T> d2 = (-s156 + s16)*spa16*spa56*t16*t17; d2 = T(1)/d2;
complex<T> d3 = spa56*t16*t17*cube(-s156 + s16); d3 = T(1)/d3;
complex<T> d4 = spa16*spa56*t16*t17; d4 = T(1)/d4;
complex<T> d5 = (s23 - s234)*spa16*spa56*t16*t17; d5 = T(1)/d5;
complex<T> d6 = spa16*spa56*t16*t17*cube(s23 - s234); d6 = T(1)/d6;
complex<T> t42 = d3*spa16; 
complex<T> t46 = d6*t13; 
complex<T> t54 = d3*t7; 
complex<T> t1 = -(d5*spa24*t12*t23) - d5*spa14*t12*t24 - d4*t30 - spa14*t22*t24*t46 - d6*spa24*t23*t47; 
complex<T> t2 = d5*t12*(spa24*t23 + spa14*t24) + spa14*t22*t24*t46 + d6*spa24*t23*t47; 
complex<T> t3 = -(d1*spa25*spb56*t12) + d2*spa15*t12*t25 + spa15*t25*t42*t48 - spa25*spb56*t10*t54; 
complex<T> t4 = d1*spa25*spb56*t12 - d2*spa15*t12*t25 - d4*t30 - spa15*t25*t42*t48 + spa25*spb56*t10*t54; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_mpmppp_G_wCI::\
C6g_mpmppp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c126, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c6, c2345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c1236));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c6, c345));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c2, c34, c5, c16));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c126));
CI_users.push_back(new Cached_Box_Integral_User(c3, c45, c6, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c6, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c34, c2, c16));
CI_users.push_back(new Cached_Box_Integral_User(c5, c6, c1, c234));
CI_users.push_back(new Cached_Box_Integral_User(c6, c1, c2, c345));
} 
  
  
template <class T> SeriesC<T> 
     C6g_mpmppp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C6g :  mpmppp G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 17:47:44 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa16 = SPA(1,6);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa26 = SPA(2,6);
complex<T> spa36 = SPA(3,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb56 = SPB(5,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb46 = SPB(4,6);
complex<T> spb14 = SPB(1,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> s16 = -(spa16*spb16);
complex<T> s126 = SS(1,2,6);
complex<T> s156 = SS(1,5,6);
complex<T> s24 = -(spa24*spb24);
complex<T> s234 = SS(2,3,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = S(3,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s123 = SS(1,2,3);
complex<T> s36 = S(3,6);
complex<T> s46 = S(4,6);
complex<T> s56 = -(spa56*spb56);
complex<T> s345 = SS(3,4,5);
complex<T> t7 = square(spa35); 
complex<T> t8 = square(spa15); 
complex<T> t10 = spa45*spa56; 
complex<T> t11 = square(spa14); 
complex<T> t13 = spa12*spa16; 
complex<T> t18 = square(spa13*spb23 + spa14*spb24); 
complex<T> t19 = square(spa34*spb24 + spa35*spb25); 
complex<T> t20 = square(spa13*spb16 + spa23*spb26); 
complex<T> t21 = square(spa13*spb14 + spa35*spb45 + spa36*spb46); 
complex<T> t36 = square(spa13); 
complex<T> t37 = square(spa36); 
complex<T> t38 = -(spa14*spb24); 
complex<T> t39 = -(spa23*T(2)); 
complex<T> t41 = cube(spa25); 
complex<T> t46 = (-s234 + s34)*spa16; 
complex<T> t47 = square(s23 + s24); 
complex<T> t49 = square(spa13*spb23 + spa14*spb24 + spa15*spb25); 
complex<T> t51 = square(spa15*spb45 + spa16*spb46); 
complex<T> t57 = cube(spa13); 
complex<T> t58 = -(spa34*spb24) - spa35*spb25; 
complex<T> t60 = -(spa13*spb16) - spa23*spb26; 
complex<T> t61 = -(spa13*spb15) - spa36*spb56; 
complex<T> t62 = spa12*spb25 - spa16*spb56; 
complex<T> t64 = -(spa13*spb14) - spa35*spb45 - spa36*spb46; 
complex<T> t65 = -(spa15*spb45) - spa16*spb46; 
complex<T> t66 = -(spa13*spb15) - spa23*spb25 - spa36*spb56; 
complex<T> t68 = s16 - s234; 
complex<T> t70 = -(spa12*spb24); 
complex<T> t71 = s12 - s345; 
complex<T> t73 = square(spa23); 
complex<T> t76 = spa36*spb26; 
complex<T> t78 = -s156 + s16; 
complex<T> t80 = cube(spa15); 
complex<T> t83 = -(spa13*spb14) - spa23*spb24 - spa35*spb45 - spa36*spb46; 
complex<T> t86 = s234*s345; 
complex<T> t97 = square(spa34); 
complex<T> t98 = cube(spa35); 
complex<T> t99 = square(spb56); 
complex<T> t102 = spa15*spa35; 
complex<T> t108 = square(spa16); 
complex<T> t109 = square(spb26); 
complex<T> t111 = spa26*spa45; 
complex<T> t112 = spa24*spa56; 
complex<T> t118 = spb16*spb34; 
complex<T> t127 = spa12*spa23; 
complex<T> t134 = square(square(spa13)); 
complex<T> t137 = square(spa12); 
complex<T> t139 = square(spb24); 
complex<T> t148 = spa13*spb15 + spa36*spb56; 
complex<T> t155 = spa23*T(2); 
complex<T> t160 = spa34*spa56; 
complex<T> t190 = s34*spa24; 
complex<T> t235 = spa14*spa34; 
complex<T> t287 = spa15*spa16; 
complex<T> t16 = -t76; 
complex<T> t17 = -(spa34*t65); 
complex<T> t23 = square(t148); 
complex<T> t24 = square(spa23*spb25 + t148); 
complex<T> t40 = spa34*t10; 
complex<T> t42 = -t102; 
complex<T> t59 = -(spa13*spb23) + t38; 
complex<T> t63 = -(spa13*spb23) - spa15*spb25 + t38; 
complex<T> t67 = spa16*t10; 
complex<T> t87 = -t118; 
complex<T> t100 = -(t36*T(2)); 
complex<T> t103 = spa34*t11; 
complex<T> t128 = t36*T(2); 
complex<T> t131 = spa36*t58; 
complex<T> t153 = square(t62); 
complex<T> t154 = t7*t8; 
complex<T> t161 = t62*t66; 
complex<T> t182 = -(t36*T(4)); 
complex<T> t183 = t13*t37; 
complex<T> t200 = spa25*t68; 
complex<T> t204 = t127*T(2); 
complex<T> t209 = spa14*t64; 
complex<T> t223 = spa15*t7; 
complex<T> t227 = spa25*t112; 
complex<T> t245 = spa16*t60; 
complex<T> t246 = spa35*t8; 
complex<T> t280 = spb12*t134; 
complex<T> t292 = spb56*t134; 
complex<T> t368 = t287*t61; 
complex<T> d12 = (-s126 + s16)*spa16*spa25*spa34*t111; d12 = T(1)/d12;
complex<T> d18 = (s16 - s345)*spa16*spa25*spa34*t111*T(6); d18 = T(1)/d18;
complex<T> d23 = spa16*spa25*spa34*t111*cube(s16 - s345)*T(3); d23 = T(1)/d23;
complex<T> d38 = t10*t46*cube(spa24); d38 = T(1)/d38;
complex<T> d39 = t10*t41*t46; d39 = T(1)/d39;
complex<T> d45 = t10*t13*t190*T(6); d45 = T(1)/d45;
complex<T> d59 = spa16*t160*square(spa25); d59 = T(1)/d59;
complex<T> d60 = spa16*t160*square(square(spa25)); d60 = T(1)/d60;
complex<T> d61 = t160*square(spa26); d61 = T(1)/d61;
complex<T> d62 = t160*square(square(spa26)); d62 = T(1)/d62;
complex<T> d65 = spa34*spa45*t13*T(2); d65 = T(1)/d65;
complex<T> d66 = t10*t13*T(2); d66 = T(1)/d66;
complex<T> d68 = spa56*t13*t155; d68 = T(1)/d68;
complex<T> d69 = spa16*t155*t160; d69 = T(1)/d69;
complex<T> d70 = spa34*t13*t155; d70 = T(1)/d70;
complex<T> d71 = t10*square(spa25); d71 = T(1)/d71;
complex<T> d72 = t10*square(square(spa25)); d72 = T(1)/d72;
complex<T> t14 = spa23*(spa15*spb25 - t59); 
complex<T> t106 = spb24*t59; 
complex<T> t119 = d59*spb45; 
complex<T> t123 = d23*t49; 
complex<T> t124 = d12*t58; 
complex<T> t129 = -t154; 
complex<T> t138 = spa26*t40; 
complex<T> t165 = d60*spb45; 
complex<T> t170 = d23*t19; 
complex<T> t186 = spa23*t40; 
complex<T> t206 = spa16*t40; 
complex<T> t210 = d61*spb45; 
complex<T> t215 = d38*t59; 
complex<T> t226 = d62*spb45; 
complex<T> t234 = t154*t204; 
complex<T> t290 = spa35*t24; 
complex<T> t299 = t36*t42; 
complex<T> d1 = t40*square(spa26); d1 = T(1)/d1;
complex<T> d4 = t40*t71*cube(spa26); d4 = T(1)/d4;
complex<T> d6 = spb12*t40*cube(spa26); d6 = T(1)/d6;
complex<T> d9 = t40*square(s36 + s46 + s56)*square(spa26); d9 = T(1)/d9;
complex<T> d17 = t40*cube(spa26); d17 = T(1)/d17;
complex<T> d20 = (s16 - s345)*t40*cube(spa26); d20 = T(1)/d20;
complex<T> d22 = t40*square(s23 + s24 + s25)*square(spa26); d22 = T(1)/d22;
complex<T> d24 = t40*t41*t68; d24 = T(1)/d24;
complex<T> d28 = t40*square(s25 + s35 + s45)*square(spa25); d28 = T(1)/d28;
complex<T> d30 = t67*square(spa24); d30 = T(1)/d30;
complex<T> d31 = spa23*spa24*t67*T(6); d31 = T(1)/d31;
complex<T> d32 = t67*cube(spa24); d32 = T(1)/d32;
complex<T> d33 = (s23 - s234)*t67*cube(spa24); d33 = T(1)/d33;
complex<T> d34 = (s23 - s234)*spa23*spa24*t67*T(6); d34 = T(1)/d34;
complex<T> d35 = t67*square(s24 + s34)*square(spa24); d35 = T(1)/d35;
complex<T> d36 = spa23*spa24*t67*cube(s23 - s234)*T(3); d36 = T(1)/d36;
complex<T> d37 = t227*t46*T(6); d37 = T(1)/d37;
complex<T> d40 = spa34*t227*t46*T(6); d40 = T(1)/d40;
complex<T> d41 = t40*t41*t46; d41 = T(1)/d41;
complex<T> d42 = t47*t67*square(spa24); d42 = T(1)/d42;
complex<T> d43 = t47*t67*square(spa25); d43 = T(1)/d43;
complex<T> d44 = spa16*t227*cube(-s234 + s34)*T(3); d44 = T(1)/d44;
complex<T> d46 = spb34*t67*cube(spa24); d46 = T(1)/d46;
complex<T> d47 = t13*t190*t40*T(6); d47 = T(1)/d47;
complex<T> d48 = (s34 - s345)*spa25*t13*t40*T(6); d48 = T(1)/d48;
complex<T> d51 = spa25*t13*t40*cube(s34 - s345)*T(3); d51 = T(1)/d51;
complex<T> d52 = t40*square(square(spa26)); d52 = T(1)/d52;
complex<T> d53 = t40*square(spa25); d53 = T(1)/d53;
complex<T> d54 = t40*square(square(spa25)); d54 = T(1)/d54;
complex<T> d55 = t67*square(square(spa24)); d55 = T(1)/d55;
complex<T> d57 = t67*square(spa25); d57 = T(1)/d57;
complex<T> d58 = t67*square(square(spa25)); d58 = T(1)/d58;
complex<T> d64 = t13*t155*t40; d64 = T(1)/d64;
complex<T> d67 = t10*t204; d67 = T(1)/d67;
complex<T> d73 = spa34*spa45*t204; d73 = T(1)/d73;
complex<T> d74 = t155*t40; d74 = T(1)/d74;
complex<T> t56 = d48*(-(spa23*spb25) + t61); 
complex<T> t88 = d64*s123; 
complex<T> t107 = -(d34*T(11)); 
complex<T> t113 = d52*s12; 
complex<T> t114 = d30*s23; 
complex<T> t117 = d55*s34; 
complex<T> t126 = d48*t62; 
complex<T> t162 = d1*s12; 
complex<T> t163 = d37*spb24; 
complex<T> t169 = d44*t18; 
complex<T> t181 = t134*(d64*t86 + d67*t87); 
complex<T> t184 = d30*spa14; 
complex<T> t187 = d34*T(11); 
complex<T> t188 = d4*spb26; 
complex<T> t192 = d43*t106; 
complex<T> t193 = d51*t153; 
complex<T> t207 = d1*spa36; 
complex<T> t212 = d44*t139; 
complex<T> t228 = d41*t59; 
complex<T> t233 = spa36*t182*t210 + t155*t183*t226 + spa12*t154*t165*t39 + t102*t119*t36*T(4); 
complex<T> d2 = (s12 - s126)*t138; d2 = T(1)/d2;
complex<T> d3 = t138*t71*T(6); d3 = T(1)/d3;
complex<T> d5 = spb12*t137*t138*T(6); d5 = T(1)/d5;
complex<T> d7 = (s12 - s126)*spa12*t138; d7 = T(1)/d7;
complex<T> d8 = spa12*t138*t71*T(6); d8 = T(1)/d8;
complex<T> d10 = t138*cube(t71)*T(3); d10 = T(1)/d10;
complex<T> d11 = spa12*spb12*t138*T(6); d11 = T(1)/d11;
complex<T> d13 = spa25*t186*t78; d13 = T(1)/d13;
complex<T> d14 = spa16*spa25*t186*t78; d14 = T(1)/d14;
complex<T> d15 = t206*square(spa25); d15 = T(1)/d15;
complex<T> d16 = spa16*spa23*t138*T(6); d16 = T(1)/d16;
complex<T> d19 = (s16 - s345)*t206*t41; d19 = T(1)/d19;
complex<T> d21 = t206*square(s23 + s24 + s25)*square(spa25); d21 = T(1)/d21;
complex<T> d25 = t186*t200*T(6); d25 = T(1)/d25;
complex<T> d26 = spa16*t186*t200*T(6); d26 = T(1)/d26;
complex<T> d27 = t206*t41*t68; d27 = T(1)/d27;
complex<T> d29 = spa25*t186*cube(t68)*T(3); d29 = T(1)/d29;
complex<T> d49 = (s34 - s345)*t206*t41; d49 = T(1)/d49;
complex<T> d50 = t206*square(s35 + s45)*square(spa25); d50 = T(1)/d50;
complex<T> d56 = t206*square(square(spa25)); d56 = T(1)/d56;
complex<T> d63 = t206*T(2); d63 = T(1)/d63;
complex<T> t2 = -(d32*spa12*spa13*t235) + t184*t36 + t187*t209*t36 - d33*spa12*t103*t64 + d33*spa23*t103*t65 + d36*t103*t21*t65 + spa34*t187*t36*t65 + d35*t103*t64*t65 + d36*t209*t51*t97 + d31*t57*T(11); 
complex<T> t35 = spa12*t124*t128 + d7*t100*t245 + d12*t36*t39*(-(spa15*spb25) + t59) + d2*t100*t76; 
complex<T> t96 = s34*spa14*t100*t114 + s23*t103*t117*t127; 
complex<T> t120 = d13*spb56; 
complex<T> t149 = spa14*t114*t182 + d55*s23*t103*t204; 
complex<T> t150 = s16*(spa36*t100*t162 + spa23*t113*t183); 
complex<T> t171 = d10*t20; 
complex<T> t180 = d53*spb16*t102*t182 + d52*s16*t155*t183 + s16*t182*t207 + d54*spb16*t234; 
complex<T> t195 = d21*t58; 
complex<T> t196 = d10*t60; 
complex<T> t197 = d14*t61; 
complex<T> t202 = spa36*t162*t182 + t113*t155*t183; 
complex<T> t203 = -(d65*spb23*t292) + s234*t134*t88; 
complex<T> t208 = d15*t102; 
complex<T> t219 = -(d69*spb45*t280) + s345*t134*t88; 
complex<T> t334 = spa34*t212; 
complex<T> t1 = -(spa23*t183*t188) + t109*t183*t196 + d7*t128*t245 + d3*t16*t36 + t207*t36 - d8*t245*t36 - d6*spa16*spa23*spa36*(t16 + t58) + d6*spa23*t37*(spa15*spb25 + spa16*spb26 - t59) + d4*t183*t60 + d9*spb26*t183*t60 + d2*t128*t76 + t108*t171*t76 + d11*t36*(t16 + t58)*T(11) + d5*spa23*t36*(spa15*spb25 + spa16*spb26 - t59)*T(11); 
complex<T> t3 = spa12*t100*t124 + d24*spa23*spb56*t129 - d20*spa23*t13*t131 + d22*spa12*t131*t14 + t102*t128*t197 + d29*spb56*t129*t23 + d15*t299 + d18*t14*t36 + t207*t36 - d17*spa12*spa13*t37 + d20*spa12*t14*t37 + d19*t127*t246*t58 + d18*spa12*t36*t58 + spa23*t137*t170*(-(spa15*spb25) + t59) + t102*t127*t195*(-(spa15*spb25) + t59) + d19*t127*t223*(-(spa15*spb25) + t59) + d12*t155*t36*(-(spa15*spb25) + t59) + d27*spa12*t129*t61 + d28*spb56*t154*t61 + d26*t299*t61 + t100*t120*t7 + d25*spb56*t36*t7 - spa12*t123*t58*t73 + d29*t368*t98*t99 + d16*spa36*t57*T(11); 
complex<T> t4 = d38*spb24*t103*t127 + d24*spa23*spb56*t154 + d33*spa23*t11*t17 + t102*t127*t192 + d36*t11*t17*t21 + t127*t223*t228 + d29*spb56*t154*t23 - t127*t215*t235 - t184*t36 + t208*t36 + t107*t209*t36 - spa23*t137*t334*t59 + d42*spa34*t127*t38*t59 + d28*spb56*t129*t61 + d27*spa12*t154*t61 + d26*t102*t36*t61 + d33*spa12*t103*t64 + d35*t11*t17*t64 + spa34*t107*t36*t65 - d25*spb56*t36*t7 + d39*spa23*t246*t70 + t169*t70*t73 - d36*t209*t51*t97 - d29*t368*t98*t99 - spa12*t163*t36*T(11) - d40*spa23*t36*t59*T(11); 
complex<T> t5 = d20*spa23*t13*t131 + t137*t14*t170 + t108*t16*t171 + spa23*t183*t188 + spa12*t102*t14*t195 - t109*t183*t196 + d19*spa12*t14*t223 - t207*t36 + t208*t36 + d8*t245*t36 - d19*t127*t246*t58 - d18*spa12*t36*t58 + d22*t127*t131*(-(spa15*spb25) + t59) + d18*spa23*t36*(-(spa15*spb25) + t59) + d20*t127*t37*(-(spa15*spb25) + t59) - d4*t183*t60 - d9*spb26*t183*t60 + d49*spa12*t129*(-(spa23*spb25) + t61) + t129*t193*(-(spa23*spb25) + t61) + d49*spa23*t154*t62 + d50*t154*(-(spa23*spb25) + t61)*t62 + spa12*t123*t58*t73 + d3*t36*t76 - d51*t290*t62*t80 - t102*t126*t36*T(11) - t36*t56*t8*T(11); 
complex<T> t6 = -(t127*t223*t228) + d42*t106*t127*t235 + t127*t215*t235 + d39*spb24*t127*t246 + d15*t299 + t184*t36 + t127*t192*t42 + spa23*t137*t334*t59 + d49*spa12*t154*(-(spa23*spb25) + t61) + t154*t193*(-(spa23*spb25) + t61) + d49*spa23*t129*t62 + d50*t129*(-(spa23*spb25) + t61)*t62 + d46*spa12*t11*(spa23*spb24 - t64) + d46*spa23*t11*(spa12*spb24 + t65) + d38*spa23*t103*t70 + spa12*spb24*t169*t73 + d51*t290*t62*t80 + t102*t126*t36*T(11) + spa12*t163*t36*T(11) + d40*spa23*t36*t59*T(11) + d47*t11*t36*(-(spa23*spb24) + t64)*T(11) + d45*spa14*t36*(spa12*spb24 + t65)*T(11) + t36*t56*t8*T(11); 
complex<T> t151 = d71*t102*t118*t128 + d56*t127*t154*t86 + t100*t208*t86 + d72*t127*t154*t87; 
complex<T> t218 = d57*spb34*t102*t182 + s34*t182*t184 + t103*t117*t204 + s156*t182*t208 + d56*s156*t234 + d58*spb34*t234 + d55*s156*spa12*t103*t39 + s156*t184*t36*T(4); 
complex<T> t241 = t102*t119*t182 + t165*t234 + t183*t226*t39 + spa36*t210*t36*T(4) + s345*(t182*t208 + d56*t234 + d52*t183*t39 + t207*t36*T(4)); 
complex<T> t251 = t100*t102*t197 + t120*t128*t7; 
complex<T> co1 = d63*spb23*t280; 
complex<T> co2 = d66*spb23*spb34*t134; 
complex<T> co3 = d68*spb34*spb45*t134; 
complex<T> co4 = d70*spb45*t292; 
complex<T> co5 = d73*spb16*t292; 
complex<T> co6 = d74*spb16*t280; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t35*(*CI_users[1]->get_value(mc,ind,mu)) + t251*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + t6*(*CI_users[6]->get_value(mc,ind,mu)) + t5*(*CI_users[7]->get_value(mc,ind,mu)) + t202*(*CI_users[8]->get_value(mc,ind,mu)) + t180*(*CI_users[9]->get_value(mc,ind,mu)) + t149*(*CI_users[10]->get_value(mc,ind,mu)) + t218*(*CI_users[11]->get_value(mc,ind,mu)) + t241*(*CI_users[12]->get_value(mc,ind,mu)) + t233*(*CI_users[13]->get_value(mc,ind,mu)) + co1*(*CI_users[14]->get_value(mc,ind,mu)) + t203*(*CI_users[15]->get_value(mc,ind,mu)) + t150*(*CI_users[16]->get_value(mc,ind,mu)) + co2*(*CI_users[17]->get_value(mc,ind,mu)) + t181*(*CI_users[18]->get_value(mc,ind,mu)) + co3*(*CI_users[19]->get_value(mc,ind,mu)) + t219*(*CI_users[20]->get_value(mc,ind,mu)) + t96*(*CI_users[21]->get_value(mc,ind,mu)) + co4*(*CI_users[22]->get_value(mc,ind,mu)) + t151*(*CI_users[23]->get_value(mc,ind,mu)) + co5*(*CI_users[24]->get_value(mc,ind,mu)) + co6*(*CI_users[25]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_mpmppp_nf_wCI::\
C6g_mpmppp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c126, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c6, c2345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c1236));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c6, c345));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
CI_users.push_back(new Cached_Box_Integral_User(c5, c34, c2, c16));
} 
  
  
template <class T> SeriesC<T> 
     C6g_mpmppp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C6g :  mpmppp nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 17:53:30 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa16 = SPA(1,6);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa56 = SPA(5,6);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa26 = SPA(2,6);
complex<T> spa36 = SPA(3,6);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb16 = SPB(1,6);
complex<T> spb45 = SPB(4,5);
complex<T> spb26 = SPB(2,6);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb56 = SPB(5,6);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb46 = SPB(4,6);
complex<T> spb14 = SPB(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> s16 = -(spa16*spb16);
complex<T> s126 = SS(1,2,6);
complex<T> s156 = SS(1,5,6);
complex<T> s24 = -(spa24*spb24);
complex<T> s234 = SS(2,3,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = S(3,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s36 = S(3,6);
complex<T> s46 = S(4,6);
complex<T> s56 = -(spa56*spb56);
complex<T> s345 = SS(3,4,5);
complex<T> t7 = square(spa35); 
complex<T> t8 = square(spa15); 
complex<T> t10 = spa45*spa56; 
complex<T> t11 = square(spa14); 
complex<T> t16 = -(spa36*spb26); 
complex<T> t18 = spa16*T(3); 
complex<T> t20 = square(spa13*spb16 + spa23*spb26); 
complex<T> t21 = square(spa13*spb14 + spa35*spb45 + spa36*spb46); 
complex<T> t36 = square(spa13); 
complex<T> t37 = square(spa36); 
complex<T> t41 = cube(spa25); 
complex<T> t46 = square(s23 + s24); 
complex<T> t49 = square(spa34*spb24 + spa35*spb25); 
complex<T> t50 = square(spa15*spb45 + spa16*spb46); 
complex<T> t53 = (s12 - s126)*T(2); 
complex<T> t63 = -(spa34*spb24) - spa35*spb25; 
complex<T> t64 = spa12*spa23; 
complex<T> t65 = -(spa13*spb15) - spa36*spb56; 
complex<T> t66 = spa12*spb25 - spa16*spb56; 
complex<T> t67 = -(spa13*spb14) - spa35*spb45 - spa36*spb46; 
complex<T> t70 = -(spa13*spb15) - spa23*spb25 - spa36*spb56; 
complex<T> t72 = -s234 + s34; 
complex<T> t73 = -(spa13*spb16) - spa23*spb26; 
complex<T> t74 = s16 - s234; 
complex<T> t76 = -(spa15*spb45) - spa16*spb46; 
complex<T> t78 = s12 - s345; 
complex<T> t84 = square(spb26); 
complex<T> t85 = -s156 + s16; 
complex<T> t100 = spb12*T(3); 
complex<T> t102 = square(spa34); 
complex<T> t103 = square(spb56); 
complex<T> t107 = spa16*spa36; 
complex<T> t113 = s12*s16; 
complex<T> t119 = spb16*spb34; 
complex<T> t127 = spa15*spa35; 
complex<T> t128 = square(spa12); 
complex<T> t134 = square(spa23); 
complex<T> t135 = spa16*T(2); 
complex<T> t136 = spa12*spb24; 
complex<T> t140 = cube(spa15); 
complex<T> t141 = cube(spa35); 
complex<T> t151 = spa13*spb15 + spa36*spb56; 
complex<T> t159 = spa14*spb24; 
complex<T> t176 = square(spa16); 
complex<T> t177 = square(spb24); 
complex<T> t210 = spa14*spa34; 
complex<T> t214 = spa25*T(6); 
complex<T> t215 = spa26*spa45; 
complex<T> t228 = spa34*spa56; 
complex<T> t278 = spa15*spa16; 
complex<T> t17 = -(spa34*t76); 
complex<T> t19 = square(spa13*spb23 + t159); 
complex<T> t23 = square(t151); 
complex<T> t24 = square(spa23*spb25 + t151); 
complex<T> t38 = -t64; 
complex<T> t39 = -t159; 
complex<T> t40 = spa34*t10; 
complex<T> t42 = -t127; 
complex<T> t48 = square(spa13*spb23 + spa15*spb25 + t159); 
complex<T> t77 = spa16*t10; 
complex<T> t105 = -(t64*T(2)); 
complex<T> t106 = spa34*t11; 
complex<T> t111 = spa24*t10; 
complex<T> t126 = spa26*t78; 
complex<T> t131 = (-(spa23*spb25) + t65)*t66; 
complex<T> t154 = square(t66); 
complex<T> t155 = t7*t8; 
complex<T> t157 = spa12*t37; 
complex<T> t174 = spa24*t72; 
complex<T> t188 = spa25*t18; 
complex<T> t192 = spa14*t67; 
complex<T> t204 = t127*square(spa13); 
complex<T> t205 = spa16*t64; 
complex<T> t209 = spa16*t73; 
complex<T> t212 = spa12*t63; 
complex<T> t225 = t135*square(spa36); 
complex<T> t256 = spa15*t7; 
complex<T> t270 = spa35*t8; 
complex<T> t295 = t141*t65; 
complex<T> t308 = t128*t177; 
complex<T> d12 = (-s126 + s16)*spa25*spa34*t135*t215; d12 = T(1)/d12;
complex<T> d18 = (s16 - s345)*spa16*spa34*t214*t215; d18 = T(1)/d18;
complex<T> d59 = spa16*t228*square(spa25); d59 = T(1)/d59;
complex<T> d60 = spa16*t228*square(square(spa25)); d60 = T(1)/d60;
complex<T> d61 = t228*square(spa26); d61 = T(1)/d61;
complex<T> d62 = t228*square(square(spa26)); d62 = T(1)/d62;
complex<T> d66 = t10*square(spa25)*T(2); d66 = T(1)/d66;
complex<T> d67 = t10*square(square(spa25)); d67 = T(1)/d67;
complex<T> t15 = -t212; 
complex<T> t44 = -t209; 
complex<T> t62 = -(spa13*spb23) + t39; 
complex<T> t68 = -(spa13*spb23) - spa15*spb25 + t39; 
complex<T> t148 = d61*spb45; 
complex<T> t156 = spa16*t40; 
complex<T> t189 = spa12*t40; 
complex<T> t190 = t36*t42; 
complex<T> t243 = spa23*t40; 
complex<T> t288 = spa35*t24; 
complex<T> d1 = t40*square(spa26); d1 = T(1)/d1;
complex<T> d2 = spa26*t40*t53; d2 = T(1)/d2;
complex<T> d3 = t126*t40*T(6); d3 = T(1)/d3;
complex<T> d4 = t40*t78*cube(spa26); d4 = T(1)/d4;
complex<T> d5 = spa26*t100*t128*t40; d5 = T(1)/d5;
complex<T> d6 = spb12*t40*cube(spa26); d6 = T(1)/d6;
complex<T> d9 = t40*square(s36 + s46 + s56)*square(spa26); d9 = T(1)/d9;
complex<T> d10 = spa26*t40*cube(t78)*T(3); d10 = T(1)/d10;
complex<T> d17 = t40*cube(spa26); d17 = T(1)/d17;
complex<T> d20 = (s16 - s345)*t40*cube(spa26); d20 = T(1)/d20;
complex<T> d22 = t40*square(s23 + s24 + s25)*square(spa26); d22 = T(1)/d22;
complex<T> d23 = spa34*t188*t215*cube(s16 - s345); d23 = T(1)/d23;
complex<T> d24 = t40*t41*t74; d24 = T(1)/d24;
complex<T> d28 = t40*square(s25 + s35 + s45)*square(spa25); d28 = T(1)/d28;
complex<T> d30 = t77*square(spa24); d30 = T(1)/d30;
complex<T> d31 = spa23*t111*t18; d31 = T(1)/d31;
complex<T> d32 = t77*cube(spa24); d32 = T(1)/d32;
complex<T> d33 = (s23 - s234)*t77*cube(spa24); d33 = T(1)/d33;
complex<T> d34 = (s23 - s234)*spa23*t111*t18; d34 = T(1)/d34;
complex<T> d35 = t77*square(s24 + s34)*square(spa24); d35 = T(1)/d35;
complex<T> d36 = spa23*t111*t18*cube(s23 - s234); d36 = T(1)/d36;
complex<T> d37 = spa56*t174*t188; d37 = T(1)/d37;
complex<T> d38 = t72*t77*cube(spa24); d38 = T(1)/d38;
complex<T> d39 = t41*t72*t77; d39 = T(1)/d39;
complex<T> d40 = t174*t188*t228; d40 = T(1)/d40;
complex<T> d42 = t46*t77*square(spa24); d42 = T(1)/d42;
complex<T> d43 = t46*t77*square(spa25); d43 = T(1)/d43;
complex<T> d44 = spa24*spa56*t188*cube(t72); d44 = T(1)/d44;
complex<T> d45 = s34*spa12*t111*t18; d45 = T(1)/d45;
complex<T> d46 = spb34*t77*cube(spa24); d46 = T(1)/d46;
complex<T> d52 = t40*square(square(spa26)); d52 = T(1)/d52;
complex<T> d53 = t40*square(spa25); d53 = T(1)/d53;
complex<T> d54 = t40*square(square(spa25)); d54 = T(1)/d54;
complex<T> d55 = t77*square(square(spa24)); d55 = T(1)/d55;
complex<T> d57 = t77*square(spa25); d57 = T(1)/d57;
complex<T> d58 = t77*square(square(spa25)); d58 = T(1)/d58;
complex<T> d63 = t40*square(spa26)*T(2); d63 = T(1)/d63;
complex<T> d64 = t77*square(spa24)*T(2); d64 = T(1)/d64;
complex<T> d65 = t135*t40*square(spa25); d65 = T(1)/d65;
complex<T> t14 = spa23*(spa15*spb25 - t62); 
complex<T> t61 = -(d59*spb45*t204) + spa36*t148*square(spa13) + d62*spa16*spb45*t105*square(spa36) + d60*spb45*t155*t64*T(2); 
complex<T> t92 = -(d52*T(2)); 
complex<T> t109 = spb24*t62; 
complex<T> t114 = d55*s23; 
complex<T> t118 = d4*spb26; 
complex<T> t137 = d30*spa14; 
complex<T> t138 = d1*spa36; 
complex<T> t169 = d23*t48; 
complex<T> t171 = d38*t62; 
complex<T> t181 = d44*spa23; 
complex<T> t186 = d10*t20; 
complex<T> t193 = d24*spa23; 
complex<T> t194 = d38*spb24; 
complex<T> t197 = d23*t49; 
complex<T> t219 = d39*spb24; 
complex<T> t237 = t63*t68; 
complex<T> t240 = t113*(d63*spa36*t36 + d52*spa16*t37*t38); 
complex<T> d7 = spa26*t189*t53; d7 = T(1)/d7;
complex<T> d8 = t126*t189*T(6); d8 = T(1)/d8;
complex<T> d11 = spa26*t100*t189; d11 = T(1)/d11;
complex<T> d13 = spa25*t243*t85*T(2); d13 = T(1)/d13;
complex<T> d14 = spa25*t135*t243*t85; d14 = T(1)/d14;
complex<T> d15 = t156*square(spa25); d15 = T(1)/d15;
complex<T> d16 = spa26*t18*t243; d16 = T(1)/d16;
complex<T> d19 = (s16 - s345)*t156*t41; d19 = T(1)/d19;
complex<T> d21 = t156*square(s23 + s24 + s25)*square(spa25); d21 = T(1)/d21;
complex<T> d25 = t214*t243*t74; d25 = T(1)/d25;
complex<T> d26 = spa23*t156*t214*t74; d26 = T(1)/d26;
complex<T> d27 = t156*t41*t74; d27 = T(1)/d27;
complex<T> d29 = spa25*t243*cube(t74)*T(3); d29 = T(1)/d29;
complex<T> d41 = t156*t41*t72; d41 = T(1)/d41;
complex<T> d47 = s34*spa24*t18*t189; d47 = T(1)/d47;
complex<T> d48 = (s34 - s345)*t188*t189; d48 = T(1)/d48;
complex<T> d49 = (s34 - s345)*t156*t41; d49 = T(1)/d49;
complex<T> d50 = t156*square(s35 + s45)*square(spa25); d50 = T(1)/d50;
complex<T> d51 = t188*t189*cube(s34 - s345); d51 = T(1)/d51;
complex<T> d56 = t156*square(square(spa25)); d56 = T(1)/d56;
complex<T> t35 = t36*(d2*spa36*spb26 + d7*t209 - d12*(spa15*spa23*spb25 + t212 - spa23*t62)); 
complex<T> t82 = -t137; 
complex<T> t83 = -t138; 
complex<T> t86 = d15*spa15; 
complex<T> t117 = d56*s345; 
complex<T> t121 = -(d13*spb56); 
complex<T> t162 = d29*t23; 
complex<T> t187 = t105*t106*t114 + s23*t137*t36; 
complex<T> t196 = d51*t154; 
complex<T> t201 = d14*t65; 
complex<T> t208 = s34*(d64*s23*spa14*t36 + t106*t114*t38); 
complex<T> t223 = d49*t66; 
complex<T> t234 = d42*t109; 
complex<T> t236 = d27*t65; 
complex<T> t238 = d49*t70; 
complex<T> t242 = t205*t92; 
complex<T> t247 = d43*t109; 
complex<T> t262 = d41*t62; 
complex<T> t1 = t16*t176*t186 + d3*spa36*spb26*t36 + d2*t16*t36 + d8*t209*t36 + t118*t205*t37 + d4*t157*t44 + d9*spb26*t157*t44 + d7*t36*t44 + d5*spa23*t36*(-(spa15*spb25) - spa16*spb26 + t62) + d6*spa23*t37*(-(spa15*spb25) - spa16*spb26 + t62) + d6*spa23*t107*(t16 + t63) - d11*t36*(t16 + t63) + t36*t83 + d10*t157*t44*t84; 
complex<T> t2 = d33*spa23*t11*t17 + d36*t11*t17*t21 + d32*spa12*spa13*t210 + d34*t17*t36 - d34*t192*t36 - d36*t102*t192*t50 + d33*spa12*t106*t67 + d35*t11*t17*t67 + t36*t82 - d31*cube(spa13); 
complex<T> t3 = d44*t134*t136*t19 + d15*t190 - spa12*t155*t236 + d29*t103*t278*t295 + d37*t136*t36 + t137*t36 + d34*t192*t36 + t106*t194*t38 + t127*t247*t38 + t256*t262*t38 + d36*t102*t192*t50 + spa34*t181*t308*t62 + d40*spa23*t36*t62 + t171*t210*t64 + t219*t270*t64 + d42*spa34*t159*t62*t64 + d26*t190*t65 - d33*spa12*t106*t67 + spb56*(-(t155*(t162 + t193 - d28*t65)) + d25*t36*t7) + d33*spa23*t106*t76 + d36*t106*t21*t76 + d34*spa34*t36*t76 + d35*t106*t67*t76; 
complex<T> t6 = d20*spa23*t107*t15 - d50*t131*t155 + d20*t14*t157 + t134*t15*t169 + spa36*spb26*t176*t186 + d15*t190 + d4*t157*t209 + d9*spb26*t157*t209 + d22*spa36*t14*t212 - spa23*t155*t223 + t138*t36 + d18*t14*t36 + d3*t16*t36 + d18*t212*t36 + spa16*t118*t37*t38 + d8*t36*t44 + spa23*t128*t197*(-(spa15*spb25) + t62) + d19*t256*(-(spa15*spb25) + t62)*t64 + d19*t270*t63*t64 + d21*t127*(-(spa15*spb25) + t62)*t63*t64 + d49*spa12*t155*(-(spa23*spb25) + t65) + t155*t196*(-(spa23*spb25) + t65) + d48*t204*t66 + d51*t140*t288*t66 + d48*t36*(-(spa23*spb25) + t65)*t8 + d10*t157*t209*t84; 
complex<T> t58 = d54*spb16*t105*t155 + d53*spb16*t204 + s16*t138*square(spa13) + s16*t242*square(spa36); 
complex<T> t180 = spa35*t86; 
complex<T> t241 = d66*t119*t190 + d65*s234*s345*t127*t36 + s234*t117*t155*t38 + d67*t119*t155*t64; 
complex<T> t253 = s12*(t138*t36 + t242*t37); 
complex<T> t254 = t201*t204 + t121*t36*t7; 
complex<T> t4 = d17*spa13*t157 + spb56*t155*t162 + spb56*t155*t193 + t128*t14*t197 + t190*t201 + d21*t127*t14*t212 + t134*t169*t212 + spa12*t155*t236 + d19*spa12*t14*t256 + d19*spa23*t15*t270 - d29*t103*t278*t295 + d12*t14*t36 + d18*t15*t36 + t180*t36 + d12*t212*t36 + d18*spa23*t36*(-(spa15*spb25) + t62) + d20*t37*(-(spa15*spb25) + t62)*t64 + d20*t107*t63*t64 + d22*spa36*(-(spa15*spb25) + t62)*t63*t64 - d28*spb56*t155*t65 + d26*t204*t65 + d13*spb56*t36*t7 - d25*spb56*t36*t7 + t36*t83 - d16*spa36*cube(spa13); 
complex<T> t5 = d50*t131*t155 - d44*t134*t136*t19 + spa23*t155*t223 - d37*t136*t36 + t180*t36 + t171*t210*t38 + t219*t270*t38 - spa34*t181*t308*t62 - d40*spa23*t36*t62 + d42*spa34*t159*t38*t62 + t106*t194*t64 + t127*t247*t64 + t256*t262*t64 + d49*spa12*t155*(spa23*spb25 - t65) + t155*t196*(spa23*spb25 - t65) + d48*t190*t66 - d51*t140*t288*t66 + d47*t11*t36*(spa23*spb24 - t67) + d46*spa12*t11*(-(spa23*spb24) + t67) - d46*spa23*t11*(t136 + t76) - d45*spa14*t36*(t136 + t76) + d48*t36*(spa23*spb25 - t65)*t8 + t36*t82; 
complex<T> t59 = d56*s156*t105*t155 + d58*spb34*t105*t155 + d57*spb34*t204 + s34*t137*square(spa13) + s156*t180*square(spa13) + s156*t82*square(spa13) + d55*t106*(s34*t105 + s156*t64*T(2)); 
complex<T> t60 = d60*spb45*t105*t155 + t105*t117*t155 + d59*spb45*t204 + d52*s345*t225*t64 + d62*spb45*t225*t64 + (-(spa36*t148) + s345*(t180 + t83))*square(spa13); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t35*(*CI_users[1]->get_value(mc,ind,mu)) + t254*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t3*(*CI_users[5]->get_value(mc,ind,mu)) + t5*(*CI_users[6]->get_value(mc,ind,mu)) + t6*(*CI_users[7]->get_value(mc,ind,mu)) + t253*(*CI_users[8]->get_value(mc,ind,mu)) + t58*(*CI_users[9]->get_value(mc,ind,mu)) + t187*(*CI_users[10]->get_value(mc,ind,mu)) + t59*(*CI_users[11]->get_value(mc,ind,mu)) + t60*(*CI_users[12]->get_value(mc,ind,mu)) + t61*(*CI_users[13]->get_value(mc,ind,mu)) + t240*(*CI_users[14]->get_value(mc,ind,mu)) + t208*(*CI_users[15]->get_value(mc,ind,mu)) + t241*(*CI_users[16]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_mppmpp_G_wCI::\
C6g_mppmpp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c126, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c126));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c1236));
CI_users.push_back(new Cached_Bubble_Integral_User(c456, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c6, c2345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c1236));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c6, c345));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c2, c34, c5, c16));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c6, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c126));
CI_users.push_back(new Cached_Box_Integral_User(c3, c45, c6, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c6, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c34, c2, c16));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c126));
CI_users.push_back(new Cached_Box_Integral_User(c5, c6, c1, c234));
CI_users.push_back(new Cached_Box_Integral_User(c6, c1, c2, c345));
} 
  
  
template <class T> SeriesC<T> 
     C6g_mppmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C6g :  mppmpp G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 18:18:54 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa16 = SPA(1,6);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa56 = SPA(5,6);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa26 = SPA(2,6);
complex<T> spa46 = SPA(4,6);
complex<T> spa36 = SPA(3,6);
complex<T> spb12 = SPB(1,2);
complex<T> spb16 = SPB(1,6);
complex<T> spa25 = SPA(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb56 = SPB(5,6);
complex<T> spb35 = SPB(3,5);
complex<T> spb36 = SPB(3,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s16 = -(spa16*spb16);
complex<T> s123 = SS(1,2,3);
complex<T> s46 = S(4,6);
complex<T> s56 = -(spa56*spb56);
complex<T> s35 = -(spa35*spb35);
complex<T> s36 = -(spa36*spb36);
complex<T> s126 = SS(1,2,6);
complex<T> s156 = SS(1,5,6);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s25 = -(spa25*spb25);
complex<T> s234 = SS(2,3,4);
complex<T> s345 = SS(3,4,5);
complex<T> t8 = square(spa46); 
complex<T> t9 = square(spa24); 
complex<T> t10 = square(spa13); 
complex<T> t11 = square(spa15); 
complex<T> t12 = spa23*spa56; 
complex<T> t15 = spa16*spa34; 
complex<T> t16 = spa12*spa45; 
complex<T> t24 = square(spa14*spb16 + spa24*spb26); 
complex<T> t25 = square(spa14*spb34 + spa15*spb35); 
complex<T> t27 = square(spa14*spb34 + spa15*spb35 + spa16*spb36); 
complex<T> t28 = square(spa14*spb16 + spa24*spb26 + spa34*spb36); 
complex<T> t29 = square(spa45*spb35 + spa46*spb36); 
complex<T> t44 = square(spa14); 
complex<T> t45 = -(spa24*spb26); 
complex<T> t46 = -(spa15*spb35); 
complex<T> t48 = cube(spa36); 
complex<T> t49 = -(spa46*T(2)); 
complex<T> t50 = cube(spa25); 
complex<T> t69 = spa34*spb23 - spa45*spb25; 
complex<T> t70 = cube(spa14); 
complex<T> t71 = spa12*spb25 - spa16*spb56; 
complex<T> t74 = -(spa14*spb15) - spa46*spb56; 
complex<T> t75 = -(spa14*spb15) - spa24*spb25 - spa46*spb56; 
complex<T> t77 = spa12*spb26 + spa13*spb36; 
complex<T> t78 = -(spa45*spb35) - spa46*spb36; 
complex<T> t79 = spa12*spa16; 
complex<T> t80 = -(spa13*spb23) - spa14*spb24; 
complex<T> t81 = spa34*spa45; 
complex<T> t83 = s12 - s345; 
complex<T> t84 = -s345 + s45; 
complex<T> t86 = -s234 + s34; 
complex<T> t89 = s12 - s123; 
complex<T> t96 = s12 - s126; 
complex<T> t97 = -s156 + s16; 
complex<T> t104 = s234*s345; 
complex<T> t105 = cube(spa15); 
complex<T> t106 = cube(spa24); 
complex<T> t116 = square(spb23); 
complex<T> t117 = square(spb56); 
complex<T> t121 = -(spa13*spa46); 
complex<T> t125 = square(spa34); 
complex<T> t126 = spa26*spa36; 
complex<T> t128 = square(spb26); 
complex<T> t129 = square(spb35); 
complex<T> t141 = spb16*spb34; 
complex<T> t157 = spa15*T(11); 
complex<T> t159 = square(square(spa14)); 
complex<T> t161 = square(spa12); 
complex<T> t165 = square(spa16); 
complex<T> t166 = square(spa45); 
complex<T> t168 = cube(spa13); 
complex<T> t177 = spa13*spb23 + spa14*spb24; 
complex<T> t182 = spa14*spb15 + spa46*spb56; 
complex<T> t189 = spa46*T(2); 
complex<T> t193 = spa36*T(3); 
complex<T> t199 = spa25*spa26; 
complex<T> t214 = cube(spa46); 
complex<T> t226 = spa45*spb56; 
complex<T> t234 = spb12*spb45; 
complex<T> t251 = spa12*spb23; 
complex<T> t270 = -(spa46*T(4)); 
complex<T> t315 = -(s345*T(2)); 
complex<T> t331 = spa12*spa35; 
complex<T> t400 = spa13*spa16; 
complex<T> t414 = spa15*spa16; 
complex<T> t21 = square(t177); 
complex<T> t22 = square(spa15*spb25 + t177); 
complex<T> t31 = square(t182); 
complex<T> t32 = square(spa24*spb25 + t182); 
complex<T> t47 = -t79; 
complex<T> t62 = (s16 - s234)*t12; 
complex<T> t72 = -(spa14*spb16) + t45; 
complex<T> t73 = -(spa15*spb25) + t80; 
complex<T> t76 = -(spa14*spb34) + t46; 
complex<T> t82 = -(spa14*spb34) - spa16*spb36 + t46; 
complex<T> t85 = -(spa14*spb16) - spa34*spb36 + t45; 
complex<T> t95 = t12*T(6); 
complex<T> t107 = -t141; 
complex<T> t108 = -t234; 
complex<T> t118 = -(t44*T(2)); 
complex<T> t120 = t12*t16; 
complex<T> t151 = spa35*t84; 
complex<T> t154 = spa24*t44; 
complex<T> t163 = spb26*t8; 
complex<T> t167 = square(t77); 
complex<T> t169 = square(t69); 
complex<T> t171 = square(t71); 
complex<T> t190 = t12*t15; 
complex<T> t223 = spa13*t44; 
complex<T> t224 = t11*t9; 
complex<T> t225 = -(t15*T(2)); 
complex<T> t227 = spa12*t69; 
complex<T> t229 = spa56*t126; 
complex<T> t243 = spa25*t86; 
complex<T> t248 = t10*t8; 
complex<T> t249 = t15*T(2); 
complex<T> t288 = -(t16*T(2)); 
complex<T> t289 = t10*t11; 
complex<T> t298 = t116*t80; 
complex<T> t303 = t16*T(2); 
complex<T> t304 = spa24*t11; 
complex<T> t305 = spa46*t10; 
complex<T> t314 = t79*t9; 
complex<T> t325 = spa13*t15; 
complex<T> t330 = spa45*t15; 
complex<T> t338 = spa13*t11; 
complex<T> t342 = spa23*t199; 
complex<T> t344 = spa15*t9; 
complex<T> t351 = t11*t75; 
complex<T> t358 = t15*t331; 
complex<T> t363 = spa16*t44; 
complex<T> t374 = spa46*t9; 
complex<T> t385 = spa15*t10; 
complex<T> t386 = spb23*t159; 
complex<T> t423 = spa13*t125; 
complex<T> t436 = t165*t28; 
complex<T> t439 = t166*t74; 
complex<T> d1 = t12*t81*square(spa26); d1 = T(1)/d1;
complex<T> d3 = t12*t81*t83*cube(spa26); d3 = T(1)/d3;
complex<T> d4 = spa45*t12*t48*t83; d4 = T(1)/d4;
complex<T> d8 = spb12*t12*t81*cube(spa26); d8 = T(1)/d8;
complex<T> d12 = t12*t81*square(s36 + s46 + s56)*square(spa26); d12 = T(1)/d12;
complex<T> d13 = spa45*t12*square(s36 + s46 + s56)*square(spa36); d13 = T(1)/d13;
complex<T> d22 = spa25*spa34*t12*t97; d22 = T(1)/d22;
complex<T> d26 = t12*t81*cube(spa26); d26 = T(1)/d26;
complex<T> d29 = (s16 - s345)*t12*t81*cube(spa26); d29 = T(1)/d29;
complex<T> d31 = t12*t81*square(s23 + s24 + s25)*square(spa26); d31 = T(1)/d31;
complex<T> d37 = spa34*t12*square(s25 + s35 + s45)*square(spa25); d37 = T(1)/d37;
complex<T> d38 = spa25*spa34*t12*cube(s16 - s234)*T(3); d38 = T(1)/d38;
complex<T> d39 = spa16*t12*t50*t86; d39 = T(1)/d39;
complex<T> d43 = spa16*t12*square(s23 + s24)*square(spa25); d43 = T(1)/d43;
complex<T> d44 = spa16*spa25*t12*cube(t86)*T(3); d44 = T(1)/d44;
complex<T> d45 = t12*t79*square(spa35); d45 = T(1)/d45;
complex<T> d47 = t12*t79*cube(spa35); d47 = T(1)/d47;
complex<T> d49 = (s34 - s345)*t12*t79*cube(spa35); d49 = T(1)/d49;
complex<T> d52 = t12*t79*square(s35 + s45)*square(spa35); d52 = T(1)/d52;
complex<T> d55 = t12*t79*t84*cube(spa35); d55 = T(1)/d55;
complex<T> d56 = spa12*t12*t48*t84; d56 = T(1)/d56;
complex<T> d59 = t12*t79*square(s34 + s35)*square(spa35); d59 = T(1)/d59;
complex<T> d60 = spa12*t12*square(s34 + s35)*square(spa36); d60 = T(1)/d60;
complex<T> d61 = spa23*spa35*t193*t79*cube(t84); d61 = T(1)/d61;
complex<T> d67 = spb45*t12*t79*cube(spa35); d67 = T(1)/d67;
complex<T> d69 = t12*t81*square(square(spa26)); d69 = T(1)/d69;
complex<T> d70 = spa45*t12*square(spa36); d70 = T(1)/d70;
complex<T> d71 = spa45*t12*square(square(spa36)); d71 = T(1)/d71;
complex<T> d73 = spa56*t16*square(spa36); d73 = T(1)/d73;
complex<T> d74 = spa56*t16*square(square(spa36)); d74 = T(1)/d74;
complex<T> d75 = spa34*t12*square(spa25); d75 = T(1)/d75;
complex<T> d76 = spa34*t12*square(square(spa25)); d76 = T(1)/d76;
complex<T> d78 = spa16*t12*square(spa25); d78 = T(1)/d78;
complex<T> d79 = spa16*t12*square(square(spa25)); d79 = T(1)/d79;
complex<T> d80 = t12*t79*square(square(spa35)); d80 = T(1)/d80;
complex<T> d81 = spa34*t12*square(spa26); d81 = T(1)/d81;
complex<T> d82 = spa12*t12*square(spa36); d82 = T(1)/d82;
complex<T> d83 = spa34*t12*square(square(spa26)); d83 = T(1)/d83;
complex<T> d84 = spa12*t12*square(square(spa36)); d84 = T(1)/d84;
complex<T> d90 = t12*square(spa36); d90 = T(1)/d90;
complex<T> d91 = t12*square(square(spa36)); d91 = T(1)/d91;
complex<T> d92 = t12*t79*T(2); d92 = T(1)/d92;
complex<T> d95 = t12*square(spa25); d95 = T(1)/d95;
complex<T> d96 = t12*square(square(spa25)); d96 = T(1)/d96;
complex<T> d98 = t12*t81*T(2); d98 = T(1)/d98;
complex<T> t131 = d1*s12; 
complex<T> t133 = d69*s16; 
complex<T> t136 = d45*s34; 
complex<T> t142 = d81*spb45; 
complex<T> t143 = d22*spb56; 
complex<T> t148 = d55*t76; 
complex<T> t175 = spa26*t95; 
complex<T> t186 = d61*t129; 
complex<T> t197 = d80*s34; 
complex<T> t202 = d83*spb45; 
complex<T> t208 = d61*t25; 
complex<T> t215 = t72*t79; 
complex<T> t232 = d73*spb23; 
complex<T> t250 = spa15*t154; 
complex<T> t254 = d60*spb35; 
complex<T> t256 = spa36*t151; 
complex<T> t262 = d59*t76; 
complex<T> t277 = d74*spb23; 
complex<T> t280 = d13*t72; 
complex<T> t294 = d49*t71; 
complex<T> t307 = spa34*t82; 
complex<T> t316 = d45*spa15; 
complex<T> t317 = d1*spa46; 
complex<T> t337 = t314*t8; 
complex<T> t353 = d29*spa12; 
complex<T> t354 = d52*spa45; 
complex<T> t355 = spa25*t62; 
complex<T> t359 = d3*spa24; 
complex<T> t362 = t224*t303; 
complex<T> t367 = t229*t83; 
complex<T> t368 = t289*t81; 
complex<T> t387 = d31*spa12; 
complex<T> t406 = spa24*t298; 
complex<T> d2 = t120*square(spa36); d2 = T(1)/d2;
complex<T> d5 = t229*t81*t96; d5 = T(1)/d5;
complex<T> d9 = spa34*t16*t229*t96; d9 = T(1)/d9;
complex<T> d11 = t120*t48*t83; d11 = T(1)/d11;
complex<T> d14 = t229*t81*cube(t83)*T(3); d14 = T(1)/d14;
complex<T> d16 = spa36*t16*t89*t95; d16 = T(1)/d16;
complex<T> d17 = t120*t48*t89; d17 = T(1)/d17;
complex<T> d18 = t120*square(s34 + s35 + s36)*square(spa36); d18 = T(1)/d18;
complex<T> d19 = t120*t193*cube(t89); d19 = T(1)/d19;
complex<T> d20 = spa36*t120*t89; d20 = T(1)/d20;
complex<T> d21 = (-s126 + s16)*t330*t342; d21 = T(1)/d21;
complex<T> d23 = spa25*t190*t97; d23 = T(1)/d23;
complex<T> d24 = t190*square(spa25); d24 = T(1)/d24;
complex<T> d27 = (s16 - s345)*t330*t342*T(6); d27 = T(1)/d27;
complex<T> d28 = (s16 - s345)*t190*t50; d28 = T(1)/d28;
complex<T> d30 = t190*square(s23 + s24 + s25)*square(spa25); d30 = T(1)/d30;
complex<T> d32 = t330*t342*cube(s16 - s345)*T(3); d32 = T(1)/d32;
complex<T> d33 = spa34*t50*t62; d33 = T(1)/d33;
complex<T> d36 = t15*t50*t62; d36 = T(1)/d36;
complex<T> d40 = spa16*t243*t95; d40 = T(1)/d40;
complex<T> d41 = t15*t243*t95; d41 = T(1)/d41;
complex<T> d42 = t190*t50*t86; d42 = T(1)/d42;
complex<T> d46 = t358*t95; d46 = T(1)/d46;
complex<T> d48 = (s34 - s345)*t190*t50; d48 = T(1)/d48;
complex<T> d50 = (s34 - s345)*spa25*spa56*t358*T(6); d50 = T(1)/d50;
complex<T> d51 = t190*square(s35 + s45)*square(spa25); d51 = T(1)/d51;
complex<T> d53 = spa25*spa56*t358*cube(s34 - s345)*T(3); d53 = T(1)/d53;
complex<T> d58 = t120*t48*t84; d58 = T(1)/d58;
complex<T> d62 = (-s123 + s45)*spa36*t16*t95; d62 = T(1)/d62;
complex<T> d63 = (-s123 + s45)*t120*t48; d63 = T(1)/d63;
complex<T> d64 = t120*square(s46 + s56)*square(spa36); d64 = T(1)/d64;
complex<T> d65 = t120*t193*cube(-s123 + s45); d65 = T(1)/d65;
complex<T> d66 = s45*spa35*t79*t95; d66 = T(1)/d66;
complex<T> d68 = s45*spa16*spa35*t16*t95; d68 = T(1)/d68;
complex<T> d72 = t120*square(square(spa36)); d72 = T(1)/d72;
complex<T> d77 = t190*square(square(spa25)); d77 = T(1)/d77;
complex<T> d85 = spa45*spa56*t249; d85 = T(1)/d85;
complex<T> d86 = t120*t249; d86 = T(1)/d86;
complex<T> d87 = t16*t249; d87 = T(1)/d87;
complex<T> d88 = spa16*spa56*t303; d88 = T(1)/d88;
complex<T> d89 = t120*T(2); d89 = T(1)/d89;
complex<T> d93 = t190*T(2); d93 = T(1)/d93;
complex<T> d94 = spa12*spa23*t249; d94 = T(1)/d94;
complex<T> d97 = spa23*spa34*t303; d97 = T(1)/d97;
complex<T> t66 = -(d16*(spa14*spb34 + spa15*spb35 + spa16*spb36)); 
complex<T> t109 = d86*s123; 
complex<T> t123 = -(d24*spa15); 
complex<T> t132 = d2*s123; 
complex<T> t139 = d77*s45; 
complex<T> t149 = d62*t77; 
complex<T> t196 = d72*s123; 
complex<T> t198 = d24*s45; 
complex<T> t205 = d65*t167; 
complex<T> t207 = d14*t24; 
complex<T> t218 = d14*t128; 
complex<T> t219 = (d86*t104 + d89*t107)*t159; 
complex<T> t222 = t197*t368*T(2) - spa15*t136*t223*T(4); 
complex<T> t235 = d53*t171; 
complex<T> t238 = d23*t74; 
complex<T> t242 = d41*t80; 
complex<T> t257 = d32*t169; 
complex<T> t265 = d17*t82; 
complex<T> t268 = t225*t248*t277 + spa46*t223*t232*T(4); 
complex<T> t269 = d24*spa15*spa24*t104*t118 + d77*t104*t16*t224 + d96*t107*t16*t224 + d95*t141*t250*T(2); 
complex<T> t282 = t186*t76; 
complex<T> t283 = d64*t77; 
complex<T> t295 = d9*t72; 
complex<T> t297 = d18*t78; 
complex<T> t301 = s45*(spa13*spa15*t118*t136 + t197*t368); 
complex<T> t302 = s12*t133*t337 + s16*t131*t154*t49; 
complex<T> t312 = d20*(t118*t307 + t223*t78*T(2)); 
complex<T> t321 = d71*spb12*t248*t249 + t131*t154*t270 + d70*spb12*t223*t270 + d69*s12*t337*T(2); 
complex<T> t329 = -(t250*T(4)); 
complex<T> t341 = d51*spa45; 
complex<T> t348 = spa23*t256; 
complex<T> t365 = d2*spa46; 
complex<T> t375 = d30*spa12; 
complex<T> t399 = d28*spa15; 
complex<T> t441 = d2*t121; 
complex<T> d6 = t367*t81*T(6); d6 = T(1)/d6;
complex<T> d7 = spb12*t161*t175*t81; d7 = T(1)/d7;
complex<T> d10 = spa34*t16*t367*T(6); d10 = T(1)/d10;
complex<T> d15 = spa34*spb12*t16*t175; d15 = T(1)/d15;
complex<T> d25 = t175*t330; d25 = T(1)/d25;
complex<T> d34 = spa34*t355*T(6); d34 = T(1)/d34;
complex<T> d35 = t15*t355*T(6); d35 = T(1)/d35;
complex<T> t1 = t223*t365 + spa34*t44*t66 + d63*spa16*t248*(spa34*spb36 - t72) + d17*spa34*t248*(spa16*spb36 - t76) + d19*spa34*t10*t29*(spa16*spb36 - t76) + spa34*t297*t305*(spa16*spb36 - t76) - d65*spa46*t436*t77 - d16*t223*t78 - d17*t15*t305*t78 + d19*t27*t423*t78 + spa16*t205*(spa34*spb36 - t72)*t8 + t283*t400*(-(spa34*spb36) + t72)*t8 + d63*t325*t77*t8 - spa46*t149*t44*T(11) + d62*t363*(spa34*spb36 - t72)*T(11); 
complex<T> t2 = t163*t165*t207 + d12*spa24*t163*t215 + t214*t215*t218 + t154*t317 + d4*t163*t325 + t189*t295*t363 - t163*t280*t400 - d6*t163*t44 + t44*t441 + d8*spa16*t374*(spa46*spb26 - t69) - d11*spa16*t248*t72 - d10*spa46*t363*t72 + d17*spa34*t248*(-(spa16*spb36) + t76) + d19*spa34*t10*t29*(-(spa16*spb36) + t76) + spa34*t297*t305*(-(spa16*spb36) + t76) + d17*t15*t305*t78 - d19*t27*t423*t78 + t215*t359*t8 + d3*t163*t47*t9 + d8*(spa16*spb26 - t73)*t8*t9 + d5*t163*t44*T(2) + spa34*t44*t66*T(11) + d15*t154*(-(spa46*spb26) + t69)*T(11) - d16*t223*t78*T(11) + d7*t44*(spa16*spb26 - t73)*t9*T(11); 
complex<T> t4 = spa45*t118*t143 + t123*t154 - d33*t224*t226 + d21*spa24*t118*t227 + d27*t154*t227 - d32*t106*t22*t227 + d28*t224*t227 - d38*t11*t226*t31 + t154*t317 + d38*t117*t414*t439 + d34*t226*t44 + d29*t374*t47*t69 + t344*t375*t69*t73 - t374*t387*t69*t73 - d36*t16*t304*t74 + d37*t226*t304*t74 - d35*spa15*t44*t74 - d26*spa12*spa14*spa24*t8 + t161*t257*t73*t9 + t16*t399*t73*t9 - d27*t44*t73*t9 - t353*t73*t8*t9 + spa15*t238*t44*T(2) + d21*t44*t73*t9*T(2) + d25*spa46*t70*T(11); 
complex<T> t5 = d33*t224*t226 + d24*t250 + d39*t224*t251 + d38*t11*t226*t31 - d44*spa34*t161*t406 - d38*t117*t414*t439 - d34*t226*t44 + d36*t16*t304*t74 - d37*t226*t304*t74 + d35*spa15*t44*t74 + d42*t16*t344*t80 - d43*t251*t344*t80 + d44*t21*t251*t9 - t154*t242*T(11) + d40*t251*t44*T(11); 
complex<T> t7 = t123*t154 - d39*t224*t251 + t223*t316 - d47*spa14*spa45*t385 + d44*spa34*t161*t406 - d48*spa45*t224*t71 + d53*spa45*t105*t32*t71 + d50*spa45*t157*t44*t71 + d49*spa45*t289*(spa24*spb25 - t74) + t304*t341*t71*(spa24*spb25 - t74) + t11*t166*t235*(-(spa24*spb25) + t74) + d48*t16*t304*(-(spa24*spb25) + t74) + t338*t354*t71*(-(spa24*spb25) + t74) - d42*t16*t344*t80 + d43*t251*t344*t80 + t294*t338*t81 - d44*t21*t251*t9 + t154*t242*T(11) - d40*t251*t44*T(11) + d46*spa13*t70*T(11) + d50*t11*t44*(-(spa24*spb25) + t74)*T(11); 
complex<T> t43 = t295*t363*t49 + d21*t118*t73*square(spa24) + d5*spb26*t118*square(spa46) + d21*t154*t227*T(2); 
complex<T> t245 = s234*t109*t159 - d87*spb56*t386; 
complex<T> t246 = d24*s156*t329 + d78*spb34*t329 + d77*s156*t362 + d79*spb34*t362; 
complex<T> t247 = d90*t189*t223*t234 + d91*t108*t15*t248 + s126*t15*t196*t248 + s126*t132*t223*t49; 
complex<T> t267 = (d93*t108 + s345*t109)*t159; 
complex<T> t287 = t196*t248*t249 + t132*t223*t270 + t223*t232*t270 + t248*t249*t277; 
complex<T> t300 = spa15*t118*t238 + spa45*t143*t44*T(2); 
complex<T> t328 = d1*s16*t154*t270 + d75*spb16*t329 + d76*spb16*t362 + t133*t337*T(2); 
complex<T> t336 = t142*t154*t270 + t198*t329 + t139*t362 + t202*t337*T(2); 
complex<T> t343 = d72*s345*t248*t249 + d84*spb45*t248*t249 + d2*s345*t223*t270 + d82*spb45*t223*t270 + t139*t224*t288 + d24*s345*t329 + d69*t315*t337 + d77*s345*t362 + d80*t315*t368 - t202*t337*T(2) + d80*s45*t368*T(2) + spa46*t142*t154*T(4) + t198*t250*T(4) + s345*t223*t316*T(4) - s45*t223*t316*T(4) + s345*t154*t317*T(4); 
complex<T> t452 = spa13*t283; 
complex<T> d54 = t348*t79*T(6); d54 = T(1)/d54;
complex<T> d57 = spa16*t16*t348*T(6); d57 = T(1)/d57;
complex<T> t201 = d54*spb35; 
complex<T> t239 = d57*t76; 
complex<T> t3 = spb35*t10*t125*t208 + d56*spb35*t15*t305 + t223*t316 - d55*spb35*t368 + t44*t441 + d67*spa34*t338*(spa13*spb35 + t71) + d66*t157*t44*(spa13*spb35 + t71) + d63*spa16*t248*(-(spa34*spb36) + t72) + d67*t289*(spa24*spb25 + spa34*spb35 - t74) - d58*spa34*t248*t76 - spa34*t254*t305*t76 + d65*spa46*t436*t77 + t283*t400*(spa34*spb36 - t72)*t8 + spa16*t205*(-(spa34*spb36) + t72)*t8 - d63*t325*t77*t8 + t168*t282*t81 + t148*t385*t81 + spb35*t262*t385*t81 + spa34*t223*t239*T(11) + spa46*t149*t44*T(11) + t10*t201*t44*T(11) + d62*t363*(-(spa34*spb36) + t72)*T(11) - d68*t11*t44*(spa24*spb25 + spa34*spb35 - t74)*T(11); 
complex<T> t6 = -(t163*t165*t207) - spb35*t10*t125*t208 - d27*t154*t227 + d32*t106*t22*t227 - d28*t224*t227 + d24*t250 - d56*spb35*t15*t305 + d3*t163*t314 - t223*t316 - t154*t317 - d4*t163*t325 + t223*t365 + d55*spb35*t368 + t163*t280*t400 + d6*t163*t44 + d29*spa46*t314*t69 + d48*spa45*t224*t71 - d53*spa45*t105*t32*t71 + d11*spa16*t248*t72 + d10*spa46*t363*t72 + t214*t218*t47*t72 - t344*t375*t69*t73 + t374*t387*t69*t73 + t11*t166*t235*(spa24*spb25 - t74) + d48*t16*t304*(spa24*spb25 - t74) + t338*t354*t71*(spa24*spb25 - t74) + d49*spa45*t289*(-(spa24*spb25) + t74) + t304*t341*t71*(-(spa24*spb25) + t74) + d58*spa34*t248*t76 + spa34*t254*t305*t76 + d12*t215*t45*t8 + t359*t47*t72*t8 - t168*t282*t81 - t294*t338*t81 - t148*t385*t81 + t10*t262*t46*t81 - t161*t257*t73*t9 - t16*t399*t73*t9 + d27*t44*t73*t9 + t353*t73*t8*t9 - spa34*t223*t239*T(11) - t10*t201*t44*T(11) - d50*spa15*spa45*t44*t71*T(11) + d50*t11*t44*(spa24*spb25 - t74)*T(11); 
complex<T> co1 = d85*spb12*t386; 
complex<T> co2 = d88*spb34*t386; 
complex<T> co3 = d92*spb34*spb45*t159; 
complex<T> co4 = d94*spb45*spb56*t159; 
complex<T> co5 = d97*spb16*spb56*t159; 
complex<T> co6 = d98*spb12*spb16*t159; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t312*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t300*(*CI_users[3]->get_value(mc,ind,mu)) + t4*(*CI_users[4]->get_value(mc,ind,mu)) + t5*(*CI_users[5]->get_value(mc,ind,mu)) + t7*(*CI_users[6]->get_value(mc,ind,mu)) + t6*(*CI_users[7]->get_value(mc,ind,mu)) + t3*(*CI_users[8]->get_value(mc,ind,mu)) + t1*(*CI_users[9]->get_value(mc,ind,mu)) + t321*(*CI_users[10]->get_value(mc,ind,mu)) + t287*(*CI_users[11]->get_value(mc,ind,mu)) + t328*(*CI_users[12]->get_value(mc,ind,mu)) + t268*(*CI_users[13]->get_value(mc,ind,mu)) + t246*(*CI_users[14]->get_value(mc,ind,mu)) + t222*(*CI_users[15]->get_value(mc,ind,mu)) + t343*(*CI_users[16]->get_value(mc,ind,mu)) + t336*(*CI_users[17]->get_value(mc,ind,mu)) + co1*(*CI_users[18]->get_value(mc,ind,mu)) + t245*(*CI_users[19]->get_value(mc,ind,mu)) + t302*(*CI_users[20]->get_value(mc,ind,mu)) + co2*(*CI_users[21]->get_value(mc,ind,mu)) + t219*(*CI_users[22]->get_value(mc,ind,mu)) + t247*(*CI_users[23]->get_value(mc,ind,mu)) + co3*(*CI_users[24]->get_value(mc,ind,mu)) + t267*(*CI_users[25]->get_value(mc,ind,mu)) + co4*(*CI_users[26]->get_value(mc,ind,mu)) + t269*(*CI_users[27]->get_value(mc,ind,mu)) + t301*(*CI_users[28]->get_value(mc,ind,mu)) + co5*(*CI_users[29]->get_value(mc,ind,mu)) + co6*(*CI_users[30]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_mppmpp_nf_wCI::\
C6g_mppmpp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c126, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c126));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c1236));
CI_users.push_back(new Cached_Bubble_Integral_User(c456, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c6, c2345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c1236));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c6, c345));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c6, c45));
CI_users.push_back(new Cached_Box_Integral_User(c5, c34, c2, c16));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c126));
} 
  
  
template <class T> SeriesC<T> 
     C6g_mppmpp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C6g :  mppmpp nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 18:37:59 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa16 = SPA(1,6);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa56 = SPA(5,6);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa26 = SPA(2,6);
complex<T> spa46 = SPA(4,6);
complex<T> spa36 = SPA(3,6);
complex<T> spb12 = SPB(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spb16 = SPB(1,6);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb36 = SPB(3,6);
complex<T> spb26 = SPB(2,6);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb56 = SPB(5,6);
complex<T> spb15 = SPB(1,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s16 = -(spa16*spb16);
complex<T> s123 = SS(1,2,3);
complex<T> s46 = S(4,6);
complex<T> s56 = -(spa56*spb56);
complex<T> s35 = -(spa35*spb35);
complex<T> s36 = -(spa36*spb36);
complex<T> s126 = SS(1,2,6);
complex<T> s156 = SS(1,5,6);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s25 = -(spa25*spb25);
complex<T> s234 = SS(2,3,4);
complex<T> s345 = SS(3,4,5);
complex<T> t8 = square(spa46); 
complex<T> t9 = square(spa24); 
complex<T> t10 = square(spa13); 
complex<T> t11 = square(spa15); 
complex<T> t12 = spa23*spa56; 
complex<T> t13 = square(spa36); 
complex<T> t16 = spa16*spa34; 
complex<T> t17 = spa12*spa45; 
complex<T> t24 = square(spa14*spb16 + spa24*spb26); 
complex<T> t25 = square(spa14*spb34 + spa15*spb35); 
complex<T> t27 = square(spa14*spb34 + spa15*spb35 + spa16*spb36); 
complex<T> t28 = square(spa14*spb16 + spa24*spb26 + spa34*spb36); 
complex<T> t29 = square(spa45*spb35 + spa46*spb36); 
complex<T> t44 = square(spa14); 
complex<T> t45 = -(spa24*spb26); 
complex<T> t46 = -(spa15*spb35); 
complex<T> t48 = cube(spa36); 
complex<T> t49 = cube(spa25); 
complex<T> t70 = spa34*spb23 - spa45*spb25; 
complex<T> t71 = spa12*spb25 - spa16*spb56; 
complex<T> t74 = -(spa14*spb15) - spa46*spb56; 
complex<T> t75 = -(spa14*spb15) - spa24*spb25 - spa46*spb56; 
complex<T> t77 = spa12*spb26 + spa13*spb36; 
complex<T> t78 = -(spa45*spb35) - spa46*spb36; 
complex<T> t79 = spa12*spa16; 
complex<T> t80 = -(spa13*spb23) - spa14*spb24; 
complex<T> t82 = spa34*spa45; 
complex<T> t83 = s12 - s345; 
complex<T> t84 = -s345 + s45; 
complex<T> t85 = -(spa15*spa24); 
complex<T> t86 = spa13*spa46; 
complex<T> t89 = s16 - s234; 
complex<T> t90 = -s234 + s34; 
complex<T> t91 = -square(spa14); 
complex<T> t97 = s12 - s126; 
complex<T> t98 = -s156 + s16; 
complex<T> t115 = square(spb23); 
complex<T> t116 = square(spb56); 
complex<T> t123 = spa36*T(3); 
complex<T> t126 = spa25*spa26; 
complex<T> t137 = spb16*spb34; 
complex<T> t153 = square(spa45); 
complex<T> t162 = square(spa16); 
complex<T> t163 = square(spa34); 
complex<T> t169 = square(spb35); 
complex<T> t174 = spa13*spb23 + spa14*spb24; 
complex<T> t178 = spa14*spb15 + spa46*spb56; 
complex<T> t189 = spa25*T(3); 
complex<T> t190 = square(spb26); 
complex<T> t214 = square(spa12); 
complex<T> t217 = cube(spa15); 
complex<T> t236 = spa26*spa36; 
complex<T> t251 = cube(spa13); 
complex<T> t252 = cube(spa24); 
complex<T> t253 = cube(spa46); 
complex<T> t254 = s45*spa15; 
complex<T> t267 = spa36*T(6); 
complex<T> t270 = spb12*spb45; 
complex<T> t305 = spa13*spa16; 
complex<T> t372 = spa13*spa34; 
complex<T> t411 = spa24*spa34; 
complex<T> t21 = square(t174); 
complex<T> t22 = square(spa15*spb25 + t174); 
complex<T> t31 = square(t178); 
complex<T> t32 = square(spa24*spb25 + t178); 
complex<T> t47 = -t79; 
complex<T> t50 = -t82; 
complex<T> t52 = t12*T(3); 
complex<T> t72 = -(spa15*spb25) + t80; 
complex<T> t73 = -(spa14*spb16) + t45; 
complex<T> t76 = -(spa14*spb34) + t46; 
complex<T> t81 = -(spa14*spb34) - spa16*spb36 + t46; 
complex<T> t87 = -(spa14*spb16) - spa34*spb36 + t45; 
complex<T> t118 = -(t10*t16); 
complex<T> t119 = t12*T(2); 
complex<T> t147 = spa35*t84; 
complex<T> t154 = t71*(-(spa24*spb25) + t74); 
complex<T> t160 = spa12*t70; 
complex<T> t164 = spb26*t8; 
complex<T> t166 = square(t77); 
complex<T> t171 = spa45*t71; 
complex<T> t184 = spa24*t44; 
complex<T> t185 = t10*t8; 
complex<T> t187 = -(t16*T(2)); 
complex<T> t191 = square(t70); 
complex<T> t192 = square(t71); 
complex<T> t211 = spa26*t83; 
complex<T> t213 = t10*t11; 
complex<T> t229 = t11*t9; 
complex<T> t230 = t12*t16; 
complex<T> t231 = -(t17*T(2)); 
complex<T> t239 = -t270; 
complex<T> t241 = spa35*t189; 
complex<T> t260 = spa45*t11; 
complex<T> t266 = spa13*t78; 
complex<T> t278 = spa25*t89; 
complex<T> t283 = t79*t8; 
complex<T> t284 = spa15*t44; 
complex<T> t285 = -(t9*T(2)); 
complex<T> t302 = spa46*t44; 
complex<T> t303 = spa12*t9; 
complex<T> t318 = t12*t17; 
complex<T> t320 = spa46*t10; 
complex<T> t329 = spa25*t90; 
complex<T> t341 = t236*t97; 
complex<T> t349 = spa25*t98; 
complex<T> t353 = spa56*t17; 
complex<T> t359 = spa13*t16; 
complex<T> t361 = t9*T(2); 
complex<T> t384 = spa46*t9; 
complex<T> t388 = spa15*t10; 
complex<T> t389 = spa12*t16; 
complex<T> t417 = t44*t85; 
complex<T> t420 = spa16*t153; 
complex<T> t422 = t10*t25; 
complex<T> t426 = spa13*t11; 
complex<T> t433 = t11*t153; 
complex<T> t435 = t214*t80; 
complex<T> t438 = spa15*t17; 
complex<T> t441 = t162*t28; 
complex<T> d1 = t12*t82*square(spa26); d1 = T(1)/d1;
complex<T> d3 = t12*t82*t83*cube(spa26); d3 = T(1)/d3;
complex<T> d4 = spa45*t12*t48*t83; d4 = T(1)/d4;
complex<T> d6 = spa56*t236*t82*t83*T(6); d6 = T(1)/d6;
complex<T> d8 = spb12*t12*t82*cube(spa26); d8 = T(1)/d8;
complex<T> d12 = t12*t82*square(s36 + s46 + s56)*square(spa26); d12 = T(1)/d12;
complex<T> d13 = spa45*t12*t13*square(s36 + s46 + s56); d13 = T(1)/d13;
complex<T> d14 = spa26*spa56*t123*t82*cube(t83); d14 = T(1)/d14;
complex<T> d21 = (-s126 + s16)*spa23*spa45*t126*t16*T(2); d21 = T(1)/d21;
complex<T> d26 = t12*t82*cube(spa26); d26 = T(1)/d26;
complex<T> d27 = (s16 - s345)*spa23*spa45*t126*t16*T(6); d27 = T(1)/d27;
complex<T> d29 = (s16 - s345)*t12*t82*cube(spa26); d29 = T(1)/d29;
complex<T> d31 = t12*t82*square(s23 + s24 + s25)*square(spa26); d31 = T(1)/d31;
complex<T> d32 = spa23*spa45*t126*t16*cube(s16 - s345)*T(3); d32 = T(1)/d32;
complex<T> d33 = spa34*t12*t49*t89; d33 = T(1)/d33;
complex<T> d37 = spa34*t12*square(s25 + s35 + s45)*square(spa25); d37 = T(1)/d37;
complex<T> d39 = spa16*t12*t49*t90; d39 = T(1)/d39;
complex<T> d43 = spa16*t12*square(s23 + s24)*square(spa25); d43 = T(1)/d43;
complex<T> d45 = t12*t79*square(spa35); d45 = T(1)/d45;
complex<T> d47 = t12*t79*cube(spa35); d47 = T(1)/d47;
complex<T> d49 = (s34 - s345)*t12*t79*cube(spa35); d49 = T(1)/d49;
complex<T> d52 = t12*t79*square(s35 + s45)*square(spa35); d52 = T(1)/d52;
complex<T> d55 = t12*t79*t84*cube(spa35); d55 = T(1)/d55;
complex<T> d56 = spa12*t12*t48*t84; d56 = T(1)/d56;
complex<T> d59 = t12*t79*square(s34 + s35)*square(spa35); d59 = T(1)/d59;
complex<T> d60 = spa12*t12*t13*square(s34 + s35); d60 = T(1)/d60;
complex<T> d61 = spa23*spa35*t123*t79*cube(t84); d61 = T(1)/d61;
complex<T> d67 = spb45*t12*t79*cube(spa35); d67 = T(1)/d67;
complex<T> d70 = t12*t82*square(square(spa26)); d70 = T(1)/d70;
complex<T> d71 = spa45*t12*t13; d71 = T(1)/d71;
complex<T> d72 = spa45*t12*square(square(spa36)); d72 = T(1)/d72;
complex<T> d76 = spa34*t12*square(spa25); d76 = T(1)/d76;
complex<T> d77 = spa34*t12*square(square(spa25)); d77 = T(1)/d77;
complex<T> d79 = spa16*t12*square(spa25); d79 = T(1)/d79;
complex<T> d80 = spa16*t12*square(square(spa25)); d80 = T(1)/d80;
complex<T> d81 = t12*t79*square(square(spa35)); d81 = T(1)/d81;
complex<T> d82 = spa34*t12*square(spa26); d82 = T(1)/d82;
complex<T> d83 = spa12*t12*t13; d83 = T(1)/d83;
complex<T> d84 = spa34*t12*square(square(spa26)); d84 = T(1)/d84;
complex<T> d85 = spa12*t12*square(square(spa36)); d85 = T(1)/d85;
complex<T> d89 = t12*square(square(spa36)); d89 = T(1)/d89;
complex<T> d92 = t12*square(square(spa25)); d92 = T(1)/d92;
complex<T> t56 = -(t70*t72); 
complex<T> t57 = -t154; 
complex<T> t114 = spb12*t52; 
complex<T> t117 = -t184; 
complex<T> t121 = -t260; 
complex<T> t127 = d82*spb45; 
complex<T> t128 = d70*s12; 
complex<T> t132 = d81*s34; 
complex<T> t144 = d55*t76; 
complex<T> t149 = -t284; 
complex<T> t165 = d45*spa13; 
complex<T> t200 = d84*spb45; 
complex<T> t215 = t72*t9; 
complex<T> t219 = t76*t82; 
complex<T> t225 = d61*t169; 
complex<T> t227 = d14*t190; 
complex<T> t228 = d32*t191; 
complex<T> t244 = d29*t70; 
complex<T> t262 = d1*spa46; 
complex<T> t263 = spa34*t81; 
complex<T> t273 = d3*t73; 
complex<T> t310 = t123*t147; 
complex<T> t311 = d49*t71; 
complex<T> t332 = t185*t187; 
complex<T> t333 = t17*t52; 
complex<T> t336 = d37*spa24; 
complex<T> t338 = d60*spb35; 
complex<T> t362 = d43*spa15; 
complex<T> t366 = d29*spa12; 
complex<T> t369 = spa56*t241; 
complex<T> t413 = d59*t10; 
complex<T> t442 = d3*t164; 
complex<T> t485 = t420*t74; 
complex<T> d2 = t13*t318; d2 = T(1)/d2;
complex<T> d5 = spa56*t341*t82*T(2); d5 = T(1)/d5;
complex<T> d9 = spa34*t341*t353*T(2); d9 = T(1)/d9;
complex<T> d10 = spa34*t236*t353*t83*T(6); d10 = T(1)/d10;
complex<T> d11 = t318*t48*t83; d11 = T(1)/d11;
complex<T> d17 = (s12 - s123)*t318*t48; d17 = T(1)/d17;
complex<T> d18 = t13*t318*square(s34 + s35 + s36); d18 = T(1)/d18;
complex<T> d20 = (s12 - s123)*spa36*t119*t17; d20 = T(1)/d20;
complex<T> d22 = spa34*t119*t349; d22 = T(1)/d22;
complex<T> d23 = t119*t16*t349; d23 = T(1)/d23;
complex<T> d24 = t230*square(spa25); d24 = T(1)/d24;
complex<T> d25 = spa26*spa45*t16*t52; d25 = T(1)/d25;
complex<T> d28 = (s16 - s345)*t230*t49; d28 = T(1)/d28;
complex<T> d30 = t230*square(s23 + s24 + s25)*square(spa25); d30 = T(1)/d30;
complex<T> d34 = spa34*t12*t278*T(6); d34 = T(1)/d34;
complex<T> d35 = t230*t278*T(6); d35 = T(1)/d35;
complex<T> d36 = t230*t49*t89; d36 = T(1)/d36;
complex<T> d38 = spa25*spa34*t52*cube(t89); d38 = T(1)/d38;
complex<T> d40 = spa16*t329*t52; d40 = T(1)/d40;
complex<T> d41 = t16*t329*t52; d41 = T(1)/d41;
complex<T> d42 = t230*t49*t90; d42 = T(1)/d42;
complex<T> d44 = spa16*spa25*t52*cube(t90); d44 = T(1)/d44;
complex<T> d46 = spa35*t389*t52; d46 = T(1)/d46;
complex<T> d48 = (s34 - s345)*t230*t49; d48 = T(1)/d48;
complex<T> d51 = t230*square(s35 + s45)*square(spa25); d51 = T(1)/d51;
complex<T> d58 = t318*t48*t84; d58 = T(1)/d58;
complex<T> d63 = (-s123 + s45)*t318*t48; d63 = T(1)/d63;
complex<T> d64 = t13*t318*square(s46 + s56); d64 = T(1)/d64;
complex<T> d66 = s45*spa35*t52*t79; d66 = T(1)/d66;
complex<T> d69 = (s12 - s123)*t267*t318; d69 = T(1)/d69;
complex<T> d73 = t318*square(square(spa36)); d73 = T(1)/d73;
complex<T> d74 = t13*t353; d74 = T(1)/d74;
complex<T> d75 = t353*square(square(spa36)); d75 = T(1)/d75;
complex<T> d78 = t230*square(square(spa25)); d78 = T(1)/d78;
complex<T> d86 = t119*t82*square(spa26); d86 = T(1)/d86;
complex<T> d87 = t119*t13*t17; d87 = T(1)/d87;
complex<T> d88 = t119*t13; d88 = T(1)/d88;
complex<T> d90 = t119*t16*square(spa25); d90 = T(1)/d90;
complex<T> d91 = t119*square(spa25); d91 = T(1)/d91;
complex<T> d93 = t119*t79*square(spa35); d93 = T(1)/d93;
complex<T> t43 = d9*spa16*t302*t73 + d21*spa24*(spa24*t44*t72 + t160*t91) + d5*spb26*t44*square(spa46); 
complex<T> t125 = d24*s45; 
complex<T> t133 = -(d78*s345); 
complex<T> t138 = -(d22*spb56); 
complex<T> t148 = -(d20*t44*(spa16*spa34*spb36 + t266 - spa34*t76)); 
complex<T> t168 = d24*spa15; 
complex<T> t193 = d73*s123; 
complex<T> t194 = d38*t31; 
complex<T> t216 = d2*t86; 
complex<T> t234 = d48*t75; 
complex<T> t235 = d78*s45; 
complex<T> t246 = d23*t74; 
complex<T> t250 = d63*t87; 
complex<T> t259 = t227*t73; 
complex<T> t293 = d36*t74; 
complex<T> t298 = d76*spa15*spb16*t184 + d77*spb16*t229*t231 + s16*t184*t262 + d70*s16*t283*t285; 
complex<T> t299 = d93*s34*spa13*t254*t44 + s45*t132*t213*t50; 
complex<T> t300 = s12*t184*t262 + t128*t283*t285 + d72*spb12*t332 + d71*spb12*t44*t86; 
complex<T> t313 = d64*t77; 
complex<T> t317 = s34*t165*t284 - t132*t213*t82*T(2); 
complex<T> t326 = d74*spb23; 
complex<T> t335 = d28*spa15; 
complex<T> t337 = d75*spb23; 
complex<T> t350 = s16*(d86*s12*spa46*t184 + t128*t47*t8*t9); 
complex<T> d7 = spa26*t114*t214*t82; d7 = T(1)/d7;
complex<T> d15 = spa26*spa34*t114*t17; d15 = T(1)/d15;
complex<T> d16 = (s12 - s123)*spa36*t333; d16 = T(1)/d16;
complex<T> d19 = spa36*t333*cube(s12 - s123); d19 = T(1)/d19;
complex<T> d50 = (s34 - s345)*t369*t389; d50 = T(1)/d50;
complex<T> d53 = t369*t389*cube(s34 - s345); d53 = T(1)/d53;
complex<T> d54 = spa23*t310*t79; d54 = T(1)/d54;
complex<T> d57 = spa16*spa23*t17*t310; d57 = T(1)/d57;
complex<T> d62 = (-s123 + s45)*spa36*t333; d62 = T(1)/d62;
complex<T> d65 = spa36*t333*cube(-s123 + s45); d65 = T(1)/d65;
complex<T> d68 = s45*spa16*spa35*t333; d68 = T(1)/d68;
complex<T> t4 = -(d39*spa12*spb23*t229) - spa24*t11*t17*t293 - d44*spb23*t21*t303 + d24*t417 + d44*t115*t411*t435 + d38*spa15*t116*t485 + d41*t184*t80 + spb23*t303*t362*t80 - d42*t438*t80*t9 + spb56*(d34*spa45*t44 + t260*t336*t74 + t121*(t194 + d33*t9)) + d40*spa12*spb23*t91 + d35*spa15*t74*t91; 
complex<T> t5 = d21*t160*t184 + t168*t184 + d31*spa46*t160*t215 - t214*t215*t228 + d33*spa45*spb56*t229 - d28*t160*t229 + d32*t160*t22*t252 + spb56*t194*t260 + spa24*t11*t17*t293 - t17*t215*t335 + d22*spa45*spb56*t44 + d27*t215*t44 - d38*spa15*t116*t485 + d30*spa15*t303*t56 + d35*t284*t74 + spb56*t121*t336*t74 + t244*t384*t79 + d26*spa12*spa14*spa24*t8 + t215*t366*t8 + d34*spa45*spb56*t91 + d27*spa24*t160*t91 + d21*t215*t91 + spa15*t246*t91 + spa24*t262*t91 - d25*spa46*cube(spa14); 
complex<T> t6 = t168*t184 + d39*spa12*spb23*t229 + d48*t171*t229 + d51*spa24*t154*t260 + d44*spb23*t21*t303 - d53*t171*t217*t32 + d47*spa14*spa45*t388 - d44*t115*t411*t435 + d40*spa12*spb23*t44 + t311*t426*t50 + d52*spa13*t260*t57 + d48*spa24*t11*t17*(spa24*spb25 - t74) + d53*t192*t433*(spa24*spb25 - t74) + d49*spa45*t213*(-(spa24*spb25) + t74) - spb23*t303*t362*t80 + d42*t438*t80*t9 + spa15*t165*t91 + d50*spa15*t171*t91 + d50*t11*(-(spa24*spb25) + t74)*t91 + d41*spa24*t80*t91 - d46*spa13*cube(spa14); 
complex<T> t7 = d27*t160*t184 + d30*spa15*t160*t215 + t214*t215*t228 + d28*t160*t229 + d14*t162*t164*t24 + t219*t225*t251 - d32*t160*t22*t252 + d52*spa13*t154*t260 + t184*t262 + spa24*t273*t283 + t165*t284 + d50*t171*t284 + d53*t171*t217*t32 + d56*spb35*t16*t320 + t17*t215*t335 + d4*t164*t359 + d59*spb35*t219*t388 + d24*t417 + d61*spb35*t163*t422 + d54*spb35*t10*t44 + t244*t384*t47 + d55*spb35*t213*t50 + d31*spa46*t303*t56 + d51*spa24*t260*t57 - d11*spa16*t185*t73 - d13*t164*t305*t73 + d49*t10*t121*(-(spa24*spb25) + t74) + d48*spa24*t11*t17*(-(spa24*spb25) + t74) + d53*t192*t433*(-(spa24*spb25) + t74) + d50*t11*t44*(-(spa24*spb25) + t74) - d58*spa34*t185*t76 - spa34*t320*t338*t76 + d57*t372*t44*t76 + t253*t259*t79 + d12*spa24*t164*t73*t79 - t215*t366*t8 + t144*t388*t82 + t311*t426*t82 + t442*t47*t9 + d48*t121*t71*t9 + d6*t164*t91 + d27*t215*t91 + t216*t91 + d10*spa16*spa46*t73*t91; 
complex<T> t69 = d85*spb45*t332 + t200*t283*t361 + spa24*spa46*t127*t91 + t165*t254*square(spa14) + t125*t85*square(spa14) + d83*spb45*t86*square(spa14) + t17*t229*t235*T(2) - d81*s45*t213*t82*T(2) + s345*(d78*t229*t231 + d73*t332 + d70*t283*t361 + spa15*t165*t91 + spa24*t262*t91 + (spa24*t168 + t216)*square(spa14) + d81*t213*t82*T(2)); 
complex<T> t204 = d65*t166; 
complex<T> t279 = d79*spa15*spb34*t184 + s156*t168*t184 + d78*s156*t229*t231 + d80*spb34*t229*t231; 
complex<T> t280 = t246*t284 + spa45*t138*t44; 
complex<T> t281 = -(t326*t44*t86) + t16*t185*t337*T(2); 
complex<T> t290 = d19*t29; 
complex<T> t316 = spa15*t125*t184 + spa46*t127*t184 + t229*t231*t235 + t200*t283*t285; 
complex<T> t330 = d90*s234*s345*spa15*t184 + s234*t133*t17*t229 + d92*t137*t17*t229 + d91*t137*t417; 
complex<T> t331 = d89*t16*t185*t270 + s126*t118*t193*t8 + d87*s123*s126*t44*t86 + d88*t239*t44*t86; 
complex<T> t342 = t193*t332 + t332*t337 + s123*t216*t44 + t326*t44*t86; 
complex<T> t1 = -(d14*t162*t164*t24) + d19*t163*t266*t27 - d4*t164*t359 + d6*t164*t44 + t216*t44 + d16*t266*t44 + t253*t259*t47 + d8*spa16*t384*(-(spa46*spb26) + t70) + d11*spa16*t185*t73 + d10*spa16*t302*t73 + d13*t164*t305*t73 + d12*t283*t45*t73 + d17*spa34*t185*(spa16*spb36 - t76) + spa34*t10*t290*(spa16*spb36 - t76) + d17*spa46*t118*t78 + d18*spa34*t320*(spa16*spb36 - t76)*t78 + spa24*t273*t47*t8 + d7*t44*(-(spa16*spb26) + t72)*t9 + t442*t79*t9 + d8*(-(spa16*spb26) + t72)*t8*t9 + d5*t164*t91 + spa24*t262*t91 + d15*spa24*(-(spa46*spb26) + t70)*t91 + d9*spa16*spa46*t73*t91 + d16*spa34*(-(spa16*spb36) + t76)*t91; 
complex<T> t2 = -(d19*t163*t266*t27) + d69*t266*t44 + d63*spa16*t185*(-(spa34*spb36) + t73) + d62*spa16*t44*(-(spa34*spb36) + t73) + d17*spa34*t185*(-(spa16*spb36) + t76) + spa34*t10*t290*(-(spa16*spb36) + t76) + d62*t302*t77 + d65*spa46*t441*t77 + d17*t16*t320*t78 + d18*spa34*t320*(-(spa16*spb36) + t76)*t78 + t305*t313*(spa34*spb36 - t73)*t8 + spa16*t204*(-(spa34*spb36) + t73)*t8 - d63*t359*t77*t8 + t216*t91 + d69*spa34*(-(spa16*spb36) + t76)*t91; 
complex<T> t3 = d56*spa46*spb35*t118 - d61*spb35*t163*t422 + t216*t44 + t219*t413*t46 + t144*t388*t50 - d67*t11*t372*(spa13*spb35 + t71) + d63*spa16*t185*(spa34*spb36 - t73) + d67*t213*(-(spa24*spb25) - spa34*spb35 + t74) + d58*spa34*t185*t76 + spa34*t320*t338*t76 + t225*t251*t50*t76 - d65*spa46*t441*t77 + spa16*t204*(spa34*spb36 - t73)*t8 + t305*t313*(-(spa34*spb36) + t73)*t8 + d63*t359*t77*t8 + d55*spb35*t213*t82 + d54*spb35*t10*t91 + spa15*t165*t91 + d66*spa15*(spa13*spb35 + t71)*t91 + d62*spa16*(-(spa34*spb36) + t73)*t91 + d68*t11*(-(spa24*spb25) - spa34*spb35 + t74)*t91 + d57*t372*t76*t91 + d62*spa46*t77*t91; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t148*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t280*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + t6*(*CI_users[6]->get_value(mc,ind,mu)) + t7*(*CI_users[7]->get_value(mc,ind,mu)) + t3*(*CI_users[8]->get_value(mc,ind,mu)) + t2*(*CI_users[9]->get_value(mc,ind,mu)) + t300*(*CI_users[10]->get_value(mc,ind,mu)) + t342*(*CI_users[11]->get_value(mc,ind,mu)) + t298*(*CI_users[12]->get_value(mc,ind,mu)) + t281*(*CI_users[13]->get_value(mc,ind,mu)) + t279*(*CI_users[14]->get_value(mc,ind,mu)) + t317*(*CI_users[15]->get_value(mc,ind,mu)) + t69*(*CI_users[16]->get_value(mc,ind,mu)) + t316*(*CI_users[17]->get_value(mc,ind,mu)) + t350*(*CI_users[18]->get_value(mc,ind,mu)) + t331*(*CI_users[19]->get_value(mc,ind,mu)) + t330*(*CI_users[20]->get_value(mc,ind,mu)) + t299*(*CI_users[21]->get_value(mc,ind,mu)));  
 return(result);
} 
  


//cc


C6g_ppmmmm_G_wCI::\
C6g_ppmmmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c2, c34, c5, c16));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c126));
CI_users.push_back(new Cached_Box_Integral_User(c3, c45, c6, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c6, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c6, c1, c234));
CI_users.push_back(new Cached_Box_Integral_User(c6, c1, c2, c345));
} 
  
  
template <class T> SeriesC<T> 
     C6g_ppmmmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C6g :  ppmmmm G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 17:38:58 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb21 = SPB(2,1);
complex<T> spb32 = SPB(3,2);
complex<T> spb43 = SPB(4,3);
complex<T> spb54 = SPB(5,4);
complex<T> spb65 = SPB(6,5);
complex<T> spa61 = SPA(6,1);
complex<T> spb61 = SPB(6,1);
complex<T> spa32 = SPA(3,2);
complex<T> spa43 = SPA(4,3);
complex<T> spa54 = SPA(5,4);
complex<T> spa65 = SPA(6,5);
complex<T> spb51 = SPB(5,1);
complex<T> spb52 = SPB(5,2);
complex<T> spa51 = SPA(5,1);
complex<T> spb62 = SPB(6,2);
complex<T> spb41 = SPB(4,1);
complex<T> spb42 = SPB(4,2);
complex<T> spa64 = SPA(6,4);
complex<T> spa41 = SPA(4,1);
complex<T> s12 = S(1,2);
complex<T> s16 = -(spb61*spa61);
complex<T> s156 = SS(1,5,6);
complex<T> s23 = -(spb32*spa32);
complex<T> s234 = SS(2,3,4);
complex<T> s123 = SS(1,2,3);
complex<T> s345 = SS(3,4,5);
complex<T> t5 = spb61*spb43; 
complex<T> t6 = spb32*spb65; 
complex<T> t9 = square(spb21*spa51 + spb62*spa65); 
complex<T> t12 = square(spb51); 
complex<T> t14 = square(spb21); 
complex<T> t16 = -(spb42*(spb51*spa54 + spb61*spa64)); 
complex<T> t22 = cube(spb21); 
complex<T> t23 = -(spb21*spa51) - spb62*spa65; 
complex<T> t24 = square(spb52); 
complex<T> t25 = square(spa65); 
complex<T> t28 = spb51*spa54 + spb61*spa64; 
complex<T> t29 = spb21*spa41 + spb52*spa54 + spb62*spa64; 
complex<T> t32 = square(spb41); 
complex<T> t33 = square(spb42); 
complex<T> t36 = spb52*spa65; 
complex<T> t41 = -(spb62*spa64); 
complex<T> d10 = spb61*spb54*spb65*T(2); d10 = T(1)/d10;
complex<T> d15 = spb32*spb43*spb54*T(2); d15 = T(1)/d15;
complex<T> t15 = spb54*t6; 
complex<T> t17 = spb41*(-(spb21*spa41) - spb52*spa54 + t41); 
complex<T> t38 = spb51*t23; 
complex<T> t39 = square(t28); 
complex<T> t40 = square(t29); 
complex<T> t55 = s12*t22; 
complex<T> t71 = spa43*t22; 
complex<T> t77 = spa65*t22; 
complex<T> d7 = spb54*spb65*t5*T(2); d7 = T(1)/d7;
complex<T> d9 = spb54*t5*T(2); d9 = T(1)/d9;
complex<T> d12 = spb61*t6*T(2); d12 = T(1)/d12;
complex<T> d13 = t5*t6*T(2); d13 = T(1)/d13;
complex<T> d14 = spb32*t5*T(2); d14 = T(1)/d14;
complex<T> d1 = (-s156 + s16)*spb43*t15*T(6); d1 = T(1)/d1;
complex<T> d2 = (-s156 + s16)*t15*t5*T(6); d2 = T(1)/d2;
complex<T> d3 = spb43*t15*cube(-s156 + s16)*T(3); d3 = T(1)/d3;
complex<T> d4 = t15*t5*T(6); d4 = T(1)/d4;
complex<T> d5 = (s23 - s234)*t15*t5*T(6); d5 = T(1)/d5;
complex<T> d6 = t15*t5*cube(s23 - s234)*T(3); d6 = T(1)/d6;
complex<T> d8 = t15*t5*T(2); d8 = T(1)/d8;
complex<T> d11 = t15*T(2); d11 = T(1)/d11;
complex<T> d16 = spb43*t15*T(2); d16 = T(1)/d16;
complex<T> t30 = d4*T(11); 
complex<T> t31 = d8*s123; 
complex<T> t42 = d6*t32; 
complex<T> t45 = d8*s234*s345*t22 - d11*spa61*t71; 
complex<T> t52 = d3*t24; 
complex<T> t61 = d6*t33; 
complex<T> t62 = d3*t9; 
complex<T> t1 = -(t16*t40*t42) - t17*t39*t61 - d5*t14*(t16 + t17)*T(11); 
complex<T> t2 = t22*t30 + t16*t40*t42 + t17*t39*t61 + d5*t14*(t16 + t17)*T(11); 
complex<T> t54 = s234*t22*t31 - d9*spa32*t77; 
complex<T> t63 = s345*t22*t31 + d13*spa54*t55; 
complex<T> t76 = spb61*t52; 
complex<T> t3 = t22*t30 - t12*t36*t62 + t25*t38*t76 - d1*t14*t36*T(11) + d2*t14*t38*T(11); 
complex<T> t4 = t12*t36*t62 - t25*t38*t76 + d1*t14*t36*T(11) - d2*t14*t38*T(11); 
complex<T> co1 = -(d7*spa32*t55); 
complex<T> co2 = d10*spa32*t71; 
complex<T> co3 = d12*spa54*t71; 
complex<T> co4 = d14*spa54*t77; 
complex<T> co5 = d15*spa61*t77; 
complex<T> co6 = -(d16*spa61*t55); 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + t54*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t45*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t63*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)) + co5*(*CI_users[11]->get_value(mc,ind,mu)) + co6*(*CI_users[12]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_ppmmmm_nf_wCI::\
C6g_ppmmmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
} 
  
  
template <class T> SeriesC<T> 
     C6g_ppmmmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C6g :  ppmmmm nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 17:39:02 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb21 = SPB(2,1);
complex<T> spb61 = SPB(6,1);
complex<T> spb32 = SPB(3,2);
complex<T> spb43 = SPB(4,3);
complex<T> spb54 = SPB(5,4);
complex<T> spb65 = SPB(6,5);
complex<T> spb51 = SPB(5,1);
complex<T> spb52 = SPB(5,2);
complex<T> spa65 = SPA(6,5);
complex<T> spa51 = SPA(5,1);
complex<T> spb62 = SPB(6,2);
complex<T> spb41 = SPB(4,1);
complex<T> spb42 = SPB(4,2);
complex<T> spa54 = SPA(5,4);
complex<T> spa64 = SPA(6,4);
complex<T> spa41 = SPA(4,1);
complex<T> s16 = S(1,6);
complex<T> s156 = SS(1,5,6);
complex<T> s23 = S(2,3);
complex<T> s234 = SS(2,3,4);
complex<T> t7 = square(spb21*spa51 + spb62*spa65); 
complex<T> t10 = square(spb51); 
complex<T> t12 = square(spb21); 
complex<T> t13 = square(spb51*spa54 + spb61*spa64); 
complex<T> t14 = square(spb21*spa41 + spb52*spa54 + spb62*spa64); 
complex<T> t16 = spb32*T(3); 
complex<T> t17 = spb43*spb54; 
complex<T> t21 = square(spb41); 
complex<T> t22 = square(spb42); 
complex<T> t23 = -(spb51*spa54) - spb61*spa64; 
complex<T> t24 = -(spb21*spa41) - spb52*spa54 - spb62*spa64; 
complex<T> t25 = -(spb21*spa51) - spb62*spa65; 
complex<T> t26 = square(spb52); 
complex<T> t27 = square(spa65); 
complex<T> t30 = cube(spb21); 
complex<T> t47 = t14*t21; 
complex<T> t48 = t26*t27; 
complex<T> d1 = (-s156 + s16)*spb65*t16*t17; d1 = T(1)/d1;
complex<T> d2 = (-s156 + s16)*spb61*spb65*t16*t17; d2 = T(1)/d2;
complex<T> d3 = spb65*t16*t17*cube(-s156 + s16); d3 = T(1)/d3;
complex<T> d4 = spb61*spb65*t16*t17; d4 = T(1)/d4;
complex<T> d5 = (s23 - s234)*spb61*spb65*t16*t17; d5 = T(1)/d5;
complex<T> d6 = spb61*spb65*t16*t17*cube(s23 - s234); d6 = T(1)/d6;
complex<T> t42 = d3*spb61; 
complex<T> t46 = d6*t13; 
complex<T> t54 = d3*t7; 
complex<T> t1 = -(d5*spb42*t12*t23) - d5*spb41*t12*t24 - d4*t30 - spb41*t22*t24*t46 - d6*spb42*t23*t47; 
complex<T> t2 = d5*t12*(spb42*t23 + spb41*t24) + spb41*t22*t24*t46 + d6*spb42*t23*t47; 
complex<T> t3 = -(d1*spb52*spa65*t12) + d2*spb51*t12*t25 + spb51*t25*t42*t48 - spb52*spa65*t10*t54; 
complex<T> t4 = d1*spb52*spa65*t12 - d2*spb51*t12*t25 - d4*t30 - spb51*t25*t42*t48 + spb52*spa65*t10*t54; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_pmpmmm_G_wCI::\
C6g_pmpmmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c126, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c6, c2345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c1236));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c6, c345));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c2, c34, c5, c16));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c126));
CI_users.push_back(new Cached_Box_Integral_User(c3, c45, c6, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c6, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c34, c2, c16));
CI_users.push_back(new Cached_Box_Integral_User(c5, c6, c1, c234));
CI_users.push_back(new Cached_Box_Integral_User(c6, c1, c2, c345));
} 
  
  
template <class T> SeriesC<T> 
     C6g_pmpmmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C6g :  pmpmmm G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 17:47:44 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb31 = SPB(3,1);
complex<T> spb41 = SPB(4,1);
complex<T> spb61 = SPB(6,1);
complex<T> spb42 = SPB(4,2);
complex<T> spb54 = SPB(5,4);
complex<T> spb65 = SPB(6,5);
complex<T> spb21 = SPB(2,1);
complex<T> spb32 = SPB(3,2);
complex<T> spb43 = SPB(4,3);
complex<T> spb62 = SPB(6,2);
complex<T> spb63 = SPB(6,3);
complex<T> spa21 = SPA(2,1);
complex<T> spa61 = SPA(6,1);
complex<T> spb51 = SPB(5,1);
complex<T> spb52 = SPB(5,2);
complex<T> spb53 = SPB(5,3);
complex<T> spa32 = SPA(3,2);
complex<T> spa43 = SPA(4,3);
complex<T> spa54 = SPA(5,4);
complex<T> spa65 = SPA(6,5);
complex<T> spa62 = SPA(6,2);
complex<T> spa42 = SPA(4,2);
complex<T> spa52 = SPA(5,2);
complex<T> spa51 = SPA(5,1);
complex<T> spa64 = SPA(6,4);
complex<T> spa41 = SPA(4,1);
complex<T> s23 = -(spb32*spa32);
complex<T> s34 = -(spb43*spa43);
complex<T> s12 = -(spb21*spa21);
complex<T> s16 = -(spb61*spa61);
complex<T> s126 = SS(1,2,6);
complex<T> s156 = SS(1,5,6);
complex<T> s24 = -(spb42*spa42);
complex<T> s234 = SS(2,3,4);
complex<T> s25 = -(spb52*spa52);
complex<T> s35 = S(3,5);
complex<T> s45 = -(spb54*spa54);
complex<T> s123 = SS(1,2,3);
complex<T> s36 = S(3,6);
complex<T> s46 = S(4,6);
complex<T> s56 = -(spb65*spa65);
complex<T> s345 = SS(3,4,5);
complex<T> t7 = square(spb53); 
complex<T> t8 = square(spb51); 
complex<T> t10 = spb54*spb65; 
complex<T> t11 = square(spb41); 
complex<T> t13 = spb21*spb61; 
complex<T> t18 = square(spb31*spa32 + spb41*spa42); 
complex<T> t19 = square(spb43*spa42 + spb53*spa52); 
complex<T> t20 = square(spb31*spa61 + spb32*spa62); 
complex<T> t21 = square(spb31*spa41 + spb53*spa54 + spb63*spa64); 
complex<T> t36 = square(spb31); 
complex<T> t37 = square(spb63); 
complex<T> t38 = -(spb41*spa42); 
complex<T> t39 = -(spb32*T(2)); 
complex<T> t41 = cube(spb52); 
complex<T> t46 = (-s234 + s34)*spb61; 
complex<T> t47 = square(s23 + s24); 
complex<T> t49 = square(spb31*spa32 + spb41*spa42 + spb51*spa52); 
complex<T> t51 = square(spb51*spa54 + spb61*spa64); 
complex<T> t57 = cube(spb31); 
complex<T> t58 = -(spb43*spa42) - spb53*spa52; 
complex<T> t60 = -(spb31*spa61) - spb32*spa62; 
complex<T> t61 = -(spb31*spa51) - spb63*spa65; 
complex<T> t62 = spb21*spa52 - spb61*spa65; 
complex<T> t64 = -(spb31*spa41) - spb53*spa54 - spb63*spa64; 
complex<T> t65 = -(spb51*spa54) - spb61*spa64; 
complex<T> t66 = -(spb31*spa51) - spb32*spa52 - spb63*spa65; 
complex<T> t68 = s16 - s234; 
complex<T> t70 = -(spb21*spa42); 
complex<T> t71 = s12 - s345; 
complex<T> t73 = square(spb32); 
complex<T> t76 = spb63*spa62; 
complex<T> t78 = -s156 + s16; 
complex<T> t80 = cube(spb51); 
complex<T> t83 = -(spb31*spa41) - spb32*spa42 - spb53*spa54 - spb63*spa64; 
complex<T> t86 = s234*s345; 
complex<T> t97 = square(spb43); 
complex<T> t98 = cube(spb53); 
complex<T> t99 = square(spa65); 
complex<T> t102 = spb51*spb53; 
complex<T> t108 = square(spb61); 
complex<T> t109 = square(spa62); 
complex<T> t111 = spb62*spb54; 
complex<T> t112 = spb42*spb65; 
complex<T> t118 = spa61*spa43; 
complex<T> t127 = spb21*spb32; 
complex<T> t134 = square(square(spb31)); 
complex<T> t137 = square(spb21); 
complex<T> t139 = square(spa42); 
complex<T> t148 = spb31*spa51 + spb63*spa65; 
complex<T> t155 = spb32*T(2); 
complex<T> t160 = spb43*spb65; 
complex<T> t190 = s34*spb42; 
complex<T> t235 = spb41*spb43; 
complex<T> t287 = spb51*spb61; 
complex<T> t16 = -t76; 
complex<T> t17 = -(spb43*t65); 
complex<T> t23 = square(t148); 
complex<T> t24 = square(spb32*spa52 + t148); 
complex<T> t40 = spb43*t10; 
complex<T> t42 = -t102; 
complex<T> t59 = -(spb31*spa32) + t38; 
complex<T> t63 = -(spb31*spa32) - spb51*spa52 + t38; 
complex<T> t67 = spb61*t10; 
complex<T> t87 = -t118; 
complex<T> t100 = -(t36*T(2)); 
complex<T> t103 = spb43*t11; 
complex<T> t128 = t36*T(2); 
complex<T> t131 = spb63*t58; 
complex<T> t153 = square(t62); 
complex<T> t154 = t7*t8; 
complex<T> t161 = t62*t66; 
complex<T> t182 = -(t36*T(4)); 
complex<T> t183 = t13*t37; 
complex<T> t200 = spb52*t68; 
complex<T> t204 = t127*T(2); 
complex<T> t209 = spb41*t64; 
complex<T> t223 = spb51*t7; 
complex<T> t227 = spb52*t112; 
complex<T> t245 = spb61*t60; 
complex<T> t246 = spb53*t8; 
complex<T> t280 = spa21*t134; 
complex<T> t292 = spa65*t134; 
complex<T> t368 = t287*t61; 
complex<T> d12 = (-s126 + s16)*spb61*spb52*spb43*t111; d12 = T(1)/d12;
complex<T> d18 = (s16 - s345)*spb61*spb52*spb43*t111*T(6); d18 = T(1)/d18;
complex<T> d23 = spb61*spb52*spb43*t111*cube(s16 - s345)*T(3); d23 = T(1)/d23;
complex<T> d38 = t10*t46*cube(spb42); d38 = T(1)/d38;
complex<T> d39 = t10*t41*t46; d39 = T(1)/d39;
complex<T> d45 = t10*t13*t190*T(6); d45 = T(1)/d45;
complex<T> d59 = spb61*t160*square(spb52); d59 = T(1)/d59;
complex<T> d60 = spb61*t160*square(square(spb52)); d60 = T(1)/d60;
complex<T> d61 = t160*square(spb62); d61 = T(1)/d61;
complex<T> d62 = t160*square(square(spb62)); d62 = T(1)/d62;
complex<T> d65 = spb43*spb54*t13*T(2); d65 = T(1)/d65;
complex<T> d66 = t10*t13*T(2); d66 = T(1)/d66;
complex<T> d68 = spb65*t13*t155; d68 = T(1)/d68;
complex<T> d69 = spb61*t155*t160; d69 = T(1)/d69;
complex<T> d70 = spb43*t13*t155; d70 = T(1)/d70;
complex<T> d71 = t10*square(spb52); d71 = T(1)/d71;
complex<T> d72 = t10*square(square(spb52)); d72 = T(1)/d72;
complex<T> t14 = spb32*(spb51*spa52 - t59); 
complex<T> t106 = spa42*t59; 
complex<T> t119 = d59*spa54; 
complex<T> t123 = d23*t49; 
complex<T> t124 = d12*t58; 
complex<T> t129 = -t154; 
complex<T> t138 = spb62*t40; 
complex<T> t165 = d60*spa54; 
complex<T> t170 = d23*t19; 
complex<T> t186 = spb32*t40; 
complex<T> t206 = spb61*t40; 
complex<T> t210 = d61*spa54; 
complex<T> t215 = d38*t59; 
complex<T> t226 = d62*spa54; 
complex<T> t234 = t154*t204; 
complex<T> t290 = spb53*t24; 
complex<T> t299 = t36*t42; 
complex<T> d1 = t40*square(spb62); d1 = T(1)/d1;
complex<T> d4 = t40*t71*cube(spb62); d4 = T(1)/d4;
complex<T> d6 = spa21*t40*cube(spb62); d6 = T(1)/d6;
complex<T> d9 = t40*square(s36 + s46 + s56)*square(spb62); d9 = T(1)/d9;
complex<T> d17 = t40*cube(spb62); d17 = T(1)/d17;
complex<T> d20 = (s16 - s345)*t40*cube(spb62); d20 = T(1)/d20;
complex<T> d22 = t40*square(s23 + s24 + s25)*square(spb62); d22 = T(1)/d22;
complex<T> d24 = t40*t41*t68; d24 = T(1)/d24;
complex<T> d28 = t40*square(s25 + s35 + s45)*square(spb52); d28 = T(1)/d28;
complex<T> d30 = t67*square(spb42); d30 = T(1)/d30;
complex<T> d31 = spb32*spb42*t67*T(6); d31 = T(1)/d31;
complex<T> d32 = t67*cube(spb42); d32 = T(1)/d32;
complex<T> d33 = (s23 - s234)*t67*cube(spb42); d33 = T(1)/d33;
complex<T> d34 = (s23 - s234)*spb32*spb42*t67*T(6); d34 = T(1)/d34;
complex<T> d35 = t67*square(s24 + s34)*square(spb42); d35 = T(1)/d35;
complex<T> d36 = spb32*spb42*t67*cube(s23 - s234)*T(3); d36 = T(1)/d36;
complex<T> d37 = t227*t46*T(6); d37 = T(1)/d37;
complex<T> d40 = spb43*t227*t46*T(6); d40 = T(1)/d40;
complex<T> d41 = t40*t41*t46; d41 = T(1)/d41;
complex<T> d42 = t47*t67*square(spb42); d42 = T(1)/d42;
complex<T> d43 = t47*t67*square(spb52); d43 = T(1)/d43;
complex<T> d44 = spb61*t227*cube(-s234 + s34)*T(3); d44 = T(1)/d44;
complex<T> d46 = spa43*t67*cube(spb42); d46 = T(1)/d46;
complex<T> d47 = t13*t190*t40*T(6); d47 = T(1)/d47;
complex<T> d48 = (s34 - s345)*spb52*t13*t40*T(6); d48 = T(1)/d48;
complex<T> d51 = spb52*t13*t40*cube(s34 - s345)*T(3); d51 = T(1)/d51;
complex<T> d52 = t40*square(square(spb62)); d52 = T(1)/d52;
complex<T> d53 = t40*square(spb52); d53 = T(1)/d53;
complex<T> d54 = t40*square(square(spb52)); d54 = T(1)/d54;
complex<T> d55 = t67*square(square(spb42)); d55 = T(1)/d55;
complex<T> d57 = t67*square(spb52); d57 = T(1)/d57;
complex<T> d58 = t67*square(square(spb52)); d58 = T(1)/d58;
complex<T> d64 = t13*t155*t40; d64 = T(1)/d64;
complex<T> d67 = t10*t204; d67 = T(1)/d67;
complex<T> d73 = spb43*spb54*t204; d73 = T(1)/d73;
complex<T> d74 = t155*t40; d74 = T(1)/d74;
complex<T> t56 = d48*(-(spb32*spa52) + t61); 
complex<T> t88 = d64*s123; 
complex<T> t107 = -(d34*T(11)); 
complex<T> t113 = d52*s12; 
complex<T> t114 = d30*s23; 
complex<T> t117 = d55*s34; 
complex<T> t126 = d48*t62; 
complex<T> t162 = d1*s12; 
complex<T> t163 = d37*spa42; 
complex<T> t169 = d44*t18; 
complex<T> t181 = t134*(d64*t86 + d67*t87); 
complex<T> t184 = d30*spb41; 
complex<T> t187 = d34*T(11); 
complex<T> t188 = d4*spa62; 
complex<T> t192 = d43*t106; 
complex<T> t193 = d51*t153; 
complex<T> t207 = d1*spb63; 
complex<T> t212 = d44*t139; 
complex<T> t228 = d41*t59; 
complex<T> t233 = spb63*t182*t210 + t155*t183*t226 + spb21*t154*t165*t39 + t102*t119*t36*T(4); 
complex<T> d2 = (s12 - s126)*t138; d2 = T(1)/d2;
complex<T> d3 = t138*t71*T(6); d3 = T(1)/d3;
complex<T> d5 = spa21*t137*t138*T(6); d5 = T(1)/d5;
complex<T> d7 = (s12 - s126)*spb21*t138; d7 = T(1)/d7;
complex<T> d8 = spb21*t138*t71*T(6); d8 = T(1)/d8;
complex<T> d10 = t138*cube(t71)*T(3); d10 = T(1)/d10;
complex<T> d11 = spb21*spa21*t138*T(6); d11 = T(1)/d11;
complex<T> d13 = spb52*t186*t78; d13 = T(1)/d13;
complex<T> d14 = spb61*spb52*t186*t78; d14 = T(1)/d14;
complex<T> d15 = t206*square(spb52); d15 = T(1)/d15;
complex<T> d16 = spb61*spb32*t138*T(6); d16 = T(1)/d16;
complex<T> d19 = (s16 - s345)*t206*t41; d19 = T(1)/d19;
complex<T> d21 = t206*square(s23 + s24 + s25)*square(spb52); d21 = T(1)/d21;
complex<T> d25 = t186*t200*T(6); d25 = T(1)/d25;
complex<T> d26 = spb61*t186*t200*T(6); d26 = T(1)/d26;
complex<T> d27 = t206*t41*t68; d27 = T(1)/d27;
complex<T> d29 = spb52*t186*cube(t68)*T(3); d29 = T(1)/d29;
complex<T> d49 = (s34 - s345)*t206*t41; d49 = T(1)/d49;
complex<T> d50 = t206*square(s35 + s45)*square(spb52); d50 = T(1)/d50;
complex<T> d56 = t206*square(square(spb52)); d56 = T(1)/d56;
complex<T> d63 = t206*T(2); d63 = T(1)/d63;
complex<T> t2 = -(d32*spb21*spb31*t235) + t184*t36 + t187*t209*t36 - d33*spb21*t103*t64 + d33*spb32*t103*t65 + d36*t103*t21*t65 + spb43*t187*t36*t65 + d35*t103*t64*t65 + d36*t209*t51*t97 + d31*t57*T(11); 
complex<T> t35 = spb21*t124*t128 + d7*t100*t245 + d12*t36*t39*(-(spb51*spa52) + t59) + d2*t100*t76; 
complex<T> t96 = s34*spb41*t100*t114 + s23*t103*t117*t127; 
complex<T> t120 = d13*spa65; 
complex<T> t149 = spb41*t114*t182 + d55*s23*t103*t204; 
complex<T> t150 = s16*(spb63*t100*t162 + spb32*t113*t183); 
complex<T> t171 = d10*t20; 
complex<T> t180 = d53*spa61*t102*t182 + d52*s16*t155*t183 + s16*t182*t207 + d54*spa61*t234; 
complex<T> t195 = d21*t58; 
complex<T> t196 = d10*t60; 
complex<T> t197 = d14*t61; 
complex<T> t202 = spb63*t162*t182 + t113*t155*t183; 
complex<T> t203 = -(d65*spa32*t292) + s234*t134*t88; 
complex<T> t208 = d15*t102; 
complex<T> t219 = -(d69*spa54*t280) + s345*t134*t88; 
complex<T> t334 = spb43*t212; 
complex<T> t1 = -(spb32*t183*t188) + t109*t183*t196 + d7*t128*t245 + d3*t16*t36 + t207*t36 - d8*t245*t36 - d6*spb61*spb32*spb63*(t16 + t58) + d6*spb32*t37*(spb51*spa52 + spb61*spa62 - t59) + d4*t183*t60 + d9*spa62*t183*t60 + d2*t128*t76 + t108*t171*t76 + d11*t36*(t16 + t58)*T(11) + d5*spb32*t36*(spb51*spa52 + spb61*spa62 - t59)*T(11); 
complex<T> t3 = spb21*t100*t124 + d24*spb32*spa65*t129 - d20*spb32*t13*t131 + d22*spb21*t131*t14 + t102*t128*t197 + d29*spa65*t129*t23 + d15*t299 + d18*t14*t36 + t207*t36 - d17*spb21*spb31*t37 + d20*spb21*t14*t37 + d19*t127*t246*t58 + d18*spb21*t36*t58 + spb32*t137*t170*(-(spb51*spa52) + t59) + t102*t127*t195*(-(spb51*spa52) + t59) + d19*t127*t223*(-(spb51*spa52) + t59) + d12*t155*t36*(-(spb51*spa52) + t59) + d27*spb21*t129*t61 + d28*spa65*t154*t61 + d26*t299*t61 + t100*t120*t7 + d25*spa65*t36*t7 - spb21*t123*t58*t73 + d29*t368*t98*t99 + d16*spb63*t57*T(11); 
complex<T> t4 = d38*spa42*t103*t127 + d24*spb32*spa65*t154 + d33*spb32*t11*t17 + t102*t127*t192 + d36*t11*t17*t21 + t127*t223*t228 + d29*spa65*t154*t23 - t127*t215*t235 - t184*t36 + t208*t36 + t107*t209*t36 - spb32*t137*t334*t59 + d42*spb43*t127*t38*t59 + d28*spa65*t129*t61 + d27*spb21*t154*t61 + d26*t102*t36*t61 + d33*spb21*t103*t64 + d35*t11*t17*t64 + spb43*t107*t36*t65 - d25*spa65*t36*t7 + d39*spb32*t246*t70 + t169*t70*t73 - d36*t209*t51*t97 - d29*t368*t98*t99 - spb21*t163*t36*T(11) - d40*spb32*t36*t59*T(11); 
complex<T> t5 = d20*spb32*t13*t131 + t137*t14*t170 + t108*t16*t171 + spb32*t183*t188 + spb21*t102*t14*t195 - t109*t183*t196 + d19*spb21*t14*t223 - t207*t36 + t208*t36 + d8*t245*t36 - d19*t127*t246*t58 - d18*spb21*t36*t58 + d22*t127*t131*(-(spb51*spa52) + t59) + d18*spb32*t36*(-(spb51*spa52) + t59) + d20*t127*t37*(-(spb51*spa52) + t59) - d4*t183*t60 - d9*spa62*t183*t60 + d49*spb21*t129*(-(spb32*spa52) + t61) + t129*t193*(-(spb32*spa52) + t61) + d49*spb32*t154*t62 + d50*t154*(-(spb32*spa52) + t61)*t62 + spb21*t123*t58*t73 + d3*t36*t76 - d51*t290*t62*t80 - t102*t126*t36*T(11) - t36*t56*t8*T(11); 
complex<T> t6 = -(t127*t223*t228) + d42*t106*t127*t235 + t127*t215*t235 + d39*spa42*t127*t246 + d15*t299 + t184*t36 + t127*t192*t42 + spb32*t137*t334*t59 + d49*spb21*t154*(-(spb32*spa52) + t61) + t154*t193*(-(spb32*spa52) + t61) + d49*spb32*t129*t62 + d50*t129*(-(spb32*spa52) + t61)*t62 + d46*spb21*t11*(spb32*spa42 - t64) + d46*spb32*t11*(spb21*spa42 + t65) + d38*spb32*t103*t70 + spb21*spa42*t169*t73 + d51*t290*t62*t80 + t102*t126*t36*T(11) + spb21*t163*t36*T(11) + d40*spb32*t36*t59*T(11) + d47*t11*t36*(-(spb32*spa42) + t64)*T(11) + d45*spb41*t36*(spb21*spa42 + t65)*T(11) + t36*t56*t8*T(11); 
complex<T> t151 = d71*t102*t118*t128 + d56*t127*t154*t86 + t100*t208*t86 + d72*t127*t154*t87; 
complex<T> t218 = d57*spa43*t102*t182 + s34*t182*t184 + t103*t117*t204 + s156*t182*t208 + d56*s156*t234 + d58*spa43*t234 + d55*s156*spb21*t103*t39 + s156*t184*t36*T(4); 
complex<T> t241 = t102*t119*t182 + t165*t234 + t183*t226*t39 + spb63*t210*t36*T(4) + s345*(t182*t208 + d56*t234 + d52*t183*t39 + t207*t36*T(4)); 
complex<T> t251 = t100*t102*t197 + t120*t128*t7; 
complex<T> co1 = d63*spa32*t280; 
complex<T> co2 = d66*spa32*spa43*t134; 
complex<T> co3 = d68*spa43*spa54*t134; 
complex<T> co4 = d70*spa54*t292; 
complex<T> co5 = d73*spa61*t292; 
complex<T> co6 = d74*spa61*t280; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t35*(*CI_users[1]->get_value(mc,ind,mu)) + t251*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + t6*(*CI_users[6]->get_value(mc,ind,mu)) + t5*(*CI_users[7]->get_value(mc,ind,mu)) + t202*(*CI_users[8]->get_value(mc,ind,mu)) + t180*(*CI_users[9]->get_value(mc,ind,mu)) + t149*(*CI_users[10]->get_value(mc,ind,mu)) + t218*(*CI_users[11]->get_value(mc,ind,mu)) + t241*(*CI_users[12]->get_value(mc,ind,mu)) + t233*(*CI_users[13]->get_value(mc,ind,mu)) + co1*(*CI_users[14]->get_value(mc,ind,mu)) + t203*(*CI_users[15]->get_value(mc,ind,mu)) + t150*(*CI_users[16]->get_value(mc,ind,mu)) + co2*(*CI_users[17]->get_value(mc,ind,mu)) + t181*(*CI_users[18]->get_value(mc,ind,mu)) + co3*(*CI_users[19]->get_value(mc,ind,mu)) + t219*(*CI_users[20]->get_value(mc,ind,mu)) + t96*(*CI_users[21]->get_value(mc,ind,mu)) + co4*(*CI_users[22]->get_value(mc,ind,mu)) + t151*(*CI_users[23]->get_value(mc,ind,mu)) + co5*(*CI_users[24]->get_value(mc,ind,mu)) + co6*(*CI_users[25]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_pmpmmm_nf_wCI::\
C6g_pmpmmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c126, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c1456));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c6, c2345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c1236));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c6, c345));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c156));
CI_users.push_back(new Cached_Box_Integral_User(c5, c34, c2, c16));
} 
  
  
template <class T> SeriesC<T> 
     C6g_pmpmmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C6g :  pmpmmm nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 17:53:30 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb31 = SPB(3,1);
complex<T> spb41 = SPB(4,1);
complex<T> spb61 = SPB(6,1);
complex<T> spb42 = SPB(4,2);
complex<T> spb54 = SPB(5,4);
complex<T> spb65 = SPB(6,5);
complex<T> spb21 = SPB(2,1);
complex<T> spb32 = SPB(3,2);
complex<T> spb43 = SPB(4,3);
complex<T> spb62 = SPB(6,2);
complex<T> spb63 = SPB(6,3);
complex<T> spb51 = SPB(5,1);
complex<T> spb52 = SPB(5,2);
complex<T> spb53 = SPB(5,3);
complex<T> spa61 = SPA(6,1);
complex<T> spa54 = SPA(5,4);
complex<T> spa62 = SPA(6,2);
complex<T> spa32 = SPA(3,2);
complex<T> spa42 = SPA(4,2);
complex<T> spa52 = SPA(5,2);
complex<T> spa65 = SPA(6,5);
complex<T> spa51 = SPA(5,1);
complex<T> spa43 = SPA(4,3);
complex<T> spa64 = SPA(6,4);
complex<T> spa41 = SPA(4,1);
complex<T> spa21 = SPA(2,1);
complex<T> s23 = -(spb32*spa32);
complex<T> s34 = -(spb43*spa43);
complex<T> s12 = -(spb21*spa21);
complex<T> s16 = -(spb61*spa61);
complex<T> s126 = SS(1,2,6);
complex<T> s156 = SS(1,5,6);
complex<T> s24 = -(spb42*spa42);
complex<T> s234 = SS(2,3,4);
complex<T> s25 = -(spb52*spa52);
complex<T> s35 = S(3,5);
complex<T> s45 = -(spb54*spa54);
complex<T> s36 = S(3,6);
complex<T> s46 = S(4,6);
complex<T> s56 = -(spb65*spa65);
complex<T> s345 = SS(3,4,5);
complex<T> t7 = square(spb53); 
complex<T> t8 = square(spb51); 
complex<T> t10 = spb54*spb65; 
complex<T> t11 = square(spb41); 
complex<T> t16 = -(spb63*spa62); 
complex<T> t18 = spb61*T(3); 
complex<T> t20 = square(spb31*spa61 + spb32*spa62); 
complex<T> t21 = square(spb31*spa41 + spb53*spa54 + spb63*spa64); 
complex<T> t36 = square(spb31); 
complex<T> t37 = square(spb63); 
complex<T> t41 = cube(spb52); 
complex<T> t46 = square(s23 + s24); 
complex<T> t49 = square(spb43*spa42 + spb53*spa52); 
complex<T> t50 = square(spb51*spa54 + spb61*spa64); 
complex<T> t53 = (s12 - s126)*T(2); 
complex<T> t63 = -(spb43*spa42) - spb53*spa52; 
complex<T> t64 = spb21*spb32; 
complex<T> t65 = -(spb31*spa51) - spb63*spa65; 
complex<T> t66 = spb21*spa52 - spb61*spa65; 
complex<T> t67 = -(spb31*spa41) - spb53*spa54 - spb63*spa64; 
complex<T> t70 = -(spb31*spa51) - spb32*spa52 - spb63*spa65; 
complex<T> t72 = -s234 + s34; 
complex<T> t73 = -(spb31*spa61) - spb32*spa62; 
complex<T> t74 = s16 - s234; 
complex<T> t76 = -(spb51*spa54) - spb61*spa64; 
complex<T> t78 = s12 - s345; 
complex<T> t84 = square(spa62); 
complex<T> t85 = -s156 + s16; 
complex<T> t100 = spa21*T(3); 
complex<T> t102 = square(spb43); 
complex<T> t103 = square(spa65); 
complex<T> t107 = spb61*spb63; 
complex<T> t113 = s12*s16; 
complex<T> t119 = spa61*spa43; 
complex<T> t127 = spb51*spb53; 
complex<T> t128 = square(spb21); 
complex<T> t134 = square(spb32); 
complex<T> t135 = spb61*T(2); 
complex<T> t136 = spb21*spa42; 
complex<T> t140 = cube(spb51); 
complex<T> t141 = cube(spb53); 
complex<T> t151 = spb31*spa51 + spb63*spa65; 
complex<T> t159 = spb41*spa42; 
complex<T> t176 = square(spb61); 
complex<T> t177 = square(spa42); 
complex<T> t210 = spb41*spb43; 
complex<T> t214 = spb52*T(6); 
complex<T> t215 = spb62*spb54; 
complex<T> t228 = spb43*spb65; 
complex<T> t278 = spb51*spb61; 
complex<T> t17 = -(spb43*t76); 
complex<T> t19 = square(spb31*spa32 + t159); 
complex<T> t23 = square(t151); 
complex<T> t24 = square(spb32*spa52 + t151); 
complex<T> t38 = -t64; 
complex<T> t39 = -t159; 
complex<T> t40 = spb43*t10; 
complex<T> t42 = -t127; 
complex<T> t48 = square(spb31*spa32 + spb51*spa52 + t159); 
complex<T> t77 = spb61*t10; 
complex<T> t105 = -(t64*T(2)); 
complex<T> t106 = spb43*t11; 
complex<T> t111 = spb42*t10; 
complex<T> t126 = spb62*t78; 
complex<T> t131 = (-(spb32*spa52) + t65)*t66; 
complex<T> t154 = square(t66); 
complex<T> t155 = t7*t8; 
complex<T> t157 = spb21*t37; 
complex<T> t174 = spb42*t72; 
complex<T> t188 = spb52*t18; 
complex<T> t192 = spb41*t67; 
complex<T> t204 = t127*square(spb31); 
complex<T> t205 = spb61*t64; 
complex<T> t209 = spb61*t73; 
complex<T> t212 = spb21*t63; 
complex<T> t225 = t135*square(spb63); 
complex<T> t256 = spb51*t7; 
complex<T> t270 = spb53*t8; 
complex<T> t295 = t141*t65; 
complex<T> t308 = t128*t177; 
complex<T> d12 = (-s126 + s16)*spb52*spb43*t135*t215; d12 = T(1)/d12;
complex<T> d18 = (s16 - s345)*spb61*spb43*t214*t215; d18 = T(1)/d18;
complex<T> d59 = spb61*t228*square(spb52); d59 = T(1)/d59;
complex<T> d60 = spb61*t228*square(square(spb52)); d60 = T(1)/d60;
complex<T> d61 = t228*square(spb62); d61 = T(1)/d61;
complex<T> d62 = t228*square(square(spb62)); d62 = T(1)/d62;
complex<T> d66 = t10*square(spb52)*T(2); d66 = T(1)/d66;
complex<T> d67 = t10*square(square(spb52)); d67 = T(1)/d67;
complex<T> t15 = -t212; 
complex<T> t44 = -t209; 
complex<T> t62 = -(spb31*spa32) + t39; 
complex<T> t68 = -(spb31*spa32) - spb51*spa52 + t39; 
complex<T> t148 = d61*spa54; 
complex<T> t156 = spb61*t40; 
complex<T> t189 = spb21*t40; 
complex<T> t190 = t36*t42; 
complex<T> t243 = spb32*t40; 
complex<T> t288 = spb53*t24; 
complex<T> d1 = t40*square(spb62); d1 = T(1)/d1;
complex<T> d2 = spb62*t40*t53; d2 = T(1)/d2;
complex<T> d3 = t126*t40*T(6); d3 = T(1)/d3;
complex<T> d4 = t40*t78*cube(spb62); d4 = T(1)/d4;
complex<T> d5 = spb62*t100*t128*t40; d5 = T(1)/d5;
complex<T> d6 = spa21*t40*cube(spb62); d6 = T(1)/d6;
complex<T> d9 = t40*square(s36 + s46 + s56)*square(spb62); d9 = T(1)/d9;
complex<T> d10 = spb62*t40*cube(t78)*T(3); d10 = T(1)/d10;
complex<T> d17 = t40*cube(spb62); d17 = T(1)/d17;
complex<T> d20 = (s16 - s345)*t40*cube(spb62); d20 = T(1)/d20;
complex<T> d22 = t40*square(s23 + s24 + s25)*square(spb62); d22 = T(1)/d22;
complex<T> d23 = spb43*t188*t215*cube(s16 - s345); d23 = T(1)/d23;
complex<T> d24 = t40*t41*t74; d24 = T(1)/d24;
complex<T> d28 = t40*square(s25 + s35 + s45)*square(spb52); d28 = T(1)/d28;
complex<T> d30 = t77*square(spb42); d30 = T(1)/d30;
complex<T> d31 = spb32*t111*t18; d31 = T(1)/d31;
complex<T> d32 = t77*cube(spb42); d32 = T(1)/d32;
complex<T> d33 = (s23 - s234)*t77*cube(spb42); d33 = T(1)/d33;
complex<T> d34 = (s23 - s234)*spb32*t111*t18; d34 = T(1)/d34;
complex<T> d35 = t77*square(s24 + s34)*square(spb42); d35 = T(1)/d35;
complex<T> d36 = spb32*t111*t18*cube(s23 - s234); d36 = T(1)/d36;
complex<T> d37 = spb65*t174*t188; d37 = T(1)/d37;
complex<T> d38 = t72*t77*cube(spb42); d38 = T(1)/d38;
complex<T> d39 = t41*t72*t77; d39 = T(1)/d39;
complex<T> d40 = t174*t188*t228; d40 = T(1)/d40;
complex<T> d42 = t46*t77*square(spb42); d42 = T(1)/d42;
complex<T> d43 = t46*t77*square(spb52); d43 = T(1)/d43;
complex<T> d44 = spb42*spb65*t188*cube(t72); d44 = T(1)/d44;
complex<T> d45 = s34*spb21*t111*t18; d45 = T(1)/d45;
complex<T> d46 = spa43*t77*cube(spb42); d46 = T(1)/d46;
complex<T> d52 = t40*square(square(spb62)); d52 = T(1)/d52;
complex<T> d53 = t40*square(spb52); d53 = T(1)/d53;
complex<T> d54 = t40*square(square(spb52)); d54 = T(1)/d54;
complex<T> d55 = t77*square(square(spb42)); d55 = T(1)/d55;
complex<T> d57 = t77*square(spb52); d57 = T(1)/d57;
complex<T> d58 = t77*square(square(spb52)); d58 = T(1)/d58;
complex<T> d63 = t40*square(spb62)*T(2); d63 = T(1)/d63;
complex<T> d64 = t77*square(spb42)*T(2); d64 = T(1)/d64;
complex<T> d65 = t135*t40*square(spb52); d65 = T(1)/d65;
complex<T> t14 = spb32*(spb51*spa52 - t62); 
complex<T> t61 = -(d59*spa54*t204) + spb63*t148*square(spb31) + d62*spb61*spa54*t105*square(spb63) + d60*spa54*t155*t64*T(2); 
complex<T> t92 = -(d52*T(2)); 
complex<T> t109 = spa42*t62; 
complex<T> t114 = d55*s23; 
complex<T> t118 = d4*spa62; 
complex<T> t137 = d30*spb41; 
complex<T> t138 = d1*spb63; 
complex<T> t169 = d23*t48; 
complex<T> t171 = d38*t62; 
complex<T> t181 = d44*spb32; 
complex<T> t186 = d10*t20; 
complex<T> t193 = d24*spb32; 
complex<T> t194 = d38*spa42; 
complex<T> t197 = d23*t49; 
complex<T> t219 = d39*spa42; 
complex<T> t237 = t63*t68; 
complex<T> t240 = t113*(d63*spb63*t36 + d52*spb61*t37*t38); 
complex<T> d7 = spb62*t189*t53; d7 = T(1)/d7;
complex<T> d8 = t126*t189*T(6); d8 = T(1)/d8;
complex<T> d11 = spb62*t100*t189; d11 = T(1)/d11;
complex<T> d13 = spb52*t243*t85*T(2); d13 = T(1)/d13;
complex<T> d14 = spb52*t135*t243*t85; d14 = T(1)/d14;
complex<T> d15 = t156*square(spb52); d15 = T(1)/d15;
complex<T> d16 = spb62*t18*t243; d16 = T(1)/d16;
complex<T> d19 = (s16 - s345)*t156*t41; d19 = T(1)/d19;
complex<T> d21 = t156*square(s23 + s24 + s25)*square(spb52); d21 = T(1)/d21;
complex<T> d25 = t214*t243*t74; d25 = T(1)/d25;
complex<T> d26 = spb32*t156*t214*t74; d26 = T(1)/d26;
complex<T> d27 = t156*t41*t74; d27 = T(1)/d27;
complex<T> d29 = spb52*t243*cube(t74)*T(3); d29 = T(1)/d29;
complex<T> d41 = t156*t41*t72; d41 = T(1)/d41;
complex<T> d47 = s34*spb42*t18*t189; d47 = T(1)/d47;
complex<T> d48 = (s34 - s345)*t188*t189; d48 = T(1)/d48;
complex<T> d49 = (s34 - s345)*t156*t41; d49 = T(1)/d49;
complex<T> d50 = t156*square(s35 + s45)*square(spb52); d50 = T(1)/d50;
complex<T> d51 = t188*t189*cube(s34 - s345); d51 = T(1)/d51;
complex<T> d56 = t156*square(square(spb52)); d56 = T(1)/d56;
complex<T> t35 = t36*(d2*spb63*spa62 + d7*t209 - d12*(spb51*spb32*spa52 + t212 - spb32*t62)); 
complex<T> t82 = -t137; 
complex<T> t83 = -t138; 
complex<T> t86 = d15*spb51; 
complex<T> t117 = d56*s345; 
complex<T> t121 = -(d13*spa65); 
complex<T> t162 = d29*t23; 
complex<T> t187 = t105*t106*t114 + s23*t137*t36; 
complex<T> t196 = d51*t154; 
complex<T> t201 = d14*t65; 
complex<T> t208 = s34*(d64*s23*spb41*t36 + t106*t114*t38); 
complex<T> t223 = d49*t66; 
complex<T> t234 = d42*t109; 
complex<T> t236 = d27*t65; 
complex<T> t238 = d49*t70; 
complex<T> t242 = t205*t92; 
complex<T> t247 = d43*t109; 
complex<T> t262 = d41*t62; 
complex<T> t1 = t16*t176*t186 + d3*spb63*spa62*t36 + d2*t16*t36 + d8*t209*t36 + t118*t205*t37 + d4*t157*t44 + d9*spa62*t157*t44 + d7*t36*t44 + d5*spb32*t36*(-(spb51*spa52) - spb61*spa62 + t62) + d6*spb32*t37*(-(spb51*spa52) - spb61*spa62 + t62) + d6*spb32*t107*(t16 + t63) - d11*t36*(t16 + t63) + t36*t83 + d10*t157*t44*t84; 
complex<T> t2 = d33*spb32*t11*t17 + d36*t11*t17*t21 + d32*spb21*spb31*t210 + d34*t17*t36 - d34*t192*t36 - d36*t102*t192*t50 + d33*spb21*t106*t67 + d35*t11*t17*t67 + t36*t82 - d31*cube(spb31); 
complex<T> t3 = d44*t134*t136*t19 + d15*t190 - spb21*t155*t236 + d29*t103*t278*t295 + d37*t136*t36 + t137*t36 + d34*t192*t36 + t106*t194*t38 + t127*t247*t38 + t256*t262*t38 + d36*t102*t192*t50 + spb43*t181*t308*t62 + d40*spb32*t36*t62 + t171*t210*t64 + t219*t270*t64 + d42*spb43*t159*t62*t64 + d26*t190*t65 - d33*spb21*t106*t67 + spa65*(-(t155*(t162 + t193 - d28*t65)) + d25*t36*t7) + d33*spb32*t106*t76 + d36*t106*t21*t76 + d34*spb43*t36*t76 + d35*t106*t67*t76; 
complex<T> t6 = d20*spb32*t107*t15 - d50*t131*t155 + d20*t14*t157 + t134*t15*t169 + spb63*spa62*t176*t186 + d15*t190 + d4*t157*t209 + d9*spa62*t157*t209 + d22*spb63*t14*t212 - spb32*t155*t223 + t138*t36 + d18*t14*t36 + d3*t16*t36 + d18*t212*t36 + spb61*t118*t37*t38 + d8*t36*t44 + spb32*t128*t197*(-(spb51*spa52) + t62) + d19*t256*(-(spb51*spa52) + t62)*t64 + d19*t270*t63*t64 + d21*t127*(-(spb51*spa52) + t62)*t63*t64 + d49*spb21*t155*(-(spb32*spa52) + t65) + t155*t196*(-(spb32*spa52) + t65) + d48*t204*t66 + d51*t140*t288*t66 + d48*t36*(-(spb32*spa52) + t65)*t8 + d10*t157*t209*t84; 
complex<T> t58 = d54*spa61*t105*t155 + d53*spa61*t204 + s16*t138*square(spb31) + s16*t242*square(spb63); 
complex<T> t180 = spb53*t86; 
complex<T> t241 = d66*t119*t190 + d65*s234*s345*t127*t36 + s234*t117*t155*t38 + d67*t119*t155*t64; 
complex<T> t253 = s12*(t138*t36 + t242*t37); 
complex<T> t254 = t201*t204 + t121*t36*t7; 
complex<T> t4 = d17*spb31*t157 + spa65*t155*t162 + spa65*t155*t193 + t128*t14*t197 + t190*t201 + d21*t127*t14*t212 + t134*t169*t212 + spb21*t155*t236 + d19*spb21*t14*t256 + d19*spb32*t15*t270 - d29*t103*t278*t295 + d12*t14*t36 + d18*t15*t36 + t180*t36 + d12*t212*t36 + d18*spb32*t36*(-(spb51*spa52) + t62) + d20*t37*(-(spb51*spa52) + t62)*t64 + d20*t107*t63*t64 + d22*spb63*(-(spb51*spa52) + t62)*t63*t64 - d28*spa65*t155*t65 + d26*t204*t65 + d13*spa65*t36*t7 - d25*spa65*t36*t7 + t36*t83 - d16*spb63*cube(spb31); 
complex<T> t5 = d50*t131*t155 - d44*t134*t136*t19 + spb32*t155*t223 - d37*t136*t36 + t180*t36 + t171*t210*t38 + t219*t270*t38 - spb43*t181*t308*t62 - d40*spb32*t36*t62 + d42*spb43*t159*t38*t62 + t106*t194*t64 + t127*t247*t64 + t256*t262*t64 + d49*spb21*t155*(spb32*spa52 - t65) + t155*t196*(spb32*spa52 - t65) + d48*t190*t66 - d51*t140*t288*t66 + d47*t11*t36*(spb32*spa42 - t67) + d46*spb21*t11*(-(spb32*spa42) + t67) - d46*spb32*t11*(t136 + t76) - d45*spb41*t36*(t136 + t76) + d48*t36*(spb32*spa52 - t65)*t8 + t36*t82; 
complex<T> t59 = d56*s156*t105*t155 + d58*spa43*t105*t155 + d57*spa43*t204 + s34*t137*square(spb31) + s156*t180*square(spb31) + s156*t82*square(spb31) + d55*t106*(s34*t105 + s156*t64*T(2)); 
complex<T> t60 = d60*spa54*t105*t155 + t105*t117*t155 + d59*spa54*t204 + d52*s345*t225*t64 + d62*spa54*t225*t64 + (-(spb63*t148) + s345*(t180 + t83))*square(spb31); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t35*(*CI_users[1]->get_value(mc,ind,mu)) + t254*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t3*(*CI_users[5]->get_value(mc,ind,mu)) + t5*(*CI_users[6]->get_value(mc,ind,mu)) + t6*(*CI_users[7]->get_value(mc,ind,mu)) + t253*(*CI_users[8]->get_value(mc,ind,mu)) + t58*(*CI_users[9]->get_value(mc,ind,mu)) + t187*(*CI_users[10]->get_value(mc,ind,mu)) + t59*(*CI_users[11]->get_value(mc,ind,mu)) + t60*(*CI_users[12]->get_value(mc,ind,mu)) + t61*(*CI_users[13]->get_value(mc,ind,mu)) + t240*(*CI_users[14]->get_value(mc,ind,mu)) + t208*(*CI_users[15]->get_value(mc,ind,mu)) + t241*(*CI_users[16]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_pmmpmm_G_wCI::\
C6g_pmmpmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c126, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c126));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c1236));
CI_users.push_back(new Cached_Bubble_Integral_User(c456, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c6, c2345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c1236));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c456));
CI_users.push_back(new Cached_Box_Integral_User(c1, c23, c4, c56));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c6, c345));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c156));
CI_users.push_back(new Cached_Box_Integral_User(c2, c34, c5, c16));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c6, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c126));
CI_users.push_back(new Cached_Box_Integral_User(c3, c45, c6, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c6, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c34, c2, c16));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c126));
CI_users.push_back(new Cached_Box_Integral_User(c5, c6, c1, c234));
CI_users.push_back(new Cached_Box_Integral_User(c6, c1, c2, c345));
} 
  
  
template <class T> SeriesC<T> 
     C6g_pmmpmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C6g :  pmmpmm G");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 18:18:54 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb21 = SPB(2,1);
complex<T> spb31 = SPB(3,1);
complex<T> spb41 = SPB(4,1);
complex<T> spb51 = SPB(5,1);
complex<T> spb61 = SPB(6,1);
complex<T> spb32 = SPB(3,2);
complex<T> spb53 = SPB(5,3);
complex<T> spb65 = SPB(6,5);
complex<T> spb43 = SPB(4,3);
complex<T> spb54 = SPB(5,4);
complex<T> spb42 = SPB(4,2);
complex<T> spb62 = SPB(6,2);
complex<T> spb64 = SPB(6,4);
complex<T> spb63 = SPB(6,3);
complex<T> spa21 = SPA(2,1);
complex<T> spa61 = SPA(6,1);
complex<T> spb52 = SPB(5,2);
complex<T> spa32 = SPA(3,2);
complex<T> spa43 = SPA(4,3);
complex<T> spa54 = SPA(5,4);
complex<T> spa65 = SPA(6,5);
complex<T> spa53 = SPA(5,3);
complex<T> spa63 = SPA(6,3);
complex<T> spa62 = SPA(6,2);
complex<T> spa42 = SPA(4,2);
complex<T> spa52 = SPA(5,2);
complex<T> spa51 = SPA(5,1);
complex<T> s34 = -(spb43*spa43);
complex<T> s45 = -(spb54*spa54);
complex<T> s12 = -(spb21*spa21);
complex<T> s16 = -(spb61*spa61);
complex<T> s123 = SS(1,2,3);
complex<T> s46 = S(4,6);
complex<T> s56 = -(spb65*spa65);
complex<T> s35 = -(spb53*spa53);
complex<T> s36 = -(spb63*spa63);
complex<T> s126 = SS(1,2,6);
complex<T> s156 = SS(1,5,6);
complex<T> s23 = -(spb32*spa32);
complex<T> s24 = -(spb42*spa42);
complex<T> s25 = -(spb52*spa52);
complex<T> s234 = SS(2,3,4);
complex<T> s345 = SS(3,4,5);
complex<T> t8 = square(spb64); 
complex<T> t9 = square(spb42); 
complex<T> t10 = square(spb31); 
complex<T> t11 = square(spb51); 
complex<T> t12 = spb32*spb65; 
complex<T> t15 = spb61*spb43; 
complex<T> t16 = spb21*spb54; 
complex<T> t24 = square(spb41*spa61 + spb42*spa62); 
complex<T> t25 = square(spb41*spa43 + spb51*spa53); 
complex<T> t27 = square(spb41*spa43 + spb51*spa53 + spb61*spa63); 
complex<T> t28 = square(spb41*spa61 + spb42*spa62 + spb43*spa63); 
complex<T> t29 = square(spb54*spa53 + spb64*spa63); 
complex<T> t44 = square(spb41); 
complex<T> t45 = -(spb42*spa62); 
complex<T> t46 = -(spb51*spa53); 
complex<T> t48 = cube(spb63); 
complex<T> t49 = -(spb64*T(2)); 
complex<T> t50 = cube(spb52); 
complex<T> t69 = spb43*spa32 - spb54*spa52; 
complex<T> t70 = cube(spb41); 
complex<T> t71 = spb21*spa52 - spb61*spa65; 
complex<T> t74 = -(spb41*spa51) - spb64*spa65; 
complex<T> t75 = -(spb41*spa51) - spb42*spa52 - spb64*spa65; 
complex<T> t77 = spb21*spa62 + spb31*spa63; 
complex<T> t78 = -(spb54*spa53) - spb64*spa63; 
complex<T> t79 = spb21*spb61; 
complex<T> t80 = -(spb31*spa32) - spb41*spa42; 
complex<T> t81 = spb43*spb54; 
complex<T> t83 = s12 - s345; 
complex<T> t84 = -s345 + s45; 
complex<T> t86 = -s234 + s34; 
complex<T> t89 = s12 - s123; 
complex<T> t96 = s12 - s126; 
complex<T> t97 = -s156 + s16; 
complex<T> t104 = s234*s345; 
complex<T> t105 = cube(spb51); 
complex<T> t106 = cube(spb42); 
complex<T> t116 = square(spa32); 
complex<T> t117 = square(spa65); 
complex<T> t121 = -(spb31*spb64); 
complex<T> t125 = square(spb43); 
complex<T> t126 = spb62*spb63; 
complex<T> t128 = square(spa62); 
complex<T> t129 = square(spa53); 
complex<T> t141 = spa61*spa43; 
complex<T> t157 = spb51*T(11); 
complex<T> t159 = square(square(spb41)); 
complex<T> t161 = square(spb21); 
complex<T> t165 = square(spb61); 
complex<T> t166 = square(spb54); 
complex<T> t168 = cube(spb31); 
complex<T> t177 = spb31*spa32 + spb41*spa42; 
complex<T> t182 = spb41*spa51 + spb64*spa65; 
complex<T> t189 = spb64*T(2); 
complex<T> t193 = spb63*T(3); 
complex<T> t199 = spb52*spb62; 
complex<T> t214 = cube(spb64); 
complex<T> t226 = spb54*spa65; 
complex<T> t234 = spa21*spa54; 
complex<T> t251 = spb21*spa32; 
complex<T> t270 = -(spb64*T(4)); 
complex<T> t315 = -(s345*T(2)); 
complex<T> t331 = spb21*spb53; 
complex<T> t400 = spb31*spb61; 
complex<T> t414 = spb51*spb61; 
complex<T> t21 = square(t177); 
complex<T> t22 = square(spb51*spa52 + t177); 
complex<T> t31 = square(t182); 
complex<T> t32 = square(spb42*spa52 + t182); 
complex<T> t47 = -t79; 
complex<T> t62 = (s16 - s234)*t12; 
complex<T> t72 = -(spb41*spa61) + t45; 
complex<T> t73 = -(spb51*spa52) + t80; 
complex<T> t76 = -(spb41*spa43) + t46; 
complex<T> t82 = -(spb41*spa43) - spb61*spa63 + t46; 
complex<T> t85 = -(spb41*spa61) - spb43*spa63 + t45; 
complex<T> t95 = t12*T(6); 
complex<T> t107 = -t141; 
complex<T> t108 = -t234; 
complex<T> t118 = -(t44*T(2)); 
complex<T> t120 = t12*t16; 
complex<T> t151 = spb53*t84; 
complex<T> t154 = spb42*t44; 
complex<T> t163 = spa62*t8; 
complex<T> t167 = square(t77); 
complex<T> t169 = square(t69); 
complex<T> t171 = square(t71); 
complex<T> t190 = t12*t15; 
complex<T> t223 = spb31*t44; 
complex<T> t224 = t11*t9; 
complex<T> t225 = -(t15*T(2)); 
complex<T> t227 = spb21*t69; 
complex<T> t229 = spb65*t126; 
complex<T> t243 = spb52*t86; 
complex<T> t248 = t10*t8; 
complex<T> t249 = t15*T(2); 
complex<T> t288 = -(t16*T(2)); 
complex<T> t289 = t10*t11; 
complex<T> t298 = t116*t80; 
complex<T> t303 = t16*T(2); 
complex<T> t304 = spb42*t11; 
complex<T> t305 = spb64*t10; 
complex<T> t314 = t79*t9; 
complex<T> t325 = spb31*t15; 
complex<T> t330 = spb54*t15; 
complex<T> t338 = spb31*t11; 
complex<T> t342 = spb32*t199; 
complex<T> t344 = spb51*t9; 
complex<T> t351 = t11*t75; 
complex<T> t358 = t15*t331; 
complex<T> t363 = spb61*t44; 
complex<T> t374 = spb64*t9; 
complex<T> t385 = spb51*t10; 
complex<T> t386 = spa32*t159; 
complex<T> t423 = spb31*t125; 
complex<T> t436 = t165*t28; 
complex<T> t439 = t166*t74; 
complex<T> d1 = t12*t81*square(spb62); d1 = T(1)/d1;
complex<T> d3 = t12*t81*t83*cube(spb62); d3 = T(1)/d3;
complex<T> d4 = spb54*t12*t48*t83; d4 = T(1)/d4;
complex<T> d8 = spa21*t12*t81*cube(spb62); d8 = T(1)/d8;
complex<T> d12 = t12*t81*square(s36 + s46 + s56)*square(spb62); d12 = T(1)/d12;
complex<T> d13 = spb54*t12*square(s36 + s46 + s56)*square(spb63); d13 = T(1)/d13;
complex<T> d22 = spb52*spb43*t12*t97; d22 = T(1)/d22;
complex<T> d26 = t12*t81*cube(spb62); d26 = T(1)/d26;
complex<T> d29 = (s16 - s345)*t12*t81*cube(spb62); d29 = T(1)/d29;
complex<T> d31 = t12*t81*square(s23 + s24 + s25)*square(spb62); d31 = T(1)/d31;
complex<T> d37 = spb43*t12*square(s25 + s35 + s45)*square(spb52); d37 = T(1)/d37;
complex<T> d38 = spb52*spb43*t12*cube(s16 - s234)*T(3); d38 = T(1)/d38;
complex<T> d39 = spb61*t12*t50*t86; d39 = T(1)/d39;
complex<T> d43 = spb61*t12*square(s23 + s24)*square(spb52); d43 = T(1)/d43;
complex<T> d44 = spb61*spb52*t12*cube(t86)*T(3); d44 = T(1)/d44;
complex<T> d45 = t12*t79*square(spb53); d45 = T(1)/d45;
complex<T> d47 = t12*t79*cube(spb53); d47 = T(1)/d47;
complex<T> d49 = (s34 - s345)*t12*t79*cube(spb53); d49 = T(1)/d49;
complex<T> d52 = t12*t79*square(s35 + s45)*square(spb53); d52 = T(1)/d52;
complex<T> d55 = t12*t79*t84*cube(spb53); d55 = T(1)/d55;
complex<T> d56 = spb21*t12*t48*t84; d56 = T(1)/d56;
complex<T> d59 = t12*t79*square(s34 + s35)*square(spb53); d59 = T(1)/d59;
complex<T> d60 = spb21*t12*square(s34 + s35)*square(spb63); d60 = T(1)/d60;
complex<T> d61 = spb32*spb53*t193*t79*cube(t84); d61 = T(1)/d61;
complex<T> d67 = spa54*t12*t79*cube(spb53); d67 = T(1)/d67;
complex<T> d69 = t12*t81*square(square(spb62)); d69 = T(1)/d69;
complex<T> d70 = spb54*t12*square(spb63); d70 = T(1)/d70;
complex<T> d71 = spb54*t12*square(square(spb63)); d71 = T(1)/d71;
complex<T> d73 = spb65*t16*square(spb63); d73 = T(1)/d73;
complex<T> d74 = spb65*t16*square(square(spb63)); d74 = T(1)/d74;
complex<T> d75 = spb43*t12*square(spb52); d75 = T(1)/d75;
complex<T> d76 = spb43*t12*square(square(spb52)); d76 = T(1)/d76;
complex<T> d78 = spb61*t12*square(spb52); d78 = T(1)/d78;
complex<T> d79 = spb61*t12*square(square(spb52)); d79 = T(1)/d79;
complex<T> d80 = t12*t79*square(square(spb53)); d80 = T(1)/d80;
complex<T> d81 = spb43*t12*square(spb62); d81 = T(1)/d81;
complex<T> d82 = spb21*t12*square(spb63); d82 = T(1)/d82;
complex<T> d83 = spb43*t12*square(square(spb62)); d83 = T(1)/d83;
complex<T> d84 = spb21*t12*square(square(spb63)); d84 = T(1)/d84;
complex<T> d90 = t12*square(spb63); d90 = T(1)/d90;
complex<T> d91 = t12*square(square(spb63)); d91 = T(1)/d91;
complex<T> d92 = t12*t79*T(2); d92 = T(1)/d92;
complex<T> d95 = t12*square(spb52); d95 = T(1)/d95;
complex<T> d96 = t12*square(square(spb52)); d96 = T(1)/d96;
complex<T> d98 = t12*t81*T(2); d98 = T(1)/d98;
complex<T> t131 = d1*s12; 
complex<T> t133 = d69*s16; 
complex<T> t136 = d45*s34; 
complex<T> t142 = d81*spa54; 
complex<T> t143 = d22*spa65; 
complex<T> t148 = d55*t76; 
complex<T> t175 = spb62*t95; 
complex<T> t186 = d61*t129; 
complex<T> t197 = d80*s34; 
complex<T> t202 = d83*spa54; 
complex<T> t208 = d61*t25; 
complex<T> t215 = t72*t79; 
complex<T> t232 = d73*spa32; 
complex<T> t250 = spb51*t154; 
complex<T> t254 = d60*spa53; 
complex<T> t256 = spb63*t151; 
complex<T> t262 = d59*t76; 
complex<T> t277 = d74*spa32; 
complex<T> t280 = d13*t72; 
complex<T> t294 = d49*t71; 
complex<T> t307 = spb43*t82; 
complex<T> t316 = d45*spb51; 
complex<T> t317 = d1*spb64; 
complex<T> t337 = t314*t8; 
complex<T> t353 = d29*spb21; 
complex<T> t354 = d52*spb54; 
complex<T> t355 = spb52*t62; 
complex<T> t359 = d3*spb42; 
complex<T> t362 = t224*t303; 
complex<T> t367 = t229*t83; 
complex<T> t368 = t289*t81; 
complex<T> t387 = d31*spb21; 
complex<T> t406 = spb42*t298; 
complex<T> d2 = t120*square(spb63); d2 = T(1)/d2;
complex<T> d5 = t229*t81*t96; d5 = T(1)/d5;
complex<T> d9 = spb43*t16*t229*t96; d9 = T(1)/d9;
complex<T> d11 = t120*t48*t83; d11 = T(1)/d11;
complex<T> d14 = t229*t81*cube(t83)*T(3); d14 = T(1)/d14;
complex<T> d16 = spb63*t16*t89*t95; d16 = T(1)/d16;
complex<T> d17 = t120*t48*t89; d17 = T(1)/d17;
complex<T> d18 = t120*square(s34 + s35 + s36)*square(spb63); d18 = T(1)/d18;
complex<T> d19 = t120*t193*cube(t89); d19 = T(1)/d19;
complex<T> d20 = spb63*t120*t89; d20 = T(1)/d20;
complex<T> d21 = (-s126 + s16)*t330*t342; d21 = T(1)/d21;
complex<T> d23 = spb52*t190*t97; d23 = T(1)/d23;
complex<T> d24 = t190*square(spb52); d24 = T(1)/d24;
complex<T> d27 = (s16 - s345)*t330*t342*T(6); d27 = T(1)/d27;
complex<T> d28 = (s16 - s345)*t190*t50; d28 = T(1)/d28;
complex<T> d30 = t190*square(s23 + s24 + s25)*square(spb52); d30 = T(1)/d30;
complex<T> d32 = t330*t342*cube(s16 - s345)*T(3); d32 = T(1)/d32;
complex<T> d33 = spb43*t50*t62; d33 = T(1)/d33;
complex<T> d36 = t15*t50*t62; d36 = T(1)/d36;
complex<T> d40 = spb61*t243*t95; d40 = T(1)/d40;
complex<T> d41 = t15*t243*t95; d41 = T(1)/d41;
complex<T> d42 = t190*t50*t86; d42 = T(1)/d42;
complex<T> d46 = t358*t95; d46 = T(1)/d46;
complex<T> d48 = (s34 - s345)*t190*t50; d48 = T(1)/d48;
complex<T> d50 = (s34 - s345)*spb52*spb65*t358*T(6); d50 = T(1)/d50;
complex<T> d51 = t190*square(s35 + s45)*square(spb52); d51 = T(1)/d51;
complex<T> d53 = spb52*spb65*t358*cube(s34 - s345)*T(3); d53 = T(1)/d53;
complex<T> d58 = t120*t48*t84; d58 = T(1)/d58;
complex<T> d62 = (-s123 + s45)*spb63*t16*t95; d62 = T(1)/d62;
complex<T> d63 = (-s123 + s45)*t120*t48; d63 = T(1)/d63;
complex<T> d64 = t120*square(s46 + s56)*square(spb63); d64 = T(1)/d64;
complex<T> d65 = t120*t193*cube(-s123 + s45); d65 = T(1)/d65;
complex<T> d66 = s45*spb53*t79*t95; d66 = T(1)/d66;
complex<T> d68 = s45*spb61*spb53*t16*t95; d68 = T(1)/d68;
complex<T> d72 = t120*square(square(spb63)); d72 = T(1)/d72;
complex<T> d77 = t190*square(square(spb52)); d77 = T(1)/d77;
complex<T> d85 = spb54*spb65*t249; d85 = T(1)/d85;
complex<T> d86 = t120*t249; d86 = T(1)/d86;
complex<T> d87 = t16*t249; d87 = T(1)/d87;
complex<T> d88 = spb61*spb65*t303; d88 = T(1)/d88;
complex<T> d89 = t120*T(2); d89 = T(1)/d89;
complex<T> d93 = t190*T(2); d93 = T(1)/d93;
complex<T> d94 = spb21*spb32*t249; d94 = T(1)/d94;
complex<T> d97 = spb32*spb43*t303; d97 = T(1)/d97;
complex<T> t66 = -(d16*(spb41*spa43 + spb51*spa53 + spb61*spa63)); 
complex<T> t109 = d86*s123; 
complex<T> t123 = -(d24*spb51); 
complex<T> t132 = d2*s123; 
complex<T> t139 = d77*s45; 
complex<T> t149 = d62*t77; 
complex<T> t196 = d72*s123; 
complex<T> t198 = d24*s45; 
complex<T> t205 = d65*t167; 
complex<T> t207 = d14*t24; 
complex<T> t218 = d14*t128; 
complex<T> t219 = (d86*t104 + d89*t107)*t159; 
complex<T> t222 = t197*t368*T(2) - spb51*t136*t223*T(4); 
complex<T> t235 = d53*t171; 
complex<T> t238 = d23*t74; 
complex<T> t242 = d41*t80; 
complex<T> t257 = d32*t169; 
complex<T> t265 = d17*t82; 
complex<T> t268 = t225*t248*t277 + spb64*t223*t232*T(4); 
complex<T> t269 = d24*spb51*spb42*t104*t118 + d77*t104*t16*t224 + d96*t107*t16*t224 + d95*t141*t250*T(2); 
complex<T> t282 = t186*t76; 
complex<T> t283 = d64*t77; 
complex<T> t295 = d9*t72; 
complex<T> t297 = d18*t78; 
complex<T> t301 = s45*(spb31*spb51*t118*t136 + t197*t368); 
complex<T> t302 = s12*t133*t337 + s16*t131*t154*t49; 
complex<T> t312 = d20*(t118*t307 + t223*t78*T(2)); 
complex<T> t321 = d71*spa21*t248*t249 + t131*t154*t270 + d70*spa21*t223*t270 + d69*s12*t337*T(2); 
complex<T> t329 = -(t250*T(4)); 
complex<T> t341 = d51*spb54; 
complex<T> t348 = spb32*t256; 
complex<T> t365 = d2*spb64; 
complex<T> t375 = d30*spb21; 
complex<T> t399 = d28*spb51; 
complex<T> t441 = d2*t121; 
complex<T> d6 = t367*t81*T(6); d6 = T(1)/d6;
complex<T> d7 = spa21*t161*t175*t81; d7 = T(1)/d7;
complex<T> d10 = spb43*t16*t367*T(6); d10 = T(1)/d10;
complex<T> d15 = spb43*spa21*t16*t175; d15 = T(1)/d15;
complex<T> d25 = t175*t330; d25 = T(1)/d25;
complex<T> d34 = spb43*t355*T(6); d34 = T(1)/d34;
complex<T> d35 = t15*t355*T(6); d35 = T(1)/d35;
complex<T> t1 = t223*t365 + spb43*t44*t66 + d63*spb61*t248*(spb43*spa63 - t72) + d17*spb43*t248*(spb61*spa63 - t76) + d19*spb43*t10*t29*(spb61*spa63 - t76) + spb43*t297*t305*(spb61*spa63 - t76) - d65*spb64*t436*t77 - d16*t223*t78 - d17*t15*t305*t78 + d19*t27*t423*t78 + spb61*t205*(spb43*spa63 - t72)*t8 + t283*t400*(-(spb43*spa63) + t72)*t8 + d63*t325*t77*t8 - spb64*t149*t44*T(11) + d62*t363*(spb43*spa63 - t72)*T(11); 
complex<T> t2 = t163*t165*t207 + d12*spb42*t163*t215 + t214*t215*t218 + t154*t317 + d4*t163*t325 + t189*t295*t363 - t163*t280*t400 - d6*t163*t44 + t44*t441 + d8*spb61*t374*(spb64*spa62 - t69) - d11*spb61*t248*t72 - d10*spb64*t363*t72 + d17*spb43*t248*(-(spb61*spa63) + t76) + d19*spb43*t10*t29*(-(spb61*spa63) + t76) + spb43*t297*t305*(-(spb61*spa63) + t76) + d17*t15*t305*t78 - d19*t27*t423*t78 + t215*t359*t8 + d3*t163*t47*t9 + d8*(spb61*spa62 - t73)*t8*t9 + d5*t163*t44*T(2) + spb43*t44*t66*T(11) + d15*t154*(-(spb64*spa62) + t69)*T(11) - d16*t223*t78*T(11) + d7*t44*(spb61*spa62 - t73)*t9*T(11); 
complex<T> t4 = spb54*t118*t143 + t123*t154 - d33*t224*t226 + d21*spb42*t118*t227 + d27*t154*t227 - d32*t106*t22*t227 + d28*t224*t227 - d38*t11*t226*t31 + t154*t317 + d38*t117*t414*t439 + d34*t226*t44 + d29*t374*t47*t69 + t344*t375*t69*t73 - t374*t387*t69*t73 - d36*t16*t304*t74 + d37*t226*t304*t74 - d35*spb51*t44*t74 - d26*spb21*spb41*spb42*t8 + t161*t257*t73*t9 + t16*t399*t73*t9 - d27*t44*t73*t9 - t353*t73*t8*t9 + spb51*t238*t44*T(2) + d21*t44*t73*t9*T(2) + d25*spb64*t70*T(11); 
complex<T> t5 = d33*t224*t226 + d24*t250 + d39*t224*t251 + d38*t11*t226*t31 - d44*spb43*t161*t406 - d38*t117*t414*t439 - d34*t226*t44 + d36*t16*t304*t74 - d37*t226*t304*t74 + d35*spb51*t44*t74 + d42*t16*t344*t80 - d43*t251*t344*t80 + d44*t21*t251*t9 - t154*t242*T(11) + d40*t251*t44*T(11); 
complex<T> t7 = t123*t154 - d39*t224*t251 + t223*t316 - d47*spb41*spb54*t385 + d44*spb43*t161*t406 - d48*spb54*t224*t71 + d53*spb54*t105*t32*t71 + d50*spb54*t157*t44*t71 + d49*spb54*t289*(spb42*spa52 - t74) + t304*t341*t71*(spb42*spa52 - t74) + t11*t166*t235*(-(spb42*spa52) + t74) + d48*t16*t304*(-(spb42*spa52) + t74) + t338*t354*t71*(-(spb42*spa52) + t74) - d42*t16*t344*t80 + d43*t251*t344*t80 + t294*t338*t81 - d44*t21*t251*t9 + t154*t242*T(11) - d40*t251*t44*T(11) + d46*spb31*t70*T(11) + d50*t11*t44*(-(spb42*spa52) + t74)*T(11); 
complex<T> t43 = t295*t363*t49 + d21*t118*t73*square(spb42) + d5*spa62*t118*square(spb64) + d21*t154*t227*T(2); 
complex<T> t245 = s234*t109*t159 - d87*spa65*t386; 
complex<T> t246 = d24*s156*t329 + d78*spa43*t329 + d77*s156*t362 + d79*spa43*t362; 
complex<T> t247 = d90*t189*t223*t234 + d91*t108*t15*t248 + s126*t15*t196*t248 + s126*t132*t223*t49; 
complex<T> t267 = (d93*t108 + s345*t109)*t159; 
complex<T> t287 = t196*t248*t249 + t132*t223*t270 + t223*t232*t270 + t248*t249*t277; 
complex<T> t300 = spb51*t118*t238 + spb54*t143*t44*T(2); 
complex<T> t328 = d1*s16*t154*t270 + d75*spa61*t329 + d76*spa61*t362 + t133*t337*T(2); 
complex<T> t336 = t142*t154*t270 + t198*t329 + t139*t362 + t202*t337*T(2); 
complex<T> t343 = d72*s345*t248*t249 + d84*spa54*t248*t249 + d2*s345*t223*t270 + d82*spa54*t223*t270 + t139*t224*t288 + d24*s345*t329 + d69*t315*t337 + d77*s345*t362 + d80*t315*t368 - t202*t337*T(2) + d80*s45*t368*T(2) + spb64*t142*t154*T(4) + t198*t250*T(4) + s345*t223*t316*T(4) - s45*t223*t316*T(4) + s345*t154*t317*T(4); 
complex<T> t452 = spb31*t283; 
complex<T> d54 = t348*t79*T(6); d54 = T(1)/d54;
complex<T> d57 = spb61*t16*t348*T(6); d57 = T(1)/d57;
complex<T> t201 = d54*spa53; 
complex<T> t239 = d57*t76; 
complex<T> t3 = spa53*t10*t125*t208 + d56*spa53*t15*t305 + t223*t316 - d55*spa53*t368 + t44*t441 + d67*spb43*t338*(spb31*spa53 + t71) + d66*t157*t44*(spb31*spa53 + t71) + d63*spb61*t248*(-(spb43*spa63) + t72) + d67*t289*(spb42*spa52 + spb43*spa53 - t74) - d58*spb43*t248*t76 - spb43*t254*t305*t76 + d65*spb64*t436*t77 + t283*t400*(spb43*spa63 - t72)*t8 + spb61*t205*(-(spb43*spa63) + t72)*t8 - d63*t325*t77*t8 + t168*t282*t81 + t148*t385*t81 + spa53*t262*t385*t81 + spb43*t223*t239*T(11) + spb64*t149*t44*T(11) + t10*t201*t44*T(11) + d62*t363*(-(spb43*spa63) + t72)*T(11) - d68*t11*t44*(spb42*spa52 + spb43*spa53 - t74)*T(11); 
complex<T> t6 = -(t163*t165*t207) - spa53*t10*t125*t208 - d27*t154*t227 + d32*t106*t22*t227 - d28*t224*t227 + d24*t250 - d56*spa53*t15*t305 + d3*t163*t314 - t223*t316 - t154*t317 - d4*t163*t325 + t223*t365 + d55*spa53*t368 + t163*t280*t400 + d6*t163*t44 + d29*spb64*t314*t69 + d48*spb54*t224*t71 - d53*spb54*t105*t32*t71 + d11*spb61*t248*t72 + d10*spb64*t363*t72 + t214*t218*t47*t72 - t344*t375*t69*t73 + t374*t387*t69*t73 + t11*t166*t235*(spb42*spa52 - t74) + d48*t16*t304*(spb42*spa52 - t74) + t338*t354*t71*(spb42*spa52 - t74) + d49*spb54*t289*(-(spb42*spa52) + t74) + t304*t341*t71*(-(spb42*spa52) + t74) + d58*spb43*t248*t76 + spb43*t254*t305*t76 + d12*t215*t45*t8 + t359*t47*t72*t8 - t168*t282*t81 - t294*t338*t81 - t148*t385*t81 + t10*t262*t46*t81 - t161*t257*t73*t9 - t16*t399*t73*t9 + d27*t44*t73*t9 + t353*t73*t8*t9 - spb43*t223*t239*T(11) - t10*t201*t44*T(11) - d50*spb51*spb54*t44*t71*T(11) + d50*t11*t44*(spb42*spa52 - t74)*T(11); 
complex<T> co1 = d85*spa21*t386; 
complex<T> co2 = d88*spa43*t386; 
complex<T> co3 = d92*spa43*spa54*t159; 
complex<T> co4 = d94*spa54*spa65*t159; 
complex<T> co5 = d97*spa61*spa65*t159; 
complex<T> co6 = d98*spa21*spa61*t159; 
complex<T> co7 = Complex(0,1); 
SeriesC<T> result = co7*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t312*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t300*(*CI_users[3]->get_value(mc,ind,mu)) + t4*(*CI_users[4]->get_value(mc,ind,mu)) + t5*(*CI_users[5]->get_value(mc,ind,mu)) + t7*(*CI_users[6]->get_value(mc,ind,mu)) + t6*(*CI_users[7]->get_value(mc,ind,mu)) + t3*(*CI_users[8]->get_value(mc,ind,mu)) + t1*(*CI_users[9]->get_value(mc,ind,mu)) + t321*(*CI_users[10]->get_value(mc,ind,mu)) + t287*(*CI_users[11]->get_value(mc,ind,mu)) + t328*(*CI_users[12]->get_value(mc,ind,mu)) + t268*(*CI_users[13]->get_value(mc,ind,mu)) + t246*(*CI_users[14]->get_value(mc,ind,mu)) + t222*(*CI_users[15]->get_value(mc,ind,mu)) + t343*(*CI_users[16]->get_value(mc,ind,mu)) + t336*(*CI_users[17]->get_value(mc,ind,mu)) + co1*(*CI_users[18]->get_value(mc,ind,mu)) + t245*(*CI_users[19]->get_value(mc,ind,mu)) + t302*(*CI_users[20]->get_value(mc,ind,mu)) + co2*(*CI_users[21]->get_value(mc,ind,mu)) + t219*(*CI_users[22]->get_value(mc,ind,mu)) + t247*(*CI_users[23]->get_value(mc,ind,mu)) + co3*(*CI_users[24]->get_value(mc,ind,mu)) + t267*(*CI_users[25]->get_value(mc,ind,mu)) + co4*(*CI_users[26]->get_value(mc,ind,mu)) + t269*(*CI_users[27]->get_value(mc,ind,mu)) + t301*(*CI_users[28]->get_value(mc,ind,mu)) + co5*(*CI_users[29]->get_value(mc,ind,mu)) + co6*(*CI_users[30]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C6g_pmmpmm_nf_wCI::\
C6g_pmmpmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c126, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c156, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c16, c2345));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c156));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c1256));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c126));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c1236));
CI_users.push_back(new Cached_Bubble_Integral_User(c456, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c3456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c456));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c6, c2345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c1456));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c34, c561));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c1256));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c126));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c1236));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c6, c345));
CI_users.push_back(new Cached_Box_Integral_User(c3, c12, c6, c45));
CI_users.push_back(new Cached_Box_Integral_User(c5, c34, c2, c16));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c126));
} 
  
  
template <class T> SeriesC<T> 
     C6g_pmmpmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C6g :  pmmpmm nf");
#endif
 
//#define TimeStamp "Mon 8 Nov 2010 18:37:59 on n2179"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb21 = SPB(2,1);
complex<T> spb31 = SPB(3,1);
complex<T> spb41 = SPB(4,1);
complex<T> spb51 = SPB(5,1);
complex<T> spb61 = SPB(6,1);
complex<T> spb32 = SPB(3,2);
complex<T> spb53 = SPB(5,3);
complex<T> spb65 = SPB(6,5);
complex<T> spb43 = SPB(4,3);
complex<T> spb54 = SPB(5,4);
complex<T> spb42 = SPB(4,2);
complex<T> spb62 = SPB(6,2);
complex<T> spb64 = SPB(6,4);
complex<T> spb63 = SPB(6,3);
complex<T> spa21 = SPA(2,1);
complex<T> spb52 = SPB(5,2);
complex<T> spa61 = SPA(6,1);
complex<T> spa32 = SPA(3,2);
complex<T> spa54 = SPA(5,4);
complex<T> spa43 = SPA(4,3);
complex<T> spa53 = SPA(5,3);
complex<T> spa63 = SPA(6,3);
complex<T> spa62 = SPA(6,2);
complex<T> spa42 = SPA(4,2);
complex<T> spa52 = SPA(5,2);
complex<T> spa65 = SPA(6,5);
complex<T> spa51 = SPA(5,1);
complex<T> s34 = -(spb43*spa43);
complex<T> s45 = -(spb54*spa54);
complex<T> s12 = -(spb21*spa21);
complex<T> s16 = -(spb61*spa61);
complex<T> s123 = SS(1,2,3);
complex<T> s46 = S(4,6);
complex<T> s56 = -(spb65*spa65);
complex<T> s35 = -(spb53*spa53);
complex<T> s36 = -(spb63*spa63);
complex<T> s126 = SS(1,2,6);
complex<T> s156 = SS(1,5,6);
complex<T> s23 = -(spb32*spa32);
complex<T> s24 = -(spb42*spa42);
complex<T> s25 = -(spb52*spa52);
complex<T> s234 = SS(2,3,4);
complex<T> s345 = SS(3,4,5);
complex<T> t8 = square(spb64); 
complex<T> t9 = square(spb42); 
complex<T> t10 = square(spb31); 
complex<T> t11 = square(spb51); 
complex<T> t12 = spb32*spb65; 
complex<T> t13 = square(spb63); 
complex<T> t16 = spb61*spb43; 
complex<T> t17 = spb21*spb54; 
complex<T> t24 = square(spb41*spa61 + spb42*spa62); 
complex<T> t25 = square(spb41*spa43 + spb51*spa53); 
complex<T> t27 = square(spb41*spa43 + spb51*spa53 + spb61*spa63); 
complex<T> t28 = square(spb41*spa61 + spb42*spa62 + spb43*spa63); 
complex<T> t29 = square(spb54*spa53 + spb64*spa63); 
complex<T> t44 = square(spb41); 
complex<T> t45 = -(spb42*spa62); 
complex<T> t46 = -(spb51*spa53); 
complex<T> t48 = cube(spb63); 
complex<T> t49 = cube(spb52); 
complex<T> t70 = spb43*spa32 - spb54*spa52; 
complex<T> t71 = spb21*spa52 - spb61*spa65; 
complex<T> t74 = -(spb41*spa51) - spb64*spa65; 
complex<T> t75 = -(spb41*spa51) - spb42*spa52 - spb64*spa65; 
complex<T> t77 = spb21*spa62 + spb31*spa63; 
complex<T> t78 = -(spb54*spa53) - spb64*spa63; 
complex<T> t79 = spb21*spb61; 
complex<T> t80 = -(spb31*spa32) - spb41*spa42; 
complex<T> t82 = spb43*spb54; 
complex<T> t83 = s12 - s345; 
complex<T> t84 = -s345 + s45; 
complex<T> t85 = -(spb51*spb42); 
complex<T> t86 = spb31*spb64; 
complex<T> t89 = s16 - s234; 
complex<T> t90 = -s234 + s34; 
complex<T> t91 = -square(spb41); 
complex<T> t97 = s12 - s126; 
complex<T> t98 = -s156 + s16; 
complex<T> t115 = square(spa32); 
complex<T> t116 = square(spa65); 
complex<T> t123 = spb63*T(3); 
complex<T> t126 = spb52*spb62; 
complex<T> t137 = spa61*spa43; 
complex<T> t153 = square(spb54); 
complex<T> t162 = square(spb61); 
complex<T> t163 = square(spb43); 
complex<T> t169 = square(spa53); 
complex<T> t174 = spb31*spa32 + spb41*spa42; 
complex<T> t178 = spb41*spa51 + spb64*spa65; 
complex<T> t189 = spb52*T(3); 
complex<T> t190 = square(spa62); 
complex<T> t214 = square(spb21); 
complex<T> t217 = cube(spb51); 
complex<T> t236 = spb62*spb63; 
complex<T> t251 = cube(spb31); 
complex<T> t252 = cube(spb42); 
complex<T> t253 = cube(spb64); 
complex<T> t254 = s45*spb51; 
complex<T> t267 = spb63*T(6); 
complex<T> t270 = spa21*spa54; 
complex<T> t305 = spb31*spb61; 
complex<T> t372 = spb31*spb43; 
complex<T> t411 = spb42*spb43; 
complex<T> t21 = square(t174); 
complex<T> t22 = square(spb51*spa52 + t174); 
complex<T> t31 = square(t178); 
complex<T> t32 = square(spb42*spa52 + t178); 
complex<T> t47 = -t79; 
complex<T> t50 = -t82; 
complex<T> t52 = t12*T(3); 
complex<T> t72 = -(spb51*spa52) + t80; 
complex<T> t73 = -(spb41*spa61) + t45; 
complex<T> t76 = -(spb41*spa43) + t46; 
complex<T> t81 = -(spb41*spa43) - spb61*spa63 + t46; 
complex<T> t87 = -(spb41*spa61) - spb43*spa63 + t45; 
complex<T> t118 = -(t10*t16); 
complex<T> t119 = t12*T(2); 
complex<T> t147 = spb53*t84; 
complex<T> t154 = t71*(-(spb42*spa52) + t74); 
complex<T> t160 = spb21*t70; 
complex<T> t164 = spa62*t8; 
complex<T> t166 = square(t77); 
complex<T> t171 = spb54*t71; 
complex<T> t184 = spb42*t44; 
complex<T> t185 = t10*t8; 
complex<T> t187 = -(t16*T(2)); 
complex<T> t191 = square(t70); 
complex<T> t192 = square(t71); 
complex<T> t211 = spb62*t83; 
complex<T> t213 = t10*t11; 
complex<T> t229 = t11*t9; 
complex<T> t230 = t12*t16; 
complex<T> t231 = -(t17*T(2)); 
complex<T> t239 = -t270; 
complex<T> t241 = spb53*t189; 
complex<T> t260 = spb54*t11; 
complex<T> t266 = spb31*t78; 
complex<T> t278 = spb52*t89; 
complex<T> t283 = t79*t8; 
complex<T> t284 = spb51*t44; 
complex<T> t285 = -(t9*T(2)); 
complex<T> t302 = spb64*t44; 
complex<T> t303 = spb21*t9; 
complex<T> t318 = t12*t17; 
complex<T> t320 = spb64*t10; 
complex<T> t329 = spb52*t90; 
complex<T> t341 = t236*t97; 
complex<T> t349 = spb52*t98; 
complex<T> t353 = spb65*t17; 
complex<T> t359 = spb31*t16; 
complex<T> t361 = t9*T(2); 
complex<T> t384 = spb64*t9; 
complex<T> t388 = spb51*t10; 
complex<T> t389 = spb21*t16; 
complex<T> t417 = t44*t85; 
complex<T> t420 = spb61*t153; 
complex<T> t422 = t10*t25; 
complex<T> t426 = spb31*t11; 
complex<T> t433 = t11*t153; 
complex<T> t435 = t214*t80; 
complex<T> t438 = spb51*t17; 
complex<T> t441 = t162*t28; 
complex<T> d1 = t12*t82*square(spb62); d1 = T(1)/d1;
complex<T> d3 = t12*t82*t83*cube(spb62); d3 = T(1)/d3;
complex<T> d4 = spb54*t12*t48*t83; d4 = T(1)/d4;
complex<T> d6 = spb65*t236*t82*t83*T(6); d6 = T(1)/d6;
complex<T> d8 = spa21*t12*t82*cube(spb62); d8 = T(1)/d8;
complex<T> d12 = t12*t82*square(s36 + s46 + s56)*square(spb62); d12 = T(1)/d12;
complex<T> d13 = spb54*t12*t13*square(s36 + s46 + s56); d13 = T(1)/d13;
complex<T> d14 = spb62*spb65*t123*t82*cube(t83); d14 = T(1)/d14;
complex<T> d21 = (-s126 + s16)*spb32*spb54*t126*t16*T(2); d21 = T(1)/d21;
complex<T> d26 = t12*t82*cube(spb62); d26 = T(1)/d26;
complex<T> d27 = (s16 - s345)*spb32*spb54*t126*t16*T(6); d27 = T(1)/d27;
complex<T> d29 = (s16 - s345)*t12*t82*cube(spb62); d29 = T(1)/d29;
complex<T> d31 = t12*t82*square(s23 + s24 + s25)*square(spb62); d31 = T(1)/d31;
complex<T> d32 = spb32*spb54*t126*t16*cube(s16 - s345)*T(3); d32 = T(1)/d32;
complex<T> d33 = spb43*t12*t49*t89; d33 = T(1)/d33;
complex<T> d37 = spb43*t12*square(s25 + s35 + s45)*square(spb52); d37 = T(1)/d37;
complex<T> d39 = spb61*t12*t49*t90; d39 = T(1)/d39;
complex<T> d43 = spb61*t12*square(s23 + s24)*square(spb52); d43 = T(1)/d43;
complex<T> d45 = t12*t79*square(spb53); d45 = T(1)/d45;
complex<T> d47 = t12*t79*cube(spb53); d47 = T(1)/d47;
complex<T> d49 = (s34 - s345)*t12*t79*cube(spb53); d49 = T(1)/d49;
complex<T> d52 = t12*t79*square(s35 + s45)*square(spb53); d52 = T(1)/d52;
complex<T> d55 = t12*t79*t84*cube(spb53); d55 = T(1)/d55;
complex<T> d56 = spb21*t12*t48*t84; d56 = T(1)/d56;
complex<T> d59 = t12*t79*square(s34 + s35)*square(spb53); d59 = T(1)/d59;
complex<T> d60 = spb21*t12*t13*square(s34 + s35); d60 = T(1)/d60;
complex<T> d61 = spb32*spb53*t123*t79*cube(t84); d61 = T(1)/d61;
complex<T> d67 = spa54*t12*t79*cube(spb53); d67 = T(1)/d67;
complex<T> d70 = t12*t82*square(square(spb62)); d70 = T(1)/d70;
complex<T> d71 = spb54*t12*t13; d71 = T(1)/d71;
complex<T> d72 = spb54*t12*square(square(spb63)); d72 = T(1)/d72;
complex<T> d76 = spb43*t12*square(spb52); d76 = T(1)/d76;
complex<T> d77 = spb43*t12*square(square(spb52)); d77 = T(1)/d77;
complex<T> d79 = spb61*t12*square(spb52); d79 = T(1)/d79;
complex<T> d80 = spb61*t12*square(square(spb52)); d80 = T(1)/d80;
complex<T> d81 = t12*t79*square(square(spb53)); d81 = T(1)/d81;
complex<T> d82 = spb43*t12*square(spb62); d82 = T(1)/d82;
complex<T> d83 = spb21*t12*t13; d83 = T(1)/d83;
complex<T> d84 = spb43*t12*square(square(spb62)); d84 = T(1)/d84;
complex<T> d85 = spb21*t12*square(square(spb63)); d85 = T(1)/d85;
complex<T> d89 = t12*square(square(spb63)); d89 = T(1)/d89;
complex<T> d92 = t12*square(square(spb52)); d92 = T(1)/d92;
complex<T> t56 = -(t70*t72); 
complex<T> t57 = -t154; 
complex<T> t114 = spa21*t52; 
complex<T> t117 = -t184; 
complex<T> t121 = -t260; 
complex<T> t127 = d82*spa54; 
complex<T> t128 = d70*s12; 
complex<T> t132 = d81*s34; 
complex<T> t144 = d55*t76; 
complex<T> t149 = -t284; 
complex<T> t165 = d45*spb31; 
complex<T> t200 = d84*spa54; 
complex<T> t215 = t72*t9; 
complex<T> t219 = t76*t82; 
complex<T> t225 = d61*t169; 
complex<T> t227 = d14*t190; 
complex<T> t228 = d32*t191; 
complex<T> t244 = d29*t70; 
complex<T> t262 = d1*spb64; 
complex<T> t263 = spb43*t81; 
complex<T> t273 = d3*t73; 
complex<T> t310 = t123*t147; 
complex<T> t311 = d49*t71; 
complex<T> t332 = t185*t187; 
complex<T> t333 = t17*t52; 
complex<T> t336 = d37*spb42; 
complex<T> t338 = d60*spa53; 
complex<T> t362 = d43*spb51; 
complex<T> t366 = d29*spb21; 
complex<T> t369 = spb65*t241; 
complex<T> t413 = d59*t10; 
complex<T> t442 = d3*t164; 
complex<T> t485 = t420*t74; 
complex<T> d2 = t13*t318; d2 = T(1)/d2;
complex<T> d5 = spb65*t341*t82*T(2); d5 = T(1)/d5;
complex<T> d9 = spb43*t341*t353*T(2); d9 = T(1)/d9;
complex<T> d10 = spb43*t236*t353*t83*T(6); d10 = T(1)/d10;
complex<T> d11 = t318*t48*t83; d11 = T(1)/d11;
complex<T> d17 = (s12 - s123)*t318*t48; d17 = T(1)/d17;
complex<T> d18 = t13*t318*square(s34 + s35 + s36); d18 = T(1)/d18;
complex<T> d20 = (s12 - s123)*spb63*t119*t17; d20 = T(1)/d20;
complex<T> d22 = spb43*t119*t349; d22 = T(1)/d22;
complex<T> d23 = t119*t16*t349; d23 = T(1)/d23;
complex<T> d24 = t230*square(spb52); d24 = T(1)/d24;
complex<T> d25 = spb62*spb54*t16*t52; d25 = T(1)/d25;
complex<T> d28 = (s16 - s345)*t230*t49; d28 = T(1)/d28;
complex<T> d30 = t230*square(s23 + s24 + s25)*square(spb52); d30 = T(1)/d30;
complex<T> d34 = spb43*t12*t278*T(6); d34 = T(1)/d34;
complex<T> d35 = t230*t278*T(6); d35 = T(1)/d35;
complex<T> d36 = t230*t49*t89; d36 = T(1)/d36;
complex<T> d38 = spb52*spb43*t52*cube(t89); d38 = T(1)/d38;
complex<T> d40 = spb61*t329*t52; d40 = T(1)/d40;
complex<T> d41 = t16*t329*t52; d41 = T(1)/d41;
complex<T> d42 = t230*t49*t90; d42 = T(1)/d42;
complex<T> d44 = spb61*spb52*t52*cube(t90); d44 = T(1)/d44;
complex<T> d46 = spb53*t389*t52; d46 = T(1)/d46;
complex<T> d48 = (s34 - s345)*t230*t49; d48 = T(1)/d48;
complex<T> d51 = t230*square(s35 + s45)*square(spb52); d51 = T(1)/d51;
complex<T> d58 = t318*t48*t84; d58 = T(1)/d58;
complex<T> d63 = (-s123 + s45)*t318*t48; d63 = T(1)/d63;
complex<T> d64 = t13*t318*square(s46 + s56); d64 = T(1)/d64;
complex<T> d66 = s45*spb53*t52*t79; d66 = T(1)/d66;
complex<T> d69 = (s12 - s123)*t267*t318; d69 = T(1)/d69;
complex<T> d73 = t318*square(square(spb63)); d73 = T(1)/d73;
complex<T> d74 = t13*t353; d74 = T(1)/d74;
complex<T> d75 = t353*square(square(spb63)); d75 = T(1)/d75;
complex<T> d78 = t230*square(square(spb52)); d78 = T(1)/d78;
complex<T> d86 = t119*t82*square(spb62); d86 = T(1)/d86;
complex<T> d87 = t119*t13*t17; d87 = T(1)/d87;
complex<T> d88 = t119*t13; d88 = T(1)/d88;
complex<T> d90 = t119*t16*square(spb52); d90 = T(1)/d90;
complex<T> d91 = t119*square(spb52); d91 = T(1)/d91;
complex<T> d93 = t119*t79*square(spb53); d93 = T(1)/d93;
complex<T> t43 = d9*spb61*t302*t73 + d21*spb42*(spb42*t44*t72 + t160*t91) + d5*spa62*t44*square(spb64); 
complex<T> t125 = d24*s45; 
complex<T> t133 = -(d78*s345); 
complex<T> t138 = -(d22*spa65); 
complex<T> t148 = -(d20*t44*(spb61*spb43*spa63 + t266 - spb43*t76)); 
complex<T> t168 = d24*spb51; 
complex<T> t193 = d73*s123; 
complex<T> t194 = d38*t31; 
complex<T> t216 = d2*t86; 
complex<T> t234 = d48*t75; 
complex<T> t235 = d78*s45; 
complex<T> t246 = d23*t74; 
complex<T> t250 = d63*t87; 
complex<T> t259 = t227*t73; 
complex<T> t293 = d36*t74; 
complex<T> t298 = d76*spb51*spa61*t184 + d77*spa61*t229*t231 + s16*t184*t262 + d70*s16*t283*t285; 
complex<T> t299 = d93*s34*spb31*t254*t44 + s45*t132*t213*t50; 
complex<T> t300 = s12*t184*t262 + t128*t283*t285 + d72*spa21*t332 + d71*spa21*t44*t86; 
complex<T> t313 = d64*t77; 
complex<T> t317 = s34*t165*t284 - t132*t213*t82*T(2); 
complex<T> t326 = d74*spa32; 
complex<T> t335 = d28*spb51; 
complex<T> t337 = d75*spa32; 
complex<T> t350 = s16*(d86*s12*spb64*t184 + t128*t47*t8*t9); 
complex<T> d7 = spb62*t114*t214*t82; d7 = T(1)/d7;
complex<T> d15 = spb62*spb43*t114*t17; d15 = T(1)/d15;
complex<T> d16 = (s12 - s123)*spb63*t333; d16 = T(1)/d16;
complex<T> d19 = spb63*t333*cube(s12 - s123); d19 = T(1)/d19;
complex<T> d50 = (s34 - s345)*t369*t389; d50 = T(1)/d50;
complex<T> d53 = t369*t389*cube(s34 - s345); d53 = T(1)/d53;
complex<T> d54 = spb32*t310*t79; d54 = T(1)/d54;
complex<T> d57 = spb61*spb32*t17*t310; d57 = T(1)/d57;
complex<T> d62 = (-s123 + s45)*spb63*t333; d62 = T(1)/d62;
complex<T> d65 = spb63*t333*cube(-s123 + s45); d65 = T(1)/d65;
complex<T> d68 = s45*spb61*spb53*t333; d68 = T(1)/d68;
complex<T> t4 = -(d39*spb21*spa32*t229) - spb42*t11*t17*t293 - d44*spa32*t21*t303 + d24*t417 + d44*t115*t411*t435 + d38*spb51*t116*t485 + d41*t184*t80 + spa32*t303*t362*t80 - d42*t438*t80*t9 + spa65*(d34*spb54*t44 + t260*t336*t74 + t121*(t194 + d33*t9)) + d40*spb21*spa32*t91 + d35*spb51*t74*t91; 
complex<T> t5 = d21*t160*t184 + t168*t184 + d31*spb64*t160*t215 - t214*t215*t228 + d33*spb54*spa65*t229 - d28*t160*t229 + d32*t160*t22*t252 + spa65*t194*t260 + spb42*t11*t17*t293 - t17*t215*t335 + d22*spb54*spa65*t44 + d27*t215*t44 - d38*spb51*t116*t485 + d30*spb51*t303*t56 + d35*t284*t74 + spa65*t121*t336*t74 + t244*t384*t79 + d26*spb21*spb41*spb42*t8 + t215*t366*t8 + d34*spb54*spa65*t91 + d27*spb42*t160*t91 + d21*t215*t91 + spb51*t246*t91 + spb42*t262*t91 - d25*spb64*cube(spb41); 
complex<T> t6 = t168*t184 + d39*spb21*spa32*t229 + d48*t171*t229 + d51*spb42*t154*t260 + d44*spa32*t21*t303 - d53*t171*t217*t32 + d47*spb41*spb54*t388 - d44*t115*t411*t435 + d40*spb21*spa32*t44 + t311*t426*t50 + d52*spb31*t260*t57 + d48*spb42*t11*t17*(spb42*spa52 - t74) + d53*t192*t433*(spb42*spa52 - t74) + d49*spb54*t213*(-(spb42*spa52) + t74) - spa32*t303*t362*t80 + d42*t438*t80*t9 + spb51*t165*t91 + d50*spb51*t171*t91 + d50*t11*(-(spb42*spa52) + t74)*t91 + d41*spb42*t80*t91 - d46*spb31*cube(spb41); 
complex<T> t7 = d27*t160*t184 + d30*spb51*t160*t215 + t214*t215*t228 + d28*t160*t229 + d14*t162*t164*t24 + t219*t225*t251 - d32*t160*t22*t252 + d52*spb31*t154*t260 + t184*t262 + spb42*t273*t283 + t165*t284 + d50*t171*t284 + d53*t171*t217*t32 + d56*spa53*t16*t320 + t17*t215*t335 + d4*t164*t359 + d59*spa53*t219*t388 + d24*t417 + d61*spa53*t163*t422 + d54*spa53*t10*t44 + t244*t384*t47 + d55*spa53*t213*t50 + d31*spb64*t303*t56 + d51*spb42*t260*t57 - d11*spb61*t185*t73 - d13*t164*t305*t73 + d49*t10*t121*(-(spb42*spa52) + t74) + d48*spb42*t11*t17*(-(spb42*spa52) + t74) + d53*t192*t433*(-(spb42*spa52) + t74) + d50*t11*t44*(-(spb42*spa52) + t74) - d58*spb43*t185*t76 - spb43*t320*t338*t76 + d57*t372*t44*t76 + t253*t259*t79 + d12*spb42*t164*t73*t79 - t215*t366*t8 + t144*t388*t82 + t311*t426*t82 + t442*t47*t9 + d48*t121*t71*t9 + d6*t164*t91 + d27*t215*t91 + t216*t91 + d10*spb61*spb64*t73*t91; 
complex<T> t69 = d85*spa54*t332 + t200*t283*t361 + spb42*spb64*t127*t91 + t165*t254*square(spb41) + t125*t85*square(spb41) + d83*spa54*t86*square(spb41) + t17*t229*t235*T(2) - d81*s45*t213*t82*T(2) + s345*(d78*t229*t231 + d73*t332 + d70*t283*t361 + spb51*t165*t91 + spb42*t262*t91 + (spb42*t168 + t216)*square(spb41) + d81*t213*t82*T(2)); 
complex<T> t204 = d65*t166; 
complex<T> t279 = d79*spb51*spa43*t184 + s156*t168*t184 + d78*s156*t229*t231 + d80*spa43*t229*t231; 
complex<T> t280 = t246*t284 + spb54*t138*t44; 
complex<T> t281 = -(t326*t44*t86) + t16*t185*t337*T(2); 
complex<T> t290 = d19*t29; 
complex<T> t316 = spb51*t125*t184 + spb64*t127*t184 + t229*t231*t235 + t200*t283*t285; 
complex<T> t330 = d90*s234*s345*spb51*t184 + s234*t133*t17*t229 + d92*t137*t17*t229 + d91*t137*t417; 
complex<T> t331 = d89*t16*t185*t270 + s126*t118*t193*t8 + d87*s123*s126*t44*t86 + d88*t239*t44*t86; 
complex<T> t342 = t193*t332 + t332*t337 + s123*t216*t44 + t326*t44*t86; 
complex<T> t1 = -(d14*t162*t164*t24) + d19*t163*t266*t27 - d4*t164*t359 + d6*t164*t44 + t216*t44 + d16*t266*t44 + t253*t259*t47 + d8*spb61*t384*(-(spb64*spa62) + t70) + d11*spb61*t185*t73 + d10*spb61*t302*t73 + d13*t164*t305*t73 + d12*t283*t45*t73 + d17*spb43*t185*(spb61*spa63 - t76) + spb43*t10*t290*(spb61*spa63 - t76) + d17*spb64*t118*t78 + d18*spb43*t320*(spb61*spa63 - t76)*t78 + spb42*t273*t47*t8 + d7*t44*(-(spb61*spa62) + t72)*t9 + t442*t79*t9 + d8*(-(spb61*spa62) + t72)*t8*t9 + d5*t164*t91 + spb42*t262*t91 + d15*spb42*(-(spb64*spa62) + t70)*t91 + d9*spb61*spb64*t73*t91 + d16*spb43*(-(spb61*spa63) + t76)*t91; 
complex<T> t2 = -(d19*t163*t266*t27) + d69*t266*t44 + d63*spb61*t185*(-(spb43*spa63) + t73) + d62*spb61*t44*(-(spb43*spa63) + t73) + d17*spb43*t185*(-(spb61*spa63) + t76) + spb43*t10*t290*(-(spb61*spa63) + t76) + d62*t302*t77 + d65*spb64*t441*t77 + d17*t16*t320*t78 + d18*spb43*t320*(-(spb61*spa63) + t76)*t78 + t305*t313*(spb43*spa63 - t73)*t8 + spb61*t204*(-(spb43*spa63) + t73)*t8 - d63*t359*t77*t8 + t216*t91 + d69*spb43*(-(spb61*spa63) + t76)*t91; 
complex<T> t3 = d56*spb64*spa53*t118 - d61*spa53*t163*t422 + t216*t44 + t219*t413*t46 + t144*t388*t50 - d67*t11*t372*(spb31*spa53 + t71) + d63*spb61*t185*(spb43*spa63 - t73) + d67*t213*(-(spb42*spa52) - spb43*spa53 + t74) + d58*spb43*t185*t76 + spb43*t320*t338*t76 + t225*t251*t50*t76 - d65*spb64*t441*t77 + spb61*t204*(spb43*spa63 - t73)*t8 + t305*t313*(-(spb43*spa63) + t73)*t8 + d63*t359*t77*t8 + d55*spa53*t213*t82 + d54*spa53*t10*t91 + spb51*t165*t91 + d66*spb51*(spb31*spa53 + t71)*t91 + d62*spb61*(-(spb43*spa63) + t73)*t91 + d68*t11*(-(spb42*spa52) - spb43*spa53 + t74)*t91 + d57*t372*t76*t91 + d62*spb64*t77*t91; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t148*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t280*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + t6*(*CI_users[6]->get_value(mc,ind,mu)) + t7*(*CI_users[7]->get_value(mc,ind,mu)) + t3*(*CI_users[8]->get_value(mc,ind,mu)) + t2*(*CI_users[9]->get_value(mc,ind,mu)) + t300*(*CI_users[10]->get_value(mc,ind,mu)) + t342*(*CI_users[11]->get_value(mc,ind,mu)) + t298*(*CI_users[12]->get_value(mc,ind,mu)) + t281*(*CI_users[13]->get_value(mc,ind,mu)) + t279*(*CI_users[14]->get_value(mc,ind,mu)) + t317*(*CI_users[15]->get_value(mc,ind,mu)) + t69*(*CI_users[16]->get_value(mc,ind,mu)) + t316*(*CI_users[17]->get_value(mc,ind,mu)) + t350*(*CI_users[18]->get_value(mc,ind,mu)) + t331*(*CI_users[19]->get_value(mc,ind,mu)) + t330*(*CI_users[20]->get_value(mc,ind,mu)) + t299*(*CI_users[21]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_mmpppp_G C6g_60_G
#define _C_mmpppp_nf C6g_60_nf
#define _C_mpmppp_G C6g_58_G
#define _C_mpmppp_nf C6g_58_nf
#define _C_mppmpp_G C6g_54_G
#define _C_mppmpp_nf C6g_54_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_mmpppp_G case 60
#define _CASE_mmpppp_nf case 60
 
#define _CASE_mpmppp_G case 58
#define _CASE_mpmppp_nf case 58
 
#define _CASE_mppmpp_G case 54
#define _CASE_mppmpp_nf case 54
 
 
#define _CASE_ppmmmm_G case 3
#define _CASE_ppmmmm_nf case 3
 
#define _CASE_pmpmmm_G case 5
#define _CASE_pmpmmm_nf case 5
 
#define _CASE_pmmpmm_G case 9
#define _CASE_pmmpmm_nf case 9
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_6g_G( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_mmpppp_G: return new 
                       C6g_mmpppp_G_wCI(ind);
    _CASE_mpmppp_G: return new 
                       C6g_mpmppp_G_wCI(ind);
    _CASE_mppmpp_G: return new 
                       C6g_mppmpp_G_wCI(ind);
 
    _CASE_ppmmmm_G: return new 
                       C6g_ppmmmm_G_wCI(ind);
    _CASE_pmpmmm_G: return new 
                       C6g_pmpmmm_G_wCI(ind);
    _CASE_pmmpmm_G: return new 
                       C6g_pmmpmm_G_wCI(ind);

      default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_6g_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_mmpppp_nf: return new 
                       C6g_mmpppp_nf_wCI(ind);
    _CASE_mpmppp_nf: return new 
                       C6g_mpmppp_nf_wCI(ind);
    _CASE_mppmpp_nf: return new 
                       C6g_mppmpp_nf_wCI(ind);
    
    _CASE_ppmmmm_nf: return new 
                       C6g_ppmmmm_nf_wCI(ind);
    _CASE_pmpmmm_nf: return new 
                       C6g_pmpmmm_nf_wCI(ind);
    _CASE_pmmpmm_nf: return new 
                       C6g_pmmpmm_nf_wCI(ind);

       default: return 0;
                   }
      }
 
 
 }
 }
