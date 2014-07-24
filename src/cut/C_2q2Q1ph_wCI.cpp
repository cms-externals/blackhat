/*
*C_2q2Q1ph_wCI.cpp
*
* Created on 12/7, 2010
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
 
#define _VERBOSE 1
 
#define mH2 mc.s(ind[2-1],ind[3-1],ind[4-1],ind[5-1])
#define SPA(i,j) mc.spa(ind[i-1],ind[j-1])
#define SPB(i,j) mc.spb(ind[i-1],ind[j-1])
#define S(i,j) mc.s(ind[i-1],ind[j-1])
#define SS(i,j,k) mc.s(ind[i-1],ind[j-1],ind[k-1])

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}
 


 
class C2q2Q1ph_phqmqpQpQm_lc_wCI : public Cut_Part_wCI {
public:
       C2q2Q1ph_phqmqpQpQm_lc_wCI
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

 
class C2q2Q1ph_phqmqpQpQm_slc_wCI : public Cut_Part_wCI {
public:
       C2q2Q1ph_phqmqpQpQm_slc_wCI
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

 
class C2q2Q1ph_phqmqpQpQm_nf_wCI : public Cut_Part_wCI {
public:
       C2q2Q1ph_phqmqpQpQm_nf_wCI
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

 
class C2q2Q1ph_phqmqpQmQp_lc_wCI : public Cut_Part_wCI {
public:
       C2q2Q1ph_phqmqpQmQp_lc_wCI
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

 
class C2q2Q1ph_phqmqpQmQp_slc_wCI : public Cut_Part_wCI {
public:
       C2q2Q1ph_phqmqpQmQp_slc_wCI
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

 
class C2q2Q1ph_phqmqpQmQp_nf_wCI : public Cut_Part_wCI {
public:
       C2q2Q1ph_phqmqpQmQp_nf_wCI
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


C2q2Q1ph_phqmqpQpQm_lc_wCI::\
C2q2Q1ph_phqmqpQpQm_lc_wCI
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
	 vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
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

         vector<int> c135; c135.push_back(ind[1-1]); c135.push_back(ind[3-1]);
                            c135.push_back(ind[5-1]);    
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1ph_phqmqpQpQm_lc_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, Qp, Qm}, lc}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1ph :  phqmqpQpQm lc");
#endif
 
//#define TimeStamp "Tue 7 Dec 2010 12:39:29 on n44"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s25 = S(2,5);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> t7 = square(spa25); 
complex<T> t8 = square(spb34); 
complex<T> d1 = spa23*spa45*T(6); d1 = T(1)/d1;
complex<T> d2 = s23 - s234; d2 = T(1)/d2;
complex<T> d3 = square(s23 - s234)*T(2); d3 = T(1)/d3;
complex<T> d4 = -s345 + s45; d4 = T(1)/d4;
complex<T> d5 = square(s345 - s45)*T(2); d5 = T(1)/d5;
complex<T> d6 = spa23*spa45; d6 = T(1)/d6;
complex<T> d7 = spa45; d7 = T(1)/d7;
complex<T> d8 = spa23; d8 = T(1)/d8;
complex<T> d9 = spa23*spa45*T(2); d9 = T(1)/d9;
complex<T> d10 = spa23*T(2); d10 = T(1)/d10;
complex<T> d11 = spa45*T(2); d11 = T(1)/d11;
complex<T> t5 = d3*spa23*spa45*t8 + d2*spa25*spb34*T(2); 
complex<T> t10 = -(d6*T(2)); 
complex<T> t11 = -(d1*T(13)); 
complex<T> t12 = -(d9*mH2); 
complex<T> t16 = d5*spa23*spa45*t8 + d4*spa25*spb34*T(2); 
complex<T> t17 = t7*T(2); 
complex<T> t23 = s25*t7; 
complex<T> t25 = s34*t7; 
complex<T> t4 = t11*t7 - d3*spa23*spa45*t8 - d2*spa25*spb34*T(2); 
complex<T> t6 = t11*t7 - d5*spa23*spa45*t8 - d4*spa25*spb34*T(2); 
complex<T> t19 = t12*t23 + d9*s235*s245*t7; 
complex<T> t22 = t12*t25 + d9*s234*s345*t7; 
complex<T> t26 = t10*t7; 
complex<T> t3 = d6*s15*t17 + mH2*t26; 
complex<T> t15 = d6*s24*t17 + s15*t26; 
complex<T> co1 = -(d7*spb23*t7*T(2)); 
complex<T> co2 = s24*t26; 
complex<T> co3 = -(d8*spb45*t7*T(2)); 
complex<T> co4 = -(d10*spb45*t23); 
complex<T> co5 = -(d11*spb23*t23); 
complex<T> co6 = -(d11*spb23*t25); 
complex<T> co7 = -(d10*spb45*t25); 
SeriesC<T> result = t4*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t16*(*CI_users[2]->get_value(mc,ind,mu)) + t6*(*CI_users[3]->get_value(mc,ind,mu)) + t3*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t19*(*CI_users[9]->get_value(mc,ind,mu)) + t22*(*CI_users[10]->get_value(mc,ind,mu)) + co4*(*CI_users[11]->get_value(mc,ind,mu)) + co5*(*CI_users[12]->get_value(mc,ind,mu)) + co6*(*CI_users[13]->get_value(mc,ind,mu)) + co7*(*CI_users[14]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1ph_phqmqpQpQm_slc_wCI::\
C2q2Q1ph_phqmqpQpQm_slc_wCI
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
	 vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
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

         vector<int> c135; c135.push_back(ind[1-1]); c135.push_back(ind[3-1]);
                            c135.push_back(ind[5-1]);    
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1ph_phqmqpQpQm_slc_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, Qp, Qm}, slc}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1ph :  phqmqpQpQm slc");
#endif
 
//#define TimeStamp "Tue 7 Dec 2010 12:39:31 on n44"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s14 = S(1,4);
complex<T> s25 = S(2,5);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> s245 = SS(2,4,5);
complex<T> t7 = square(spa25); 
complex<T> t8 = square(spb34); 
complex<T> t9 = -(spa23*spa45); 
complex<T> d1 = spa23*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = s23 - s234; d2 = T(1)/d2;
complex<T> d3 = square(s23 - s234)*T(2); d3 = T(1)/d3;
complex<T> d4 = -s345 + s45; d4 = T(1)/d4;
complex<T> d5 = square(s345 - s45)*T(2); d5 = T(1)/d5;
complex<T> d6 = spa23*spa45; d6 = T(1)/d6;
complex<T> d7 = spa45; d7 = T(1)/d7;
complex<T> d8 = spa23; d8 = T(1)/d8;
complex<T> d9 = spa23*T(2); d9 = T(1)/d9;
complex<T> d10 = spa45*T(2); d10 = T(1)/d10;
complex<T> t4 = d2*spa25*spb34 + d3*t8*t9; 
complex<T> t10 = d1*T(3); 
complex<T> t12 = d8*spb45; 
complex<T> t16 = d4*spa25*spb34 + d5*t8*t9; 
complex<T> t17 = d6*t7; 
complex<T> t20 = mH2*t7; 
complex<T> t24 = d1*t7; 
complex<T> t3 = (-mH2 + s15)*t17*T(2); 
complex<T> t5 = -(d2*spa25*spb34) + t10*t7 + d3*spa23*spa45*t8; 
complex<T> t6 = -(d4*spa25*spb34) + t10*t7 + d5*spa23*spa45*t8; 
complex<T> t13 = -t17; 
complex<T> t22 = d10*spb23*t20 + s234*s235*t24; 
complex<T> t25 = d9*spb45*t20 + s245*s345*t24; 
complex<T> t29 = t12*t7; 
complex<T> t15 = s25*t13 + s14*t17; 
complex<T> t19 = s15*t13 + s24*t17; 
complex<T> t28 = s12*t17 + t29; 
complex<T> t31 = s13*t17 + t29; 
complex<T> co1 = d7*spb23*t7*T(2); 
complex<T> co2 = s24*t13; 
complex<T> co3 = s25*t17; 
SeriesC<T> result = t5*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t16*(*CI_users[2]->get_value(mc,ind,mu)) + t6*(*CI_users[3]->get_value(mc,ind,mu)) + t3*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t31*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t15*(*CI_users[9]->get_value(mc,ind,mu)) + t19*(*CI_users[10]->get_value(mc,ind,mu)) + t28*(*CI_users[11]->get_value(mc,ind,mu)) + t25*(*CI_users[12]->get_value(mc,ind,mu)) + t22*(*CI_users[13]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1ph_phqmqpQpQm_nf_wCI::\
C2q2Q1ph_phqmqpQpQm_nf_wCI
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
	 vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
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

         vector<int> c135; c135.push_back(ind[1-1]); c135.push_back(ind[3-1]);
                            c135.push_back(ind[5-1]);    
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1ph_phqmqpQpQm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, Qp, Qm}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1ph :  phqmqpQpQm nf");
#endif
 
//#define TimeStamp "Tue 7 Dec 2010 12:39:32 on n44"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa25); 
complex<T> d1 = spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = d1*t1*T(2); 
SeriesC<T> result = co1*((*CI_users[0]->get_value(mc,ind,mu)) + (*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1ph_phqmqpQmQp_lc_wCI::\
C2q2Q1ph_phqmqpQmQp_lc_wCI
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
	 vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
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

         vector<int> c135; c135.push_back(ind[1-1]); c135.push_back(ind[3-1]);
                            c135.push_back(ind[5-1]);    
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1ph_phqmqpQmQp_lc_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, Qm, Qp}, lc}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1ph :  phqmqpQmQp lc");
#endif
 
//#define TimeStamp "Tue 7 Dec 2010 12:39:40 on n44"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s14 = S(1,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> s234 = SS(2,3,4);
complex<T> t12 = square(spa24); 
complex<T> t13 = square(spa25); 
complex<T> t14 = square(spa34); 
complex<T> t15 = -(spa23*spa45); 
complex<T> t16 = square(spb35); 
complex<T> t19 = square(s25); 
complex<T> t20 = square(s34); 
complex<T> t40 = spa25*spb35; 
complex<T> d1 = spa23*spa45*T(6); d1 = T(1)/d1;
complex<T> d2 = square(spa35)*T(2); d2 = T(1)/d2;
complex<T> d3 = square(s23 - s235)*square(spa35)*T(2); d3 = T(1)/d3;
complex<T> d4 = spa35*spb25; d4 = T(1)/d4;
complex<T> d5 = (s23 - s235)*spb25; d5 = T(1)/d5;
complex<T> d6 = (-s235 + s25)*spa35; d6 = T(1)/d6;
complex<T> d7 = (s34 - s345)*spa35; d7 = T(1)/d7;
complex<T> d8 = square(s345 - s45)*square(spa35)*T(2); d8 = T(1)/d8;
complex<T> d9 = spa35*spb34; d9 = T(1)/d9;
complex<T> d10 = (-s345 + s45)*spb34; d10 = T(1)/d10;
complex<T> d11 = spa23*spa45; d11 = T(1)/d11;
complex<T> d12 = spa45; d12 = T(1)/d12;
complex<T> d13 = spa45*square(spa35); d13 = T(1)/d13;
complex<T> d14 = spa23*spa45*square(spa35); d14 = T(1)/d14;
complex<T> d15 = spa23; d15 = T(1)/d15;
complex<T> d16 = spa23*square(spa35); d16 = T(1)/d16;
complex<T> d17 = spa23*spa45*T(2); d17 = T(1)/d17;
complex<T> d18 = spa23*T(2); d18 = T(1)/d18;
complex<T> d19 = spa45*square(spa35)*T(2); d19 = T(1)/d19;
complex<T> d20 = spa45*T(2); d20 = T(1)/d20;
complex<T> d21 = spa23*square(spa35)*T(2); d21 = T(1)/d21;
complex<T> t17 = -(d11*T(2)); 
complex<T> t22 = -(d1*T(13)); 
complex<T> t23 = -(d17*mH2); 
complex<T> t26 = d15*spb45; 
complex<T> t30 = d2*spa23; 
complex<T> t36 = d11*t12; 
complex<T> t42 = d17*t12; 
complex<T> t45 = t13*t14; 
complex<T> t49 = d6*spa34; 
complex<T> t50 = s25*t12; 
complex<T> t53 = d7*spa34; 
complex<T> t54 = s34*t12; 
complex<T> t7 = d8*t15*t20 + t12*t22 + spa45*t30 - d10*spa25*t16*T(2) + d9*t40*T(2); 
complex<T> t8 = d2*t15 + d8*spa23*spa45*t20 - t40*t53 + d10*spa25*t16*T(2) - d9*t40*T(2); 
complex<T> t9 = d3*t15*t19 + t12*t22 + spa45*t30 + d4*spa34*spb35*T(2) - d5*spa34*t16*T(2); 
complex<T> t10 = d2*t15 + d3*spa23*spa45*t19 - t40*t49 - d4*spa34*spb35*T(2) + d5*spa34*t16*T(2); 
complex<T> t27 = -t45; 
complex<T> t28 = -t36; 
complex<T> t39 = s235*s245*t42 + t23*t50; 
complex<T> t44 = s234*s345*t42 + t23*t54; 
complex<T> t46 = t12*t17; 
complex<T> t57 = t12*t26; 
complex<T> t4 = mH2*t46 + s15*t36*T(2); 
complex<T> t5 = d14*s14*t27 + s25*t28 + s14*t36 + d14*s25*t45; 
complex<T> t33 = s15*t46 + s24*t36*T(2); 
complex<T> t34 = s34*(t28 + d14*t45); 
complex<T> t58 = spb23*t27; 
complex<T> t59 = spb45*t27; 
complex<T> t6 = -(d12*spb23*t12) + d13*t58; 
complex<T> t11 = d14*s12*t27 + s12*t36 + t57 + d16*t59; 
complex<T> co1 = t40*t49; 
complex<T> co2 = t40*t53; 
complex<T> co3 = s24*t46; 
complex<T> co4 = -(t57*T(2)); 
complex<T> co5 = -(d18*spb45*t50); 
complex<T> co6 = d19*s25*t58; 
complex<T> co7 = -(d20*spb23*t54); 
complex<T> co8 = d21*s34*t59; 
SeriesC<T> result = t9*(*CI_users[0]->get_value(mc,ind,mu)) + t10*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + t7*(*CI_users[5]->get_value(mc,ind,mu)) + t4*(*CI_users[6]->get_value(mc,ind,mu)) + t6*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t5*(*CI_users[9]->get_value(mc,ind,mu)) + t33*(*CI_users[10]->get_value(mc,ind,mu)) + t34*(*CI_users[11]->get_value(mc,ind,mu)) + t11*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)) + t39*(*CI_users[14]->get_value(mc,ind,mu)) + t44*(*CI_users[15]->get_value(mc,ind,mu)) + co5*(*CI_users[16]->get_value(mc,ind,mu)) + co6*(*CI_users[17]->get_value(mc,ind,mu)) + co7*(*CI_users[18]->get_value(mc,ind,mu)) + co8*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1ph_phqmqpQmQp_slc_wCI::\
C2q2Q1ph_phqmqpQmQp_slc_wCI
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
	 vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
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

         vector<int> c135; c135.push_back(ind[1-1]); c135.push_back(ind[3-1]);
                            c135.push_back(ind[5-1]);    
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c23, c4));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1ph_phqmqpQmQp_slc_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, Qm, Qp}, slc}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1ph :  phqmqpQmQp slc");
#endif
 
//#define TimeStamp "Tue 7 Dec 2010 12:39:43 on n44"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s14 = S(1,4);
complex<T> s25 = S(2,5);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s235 = SS(2,3,5);
complex<T> s234 = SS(2,3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> s245 = SS(2,4,5);
complex<T> t7 = square(spa24); 
complex<T> t8 = square(spb35); 
complex<T> t9 = -(spa23*spa45); 
complex<T> d1 = spa23*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = s23 - s235; d2 = T(1)/d2;
complex<T> d3 = square(s23 - s235)*T(2); d3 = T(1)/d3;
complex<T> d4 = -s345 + s45; d4 = T(1)/d4;
complex<T> d5 = square(s345 - s45)*T(2); d5 = T(1)/d5;
complex<T> d6 = spa23*spa45; d6 = T(1)/d6;
complex<T> d7 = spa45; d7 = T(1)/d7;
complex<T> d8 = spa23; d8 = T(1)/d8;
complex<T> d9 = spa23*T(2); d9 = T(1)/d9;
complex<T> d10 = spa45*T(2); d10 = T(1)/d10;
complex<T> t4 = -(d2*spa24*spb35) + d3*t8*t9; 
complex<T> t10 = d1*T(3); 
complex<T> t12 = d8*spb45; 
complex<T> t16 = -(d4*spa24*spb35) + d5*t8*t9; 
complex<T> t17 = d6*t7; 
complex<T> t20 = mH2*t7; 
complex<T> t24 = d1*t7; 
complex<T> t3 = (-mH2 + s15)*t17*T(2); 
complex<T> t5 = d2*spa24*spb35 + t10*t7 + d3*spa23*spa45*t8; 
complex<T> t6 = d4*spa24*spb35 + t10*t7 + d5*spa23*spa45*t8; 
complex<T> t13 = -t17; 
complex<T> t22 = d10*spb23*t20 + s234*s235*t24; 
complex<T> t25 = d9*spb45*t20 + s245*s345*t24; 
complex<T> t29 = t12*t7; 
complex<T> t15 = s25*t13 + s14*t17; 
complex<T> t19 = s15*t13 + s24*t17; 
complex<T> t28 = s12*t17 + t29; 
complex<T> t31 = s13*t17 + t29; 
complex<T> co1 = d7*spb23*t7*T(2); 
complex<T> co2 = s24*t13; 
complex<T> co3 = s25*t17; 
SeriesC<T> result = t5*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t16*(*CI_users[2]->get_value(mc,ind,mu)) + t6*(*CI_users[3]->get_value(mc,ind,mu)) + t3*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t31*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + t15*(*CI_users[9]->get_value(mc,ind,mu)) + t19*(*CI_users[10]->get_value(mc,ind,mu)) + t28*(*CI_users[11]->get_value(mc,ind,mu)) + t25*(*CI_users[12]->get_value(mc,ind,mu)) + t22*(*CI_users[13]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1ph_phqmqpQmQp_nf_wCI::\
C2q2Q1ph_phqmqpQmQp_nf_wCI
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
	 vector<int> c24;  c24.push_back(ind[2-1]); c24.push_back(ind[4-1]);
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

         vector<int> c135; c135.push_back(ind[1-1]); c135.push_back(ind[3-1]);
                            c135.push_back(ind[5-1]);    
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1ph_phqmqpQmQp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, Qm, Qp}, nf}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1ph :  phqmqpQmQp nf");
#endif
 
//#define TimeStamp "Tue 7 Dec 2010 12:39:43 on n44"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa24); 
complex<T> d1 = spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = d1*t1*T(2); 
SeriesC<T> result = co1*((*CI_users[0]->get_value(mc,ind,mu)) + (*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_phqmqpQpQm_lc C2q2Q1ph_6356_lc
#define _C_phqmqpQpQm_slc C2q2Q1ph_6356_slc
#define _C_phqmqpQpQm_nf C2q2Q1ph_6356_nf
#define _C_phqmqpQmQp_lc C2q2Q1ph_7436_lc
#define _C_phqmqpQmQp_slc C2q2Q1ph_7436_slc
#define _C_phqmqpQmQp_nf C2q2Q1ph_7436_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_phqmqpQpQm_lc case 6356
 
#define _CASE_phqmqpQpQm_slc case 6356
 
#define _CASE_phqmqpQpQm_nf case 6356
 
#define _CASE_phqmqpQmQp_lc case 7436
 
#define _CASE_phqmqpQmQp_slc case 7436
 
#define _CASE_phqmqpQmQp_nf case 7436
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2Q1ph_lc( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_phqmqpQpQm_lc: return new 
                       C2q2Q1ph_phqmqpQpQm_lc_wCI(ind);
    _CASE_phqmqpQmQp_lc: return new 
                       C2q2Q1ph_phqmqpQmQp_lc_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q1ph_nf( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_phqmqpQpQm_nf: return new 
                       C2q2Q1ph_phqmqpQpQm_nf_wCI(ind);
    _CASE_phqmqpQmQp_nf: return new 
                       C2q2Q1ph_phqmqpQmQp_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q1ph_slc( int hc,const std::vector<int>& ind) 
{ 
    switch (hc) {
    _CASE_phqmqpQpQm_slc: return new 
                       C2q2Q1ph_phqmqpQpQm_slc_wCI(ind);
    _CASE_phqmqpQmQp_slc: return new 
                       C2q2Q1ph_phqmqpQmQp_slc_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
