/*
*C_2q2Q_wCI.cpp
*
* Created on 11/15, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*
* convention: first fermion is always the antifermion
*            
*
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
 



 
class C2q2Q_qmQpQmqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmQpQmqp_LLT_wCI
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

 
class C2q2Q_qmQmQpqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmQmQpqp_LLT_wCI
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

 
class C2q2Q_qmqpQpQm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmqpQpQm_LLT_wCI
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

 
class C2q2Q_qmqpQmQp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmqpQmQp_LLT_wCI
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

 
class C2q2Q_qmQpQmqp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmQpQmqp_LRT_wCI
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

 
class C2q2Q_qmQmQpqp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmQmQpqp_LRT_wCI
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

 
class C2q2Q_qmqpQpQm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmqpQpQm_LRT_wCI
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

 
class C2q2Q_qmqpQmQp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmqpQmQp_LRT_wCI
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

 
class C2q2Q_qmqpQpQm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmqpQpQm_nfLLT_wCI
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

 
class C2q2Q_qmqpQmQp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmqpQmQp_nfLLT_wCI
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

 
class C2q2Q_qmqpQpQm_nfLRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmqpQpQm_nfLRT_wCI
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

 
class C2q2Q_qmqpQmQp_nfLRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q_qmqpQmQp_nfLRT_wCI
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


C2q2Q_qmQpQmqp_LLT_wCI::\
C2q2Q_qmQpQmqp_LLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q_qmQpQmqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmQpQmqp LLT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:12 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb14 = SPB(1,4);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa14*spa23*T(2); d1 = T(1)/d1;
complex<T> d2 = spa23; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d2*spb14*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q_qmQmQpqp_LLT_wCI::\
C2q2Q_qmQmQpqp_LLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q_qmQmQpqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmQmQpqp LLT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:13 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb14 = SPB(1,4);
complex<T> t1 = square(spa12); 
complex<T> d1 = spa14*spa23*T(2); d1 = T(1)/d1;
complex<T> d2 = spa23; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d2*spb14*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q_qmqpQpQm_LLT_wCI::\
C2q2Q_qmqpQpQm_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q_qmqpQpQm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQpQm LLT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:13 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> s23 = S(2,3);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(2)); 
complex<T> co2 = -(d2*spb12*t1*T(2)); 
complex<T> co3 = -(d2*s23*spb12*t1); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)) + co3*(*CI_users[2]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q_qmqpQmQp_LLT_wCI::\
C2q2Q_qmqpQmQp_LLT_wCI
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
     C2q2Q_qmqpQmQp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQmQp LLT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:14 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb12 = SPB(1,2);
complex<T> s23 = S(2,3);
complex<T> t4 = square(spa13); 
complex<T> t7 = cube(s23); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> d2 = spa34*spb13; d2 = T(1)/d2;
complex<T> d3 = spa34; d3 = T(1)/d3;
complex<T> d4 = spa34*square(spb13); d4 = T(1)/d4;
complex<T> d5 = spa12*spa34; d5 = T(1)/d5;
complex<T> d6 = spa12*spa34*square(spb13); d6 = T(1)/d6;
complex<T> d7 = spa34*T(2); d7 = T(1)/d7;
complex<T> d8 = spa34*square(spb13)*T(2); d8 = T(1)/d8;
complex<T> t1 = -(d5*s23*t4) + d6*t7; 
complex<T> t3 = -(spb12*(d3*t4 + d4*square(s23))); 
complex<T> t8 = d2*spa13; 
complex<T> t9 = -(spb12*(d7*s23*t4 + d8*t7)); 
complex<T> t11 = spb12*t8; 
complex<T> t2 = t11 - d1*t4*T(2); 
complex<T> co1 = -t11; 
complex<T> co2 = Complex(0,1); 
SeriesC<T> result = co2*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t9*(*CI_users[4]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  


C2q2Q_qmQpQmqp_LRT_wCI::\
C2q2Q_qmQpQmqp_LRT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2Q_qmQpQmqp_LRT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmQpQmqp LRT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2Q_qmQmQpqp_LRT_wCI::\
C2q2Q_qmQmQpqp_LRT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2Q_qmQmQpqp_LRT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmQmQpqp LRT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2Q_qmqpQpQm_LRT_wCI::\
C2q2Q_qmqpQpQm_LRT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q_qmqpQpQm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQpQm LRT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:14 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb12*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q_qmqpQmQp_LRT_wCI::\
C2q2Q_qmqpQmQp_LRT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q_qmqpQmQp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQmQp LRT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:15 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb12*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q_qmqpQpQm_nfLLT_wCI::\
C2q2Q_qmqpQpQm_nfLLT_wCI
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
     C2q2Q_qmqpQpQm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQpQm nfLLT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:15 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q_qmqpQmQp_nfLLT_wCI::\
C2q2Q_qmqpQmQp_nfLLT_wCI
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
     C2q2Q_qmqpQmQp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQmQp nfLLT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:16 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q_qmqpQpQm_nfLRT_wCI::\
C2q2Q_qmqpQpQm_nfLRT_wCI
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
     C2q2Q_qmqpQpQm_nfLRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm}, nfLRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQpQm nfLRT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:16 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q_qmqpQmQp_nfLRT_wCI::\
C2q2Q_qmqpQmQp_nfLRT_wCI
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
     C2q2Q_qmqpQmQp_nfLRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp}, nfLRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q :  qmqpQmQp nfLRT");
#endif
 
//#define TimeStamp "Mon 15 Nov 2010 09:04:16 on n2174"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qmQpQmqp_LLT C2q2Q_607_LLT
#define _C_qmQmQpqp_LLT C2q2Q_637_LLT
#define _C_qmqpQpQm_LLT C2q2Q_1057_LLT
#define _C_qmqpQmQp_LLT C2q2Q_1237_LLT
#define _C_qmQpQmqp_LRT C2q2Q_607_LRT
#define _C_qmQmQpqp_LRT C2q2Q_637_LRT
#define _C_qmqpQpQm_LRT C2q2Q_1057_LRT
#define _C_qmqpQmQp_LRT C2q2Q_1237_LRT
#define _C_qmqpQpQm_nfLLT C2q2Q_1057_nfLLT
#define _C_qmqpQmQp_nfLLT C2q2Q_1237_nfLLT
#define _C_qmqpQpQm_nfLRT C2q2Q_1057_nfLRT
#define _C_qmqpQmQp_nfLRT C2q2Q_1237_nfLRT
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmQpQmqp_LLT case 607
 
#define _CASE_qmQmQpqp_LLT case 637
 
#define _CASE_qmqpQpQm_LLT case 1057
 
#define _CASE_qmqpQmQp_LLT case 1237
 
#define _CASE_qmQpQmqp_LRT case 607
 
#define _CASE_qmQmQpqp_LRT case 637
 
#define _CASE_qmqpQpQm_LRT case 1057
 
#define _CASE_qmqpQmQp_LRT case 1237
 
#define _CASE_qmqpQpQm_nfLLT case 1057
 
#define _CASE_qmqpQmQp_nfLLT case 1237
 
#define _CASE_qmqpQpQm_nfLRT case 1057
 
#define _CASE_qmqpQmQp_nfLRT case 1237
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2Q_LLT( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmQpQmqp_LLT: return new 
                       C2q2Q_qmQpQmqp_LLT_wCI(ind);
    _CASE_qmQmQpqp_LLT: return new 
                       C2q2Q_qmQmQpqp_LLT_wCI(ind);
    _CASE_qmqpQpQm_LLT: return new 
                       C2q2Q_qmqpQpQm_LLT_wCI(ind);
    _CASE_qmqpQmQp_LLT: return new 
                       C2q2Q_qmqpQmQp_LLT_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q_LRT( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qmQpQmqp_LRT: return new 
                       C2q2Q_qmQpQmqp_LRT_wCI(ind);
    _CASE_qmQmQpqp_LRT: return new 
                       C2q2Q_qmQmQpqp_LRT_wCI(ind);
    _CASE_qmqpQpQm_LRT: return new 
                       C2q2Q_qmqpQpQm_LRT_wCI(ind);
    _CASE_qmqpQmQp_LRT: return new 
                       C2q2Q_qmqpQmQp_LRT_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q_nfLLT( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_qmqpQpQm_nfLLT: return new 
                       C2q2Q_qmqpQpQm_nfLLT_wCI(ind);
    _CASE_qmqpQmQp_nfLLT: return new 
                       C2q2Q_qmqpQmQp_nfLLT_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q_nfLRT( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_qmqpQpQm_nfLRT: return new 
                       C2q2Q_qmqpQpQm_nfLRT_wCI(ind);
    _CASE_qmqpQmQp_nfLRT: return new 
                       C2q2Q_qmqpQmQp_nfLRT_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
