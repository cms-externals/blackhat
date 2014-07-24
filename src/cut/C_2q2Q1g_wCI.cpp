/*
*C_2q2Q1g_wCI.cpp
*
* Created on 12/22, 2010
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
 
 
class C2q2Q1g_qmqpmQmQp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpmQmQp_LLT_wCI
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

 
class C2q2Q1g_qmqppQmQp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqppQmQp_LLT_wCI
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

 
class C2q2Q1g_qmqpmQpQm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpmQpQm_LLT_wCI
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

 
class C2q2Q1g_qmqppQpQm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqppQpQm_LLT_wCI
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

 
class C2q2Q1g_qmqpQmmQp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmmQp_LLT_wCI
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

 
class C2q2Q1g_qmqpQmpQp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmpQp_LLT_wCI
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

 
class C2q2Q1g_qmqpQpmQm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpmQm_LLT_wCI
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

 
class C2q2Q1g_qmqpQppQm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQppQm_LLT_wCI
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

 
class C2q2Q1g_qmqpQmQpm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmQpm_LLT_wCI
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

 
class C2q2Q1g_qmqpQmQpp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmQpp_LLT_wCI
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

 
class C2q2Q1g_qmqpQpQmm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpQmm_LLT_wCI
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

 
class C2q2Q1g_qmqpQpQmp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpQmp_LLT_wCI
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

 
class C2q2Q1g_qmmqpQmQp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmmqpQmQp_LLT_wCI
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

 
class C2q2Q1g_qmpqpQmQp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmpqpQmQp_LLT_wCI
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

 
class C2q2Q1g_qmmqpQpQm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmmqpQpQm_LLT_wCI
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

 
class C2q2Q1g_qmpqpQpQm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmpqpQpQm_LLT_wCI
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

 
class C2q2Q1g_qmqpmQmQp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpmQmQp_LRT_wCI
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

 
class C2q2Q1g_qmqppQmQp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqppQmQp_LRT_wCI
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

 
class C2q2Q1g_qmqpmQpQm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpmQpQm_LRT_wCI
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

 
class C2q2Q1g_qmqppQpQm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqppQpQm_LRT_wCI
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

 
class C2q2Q1g_qmqpQmmQp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmmQp_LRT_wCI
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

 
class C2q2Q1g_qmqpQmpQp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmpQp_LRT_wCI
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

 
class C2q2Q1g_qmqpQpmQm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpmQm_LRT_wCI
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

 
class C2q2Q1g_qmqpQppQm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQppQm_LRT_wCI
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

 
class C2q2Q1g_qmqpQmQpm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmQpm_LRT_wCI
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

 
class C2q2Q1g_qmqpQmQpp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmQpp_LRT_wCI
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

 
class C2q2Q1g_qmqpQpQmm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpQmm_LRT_wCI
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

 
class C2q2Q1g_qmqpQpQmp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpQmp_LRT_wCI
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

 
class C2q2Q1g_qmmqpQmQp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmmqpQmQp_LRT_wCI
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

 
class C2q2Q1g_qmpqpQmQp_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmpqpQmQp_LRT_wCI
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

 
class C2q2Q1g_qmmqpQpQm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmmqpQpQm_LRT_wCI
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

 
class C2q2Q1g_qmpqpQpQm_LRT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmpqpQpQm_LRT_wCI
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

 
class C2q2Q1g_qmmQmQpqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmmQmQpqp_LLT_wCI
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

 
class C2q2Q1g_qmmQpQmqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmmQpQmqp_LLT_wCI
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

 
class C2q2Q1g_qmpQmQpqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmpQmQpqp_LLT_wCI
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

 
class C2q2Q1g_qmpQpQmqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmpQpQmqp_LLT_wCI
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

 
class C2q2Q1g_qmQmmQpqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQmmQpqp_LLT_wCI
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

 
class C2q2Q1g_qmQpmQmqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQpmQmqp_LLT_wCI
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

 
class C2q2Q1g_qmQmpQpqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQmpQpqp_LLT_wCI
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

 
class C2q2Q1g_qmQppQmqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQppQmqp_LLT_wCI
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

 
class C2q2Q1g_qmQmQpqpm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQmQpqpm_LLT_wCI
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

 
class C2q2Q1g_qmQmQpqpp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQmQpqpp_LLT_wCI
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

 
class C2q2Q1g_qmQpQmqpm_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQpQmqpm_LLT_wCI
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

 
class C2q2Q1g_qmQpQmqpp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQpQmqpp_LLT_wCI
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

 
class C2q2Q1g_qmQmQpmqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQmQpmqp_LLT_wCI
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

 
class C2q2Q1g_qmQmQppqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQmQppqp_LLT_wCI
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

 
class C2q2Q1g_qmQpQmmqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQpQmmqp_LLT_wCI
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

 
class C2q2Q1g_qmQpQmpqp_LLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmQpQmpqp_LLT_wCI
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

 
class C2q2Q1g_qmqpmQmQp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpmQmQp_nfLLT_wCI
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

 
class C2q2Q1g_qmqppQmQp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqppQmQp_nfLLT_wCI
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

 
class C2q2Q1g_qmqpmQpQm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpmQpQm_nfLLT_wCI
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

 
class C2q2Q1g_qmqppQpQm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqppQpQm_nfLLT_wCI
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

 
class C2q2Q1g_qmqpQmmQp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmmQp_nfLLT_wCI
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

 
class C2q2Q1g_qmqpQmpQp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmpQp_nfLLT_wCI
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

 
class C2q2Q1g_qmqpQpmQm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpmQm_nfLLT_wCI
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

 
class C2q2Q1g_qmqpQppQm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQppQm_nfLLT_wCI
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

 
class C2q2Q1g_qmqpQmQpm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmQpm_nfLLT_wCI
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

 
class C2q2Q1g_qmqpQmQpp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQmQpp_nfLLT_wCI
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

 
class C2q2Q1g_qmqpQpQmm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpQmm_nfLLT_wCI
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

 
class C2q2Q1g_qmqpQpQmp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmqpQpQmp_nfLLT_wCI
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

 
class C2q2Q1g_qmmqpQmQp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmmqpQmQp_nfLLT_wCI
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

 
class C2q2Q1g_qmpqpQmQp_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmpqpQmQp_nfLLT_wCI
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

 
class C2q2Q1g_qmmqpQpQm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmmqpQpQm_nfLLT_wCI
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

 
class C2q2Q1g_qmpqpQpQm_nfLLT_wCI : public Cut_Part_wCI {
public:
       C2q2Q1g_qmpqpQpQm_nfLLT_wCI
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


C2q2Q1g_qmqpmQmQp_LLT_wCI::\
C2q2Q1g_qmqpmQmQp_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpmQmQp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpmQmQp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:28:04 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa35 = SPA(3,5);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = -(spa14*spb14);
complex<T> s13 = -(spa13*spb13);
complex<T> s12 = -(spa12*spb12);
complex<T> t1 = square(s12 - s45); 
complex<T> t2 = spb12*T(2); 
complex<T> t3 = spb34*spb45; 
complex<T> t5 = square(s23 - s45); 
complex<T> t16 = square(spb25); 
complex<T> t17 = square(spa13); 
complex<T> t18 = square(spb15); 
complex<T> t19 = square(spb23); 
complex<T> t20 = square(spa14); 
complex<T> t21 = square(spb24); 
complex<T> t25 = s15 - s23 + s45; 
complex<T> t56 = s14*spa14; 
complex<T> t67 = s13*spb12; 
complex<T> t70 = spb24*spb25; 
complex<T> t90 = spb14*spb23; 
complex<T> t93 = spa13*spb25; 
complex<T> d3 = cube(s12 - s45)*T(3); d3 = T(1)/d3;
complex<T> d5 = (s12 - s45)*spb23*spb45*T(2); d5 = T(1)/d5;
complex<T> d13 = (s15 - s23)*spb23*spb34; d13 = T(1)/d13;
complex<T> d15 = (s23 - s45)*spb23*spb34; d15 = T(1)/d15;
complex<T> d22 = spb34*spb35; d22 = T(1)/d22;
complex<T> d28 = spb12*spb35; d28 = T(1)/d28;
complex<T> d29 = spb12*spb34*cube(spb13); d29 = T(1)/d29;
complex<T> d30 = spb12*spb23*spb34; d30 = T(1)/d30;
complex<T> d31 = spb12*spb34*spb35; d31 = T(1)/d31;
complex<T> t6 = spb13*t3; 
complex<T> t35 = d3*spa34; 
complex<T> t39 = -t56; 
complex<T> t52 = spb24*t16; 
complex<T> t53 = spb14*t17; 
complex<T> t54 = t18*t19; 
complex<T> t66 = spb35*t17; 
complex<T> t71 = spb15*t20; 
complex<T> t72 = spb23*t3; 
complex<T> t88 = s45*t16; 
complex<T> t92 = spb24*t18; 
complex<T> d1 = spb12*spb34*t1; d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spb34*t2; d2 = T(1)/d2;
complex<T> d4 = spb45*t1*T(2); d4 = T(1)/d4;
complex<T> d8 = t1*t2*t3; d8 = T(1)/d8;
complex<T> d11 = spb23*spb45*t1*T(3); d11 = T(1)/d11;
complex<T> d12 = (s15 - s23)*spb23*spb34*t25; d12 = T(1)/d12;
complex<T> d14 = (s23 - s45)*spb23*spb34*t25; d14 = T(1)/d14;
complex<T> d18 = t2*t3*t5; d18 = T(1)/d18;
complex<T> d23 = t3*cube(spb13); d23 = T(1)/d23;
complex<T> d24 = spb23*spb34*square(t25); d24 = T(1)/d24;
complex<T> d25 = spb34*square(t25); d25 = T(1)/d25;
complex<T> d26 = spb12*t3*cube(spb13); d26 = T(1)/d26;
complex<T> d27 = spb12*t3; d27 = T(1)/d27;
complex<T> d32 = t3*T(2); d32 = T(1)/d32;
complex<T> d33 = spb23*spb34*square(t25)*T(2); d33 = T(1)/d33;
complex<T> d34 = spb45*t2; d34 = T(1)/d34;
complex<T> d35 = t3*cube(spb13)*T(2); d35 = T(1)/d35;
complex<T> d36 = spb23*t2; d36 = T(1)/d36;
complex<T> d37 = spb23*spb34*t2; d37 = T(1)/d37;
complex<T> d39 = spb35*t2; d39 = T(1)/d39;
complex<T> t10 = spa12*(d22*t16 + d23*spb14*t54); 
complex<T> t30 = d24*s15; 
complex<T> t32 = d33*s45; 
complex<T> t37 = -t52; 
complex<T> t38 = -t54; 
complex<T> t47 = d32*spa23; 
complex<T> t64 = -t71; 
complex<T> t99 = spa45*t52; 
complex<T> d6 = t1*t2*t6; d6 = T(1)/d6;
complex<T> d7 = (s12 - s45)*t6*t67; d7 = T(1)/d7;
complex<T> d9 = spb12*t72*T(6); d9 = T(1)/d9;
complex<T> d10 = (s12 - s45)*spb12*t72; d10 = T(1)/d10;
complex<T> d16 = t2*t5*t6; d16 = T(1)/d16;
complex<T> d17 = (s23 - s45)*t6*t67; d17 = T(1)/d17;
complex<T> d19 = (s23 - s45)*t72; d19 = T(1)/d19;
complex<T> d20 = (s23 - s45)*spb12*t72; d20 = T(1)/d20;
complex<T> d21 = t2*t72; d21 = T(1)/d21;
complex<T> d38 = t72*T(2); d38 = T(1)/d38;
complex<T> t9 = d12*t21*t64 + d13*spa14*t70 + d12*t56*t70; 
complex<T> t11 = d19*spa12*t37 + d20*s13*t52 + d17*t38*t53 + d16*t53*t54 + d14*t21*t64 - d13*spa14*t70 + d15*spa14*t70 + d12*t39*t70 + d14*t56*t70 + d12*t21*t71 - d18*spb23*t17*t92; 
complex<T> t13 = d29*spa45*spb14*t38 + d24*s45*(t39*t70 + t21*t71) + d31*t88 + d30*t99*T(2); 
complex<T> t15 = d5*spa13*t16 + d4*spb15*spb25*t17 + d9*t52 + d10*s13*t52 + d8*s23*spa23*t52 + d7*t38*t53 + d6*t53*t54 + d11*spb12*spb25*t66 + t2*t35*t66 - d1*spa34*t90*t93 - d2*spa35*t16*T(5); 
complex<T> t50 = t30*(t39*t70 + t21*t71); 
complex<T> t61 = spa12*t47; 
complex<T> t62 = s15*t32*(t39*t70 + t21*t71); 
complex<T> t95 = spa23*t37; 
complex<T> t12 = d26*s23*spb14*t38 + d25*spa23*(t39*t70 + t21*t71) + d27*t95; 
complex<T> t14 = -(d5*spa13*t16) - d4*spb15*spb25*t17 + d10*s13*t37 + d20*s13*t37 + d21*t52 + d19*spa12*t52 + d16*t38*t53 + d6*t38*t53 + d17*t53*t54 + d7*t53*t54 - d11*spb12*spb25*t66 - d15*spa14*t70 + d14*t39*t70 + d14*t21*t71 + d18*spb23*t17*t92 + d1*spa34*t90*t93 + d8*s23*t95 - spb12*t35*t66*T(2) + d2*spa35*t16*T(5); 
complex<T> t51 = d35*s23*spa12*spb14*t54 + t52*t61; 
complex<T> co1 = -(d28*spa34*t16); 
complex<T> co2 = t37*t61; 
complex<T> co3 = d34*spa34*t95; 
complex<T> co4 = d36*spa34*spa45*t37; 
complex<T> co5 = d37*s15*t99; 
complex<T> co6 = d38*s15*spa12*t52; 
complex<T> co7 = -(d39*spa34*t88); 
complex<T> co8 = Complex(0,1); 
SeriesC<T> result = co8*(t15*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t11*(*CI_users[2]->get_value(mc,ind,mu)) + t14*(*CI_users[4]->get_value(mc,ind,mu)) + t10*(*CI_users[5]->get_value(mc,ind,mu)) + t50*(*CI_users[6]->get_value(mc,ind,mu)) + t12*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t13*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t62*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + t51*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[14]->get_value(mc,ind,mu)) + co5*(*CI_users[15]->get_value(mc,ind,mu)) + co6*(*CI_users[16]->get_value(mc,ind,mu)) + co7*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqppQmQp_LLT_wCI::\
C2q2Q1g_qmqppQmQp_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqppQmQp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqppQmQp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:28:56 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s25 = -(spa25*spb25);
complex<T> s15 = S(1,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = square(s12 - s45); 
complex<T> t2 = spa45*T(2); 
complex<T> t3 = spa12*spa23; 
complex<T> t5 = square(s12 - s34); 
complex<T> t17 = square(spa14); 
complex<T> t18 = square(spb35); 
complex<T> t19 = square(spa15); 
complex<T> t20 = square(spa34); 
complex<T> t21 = square(spa24); 
complex<T> t22 = square(spb25); 
complex<T> t26 = s12 + s15 - s34; 
complex<T> t58 = s25*spb25; 
complex<T> t70 = s34*spa25; 
complex<T> t71 = s35*spa45; 
complex<T> t72 = spa14*spa24; 
complex<T> t91 = spa14*spa25; 
complex<T> t95 = spa34*spb35; 
complex<T> d5 = (s12 - s34)*spa23*spa34; d5 = T(1)/d5;
complex<T> d8 = (s12 - s45)*spa12*spa34*T(2); d8 = T(1)/d8;
complex<T> d17 = cube(s12 - s45)*T(3); d17 = T(1)/d17;
complex<T> d19 = (s15 - s34)*spa23*spa34; d19 = T(1)/d19;
complex<T> d22 = spa13*spa23*spa45; d22 = T(1)/d22;
complex<T> d23 = spa23*spa34*spa45; d23 = T(1)/d23;
complex<T> d24 = spa23*spa45*cube(spa35); d24 = T(1)/d24;
complex<T> d26 = spa13*spa45; d26 = T(1)/d26;
complex<T> d30 = spa13*spa23; d30 = T(1)/d30;
complex<T> t6 = spa35*t3; 
complex<T> t36 = d17*spb23; 
complex<T> t41 = -t58; 
complex<T> t56 = spa25*t18; 
complex<T> t57 = t19*t20; 
complex<T> t65 = spa24*t17; 
complex<T> t68 = spa13*t18; 
complex<T> t73 = spa15*t22; 
complex<T> t74 = spa34*t3; 
complex<T> t93 = spa24*t19; 
complex<T> d4 = (s12 - s45)*spa23*t2; d4 = T(1)/d4;
complex<T> d6 = (s12 - s34)*spa23*spa34*t26; d6 = T(1)/d6;
complex<T> d7 = t1*t2*t3; d7 = T(1)/d7;
complex<T> d9 = spa23*spa45*t1; d9 = T(1)/d9;
complex<T> d10 = spa12*t1*T(2); d10 = T(1)/d10;
complex<T> d11 = t2*t3*t5; d11 = T(1)/d11;
complex<T> d16 = spa12*spa34*t1*T(3); d16 = T(1)/d16;
complex<T> d20 = (s15 - s34)*spa23*spa34*t26; d20 = T(1)/d20;
complex<T> d25 = spa23*spa34*square(t26); d25 = T(1)/d25;
complex<T> d27 = spa45*t3*cube(spa35); d27 = T(1)/d27;
complex<T> d28 = spa45*t3; d28 = T(1)/d28;
complex<T> d29 = spa23*square(t26); d29 = T(1)/d29;
complex<T> d31 = t3*cube(spa35); d31 = T(1)/d31;
complex<T> d32 = spa13*t2; d32 = T(1)/d32;
complex<T> d33 = spa34*t2; d33 = T(1)/d33;
complex<T> d34 = spa12*t2; d34 = T(1)/d34;
complex<T> d35 = t3*cube(spa35)*T(2); d35 = T(1)/d35;
complex<T> d37 = spa23*spa34*t2; d37 = T(1)/d37;
complex<T> d38 = spa23*spa34*square(t26)*T(2); d38 = T(1)/d38;
complex<T> t9 = spb23*(d32*s12*t17 + d33*spb12*t65); 
complex<T> t10 = -(d19*spb25*t72) + d20*t41*t72 + d20*t21*t73; 
complex<T> t39 = -t65; 
complex<T> t40 = -t57; 
complex<T> t47 = d38*s15; 
complex<T> t66 = -t73; 
complex<T> t84 = spb34*t65; 
complex<T> d1 = t2*t74; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa45*t74; d2 = T(1)/d2;
complex<T> d3 = (s12 - s45)*spa45*t74; d3 = T(1)/d3;
complex<T> d12 = t2*t5*t6; d12 = T(1)/d12;
complex<T> d13 = (s12 - s34)*t6*t71; d13 = T(1)/d13;
complex<T> d14 = t1*t2*t6; d14 = T(1)/d14;
complex<T> d15 = (s12 - s45)*t6*t71; d15 = T(1)/d15;
complex<T> d18 = (s12 - s34)*t74; d18 = T(1)/d18;
complex<T> d21 = spa45*t74*T(6); d21 = T(1)/d21;
complex<T> d36 = t74*T(2); d36 = T(1)/d36;
complex<T> t11 = -(d22*s12*t17) + d24*spa25*spb12*t57 + d25*s12*t21*t66 + d25*s12*t58*t72 - d23*spb12*t65*T(2); 
complex<T> t12 = d27*t57*t70 + d29*spb34*(t21*t66 + t58*t72) + d28*t84; 
complex<T> t13 = d18*spb45*t39 + d12*t40*t56 + d13*t40*t56 + d2*s35*t65 + d20*t21*t66 + d6*t21*t66 + d19*spb25*t72 + d5*spb25*t72 + d20*t58*t72 + d6*t58*t72 + d11*spa34*t18*t93; 
complex<T> t14 = -(d8*spb35*t17) + d10*spa14*spa15*t18 + d1*t39 + d2*s35*t39 + d3*s35*t39 + d12*t56*t57 + d13*t56*t57 + d14*t56*t57 + d15*t56*t57 + d18*spb45*t65 + d16*spa14*spa45*t68 - d5*spb25*t72 + d6*t41*t72 + d6*t21*t73 + d7*s34*t84 - d11*spa34*t18*t93 - d9*spb23*t91*t95 - spa45*t36*t68*T(2) + d4*spb13*t17*T(5); 
complex<T> t16 = d8*spb35*t17 - d10*spa14*spa15*t18 + d21*t39 + d7*s34*spb34*t39 + d14*t40*t56 + d15*t40*t56 + d3*s35*t65 - d16*spa14*spa45*t68 + t2*t36*t68 + d9*spb23*t91*t95 - d4*spb13*t17*T(5); 
complex<T> t32 = -t47; 
complex<T> t54 = d25*s15*(t21*t66 + t58*t72); 
complex<T> t99 = spb45*t40; 
complex<T> t103 = s15*t39; 
complex<T> t15 = -(d30*spb45*t17) + d31*spa25*t99; 
complex<T> t64 = d37*spb12*t103 + s12*t47*t58*t72 + s12*t21*t32*t73; 
complex<T> co1 = d26*spb23*t17; 
complex<T> co2 = d34*spb23*t84; 
complex<T> co3 = d35*t70*t99; 
complex<T> co4 = d36*spb45*t103; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t14*(*CI_users[0]->get_value(mc,ind,mu)) + t10*(*CI_users[1]->get_value(mc,ind,mu)) + t13*(*CI_users[3]->get_value(mc,ind,mu)) + t16*(*CI_users[4]->get_value(mc,ind,mu)) + t11*(*CI_users[5]->get_value(mc,ind,mu)) + t54*(*CI_users[6]->get_value(mc,ind,mu)) + co1*(*CI_users[7]->get_value(mc,ind,mu)) + t12*(*CI_users[8]->get_value(mc,ind,mu)) + t15*(*CI_users[9]->get_value(mc,ind,mu)) + t9*(*CI_users[10]->get_value(mc,ind,mu)) + co2*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)) + t64*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpmQpQm_LLT_wCI::\
C2q2Q1g_qmqpmQpQm_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpmQpQm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpmQpQm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:29:43 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa35 = SPA(3,5);
complex<T> s15 = S(1,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s12 = -(spa12*spb12);
complex<T> s35 = -(spa35*spb35);
complex<T> t2 = spb12*spb23; 
complex<T> t3 = spb34*spb45; 
complex<T> t4 = square(s12 - s34); 
complex<T> t5 = square(s23 - s45); 
complex<T> t6 = spb12*T(2); 
complex<T> t17 = square(spb24); 
complex<T> t18 = square(spb14); 
complex<T> t19 = square(spb25); 
complex<T> t20 = square(spa13); 
complex<T> t21 = square(spa35); 
complex<T> t22 = square(spb23); 
complex<T> t23 = square(spb34); 
complex<T> t24 = s12 - s45; 
complex<T> t35 = cube(spb24); 
complex<T> t36 = cube(spb14); 
complex<T> t37 = cube(spb25); 
complex<T> t38 = spb45*T(2); 
complex<T> t78 = s35*spb45; 
complex<T> t86 = s13*spb12; 
complex<T> d1 = (s12 - s34)*spb12*spb34; d1 = T(1)/d1;
complex<T> d21 = (s23 - s45)*spb23*spb45; d21 = T(1)/d21;
complex<T> d28 = spb23*spb45*cube(spb35); d28 = T(1)/d28;
complex<T> d33 = spb12*spb34*cube(spb13); d33 = T(1)/d33;
complex<T> t1 = square(t24); 
complex<T> t7 = spb13*t3; 
complex<T> t57 = t22*t36; 
complex<T> t58 = t23*t37; 
complex<T> t74 = t19*t21; 
complex<T> t75 = spb35*t2; 
complex<T> t85 = spb23*t20; 
complex<T> t91 = spa12*t35; 
complex<T> t94 = -(spa23*t35); 
complex<T> t98 = spa45*t35; 
complex<T> d3 = spb12*spb34*t24*T(6); d3 = T(1)/d3;
complex<T> d4 = cube(t24)*T(3); d4 = T(1)/d4;
complex<T> d6 = spb23*spb45*t24*T(6); d6 = T(1)/d6;
complex<T> d7 = t2*t24*t38; d7 = T(1)/d7;
complex<T> d11 = t24*t3*t6; d11 = T(1)/d11;
complex<T> d12 = t2*t3*T(3); d12 = T(1)/d12;
complex<T> d13 = (s12 - s34)*t2*t3; d13 = T(1)/d13;
complex<T> d14 = t2*t24*t3; d14 = T(1)/d14;
complex<T> d15 = t2*t38*t4; d15 = T(1)/d15;
complex<T> d24 = t3*t5*t6; d24 = T(1)/d24;
complex<T> d25 = (s23 - s45)*t2*t3; d25 = T(1)/d25;
complex<T> d26 = t3*cube(spb13); d26 = T(1)/d26;
complex<T> d27 = spb23*t3; d27 = T(1)/d27;
complex<T> d29 = spb12*t3*cube(spb13); d29 = T(1)/d29;
complex<T> d30 = spb12*t3; d30 = T(1)/d30;
complex<T> d31 = spb45*t2; d31 = T(1)/d31;
complex<T> d32 = spb45*t2*cube(spb35); d32 = T(1)/d32;
complex<T> d34 = spb34*t2; d34 = T(1)/d34;
complex<T> d35 = t2*cube(spb35); d35 = T(1)/d35;
complex<T> d36 = t3*T(2); d36 = T(1)/d36;
complex<T> d37 = spb45*t6; d37 = T(1)/d37;
complex<T> d38 = t3*cube(spb13)*T(2); d38 = T(1)/d38;
complex<T> d39 = t2*T(2); d39 = T(1)/d39;
complex<T> d40 = spb34*t2*T(2); d40 = T(1)/d40;
complex<T> d41 = spb23*t3*T(2); d41 = T(1)/d41;
complex<T> d42 = t2*cube(spb35)*T(2); d42 = T(1)/d42;
complex<T> t10 = -(d29*s23*t57) + d30*t94; 
complex<T> t14 = -(d33*spa45*t57) + d35*spa45*t58 + d34*t98; 
complex<T> t31 = d36*spa12; 
complex<T> t33 = d39*spa34; 
complex<T> t40 = -t85; 
complex<T> t49 = -t58; 
complex<T> d2 = spb34*t1*t6; d2 = T(1)/d2;
complex<T> d5 = spb23*t1*t38; d5 = T(1)/d5;
complex<T> d8 = t1*t6*t7; d8 = T(1)/d8;
complex<T> d9 = t24*t7*t86; d9 = T(1)/d9;
complex<T> d10 = t1*t3*t6; d10 = T(1)/d10;
complex<T> d16 = t1*t2*t38; d16 = T(1)/d16;
complex<T> d17 = t38*t4*t75; d17 = T(1)/d17;
complex<T> d18 = (s12 - s34)*t75*t78; d18 = T(1)/d18;
complex<T> d19 = t1*t38*t75; d19 = T(1)/d19;
complex<T> d20 = t24*t75*t78; d20 = T(1)/d20;
complex<T> d22 = t5*t6*t7; d22 = T(1)/d22;
complex<T> d23 = (s23 - s45)*t7*t86; d23 = T(1)/d23;
complex<T> t9 = d21*spa13*t17 + d25*s13*t35 + d24*spb24*t18*t40 + d22*t20*t57 - d23*t20*t57; 
complex<T> t12 = -(d31*spa34*t35) + d32*s34*t49; 
complex<T> t13 = -(d1*spa35*t17) - d13*s35*t35 + d17*t21*t58 + d18*t21*t58 - d15*spb24*spb34*t74; 
complex<T> t15 = d6*spa13*t17 + d5*s35*spa13*t17 + d1*spa35*t17 + d3*spa35*t17 - d2*s13*spa35*t17 - d4*spa35*spb12*spb34*t20 - d4*spa13*spb23*spb45*t21 + d12*t35 + d14*s13*t35 + d13*s35*t35 + d14*s35*t35 - d7*spa34*t35 + d10*spb24*t18*t40 + d17*t21*t49 + d18*t21*t49 + d19*t21*t49 + d20*t21*t49 + d8*t20*t57 - d9*t20*t57 + d15*spb24*spb34*t74 + d16*spb24*spb34*t74 + d11*t94; 
complex<T> t16 = -(d21*spa13*t17) - d6*spa13*t17 - d5*s35*spa13*t17 - d3*spa35*t17 + d2*s13*spa35*t17 + d4*spa35*spb12*spb34*t20 + d4*spa13*spb23*spb45*t21 + d12*t35 - d14*s13*t35 - d25*s13*t35 - d14*s35*t35 + d11*spa23*t35 + d7*spa34*t35 - d22*t20*t57 + d23*t20*t57 - d8*t20*t57 + d9*t20*t57 + d19*t21*t58 + d20*t21*t58 - d16*spb24*spb34*t74 + d10*spb24*t18*t85 + d24*spb24*t18*t85; 
complex<T> t46 = d28*spa12*t49 + d26*spa12*t57 + d27*t91; 
complex<T> t56 = spa23*t31; 
complex<T> t65 = spa45*t33; 
complex<T> t103 = t35*t56; 
complex<T> t105 = t35*t65; 
complex<T> t11 = t103 + d38*s23*spa12*t57; 
complex<T> t47 = t105 + d42*s34*spa45*t58; 
complex<T> co1 = -t103; 
complex<T> co2 = d37*spa34*t94; 
complex<T> co3 = -t105; 
complex<T> co4 = d40*s15*t98; 
complex<T> co5 = d41*s15*t91; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t15*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t13*(*CI_users[2]->get_value(mc,ind,mu)) + t16*(*CI_users[3]->get_value(mc,ind,mu)) + t46*(*CI_users[4]->get_value(mc,ind,mu)) + t10*(*CI_users[5]->get_value(mc,ind,mu)) + t12*(*CI_users[6]->get_value(mc,ind,mu)) + t14*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + t11*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + co5*(*CI_users[13]->get_value(mc,ind,mu)) + t47*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqppQpQm_LLT_wCI::\
C2q2Q1g_qmqppQpQm_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqppQpQm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqppQpQm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:30:49 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t2 = spa23*spa34; 
complex<T> t7 = square(spa15); 
complex<T> t8 = s12 - s45; 
complex<T> t11 = square(spb23); 
complex<T> t12 = square(spb34); 
complex<T> t13 = -(spa24*spb12); 
complex<T> t29 = spb23*spb34; 
complex<T> t32 = spa24*spa35; 
complex<T> t36 = spa45*spb15; 
complex<T> d13 = spa13*spa23*spa45; d13 = T(1)/d13;
complex<T> d14 = spa34*spa35; d14 = T(1)/d14;
complex<T> d16 = spa13*spa45; d16 = T(1)/d16;
complex<T> d17 = spa12*spa35; d17 = T(1)/d17;
complex<T> d18 = spa12*spa34*spa35; d18 = T(1)/d18;
complex<T> d19 = spa13*spa23; d19 = T(1)/d19;
complex<T> d21 = spa34*spa45*T(2); d21 = T(1)/d21;
complex<T> d22 = spa12*spa45*T(2); d22 = T(1)/d22;
complex<T> d23 = spa13*spa45*T(2); d23 = T(1)/d23;
complex<T> d24 = spa12*spa23*T(2); d24 = T(1)/d24;
complex<T> d27 = spa12*spa35*T(2); d27 = T(1)/d27;
complex<T> t1 = square(t8); 
complex<T> t46 = spa35*t29; 
complex<T> t47 = spa24*t7; 
complex<T> t66 = s12*t7; 
complex<T> t68 = s45*t7; 
complex<T> t70 = t13*t7; 
complex<T> d1 = spa12*spa45*t2*T(3); d1 = T(1)/d1;
complex<T> d2 = spa23*spa45*t8; d2 = T(1)/d2;
complex<T> d3 = spa34*t8*T(2); d3 = T(1)/d3;
complex<T> d4 = spa23*spa45*t8*T(6); d4 = T(1)/d4;
complex<T> d6 = spa23*t8*T(2); d6 = T(1)/d6;
complex<T> d7 = spa12*spa34*t8*T(6); d7 = T(1)/d7;
complex<T> d11 = t2*cube(t8)*T(3); d11 = T(1)/d11;
complex<T> d12 = spa12*spa34*t8; d12 = T(1)/d12;
complex<T> d15 = spa45*t2; d15 = T(1)/d15;
complex<T> d20 = spa12*t2; d20 = T(1)/d20;
complex<T> d25 = spa12*t2*T(2); d25 = T(1)/d25;
complex<T> d26 = spa45*t2*T(2); d26 = T(1)/d26;
complex<T> t4 = -(d13*t66) - d14*spb12*t7 + d15*t70; 
complex<T> t15 = d2*spb13; 
complex<T> t19 = -t47; 
complex<T> t31 = d11*spa13; 
complex<T> d5 = spa34*t1*T(6); d5 = T(1)/d5;
complex<T> d8 = spa23*t1*T(2); d8 = T(1)/d8;
complex<T> d9 = spa23*t1*T(6); d9 = T(1)/d9;
complex<T> d10 = spa34*t1*T(2); d10 = T(1)/d10;
complex<T> t39 = t31*t32; 
complex<T> t44 = d8*spa25; 
complex<T> t50 = t31*t36; 
complex<T> t63 = spb45*t19; 
complex<T> t5 = -(d4*spa15*spa25*spb23) - d7*spa14*spa15*spb34 - d5*spa12*spa35*t11 + d9*spa13*spa45*t12 + d1*t19 - d5*spa13*spa45*t29 - s12*t29*t39 - s45*t29*t39 - spa13*t29*t44 + d9*spa12*t46 + d10*spa14*t46 - spa12*t46*t50*T(2) - d12*spb35*t7*T(2) - t15*t7*T(2) + d3*spa15*spb23*T(3) + d6*spa15*spb34*T(3); 
complex<T> t6 = d4*spa15*spa25*spb23 + d7*spa14*spa15*spb34 + d5*spa12*spa35*t11 - d9*spa13*spa45*t12 + d1*t19 + d5*spa13*spa45*t29 + s12*t29*t39 + s45*t29*t39 + spa13*t29*t44 - d9*spa12*t46 - d10*spa14*t46 + spa12*t46*t50*T(2) + d12*spb35*t7*T(2) + t15*t7*T(2) - d3*spa15*spb23*T(3) - d6*spa15*spb34*T(3); 
complex<T> t27 = d20*t63 - d18*t68 - d19*spb45*t7; 
complex<T> co1 = d16*spb23*t7; 
complex<T> co2 = d17*spb34*t7; 
complex<T> co3 = d21*spb12*spb23*t47; 
complex<T> co4 = d22*t29*t47; 
complex<T> co5 = d23*spb23*t66; 
complex<T> co6 = d24*spb34*spb45*t47; 
complex<T> co7 = d25*s15*t63; 
complex<T> co8 = d26*s15*t70; 
complex<T> co9 = d27*spb34*t68; 
complex<T> co10 = Complex(0,1); 
SeriesC<T> result = co10*(t6*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[3]->get_value(mc,ind,mu)) + t4*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t27*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)) + co5*(*CI_users[10]->get_value(mc,ind,mu)) + co6*(*CI_users[11]->get_value(mc,ind,mu)) + co7*(*CI_users[12]->get_value(mc,ind,mu)) + co8*(*CI_users[13]->get_value(mc,ind,mu)) + co9*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmmQp_LLT_wCI::\
C2q2Q1g_qmqpQmmQp_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmmQp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, m, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmmQp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:37:19 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa13 = SPA(1,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s15 = -(spa15*spb15);
complex<T> s14 = -(spa14*spb14);
complex<T> t2 = spb34*T(2); 
complex<T> t4 = square(s23 - s45); 
complex<T> t6 = spb15*spb45; 
complex<T> t9 = spb14*spb15; 
complex<T> t10 = spb45*T(2); 
complex<T> t22 = square(spb25); 
complex<T> t23 = square(spa14); 
complex<T> t24 = square(spb12); 
complex<T> t25 = square(spb15); 
complex<T> t26 = square(spb24); 
complex<T> t27 = square(spb45); 
complex<T> t28 = square(spb35); 
complex<T> t29 = square(s14); 
complex<T> t30 = s15 - s23; 
complex<T> t31 = -(spb13*T(3)); 
complex<T> t32 = s15 - s23 + s45; 
complex<T> t33 = square(spa13); 
complex<T> t43 = cube(spa14); 
complex<T> t44 = -(spb25*T(3)); 
complex<T> t49 = s45*spa45; 
complex<T> t54 = s12*s23; 
complex<T> t58 = cube(spb25); 
complex<T> t59 = spb23*spb45; 
complex<T> t63 = -(spa23*spa34); 
complex<T> t68 = spa13*T(3); 
complex<T> t81 = spa45*T(2); 
complex<T> t83 = -(spa34*spb13); 
complex<T> t88 = spb25*spb35; 
complex<T> t111 = spa13*spb35; 
complex<T> t117 = spa12*spb13; 
complex<T> t120 = spa14*spb25; 
complex<T> t130 = spb34*T(4); 
complex<T> d2 = (s12 - s45)*spb34; d2 = T(1)/d2;
complex<T> d4 = spb34*square(s12 - s45); d4 = T(1)/d4;
complex<T> d6 = spb12*spb34*spb45*T(6); d6 = T(1)/d6;
complex<T> d8 = (s12 - s45)*spb34*spb45; d8 = T(1)/d8;
complex<T> d9 = s13*(s12 - s45)*spb34*spb45; d9 = T(1)/d9;
complex<T> d33 = (s23 - s45)*spb34*spb45; d33 = T(1)/d33;
complex<T> d36 = s13*(s23 - s45)*spb34*spb45; d36 = T(1)/d36;
complex<T> d43 = spb34; d43 = T(1)/d43;
complex<T> d44 = spb34*spb45; d44 = T(1)/d44;
complex<T> d45 = spb13*spb34*spb45; d45 = T(1)/d45;
complex<T> d46 = spb34*spb45*square(spb13); d46 = T(1)/d46;
complex<T> d55 = s12; d55 = T(1)/d55;
complex<T> d56 = spb12*spb45; d56 = T(1)/d56;
complex<T> d57 = s12*spb45; d57 = T(1)/d57;
complex<T> d59 = s12*spb34; d59 = T(1)/d59;
complex<T> d60 = spb12*spb34; d60 = T(1)/d60;
complex<T> d61 = spb13*spb34; d61 = T(1)/d61;
complex<T> d63 = spb34*square(spb13); d63 = T(1)/d63;
complex<T> d68 = spb15*spb24*T(2); d68 = T(1)/d68;
complex<T> d71 = spb12*spb15*spb23*T(2); d71 = T(1)/d71;
complex<T> d72 = spb15*spb23*T(2); d72 = T(1)/d72;
complex<T> d79 = s12*T(2); d79 = T(1)/d79;
complex<T> t1 = square(s45 + t30); 
complex<T> t3 = square(t30); 
complex<T> t7 = spb23*(s45 + t30); 
complex<T> t42 = d68*spa34; 
complex<T> t45 = -t111; 
complex<T> t46 = -(spb13*t24); 
complex<T> t48 = -t59; 
complex<T> t65 = spb13*t23; 
complex<T> t67 = t24*t27; 
complex<T> t76 = d79*spa34*(-(s45*t120) + spa13*spa45*t88); 
complex<T> t78 = t25*t26; 
complex<T> t79 = t23*t27; 
complex<T> t90 = spb12*t33; 
complex<T> t99 = spb35*t22; 
complex<T> t102 = d4*spa13; 
complex<T> t121 = -(spb12*t28); 
complex<T> t122 = d2*T(2); 
complex<T> t128 = spa14*t22; 
complex<T> t137 = t111*t44; 
complex<T> t138 = d45*T(2); 
complex<T> t142 = t68*t88; 
complex<T> t145 = spb35*t24; 
complex<T> t153 = t22*T(3); 
complex<T> t156 = spb25*t23; 
complex<T> t161 = spa14*t44; 
complex<T> d1 = s12*t2; d1 = T(1)/d1;
complex<T> d3 = s12*(s12 - s45)*t2; d3 = T(1)/d3;
complex<T> d5 = s12*t2*square(s12 - s45); d5 = T(1)/d5;
complex<T> d7 = s12*spb45*t2; d7 = T(1)/d7;
complex<T> d10 = t2*square(s12 - s45); d10 = T(1)/d10;
complex<T> d14 = spb15*t2*t30; d14 = T(1)/d14;
complex<T> d15 = spb15*spb34*t30; d15 = T(1)/d15;
complex<T> d19 = spb23*t10*t30; d19 = T(1)/d19;
complex<T> d20 = spb23*spb34*t6; d20 = T(1)/d20;
complex<T> d26 = (s23 - s45)*spb15*spb23*t2; d26 = T(1)/d26;
complex<T> d29 = spb45*t4; d29 = T(1)/d29;
complex<T> d30 = spb23*t10*t4; d30 = T(1)/d30;
complex<T> d31 = (s23 - s45)*spb12*spb45*t2; d31 = T(1)/d31;
complex<T> d32 = spb12*spb23*t2*t6; d32 = T(1)/d32;
complex<T> d34 = (s23 - s45)*spb34*t59; d34 = T(1)/d34;
complex<T> d37 = spb15*spb23*t2*t4; d37 = T(1)/d37;
complex<T> d38 = spb23*t2*t4*t9; d38 = T(1)/d38;
complex<T> d40 = s12*t130; d40 = T(1)/d40;
complex<T> d41 = spb12*spb45*t130; d41 = T(1)/d41;
complex<T> d42 = s12*spb45*t130; d42 = T(1)/d42;
complex<T> d48 = t59*cube(s45 + t30); d48 = T(1)/d48;
complex<T> d52 = spb45*cube(s45 + t30); d52 = T(1)/d52;
complex<T> d58 = spb23*cube(s45 + t30); d58 = T(1)/d58;
complex<T> d65 = t2*t6; d65 = T(1)/d65;
complex<T> d66 = spb13*spb45*t2; d66 = T(1)/d66;
complex<T> d67 = spb45*t2*square(spb13); d67 = T(1)/d67;
complex<T> d69 = spb12*t6*T(2); d69 = T(1)/d69;
complex<T> d70 = t6*T(2); d70 = T(1)/d70;
complex<T> d73 = spb23*cube(s45 + t30)*T(2); d73 = T(1)/d73;
complex<T> d75 = spb12*spb23*t2; d75 = T(1)/d75;
complex<T> d76 = spb23*t2; d76 = T(1)/d76;
complex<T> d78 = t2*t59; d78 = T(1)/d78;
complex<T> t11 = d42*t137 + d41*t153 + d40*t161; 
complex<T> t12 = -(d43*t120) + d46*s12*t121 + d44*spa12*t22 + d44*spb25*t45 + s12*t138*t88; 
complex<T> t37 = d31*s13; 
complex<T> t39 = -(d34*spa15); 
complex<T> t51 = -(d30*s14); 
complex<T> t64 = -t99; 
complex<T> t66 = -t78; 
complex<T> t69 = spb23*t1; 
complex<T> t70 = d10*spa34; 
complex<T> t73 = d3*s45; 
complex<T> t80 = -t90; 
complex<T> t87 = spa34*(-(d55*t120) + d56*t22 + d57*spb25*t45); 
complex<T> t89 = t65*T(3); 
complex<T> t92 = -(d29*spb15); 
complex<T> t94 = d5*t49; 
complex<T> t104 = d32*t31; 
complex<T> t109 = spa45*(d71*t58*t83 + d72*spa34*t99); 
complex<T> t113 = d37*spb45; 
complex<T> t129 = t28*t90; 
complex<T> t136 = t65*t67; 
complex<T> t139 = s23*t42; 
complex<T> t141 = t46*t79; 
complex<T> t147 = t102*t81; 
complex<T> t154 = t145*t23; 
complex<T> t164 = d26*t128; 
complex<T> t204 = d66*t88; 
complex<T> d11 = spb23*t3*T(2); d11 = T(1)/d11;
complex<T> d12 = t2*t3; d12 = T(1)/d12;
complex<T> d13 = spb34*t30*t7; d13 = T(1)/d13;
complex<T> d16 = spb15*spb23*t2*t3; d16 = T(1)/d16;
complex<T> d17 = t1*t30*t59; d17 = T(1)/d17;
complex<T> d18 = t10*t3*t7; d18 = T(1)/d18;
complex<T> d21 = spb34*t30*t6*t7; d21 = T(1)/d21;
complex<T> d22 = spb15*spb34*t3; d22 = T(1)/d22;
complex<T> d23 = spb23*t2*t3*t9; d23 = T(1)/d23;
complex<T> d24 = spb34*t30*t7*t9; d24 = T(1)/d24;
complex<T> d25 = (s23 - s45)*spb34*t7; d25 = T(1)/d25;
complex<T> d27 = (s23 - s45)*t1*t59; d27 = T(1)/d27;
complex<T> d28 = t10*t4*t7; d28 = T(1)/d28;
complex<T> d35 = (s23 - s45)*spb34*t6*t7; d35 = T(1)/d35;
complex<T> d39 = (s23 - s45)*spb34*t7*t9; d39 = T(1)/d39;
complex<T> d49 = spb34*t1*t59; d49 = T(1)/d49;
complex<T> d51 = spb34*t1; d51 = T(1)/d51;
complex<T> d53 = spb34*t1*t6; d53 = T(1)/d53;
complex<T> d54 = spb34*t1*t9; d54 = T(1)/d54;
complex<T> t13 = d67*t121*t54 - t204*t54 + d65*spa23*(-(t117*t58) + s12*t64); 
complex<T> t14 = d8*t137 + d7*t142 + spa14*t59*t70 + d9*t28*t80 + t142*t94 - d2*t120*T(2) - spa45*t102*t88*T(2) + d1*t120*T(3) + t120*t73*T(3) + d6*t22*T(13); 
complex<T> t53 = d25*spb24; 
complex<T> t56 = d78*spa15*(-(t117*t58) + s12*t64); 
complex<T> t82 = d35*t29; 
complex<T> t100 = t43*t66; 
complex<T> t114 = t37*T(3); 
complex<T> t135 = d22*spb45; 
complex<T> t140 = d27*t78; 
complex<T> t144 = spb25*t89; 
complex<T> t160 = spa13*t92; 
complex<T> t166 = t48*t70; 
complex<T> t194 = t139*t22; 
complex<T> d47 = spb34*t69; d47 = T(1)/d47;
complex<T> d50 = spb14*spb34*t69; d50 = T(1)/d50;
complex<T> d62 = spb15*spb34*t69; d62 = T(1)/d62;
complex<T> d64 = spb34*t69*t9; d64 = T(1)/d64;
complex<T> d74 = t2*t69; d74 = T(1)/d74;
complex<T> d77 = spb14*t2*t69; d77 = T(1)/d77;
complex<T> t15 = d17*t100 + d18*t100 + spa34*spb13*t120*t135 + d23*t141 + d24*t141 + d14*spa24*t153 - d11*spb12*t156 + d12*spb15*spb24*t23 + d13*spb24*t156*t31 + d16*t49*t64 - d15*spa34*t88 + d21*t29*t99 + d20*t99*T(2) - d19*t128*T(3); 
complex<T> t16 = d27*t100 + d29*spa13*spb15*t120 + d30*s14*t128 + d23*t136 + d24*t136 + d38*t136 + d33*t137 + d39*t141 + d13*spb24*t144 - t113*t154 + d11*spb12*t156 - d12*spb15*spb24*t23 + t156*t31*t53 + t104*t58 + d21*t29*t64 + d17*t43*t78 + d18*t43*t78 + d28*t43*t78 + d36*t28*t80 + t120*t135*t83 + d15*spa34*t88 + d20*t99 + d34*spa15*t99 + d16*t49*t99 + t82*t99 + d19*t128*T(3) + spb13*t164*T(3) - d14*spa24*t22*T(3) - t22*t37*T(3); 
complex<T> t17 = d28*t100 + d36*t129 + d9*t129 + d39*t136 + d42*t137 + d38*t141 + d33*t142 + d8*t142 + d41*t153 + t113*t154 + t120*(t122 + t160) + d40*t161 + spa14*t166 + t114*t22 + t164*t31 + t140*t43 + t128*t51 + t144*t53 + t104*t58 + t161*t73 + t64*t82 + t147*t88 + t137*t94 + t39*t99; 
complex<T> t18 = d58*spa45*t100 + d59*s45*t120 + d63*spa45*t121 + d64*s45*t136 + d47*s45*spb24*t144 + d60*spa45*t22 + d59*spa45*spb25*t45 + d61*t81*t88 + d62*spa45*t29*t99; 
complex<T> t19 = d50*spa15*t141 + d47*s15*spb24*t144 + d48*s15*t43*t78 + d49*spa15*t29*t99; 
complex<T> t20 = d73*s15*spa45*t100 + d77*s45*spa15*t141 + d74*s15*s45*spb24*t144 - d75*spa15*spa45*spb13*t58 + d74*spa15*spa45*t29*t64 + d76*spa15*spa45*t99; 
complex<T> t21 = d46*s23*t121 + d54*spa23*t136 + d51*spa23*spb24*t144 + d53*spa23*t29*t64 + d52*spa23*t43*t78 + s23*t138*t88; 
complex<T> t98 = t194 + d69*spb13*t58*t63 + d70*spa23*spa34*t99; 
complex<T> co1 = t204*t54*T(3); 
complex<T> co2 = -t194; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t14*(*CI_users[0]->get_value(mc,ind,mu)) + t15*(*CI_users[1]->get_value(mc,ind,mu)) + t16*(*CI_users[2]->get_value(mc,ind,mu)) + t11*(*CI_users[3]->get_value(mc,ind,mu)) + t17*(*CI_users[4]->get_value(mc,ind,mu)) + t12*(*CI_users[5]->get_value(mc,ind,mu)) + t19*(*CI_users[6]->get_value(mc,ind,mu)) + t21*(*CI_users[7]->get_value(mc,ind,mu)) + t87*(*CI_users[8]->get_value(mc,ind,mu)) + t18*(*CI_users[9]->get_value(mc,ind,mu)) + t13*(*CI_users[10]->get_value(mc,ind,mu)) + t98*(*CI_users[11]->get_value(mc,ind,mu)) + co1*(*CI_users[12]->get_value(mc,ind,mu)) + t109*(*CI_users[13]->get_value(mc,ind,mu)) + co2*(*CI_users[14]->get_value(mc,ind,mu)) + t20*(*CI_users[15]->get_value(mc,ind,mu)) + t56*(*CI_users[16]->get_value(mc,ind,mu)) + t76*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmpQp_LLT_wCI::\
C2q2Q1g_qmqpQmpQp_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmpQp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, p, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmpQp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:41:36 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa13 = SPA(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s14 = -(spa14*spb14);
complex<T> s45 = -(spa45*spb45);
complex<T> s24 = -(spa24*spb24);
complex<T> t1 = spa45*T(2); 
complex<T> t4 = spa15*spa24; 
complex<T> t8 = spa15*spa34; 
complex<T> t9 = (s15 - s34)*spa45; 
complex<T> t11 = square(s24); 
complex<T> t23 = square(spa13); 
complex<T> t24 = square(spb24); 
complex<T> t25 = square(spa12); 
complex<T> t26 = square(spa14); 
complex<T> t27 = square(spa23); 
complex<T> t28 = square(spa34); 
complex<T> t30 = square(spa35); 
complex<T> t31 = square(spb25); 
complex<T> t32 = -(spa25*T(3)); 
complex<T> t34 = s12 - s34; 
complex<T> t50 = cube(spb24); 
complex<T> t53 = -(spb25*T(2)); 
complex<T> t65 = cube(spa13); 
complex<T> t66 = spa25*T(3); 
complex<T> t67 = spb23*spb34; 
complex<T> t69 = spa13*spa35; 
complex<T> t83 = spb25*T(2); 
complex<T> t90 = s15*spb14; 
complex<T> t97 = spb24*T(3); 
complex<T> t112 = spa25*spb12; 
complex<T> t118 = spa45*T(4); 
complex<T> t124 = spa13*spb24; 
complex<T> d1 = spa12*spa34*spa45*T(6); d1 = T(1)/d1;
complex<T> d21 = spa15*square(s15 - s23)*T(2); d21 = T(1)/d21;
complex<T> d23 = (s15 - s23)*s24*spa15*spa45; d23 = T(1)/d23;
complex<T> d36 = spa34*square(s15 - s34); d36 = T(1)/d36;
complex<T> d38 = (s15 - s23)*spa23*spa45; d38 = T(1)/d38;
complex<T> d39 = spa23*spa45*square(s15 - s23); d39 = T(1)/d39;
complex<T> d43 = spa34*spa45; d43 = T(1)/d43;
complex<T> d44 = spa45; d44 = T(1)/d44;
complex<T> d46 = spa34*cube(spa24); d46 = T(1)/d46;
complex<T> d47 = spa45*square(spa24); d47 = T(1)/d47;
complex<T> d48 = spa23*spa45*cube(spa24); d48 = T(1)/d48;
complex<T> d49 = spa23*spa34*spa45; d49 = T(1)/d49;
complex<T> d51 = spa15*spa45*square(spa24); d51 = T(1)/d51;
complex<T> d52 = spa15*spa45*cube(spa24); d52 = T(1)/d52;
complex<T> d54 = spa15*spa23*spa45*cube(spa24); d54 = T(1)/d54;
complex<T> d55 = s12*spa45; d55 = T(1)/d55;
complex<T> d56 = spa15*cube(spa24); d56 = T(1)/d56;
complex<T> d57 = spa12*spa45; d57 = T(1)/d57;
complex<T> d58 = spa15*spa23*spa45; d58 = T(1)/d58;
complex<T> d60 = spa12*spa34; d60 = T(1)/d60;
complex<T> d61 = s12; d61 = T(1)/d61;
complex<T> d62 = s12*spa34; d62 = T(1)/d62;
complex<T> d64 = spa23*square(s15 - s23 + s45)*T(2); d64 = T(1)/d64;
complex<T> d68 = s12*T(2); d68 = T(1)/d68;
complex<T> d69 = spa12*spa15*spa23*T(2); d69 = T(1)/d69;
complex<T> d70 = spa15*spa23*T(2); d70 = T(1)/d70;
complex<T> d73 = spa15*cube(spa24)*T(2); d73 = T(1)/d73;
complex<T> d74 = spa12*spa23*spa34*T(2); d74 = T(1)/d74;
complex<T> d75 = spa23*spa34*T(2); d75 = T(1)/d75;
complex<T> t48 = -t69; 
complex<T> t58 = -t67; 
complex<T> t60 = d64*s14; 
complex<T> t70 = t24*t28; 
complex<T> t72 = spa23*t1; 
complex<T> t73 = spa14*t32; 
complex<T> t84 = spa25*t25; 
complex<T> t85 = t26*t27; 
complex<T> t94 = spa35*t23; 
complex<T> t95 = spa12*t31; 
complex<T> t102 = d39*spa25; 
complex<T> t107 = t8*T(2); 
complex<T> t114 = spb45*t90; 
complex<T> t116 = spb25*t69; 
complex<T> t125 = d36*spa23; 
complex<T> t129 = spa13*t24; 
complex<T> t132 = spa35*t25; 
complex<T> t134 = t69*t83; 
complex<T> t138 = spb24*t23; 
complex<T> d2 = s12*t1; d2 = T(1)/d2;
complex<T> d3 = spa45*t34; d3 = T(1)/d3;
complex<T> d4 = s12*t1*t34; d4 = T(1)/d4;
complex<T> d5 = s12*spa34*t1; d5 = T(1)/d5;
complex<T> d6 = spa34*spa45*t34; d6 = T(1)/d6;
complex<T> d7 = spa34*spa45*t34*(s15 + t34); d7 = T(1)/d7;
complex<T> d8 = spa45*square(t34); d8 = T(1)/d8;
complex<T> d9 = s12*t1*square(t34); d9 = T(1)/d9;
complex<T> d10 = t1*square(t34); d10 = T(1)/d10;
complex<T> d11 = (s15 - s34)*spa12*spa34*t1; d11 = T(1)/d11;
complex<T> d13 = spa23*spa45*t8; d13 = T(1)/d13;
complex<T> d14 = (s15 - s23)*spa23*spa45*t8; d14 = T(1)/d14;
complex<T> d15 = spa23*t8*t9; d15 = T(1)/d15;
complex<T> d17 = t8*t9; d17 = T(1)/d17;
complex<T> d22 = t1*square(s15 - s23); d22 = T(1)/d22;
complex<T> d24 = s24*spa15*t9; d24 = T(1)/d24;
complex<T> d26 = (s15 - s23)*s24*spa23*spa45*t4; d26 = T(1)/d26;
complex<T> d28 = s24*spa23*t4*t9; d28 = T(1)/d28;
complex<T> d30 = (s15 - s23)*t11*t8; d30 = T(1)/d30;
complex<T> d33 = (s15 - s34)*t11*t8; d33 = T(1)/d33;
complex<T> d34 = spa34*t9; d34 = T(1)/d34;
complex<T> d35 = spa34*(s15 + t34)*t9; d35 = T(1)/d35;
complex<T> d40 = spa12*spa34*t118; d40 = T(1)/d40;
complex<T> d41 = s12*t118; d41 = T(1)/d41;
complex<T> d42 = s12*spa34*t118; d42 = T(1)/d42;
complex<T> d45 = spa34*spa45*square(s15 + t34); d45 = T(1)/d45;
complex<T> d50 = t8*cube(spa24); d50 = T(1)/d50;
complex<T> d53 = spa45*t8; d53 = T(1)/d53;
complex<T> d59 = spa45*square(s15 + t34); d59 = T(1)/d59;
complex<T> d63 = t1*t8; d63 = T(1)/d63;
complex<T> d65 = spa34*t1*square(s15 + t34); d65 = T(1)/d65;
complex<T> d66 = spa12*spa15*t1; d66 = T(1)/d66;
complex<T> d67 = spa15*t1; d67 = T(1)/d67;
complex<T> d71 = spa15*t1*square(spa24); d71 = T(1)/d71;
complex<T> d72 = spa15*t1*cube(spa24); d72 = T(1)/d72;
complex<T> t13 = d63*spb23*(t112*t65 + s12*t94); 
complex<T> t14 = d44*t124 + d43*(t116 - spb12*t23) + d45*s12*(s25*t134 + t30*t95); 
complex<T> t39 = d65*s12; 
complex<T> t40 = -t60; 
complex<T> t43 = d11*s25; 
complex<T> t51 = -t84; 
complex<T> t52 = -t85; 
complex<T> t71 = -t95; 
complex<T> t78 = d10*spa15; 
complex<T> t81 = spb45*(d68*s34*t124 + d68*spb25*spb34*t48 + d69*spa25*spb34*t65 - d70*spb34*t94); 
complex<T> t93 = d66*spa25*t65*t67 + d67*t58*t94; 
complex<T> t104 = spb45*(d62*t116 + d61*t124 - d60*t23); 
complex<T> t108 = -(d40*T(3)); 
complex<T> t119 = d42*T(3); 
complex<T> t147 = spa13*t73; 
complex<T> t175 = t114*t23; 
complex<T> d12 = spa12*t72*t8; d12 = T(1)/d12;
complex<T> d16 = (s15 - s23)*t72; d16 = T(1)/d16;
complex<T> d18 = (s15 - s23)*t107; d18 = T(1)/d18;
complex<T> d19 = t107*square(s15 - s34); d19 = T(1)/d19;
complex<T> d20 = (s15 - s34)*spa15*t72; d20 = T(1)/d20;
complex<T> d25 = t4*t72*square(s15 - s23); d25 = T(1)/d25;
complex<T> d27 = t4*t72*square(s15 - s34); d27 = T(1)/d27;
complex<T> d29 = spa15*t72*square(s15 - s34); d29 = T(1)/d29;
complex<T> d31 = s24*t107*square(s15 - s23); d31 = T(1)/d31;
complex<T> d32 = s24*t107*square(s15 - s34); d32 = T(1)/d32;
complex<T> d37 = spa15*t72*square(s15 - s23); d37 = T(1)/d37;
complex<T> d76 = spa34*t72; d76 = T(1)/d76;
complex<T> t12 = d51*s23*t147 + d50*s23*t52 + d52*spb23*t28*t84 - d53*spb23*t94; 
complex<T> t15 = spb15*(d47*t147 + d48*t28*t51 + d46*t52 + d49*t94) + d45*s15*(s25*t134 + t30*t95); 
complex<T> t17 = d71*s23*s34*t147 + d72*s34*spb23*t28*t84 + d73*s23*spb34*t85 + d67*t67*t94; 
complex<T> t18 = d55*spb34*t116 - d55*s34*t124 + d59*s25*spb34*t134 + d51*s34*t147 - d57*spb34*t23 + d54*s34*t28*t51 + d56*spb34*t85 - d58*spb34*t94 + d59*spb34*t30*t95; 
complex<T> t19 = d6*t116 + d8*spb34*t134 + d7*s25*t53*t69 + d7*t30*t71 - spa34*spb24*spb45*t78 + d3*t124*T(2) - d5*t116*T(3) - d9*s34*spb34*t116*T(3) - d2*t124*T(3) - d4*s34*t124*T(3) - d1*t23*T(13); 
complex<T> t45 = d16*spb14; 
complex<T> t46 = d76*spb15; 
complex<T> t64 = t175*t40 + spb15*spb45*(d74*spa25*t65 - d75*t94); 
complex<T> t74 = s15*t39; 
complex<T> t123 = t116*t119 + t108*t23 + d41*spa13*t97; 
complex<T> t159 = d12*t66; 
complex<T> t20 = t123 + spb25*t124*t125 + d35*s25*t134 + d7*s25*t134 + d19*s24*t138 - d29*spa34*t132*t24 + d34*spb25*t48 + d6*spb25*t48 + d32*t50*t52 + d33*t50*t52 + t159*t65 + d24*spa14*t129*t66 + d20*t138*t66 + d8*spb34*t53*t69 + spa34*spb24*spb45*t78 + d27*t70*t84 + d28*t70*t84 - d15*s24*t94 + d17*spb23*t94 + d35*t30*t95 + d7*t30*t95 + d4*s34*spa13*t97 - d3*t124*T(2) + d9*s34*spb34*t116*T(3) - t23*t43*T(3); 
complex<T> t21 = -(spa34*spb45*t102*t124) + d21*spa12*t129 - d22*spa14*spa23*t24 + d38*spb45*t48 + d30*t50*t52 + d31*t50*t52 + d23*spa14*t129*t66 + d25*t70*t84 + d26*t70*t84 - d14*s24*t94 + d37*s34*spb34*t94 - d13*t94*T(2) - d18*t138*T(3) + t23*t45*T(3); 
complex<T> t22 = d34*t116 + spa34*spb45*t102*t124 - spb25*t124*t125 - d21*spa12*t129 - d19*s24*t138 + d22*spa14*spa23*t24 + d29*spa34*t132*t24 + d23*t147*t24 + d24*t147*t24 + d20*t138*t32 + t159*t65 + d38*spb45*t69 + d35*s25*t53*t69 + d25*t51*t70 + d26*t51*t70 + d27*t51*t70 + d28*t51*t70 + d35*t30*t71 + d30*t50*t85 + d31*t50*t85 + d32*t50*t85 + d33*t50*t85 - d13*t94 + d14*s24*t94 + d15*s24*t94 - d17*spb23*t94 - d37*s34*spb34*t94 + d18*t23*t97 + t23*t43*T(3) - t23*t45*T(3); 
complex<T> t103 = s25*t74; 
complex<T> t16 = spb25*t103*t48 + t30*t74*t95; 
complex<T> t115 = t112*t46*t65 + s12*t46*t94 + t103*t116*T(3); 
complex<T> co1 = t175*t60; 
complex<T> co2 = Complex(0,1); 
SeriesC<T> result = co2*(t19*(*CI_users[0]->get_value(mc,ind,mu)) + t22*(*CI_users[1]->get_value(mc,ind,mu)) + t21*(*CI_users[2]->get_value(mc,ind,mu)) + t20*(*CI_users[3]->get_value(mc,ind,mu)) + t123*(*CI_users[4]->get_value(mc,ind,mu)) + t14*(*CI_users[5]->get_value(mc,ind,mu)) + t15*(*CI_users[6]->get_value(mc,ind,mu)) + t12*(*CI_users[7]->get_value(mc,ind,mu)) + t18*(*CI_users[8]->get_value(mc,ind,mu)) + t104*(*CI_users[9]->get_value(mc,ind,mu)) + t13*(*CI_users[10]->get_value(mc,ind,mu)) + co1*(*CI_users[11]->get_value(mc,ind,mu)) + t16*(*CI_users[12]->get_value(mc,ind,mu)) + t93*(*CI_users[13]->get_value(mc,ind,mu)) + t81*(*CI_users[14]->get_value(mc,ind,mu)) + t17*(*CI_users[15]->get_value(mc,ind,mu)) + t64*(*CI_users[16]->get_value(mc,ind,mu)) + t115*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpmQm_LLT_wCI::\
C2q2Q1g_qmqpQpmQm_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpmQm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, m, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpmQm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:41:45 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa14 = SPA(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> s15 = S(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s23 = S(2,3);
complex<T> t1 = spb45*T(2); 
complex<T> t7 = square(spb23); 
complex<T> t10 = -(spa34*spa45); 
complex<T> t11 = -(spb23*T(3)); 
complex<T> t12 = spa15*spb35; 
complex<T> t19 = s34*spa14; 
complex<T> t20 = spa45*spb25; 
complex<T> t21 = spb45*T(4); 
complex<T> t26 = spa14*spb23; 
complex<T> d2 = (s12 - s34)*spb45; d2 = T(1)/d2;
complex<T> d5 = spb12*spb34*spb45*T(3); d5 = T(1)/d5;
complex<T> d7 = spb45*square(s12 - s34); d7 = T(1)/d7;
complex<T> d14 = spb45; d14 = T(1)/d14;
complex<T> d15 = spb34*spb45; d15 = T(1)/d15;
complex<T> d16 = s12*spb45; d16 = T(1)/d16;
complex<T> d17 = spb12*spb45; d17 = T(1)/d17;
complex<T> d18 = s12; d18 = T(1)/d18;
complex<T> d19 = spb12*spb34; d19 = T(1)/d19;
complex<T> d20 = s12*spb34; d20 = T(1)/d20;
complex<T> d23 = s12*T(2); d23 = T(1)/d23;
complex<T> d24 = spb12*T(2); d24 = T(1)/d24;
complex<T> d25 = spb12*spb34*T(2); d25 = T(1)/d25;
complex<T> t13 = -t19; 
complex<T> t16 = d7*spa34; 
complex<T> t18 = spb23*t12; 
complex<T> t41 = spa12*t7; 
complex<T> t45 = spa34*t7; 
complex<T> t46 = spa45*t7; 
complex<T> d1 = s12*t1; d1 = T(1)/d1;
complex<T> d3 = s12*(s12 - s34)*t1; d3 = T(1)/d3;
complex<T> d4 = (s12 - s34)*spb12*t1; d4 = T(1)/d4;
complex<T> d6 = t1*square(s12 - s34); d6 = T(1)/d6;
complex<T> d8 = s12*t1*square(s12 - s34); d8 = T(1)/d8;
complex<T> d9 = s12*spb34*t1; d9 = T(1)/d9;
complex<T> d10 = (s12 - s34)*spb34*t1; d10 = T(1)/d10;
complex<T> d11 = s12*t21; d11 = T(1)/d11;
complex<T> d12 = spb12*spb34*t21; d12 = T(1)/d12;
complex<T> d13 = s12*spb34*t21; d13 = T(1)/d13;
complex<T> d21 = spb34*t1; d21 = T(1)/d21;
complex<T> d22 = spb12*t1; d22 = T(1)/d22;
complex<T> t17 = d14*t26 + d15*(t18 + t41); 
complex<T> t22 = d11*T(3); 
complex<T> t23 = d8*spa34; 
complex<T> t24 = d6*spb34; 
complex<T> t25 = d16*spb23*t13 + d16*spa34*t18 + d17*t45; 
complex<T> t31 = d20*spa45*t18 + d18*spa45*t26 + d19*t46; 
complex<T> t32 = t18*T(3); 
complex<T> t35 = d23*t10*t18 + d23*spa45*spb23*t19 + d24*t10*t7; 
complex<T> t4 = t22*t26 + d13*t32 + d12*t7*T(3); 
complex<T> t5 = d1*spa14*t11 + d9*t11*t12 + d3*t11*t19 + s34*t11*t12*t23 - spa14*t20*t24 + d10*t32 + t16*t18*T(2) + d2*t26*T(2) + d5*t7*T(2) + d4*spb23*t20*T(3); 
complex<T> t6 = d10*t11*t12 + d4*t11*t20 + spa14*t20*t24 + t22*t26 + d13*t32 + s34*t23*t32 - t16*t18*T(2) - d2*t26*T(2) + d3*spb23*t19*T(3) - d12*t7*T(3); 
complex<T> co1 = d21*s23*t41; 
complex<T> co2 = d22*s23*t45; 
complex<T> co3 = d25*s15*t46; 
complex<T> co4 = d21*s15*t41; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t4*(*CI_users[2]->get_value(mc,ind,mu)) + t17*(*CI_users[3]->get_value(mc,ind,mu)) + t25*(*CI_users[4]->get_value(mc,ind,mu)) + t31*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t35*(*CI_users[8]->get_value(mc,ind,mu)) + co3*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQppQm_LLT_wCI::\
C2q2Q1g_qmqpQppQm_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQppQm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, p, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQppQm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:41:56 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> s15 = S(1,5);
complex<T> s23 = S(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = spa34*T(2); 
complex<T> t7 = square(spa15); 
complex<T> t10 = -(spa15*T(3)); 
complex<T> t14 = spa34*T(4); 
complex<T> t18 = spa35*spb23; 
complex<T> t20 = s45*spb24; 
complex<T> t21 = spa13*spb34; 
complex<T> t26 = -(spa15*spb24); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> d5 = (s12 - s45)*spa34; d5 = T(1)/d5;
complex<T> d9 = spa34*square(s12 - s45); d9 = T(1)/d9;
complex<T> d14 = spa34*spa45; d14 = T(1)/d14;
complex<T> d15 = spa34; d15 = T(1)/d15;
complex<T> d16 = spa12*spa45; d16 = T(1)/d16;
complex<T> d17 = s12*spa45; d17 = T(1)/d17;
complex<T> d18 = s12; d18 = T(1)/d18;
complex<T> d19 = s12*spa34; d19 = T(1)/d19;
complex<T> d20 = spa12*spa34; d20 = T(1)/d20;
complex<T> d22 = spa12*spa45*T(2); d22 = T(1)/d22;
complex<T> d23 = spa12*T(2); d23 = T(1)/d23;
complex<T> d25 = s12*T(2); d25 = T(1)/d25;
complex<T> t11 = -t18; 
complex<T> t12 = -t20; 
complex<T> t24 = d9*spb45; 
complex<T> t44 = -(spb12*t7); 
complex<T> t48 = -(spb45*t7); 
complex<T> t51 = spb34*t7; 
complex<T> d2 = s12*spa45*t1; d2 = T(1)/d2;
complex<T> d3 = (s12 - s45)*spa45*t1; d3 = T(1)/d3;
complex<T> d4 = s12*t1; d4 = T(1)/d4;
complex<T> d6 = s12*(s12 - s45)*t1; d6 = T(1)/d6;
complex<T> d7 = (s12 - s45)*spa12*t1; d7 = T(1)/d7;
complex<T> d8 = t1*square(s12 - s45); d8 = T(1)/d8;
complex<T> d10 = s12*t1*square(s12 - s45); d10 = T(1)/d10;
complex<T> d11 = spa12*spa45*t14; d11 = T(1)/d11;
complex<T> d12 = s12*spa45*t14; d12 = T(1)/d12;
complex<T> d13 = s12*t14; d13 = T(1)/d13;
complex<T> d21 = spa45*t1; d21 = T(1)/d21;
complex<T> d24 = spa12*t1; d24 = T(1)/d24;
complex<T> t16 = d10*spb45; 
complex<T> t17 = d25*spa15*spb34*(t12 + spb45*t18); 
complex<T> t22 = d8*spa45; 
complex<T> t25 = d15*t26 + d14*(spa15*t11 + t44); 
complex<T> t31 = d12*t18; 
complex<T> t37 = d19*spa15*(spb45*t11 + t20) + d20*t48; 
complex<T> t39 = d13*spb24; 
complex<T> t47 = -t51; 
complex<T> t4 = t10*(t31 + t39) - d11*t7*T(3); 
complex<T> t30 = s45*t16; 
complex<T> t32 = d17*spa15*spb34*t11 + d18*spb34*t26 + d16*t47; 
complex<T> t5 = d6*t10*t20 - spb24*t21*t22 + t10*t18*t30 + t10*t31 + t10*t39 + d5*spa15*spb24*T(2) + spa15*t18*t24*T(2) + d3*spa15*t18*T(3) + d7*spa15*t21*T(3) + d11*t7*T(3); 
complex<T> t6 = d3*t10*t18 + d7*t10*t21 + spb24*t21*t22 - d5*spa15*spb24*T(2) - spa15*t18*t24*T(2) - d1*t7*T(2) + d4*spa15*spb24*T(3) + d2*spa15*t18*T(3) + d6*spa15*t20*T(3) + spa15*t18*t30*T(3); 
complex<T> co1 = d21*s23*t44; 
complex<T> co2 = d22*s23*t47; 
complex<T> co3 = d23*spb45*t51; 
complex<T> co4 = d24*s15*t48; 
complex<T> co5 = d21*s15*t44; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t6*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + t25*(*CI_users[3]->get_value(mc,ind,mu)) + t32*(*CI_users[4]->get_value(mc,ind,mu)) + t37*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)) + co5*(*CI_users[10]->get_value(mc,ind,mu)) + t17*(*CI_users[11]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmQpm_LLT_wCI::\
C2q2Q1g_qmqpQmQpm_LLT_wCI
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
     C2q2Q1g_qmqpQmQpm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, m}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmQpm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:42:39 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = S(2,3);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = -(spa35*spb35);
complex<T> t1 = square(s12 - s34); 
complex<T> t3 = spb15*T(2); 
complex<T> t4 = spb12*spb34; 
complex<T> t6 = (s12 - s45)*spb15; 
complex<T> t18 = square(spb24); 
complex<T> t19 = square(spa35); 
complex<T> t20 = square(spb23); 
complex<T> t21 = square(spb45); 
complex<T> t22 = square(spb14); 
complex<T> t24 = square(spa13); 
complex<T> t26 = cube(spb35); 
complex<T> t27 = s12 + s15 - s34; 
complex<T> t31 = spb15*spb34; 
complex<T> t34 = s12*s23; 
complex<T> t44 = s25*spa25; 
complex<T> t45 = spb25*spb34; 
complex<T> t67 = spb15*spb45; 
complex<T> t78 = spb14*spb24; 
complex<T> t91 = spb13*spb24; 
complex<T> t93 = spa35*spb45; 
complex<T> d5 = cube(s12 - s34)*T(3); d5 = T(1)/d5;
complex<T> d7 = (s12 - s34)*spb12*spb45*T(2); d7 = T(1)/d7;
complex<T> d31 = spb15*square(spb13); d31 = T(1)/d31;
complex<T> d32 = spb13*spb15; d32 = T(1)/d32;
complex<T> d42 = spb34*spb45*T(2); d42 = T(1)/d42;
complex<T> t7 = spb35*t4; 
complex<T> t25 = -t44; 
complex<T> t57 = spb14*t18; 
complex<T> t58 = t19*t21; 
complex<T> t59 = spb13*t20; 
complex<T> t72 = spa15*t45; 
complex<T> t73 = spb23*t22; 
complex<T> t85 = t18*t44; 
complex<T> t92 = spb45*t20; 
complex<T> d1 = spb12*t1*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spb34*t3; d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*t27*t31; d3 = T(1)/d3;
complex<T> d4 = t1*t3*t4; d4 = T(1)/d4;
complex<T> d6 = s13*spb45*t6; d6 = T(1)/d6;
complex<T> d8 = spb12*spb45*t6; d8 = T(1)/d8;
complex<T> d9 = spb45*t3*t4; d9 = T(1)/d9;
complex<T> d10 = (s12 - s34)*t4*t67; d10 = T(1)/d10;
complex<T> d11 = spb45*t4*t6; d11 = T(1)/d11;
complex<T> d12 = spb12*spb45*t1*T(3); d12 = T(1)/d12;
complex<T> d13 = t3*t4*square(s12 - s45); d13 = T(1)/d13;
complex<T> d14 = t1*t31; d14 = T(1)/d14;
complex<T> d19 = (s15 - s34)*t31; d19 = T(1)/d19;
complex<T> d20 = (s15 - s34)*t27*t31; d20 = T(1)/d20;
complex<T> d21 = s13*(s23 - s45)*t67; d21 = T(1)/d21;
complex<T> d22 = t4*t67*T(6); d22 = T(1)/d22;
complex<T> d23 = t31*square(t27); d23 = T(1)/d23;
complex<T> d24 = t67*square(spb13); d24 = T(1)/d24;
complex<T> d25 = spb13*t67; d25 = T(1)/d25;
complex<T> d26 = spb45*t31; d26 = T(1)/d26;
complex<T> d27 = t26*t31; d27 = T(1)/d27;
complex<T> d28 = spb34*square(t27); d28 = T(1)/d28;
complex<T> d29 = spb15*square(t27); d29 = T(1)/d29;
complex<T> d30 = spb12*spb15*t26; d30 = T(1)/d30;
complex<T> d33 = spb15*t4; d33 = T(1)/d33;
complex<T> d34 = spb15*t26*t4; d34 = T(1)/d34;
complex<T> d35 = spb34*spb45*t3; d35 = T(1)/d35;
complex<T> d36 = spb34*square(t27)*T(2); d36 = T(1)/d36;
complex<T> d37 = spb12*spb45*t3; d37 = T(1)/d37;
complex<T> d38 = spb45*t3*square(spb13); d38 = T(1)/d38;
complex<T> d39 = spb13*spb45*t3; d39 = T(1)/d39;
complex<T> d40 = spb12*t3; d40 = T(1)/d40;
complex<T> d41 = t4*T(2); d41 = T(1)/d41;
complex<T> d43 = spb12*t26*t3; d43 = T(1)/d43;
complex<T> t9 = s23*(d24*t73 + d25*t78); 
complex<T> t10 = d19*spa25*t18 + d20*t85; 
complex<T> t14 = spa34*(d29*t18*t25 + d30*t21*t59); 
complex<T> t38 = d40*spa34; 
complex<T> t40 = -t57; 
complex<T> t41 = -t59; 
complex<T> t55 = t34*(d38*t73 + d39*t78); 
complex<T> t65 = d21*t24; 
complex<T> t99 = spa12*t57; 
complex<T> d15 = t1*t3*t7; d15 = T(1)/d15;
complex<T> d16 = (s12 - s34)*s35*spb15*t7; d16 = T(1)/d16;
complex<T> d17 = t3*t7*square(s12 - s45); d17 = T(1)/d17;
complex<T> d18 = s35*t6*t7; d18 = T(1)/d18;
complex<T> t13 = d23*s12*t18*t25 + d27*spa12*t21*t41 + d24*s12*t73 + d25*s12*t78 + d26*t99*T(2); 
complex<T> t15 = -(d19*spa25*t18) - d7*spa35*t18 + d1*spb23*spb24*t19 + d20*t18*t25 + d3*t18*t25 + d10*s35*t40 + d12*spb24*t19*t45 + d22*t57 + d4*s45*spa45*t57 + d15*t58*t59 + d16*t58*t59 - d14*spa15*t91*t93 - d5*t19*t72*T(2) + d2*spa25*t18*T(3); 
complex<T> t17 = d11*s35*t40 + d8*spa34*t57 + d17*t58*t59 + d18*t58*t59 - d6*t24*t73 - t65*t73 - d13*spb14*t19*t92; 
complex<T> t71 = spa45*t38; 
complex<T> t97 = spa45*t40; 
complex<T> t12 = d34*s45*t21*t41 + d31*spa45*t73 + d32*spa45*t78 + d33*t97; 
complex<T> t16 = d7*spa35*t18 - d1*spb23*spb24*t19 + d8*spa34*t40 - d12*spb24*t19*t45 + d9*t57 + d10*s35*t57 + d11*s35*t57 + d15*t41*t58 + d16*t41*t58 + d17*t41*t58 + d18*t41*t58 + d6*t24*t73 + d3*t85 + d13*spb14*t19*t92 + d14*spa15*t91*t93 + d4*s45*t97 + d5*t19*t72*T(2) - d2*spa25*t18*T(3); 
complex<T> t56 = d43*s45*spa34*t21*t59 + t57*t71; 
complex<T> co1 = t65*t73; 
complex<T> co2 = d28*spa15*t85; 
complex<T> co3 = d35*s23*t99; 
complex<T> co4 = d36*s12*spa15*t85; 
complex<T> co5 = d37*s23*spa34*t57; 
complex<T> co6 = t40*t71; 
complex<T> co7 = d41*spa15*t97; 
complex<T> co8 = d42*spa12*spa15*t40; 
complex<T> co9 = Complex(0,1); 
SeriesC<T> result = co9*(t16*(*CI_users[0]->get_value(mc,ind,mu)) + t10*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + t15*(*CI_users[3]->get_value(mc,ind,mu)) + t17*(*CI_users[4]->get_value(mc,ind,mu)) + t13*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t9*(*CI_users[7]->get_value(mc,ind,mu)) + t14*(*CI_users[8]->get_value(mc,ind,mu)) + t12*(*CI_users[9]->get_value(mc,ind,mu)) + co3*(*CI_users[10]->get_value(mc,ind,mu)) + co4*(*CI_users[11]->get_value(mc,ind,mu)) + co5*(*CI_users[12]->get_value(mc,ind,mu)) + t55*(*CI_users[13]->get_value(mc,ind,mu)) + co6*(*CI_users[14]->get_value(mc,ind,mu)) + co7*(*CI_users[15]->get_value(mc,ind,mu)) + co8*(*CI_users[16]->get_value(mc,ind,mu)) + t56*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmQpp_LLT_wCI::\
C2q2Q1g_qmqpQmQpp_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmQpp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, p}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmQpp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:43:31 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> s23 = S(2,3);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> t2 = square(s12 - s34); 
complex<T> t3 = spa15*spa45; 
complex<T> t4 = spa12*T(2); 
complex<T> t6 = spa34*spa45; 
complex<T> t17 = square(spa13); 
complex<T> t18 = square(spb25); 
complex<T> t19 = square(spa23); 
complex<T> t20 = square(spa15); 
complex<T> t21 = square(s25); 
complex<T> t22 = square(spa14); 
complex<T> t24 = s15 - s34; 
complex<T> t25 = square(spb24); 
complex<T> t32 = -(s23*s34); 
complex<T> t67 = -(spa13*spa14); 
complex<T> t83 = spa15*spa24; 
complex<T> t86 = spa13*spb25; 
complex<T> d4 = (s12 - s34)*spa15*spa34*T(2); d4 = T(1)/d4;
complex<T> d11 = cube(s12 - s34)*T(3); d11 = T(1)/d11;
complex<T> d22 = spa35*spa45; d22 = T(1)/d22;
complex<T> d24 = spa45*square(spa24); d24 = T(1)/d24;
complex<T> d25 = spa24*spa45; d25 = T(1)/d25;
complex<T> d30 = spa12*spa35*spa45; d30 = T(1)/d30;
complex<T> d34 = spa12*spa35; d34 = T(1)/d34;
complex<T> t1 = square(s12 + t24); 
complex<T> t5 = square(t24); 
complex<T> t7 = spa12*(s12 + t24); 
complex<T> t34 = d11*spb45; 
complex<T> t38 = spa25*t6; 
complex<T> t48 = spa24*t20; 
complex<T> t49 = spa14*t17; 
complex<T> t50 = t18*t19; 
complex<T> t73 = spa23*t22; 
complex<T> t74 = spa13*t18; 
complex<T> d1 = spa12*spa34*t3*T(6); d1 = T(1)/d1;
complex<T> d3 = t2*t4*t6; d3 = T(1)/d3;
complex<T> d5 = spa34*t2*T(2); d5 = T(1)/d5;
complex<T> d6 = spa15*spa34*t2*T(3); d6 = T(1)/d6;
complex<T> d9 = (s12 - s34)*spa45*t4; d9 = T(1)/d9;
complex<T> d10 = spa12*spa45*t2; d10 = T(1)/d10;
complex<T> d13 = spa34*t24*t3; d13 = T(1)/d13;
complex<T> d14 = (s15 - s23)*s24*t3; d14 = T(1)/d14;
complex<T> d15 = s24*t24*t3; d15 = T(1)/d15;
complex<T> d19 = spa34*t3*t4; d19 = T(1)/d19;
complex<T> d20 = spa34*t3; d20 = T(1)/d20;
complex<T> d28 = t3*square(spa24); d28 = T(1)/d28;
complex<T> d29 = spa24*t3; d29 = T(1)/d29;
complex<T> d31 = spa12*t3; d31 = T(1)/d31;
complex<T> d35 = spa34*t3*T(2); d35 = T(1)/d35;
complex<T> d36 = t3*square(spa24)*T(2); d36 = T(1)/d36;
complex<T> d37 = spa24*t3*T(2); d37 = T(1)/d37;
complex<T> d38 = t3*t4; d38 = T(1)/d38;
complex<T> d39 = spa35*t4; d39 = T(1)/d39;
complex<T> d40 = spa15*t4; d40 = T(1)/d40;
complex<T> d41 = spa34*t4; d41 = T(1)/d41;
complex<T> d42 = t6*T(2); d42 = T(1)/d42;
complex<T> t15 = spb45*(d39*s34*t17 + d40*spb34*t49); 
complex<T> t35 = -t49; 
complex<T> t36 = -t50; 
complex<T> t57 = d14*t25; 
complex<T> t58 = -t73; 
complex<T> t65 = spa35*t34; 
complex<T> t71 = spb15*t49; 
complex<T> d2 = (s12 - s34)*spa34*t3*t7; d2 = T(1)/d2;
complex<T> d7 = t2*t38*t4; d7 = T(1)/d7;
complex<T> d8 = (s12 - s34)*t38*t7; d8 = T(1)/d8;
complex<T> d12 = spa34*t24*t3*t7; d12 = T(1)/d12;
complex<T> d16 = t4*t5*t6; d16 = T(1)/d16;
complex<T> d17 = t38*t4*t5; d17 = T(1)/d17;
complex<T> d18 = t24*t38*t7; d18 = T(1)/d18;
complex<T> d21 = spa34*t1*t3; d21 = T(1)/d21;
complex<T> d23 = t1*t38; d23 = T(1)/d23;
complex<T> d26 = spa12*t1*t6; d26 = T(1)/d26;
complex<T> d27 = spa12*t1*t38; d27 = T(1)/d27;
complex<T> d32 = spa12*t1*t3; d32 = T(1)/d32;
complex<T> d33 = spa12*spa25*spa45*t1; d33 = T(1)/d33;
complex<T> d43 = t1*t6*T(2); d43 = T(1)/d43;
complex<T> d44 = t1*t38*T(2); d44 = T(1)/d44;
complex<T> t9 = s23*(d28*t58 + d29*t67); 
complex<T> t11 = d27*s15*t48*t50 + d24*spb15*t58 + d25*spb15*t67 + d26*t21*t71; 
complex<T> t13 = -(d30*s34*t17) + d31*spb34*t35 + d32*spb34*t21*t35 + d33*spb34*t48*t50 + d28*s34*t58 + d29*s34*t67; 
complex<T> t46 = d37*spa13*spa14*t32 + d38*s23*spb34*t35 + d36*t32*t73; 
complex<T> t66 = t36*t48; 
complex<T> t78 = d16*spa15; 
complex<T> t94 = spb12*t35; 
complex<T> t10 = d13*spb12*t49 + d12*t21*t49 + d15*t25*t58 + t57*t58 + d17*t66 + d18*t66 + spa14*t50*t78; 
complex<T> t12 = -(d22*spb12*t17) + d21*spb12*t21*t49 + d23*spb12*t66 + d20*t94; 
complex<T> t14 = -(d4*spb25*t17) + d1*t35 + d3*s15*spb15*t35 + d2*t21*t49 + d7*t66 + d8*t66 - d5*spa23*t74 - d6*spa12*spa35*t74 + d10*spb45*t83*t86 - spa12*t18*t65*T(2) + d9*spb35*t17*T(5); 
complex<T> t16 = d4*spb25*t17 + d19*t35 + d12*t21*t35 + d2*t21*t35 + d17*t48*t50 + d18*t48*t50 + d7*t48*t50 + d8*t48*t50 + t18*t4*t65 + d3*s15*t71 + d15*t25*t73 + d5*spa23*t74 + d6*spa12*spa35*t74 + spa14*t36*t78 - d10*spb45*t83*t86 + d13*t94 - d9*spb35*t17*T(5); 
complex<T> t47 = d44*s15*spb12*t66 + d42*spb12*t71 + d43*spb15*t21*t94; 
complex<T> co1 = t57*t73; 
complex<T> co2 = d34*spb45*t17; 
complex<T> co3 = d35*s23*t94; 
complex<T> co4 = d41*spb45*t71; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t14*(*CI_users[0]->get_value(mc,ind,mu)) + t10*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + t16*(*CI_users[3]->get_value(mc,ind,mu)) + t12*(*CI_users[5]->get_value(mc,ind,mu)) + t11*(*CI_users[6]->get_value(mc,ind,mu)) + t9*(*CI_users[7]->get_value(mc,ind,mu)) + t13*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + co3*(*CI_users[10]->get_value(mc,ind,mu)) + t46*(*CI_users[11]->get_value(mc,ind,mu)) + t15*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)) + t47*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpQmm_LLT_wCI::\
C2q2Q1g_qmqpQpQmm_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpQmm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, m}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpQmm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:44:12 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb45 = SPB(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> t3 = spb34*T(2); 
complex<T> t4 = spb15*spb45; 
complex<T> t11 = square(spb23); 
complex<T> t12 = s12 - s34; 
complex<T> t14 = -(s25*spa25); 
complex<T> t17 = square(spa15); 
complex<T> t18 = square(spa45); 
complex<T> t19 = spa12*spb14; 
complex<T> t36 = spa15*spa45; 
complex<T> t40 = spb14*spb35; 
complex<T> d14 = (s15 - s34)*spb15*spb34; d14 = T(1)/d14;
complex<T> d18 = spb35*spb45; d18 = T(1)/d18;
complex<T> d22 = spb12*spb35*spb45; d22 = T(1)/d22;
complex<T> d23 = spb12*spb35; d23 = T(1)/d23;
complex<T> d27 = spb12*spb15*T(2); d27 = T(1)/d27;
complex<T> d30 = spb12*spb35*T(2); d30 = T(1)/d30;
complex<T> t1 = square(t12); 
complex<T> t23 = -t36; 
complex<T> t38 = spa25*t11; 
complex<T> t46 = spb14*t11; 
complex<T> t47 = spb35*t36; 
complex<T> t69 = t11*t19; 
complex<T> t71 = s34*t11; 
complex<T> d1 = spb15*t12*T(2); d1 = T(1)/d1;
complex<T> d3 = spb15*spb34*t12*T(6); d3 = T(1)/d3;
complex<T> d4 = spb15*spb34*t12; d4 = T(1)/d4;
complex<T> d5 = spb15*spb34*t12*(s15 + t12); d5 = T(1)/d5;
complex<T> d7 = spb45*t12*T(2); d7 = T(1)/d7;
complex<T> d8 = spb12*spb45*t12; d8 = T(1)/d8;
complex<T> d9 = spb12*spb45*t12*T(6); d9 = T(1)/d9;
complex<T> d10 = spb12*spb34*t4*T(3); d10 = T(1)/d10;
complex<T> d13 = t4*cube(t12)*T(3); d13 = T(1)/d13;
complex<T> d15 = (s15 - s34)*spb15*spb34*(s15 + t12); d15 = T(1)/d15;
complex<T> d16 = spb15*spb34*square(s15 + t12); d16 = T(1)/d16;
complex<T> d17 = spb34*t4; d17 = T(1)/d17;
complex<T> d19 = spb34*square(s15 + t12); d19 = T(1)/d19;
complex<T> d20 = spb15*square(s15 + t12); d20 = T(1)/d20;
complex<T> d21 = spb12*t4; d21 = T(1)/d21;
complex<T> d24 = t3*t4; d24 = T(1)/d24;
complex<T> d25 = t3*square(s15 + t12); d25 = T(1)/d25;
complex<T> d26 = spb12*t4*T(2); d26 = T(1)/d26;
complex<T> d28 = spb12*t3; d28 = T(1)/d28;
complex<T> d29 = spb45*t3; d29 = T(1)/d29;
complex<T> t37 = d13*spb25; 
complex<T> t45 = d18*spa12*t11 + d16*s12*t11*t14 + d17*t69; 
complex<T> t67 = s25*t38; 
complex<T> t70 = spa34*t46; 
complex<T> d2 = spb15*t1*T(2); d2 = T(1)/d2;
complex<T> d6 = spb15*t1*T(6); d6 = T(1)/d6;
complex<T> d11 = spb45*t1*T(6); d11 = T(1)/d11;
complex<T> d12 = spb45*t1*T(2); d12 = T(1)/d12;
complex<T> t7 = d14*t38 + d15*t67; 
complex<T> t35 = d20*spa34*t11*t14 + d21*t70 + d22*t71; 
complex<T> t41 = d2*spb13; 
complex<T> t55 = t37*t40; 
complex<T> t9 = -(d3*spa15*spb13*spb23) - d9*spa45*spb23*spb24 - d11*spb12*spb35*t17 + d6*spb25*spb34*t18 + d11*spb25*spb34*t23 - d4*t38 + spb25*t23*t41 + d10*t46 + d6*spb12*t47 + d12*spb24*t47 + s12*t23*t55 + s34*t23*t55 + d5*t67 - d8*spa35*t11*T(2) - spa23*spb12*spb34*t37*t47*T(2) + d7*spa15*spb23*T(3) + d1*spa45*spb23*T(3); 
complex<T> t10 = d3*spa15*spb13*spb23 + d9*spa45*spb23*spb24 + d15*t11*t14 + d5*t11*t14 + d11*spb12*spb35*t17 - d6*spb25*spb34*t18 + d6*spb12*spb35*t23 + d12*spb24*spb35*t23 + d11*spb25*spb34*t36 - d14*t38 + d4*t38 + spb25*t36*t41 + d10*t46 + spa23*spb12*t3*t37*t47 + s12*t36*t55 + s34*t36*t55 + d8*spa35*t11*T(2) - d7*spa15*spb23*T(3) - d1*spa45*spb23*T(3); 
complex<T> co1 = d19*spa15*t67; 
complex<T> co2 = -(d23*spa45*t11); 
complex<T> co3 = d24*s23*t69; 
complex<T> co4 = d25*s12*spa15*t67; 
complex<T> co5 = d26*s23*t70; 
complex<T> co6 = -(d27*spa45*t70); 
complex<T> co7 = d28*t23*t46; 
complex<T> co8 = -(d29*spa15*t69); 
complex<T> co9 = -(d30*spa45*t71); 
complex<T> co10 = Complex(0,1); 
SeriesC<T> result = co10*(t9*(*CI_users[0]->get_value(mc,ind,mu)) + t7*(*CI_users[1]->get_value(mc,ind,mu)) + t10*(*CI_users[2]->get_value(mc,ind,mu)) + t45*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t35*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + co3*(*CI_users[8]->get_value(mc,ind,mu)) + co4*(*CI_users[9]->get_value(mc,ind,mu)) + co5*(*CI_users[10]->get_value(mc,ind,mu)) + co6*(*CI_users[11]->get_value(mc,ind,mu)) + co7*(*CI_users[12]->get_value(mc,ind,mu)) + co8*(*CI_users[13]->get_value(mc,ind,mu)) + co9*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpQmp_LLT_wCI::\
C2q2Q1g_qmqpQpQmp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpQmp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, p}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpQmp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:45:03 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa35 = SPA(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = S(2,3);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> t3 = spa15*spa34; 
complex<T> t4 = square(s15 - s34); 
complex<T> t5 = square(s12 - s45); 
complex<T> t6 = spa12*T(2); 
complex<T> t7 = spa34*spa45; 
complex<T> t18 = square(spa14); 
complex<T> t19 = square(spb25); 
complex<T> t20 = square(spa13); 
complex<T> t21 = square(spa24); 
complex<T> t22 = square(spb35); 
complex<T> t23 = square(spa15); 
complex<T> t24 = square(spa45); 
complex<T> t25 = square(s25); 
complex<T> t26 = s12 - s34; 
complex<T> t27 = s12 + s15 - s34; 
complex<T> t29 = cube(spa35); 
complex<T> t36 = cube(spa14); 
complex<T> t37 = cube(spa13); 
complex<T> t38 = cube(spa24); 
complex<T> t40 = spa12*spa45; 
complex<T> t63 = spb12*spb15; 
complex<T> t1 = square(s15 + t26); 
complex<T> t2 = square(t26); 
complex<T> t8 = spa12*t26; 
complex<T> t60 = t19*t23; 
complex<T> t61 = t24*t37; 
complex<T> t70 = spa35*t3; 
complex<T> t77 = spa14*t22; 
complex<T> t82 = spb45*t36; 
complex<T> t97 = spb12*t36; 
complex<T> d1 = t3*t40*T(3); d1 = T(1)/d1;
complex<T> d4 = (s12 - s45)*t3*t40; d4 = T(1)/d4;
complex<T> d5 = t26*t6*t7; d5 = T(1)/d5;
complex<T> d6 = t26*t3*T(6); d6 = T(1)/d6;
complex<T> d13 = (s12 - s45)*t40; d13 = T(1)/d13;
complex<T> d14 = cube(t26)*T(3); d14 = T(1)/d14;
complex<T> d16 = t3*t5*t6; d16 = T(1)/d16;
complex<T> d21 = t26*t3*t6; d21 = T(1)/d21;
complex<T> d22 = (s15 - s34)*(s15 + t26)*t3*t40; d22 = T(1)/d22;
complex<T> d23 = (s15 - s34)*t3; d23 = T(1)/d23;
complex<T> d24 = t4*t6*t7; d24 = T(1)/d24;
complex<T> d25 = spa25*t4*t6*t7; d25 = T(1)/d25;
complex<T> d26 = (s15 - s34)*spa12*spa25*(s15 + t26)*t7; d26 = T(1)/d26;
complex<T> d27 = spa45*t3; d27 = T(1)/d27;
complex<T> d29 = t29*t3; d29 = T(1)/d29;
complex<T> d34 = spa12*spa15*t29; d34 = T(1)/d34;
complex<T> d36 = spa12*t29*t3; d36 = T(1)/d36;
complex<T> d37 = spa12*t3; d37 = T(1)/d37;
complex<T> d38 = spa45*t3*T(2); d38 = T(1)/d38;
complex<T> d41 = spa15*spa45*t6; d41 = T(1)/d41;
complex<T> d42 = spa15*t6; d42 = T(1)/d42;
complex<T> d43 = spa34*t6; d43 = T(1)/d43;
complex<T> d44 = t7*T(2); d44 = T(1)/d44;
complex<T> d45 = spa15*t29*t6; d45 = T(1)/d45;
complex<T> t15 = d36*s45*t61 + d37*t82; 
complex<T> t45 = d42*spb34; 
complex<T> t51 = -t60; 
complex<T> t52 = -t61; 
complex<T> t59 = d24*t21; 
complex<T> d2 = spa45*(s15 + t26)*t3*t8; d2 = T(1)/d2;
complex<T> d3 = spa45*t3*t8; d3 = T(1)/d3;
complex<T> d7 = t2*t3*T(2); d7 = T(1)/d7;
complex<T> d8 = t2*t6*t7; d8 = T(1)/d8;
complex<T> d9 = spa25*t2*t6*t7; d9 = T(1)/d9;
complex<T> d10 = spa25*(s15 + t26)*t7*t8; d10 = T(1)/d10;
complex<T> d11 = spa45*t2*t6; d11 = T(1)/d11;
complex<T> d12 = spa45*t8*T(6); d12 = T(1)/d12;
complex<T> d15 = t2*t3*t6; d15 = T(1)/d15;
complex<T> d17 = t2*t6*t70; d17 = T(1)/d17;
complex<T> d18 = s35*t70*t8; d18 = T(1)/d18;
complex<T> d19 = t5*t6*t70; d19 = T(1)/d19;
complex<T> d20 = s35*(s12 - s45)*spa12*t70; d20 = T(1)/d20;
complex<T> d28 = spa45*t1*t3; d28 = T(1)/d28;
complex<T> d30 = spa25*t1*t7; d30 = T(1)/d30;
complex<T> d31 = spa12*t1*t7; d31 = T(1)/d31;
complex<T> d32 = spa12*spa25*t1*t7; d32 = T(1)/d32;
complex<T> d33 = spa15*t1*t40; d33 = T(1)/d33;
complex<T> d35 = spa25*t1*t40; d35 = T(1)/d35;
complex<T> d39 = t1*t7*T(2); d39 = T(1)/d39;
complex<T> d40 = spa25*t1*t7*T(2); d40 = T(1)/d40;
complex<T> t11 = d31*spb15*t25*t36 + d32*s15*t38*t60; 
complex<T> t13 = spb34*(-(d33*t25*t36) + d34*t52 + d35*t38*t60); 
complex<T> t14 = d13*spb35*t18 + d4*s35*t36 + d19*t22*t52 + d20*t22*t52 + d16*spa45*t20*t77; 
complex<T> t34 = -t45; 
complex<T> t66 = d8*t21; 
complex<T> t76 = t38*t51; 
complex<T> t10 = -(d23*spb25*t18) + d22*t25*t36 + spa14*spa15*t19*t59 + d25*t76 + d26*t76; 
complex<T> t12 = d29*spb12*t61 + d30*spb12*t76 + d28*t25*t97 - d27*t97*T(2); 
complex<T> t16 = d23*spb25*t18 + d6*spb25*t18 + d7*s35*spb25*t18 + d12*spb35*t18 - d11*s25*spb35*t18 - d14*spb25*t22*t3 - d1*t36 + d3*s35*t36 - d5*spb15*t36 - d2*t25*t36 - d22*t25*t36 - d14*spb35*t19*t40 + d17*t22*t52 + d18*t22*t52 - spa14*spa15*t19*t59 + d10*t38*t60 + d25*t38*t60 + d26*t38*t60 + d9*t38*t60 - spa14*spa15*t19*t66 + d15*spa45*t20*t77 - d21*t82; 
complex<T> t17 = -(d6*spb25*t18) - d7*s35*spb25*t18 - d12*spb35*t18 - d13*spb35*t18 + d11*s25*spb35*t18 + d14*spb25*t22*t3 - d1*t36 - d3*s35*t36 - d4*s35*t36 + d5*spb15*t36 + d2*t25*t36 + d14*spb35*t19*t40 + d17*t22*t61 + d18*t22*t61 + d19*t22*t61 + d20*t22*t61 + spa14*spa15*t19*t66 + d10*t76 + d9*t76 - d15*spa45*t20*t77 - d16*spa45*t20*t77 + d21*t82; 
complex<T> t49 = d45*s45*spb34*t52 + t34*t82; 
complex<T> t50 = -(d39*t25*t36*t63) + d40*s15*spb12*t76; 
complex<T> co1 = -(d38*s23*t97); 
complex<T> co2 = -(d41*s23*spb34*t36); 
complex<T> co3 = t45*t82; 
complex<T> co4 = d43*spb15*t82; 
complex<T> co5 = d44*t36*t63; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t17*(*CI_users[0]->get_value(mc,ind,mu)) + t10*(*CI_users[1]->get_value(mc,ind,mu)) + t16*(*CI_users[2]->get_value(mc,ind,mu)) + t14*(*CI_users[3]->get_value(mc,ind,mu)) + t12*(*CI_users[4]->get_value(mc,ind,mu)) + t11*(*CI_users[5]->get_value(mc,ind,mu)) + t13*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t50*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + co5*(*CI_users[13]->get_value(mc,ind,mu)) + t49*(*CI_users[14]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmmqpQmQp_LLT_wCI::\
C2q2Q1g_qmmqpQmQp_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmmqpQmQp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, qp, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmmqpQmQp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:51:32 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> s14 = -(spa14*spb14);
complex<T> s15 = -(spa15*spb15);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s24 = -(spa24*spb24);
complex<T> t2 = spb12*T(2); 
complex<T> t4 = spb15*spb24; 
complex<T> t8 = spb15*spb23; 
complex<T> t9 = (s15 - s23)*spb12; 
complex<T> t11 = square(s24); 
complex<T> t23 = square(spb35); 
complex<T> t24 = square(spa24); 
complex<T> t25 = square(spb45); 
complex<T> t26 = square(spb23); 
complex<T> t27 = square(spb25); 
complex<T> t28 = square(spb34); 
complex<T> t30 = square(spa14); 
complex<T> t31 = square(spb13); 
complex<T> t32 = -(spb14*T(3)); 
complex<T> t36 = s15 - s23 + s45; 
complex<T> t49 = cube(spa24); 
complex<T> t53 = -(spb13*T(2)); 
complex<T> t54 = -(s23*T(3)); 
complex<T> t69 = cube(spb35); 
complex<T> t71 = -(spa15*spb14); 
complex<T> t74 = spa14*spb35; 
complex<T> t78 = spb14*T(3); 
complex<T> t93 = spa12*spa23; 
complex<T> t96 = s14*s45; 
complex<T> t108 = s15*spa25; 
complex<T> t126 = spa24*spb35; 
complex<T> t132 = spa45*spb14; 
complex<T> t135 = spb12*T(4); 
complex<T> d12 = s24*(s15 - s34)*spb12*spb15; d12 = T(1)/d12;
complex<T> d13 = (s15 - s34)*spb12*spb34; d13 = T(1)/d13;
complex<T> d14 = spb12*spb34*square(s15 - s34); d14 = T(1)/d14;
complex<T> d15 = spb23*square(s15 - s23); d15 = T(1)/d15;
complex<T> d27 = spb15*square(s15 - s34)*T(2); d27 = T(1)/d27;
complex<T> d34 = (-s23 + s45)*spb12; d34 = T(1)/d34;
complex<T> d36 = spb12*square(s23 - s45); d36 = T(1)/d36;
complex<T> d38 = (s23 - s45)*spb12*spb23; d38 = T(1)/d38;
complex<T> d42 = spb12*spb23*spb45*T(6); d42 = T(1)/d42;
complex<T> d43 = s45; d43 = T(1)/d43;
complex<T> d44 = s45*spb23; d44 = T(1)/d44;
complex<T> d45 = spb23*spb45; d45 = T(1)/d45;
complex<T> d46 = spb23*cube(spb24); d46 = T(1)/d46;
complex<T> d48 = spb12*square(spb24); d48 = T(1)/d48;
complex<T> d49 = spb12*spb23*spb34; d49 = T(1)/d49;
complex<T> d50 = spb12*spb34*cube(spb24); d50 = T(1)/d50;
complex<T> d51 = spb15*cube(spb24); d51 = T(1)/d51;
complex<T> d52 = s45*spb12; d52 = T(1)/d52;
complex<T> d54 = spb12*spb15*square(spb24); d54 = T(1)/d54;
complex<T> d55 = spb12*spb15*spb34; d55 = T(1)/d55;
complex<T> d56 = spb12*spb45; d56 = T(1)/d56;
complex<T> d57 = spb12*spb15*spb34*cube(spb24); d57 = T(1)/d57;
complex<T> d60 = spb12*spb15*cube(spb24); d60 = T(1)/d60;
complex<T> d61 = spb12; d61 = T(1)/d61;
complex<T> d62 = spb12*spb23; d62 = T(1)/d62;
complex<T> d63 = spb15*spb34*T(2); d63 = T(1)/d63;
complex<T> d64 = spb15*spb34*spb45*T(2); d64 = T(1)/d64;
complex<T> d66 = spb34*square(s12 + s15 - s34)*T(2); d66 = T(1)/d66;
complex<T> d67 = spb15*cube(spb24)*T(2); d67 = T(1)/d67;
complex<T> d71 = s45*T(2); d71 = T(1)/d71;
complex<T> d74 = spb23*spb34*T(2); d74 = T(1)/d74;
complex<T> d75 = spb23*spb34*spb45*T(2); d75 = T(1)/d75;
complex<T> t38 = spb15*t2; 
complex<T> t48 = -t74; 
complex<T> t61 = d36*spa23; 
complex<T> t62 = d66*s25; 
complex<T> t75 = t24*t26; 
complex<T> t76 = -(spb45*t31); 
complex<T> t77 = t27*t28; 
complex<T> t88 = spb14*t25; 
complex<T> t99 = spb13*t23; 
complex<T> t100 = spb25*t24; 
complex<T> t102 = spb45*t30; 
complex<T> t103 = d15*spb34; 
complex<T> t110 = spb13*t74; 
complex<T> t114 = spb34*t4; 
complex<T> t117 = spa12*t108; 
complex<T> t120 = -t126; 
complex<T> t122 = t8*T(2); 
complex<T> t124 = d14*spb14; 
complex<T> t127 = d34*T(2); 
complex<T> t134 = spa24*t23; 
complex<T> t137 = spb23*t25; 
complex<T> t140 = t53*t74; 
complex<T> t150 = spb35*t78; 
complex<T> t167 = spb35*t32; 
complex<T> d1 = s45*t135; d1 = T(1)/d1;
complex<T> d2 = s45*spb23*t135; d2 = T(1)/d2;
complex<T> d3 = spb23*spb45*t135; d3 = T(1)/d3;
complex<T> d4 = t2*square(s15 - s34); d4 = T(1)/d4;
complex<T> d5 = (s15 - s23)*t11*t8; d5 = T(1)/d5;
complex<T> d8 = (s15 - s34)*t11*t8; d8 = T(1)/d8;
complex<T> d9 = spb23*t9; d9 = T(1)/d9;
complex<T> d10 = spb23*t36*t9; d10 = T(1)/d10;
complex<T> d11 = s24*spb15*t9; d11 = T(1)/d11;
complex<T> d18 = t8*t9; d18 = T(1)/d18;
complex<T> d19 = (s15 - s34)*spb34*t2; d19 = T(1)/d19;
complex<T> d22 = spb12*spb34*t8; d22 = T(1)/d22;
complex<T> d23 = spb34*t8*t9; d23 = T(1)/d23;
complex<T> d24 = (s15 - s34)*spb12*spb34*t8; d24 = T(1)/d24;
complex<T> d25 = (s15 - s23)*spb23*spb45*t2; d25 = T(1)/d25;
complex<T> d26 = spb34*spb45*t2*t8; d26 = T(1)/d26;
complex<T> d33 = t2*square(s23 - s45); d33 = T(1)/d33;
complex<T> d35 = s45*(-s23 + s45)*t2; d35 = T(1)/d35;
complex<T> d37 = s45*t2*square(s23 - s45); d37 = T(1)/d37;
complex<T> d39 = (s23 - s45)*spb12*spb23*t36; d39 = T(1)/d39;
complex<T> d40 = s45*t2; d40 = T(1)/d40;
complex<T> d41 = s45*spb23*t2; d41 = T(1)/d41;
complex<T> d47 = spb12*spb23*square(t36); d47 = T(1)/d47;
complex<T> d53 = spb12*square(t36); d53 = T(1)/d53;
complex<T> d58 = t8*cube(spb24); d58 = T(1)/d58;
complex<T> d59 = spb12*t8; d59 = T(1)/d59;
complex<T> d65 = spb23*t2*square(t36); d65 = T(1)/d65;
complex<T> d72 = t2*t8; d72 = T(1)/d72;
complex<T> d73 = spb23*spb34*t2; d73 = T(1)/d73;
complex<T> t17 = d47*s15*(s14*t140 + t30*t76) + spa15*(d48*spb25*t150 + d46*t77 + d50*t26*t88 - d49*t99); 
complex<T> t39 = d25*s14; 
complex<T> t42 = -(d18*spa34); 
complex<T> t43 = -t62; 
complex<T> t47 = d19*spa25; 
complex<T> t51 = -t88; 
complex<T> t52 = -t77; 
complex<T> t58 = d10*s14; 
complex<T> t59 = d65*s15; 
complex<T> t67 = -(d64*spb14*t69*t93) + d63*t93*t99; 
complex<T> t73 = -(d72*spa34); 
complex<T> t80 = d39*s14; 
complex<T> t81 = d23*s24; 
complex<T> t82 = d37*spa23; 
complex<T> t94 = d10*t30; 
complex<T> t97 = d33*spb15; 
complex<T> t104 = d26*t32; 
complex<T> t109 = d71*(s23*spa12*t120 + t110*t93); 
complex<T> t121 = t102*t31; 
complex<T> t131 = -(d1*T(3)); 
complex<T> t141 = t49*t77; 
complex<T> t145 = spb13*t48; 
complex<T> t147 = -(d2*T(3)); 
complex<T> t151 = d11*t100; 
complex<T> t155 = t110*T(2); 
complex<T> t157 = d3*T(3); 
complex<T> t158 = spa24*t103; 
complex<T> t162 = d35*t126; 
complex<T> t166 = t137*t24; 
complex<T> t195 = t117*t23; 
complex<T> d6 = s24*t122*square(s15 - s23); d6 = T(1)/d6;
complex<T> d7 = s24*t122*square(s15 - s34); d7 = T(1)/d7;
complex<T> d16 = t122*square(s15 - s23); d16 = T(1)/d16;
complex<T> d17 = (s15 - s34)*t122; d17 = T(1)/d17;
complex<T> d20 = spb34*t38*square(s15 - s34); d20 = T(1)/d20;
complex<T> d21 = (s15 - s23)*spb34*t38; d21 = T(1)/d21;
complex<T> d28 = spb34*t38*square(s15 - s23); d28 = T(1)/d28;
complex<T> d29 = t114*t2*square(s15 - s23); d29 = T(1)/d29;
complex<T> d30 = s24*t114*t9; d30 = T(1)/d30;
complex<T> d31 = t114*t2*square(s15 - s34); d31 = T(1)/d31;
complex<T> d32 = s24*(s15 - s34)*spb12*t114; d32 = T(1)/d32;
complex<T> d68 = t38*square(spb24); d68 = T(1)/d68;
complex<T> d69 = spb45*t38; d69 = T(1)/d69;
complex<T> d70 = t38*cube(spb24); d70 = T(1)/d70;
complex<T> t12 = t126*t131 + t110*t147 + t157*t23; 
complex<T> t14 = d61*t120 + d62*t145 + d62*spa45*t23 + d47*s45*t30*t76 + d47*t140*t96; 
complex<T> t16 = d52*s23*t126 + d53*s14*spa23*t140 + d52*spa23*t145 + d54*s23*spb25*t150 + d56*spa23*t23 + d51*spa23*t52 + d53*spa23*t30*t76 + d57*s23*t26*t88 + d55*spa23*t99; 
complex<T> t21 = d54*s34*spb25*t150 + d60*spa34*t26*t51 + d58*s34*t77 + d59*spa34*t99; 
complex<T> t22 = d68*s23*s34*spb25*t150 + d70*s23*spa34*t26*t51 + d67*s34*spa23*t52 - d69*spa23*spa34*spb14*t69; 
complex<T> t41 = -(d16*s24); 
complex<T> t85 = spa12*(d43*t120 + d44*t145 + d45*t23); 
complex<T> t98 = t195*t43 + d75*spa12*t69*t71 + d74*spa12*spa15*t99; 
complex<T> t111 = t51*t75; 
complex<T> t119 = t73*(t132*t69 + s45*t99); 
complex<T> t123 = t39*T(3); 
complex<T> t133 = t59*t96; 
complex<T> t136 = d28*spb13; 
complex<T> t160 = d21*t32; 
complex<T> t172 = t110*t82; 
complex<T> t13 = d38*t110 + t140*t61 + d39*t30*t76 + t140*t80 + spa12*spa24*spb23*t97 - d34*t126*T(2) + d41*t110*T(3) + d40*t126*T(3) + s23*t162*T(3) + s23*t172*T(3) + d42*t23*T(13); 
complex<T> t18 = d9*t110 + d29*t111 + d30*t111 + d39*t121 + t126*t127 + t126*t131 + d5*t141 + d6*t141 + d38*t145 + t110*t147 + t134*t160 + t136*t166 + t151*t167 + t123*t23 + t157*t23 + t134*t41 + t158*t48 + t162*t54 + t172*t54 + t140*t58 + t155*t61 + t104*t69 + t155*t80 + t76*t94 - spa12*spa24*spb23*t97 + t42*t99 + t81*t99; 
complex<T> t19 = d13*spa12*spb13*spb35 + d4*spb34*t100 + d31*t111 + d32*t111 + spa12*spb23*t124*t126 + d7*t141 + d8*t141 + d12*t100*t167 - d27*spb35*spb45*t24 + d24*s24*t99 - d20*s23*spa23*t99 + d22*t99*T(2) + d17*t134*T(3) - t23*t47*T(3); 
complex<T> t20 = -(d13*spa12*spb13*spb35) - d4*spb34*t100 + spa12*spb23*t120*t124 + d16*s24*t134 + d9*t145 + d12*t100*t150 + t150*t151 - t136*t166 + d27*spb35*spb45*t24 + d5*t49*t52 + d6*t49*t52 + d7*t49*t52 + d8*t49*t52 + t155*t58 + t104*t69 + t158*t74 + d21*t134*t78 + d29*t75*t88 + d30*t75*t88 + d31*t75*t88 + d32*t75*t88 + spb45*t31*t94 + d22*t99 - d24*s24*t99 + d20*s23*spa23*t99 + d18*spa34*t99 - t81*t99 - d17*t134*T(3) - t23*t39*T(3) + t23*t47*T(3); 
complex<T> t187 = t110*t133; 
complex<T> t15 = t187 + d73*spa45*t69*t71 + s45*t30*t59*t76 - d73*s45*spa15*t99; 
complex<T> co1 = -(t187*T(3)); 
complex<T> co2 = t195*t62; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t12*(*CI_users[0]->get_value(mc,ind,mu)) + t20*(*CI_users[1]->get_value(mc,ind,mu)) + t18*(*CI_users[2]->get_value(mc,ind,mu)) + t19*(*CI_users[3]->get_value(mc,ind,mu)) + t13*(*CI_users[4]->get_value(mc,ind,mu)) + t85*(*CI_users[5]->get_value(mc,ind,mu)) + t17*(*CI_users[6]->get_value(mc,ind,mu)) + t16*(*CI_users[7]->get_value(mc,ind,mu)) + t21*(*CI_users[8]->get_value(mc,ind,mu)) + t14*(*CI_users[9]->get_value(mc,ind,mu)) + t67*(*CI_users[10]->get_value(mc,ind,mu)) + co1*(*CI_users[11]->get_value(mc,ind,mu)) + co2*(*CI_users[12]->get_value(mc,ind,mu)) + t22*(*CI_users[13]->get_value(mc,ind,mu)) + t109*(*CI_users[14]->get_value(mc,ind,mu)) + t119*(*CI_users[15]->get_value(mc,ind,mu)) + t15*(*CI_users[16]->get_value(mc,ind,mu)) + t98*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmpqpQmQp_LLT_wCI::\
C2q2Q1g_qmpqpQmQp_LLT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmpqpQmQp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpqpQmQp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:03 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s35 = -(spa35*spb35);
complex<T> t2 = spa23*T(2); 
complex<T> t4 = square(s12 - s34); 
complex<T> t6 = spa12*spa15; 
complex<T> t9 = spa12*T(2); 
complex<T> t10 = spa15*spa25; 
complex<T> t22 = square(spa14); 
complex<T> t23 = square(spb25); 
complex<T> t24 = square(spa45); 
complex<T> t25 = square(spa12); 
complex<T> t26 = square(spa15); 
complex<T> t27 = square(spa24); 
complex<T> t28 = square(s25); 
complex<T> t29 = s15 - s34; 
complex<T> t30 = square(spa13); 
complex<T> t31 = s12 + s15 - s34; 
complex<T> t32 = -(spa35*T(3)); 
complex<T> t34 = square(spb35); 
complex<T> t37 = s34*s45; 
complex<T> t43 = cube(spb25); 
complex<T> t47 = -(spb35*T(3)); 
complex<T> t50 = s12*spb12; 
complex<T> t58 = cube(spa14); 
complex<T> t59 = spa12*spa34; 
complex<T> t61 = spa35*T(3); 
complex<T> t64 = spb23*spb34; 
complex<T> t65 = spa13*spa14; 
complex<T> t69 = spb35*T(3); 
complex<T> t71 = spa23*T(4); 
complex<T> t83 = spa15*spa34; 
complex<T> t102 = spa14*spb25; 
complex<T> t110 = spa35*spb45; 
complex<T> d7 = (-s12 + s45)*spa23; d7 = T(1)/d7;
complex<T> d18 = (s12 - s34)*spa12*spa23; d18 = T(1)/d18;
complex<T> d19 = (s12 - s45)*spa12*spa23; d19 = T(1)/d19;
complex<T> d21 = spa23*square(s12 - s45); d21 = T(1)/d21;
complex<T> d24 = (s12 - s34)*s35*spa12*spa23; d24 = T(1)/d24;
complex<T> d25 = s35*(s12 - s45)*spa12*spa23; d25 = T(1)/d25;
complex<T> d40 = spa12*spa23*spa45*T(6); d40 = T(1)/d40;
complex<T> d44 = spa23*spa35; d44 = T(1)/d44;
complex<T> d45 = spa23*spa45; d45 = T(1)/d45;
complex<T> d46 = spa23*square(spa35); d46 = T(1)/d46;
complex<T> d47 = s45*spa23; d47 = T(1)/d47;
complex<T> d54 = spa12*spa45; d54 = T(1)/d54;
complex<T> d55 = s45; d55 = T(1)/d55;
complex<T> d56 = s45*spa12; d56 = T(1)/d56;
complex<T> d57 = spa12*spa23*spa35; d57 = T(1)/d57;
complex<T> d58 = spa12*spa23*square(spa35); d58 = T(1)/d58;
complex<T> d63 = spa23; d63 = T(1)/d63;
complex<T> d64 = spa12*spa23; d64 = T(1)/d64;
complex<T> d67 = s45*T(2); d67 = T(1)/d67;
complex<T> d71 = spa15*spa24*T(2); d71 = T(1)/d71;
complex<T> t1 = square(s12 + t29); 
complex<T> t3 = square(t29); 
complex<T> t7 = spa34*(s12 + t29); 
complex<T> t40 = d71*spb23; 
complex<T> t44 = -t65; 
complex<T> t66 = t23*t25; 
complex<T> t68 = spa35*t24; 
complex<T> t75 = d21*spb12; 
complex<T> t82 = spa24*t23; 
complex<T> t84 = -(spa45*t34); 
complex<T> t90 = spb23*(d55*t102 - d54*t22 + d56*spb35*t65); 
complex<T> t92 = t26*t27; 
complex<T> t93 = spa45*t30; 
complex<T> t105 = spa13*t28; 
complex<T> t113 = spa13*t22; 
complex<T> t128 = spb25*t22; 
complex<T> t144 = spa13*t24; 
complex<T> d2 = spa12*spa45*t71; d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*spa12*spa45*t2; d3 = T(1)/d3;
complex<T> d4 = spa34*spa45*t2*t6; d4 = T(1)/d4;
complex<T> d5 = (s12 - s34)*spa23*t59; d5 = T(1)/d5;
complex<T> d6 = s45*t71; d6 = T(1)/d6;
complex<T> d8 = s45*(-s12 + s45)*t2; d8 = T(1)/d8;
complex<T> d9 = spa34*t4*t9; d9 = T(1)/d9;
complex<T> d10 = (s12 - s34)*t2*t83; d10 = T(1)/d10;
complex<T> d11 = t2*square(s12 - s45); d11 = T(1)/d11;
complex<T> d13 = t2*t4*t83; d13 = T(1)/d13;
complex<T> d14 = spa34*t10*t2*t4; d14 = T(1)/d14;
complex<T> d20 = s45*spa12*t71; d20 = T(1)/d20;
complex<T> d22 = s45*t2*square(s12 - s45); d22 = T(1)/d22;
complex<T> d23 = spa12*t4; d23 = T(1)/d23;
complex<T> d26 = spa23*spa34*t6; d26 = T(1)/d26;
complex<T> d29 = spa15*spa23*t29; d29 = T(1)/d29;
complex<T> d30 = spa15*t2*t29; d30 = T(1)/d30;
complex<T> d31 = spa34*t29*t9; d31 = T(1)/d31;
complex<T> d41 = s45*t2; d41 = T(1)/d41;
complex<T> d42 = s45*spa12*t2; d42 = T(1)/d42;
complex<T> d50 = spa34*cube(s12 + t29); d50 = T(1)/d50;
complex<T> d53 = t59*cube(s12 + t29); d53 = T(1)/d53;
complex<T> d62 = spa12*cube(s12 + t29); d62 = T(1)/d62;
complex<T> d65 = t83*T(2); d65 = T(1)/d65;
complex<T> d66 = spa45*t83*T(2); d66 = T(1)/d66;
complex<T> d70 = spa34*cube(s12 + t29)*T(2); d70 = T(1)/d70;
complex<T> d72 = t6*T(2); d72 = T(1)/d72;
complex<T> d73 = spa45*t6*T(2); d73 = T(1)/d73;
complex<T> d74 = spa12*spa35*t2; d74 = T(1)/d74;
complex<T> d75 = t2*t6; d75 = T(1)/d75;
complex<T> d76 = t2*t59; d76 = T(1)/d76;
complex<T> d77 = spa34*t2; d77 = T(1)/d77;
complex<T> d78 = spa34*spa45*t2; d78 = T(1)/d78;
complex<T> d79 = spa12*t2*square(spa35); d79 = T(1)/d79;
complex<T> t21 = d63*t102 - d64*spb45*t22 + d64*spb35*t65 + d58*s45*t93 - d57*s45*t65*T(2); 
complex<T> t38 = d3*s35; 
complex<T> t42 = d75*spb34; 
complex<T> t45 = -t68; 
complex<T> t57 = d76*spb15*(s45*t113 + t110*t58); 
complex<T> t67 = -t92; 
complex<T> t70 = spa34*t1; 
complex<T> t72 = d11*spb23; 
complex<T> t77 = d74*t37; 
complex<T> t81 = -t113; 
complex<T> t87 = -(d2*T(3)); 
complex<T> t104 = d23*spb35; 
complex<T> t116 = d6*T(3); 
complex<T> t118 = spb35*t75; 
complex<T> t120 = spa14*t82; 
complex<T> t127 = s34*t40; 
complex<T> t139 = t66*t68; 
complex<T> t142 = t43*t92; 
complex<T> t157 = d4*t61; 
complex<T> d1 = (s12 - s34)*spa23*t6*t7; d1 = T(1)/d1;
complex<T> d12 = (s12 - s34)*spa23*t7; d12 = T(1)/d12;
complex<T> d15 = (s12 - s34)*spa23*t10*t7; d15 = T(1)/d15;
complex<T> d16 = (s12 - s34)*t1*t59; d16 = T(1)/d16;
complex<T> d17 = t4*t7*t9; d17 = T(1)/d17;
complex<T> d27 = spa23*t29*t6*t7; d27 = T(1)/d27;
complex<T> d28 = t2*t3*t83; d28 = T(1)/d28;
complex<T> d32 = spa15*spa23*t3; d32 = T(1)/d32;
complex<T> d33 = t2*t3; d33 = T(1)/d33;
complex<T> d34 = spa23*t29*t7; d34 = T(1)/d34;
complex<T> d35 = spa34*t3*T(2); d35 = T(1)/d35;
complex<T> d36 = spa34*t10*t2*t3; d36 = T(1)/d36;
complex<T> d37 = spa23*t10*t29*t7; d37 = T(1)/d37;
complex<T> d38 = t1*t29*t59; d38 = T(1)/d38;
complex<T> d39 = t3*t7*t9; d39 = T(1)/d39;
complex<T> d43 = spa23*t1*t83; d43 = T(1)/d43;
complex<T> d51 = spa23*t1*t59; d51 = T(1)/d51;
complex<T> d59 = spa23*t1*t6; d59 = T(1)/d59;
complex<T> d60 = spa23*t1; d60 = T(1)/d60;
complex<T> d61 = spa23*t1*t10; d61 = T(1)/d61;
complex<T> t11 = t65*t77 + d79*t37*t93; 
complex<T> t12 = spb12*spb15*(d78*spa35*t58 + d77*t81); 
complex<T> t16 = d59*spb34*t105*t22 + d60*spb34*t120*t32 + d61*spb34*t45*t66 + d62*spb34*t43*t67 + d58*s34*t93 - d57*s34*t65*T(2); 
complex<T> t19 = d19*t47*t65 + d42*t47*t65 + d22*t47*t50*t65 - spb25*t59*t72 + d25*t30*t84 + d7*t102*T(2) + t118*t65*T(2) - d41*t102*T(3) - d8*s12*t102*T(3) - d40*t22*T(13); 
complex<T> t80 = spb23*(d67*s12*t102 + d67*spb12*spb35*t44 + d66*spa35*spb12*t58 + d65*spb12*t81); 
complex<T> t101 = s45*t113*t42 + t110*t42*t58 - t65*t77*T(3); 
complex<T> t117 = d32*spb23; 
complex<T> t119 = t102*t116 + d20*t65*t69 + t22*t87; 
complex<T> t172 = t127*t22; 
complex<T> d48 = spa23*t70; d48 = T(1)/d48;
complex<T> d49 = spa23*t10*t70; d49 = T(1)/d49;
complex<T> d52 = spa23*spa25*t70; d52 = T(1)/d52;
complex<T> d68 = t2*t70; d68 = T(1)/d68;
complex<T> d69 = spa25*t2*t70; d69 = T(1)/d69;
complex<T> t13 = -(spa12*spa35*t102*t117) + d36*t139 + d37*t139 + d38*t142 + d39*t142 + d35*spa14*spa45*t23 + d28*t113*t50 + d34*t120*t61 + d29*spb23*t65 + d27*t28*t81 - d33*spa15*t82 - d26*t113*T(2) + d31*t128*T(3) - d30*spb24*t22*T(3); 
complex<T> t15 = d69*s12*spb15*t139 + d70*s15*spb12*t142 + d68*spb12*spb15*t105*t22 + d68*s12*s15*t120*t32; 
complex<T> t18 = -(spa15*t102*t104) + d5*spb15*t113 + spa12*spa35*t102*t117 - d9*s25*t128 + d1*t105*t22 + d27*t105*t22 - d35*spa14*spa45*t23 + d13*spa12*t144*t23 + d12*t120*t32 + d34*t120*t32 + d29*spb23*t44 + t157*t58 + d10*t128*t61 + d18*t47*t65 + d14*t45*t66 + d15*t45*t66 + d36*t45*t66 + d37*t45*t66 + d16*t43*t67 + d17*t43*t67 + d38*t43*t67 + d39*t43*t67 + d26*t81 + d28*t50*t81 + d33*spa15*t82 + d24*t30*t84 - d31*t128*T(3) + d30*spb24*t22*T(3) - t22*t38*T(3); 
complex<T> t20 = spa15*t102*t104 + t119 + d9*s25*t128 + d14*t139 + d15*t139 + d16*t142 + d17*t142 - d13*spa12*t144*t23 + d10*t128*t32 + t157*t58 + d12*t120*t61 + d18*t65*t69 + d19*t65*t69 + d22*t50*t65*t69 + spb25*t59*t72 + d5*spb15*t81 + d1*t28*t81 + d24*t34*t93 + d25*t34*t93 - d7*t102*T(2) - t118*t65*T(2) + d8*s12*t102*T(3) + t22*t38*T(3); 
complex<T> t112 = t172 + d73*spa35*t58*t64 + d72*t64*t81; 
complex<T> t125 = d48*t32; 
complex<T> t14 = s15*t120*t125 + d52*spb15*t139 + d53*s15*t43*t67 + d51*spb15*t28*t81; 
complex<T> t17 = d47*(-(s12*t102) + spb12*spb35*t65) + s12*(t120*t125 + d49*t45*t66) + spb12*(d50*t142 - d45*t22 + d43*t28*t81 + d46*t93 - d44*t65*T(2)); 
complex<T> co1 = -t172; 
complex<T> co2 = Complex(0,1); 
SeriesC<T> result = co2*(t20*(*CI_users[0]->get_value(mc,ind,mu)) + t13*(*CI_users[1]->get_value(mc,ind,mu)) + t119*(*CI_users[2]->get_value(mc,ind,mu)) + t18*(*CI_users[3]->get_value(mc,ind,mu)) + t19*(*CI_users[4]->get_value(mc,ind,mu)) + t17*(*CI_users[5]->get_value(mc,ind,mu)) + t14*(*CI_users[6]->get_value(mc,ind,mu)) + t90*(*CI_users[7]->get_value(mc,ind,mu)) + t16*(*CI_users[8]->get_value(mc,ind,mu)) + t21*(*CI_users[9]->get_value(mc,ind,mu)) + t80*(*CI_users[10]->get_value(mc,ind,mu)) + t15*(*CI_users[11]->get_value(mc,ind,mu)) + t112*(*CI_users[12]->get_value(mc,ind,mu)) + t101*(*CI_users[13]->get_value(mc,ind,mu)) + co1*(*CI_users[14]->get_value(mc,ind,mu)) + t57*(*CI_users[15]->get_value(mc,ind,mu)) + t12*(*CI_users[16]->get_value(mc,ind,mu)) + t11*(*CI_users[17]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmmqpQpQm_LLT_wCI::\
C2q2Q1g_qmmqpQpQm_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmmqpQpQm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, qp, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmmqpQpQm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:14 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb13 = SPB(1,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> s15 = S(1,5);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = spb12*T(2); 
complex<T> t7 = square(spb34); 
complex<T> t10 = -(spb34*T(3)); 
complex<T> t14 = spb12*T(4); 
complex<T> t19 = spa15*spb13; 
complex<T> t20 = spb34*T(3); 
complex<T> t21 = s23*spa25; 
complex<T> t22 = spa12*spb14; 
complex<T> t34 = spa25*spb34; 
complex<T> d5 = (-s23 + s45)*spb12; d5 = T(1)/d5;
complex<T> d7 = spb12*square(s23 - s45); d7 = T(1)/d7;
complex<T> d13 = spb12*spb23*spb45*T(3); d13 = T(1)/d13;
complex<T> d14 = s45; d14 = T(1)/d14;
complex<T> d15 = s45*spb23; d15 = T(1)/d15;
complex<T> d16 = spb23*spb45; d16 = T(1)/d16;
complex<T> d17 = s45*spb12; d17 = T(1)/d17;
complex<T> d18 = spb12*spb45; d18 = T(1)/d18;
complex<T> d19 = spb12; d19 = T(1)/d19;
complex<T> d20 = spb12*spb23; d20 = T(1)/d20;
complex<T> d21 = spb45*T(2); d21 = T(1)/d21;
complex<T> d23 = s45*T(2); d23 = T(1)/d23;
complex<T> d25 = spb23*spb45*T(2); d25 = T(1)/d25;
complex<T> t11 = -t19; 
complex<T> t12 = -t21; 
complex<T> t15 = d7*spa23; 
complex<T> t27 = spb34*t19; 
complex<T> t30 = d23*spa12; 
complex<T> t42 = spa12*t7; 
complex<T> t45 = spa45*t7; 
complex<T> t50 = spa23*t7; 
complex<T> d1 = s45*t14; d1 = T(1)/d1;
complex<T> d2 = s45*spb23*t14; d2 = T(1)/d2;
complex<T> d3 = spb23*spb45*t14; d3 = T(1)/d3;
complex<T> d4 = t1*square(s23 - s45); d4 = T(1)/d4;
complex<T> d6 = s45*(-s23 + s45)*t1; d6 = T(1)/d6;
complex<T> d8 = s45*t1*square(s23 - s45); d8 = T(1)/d8;
complex<T> d9 = (s23 - s45)*spb23*t1; d9 = T(1)/d9;
complex<T> d10 = (s23 - s45)*spb45*t1; d10 = T(1)/d10;
complex<T> d11 = s45*t1; d11 = T(1)/d11;
complex<T> d12 = s45*spb23*t1; d12 = T(1)/d12;
complex<T> d22 = spb45*t1; d22 = T(1)/d22;
complex<T> d24 = spb23*t1; d24 = T(1)/d24;
complex<T> t18 = spb34*(spa23*t11 + t21)*t30; 
complex<T> t23 = d4*spb23; 
complex<T> t25 = d8*spa23; 
complex<T> t26 = d15*spa12*t27 + d14*spa12*t34 + d16*t42; 
complex<T> t31 = d1*spa25; 
complex<T> t33 = d17*spb34*t12 + d17*spa23*t27 + d18*t50; 
complex<T> t36 = d2*t19; 
complex<T> t38 = d19*t34 + d20*(t27 + t45); 
complex<T> t3 = d11*spa25*t10 + d12*t10*t19 + d9*t10*t19 + d6*t10*t21 + d10*t10*t22 - spa25*t22*t23 + s23*t10*t19*t25 + t15*t27*T(2) + d5*t34*T(2) + d13*t7*T(2); 
complex<T> t4 = d9*t19*t20 + d6*t20*t21 + d10*t20*t22 + spa25*t22*t23 + s23*t19*t20*t25 + t20*t31 + t20*t36 - t15*t27*T(2) - d5*t34*T(2) - d3*t7*T(3); 
complex<T> t5 = t20*(t31 + t36) + d3*t7*T(3); 
complex<T> co1 = -(d21*spa23*t42); 
complex<T> co2 = d22*s34*t50; 
complex<T> co3 = d24*s34*t45; 
complex<T> co4 = d24*s15*t45; 
complex<T> co5 = d25*s15*t42; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)) + t26*(*CI_users[3]->get_value(mc,ind,mu)) + t33*(*CI_users[4]->get_value(mc,ind,mu)) + t38*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t18*(*CI_users[8]->get_value(mc,ind,mu)) + co3*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)) + co5*(*CI_users[11]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmpqpQpQm_LLT_wCI::\
C2q2Q1g_qmpqpQpQm_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmpqpQpQm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpqpQpQm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:26 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> t1 = spa23*T(2); 
complex<T> t7 = square(spa15); 
complex<T> t10 = spb12*spb23; 
complex<T> t11 = -(spa15*T(3)); 
complex<T> t19 = spa13*spb34; 
complex<T> t21 = spa35*spb23; 
complex<T> t22 = s12*spb24; 
complex<T> t23 = spa23*T(4); 
complex<T> d4 = (-s12 + s45)*spa23; d4 = T(1)/d4;
complex<T> d9 = spa23*square(s12 - s45); d9 = T(1)/d9;
complex<T> d11 = spa12*spa23*spa45*T(3); d11 = T(1)/d11;
complex<T> d14 = spa23*spa45; d14 = T(1)/d14;
complex<T> d15 = s45*spa23; d15 = T(1)/d15;
complex<T> d16 = spa12*spa45; d16 = T(1)/d16;
complex<T> d17 = s45; d17 = T(1)/d17;
complex<T> d18 = s45*spa12; d18 = T(1)/d18;
complex<T> d19 = spa23; d19 = T(1)/d19;
complex<T> d20 = spa12*spa23; d20 = T(1)/d20;
complex<T> d21 = spa45*T(2); d21 = T(1)/d21;
complex<T> d22 = s45*T(2); d22 = T(1)/d22;
complex<T> d23 = spa12*spa45*T(2); d23 = T(1)/d23;
complex<T> t12 = -t19; 
complex<T> t14 = -t22; 
complex<T> t25 = d9*spb12; 
complex<T> t45 = -(spb45*t7); 
complex<T> t49 = -(spb12*t7); 
complex<T> t50 = -(spb23*t7); 
complex<T> d1 = spa12*spa45*t23; d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa45*t1; d2 = T(1)/d2;
complex<T> d3 = s45*t23; d3 = T(1)/d3;
complex<T> d5 = s45*(-s12 + s45)*t1; d5 = T(1)/d5;
complex<T> d6 = t1*square(s12 - s45); d6 = T(1)/d6;
complex<T> d7 = (s12 - s45)*spa12*t1; d7 = T(1)/d7;
complex<T> d8 = s45*spa12*t23; d8 = T(1)/d8;
complex<T> d10 = s45*t1*square(s12 - s45); d10 = T(1)/d10;
complex<T> d12 = s45*t1; d12 = T(1)/d12;
complex<T> d13 = s45*spa12*t1; d13 = T(1)/d13;
complex<T> d24 = spa12*t1; d24 = T(1)/d24;
complex<T> d25 = spa45*t1; d25 = T(1)/d25;
complex<T> t16 = d10*spb12; 
complex<T> t18 = d15*spa15*(spb12*t12 + t22) + d14*t49; 
complex<T> t24 = d6*spa12; 
complex<T> t26 = -(d17*spa15*spb23*spb24) + d18*spa15*spb23*t12 + d16*t50; 
complex<T> t29 = d3*spb24; 
complex<T> t33 = -(d19*spa15*spb24) + d20*(spa15*t12 + t45); 
complex<T> t36 = d8*t19; 
complex<T> t37 = d22*spa15*(spb23*t14 + t10*t19) + d21*t10*t7; 
complex<T> t3 = t11*(t29 + t36) - d1*t7*T(3); 
complex<T> t30 = s12*t16; 
complex<T> t4 = d7*t11*t19 + d2*t11*t21 + d5*t11*t22 - spb24*t21*t24 + t11*t29 + t11*t19*t30 + t11*t36 + d4*spa15*spb24*T(2) + spa15*t19*t25*T(2) + d1*t7*T(3); 
complex<T> t5 = spb24*t21*t24 - d4*spa15*spb24*T(2) - spa15*t19*t25*T(2) - d11*t7*T(2) + d12*spa15*spb24*T(3) + d13*spa15*t19*T(3) + d7*spa15*t19*T(3) + d2*spa15*t21*T(3) + d5*spa15*t22*T(3) + spa15*t19*t30*T(3); 
complex<T> co1 = d23*s34*t50; 
complex<T> co2 = d24*s34*t45; 
complex<T> co3 = d24*s15*t45; 
complex<T> co4 = d25*s15*t49; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + t18*(*CI_users[3]->get_value(mc,ind,mu)) + t26*(*CI_users[4]->get_value(mc,ind,mu)) + t33*(*CI_users[5]->get_value(mc,ind,mu)) + t37*(*CI_users[6]->get_value(mc,ind,mu)) + co1*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu)) + co3*(*CI_users[9]->get_value(mc,ind,mu)) + co4*(*CI_users[10]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpmQmQp_LRT_wCI::\
C2q2Q1g_qmqpmQmQp_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpmQmQp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpmQmQp LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:28 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb45 = SPB(4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> t6 = square(spb25); 
complex<T> t8 = square(spa13); 
complex<T> d1 = (s12 - s45)*spb12*spb34; d1 = T(1)/d1;
complex<T> d2 = spb34*square(s12 - s45)*T(2); d2 = T(1)/d2;
complex<T> d3 = spb12*spb23*spb34*spb45*T(2); d3 = T(1)/d3;
complex<T> d4 = spb34*spb35; d4 = T(1)/d4;
complex<T> d5 = spb12*spb35; d5 = T(1)/d5;
complex<T> d6 = spb12*spb23*spb34; d6 = T(1)/d6;
complex<T> d7 = spb12*spb34*spb35; d7 = T(1)/d7;
complex<T> d8 = spb12*spb35*T(2); d8 = T(1)/d8;
complex<T> t4 = d1*spa34*spb24*spb25 + d1*spa35*t6 - d2*spb12*spb35*t8; 
complex<T> t5 = -(d1*(spa34*spb24*spb25 + spa35*t6)) + d2*spb12*spb35*t8 - d3*spb24*t6*T(3); 
complex<T> t18 = s45*t6; 
complex<T> t3 = -(d7*t18) - d6*spa45*spb24*t6; 
complex<T> co1 = -(d4*spa12*t6); 
complex<T> co2 = d5*spa34*t6; 
complex<T> co3 = d8*spa34*t18; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)) + t3*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqppQmQp_LRT_wCI::\
C2q2Q1g_qmqppQmQp_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqppQmQp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqppQmQp LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:31 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t2 = spa12*spa34; 
complex<T> t4 = square(s35); 
complex<T> t8 = square(spa13); 
complex<T> t9 = square(spa45); 
complex<T> t10 = cube(spb35); 
complex<T> t11 = square(spa14); 
complex<T> t23 = spb25*spb35; 
complex<T> d2 = spa34*square(s12 - s34); d2 = T(1)/d2;
complex<T> d8 = spa34*cube(spa35); d8 = T(1)/d8;
complex<T> d9 = spa12*cube(spa35); d9 = T(1)/d9;
complex<T> d12 = spa12*cube(spa35)*T(2); d12 = T(1)/d12;
complex<T> t17 = t2*T(2); 
complex<T> t22 = t8*t9; 
complex<T> d3 = (s12 - s34)*t2*t4; d3 = T(1)/d3;
complex<T> d6 = (s12 - s45)*t2*t4; d6 = T(1)/d6;
complex<T> d10 = t2*cube(spa35); d10 = T(1)/d10;
complex<T> d11 = spa23*t2; d11 = T(1)/d11;
complex<T> t7 = d11*spa24*spb45*t11 + d10*s45*t22; 
complex<T> t16 = -t22; 
complex<T> d1 = t17*square(s12 - s34); d1 = T(1)/d1;
complex<T> d4 = s35*t17*square(s12 - s34); d4 = T(1)/d4;
complex<T> d5 = s35*t17*square(s12 - s45); d5 = T(1)/d5;
complex<T> d7 = spa23*spa45*t17; d7 = T(1)/d7;
complex<T> t6 = (d5 + d6)*t10*t22 + d7*spa24*t11*T(3); 
complex<T> t18 = d1*s35; 
complex<T> t5 = -(spb35*t11*t18) + d3*t10*t22 + d4*t10*t22 - d2*spa14*spa45*t23; 
complex<T> t21 = d3*t10*t16 + d4*t10*t16 + d5*t10*t16 + d6*t10*t16 + spb35*t11*t18 + d2*spa14*spa45*t23; 
complex<T> co1 = d8*spb12*t22; 
complex<T> co2 = d9*spb34*t16; 
complex<T> co3 = d12*s45*spb34*t16; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t21*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)) + t7*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpmQpQm_LRT_wCI::\
C2q2Q1g_qmqpmQpQm_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpmQpQm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpmQpQm LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:34 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t2 = spb12*spb34; 
complex<T> t4 = square(s35); 
complex<T> t8 = square(spb23); 
complex<T> t9 = square(spb45); 
complex<T> t10 = cube(spa35); 
complex<T> t11 = square(spb24); 
complex<T> t30 = spa35*spb24; 
complex<T> d2 = spb34*square(s12 - s34); d2 = T(1)/d2;
complex<T> d8 = spb34*cube(spb35); d8 = T(1)/d8;
complex<T> d9 = spb12*cube(spb35); d9 = T(1)/d9;
complex<T> d12 = spb12*cube(spb35)*T(2); d12 = T(1)/d12;
complex<T> t18 = t2*T(2); 
complex<T> t23 = t8*t9; 
complex<T> d3 = (s12 - s34)*t2*t4; d3 = T(1)/d3;
complex<T> d6 = (s12 - s45)*t2*t4; d6 = T(1)/d6;
complex<T> d10 = spb23*t2; d10 = T(1)/d10;
complex<T> d11 = t2*cube(spb35); d11 = T(1)/d11;
complex<T> t16 = -t23; 
complex<T> d1 = t18*square(s12 - s34); d1 = T(1)/d1;
complex<T> d4 = s35*t18*square(s12 - s34); d4 = T(1)/d4;
complex<T> d5 = s35*t18*square(s12 - s45); d5 = T(1)/d5;
complex<T> d7 = spb23*spb45*t18; d7 = T(1)/d7;
complex<T> t5 = d11*s45*t16 - d10*spa45*cube(spb24); 
complex<T> t6 = d1*s35*spa35*t11 + d3*t10*t16 + d4*t10*t16 - d2*spa15*spb45*t30; 
complex<T> t7 = (d5 + d6)*t10*t16 - d7*cube(spb24)*T(3); 
complex<T> t15 = -(d1*s35); 
complex<T> t22 = spa35*t11*t15 + d3*t10*t23 + d4*t10*t23 + d5*t10*t23 + d6*t10*t23 + d2*spa15*spb45*t30; 
complex<T> co1 = d8*spa12*t16; 
complex<T> co2 = d9*spa34*t23; 
complex<T> co3 = d12*s45*spa34*t23; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t22*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t7*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)) + t5*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqppQpQm_LRT_wCI::\
C2q2Q1g_qmqppQpQm_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqppQpQm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqppQpQm LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:35 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t5 = square(spa15); 
complex<T> t8 = square(spb23); 
complex<T> d1 = (s12 - s45)*spa34; d1 = T(1)/d1;
complex<T> d2 = spa34*square(s12 - s45)*T(2); d2 = T(1)/d2;
complex<T> d3 = spa12*spa23*spa34*spa45*T(2); d3 = T(1)/d3;
complex<T> d4 = spa34*spa35; d4 = T(1)/d4;
complex<T> d5 = spa12*spa35; d5 = T(1)/d5;
complex<T> d6 = spa12*spa34*spa35; d6 = T(1)/d6;
complex<T> d7 = spa12*spa23*spa34; d7 = T(1)/d7;
complex<T> d8 = spa12*spa35*T(2); d8 = T(1)/d8;
complex<T> t2 = d1*spa15*spb23 - d2*spa12*spa35*t8 + d3*spa24*t5*T(3); 
complex<T> t3 = -(d1*spa15*spb23) + d2*spa12*spa35*t8; 
complex<T> t13 = s45*t5; 
complex<T> t4 = d6*t13 + d7*spa24*spb45*t5; 
complex<T> co1 = d4*spb12*t5; 
complex<T> co2 = -(d5*spb34*t5); 
complex<T> co3 = -(d8*spb34*t13); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + t4*(*CI_users[4]->get_value(mc,ind,mu)) + co3*(*CI_users[5]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmmQp_LRT_wCI::\
C2q2Q1g_qmqpQmmQp_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmmQp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, m, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmmQp LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:37 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = square(s12 - s45); 
complex<T> t4 = square(spb25); 
complex<T> t16 = spa34*spb23; 
complex<T> t17 = s12*spb25; 
complex<T> d4 = spb12*spb34*spb45*T(2); d4 = T(1)/d4;
complex<T> d5 = spb34*spb45; d5 = T(1)/d5;
complex<T> d6 = spb12*T(2); d6 = T(1)/d6;
complex<T> t12 = spa12*t4; 
complex<T> d1 = spb34*t1; d1 = T(1)/d1;
complex<T> d2 = spb34*t1*T(2); d2 = T(1)/d2;
complex<T> d3 = spb12*spb34*t1*T(2); d3 = T(1)/d3;
complex<T> t8 = d3*spa45; 
complex<T> t15 = d2*spa14; 
complex<T> t2 = -(spb45*t15*t16) - d1*spa45*t12*T(2) + d1*spa14*t17*T(2) - s45*spb25*t15*T(3) - d4*t4*T(3) - s45*t4*t8*T(3); 
complex<T> t3 = spb45*t15*t16 + d1*spa45*t12*T(2) - d1*spa14*t17*T(2) + s45*spb25*t15*T(3) + s45*t4*t8*T(3); 
complex<T> co1 = -(d5*t12); 
complex<T> co2 = -(d6*spa34*spa45*t4); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmpQp_LRT_wCI::\
C2q2Q1g_qmqpQmpQp_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmpQp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, p, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmpQp LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:39 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa15 = SPA(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> t1 = square(s12 - s34); 
complex<T> t4 = square(spa13); 
complex<T> t16 = spa15*spa34; 
complex<T> t17 = s12*spa13; 
complex<T> d1 = spa12*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d5 = spa34*spa45; d5 = T(1)/d5;
complex<T> d6 = spa12*T(2); d6 = T(1)/d6;
complex<T> t12 = spb12*t4; 
complex<T> d2 = spa45*t1; d2 = T(1)/d2;
complex<T> d3 = spa45*t1*T(2); d3 = T(1)/d3;
complex<T> d4 = spa12*spa45*t1*T(2); d4 = T(1)/d4;
complex<T> t8 = d4*spb34; 
complex<T> t15 = d3*spb24; 
complex<T> t2 = -(spb45*t15*t16) - d2*spb34*t12*T(2) + d2*spb24*t17*T(2) - s34*spa13*t15*T(3) - s34*t4*t8*T(3); 
complex<T> t3 = spb45*t15*t16 + d2*spb34*t12*T(2) - d2*spb24*t17*T(2) + s34*spa13*t15*T(3) + d1*t4*T(3) + s34*t4*t8*T(3); 
complex<T> co1 = d5*t12; 
complex<T> co2 = d6*spb34*spb45*t4; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpmQm_LRT_wCI::\
C2q2Q1g_qmqpQpmQm_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpmQm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, m, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpmQm LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:41 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa14 = SPA(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> t1 = square(s12 - s34); 
complex<T> t4 = square(spb23); 
complex<T> t16 = spa45*spb25; 
complex<T> t17 = s12*spb23; 
complex<T> d4 = spb12*spb34*spb45*T(2); d4 = T(1)/d4;
complex<T> d5 = spb34*spb45; d5 = T(1)/d5;
complex<T> d6 = spb12*T(2); d6 = T(1)/d6;
complex<T> t12 = spa12*t4; 
complex<T> d1 = spb45*t1; d1 = T(1)/d1;
complex<T> d2 = spb45*t1*T(2); d2 = T(1)/d2;
complex<T> d3 = spb12*spb45*t1*T(2); d3 = T(1)/d3;
complex<T> t8 = d3*spa34; 
complex<T> t15 = d2*spa14; 
complex<T> t2 = -(spb34*t15*t16) + d1*(spa34*t12 + spa14*t17)*T(2) - s34*spb23*t15*T(3) + s34*t4*t8*T(3); 
complex<T> t3 = spb34*t15*t16 - d1*(spa34*t12 + spa14*t17)*T(2) + s34*spb23*t15*T(3) - d4*t4*T(3) - s34*t4*t8*T(3); 
complex<T> co1 = -(d5*t12); 
complex<T> co2 = -(d6*spa34*spa45*t4); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQppQm_LRT_wCI::\
C2q2Q1g_qmqpQppQm_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQppQm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, p, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQppQm LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:43 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa13 = SPA(1,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = square(s12 - s45); 
complex<T> t4 = square(spa15); 
complex<T> t16 = spa13*spa45; 
complex<T> t17 = s12*spa15; 
complex<T> d1 = spa12*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d5 = spa34*spa45; d5 = T(1)/d5;
complex<T> d6 = spa12*T(2); d6 = T(1)/d6;
complex<T> t12 = spb12*t4; 
complex<T> d2 = spa34*t1; d2 = T(1)/d2;
complex<T> d3 = spa34*t1*T(2); d3 = T(1)/d3;
complex<T> d4 = spa12*spa34*t1*T(2); d4 = T(1)/d4;
complex<T> t8 = d4*spb45; 
complex<T> t15 = d3*spb24; 
complex<T> t2 = spb34*t15*t16 - d2*(spb45*t12 + spb24*t17)*T(2) + s45*spa15*t15*T(3) - s45*t4*t8*T(3); 
complex<T> t3 = -(spb34*t15*t16) + d2*(spb45*t12 + spb24*t17)*T(2) - s45*spa15*t15*T(3) + d1*t4*T(3) + s45*t4*t8*T(3); 
complex<T> co1 = d5*t12; 
complex<T> co2 = d6*spb34*spb45*t4; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmQpm_LRT_wCI::\
C2q2Q1g_qmqpQmQpm_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmQpm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, m}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmQpm LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:45 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> t2 = spb12*spb45; 
complex<T> t4 = square(s35); 
complex<T> t8 = square(spb25); 
complex<T> t9 = square(spb34); 
complex<T> t10 = cube(spa35); 
complex<T> t11 = square(spb24); 
complex<T> t23 = spa35*spb24; 
complex<T> d2 = spb45*square(s12 - s45); d2 = T(1)/d2;
complex<T> d8 = spb45*cube(spb35); d8 = T(1)/d8;
complex<T> d11 = spb12*cube(spb35); d11 = T(1)/d11;
complex<T> d12 = spb12*cube(spb35)*T(2); d12 = T(1)/d12;
complex<T> t17 = t2*T(2); 
complex<T> t19 = -(d2*spa13); 
complex<T> t22 = t8*t9; 
complex<T> d3 = (s12 - s34)*t2*t4; d3 = T(1)/d3;
complex<T> d6 = (s12 - s45)*t2*t4; d6 = T(1)/d6;
complex<T> d9 = spb15*t2; d9 = T(1)/d9;
complex<T> d10 = t2*cube(spb35); d10 = T(1)/d10;
complex<T> t16 = -t22; 
complex<T> d1 = t17*square(s12 - s45); d1 = T(1)/d1;
complex<T> d4 = s35*t17*square(s12 - s34); d4 = T(1)/d4;
complex<T> d5 = s35*t17*square(s12 - s45); d5 = T(1)/d5;
complex<T> d7 = spb15*spb34*t17; d7 = T(1)/d7;
complex<T> t5 = -(d9*spa34*spb14*t11) + d10*s34*t16; 
complex<T> t6 = (d3 + d4)*t10*t16 - d7*spb14*t11*T(3); 
complex<T> t7 = d1*s35*spa35*t11 + d5*t10*t16 + d6*t10*t16 + d2*spa13*spb34*t23; 
complex<T> t15 = -(d1*s35); 
complex<T> t21 = spa35*t11*t15 + d3*t10*t22 + d4*t10*t22 + d5*t10*t22 + d6*t10*t22 + spb34*t19*t23; 
complex<T> co1 = d8*spa12*t16; 
complex<T> co2 = d11*spa45*t22; 
complex<T> co3 = d12*s34*spa45*t22; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t21*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t7*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmQpp_LRT_wCI::\
C2q2Q1g_qmqpQmQpp_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmQpp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, p}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmQpp LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:47 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> t6 = square(spa13); 
complex<T> t9 = square(spb25); 
complex<T> d1 = spa45*square(s12 - s34)*T(2); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa12*spa45; d2 = T(1)/d2;
complex<T> d3 = spa12*spa15*spa34*spa45*T(2); d3 = T(1)/d3;
complex<T> d4 = spa35*spa45; d4 = T(1)/d4;
complex<T> d5 = spa12*spa35*spa45; d5 = T(1)/d5;
complex<T> d6 = spa12*spa15*spa45; d6 = T(1)/d6;
complex<T> d7 = spa12*spa35; d7 = T(1)/d7;
complex<T> d8 = spa12*spa35*T(2); d8 = T(1)/d8;
complex<T> t4 = -(d2*(spa13*spa14*spb45 + spb35*t6)) + d1*spa12*spa35*t9; 
complex<T> t5 = d2*spa13*spa14*spb45 + d2*spb35*t6 - d1*spa12*spa35*t9 + d3*spa14*t6*T(3); 
complex<T> t17 = s34*t6; 
complex<T> t3 = d5*t17 + d6*spa14*spb34*t6; 
complex<T> co1 = d4*spb12*t6; 
complex<T> co2 = -(d7*spb45*t6); 
complex<T> co3 = -(d8*spb45*t17); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + t3*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpQmm_LRT_wCI::\
C2q2Q1g_qmqpQpQmm_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpQmm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, m}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpQmm LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:49 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> t5 = square(spb23); 
complex<T> t8 = square(spa15); 
complex<T> d1 = (s12 - s34)*spb45; d1 = T(1)/d1;
complex<T> d2 = spb45*square(s12 - s34)*T(2); d2 = T(1)/d2;
complex<T> d3 = spb12*spb15*spb34*spb45*T(2); d3 = T(1)/d3;
complex<T> d4 = spb35*spb45; d4 = T(1)/d4;
complex<T> d5 = spb12*spb15*spb45; d5 = T(1)/d5;
complex<T> d6 = spb12*spb35*spb45; d6 = T(1)/d6;
complex<T> d7 = spb12*spb35; d7 = T(1)/d7;
complex<T> d8 = spb12*spb35*T(2); d8 = T(1)/d8;
complex<T> t3 = d1*spa15*spb23 - d2*spb12*spb35*t8; 
complex<T> t4 = -(d1*spa15*spb23) + d2*spb12*spb35*t8 - d3*spb14*t5*T(3); 
complex<T> t13 = s34*t5; 
complex<T> t2 = -(d6*t13) - d5*spa34*spb14*t5; 
complex<T> co1 = -(d4*spa12*t5); 
complex<T> co2 = d7*spa45*t5; 
complex<T> co3 = d8*spa45*t13; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)) + co3*(*CI_users[5]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpQmp_LRT_wCI::\
C2q2Q1g_qmqpQpQmp_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpQmp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, p}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpQmp LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:52 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> t2 = spa12*spa45; 
complex<T> t4 = square(s35); 
complex<T> t8 = square(spa15); 
complex<T> t9 = square(spa34); 
complex<T> t10 = cube(spb35); 
complex<T> t11 = square(spa14); 
complex<T> t30 = spa14*spa34; 
complex<T> d2 = spa45*square(s12 - s45); d2 = T(1)/d2;
complex<T> d8 = spa45*cube(spa35); d8 = T(1)/d8;
complex<T> d11 = spa12*cube(spa35); d11 = T(1)/d11;
complex<T> d12 = spa12*cube(spa35)*T(2); d12 = T(1)/d12;
complex<T> t18 = t2*T(2); 
complex<T> t23 = t8*t9; 
complex<T> d3 = (s12 - s34)*t2*t4; d3 = T(1)/d3;
complex<T> d6 = (s12 - s45)*t2*t4; d6 = T(1)/d6;
complex<T> d9 = t2*cube(spa35); d9 = T(1)/d9;
complex<T> d10 = spa15*t2; d10 = T(1)/d10;
complex<T> t5 = d9*s34*t23 + d10*spb34*cube(spa14); 
complex<T> t16 = -t23; 
complex<T> d1 = t18*square(s12 - s45); d1 = T(1)/d1;
complex<T> d4 = s35*t18*square(s12 - s34); d4 = T(1)/d4;
complex<T> d5 = s35*t18*square(s12 - s45); d5 = T(1)/d5;
complex<T> d7 = spa15*spa34*t18; d7 = T(1)/d7;
complex<T> t6 = (d3 + d4)*t10*t23 + d7*cube(spa14)*T(3); 
complex<T> t7 = d1*s35*spb35*t11 + d3*t10*t16 + d4*t10*t16 + d5*t10*t16 + d6*t10*t16 - d2*spb23*spb35*t30; 
complex<T> t15 = -(d1*s35); 
complex<T> t22 = (d5 + d6)*t10*t23 + spb35*(t11*t15 + d2*spb23*t30); 
complex<T> co1 = d8*spb12*t23; 
complex<T> co2 = d11*spb45*t16; 
complex<T> co3 = d12*s34*spb45*t16; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t7*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t22*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmmqpQmQp_LRT_wCI::\
C2q2Q1g_qmmqpQmQp_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmmqpQmQp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, qp, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmmqpQmQp LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:53 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb35); 
complex<T> d1 = spb12*spb23*spb45*T(2); d1 = T(1)/d1;
complex<T> d2 = spb12*spb23; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d2*spa45*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmpqpQmQp_LRT_wCI::\
C2q2Q1g_qmpqpQmQp_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmpqpQmQp_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpqpQmQp LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:54 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa23*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = spa12*spa23; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb45*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmmqpQpQm_LRT_wCI::\
C2q2Q1g_qmmqpQpQm_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmmqpQpQm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, qp, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmmqpQpQm LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:55 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb34); 
complex<T> d1 = spb12*spb23*spb45*T(2); d1 = T(1)/d1;
complex<T> d2 = spb12*spb23; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d2*spa45*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmpqpQpQm_LRT_wCI::\
C2q2Q1g_qmpqpQpQm_LRT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmpqpQpQm_LRT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpqpQpQm LRT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:56 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spa15); 
complex<T> d1 = spa12*spa23*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = spa12*spa23; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb45*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmmQmQpqp_LLT_wCI::\
C2q2Q1g_qmmQmQpqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmmQmQpqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, Qm, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmmQmQpqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:56:58 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> t6 = square(spb45); 
complex<T> t9 = square(spa23); 
complex<T> t10 = s15*spa25; 
complex<T> d1 = spb12*square(s15 - s34)*T(2); d1 = T(1)/d1;
complex<T> d2 = (s15 - s34)*spb12; d2 = T(1)/d2;
complex<T> d3 = spb12*spb15*spb23*spb34*T(2); d3 = T(1)/d3;
complex<T> d4 = (s12 + s15 - s34)*spb34; d4 = T(1)/d4;
complex<T> d5 = (s12 + s15 - s34)*spb12*spb34; d5 = T(1)/d5;
complex<T> d6 = spb12*spb23*spb34; d6 = T(1)/d6;
complex<T> d7 = (s12 + s15 - s34)*spb12; d7 = T(1)/d7;
complex<T> d8 = (s12 + s15 - s34)*spb34*T(2); d8 = T(1)/d8;
complex<T> t3 = d2*spa23*spb45 + d1*spb25*spb34*t9; 
complex<T> t4 = -(d2*spa23*spb45) - d1*spb25*spb34*t9 + d3*spb13*t6*T(3); 
complex<T> t16 = t10*t6; 
complex<T> t5 = d5*t16 + d6*spa15*spb13*t6; 
complex<T> co1 = -(d4*spa12*spa25*t6); 
complex<T> co2 = d7*spa25*spa34*t6; 
complex<T> co3 = -(d8*spa12*t16); 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + t5*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)) + co3*(*CI_users[5]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmmQpQmqp_LLT_wCI::\
C2q2Q1g_qmmQpQmqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmmQpQmqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, Qp, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmmQpQmqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:00 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> t3 = spb12*spb34; 
complex<T> t8 = square(spb35); 
complex<T> t9 = -(s25*spa25); 
complex<T> t13 = square(spa24); 
complex<T> d3 = spb12*square(s15 - s34)*T(2); d3 = T(1)/d3;
complex<T> d7 = spb34*square(s12 + s15 - s34); d7 = T(1)/d7;
complex<T> d10 = spb12*square(s12 + s15 - s34); d10 = T(1)/d10;
complex<T> d11 = spb34*square(s12 + s15 - s34)*T(2); d11 = T(1)/d11;
complex<T> t21 = spa25*t8; 
complex<T> t24 = t8*t9; 
complex<T> d1 = (s12 - s34)*t3; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*(s12 + s15 - s34)*t3; d2 = T(1)/d2;
complex<T> d4 = (s15 - s34)*t3; d4 = T(1)/d4;
complex<T> d5 = (s15 - s34)*(s12 + s15 - s34)*t3; d5 = T(1)/d5;
complex<T> d6 = spb15*spb23*t3*T(2); d6 = T(1)/d6;
complex<T> d8 = t3*square(s12 + s15 - s34); d8 = T(1)/d8;
complex<T> d9 = spb23*t3; d9 = T(1)/d9;
complex<T> t5 = d4*spa12*spb13*spb35 + d3*spb25*spb34*t13 - d1*t21 + d2*t24 + d5*t24; 
complex<T> t20 = d8*s15*t24 + d9*spa15*spb13*t8; 
complex<T> t23 = s25*t21; 
complex<T> t4 = d1*t21 + d2*t23; 
complex<T> t6 = -(d4*spa12*spb13*spb35) - d3*spb25*spb34*t13 + d5*t23 + d6*spb13*t8*T(3); 
complex<T> co1 = d7*spa12*t23; 
complex<T> co2 = d10*spa34*t24; 
complex<T> co3 = d11*s15*spa12*t23; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + t20*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmpQmQpqp_LLT_wCI::\
C2q2Q1g_qmpQmQpqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmpQmQpqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, Qm, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpQmQpqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:03 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = S(2,5);
complex<T> t4 = spa12*spa34; 
complex<T> t10 = square(spa15); 
complex<T> t11 = square(spa23); 
complex<T> t12 = cube(spb25); 
complex<T> t13 = s12 + s15 - s34; 
complex<T> t14 = square(spa13); 
complex<T> t26 = spa13*spb45; 
complex<T> d4 = spa12*square(s12 - s34); d4 = T(1)/d4;
complex<T> t1 = square(t13); 
complex<T> t9 = t4*T(2); 
complex<T> t24 = t10*t11; 
complex<T> d8 = spa34*cube(t13); d8 = T(1)/d8;
complex<T> d9 = spa23*t4; d9 = T(1)/d9;
complex<T> d10 = t4*cube(t13); d10 = T(1)/d10;
complex<T> d11 = spa12*cube(t13); d11 = T(1)/d11;
complex<T> d12 = spa34*cube(t13)*T(2); d12 = T(1)/d12;
complex<T> t18 = -t24; 
complex<T> t29 = t12*t24; 
complex<T> d1 = t9*square(s12 - s34); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*t1*t4; d2 = T(1)/d2;
complex<T> d3 = t13*t9*square(s12 - s34); d3 = T(1)/d3;
complex<T> d5 = spa15*spa23*t9; d5 = T(1)/d5;
complex<T> d6 = (s15 - s34)*t1*t4; d6 = T(1)/d6;
complex<T> d7 = t13*t9*square(s15 - s34); d7 = T(1)/d7;
complex<T> t6 = (d6 + d7)*t29 - d5*cube(spa13)*T(3); 
complex<T> t8 = d1*s25*spb25*t14 - d4*spa15*spb25*t26 + (d2 + d3)*t29; 
complex<T> t17 = -(d1*s25); 
complex<T> t28 = t12*t18; 
complex<T> t7 = d10*s15*t28 - d9*spb15*cube(spa13); 
complex<T> t23 = spb25*(t14*t17 + d4*spa15*t26) + (d2 + d3 + d6 + d7)*t28; 
complex<T> co1 = d8*spb12*t29; 
complex<T> co2 = d11*spb34*t28; 
complex<T> co3 = d12*s15*spb12*t29; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t8*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t23*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmpQpQmqp_LLT_wCI::\
C2q2Q1g_qmpQpQmqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmpQpQmqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, Qp, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpQpQmqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:06 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = S(2,5);
complex<T> t4 = spa12*spa34; 
complex<T> t10 = square(spa15); 
complex<T> t11 = square(spa24); 
complex<T> t12 = cube(spb25); 
complex<T> t13 = s12 + s15 - s34; 
complex<T> t14 = square(spa14); 
complex<T> t25 = spb25*spb35; 
complex<T> d4 = spa12*square(s12 - s34); d4 = T(1)/d4;
complex<T> t1 = square(t13); 
complex<T> t9 = t4*T(2); 
complex<T> t24 = t10*t11; 
complex<T> d8 = spa34*cube(t13); d8 = T(1)/d8;
complex<T> d9 = spa23*t4; d9 = T(1)/d9;
complex<T> d10 = t4*cube(t13); d10 = T(1)/d10;
complex<T> d11 = spa12*cube(t13); d11 = T(1)/d11;
complex<T> d12 = spa34*cube(t13)*T(2); d12 = T(1)/d12;
complex<T> t18 = -t24; 
complex<T> t30 = t12*t24; 
complex<T> d1 = t9*square(s12 - s34); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*t1*t4; d2 = T(1)/d2;
complex<T> d3 = t13*t9*square(s12 - s34); d3 = T(1)/d3;
complex<T> d5 = spa15*spa23*t9; d5 = T(1)/d5;
complex<T> d6 = (s15 - s34)*t1*t4; d6 = T(1)/d6;
complex<T> d7 = t13*t9*square(s15 - s34); d7 = T(1)/d7;
complex<T> t6 = (d6 + d7)*t30 - d5*spa13*t14*T(3); 
complex<T> t19 = d1*s25; 
complex<T> t28 = t12*t18; 
complex<T> t7 = -(d9*spa13*spb15*t14) + d10*s15*t28; 
complex<T> t8 = -(spb25*t14*t19) - d4*spa14*spa15*t25 + (d2 + d3 + d6 + d7)*t28; 
complex<T> t23 = spb25*t14*t19 + d4*spa14*spa15*t25 + (d2 + d3)*t30; 
complex<T> co1 = d8*spb12*t30; 
complex<T> co2 = d11*spb34*t28; 
complex<T> co3 = d12*s15*spb12*t30; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t23*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t8*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQmmQpqp_LLT_wCI::\
C2q2Q1g_qmQmmQpqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQmmQpqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, m, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQmmQpqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:07 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> t1 = square(spb45); 
complex<T> d1 = spb15*spb23*spb34*T(2); d1 = T(1)/d1;
complex<T> d2 = spb23*spb34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spa15*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQpmQmqp_LLT_wCI::\
C2q2Q1g_qmQpmQmqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQpmQmqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, m, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQpmQmqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:08 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> t1 = square(spb25); 
complex<T> d1 = spb15*spb23*spb34*T(2); d1 = T(1)/d1;
complex<T> d2 = spb23*spb34; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spa15*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQmpQpqp_LLT_wCI::\
C2q2Q1g_qmQmpQpqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQmpQpqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, p, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQmpQpqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:09 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> t1 = square(spa12); 
complex<T> d1 = spa15*spa23*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = spa23*spa34; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d2*spb15*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQppQmqp_LLT_wCI::\
C2q2Q1g_qmQppQmqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQppQmqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, p, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQppQmqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:11 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa15*spa23*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = spa23*spa34; d2 = T(1)/d2;
complex<T> co1 = -(d1*t1*T(3)); 
complex<T> co2 = -(d2*spb15*t1); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*(*CI_users[0]->get_value(mc,ind,mu)) + co2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQmQpqpm_LLT_wCI::\
C2q2Q1g_qmQmQpqpm_LLT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQmQpqpm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp, m}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQmQpqpm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:13 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb13 = SPB(1,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = square(s23 - s45); 
complex<T> t4 = square(spb34); 
complex<T> t16 = spa15*spb13; 
complex<T> t17 = s23*spb34; 
complex<T> d4 = spb15*spb23*spb45*T(2); d4 = T(1)/d4;
complex<T> d5 = spb15*spb45; d5 = T(1)/d5;
complex<T> d6 = spb23*T(2); d6 = T(1)/d6;
complex<T> t12 = spa23*t4; 
complex<T> d1 = spb15*t1; d1 = T(1)/d1;
complex<T> d2 = spb15*t1*T(2); d2 = T(1)/d2;
complex<T> d3 = spb15*spb23*t1*T(2); d3 = T(1)/d3;
complex<T> t8 = d3*spa45; 
complex<T> t15 = d2*spa25; 
complex<T> t2 = -(spb45*t15*t16) + d1*(spa45*t12 + spa25*t17)*T(2) - s45*spb34*t15*T(3) + d4*t4*T(3) + s45*t4*t8*T(3); 
complex<T> t3 = spb45*t15*t16 - d1*(spa45*t12 + spa25*t17)*T(2) + s45*spb34*t15*T(3) - s45*t4*t8*T(3); 
complex<T> co1 = d5*t12; 
complex<T> co2 = d6*spa15*spa45*t4; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQmQpqpp_LLT_wCI::\
C2q2Q1g_qmQmQpqpp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQmQpqpp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, qp, p}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQmQpqpp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:14 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa24 = SPA(2,4);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> t1 = square(s15 - s23); 
complex<T> t4 = square(spa12); 
complex<T> t16 = spa15*spa24; 
complex<T> t17 = s23*spa12; 
complex<T> d4 = spa15*spa23*spa45*T(2); d4 = T(1)/d4;
complex<T> d5 = spa15*spa45; d5 = T(1)/d5;
complex<T> d6 = spa23*T(2); d6 = T(1)/d6;
complex<T> t12 = spb15*t4; 
complex<T> d1 = spa23*spa45*t1*T(2); d1 = T(1)/d1;
complex<T> d2 = spa45*t1; d2 = T(1)/d2;
complex<T> d3 = spa45*t1*T(2); d3 = T(1)/d3;
complex<T> t8 = d1*spb15; 
complex<T> t15 = d3*spb35; 
complex<T> t2 = -(spb45*t15*t16) + d2*(spb23*t12 + spb35*t17)*T(2) - s15*spa12*t15*T(3) + s15*t4*t8*T(3); 
complex<T> t3 = spb45*t15*t16 - d2*(spb23*t12 + spb35*t17)*T(2) + s15*spa12*t15*T(3) - d4*t4*T(3) - s15*t4*t8*T(3); 
complex<T> co1 = -(d5*spb23*t4); 
complex<T> co2 = -(d6*spb45*t12); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQpQmqpm_LLT_wCI::\
C2q2Q1g_qmQpQmqpm_LLT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQpQmqpm_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp, m}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQpQmqpm LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:16 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = square(s23 - s45); 
complex<T> t4 = square(spb24); 
complex<T> t16 = spa15*spb12; 
complex<T> t17 = s23*spb24; 
complex<T> d4 = spb15*spb23*spb45*T(2); d4 = T(1)/d4;
complex<T> d5 = spb15*spb45; d5 = T(1)/d5;
complex<T> d6 = spb23*T(2); d6 = T(1)/d6;
complex<T> t12 = spa23*t4; 
complex<T> d1 = spb15*t1; d1 = T(1)/d1;
complex<T> d2 = spb15*t1*T(2); d2 = T(1)/d2;
complex<T> d3 = spb15*spb23*t1*T(2); d3 = T(1)/d3;
complex<T> t8 = d3*spa45; 
complex<T> t15 = d2*spa35; 
complex<T> t2 = -(spb45*t15*t16) - d1*spa45*t12*T(2) + d1*spa35*t17*T(2) - s45*spb24*t15*T(3) - s45*t4*t8*T(3); 
complex<T> t3 = spb45*t15*t16 + d1*spa45*t12*T(2) - d1*spa35*t17*T(2) + s45*spb24*t15*T(3) + d4*t4*T(3) + s45*t4*t8*T(3); 
complex<T> co1 = d5*t12; 
complex<T> co2 = d6*spa15*spa45*t4; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQpQmqpp_LLT_wCI::\
C2q2Q1g_qmQpQmqpp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQpQmqpp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, qp, p}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQpQmqpp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:18 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> t1 = square(s15 - s23); 
complex<T> t4 = square(spa13); 
complex<T> t16 = spa15*spa34; 
complex<T> t17 = s23*spa13; 
complex<T> d4 = spa15*spa23*spa45*T(2); d4 = T(1)/d4;
complex<T> d5 = spa15*spa45; d5 = T(1)/d5;
complex<T> d6 = spa23*T(2); d6 = T(1)/d6;
complex<T> t12 = spb15*t4; 
complex<T> d1 = spa23*spa45*t1*T(2); d1 = T(1)/d1;
complex<T> d2 = spa45*t1; d2 = T(1)/d2;
complex<T> d3 = spa45*t1*T(2); d3 = T(1)/d3;
complex<T> t8 = d1*spb15; 
complex<T> t15 = d3*spb25; 
complex<T> t2 = -(spb45*t15*t16) - d2*spb23*t12*T(2) + d2*spb25*t17*T(2) - s15*spa13*t15*T(3) - d4*t4*T(3) - s15*t4*t8*T(3); 
complex<T> t3 = spb45*t15*t16 + d2*spb23*t12*T(2) - d2*spb25*t17*T(2) + s15*spa13*t15*T(3) + s15*t4*t8*T(3); 
complex<T> co1 = -(d5*spb23*t4); 
complex<T> co2 = -(d6*spb45*t12); 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQmQpmqp_LLT_wCI::\
C2q2Q1g_qmQmQpmqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQmQpmqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, m, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQmQpmqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:21 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = S(1,4);
complex<T> t4 = spb23*spb45; 
complex<T> t10 = cube(spa14); 
complex<T> t11 = square(spb15); 
complex<T> t12 = square(spb34); 
complex<T> t13 = s15 - s23 + s45; 
complex<T> t14 = square(spb35); 
complex<T> t31 = spa14*spb15; 
complex<T> d6 = spb45*square(s23 - s45); d6 = T(1)/d6;
complex<T> t1 = square(t13); 
complex<T> t9 = t4*T(2); 
complex<T> t21 = -(d6*spa12); 
complex<T> t25 = t10*t11; 
complex<T> d8 = t4*cube(t13); d8 = T(1)/d8;
complex<T> d9 = spb34*t4; d9 = T(1)/d9;
complex<T> d10 = spb45*cube(t13); d10 = T(1)/d10;
complex<T> d11 = spb23*cube(t13); d11 = T(1)/d11;
complex<T> d12 = spb23*cube(t13)*T(2); d12 = T(1)/d12;
complex<T> t18 = -t25; 
complex<T> t30 = t12*t25; 
complex<T> d1 = (s15 - s23)*t1*t4; d1 = T(1)/d1;
complex<T> d2 = t13*t9*square(s15 - s23); d2 = T(1)/d2;
complex<T> d3 = spb15*spb34*t9; d3 = T(1)/d3;
complex<T> d4 = (s23 - s45)*t1*t4; d4 = T(1)/d4;
complex<T> d5 = t13*t9*square(s23 - s45); d5 = T(1)/d5;
complex<T> d7 = t9*square(s23 - s45); d7 = T(1)/d7;
complex<T> t8 = d8*s15*t30 + d9*spa15*cube(spb35); 
complex<T> t20 = d7*s14; 
complex<T> t29 = t12*t18; 
complex<T> t6 = -(spa14*t14*t20) + d5*t29 + d4*t30 + d6*spa12*spb35*t31; 
complex<T> t7 = (d1 + d2)*t29 + d3*cube(spb35)*T(3); 
complex<T> t24 = spa14*t14*t20 + d4*t29 + d1*t30 + d2*t30 + d5*t30 + spb35*t21*t31; 
complex<T> co1 = d10*spa23*t30; 
complex<T> co2 = d11*spa45*t29; 
complex<T> co3 = d12*s15*spa45*t29; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t7*(*CI_users[0]->get_value(mc,ind,mu)) + t24*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t8*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQmQppqp_LLT_wCI::\
C2q2Q1g_qmQmQppqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQmQppqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qm, Qp, p, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQmQppqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:23 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t6 = square(spa12); 
complex<T> t9 = square(spb34); 
complex<T> t14 = s15*spb14; 
complex<T> d1 = spa15*spa23*spa34*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*spa45; d2 = T(1)/d2;
complex<T> d3 = spa45*square(s15 - s23)*T(2); d3 = T(1)/d3;
complex<T> d4 = (s15 - s23 + s45)*spa23*spa45; d4 = T(1)/d4;
complex<T> d5 = spa23*spa34*spa45; d5 = T(1)/d5;
complex<T> d6 = (s15 - s23 + s45)*spa45; d6 = T(1)/d6;
complex<T> d7 = (s15 - s23 + s45)*spa23; d7 = T(1)/d7;
complex<T> d8 = (s15 - s23 + s45)*spa23*T(2); d8 = T(1)/d8;
complex<T> t3 = -((d5*spa35*spb15 + d4*t14)*t6); 
complex<T> t4 = -(d2*spa12*spb34) - d3*spa14*spa23*t9; 
complex<T> t5 = d2*spa12*spb34 + d3*spa14*spa23*t9 - d1*spa35*t6*T(3); 
complex<T> co1 = -(d6*spb14*spb23*t6); 
complex<T> co2 = d7*spb14*spb45*t6; 
complex<T> co3 = d8*spb45*t14*t6; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)) + co1*(*CI_users[3]->get_value(mc,ind,mu)) + co2*(*CI_users[4]->get_value(mc,ind,mu)) + co3*(*CI_users[5]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQpQmmqp_LLT_wCI::\
C2q2Q1g_qmQpQmmqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQpQmmqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, m, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQpQmmqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:26 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = S(1,4);
complex<T> t4 = spb23*spb45; 
complex<T> t10 = cube(spa14); 
complex<T> t11 = square(spb15); 
complex<T> t12 = square(spb24); 
complex<T> t13 = s15 - s23 + s45; 
complex<T> t14 = square(spb25); 
complex<T> t31 = spa14*spb15; 
complex<T> d6 = spb45*square(s23 - s45); d6 = T(1)/d6;
complex<T> t1 = square(t13); 
complex<T> t9 = t4*T(2); 
complex<T> t25 = t10*t11; 
complex<T> d8 = t4*cube(t13); d8 = T(1)/d8;
complex<T> d9 = spb34*t4; d9 = T(1)/d9;
complex<T> d10 = spb45*cube(t13); d10 = T(1)/d10;
complex<T> d11 = spb23*cube(t13); d11 = T(1)/d11;
complex<T> d12 = spb23*cube(t13)*T(2); d12 = T(1)/d12;
complex<T> t18 = -t25; 
complex<T> t30 = t12*t25; 
complex<T> d1 = (s15 - s23)*t1*t4; d1 = T(1)/d1;
complex<T> d2 = t13*t9*square(s15 - s23); d2 = T(1)/d2;
complex<T> d3 = spb15*spb34*t9; d3 = T(1)/d3;
complex<T> d4 = (s23 - s45)*t1*t4; d4 = T(1)/d4;
complex<T> d5 = t13*t9*square(s23 - s45); d5 = T(1)/d5;
complex<T> d7 = t9*square(s23 - s45); d7 = T(1)/d7;
complex<T> t8 = d9*spa15*spb35*t14 + d8*s15*t30; 
complex<T> t19 = d7*s14; 
complex<T> t29 = t12*t18; 
complex<T> t6 = -(spa14*t14*t19) + d5*t29 + d4*t30 - d6*spa13*spb25*t31; 
complex<T> t7 = (d1 + d2)*t29 + d3*spb35*t14*T(3); 
complex<T> t24 = spa14*t14*t19 + d4*t29 + d1*t30 + d2*t30 + d5*t30 + d6*spa13*spb25*t31; 
complex<T> co1 = d10*spa23*t30; 
complex<T> co2 = d11*spa45*t29; 
complex<T> co3 = d12*s15*spa45*t29; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t7*(*CI_users[0]->get_value(mc,ind,mu)) + t24*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t8*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmQpQmpqp_LLT_wCI::\
C2q2Q1g_qmQpQmpqp_LLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c15, c234));
CI_users.push_back(new Cached_Bubble_Integral_User(c23, c145));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmQpQmpqp_LLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, Qp, Qm, p, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmQpQmpqp LLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:28 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa15 = SPA(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = -(spa14*spb14);
complex<T> s15 = -(spa15*spb15);
complex<T> t3 = spa23*spa45; 
complex<T> t8 = square(spa13); 
complex<T> t9 = -(s14*spb14); 
complex<T> t13 = square(spb24); 
complex<T> d3 = spa45*square(s15 - s23)*T(2); d3 = T(1)/d3;
complex<T> d9 = spa45*square(s15 - s23 + s45); d9 = T(1)/d9;
complex<T> d10 = spa23*square(s15 - s23 + s45); d10 = T(1)/d10;
complex<T> d11 = spa23*square(s15 - s23 + s45)*T(2); d11 = T(1)/d11;
complex<T> t21 = spb14*t8; 
complex<T> t24 = t8*t9; 
complex<T> d1 = spa15*spa34*t3*T(2); d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*(s15 - s23 + s45)*t3; d2 = T(1)/d2;
complex<T> d4 = (s15 - s23)*t3; d4 = T(1)/d4;
complex<T> d5 = (s23 - s45)*t3; d5 = T(1)/d5;
complex<T> d6 = (s23 - s45)*(s15 - s23 + s45)*t3; d6 = T(1)/d6;
complex<T> d7 = t3*square(s15 - s23 + s45); d7 = T(1)/d7;
complex<T> d8 = spa34*t3; d8 = T(1)/d8;
complex<T> t6 = d4*spa13*spa35*spb45 + d3*spa14*spa23*t13 + d2*t24 - d1*spa35*t8*T(3); 
complex<T> t23 = s14*t21; 
complex<T> t4 = d5*t21 + d6*t23; 
complex<T> t5 = -(d4*spa13*spa35*spb45) - d3*spa14*spa23*t13 - d5*t21 + d2*t23 + d6*t24; 
complex<T> t20 = d7*s15*t23 - d8*spa35*spb15*t8; 
complex<T> co1 = d9*spb23*t23; 
complex<T> co2 = d10*spb45*t24; 
complex<T> co3 = d11*s15*spb45*t24; 
complex<T> co4 = Complex(0,1); 
SeriesC<T> result = co4*(t6*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t4*(*CI_users[2]->get_value(mc,ind,mu)) + t20*(*CI_users[3]->get_value(mc,ind,mu)) + co1*(*CI_users[4]->get_value(mc,ind,mu)) + co2*(*CI_users[5]->get_value(mc,ind,mu)) + co3*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpmQmQp_nfLLT_wCI::\
C2q2Q1g_qmqpmQmQp_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpmQmQp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpmQmQp nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:30 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> t5 = spb23*spb45; 
complex<T> t6 = square(spa13); 
complex<T> t11 = spb12*spb35; 
complex<T> d1 = cube(s12 - s45)*T(3); d1 = T(1)/d1;
complex<T> t16 = t11*t6; 
complex<T> d2 = (s12 - s45)*t5; d2 = T(1)/d2;
complex<T> d3 = spb12*spb34*t5*T(3); d3 = T(1)/d3;
complex<T> d4 = t5*square(s12 - s45)*T(3); d4 = T(1)/d4;
complex<T> t3 = -(d4*spb25*t16) + d2*spa13*square(spb25) - d1*spa34*t16*T(2) - d3*spb24*square(spb25)*T(2); 
complex<T> t4 = d4*spb25*t16 - d2*spa13*square(spb25) + d1*spa34*t16*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqppQmQp_nfLLT_wCI::\
C2q2Q1g_qmqppQmQp_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqppQmQp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqppQmQp nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:32 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> t3 = spa12*spa34; 
complex<T> t6 = square(spb35); 
complex<T> t11 = spa13*spa45; 
complex<T> d3 = cube(s12 - s45)*T(3); d3 = T(1)/d3;
complex<T> t16 = t11*t6; 
complex<T> d1 = (s12 - s45)*t3; d1 = T(1)/d1;
complex<T> d2 = t3*square(s12 - s45)*T(3); d2 = T(1)/d2;
complex<T> d4 = spa23*spa45*t3*T(3); d4 = T(1)/d4;
complex<T> t4 = d2*spa14*t16 + d1*spb35*square(spa14) - d3*spb23*t16*T(2) + d4*spa24*square(spa14)*T(2); 
complex<T> t5 = -(d2*spa14*t16) - d1*spb35*square(spa14) + d3*spb23*t16*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpmQpQm_nfLLT_wCI::\
C2q2Q1g_qmqpmQpQm_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpmQpQm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, m, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpmQpQm nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:57:34 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> t10 = square(spa35); 
complex<T> t14 = cube(spb24); 
complex<T> t16 = square(spa13); 
complex<T> d1 = (s12 - s45)*spb12*spb34*T(3); d1 = T(1)/d1;
complex<T> d2 = cube(s12 - s45)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s12 - s45)*spb23*spb45*T(3); d3 = T(1)/d3;
complex<T> d4 = spb12*spb23*spb34*spb45*T(3); d4 = T(1)/d4;
complex<T> t17 = d2*spa35; 
complex<T> t22 = spb23*t10; 
complex<T> t5 = -(d4*t14) - spb12*spb34*t16*t17 - d2*spa13*spb45*t22 - d3*spa13*square(spb24) - d1*spa35*square(spb24); 
complex<T> t6 = -(d4*t14) + spb12*spb34*t16*t17 + d2*spa13*spb45*t22 + d3*spa13*square(spb24) + d1*spa35*square(spb24); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t6*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqppQpQm_nfLLT_wCI::\
C2q2Q1g_qmqppQpQm_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqppQpQm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, p, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqppQpQm nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:06 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> t7 = square(spa15); 
complex<T> t8 = s12 - s45; 
complex<T> t11 = square(spb23); 
complex<T> t12 = square(spb34); 
complex<T> t18 = spa12*spa35; 
complex<T> t26 = spb23*spb34; 
complex<T> t28 = spa24*spa35; 
complex<T> d1 = spa12*spa23*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> t13 = d1*spa24; 
complex<T> t16 = -t26; 
complex<T> t34 = spb15*t18; 
complex<T> d2 = spa23*spa45*t8*T(2); d2 = T(1)/d2;
complex<T> d3 = spa23*spa45*t8*T(6); d3 = T(1)/d3;
complex<T> d4 = spa34*square(t8)*T(6); d4 = T(1)/d4;
complex<T> d5 = spa12*spa34*t8*T(6); d5 = T(1)/d5;
complex<T> d6 = spa23*square(t8)*T(6); d6 = T(1)/d6;
complex<T> d7 = spa23*spa34*cube(t8)*T(3); d7 = T(1)/d7;
complex<T> d8 = spa12*spa34*t8*T(2); d8 = T(1)/d8;
complex<T> t27 = d7*spa13; 
complex<T> t32 = t26*t27; 
complex<T> t4 = -(d3*spa15*spa25*spb23) - d5*spa14*spa15*spb34 + d6*spa13*spa45*t12 + d4*spa13*spa45*t16 - d4*t11*t18 + d6*t18*t26 + s12*t16*t27*t28 + s45*t16*t27*t28 - d2*spb13*t7 - d8*spb35*t7 + t13*t7 - spa45*t32*t34*T(2); 
complex<T> t5 = d3*spa15*spa25*spb23 + d5*spa14*spa15*spb34 - d6*spa13*spa45*t12 + d4*t11*t18 + d6*t16*t18 + d4*spa13*spa45*t26 + s12*t28*t32 + s45*t28*t32 + d2*spb13*t7 + d8*spb35*t7 + t13*t7 + spa45*t32*t34*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmmQp_nfLLT_wCI::\
C2q2Q1g_qmqpQmmQp_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmmQp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, m, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmmQp nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:07 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb25); 
complex<T> d1 = spb12*spb34*spb45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmpQp_nfLLT_wCI::\
C2q2Q1g_qmqpQmpQp_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmpQp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, p, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmpQp nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:08 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa13); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpmQm_nfLLT_wCI::\
C2q2Q1g_qmqpQpmQm_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpmQm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, m, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpmQm nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:09 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb23); 
complex<T> d1 = spb12*spb34*spb45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQppQm_nfLLT_wCI::\
C2q2Q1g_qmqpQppQm_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQppQm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, p, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQppQm nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:11 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa15); 
complex<T> d1 = spa12*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmQpm_nfLLT_wCI::\
C2q2Q1g_qmqpQmQpm_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmQpm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, m}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmQpm nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:12 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t5 = spb12*spb45; 
complex<T> t6 = square(spa35); 
complex<T> t11 = spb25*spb34; 
complex<T> d1 = cube(s12 - s34)*T(3); d1 = T(1)/d1;
complex<T> t16 = t11*t6; 
complex<T> d2 = (s12 - s34)*t5; d2 = T(1)/d2;
complex<T> d3 = t5*square(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> d4 = spb15*spb34*t5*T(3); d4 = T(1)/d4;
complex<T> t3 = -(d3*spb24*t16) - d2*spa35*square(spb24) + d1*spa15*t16*T(2) - d4*spb14*square(spb24)*T(2); 
complex<T> t4 = d3*spb24*t16 + d2*spa35*square(spb24) - d1*spa15*t16*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQmQpp_nfLLT_wCI::\
C2q2Q1g_qmqpQmQpp_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQmQpp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qm, Qp, p}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQmQpp nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:14 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t3 = spa15*spa34; 
complex<T> t6 = square(spb25); 
complex<T> t11 = spa12*spa35; 
complex<T> d4 = cube(s12 - s34)*T(3); d4 = T(1)/d4;
complex<T> t16 = t11*t6; 
complex<T> d1 = spa12*spa45*t3*T(3); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*t3; d2 = T(1)/d2;
complex<T> d3 = t3*square(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> t4 = d2*spb25*square(spa13) - t16*(d3*spa13 + d4*spb45*T(2)); 
complex<T> t5 = d3*spa13*t16 - d2*spb25*square(spa13) + d4*spb45*t16*T(2) + d1*spa14*square(spa13)*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpQmm_nfLLT_wCI::\
C2q2Q1g_qmqpQpQmm_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpQmm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, m}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpQmm nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:47 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb35 = SPB(3,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t7 = square(spb23); 
complex<T> t8 = s12 - s34; 
complex<T> t11 = square(spa15); 
complex<T> t12 = square(spa45); 
complex<T> t18 = spb12*spb35; 
complex<T> t26 = spa15*spa45; 
complex<T> t28 = spb14*spb35; 
complex<T> d6 = spb12*spb15*spb34*spb45*T(3); d6 = T(1)/d6;
complex<T> t15 = -(d6*spb14); 
complex<T> t16 = -t26; 
complex<T> t34 = spa23*t18; 
complex<T> d1 = spb15*spb34*t8*T(6); d1 = T(1)/d1;
complex<T> d2 = spb15*spb34*t8*T(2); d2 = T(1)/d2;
complex<T> d3 = spb15*square(t8)*T(6); d3 = T(1)/d3;
complex<T> d4 = spb12*spb45*t8*T(2); d4 = T(1)/d4;
complex<T> d5 = spb12*spb45*t8*T(6); d5 = T(1)/d5;
complex<T> d7 = spb45*square(t8)*T(6); d7 = T(1)/d7;
complex<T> d8 = spb15*spb45*cube(t8)*T(3); d8 = T(1)/d8;
complex<T> t27 = d8*spb25; 
complex<T> t32 = t26*t27; 
complex<T> t4 = -(d1*spa15*spb13*spb23) - d5*spa45*spb23*spb24 + d3*spb25*spb34*t12 + d7*spb25*spb34*t16 - d7*t11*t18 + d3*t18*t26 + s12*t16*t27*t28 + s34*t16*t27*t28 - d2*spa25*t7 - d4*spa35*t7 + t15*t7 - spb34*t32*t34*T(2); 
complex<T> t5 = d1*spa15*spb13*spb23 + d5*spa45*spb23*spb24 - d3*spb25*spb34*t12 + d7*t11*t18 + d3*t16*t18 + d7*spb25*spb34*t26 + s12*t28*t32 + s34*t28*t32 + d2*spa25*t7 + d4*spa35*t7 + t15*t7 + spb34*t32*t34*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmqpQpQmp_nfLLT_wCI::\
C2q2Q1g_qmqpQpQmp_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c12, c345));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmqpQpQmp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, qp, Qp, Qm, p}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmqpQpQmp nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:48 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t10 = square(spb35); 
complex<T> t13 = cube(spa14); 
complex<T> t16 = square(spb25); 
complex<T> d1 = spa12*spa15*spa34*spa45*T(3); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa15*spa34*T(3); d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*spa12*spa45*T(3); d3 = T(1)/d3;
complex<T> d4 = cube(s12 - s34)*T(3); d4 = T(1)/d4;
complex<T> t17 = d4*spb35; 
complex<T> t22 = spa15*t10; 
complex<T> t5 = d1*t13 - spa12*spa45*t16*t17 - d4*spa34*spb25*t22 - d2*spb25*square(spa14) - d3*spb35*square(spa14); 
complex<T> t6 = d1*t13 + spa12*spa45*t16*t17 + d4*spa34*spb25*t22 + d2*spb25*square(spa14) + d3*spb35*square(spa14); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t5*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C2q2Q1g_qmmqpQmQp_nfLLT_wCI::\
C2q2Q1g_qmmqpQmQp_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmmqpQmQp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, qp, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmmqpQmQp nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:49 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb35); 
complex<T> d1 = spb12*spb23*spb45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1g_qmpqpQmQp_nfLLT_wCI::\
C2q2Q1g_qmpqpQmQp_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmpqpQmQp_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpqpQmQp nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:51 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa14); 
complex<T> d1 = spa12*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1g_qmmqpQpQm_nfLLT_wCI::\
C2q2Q1g_qmmqpQpQm_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmmqpQpQm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, m, qp, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmmqpQpQm nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:52 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb34); 
complex<T> d1 = spb12*spb23*spb45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,-2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2Q1g_qmpqpQpQm_nfLLT_wCI::\
C2q2Q1g_qmpqpQpQm_nfLLT_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C2q2Q1g_qmpqpQpQm_nfLLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{qm, p, qp, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("C2q2Q1g :  qmpqpQpQm nfLLT");
#endif
 
//#define TimeStamp "Wed 22 Dec 2010 16:58:53 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> t1 = square(spa15); 
complex<T> d1 = spa12*spa23*spa45*T(3); d1 = T(1)/d1;
complex<T> co1 = Complex(0,2)*d1*t1; 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qmqpmQmQp_LLT C2q2Q1g_7357_LLT
#define _C_qmqppQmQp_LLT C2q2Q1g_7465_LLT
#define _C_qmqpmQpQm_LLT C2q2Q1g_6277_LLT
#define _C_qmqppQpQm_LLT C2q2Q1g_6385_LLT
#define _C_qmqpQmmQp_LLT C2q2Q1g_6637_LLT
#define _C_qmqpQmpQp_LLT C2q2Q1g_7285_LLT
#define _C_qmqpQpmQm_LLT C2q2Q1g_5377_LLT
#define _C_qmqpQppQm_LLT C2q2Q1g_6025_LLT
#define _C_qmqpQmQpm_LLT C2q2Q1g_1237_LLT
#define _C_qmqpQmQpp_LLT C2q2Q1g_5125_LLT
#define _C_qmqpQpQmm_LLT C2q2Q1g_1057_LLT
#define _C_qmqpQpQmp_LLT C2q2Q1g_4945_LLT
#define _C_qmmqpQmQp_LLT C2q2Q1g_7417_LLT
#define _C_qmpqpQmQp_LLT C2q2Q1g_7435_LLT
#define _C_qmmqpQpQm_LLT C2q2Q1g_6337_LLT
#define _C_qmpqpQpQm_LLT C2q2Q1g_6355_LLT
#define _C_qmqpmQmQp_LRT C2q2Q1g_7357_LRT
#define _C_qmqppQmQp_LRT C2q2Q1g_7465_LRT
#define _C_qmqpmQpQm_LRT C2q2Q1g_6277_LRT
#define _C_qmqppQpQm_LRT C2q2Q1g_6385_LRT
#define _C_qmqpQmmQp_LRT C2q2Q1g_6637_LRT
#define _C_qmqpQmpQp_LRT C2q2Q1g_7285_LRT
#define _C_qmqpQpmQm_LRT C2q2Q1g_5377_LRT
#define _C_qmqpQppQm_LRT C2q2Q1g_6025_LRT
#define _C_qmqpQmQpm_LRT C2q2Q1g_1237_LRT
#define _C_qmqpQmQpp_LRT C2q2Q1g_5125_LRT
#define _C_qmqpQpQmm_LRT C2q2Q1g_1057_LRT
#define _C_qmqpQpQmp_LRT C2q2Q1g_4945_LRT
#define _C_qmmqpQmQp_LRT C2q2Q1g_7417_LRT
#define _C_qmpqpQmQp_LRT C2q2Q1g_7435_LRT
#define _C_qmmqpQpQm_LRT C2q2Q1g_6337_LRT
#define _C_qmpqpQpQm_LRT C2q2Q1g_6355_LRT
#define _C_qmmQmQpqp_LLT C2q2Q1g_3817_LLT
#define _C_qmmQpQmqp_LLT C2q2Q1g_3637_LLT
#define _C_qmpQmQpqp_LLT C2q2Q1g_3835_LLT
#define _C_qmpQpQmqp_LLT C2q2Q1g_3655_LLT
#define _C_qmQmmQpqp_LLT C2q2Q1g_3697_LLT
#define _C_qmQpmQmqp_LLT C2q2Q1g_3487_LLT
#define _C_qmQmpQpqp_LLT C2q2Q1g_3805_LLT
#define _C_qmQppQmqp_LLT C2q2Q1g_3595_LLT
#define _C_qmQmQpqpm_LLT C2q2Q1g_637_LLT
#define _C_qmQmQpqpp_LLT C2q2Q1g_4525_LLT
#define _C_qmQpQmqpm_LLT C2q2Q1g_607_LLT
#define _C_qmQpQmqpp_LLT C2q2Q1g_4495_LLT
#define _C_qmQmQpmqp_LLT C2q2Q1g_2797_LLT
#define _C_qmQmQppqp_LLT C2q2Q1g_3445_LLT
#define _C_qmQpQmmqp_LLT C2q2Q1g_2767_LLT
#define _C_qmQpQmpqp_LLT C2q2Q1g_3415_LLT
#define _C_qmqpmQmQp_nfLLT C2q2Q1g_7357_nfLLT
#define _C_qmqppQmQp_nfLLT C2q2Q1g_7465_nfLLT
#define _C_qmqpmQpQm_nfLLT C2q2Q1g_6277_nfLLT
#define _C_qmqppQpQm_nfLLT C2q2Q1g_6385_nfLLT
#define _C_qmqpQmmQp_nfLLT C2q2Q1g_6637_nfLLT
#define _C_qmqpQmpQp_nfLLT C2q2Q1g_7285_nfLLT
#define _C_qmqpQpmQm_nfLLT C2q2Q1g_5377_nfLLT
#define _C_qmqpQppQm_nfLLT C2q2Q1g_6025_nfLLT
#define _C_qmqpQmQpm_nfLLT C2q2Q1g_1237_nfLLT
#define _C_qmqpQmQpp_nfLLT C2q2Q1g_5125_nfLLT
#define _C_qmqpQpQmm_nfLLT C2q2Q1g_1057_nfLLT
#define _C_qmqpQpQmp_nfLLT C2q2Q1g_4945_nfLLT
#define _C_qmmqpQmQp_nfLLT C2q2Q1g_7417_nfLLT
#define _C_qmpqpQmQp_nfLLT C2q2Q1g_7435_nfLLT
#define _C_qmmqpQpQm_nfLLT C2q2Q1g_6337_nfLLT
#define _C_qmpqpQpQm_nfLLT C2q2Q1g_6355_nfLLT
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmqpmQmQp_LLT case 7357
 
#define _CASE_qmqppQmQp_LLT case 7465
 
#define _CASE_qmqpmQpQm_LLT case 6277
 
#define _CASE_qmqppQpQm_LLT case 6385
 
#define _CASE_qmqpQmmQp_LLT case 6637
 
#define _CASE_qmqpQmpQp_LLT case 7285
 
#define _CASE_qmqpQpmQm_LLT case 5377
 
#define _CASE_qmqpQppQm_LLT case 6025
 
#define _CASE_qmqpQmQpm_LLT case 1237
 
#define _CASE_qmqpQmQpp_LLT case 5125
 
#define _CASE_qmqpQpQmm_LLT case 1057
 
#define _CASE_qmqpQpQmp_LLT case 4945
 
#define _CASE_qmmqpQmQp_LLT case 7417
 
#define _CASE_qmpqpQmQp_LLT case 7435
 
#define _CASE_qmmqpQpQm_LLT case 6337
 
#define _CASE_qmpqpQpQm_LLT case 6355
 
#define _CASE_qmqpmQmQp_LRT case 7357
 
#define _CASE_qmqppQmQp_LRT case 7465
 
#define _CASE_qmqpmQpQm_LRT case 6277
 
#define _CASE_qmqppQpQm_LRT case 6385
 
#define _CASE_qmqpQmmQp_LRT case 6637
 
#define _CASE_qmqpQmpQp_LRT case 7285
 
#define _CASE_qmqpQpmQm_LRT case 5377
 
#define _CASE_qmqpQppQm_LRT case 6025
 
#define _CASE_qmqpQmQpm_LRT case 1237
 
#define _CASE_qmqpQmQpp_LRT case 5125
 
#define _CASE_qmqpQpQmm_LRT case 1057
 
#define _CASE_qmqpQpQmp_LRT case 4945
 
#define _CASE_qmmqpQmQp_LRT case 7417
 
#define _CASE_qmpqpQmQp_LRT case 7435
 
#define _CASE_qmmqpQpQm_LRT case 6337
 
#define _CASE_qmpqpQpQm_LRT case 6355
 
#define _CASE_qmmQmQpqp_LLT case 3817
 
#define _CASE_qmmQpQmqp_LLT case 3637
 
#define _CASE_qmpQmQpqp_LLT case 3835
 
#define _CASE_qmpQpQmqp_LLT case 3655
 
#define _CASE_qmQmmQpqp_LLT case 3697
 
#define _CASE_qmQpmQmqp_LLT case 3487
 
#define _CASE_qmQmpQpqp_LLT case 3805
 
#define _CASE_qmQppQmqp_LLT case 3595
 
#define _CASE_qmQmQpqpm_LLT case 637
 
#define _CASE_qmQmQpqpp_LLT case 4525
 
#define _CASE_qmQpQmqpm_LLT case 607
 
#define _CASE_qmQpQmqpp_LLT case 4495
 
#define _CASE_qmQmQpmqp_LLT case 2797
 
#define _CASE_qmQmQppqp_LLT case 3445
 
#define _CASE_qmQpQmmqp_LLT case 2767
 
#define _CASE_qmQpQmpqp_LLT case 3415
 
#define _CASE_qmqpmQmQp_nfLLT case 7357
 
#define _CASE_qmqppQmQp_nfLLT case 7465
 
#define _CASE_qmqpmQpQm_nfLLT case 6277
 
#define _CASE_qmqppQpQm_nfLLT case 6385
 
#define _CASE_qmqpQmmQp_nfLLT case 6637
 
#define _CASE_qmqpQmpQp_nfLLT case 7285
 
#define _CASE_qmqpQpmQm_nfLLT case 5377
 
#define _CASE_qmqpQppQm_nfLLT case 6025
 
#define _CASE_qmqpQmQpm_nfLLT case 1237
 
#define _CASE_qmqpQmQpp_nfLLT case 5125
 
#define _CASE_qmqpQpQmm_nfLLT case 1057
 
#define _CASE_qmqpQpQmp_nfLLT case 4945
 
#define _CASE_qmmqpQmQp_nfLLT case 7417
 
#define _CASE_qmpqpQmQp_nfLLT case 7435
 
#define _CASE_qmmqpQpQm_nfLLT case 6337
 
#define _CASE_qmpqpQpQm_nfLLT case 6355
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2Q1g_LLT( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_qmqpmQmQp_LLT: return new 
                       C2q2Q1g_qmqpmQmQp_LLT_wCI(ind);
    _CASE_qmqppQmQp_LLT: return new 
                       C2q2Q1g_qmqppQmQp_LLT_wCI(ind);
    _CASE_qmqpmQpQm_LLT: return new 
                       C2q2Q1g_qmqpmQpQm_LLT_wCI(ind);
    _CASE_qmqppQpQm_LLT: return new 
                       C2q2Q1g_qmqppQpQm_LLT_wCI(ind);
    _CASE_qmqpQmmQp_LLT: return new 
                       C2q2Q1g_qmqpQmmQp_LLT_wCI(ind);
    _CASE_qmqpQmpQp_LLT: return new 
                       C2q2Q1g_qmqpQmpQp_LLT_wCI(ind);
    _CASE_qmqpQpmQm_LLT: return new 
                       C2q2Q1g_qmqpQpmQm_LLT_wCI(ind);
    _CASE_qmqpQppQm_LLT: return new 
                       C2q2Q1g_qmqpQppQm_LLT_wCI(ind);
    _CASE_qmqpQmQpm_LLT: return new 
                       C2q2Q1g_qmqpQmQpm_LLT_wCI(ind);
    _CASE_qmqpQmQpp_LLT: return new 
                       C2q2Q1g_qmqpQmQpp_LLT_wCI(ind);
    _CASE_qmqpQpQmm_LLT: return new 
                       C2q2Q1g_qmqpQpQmm_LLT_wCI(ind);
    _CASE_qmqpQpQmp_LLT: return new 
                       C2q2Q1g_qmqpQpQmp_LLT_wCI(ind);
    _CASE_qmmqpQmQp_LLT: return new 
                       C2q2Q1g_qmmqpQmQp_LLT_wCI(ind);
    _CASE_qmpqpQmQp_LLT: return new 
                       C2q2Q1g_qmpqpQmQp_LLT_wCI(ind);
    _CASE_qmmqpQpQm_LLT: return new 
                       C2q2Q1g_qmmqpQpQm_LLT_wCI(ind);
    _CASE_qmpqpQpQm_LLT: return new 
                       C2q2Q1g_qmpqpQpQm_LLT_wCI(ind);
    _CASE_qmmQmQpqp_LLT: return new 
                       C2q2Q1g_qmmQmQpqp_LLT_wCI(ind);
    _CASE_qmmQpQmqp_LLT: return new 
                       C2q2Q1g_qmmQpQmqp_LLT_wCI(ind);
    _CASE_qmpQmQpqp_LLT: return new 
                       C2q2Q1g_qmpQmQpqp_LLT_wCI(ind);
    _CASE_qmpQpQmqp_LLT: return new 
                       C2q2Q1g_qmpQpQmqp_LLT_wCI(ind);
    _CASE_qmQmmQpqp_LLT: return new 
                       C2q2Q1g_qmQmmQpqp_LLT_wCI(ind);
    _CASE_qmQpmQmqp_LLT: return new 
                       C2q2Q1g_qmQpmQmqp_LLT_wCI(ind);
    _CASE_qmQmpQpqp_LLT: return new 
                       C2q2Q1g_qmQmpQpqp_LLT_wCI(ind);
    _CASE_qmQppQmqp_LLT: return new 
                       C2q2Q1g_qmQppQmqp_LLT_wCI(ind);
    _CASE_qmQmQpqpm_LLT: return new 
                       C2q2Q1g_qmQmQpqpm_LLT_wCI(ind);
    _CASE_qmQmQpqpp_LLT: return new 
                       C2q2Q1g_qmQmQpqpp_LLT_wCI(ind);
    _CASE_qmQpQmqpm_LLT: return new 
                       C2q2Q1g_qmQpQmqpm_LLT_wCI(ind);
    _CASE_qmQpQmqpp_LLT: return new 
                       C2q2Q1g_qmQpQmqpp_LLT_wCI(ind);
    _CASE_qmQmQpmqp_LLT: return new 
                       C2q2Q1g_qmQmQpmqp_LLT_wCI(ind);
    _CASE_qmQmQppqp_LLT: return new 
                       C2q2Q1g_qmQmQppqp_LLT_wCI(ind);
    _CASE_qmQpQmmqp_LLT: return new 
                       C2q2Q1g_qmQpQmmqp_LLT_wCI(ind);
    _CASE_qmQpQmpqp_LLT: return new 
                       C2q2Q1g_qmQpQmpqp_LLT_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q1g_LRT( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_qmqpmQmQp_LRT: return new 
                       C2q2Q1g_qmqpmQmQp_LRT_wCI(ind);
    _CASE_qmqppQmQp_LRT: return new 
                       C2q2Q1g_qmqppQmQp_LRT_wCI(ind);
    _CASE_qmqpmQpQm_LRT: return new 
                       C2q2Q1g_qmqpmQpQm_LRT_wCI(ind);
    _CASE_qmqppQpQm_LRT: return new 
                       C2q2Q1g_qmqppQpQm_LRT_wCI(ind);
    _CASE_qmqpQmmQp_LRT: return new 
                       C2q2Q1g_qmqpQmmQp_LRT_wCI(ind);
    _CASE_qmqpQmpQp_LRT: return new 
                       C2q2Q1g_qmqpQmpQp_LRT_wCI(ind);
    _CASE_qmqpQpmQm_LRT: return new 
                       C2q2Q1g_qmqpQpmQm_LRT_wCI(ind);
    _CASE_qmqpQppQm_LRT: return new 
                       C2q2Q1g_qmqpQppQm_LRT_wCI(ind);
    _CASE_qmqpQmQpm_LRT: return new 
                       C2q2Q1g_qmqpQmQpm_LRT_wCI(ind);
    _CASE_qmqpQmQpp_LRT: return new 
                       C2q2Q1g_qmqpQmQpp_LRT_wCI(ind);
    _CASE_qmqpQpQmm_LRT: return new 
                       C2q2Q1g_qmqpQpQmm_LRT_wCI(ind);
    _CASE_qmqpQpQmp_LRT: return new 
                       C2q2Q1g_qmqpQpQmp_LRT_wCI(ind);
    _CASE_qmmqpQmQp_LRT: return new 
                       C2q2Q1g_qmmqpQmQp_LRT_wCI(ind);
    _CASE_qmpqpQmQp_LRT: return new 
                       C2q2Q1g_qmpqpQmQp_LRT_wCI(ind);
    _CASE_qmmqpQpQm_LRT: return new 
                       C2q2Q1g_qmmqpQpQm_LRT_wCI(ind);
    _CASE_qmpqpQpQm_LRT: return new 
                       C2q2Q1g_qmpqpQpQm_LRT_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2Q1g_nfLLT( int hc,const std::vector<int>& ind) 
{ 
    switch (hc) {
    _CASE_qmqpmQmQp_nfLLT: return new 
                       C2q2Q1g_qmqpmQmQp_nfLLT_wCI(ind);
    _CASE_qmqppQmQp_nfLLT: return new 
                       C2q2Q1g_qmqppQmQp_nfLLT_wCI(ind);
    _CASE_qmqpmQpQm_nfLLT: return new 
                       C2q2Q1g_qmqpmQpQm_nfLLT_wCI(ind);
    _CASE_qmqppQpQm_nfLLT: return new 
                       C2q2Q1g_qmqppQpQm_nfLLT_wCI(ind);
    _CASE_qmqpQmmQp_nfLLT: return new 
                       C2q2Q1g_qmqpQmmQp_nfLLT_wCI(ind);
    _CASE_qmqpQmpQp_nfLLT: return new 
                       C2q2Q1g_qmqpQmpQp_nfLLT_wCI(ind);
    _CASE_qmqpQpmQm_nfLLT: return new 
                       C2q2Q1g_qmqpQpmQm_nfLLT_wCI(ind);
    _CASE_qmqpQppQm_nfLLT: return new 
                       C2q2Q1g_qmqpQppQm_nfLLT_wCI(ind);
    _CASE_qmqpQmQpm_nfLLT: return new 
                       C2q2Q1g_qmqpQmQpm_nfLLT_wCI(ind);
    _CASE_qmqpQmQpp_nfLLT: return new 
                       C2q2Q1g_qmqpQmQpp_nfLLT_wCI(ind);
    _CASE_qmqpQpQmm_nfLLT: return new 
                       C2q2Q1g_qmqpQpQmm_nfLLT_wCI(ind);
    _CASE_qmqpQpQmp_nfLLT: return new 
                       C2q2Q1g_qmqpQpQmp_nfLLT_wCI(ind);
    _CASE_qmmqpQmQp_nfLLT: return new 
                       C2q2Q1g_qmmqpQmQp_nfLLT_wCI(ind);
    _CASE_qmpqpQmQp_nfLLT: return new 
                       C2q2Q1g_qmpqpQmQp_nfLLT_wCI(ind);
    _CASE_qmmqpQpQm_nfLLT: return new 
                       C2q2Q1g_qmmqpQpQm_nfLLT_wCI(ind);
    _CASE_qmpqpQpQm_nfLLT: return new 
                       C2q2Q1g_qmpqpQpQm_nfLLT_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
