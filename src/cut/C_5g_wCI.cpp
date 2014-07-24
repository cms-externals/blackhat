/*
*C_5g_wCI.cpp
*
* Created on 10/19, 2010
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
 

 
class C5g_mmppp_G_wCI : public Cut_Part_wCI {
public:
       C5g_mmppp_G_wCI
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

 
class C5g_mpmpp_G_wCI : public Cut_Part_wCI {
public:
       C5g_mpmpp_G_wCI
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

 
class C5g_mppmp_G_wCI : public Cut_Part_wCI {
public:
       C5g_mppmp_G_wCI
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

 
class C5g_mpppm_G_wCI : public Cut_Part_wCI {
public:
       C5g_mpppm_G_wCI
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

 
class C5g_pmmpp_G_wCI : public Cut_Part_wCI {
public:
       C5g_pmmpp_G_wCI
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

 
class C5g_pmpmp_G_wCI : public Cut_Part_wCI {
public:
       C5g_pmpmp_G_wCI
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

 
class C5g_pmppm_G_wCI : public Cut_Part_wCI {
public:
       C5g_pmppm_G_wCI
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

 
class C5g_ppmmp_G_wCI : public Cut_Part_wCI {
public:
       C5g_ppmmp_G_wCI
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

 
class C5g_ppmpm_G_wCI : public Cut_Part_wCI {
public:
       C5g_ppmpm_G_wCI
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

 
class C5g_pppmm_G_wCI : public Cut_Part_wCI {
public:
       C5g_pppmm_G_wCI
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

 
class C5g_mmppp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mmppp_nf_wCI
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

 
class C5g_mpmpp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mpmpp_nf_wCI
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

 
class C5g_mppmp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mppmp_nf_wCI
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

 
class C5g_mpppm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mpppm_nf_wCI
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

 
class C5g_pmmpp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_pmmpp_nf_wCI
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

 
class C5g_pmpmp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_pmpmp_nf_wCI
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

 
class C5g_pmppm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_pmppm_nf_wCI
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

 
class C5g_ppmmp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_ppmmp_nf_wCI
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

 
class C5g_ppmpm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_ppmpm_nf_wCI
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

 
class C5g_pppmm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_pppmm_nf_wCI
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

 
class C5g_ppmmm_G_wCI : public Cut_Part_wCI {
public:
       C5g_ppmmm_G_wCI
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

 
class C5g_pmpmm_G_wCI : public Cut_Part_wCI {
public:
       C5g_pmpmm_G_wCI
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

 
class C5g_pmmpm_G_wCI : public Cut_Part_wCI {
public:
       C5g_pmmpm_G_wCI
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

 
class C5g_pmmmp_G_wCI : public Cut_Part_wCI {
public:
       C5g_pmmmp_G_wCI
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

 
class C5g_mppmm_G_wCI : public Cut_Part_wCI {
public:
       C5g_mppmm_G_wCI
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

 
class C5g_mpmpm_G_wCI : public Cut_Part_wCI {
public:
       C5g_mpmpm_G_wCI
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

 
class C5g_mpmmp_G_wCI : public Cut_Part_wCI {
public:
       C5g_mpmmp_G_wCI
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

 
class C5g_mmppm_G_wCI : public Cut_Part_wCI {
public:
       C5g_mmppm_G_wCI
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

 
class C5g_mmpmp_G_wCI : public Cut_Part_wCI {
public:
       C5g_mmpmp_G_wCI
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

 
class C5g_mmmpp_G_wCI : public Cut_Part_wCI {
public:
       C5g_mmmpp_G_wCI
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

 
class C5g_ppmmm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_ppmmm_nf_wCI
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

 
class C5g_pmpmm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_pmpmm_nf_wCI
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

 
class C5g_pmmpm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_pmmpm_nf_wCI
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

 
class C5g_pmmmp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_pmmmp_nf_wCI
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

 
class C5g_mppmm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mppmm_nf_wCI
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

 
class C5g_mpmpm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mpmpm_nf_wCI
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

 
class C5g_mpmmp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mpmmp_nf_wCI
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

 
class C5g_mmppm_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mmppm_nf_wCI
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

 
class C5g_mmpmp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mmpmp_nf_wCI
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

 
class C5g_mmmpp_nf_wCI : public Cut_Part_wCI {
public:
       C5g_mmmpp_nf_wCI
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


C5g_mmppp_G_wCI::\
C5g_mmppp_G_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mmppp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mmppp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:45:40 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa24 = SPA(2,4);
complex<T> s12 = S(1,2);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> t1 = spa34*spa45; 
complex<T> t5 = square(spa12); 
complex<T> t6 = -(spa14*spb34); 
complex<T> t10 = cube(spa12); 
complex<T> t14 = spa24*spb45; 
complex<T> d5 = spa15*spa45*T(2); d5 = T(1)/d5;
complex<T> d6 = spa15*spa23*T(2); d6 = T(1)/d6;
complex<T> d7 = spa23*spa34*T(2); d7 = T(1)/d7;
complex<T> t7 = -t14; 
complex<T> d1 = spa15*spa23*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*spa15*spa23*t1*T(6); d2 = T(1)/d2;
complex<T> d3 = t1*cube(s15 - s23)*T(3); d3 = T(1)/d3;
complex<T> d4 = spa15*t1*T(2); d4 = T(1)/d4;
complex<T> d8 = spa23*t1*T(2); d8 = T(1)/d8;
complex<T> t11 = spa23*t6 + spa15*t7; 
complex<T> t13 = d1*T(11); 
complex<T> t17 = d3*t11; 
complex<T> t3 = t10*t13 + t14*t17*t6 - d2*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa14*spb34*t14*t17 + d2*t11*t5*T(11); 
complex<T> co1 = -(d4*s12*spb23*t10); 
complex<T> co2 = d5*spb23*spb34*t10; 
complex<T> co3 = d6*spb34*spb45*t10; 
complex<T> co4 = d7*spb15*spb45*t10; 
complex<T> co5 = -(d8*s12*spb15*t10); 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mpmpp_G_wCI::\
C5g_mpmpp_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mpmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mpmpp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:47:06 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = -(spa12*spb12);
complex<T> t5 = square(spa25); 
complex<T> t6 = spa15*spa45; 
complex<T> t12 = spa12*spa23; 
complex<T> t13 = cube(spa13); 
complex<T> t14 = spa34*spa45; 
complex<T> t23 = square(spb25); 
complex<T> t24 = square(spa13); 
complex<T> t25 = square(spb24); 
complex<T> t26 = square(spa14); 
complex<T> t27 = square(spa35); 
complex<T> t30 = s15 - s23; 
complex<T> t32 = s12 + s15 - s34; 
complex<T> t33 = -(spa12*spa34*spb24) - spa15*spa23*spb25; 
complex<T> t34 = square(spa15); 
complex<T> t36 = square(spa34); 
complex<T> t45 = cube(spb24); 
complex<T> t47 = cube(spa14); 
complex<T> t48 = cube(spa35); 
complex<T> t55 = square(square(spa13)); 
complex<T> t57 = square(spa12); 
complex<T> t58 = square(spa23); 
complex<T> t59 = cube(spb25); 
complex<T> t76 = -(spa14*T(4)); 
complex<T> d19 = spa45*T(3); d19 = T(1)/d19;
complex<T> d20 = cube(s15 - s34); d20 = T(1)/d20;
complex<T> d23 = spa45*square(spa24); d23 = T(1)/d23;
complex<T> d24 = spa45*square(square(spa24)); d24 = T(1)/d24;
complex<T> t4 = square(t32); 
complex<T> t22 = t32*t5; 
complex<T> t28 = -(t12*T(2)); 
complex<T> t42 = -(t24*T(2)); 
complex<T> t63 = t12*T(2); 
complex<T> t64 = spa15*t23; 
complex<T> t65 = spa34*t26; 
complex<T> t66 = spa35*t24; 
complex<T> t73 = t24*t25; 
complex<T> t74 = spa34*t12; 
complex<T> t85 = t34*t58; 
complex<T> t88 = t25*t26; 
complex<T> t91 = spa13*t33; 
complex<T> t93 = t45*t57; 
complex<T> t96 = spb12*t55; 
complex<T> d1 = (s12 - s34)*t14*t32; d1 = T(1)/d1;
complex<T> d2 = t14*t5*square(s12 - s34); d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*t14*t32*t5; d3 = T(1)/d3;
complex<T> d4 = spa25*t14*cube(s12 - s34)*T(3); d4 = T(1)/d4;
complex<T> d6 = s24*t30*t6; d6 = T(1)/d6;
complex<T> d7 = s24*(s15 - s34)*t6; d7 = T(1)/d7;
complex<T> d8 = t6*square(spa24)*square(t30); d8 = T(1)/d8;
complex<T> d9 = s24*t30*t6*square(spa24); d9 = T(1)/d9;
complex<T> d10 = t6*square(s15 - s34)*square(spa24); d10 = T(1)/d10;
complex<T> d11 = s24*(s15 - s34)*t6*square(spa24); d11 = T(1)/d11;
complex<T> d12 = spa24*t6*cube(t30)*T(3); d12 = T(1)/d12;
complex<T> d13 = (s15 - s34)*t14*t32; d13 = T(1)/d13;
complex<T> d14 = t14*t5*square(s15 - s34); d14 = T(1)/d14;
complex<T> d15 = (s15 - s34)*t14*t32*t5; d15 = T(1)/d15;
complex<T> d17 = spa24*t6*T(3); d17 = T(1)/d17;
complex<T> d18 = spa25*t14*T(3); d18 = T(1)/d18;
complex<T> d25 = t6*square(spa24); d25 = T(1)/d25;
complex<T> d26 = t6*square(square(spa24)); d26 = T(1)/d26;
complex<T> d29 = spa34*t6*T(2); d29 = T(1)/d29;
complex<T> d30 = spa12*t6*T(2); d30 = T(1)/d30;
complex<T> d33 = spa23*t14*T(2); d33 = T(1)/d33;
complex<T> t40 = d25*s23; 
complex<T> t41 = d26*s34; 
complex<T> t43 = -t74; 
complex<T> t44 = -t64; 
complex<T> t60 = d18*spa35; 
complex<T> t61 = d17*t36; 
complex<T> t62 = d12*t45; 
complex<T> t71 = t27*t64; 
complex<T> t77 = t23*t66; 
complex<T> t83 = t25*t65; 
complex<T> t84 = d4*t57; 
complex<T> d5 = t6*t74*T(6); d5 = T(1)/d5;
complex<T> d16 = (s15 - s34)*t6*t74*T(6); d16 = T(1)/d16;
complex<T> d21 = t14*t4; d21 = T(1)/d21;
complex<T> d22 = t14*t4*t5; d22 = T(1)/d22;
complex<T> d27 = spa45*t4; d27 = T(1)/d27;
complex<T> d28 = spa45*t4*t5; d28 = T(1)/d28;
complex<T> d31 = spa15*t63; d31 = T(1)/d31;
complex<T> d32 = spa34*t63; d32 = T(1)/d32;
complex<T> t18 = t24*t40*t76 + d26*s23*spa12*spa23*t65*T(2); 
complex<T> t21 = d2*spa12*spa23*t27*t44 - d3*spa12*spa23*t71*T(2) + t48*t59*t84*T(2) + d1*t77*T(4); 
complex<T> t37 = d5*T(11); 
complex<T> t38 = d21*s12; 
complex<T> t39 = d22*s15; 
complex<T> t54 = d9*t28*t83 + d8*t43*t88 - t47*t58*t62*T(2) + d6*spa14*t73*T(4); 
complex<T> t69 = s34*spa14*t40*t42 + s23*t12*t41*t65; 
complex<T> t81 = t41*t63*t65 + d28*spb34*t63*t71 + d25*s34*t24*t76 - d27*spb34*t77*T(4); 
complex<T> t1 = t37*t55 + d14*t12*t71 + d2*t12*t71 + d15*t63*t71 + d3*t63*t71 + d11*t28*t83 + d10*t43*t88 - d19*d20*spb24*spb25*t91 - t48*t59*t84*T(2) - d20*t59*t60*t85*T(2) - d20*spa14*t61*t93*T(2) + d7*spa14*t73*T(4) - d1*t77*T(4) - d13*t77*T(4) - d16*t13*t33*T(11); 
complex<T> t2 = d14*t12*t27*t44 + t37*t55 + d15*t28*t71 + d6*t73*t76 + d7*t73*t76 + d11*t63*t83 + d9*t63*t83 + d10*t74*t88 + d8*t74*t88 + d19*d20*spb24*spb25*t91 + t47*t58*t62*T(2) + d20*t59*t60*t85*T(2) + d20*spa14*t61*t93*T(2) + d13*t77*T(4) + d16*t13*t33*T(11); 
complex<T> t20 = s15*spa35*t23*t38*t42 + s12*spa12*spa23*t39*t71 + d33*spb15*t96; 
complex<T> t70 = d24*spb15*t63*t65 + t39*t63*t71 + d23*spb15*t24*t76 - d21*s15*t77*T(4); 
complex<T> t75 = d22*s12*t63*t71 - t38*t77*T(4); 
complex<T> co1 = d29*spb23*t96; 
complex<T> co2 = d30*spb23*spb34*t55; 
complex<T> co3 = d31*spb34*spb45*t55; 
complex<T> co4 = d32*spb15*spb45*t55; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t21*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t54*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t75*(*CI_users[4]->get_value(mc,ind,mu)) + t70*(*CI_users[5]->get_value(mc,ind,mu)) + t18*(*CI_users[6]->get_value(mc,ind,mu)) + t81*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + co3*(*CI_users[10]->get_value(mc,ind,mu)) + t69*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + t20*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mppmp_G_wCI::\
C5g_mppmp_G_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C5g_mppmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mppmp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:48:20 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t5 = spa12*spa23; 
complex<T> t6 = square(spa25); 
complex<T> t12 = spa15*spa45; 
complex<T> t13 = cube(spa14); 
complex<T> t14 = spa23*spa34; 
complex<T> t22 = square(spb25); 
complex<T> t23 = square(spa14); 
complex<T> t24 = square(spb35); 
complex<T> t25 = square(spa13); 
complex<T> t26 = square(spa24); 
complex<T> t28 = s12 - s34; 
complex<T> t29 = s15 - s34; 
complex<T> t30 = s12 + s15 - s34; 
complex<T> t31 = spa12*spa45*spb25 + spa15*spa34*spb35; 
complex<T> t32 = s12 - s45; 
complex<T> t33 = square(spa12); 
complex<T> t34 = square(spa34); 
complex<T> t45 = cube(spb35); 
complex<T> t55 = square(square(spa14)); 
complex<T> t57 = square(spa15); 
complex<T> t58 = square(spa45); 
complex<T> t59 = cube(spb25); 
complex<T> t60 = cube(spa24); 
complex<T> t70 = cube(spa13); 
complex<T> t78 = -(spa13*T(4)); 
complex<T> t92 = spa14*spb25; 
complex<T> d15 = spa23*T(3); d15 = T(1)/d15;
complex<T> d21 = spa23*square(spa35); d21 = T(1)/d21;
complex<T> d22 = spa23*square(square(spa35)); d22 = T(1)/d22;
complex<T> t27 = -(t12*T(2)); 
complex<T> t64 = t12*T(2); 
complex<T> t65 = spa12*t22; 
complex<T> t66 = spa34*t25; 
complex<T> t67 = spa24*t23; 
complex<T> t75 = t12*t24; 
complex<T> t89 = t23*t24; 
complex<T> t98 = spb34*t55; 
complex<T> t101 = spa24*t59; 
complex<T> t107 = t57*t59; 
complex<T> d1 = spa34*t12*t5*T(6); d1 = T(1)/d1;
complex<T> d2 = t14*t28*(s15 + t28); d2 = T(1)/d2;
complex<T> d3 = t14*t6*square(t28); d3 = T(1)/d3;
complex<T> d4 = t14*t28*(s15 + t28)*t6; d4 = T(1)/d4;
complex<T> d5 = s35*t28*t5; d5 = T(1)/d5;
complex<T> d6 = s35*t32*t5; d6 = T(1)/d6;
complex<T> d7 = t5*square(spa35)*square(t28); d7 = T(1)/d7;
complex<T> d8 = s35*t28*t5*square(spa35); d8 = T(1)/d8;
complex<T> d9 = t5*square(spa35)*square(t32); d9 = T(1)/d9;
complex<T> d10 = s35*t32*t5*square(spa35); d10 = T(1)/d10;
complex<T> d11 = spa35*t5*cube(t32)*T(3); d11 = T(1)/d11;
complex<T> d12 = spa34*t12*t28*t5*T(6); d12 = T(1)/d12;
complex<T> d13 = spa25*t14*T(3); d13 = T(1)/d13;
complex<T> d14 = spa35*t5*T(3); d14 = T(1)/d14;
complex<T> d16 = cube(t28); d16 = T(1)/d16;
complex<T> d17 = t14*(s15 + t28)*t29; d17 = T(1)/d17;
complex<T> d18 = t14*t6*square(t29); d18 = T(1)/d18;
complex<T> d19 = t14*(s15 + t28)*t29*t6; d19 = T(1)/d19;
complex<T> d20 = spa25*t14*cube(t29)*T(3); d20 = T(1)/d20;
complex<T> d23 = t14*square(s15 + t28); d23 = T(1)/d23;
complex<T> d24 = t14*t6*square(s15 + t28); d24 = T(1)/d24;
complex<T> d25 = t5*square(spa35); d25 = T(1)/d25;
complex<T> d26 = t5*square(square(spa35)); d26 = T(1)/d26;
complex<T> d27 = spa23*square(s15 + t28); d27 = T(1)/d27;
complex<T> d28 = spa23*t6*square(s15 + t28); d28 = T(1)/d28;
complex<T> d31 = spa15*t5*T(2); d31 = T(1)/d31;
complex<T> d32 = spa34*t5*T(2); d32 = T(1)/d32;
complex<T> d33 = spa45*t14*T(2); d33 = T(1)/d33;
complex<T> t36 = d1*T(11); 
complex<T> t37 = d23*s12; 
complex<T> t38 = d24*s15; 
complex<T> t39 = d25*s34; 
complex<T> t40 = d26*s45; 
complex<T> t42 = -t65; 
complex<T> t43 = -t66; 
complex<T> t61 = d13*t33; 
complex<T> t62 = d14*t34; 
complex<T> t63 = d11*t45; 
complex<T> t72 = t26*t65; 
complex<T> t79 = t22*t67; 
complex<T> t85 = t24*t66; 
complex<T> d29 = spa34*t64; d29 = T(1)/d29;
complex<T> d30 = spa12*t64; d30 = T(1)/d30;
complex<T> t18 = d25*s45*t23*t78 + spa15*spa45*t40*t66*T(2); 
complex<T> t19 = s12*spa15*spa45*t38*t72 - s15*t37*t79*T(2); 
complex<T> t21 = s34*spa15*spa45*t40*t66 + d31*spb45*t98 - s45*spa13*t23*t39*T(2); 
complex<T> t52 = d18*t12*t26*t42 + d19*t27*t72 + d20*t107*t60*T(2) + d17*t79*T(4); 
complex<T> t53 = d9*t43*t75 + d10*t27*t85 - t58*t63*t70*T(2) + d6*spa13*t89*T(4); 
complex<T> t84 = -(t79*T(4)); 
complex<T> t95 = t45*t62; 
complex<T> t1 = d3*t12*t26*t42 + t36*t55 + d4*t27*t72 + d7*t66*t75 + d9*t66*t75 + d10*t64*t85 + d8*t64*t85 + d5*t78*t89 + d6*t78*t89 - d15*d16*spb35*t31*t92 + d16*t101*t58*t61*T(2) + t58*t63*t70*T(2) + d16*spa13*t57*t95*T(2) + d2*t79*T(4) - d12*t13*t31*T(11); 
complex<T> t2 = t36*t55 + d18*t12*t72 + d3*t12*t72 + d19*t64*t72 + d4*t64*t72 + d7*t43*t75 + d17*t84 + d2*t84 + d8*t27*t85 + d15*d16*spb35*t31*t92 - d20*t107*t60*T(2) - d16*t101*t58*t61*T(2) - d16*spa13*t57*t95*T(2) + d5*spa13*t89*T(4) + d12*t13*t31*T(11); 
complex<T> t77 = d22*spb12*t64*t66 + d24*s12*t64*t72 + d21*spb12*t23*t78 + t37*t84; 
complex<T> t83 = t38*t64*t72 + d23*s15*t84; 
complex<T> t87 = d26*s34*t64*t66 + d28*spb34*t64*t72 + t23*t39*t78 + d27*spb34*t84; 
complex<T> co1 = d29*spb12*spb23*t55; 
complex<T> co2 = d30*spb23*t98; 
complex<T> co3 = d32*spb15*spb45*t55; 
complex<T> co4 = d33*spb12*spb15*t55; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t52*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t53*(*CI_users[3]->get_value(mc,ind,mu)) + t77*(*CI_users[4]->get_value(mc,ind,mu)) + t83*(*CI_users[5]->get_value(mc,ind,mu)) + t87*(*CI_users[6]->get_value(mc,ind,mu)) + t18*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t19*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t21*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mpppm_G_wCI::\
C5g_mpppm_G_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mpppm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mpppm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:48:22 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = spa23*spa34; 
complex<T> t5 = square(spa15); 
complex<T> t10 = cube(spa15); 
complex<T> t14 = spa13*spb34; 
complex<T> d4 = spa34*spa45*T(2); d4 = T(1)/d4;
complex<T> d5 = spa12*spa45*T(2); d5 = T(1)/d5;
complex<T> d6 = spa12*spa23*T(2); d6 = T(1)/d6;
complex<T> t11 = spa12*spa35*spb23 + spa45*t14; 
complex<T> d1 = spa12*spa45*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa12*spa45*t1*T(6); d2 = T(1)/d2;
complex<T> d3 = t1*cube(s12 - s45)*T(3); d3 = T(1)/d3;
complex<T> d7 = spa12*t1*T(2); d7 = T(1)/d7;
complex<T> d8 = spa45*t1*T(2); d8 = T(1)/d8;
complex<T> t13 = d1*T(11); 
complex<T> t20 = d3*t11; 
complex<T> t3 = t10*t13 - spa35*spb23*t14*t20 - d2*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa35*spb23*t14*t20 + d2*t11*t5*T(11); 
complex<T> co1 = d4*spb12*spb23*t10; 
complex<T> co2 = d5*spb23*spb34*t10; 
complex<T> co3 = d6*spb34*spb45*t10; 
complex<T> co4 = -(d7*s15*spb45*t10); 
complex<T> co5 = -(d8*s15*spb12*t10); 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmmpp_G_wCI::\
C5g_pmmpp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmmpp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:48:24 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> s23 = S(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> t1 = spa15*spa45; 
complex<T> t5 = square(spa23); 
complex<T> t6 = -(spa35*spb15); 
complex<T> t10 = cube(spa23); 
complex<T> t14 = spa25*spb45; 
complex<T> d6 = spa12*spa15*T(2); d6 = T(1)/d6;
complex<T> d7 = spa12*spa34*T(2); d7 = T(1)/d7;
complex<T> d8 = spa34*spa45*T(2); d8 = T(1)/d8;
complex<T> t7 = -t14; 
complex<T> d1 = spa12*spa34*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa12*spa34*t1*T(6); d2 = T(1)/d2;
complex<T> d3 = t1*cube(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> d4 = spa34*t1*T(2); d4 = T(1)/d4;
complex<T> d5 = spa12*t1*T(2); d5 = T(1)/d5;
complex<T> t11 = spa12*t6 + spa34*t7; 
complex<T> t13 = d1*T(11); 
complex<T> t17 = d3*t11; 
complex<T> t3 = t10*t13 + t14*t17*t6 - d2*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa35*spb15*t14*t17 + d2*t11*t5*T(11); 
complex<T> co1 = -(d4*s23*spb12*t10); 
complex<T> co2 = -(d5*s23*spb34*t10); 
complex<T> co3 = d6*spb34*spb45*t10; 
complex<T> co4 = d7*spb15*spb45*t10; 
complex<T> co5 = d8*spb12*spb15*t10; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmpmp_G_wCI::\
C5g_pmpmp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
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
     C5g_pmpmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmpmp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:49:27 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = -(spa12*spb12);
complex<T> s13 = -(spa13*spb13);
complex<T> s35 = -(spa35*spb35);
complex<T> t5 = spa12*spa15; 
complex<T> t11 = spa23*spa34; 
complex<T> t12 = cube(spa24); 
complex<T> t21 = square(spa24); 
complex<T> t22 = square(spb13); 
complex<T> t23 = square(spb35); 
complex<T> t24 = square(spa14); 
complex<T> t25 = square(spa25); 
complex<T> t27 = s12 - s45; 
complex<T> t28 = s12 - s34; 
complex<T> t29 = s23 - s45; 
complex<T> t30 = -(spa12*spa34*spb13) - spa23*spa45*spb35; 
complex<T> t31 = square(spa12); 
complex<T> t34 = square(spa45); 
complex<T> t43 = cube(spb13); 
complex<T> t44 = cube(spb35); 
complex<T> t45 = cube(spa25); 
complex<T> t54 = square(square(spa24)); 
complex<T> t56 = square(spa23); 
complex<T> t57 = square(spa34); 
complex<T> t68 = cube(spa14); 
complex<T> t78 = spa24*spb13; 
complex<T> d13 = spa13*spa15*spa45*T(3); d13 = T(1)/d13;
complex<T> d15 = spa15*T(3); d15 = T(1)/d15;
complex<T> d21 = spa15*spa45*square(spa13); d21 = T(1)/d21;
complex<T> d22 = spa15*spa45*square(square(spa13)); d22 = T(1)/d22;
complex<T> d23 = spa15*square(spa35); d23 = T(1)/d23;
complex<T> d24 = spa15*square(square(spa35)); d24 = T(1)/d24;
complex<T> d27 = spa15*square(spa13); d27 = T(1)/d27;
complex<T> d28 = spa15*square(square(spa13)); d28 = T(1)/d28;
complex<T> d29 = spa15*spa34*spa45*T(2); d29 = T(1)/d29;
complex<T> t26 = -(t11*T(2)); 
complex<T> t36 = d21*s12; 
complex<T> t37 = d22*s23; 
complex<T> t41 = -(spa45*t11); 
complex<T> t59 = d13*t31; 
complex<T> t62 = t11*T(2); 
complex<T> t63 = spa12*t24; 
complex<T> t64 = spa45*t25; 
complex<T> t65 = t21*t22; 
complex<T> t70 = -(t21*T(4)); 
complex<T> t76 = spa25*t23; 
complex<T> t89 = t23*t25; 
complex<T> t90 = spb12*t54; 
complex<T> t93 = t34*t44; 
complex<T> d1 = spa45*t11*t5*T(6); d1 = T(1)/d1;
complex<T> d2 = s13*spa15*spa45*t27; d2 = T(1)/d2;
complex<T> d3 = spa15*spa45*square(spa13)*square(t27); d3 = T(1)/d3;
complex<T> d4 = s13*spa15*spa45*t27*square(spa13); d4 = T(1)/d4;
complex<T> d5 = s35*t28*t5; d5 = T(1)/d5;
complex<T> d6 = s35*t27*t5; d6 = T(1)/d6;
complex<T> d7 = t5*square(spa35)*square(t28); d7 = T(1)/d7;
complex<T> d8 = s35*t28*t5*square(spa35); d8 = T(1)/d8;
complex<T> d9 = t5*square(spa35)*square(t27); d9 = T(1)/d9;
complex<T> d10 = s35*t27*t5*square(spa35); d10 = T(1)/d10;
complex<T> d11 = spa35*t5*cube(t28)*T(3); d11 = T(1)/d11;
complex<T> d12 = spa45*t11*t27*t5*T(6); d12 = T(1)/d12;
complex<T> d14 = spa35*t5*T(3); d14 = T(1)/d14;
complex<T> d16 = cube(t27); d16 = T(1)/d16;
complex<T> d17 = s13*spa15*spa45*t29; d17 = T(1)/d17;
complex<T> d18 = spa15*spa45*square(spa13)*square(t29); d18 = T(1)/d18;
complex<T> d19 = s13*spa15*spa45*t29*square(spa13); d19 = T(1)/d19;
complex<T> d20 = spa13*spa15*spa45*cube(t29)*T(3); d20 = T(1)/d20;
complex<T> d25 = t5*square(spa35); d25 = T(1)/d25;
complex<T> d26 = t5*square(square(spa35)); d26 = T(1)/d26;
complex<T> d30 = spa45*t5*T(2); d30 = T(1)/d30;
complex<T> d31 = spa23*t5*T(2); d31 = T(1)/d31;
complex<T> t35 = d1*T(11); 
complex<T> t38 = d25*s34; 
complex<T> t39 = d26*s45; 
complex<T> t42 = -t63; 
complex<T> t58 = d14*spa25; 
complex<T> t60 = d20*t43; 
complex<T> t61 = d11*t44; 
complex<T> t71 = t22*t63; 
complex<T> t72 = t23*t64; 
complex<T> t79 = d22*s12*t62*t63 + d24*spb12*t62*t64 + d23*spa25*spb12*t70 + spa14*t36*t70; 
complex<T> t85 = t43*t59; 
complex<T> t92 = spa23*t37; 
complex<T> d32 = spa12*t62; d32 = T(1)/d32;
complex<T> d33 = spa45*t62; d33 = T(1)/d33;
complex<T> t1 = t35*t54 + d18*t11*t71 + d3*t11*t71 + d19*t26*t71 + d4*t26*t71 + d10*t26*t72 - d15*d16*spb35*t30*t78 + d9*t41*t89 - t56*t60*t68*T(2) - d16*spa14*t57*t85*T(2) - d16*t56*t58*t93*T(2) + d17*spa14*t65*T(4) + d2*spa14*t65*T(4) + d6*t21*t76*T(4) - d12*t12*t30*T(11); 
complex<T> t2 = d3*t11*t22*t42 + t35*t54 + d4*t62*t71 + d7*t11*t72 + d9*t11*t72 + d10*t62*t72 + d8*t62*t72 + d5*t70*t76 + d6*t70*t76 + d15*d16*spb35*t30*t78 + t45*t57*t61*T(2) + d16*spa14*t57*t85*T(2) + d16*t56*t58*t93*T(2) - d2*spa14*t65*T(4) + d12*t12*t30*T(11); 
complex<T> t17 = d21*s23*spa14*t70 + spa34*t63*t92*T(2); 
complex<T> t18 = s34*spa23*spa34*t39*t64 - s45*spa25*t21*t38*T(2); 
complex<T> t19 = d29*spb23*t90 + s12*spa34*t63*t92 - s23*spa14*t21*t36*T(2); 
complex<T> t51 = d18*t11*t22*t42 + d19*t62*t71 + t56*t60*t68*T(2) - d17*spa14*t65*T(4); 
complex<T> t53 = d8*t26*t72 + d7*t41*t89 - t45*t57*t61*T(2) + d5*t21*t76*T(4); 
complex<T> t74 = d28*spb45*t62*t63 + t39*t62*t64 + d25*s45*spa25*t70 + d27*spa14*spb45*t70; 
complex<T> t86 = d26*s34*t62*t64 + spa25*t38*t70; 
complex<T> co1 = d30*spb23*spb34*t54; 
complex<T> co2 = d31*spb34*spb45*t54; 
complex<T> co3 = d32*spb15*spb45*t54; 
complex<T> co4 = d33*spb15*t90; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t51*(*CI_users[1]->get_value(mc,ind,mu)) + t53*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t79*(*CI_users[4]->get_value(mc,ind,mu)) + t17*(*CI_users[5]->get_value(mc,ind,mu)) + t86*(*CI_users[6]->get_value(mc,ind,mu)) + t74*(*CI_users[7]->get_value(mc,ind,mu)) + t19*(*CI_users[8]->get_value(mc,ind,mu)) + co1*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + t18*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmppm_G_wCI::\
C5g_pmppm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
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
     C5g_pmppm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmppm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:50:51 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = -(spa15*spb15);
complex<T> t5 = square(spa14); 
complex<T> t6 = spa34*spa45; 
complex<T> t12 = spa12*spa15; 
complex<T> t13 = cube(spa25); 
complex<T> t14 = spa23*spa34; 
complex<T> t22 = square(spb14); 
complex<T> t23 = square(spa25); 
complex<T> t24 = square(spb13); 
complex<T> t25 = square(spa24); 
complex<T> t26 = square(spa35); 
complex<T> t28 = s23 - s45; 
complex<T> t29 = s15 - s23; 
complex<T> t30 = s15 - s23 + s45; 
complex<T> t31 = spa15*spa23*spb13 + spa12*spa45*spb14; 
complex<T> t32 = s12 - s45; 
complex<T> t34 = square(spa23); 
complex<T> t35 = square(spa45); 
complex<T> t44 = cube(spb13); 
complex<T> t46 = cube(spa24); 
complex<T> t47 = cube(spa35); 
complex<T> t55 = square(square(spa25)); 
complex<T> t57 = square(spa12); 
complex<T> t58 = square(spa15); 
complex<T> t59 = cube(spb14); 
complex<T> t72 = -(spa35*T(4)); 
complex<T> d19 = spa34*T(3); d19 = T(1)/d19;
complex<T> d27 = spa34*square(spa13); d27 = T(1)/d27;
complex<T> d28 = spa34*square(square(spa13)); d28 = T(1)/d28;
complex<T> t41 = -(t23*T(2)); 
complex<T> t42 = -(spa23*t12); 
complex<T> t63 = t12*T(2); 
complex<T> t64 = spa45*t22; 
complex<T> t65 = spa23*t26; 
complex<T> t66 = spa24*t23; 
complex<T> t74 = t23*t24; 
complex<T> t77 = -(t22*T(4)); 
complex<T> t81 = t12*t6; 
complex<T> t85 = spa25*t31; 
complex<T> t88 = t24*t26; 
complex<T> t92 = spb15*t55; 
complex<T> d1 = s13*t32*t6; d1 = T(1)/d1;
complex<T> d2 = t6*square(spa13)*square(t32); d2 = T(1)/d2;
complex<T> d3 = s13*t32*t6*square(spa13); d3 = T(1)/d3;
complex<T> d4 = spa13*t6*cube(t32)*T(3); d4 = T(1)/d4;
complex<T> d5 = t14*t29*(s45 + t29); d5 = T(1)/d5;
complex<T> d6 = t14*t5*square(t29); d6 = T(1)/d6;
complex<T> d7 = t14*t29*(s45 + t29)*t5; d7 = T(1)/d7;
complex<T> d8 = spa14*t14*cube(t29)*T(3); d8 = T(1)/d8;
complex<T> d10 = s13*t28*t6; d10 = T(1)/d10;
complex<T> d11 = t6*square(spa13)*square(t28); d11 = T(1)/d11;
complex<T> d12 = s13*t28*t6*square(spa13); d12 = T(1)/d12;
complex<T> d13 = t14*t28*(s45 + t29); d13 = T(1)/d13;
complex<T> d14 = t14*t5*square(t28); d14 = T(1)/d14;
complex<T> d15 = t14*t28*(s45 + t29)*t5; d15 = T(1)/d15;
complex<T> d17 = spa13*t6*T(3); d17 = T(1)/d17;
complex<T> d18 = spa14*t14*T(3); d18 = T(1)/d18;
complex<T> d20 = cube(t28); d20 = T(1)/d20;
complex<T> d21 = t6*square(spa13); d21 = T(1)/d21;
complex<T> d22 = t6*square(square(spa13)); d22 = T(1)/d22;
complex<T> d23 = t14*square(s45 + t29); d23 = T(1)/d23;
complex<T> d24 = t14*t5*square(s45 + t29); d24 = T(1)/d24;
complex<T> d25 = spa34*square(s45 + t29); d25 = T(1)/d25;
complex<T> d26 = spa34*t5*square(s45 + t29); d26 = T(1)/d26;
complex<T> d29 = spa15*t6*T(2); d29 = T(1)/d29;
complex<T> d32 = spa12*t14*T(2); d32 = T(1)/d32;
complex<T> d33 = spa23*t6*T(2); d33 = T(1)/d33;
complex<T> t37 = d21*s12; 
complex<T> t38 = d23*s15; 
complex<T> t39 = d22*s23; 
complex<T> t40 = d24*s45; 
complex<T> t43 = -t64; 
complex<T> t60 = d17*spa35; 
complex<T> t61 = d18*t35; 
complex<T> t62 = d4*t44; 
complex<T> t71 = t25*t64; 
complex<T> t78 = t24*t65; 
complex<T> t80 = d8*t46; 
complex<T> t99 = t22*t66; 
complex<T> d9 = spa23*t81*T(6); d9 = T(1)/d9;
complex<T> d16 = spa23*t28*t81*T(6); d16 = T(1)/d16;
complex<T> d30 = spa45*t63; d30 = T(1)/d30;
complex<T> d31 = spa23*t63; d31 = T(1)/d31;
complex<T> t18 = t23*t37*t72 + d22*s12*spa12*spa15*t65*T(2); 
complex<T> t20 = d6*spa12*spa15*t25*t43 - d7*spa12*spa15*t71*T(2) + t58*t59*t80*T(2) + d5*t99*T(4); 
complex<T> t21 = s45*spa24*t22*t38*t41 + s15*spa12*spa15*t40*t71 + d32*spb45*t92; 
complex<T> t36 = d9*T(11); 
complex<T> t54 = d1*t72*t74 + d3*t63*t78 + d2*t42*t88 + t47*t57*t62*T(2); 
complex<T> t69 = s23*spa35*t37*t41 + s12*t12*t39*t65; 
complex<T> t70 = d28*spb45*t63*t65 + t40*t63*t71 + d27*spb45*t23*t72 + d23*s45*t66*t77; 
complex<T> t76 = d24*s15*t63*t71 + t38*t66*t77; 
complex<T> t82 = t39*t63*t65 + d26*spb23*t63*t71 + d21*s23*t23*t72 + d25*spb23*t66*t77; 
complex<T> t86 = t34*t60; 
complex<T> t89 = t57*t61; 
complex<T> t1 = t36*t55 + d14*t12*t71 + d6*t12*t71 + d7*t63*t71 + d10*t72*t74 + d5*t66*t77 + d12*t63*t78 - d19*d20*spb13*spb14*t85 + d11*t42*t88 - d15*t12*t71*T(2) - t58*t59*t80*T(2) + d20*t44*t58*t86*T(2) + d20*spa24*t59*t89*T(2) + d13*t99*T(4) - d16*t13*t31*T(11); 
complex<T> t2 = d14*t12*t25*t43 + t36*t55 + d15*t63*t71 + d13*t66*t77 + d11*t12*t78 + d2*t12*t78 + d19*d20*spb13*spb14*t85 - t47*t57*t62*T(2) - d12*t12*t78*T(2) - d3*t12*t78*T(2) - d20*t44*t58*t86*T(2) - d20*spa24*t59*t89*T(2) + d1*spa35*t74*T(4) + d10*spa35*t74*T(4) + d16*t13*t31*T(11); 
complex<T> co1 = d29*spb12*spb23*t55; 
complex<T> co2 = d30*spb23*spb34*t55; 
complex<T> co3 = d31*spb34*spb45*t55; 
complex<T> co4 = d33*spb12*t92; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t54*(*CI_users[0]->get_value(mc,ind,mu)) + t20*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t18*(*CI_users[4]->get_value(mc,ind,mu)) + t76*(*CI_users[5]->get_value(mc,ind,mu)) + t82*(*CI_users[6]->get_value(mc,ind,mu)) + t70*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + t69*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t21*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_ppmmp_G_wCI::\
C5g_ppmmp_G_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_ppmmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  ppmmp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:50:52 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa14 = SPA(1,4);
complex<T> spa13 = SPA(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = S(3,4);
complex<T> t1 = spa12*spa15; 
complex<T> t5 = square(spa34); 
complex<T> t6 = -(spa14*spb12); 
complex<T> t10 = cube(spa34); 
complex<T> t14 = spa13*spb15; 
complex<T> d4 = spa15*spa45*T(2); d4 = T(1)/d4;
complex<T> d7 = spa12*spa23*T(2); d7 = T(1)/d7;
complex<T> d8 = spa23*spa45*T(2); d8 = T(1)/d8;
complex<T> t7 = -t14; 
complex<T> d1 = spa23*spa45*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = (s23 - s45)*spa23*spa45*t1*T(6); d2 = T(1)/d2;
complex<T> d3 = t1*cube(s23 - s45)*T(3); d3 = T(1)/d3;
complex<T> d5 = spa45*t1*T(2); d5 = T(1)/d5;
complex<T> d6 = spa23*t1*T(2); d6 = T(1)/d6;
complex<T> t11 = spa23*t6 + spa45*t7; 
complex<T> t13 = d1*T(11); 
complex<T> t17 = d3*t11; 
complex<T> t3 = t10*t13 + t14*t17*t6 - d2*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa14*spb12*t14*t17 + d2*t11*t5*T(11); 
complex<T> co1 = d4*spb12*spb23*t10; 
complex<T> co2 = -(d5*s34*spb23*t10); 
complex<T> co3 = -(d6*s34*spb45*t10); 
complex<T> co4 = d7*spb15*spb45*t10; 
complex<T> co5 = d8*spb12*spb15*t10; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_ppmpm_G_wCI::\
C5g_ppmpm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_ppmpm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  ppmpm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:52:06 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> t5 = square(spa14); 
complex<T> t6 = spa12*spa15; 
complex<T> t12 = spa34*spa45; 
complex<T> t13 = cube(spa35); 
complex<T> t14 = spa12*spa23; 
complex<T> t23 = square(spb14); 
complex<T> t24 = square(spa35); 
complex<T> t25 = square(spb24); 
complex<T> t26 = square(spa13); 
complex<T> t27 = square(spa25); 
complex<T> t30 = s15 - s34; 
complex<T> t32 = s15 - s23 + s45; 
complex<T> t33 = -(spa15*spa34*spb14) - spa23*spa45*spb24; 
complex<T> t34 = square(spa15); 
complex<T> t35 = square(spa23); 
complex<T> t46 = cube(spb24); 
complex<T> t55 = square(square(spa35)); 
complex<T> t57 = square(spa34); 
complex<T> t58 = square(spa45); 
complex<T> t59 = cube(spb14); 
complex<T> t60 = cube(spa25); 
complex<T> t70 = cube(spa13); 
complex<T> t79 = -(spa25*T(4)); 
complex<T> t93 = spa35*spb14; 
complex<T> d15 = spa12*T(3); d15 = T(1)/d15;
complex<T> d16 = cube(s15 - s23); d16 = T(1)/d16;
complex<T> d23 = spa12*square(spa24); d23 = T(1)/d23;
complex<T> d24 = spa12*square(square(spa24)); d24 = T(1)/d24;
complex<T> t4 = square(t32); 
complex<T> t22 = t32*t5; 
complex<T> t28 = -(t12*T(2)); 
complex<T> t64 = t12*T(2); 
complex<T> t65 = spa15*t23; 
complex<T> t66 = spa23*t27; 
complex<T> t67 = spa13*t24; 
complex<T> t75 = t12*t25; 
complex<T> t90 = t24*t25; 
complex<T> t98 = spb23*t55; 
complex<T> d1 = spa23*t12*t6*T(6); d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*t14*t32; d2 = T(1)/d2;
complex<T> d3 = t14*t5*square(s15 - s23); d3 = T(1)/d3;
complex<T> d4 = (s15 - s23)*t14*t32*t5; d4 = T(1)/d4;
complex<T> d5 = (s15 - s23)*s24*t6; d5 = T(1)/d5;
complex<T> d6 = s24*t30*t6; d6 = T(1)/d6;
complex<T> d7 = t6*square(s15 - s23)*square(spa24); d7 = T(1)/d7;
complex<T> d8 = (s15 - s23)*s24*t6*square(spa24); d8 = T(1)/d8;
complex<T> d9 = t6*square(spa24)*square(t30); d9 = T(1)/d9;
complex<T> d10 = s24*t30*t6*square(spa24); d10 = T(1)/d10;
complex<T> d11 = spa24*t6*cube(t30)*T(3); d11 = T(1)/d11;
complex<T> d12 = (s15 - s23)*spa23*t12*t6*T(6); d12 = T(1)/d12;
complex<T> d13 = spa14*t14*T(3); d13 = T(1)/d13;
complex<T> d14 = spa24*t6*T(3); d14 = T(1)/d14;
complex<T> d17 = (s23 - s45)*t14*t32; d17 = T(1)/d17;
complex<T> d18 = t14*t5*square(s23 - s45); d18 = T(1)/d18;
complex<T> d19 = (s23 - s45)*t14*t32*t5; d19 = T(1)/d19;
complex<T> d20 = spa14*t14*cube(s23 - s45)*T(3); d20 = T(1)/d20;
complex<T> d25 = t6*square(spa24); d25 = T(1)/d25;
complex<T> d26 = t6*square(square(spa24)); d26 = T(1)/d26;
complex<T> d30 = spa45*t6*T(2); d30 = T(1)/d30;
complex<T> d31 = spa23*t6*T(2); d31 = T(1)/d31;
complex<T> d32 = spa34*t14*T(2); d32 = T(1)/d32;
complex<T> t37 = d1*T(11); 
complex<T> t39 = d25*s23; 
complex<T> t40 = d26*s34; 
complex<T> t43 = -t65; 
complex<T> t44 = -t66; 
complex<T> t61 = d13*t34; 
complex<T> t62 = d14*t35; 
complex<T> t63 = d11*t46; 
complex<T> t72 = d20*t59; 
complex<T> t73 = t26*t65; 
complex<T> t80 = t23*t67; 
complex<T> t86 = t25*t66; 
complex<T> d21 = t14*t4; d21 = T(1)/d21;
complex<T> d22 = t14*t4*t5; d22 = T(1)/d22;
complex<T> d27 = spa12*t4; d27 = T(1)/d27;
complex<T> d28 = spa12*t4*t5; d28 = T(1)/d28;
complex<T> d29 = spa15*t64; d29 = T(1)/d29;
complex<T> d33 = spa23*t64; d33 = T(1)/d33;
complex<T> t2 = d3*t12*t26*t43 + t37*t55 + d4*t28*t73 + d7*t66*t75 + d9*t66*t75 + d10*t64*t86 + d8*t64*t86 + d5*t79*t90 + d6*t79*t90 + d15*d16*spb24*t33*t93 + d16*spa13*t57*t59*t61*T(2) + d16*spa25*t46*t58*t62*T(2) + t57*t60*t63*T(2) + d2*t80*T(4) + d12*t13*t33*T(11); 
complex<T> t18 = d25*s34*t24*t79 + spa34*spa45*t40*t66*T(2); 
complex<T> t21 = s23*spa34*spa45*t40*t66 + d30*spb34*t98 - s34*spa25*t24*t39*T(2); 
complex<T> t38 = d21*s15; 
complex<T> t41 = d22*s45; 
complex<T> t53 = d9*t44*t75 + d10*t28*t86 - t57*t60*t63*T(2) + d6*spa25*t90*T(4); 
complex<T> t85 = -(t80*T(4)); 
complex<T> t1 = t37*t55 + d18*t12*t73 + d3*t12*t73 + d19*t28*t73 + d4*t64*t73 + d7*t44*t75 + d2*t85 + d8*t28*t86 - d15*d16*spb24*t33*t93 - d16*spa13*t57*t59*t61*T(2) - d16*spa25*t46*t58*t62*T(2) + t58*t70*t72*T(2) + d17*t80*T(4) + d5*spa25*t90*T(4) - d12*t13*t33*T(11); 
complex<T> t19 = s15*spa34*spa45*t41*t73 - s45*t38*t80*T(2); 
complex<T> t52 = d18*t12*t26*t43 + d19*t64*t73 + d17*t85 - t58*t70*t72*T(2); 
complex<T> t78 = d24*spb15*t64*t66 + d22*s15*t64*t73 + d23*spb15*t24*t79 + t38*t85; 
complex<T> t84 = t41*t64*t73 + d21*s45*t85; 
complex<T> t88 = d26*s23*t64*t66 + d28*spb23*t64*t73 + t24*t39*t79 + d27*spb23*t85; 
complex<T> co1 = d29*spb12*t98; 
complex<T> co2 = d31*spb34*spb45*t55; 
complex<T> co3 = d32*spb15*spb45*t55; 
complex<T> co4 = d33*spb12*spb15*t55; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t53*(*CI_users[2]->get_value(mc,ind,mu)) + t52*(*CI_users[3]->get_value(mc,ind,mu)) + t78*(*CI_users[4]->get_value(mc,ind,mu)) + t88*(*CI_users[5]->get_value(mc,ind,mu)) + t18*(*CI_users[6]->get_value(mc,ind,mu)) + t84*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t19*(*CI_users[9]->get_value(mc,ind,mu)) + t21*(*CI_users[10]->get_value(mc,ind,mu)) + co2*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pppmm_G_wCI::\
C5g_pppmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pppmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  pppmm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:52:08 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> s45 = S(4,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> t1 = spa12*spa23; 
complex<T> t5 = square(spa45); 
complex<T> t6 = -(spa24*spb12); 
complex<T> t10 = cube(spa45); 
complex<T> t14 = spa25*spb23; 
complex<T> d4 = spa15*spa34*T(2); d4 = T(1)/d4;
complex<T> d5 = spa12*spa15*T(2); d5 = T(1)/d5;
complex<T> d8 = spa23*spa34*T(2); d8 = T(1)/d8;
complex<T> t7 = -t14; 
complex<T> d1 = spa15*spa34*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = (s15 - s34)*spa15*spa34*t1*T(6); d2 = T(1)/d2;
complex<T> d3 = t1*cube(s15 - s34)*T(3); d3 = T(1)/d3;
complex<T> d6 = spa15*t1*T(2); d6 = T(1)/d6;
complex<T> d7 = spa34*t1*T(2); d7 = T(1)/d7;
complex<T> t11 = spa15*t6 + spa34*t7; 
complex<T> t13 = d1*T(11); 
complex<T> t17 = d3*t11; 
complex<T> t3 = t10*t13 + t14*t17*t6 - d2*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa24*spb12*t14*t17 + d2*t11*t5*T(11); 
complex<T> co1 = d4*spb12*spb23*t10; 
complex<T> co2 = d5*spb23*spb34*t10; 
complex<T> co3 = -(d6*s45*spb34*t10); 
complex<T> co4 = -(d7*s45*spb15*t10); 
complex<T> co5 = d8*spb12*spb15*t10; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mmppp_nf_wCI::\
C5g_mmppp_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C5g_mmppp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mmppp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:52:09 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s23 = S(2,3);
complex<T> t5 = -(spa14*spb34); 
complex<T> t7 = spa34*T(3); 
complex<T> t11 = cube(spa12); 
complex<T> t12 = spa24*spb45; 
complex<T> t6 = -t12; 
complex<T> d1 = spa15*spa23*spa45*t7; d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*spa15*spa23*spa45*t7; d2 = T(1)/d2;
complex<T> d3 = spa45*t7*cube(s15 - s23); d3 = T(1)/d3;
complex<T> t9 = spa23*t5 + spa15*t6; 
complex<T> t17 = d3*t9; 
complex<T> t2 = -(d1*t11) + t12*t17*t5 - d2*t9*square(spa12); 
complex<T> t3 = -(d1*t11) + spa14*spb34*t12*t17 + d2*t9*square(spa12); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mpmpp_nf_wCI::\
C5g_mpmpp_nf_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mpmpp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mpmpp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:53:22 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> s23 = S(2,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = S(1,2);
complex<T> t5 = spa34*T(2); 
complex<T> t6 = square(spa25); 
complex<T> t8 = spa15*spa45; 
complex<T> t12 = spa12*spa23; 
complex<T> t14 = spa34*spa45; 
complex<T> t23 = square(spb25); 
complex<T> t24 = square(spa13); 
complex<T> t25 = square(spb24); 
complex<T> t26 = square(spa14); 
complex<T> t27 = square(spa35); 
complex<T> t30 = s15 - s23; 
complex<T> t32 = s12 + s15 - s34; 
complex<T> t33 = -(spa34*T(2)); 
complex<T> t34 = -(spa12*spa34*spb24) - spa15*spa23*spb25; 
complex<T> t39 = -(spa15*T(2)); 
complex<T> t42 = cube(spb24); 
complex<T> t44 = square(spa34); 
complex<T> t53 = square(spa12); 
complex<T> t54 = square(spa23); 
complex<T> t55 = cube(spb25); 
complex<T> t56 = square(square(spa13)); 
complex<T> t58 = square(spa15); 
complex<T> t59 = cube(spa35); 
complex<T> t65 = spa15*T(2); 
complex<T> t68 = cube(spa14); 
complex<T> t83 = spa34*T(3); 
complex<T> d19 = spa45*T(3); d19 = T(1)/d19;
complex<T> d20 = cube(s15 - s34); d20 = T(1)/d20;
complex<T> d23 = spa45*square(spa24); d23 = T(1)/d23;
complex<T> d24 = spa45*square(square(spa24)); d24 = T(1)/d24;
complex<T> t22 = t32*t6; 
complex<T> t38 = t23*t27; 
complex<T> t40 = -(t12*t26); 
complex<T> t64 = t12*t25; 
complex<T> t67 = d19*t34; 
complex<T> t71 = spa35*t23; 
complex<T> t72 = spa14*t24; 
complex<T> t73 = spa15*t12; 
complex<T> t86 = t26*t33; 
complex<T> t96 = t42*t53; 
complex<T> t104 = t55*t59; 
complex<T> d1 = (s12 - s34)*t14*t32; d1 = T(1)/d1;
complex<T> d2 = t14*t6*square(s12 - s34); d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*t14*t32*t6; d3 = T(1)/d3;
complex<T> d4 = spa25*t14*cube(s12 - s34)*T(3); d4 = T(1)/d4;
complex<T> d5 = t12*t8*t83; d5 = T(1)/d5;
complex<T> d6 = s24*t30*t8; d6 = T(1)/d6;
complex<T> d7 = s24*(s15 - s34)*t8; d7 = T(1)/d7;
complex<T> d8 = t8*square(spa24)*square(t30); d8 = T(1)/d8;
complex<T> d9 = s24*t30*t8*square(spa24); d9 = T(1)/d9;
complex<T> d10 = t8*square(s15 - s34)*square(spa24); d10 = T(1)/d10;
complex<T> d11 = s24*(s15 - s34)*t8*square(spa24); d11 = T(1)/d11;
complex<T> d12 = spa24*t8*cube(t30)*T(3); d12 = T(1)/d12;
complex<T> d13 = (s15 - s34)*t14*t32; d13 = T(1)/d13;
complex<T> d14 = t14*t6*square(s15 - s34); d14 = T(1)/d14;
complex<T> d15 = (s15 - s34)*t14*t32*t6; d15 = T(1)/d15;
complex<T> d16 = (s15 - s34)*t12*t8*t83; d16 = T(1)/d16;
complex<T> d17 = spa24*t8*T(3); d17 = T(1)/d17;
complex<T> d18 = spa25*t14*T(3); d18 = T(1)/d18;
complex<T> d21 = t14*square(t32); d21 = T(1)/d21;
complex<T> d22 = t14*t6*square(t32); d22 = T(1)/d22;
complex<T> d25 = t8*square(spa24); d25 = T(1)/d25;
complex<T> d26 = t8*square(square(spa24)); d26 = T(1)/d26;
complex<T> d27 = spa45*square(t32); d27 = T(1)/d27;
complex<T> d28 = spa45*t6*square(t32); d28 = T(1)/d28;
complex<T> d29 = t8*square(spa24)*T(2); d29 = T(1)/d29;
complex<T> d30 = spa45*t5*square(t32); d30 = T(1)/d30;
complex<T> t29 = -t73; 
complex<T> t37 = -(d26*s23); 
complex<T> t45 = d22*s12; 
complex<T> t60 = d12*t42; 
complex<T> t61 = d17*t44; 
complex<T> t62 = -t71; 
complex<T> t80 = t26*t64; 
complex<T> t81 = t24*t71; 
complex<T> t84 = d18*t58; 
complex<T> t85 = t12*t38; 
complex<T> t101 = spa12*t86; 
complex<T> t18 = s23*(d26*spa23*t101 + d25*t72); 
complex<T> t20 = d24*spa23*spb15*t101 + d22*s15*spa12*spa23*t38*t39 + d23*spb15*t72 + d21*s15*t81; 
complex<T> t21 = d26*s34*spa23*t101 + d28*spa12*spa23*spb34*t38*t39 + d25*s34*t72 + d27*spb34*t81; 
complex<T> t50 = d1*t24*t62 + d2*t38*t73 + d3*t65*t85 - d4*t104*t53*T(2); 
complex<T> t52 = -(d6*t25*t72) + d8*spa34*t80 + d9*t5*t80 + t54*t60*t68*T(2); 
complex<T> t70 = s15*(t29*t38*t45 + d30*s12*t81); 
complex<T> t78 = d21*s12*t81 + t39*t45*t85; 
complex<T> t79 = s34*(spa34*t12*t26*t37 + d29*s23*t72); 
complex<T> t97 = spa35*t84; 
complex<T> t1 = d10*spa34*t25*t40 + d8*spa34*t25*t40 - d5*t56 + d13*t24*t62 - d20*spa13*spb24*spb25*t67 + d6*t25*t72 + d7*t25*t72 + d14*t38*t73 + d11*t33*t80 + d9*t33*t80 + d15*t65*t85 - d16*t34*cube(spa13) - t54*t60*t68*T(2) - d20*spa14*t61*t96*T(2) - d20*t54*t55*t97*T(2); 
complex<T> t2 = d14*t29*t38 + d2*t29*t38 - d5*t56 + d20*spa13*spb24*spb25*t67 - d7*t25*t72 + d10*spa34*t80 + d11*t5*t80 + d1*t81 + d13*t81 + d15*t39*t85 + d3*t39*t85 + d16*t34*cube(spa13) + d4*t104*t53*T(2) + d20*spa14*t61*t96*T(2) + d20*t54*t55*t97*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t50*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t52*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t78*(*CI_users[4]->get_value(mc,ind,mu)) + t20*(*CI_users[5]->get_value(mc,ind,mu)) + t18*(*CI_users[6]->get_value(mc,ind,mu)) + t21*(*CI_users[7]->get_value(mc,ind,mu)) + t79*(*CI_users[8]->get_value(mc,ind,mu)) + t70*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mppmp_nf_wCI::\
C5g_mppmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mppmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mppmp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:54:32 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa15 = SPA(1,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> s45 = S(4,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = -(spa12*spb12);
complex<T> s15 = S(1,5);
complex<T> s35 = -(spa35*spb35);
complex<T> t5 = spa34*T(2); 
complex<T> t6 = square(spa25); 
complex<T> t8 = spa12*spa23; 
complex<T> t12 = spa15*spa45; 
complex<T> t14 = spa23*spa34; 
complex<T> t22 = square(spb25); 
complex<T> t23 = square(spa14); 
complex<T> t24 = square(spb35); 
complex<T> t25 = square(spa13); 
complex<T> t26 = square(spa24); 
complex<T> t27 = s12 - s34; 
complex<T> t29 = s15 - s34; 
complex<T> t30 = s12 + s15 - s34; 
complex<T> t31 = -(spa34*T(2)); 
complex<T> t32 = spa12*spa45*spb25 + spa15*spa34*spb35; 
complex<T> t33 = s12 - s45; 
complex<T> t38 = -(spa12*T(2)); 
complex<T> t41 = cube(spb35); 
complex<T> t53 = square(spa15); 
complex<T> t54 = square(spa45); 
complex<T> t55 = cube(spb25); 
complex<T> t56 = square(spa12); 
complex<T> t58 = square(square(spa14)); 
complex<T> t59 = cube(spa24); 
complex<T> t60 = square(spa34); 
complex<T> t65 = spa12*T(2); 
complex<T> t70 = cube(spa13); 
complex<T> t77 = spa34*T(3); 
complex<T> d15 = spa23*T(3); d15 = T(1)/d15;
complex<T> d21 = spa23*square(spa35); d21 = T(1)/d21;
complex<T> d22 = spa23*square(square(spa35)); d22 = T(1)/d22;
complex<T> t37 = t22*t26; 
complex<T> t39 = -(t12*t25); 
complex<T> t64 = t12*t24; 
complex<T> t69 = d15*t32; 
complex<T> t73 = spa24*t22; 
complex<T> t74 = spa13*t23; 
complex<T> t75 = spa12*t12; 
complex<T> t87 = t25*t31; 
complex<T> t95 = t53*t59; 
complex<T> t103 = spa24*t55; 
complex<T> d1 = t12*t77*t8; d1 = T(1)/d1;
complex<T> d2 = t14*t27*(s15 + t27); d2 = T(1)/d2;
complex<T> d3 = t14*t6*square(t27); d3 = T(1)/d3;
complex<T> d4 = t14*t27*(s15 + t27)*t6; d4 = T(1)/d4;
complex<T> d5 = s35*t27*t8; d5 = T(1)/d5;
complex<T> d6 = s35*t33*t8; d6 = T(1)/d6;
complex<T> d7 = t8*square(spa35)*square(t27); d7 = T(1)/d7;
complex<T> d8 = s35*t27*t8*square(spa35); d8 = T(1)/d8;
complex<T> d9 = t8*square(spa35)*square(t33); d9 = T(1)/d9;
complex<T> d10 = s35*t33*t8*square(spa35); d10 = T(1)/d10;
complex<T> d11 = spa35*t8*cube(t33)*T(3); d11 = T(1)/d11;
complex<T> d12 = t12*t27*t77*t8; d12 = T(1)/d12;
complex<T> d13 = spa25*t14*T(3); d13 = T(1)/d13;
complex<T> d14 = spa35*t8*T(3); d14 = T(1)/d14;
complex<T> d16 = cube(t27); d16 = T(1)/d16;
complex<T> d17 = t14*(s15 + t27)*t29; d17 = T(1)/d17;
complex<T> d18 = t14*t6*square(t29); d18 = T(1)/d18;
complex<T> d19 = t14*(s15 + t27)*t29*t6; d19 = T(1)/d19;
complex<T> d20 = spa25*t14*cube(t29)*T(3); d20 = T(1)/d20;
complex<T> d23 = t14*square(s15 + t27); d23 = T(1)/d23;
complex<T> d24 = t14*t6*square(s15 + t27); d24 = T(1)/d24;
complex<T> d25 = t8*square(spa35); d25 = T(1)/d25;
complex<T> d26 = t8*square(square(spa35)); d26 = T(1)/d26;
complex<T> d27 = spa23*square(s15 + t27); d27 = T(1)/d27;
complex<T> d28 = spa23*t6*square(s15 + t27); d28 = T(1)/d28;
complex<T> d29 = spa23*t5*square(s15 + t27); d29 = T(1)/d29;
complex<T> d30 = t8*square(spa35)*T(2); d30 = T(1)/d30;
complex<T> t28 = -t75; 
complex<T> t36 = -(d26*s34); 
complex<T> t44 = d24*s12; 
complex<T> t61 = d11*t41; 
complex<T> t62 = -t73; 
complex<T> t81 = t25*t64; 
complex<T> t82 = t23*t73; 
complex<T> t85 = d13*t56; 
complex<T> t86 = t12*t37; 
complex<T> t92 = d14*t41; 
complex<T> t100 = spa15*t87; 
complex<T> t18 = s45*(d26*spa45*t100 + d25*t74); 
complex<T> t19 = d22*spa45*spb12*t100 + spa15*spa45*t37*t38*t44 + d21*spb12*t74 + d23*s12*t82; 
complex<T> t20 = d26*s34*spa45*t100 + d28*spa15*spa45*spb34*t37*t38 + d25*s34*t74 + d27*spb34*t82; 
complex<T> t50 = d17*t23*t62 + d18*t37*t75 + d19*t65*t86 - d20*t55*t95*T(2); 
complex<T> t52 = -(d6*t24*t74) + d9*spa34*t81 + d10*t5*t81 + t54*t61*t70*T(2); 
complex<T> t72 = s15*(d23*t82 + d24*t38*t86); 
complex<T> t79 = s15*(t28*t37*t44 + d29*s12*t82); 
complex<T> t80 = s45*(spa34*t12*t25*t36 + d30*s34*t74); 
complex<T> t106 = t60*t92; 
complex<T> t1 = d18*t28*t37 + d3*t28*t37 - d1*t58 - d16*spa14*spb25*spb35*t69 - d5*t24*t74 + d7*spa34*t81 + d8*t5*t81 + d17*t82 + d2*t82 + d19*t38*t86 + d4*t38*t86 - d12*t32*cube(spa14) + d16*spa13*t106*t53*T(2) + d16*t103*t54*t85*T(2) + d20*t55*t95*T(2); 
complex<T> t2 = d7*spa34*t24*t39 + d9*spa34*t24*t39 - d1*t58 + d2*t23*t62 + d16*spa14*spb25*spb35*t69 + d5*t24*t74 + d6*t24*t74 + d3*t37*t75 + d10*t31*t81 + d8*t31*t81 + d4*t65*t86 + d12*t32*cube(spa14) - d16*spa13*t106*t53*T(2) - t54*t61*t70*T(2) - d16*t103*t54*t85*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t50*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t52*(*CI_users[3]->get_value(mc,ind,mu)) + t19*(*CI_users[4]->get_value(mc,ind,mu)) + t72*(*CI_users[5]->get_value(mc,ind,mu)) + t20*(*CI_users[6]->get_value(mc,ind,mu)) + t18*(*CI_users[7]->get_value(mc,ind,mu)) + t79*(*CI_users[8]->get_value(mc,ind,mu)) + t80*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mpppm_nf_wCI::\
C5g_mpppm_nf_wCI
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
     C5g_mpppm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mpppm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:54:34 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> t7 = spa23*T(3); 
complex<T> t11 = cube(spa15); 
complex<T> t12 = spa13*spb34; 
complex<T> t9 = spa12*spa35*spb23 + spa45*t12; 
complex<T> d1 = spa12*spa34*spa45*t7; d1 = T(1)/d1;
complex<T> d2 = (s12 - s45)*spa12*spa34*spa45*t7; d2 = T(1)/d2;
complex<T> d3 = spa34*t7*cube(s12 - s45); d3 = T(1)/d3;
complex<T> t18 = d3*t12; 
complex<T> t2 = -(d1*t11) - spa35*spb23*t18*t9 - d2*t9*square(spa15); 
complex<T> t3 = -(d1*t11) + spa35*spb23*t18*t9 + d2*t9*square(spa15); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmmpp_nf_wCI::\
C5g_pmmpp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmmpp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmmpp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:54:35 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb15 = SPB(1,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t5 = -(spa35*spb15); 
complex<T> t7 = spa15*T(3); 
complex<T> t11 = cube(spa23); 
complex<T> t12 = spa25*spb45; 
complex<T> t6 = -t12; 
complex<T> d1 = spa12*spa34*spa45*t7; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spa12*spa34*spa45*t7; d2 = T(1)/d2;
complex<T> d3 = spa45*t7*cube(s12 - s34); d3 = T(1)/d3;
complex<T> t9 = spa12*t5 + spa34*t6; 
complex<T> t17 = d3*t9; 
complex<T> t2 = -(d1*t11) + t12*t17*t5 - d2*t9*square(spa23); 
complex<T> t3 = -(d1*t11) + spa35*spb15*t12*t17 + d2*t9*square(spa23); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmpmp_nf_wCI::\
C5g_pmpmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmpmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmpmp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:55:26 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = S(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s35 = -(spa35*spb35);
complex<T> t6 = spa12*spa15; 
complex<T> t9 = spa45*T(3); 
complex<T> t11 = spa23*spa34; 
complex<T> t22 = square(spa24); 
complex<T> t23 = square(spb13); 
complex<T> t24 = square(spb35); 
complex<T> t25 = square(spa14); 
complex<T> t26 = square(spa25); 
complex<T> t27 = s12 - s45; 
complex<T> t28 = -(spa23*T(2)); 
complex<T> t29 = spa12*spa34; 
complex<T> t30 = s12 - s34; 
complex<T> t31 = s23 - s45; 
complex<T> t32 = -(spa12*spa34*spb13) - spa23*spa45*spb35; 
complex<T> t35 = square(spa45); 
complex<T> t40 = cube(spb13); 
complex<T> t41 = cube(spb35); 
complex<T> t54 = square(spa23); 
complex<T> t55 = square(spa34); 
complex<T> t56 = square(spa12); 
complex<T> t58 = square(square(spa24)); 
complex<T> t59 = cube(spa25); 
complex<T> t70 = cube(spa14); 
complex<T> t88 = spa24*spb13; 
complex<T> d15 = spa15*T(3); d15 = T(1)/d15;
complex<T> d21 = spa15*spa45*square(spa13); d21 = T(1)/d21;
complex<T> d22 = spa15*spa45*square(square(spa13)); d22 = T(1)/d22;
complex<T> d23 = spa15*square(spa35); d23 = T(1)/d23;
complex<T> d24 = spa15*square(square(spa35)); d24 = T(1)/d24;
complex<T> d27 = spa15*square(spa13); d27 = T(1)/d27;
complex<T> d28 = spa15*square(square(spa13)); d28 = T(1)/d28;
complex<T> d29 = spa15*spa45*square(spa13)*T(2); d29 = T(1)/d29;
complex<T> t36 = -(d22*s12); 
complex<T> t63 = spa45*t26; 
complex<T> t64 = spa12*t11; 
complex<T> t65 = t23*t25; 
complex<T> t66 = -(t22*t24); 
complex<T> t74 = t11*t24; 
complex<T> t75 = spa14*t22; 
complex<T> t83 = spa25*t22; 
complex<T> t90 = spa34*t28; 
complex<T> t94 = t28*t29; 
complex<T> d1 = t11*t6*t9; d1 = T(1)/d1;
complex<T> d2 = s13*spa15*spa45*t27; d2 = T(1)/d2;
complex<T> d3 = spa15*spa45*square(spa13)*square(t27); d3 = T(1)/d3;
complex<T> d4 = s13*spa15*spa45*t27*square(spa13); d4 = T(1)/d4;
complex<T> d5 = s35*t30*t6; d5 = T(1)/d5;
complex<T> d6 = s35*t27*t6; d6 = T(1)/d6;
complex<T> d7 = t6*square(spa35)*square(t30); d7 = T(1)/d7;
complex<T> d8 = s35*t30*t6*square(spa35); d8 = T(1)/d8;
complex<T> d9 = t6*square(spa35)*square(t27); d9 = T(1)/d9;
complex<T> d10 = s35*t27*t6*square(spa35); d10 = T(1)/d10;
complex<T> d11 = spa35*t6*cube(t30)*T(3); d11 = T(1)/d11;
complex<T> d12 = t11*t27*t6*t9; d12 = T(1)/d12;
complex<T> d13 = spa13*spa15*t9; d13 = T(1)/d13;
complex<T> d14 = spa35*t6*T(3); d14 = T(1)/d14;
complex<T> d16 = cube(t27); d16 = T(1)/d16;
complex<T> d17 = s13*spa15*spa45*t31; d17 = T(1)/d17;
complex<T> d18 = spa15*spa45*square(spa13)*square(t31); d18 = T(1)/d18;
complex<T> d19 = s13*spa15*spa45*t31*square(spa13); d19 = T(1)/d19;
complex<T> d20 = spa13*spa15*t9*cube(t31); d20 = T(1)/d20;
complex<T> d25 = t6*square(spa35); d25 = T(1)/d25;
complex<T> d26 = t6*square(square(spa35)); d26 = T(1)/d26;
complex<T> d30 = t6*square(spa35)*T(2); d30 = T(1)/d30;
complex<T> t17 = s23*(d21*t75 + d22*t25*t94); 
complex<T> t19 = d21*s12*t75 + d23*spb12*t83 + d24*spb12*t63*t90 + d22*s12*t25*t94; 
complex<T> t37 = -(d26*s34); 
complex<T> t38 = -t63; 
complex<T> t39 = -t64; 
complex<T> t60 = d14*t35; 
complex<T> t61 = d13*t40; 
complex<T> t62 = d11*t41; 
complex<T> t72 = d20*t40; 
complex<T> t77 = d25*spa25; 
complex<T> t81 = s23*(t25*t36*t64 + d29*s12*t75); 
complex<T> t18 = s34*(t22*t77 + d26*t63*t90); 
complex<T> t21 = d27*spb45*t75 + s45*t22*t77 + d26*s45*t63*t90 + d28*spb45*t25*t94; 
complex<T> t51 = d18*t64*t65 + d17*t23*t75 - d19*t64*t65*T(2) - t54*t70*t72*T(2); 
complex<T> t52 = d5*spa25*t66 + t55*t59*t62*T(2) + t63*t74*(d7 + d8*T(2)); 
complex<T> t73 = s45*(t11*t37*t63 + d30*s34*t83); 
complex<T> t85 = t41*t60; 
complex<T> t1 = -(d1*t58) + d3*t64*t65 + d7*t38*t74 + d9*t38*t74 + d2*t23*t75 + d5*t24*t83 + d6*t24*t83 - d15*d16*spb35*t32*t88 - d12*t32*cube(spa24) - d16*spa14*t55*t56*t61*T(2) - t55*t59*t62*T(2) - d4*t64*t65*T(2) - d10*t63*t74*T(2) - d8*t63*t74*T(2) - d16*spa25*t54*t85*T(2); 
complex<T> t2 = -(d1*t58) + d18*t39*t65 + d3*t39*t65 + d6*spa25*t66 + d9*t63*t74 - d17*t23*t75 - d2*t23*t75 + d15*d16*spb35*t32*t88 + d12*t32*cube(spa24) + d16*spa14*t55*t56*t61*T(2) + d19*t64*t65*T(2) + d4*t64*t65*T(2) + t54*t70*t72*T(2) + d10*t63*t74*T(2) + d16*spa25*t54*t85*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t51*(*CI_users[1]->get_value(mc,ind,mu)) + t52*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t19*(*CI_users[4]->get_value(mc,ind,mu)) + t17*(*CI_users[5]->get_value(mc,ind,mu)) + t18*(*CI_users[6]->get_value(mc,ind,mu)) + t21*(*CI_users[7]->get_value(mc,ind,mu)) + t81*(*CI_users[8]->get_value(mc,ind,mu)) + t73*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmppm_nf_wCI::\
C5g_pmppm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmppm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmppm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:56:35 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa13 = SPA(1,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb13 = SPB(1,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = S(1,5);
complex<T> t5 = spa23*T(2); 
complex<T> t6 = square(spa14); 
complex<T> t8 = spa34*spa45; 
complex<T> t12 = spa12*spa15; 
complex<T> t14 = spa23*spa34; 
complex<T> t22 = square(spb14); 
complex<T> t23 = square(spa25); 
complex<T> t24 = square(spb13); 
complex<T> t25 = square(spa24); 
complex<T> t26 = square(spa35); 
complex<T> t27 = s23 - s45; 
complex<T> t29 = s15 - s23; 
complex<T> t30 = s15 - s23 + s45; 
complex<T> t31 = -(spa23*T(2)); 
complex<T> t32 = spa15*spa23*spb13 + spa12*spa45*spb14; 
complex<T> t33 = s12 - s45; 
complex<T> t38 = -(spa45*T(2)); 
complex<T> t41 = cube(spb13); 
complex<T> t53 = square(spa12); 
complex<T> t54 = square(spa15); 
complex<T> t55 = cube(spb14); 
complex<T> t56 = square(spa23); 
complex<T> t58 = square(square(spa25)); 
complex<T> t59 = cube(spa35); 
complex<T> t60 = square(spa45); 
complex<T> t65 = spa45*T(2); 
complex<T> t70 = cube(spa24); 
complex<T> d19 = spa34*T(3); d19 = T(1)/d19;
complex<T> d27 = spa34*square(spa13); d27 = T(1)/d27;
complex<T> d28 = spa34*square(square(spa13)); d28 = T(1)/d28;
complex<T> t28 = -(spa45*t12); 
complex<T> t37 = t22*t25; 
complex<T> t39 = -(t12*t26); 
complex<T> t64 = t12*t24; 
complex<T> t69 = d19*t32; 
complex<T> t74 = spa24*t22; 
complex<T> t75 = spa35*t23; 
complex<T> t82 = t26*t31; 
complex<T> t92 = spa35*t56; 
complex<T> t97 = t41*t53; 
complex<T> d1 = s13*t33*t8; d1 = T(1)/d1;
complex<T> d2 = t8*square(spa13)*square(t33); d2 = T(1)/d2;
complex<T> d3 = s13*t33*t8*square(spa13); d3 = T(1)/d3;
complex<T> d4 = spa13*t8*cube(t33)*T(3); d4 = T(1)/d4;
complex<T> d5 = t14*t29*(s45 + t29); d5 = T(1)/d5;
complex<T> d6 = t14*t6*square(t29); d6 = T(1)/d6;
complex<T> d7 = t14*t29*(s45 + t29)*t6; d7 = T(1)/d7;
complex<T> d8 = spa14*t14*cube(t29)*T(3); d8 = T(1)/d8;
complex<T> d9 = spa23*t12*t8*T(3); d9 = T(1)/d9;
complex<T> d10 = s13*t27*t8; d10 = T(1)/d10;
complex<T> d11 = t8*square(spa13)*square(t27); d11 = T(1)/d11;
complex<T> d12 = s13*t27*t8*square(spa13); d12 = T(1)/d12;
complex<T> d13 = t14*t27*(s45 + t29); d13 = T(1)/d13;
complex<T> d14 = t14*t6*square(t27); d14 = T(1)/d14;
complex<T> d15 = t14*t27*(s45 + t29)*t6; d15 = T(1)/d15;
complex<T> d16 = spa23*t12*t27*t8*T(3); d16 = T(1)/d16;
complex<T> d17 = spa13*t8*T(3); d17 = T(1)/d17;
complex<T> d18 = spa14*t14*T(3); d18 = T(1)/d18;
complex<T> d20 = cube(t27); d20 = T(1)/d20;
complex<T> d21 = t8*square(spa13); d21 = T(1)/d21;
complex<T> d22 = t8*square(square(spa13)); d22 = T(1)/d22;
complex<T> d23 = t14*square(s45 + t29); d23 = T(1)/d23;
complex<T> d24 = t14*t6*square(s45 + t29); d24 = T(1)/d24;
complex<T> d25 = spa34*square(s45 + t29); d25 = T(1)/d25;
complex<T> d26 = spa34*t6*square(s45 + t29); d26 = T(1)/d26;
complex<T> d29 = t8*square(spa13)*T(2); d29 = T(1)/d29;
complex<T> d30 = spa34*t5*square(s45 + t29); d30 = T(1)/d30;
complex<T> t35 = -(d22*s12); 
complex<T> t44 = d24*s15; 
complex<T> t61 = d17*t41; 
complex<T> t73 = d8*t55; 
complex<T> t79 = d18*t60; 
complex<T> t81 = t12*t37; 
complex<T> t87 = t23*t74; 
complex<T> t88 = t26*t64; 
complex<T> t100 = spa12*t82; 
complex<T> t18 = s12*(d22*spa15*t100 + d21*t75); 
complex<T> t20 = d22*s23*spa15*t100 + d26*spa12*spa15*spb23*t37*t38 + d21*s23*t75 + d25*spb23*t87; 
complex<T> t21 = d28*spa15*spb45*t100 + d24*s45*spa12*spa15*t37*t38 + d27*spb45*t75 + d23*s45*t87; 
complex<T> t50 = d6*spa45*t81 + d7*t65*t81 - d5*t87 - t54*t70*t73*T(2); 
complex<T> t52 = d1*t24*t75 + d3*t64*t82 + d2*spa23*t88 - d4*t59*t97*T(2); 
complex<T> t72 = s23*(spa23*t12*t26*t35 + d29*s12*t75); 
complex<T> t80 = t38*t44*t81 + d23*s15*t87; 
complex<T> t86 = s45*(t28*t37*t44 + d30*s15*t87); 
complex<T> t106 = t55*t79; 
complex<T> t1 = d11*spa23*t24*t39 + d2*spa23*t24*t39 - d9*t58 - d20*spa25*spb13*spb14*t69 - d1*t24*t75 - d10*t24*t75 + d14*spa45*t81 + d15*t38*t81 + d13*t87 + d12*t5*t88 + d3*t5*t88 - d16*t32*cube(spa25) + d20*spa24*t106*t53*T(2) + d20*t54*t61*t92*T(2) + d4*t59*t97*T(2); 
complex<T> t2 = d14*t28*t37 + d6*t28*t37 - d9*t58 + d20*spa25*spb13*spb14*t69 + d10*t24*t75 + d7*t38*t81 + d15*t65*t81 + d12*t64*t82 - d13*t87 + d5*t87 + d11*spa23*t88 + d16*t32*cube(spa25) - d20*spa24*t106*t53*T(2) + t54*t70*t73*T(2) - d20*t54*t61*t92*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t52*(*CI_users[0]->get_value(mc,ind,mu)) + t50*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t18*(*CI_users[4]->get_value(mc,ind,mu)) + t80*(*CI_users[5]->get_value(mc,ind,mu)) + t20*(*CI_users[6]->get_value(mc,ind,mu)) + t21*(*CI_users[7]->get_value(mc,ind,mu)) + t72*(*CI_users[8]->get_value(mc,ind,mu)) + t86*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_ppmmp_nf_wCI::\
C5g_ppmmp_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C5g_ppmmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  ppmmp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:56:36 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa13 = SPA(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> s23 = S(2,3);
complex<T> s45 = S(4,5);
complex<T> t5 = -(spa14*spb12); 
complex<T> t7 = spa12*T(3); 
complex<T> t11 = cube(spa34); 
complex<T> t12 = spa13*spb15; 
complex<T> t6 = -t12; 
complex<T> d1 = spa15*spa23*spa45*t7; d1 = T(1)/d1;
complex<T> d2 = (s23 - s45)*spa15*spa23*spa45*t7; d2 = T(1)/d2;
complex<T> d3 = spa15*t7*cube(s23 - s45); d3 = T(1)/d3;
complex<T> t9 = spa23*t5 + spa45*t6; 
complex<T> t17 = d3*t9; 
complex<T> t2 = -(d1*t11) + t12*t17*t5 - d2*t9*square(spa34); 
complex<T> t3 = -(d1*t11) + spa14*spb12*t12*t17 + d2*t9*square(spa34); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_ppmpm_nf_wCI::\
C5g_ppmpm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
} 
  
  
template <class T> SeriesC<T> 
     C5g_ppmpm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  ppmpm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:57:43 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spa14 = SPA(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = S(4,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> t5 = spa23*T(2); 
complex<T> t6 = square(spa14); 
complex<T> t8 = spa12*spa15; 
complex<T> t12 = spa34*spa45; 
complex<T> t14 = spa12*spa23; 
complex<T> t23 = square(spb14); 
complex<T> t24 = square(spa35); 
complex<T> t25 = square(spb24); 
complex<T> t26 = square(spa13); 
complex<T> t27 = square(spa25); 
complex<T> t30 = s15 - s34; 
complex<T> t32 = s15 - s23 + s45; 
complex<T> t33 = -(spa23*T(2)); 
complex<T> t34 = -(spa15*spa34*spb14) - spa23*spa45*spb24; 
complex<T> t39 = -(spa15*T(2)); 
complex<T> t42 = cube(spb24); 
complex<T> t44 = square(spa23); 
complex<T> t53 = square(spa34); 
complex<T> t54 = square(spa45); 
complex<T> t55 = cube(spb14); 
complex<T> t57 = square(spa15); 
complex<T> t58 = cube(spa25); 
complex<T> t59 = square(square(spa35)); 
complex<T> t69 = cube(spa13); 
complex<T> d15 = spa12*T(3); d15 = T(1)/d15;
complex<T> d16 = cube(s15 - s23); d16 = T(1)/d16;
complex<T> d23 = spa12*square(spa24); d23 = T(1)/d23;
complex<T> d24 = spa12*square(square(spa24)); d24 = T(1)/d24;
complex<T> t4 = square(t32); 
complex<T> t22 = t32*t6; 
complex<T> t38 = t23*t26; 
complex<T> t40 = -(t12*t27); 
complex<T> t64 = t12*t25; 
complex<T> t67 = d15*t34; 
complex<T> t73 = spa13*t23; 
complex<T> t74 = spa25*t24; 
complex<T> t75 = spa15*t12; 
complex<T> t91 = spa25*t42; 
complex<T> t100 = t53*t55; 
complex<T> d1 = spa23*t12*t8*T(3); d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*t14*t32; d2 = T(1)/d2;
complex<T> d3 = t14*t6*square(s15 - s23); d3 = T(1)/d3;
complex<T> d4 = (s15 - s23)*t14*t32*t6; d4 = T(1)/d4;
complex<T> d5 = (s15 - s23)*s24*t8; d5 = T(1)/d5;
complex<T> d6 = s24*t30*t8; d6 = T(1)/d6;
complex<T> d7 = t8*square(s15 - s23)*square(spa24); d7 = T(1)/d7;
complex<T> d8 = (s15 - s23)*s24*t8*square(spa24); d8 = T(1)/d8;
complex<T> d9 = t8*square(spa24)*square(t30); d9 = T(1)/d9;
complex<T> d10 = s24*t30*t8*square(spa24); d10 = T(1)/d10;
complex<T> d11 = spa24*t8*cube(t30)*T(3); d11 = T(1)/d11;
complex<T> d12 = (s15 - s23)*spa23*t12*t8*T(3); d12 = T(1)/d12;
complex<T> d13 = spa14*t14*T(3); d13 = T(1)/d13;
complex<T> d14 = spa24*t8*T(3); d14 = T(1)/d14;
complex<T> d17 = (s23 - s45)*t14*t32; d17 = T(1)/d17;
complex<T> d18 = t14*t6*square(s23 - s45); d18 = T(1)/d18;
complex<T> d19 = (s23 - s45)*t14*t32*t6; d19 = T(1)/d19;
complex<T> d20 = spa14*t14*cube(s23 - s45)*T(3); d20 = T(1)/d20;
complex<T> d25 = t8*square(spa24); d25 = T(1)/d25;
complex<T> d26 = t8*square(square(spa24)); d26 = T(1)/d26;
complex<T> d30 = t8*square(spa24)*T(2); d30 = T(1)/d30;
complex<T> t18 = s34*(d26*spa34*spa45*t27*t33 + d25*t74); 
complex<T> t29 = -t75; 
complex<T> t37 = -(d26*s23); 
complex<T> t60 = d11*t42; 
complex<T> t61 = d14*t44; 
complex<T> t62 = -t73; 
complex<T> t72 = d20*t55; 
complex<T> t79 = t27*t64; 
complex<T> t80 = t24*t73; 
complex<T> t82 = d13*t57; 
complex<T> t84 = t38*t39; 
complex<T> d21 = t14*t4; d21 = T(1)/d21;
complex<T> d22 = t14*t4*t6; d22 = T(1)/d22;
complex<T> d27 = spa12*t4; d27 = T(1)/d27;
complex<T> d28 = spa12*t4*t6; d28 = T(1)/d28;
complex<T> d29 = spa12*t4*t5; d29 = T(1)/d29;
complex<T> t1 = d7*spa23*t25*t40 + d9*spa23*t25*t40 - d1*t59 + d2*t24*t62 - d16*spa35*spb14*spb24*t67 + d5*t25*t74 + d6*t25*t74 + d3*t38*t75 + d10*t33*t79 + d8*t33*t79 - d12*t34*cube(spa35) - t53*t58*t60*T(2) + d4*t38*t75*T(2) - d16*spa13*t100*t82*T(2) - d16*t54*t61*t91*T(2); 
complex<T> t2 = d18*t29*t38 + d3*t29*t38 - d1*t59 + d17*t24*t62 + d16*spa35*spb14*spb24*t67 - d5*t25*t74 + d7*spa23*t79 + d8*t5*t79 + d2*t80 + d4*t12*t84 + d12*t34*cube(spa35) - t54*t69*t72*T(2) + d19*t38*t75*T(2) + d16*spa13*t100*t82*T(2) + d16*t54*t61*t91*T(2); 
complex<T> t20 = d26*s23*spa34*spa45*t27*t33 + d25*s23*t74 + d27*spb23*t80 + d28*spa34*spa45*spb23*t84; 
complex<T> t45 = d22*s15; 
complex<T> t50 = d18*t38*t75 + d17*t80 + d19*t12*t84 + t54*t69*t72*T(2); 
complex<T> t52 = -(d6*t25*t74) + d9*spa23*t79 + d10*t5*t79 + t53*t58*t60*T(2); 
complex<T> t71 = s34*(spa23*t12*t27*t37 + d30*s23*t74); 
complex<T> t78 = s45*(d21*t80 + d22*t12*t84); 
complex<T> t19 = d24*spa34*spa45*spb15*t27*t33 + d23*spb15*t74 + d21*s15*t80 + spa34*spa45*t45*t84; 
complex<T> t83 = s45*(t29*t38*t45 + d29*s15*t80); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t52*(*CI_users[2]->get_value(mc,ind,mu)) + t50*(*CI_users[3]->get_value(mc,ind,mu)) + t19*(*CI_users[4]->get_value(mc,ind,mu)) + t20*(*CI_users[5]->get_value(mc,ind,mu)) + t18*(*CI_users[6]->get_value(mc,ind,mu)) + t78*(*CI_users[7]->get_value(mc,ind,mu)) + t83*(*CI_users[8]->get_value(mc,ind,mu)) + t71*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pppmm_nf_wCI::\
C5g_pppmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pppmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  pppmm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:57:45 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb12 = SPB(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> s15 = S(1,5);
complex<T> s34 = S(3,4);
complex<T> t5 = -(spa24*spb12); 
complex<T> t7 = spa12*T(3); 
complex<T> t11 = cube(spa45); 
complex<T> t12 = spa25*spb23; 
complex<T> t6 = -t12; 
complex<T> d1 = spa15*spa23*spa34*t7; d1 = T(1)/d1;
complex<T> d2 = (s15 - s34)*spa15*spa23*spa34*t7; d2 = T(1)/d2;
complex<T> d3 = spa23*t7*cube(s15 - s34); d3 = T(1)/d3;
complex<T> t9 = spa15*t5 + spa34*t6; 
complex<T> t17 = d3*t9; 
complex<T> t2 = -(d1*t11) + t12*t17*t5 - d2*t9*square(spa45); 
complex<T> t3 = -(d1*t11) + spa24*spb12*t12*t17 + d2*t9*square(spa45); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_ppmmm_G_wCI::\
C5g_ppmmm_G_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_ppmmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  ppmmm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:57:46 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = S(1,2);
complex<T> t1 = spb34*spb45; 
complex<T> t5 = square(spb12); 
complex<T> t10 = cube(spb12); 
complex<T> t14 = spa45*spb24; 
complex<T> d5 = spb15*spb45*T(2); d5 = T(1)/d5;
complex<T> d6 = spb15*spb23*T(2); d6 = T(1)/d6;
complex<T> d7 = spb23*spb34*T(2); d7 = T(1)/d7;
complex<T> t11 = spa34*spb14*spb23 + spb15*t14; 
complex<T> d1 = spb15*spb23*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*spb15*spb23*t1*T(6); d2 = T(1)/d2;
complex<T> d3 = t1*cube(s15 - s23)*T(3); d3 = T(1)/d3;
complex<T> d4 = spb15*t1*T(2); d4 = T(1)/d4;
complex<T> d8 = spb23*t1*T(2); d8 = T(1)/d8;
complex<T> t13 = -(d1*T(11)); 
complex<T> t20 = d3*t11; 
complex<T> t3 = t10*t13 - spa34*spb14*t14*t20 - d2*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa34*spb14*t14*t20 + d2*t11*t5*T(11); 
complex<T> co1 = d4*s12*spa23*t10; 
complex<T> co2 = -(d5*spa23*spa34*t10); 
complex<T> co3 = -(d6*spa34*spa45*t10); 
complex<T> co4 = -(d7*spa15*spa45*t10); 
complex<T> co5 = d8*s12*spa15*t10; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmpmm_G_wCI::\
C5g_pmpmm_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmpmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmpmm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 19:59:22 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = -(spa12*spb12);
complex<T> t5 = square(spb25); 
complex<T> t6 = spb15*spb45; 
complex<T> t12 = spb12*spb23; 
complex<T> t13 = cube(spb13); 
complex<T> t14 = spb34*spb45; 
complex<T> t23 = square(spa25); 
complex<T> t24 = square(spb13); 
complex<T> t25 = square(spa24); 
complex<T> t26 = square(spb14); 
complex<T> t27 = square(spb35); 
complex<T> t30 = s15 - s23; 
complex<T> t32 = s12 + s15 - s34; 
complex<T> t33 = spa25*spb15*spb23 + spa24*spb12*spb34; 
complex<T> t34 = square(spb15); 
complex<T> t36 = square(spb34); 
complex<T> t45 = cube(spa24); 
complex<T> t46 = square(spb12); 
complex<T> t47 = square(spb23); 
complex<T> t49 = cube(spb35); 
complex<T> t58 = square(square(spb13)); 
complex<T> t60 = cube(spa25); 
complex<T> t61 = cube(spb14); 
complex<T> t86 = spa24*spb13; 
complex<T> d6 = spb45*T(3); d6 = T(1)/d6;
complex<T> d8 = cube(s15 - s34); d8 = T(1)/d8;
complex<T> d23 = spb45*square(spb24); d23 = T(1)/d23;
complex<T> d24 = spb45*square(square(spb24)); d24 = T(1)/d24;
complex<T> t4 = square(t32); 
complex<T> t22 = t32*t5; 
complex<T> t28 = -(t12*T(2)); 
complex<T> t42 = t24*T(2); 
complex<T> t66 = spb15*t23; 
complex<T> t67 = spb34*t26; 
complex<T> t68 = t12*t25; 
complex<T> t74 = t24*T(4); 
complex<T> t76 = spb14*t25; 
complex<T> t80 = spb35*t23; 
complex<T> t91 = -(spa12*t58); 
complex<T> t98 = t36*t45; 
complex<T> t102 = t49*t60; 
complex<T> d1 = (s12 - s34)*t14*t32; d1 = T(1)/d1;
complex<T> d2 = t14*t5*square(s12 - s34); d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*t14*t32*t5; d3 = T(1)/d3;
complex<T> d4 = spb25*t14*cube(s12 - s34)*T(3); d4 = T(1)/d4;
complex<T> d5 = spb24*t6*T(3); d5 = T(1)/d5;
complex<T> d7 = spb25*t14*T(3); d7 = T(1)/d7;
complex<T> d9 = s24*t30*t6; d9 = T(1)/d9;
complex<T> d10 = s24*(s15 - s34)*t6; d10 = T(1)/d10;
complex<T> d11 = spb24*t6*cube(t30)*T(3); d11 = T(1)/d11;
complex<T> d12 = spb34*t12*t6*T(6); d12 = T(1)/d12;
complex<T> d13 = t6*square(spb24)*square(t30); d13 = T(1)/d13;
complex<T> d14 = s24*t30*t6*square(spb24); d14 = T(1)/d14;
complex<T> d15 = t6*square(s15 - s34)*square(spb24); d15 = T(1)/d15;
complex<T> d16 = s24*(s15 - s34)*t6*square(spb24); d16 = T(1)/d16;
complex<T> d17 = (s15 - s34)*spb34*t12*t6*T(6); d17 = T(1)/d17;
complex<T> d18 = (s15 - s34)*t14*t32; d18 = T(1)/d18;
complex<T> d19 = t14*t5*square(s15 - s34); d19 = T(1)/d19;
complex<T> d20 = (s15 - s34)*t14*t32*t5; d20 = T(1)/d20;
complex<T> d25 = t6*square(spb24); d25 = T(1)/d25;
complex<T> d26 = t6*square(square(spb24)); d26 = T(1)/d26;
complex<T> d29 = spb34*t6*T(2); d29 = T(1)/d29;
complex<T> d30 = spb12*t6*T(2); d30 = T(1)/d30;
complex<T> d31 = spb15*t12*T(2); d31 = T(1)/d31;
complex<T> d32 = spb34*t12*T(2); d32 = T(1)/d32;
complex<T> d33 = spb23*t14*T(2); d33 = T(1)/d33;
complex<T> t37 = -(d12*T(11)); 
complex<T> t40 = d25*s23; 
complex<T> t43 = -t66; 
complex<T> t44 = -t67; 
complex<T> t51 = d26*s34; 
complex<T> t62 = d7*t34; 
complex<T> t63 = d11*t45; 
complex<T> t64 = d4*t46; 
complex<T> t69 = d1*spb35; 
complex<T> t75 = t27*t66; 
complex<T> t82 = t67*t68; 
complex<T> t87 = d5*spb14; 
complex<T> t89 = t28*t67; 
complex<T> d21 = t14*t4; d21 = T(1)/d21;
complex<T> d22 = t14*t4*t5; d22 = T(1)/d22;
complex<T> d27 = spb45*t4; d27 = T(1)/d27;
complex<T> d28 = spb45*t4*t5; d28 = T(1)/d28;
complex<T> t19 = spb14*t40*t74 - d26*s23*spb12*spb23*t67*T(2); 
complex<T> t38 = d21*s12; 
complex<T> t41 = -t51; 
complex<T> t50 = d22*s15; 
complex<T> t57 = t47*t61*t63*T(2) + t82*(d13 + d14*T(2)) - d9*t24*t76*T(4); 
complex<T> t79 = d25*s34*spb14*t74 + d28*spa34*t28*t75 + d27*spa34*t74*t80 + t51*t89; 
complex<T> t97 = spb12*t75; 
complex<T> t108 = t102*t64; 
complex<T> t114 = t87*t98; 
complex<T> t1 = t37*t58 + d13*t44*t68 + d15*t44*t68 + d19*t12*t75 + d10*t74*t76 + d9*t74*t76 + d14*t25*t89 + d16*t25*t89 - t47*t61*t63*T(2) + d20*t12*t75*T(2) + d8*(d6*spa25*t33*t86 - t114*t46*T(2) - spb35*t47*t60*t62*T(2)) - d18*t24*t80*T(4) + d17*t13*t33*T(11); 
complex<T> t2 = d19*t12*t27*t43 + d2*t12*t27*t43 + t37*t58 + t23*t69*t74 + d20*t28*t75 + d3*t28*t75 + d18*t74*t80 + d15*t82 - d6*d8*spa25*t33*t86 + t108*T(2) + d8*t114*t46*T(2) + d8*spb35*t47*t60*t62*T(2) + d16*t82*T(2) - d10*t24*t76*T(4) - d17*t13*t33*T(11); 
complex<T> t20 = s12*spb12*spb23*t27*t43*t50 + s15*t38*t42*t80 + d33*spa15*t91; 
complex<T> t21 = -(t108*T(2)) + spb23*t97*(d2 + d3*T(2)) - t23*t24*t69*T(4); 
complex<T> t72 = d22*s12*t28*t75 + t38*t74*t80; 
complex<T> t73 = s34*spb14*t40*t42 + s23*t12*t41*t67; 
complex<T> t84 = d23*spa15*spb14*t74 + t28*t50*t75 + d21*s15*t74*t80 + d24*spa15*t89; 
complex<T> co1 = d29*spa23*t91; 
complex<T> co2 = -(d30*spa23*spa34*t58); 
complex<T> co3 = -(d31*spa34*spa45*t58); 
complex<T> co4 = -(d32*spa15*spa45*t58); 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t21*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t57*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t72*(*CI_users[4]->get_value(mc,ind,mu)) + t84*(*CI_users[5]->get_value(mc,ind,mu)) + t19*(*CI_users[6]->get_value(mc,ind,mu)) + t79*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + co3*(*CI_users[10]->get_value(mc,ind,mu)) + t73*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + t20*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmmpm_G_wCI::\
C5g_pmmpm_G_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmmpm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmmpm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:00:58 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa25 = SPA(2,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa35 = SPA(3,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = -(spa35*spb35);
complex<T> t5 = spb12*spb23; 
complex<T> t6 = square(spb25); 
complex<T> t12 = spb15*spb45; 
complex<T> t13 = cube(spb14); 
complex<T> t14 = spb23*spb34; 
complex<T> t23 = square(spa25); 
complex<T> t24 = square(spb14); 
complex<T> t25 = square(spa35); 
complex<T> t26 = square(spb13); 
complex<T> t27 = square(spb24); 
complex<T> t29 = s12 - s34; 
complex<T> t30 = s15 - s34; 
complex<T> t31 = s12 + s15 - s34; 
complex<T> t32 = -(spa35*spb15*spb34) - spa25*spb12*spb45; 
complex<T> t33 = s12 - s45; 
complex<T> t34 = square(spb12); 
complex<T> t35 = square(spb34); 
complex<T> t45 = cube(spa35); 
complex<T> t46 = square(spb15); 
complex<T> t47 = square(spb45); 
complex<T> t49 = cube(spb24); 
complex<T> t58 = square(square(spb14)); 
complex<T> t60 = cube(spa25); 
complex<T> t61 = cube(spb13); 
complex<T> t75 = s34*spb15; 
complex<T> t90 = spa25*spb14; 
complex<T> d15 = spb23*T(3); d15 = T(1)/d15;
complex<T> d22 = spb23*square(spb35); d22 = T(1)/d22;
complex<T> d24 = spb23*square(square(spb35)); d24 = T(1)/d24;
complex<T> t4 = square(s15 + t29); 
complex<T> t28 = -(t12*T(2)); 
complex<T> t42 = t24*T(2); 
complex<T> t65 = -(t24*T(4)); 
complex<T> t66 = spb12*t23; 
complex<T> t67 = spb34*t26; 
complex<T> t68 = t12*t27; 
complex<T> t72 = t24*T(4); 
complex<T> t73 = spb24*t23; 
complex<T> t82 = spb13*t25; 
complex<T> t95 = t35*t45; 
complex<T> t99 = -(spa34*t58); 
complex<T> d1 = s35*t29*t5; d1 = T(1)/d1;
complex<T> d2 = s35*t33*t5; d2 = T(1)/d2;
complex<T> d3 = t14*t29*(s15 + t29); d3 = T(1)/d3;
complex<T> d4 = spb34*t12*t5*T(6); d4 = T(1)/d4;
complex<T> d5 = t14*t6*square(t29); d5 = T(1)/d5;
complex<T> d6 = t14*t29*(s15 + t29)*t6; d6 = T(1)/d6;
complex<T> d7 = t5*square(spb35)*square(t29); d7 = T(1)/d7;
complex<T> d8 = s35*t29*t5*square(spb35); d8 = T(1)/d8;
complex<T> d9 = t5*square(spb35)*square(t33); d9 = T(1)/d9;
complex<T> d10 = s35*t33*t5*square(spb35); d10 = T(1)/d10;
complex<T> d11 = spb35*t5*cube(t33)*T(3); d11 = T(1)/d11;
complex<T> d12 = spb34*t12*t29*t5*T(6); d12 = T(1)/d12;
complex<T> d13 = spb35*t5*T(3); d13 = T(1)/d13;
complex<T> d14 = spb25*t14*T(3); d14 = T(1)/d14;
complex<T> d16 = cube(t29); d16 = T(1)/d16;
complex<T> d17 = t14*(s15 + t29)*t30; d17 = T(1)/d17;
complex<T> d18 = spb25*t14*cube(t30)*T(3); d18 = T(1)/d18;
complex<T> d19 = t14*t6*square(t30); d19 = T(1)/d19;
complex<T> d20 = t14*(s15 + t29)*t30*t6; d20 = T(1)/d20;
complex<T> d26 = t5*square(spb35); d26 = T(1)/d26;
complex<T> d28 = t5*square(square(spb35)); d28 = T(1)/d28;
complex<T> d29 = spb34*t12*T(2); d29 = T(1)/d29;
complex<T> d30 = spb12*t12*T(2); d30 = T(1)/d30;
complex<T> d31 = spb15*t5*T(2); d31 = T(1)/d31;
complex<T> d32 = spb34*t5*T(2); d32 = T(1)/d32;
complex<T> d33 = spb45*t14*T(2); d33 = T(1)/d33;
complex<T> t37 = -(d4*T(11)); 
complex<T> t44 = -t67; 
complex<T> t51 = d28*s45; 
complex<T> t62 = d14*t34; 
complex<T> t63 = d11*t45; 
complex<T> t64 = d18*t46; 
complex<T> t74 = t25*t67; 
complex<T> t80 = t27*t66; 
complex<T> t89 = d13*spb13; 
complex<T> d21 = t14*t4; d21 = T(1)/d21;
complex<T> d23 = t14*t4*t6; d23 = T(1)/d23;
complex<T> d25 = spb23*t4; d25 = T(1)/d25;
complex<T> d27 = spb23*t4*t6; d27 = T(1)/d27;
complex<T> t2 = t37*t58 - d19*t66*t68 - d5*t66*t68 + d17*t72*t73 + d3*t72*t73 + d7*t12*t74 + d20*t28*t80 + d6*t28*t80 + d1*t65*t82 + d15*d16*spa35*t32*t90 + d16*spb24*t47*t60*t62*T(2) + t49*t60*t64*T(2) + d8*t12*t74*T(2) + d16*t46*t89*t95*T(2) + d12*t13*t32*T(11); 
complex<T> t19 = d17*t65*t73 - t49*t60*t64*T(2) + spb15*spb45*t80*(d19 + d20*T(2)); 
complex<T> t20 = d26*s34*spb13*t72 + d25*spa34*t72*t73 - spb45*(d28*t67*t75 + d27*spa34*spb15*t80)*T(2); 
complex<T> t21 = d26*s34*s45*spb13*t42 + spb45*t44*t51*t75 + d31*spa45*t99; 
complex<T> t38 = d21*s12; 
complex<T> t50 = d23*s15; 
complex<T> t78 = t28*t51*t67 + d26*s45*spb13*t72; 
complex<T> t102 = t61*t63; 
complex<T> t39 = -t50; 
complex<T> t71 = d21*s15*t72*t73 + t28*t50*t80; 
complex<T> t85 = d24*spa12*t28*t67 + d22*spa12*spb13*t72 + t38*t72*t73 + d23*s12*t28*t80; 
complex<T> t104 = t102*t47; 
complex<T> t1 = d7*t12*t25*t44 + d9*t12*t25*t44 + t37*t58 + d5*t66*t68 + d3*t65*t73 + d10*t28*t74 + d8*t28*t74 + d1*t72*t82 + d2*t72*t82 - d15*d16*spa35*t32*t90 - t104*T(2) - d16*spb24*t47*t60*t62*T(2) + d6*t66*t68*T(2) - d16*t46*t89*t95*T(2) - d12*t13*t32*T(11); 
complex<T> t22 = d9*spb15*spb45*t74 + d2*t65*t82 + t104*T(2) + d10*spb15*spb45*t74*T(2); 
complex<T> t79 = s12*t39*t66*t68 + s15*t38*t42*t73; 
complex<T> co1 = -(d29*spa12*spa23*t58); 
complex<T> co2 = d30*spa23*t99; 
complex<T> co3 = -(d32*spa15*spa45*t58); 
complex<T> co4 = -(d33*spa12*spa15*t58); 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t19*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t22*(*CI_users[3]->get_value(mc,ind,mu)) + t85*(*CI_users[4]->get_value(mc,ind,mu)) + t71*(*CI_users[5]->get_value(mc,ind,mu)) + t20*(*CI_users[6]->get_value(mc,ind,mu)) + t78*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t79*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t21*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmmmp_G_wCI::\
C5g_pmmmp_G_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmmmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmmmp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:01:00 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> s15 = S(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = spb23*spb34; 
complex<T> t5 = square(spb15); 
complex<T> t6 = -(spa34*spb13); 
complex<T> t10 = cube(spb15); 
complex<T> t14 = spa23*spb35; 
complex<T> d4 = spb34*spb45*T(2); d4 = T(1)/d4;
complex<T> d5 = spb12*spb45*T(2); d5 = T(1)/d5;
complex<T> d6 = spb12*spb23*T(2); d6 = T(1)/d6;
complex<T> t7 = -t14; 
complex<T> d1 = spb12*spb45*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = t1*cube(s12 - s45)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s12 - s45)*spb12*spb45*t1*T(6); d3 = T(1)/d3;
complex<T> d7 = spb12*t1*T(2); d7 = T(1)/d7;
complex<T> d8 = spb45*t1*T(2); d8 = T(1)/d8;
complex<T> t11 = spb45*t6 + spb12*t7; 
complex<T> t13 = -(d1*T(11)); 
complex<T> t17 = d2*t11; 
complex<T> t3 = t10*t13 + t14*t17*t6 - d3*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa34*spb13*t14*t17 + d3*t11*t5*T(11); 
complex<T> co1 = -(d4*spa12*spa23*t10); 
complex<T> co2 = -(d5*spa23*spa34*t10); 
complex<T> co3 = -(d6*spa34*spa45*t10); 
complex<T> co4 = d7*s15*spa45*t10; 
complex<T> co5 = d8*s15*spa12*t10; 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mppmm_G_wCI::\
C5g_mppmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mppmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mppmm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:01:02 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa12 = SPA(1,2);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s23 = S(2,3);
complex<T> t1 = spb15*spb45; 
complex<T> t5 = square(spb23); 
complex<T> t10 = cube(spb23); 
complex<T> t14 = spa15*spb35; 
complex<T> d6 = spb12*spb15*T(2); d6 = T(1)/d6;
complex<T> d7 = spb12*spb34*T(2); d7 = T(1)/d7;
complex<T> d8 = spb34*spb45*T(2); d8 = T(1)/d8;
complex<T> t11 = spa45*spb25*spb34 + spb12*t14; 
complex<T> d1 = spb12*spb34*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spb12*spb34*t1*T(6); d2 = T(1)/d2;
complex<T> d3 = t1*cube(s12 - s34)*T(3); d3 = T(1)/d3;
complex<T> d4 = spb34*t1*T(2); d4 = T(1)/d4;
complex<T> d5 = spb12*t1*T(2); d5 = T(1)/d5;
complex<T> t13 = -(d1*T(11)); 
complex<T> t20 = d3*t11; 
complex<T> t3 = t10*t13 - spa45*spb25*t14*t20 - d2*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa45*spb25*t14*t20 + d2*t11*t5*T(11); 
complex<T> co1 = d4*s23*spa12*t10; 
complex<T> co2 = d5*s23*spa34*t10; 
complex<T> co3 = -(d6*spa34*spa45*t10); 
complex<T> co4 = -(d7*spa15*spa45*t10); 
complex<T> co5 = -(d8*spa12*spa15*t10); 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mpmpm_G_wCI::\
C5g_mpmpm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
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
     C5g_mpmpm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mpmpm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:02:18 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa35 = SPA(3,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = -(spa35*spb35);
complex<T> t5 = spb12*spb15; 
complex<T> t11 = spb23*spb34; 
complex<T> t12 = cube(spb24); 
complex<T> t21 = square(spb24); 
complex<T> t22 = square(spa13); 
complex<T> t23 = square(spa35); 
complex<T> t24 = square(spb14); 
complex<T> t25 = square(spb25); 
complex<T> t27 = s12 - s45; 
complex<T> t28 = s12 - s34; 
complex<T> t29 = s23 - s45; 
complex<T> t30 = spa13*spb12*spb34 + spa35*spb23*spb45; 
complex<T> t31 = square(spb12); 
complex<T> t34 = square(spb45); 
complex<T> t43 = cube(spa13); 
complex<T> t44 = cube(spa35); 
complex<T> t45 = square(spb23); 
complex<T> t46 = square(spb34); 
complex<T> t56 = square(square(spb24)); 
complex<T> t57 = cube(spb14); 
complex<T> t58 = cube(spb25); 
complex<T> t79 = spa13*spb24; 
complex<T> d13 = spb13*spb15*spb45*T(3); d13 = T(1)/d13;
complex<T> d15 = spb15*T(3); d15 = T(1)/d15;
complex<T> d21 = spb15*square(spb35); d21 = T(1)/d21;
complex<T> d22 = spb15*spb45*square(spb13); d22 = T(1)/d22;
complex<T> d23 = spb15*spb45*square(square(spb13)); d23 = T(1)/d23;
complex<T> d24 = spb15*square(square(spb35)); d24 = T(1)/d24;
complex<T> d27 = spb15*square(spb13); d27 = T(1)/d27;
complex<T> d28 = spb15*square(square(spb13)); d28 = T(1)/d28;
complex<T> d29 = spb15*spb34*spb45*T(2); d29 = T(1)/d29;
complex<T> t26 = -(t11*T(2)); 
complex<T> t36 = d22*s12; 
complex<T> t64 = spb12*t24; 
complex<T> t65 = spb45*t25; 
complex<T> t66 = t11*t22; 
complex<T> t70 = t21*T(4); 
complex<T> t71 = t11*t23; 
complex<T> t72 = spb14*t22; 
complex<T> t78 = d13*spb14; 
complex<T> t86 = t31*t43; 
complex<T> t89 = t34*t44; 
complex<T> t90 = -(spa12*t56); 
complex<T> d1 = s35*t28*t5; d1 = T(1)/d1;
complex<T> d2 = s35*t27*t5; d2 = T(1)/d2;
complex<T> d3 = spb35*t5*cube(t28)*T(3); d3 = T(1)/d3;
complex<T> d4 = s13*spb15*spb45*t27; d4 = T(1)/d4;
complex<T> d5 = spb45*t11*t5*T(6); d5 = T(1)/d5;
complex<T> d6 = spb15*spb45*square(spb13)*square(t27); d6 = T(1)/d6;
complex<T> d7 = s13*spb15*spb45*t27*square(spb13); d7 = T(1)/d7;
complex<T> d8 = t5*square(spb35)*square(t28); d8 = T(1)/d8;
complex<T> d9 = s35*t28*t5*square(spb35); d9 = T(1)/d9;
complex<T> d10 = t5*square(spb35)*square(t27); d10 = T(1)/d10;
complex<T> d11 = s35*t27*t5*square(spb35); d11 = T(1)/d11;
complex<T> d12 = spb45*t11*t27*t5*T(6); d12 = T(1)/d12;
complex<T> d14 = spb35*t5*T(3); d14 = T(1)/d14;
complex<T> d16 = cube(t27); d16 = T(1)/d16;
complex<T> d17 = spb13*spb15*spb45*cube(t29)*T(3); d17 = T(1)/d17;
complex<T> d18 = s13*spb15*spb45*t29; d18 = T(1)/d18;
complex<T> d19 = spb15*spb45*square(spb13)*square(t29); d19 = T(1)/d19;
complex<T> d20 = s13*spb15*spb45*t29*square(spb13); d20 = T(1)/d20;
complex<T> d25 = t5*square(spb35); d25 = T(1)/d25;
complex<T> d26 = t5*square(square(spb35)); d26 = T(1)/d26;
complex<T> d30 = spb45*t5*T(2); d30 = T(1)/d30;
complex<T> d31 = spb23*t5*T(2); d31 = T(1)/d31;
complex<T> d32 = spb12*t11*T(2); d32 = T(1)/d32;
complex<T> d33 = spb45*t11*T(2); d33 = T(1)/d33;
complex<T> t17 = d22*s23*spb14*t70 - d23*s23*spb23*spb34*t64*T(2); 
complex<T> t35 = -(d5*T(11)); 
complex<T> t38 = d25*s34; 
complex<T> t41 = -t64; 
complex<T> t42 = -t65; 
complex<T> t48 = d26*s45; 
complex<T> t60 = d3*t44; 
complex<T> t61 = d17*t45; 
complex<T> t74 = d1*spb25; 
complex<T> t80 = d14*spb25; 
complex<T> t82 = t26*t65; 
complex<T> t92 = t26*t64; 
complex<T> t95 = t64*t66; 
complex<T> t97 = t65*t71; 
complex<T> t18 = d23*s12*s23*spb23*spb34*t41 + d29*spa23*t90 + s23*spb14*t21*t36*T(2); 
complex<T> t20 = s34*spb23*spb34*t42*t48 + s45*spb25*t21*t38*T(2); 
complex<T> t54 = d18*t70*t72 + d20*t22*t92 + d19*t95 - t43*t57*t61*T(2); 
complex<T> t55 = t46*t58*t60*T(2) + t97*(d8 + d9*T(2)) - t21*t23*t74*T(4); 
complex<T> t69 = spb25*t38*t70 + d26*s34*t82; 
complex<T> t75 = d27*spa45*spb14*t70 + d25*s45*spb25*t70 + t48*t82 + d28*spa45*t92; 
complex<T> t81 = d21*spa12*spb25*t70 + spb14*t36*t70 + d24*spa12*t82 + d23*s12*t92; 
complex<T> t105 = t45*t80; 
complex<T> t108 = t105*t89; 
complex<T> t1 = t35*t56 + d19*t41*t66 + d6*t41*t66 + d10*t97 + t43*t57*t61*T(2) + d20*t95*T(2) + d7*t95*T(2) + d11*t97*T(2) + d16*(-(d15*spa35*t30*t79) + t108*T(2) + t46*t78*t86*T(2)) - d2*spb25*t21*t23*T(4) - d18*t21*t72*T(4) - d4*t21*t72*T(4) - d12*t12*t30*T(11); 
complex<T> t2 = t35*t56 + d2*spb25*t23*t70 + d10*t42*t71 + d8*t42*t71 + d4*t70*t72 + t23*t70*t74 + d11*t23*t82 + d9*t23*t82 + d7*t22*t92 + d6*t95 - t46*t58*t60*T(2) + d16*(d15*spa35*t30*t79 - t108*T(2) - t46*t78*t86*T(2)) + d12*t12*t30*T(11); 
complex<T> co1 = -(d30*spa23*spa34*t56); 
complex<T> co2 = -(d31*spa34*spa45*t56); 
complex<T> co3 = -(d32*spa15*spa45*t56); 
complex<T> co4 = d33*spa15*t90; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t54*(*CI_users[1]->get_value(mc,ind,mu)) + t55*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t81*(*CI_users[4]->get_value(mc,ind,mu)) + t17*(*CI_users[5]->get_value(mc,ind,mu)) + t69*(*CI_users[6]->get_value(mc,ind,mu)) + t75*(*CI_users[7]->get_value(mc,ind,mu)) + t18*(*CI_users[8]->get_value(mc,ind,mu)) + co1*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + t20*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mpmmp_G_wCI::\
C5g_mpmmp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
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
     C5g_mpmmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mpmmp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:03:50 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa13 = SPA(1,3);
complex<T> spa23 = SPA(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> s12 = -(spa12*spb12);
complex<T> s23 = -(spa23*spb23);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = -(spa15*spb15);
complex<T> t5 = square(spb14); 
complex<T> t6 = spb34*spb45; 
complex<T> t12 = spb12*spb15; 
complex<T> t13 = cube(spb25); 
complex<T> t14 = spb23*spb34; 
complex<T> t22 = square(spa14); 
complex<T> t23 = square(spb25); 
complex<T> t24 = square(spa13); 
complex<T> t25 = square(spb24); 
complex<T> t26 = square(spb35); 
complex<T> t28 = s23 - s45; 
complex<T> t29 = s15 - s23; 
complex<T> t30 = s15 - s23 + s45; 
complex<T> t31 = -(spa13*spb15*spb23) - spa14*spb12*spb45; 
complex<T> t32 = s12 - s45; 
complex<T> t34 = square(spb23); 
complex<T> t35 = square(spb45); 
complex<T> t44 = cube(spa13); 
complex<T> t45 = square(spb12); 
complex<T> t46 = square(spb15); 
complex<T> t48 = cube(spb24); 
complex<T> t58 = square(square(spb25)); 
complex<T> t60 = cube(spa14); 
complex<T> t61 = cube(spb35); 
complex<T> t84 = spa13*spb25; 
complex<T> d19 = spb34*T(3); d19 = T(1)/d19;
complex<T> d27 = spb34*square(spb13); d27 = T(1)/d27;
complex<T> d28 = spb34*square(square(spb13)); d28 = T(1)/d28;
complex<T> t27 = -(t12*T(2)); 
complex<T> t41 = t23*T(2); 
complex<T> t66 = spb45*t22; 
complex<T> t67 = spb23*t26; 
complex<T> t68 = t12*t24; 
complex<T> t74 = t23*T(4); 
complex<T> t77 = spb35*t24; 
complex<T> t82 = spb24*t22; 
complex<T> t90 = t34*t44; 
complex<T> t93 = -(spa15*t58); 
complex<T> t101 = t48*t60; 
complex<T> d1 = s13*t32*t6; d1 = T(1)/d1;
complex<T> d2 = t6*square(spb13)*square(t32); d2 = T(1)/d2;
complex<T> d3 = s13*t32*t6*square(spb13); d3 = T(1)/d3;
complex<T> d4 = spb13*t6*cube(t32)*T(3); d4 = T(1)/d4;
complex<T> d5 = spb14*t14*cube(t29)*T(3); d5 = T(1)/d5;
complex<T> d6 = t14*t29*(s45 + t29); d6 = T(1)/d6;
complex<T> d7 = t14*t5*square(t29); d7 = T(1)/d7;
complex<T> d8 = t14*t29*(s45 + t29)*t5; d8 = T(1)/d8;
complex<T> d9 = t14*t28*(s45 + t29); d9 = T(1)/d9;
complex<T> d10 = spb23*t12*t6*T(6); d10 = T(1)/d10;
complex<T> d11 = s13*t28*t6; d11 = T(1)/d11;
complex<T> d12 = t6*square(spb13)*square(t28); d12 = T(1)/d12;
complex<T> d13 = s13*t28*t6*square(spb13); d13 = T(1)/d13;
complex<T> d14 = t14*t5*square(t28); d14 = T(1)/d14;
complex<T> d15 = t14*t28*(s45 + t29)*t5; d15 = T(1)/d15;
complex<T> d16 = spb23*t12*t28*t6*T(6); d16 = T(1)/d16;
complex<T> d17 = spb13*t6*T(3); d17 = T(1)/d17;
complex<T> d18 = spb14*t14*T(3); d18 = T(1)/d18;
complex<T> d20 = cube(t28); d20 = T(1)/d20;
complex<T> d21 = t6*square(spb13); d21 = T(1)/d21;
complex<T> d22 = t6*square(square(spb13)); d22 = T(1)/d22;
complex<T> d23 = t14*square(s45 + t29); d23 = T(1)/d23;
complex<T> d24 = t14*t5*square(s45 + t29); d24 = T(1)/d24;
complex<T> d25 = spb34*square(s45 + t29); d25 = T(1)/d25;
complex<T> d26 = spb34*t5*square(s45 + t29); d26 = T(1)/d26;
complex<T> d29 = spb15*t6*T(2); d29 = T(1)/d29;
complex<T> d30 = spb45*t12*T(2); d30 = T(1)/d30;
complex<T> d31 = spb23*t12*T(2); d31 = T(1)/d31;
complex<T> d32 = spb12*t14*T(2); d32 = T(1)/d32;
complex<T> d33 = spb23*t6*T(2); d33 = T(1)/d33;
complex<T> t36 = -(d10*T(11)); 
complex<T> t37 = d21*s12; 
complex<T> t38 = d23*s15; 
complex<T> t42 = -t66; 
complex<T> t43 = -t67; 
complex<T> t49 = d22*s23; 
complex<T> t50 = d24*s45; 
complex<T> t62 = d18*t35; 
complex<T> t63 = d4*t44; 
complex<T> t64 = d5*t46; 
complex<T> t69 = d6*spb24; 
complex<T> t75 = t25*t66; 
complex<T> t85 = d17*spb35; 
complex<T> t91 = t27*t67; 
complex<T> t94 = t67*t68; 
complex<T> t18 = spb35*t37*t74 - d22*s12*spb12*spb15*t67*T(2); 
complex<T> t21 = s15*spb12*spb15*t25*t42*t50 + s45*t38*t41*t82 + d32*spa45*t93; 
complex<T> t39 = -t49; 
complex<T> t56 = d1*t74*t77 + d3*t24*t91 + d2*t94 - t45*t61*t63*T(2); 
complex<T> t72 = d24*s15*t27*t75 + t38*t74*t82; 
complex<T> t80 = d21*s23*spb35*t74 + d26*spa23*t27*t75 + d25*spa23*t74*t82 + t49*t91; 
complex<T> t86 = d27*spa45*spb35*t74 + t27*t50*t75 + d23*s45*t74*t82 + d28*spa45*t91; 
complex<T> t97 = spb15*t75; 
complex<T> t107 = t101*t64; 
complex<T> t110 = t85*t90; 
complex<T> t1 = d14*t12*t25*t42 + d7*t12*t25*t42 + t36*t58 + t22*t69*t74 + d8*t27*t75 + d11*t74*t77 - d19*d20*spa14*t31*t84 + d13*t24*t91 + d12*t94 + t107*T(2) - d20*t110*t46*T(2) - d20*spb24*t45*t60*t62*T(2) + d15*t12*t75*T(2) - d9*t23*t82*T(4) - d16*t13*t31*T(11); 
complex<T> t2 = t36*t58 + d12*t43*t68 + d2*t43*t68 + d14*t12*t75 + d15*t27*t75 + d9*t74*t82 + t45*t61*t63*T(2) + d13*t94*T(2) + d3*t94*T(2) + d20*(d19*spa14*t31*t84 + t110*t46*T(2) + spb24*t45*t60*t62*T(2)) - d1*t23*t77*T(4) - d11*t23*t77*T(4) + d16*t13*t31*T(11); 
complex<T> t20 = -(t107*T(2)) + spb12*t97*(d7 + d8*T(2)) - t22*t23*t69*T(4); 
complex<T> t73 = s23*spb35*t37*t41 + s12*t12*t39*t67; 
complex<T> co1 = -(d29*spa12*spa23*t58); 
complex<T> co2 = -(d30*spa23*spa34*t58); 
complex<T> co3 = -(d31*spa34*spa45*t58); 
complex<T> co4 = d33*spa12*t93; 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t56*(*CI_users[0]->get_value(mc,ind,mu)) + t20*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t18*(*CI_users[4]->get_value(mc,ind,mu)) + t72*(*CI_users[5]->get_value(mc,ind,mu)) + t80*(*CI_users[6]->get_value(mc,ind,mu)) + t86*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + co2*(*CI_users[9]->get_value(mc,ind,mu)) + t73*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t21*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mmppm_G_wCI::\
C5g_mmppm_G_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mmppm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mmppm G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:03:52 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa15 = SPA(1,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb15 = SPB(1,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> s34 = S(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s45 = -(spa45*spb45);
complex<T> t1 = spb12*spb15; 
complex<T> t5 = square(spb34); 
complex<T> t10 = cube(spb34); 
complex<T> t14 = spa12*spb14; 
complex<T> d4 = spb15*spb45*T(2); d4 = T(1)/d4;
complex<T> d7 = spb12*spb23*T(2); d7 = T(1)/d7;
complex<T> d8 = spb23*spb45*T(2); d8 = T(1)/d8;
complex<T> t11 = spa15*spb13*spb45 + spb23*t14; 
complex<T> d1 = spb23*spb45*t1*T(6); d1 = T(1)/d1;
complex<T> d2 = t1*cube(s23 - s45)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s23 - s45)*spb23*spb45*t1*T(6); d3 = T(1)/d3;
complex<T> d5 = spb45*t1*T(2); d5 = T(1)/d5;
complex<T> d6 = spb23*t1*T(2); d6 = T(1)/d6;
complex<T> t13 = -(d1*T(11)); 
complex<T> t20 = d2*t11; 
complex<T> t3 = t10*t13 - spa15*spb13*t14*t20 - d3*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa15*spb13*t14*t20 + d3*t11*t5*T(11); 
complex<T> co1 = -(d4*spa12*spa23*t10); 
complex<T> co2 = d5*s34*spa23*t10; 
complex<T> co3 = d6*s34*spa45*t10; 
complex<T> co4 = -(d7*spa15*spa45*t10); 
complex<T> co5 = -(d8*spa12*spa15*t10); 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mmpmp_G_wCI::\
C5g_mmpmp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mmpmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mmpmp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:05:21 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa14 = SPA(1,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa24 = SPA(2,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s24 = -(spa24*spb24);
complex<T> t5 = square(spb14); 
complex<T> t6 = spb12*spb15; 
complex<T> t12 = spb34*spb45; 
complex<T> t13 = cube(spb35); 
complex<T> t14 = spb12*spb23; 
complex<T> t23 = square(spa14); 
complex<T> t24 = square(spb35); 
complex<T> t25 = square(spa24); 
complex<T> t26 = square(spb13); 
complex<T> t27 = square(spb25); 
complex<T> t30 = s15 - s34; 
complex<T> t32 = s15 - s23 + s45; 
complex<T> t33 = spa14*spb15*spb34 + spa24*spb23*spb45; 
complex<T> t34 = square(spb15); 
complex<T> t35 = square(spb23); 
complex<T> t45 = cube(spa24); 
complex<T> t46 = square(spb34); 
complex<T> t47 = square(spb45); 
complex<T> t49 = cube(spb13); 
complex<T> t58 = square(square(spb35)); 
complex<T> t60 = cube(spa14); 
complex<T> t61 = cube(spb25); 
complex<T> d15 = spb12*T(3); d15 = T(1)/d15;
complex<T> d16 = cube(s15 - s23); d16 = T(1)/d16;
complex<T> d22 = spb12*square(spb24); d22 = T(1)/d22;
complex<T> d24 = spb12*square(square(spb24)); d24 = T(1)/d24;
complex<T> t4 = square(t32); 
complex<T> t22 = t32*t5; 
complex<T> t28 = -(t12*T(2)); 
complex<T> t42 = t24*T(2); 
complex<T> t65 = -(t24*T(4)); 
complex<T> t66 = spb15*t23; 
complex<T> t67 = spb23*t27; 
complex<T> t68 = t12*t26; 
complex<T> t73 = t24*T(4); 
complex<T> t74 = spb13*t23; 
complex<T> t80 = spb25*t25; 
complex<T> t88 = spb35*t33; 
complex<T> t94 = -(spa23*t58); 
complex<T> t99 = t49*t60; 
complex<T> t101 = t35*t45; 
complex<T> d1 = spb24*t6*cube(t30)*T(3); d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*t14*t32; d2 = T(1)/d2;
complex<T> d3 = (s15 - s23)*s24*t6; d3 = T(1)/d3;
complex<T> d4 = s24*t30*t6; d4 = T(1)/d4;
complex<T> d5 = spb23*t12*t6*T(6); d5 = T(1)/d5;
complex<T> d6 = t14*t5*square(s15 - s23); d6 = T(1)/d6;
complex<T> d7 = (s15 - s23)*t14*t32*t5; d7 = T(1)/d7;
complex<T> d8 = t6*square(s15 - s23)*square(spb24); d8 = T(1)/d8;
complex<T> d9 = (s15 - s23)*s24*t6*square(spb24); d9 = T(1)/d9;
complex<T> d10 = t6*square(spb24)*square(t30); d10 = T(1)/d10;
complex<T> d11 = s24*t30*t6*square(spb24); d11 = T(1)/d11;
complex<T> d12 = (s15 - s23)*spb23*t12*t6*T(6); d12 = T(1)/d12;
complex<T> d13 = spb14*t14*T(3); d13 = T(1)/d13;
complex<T> d14 = spb24*t6*T(3); d14 = T(1)/d14;
complex<T> d17 = (s23 - s45)*t14*t32; d17 = T(1)/d17;
complex<T> d18 = t14*t5*square(s23 - s45); d18 = T(1)/d18;
complex<T> d19 = (s23 - s45)*t14*t32*t5; d19 = T(1)/d19;
complex<T> d20 = spb14*t14*cube(s23 - s45)*T(3); d20 = T(1)/d20;
complex<T> d26 = t6*square(spb24); d26 = T(1)/d26;
complex<T> d28 = t6*square(square(spb24)); d28 = T(1)/d28;
complex<T> d29 = spb15*t12*T(2); d29 = T(1)/d29;
complex<T> d30 = spb45*t6*T(2); d30 = T(1)/d30;
complex<T> d31 = spb23*t6*T(2); d31 = T(1)/d31;
complex<T> d32 = spb34*t14*T(2); d32 = T(1)/d32;
complex<T> d33 = spb23*t12*T(2); d33 = T(1)/d33;
complex<T> t37 = -(d5*T(11)); 
complex<T> t39 = d26*s23; 
complex<T> t43 = -t66; 
complex<T> t44 = -t67; 
complex<T> t50 = d28*s34; 
complex<T> t62 = d13*t34; 
complex<T> t63 = d1*t45; 
complex<T> t64 = d20*t47; 
complex<T> t75 = t25*t67; 
complex<T> t79 = t26*t66; 
complex<T> t84 = t66*t68; 
complex<T> t87 = d14*spb25; 
complex<T> d21 = t14*t4; d21 = T(1)/d21;
complex<T> d23 = t14*t4*t5; d23 = T(1)/d23;
complex<T> d25 = spb12*t4; d25 = T(1)/d25;
complex<T> d27 = spb12*t4*t5; d27 = T(1)/d27;
complex<T> t19 = s45*(d21*t73*t74 - d23*spb34*spb45*t79*T(2)); 
complex<T> t20 = s34*spb25*t39*t42 + s23*spb34*spb45*t44*t50 + d30*spa34*t94; 
complex<T> t21 = d10*spb34*spb45*t75 + d4*t65*t80 + t46*t61*t63*T(2) + d11*spb34*spb45*t75*T(2); 
complex<T> t38 = d21*s15; 
complex<T> t41 = -(d23*s45); 
complex<T> t71 = d28*s23*t28*t67 + spb25*t39*t73 + d25*spa23*t73*t74 + d27*spa23*t28*t79; 
complex<T> t77 = t28*t50*t67 + d26*s34*spb25*t73; 
complex<T> t103 = t47*t87; 
complex<T> t110 = t64*t99; 
complex<T> t57 = d17*t73*t74 + d19*t28*t79 + d18*t84 + t110*T(2); 
complex<T> t72 = s45*t38*t42*t74 + s15*t41*t84; 
complex<T> t83 = d24*spa15*t28*t67 + d22*spa15*spb25*t73 + t38*t73*t74 + d23*s15*t28*t79; 
complex<T> t105 = t101*t103; 
complex<T> t1 = t37*t58 + d18*t43*t68 + d6*t43*t68 + d17*t65*t74 + d2*t73*t74 + d8*t12*t75 + d7*t28*t79 + d3*t65*t80 - t110*T(2) + d9*t12*t75*T(2) + d19*t84*T(2) + d16*(-(d15*spa14*spa24*t88) + t105*T(2) + spb13*t46*t60*t62*T(2)) - d12*t13*t33*T(11); 
complex<T> t2 = d10*t12*t25*t44 + d8*t12*t25*t44 + t37*t58 + d2*t65*t74 + d11*t28*t75 + d9*t28*t75 + d3*t73*t80 + d4*t73*t80 + d6*t84 - t46*t61*t63*T(2) + d7*t84*T(2) + d16*(d15*spa14*spa24*t88 - t105*T(2) - spb13*t46*t60*t62*T(2)) + d12*t13*t33*T(11); 
complex<T> co1 = d29*spa12*t94; 
complex<T> co2 = -(d31*spa34*spa45*t58); 
complex<T> co3 = -(d32*spa15*spa45*t58); 
complex<T> co4 = -(d33*spa12*spa15*t58); 
complex<T> co5 = Complex(0,1); 
SeriesC<T> result = co5*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t21*(*CI_users[2]->get_value(mc,ind,mu)) + t57*(*CI_users[3]->get_value(mc,ind,mu)) + t83*(*CI_users[4]->get_value(mc,ind,mu)) + t71*(*CI_users[5]->get_value(mc,ind,mu)) + t77*(*CI_users[6]->get_value(mc,ind,mu)) + t19*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t72*(*CI_users[9]->get_value(mc,ind,mu)) + t20*(*CI_users[10]->get_value(mc,ind,mu)) + co2*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + co4*(*CI_users[13]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mmmpp_G_wCI::\
C5g_mmmpp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mmmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C5g :  mmmpp G");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:05:23 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb34 = SPB(3,4);
complex<T> spa15 = SPA(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> s45 = S(4,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s34 = -(spa34*spb34);
complex<T> t1 = spb12*spb23; 
complex<T> t5 = square(spb45); 
complex<T> t10 = cube(spb45); 
complex<T> t14 = spa23*spb25; 
complex<T> d4 = spb15*spb34*T(2); d4 = T(1)/d4;
complex<T> d5 = spb12*spb15*T(2); d5 = T(1)/d5;
complex<T> d8 = spb23*spb34*T(2); d8 = T(1)/d8;
complex<T> t11 = spa12*spb15*spb24 + spb34*t14; 
complex<T> d1 = t1*cube(s15 - s34)*T(3); d1 = T(1)/d1;
complex<T> d2 = (s15 - s34)*spb15*spb34*t1*T(6); d2 = T(1)/d2;
complex<T> d3 = spb15*spb34*t1*T(6); d3 = T(1)/d3;
complex<T> d6 = spb15*t1*T(2); d6 = T(1)/d6;
complex<T> d7 = spb34*t1*T(2); d7 = T(1)/d7;
complex<T> t13 = -(d3*T(11)); 
complex<T> t20 = d1*t11; 
complex<T> t3 = t10*t13 - spa12*spb24*t14*t20 - d2*t11*t5*T(11); 
complex<T> t4 = t10*t13 + spa12*spb24*t14*t20 + d2*t11*t5*T(11); 
complex<T> co1 = -(d4*spa12*spa23*t10); 
complex<T> co2 = -(d5*spa23*spa34*t10); 
complex<T> co3 = d6*s45*spa34*t10; 
complex<T> co4 = d7*s45*spa15*t10; 
complex<T> co5 = -(d8*spa12*spa15*t10); 
complex<T> co6 = Complex(0,1); 
SeriesC<T> result = co6*(t4*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + co1*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + co3*(*CI_users[4]->get_value(mc,ind,mu)) + co4*(*CI_users[5]->get_value(mc,ind,mu)) + co5*(*CI_users[6]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_ppmmm_nf_wCI::\
C5g_ppmmm_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C5g_ppmmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, p, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  ppmmm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:05:24 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> s15 = S(1,5);
complex<T> s23 = S(2,3);
complex<T> t8 = spb34*T(3); 
complex<T> t11 = cube(spb12); 
complex<T> t12 = spa45*spb24; 
complex<T> t9 = spa34*spb14*spb23 + spb15*t12; 
complex<T> d1 = spb15*spb23*spb45*t8; d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*spb15*spb23*spb45*t8; d2 = T(1)/d2;
complex<T> d3 = spb45*t8*cube(s15 - s23); d3 = T(1)/d3;
complex<T> t18 = d3*t12; 
complex<T> t2 = d1*t11 - t9*(spa34*spb14*t18 + d2*square(spb12)); 
complex<T> t3 = d1*t11 + spa34*spb14*t18*t9 + d2*t9*square(spb12); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmpmm_nf_wCI::\
C5g_pmpmm_nf_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c1, c2, c34));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmpmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmpmm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:06:34 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> s23 = S(2,3);
complex<T> s34 = -(spa34*spb34);
complex<T> s15 = -(spa15*spb15);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = S(1,2);
complex<T> t5 = spb34*T(2); 
complex<T> t6 = square(spb25); 
complex<T> t8 = spb15*spb45; 
complex<T> t12 = spb12*spb23; 
complex<T> t13 = -(spa25*spb13); 
complex<T> t14 = spb34*spb45; 
complex<T> t23 = square(spa25); 
complex<T> t24 = square(spb13); 
complex<T> t25 = square(spa24); 
complex<T> t26 = square(spb14); 
complex<T> t27 = square(spb35); 
complex<T> t30 = s15 - s23; 
complex<T> t32 = s12 + s15 - s34; 
complex<T> t33 = spa25*spb15*spb23 + spa24*spb12*spb34; 
complex<T> t34 = -(spb34*T(2)); 
complex<T> t39 = -(spb15*T(2)); 
complex<T> t41 = cube(spa24); 
complex<T> t44 = square(spb34); 
complex<T> t51 = cube(spa25); 
complex<T> t52 = square(spb12); 
complex<T> t53 = square(spb23); 
complex<T> t54 = square(square(spb13)); 
complex<T> t56 = square(spb15); 
complex<T> t57 = cube(spb35); 
complex<T> t67 = cube(spb14); 
complex<T> t80 = spb34*T(3); 
complex<T> d6 = spb45*T(3); d6 = T(1)/d6;
complex<T> d8 = cube(s15 - s34); d8 = T(1)/d8;
complex<T> d23 = spb45*square(spb24); d23 = T(1)/d23;
complex<T> d24 = spb45*square(square(spb24)); d24 = T(1)/d24;
complex<T> t22 = t32*t6; 
complex<T> t38 = t23*t27; 
complex<T> t61 = -(spb14*t24); 
complex<T> t62 = spb15*t12; 
complex<T> t63 = spb34*t25; 
complex<T> t66 = d6*t33; 
complex<T> t70 = t12*t26; 
complex<T> t82 = spb35*t23; 
complex<T> t107 = t51*t52; 
complex<T> d1 = (s12 - s34)*t14*t32; d1 = T(1)/d1;
complex<T> d2 = t14*t6*square(s12 - s34); d2 = T(1)/d2;
complex<T> d3 = (s12 - s34)*t14*t32*t6; d3 = T(1)/d3;
complex<T> d4 = spb25*t14*cube(s12 - s34)*T(3); d4 = T(1)/d4;
complex<T> d5 = spb24*t8*T(3); d5 = T(1)/d5;
complex<T> d7 = spb25*t14*T(3); d7 = T(1)/d7;
complex<T> d9 = s24*t30*t8; d9 = T(1)/d9;
complex<T> d10 = s24*(s15 - s34)*t8; d10 = T(1)/d10;
complex<T> d11 = spb24*t8*cube(t30)*T(3); d11 = T(1)/d11;
complex<T> d12 = t12*t8*t80; d12 = T(1)/d12;
complex<T> d13 = t8*square(spb24)*square(t30); d13 = T(1)/d13;
complex<T> d14 = s24*t30*t8*square(spb24); d14 = T(1)/d14;
complex<T> d15 = t8*square(s15 - s34)*square(spb24); d15 = T(1)/d15;
complex<T> d16 = s24*(s15 - s34)*t8*square(spb24); d16 = T(1)/d16;
complex<T> d17 = (s15 - s34)*t12*t8*t80; d17 = T(1)/d17;
complex<T> d18 = (s15 - s34)*t14*t32; d18 = T(1)/d18;
complex<T> d19 = t14*t6*square(s15 - s34); d19 = T(1)/d19;
complex<T> d20 = (s15 - s34)*t14*t32*t6; d20 = T(1)/d20;
complex<T> d21 = t14*square(t32); d21 = T(1)/d21;
complex<T> d22 = t14*t6*square(t32); d22 = T(1)/d22;
complex<T> d25 = t8*square(spb24); d25 = T(1)/d25;
complex<T> d26 = t8*square(square(spb24)); d26 = T(1)/d26;
complex<T> d27 = spb45*square(t32); d27 = T(1)/d27;
complex<T> d28 = spb45*t6*square(t32); d28 = T(1)/d28;
complex<T> d29 = t8*square(spb24)*T(2); d29 = T(1)/d29;
complex<T> d30 = spb45*t5*square(t32); d30 = T(1)/d30;
complex<T> t29 = -t62; 
complex<T> t36 = d22*s12; 
complex<T> t37 = d26*s23; 
complex<T> t40 = -t70; 
complex<T> t58 = d11*t41; 
complex<T> t59 = d5*t44; 
complex<T> t60 = -t82; 
complex<T> t77 = t38*t62; 
complex<T> t88 = d7*t56; 
complex<T> t101 = spb12*t38; 
complex<T> t2 = d19*t29*t38 + d20*t12*t38*t39 + d12*t54 + d10*t25*t61 + d9*t25*t61 + d8*spa24*t13*t66 + d14*t25*t5*t70 + d16*t25*t5*t70 + d13*t63*t70 + d15*t63*t70 + d18*t24*t82 - d17*t33*cube(spb13) + d8*spb14*t41*t52*t59*T(2) + t53*t58*t67*T(2) + d8*spb35*t51*t53*t88*T(2); 
complex<T> t19 = s34*(spb12*spb23*spb34*t26*t37 + d29*s23*t61); 
complex<T> t48 = d2*t29*t38 + d3*t12*t38*t39 + d1*t24*t82 + d4*t107*t57*T(2); 
complex<T> t49 = d9*spb14*t24*t25 + d13*t40*t63 + d14*t25*t34*t70 - t53*t58*t67*T(2); 
complex<T> t71 = t24*t60; 
complex<T> t75 = d25*s23*t61 + t37*t5*t70; 
complex<T> t1 = d10*spb14*t24*t25 + d12*t54 + d15*t40*t63 + d8*spa24*spa25*spb13*t66 + d16*t25*t34*t70 + d1*t71 + d18*t71 + d19*t77 + d2*t77 + d17*t33*cube(spb13) - d4*t107*t57*T(2) - d8*spb14*t41*t52*t59*T(2) + d20*t77*T(2) + d3*t77*T(2) - d8*spb35*t51*t53*t88*T(2); 
complex<T> t20 = d26*s34*spb12*spb23*t26*t5 + d25*s34*t61 + d27*spa34*t71 + d28*spa34*spb15*spb23*t101*T(2); 
complex<T> t21 = d24*spa15*spb12*spb23*t26*t5 + d23*spa15*t61 + d21*s15*t71 + d22*s15*spb15*spb23*t101*T(2); 
complex<T> t69 = d21*s12*t71 + t36*t77*T(2); 
complex<T> t76 = s15*(d30*s12*t71 + t36*t77); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t48*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t49*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t69*(*CI_users[4]->get_value(mc,ind,mu)) + t21*(*CI_users[5]->get_value(mc,ind,mu)) + t75*(*CI_users[6]->get_value(mc,ind,mu)) + t20*(*CI_users[7]->get_value(mc,ind,mu)) + t19*(*CI_users[8]->get_value(mc,ind,mu)) + t76*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmmpm_nf_wCI::\
C5g_pmmpm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c2, c1, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
} 
  
  
template <class T> SeriesC<T> 
     C5g_pmmpm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmmpm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:07:45 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> s15 = S(1,5);
complex<T> s12 = -(spa12*spb12);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = S(4,5);
complex<T> s35 = -(spa35*spb35);
complex<T> t5 = spb34*T(2); 
complex<T> t6 = square(spb25); 
complex<T> t8 = spb12*spb23; 
complex<T> t12 = spb15*spb45; 
complex<T> t14 = spb23*spb34; 
complex<T> t22 = square(spa25); 
complex<T> t23 = square(spb14); 
complex<T> t24 = square(spa35); 
complex<T> t25 = square(spb13); 
complex<T> t26 = square(spb24); 
complex<T> t27 = s12 - s34; 
complex<T> t29 = s15 - s34; 
complex<T> t30 = s12 + s15 - s34; 
complex<T> t31 = -(spa35*spb15*spb34) - spa25*spb12*spb45; 
complex<T> t32 = s12 - s45; 
complex<T> t33 = -(spb34*T(2)); 
complex<T> t38 = -(spb12*T(2)); 
complex<T> t40 = cube(spa35); 
complex<T> t51 = cube(spa25); 
complex<T> t52 = square(spb15); 
complex<T> t53 = square(spb45); 
complex<T> t54 = square(spb12); 
complex<T> t56 = square(square(spb14)); 
complex<T> t57 = cube(spb24); 
complex<T> t58 = square(spb34); 
complex<T> t68 = cube(spb13); 
complex<T> t73 = spb12*T(2); 
complex<T> t74 = spb34*T(3); 
complex<T> d15 = spb23*T(3); d15 = T(1)/d15;
complex<T> d22 = spb23*square(spb35); d22 = T(1)/d22;
complex<T> d24 = spb23*square(square(spb35)); d24 = T(1)/d24;
complex<T> t37 = t22*t26; 
complex<T> t61 = -(spb13*t23); 
complex<T> t62 = spb12*t12; 
complex<T> t63 = spb34*t24; 
complex<T> t67 = d15*t31; 
complex<T> t71 = t12*t25; 
complex<T> t83 = spb24*t22; 
complex<T> t84 = spb15*t25; 
complex<T> t93 = spb24*t54; 
complex<T> t109 = t51*t57; 
complex<T> t112 = t51*t53; 
complex<T> d1 = s35*t27*t8; d1 = T(1)/d1;
complex<T> d2 = s35*t32*t8; d2 = T(1)/d2;
complex<T> d3 = t14*t27*(s15 + t27); d3 = T(1)/d3;
complex<T> d4 = t12*t74*t8; d4 = T(1)/d4;
complex<T> d5 = t14*t6*square(t27); d5 = T(1)/d5;
complex<T> d6 = t14*t27*(s15 + t27)*t6; d6 = T(1)/d6;
complex<T> d7 = t8*square(spb35)*square(t27); d7 = T(1)/d7;
complex<T> d8 = s35*t27*t8*square(spb35); d8 = T(1)/d8;
complex<T> d9 = t8*square(spb35)*square(t32); d9 = T(1)/d9;
complex<T> d10 = s35*t32*t8*square(spb35); d10 = T(1)/d10;
complex<T> d11 = spb35*t8*cube(t32)*T(3); d11 = T(1)/d11;
complex<T> d12 = t12*t27*t74*t8; d12 = T(1)/d12;
complex<T> d13 = spb35*t8*T(3); d13 = T(1)/d13;
complex<T> d14 = spb25*t14*T(3); d14 = T(1)/d14;
complex<T> d16 = cube(t27); d16 = T(1)/d16;
complex<T> d17 = t14*(s15 + t27)*t29; d17 = T(1)/d17;
complex<T> d18 = spb25*t14*cube(t29)*T(3); d18 = T(1)/d18;
complex<T> d19 = t14*t6*square(t29); d19 = T(1)/d19;
complex<T> d20 = t14*(s15 + t27)*t29*t6; d20 = T(1)/d20;
complex<T> d21 = t14*square(s15 + t27); d21 = T(1)/d21;
complex<T> d23 = t14*t6*square(s15 + t27); d23 = T(1)/d23;
complex<T> d25 = spb23*square(s15 + t27); d25 = T(1)/d25;
complex<T> d26 = t8*square(spb35); d26 = T(1)/d26;
complex<T> d27 = spb23*t6*square(s15 + t27); d27 = T(1)/d27;
complex<T> d28 = t8*square(square(spb35)); d28 = T(1)/d28;
complex<T> d29 = spb23*t5*square(s15 + t27); d29 = T(1)/d29;
complex<T> d30 = t8*square(spb35)*T(2); d30 = T(1)/d30;
complex<T> t28 = -t62; 
complex<T> t35 = d23*s12; 
complex<T> t36 = d28*s34; 
complex<T> t39 = -t71; 
complex<T> t59 = d11*t40; 
complex<T> t60 = -t83; 
complex<T> t76 = s45*(d26*t61 + d28*t5*t71); 
complex<T> t78 = t37*t62; 
complex<T> t82 = d13*t58; 
complex<T> t101 = t12*t37; 
complex<T> t19 = s45*(d30*s34*t61 + spb34*spb45*t36*t84); 
complex<T> t48 = d19*t28*t37 + d20*t101*t38 + d17*t23*t83 + d18*t109*t52*T(2); 
complex<T> t50 = d2*spb13*t23*t24 + d9*t39*t63 + d10*t24*t33*t71 - t53*t59*t68*T(2); 
complex<T> t72 = t23*t60; 
complex<T> t106 = t40*t82; 
complex<T> t1 = d1*spb13*t23*t24 + d4*t56 + d7*t39*t63 - d16*spa25*spa35*spb14*t67 + d8*t24*t33*t71 + d17*t72 + d3*t72 + d19*t78 + d5*t78 - d12*t31*cube(spb14) - d16*spb13*t106*t52*T(2) - d18*t109*t52*T(2) + d20*t78*T(2) + d6*t78*T(2) - d14*d16*t112*t93*T(2); 
complex<T> t2 = d5*t28*t37 + d6*t101*t38 + d4*t56 + d1*t24*t61 + d2*t24*t61 + d16*spa25*spa35*spb14*t67 + d10*t24*t5*t71 + d8*t24*t5*t71 + d7*t63*t71 + d9*t63*t71 + d3*t23*t83 + d12*t31*cube(spb14) + d16*spb13*t106*t52*T(2) + t53*t59*t68*T(2) + d14*d16*t112*t93*T(2); 
complex<T> t20 = d26*s34*t61 + d25*spa34*t72 + d27*spa34*spb15*spb45*t37*t73 + spb45*t36*t5*t84; 
complex<T> t21 = d22*spa12*t61 + d21*s12*t72 + spb15*spb45*t35*t37*t73 + d24*spa12*spb45*t5*t84; 
complex<T> t70 = s15*(d21*t72 + d23*t78*T(2)); 
complex<T> t77 = s15*(d29*s12*t72 + t35*t78); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t48*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t50*(*CI_users[3]->get_value(mc,ind,mu)) + t21*(*CI_users[4]->get_value(mc,ind,mu)) + t70*(*CI_users[5]->get_value(mc,ind,mu)) + t20*(*CI_users[6]->get_value(mc,ind,mu)) + t76*(*CI_users[7]->get_value(mc,ind,mu)) + t77*(*CI_users[8]->get_value(mc,ind,mu)) + t19*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_pmmmp_nf_wCI::\
C5g_pmmmp_nf_wCI
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
     C5g_pmmmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{p, m, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  pmmmp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:07:47 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb13 = SPB(1,3);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> t4 = square(spb15); 
complex<T> t5 = -(spa34*spb13); 
complex<T> t7 = spb23*T(3); 
complex<T> t11 = cube(spb15); 
complex<T> t12 = spa23*spb35; 
complex<T> t6 = -t12; 
complex<T> d1 = spb12*spb34*spb45*t7; d1 = T(1)/d1;
complex<T> d2 = spb34*t7*cube(s12 - s45); d2 = T(1)/d2;
complex<T> d3 = (s12 - s45)*spb12*spb34*spb45*t7; d3 = T(1)/d3;
complex<T> t9 = spb45*t5 + spb12*t6; 
complex<T> t17 = d2*t9; 
complex<T> t2 = d1*t11 + t12*t17*t5 - d3*t4*t9; 
complex<T> t3 = d1*t11 + spa34*spb13*t12*t17 + d3*t4*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t3*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mppmm_nf_wCI::\
C5g_mppmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mppmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mppmm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:07:48 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = S(1,2);
complex<T> s34 = S(3,4);
complex<T> t7 = spb15*T(3); 
complex<T> t11 = cube(spb23); 
complex<T> t12 = spa15*spb35; 
complex<T> t9 = spa45*spb25*spb34 + spb12*t12; 
complex<T> d1 = spb12*spb34*spb45*t7; d1 = T(1)/d1;
complex<T> d2 = (s12 - s34)*spb12*spb34*spb45*t7; d2 = T(1)/d2;
complex<T> d3 = spb45*t7*cube(s12 - s34); d3 = T(1)/d3;
complex<T> t18 = d3*t12; 
complex<T> t2 = d1*t11 - t9*(spa45*spb25*t18 + d2*square(spb23)); 
complex<T> t3 = d1*t11 + spa45*spb25*t18*t9 + d2*t9*square(spb23); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mpmpm_nf_wCI::\
C5g_mpmpm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mpmpm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mpmpm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:08:40 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa13 = SPA(1,3);
complex<T> spa12 = SPA(1,2);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> s23 = S(2,3);
complex<T> s12 = -(spa12*spb12);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = -(spa13*spb13);
complex<T> s34 = S(3,4);
complex<T> s35 = -(spa35*spb35);
complex<T> t6 = spb12*spb15; 
complex<T> t9 = spb45*T(3); 
complex<T> t11 = spb23*spb34; 
complex<T> t22 = square(spb24); 
complex<T> t23 = square(spa13); 
complex<T> t24 = square(spa35); 
complex<T> t25 = square(spb14); 
complex<T> t26 = square(spb25); 
complex<T> t27 = s12 - s45; 
complex<T> t28 = spb23*T(2); 
complex<T> t29 = spb12*spb34; 
complex<T> t30 = s12 - s34; 
complex<T> t31 = s23 - s45; 
complex<T> t32 = spa13*spb12*spb34 + spa35*spb23*spb45; 
complex<T> t35 = square(spb45); 
complex<T> t40 = cube(spa13); 
complex<T> t41 = cube(spa35); 
complex<T> t52 = square(spb23); 
complex<T> t53 = square(spb34); 
complex<T> t54 = square(spb12); 
complex<T> t56 = square(square(spb24)); 
complex<T> t57 = cube(spb25); 
complex<T> t68 = cube(spb14); 
complex<T> t87 = spa13*spa35; 
complex<T> d15 = spb15*T(3); d15 = T(1)/d15;
complex<T> d21 = spb15*square(spb35); d21 = T(1)/d21;
complex<T> d22 = spb15*spb45*square(spb13); d22 = T(1)/d22;
complex<T> d23 = spb15*spb45*square(square(spb13)); d23 = T(1)/d23;
complex<T> d24 = spb15*square(square(spb35)); d24 = T(1)/d24;
complex<T> d27 = spb15*square(spb13); d27 = T(1)/d27;
complex<T> d28 = spb15*square(square(spb13)); d28 = T(1)/d28;
complex<T> d29 = spb15*spb45*square(spb13)*T(2); d29 = T(1)/d29;
complex<T> t36 = d23*s12; 
complex<T> t42 = -(spb14*t22); 
complex<T> t61 = spb45*t26; 
complex<T> t62 = spb12*t11; 
complex<T> t63 = t23*t25; 
complex<T> t64 = -(t22*t24); 
complex<T> t72 = t11*t24; 
complex<T> t83 = spb25*t22; 
complex<T> t86 = spb34*t28; 
complex<T> t91 = t28*t29; 
complex<T> t96 = t53*t54; 
complex<T> d1 = s35*t30*t6; d1 = T(1)/d1;
complex<T> d2 = s35*t27*t6; d2 = T(1)/d2;
complex<T> d3 = spb35*t6*cube(t30)*T(3); d3 = T(1)/d3;
complex<T> d4 = s13*spb15*spb45*t27; d4 = T(1)/d4;
complex<T> d5 = t11*t6*t9; d5 = T(1)/d5;
complex<T> d6 = spb15*spb45*square(spb13)*square(t27); d6 = T(1)/d6;
complex<T> d7 = s13*spb15*spb45*t27*square(spb13); d7 = T(1)/d7;
complex<T> d8 = t6*square(spb35)*square(t30); d8 = T(1)/d8;
complex<T> d9 = s35*t30*t6*square(spb35); d9 = T(1)/d9;
complex<T> d10 = t6*square(spb35)*square(t27); d10 = T(1)/d10;
complex<T> d11 = s35*t27*t6*square(spb35); d11 = T(1)/d11;
complex<T> d12 = t11*t27*t6*t9; d12 = T(1)/d12;
complex<T> d13 = spb13*spb15*t9; d13 = T(1)/d13;
complex<T> d14 = spb35*t6*T(3); d14 = T(1)/d14;
complex<T> d16 = cube(t27); d16 = T(1)/d16;
complex<T> d17 = spb13*spb15*t9*cube(t31); d17 = T(1)/d17;
complex<T> d18 = s13*spb15*spb45*t31; d18 = T(1)/d18;
complex<T> d19 = spb15*spb45*square(spb13)*square(t31); d19 = T(1)/d19;
complex<T> d20 = s13*spb15*spb45*t31*square(spb13); d20 = T(1)/d20;
complex<T> d25 = t6*square(spb35); d25 = T(1)/d25;
complex<T> d26 = t6*square(square(spb35)); d26 = T(1)/d26;
complex<T> d30 = t6*square(spb35)*T(2); d30 = T(1)/d30;
complex<T> t17 = s23*(d22*t42 + d23*t25*t91); 
complex<T> t37 = d26*s34; 
complex<T> t38 = -t61; 
complex<T> t39 = -t62; 
complex<T> t58 = d14*t35; 
complex<T> t59 = d13*t40; 
complex<T> t60 = d3*t41; 
complex<T> t70 = d17*t40; 
complex<T> t73 = -t83; 
complex<T> t74 = -(t63*T(2)); 
complex<T> t77 = s23*(d29*s12*t42 + t25*t36*t62); 
complex<T> t82 = t62*t63; 
complex<T> t19 = d25*s34*t73 + t37*t61*t86; 
complex<T> t20 = d27*spa45*t42 + d25*s45*t73 + d26*s45*t61*t86 + d28*spa45*t25*t91; 
complex<T> t21 = d22*s12*t42 + d21*spa12*t73 + d24*spa12*t61*t86 + t25*t36*t91; 
complex<T> t49 = d8*t38*t72 + d1*t24*t83 - t53*t57*t60*T(2) - d9*t61*t72*T(2); 
complex<T> t51 = d18*t23*t42 + d19*t39*t63 + t52*t68*t70*T(2) + d20*t82*T(2); 
complex<T> t71 = s45*(t11*t37*t61 + d30*s34*t73); 
complex<T> t85 = t41*t58; 
complex<T> t1 = d4*t23*t42 + d5*t56 + d6*t39*t63 + d1*spb25*t64 + d2*spb25*t64 + d10*t61*t72 + d8*t61*t72 - d15*d16*spb24*t32*t87 - d12*t32*cube(spb24) + t53*t57*t60*T(2) + d11*t61*t72*T(2) + d9*t61*t72*T(2) + d7*t82*T(2) + d16*spb25*t52*t85*T(2) + d16*spb14*t59*t96*T(2); 
complex<T> t2 = d18*spb14*t22*t23 + d4*spb14*t22*t23 + d5*t56 + d10*t38*t72 + d20*t62*t74 + d7*t62*t74 + d19*t82 + d6*t82 + d2*t24*t83 + d15*d16*spb24*t32*t87 + d12*t32*cube(spb24) - t52*t68*t70*T(2) - d11*t61*t72*T(2) - d16*spb25*t52*t85*T(2) - d16*spb14*t59*t96*T(2); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t51*(*CI_users[1]->get_value(mc,ind,mu)) + t49*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t21*(*CI_users[4]->get_value(mc,ind,mu)) + t17*(*CI_users[5]->get_value(mc,ind,mu)) + t19*(*CI_users[6]->get_value(mc,ind,mu)) + t20*(*CI_users[7]->get_value(mc,ind,mu)) + t77*(*CI_users[8]->get_value(mc,ind,mu)) + t71*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mpmmp_nf_wCI::\
C5g_mpmmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c2, c345));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c1, c45));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c1, c23));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mpmmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, p, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mpmmp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:09:53 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb13 = SPB(1,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa13 = SPA(1,3);
complex<T> spa14 = SPA(1,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb14 = SPB(1,4);
complex<T> spa45 = SPA(4,5);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s13 = -(spa13*spb13);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = S(1,5);
complex<T> t5 = spb23*T(2); 
complex<T> t6 = square(spb14); 
complex<T> t8 = spb34*spb45; 
complex<T> t12 = spb12*spb15; 
complex<T> t14 = spb23*spb34; 
complex<T> t22 = square(spa14); 
complex<T> t23 = square(spb25); 
complex<T> t24 = square(spa13); 
complex<T> t25 = square(spb24); 
complex<T> t26 = square(spb35); 
complex<T> t27 = s23 - s45; 
complex<T> t29 = s15 - s23; 
complex<T> t30 = s15 - s23 + s45; 
complex<T> t31 = -(spa13*spb15*spb23) - spa14*spb12*spb45; 
complex<T> t32 = s12 - s45; 
complex<T> t38 = -(spb45*T(2)); 
complex<T> t40 = cube(spa13); 
complex<T> t51 = cube(spa14); 
complex<T> t52 = square(spb12); 
complex<T> t53 = square(spb15); 
complex<T> t54 = square(spb23); 
complex<T> t56 = square(square(spb25)); 
complex<T> t57 = cube(spb35); 
complex<T> t58 = square(spb45); 
complex<T> t68 = cube(spb24); 
complex<T> t74 = spb45*T(2); 
complex<T> t75 = spb23*T(3); 
complex<T> d19 = spb34*T(3); d19 = T(1)/d19;
complex<T> d27 = spb34*square(spb13); d27 = T(1)/d27;
complex<T> d28 = spb34*square(square(spb13)); d28 = T(1)/d28;
complex<T> t37 = t22*t25; 
complex<T> t61 = -(spb35*t23); 
complex<T> t62 = spb45*t12; 
complex<T> t63 = spb23*t24; 
complex<T> t67 = d19*t31; 
complex<T> t72 = t12*t26; 
complex<T> t84 = spb24*t22; 
complex<T> t92 = t40*t57; 
complex<T> d1 = s13*t32*t8; d1 = T(1)/d1;
complex<T> d2 = t8*square(spb13)*square(t32); d2 = T(1)/d2;
complex<T> d3 = s13*t32*t8*square(spb13); d3 = T(1)/d3;
complex<T> d4 = spb13*t8*cube(t32)*T(3); d4 = T(1)/d4;
complex<T> d5 = spb14*t14*cube(t29)*T(3); d5 = T(1)/d5;
complex<T> d6 = t14*t29*(s45 + t29); d6 = T(1)/d6;
complex<T> d7 = t14*t6*square(t29); d7 = T(1)/d7;
complex<T> d8 = t14*t29*(s45 + t29)*t6; d8 = T(1)/d8;
complex<T> d9 = t14*t27*(s45 + t29); d9 = T(1)/d9;
complex<T> d10 = t12*t75*t8; d10 = T(1)/d10;
complex<T> d11 = s13*t27*t8; d11 = T(1)/d11;
complex<T> d12 = t8*square(spb13)*square(t27); d12 = T(1)/d12;
complex<T> d13 = s13*t27*t8*square(spb13); d13 = T(1)/d13;
complex<T> d14 = t14*t6*square(t27); d14 = T(1)/d14;
complex<T> d15 = t14*t27*(s45 + t29)*t6; d15 = T(1)/d15;
complex<T> d16 = t12*t27*t75*t8; d16 = T(1)/d16;
complex<T> d17 = spb13*t8*T(3); d17 = T(1)/d17;
complex<T> d18 = spb14*t14*T(3); d18 = T(1)/d18;
complex<T> d20 = cube(t27); d20 = T(1)/d20;
complex<T> d21 = t8*square(spb13); d21 = T(1)/d21;
complex<T> d22 = t8*square(square(spb13)); d22 = T(1)/d22;
complex<T> d23 = t14*square(s45 + t29); d23 = T(1)/d23;
complex<T> d24 = t14*t6*square(s45 + t29); d24 = T(1)/d24;
complex<T> d25 = spb34*square(s45 + t29); d25 = T(1)/d25;
complex<T> d26 = spb34*t6*square(s45 + t29); d26 = T(1)/d26;
complex<T> d29 = t8*square(spb13)*T(2); d29 = T(1)/d29;
complex<T> d30 = spb34*t5*square(s45 + t29); d30 = T(1)/d30;
complex<T> t28 = -t62; 
complex<T> t35 = d22*s12; 
complex<T> t36 = d24*s15; 
complex<T> t39 = -t72; 
complex<T> t59 = d17*t40; 
complex<T> t60 = -t84; 
complex<T> t70 = d5*t51; 
complex<T> t78 = t37*t62; 
complex<T> t83 = t24*t72; 
complex<T> t86 = d18*t58; 
complex<T> t106 = t12*t37; 
complex<T> t18 = s23*(spb12*spb15*spb23*t26*t35 + d29*s12*t61); 
complex<T> t48 = d7*t28*t37 + d8*t106*t38 + d6*t23*t84 + t53*t68*t70*T(2); 
complex<T> t50 = d1*t24*t61 + d2*t39*t63 + d3*t5*t83 + d4*t52*t92*T(2); 
complex<T> t71 = d21*s12*t61 + t35*t5*t72; 
complex<T> t73 = t23*t60; 
complex<T> t90 = t54*t59; 
complex<T> t96 = t51*t86; 
complex<T> t1 = d1*spb35*t23*t24 + d11*spb35*t23*t24 + d14*t28*t37 + d10*t56 - d20*spa13*spa14*spb25*t67 + d12*t63*t72 + d2*t63*t72 + d9*t73 - d16*t31*cube(spb25) - d13*t63*t72*T(2) - d3*t63*t72*T(2) + d15*t78*T(2) - d20*spb35*t53*t90*T(2) - d4*t52*t92*T(2) - d20*spb24*t52*t96*T(2); 
complex<T> t2 = d15*t106*t38 + d10*t56 + d11*t24*t61 + d12*t39*t63 + d20*spa13*spa14*spb25*t67 + d6*t73 + d14*t78 + d7*t78 + d13*t5*t83 + d9*t23*t84 + d16*t31*cube(spb25) - t53*t68*t70*T(2) + d8*t78*T(2) + d20*spb35*t53*t90*T(2) + d20*spb24*t52*t96*T(2); 
complex<T> t20 = d22*s23*spb12*spb15*t26*t5 + d21*s23*t61 + d25*spa23*t73 + d26*spa23*spb12*spb15*t37*t74; 
complex<T> t21 = d28*spa45*spb12*spb15*t26*t5 + d27*spa45*t61 + d23*s45*t73 + d24*s45*spb12*spb15*t37*t74; 
complex<T> t77 = d23*s15*t73 + t36*t78*T(2); 
complex<T> t82 = s45*(d30*s15*t73 + t36*t78); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t50*(*CI_users[0]->get_value(mc,ind,mu)) + t48*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t71*(*CI_users[4]->get_value(mc,ind,mu)) + t77*(*CI_users[5]->get_value(mc,ind,mu)) + t20*(*CI_users[6]->get_value(mc,ind,mu)) + t21*(*CI_users[7]->get_value(mc,ind,mu)) + t18*(*CI_users[8]->get_value(mc,ind,mu)) + t82*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mmppm_nf_wCI::\
C5g_mmppm_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C5g_mmppm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mmppm nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:09:54 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb12 = SPB(1,2);
complex<T> spb15 = SPB(1,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa12 = SPA(1,2);
complex<T> spa15 = SPA(1,5);
complex<T> spb13 = SPB(1,3);
complex<T> spb14 = SPB(1,4);
complex<T> s23 = S(2,3);
complex<T> s45 = S(4,5);
complex<T> t4 = square(spb34); 
complex<T> t7 = spb12*T(3); 
complex<T> t11 = cube(spb34); 
complex<T> t12 = spa12*spb14; 
complex<T> t9 = spa15*spb13*spb45 + spb23*t12; 
complex<T> d1 = spb15*spb23*spb45*t7; d1 = T(1)/d1;
complex<T> d2 = spb15*t7*cube(s23 - s45); d2 = T(1)/d2;
complex<T> d3 = (s23 - s45)*spb15*spb23*spb45*t7; d3 = T(1)/d3;
complex<T> t15 = d2*t9; 
complex<T> t2 = d1*t11 - spa15*spb13*t12*t15 - d3*t4*t9; 
complex<T> t3 = d1*t11 + spa15*spb13*t12*t15 + d3*t4*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mmpmp_nf_wCI::\
C5g_mmpmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mmpmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mmpmp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:11:01 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa14 = SPA(1,4);
complex<T> spb12 = SPB(1,2);
complex<T> spb13 = SPB(1,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb35 = SPB(3,5);
complex<T> spb14 = SPB(1,4);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa15 = SPA(1,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> s45 = S(4,5);
complex<T> s15 = -(spa15*spb15);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = S(3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> t5 = spb23*T(2); 
complex<T> t6 = square(spb14); 
complex<T> t12 = spb34*spb45; 
complex<T> t14 = spb12*spb23; 
complex<T> t23 = square(spa14); 
complex<T> t24 = square(spb35); 
complex<T> t25 = square(spa24); 
complex<T> t26 = square(spb13); 
complex<T> t27 = square(spb25); 
complex<T> t30 = s15 - s34; 
complex<T> t32 = s15 - s23 + s45; 
complex<T> t33 = spa14*spb15*spb34 + spa24*spb23*spb45; 
complex<T> t34 = -(spb23*T(2)); 
complex<T> t39 = -(spb15*T(2)); 
complex<T> t41 = cube(spa24); 
complex<T> t44 = square(spb23); 
complex<T> t51 = cube(spa14); 
complex<T> t52 = square(spb34); 
complex<T> t53 = square(spb45); 
complex<T> t55 = square(spb15); 
complex<T> t56 = cube(spb25); 
complex<T> t57 = square(square(spb35)); 
complex<T> t67 = cube(spb13); 
complex<T> t73 = spb15*T(2); 
complex<T> t74 = spb23*T(3); 
complex<T> d3 = (s15 - s23)*s24*spb12*spb15; d3 = T(1)/d3;
complex<T> d8 = spb12*spb15*square(s15 - s23)*square(spb24); d8 = T(1)/d8;
complex<T> d9 = (s15 - s23)*s24*spb12*spb15*square(spb24); d9 = T(1)/d9;
complex<T> d14 = spb12*spb15*spb24*T(3); d14 = T(1)/d14;
complex<T> d15 = spb12*T(3); d15 = T(1)/d15;
complex<T> d16 = cube(s15 - s23); d16 = T(1)/d16;
complex<T> d22 = spb12*square(spb24); d22 = T(1)/d22;
complex<T> d24 = spb12*square(square(spb24)); d24 = T(1)/d24;
complex<T> d26 = spb12*spb15*square(spb24); d26 = T(1)/d26;
complex<T> d28 = spb12*spb15*square(square(spb24)); d28 = T(1)/d28;
complex<T> t4 = square(t32); 
complex<T> t22 = t32*t6; 
complex<T> t37 = d28*s23; 
complex<T> t38 = t23*t26; 
complex<T> t59 = d14*t44; 
complex<T> t61 = -(spb25*t24); 
complex<T> t62 = spb15*t12; 
complex<T> t63 = spb23*t25; 
complex<T> t65 = d15*t33; 
complex<T> t71 = t12*t27; 
complex<T> t82 = spb13*t23; 
complex<T> t83 = spb34*t27; 
complex<T> t95 = spb25*t41; 
complex<T> d1 = spb12*spb15*spb24*cube(t30)*T(3); d1 = T(1)/d1;
complex<T> d2 = (s15 - s23)*t14*t32; d2 = T(1)/d2;
complex<T> d4 = s24*spb12*spb15*t30; d4 = T(1)/d4;
complex<T> d6 = t14*t6*square(s15 - s23); d6 = T(1)/d6;
complex<T> d7 = (s15 - s23)*t14*t32*t6; d7 = T(1)/d7;
complex<T> d10 = spb12*spb15*square(spb24)*square(t30); d10 = T(1)/d10;
complex<T> d11 = s24*spb12*spb15*t30*square(spb24); d11 = T(1)/d11;
complex<T> d13 = spb14*t14*T(3); d13 = T(1)/d13;
complex<T> d17 = (s23 - s45)*t14*t32; d17 = T(1)/d17;
complex<T> d18 = t14*t6*square(s23 - s45); d18 = T(1)/d18;
complex<T> d19 = (s23 - s45)*t14*t32*t6; d19 = T(1)/d19;
complex<T> d20 = spb14*t14*cube(s23 - s45)*T(3); d20 = T(1)/d20;
complex<T> d30 = spb12*t73*square(spb24); d30 = T(1)/d30;
complex<T> t19 = s34*(d30*s23*t61 + spb23*spb45*t37*t83); 
complex<T> t29 = -t62; 
complex<T> t40 = -t71; 
complex<T> t58 = d1*t41; 
complex<T> t60 = -t82; 
complex<T> t69 = d20*t51; 
complex<T> t70 = s34*(d26*t61 + d28*t5*t71); 
complex<T> t77 = t38*t62; 
complex<T> t87 = d13*t55; 
complex<T> t102 = t12*t38; 
complex<T> d5 = t14*t62*T(3); d5 = T(1)/d5;
complex<T> d12 = (s15 - s23)*t14*t62*T(3); d12 = T(1)/d12;
complex<T> d21 = t14*t4; d21 = T(1)/d21;
complex<T> d23 = t14*t4*t6; d23 = T(1)/d23;
complex<T> d25 = spb12*t4; d25 = T(1)/d25;
complex<T> d27 = spb12*t4*t6; d27 = T(1)/d27;
complex<T> d29 = spb12*t4*t5; d29 = T(1)/d29;
complex<T> t1 = d6*t29*t38 + d7*t102*t39 + d5*t57 + d3*t25*t61 + d4*t25*t61 - d16*spa14*spa24*spb35*t65 + d11*t25*t5*t71 + d9*t25*t5*t71 + d10*t63*t71 + d8*t63*t71 + d2*t24*t82 - d12*t33*cube(spb35) + t52*t56*t58*T(2) + d16*spb13*t51*t52*t87*T(2) + d16*t53*t59*t95*T(2); 
complex<T> t36 = d23*s15; 
complex<T> t49 = d4*spb25*t24*t25 + d10*t40*t63 + d11*t25*t34*t71 - t52*t56*t58*T(2); 
complex<T> t72 = t24*t60; 
complex<T> t2 = d3*spb25*t24*t25 + d19*t102*t39 + d5*t57 + d8*t40*t63 + d16*spa14*spa24*spb35*t65 + d9*t25*t34*t71 + d2*t72 + d18*t77 + d6*t77 + d17*t24*t82 + d12*t33*cube(spb35) + t53*t67*t69*T(2) + d7*t77*T(2) - d16*spb13*t51*t52*t87*T(2) - d16*t53*t59*t95*T(2); 
complex<T> t20 = d26*s23*t61 + d25*spa23*t72 + d27*spa23*spb34*spb45*t38*t73 + spb45*t37*t5*t83; 
complex<T> t21 = d22*spa15*t61 + d21*s15*t72 + spb34*spb45*t36*t38*t73 + d24*spa15*spb45*t5*t83; 
complex<T> t48 = d18*t29*t38 + d17*t72 - t53*t67*t69*T(2) + d19*t77*T(2); 
complex<T> t76 = s45*(d21*t72 + d23*t77*T(2)); 
complex<T> t81 = s45*(d29*s15*t72 + t36*t77); 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t1*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t49*(*CI_users[2]->get_value(mc,ind,mu)) + t48*(*CI_users[3]->get_value(mc,ind,mu)) + t21*(*CI_users[4]->get_value(mc,ind,mu)) + t20*(*CI_users[5]->get_value(mc,ind,mu)) + t70*(*CI_users[6]->get_value(mc,ind,mu)) + t76*(*CI_users[7]->get_value(mc,ind,mu)) + t81*(*CI_users[8]->get_value(mc,ind,mu)) + t19*(*CI_users[9]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  

C5g_mmmpp_nf_wCI::\
C5g_mmmpp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C5g_mmmpp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{m, m, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C5g :  mmmpp nf");
#endif
 
//#define TimeStamp "Tue 19 Oct 2010 20:11:03 on n45"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spb12 = SPB(1,2);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb15 = SPB(1,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s34 = S(3,4);
complex<T> t4 = square(spb45); 
complex<T> t7 = spb12*T(3); 
complex<T> t11 = cube(spb45); 
complex<T> t12 = spa23*spb25; 
complex<T> t9 = spa12*spb15*spb24 + spb34*t12; 
complex<T> d1 = spb23*t7*cube(s15 - s34); d1 = T(1)/d1;
complex<T> d2 = (s15 - s34)*spb15*spb23*spb34*t7; d2 = T(1)/d2;
complex<T> d3 = spb15*spb23*spb34*t7; d3 = T(1)/d3;
complex<T> t15 = d1*t9; 
complex<T> t2 = d3*t11 - spa12*spb24*t12*t15 - d2*t4*t9; 
complex<T> t3 = d3*t11 + spa12*spb24*t12*t15 + d2*t4*t9; 
complex<T> co1 = Complex(0,1); 
SeriesC<T> result = co1*(t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_mmppp_G C5g_28_G
#define _C_mpmpp_G C5g_26_G
#define _C_mppmp_G C5g_22_G
#define _C_mpppm_G C5g_14_G
#define _C_pmmpp_G C5g_25_G
#define _C_pmpmp_G C5g_21_G
#define _C_pmppm_G C5g_13_G
#define _C_ppmmp_G C5g_19_G
#define _C_ppmpm_G C5g_11_G
#define _C_pppmm_G C5g_7_G
#define _C_mmppp_nf C5g_28_nf
#define _C_mpmpp_nf C5g_26_nf
#define _C_mppmp_nf C5g_22_nf
#define _C_mpppm_nf C5g_14_nf
#define _C_pmmpp_nf C5g_25_nf
#define _C_pmpmp_nf C5g_21_nf
#define _C_pmppm_nf C5g_13_nf
#define _C_ppmmp_nf C5g_19_nf
#define _C_ppmpm_nf C5g_11_nf
#define _C_pppmm_nf C5g_7_nf
#define _C_ppmmm_G C5g_3_G
#define _C_pmpmm_G C5g_5_G
#define _C_pmmpm_G C5g_9_G
#define _C_pmmmp_G C5g_17_G
#define _C_mppmm_G C5g_6_G
#define _C_mpmpm_G C5g_10_G
#define _C_mpmmp_G C5g_18_G
#define _C_mmppm_G C5g_12_G
#define _C_mmpmp_G C5g_20_G
#define _C_mmmpp_G C5g_24_G
#define _C_ppmmm_nf C5g_3_nf
#define _C_pmpmm_nf C5g_5_nf
#define _C_pmmpm_nf C5g_9_nf
#define _C_pmmmp_nf C5g_17_nf
#define _C_mppmm_nf C5g_6_nf
#define _C_mpmpm_nf C5g_10_nf
#define _C_mpmmp_nf C5g_18_nf
#define _C_mmppm_nf C5g_12_nf
#define _C_mmpmp_nf C5g_20_nf
#define _C_mmmpp_nf C5g_24_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_mmppp_G case 28
 
#define _CASE_mpmpp_G case 26
 
#define _CASE_mppmp_G case 22
 
#define _CASE_mpppm_G case 14
 
#define _CASE_pmmpp_G case 25
 
#define _CASE_pmpmp_G case 21
 
#define _CASE_pmppm_G case 13
 
#define _CASE_ppmmp_G case 19
 
#define _CASE_ppmpm_G case 11
 
#define _CASE_pppmm_G case 7
 
#define _CASE_mmppp_nf case 28
 
#define _CASE_mpmpp_nf case 26
 
#define _CASE_mppmp_nf case 22
 
#define _CASE_mpppm_nf case 14
 
#define _CASE_pmmpp_nf case 25
 
#define _CASE_pmpmp_nf case 21
 
#define _CASE_pmppm_nf case 13
 
#define _CASE_ppmmp_nf case 19
 
#define _CASE_ppmpm_nf case 11
 
#define _CASE_pppmm_nf case 7
 
#define _CASE_ppmmm_G case 3
 
#define _CASE_pmpmm_G case 5
 
#define _CASE_pmmpm_G case 9
 
#define _CASE_pmmmp_G case 17
 
#define _CASE_mppmm_G case 6
 
#define _CASE_mpmpm_G case 10
 
#define _CASE_mpmmp_G case 18
 
#define _CASE_mmppm_G case 12
 
#define _CASE_mmpmp_G case 20
 
#define _CASE_mmmpp_G case 24
 
#define _CASE_ppmmm_nf case 3
 
#define _CASE_pmpmm_nf case 5
 
#define _CASE_pmmpm_nf case 9
 
#define _CASE_pmmmp_nf case 17
 
#define _CASE_mppmm_nf case 6
 
#define _CASE_mpmpm_nf case 10
 
#define _CASE_mpmmp_nf case 18
 
#define _CASE_mmppm_nf case 12
 
#define _CASE_mmpmp_nf case 20
 
#define _CASE_mmmpp_nf case 24
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_5g_G( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_mmppp_G: return new 
                       C5g_mmppp_G_wCI(ind);
    _CASE_mpmpp_G: return new 
                       C5g_mpmpp_G_wCI(ind);
    _CASE_mppmp_G: return new 
                       C5g_mppmp_G_wCI(ind);
    _CASE_mpppm_G: return new 
                       C5g_mpppm_G_wCI(ind);
    _CASE_pmmpp_G: return new 
                       C5g_pmmpp_G_wCI(ind);
    _CASE_pmpmp_G: return new 
                       C5g_pmpmp_G_wCI(ind);
    _CASE_pmppm_G: return new 
                       C5g_pmppm_G_wCI(ind);
    _CASE_ppmmp_G: return new 
                       C5g_ppmmp_G_wCI(ind);
    _CASE_ppmpm_G: return new 
                       C5g_ppmpm_G_wCI(ind);
    _CASE_pppmm_G: return new 
                       C5g_pppmm_G_wCI(ind);
    _CASE_ppmmm_G: return new 
                       C5g_ppmmm_G_wCI(ind);
    _CASE_pmpmm_G: return new 
                       C5g_pmpmm_G_wCI(ind);
    _CASE_pmmpm_G: return new 
                       C5g_pmmpm_G_wCI(ind);
    _CASE_pmmmp_G: return new 
                       C5g_pmmmp_G_wCI(ind);
    _CASE_mppmm_G: return new 
                       C5g_mppmm_G_wCI(ind);
    _CASE_mpmpm_G: return new 
                       C5g_mpmpm_G_wCI(ind);
    _CASE_mpmmp_G: return new 
                       C5g_mpmmp_G_wCI(ind);
    _CASE_mmppm_G: return new 
                       C5g_mmppm_G_wCI(ind);
    _CASE_mmpmp_G: return new 
                       C5g_mmpmp_G_wCI(ind);
    _CASE_mmmpp_G: return new 
                       C5g_mmmpp_G_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_5g_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_mmppp_nf: return new 
                       C5g_mmppp_nf_wCI(ind);
    _CASE_mpmpp_nf: return new 
                       C5g_mpmpp_nf_wCI(ind);
    _CASE_mppmp_nf: return new 
                       C5g_mppmp_nf_wCI(ind);
    _CASE_mpppm_nf: return new 
                       C5g_mpppm_nf_wCI(ind);
    _CASE_pmmpp_nf: return new 
                       C5g_pmmpp_nf_wCI(ind);
    _CASE_pmpmp_nf: return new 
                       C5g_pmpmp_nf_wCI(ind);
    _CASE_pmppm_nf: return new 
                       C5g_pmppm_nf_wCI(ind);
    _CASE_ppmmp_nf: return new 
                       C5g_ppmmp_nf_wCI(ind);
    _CASE_ppmpm_nf: return new 
                       C5g_ppmpm_nf_wCI(ind);
    _CASE_pppmm_nf: return new 
                       C5g_pppmm_nf_wCI(ind);
    _CASE_ppmmm_nf: return new 
                       C5g_ppmmm_nf_wCI(ind);
    _CASE_pmpmm_nf: return new 
                       C5g_pmpmm_nf_wCI(ind);
    _CASE_pmmpm_nf: return new 
                       C5g_pmmpm_nf_wCI(ind);
    _CASE_pmmmp_nf: return new 
                       C5g_pmmmp_nf_wCI(ind);
    _CASE_mppmm_nf: return new 
                       C5g_mppmm_nf_wCI(ind);
    _CASE_mpmpm_nf: return new 
                       C5g_mpmpm_nf_wCI(ind);
    _CASE_mpmmp_nf: return new 
                       C5g_mpmmp_nf_wCI(ind);
    _CASE_mmppm_nf: return new 
                       C5g_mmppm_nf_wCI(ind);
    _CASE_mmpmp_nf: return new 
                       C5g_mmpmp_nf_wCI(ind);
    _CASE_mmmpp_nf: return new 
                       C5g_mmmpp_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
