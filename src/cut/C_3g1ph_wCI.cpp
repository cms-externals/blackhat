/*
*C_3g1ph_wCI.cpp
*
* Created on 12/9, 2010
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
 
#define mH2 mc.s(ind[2-1],ind[3-1],ind[4-1])
#define SPA(i,j) mc.spa(ind[i-1],ind[j-1])
#define SPB(i,j) mc.spb(ind[i-1],ind[j-1])
#define S(i,j) mc.s(ind[i-1],ind[j-1])
#define SS(i,j,k) mc.s(ind[i-1],ind[j-1],ind[k-1])

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}
 


 
class C3g1ph_phppp_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phppp_G_wCI
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

 
class C3g1ph_phppp_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phppp_nf_wCI
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

 
class C3g1ph_phmpp_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phmpp_G_wCI
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

 
class C3g1ph_phmpp_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phmpp_nf_wCI
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

 
class C3g1ph_phpmp_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phpmp_G_wCI
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

 
class C3g1ph_phpmp_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phpmp_nf_wCI
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

 
class C3g1ph_phppm_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phppm_G_wCI
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

 
class C3g1ph_phppm_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phppm_nf_wCI
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

 
class C3g1ph_phmmp_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phmmp_G_wCI
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

 
class C3g1ph_phmmp_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phmmp_nf_wCI
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

 
class C3g1ph_phmpm_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phmpm_G_wCI
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

 
class C3g1ph_phmpm_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phmpm_nf_wCI
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

 
class C3g1ph_phpmm_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phpmm_G_wCI
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

 
class C3g1ph_phpmm_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phpmm_nf_wCI
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

 
class C3g1ph_phmmm_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phmmm_G_wCI
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

 
class C3g1ph_phmmm_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phmmm_nf_wCI
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

 
class C3g1ph_phdppp_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdppp_G_wCI
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

 
class C3g1ph_phdppp_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdppp_nf_wCI
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

 
class C3g1ph_phdmpp_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdmpp_G_wCI
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

 
class C3g1ph_phdmpp_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdmpp_nf_wCI
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

 
class C3g1ph_phdpmp_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdpmp_G_wCI
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

 
class C3g1ph_phdpmp_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdpmp_nf_wCI
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

 
class C3g1ph_phdppm_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdppm_G_wCI
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

 
class C3g1ph_phdppm_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdppm_nf_wCI
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

 
class C3g1ph_phdmmp_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdmmp_G_wCI
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

 
class C3g1ph_phdmmp_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdmmp_nf_wCI
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

 
class C3g1ph_phdmpm_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdmpm_G_wCI
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

 
class C3g1ph_phdmpm_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdmpm_nf_wCI
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

 
class C3g1ph_phdpmm_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdpmm_G_wCI
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

 
class C3g1ph_phdpmm_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdpmm_nf_wCI
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

 
class C3g1ph_phdmmm_G_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdmmm_G_wCI
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

 
class C3g1ph_phdmmm_nf_wCI : public Cut_Part_wCI {
public:
       C3g1ph_phdmmm_nf_wCI
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



C3g1ph_phppp_G_wCI::\
C3g1ph_phppp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phppp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phppp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phppp_nf_wCI::\
C3g1ph_phppp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phppp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phppp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phmpp_G_wCI::\
C3g1ph_phmpp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phmpp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phmpp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phmpp_nf_wCI::\
C3g1ph_phmpp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phmpp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phmpp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phpmp_G_wCI::\
C3g1ph_phpmp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phpmp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phpmp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phpmp_nf_wCI::\
C3g1ph_phpmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phpmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phpmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phppm_G_wCI::\
C3g1ph_phppm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phppm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phppm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phppm_nf_wCI::\
C3g1ph_phppm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phppm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phppm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C3g1ph_phmmp_G_wCI::\
C3g1ph_phmmp_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c4));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c4, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C3g1ph_phmmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C3g1ph :  phmmp G");
#endif
 
//#define TimeStamp "Thu 9 Dec 2010 10:13:26 on n39"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> s23 = S(2,3);
complex<T> t2 = cube(spa23); 
complex<T> d1 = spa24*spa34; d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> d3 = spa24; d3 = T(1)/d3;
complex<T> d4 = spa24*T(2); d4 = T(1)/d4;
complex<T> d5 = T(2); d5 = T(1)/d5;
complex<T> d6 = spa34*T(2); d6 = T(1)/d6;
complex<T> t4 = s23*t2; 
complex<T> t1 = d1*(-(mH2*t2) + t4)*T(3); 
complex<T> co1 = -(d1*t4*T(2)); 
complex<T> co2 = -(d2*spb24*t2); 
complex<T> co3 = -(d3*spb34*t2); 
complex<T> co4 = -(d4*spb34*t4); 
complex<T> co5 = d5*spb24*spb34*t2; 
complex<T> co6 = -(d6*spb24*t4); 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)) + co4*(*CI_users[4]->get_value(mc,ind,mu)) + co5*(*CI_users[5]->get_value(mc,ind,mu)) + co6*(*CI_users[6]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C3g1ph_phmmp_nf_wCI::\
C3g1ph_phmmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phmmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phmmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C3g1ph_phmpm_G_wCI::\
C3g1ph_phmpm_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c4));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c4, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C3g1ph_phmpm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C3g1ph :  phmpm G");
#endif
 
//#define TimeStamp "Thu 9 Dec 2010 10:13:26 on n39"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> s24 = S(2,4);
complex<T> t2 = cube(spa24); 
complex<T> d1 = spa23*spa34; d1 = T(1)/d1;
complex<T> d2 = spa34; d2 = T(1)/d2;
complex<T> d3 = spa23; d3 = T(1)/d3;
complex<T> d4 = T(2); d4 = T(1)/d4;
complex<T> d5 = spa23*T(2); d5 = T(1)/d5;
complex<T> d6 = spa34*T(2); d6 = T(1)/d6;
complex<T> t3 = d2*spb23; 
complex<T> t1 = -(t2*(d1*mH2 + t3)*T(3)); 
complex<T> co1 = t2*t3*T(2); 
complex<T> co2 = d1*s24*t2; 
complex<T> co3 = -(d3*spb34*t2); 
complex<T> co4 = d4*spb23*spb34*t2; 
complex<T> co5 = -(d5*s24*spb34*t2); 
complex<T> co6 = -(d6*s24*spb23*t2); 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)) + co4*(*CI_users[4]->get_value(mc,ind,mu)) + co5*(*CI_users[5]->get_value(mc,ind,mu)) + co6*(*CI_users[6]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C3g1ph_phmpm_nf_wCI::\
C3g1ph_phmpm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phmpm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phmpm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C3g1ph_phpmm_G_wCI::\
C3g1ph_phpmm_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c4));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c4, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C3g1ph_phpmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C3g1ph :  phpmm G");
#endif
 
//#define TimeStamp "Thu 9 Dec 2010 10:13:27 on n39"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> s34 = S(3,4);
complex<T> t2 = cube(spa34); 
complex<T> d1 = spa23*spa24; d1 = T(1)/d1;
complex<T> d2 = spa24; d2 = T(1)/d2;
complex<T> d3 = spa23; d3 = T(1)/d3;
complex<T> d4 = spa24*T(2); d4 = T(1)/d4;
complex<T> d5 = spa23*T(2); d5 = T(1)/d5;
complex<T> d6 = T(2); d6 = T(1)/d6;
complex<T> t3 = d2*spb23; 
complex<T> t1 = -(t2*(d1*mH2 + t3)*T(3)); 
complex<T> co1 = t2*t3*T(2); 
complex<T> co2 = -(d3*spb24*t2); 
complex<T> co3 = d1*s34*t2; 
complex<T> co4 = -(d4*s34*spb23*t2); 
complex<T> co5 = -(d5*s34*spb24*t2); 
complex<T> co6 = d6*spb23*spb24*t2; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)) + co4*(*CI_users[4]->get_value(mc,ind,mu)) + co5*(*CI_users[5]->get_value(mc,ind,mu)) + co6*(*CI_users[6]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C3g1ph_phpmm_nf_wCI::\
C3g1ph_phpmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phpmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phpmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C3g1ph_phmmm_G_wCI::\
C3g1ph_phmmm_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c4));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c4, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C3g1ph_phmmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C3g1ph :  phmmm G");
#endif
 
//#define TimeStamp "Thu 9 Dec 2010 10:13:27 on n39"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> t2 = square(mH2); 
complex<T> d1 = spb24*spb34; d1 = T(1)/d1;
complex<T> d2 = spb23*spb24*spb34; d2 = T(1)/d2;
complex<T> d3 = spb23*spb34; d3 = T(1)/d3;
complex<T> d4 = spb23*spb24; d4 = T(1)/d4;
complex<T> d5 = spb24*T(2); d5 = T(1)/d5;
complex<T> d6 = spb23*T(2); d6 = T(1)/d6;
complex<T> d7 = spb34*T(2); d7 = T(1)/d7;
complex<T> t3 = d1*spa23; 
complex<T> t5 = t2*t3; 
complex<T> t1 = (t5 + d2*cube(mH2))*T(3); 
complex<T> co1 = -(t5*T(2)); 
complex<T> co2 = d3*spa24*t2; 
complex<T> co3 = d4*spa34*t2; 
complex<T> co4 = -(d5*spa23*spa34*t2); 
complex<T> co5 = -(d6*spa24*spa34*t2); 
complex<T> co6 = -(d7*spa23*spa24*t2); 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)) + co4*(*CI_users[4]->get_value(mc,ind,mu)) + co5*(*CI_users[5]->get_value(mc,ind,mu)) + co6*(*CI_users[6]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C3g1ph_phmmm_nf_wCI::\
C3g1ph_phmmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phmmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phmmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C3g1ph_phdppp_G_wCI::\
C3g1ph_phdppp_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c4));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c4, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C3g1ph_phdppp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdppp G");
#endif
 
//#define TimeStamp "Thu 9 Dec 2010 10:13:28 on n39"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> t2 = square(mH2); 
complex<T> d1 = spa23*spa24*spa34; d1 = T(1)/d1;
complex<T> d2 = spa24*spa34; d2 = T(1)/d2;
complex<T> d3 = spa23*spa34; d3 = T(1)/d3;
complex<T> d4 = spa23*spa24; d4 = T(1)/d4;
complex<T> d5 = spa24*T(2); d5 = T(1)/d5;
complex<T> d6 = spa23*T(2); d6 = T(1)/d6;
complex<T> d7 = spa34*T(2); d7 = T(1)/d7;
complex<T> t3 = d2*spb23; 
complex<T> t5 = t2*t3; 
complex<T> t1 = -((t5 + d1*cube(mH2))*T(3)); 
complex<T> co1 = t5*T(2); 
complex<T> co2 = -(d3*spb24*t2); 
complex<T> co3 = -(d4*spb34*t2); 
complex<T> co4 = d5*spb23*spb34*t2; 
complex<T> co5 = d6*spb24*spb34*t2; 
complex<T> co6 = d7*spb23*spb24*t2; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)) + co4*(*CI_users[4]->get_value(mc,ind,mu)) + co5*(*CI_users[5]->get_value(mc,ind,mu)) + co6*(*CI_users[6]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C3g1ph_phdppp_nf_wCI::\
C3g1ph_phdppp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdppp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdppp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C3g1ph_phdmpp_G_wCI::\
C3g1ph_phdmpp_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c4));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c4, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C3g1ph_phdmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdmpp G");
#endif
 
//#define TimeStamp "Thu 9 Dec 2010 10:13:29 on n39"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> s34 = S(3,4);
complex<T> t2 = cube(spb34); 
complex<T> d1 = spb24; d1 = T(1)/d1;
complex<T> d2 = spb23*spb24; d2 = T(1)/d2;
complex<T> d3 = spb23; d3 = T(1)/d3;
complex<T> d4 = spb24*T(2); d4 = T(1)/d4;
complex<T> d5 = spb23*T(2); d5 = T(1)/d5;
complex<T> d6 = T(2); d6 = T(1)/d6;
complex<T> t3 = d1*spa23; 
complex<T> t1 = t2*(d2*mH2 + t3)*T(3); 
complex<T> co1 = -(t2*t3*T(2)); 
complex<T> co2 = d3*spa24*t2; 
complex<T> co3 = -(d2*s34*t2); 
complex<T> co4 = d4*s34*spa23*t2; 
complex<T> co5 = d5*s34*spa24*t2; 
complex<T> co6 = -(d6*spa23*spa24*t2); 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)) + co4*(*CI_users[4]->get_value(mc,ind,mu)) + co5*(*CI_users[5]->get_value(mc,ind,mu)) + co6*(*CI_users[6]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C3g1ph_phdmpp_nf_wCI::\
C3g1ph_phdmpp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdmpp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdmpp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C3g1ph_phdpmp_G_wCI::\
C3g1ph_phdpmp_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c4));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c4, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C3g1ph_phdpmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdpmp G");
#endif
 
//#define TimeStamp "Thu 9 Dec 2010 10:13:29 on n39"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> s24 = S(2,4);
complex<T> t2 = cube(spb24); 
complex<T> d1 = spb34; d1 = T(1)/d1;
complex<T> d2 = spb23*spb34; d2 = T(1)/d2;
complex<T> d3 = spb23; d3 = T(1)/d3;
complex<T> d4 = T(2); d4 = T(1)/d4;
complex<T> d5 = spb23*T(2); d5 = T(1)/d5;
complex<T> d6 = spb34*T(2); d6 = T(1)/d6;
complex<T> t3 = d1*spa23; 
complex<T> t1 = t2*(d2*mH2 + t3)*T(3); 
complex<T> co1 = -(t2*t3*T(2)); 
complex<T> co2 = -(d2*s24*t2); 
complex<T> co3 = d3*spa34*t2; 
complex<T> co4 = -(d4*spa23*spa34*t2); 
complex<T> co5 = d5*s24*spa34*t2; 
complex<T> co6 = d6*s24*spa23*t2; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)) + co4*(*CI_users[4]->get_value(mc,ind,mu)) + co5*(*CI_users[5]->get_value(mc,ind,mu)) + co6*(*CI_users[6]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C3g1ph_phdpmp_nf_wCI::\
C3g1ph_phdpmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdpmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdpmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C3g1ph_phdppm_G_wCI::\
C3g1ph_phdppm_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c4));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c14));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c12));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c4, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c4, c1));
} 
  
  
template <class T> SeriesC<T> 
     C3g1ph_phdppm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdppm G");
#endif
 
//#define TimeStamp "Thu 9 Dec 2010 10:13:30 on n39"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> s23 = S(2,3);
complex<T> t2 = cube(spb23); 
complex<T> d1 = spb24*spb34; d1 = T(1)/d1;
complex<T> d2 = spb34; d2 = T(1)/d2;
complex<T> d3 = spb24; d3 = T(1)/d3;
complex<T> d4 = spb24*T(2); d4 = T(1)/d4;
complex<T> d5 = T(2); d5 = T(1)/d5;
complex<T> d6 = spb34*T(2); d6 = T(1)/d6;
complex<T> t4 = s23*t2; 
complex<T> t1 = d1*(mH2*t2 - t4)*T(3); 
complex<T> co1 = d1*t4*T(2); 
complex<T> co2 = d2*spa24*t2; 
complex<T> co3 = d3*spa34*t2; 
complex<T> co4 = d4*spa34*t4; 
complex<T> co5 = -(d5*spa24*spa34*t2); 
complex<T> co6 = d6*spa24*t4; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + co2*(*CI_users[2]->get_value(mc,ind,mu)) + co3*(*CI_users[3]->get_value(mc,ind,mu)) + co4*(*CI_users[4]->get_value(mc,ind,mu)) + co5*(*CI_users[5]->get_value(mc,ind,mu)) + co6*(*CI_users[6]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C3g1ph_phdppm_nf_wCI::\
C3g1ph_phdppm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdppm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdppm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phdmmp_G_wCI::\
C3g1ph_phdmmp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdmmp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdmmp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phdmmp_nf_wCI::\
C3g1ph_phdmmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdmmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdmmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phdmpm_G_wCI::\
C3g1ph_phdmpm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdmpm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdmpm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phdmpm_nf_wCI::\
C3g1ph_phdmpm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdmpm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdmpm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phdpmm_G_wCI::\
C3g1ph_phdpmm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdpmm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdpmm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phdpmm_nf_wCI::\
C3g1ph_phdpmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdpmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdpmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phdmmm_G_wCI::\
C3g1ph_phdmmm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdmmm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdmmm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C3g1ph_phdmmm_nf_wCI::\
C3g1ph_phdmmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C3g1ph_phdmmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C3g1ph :  phdmmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
 // *************** table of switch values ************* 
 
#define _C_phppp_G C3g1ph_253_G
#define _C_phppp_nf C3g1ph_253_nf
#define _C_phmpp_G C3g1ph_241_G
#define _C_phmpp_nf C3g1ph_241_nf
#define _C_phpmp_G C3g1ph_205_G
#define _C_phpmp_nf C3g1ph_205_nf
#define _C_phppm_G C3g1ph_61_G
#define _C_phppm_nf C3g1ph_61_nf
#define _C_phmmp_G C3g1ph_193_G
#define _C_phmmp_nf C3g1ph_193_nf
#define _C_phmpm_G C3g1ph_49_G
#define _C_phmpm_nf C3g1ph_49_nf
#define _C_phpmm_G C3g1ph_13_G
#define _C_phpmm_nf C3g1ph_13_nf
#define _C_phmmm_G C3g1ph_1_G
#define _C_phmmm_nf C3g1ph_1_nf
#define _C_phdppp_G C3g1ph_254_G
#define _C_phdppp_nf C3g1ph_254_nf
#define _C_phdmpp_G C3g1ph_242_G
#define _C_phdmpp_nf C3g1ph_242_nf
#define _C_phdpmp_G C3g1ph_206_G
#define _C_phdpmp_nf C3g1ph_206_nf
#define _C_phdppm_G C3g1ph_62_G
#define _C_phdppm_nf C3g1ph_62_nf
#define _C_phdmmp_G C3g1ph_194_G
#define _C_phdmmp_nf C3g1ph_194_nf
#define _C_phdmpm_G C3g1ph_50_G
#define _C_phdmpm_nf C3g1ph_50_nf
#define _C_phdpmm_G C3g1ph_14_G
#define _C_phdpmm_nf C3g1ph_14_nf
#define _C_phdmmm_G C3g1ph_2_G
#define _C_phdmmm_nf C3g1ph_2_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_phppp_G case 253
 
#define _CASE_phppp_nf case 253
 
#define _CASE_phmpp_G case 241
 
#define _CASE_phmpp_nf case 241
 
#define _CASE_phpmp_G case 205
 
#define _CASE_phpmp_nf case 205
 
#define _CASE_phppm_G case 61
 
#define _CASE_phppm_nf case 61
 
#define _CASE_phmmp_G case 193
 
#define _CASE_phmmp_nf case 193
 
#define _CASE_phmpm_G case 49
 
#define _CASE_phmpm_nf case 49
 
#define _CASE_phpmm_G case 13
 
#define _CASE_phpmm_nf case 13
 
#define _CASE_phmmm_G case 1
 
#define _CASE_phmmm_nf case 1
 
#define _CASE_phdppp_G case 254
 
#define _CASE_phdppp_nf case 254
 
#define _CASE_phdmpp_G case 242
 
#define _CASE_phdmpp_nf case 242
 
#define _CASE_phdpmp_G case 206
 
#define _CASE_phdpmp_nf case 206
 
#define _CASE_phdppm_G case 62
 
#define _CASE_phdppm_nf case 62
 
#define _CASE_phdmmp_G case 194
 
#define _CASE_phdmmp_nf case 194
 
#define _CASE_phdmpm_G case 50
 
#define _CASE_phdmpm_nf case 50
 
#define _CASE_phdpmm_G case 14
 
#define _CASE_phdpmm_nf case 14
 
#define _CASE_phdmmm_G case 2
 
#define _CASE_phdmmm_nf case 2
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_3g1ph_G( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_phppp_G: return new 
                       C3g1ph_phppp_G_wCI(ind);
    _CASE_phmpp_G: return new 
                       C3g1ph_phmpp_G_wCI(ind);
    _CASE_phpmp_G: return new 
                       C3g1ph_phpmp_G_wCI(ind);
    _CASE_phppm_G: return new 
                       C3g1ph_phppm_G_wCI(ind);
    _CASE_phmmp_G: return new 
                       C3g1ph_phmmp_G_wCI(ind);
    _CASE_phmpm_G: return new 
                       C3g1ph_phmpm_G_wCI(ind);
    _CASE_phpmm_G: return new 
                       C3g1ph_phpmm_G_wCI(ind);
    _CASE_phmmm_G: return new 
                       C3g1ph_phmmm_G_wCI(ind);
    _CASE_phdppp_G: return new 
                       C3g1ph_phdppp_G_wCI(ind);
    _CASE_phdmpp_G: return new 
                       C3g1ph_phdmpp_G_wCI(ind);
    _CASE_phdpmp_G: return new 
                       C3g1ph_phdpmp_G_wCI(ind);
    _CASE_phdppm_G: return new 
                       C3g1ph_phdppm_G_wCI(ind);
    _CASE_phdmmp_G: return new 
                       C3g1ph_phdmmp_G_wCI(ind);
    _CASE_phdmpm_G: return new 
                       C3g1ph_phdmpm_G_wCI(ind);
    _CASE_phdpmm_G: return new 
                       C3g1ph_phdpmm_G_wCI(ind);
    _CASE_phdmmm_G: return new 
                       C3g1ph_phdmmm_G_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_3g1ph_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_phppp_nf: return new 
                       C3g1ph_phppp_nf_wCI(ind);
    _CASE_phmpp_nf: return new 
                       C3g1ph_phmpp_nf_wCI(ind);
    _CASE_phpmp_nf: return new 
                       C3g1ph_phpmp_nf_wCI(ind);
    _CASE_phppm_nf: return new 
                       C3g1ph_phppm_nf_wCI(ind);
    _CASE_phmmp_nf: return new 
                       C3g1ph_phmmp_nf_wCI(ind);
    _CASE_phmpm_nf: return new 
                       C3g1ph_phmpm_nf_wCI(ind);
    _CASE_phpmm_nf: return new 
                       C3g1ph_phpmm_nf_wCI(ind);
    _CASE_phmmm_nf: return new 
                       C3g1ph_phmmm_nf_wCI(ind);
    _CASE_phdppp_nf: return new 
                       C3g1ph_phdppp_nf_wCI(ind);
    _CASE_phdmpp_nf: return new 
                       C3g1ph_phdmpp_nf_wCI(ind);
    _CASE_phdpmp_nf: return new 
                       C3g1ph_phdpmp_nf_wCI(ind);
    _CASE_phdppm_nf: return new 
                       C3g1ph_phdppm_nf_wCI(ind);
    _CASE_phdmmp_nf: return new 
                       C3g1ph_phdmmp_nf_wCI(ind);
    _CASE_phdmpm_nf: return new 
                       C3g1ph_phdmpm_nf_wCI(ind);
    _CASE_phdpmm_nf: return new 
                       C3g1ph_phdpmm_nf_wCI(ind);
    _CASE_phdmmm_nf: return new 
                       C3g1ph_phdmmm_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
