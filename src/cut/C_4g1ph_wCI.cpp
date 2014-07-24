/*
*C_4g1ph_wCI.cpp
*
* Created on 12/8, 2010
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
 


template<class T> complex<T> coeff3mass(momentum_configuration<T>& mc,
							const Index_Vector& ind, int i1, int i2, int i3, int i4){
    complex<T> K1sq = mH2;
    complex<T> K2sq = S(i1,i2);
    complex<T> K1dK2 = -S(i1,i2) - (S(i1,i3)+S(i1,i4)+S(i2,i3)+S(i2,i4))/T(2);
    complex<T> gammaP = K1dK2 + sqrt(pow(K1dK2,2)-K1sq*K2sq);
    complex<T> gammaM = K1dK2 - sqrt(pow(K1dK2,2)-K1sq*K2sq);
    
    complex<T> a1 = -gammaP-mH2;
    complex<T> a2 = a1;
    complex<T> a3 = -gammaP;
    complex<T> a4 = a3;

    complex<T> resultP = pow(mH2,2)*pow(SPA(i3,i4),3)*
(a3*SPA(i2,i3)*SPB(i3,i1)+a4*SPA(i2,i4)*SPB(i4,i1))*
(a1*SPA(i2,i1)*SPB(i1,i3)+a4*SPA(i2,i4)*SPB(i4,i3))*
(a1*SPA(i2,i1)*SPB(i1,i4)+a3*SPA(i2,i3)*SPB(i3,i4))/
(
(a2*S(i1,i2)+a3*S(i1,i3)+a4*S(i1,i4))*
(a1*S(i1,i3)+a2*S(i2,i3)+a4*S(i3,i4))*
(a1*S(i1,i4)+a2*S(i2,i4)+a3*S(i3,i4))
)/gammaP/(gammaP+mH2)/SPA(i1,i2);

    a1 = -gammaM-mH2;
    a2 = a1;
    a3 = -gammaM;
    a4 = a3;

    complex<T> resultM = pow(mH2,2)*pow(SPA(i3,i4),3)*
(a3*SPA(i2,i3)*SPB(i3,i1)+a4*SPA(i2,i4)*SPB(i4,i1))*
(a1*SPA(i2,i1)*SPB(i1,i3)+a4*SPA(i2,i4)*SPB(i4,i3))*
(a1*SPA(i2,i1)*SPB(i1,i4)+a3*SPA(i2,i3)*SPB(i3,i4))/
(
(a2*S(i1,i2)+a3*S(i1,i3)+a4*S(i1,i4))*
(a1*S(i1,i3)+a2*S(i2,i3)+a4*S(i3,i4))*
(a1*S(i1,i4)+a2*S(i2,i4)+a3*S(i3,i4))
)/gammaM/(gammaM+mH2)/SPA(i1,i2);

    return(resultP+resultM);}


  template complex<R> coeff3mass(momentum_configuration<R>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);
  template complex<RHP> coeff3mass(momentum_configuration<RHP>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);
  template  complex<RVHP> coeff3mass(momentum_configuration<RVHP>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);


 
class C4g1ph_phpppp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpppp_G_wCI
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

 
class C4g1ph_phpppp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpppp_nf_wCI
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

 
class C4g1ph_phmppp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmppp_G_wCI
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

 
class C4g1ph_phmppp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmppp_nf_wCI
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

 
class C4g1ph_phpmpp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpmpp_G_wCI
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

 
class C4g1ph_phpmpp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpmpp_nf_wCI
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

 
class C4g1ph_phppmp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phppmp_G_wCI
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

 
class C4g1ph_phppmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phppmp_nf_wCI
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

 
class C4g1ph_phpppm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpppm_G_wCI
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

 
class C4g1ph_phpppm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpppm_nf_wCI
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

 
class C4g1ph_phmmpp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmmpp_G_wCI
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

 
class C4g1ph_phmmpp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmmpp_nf_wCI
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

 
class C4g1ph_phpmmp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpmmp_G_wCI
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

 
class C4g1ph_phpmmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpmmp_nf_wCI
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

 
class C4g1ph_phppmm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phppmm_G_wCI
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

 
class C4g1ph_phppmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phppmm_nf_wCI
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

 
class C4g1ph_phmppm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmppm_G_wCI
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

 
class C4g1ph_phmppm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmppm_nf_wCI
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

 
class C4g1ph_phmpmp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmpmp_G_wCI
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

 
class C4g1ph_phmpmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmpmp_nf_wCI
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

 
class C4g1ph_phpmpm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpmpm_G_wCI
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

 
class C4g1ph_phpmpm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpmpm_nf_wCI
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

 
class C4g1ph_phpmmm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpmmm_G_wCI
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

 
class C4g1ph_phpmmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phpmmm_nf_wCI
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

 
class C4g1ph_phmpmm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmpmm_G_wCI
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

 
class C4g1ph_phmpmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmpmm_nf_wCI
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

 
class C4g1ph_phmmpm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmmpm_G_wCI
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

 
class C4g1ph_phmmpm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmmpm_nf_wCI
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

 
class C4g1ph_phmmmp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmmmp_G_wCI
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

 
class C4g1ph_phmmmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmmmp_nf_wCI
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

 
class C4g1ph_phmmmm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmmmm_G_wCI
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

 
class C4g1ph_phmmmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phmmmm_nf_wCI
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

 
class C4g1ph_phdpppp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpppp_G_wCI
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

 
class C4g1ph_phdpppp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpppp_nf_wCI
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

 
class C4g1ph_phdmppp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmppp_G_wCI
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

 
class C4g1ph_phdmppp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmppp_nf_wCI
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

 
class C4g1ph_phdpmpp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpmpp_G_wCI
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

 
class C4g1ph_phdpmpp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpmpp_nf_wCI
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

 
class C4g1ph_phdppmp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdppmp_G_wCI
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

 
class C4g1ph_phdppmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdppmp_nf_wCI
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

 
class C4g1ph_phdpppm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpppm_G_wCI
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

 
class C4g1ph_phdpppm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpppm_nf_wCI
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

 
class C4g1ph_phdmmpp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmmpp_G_wCI
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

 
class C4g1ph_phdmmpp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmmpp_nf_wCI
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

 
class C4g1ph_phdpmmp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpmmp_G_wCI
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

 
class C4g1ph_phdpmmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpmmp_nf_wCI
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

 
class C4g1ph_phdppmm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdppmm_G_wCI
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

 
class C4g1ph_phdppmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdppmm_nf_wCI
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

 
class C4g1ph_phdmppm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmppm_G_wCI
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

 
class C4g1ph_phdmppm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmppm_nf_wCI
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

 
class C4g1ph_phdmpmp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmpmp_G_wCI
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

 
class C4g1ph_phdmpmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmpmp_nf_wCI
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

 
class C4g1ph_phdpmpm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpmpm_G_wCI
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

 
class C4g1ph_phdpmpm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpmpm_nf_wCI
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

 
class C4g1ph_phdpmmm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpmmm_G_wCI
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

 
class C4g1ph_phdpmmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdpmmm_nf_wCI
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

 
class C4g1ph_phdmpmm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmpmm_G_wCI
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

 
class C4g1ph_phdmpmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmpmm_nf_wCI
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

 
class C4g1ph_phdmmpm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmmpm_G_wCI
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

 
class C4g1ph_phdmmpm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmmpm_nf_wCI
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

 
class C4g1ph_phdmmmp_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmmmp_G_wCI
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

 
class C4g1ph_phdmmmp_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmmmp_nf_wCI
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

 
class C4g1ph_phdmmmm_G_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmmmm_G_wCI
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

 
class C4g1ph_phdmmmm_nf_wCI : public Cut_Part_wCI {
public:
       C4g1ph_phdmmmm_nf_wCI
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



C4g1ph_phpppp_G_wCI::\
C4g1ph_phpppp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phpppp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpppp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phpppp_nf_wCI::\
C4g1ph_phpppp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phpppp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpppp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phmppp_G_wCI::\
C4g1ph_phmppp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phmppp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmppp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phmppp_nf_wCI::\
C4g1ph_phmppp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phmppp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmppp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phpmpp_G_wCI::\
C4g1ph_phpmpp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phpmpp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpmpp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phpmpp_nf_wCI::\
C4g1ph_phpmpp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phpmpp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpmpp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phppmp_G_wCI::\
C4g1ph_phppmp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phppmp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phppmp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phppmp_nf_wCI::\
C4g1ph_phppmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phppmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phppmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phpppm_G_wCI::\
C4g1ph_phpppm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phpppm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpppm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phpppm_nf_wCI::\
C4g1ph_phpppm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phpppm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpppm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C4g1ph_phmmpp_G_wCI::\
C4g1ph_phmmpp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmmpp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:52:52 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s23 = S(2,3);
complex<T> s14 = S(1,4);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = S(3,5);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spa23); 
complex<T> t11 = square(spb45); 
complex<T> t12 = -(spa25*spa34); 
complex<T> t23 = cube(spa23); 
complex<T> t24 = cube(spb45); 
complex<T> t25 = spa24*spa35; 
complex<T> d1 = (s24 + s45)*spa45*T(3); d1 = T(1)/d1;
complex<T> d2 = spa45*cube(s24 + s45)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s35 + s45)*spa45*T(3); d3 = T(1)/d3;
complex<T> d4 = spa45*cube(s35 + s45)*T(3); d4 = T(1)/d4;
complex<T> d5 = spa25*spa34*spa45; d5 = T(1)/d5;
complex<T> d6 = spa25*spa34; d6 = T(1)/d6;
complex<T> d7 = spa34*spa45; d7 = T(1)/d7;
complex<T> d8 = spa25*spa34*spa45*T(2); d8 = T(1)/d8;
complex<T> d9 = spa25*spa34*T(2); d9 = T(1)/d9;
complex<T> d10 = spa34*spa45*T(2); d10 = T(1)/d10;
complex<T> d11 = spa25*spa45*T(2); d11 = T(1)/d11;
complex<T> d12 = spa34*T(2); d12 = T(1)/d12;
complex<T> d13 = spa25*T(2); d13 = T(1)/d13;
complex<T> t3 = square(s24*spa23 - spb45*t25); 
complex<T> t4 = square(s35*spa23 - spb45*t25); 
complex<T> t14 = -(d5*T(3)); 
complex<T> t16 = d8*s234; 
complex<T> t19 = d10*spb25; 
complex<T> t20 = d11*spb34; 
complex<T> t29 = d4*s35; 
complex<T> t32 = d5*T(3); 
complex<T> t33 = d8*s245; 
complex<T> t34 = d7*spb25; 
complex<T> t37 = spb45*t23; 
complex<T> t49 = spa23*t11; 
complex<T> t56 = spa25*t24; 
complex<T> t58 = d5*t23; 
complex<T> t67 = s24*t23; 
complex<T> t7 = (-mH2 + s15)*t58*T(4); 
complex<T> t9 = d2*(spb45*t3 + s24*t12*t49 + spa34*t25*t56) + d1*spb45*t10*T(11); 
complex<T> t42 = -(d8*mH2*s23*t23) + s235*t16*t23; 
complex<T> t43 = d4*spb45*t4 + t12*t29*t49 + d4*spa34*t25*t56 + d3*spb45*t10*T(11); 
complex<T> t48 = (s345*t16 + mH2*t20)*t23; 
complex<T> t51 = s15*t14*t23 + t32*t67; 
complex<T> t54 = t23*(mH2*t19 + s235*t33); 
complex<T> t60 = s345*t23*t33 + d9*mH2*t37; 
complex<T> t61 = d6*t37; 
complex<T> t69 = t23*t34; 
complex<T> t70 = spa25*t49; 
complex<T> t1 = d2*(t12*t24*t25 - spb45*t3 + s24*spa34*t70) - d1*spb45*t10*T(11); 
complex<T> t35 = s12*t58 + t61; 
complex<T> t36 = d4*t12*t24*t25 - d4*spb45*t4 + spa34*t29*t70 - d3*spb45*t10*T(11); 
complex<T> t57 = s14*t58 + t69; 
complex<T> t62 = s13*t58 + t61; 
complex<T> co1 = t14*t67; 
complex<T> co2 = -t69; 
complex<T> co3 = -(t61*T(2)); 
complex<T> co4 = d12*spb25*t37; 
complex<T> co5 = -(s23*t19*t23); 
complex<T> co6 = -(s23*t20*t23); 
complex<T> co7 = d13*spb34*t37; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t36*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t62*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t57*(*CI_users[8]->get_value(mc,ind,mu)) + t51*(*CI_users[9]->get_value(mc,ind,mu)) + t35*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t60*(*CI_users[12]->get_value(mc,ind,mu)) + t54*(*CI_users[13]->get_value(mc,ind,mu)) + t42*(*CI_users[14]->get_value(mc,ind,mu)) + t48*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmmpp_nf_wCI::\
C4g1ph_phmmpp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmmpp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmmpp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:52:55 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> s24 = S(2,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = S(3,5);
complex<T> t9 = square(spb45); 
complex<T> t10 = square(spa23); 
complex<T> t18 = cube(spb45); 
complex<T> t19 = spa24*spa35; 
complex<T> t20 = spa25*spa34; 
complex<T> d1 = (s24 + s45)*spa45*T(9); d1 = T(1)/d1;
complex<T> d2 = spa45*cube(s24 + s45)*T(9); d2 = T(1)/d2;
complex<T> d3 = (s35 + s45)*spa45*T(9); d3 = T(1)/d3;
complex<T> d4 = spa45*cube(s35 + s45)*T(9); d4 = T(1)/d4;
complex<T> t3 = square(s24*spa23 - spb45*t19); 
complex<T> t4 = square(s35*spa23 - spb45*t19); 
complex<T> t11 = -t20; 
complex<T> t24 = d4*s35; 
complex<T> t25 = -(t10*T(2)); 
complex<T> t27 = t18*t19; 
complex<T> t31 = spa23*t9; 
complex<T> t1 = d1*spb45*t25 + d2*(t11*t27 - spb45*t3 + s24*t20*t31); 
complex<T> t8 = d2*(t20*t27 + spb45*t3 + s24*t11*t31) + d1*spb45*t10*T(2); 
complex<T> t26 = d4*t20*t27 + t11*t24*t31 + d4*spb45*t4 + d3*spb45*t10*T(2); 
complex<T> t29 = -(d4*t4); 
complex<T> t30 = d3*spb45*t25 + d4*t11*t27 + spb45*t29 + t20*t24*t31; 
SeriesC<T> result = t8*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t30*(*CI_users[2]->get_value(mc,ind,mu)) + t26*(*CI_users[3]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phpmmp_G_wCI::\
C4g1ph_phpmmp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phpmmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpmmp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:53:03 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s34 = S(3,4);
complex<T> s14 = S(1,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = S(3,5);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spa34); 
complex<T> t11 = square(spb25); 
complex<T> t12 = -(spa23*spa45); 
complex<T> t23 = cube(spa34); 
complex<T> t24 = cube(spb25); 
complex<T> t25 = spa24*spa35; 
complex<T> d1 = (s25 + s35)*spa25*T(3); d1 = T(1)/d1;
complex<T> d2 = spa25*cube(s25 + s35)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s24 + s25)*spa25*T(3); d3 = T(1)/d3;
complex<T> d4 = spa25*cube(s24 + s25)*T(3); d4 = T(1)/d4;
complex<T> d5 = spa23*spa25*spa45; d5 = T(1)/d5;
complex<T> d6 = spa23*spa25; d6 = T(1)/d6;
complex<T> d7 = spa23*spa45; d7 = T(1)/d7;
complex<T> d8 = spa23*spa25*spa45*T(2); d8 = T(1)/d8;
complex<T> d9 = spa23*spa25*T(2); d9 = T(1)/d9;
complex<T> d10 = spa23*spa45*T(2); d10 = T(1)/d10;
complex<T> d11 = spa25*spa45*T(2); d11 = T(1)/d11;
complex<T> d12 = spa23*T(2); d12 = T(1)/d12;
complex<T> d13 = spa45*T(2); d13 = T(1)/d13;
complex<T> t2 = -(d2*spb25); 
complex<T> t3 = square(s24*spa34 - spb25*t25); 
complex<T> t4 = square(s35*spa34 - spb25*t25); 
complex<T> t14 = -(d5*T(3)); 
complex<T> t15 = d6*spb45; 
complex<T> t16 = -(d1*T(11)); 
complex<T> t17 = d8*s234; 
complex<T> t20 = d11*spb23; 
complex<T> t29 = d2*s35; 
complex<T> t31 = d1*T(11); 
complex<T> t32 = d5*T(3); 
complex<T> t33 = d8*s245; 
complex<T> t34 = d9*spb45; 
complex<T> t37 = d5*t23; 
complex<T> t38 = t24*t25; 
complex<T> t43 = spb25*t23; 
complex<T> t49 = spa34*t11; 
complex<T> t62 = s24*t23; 
complex<T> t8 = (-mH2 + s15)*t37*T(4); 
complex<T> t9 = d4*(spb25*t3 + spa23*spa45*t38 + s24*t12*t49) + d3*spb25*t10*T(11); 
complex<T> t36 = spb25*t10*t31 + d2*spa23*spa45*t38 + d2*spb25*t4 + t12*t29*t49; 
complex<T> t41 = -(d8*mH2*s34*t23) + s345*t17*t23; 
complex<T> t46 = (s235*t17 + mH2*t20)*t23; 
complex<T> t48 = s15*t14*t23 + t32*t62; 
complex<T> t50 = t23*(s345*t33 + mH2*t34); 
complex<T> t57 = t15*t23; 
complex<T> t59 = s235*t23*t33 + d10*mH2*t43; 
complex<T> t66 = d7*t43; 
complex<T> t67 = spa23*t49; 
complex<T> t1 = d4*(-(spb25*t3) + t12*t38 + s24*spa45*t67) - d3*spb25*t10*T(11); 
complex<T> t35 = s14*t37 + t66; 
complex<T> t42 = spb25*t10*t16 + d2*t12*t38 + t2*t4 + spa45*t29*t67; 
complex<T> t53 = s12*t37 + t57; 
complex<T> t56 = s13*t37 + t57; 
complex<T> co1 = t14*t62; 
complex<T> co2 = -t66; 
complex<T> co3 = -(t57*T(2)); 
complex<T> co4 = d12*spb45*t43; 
complex<T> co5 = d13*spb23*t43; 
complex<T> co6 = -(s34*t20*t23); 
complex<T> co7 = -(s34*t23*t34); 
SeriesC<T> result = t36*(*CI_users[0]->get_value(mc,ind,mu)) + t42*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t56*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t35*(*CI_users[8]->get_value(mc,ind,mu)) + t48*(*CI_users[9]->get_value(mc,ind,mu)) + t53*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t50*(*CI_users[12]->get_value(mc,ind,mu)) + t59*(*CI_users[13]->get_value(mc,ind,mu)) + t46*(*CI_users[14]->get_value(mc,ind,mu)) + t41*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phpmmp_nf_wCI::\
C4g1ph_phpmmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phpmmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpmmp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:53:06 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> s24 = S(2,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = S(3,5);
complex<T> t9 = square(spb25); 
complex<T> t10 = square(spa34); 
complex<T> t18 = cube(spb25); 
complex<T> t19 = spa24*spa35; 
complex<T> t20 = spa23*spa45; 
complex<T> d1 = (s25 + s35)*spa25*T(9); d1 = T(1)/d1;
complex<T> d2 = spa25*cube(s25 + s35)*T(9); d2 = T(1)/d2;
complex<T> d3 = (s24 + s25)*spa25*T(9); d3 = T(1)/d3;
complex<T> d4 = spa25*cube(s24 + s25)*T(9); d4 = T(1)/d4;
complex<T> t2 = -(d2*spb25); 
complex<T> t3 = square(s24*spa34 - spb25*t19); 
complex<T> t4 = square(s35*spa34 - spb25*t19); 
complex<T> t11 = -t20; 
complex<T> t13 = -(d1*T(2)); 
complex<T> t24 = d2*s35; 
complex<T> t27 = t18*t19; 
complex<T> t31 = spa34*t9; 
complex<T> t1 = d4*(t11*t27 - spb25*t3 + s24*t20*t31) - d3*spb25*t10*T(2); 
complex<T> t8 = d4*(t20*t27 + spb25*t3 + s24*t11*t31) + d3*spb25*t10*T(2); 
complex<T> t26 = d2*t20*t27 + t11*t24*t31 + d2*spb25*t4 + d1*spb25*t10*T(2); 
complex<T> t30 = spb25*t10*t13 + d2*t11*t27 + t20*t24*t31 + t2*t4; 
SeriesC<T> result = t30*(*CI_users[0]->get_value(mc,ind,mu)) + t26*(*CI_users[1]->get_value(mc,ind,mu)) + t8*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phppmm_G_wCI::\
C4g1ph_phppmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phppmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phppmm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:53:12 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> s24 = S(2,4);
complex<T> s45 = S(4,5);
complex<T> s15 = S(1,5);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s35 = S(3,5);
complex<T> s14 = S(1,4);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spa45); 
complex<T> t11 = square(spb23); 
complex<T> t23 = cube(spa45); 
complex<T> t24 = spa25*spa34; 
complex<T> t25 = cube(spb23); 
complex<T> t26 = spa24*spa35; 
complex<T> d1 = (s23 + s24)*spa23*T(3); d1 = T(1)/d1;
complex<T> d2 = spa23*cube(s23 + s24)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s23 + s35)*spa23*T(3); d3 = T(1)/d3;
complex<T> d4 = spa23*cube(s23 + s35)*T(3); d4 = T(1)/d4;
complex<T> d6 = spa23*spa34; d6 = T(1)/d6;
complex<T> d8 = spa23*spa34*T(2); d8 = T(1)/d8;
complex<T> d10 = spa23*spa25*T(2); d10 = T(1)/d10;
complex<T> d11 = spa34*T(2); d11 = T(1)/d11;
complex<T> d12 = spa25*T(2); d12 = T(1)/d12;
complex<T> t3 = square(s24*spa45 - spb23*t26); 
complex<T> t4 = square(s35*spa45 - spb23*t26); 
complex<T> t12 = -t24; 
complex<T> t21 = d6*spb25; 
complex<T> t22 = d10*spb34; 
complex<T> t30 = d4*s35; 
complex<T> t35 = d8*spb25; 
complex<T> t39 = t25*t26; 
complex<T> t48 = spa45*t11; 
complex<T> d5 = spa23*t24; d5 = T(1)/d5;
complex<T> d7 = spa23*t24*T(2); d7 = T(1)/d7;
complex<T> d9 = t24*T(2); d9 = T(1)/d9;
complex<T> t1 = d2*(-(spb23*t3) + t12*t39 + s24*t24*t48) - d1*spb23*t10*T(11); 
complex<T> t9 = d2*(spb23*t3 + t24*t39 + s24*t12*t48) + d1*spb23*t10*T(11); 
complex<T> t14 = -(d5*T(3)); 
complex<T> t16 = d7*s234; 
complex<T> t33 = -(d5*s45); 
complex<T> t34 = d7*s245; 
complex<T> t37 = d4*t12*t39 - d4*spb23*t4 + t24*t30*t48 - d3*spb23*t10*T(11); 
complex<T> t38 = d5*t23; 
complex<T> t43 = d4*t24*t39 + d4*spb23*t4 + t12*t30*t48 + d3*spb23*t10*T(11); 
complex<T> t63 = t21*t23; 
complex<T> t8 = (-mH2 + s15)*t38*T(4); 
complex<T> t42 = (d9*mH2*spb23 + s235*t16)*t23; 
complex<T> t47 = (s345*t16 + mH2*t22)*t23; 
complex<T> t49 = -(d7*mH2*s45*t23) + s345*t23*t34; 
complex<T> t52 = t23*(s235*t34 + mH2*t35); 
complex<T> t55 = t23*t33 + s12*t38; 
complex<T> t58 = t23*t33 + s13*t38; 
complex<T> t61 = s14*t38 + t63; 
complex<T> t62 = t14*t23; 
complex<T> t36 = s15*t62 + s24*t38*T(3); 
complex<T> co1 = s24*t62; 
complex<T> co2 = -t63; 
complex<T> co3 = s45*t38*T(2); 
complex<T> co4 = -(s45*t23*t35); 
complex<T> co5 = d11*spb23*spb25*t23; 
complex<T> co6 = d12*spb23*spb34*t23; 
complex<T> co7 = -(s45*t22*t23); 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t37*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t58*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t61*(*CI_users[8]->get_value(mc,ind,mu)) + t36*(*CI_users[9]->get_value(mc,ind,mu)) + t55*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t49*(*CI_users[12]->get_value(mc,ind,mu)) + t52*(*CI_users[13]->get_value(mc,ind,mu)) + t42*(*CI_users[14]->get_value(mc,ind,mu)) + t47*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phppmm_nf_wCI::\
C4g1ph_phppmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phppmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phppmm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:53:15 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = S(2,4);
complex<T> s35 = S(3,5);
complex<T> t9 = square(spb23); 
complex<T> t11 = square(spa45); 
complex<T> t18 = cube(spb23); 
complex<T> t19 = spa24*spa35; 
complex<T> t20 = spa25*spa34; 
complex<T> d1 = (s23 + s24)*spa23*T(9); d1 = T(1)/d1;
complex<T> d2 = spa23*cube(s23 + s24)*T(9); d2 = T(1)/d2;
complex<T> d3 = (s23 + s35)*spa23*T(9); d3 = T(1)/d3;
complex<T> d4 = spa23*cube(s23 + s35)*T(9); d4 = T(1)/d4;
complex<T> t3 = square(s24*spa45 - spb23*t19); 
complex<T> t4 = square(s35*spa45 - spb23*t19); 
complex<T> t10 = -t20; 
complex<T> t24 = d4*s35; 
complex<T> t25 = -(t11*T(2)); 
complex<T> t27 = t18*t19; 
complex<T> t31 = spa45*t9; 
complex<T> t1 = d1*spb23*t25 + d2*(t10*t27 - spb23*t3 + s24*t20*t31); 
complex<T> t8 = d2*(t20*t27 + spb23*t3 + s24*t10*t31) + d1*spb23*t11*T(2); 
complex<T> t26 = d4*t20*t27 + t10*t24*t31 + d4*spb23*t4 + d3*spb23*t11*T(2); 
complex<T> t29 = -(d4*t4); 
complex<T> t30 = d3*spb23*t25 + d4*t10*t27 + spb23*t29 + t20*t24*t31; 
SeriesC<T> result = t8*(*CI_users[0]->get_value(mc,ind,mu)) + t26*(*CI_users[1]->get_value(mc,ind,mu)) + t30*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmppm_G_wCI::\
C4g1ph_phmppm_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmppm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmppm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:53:22 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s14 = S(1,4);
complex<T> s25 = S(2,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = S(3,5);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spa25); 
complex<T> t11 = square(spb34); 
complex<T> t23 = cube(spa25); 
complex<T> t24 = spa23*spa45; 
complex<T> t25 = cube(spb34); 
complex<T> t26 = spa24*spa35; 
complex<T> d1 = (s24 + s34)*spa34*T(3); d1 = T(1)/d1;
complex<T> d2 = spa34*cube(s24 + s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s34 + s35)*spa34*T(3); d3 = T(1)/d3;
complex<T> d4 = spa34*cube(s34 + s35)*T(3); d4 = T(1)/d4;
complex<T> d6 = spa23*spa34; d6 = T(1)/d6;
complex<T> d8 = spa23*spa34*T(2); d8 = T(1)/d8;
complex<T> d9 = spa34*spa45*T(2); d9 = T(1)/d9;
complex<T> d11 = spa45*T(2); d11 = T(1)/d11;
complex<T> d12 = spa23*T(2); d12 = T(1)/d12;
complex<T> t3 = square(s24*spa25 - spb34*t26); 
complex<T> t4 = square(s35*spa25 - spb34*t26); 
complex<T> t12 = -t24; 
complex<T> t15 = d6*spb45; 
complex<T> t20 = d9*spb23; 
complex<T> t30 = d4*s35; 
complex<T> t35 = d8*spb45; 
complex<T> t39 = t25*t26; 
complex<T> t51 = spa25*t11; 
complex<T> t53 = -(mH2*t23); 
complex<T> t63 = s24*t23; 
complex<T> d5 = spa34*t24; d5 = T(1)/d5;
complex<T> d7 = spa34*t24*T(2); d7 = T(1)/d7;
complex<T> d10 = t24*T(2); d10 = T(1)/d10;
complex<T> t1 = d2*(-(spb34*t3) + t12*t39 + s24*t24*t51) - d1*spb34*t10*T(11); 
complex<T> t9 = d2*(spb34*t3 + t24*t39 + s24*t12*t51) + d1*spb34*t10*T(11); 
complex<T> t14 = -(d5*T(3)); 
complex<T> t17 = d7*s234; 
complex<T> t33 = d5*T(3); 
complex<T> t34 = d7*s245; 
complex<T> t37 = d4*t12*t39 - d4*spb34*t4 + t24*t30*t51 - d3*spb34*t10*T(11); 
complex<T> t43 = d4*t24*t39 + d4*spb34*t4 + t12*t30*t51 + d3*spb34*t10*T(11); 
complex<T> t44 = d5*t23; 
complex<T> t59 = t15*t23; 
complex<T> t7 = (-mH2 + s15)*t44*T(4); 
complex<T> t36 = (s14 - s25)*t44; 
complex<T> t42 = s12*t44 + t59; 
complex<T> t48 = (d10*mH2*spb34 + s345*t17)*t23; 
complex<T> t50 = (s235*t17 + mH2*t20)*t23; 
complex<T> t52 = s15*t14*t23 + t33*t63; 
complex<T> t55 = t23*(s345*t34 + mH2*t35); 
complex<T> t58 = s13*t44 + t59; 
complex<T> t61 = s235*t23*t34 + d7*s25*t53; 
complex<T> co1 = t14*t63; 
complex<T> co2 = s25*t44; 
complex<T> co3 = -(t59*T(2)); 
complex<T> co4 = -(s25*t23*t35); 
complex<T> co5 = -(s25*t20*t23); 
complex<T> co6 = d11*spb23*spb34*t23; 
complex<T> co7 = d12*spb34*spb45*t23; 
SeriesC<T> result = t9*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t37*(*CI_users[2]->get_value(mc,ind,mu)) + t43*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t58*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t36*(*CI_users[8]->get_value(mc,ind,mu)) + t52*(*CI_users[9]->get_value(mc,ind,mu)) + t42*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t55*(*CI_users[12]->get_value(mc,ind,mu)) + t61*(*CI_users[13]->get_value(mc,ind,mu)) + t50*(*CI_users[14]->get_value(mc,ind,mu)) + t48*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmppm_nf_wCI::\
C4g1ph_phmppm_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmppm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmppm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:53:25 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> s24 = S(2,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = S(3,5);
complex<T> t9 = square(spb34); 
complex<T> t10 = square(spa25); 
complex<T> t18 = cube(spb34); 
complex<T> t19 = spa24*spa35; 
complex<T> t20 = spa23*spa45; 
complex<T> d1 = (s24 + s34)*spa34*T(9); d1 = T(1)/d1;
complex<T> d2 = spa34*cube(s24 + s34)*T(9); d2 = T(1)/d2;
complex<T> d3 = (s34 + s35)*spa34*T(9); d3 = T(1)/d3;
complex<T> d4 = spa34*cube(s34 + s35)*T(9); d4 = T(1)/d4;
complex<T> t3 = square(s24*spa25 - spb34*t19); 
complex<T> t4 = square(s35*spa25 - spb34*t19); 
complex<T> t11 = -t20; 
complex<T> t24 = d4*s35; 
complex<T> t25 = -(t10*T(2)); 
complex<T> t27 = t18*t19; 
complex<T> t31 = spa25*t9; 
complex<T> t1 = d1*spb34*t25 + d2*(t11*t27 - spb34*t3 + s24*t20*t31); 
complex<T> t8 = d2*(t20*t27 + spb34*t3 + s24*t11*t31) + d1*spb34*t10*T(2); 
complex<T> t26 = d4*t20*t27 + t11*t24*t31 + d4*spb34*t4 + d3*spb34*t10*T(2); 
complex<T> t29 = -(d4*t4); 
complex<T> t30 = d3*spb34*t25 + d4*t11*t27 + spb34*t29 + t20*t24*t31; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t8*(*CI_users[1]->get_value(mc,ind,mu)) + t26*(*CI_users[2]->get_value(mc,ind,mu)) + t30*(*CI_users[3]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmpmp_G_wCI::\
C4g1ph_phmpmp_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmpmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmpmp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:53:59 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s14 = S(1,4);
complex<T> s35 = -(spa35*spb35);
complex<T> s13 = S(1,3);
complex<T> s12 = S(1,2);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t17 = square(spa24); 
complex<T> t20 = square(spb35); 
complex<T> t33 = square(square(spa24)); 
complex<T> t34 = spa23*spa45; 
complex<T> t35 = spa25*spa34; 
complex<T> t36 = -(spb35*T(2)); 
complex<T> t38 = s23*s25; 
complex<T> t39 = s34*s45; 
complex<T> t43 = cube(spb35); 
complex<T> t44 = square(spa23); 
complex<T> t45 = square(spa25); 
complex<T> t46 = square(spa34); 
complex<T> t47 = square(spa45); 
complex<T> t52 = spb35*T(2); 
complex<T> d1 = -((s25 + s35)*spa35); d1 = T(1)/d1;
complex<T> d2 = -((s25 + s35)*cube(spa35)); d2 = T(1)/d2;
complex<T> d3 = square(s25 + s35)*square(spa35); d3 = T(1)/d3;
complex<T> d4 = -(spa35*cube(s25 + s35)*T(3)); d4 = T(1)/d4;
complex<T> d5 = -((s23 + s35)*spa35); d5 = T(1)/d5;
complex<T> d6 = -((s23 + s35)*cube(spa35)); d6 = T(1)/d6;
complex<T> d7 = square(s23 + s35)*square(spa35); d7 = T(1)/d7;
complex<T> d8 = -(spa35*cube(s23 + s35)*T(3)); d8 = T(1)/d8;
complex<T> d9 = -((s35 + s45)*spa35); d9 = T(1)/d9;
complex<T> d10 = -((s35 + s45)*cube(spa35)); d10 = T(1)/d10;
complex<T> d11 = square(s35 + s45)*square(spa35); d11 = T(1)/d11;
complex<T> d12 = -(spa35*cube(s35 + s45)*T(3)); d12 = T(1)/d12;
complex<T> d13 = -((s34 + s35)*spa35); d13 = T(1)/d13;
complex<T> d14 = -((s34 + s35)*cube(spa35)); d14 = T(1)/d14;
complex<T> d15 = square(s34 + s35)*square(spa35); d15 = T(1)/d15;
complex<T> d16 = -(spa35*cube(s34 + s35)*T(3)); d16 = T(1)/d16;
complex<T> d18 = square(spa35); d18 = T(1)/d18;
complex<T> d19 = square(square(spa35)); d19 = T(1)/d19;
complex<T> d27 = spa23*spa34*T(2); d27 = T(1)/d27;
complex<T> d28 = square(square(spa35))*T(2); d28 = T(1)/d28;
complex<T> d29 = spa34*spa45*T(2); d29 = T(1)/d29;
complex<T> d30 = spa25*spa45*T(2); d30 = T(1)/d30;
complex<T> d31 = spa23*spa25*T(2); d31 = T(1)/d31;
complex<T> t18 = -t35; 
complex<T> t25 = -(d18*T(4)); 
complex<T> t26 = -t38; 
complex<T> t27 = -t39; 
complex<T> t31 = d29*spb23; 
complex<T> t54 = d19*t34; 
complex<T> t55 = t45*t46; 
complex<T> t60 = -(t17*T(4)); 
complex<T> t61 = t20*t34; 
complex<T> t62 = d28*t35; 
complex<T> t66 = -(t35*T(2)); 
complex<T> t67 = spb35*t17; 
complex<T> t68 = t44*t47; 
complex<T> t71 = t35*T(2); 
complex<T> t72 = mH2*t33; 
complex<T> t81 = d18*t17; 
complex<T> t97 = spb34*t33; 
complex<T> d17 = t34*t35; d17 = T(1)/d17;
complex<T> d20 = spa23*t35; d20 = T(1)/d20;
complex<T> d21 = spa34*t34; d21 = T(1)/d21;
complex<T> d24 = spa34*t34*T(2); d24 = T(1)/d24;
complex<T> d26 = spa25*t34*T(2); d26 = T(1)/d26;
complex<T> t28 = -(d17*T(3)); 
complex<T> t29 = d20*spb45; 
complex<T> t32 = d21*spb25; 
complex<T> t41 = spb25*t31; 
complex<T> t53 = d17*t33; 
complex<T> t56 = d18*t26; 
complex<T> t59 = s23*(t17*t25 + t54*t71); 
complex<T> t64 = s34*(t17*t25 + t54*t71); 
complex<T> t73 = t43*t68; 
complex<T> t77 = t43*t55; 
complex<T> t87 = t34*t62; 
complex<T> d22 = t34*t71; d22 = T(1)/d22;
complex<T> d23 = spa23*t71; d23 = T(1)/d23;
complex<T> d25 = spa45*t71; d25 = T(1)/d25;
complex<T> t10 = (-mH2 + s15)*t53*T(4); 
complex<T> t11 = t17*t56 + t38*t87; 
complex<T> t13 = t33*t41 + t17*t56 + t38*t87; 
complex<T> t14 = d10*t34*t35*t36 + d11*t18*t61 + d12*t77*T(2) + d9*t67*T(4); 
complex<T> t16 = t27*t81 + t39*t87 + d31*spb45*t97; 
complex<T> t30 = d22*s234; 
complex<T> t40 = d22*s245; 
complex<T> t42 = d10*t34*t35*t52 + d14*t34*t35*t52 + d13*spb35*t60 + d9*spb35*t60 + d11*t35*t61 + d15*t35*t61 - d16*t73*T(2) - d12*t77*T(2); 
complex<T> t50 = t27*t81 + t39*t87; 
complex<T> t58 = d14*t34*t35*t36 + d15*t18*t61 + d16*t73*T(2) + d13*t67*T(4); 
complex<T> t84 = t29*t33; 
complex<T> t99 = t28*t33; 
complex<T> t100 = d4*t73; 
complex<T> t102 = t32*t33; 
complex<T> t105 = d8*t77; 
complex<T> t12 = t102 + s25*(t17*t25 + t54*t71) + s14*(t53 + t54*t66 + t81*T(4)); 
complex<T> t15 = s45*(t17*t25 + t54*t71) + t84 + s12*(t53 + t54*t66 + t81*T(4)); 
complex<T> t49 = s15*t99 + s24*t53*T(3); 
complex<T> t57 = s235*t30*t33 + d25*spb23*t72; 
complex<T> t63 = s345*t30*t33 + d26*spb34*t72; 
complex<T> t65 = d2*t34*t35*t36 + d3*t18*t61 + t100*T(2) + d1*t67*T(4); 
complex<T> t69 = s235*t33*t40 + d24*spb25*t72; 
complex<T> t70 = d6*t34*t35*t36 + d7*t18*t61 + t105*T(2) + d5*t67*T(4); 
complex<T> t74 = s13*t53 + t84; 
complex<T> t75 = d2*t34*t35*t52 + d6*t34*t35*t52 + d1*spb35*t60 + d5*spb35*t60 + d3*t35*t61 + d7*t35*t61 - t100*T(2) - t105*T(2); 
complex<T> t79 = s345*t33*t40 + d23*spb45*t72; 
complex<T> co1 = s24*t99; 
complex<T> co2 = -t102; 
complex<T> co3 = -(t84*T(2)); 
complex<T> co4 = d27*spb25*spb45*t33; 
complex<T> co5 = d30*spb23*t97; 
SeriesC<T> result = t65*(*CI_users[0]->get_value(mc,ind,mu)) + t75*(*CI_users[1]->get_value(mc,ind,mu)) + t70*(*CI_users[2]->get_value(mc,ind,mu)) + t14*(*CI_users[3]->get_value(mc,ind,mu)) + t42*(*CI_users[4]->get_value(mc,ind,mu)) + t58*(*CI_users[5]->get_value(mc,ind,mu)) + t10*(*CI_users[6]->get_value(mc,ind,mu)) + t59*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t74*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t12*(*CI_users[11]->get_value(mc,ind,mu)) + t49*(*CI_users[12]->get_value(mc,ind,mu)) + t64*(*CI_users[13]->get_value(mc,ind,mu)) + t15*(*CI_users[14]->get_value(mc,ind,mu)) + co3*(*CI_users[15]->get_value(mc,ind,mu)) + t79*(*CI_users[16]->get_value(mc,ind,mu)) + t69*(*CI_users[17]->get_value(mc,ind,mu)) + t57*(*CI_users[18]->get_value(mc,ind,mu)) + t63*(*CI_users[19]->get_value(mc,ind,mu)) + co4*(*CI_users[20]->get_value(mc,ind,mu)) + t13*(*CI_users[21]->get_value(mc,ind,mu)) + t50*(*CI_users[22]->get_value(mc,ind,mu)) + co5*(*CI_users[23]->get_value(mc,ind,mu)) + t11*(*CI_users[24]->get_value(mc,ind,mu)) + t16*(*CI_users[25]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmpmp_nf_wCI::\
C4g1ph_phmpmp_nf_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmpmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmpmp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:54:13 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> s23 = S(2,3);
complex<T> s14 = S(1,4);
complex<T> s25 = S(2,5);
complex<T> s34 = S(3,4);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> s35 = -(spa35*spb35);
complex<T> t6 = cube(spa35); 
complex<T> t14 = square(spa24); 
complex<T> t15 = square(spb35); 
complex<T> t21 = -s23 - s35; 
complex<T> t22 = s34*s45; 
complex<T> t27 = spa23*spa45; 
complex<T> t28 = spa25*spa34; 
complex<T> t29 = cube(spb35); 
complex<T> t36 = square(spa23); 
complex<T> t37 = square(spa25); 
complex<T> t38 = square(spa34); 
complex<T> t39 = square(spa45); 
complex<T> t44 = -(spb35*T(2)); 
complex<T> t54 = spb35*T(2); 
complex<T> d1 = -((s25 + s35)*spa35*T(3)); d1 = T(1)/d1;
complex<T> d3 = square(s25 + s35)*square(spa35)*T(3); d3 = T(1)/d3;
complex<T> d4 = -(spa35*cube(s25 + s35)*T(9)); d4 = T(1)/d4;
complex<T> d9 = -((s35 + s45)*spa35*T(3)); d9 = T(1)/d9;
complex<T> d11 = square(s35 + s45)*square(spa35)*T(3); d11 = T(1)/d11;
complex<T> d12 = -(spa35*cube(s35 + s45)*T(9)); d12 = T(1)/d12;
complex<T> d13 = -((s34 + s35)*spa35*T(3)); d13 = T(1)/d13;
complex<T> d15 = square(s34 + s35)*square(spa35)*T(3); d15 = T(1)/d15;
complex<T> d16 = -(spa35*cube(s34 + s35)*T(9)); d16 = T(1)/d16;
complex<T> d17 = square(spa35)*T(3); d17 = T(1)/d17;
complex<T> d18 = square(square(spa35))*T(3); d18 = T(1)/d18;
complex<T> d19 = square(spa35)*T(12); d19 = T(1)/d19;
complex<T> d20 = square(square(spa35))*T(6); d20 = T(1)/d20;
complex<T> t16 = -t28; 
complex<T> t31 = -(d18*T(2)); 
complex<T> t34 = d19*t22; 
complex<T> t40 = -(t29*T(2)); 
complex<T> t43 = t27*t28; 
complex<T> t45 = d17*t14; 
complex<T> t51 = d18*T(2); 
complex<T> t55 = t37*t38; 
complex<T> t56 = t36*t39; 
complex<T> t58 = spb35*t14; 
complex<T> d2 = -((s25 + s35)*t6*T(3)); d2 = T(1)/d2;
complex<T> d5 = spa35*t21*T(3); d5 = T(1)/d5;
complex<T> d6 = t21*t6*T(3); d6 = T(1)/d6;
complex<T> d7 = square(spa35)*square(t21)*T(3); d7 = T(1)/d7;
complex<T> d8 = spa35*cube(t21)*T(9); d8 = T(1)/d8;
complex<T> d10 = -((s35 + s45)*t6*T(3)); d10 = T(1)/d10;
complex<T> d14 = -((s34 + s35)*t6*T(3)); d14 = T(1)/d14;
complex<T> t12 = s25*t31*t43 - s14*t45 + s25*t45 + s14*t43*t51; 
complex<T> t26 = s45*t31*t43 - s12*t45 + s45*t45 + s12*t43*t51; 
complex<T> t42 = s23*(t31*t43 + t45); 
complex<T> t48 = s34*(t31*t43 + t45); 
complex<T> t49 = -t58; 
complex<T> t50 = t16*t27; 
complex<T> t7 = s23*s25*(d19*t14 + d20*t50); 
complex<T> t13 = d11*t15*t43 + d9*t49 + d10*t43*t54 + d12*t40*t55; 
complex<T> t35 = d7*t15*t43 + d5*t49 + d6*t43*t54 + d8*t40*t55; 
complex<T> t41 = t14*t34 + d20*t22*t50; 
complex<T> t47 = d15*t15*t43 + d13*t49 + d14*t43*t54 + d16*t40*t56; 
complex<T> t53 = d3*t15*t43 + d1*t49 + d2*t43*t54 + d4*t40*t56; 
complex<T> t57 = d10*t43*t44 + d14*t43*t44 + d11*t15*t50 + d15*t15*t50 + d13*t58 + d9*t58 + d12*t29*t55*T(2) + d16*t29*t56*T(2); 
complex<T> t60 = d2*t43*t44 + d6*t43*t44 + d3*t15*t50 + d7*t15*t50 + d1*t58 + d5*t58 + d8*t29*t55*T(2) + d4*t29*t56*T(2); 
SeriesC<T> result = t53*(*CI_users[0]->get_value(mc,ind,mu)) + t60*(*CI_users[1]->get_value(mc,ind,mu)) + t35*(*CI_users[2]->get_value(mc,ind,mu)) + t13*(*CI_users[3]->get_value(mc,ind,mu)) + t57*(*CI_users[4]->get_value(mc,ind,mu)) + t47*(*CI_users[5]->get_value(mc,ind,mu)) + t42*(*CI_users[6]->get_value(mc,ind,mu)) + t12*(*CI_users[7]->get_value(mc,ind,mu)) + t48*(*CI_users[8]->get_value(mc,ind,mu)) + t26*(*CI_users[9]->get_value(mc,ind,mu)) + t7*(*CI_users[10]->get_value(mc,ind,mu)) + t41*(*CI_users[11]->get_value(mc,ind,mu)) + t7*(*CI_users[12]->get_value(mc,ind,mu)) + t41*(*CI_users[13]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phpmpm_G_wCI::\
C4g1ph_phpmpm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phpmpm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpmpm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:54:51 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s24 = -(spa24*spb24);
complex<T> s14 = S(1,4);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t18 = square(spa35); 
complex<T> t21 = square(spb24); 
complex<T> t34 = square(square(spa35)); 
complex<T> t35 = spa23*spa45; 
complex<T> t36 = spa25*spa34; 
complex<T> t37 = -(spb24*T(2)); 
complex<T> t39 = s23*s34; 
complex<T> t40 = s25*s45; 
complex<T> t45 = cube(spb24); 
complex<T> t46 = square(spa23); 
complex<T> t47 = square(spa25); 
complex<T> t48 = square(spa34); 
complex<T> t49 = square(spa45); 
complex<T> t54 = spb24*T(2); 
complex<T> d1 = -((s24 + s34)*spa24); d1 = T(1)/d1;
complex<T> d2 = -((s24 + s34)*cube(spa24)); d2 = T(1)/d2;
complex<T> d3 = square(s24 + s34)*square(spa24); d3 = T(1)/d3;
complex<T> d4 = -(spa24*cube(s24 + s34)*T(3)); d4 = T(1)/d4;
complex<T> d5 = -((s23 + s24)*spa24); d5 = T(1)/d5;
complex<T> d6 = -((s23 + s24)*cube(spa24)); d6 = T(1)/d6;
complex<T> d7 = square(s23 + s24)*square(spa24); d7 = T(1)/d7;
complex<T> d8 = -(spa24*cube(s23 + s24)*T(3)); d8 = T(1)/d8;
complex<T> d9 = -((s24 + s25)*spa24); d9 = T(1)/d9;
complex<T> d10 = -((s24 + s45)*spa24); d10 = T(1)/d10;
complex<T> d11 = -((s24 + s25)*cube(spa24)); d11 = T(1)/d11;
complex<T> d12 = -((s24 + s45)*cube(spa24)); d12 = T(1)/d12;
complex<T> d13 = square(s24 + s25)*square(spa24); d13 = T(1)/d13;
complex<T> d14 = square(s24 + s45)*square(spa24); d14 = T(1)/d14;
complex<T> d15 = -(spa24*cube(s24 + s45)*T(3)); d15 = T(1)/d15;
complex<T> d16 = -(spa24*cube(s24 + s25)*T(3)); d16 = T(1)/d16;
complex<T> d18 = square(spa24); d18 = T(1)/d18;
complex<T> d19 = square(square(spa24)); d19 = T(1)/d19;
complex<T> d20 = spa24; d20 = T(1)/d20;
complex<T> d21 = cube(spa24); d21 = T(1)/d21;
complex<T> d29 = square(square(spa24))*T(2); d29 = T(1)/d29;
complex<T> d30 = spa23*spa34*T(2); d30 = T(1)/d30;
complex<T> d31 = spa34*spa45*T(2); d31 = T(1)/d31;
complex<T> d32 = spa25*spa45*T(2); d32 = T(1)/d32;
complex<T> d33 = spa23*spa25*T(2); d33 = T(1)/d33;
complex<T> t19 = -t36; 
complex<T> t26 = -(d18*T(4)); 
complex<T> t27 = -t39; 
complex<T> t28 = -t40; 
complex<T> t32 = d32*spb23; 
complex<T> t56 = d19*t35; 
complex<T> t57 = t47*t48; 
complex<T> t62 = -(t18*T(4)); 
complex<T> t63 = t21*t35; 
complex<T> t64 = d29*t36; 
complex<T> t68 = -(t36*T(2)); 
complex<T> t69 = spb24*t18; 
complex<T> t70 = t46*t49; 
complex<T> t73 = t36*T(2); 
complex<T> t74 = mH2*t34; 
complex<T> t98 = spb25*t34; 
complex<T> d17 = t35*t36; d17 = T(1)/d17;
complex<T> d22 = spa23*t36; d22 = T(1)/d22;
complex<T> d23 = spa34*t35; d23 = T(1)/d23;
complex<T> d26 = spa34*t35*T(2); d26 = T(1)/d26;
complex<T> d28 = spa25*t35*T(2); d28 = T(1)/d28;
complex<T> t13 = d12*t35*t36*t37 + d14*t19*t63 + d15*t45*t57*T(2) + d10*t69*T(4); 
complex<T> t29 = -(d17*T(3)); 
complex<T> t30 = d22*spb45; 
complex<T> t41 = d23*spb25; 
complex<T> t43 = spb34*t32; 
complex<T> t55 = d17*t34; 
complex<T> t58 = d18*t27; 
complex<T> t61 = s23*(t18*t26 + t56*t73); 
complex<T> t66 = s34*(t18*t26 + t56*t73); 
complex<T> t72 = d6*t35*t36*t37 + d7*t19*t63 + d8*t45*t57*T(2) + d5*t69*T(4); 
complex<T> t75 = t45*t70; 
complex<T> t86 = t35*t64; 
complex<T> d24 = t35*t73; d24 = T(1)/d24;
complex<T> d25 = spa23*t73; d25 = T(1)/d25;
complex<T> d27 = spa45*t73; d27 = T(1)/d27;
complex<T> t9 = (-mH2 + s15)*t55*T(4); 
complex<T> t10 = t18*t58 + t39*t86; 
complex<T> t11 = s15*t29*t34 + d21*t35*t36*t37 + s15*t56*t68 + s24*t55*T(3) + d18*s15*t18*T(4) + d20*t69*T(4); 
complex<T> t12 = s24*t29*t34 + d21*t35*t36*t54 + d20*spb24*t62; 
complex<T> t14 = s25*t18*t26 - t34*t41 + s25*t56*t73; 
complex<T> t15 = t34*t43 + t18*t58 + t39*t86; 
complex<T> t31 = d24*s234; 
complex<T> t42 = d24*s245; 
complex<T> t44 = d11*t35*t36*t54 + d12*t35*t36*t54 + d10*spb24*t62 + d9*spb24*t62 + d13*t36*t63 + d14*t36*t63 - d15*t45*t57*T(2) - d16*t75*T(2); 
complex<T> t52 = d18*t18*t28 + t40*t86; 
complex<T> t60 = d11*t35*t36*t37 + d13*t19*t63 + d16*t75*T(2) + d9*t69*T(4); 
complex<T> t76 = t34*t41 + s14*t55; 
complex<T> t83 = t30*t34; 
complex<T> t101 = d4*t75; 
complex<T> t16 = s45*t18*t26 + s13*t55 + s13*t56*t68 + s45*t56*t73 + t83 + d18*s13*t18*T(4); 
complex<T> t17 = t52 + d30*spb45*t98; 
complex<T> t51 = s12*t55 + t83; 
complex<T> t59 = s235*t31*t34 + d27*spb23*t74; 
complex<T> t65 = s345*t31*t34 + d28*spb34*t74; 
complex<T> t67 = d2*t35*t36*t37 + d3*t19*t63 + t101*T(2) + d1*t69*T(4); 
complex<T> t71 = s235*t34*t42 + d26*spb25*t74; 
complex<T> t77 = d2*t35*t36*t54 + d6*t35*t36*t54 + d1*spb24*t62 + d5*spb24*t62 + d3*t36*t63 + d7*t36*t63 - t101*T(2) - d8*t45*t57*T(2); 
complex<T> t81 = s345*t34*t42 + d25*spb45*t74; 
complex<T> co1 = -(t83*T(2)); 
complex<T> co2 = d31*spb23*t98; 
complex<T> co3 = d33*spb34*spb45*t34; 
SeriesC<T> result = t67*(*CI_users[0]->get_value(mc,ind,mu)) + t77*(*CI_users[1]->get_value(mc,ind,mu)) + t44*(*CI_users[2]->get_value(mc,ind,mu)) + t13*(*CI_users[3]->get_value(mc,ind,mu)) + t72*(*CI_users[4]->get_value(mc,ind,mu)) + t60*(*CI_users[5]->get_value(mc,ind,mu)) + t9*(*CI_users[6]->get_value(mc,ind,mu)) + t61*(*CI_users[7]->get_value(mc,ind,mu)) + t12*(*CI_users[8]->get_value(mc,ind,mu)) + t16*(*CI_users[9]->get_value(mc,ind,mu)) + t14*(*CI_users[10]->get_value(mc,ind,mu)) + t76*(*CI_users[11]->get_value(mc,ind,mu)) + t11*(*CI_users[12]->get_value(mc,ind,mu)) + t66*(*CI_users[13]->get_value(mc,ind,mu)) + t51*(*CI_users[14]->get_value(mc,ind,mu)) + co1*(*CI_users[15]->get_value(mc,ind,mu)) + t81*(*CI_users[16]->get_value(mc,ind,mu)) + t71*(*CI_users[17]->get_value(mc,ind,mu)) + t59*(*CI_users[18]->get_value(mc,ind,mu)) + t65*(*CI_users[19]->get_value(mc,ind,mu)) + t10*(*CI_users[20]->get_value(mc,ind,mu)) + t17*(*CI_users[21]->get_value(mc,ind,mu)) + co2*(*CI_users[22]->get_value(mc,ind,mu)) + t15*(*CI_users[23]->get_value(mc,ind,mu)) + t52*(*CI_users[24]->get_value(mc,ind,mu)) + co3*(*CI_users[25]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phpmpm_nf_wCI::\
C4g1ph_phpmpm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phpmpm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpmpm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:55:08 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> s23 = S(2,3);
complex<T> s25 = S(2,5);
complex<T> s34 = S(3,4);
complex<T> s13 = S(1,3);
complex<T> s45 = S(4,5);
complex<T> s15 = S(1,5);
complex<T> s24 = -(spa24*spb24);
complex<T> t15 = square(spa35); 
complex<T> t16 = square(spb24); 
complex<T> t19 = -s24 - s25; 
complex<T> t20 = -s24 - s34; 
complex<T> t21 = -s24 - s45; 
complex<T> t22 = -s23 - s24; 
complex<T> t23 = s25*s45; 
complex<T> t28 = spa23*spa45; 
complex<T> t29 = spa25*spa34; 
complex<T> t30 = -(spb24*T(2)); 
complex<T> t36 = cube(spb24); 
complex<T> t37 = square(spa23); 
complex<T> t38 = square(spa25); 
complex<T> t39 = square(spa34); 
complex<T> t40 = square(spa45); 
complex<T> d17 = square(spa24)*T(3); d17 = T(1)/d17;
complex<T> d18 = square(square(spa24))*T(3); d18 = T(1)/d18;
complex<T> d19 = spa24*T(3); d19 = T(1)/d19;
complex<T> d20 = cube(spa24)*T(3); d20 = T(1)/d20;
complex<T> d21 = square(spa24)*T(12); d21 = T(1)/d21;
complex<T> d22 = square(square(spa24))*T(6); d22 = T(1)/d22;
complex<T> t17 = -t29; 
complex<T> t34 = d21*t23; 
complex<T> t44 = t28*t29; 
complex<T> t46 = -(d18*T(2)); 
complex<T> t51 = spb24*t15; 
complex<T> t55 = -(t36*T(2)); 
complex<T> t64 = t38*t39; 
complex<T> t65 = t37*t40; 
complex<T> d1 = spa24*t20*T(3); d1 = T(1)/d1;
complex<T> d2 = t20*cube(spa24)*T(3); d2 = T(1)/d2;
complex<T> d3 = square(spa24)*square(t20)*T(3); d3 = T(1)/d3;
complex<T> d4 = spa24*cube(t20)*T(9); d4 = T(1)/d4;
complex<T> d5 = spa24*t22*T(3); d5 = T(1)/d5;
complex<T> d6 = t22*cube(spa24)*T(3); d6 = T(1)/d6;
complex<T> d7 = square(spa24)*square(t22)*T(3); d7 = T(1)/d7;
complex<T> d8 = spa24*cube(t22)*T(9); d8 = T(1)/d8;
complex<T> d9 = spa24*t19*T(3); d9 = T(1)/d9;
complex<T> d10 = spa24*t21*T(3); d10 = T(1)/d10;
complex<T> d11 = t19*cube(spa24)*T(3); d11 = T(1)/d11;
complex<T> d12 = t21*cube(spa24)*T(3); d12 = T(1)/d12;
complex<T> d13 = square(spa24)*square(t19)*T(3); d13 = T(1)/d13;
complex<T> d14 = square(spa24)*square(t21)*T(3); d14 = T(1)/d14;
complex<T> d15 = spa24*cube(t21)*T(9); d15 = T(1)/d15;
complex<T> d16 = spa24*cube(t19)*T(9); d16 = T(1)/d16;
complex<T> t12 = s25*(d17*t15 + t44*t46); 
complex<T> t13 = d20*t30*t44 + d19*t51; 
complex<T> t25 = -(d15*T(2)); 
complex<T> t42 = d15*T(2); 
complex<T> t43 = s23*(d17*t15 + t44*t46); 
complex<T> t45 = -t51; 
complex<T> t49 = s34*(d17*t15 + t44*t46); 
complex<T> t50 = t44*T(2); 
complex<T> t54 = t17*t28; 
complex<T> t7 = s23*s34*(d21*t15 + d22*t54); 
complex<T> t27 = d17*(-s13 + s45)*t15 + s45*t44*t46 + d18*s13*t50; 
complex<T> t41 = t15*t34 + d22*t23*t54; 
complex<T> t48 = d11*t30*t44 + d12*t30*t44 + d10*t51 + d9*t51 + d13*t16*t54 + d14*t16*t54 + t36*t42*t64 + d16*t36*t65*T(2); 
complex<T> t58 = d2*t30*t44 + d6*t30*t44 + d1*t51 + d5*t51 + d3*t16*t54 + d7*t16*t54 + d8*t36*t64*T(2) + d4*t36*t65*T(2); 
complex<T> t61 = spb24*t50; 
complex<T> t14 = -(d17*s15*t15) + d19*t45 + d18*s15*t50 + d20*t61; 
complex<T> t35 = d13*t16*t44 + d9*t45 + d11*t61 + d16*t55*t65; 
complex<T> t53 = d14*t16*t44 + d10*t45 + d12*t61 + t25*t36*t64; 
complex<T> t60 = d7*t16*t44 + d5*t45 + d6*t61 + d8*t55*t64; 
complex<T> t62 = d3*t16*t44 + d1*t45 + d2*t61 + d4*t55*t65; 
SeriesC<T> result = t62*(*CI_users[0]->get_value(mc,ind,mu)) + t58*(*CI_users[1]->get_value(mc,ind,mu)) + t48*(*CI_users[2]->get_value(mc,ind,mu)) + t53*(*CI_users[3]->get_value(mc,ind,mu)) + t60*(*CI_users[4]->get_value(mc,ind,mu)) + t35*(*CI_users[5]->get_value(mc,ind,mu)) + t43*(*CI_users[6]->get_value(mc,ind,mu)) + t13*(*CI_users[7]->get_value(mc,ind,mu)) + t27*(*CI_users[8]->get_value(mc,ind,mu)) + t12*(*CI_users[9]->get_value(mc,ind,mu)) + t14*(*CI_users[10]->get_value(mc,ind,mu)) + t49*(*CI_users[11]->get_value(mc,ind,mu)) + t7*(*CI_users[12]->get_value(mc,ind,mu)) + t41*(*CI_users[13]->get_value(mc,ind,mu)) + t7*(*CI_users[14]->get_value(mc,ind,mu)) + t41*(*CI_users[15]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phpmmm_G_wCI::\
C4g1ph_phpmmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c45));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c25, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phpmmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpmmm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:58:29 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = -(spa34*spb34);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> s14 = S(1,4);
complex<T> s235 = SS(2,3,5);
complex<T> s35 = -(spa35*spb35);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s345 = SS(3,4,5);
complex<T> s12 = S(1,2);
complex<T> s15 = S(1,5);
complex<T> t18 = square(spa34*spb23 - spa45*spb25); 
complex<T> t23 = square(square(spa34*spb23 - spa45*spb25)); 
complex<T> t24 = spb23*spb25; 
complex<T> t27 = s235*T(2); 
complex<T> t42 = square(mH2); 
complex<T> t43 = square(spa45); 
complex<T> t44 = square(spa34); 
complex<T> t45 = square(spb23); 
complex<T> t46 = square(spb25); 
complex<T> t47 = cube(s345); 
complex<T> t48 = square(square(spa35)); 
complex<T> t49 = square(spb24); 
complex<T> t51 = spa23*spb24 + spa35*spb45; 
complex<T> t52 = spa25*spb24 + spa35*spb34; 
complex<T> t53 = spa24*spb34 + spa25*spb35; 
complex<T> t56 = spa23*spb35 + spa24*spb45; 
complex<T> t57 = -(spa24*T(2)); 
complex<T> t58 = s23*s25; 
complex<T> t65 = s235*T(3); 
complex<T> t72 = cube(spa45); 
complex<T> t73 = cube(spa34); 
complex<T> t78 = square(spa24); 
complex<T> t86 = spa23*spa25; 
 complex<T> t99 = coeff3mass<T>(mc,ind,2,3,4,5); 
 complex<T> t100 = coeff3mass<T>(mc,ind,2,5,4,3); 
complex<T> t111 = square(square(s345)); 
complex<T> t112 = spa24*spb23; 
complex<T> t125 = spa34*spa45; 
complex<T> t129 = spa24*spb25; 
complex<T> t189 = spa34*T(2); 
complex<T> t190 = spa35*spb23; 
complex<T> t197 = spa35*spb25; 
complex<T> d1 = spb34*square(s24 + s34)*T(3); d1 = T(1)/d1;
complex<T> d2 = spb34*cube(s24 + s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = s234*(s24 + s34)*spb34*T(3); d3 = T(1)/d3;
complex<T> d4 = spa35*cube(spb35); d4 = T(1)/d4;
complex<T> d5 = (s25 + s35)*spa35*cube(spb35); d5 = T(1)/d5;
complex<T> d6 = s235*spa35*cube(spb35); d6 = T(1)/d6;
complex<T> d7 = square(s25 + s35)*square(spb35); d7 = T(1)/d7;
complex<T> d8 = s235*square(spb35); d8 = T(1)/d8;
complex<T> d9 = spb35*square(s25 + s35)*T(3); d9 = T(1)/d9;
complex<T> d10 = spb35*cube(s25 + s35)*T(3); d10 = T(1)/d10;
complex<T> d12 = (s23 + s35)*spa35*cube(spb35); d12 = T(1)/d12;
complex<T> d13 = square(s23 + s35)*square(spb35); d13 = T(1)/d13;
complex<T> d14 = spb35*cube(s23 + s35)*T(3); d14 = T(1)/d14;
complex<T> d15 = spb35*square(s23 + s35)*T(3); d15 = T(1)/d15;
complex<T> d17 = spb45*cube(s24 + s45)*T(3); d17 = T(1)/d17;
complex<T> d18 = spb45*square(s24 + s45)*T(3); d18 = T(1)/d18;
complex<T> d19 = s245*(s24 + s45)*spb45*T(3); d19 = T(1)/d19;
complex<T> d28 = s235*square(square(spb35)); d28 = T(1)/d28;
complex<T> d29 = s235*cube(spb35); d29 = T(1)/d29;
complex<T> t19 = square(spa34*spb24 + t197); 
complex<T> t20 = square(spa45*spb24 + t190); 
complex<T> t22 = cube(spa34*spb24 + t197); 
complex<T> t26 = cube(spa45*spb24 + t190); 
complex<T> t29 = t51*t52; 
complex<T> t32 = t53*t56; 
complex<T> t54 = spa34*spb35 + t129; 
complex<T> t55 = spa45*spb35 + t112; 
complex<T> t62 = d29*s14; 
complex<T> t64 = d2*s234; 
complex<T> t67 = d17*s245; 
complex<T> t75 = -(t45*T(2)); 
complex<T> t83 = s245*t51; 
complex<T> t85 = s234*t52; 
complex<T> t107 = -t125; 
complex<T> t108 = t42*t48; 
complex<T> t109 = -(t46*T(2)); 
complex<T> t110 = d8*t18; 
complex<T> t138 = d18*spb24; 
complex<T> t144 = d4*T(4); 
complex<T> t146 = t46*t78; 
complex<T> t155 = t42*t72; 
complex<T> t166 = t125*t24; 
complex<T> t169 = s45*t42; 
complex<T> t176 = t42*t73; 
complex<T> t183 = spa45*t46; 
complex<T> t196 = spa45*t44; 
complex<T> t217 = t129*t45; 
complex<T> d11 = (s25 + s35)*spb35*t65; d11 = T(1)/d11;
complex<T> d16 = (s23 + s35)*spb35*t65; d16 = T(1)/d16;
complex<T> d53 = spb25*spb45*t51*T(2); d53 = T(1)/d53;
complex<T> d56 = spb23*spb34*t52*T(2); d56 = T(1)/d56;
complex<T> d58 = t27*square(square(spb35)); d58 = T(1)/d58;
complex<T> d59 = t27*cube(spb35); d59 = T(1)/d59;
complex<T> d60 = t27*square(spb35); d60 = T(1)/d60;
complex<T> t8 = d15*spa35*t107*t24 + d9*spa35*t107*t24 + d5*s235*t109*t43 + d10*s235*spa35*t109*t43 + d13*s235*t44*t45 + d15*spa35*t44*t45 + t144*t44*t45 + d7*s235*t43*t46 + d9*spa35*t43*t46 + t144*t43*t46 + d12*s235*t44*t75 + d14*s235*spa35*t44*t75 + t110*T(2) - d6*s23*t18*T(2) - d6*s25*t18*T(2) - d4*t166*T(8) - d11*spa35*t18*T(11) - d16*spa35*t18*T(11); 
complex<T> t28 = t54*t55; 
complex<T> t74 = -t108; 
complex<T> t89 = t53*t55; 
complex<T> t90 = t54*t56; 
complex<T> t98 = d3*t20; 
complex<T> t139 = d19*t19; 
complex<T> t160 = d28*t146; 
complex<T> t187 = t49*t64; 
complex<T> t194 = t49*t67; 
complex<T> t208 = t112*t183; 
complex<T> d20 = spb23*spb34*t85; d20 = T(1)/d20;
complex<T> d24 = spb34*spb45*t32; d24 = T(1)/d24;
complex<T> d25 = s235*t29*t86; d25 = T(1)/d25;
complex<T> d26 = spb25*spb45*t83; d26 = T(1)/d26;
complex<T> d27 = spb34*t85*T(2); d27 = T(1)/d27;
complex<T> d33 = spb34*spb45*t32*T(2); d33 = T(1)/d33;
complex<T> d34 = spa25*t27*t29; d34 = T(1)/d34;
complex<T> d35 = spb25*spb45*t83*T(2); d35 = T(1)/d35;
complex<T> d36 = spb25*t83; d36 = T(1)/d36;
complex<T> d37 = spb23*spb34*t85*T(2); d37 = T(1)/d37;
complex<T> d41 = spa23*t27*t29; d41 = T(1)/d41;
complex<T> d42 = spb45*t83*T(2); d42 = T(1)/d42;
complex<T> d44 = s235*spa23*t29; d44 = T(1)/d44;
complex<T> d45 = spb23*t85*T(2); d45 = T(1)/d45;
complex<T> d47 = spb45*t32*T(2); d47 = T(1)/d47;
complex<T> d48 = t27*t29*t86; d48 = T(1)/d48;
complex<T> d49 = spb34*t32; d49 = T(1)/d49;
complex<T> d50 = spb34*t32*T(2); d50 = T(1)/d50;
complex<T> d51 = spb25*t83*T(2); d51 = T(1)/d51;
complex<T> d55 = t29*t86*T(2); d55 = T(1)/d55;
complex<T> d62 = s235*t29*T(4); d62 = T(1)/d62;
complex<T> d63 = t32*T(4); d63 = T(1)/d63;
complex<T> d64 = t85*T(2); d64 = T(1)/d64;
complex<T> d65 = t83*T(2); d65 = T(1)/d65;
complex<T> t2 = d1*spb24*t107*t190 + t187*t189*t43 - d1*spa34*t43*t49 + spa34*t98*T(11); 
complex<T> t7 = -t110 + t144*t166 + d1*spb24*t125*t190 - d7*s235*t43*t46 + d5*t27*t43*t46 + d10*spa35*t27*t43*t46 + d9*spa35*(t166 - t43*t46) + d1*spa34*t43*t49 + d6*s25*t18*T(2) - spa34*t187*t43*T(2) - d4*t43*t46*T(4) + d11*spa35*t18*T(11) - spa34*t98*T(11); 
complex<T> t35 = (d24*s12 + d49*spa45)*t47; 
complex<T> t95 = d33*s25; 
complex<T> t119 = d33*s23; 
complex<T> t216 = t194*t196; 
complex<T> d21 = s245*spa25*t89; d21 = T(1)/d21;
complex<T> d22 = s235*t24*t28; d22 = T(1)/d22;
complex<T> d23 = s234*spa23*t90; d23 = T(1)/d23;
complex<T> d30 = s245*spa25*t89*T(2); d30 = T(1)/d30;
complex<T> d31 = spb25*t27*t28; d31 = T(1)/d31;
complex<T> d32 = s234*t90*T(2); d32 = T(1)/d32;
complex<T> d38 = s245*t89*T(2); d38 = T(1)/d38;
complex<T> d39 = spb23*t27*t28; d39 = T(1)/d39;
complex<T> d40 = s234*spa23*t90*T(2); d40 = T(1)/d40;
complex<T> d43 = s235*spb23*t28; d43 = T(1)/d43;
complex<T> d46 = t24*t27*t28; d46 = T(1)/d46;
complex<T> d52 = spa25*t89*T(2); d52 = T(1)/d52;
complex<T> d54 = t24*t28*T(2); d54 = T(1)/d54;
complex<T> d57 = spa23*t90*T(2); d57 = T(1)/d57;
complex<T> d61 = s235*t28*T(4); d61 = T(1)/d61;
complex<T> t1 = d62*t108*t24 + t110*t58 + d60*t166*t58 + d59*t208*t58 + d59*spa34*t217*t58 + d58*t146*t45*t58 + d61*t23*t86; 
complex<T> t3 = s25*(d57*t176 + d56*t26); 
complex<T> t4 = t107*t138*t197 - d18*t196*t49 + t216*T(2) + spa45*t139*T(11); 
complex<T> t5 = s23*(d52*t155 + d53*t22); 
complex<T> t6 = d21*(s13 - s45)*t155 + (d26*s13 + d36*spa45)*t22; 
complex<T> t9 = -t110 + d15*spa35*t166 + t144*t166 + t125*t138*t197 - d13*s235*t44*t45 - d15*spa35*t44*t45 + d12*t27*t44*t45 + d14*spa35*t27*t44*t45 + d18*t196*t49 + d6*s23*t18*T(2) - t216*T(2) - d4*t44*t45*T(4) - spa45*t139*T(11) + d16*spa35*t18*T(11); 
complex<T> t10 = d25*s15*t108 + d21*s15*t155 + d23*s15*t176 - d26*mH2*t22 + d26*s15*t22 - d22*mH2*t23 + d22*s15*t23 - d20*mH2*t26 + d20*s15*t26 - d24*mH2*t47 + d24*s15*t47 - d25*t48*cube(mH2) - d21*t72*cube(mH2) - d23*t73*cube(mH2); 
complex<T> t11 = -(s24*(d21*t155 + d26*t22 + d22*t23 + d24*t47 - d25*t74)); 
complex<T> t12 = d30*s34*t155 + d40*s34*t176 - d35*s34*t22 + d46*s34*t23 + d45*spa34*t26 - d47*spa34*t47 + d48*s34*t74; 
complex<T> t13 = s34*(d55*t108 + d54*t23); 
complex<T> t14 = -(d51*spa45*t22) + d46*s45*t23 - d37*s45*t26 + d40*t169*t73 + d48*s45*t74 + d30*s45*t155*T(3) - d50*spa45*t47*T(3); 
complex<T> t16 = -(d38*spb25*t155) + d40*s25*t176 + d42*spa25*t22 - d37*s25*t26 + d41*spb25*t74 + t47*t95 - d39*spa25*t23*T(3); 
complex<T> t17 = d25*s14*t108 + d44*spb25*t108 + d29*s25*t189*t217 + d22*s14*t23 + d43*spa25*t23 + spb23*t183*t57*t62 + spa34*spb25*t45*t57*t62 + s14*t160*t75 - d8*s14*t166*T(2) + d8*s25*t166*T(2) + d29*s25*t208*T(2) + s25*t160*t45*T(2) - s14*t110*T(4) + s25*t110*T(4); 
complex<T> t37 = s34*(d52*t155 + d53*t22); 
complex<T> t38 = -(d38*s45*spb25*t155) + d65*spa25*spa45*t22; 
complex<T> t39 = -((s15 - s24)*(d21*t155 + d26*t22 + d22*t23 + d24*t47)) + d25*(s24*t108 + s15*t74); 
complex<T> t40 = s45*(d55*t108 + d54*t23); 
complex<T> t106 = d56*s45*t26 + d57*t169*t73; 
complex<T> t137 = -(d32*spb23); 
complex<T> t15 = d34*spb23*t108 + t137*t176 - d31*spa23*t23 + d27*spa23*t26 + t119*t47 + s23*(d30*t155 + d29*t189*t217 - d35*t22 + d8*t166*T(2) + d29*t208*T(2) + t160*t45*T(2) + t110*T(4)); 
complex<T> t41 = s34*t137*t176 + d64*spa23*spa34*t26; 
complex<T> co1 = -t99; 
complex<T> co2 = -t100; 
complex<T> co3 = t111*t119; 
complex<T> co4 = t111*t95; 
complex<T> co5 = d63*t125*t47; 
SeriesC<T> result = t7*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t8*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + t9*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t10*(*CI_users[7]->get_value(mc,ind,mu)) + t15*(*CI_users[8]->get_value(mc,ind,mu)) + t11*(*CI_users[9]->get_value(mc,ind,mu)) + t6*(*CI_users[10]->get_value(mc,ind,mu)) + t16*(*CI_users[11]->get_value(mc,ind,mu)) + t17*(*CI_users[12]->get_value(mc,ind,mu)) + t39*(*CI_users[13]->get_value(mc,ind,mu)) + t12*(*CI_users[14]->get_value(mc,ind,mu)) + t35*(*CI_users[15]->get_value(mc,ind,mu)) + t14*(*CI_users[16]->get_value(mc,ind,mu)) + co3*(*CI_users[17]->get_value(mc,ind,mu)) + co4*(*CI_users[18]->get_value(mc,ind,mu)) + t5*(*CI_users[19]->get_value(mc,ind,mu)) + t37*(*CI_users[20]->get_value(mc,ind,mu)) + t13*(*CI_users[21]->get_value(mc,ind,mu)) + t40*(*CI_users[22]->get_value(mc,ind,mu)) + t3*(*CI_users[23]->get_value(mc,ind,mu)) + t106*(*CI_users[24]->get_value(mc,ind,mu)) + t1*(*CI_users[25]->get_value(mc,ind,mu)) + co5*(*CI_users[26]->get_value(mc,ind,mu)) + t41*(*CI_users[27]->get_value(mc,ind,mu)) + t38*(*CI_users[28]->get_value(mc,ind,mu)) + t1*(*CI_users[29]->get_value(mc,ind,mu)) + co5*(*CI_users[30]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phpmmm_nf_wCI::\
C4g1ph_phpmmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phpmmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, p, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phpmmm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 20:59:12 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s24 = -(spa24*spb24);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s23 = S(2,3);
complex<T> s235 = SS(2,3,5);
complex<T> s14 = S(1,4);
complex<T> s25 = S(2,5);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> s245 = SS(2,4,5);
complex<T> t8 = square(spa34*spb23 - spa45*spb25); 
complex<T> t9 = -(spa34*T(2)); 
complex<T> t10 = -(s23*s25); 
complex<T> t11 = square(spa35*spb23 + spa45*spb24); 
complex<T> t12 = square(spa34*spb24 + spa35*spb25); 
complex<T> t13 = square(spb35); 
complex<T> t20 = square(spb23); 
complex<T> t21 = square(spb25); 
complex<T> t22 = square(spa34); 
complex<T> t23 = square(spa45); 
complex<T> t24 = square(spb24); 
complex<T> t25 = cube(spb35); 
complex<T> t27 = square(spa24); 
complex<T> t47 = spb35*T(9); 
complex<T> t57 = spa45*spb23; 
complex<T> t58 = spa34*spb25; 
complex<T> t106 = spa34*spa35; 
complex<T> t108 = spa35*spa45; 
complex<T> d1 = spb34*square(s24 + s34)*T(9); d1 = T(1)/d1;
complex<T> d2 = spb34*cube(s24 + s34)*T(9); d2 = T(1)/d2;
complex<T> d3 = s234*(s24 + s34)*spb34*T(9); d3 = T(1)/d3;
complex<T> d17 = spb45*cube(s24 + s45)*T(9); d17 = T(1)/d17;
complex<T> d18 = spb45*square(s24 + s45)*T(9); d18 = T(1)/d18;
complex<T> d19 = s245*(s24 + s45)*spb45*T(9); d19 = T(1)/d19;
complex<T> d20 = s235*square(square(spb35))*T(3); d20 = T(1)/d20;
complex<T> d22 = s235*square(square(spb35))*T(6); d22 = T(1)/d22;
complex<T> t30 = d2*s234; 
complex<T> t31 = d17*s245; 
complex<T> t40 = -t57; 
complex<T> t62 = d1*spb24; 
complex<T> t63 = d20*t27; 
complex<T> t67 = t20*t22; 
complex<T> t74 = d18*spb24; 
complex<T> t76 = t21*t23; 
complex<T> t82 = t57*t58; 
complex<T> t88 = -(t21*T(2)); 
complex<T> t89 = t22*t24; 
complex<T> t103 = spb25*t20; 
complex<T> t104 = t23*t24; 
complex<T> t123 = t57*t9; 
complex<T> d4 = spa35*t25*T(3); d4 = T(1)/d4;
complex<T> d5 = (s25 + s35)*spa35*t25*T(3); d5 = T(1)/d5;
complex<T> d6 = s235*spa35*t25*T(3); d6 = T(1)/d6;
complex<T> d7 = t13*square(s25 + s35)*T(3); d7 = T(1)/d7;
complex<T> d8 = s235*t13*T(3); d8 = T(1)/d8;
complex<T> d9 = t47*square(s25 + s35); d9 = T(1)/d9;
complex<T> d10 = t47*cube(s25 + s35); d10 = T(1)/d10;
complex<T> d11 = s235*(s25 + s35)*t47; d11 = T(1)/d11;
complex<T> d12 = (s23 + s35)*spa35*t25*T(3); d12 = T(1)/d12;
complex<T> d13 = t13*square(s23 + s35)*T(3); d13 = T(1)/d13;
complex<T> d14 = t47*cube(s23 + s35); d14 = T(1)/d14;
complex<T> d15 = t47*square(s23 + s35); d15 = T(1)/d15;
complex<T> d16 = s235*(s23 + s35)*t47; d16 = T(1)/d16;
complex<T> d21 = s235*t25*T(3); d21 = T(1)/d21;
complex<T> d23 = s235*t25*T(6); d23 = T(1)/d23;
complex<T> d24 = s235*t13*T(6); d24 = T(1)/d24;
complex<T> d25 = s235*t13*T(12); d25 = T(1)/d25;
complex<T> t2 = d1*spa34*t104 + t106*t57*t62 + d3*t11*t9 + t104*t30*t9; 
complex<T> t3 = t108*t58*t74 - d19*spa45*t12*T(2) + spa45*t89*(d18 - t31*T(2)); 
complex<T> t29 = d21*s14; 
complex<T> t39 = d23*t10; 
complex<T> t49 = d21*spa24; 
complex<T> t69 = -(d4*T(4)); 
complex<T> t80 = t20*t63; 
complex<T> t84 = d10*spa35; 
complex<T> t95 = d14*spa35; 
complex<T> t96 = d8*spb25; 
complex<T> t1 = d22*t10*t20*t21*t27 + spa24*t21*t39*t57 + spa24*t20*t39*t58 + d25*t10*t8 + d24*t10*t82; 
complex<T> t4 = s25*t49*t57*t88 + s25*t80*t88 + s25*t103*t49*t9 + s25*t123*t96 + spa24*t29*(t21*t57 + t20*t58)*T(2) + s14*t21*t80*T(2) + d8*(s14*t8 - s25*t8 + s14*t82*T(2)); 
complex<T> t5 = -(d1*spa34*t104) + d9*spa35*t40*t58 + t106*t40*t62 + d7*s235*t76 + d9*spa35*t76 + d8*t8 + t69*t82 + d3*spa34*t11*T(2) + spa34*t104*t30*T(2) - d5*s235*t76*T(2) - d6*s25*t8*T(2) - d11*spa35*t8*T(2) - s235*t76*t84*T(2) + d4*t76*T(4); 
complex<T> t6 = d13*s235*t67 + d15*spa35*(t40*t58 + t67) - t108*t58*t74 + d8*t8 + t69*t82 - d18*spa45*t89 + d19*spa45*t12*T(2) - d12*s235*t67*T(2) - d6*s23*t8*T(2) - d16*spa35*t8*T(2) + spa45*t31*t89*T(2) - s235*t67*t95*T(2) + d4*t67*T(4); 
complex<T> t7 = -(d13*s235*t67) - d15*spa35*t67 + t67*t69 - d7*s235*t76 - d9*spa35*t76 + t69*t76 + d15*spa35*t82 + d9*spa35*t82 + d12*s235*t67*T(2) + d5*s235*t76*T(2) - d8*t8*T(2) + d6*s23*t8*T(2) + d6*s25*t8*T(2) + d11*spa35*t8*T(2) + d16*spa35*t8*T(2) + s235*t76*t84*T(2) + s235*t67*t95*T(2) + d4*t82*T(8); 
complex<T> t18 = s23*(-(d8*t8) + t49*t57*t88 + t80*t88 + t103*t49*t9 + t123*t96); 
SeriesC<T> result = t5*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t7*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t6*(*CI_users[4]->get_value(mc,ind,mu)) + t18*(*CI_users[5]->get_value(mc,ind,mu)) + t4*(*CI_users[6]->get_value(mc,ind,mu)) + t1*(*CI_users[7]->get_value(mc,ind,mu)) + t1*(*CI_users[8]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmpmm_G_wCI::\
C4g1ph_phmpmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c45));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmpmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmpmm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:03:17 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> s235 = SS(2,3,5);
complex<T> s14 = S(1,4);
complex<T> s35 = -(spa35*spb35);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s345 = SS(3,4,5);
complex<T> s12 = S(1,2);
complex<T> s15 = S(1,5);
complex<T> t19 = square(spa25*spb23 - spa45*spb34); 
complex<T> t23 = s23*s34; 
complex<T> t25 = spb23*spb34; 
complex<T> t29 = square(square(spa25*spb23 - spa45*spb34)); 
complex<T> t30 = s234*T(2); 
complex<T> t45 = square(mH2); 
complex<T> t46 = square(spa25); 
complex<T> t47 = square(spa45); 
complex<T> t48 = square(spb23); 
complex<T> t49 = square(spb34); 
complex<T> t50 = cube(s245); 
complex<T> t51 = square(square(spa24)); 
complex<T> t53 = square(spb35); 
complex<T> t54 = spa24*spb25 + spa34*spb35; 
complex<T> t55 = spa23*spb35 + spa24*spb45; 
complex<T> t56 = -(spa35*T(2)); 
complex<T> t58 = spa34*spb24 + spa35*spb25; 
complex<T> t59 = spa25*spb24 + spa35*spb34; 
complex<T> t60 = spa23*spb24 + spa35*spb45; 
complex<T> t66 = s234*T(3); 
complex<T> t76 = cube(spa25); 
complex<T> t77 = cube(spa45); 
complex<T> t78 = square(spa35); 
complex<T> t91 = spa23*spa34; 
complex<T> t103 = coeff3mass<T>(mc,ind,3,2,5,4); 
complex<T> t104 = coeff3mass<T>(mc,ind,3,4,5,2); 
complex<T> t114 = spa35*spb23; 
complex<T> t115 = square(square(s245)); 
complex<T> t127 = spa25*spa45; 
complex<T> t182 = spa24*spb34; 
complex<T> t186 = spa25*T(2); 
complex<T> t196 = spa24*spb23; 
complex<T> d1 = spa24*cube(spb24); d1 = T(1)/d1;
complex<T> d2 = spb24*square(s24 + s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s24 + s34)*spa24*cube(spb24); d3 = T(1)/d3;
complex<T> d4 = square(s24 + s34)*square(spb24); d4 = T(1)/d4;
complex<T> d5 = spb24*cube(s24 + s34)*T(3); d5 = T(1)/d5;
complex<T> d6 = s234*spa24*cube(spb24); d6 = T(1)/d6;
complex<T> d7 = s234*square(spb24); d7 = T(1)/d7;
complex<T> d9 = spb25*square(s25 + s35)*T(3); d9 = T(1)/d9;
complex<T> d10 = spb25*cube(s25 + s35)*T(3); d10 = T(1)/d10;
complex<T> d11 = s235*(s25 + s35)*spb25*T(3); d11 = T(1)/d11;
complex<T> d12 = (s23 + s24)*spa24*cube(spb24); d12 = T(1)/d12;
complex<T> d13 = square(s23 + s24)*square(spb24); d13 = T(1)/d13;
complex<T> d14 = spb24*cube(s23 + s24)*T(3); d14 = T(1)/d14;
complex<T> d15 = spb24*square(s23 + s24)*T(3); d15 = T(1)/d15;
complex<T> d17 = spb45*square(s35 + s45)*T(3); d17 = T(1)/d17;
complex<T> d18 = spb45*cube(s35 + s45)*T(3); d18 = T(1)/d18;
complex<T> d19 = s345*(s35 + s45)*spb45*T(3); d19 = T(1)/d19;
complex<T> d28 = s234*cube(spb24); d28 = T(1)/d28;
complex<T> d29 = s234*square(square(spb24)); d29 = T(1)/d29;
complex<T> d36 = s234*spb24; d36 = T(1)/d36;
complex<T> t20 = square(spa25*spb35 + t182); 
complex<T> t21 = square(spa45*spb35 + t196); 
complex<T> t22 = d7*T(2); 
complex<T> t24 = d36*spa24; 
complex<T> t26 = cube(spa25*spb35 + t182); 
complex<T> t27 = cube(spa45*spb35 + t196); 
complex<T> t34 = s235*t54; 
complex<T> t35 = s345*t55; 
complex<T> t36 = t58*t60; 
complex<T> t57 = spa45*spb24 + t114; 
complex<T> t67 = d10*s235; 
complex<T> t72 = d18*s345; 
complex<T> t79 = -(t48*T(2)); 
complex<T> t90 = t54*t55; 
complex<T> t94 = t59*t60; 
complex<T> t110 = -t127; 
complex<T> t112 = spa24*t49; 
complex<T> t123 = d28*s15; 
complex<T> t130 = spb34*t48; 
complex<T> t134 = d29*t78; 
complex<T> t138 = d9*spb35; 
complex<T> t141 = t127*t25; 
complex<T> t145 = d7*t19; 
complex<T> t152 = d28*t78; 
complex<T> t153 = t46*t48; 
complex<T> t155 = d1*T(4); 
complex<T> t164 = t47*t49; 
complex<T> t171 = t45*t76; 
complex<T> t179 = t45*t77; 
complex<T> t187 = spa45*t114; 
complex<T> t198 = t48*T(2); 
complex<T> t229 = s234*t46; 
complex<T> d8 = (s24 + s34)*spb24*t66; d8 = T(1)/d8;
complex<T> d16 = (s23 + s24)*spb24*t66; d16 = T(1)/d16;
complex<T> d53 = spb34*spb45*t55*T(2); d53 = T(1)/d53;
complex<T> d54 = spb23*spb25*t54*T(2); d54 = T(1)/d54;
complex<T> d58 = t30*cube(spb24); d58 = T(1)/d58;
complex<T> d59 = t30*square(spb24); d59 = T(1)/d59;
complex<T> d60 = t30*square(square(spb24)); d60 = T(1)/d60;
complex<T> t33 = t57*t59; 
complex<T> t93 = t57*t58; 
complex<T> t139 = d19*t20; 
complex<T> t157 = t134*t49; 
complex<T> t173 = spa35*t130; 
complex<T> t185 = t53*t67; 
complex<T> t199 = t112*t47; 
complex<T> t207 = t19*t24; 
complex<T> t211 = t187*t49; 
complex<T> d22 = spb23*spb25*t34; d22 = T(1)/d22;
complex<T> d23 = s234*t90*t91; d23 = T(1)/d23;
complex<T> d24 = spb34*spb45*t35; d24 = T(1)/d24;
complex<T> d25 = s235*spa23*t94; d25 = T(1)/d25;
complex<T> d26 = spb25*spb45*t36; d26 = T(1)/d26;
complex<T> d31 = spb25*t34*T(2); d31 = T(1)/d31;
complex<T> d32 = spa34*t30*t90; d32 = T(1)/d32;
complex<T> d33 = spb34*spb45*t35*T(2); d33 = T(1)/d33;
complex<T> d34 = s235*t94*T(2); d34 = T(1)/d34;
complex<T> d35 = spb25*spb45*t36*T(2); d35 = T(1)/d35;
complex<T> d37 = spb25*t36; d37 = T(1)/d37;
complex<T> d39 = spb23*t34*T(2); d39 = T(1)/d39;
complex<T> d40 = t30*t90*t91; d40 = T(1)/d40;
complex<T> d41 = s235*spa23*t94*T(2); d41 = T(1)/d41;
complex<T> d42 = spb45*t36*T(2); d42 = T(1)/d42;
complex<T> d43 = spb23*t34; d43 = T(1)/d43;
complex<T> d46 = spb23*spb25*t34*T(2); d46 = T(1)/d46;
complex<T> d47 = spa23*t30*t90; d47 = T(1)/d47;
complex<T> d48 = spb45*t35*T(2); d48 = T(1)/d48;
complex<T> d49 = spb34*t35; d49 = T(1)/d49;
complex<T> d50 = spb34*t35*T(2); d50 = T(1)/d50;
complex<T> d51 = spb25*t36*T(2); d51 = T(1)/d51;
complex<T> d55 = spa23*t94*T(2); d55 = T(1)/d55;
complex<T> d57 = t90*t91*T(2); d57 = T(1)/d57;
complex<T> d62 = s234*t90*T(4); d62 = T(1)/d62;
complex<T> d63 = t36*T(4); d63 = T(1)/d63;
complex<T> d64 = t34*T(2); d64 = T(1)/d64;
complex<T> d65 = t35*T(2); d65 = T(1)/d65;
complex<T> t2 = d13*s234*t153 + d15*spa24*t153 + t153*t155 + d4*s234*t164 + t155*t164 + d2*t199 + t19*t22 + d15*spa24*t110*t25 + d2*spa24*t110*t25 + d12*t229*t79 + d14*spa24*t229*t79 - d3*s234*t164*T(2) - d6*s23*t19*T(2) - d6*s34*t19*T(2) - d5*s234*t199*T(2) - d1*t141*T(8) - d16*spa24*t19*T(11) - d8*spa24*t19*T(11); 
complex<T> t4 = -t145 + d15*spa24*(t141 - t153) - d13*s234*t153 + t141*t155 + d17*spb35*t127*t182 + d12*t153*t30 + d14*spa24*t153*t30 + d6*s23*t19*T(2) - d1*t153*T(4) + d16*spa24*t19*T(11) + spa45*(t46*t53*(d17 - t72*T(2)) - t139*T(11)); 
complex<T> t5 = d17*spb35*t110*t182 + spa45*(-(t46*t53*(d17 - t72*T(2))) + t139*T(11)); 
complex<T> t8 = d2*spa24*t141 - t145 + t141*t155 - d4*s234*t164 + t127*t138*t196 - d2*t199 + d3*t164*t30 + d5*t199*t30 + d9*spa25*t47*t53 + d6*s34*t19*T(2) - spa25*t185*t47*T(2) - d1*t164*T(4) + d8*spa24*t19*T(11) - d11*spa25*t21*T(11); 
complex<T> t9 = t110*t138*t196 + t185*t186*t47 - d9*spa25*t47*t53 + d11*spa25*t21*T(11); 
complex<T> t14 = s34*(d55*t171 + d54*t27); 
complex<T> t18 = d25*(s14 - s25)*t171 + (d22*s14 + d43*spa25)*t27; 
complex<T> t39 = (d26*s13 + d37*spa45)*t50; 
complex<T> t43 = s45*(d55*t171 + d54*t27); 
complex<T> t44 = -(d34*s25*spb23*t171) + d64*spa23*spa25*t27; 
complex<T> t70 = d35*s34; 
complex<T> t97 = d35*s23; 
complex<T> d20 = s345*spa34*t93; d20 = T(1)/d20;
complex<T> d21 = s234*t25*t33; d21 = T(1)/d21;
complex<T> d27 = s345*spa34*t93*T(2); d27 = T(1)/d27;
complex<T> d30 = spb34*t30*t33; d30 = T(1)/d30;
complex<T> d38 = t25*t30*t33; d38 = T(1)/d38;
complex<T> d44 = s345*t93*T(2); d44 = T(1)/d44;
complex<T> d45 = spb23*t30*t33; d45 = T(1)/d45;
complex<T> d52 = spa34*t93*T(2); d52 = T(1)/d52;
complex<T> d56 = t25*t33*T(2); d56 = T(1)/d56;
complex<T> d61 = s234*t33*T(4); d61 = T(1)/d61;
complex<T> t1 = d59*t141*t23 + t145*t23 + d58*spa25*t173*t23 + d58*t211*t23 + d62*t25*t45*t51 + d60*t23*t48*t49*t78 + d61*t29*t91; 
complex<T> t3 = s25*(d56*t29 + d57*t45*t51); 
complex<T> t6 = s23*(d52*t179 + d53*t26); 
complex<T> t7 = d20*(s12 - s45)*t179 + (d24*s12 + d49*spa45)*t26; 
complex<T> t10 = d25*s15*t171 + d20*s15*t179 - d24*mH2*t26 + d24*s15*t26 - d22*mH2*t27 + d22*s15*t27 - d21*mH2*t29 + d21*s15*t29 - d26*mH2*t50 + d26*s15*t50 + d23*s15*t45*t51 - d23*t51*cube(mH2) - d25*t76*cube(mH2) - d20*t77*cube(mH2); 
complex<T> t11 = -(d25*s24*t171) - d20*s24*t179 + t112*t152*t198 + t112*t187*t22 - d24*s24*t26 - d22*s24*t27 + spa25*spa35*t182*t22*t48 - d26*s24*t50 + t141*t24*T(2) + t207*T(4); 
complex<T> t12 = d25*s24*t171 + d20*s24*t179 + d24*s24*t26 + d22*s24*t27 + d26*s24*t50 + spa25*t123*t130*t56 + spa45*spb23*t123*t49*t56 + t112*t152*t79 - t141*t24*T(2) + d7*(spa25*t182*t48*t56 + spa45*t196*t49*t56 - s15*t141*T(2)) - t207*T(4) - s15*(d25*t171 + d20*t179 + d24*t26 + d22*t27 + d26*t50 - t157*t79 + t145*T(4)); 
complex<T> t13 = -(d44*spb34*t179) + d48*spa34*t26 - d45*spa34*t29 + d47*spb34*t45*t51 + t50*t70 + s34*(d41*t171 + d28*t173*t186 + t157*t198 + t141*t22 - d46*t27 + d28*t211*T(2) + t145*T(4)); 
complex<T> t15 = d41*s45*t171 - d50*spa45*t26 - d46*s45*t27 + d38*s45*t29 - d40*s45*t45*t51 + d27*s45*t179*T(3) - d51*spa45*t50*T(3); 
complex<T> t16 = -(d34*spb23*t171) + d31*spa23*t27 - d30*spa23*t29 + d32*spb23*t45*t51 + t50*t97 + s23*(d27*t179 + d28*t173*t186 + t157*t198 + t141*t22 - d33*t26 + d28*t211*T(2) + t145*T(4)); 
complex<T> t17 = d27*s25*t179 - d33*s25*t26 - d39*spa25*t27 + d38*s25*t29 - d42*spa25*t50 - d40*s25*t45*t51 + d41*s25*t171*T(3); 
complex<T> t41 = s25*(d52*t179 + d53*t26); 
complex<T> t42 = -(d44*s45*spb34*t179) + d65*spa34*spa45*t26; 
complex<T> t75 = s45*(d56*t29 + d57*t45*t51); 
complex<T> co1 = -t103; 
complex<T> co2 = -t104; 
complex<T> co3 = t115*t97; 
complex<T> co4 = t115*t70; 
complex<T> co5 = d63*t127*t50; 
SeriesC<T> result = t8*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t9*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t10*(*CI_users[7]->get_value(mc,ind,mu)) + t16*(*CI_users[8]->get_value(mc,ind,mu)) + t11*(*CI_users[9]->get_value(mc,ind,mu)) + t39*(*CI_users[10]->get_value(mc,ind,mu)) + t17*(*CI_users[11]->get_value(mc,ind,mu)) + t18*(*CI_users[12]->get_value(mc,ind,mu)) + t12*(*CI_users[13]->get_value(mc,ind,mu)) + t13*(*CI_users[14]->get_value(mc,ind,mu)) + t7*(*CI_users[15]->get_value(mc,ind,mu)) + t15*(*CI_users[16]->get_value(mc,ind,mu)) + t6*(*CI_users[17]->get_value(mc,ind,mu)) + t41*(*CI_users[18]->get_value(mc,ind,mu)) + co3*(*CI_users[19]->get_value(mc,ind,mu)) + co4*(*CI_users[20]->get_value(mc,ind,mu)) + t14*(*CI_users[21]->get_value(mc,ind,mu)) + t43*(*CI_users[22]->get_value(mc,ind,mu)) + t3*(*CI_users[23]->get_value(mc,ind,mu)) + t75*(*CI_users[24]->get_value(mc,ind,mu)) + t1*(*CI_users[25]->get_value(mc,ind,mu)) + co5*(*CI_users[26]->get_value(mc,ind,mu)) + t1*(*CI_users[27]->get_value(mc,ind,mu)) + co5*(*CI_users[28]->get_value(mc,ind,mu)) + t44*(*CI_users[29]->get_value(mc,ind,mu)) + t42*(*CI_users[30]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmpmm_nf_wCI::\
C4g1ph_phmpmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmpmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmpmm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:04:15 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = S(2,3);
complex<T> s234 = SS(2,3,4);
complex<T> s34 = S(3,4);
complex<T> s15 = S(1,5);
complex<T> s24 = -(spa24*spb24);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = -(spa35*spb35);
complex<T> s235 = SS(2,3,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> t9 = square(spa25*spb23 - spa45*spb34); 
complex<T> t10 = -(spa25*T(2)); 
complex<T> t11 = -(s23*s34); 
complex<T> t12 = square(spa24*spb34 + spa25*spb35); 
complex<T> t13 = square(spa24*spb23 + spa45*spb35); 
complex<T> t15 = square(spb24); 
complex<T> t22 = square(spb23); 
complex<T> t23 = square(spb34); 
complex<T> t24 = square(spa25); 
complex<T> t25 = square(spa45); 
complex<T> t26 = cube(spb24); 
complex<T> t27 = square(spb35); 
complex<T> t29 = square(spa35); 
complex<T> t49 = spb24*T(9); 
complex<T> t60 = spa45*spb23; 
complex<T> t61 = spa25*spb34; 
complex<T> t104 = spa24*spa25; 
complex<T> d9 = spb25*square(s25 + s35)*T(9); d9 = T(1)/d9;
complex<T> d10 = spb25*cube(s25 + s35)*T(9); d10 = T(1)/d10;
complex<T> d11 = s235*(s25 + s35)*spb25*T(9); d11 = T(1)/d11;
complex<T> d17 = spb45*square(s35 + s45)*T(9); d17 = T(1)/d17;
complex<T> d18 = spb45*cube(s35 + s45)*T(9); d18 = T(1)/d18;
complex<T> d19 = s345*(s35 + s45)*spb45*T(9); d19 = T(1)/d19;
complex<T> d21 = s234*square(square(spb24))*T(3); d21 = T(1)/d21;
complex<T> d22 = s234*spb24*T(3); d22 = T(1)/d22;
complex<T> d25 = s234*square(square(spb24))*T(6); d25 = T(1)/d25;
complex<T> t32 = d10*s235; 
complex<T> t35 = d18*s345; 
complex<T> t42 = -t60; 
complex<T> t64 = d17*spb35; 
complex<T> t66 = d21*t29; 
complex<T> t69 = t22*t24; 
complex<T> t76 = d9*spb35; 
complex<T> t78 = t23*t25; 
complex<T> t83 = t22*t29; 
complex<T> t85 = t60*t61; 
complex<T> t92 = spb34*t10; 
complex<T> t93 = spa45*t27; 
complex<T> t97 = -(t23*T(2)); 
complex<T> t102 = t23*t60; 
complex<T> d1 = spa24*t26*T(3); d1 = T(1)/d1;
complex<T> d2 = t49*square(s24 + s34); d2 = T(1)/d2;
complex<T> d3 = (s24 + s34)*spa24*t26*T(3); d3 = T(1)/d3;
complex<T> d4 = t15*square(s24 + s34)*T(3); d4 = T(1)/d4;
complex<T> d5 = t49*cube(s24 + s34); d5 = T(1)/d5;
complex<T> d6 = s234*spa24*t26*T(3); d6 = T(1)/d6;
complex<T> d7 = s234*t15*T(3); d7 = T(1)/d7;
complex<T> d8 = s234*(s24 + s34)*t49; d8 = T(1)/d8;
complex<T> d12 = (s23 + s24)*spa24*t26*T(3); d12 = T(1)/d12;
complex<T> d13 = t15*square(s23 + s24)*T(3); d13 = T(1)/d13;
complex<T> d14 = t49*cube(s23 + s24); d14 = T(1)/d14;
complex<T> d15 = t49*square(s23 + s24); d15 = T(1)/d15;
complex<T> d16 = s234*(s23 + s24)*t49; d16 = T(1)/d16;
complex<T> d20 = s234*t26*T(3); d20 = T(1)/d20;
complex<T> d23 = s234*t26*T(6); d23 = T(1)/d23;
complex<T> d24 = s234*t15*T(6); d24 = T(1)/d24;
complex<T> d26 = s234*t15*T(12); d26 = T(1)/d26;
complex<T> t7 = d11*t10*t13 + d9*spa25*t25*t27 + t10*t25*t27*t32 + t104*t60*t76; 
complex<T> t31 = d20*s15; 
complex<T> t41 = d23*t11; 
complex<T> t50 = -(d7*t9); 
complex<T> t52 = d20*spa35; 
complex<T> t53 = d6*s23; 
complex<T> t68 = d6*s34; 
complex<T> t71 = -(d1*T(4)); 
complex<T> t73 = d7*spa35; 
complex<T> t88 = d14*s234; 
complex<T> t94 = t22*t66; 
complex<T> t99 = d5*s234; 
complex<T> t111 = spa24*t64; 
complex<T> t112 = t23*t83; 
complex<T> t117 = t60*t92; 
complex<T> t1 = d25*t11*t112 + spa35*t41*(t102 + t22*t61) + t11*(d24*t85 + d26*t9); 
complex<T> t2 = spa24*(d22*(t117 - t9) + t22*t73*t92 + t60*t73*t97 + d20*t83*t97); 
complex<T> t5 = spa45*t111*t61 - d19*spa45*t12*T(2) + t24*t93*(d17 - t35*T(2)); 
complex<T> t20 = s23*(d7*t117 + t50 + t22*t52*t92 + t52*t60*t97 + t94*t97); 
complex<T> t59 = s34*(d7*t117 + t50 + t22*t52*t92 + t52*t60*t97 + t94*t97); 
complex<T> t90 = spa35*t31; 
complex<T> t91 = t53*t9; 
complex<T> t96 = t68*t9; 
complex<T> t3 = d7*s15*t9 + d20*spa24*t112*T(2) + d7*s15*t85*T(2) + t102*t90*T(2) + t22*t61*t90*T(2) + s15*t23*t94*T(2) + spa24*(d22*t9 + t102*t73*T(2) + t22*t61*t73*T(2) + d22*t85*T(2)); 
complex<T> t4 = -(d13*s234*t69) - d15*spa24*t69 + t69*t71 - d4*s234*t78 - d2*spa24*t78 + t71*t78 + d15*spa24*t85 + d2*spa24*t85 + d12*s234*t69*T(2) + d3*s234*t78*T(2) + spa24*t69*t88*T(2) - d7*t9*T(2) + d16*spa24*t9*T(2) + d8*spa24*t9*T(2) + t91*T(2) + t96*T(2) + spa24*t78*t99*T(2) + d1*t85*T(8); 
complex<T> t6 = -(spa45*t111*t61) + d15*spa24*t42*t61 + d13*s234*t69 + d15*spa24*t69 + t71*t85 + d7*t9 - d17*t24*t93 + d19*spa45*t12*T(2) - d12*s234*t69*T(2) - spa24*t69*t88*T(2) - d16*spa24*t9*T(2) - t91*T(2) + t24*t35*t93*T(2) + d1*t69*T(4); 
complex<T> t8 = -(d9*spa25*t25*t27) + d2*spa24*t42*t61 + t104*t42*t76 + d4*s234*t78 + d2*spa24*t78 + t71*t85 + d7*t9 + d11*spa25*t13*T(2) + spa25*t25*t27*t32*T(2) - d3*s234*t78*T(2) - d8*spa24*t9*T(2) - t96*T(2) - spa24*t78*t99*T(2) + d1*t78*T(4); 
SeriesC<T> result = t8*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t7*(*CI_users[2]->get_value(mc,ind,mu)) + t6*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + t20*(*CI_users[5]->get_value(mc,ind,mu)) + t2*(*CI_users[6]->get_value(mc,ind,mu)) + t3*(*CI_users[7]->get_value(mc,ind,mu)) + t59*(*CI_users[8]->get_value(mc,ind,mu)) + t1*(*CI_users[9]->get_value(mc,ind,mu)) + t1*(*CI_users[10]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmmpm_G_wCI::\
C4g1ph_phmmpm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c45, c23));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmmpm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmmpm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:07:35 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> s235 = SS(2,3,5);
complex<T> s14 = S(1,4);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s345 = SS(3,4,5);
complex<T> s12 = S(1,2);
complex<T> s15 = S(1,5);
complex<T> s35 = -(spa35*spb35);
complex<T> t18 = square(spa23*spb34 - spa25*spb45); 
complex<T> t19 = square(spa23*spb24 + spa35*spb45); 
complex<T> t22 = square(square(spa23*spb34 - spa25*spb45)); 
complex<T> t23 = cube(spa23*spb24 + spa35*spb45); 
complex<T> t24 = spb34*spb45; 
complex<T> t27 = s345*T(2); 
complex<T> t43 = square(mH2); 
complex<T> t44 = square(spa25); 
complex<T> t45 = square(spa23); 
complex<T> t46 = square(spb34); 
complex<T> t47 = square(spb45); 
complex<T> t48 = cube(s235); 
complex<T> t49 = square(square(spa35)); 
complex<T> t50 = square(spb24); 
complex<T> t51 = spa34*spb24 + spa35*spb25; 
complex<T> t53 = spa35*spb23 + spa45*spb24; 
complex<T> t54 = spa24*spb34 + spa25*spb35; 
complex<T> t55 = spa24*spb23 + spa45*spb35; 
complex<T> t56 = spa24*spb25 + spa34*spb35; 
complex<T> t57 = spa23*spb35 + spa24*spb45; 
complex<T> t58 = -(spa24*T(2)); 
complex<T> t59 = s34*s45; 
complex<T> t73 = cube(spa25); 
complex<T> t74 = cube(spa23); 
complex<T> t79 = square(spa24); 
complex<T> t82 = s24 + s25; 
complex<T> t83 = spb35*T(3); 
complex<T> t89 = spa34*spa45; 
complex<T> t91 = spb25*spb45; 
complex<T> t100 = coeff3mass<T>(mc,ind,4,3,2,5); 
complex<T> t101 = coeff3mass<T>(mc,ind,4,5,2,3); 
complex<T> t113 = square(square(s235)); 
complex<T> t116 = spb23*spb34; 
complex<T> t127 = spa23*spa25; 
complex<T> t143 = s34*T(2); 
complex<T> t181 = spa25*spb34; 
complex<T> t198 = spa35*spb34; 
complex<T> d1 = spb23*cube(s23 + s24)*T(3); d1 = T(1)/d1;
complex<T> d2 = spb23*square(s23 + s24)*T(3); d2 = T(1)/d2;
complex<T> d3 = s234*(s23 + s24)*spb23*T(3); d3 = T(1)/d3;
complex<T> d7 = spa35*cube(spb35); d7 = T(1)/d7;
complex<T> d9 = (s35 + s45)*spa35*cube(spb35); d9 = T(1)/d9;
complex<T> d10 = square(s35 + s45)*square(spb35); d10 = T(1)/d10;
complex<T> d12 = s345*spa35*cube(spb35); d12 = T(1)/d12;
complex<T> d13 = s345*square(spb35); d13 = T(1)/d13;
complex<T> d15 = (s34 + s35)*spa35*cube(spb35); d15 = T(1)/d15;
complex<T> d16 = square(s34 + s35)*square(spb35); d16 = T(1)/d16;
complex<T> d43 = s345*cube(spb35); d43 = T(1)/d43;
complex<T> d44 = s345*square(square(spb35)); d44 = T(1)/d44;
complex<T> t20 = square(spa25*spb24 + t198); 
complex<T> t26 = cube(spa25*spb24 + t198); 
complex<T> t28 = t51*t53; 
complex<T> t29 = t54*t57; 
complex<T> t32 = t55*t56; 
complex<T> t63 = d43*s12; 
complex<T> t65 = d1*s234; 
complex<T> t76 = -(t46*T(2)); 
complex<T> t88 = t54*t55; 
complex<T> t90 = t56*t57; 
complex<T> t97 = d12*s45; 
complex<T> t109 = -t127; 
complex<T> t110 = t43*t49; 
complex<T> t111 = -(t47*T(2)); 
complex<T> t112 = d13*t18; 
complex<T> t114 = d43*spa24; 
complex<T> t132 = spb45*t46; 
complex<T> t136 = d44*t79; 
complex<T> t152 = t43*t73; 
complex<T> t154 = d7*T(4); 
complex<T> t165 = t127*t24; 
complex<T> t182 = t43*t74; 
complex<T> t193 = spa23*t44; 
complex<T> t197 = spa25*t45; 
complex<T> d4 = spb25*cube(t82)*T(3); d4 = T(1)/d4;
complex<T> d5 = spb25*square(t82)*T(3); d5 = T(1)/d5;
complex<T> d6 = s245*spb25*t82*T(3); d6 = T(1)/d6;
complex<T> d8 = t83*square(s35 + s45); d8 = T(1)/d8;
complex<T> d11 = t83*cube(s35 + s45); d11 = T(1)/d11;
complex<T> d14 = s345*(s35 + s45)*t83; d14 = T(1)/d14;
complex<T> d17 = t83*cube(s34 + s35); d17 = T(1)/d17;
complex<T> d18 = t83*square(s34 + s35); d18 = T(1)/d18;
complex<T> d19 = s345*(s34 + s35)*t83; d19 = T(1)/d19;
complex<T> d21 = s234*t116*t53; d21 = T(1)/d21;
complex<T> d26 = s245*t51*t91; d26 = T(1)/d26;
complex<T> d28 = s234*spb34*t53*T(2); d28 = T(1)/d28;
complex<T> d33 = s245*t51*t91*T(2); d33 = T(1)/d33;
complex<T> d35 = s245*spb25*t51; d35 = T(1)/d35;
complex<T> d36 = s234*t116*t53*T(2); d36 = T(1)/d36;
complex<T> d38 = s245*spb45*t51*T(2); d38 = T(1)/d38;
complex<T> d41 = s234*spb23*t53*T(2); d41 = T(1)/d41;
complex<T> d52 = s245*spb25*t51*T(2); d52 = T(1)/d52;
complex<T> d56 = t51*t91*T(2); d56 = T(1)/d56;
complex<T> d57 = t116*t53*T(2); d57 = T(1)/d57;
complex<T> d59 = s234*t53*T(2); d59 = T(1)/d59;
complex<T> d60 = s245*t51*T(2); d60 = T(1)/d60;
complex<T> d63 = t27*cube(spb35); d63 = T(1)/d63;
complex<T> d64 = t27*square(spb35); d64 = T(1)/d64;
complex<T> d65 = t27*square(square(spb35)); d65 = T(1)/d65;
complex<T> t68 = d4*s245; 
complex<T> t75 = -t110; 
complex<T> t99 = d3*t20; 
complex<T> t138 = d5*spb24; 
complex<T> t140 = d63*t59; 
complex<T> t147 = t136*t47; 
complex<T> t159 = s45*t114; 
complex<T> t207 = t18*t97; 
complex<T> d20 = s345*t28*t89; d20 = T(1)/d20;
complex<T> d22 = s245*spa45*t88; d22 = T(1)/d22;
complex<T> d23 = spb23*spb25*t32; d23 = T(1)/d23;
complex<T> d24 = s234*spa34*t90; d24 = T(1)/d24;
complex<T> d25 = s345*t24*t29; d25 = T(1)/d25;
complex<T> d27 = t27*t28*t89; d27 = T(1)/d27;
complex<T> d29 = s245*spa45*t88*T(2); d29 = T(1)/d29;
complex<T> d30 = spb25*t32*T(2); d30 = T(1)/d30;
complex<T> d31 = s234*spa34*t90*T(2); d31 = T(1)/d31;
complex<T> d32 = t24*t27*t29; d32 = T(1)/d32;
complex<T> d34 = s245*t88; d34 = T(1)/d34;
complex<T> d37 = spb23*t32*T(2); d37 = T(1)/d37;
complex<T> d39 = spb23*t32; d39 = T(1)/d39;
complex<T> d40 = spa45*t27*t28; d40 = T(1)/d40;
complex<T> d42 = spb23*spb25*t32*T(2); d42 = T(1)/d42;
complex<T> d45 = s234*t90*T(2); d45 = T(1)/d45;
complex<T> d46 = spb45*t27*t29; d46 = T(1)/d46;
complex<T> d47 = s345*spa34*t28; d47 = T(1)/d47;
complex<T> d48 = s345*spb34*t29; d48 = T(1)/d48;
complex<T> d49 = spa34*t27*t28; d49 = T(1)/d49;
complex<T> d50 = s245*t88*T(2); d50 = T(1)/d50;
complex<T> d51 = spb34*t27*t29; d51 = T(1)/d51;
complex<T> d53 = t28*t89*T(2); d53 = T(1)/d53;
complex<T> d54 = t24*t29*T(2); d54 = T(1)/d54;
complex<T> d55 = spa45*t88*T(2); d55 = T(1)/d55;
complex<T> d58 = spa34*t90*T(2); d58 = T(1)/d58;
complex<T> d61 = t32*T(4); d61 = T(1)/d61;
complex<T> d62 = s345*t28*T(4); d62 = T(1)/d62;
complex<T> d66 = s345*t29*T(4); d66 = T(1)/d66;
complex<T> t1 = spa23*spa24*t132*t140 + d62*t110*t24 + spa24*t140*t181*t47 + t112*t59 + d64*t165*t59 + d65*t46*t47*t59*t79 + d66*t22*t89; 
complex<T> t2 = d2*spb24*t109*t198 - d2*t193*t50 + t193*t50*t65*T(2) + spa23*t99*T(11); 
complex<T> t3 = s25*(d58*t182 + d57*t26); 
complex<T> t4 = -t112 + t154*t165 + d2*spb24*t127*t198 - d10*s345*t44*t47 + d9*t27*t44*t47 + d11*spa35*t27*t44*t47 + d8*spa35*(t165 - t44*t47) + d2*t193*t50 + t207*T(2) - t193*t50*t65*T(2) - d7*t44*t47*T(4) + d14*spa35*t18*T(11) - spa23*t99*T(11); 
complex<T> t5 = d18*spa35*t109*t24 + d8*spa35*t109*t24 + d9*s345*t111*t44 + d11*s345*spa35*t111*t44 + d16*s345*t45*t46 + d18*spa35*t45*t46 + t154*t45*t46 + d10*s345*t44*t47 + d8*spa35*t44*t47 + t154*t44*t47 + d15*s345*t45*t76 + d17*s345*spa35*t45*t76 + t112*T(2) - d12*s34*t18*T(2) - t207*T(2) - d7*t165*T(8) - d14*spa35*t18*T(11) - d19*spa35*t18*T(11); 
complex<T> t6 = s23*(d53*t110 + d54*t22); 
complex<T> t7 = d20*s12*t110 + d47*spb45*t110 + d25*s12*t22 + d48*spa45*t22 + spa23*t132*t58*t63 + t181*t47*t58*t63 + s12*t147*t76 + spa23*t132*t159*T(2) - d13*s12*t165*T(2) + d13*s45*t165*T(2) + s45*t147*t46*T(2) + t159*t181*t47*T(2) - s12*t112*T(4) + s45*t112*T(4); 
complex<T> t10 = d20*s15*t110 + d22*s15*t152 + d24*s15*t182 - d25*mH2*t22 + d25*s15*t22 - d26*mH2*t23 + d26*s15*t23 - d21*mH2*t26 + d21*s15*t26 - d23*mH2*t48 + d23*s15*t48 - d20*t49*cube(mH2) - d22*t73*cube(mH2) - d24*t74*cube(mH2); 
complex<T> t12 = -(s24*(d22*t152 + d25*t22 + d26*t23 + d23*t48 - d20*t75)); 
complex<T> t14 = s34*(d55*t152 + d56*t23); 
complex<T> t16 = d22*s13*t152 + d34*spb45*t152 + d26*s13*t23 + d35*spa45*t23; 
complex<T> t35 = (d23*s14 + d39*spa25)*t48; 
complex<T> t37 = s25*(d53*t110 + d54*t22); 
complex<T> t38 = s23*(d55*t152 + d56*t23); 
complex<T> t39 = -((s15 - s24)*(d22*t152 + d25*t22 + d26*t23 + d23*t48)) + d20*(s24*t110 + s15*t75); 
complex<T> t40 = -(d50*s25*spb45*t152) + d60*spa25*spa45*t23; 
complex<T> t41 = s45*(d58*t182 + d57*t26); 
complex<T> t42 = -(d45*s23*spb34*t182) + d59*spa23*spa34*t26; 
complex<T> t95 = d42*s34; 
complex<T> t122 = d42*s45; 
complex<T> t168 = d31*t43; 
complex<T> t185 = spb45*t138; 
complex<T> t8 = -t112 + t154*t165 + d12*t143*t18 + spa35*t127*t185 - d16*s345*t45*t46 + d15*t27*t45*t46 + d17*spa35*t27*t45*t46 + d18*spa35*(t165 - t45*t46) + d5*t197*t50 - t197*t50*t68*T(2) - d7*t45*t46*T(4) + d19*spa35*t18*T(11) - d6*spa25*t19*T(11); 
complex<T> t9 = spa35*t109*t185 - t197*t50*(d5 - t68*T(2)) + d6*spa25*t19*T(11); 
complex<T> t11 = d29*s23*t152 + d32*s23*t22 - d33*s23*t23 + d28*spa23*t26 - d30*spa23*t48 + s23*t168*t74 + d27*s23*t75; 
complex<T> t13 = d40*spb34*t110 + spa23*t114*t132*t143 + d13*t143*t165 - d45*spb34*t182 - d46*spa34*t22 + d41*spa34*t26 + t143*t147*t46 + t114*t143*t181*t47 + t48*t95 + s34*(d29*t152 - d33*t23 + t112*T(4)); 
complex<T> t15 = d29*s25*t152 + d32*s25*t22 + d38*spa25*t23 - d36*s25*t26 + s25*t168*t74 + d27*s25*t75 - d37*spa25*t48*T(3); 
complex<T> t17 = -(d52*spa45*t23) - d36*s45*t26 + t122*t48 + s45*t168*t74 + d49*spb45*t75 - d50*spb45*t152*T(3) - d51*spa45*t22*T(3); 
complex<T> co1 = -t100; 
complex<T> co2 = -t101; 
complex<T> co3 = t113*t95; 
complex<T> co4 = t113*t122; 
complex<T> co5 = d61*t127*t48; 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t4*(*CI_users[2]->get_value(mc,ind,mu)) + t5*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t10*(*CI_users[7]->get_value(mc,ind,mu)) + t11*(*CI_users[8]->get_value(mc,ind,mu)) + t12*(*CI_users[9]->get_value(mc,ind,mu)) + t16*(*CI_users[10]->get_value(mc,ind,mu)) + t15*(*CI_users[11]->get_value(mc,ind,mu)) + t35*(*CI_users[12]->get_value(mc,ind,mu)) + t39*(*CI_users[13]->get_value(mc,ind,mu)) + t13*(*CI_users[14]->get_value(mc,ind,mu)) + t7*(*CI_users[15]->get_value(mc,ind,mu)) + t17*(*CI_users[16]->get_value(mc,ind,mu)) + t6*(*CI_users[17]->get_value(mc,ind,mu)) + t37*(*CI_users[18]->get_value(mc,ind,mu)) + t38*(*CI_users[19]->get_value(mc,ind,mu)) + t14*(*CI_users[20]->get_value(mc,ind,mu)) + co3*(*CI_users[21]->get_value(mc,ind,mu)) + co4*(*CI_users[22]->get_value(mc,ind,mu)) + t3*(*CI_users[23]->get_value(mc,ind,mu)) + t41*(*CI_users[24]->get_value(mc,ind,mu)) + t42*(*CI_users[25]->get_value(mc,ind,mu)) + t40*(*CI_users[26]->get_value(mc,ind,mu)) + co5*(*CI_users[27]->get_value(mc,ind,mu)) + t1*(*CI_users[28]->get_value(mc,ind,mu)) + co5*(*CI_users[29]->get_value(mc,ind,mu)) + t1*(*CI_users[30]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmmpm_nf_wCI::\
C4g1ph_phmmpm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmmpm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmmpm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:08:20 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s234 = SS(2,3,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s245 = SS(2,4,5);
complex<T> s34 = S(3,4);
complex<T> s345 = SS(3,4,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> s35 = -(spa35*spb35);
complex<T> t8 = square(spa23*spb34 - spa25*spb45); 
complex<T> t9 = -(spa23*T(2)); 
complex<T> t10 = -(s34*s45); 
complex<T> t11 = square(spa25*spb24 + spa35*spb34); 
complex<T> t12 = square(spa23*spb24 + spa35*spb45); 
complex<T> t13 = square(spb35); 
complex<T> t20 = square(spb34); 
complex<T> t21 = square(spb45); 
complex<T> t22 = square(spa23); 
complex<T> t23 = square(spa25); 
complex<T> t24 = square(spb24); 
complex<T> t25 = cube(spb35); 
complex<T> t27 = square(spa24); 
complex<T> t47 = spb35*T(9); 
complex<T> t57 = spa25*spb34; 
complex<T> t58 = spa23*spb45; 
complex<T> t107 = spa23*spa35; 
complex<T> t109 = spa25*spa35; 
complex<T> d1 = spb23*cube(s23 + s24)*T(9); d1 = T(1)/d1;
complex<T> d2 = spb23*square(s23 + s24)*T(9); d2 = T(1)/d2;
complex<T> d3 = s234*(s23 + s24)*spb23*T(9); d3 = T(1)/d3;
complex<T> d4 = spb25*cube(s24 + s25)*T(9); d4 = T(1)/d4;
complex<T> d5 = spb25*square(s24 + s25)*T(9); d5 = T(1)/d5;
complex<T> d6 = s245*(s24 + s25)*spb25*T(9); d6 = T(1)/d6;
complex<T> d21 = s345*square(square(spb35))*T(3); d21 = T(1)/d21;
complex<T> d24 = s345*square(square(spb35))*T(6); d24 = T(1)/d24;
complex<T> t30 = d1*s234; 
complex<T> t32 = d4*s245; 
complex<T> t40 = -t57; 
complex<T> t62 = d2*spb24; 
complex<T> t63 = d21*t27; 
complex<T> t74 = d5*spb24; 
complex<T> t76 = t20*t22; 
complex<T> t77 = t21*t23; 
complex<T> t83 = t57*t58; 
complex<T> t89 = -(t21*T(2)); 
complex<T> t104 = spb45*t20; 
complex<T> t105 = t23*t24; 
complex<T> t126 = t57*t9; 
complex<T> d7 = spa35*t25*T(3); d7 = T(1)/d7;
complex<T> d8 = t47*square(s35 + s45); d8 = T(1)/d8;
complex<T> d9 = (s35 + s45)*spa35*t25*T(3); d9 = T(1)/d9;
complex<T> d10 = t13*square(s35 + s45)*T(3); d10 = T(1)/d10;
complex<T> d11 = t47*cube(s35 + s45); d11 = T(1)/d11;
complex<T> d12 = s345*spa35*t25*T(3); d12 = T(1)/d12;
complex<T> d13 = s345*t13*T(3); d13 = T(1)/d13;
complex<T> d14 = s345*(s35 + s45)*t47; d14 = T(1)/d14;
complex<T> d15 = (s34 + s35)*spa35*t25*T(3); d15 = T(1)/d15;
complex<T> d16 = t13*square(s34 + s35)*T(3); d16 = T(1)/d16;
complex<T> d17 = t47*cube(s34 + s35); d17 = T(1)/d17;
complex<T> d18 = t47*square(s34 + s35); d18 = T(1)/d18;
complex<T> d19 = s345*(s34 + s35)*t47; d19 = T(1)/d19;
complex<T> d20 = s345*t25*T(3); d20 = T(1)/d20;
complex<T> d22 = s345*t25*T(6); d22 = T(1)/d22;
complex<T> d23 = s345*t13*T(6); d23 = T(1)/d23;
complex<T> d25 = s345*t13*T(12); d25 = T(1)/d25;
complex<T> t2 = d2*spa23*t105 + t107*t57*t62 + d3*t11*t9 + t105*t30*t9; 
complex<T> t29 = d20*s12; 
complex<T> t39 = d22*t10; 
complex<T> t49 = d20*spa24; 
complex<T> t69 = -(d7*T(4)); 
complex<T> t85 = d11*spa35; 
complex<T> t91 = t20*t63; 
complex<T> t97 = d17*spa35; 
complex<T> t98 = d13*spb45; 
complex<T> t99 = t22*t32; 
complex<T> t1 = d24*t10*t20*t21*t27 + spa24*t21*t39*t57 + spa24*t20*t39*t58 + d25*t10*t8 + d23*t10*t83; 
complex<T> t3 = s45*t49*t57*t89 + s45*t104*t49*t9 + s45*t89*t91 + s45*t126*t98 + spa24*t29*(t21*t57 + t20*t58)*T(2) + s12*t21*t91*T(2) + d13*(s12*t8 - s45*t8 + s12*t83*T(2)); 
complex<T> t4 = -(d2*spa23*t105) + d8*spa35*t40*t58 + t107*t40*t62 + d10*s345*t77 + d8*spa35*t77 + d13*t8 + t69*t83 + d3*spa23*t11*T(2) + spa23*t105*t30*T(2) - d9*s345*t77*T(2) - d12*s45*t8*T(2) - d14*spa35*t8*T(2) - s345*t77*t85*T(2) + d7*t77*T(4); 
complex<T> t5 = -(d16*s345*t76) - d18*spa35*t76 + t69*t76 - d10*s345*t77 - d8*spa35*t77 + t69*t77 + d18*spa35*t83 + d8*spa35*t83 + d15*s345*t76*T(2) + d9*s345*t77*T(2) - d13*t8*T(2) + d12*s34*t8*T(2) + d12*s45*t8*T(2) + d14*spa35*t8*T(2) + d19*spa35*t8*T(2) + s345*t77*t85*T(2) + s345*t76*t97*T(2) + d7*t83*T(8); 
complex<T> t6 = d5*spa25*t22*t24 + t109*t58*t74 - d6*spa25*t12*T(2) - spa25*t24*t99*T(2); 
complex<T> t7 = -(d5*spa25*t22*t24) + d18*spa35*t40*t58 - t109*t58*t74 + d16*s345*t76 + d18*spa35*t76 + d13*t8 + t69*t83 + d6*spa25*t12*T(2) - d15*s345*t76*T(2) - d12*s34*t8*T(2) - d19*spa35*t8*T(2) - s345*t76*t97*T(2) + spa25*t24*t99*T(2) + d7*t76*T(4); 
complex<T> t18 = s34*(-(d13*t8) + t49*t57*t89 + t104*t49*t9 + t89*t91 + t126*t98); 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t4*(*CI_users[2]->get_value(mc,ind,mu)) + t5*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + t18*(*CI_users[5]->get_value(mc,ind,mu)) + t3*(*CI_users[6]->get_value(mc,ind,mu)) + t1*(*CI_users[7]->get_value(mc,ind,mu)) + t1*(*CI_users[8]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmmmp_G_wCI::\
C4g1ph_phmmmp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c25, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c45, c23));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmmmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmmmp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:11:58 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s14 = S(1,4);
complex<T> s35 = -(spa35*spb35);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s24 = -(spa24*spb24);
complex<T> s345 = SS(3,4,5);
complex<T> s12 = S(1,2);
complex<T> s15 = S(1,5);
complex<T> t18 = square(spa23*spb25 - spa34*spb45); 
complex<T> t24 = square(square(spa23*spb25 - spa34*spb45)); 
complex<T> t25 = spb25*spb45; 
complex<T> t27 = s245*T(2); 
complex<T> t42 = square(mH2); 
complex<T> t43 = square(spa23); 
complex<T> t44 = square(spa34); 
complex<T> t45 = square(spb25); 
complex<T> t46 = square(spb45); 
complex<T> t47 = square(square(spa24)); 
complex<T> t48 = cube(s234); 
complex<T> t50 = square(spb35); 
complex<T> t51 = spa24*spb34 + spa25*spb35; 
complex<T> t52 = spa24*spb23 + spa45*spb35; 
complex<T> t53 = spa35*spb23 + spa45*spb24; 
complex<T> t55 = spa25*spb24 + spa35*spb34; 
complex<T> t57 = -(spa35*T(2)); 
complex<T> t58 = s25*s45; 
complex<T> t66 = s245*T(3); 
complex<T> t73 = cube(spa23); 
complex<T> t74 = cube(spa34); 
complex<T> t79 = square(spa35); 
complex<T> t89 = square(square(s234)); 
complex<T> t90 = spa25*spa45; 
complex<T> t100 = coeff3mass<T>(mc,ind,5,2,3,4); 
complex<T> t101 = coeff3mass<T>(mc,ind,5,4,3,2); 
complex<T> t111 = spa35*spb25; 
complex<T> t125 = spa23*spa34; 
complex<T> t129 = spa35*spb45; 
complex<T> t131 = spb23*spb25; 
complex<T> t156 = spb34*spb45; 
complex<T> t186 = spa23*T(2); 
complex<T> t187 = spa24*spb25; 
complex<T> t195 = spa24*spb45; 
complex<T> d1 = spb23*square(s23 + s35)*T(3); d1 = T(1)/d1;
complex<T> d2 = spb23*cube(s23 + s35)*T(3); d2 = T(1)/d2;
complex<T> d3 = s235*(s23 + s35)*spb23*T(3); d3 = T(1)/d3;
complex<T> d4 = spa24*cube(spb24); d4 = T(1)/d4;
complex<T> d5 = (s24 + s25)*spa24*cube(spb24); d5 = T(1)/d5;
complex<T> d6 = square(s24 + s25)*square(spb24); d6 = T(1)/d6;
complex<T> d7 = spb24*cube(s24 + s25)*T(3); d7 = T(1)/d7;
complex<T> d8 = spb24*square(s24 + s25)*T(3); d8 = T(1)/d8;
complex<T> d9 = spb24*square(s24 + s45)*T(3); d9 = T(1)/d9;
complex<T> d10 = (s24 + s45)*spa24*cube(spb24); d10 = T(1)/d10;
complex<T> d11 = square(s24 + s45)*square(spb24); d11 = T(1)/d11;
complex<T> d12 = spb24*cube(s24 + s45)*T(3); d12 = T(1)/d12;
complex<T> d13 = s245*spa24*cube(spb24); d13 = T(1)/d13;
complex<T> d14 = s245*square(spb24); d14 = T(1)/d14;
complex<T> d17 = spb34*cube(s34 + s35)*T(3); d17 = T(1)/d17;
complex<T> d18 = spb34*square(s34 + s35)*T(3); d18 = T(1)/d18;
complex<T> d19 = s345*(s34 + s35)*spb34*T(3); d19 = T(1)/d19;
complex<T> d34 = s245*cube(spb24); d34 = T(1)/d34;
complex<T> d36 = s245*square(square(spb24)); d36 = T(1)/d36;
complex<T> t19 = square(spa34*spb35 + t187); 
complex<T> t20 = square(spa23*spb35 + t195); 
complex<T> t22 = cube(spa34*spb35 + t187); 
complex<T> t23 = cube(spa23*spb35 + t195); 
complex<T> t28 = t51*t52; 
complex<T> t32 = t53*t55; 
complex<T> t54 = spa34*spb24 + t111; 
complex<T> t56 = spa23*spb24 + t129; 
complex<T> t62 = d34*s13; 
complex<T> t64 = d2*s235; 
complex<T> t68 = d17*s345; 
complex<T> t76 = -(t45*T(2)); 
complex<T> t97 = d13*s45; 
complex<T> t107 = -t125; 
complex<T> t108 = t42*t47; 
complex<T> t109 = -(t46*T(2)); 
complex<T> t110 = d14*t18; 
complex<T> t137 = d18*spb35; 
complex<T> t141 = d4*T(4); 
complex<T> t144 = t46*t79; 
complex<T> t151 = t42*t74; 
complex<T> t163 = t42*t73; 
complex<T> t173 = t125*t25; 
complex<T> t180 = spa34*t46; 
complex<T> t184 = d36*t45; 
complex<T> t202 = t129*t45; 
complex<T> d15 = (s24 + s25)*spb24*t66; d15 = T(1)/d15;
complex<T> d16 = (s24 + s45)*spb24*t66; d16 = T(1)/d16;
complex<T> d23 = s235*t131*t52; d23 = T(1)/d23;
complex<T> d24 = s345*t156*t51; d24 = T(1)/d24;
complex<T> d30 = s235*spb25*t52*T(2); d30 = T(1)/d30;
complex<T> d31 = s345*t156*t51*T(2); d31 = T(1)/d31;
complex<T> d40 = s235*spb23*t52*T(2); d40 = T(1)/d40;
complex<T> d43 = s235*spb23*t52; d43 = T(1)/d43;
complex<T> d46 = s235*t131*t52*T(2); d46 = T(1)/d46;
complex<T> d47 = s345*spb45*t51*T(2); d47 = T(1)/d47;
complex<T> d49 = s345*spb34*t51; d49 = T(1)/d49;
complex<T> d52 = s345*spb34*t51*T(2); d52 = T(1)/d52;
complex<T> d55 = t156*t51*T(2); d55 = T(1)/d55;
complex<T> d58 = t131*t52*T(2); d58 = T(1)/d58;
complex<T> d61 = t27*square(spb24); d61 = T(1)/d61;
complex<T> d62 = t27*cube(spb24); d62 = T(1)/d62;
complex<T> d64 = t27*square(square(spb24)); d64 = T(1)/d64;
complex<T> d66 = s235*t52*T(2); d66 = T(1)/d66;
complex<T> d67 = s345*t51*T(2); d67 = T(1)/d67;
complex<T> t29 = t54*t56; 
complex<T> t75 = -t108; 
complex<T> t87 = t53*t54; 
complex<T> t88 = t55*t56; 
complex<T> t159 = d3*t19; 
complex<T> t185 = t50*t64; 
complex<T> t192 = t50*t68; 
complex<T> t199 = t111*t180; 
complex<T> t207 = t18*t97; 
complex<T> d21 = spb23*spb34*t32; d21 = T(1)/d21;
complex<T> d22 = s245*t28*t90; d22 = T(1)/d22;
complex<T> d28 = spb34*t32*T(2); d28 = T(1)/d28;
complex<T> d29 = t27*t28*t90; d29 = T(1)/d29;
complex<T> d35 = s245*spa25*t28; d35 = T(1)/d35;
complex<T> d38 = spb23*spb34*t32*T(2); d38 = T(1)/d38;
complex<T> d39 = spa45*t27*t28; d39 = T(1)/d39;
complex<T> d45 = spb23*t32*T(2); d45 = T(1)/d45;
complex<T> d51 = spa25*t27*t28; d51 = T(1)/d51;
complex<T> d56 = t28*t90*T(2); d56 = T(1)/d56;
complex<T> d60 = t32*T(4); d60 = T(1)/d60;
complex<T> d63 = s245*t28*T(4); d63 = T(1)/d63;
complex<T> t2 = d1*spb35*t107*t187 + t185*t186*t44 - d1*spa23*t44*t50 + spa23*t159*T(11); 
complex<T> t5 = t107*t137*t195 + spa34*(-(d18*t43*t50) + t192*t43*T(2) + d19*t20*T(11)); 
complex<T> t8 = -t110 + t141*t173 + t125*t137*t195 - d6*s245*t43*t45 + d5*t27*t43*t45 + d7*spa24*t27*t43*t45 + d8*spa24*(t173 - t43*t45) + d18*spa34*t43*t50 + d13*s25*t18*T(2) - spa34*t192*t43*T(2) - d4*t43*t45*T(4) + d15*spa24*t18*T(11) - d19*spa34*t20*T(11); 
complex<T> t9 = d8*spa24*t107*t25 + d9*spa24*t107*t25 + d10*s245*t109*t44 + d12*s245*spa24*t109*t44 + d6*s245*t43*t45 + d8*spa24*t43*t45 + t141*t43*t45 + d11*s245*t44*t46 + d9*spa24*t44*t46 + t141*t44*t46 + d5*s245*t43*t76 + d7*s245*spa24*t43*t76 + t110*T(2) - d13*s25*t18*T(2) - t207*T(2) - d4*t173*T(8) - d15*spa24*t18*T(11) - d16*spa24*t18*T(11); 
complex<T> t10 = -t110 + d9*spa24*t173 + t141*t173 + d1*spb35*t125*t187 - d11*s245*t44*t46 - d9*spa24*t44*t46 + d10*t27*t44*t46 + d12*spa24*t27*t44*t46 + t207*T(2) - d4*t44*t46*T(4) + d16*spa24*t18*T(11) + spa23*(d1*t44*t50 - t185*t44*T(2) - t159*T(11)); 
complex<T> t121 = d38*s25; 
complex<T> t136 = d38*s45; 
complex<T> d20 = s345*spa45*t87; d20 = T(1)/d20;
complex<T> d25 = s235*spa25*t88; d25 = T(1)/d25;
complex<T> d26 = s245*t25*t29; d26 = T(1)/d26;
complex<T> d27 = s345*spa45*t87*T(2); d27 = T(1)/d27;
complex<T> d32 = s235*spa25*t88*T(2); d32 = T(1)/d32;
complex<T> d33 = t25*t27*t29; d33 = T(1)/d33;
complex<T> d37 = s245*spb25*t29; d37 = T(1)/d37;
complex<T> d41 = s235*t88*T(2); d41 = T(1)/d41;
complex<T> d42 = spb45*t27*t29; d42 = T(1)/d42;
complex<T> d44 = s235*t88; d44 = T(1)/d44;
complex<T> d48 = s345*t87; d48 = T(1)/d48;
complex<T> d50 = s345*t87*T(2); d50 = T(1)/d50;
complex<T> d53 = spb25*t27*t29; d53 = T(1)/d53;
complex<T> d54 = spa45*t87*T(2); d54 = T(1)/d54;
complex<T> d57 = t25*t29*T(2); d57 = T(1)/d57;
complex<T> d59 = spa25*t88*T(2); d59 = T(1)/d59;
complex<T> d65 = s245*t29*T(4); d65 = T(1)/d65;
complex<T> t1 = d63*t108*t25 + t110*t58 + d61*t173*t58 + d62*t199*t58 + d62*spa23*t202*t58 + d64*t144*t45*t58 + d65*t24*t90; 
complex<T> t3 = s34*(d59*t163 + d58*t22); 
complex<T> t4 = d25*s14*t163 + d44*spb25*t163 + d23*s14*t22 + d43*spa25*t22; 
complex<T> t6 = s23*(d54*t151 + d55*t23); 
complex<T> t7 = d20*s12*t151 + d48*spb45*t151 + d24*s12*t23 + d49*spa45*t23; 
complex<T> t11 = d22*s15*t108 + d20*s15*t151 + d25*s15*t163 - d23*mH2*t22 + d23*s15*t22 - d24*mH2*t23 + d24*s15*t23 - d26*mH2*t24 + d26*s15*t24 - d21*mH2*t48 + d21*s15*t48 - d22*t47*cube(mH2) - d25*t73*cube(mH2) - d20*t74*cube(mH2); 
complex<T> t12 = d27*s23*t151 + d32*s23*t163 + d30*spa23*t22 - d31*s23*t23 + d33*s23*t24 - d28*spa23*t48 + d29*s23*t75; 
complex<T> t13 = -(s24*(d20*t151 + d25*t163 + d23*t22 + d24*t23 + d26*t24 - d22*t75)); 
complex<T> t14 = d27*s34*t151 + d32*s34*t163 - d46*s34*t22 + d47*spa34*t23 + d33*s34*t24 - d45*spa34*t48 + d29*s34*t75; 
complex<T> t15 = d39*spb25*t108 - d40*spa25*t22 - d42*spa25*t24 + t121*t48 - d41*spb25*t163*T(3) + s25*(d27*t151 + d34*t186*t202 - d31*t23 + d14*t173*T(2) + t144*t184*T(2) + d34*t199*T(2) + t110*T(4)); 
complex<T> t16 = d22*s13*t108 + d35*spb45*t108 + d34*s45*t186*t202 + d26*s13*t24 + d37*spa45*t24 + spb25*t180*t57*t62 + spa23*spb45*t45*t57*t62 + d36*s13*t144*t76 - d14*s13*t173*T(2) + d14*s45*t173*T(2) + s45*t144*t184*T(2) + d34*s45*t199*T(2) - s13*t110*T(4) + s45*t110*T(4); 
complex<T> t17 = d32*s45*t163 - d46*s45*t22 - d52*spa45*t23 + t136*t48 + d51*spb45*t75 - d50*spb45*t151*T(3) - d53*spa45*t24*T(3); 
complex<T> t36 = s45*(d59*t163 + d58*t22); 
complex<T> t37 = -(d41*s23*spb25*t163) + d66*spa23*spa25*t22; 
complex<T> t38 = s25*(d54*t151 + d55*t23); 
complex<T> t39 = -(d50*s34*spb45*t151) + d67*spa34*spa45*t23; 
complex<T> t40 = s23*(d56*t108 + d57*t24); 
complex<T> t41 = -((s15 - s24)*(d20*t151 + d25*t163 + d23*t22 + d24*t23 + d26*t24)) + d22*(s24*t108 + s15*t75); 
complex<T> t72 = s34*(d56*t108 + d57*t24); 
complex<T> co1 = -t100; 
complex<T> co2 = -t101; 
complex<T> co3 = t121*t89; 
complex<T> co4 = t136*t89; 
complex<T> co5 = d60*t125*t48; 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t10*(*CI_users[2]->get_value(mc,ind,mu)) + t5*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t11*(*CI_users[7]->get_value(mc,ind,mu)) + t12*(*CI_users[8]->get_value(mc,ind,mu)) + t13*(*CI_users[9]->get_value(mc,ind,mu)) + t16*(*CI_users[10]->get_value(mc,ind,mu)) + t15*(*CI_users[11]->get_value(mc,ind,mu)) + t4*(*CI_users[12]->get_value(mc,ind,mu)) + t41*(*CI_users[13]->get_value(mc,ind,mu)) + t14*(*CI_users[14]->get_value(mc,ind,mu)) + t7*(*CI_users[15]->get_value(mc,ind,mu)) + t17*(*CI_users[16]->get_value(mc,ind,mu)) + t6*(*CI_users[17]->get_value(mc,ind,mu)) + t38*(*CI_users[18]->get_value(mc,ind,mu)) + t40*(*CI_users[19]->get_value(mc,ind,mu)) + t72*(*CI_users[20]->get_value(mc,ind,mu)) + t3*(*CI_users[21]->get_value(mc,ind,mu)) + t36*(*CI_users[22]->get_value(mc,ind,mu)) + co3*(*CI_users[23]->get_value(mc,ind,mu)) + co4*(*CI_users[24]->get_value(mc,ind,mu)) + co5*(*CI_users[25]->get_value(mc,ind,mu)) + t1*(*CI_users[26]->get_value(mc,ind,mu)) + t37*(*CI_users[27]->get_value(mc,ind,mu)) + t39*(*CI_users[28]->get_value(mc,ind,mu)) + co5*(*CI_users[29]->get_value(mc,ind,mu)) + t1*(*CI_users[30]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmmmp_nf_wCI::\
C4g1ph_phmmmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmmmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmmmp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:12:47 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s35 = -(spa35*spb35);
complex<T> s235 = SS(2,3,5);
complex<T> s25 = S(2,5);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s45 = S(4,5);
complex<T> s24 = -(spa24*spb24);
complex<T> s34 = -(spa34*spb34);
complex<T> s345 = SS(3,4,5);
complex<T> t8 = square(spa23*spb25 - spa34*spb45); 
complex<T> t9 = -(spa23*T(2)); 
complex<T> t10 = -(s25*s45); 
complex<T> t11 = square(spa24*spb25 + spa34*spb35); 
complex<T> t12 = square(spa23*spb35 + spa24*spb45); 
complex<T> t13 = square(spb24); 
complex<T> t20 = square(spb25); 
complex<T> t21 = square(spb45); 
complex<T> t22 = square(spa23); 
complex<T> t23 = square(spa34); 
complex<T> t24 = cube(spb24); 
complex<T> t25 = square(spb35); 
complex<T> t27 = square(spa35); 
complex<T> t47 = spb24*T(9); 
complex<T> t57 = spa34*spb25; 
complex<T> t58 = spa23*spb45; 
complex<T> t107 = spa23*spa24; 
complex<T> t109 = spa24*spa34; 
complex<T> d1 = spb23*square(s23 + s35)*T(9); d1 = T(1)/d1;
complex<T> d2 = spb23*cube(s23 + s35)*T(9); d2 = T(1)/d2;
complex<T> d3 = s235*(s23 + s35)*spb23*T(9); d3 = T(1)/d3;
complex<T> d17 = spb34*cube(s34 + s35)*T(9); d17 = T(1)/d17;
complex<T> d18 = spb34*square(s34 + s35)*T(9); d18 = T(1)/d18;
complex<T> d19 = s345*(s34 + s35)*spb34*T(9); d19 = T(1)/d19;
complex<T> d21 = s245*square(square(spb24))*T(3); d21 = T(1)/d21;
complex<T> d24 = s245*square(square(spb24))*T(6); d24 = T(1)/d24;
complex<T> t30 = d2*s235; 
complex<T> t32 = d17*s345; 
complex<T> t40 = -t57; 
complex<T> t62 = d1*spb35; 
complex<T> t63 = d21*t27; 
complex<T> t74 = d18*spb35; 
complex<T> t76 = t20*t22; 
complex<T> t77 = t21*t23; 
complex<T> t83 = t57*t58; 
complex<T> t89 = -(t21*T(2)); 
complex<T> t104 = spb45*t20; 
complex<T> t105 = t23*t25; 
complex<T> t127 = t57*t9; 
complex<T> d4 = spa24*t24*T(3); d4 = T(1)/d4;
complex<T> d5 = (s24 + s25)*spa24*t24*T(3); d5 = T(1)/d5;
complex<T> d6 = t13*square(s24 + s25)*T(3); d6 = T(1)/d6;
complex<T> d7 = t47*cube(s24 + s25); d7 = T(1)/d7;
complex<T> d8 = t47*square(s24 + s25); d8 = T(1)/d8;
complex<T> d9 = t47*square(s24 + s45); d9 = T(1)/d9;
complex<T> d10 = (s24 + s45)*spa24*t24*T(3); d10 = T(1)/d10;
complex<T> d11 = t13*square(s24 + s45)*T(3); d11 = T(1)/d11;
complex<T> d12 = t47*cube(s24 + s45); d12 = T(1)/d12;
complex<T> d13 = s245*spa24*t24*T(3); d13 = T(1)/d13;
complex<T> d14 = s245*t13*T(3); d14 = T(1)/d14;
complex<T> d15 = s245*(s24 + s25)*t47; d15 = T(1)/d15;
complex<T> d16 = s245*(s24 + s45)*t47; d16 = T(1)/d16;
complex<T> d20 = s245*t24*T(3); d20 = T(1)/d20;
complex<T> d22 = s245*t13*T(6); d22 = T(1)/d22;
complex<T> d23 = s245*t24*T(6); d23 = T(1)/d23;
complex<T> d25 = s245*t13*T(12); d25 = T(1)/d25;
complex<T> t2 = d1*spa23*t105 + t107*t57*t62 + d3*t11*t9 + t105*t30*t9; 
complex<T> t29 = d20*s13; 
complex<T> t39 = d23*t10; 
complex<T> t50 = d20*spa35; 
complex<T> t69 = -(d4*T(4)); 
complex<T> t85 = d12*spa24; 
complex<T> t91 = t20*t63; 
complex<T> t98 = d14*spb45; 
complex<T> t99 = t22*t32; 
complex<T> t103 = d7*spa24; 
complex<T> t1 = d24*t10*t20*t21*t27 + spa35*t21*t39*t57 + spa35*t20*t39*t58 + d25*t10*t8 + d22*t10*t83; 
complex<T> t3 = d18*spa34*t22*t25 + t109*t58*t74 - d19*spa34*t12*T(2) - spa34*t25*t99*T(2); 
complex<T> t4 = s45*t50*t57*t89 + s45*t104*t50*t9 + s45*t89*t91 + s45*t127*t98 + spa35*t29*(t21*t57 + t20*t58)*T(2) + s13*t21*t91*T(2) + d14*(s13*t8 - s45*t8 + s13*t83*T(2)); 
complex<T> t5 = -(d18*spa34*t22*t25) + d8*spa24*t40*t58 - t109*t58*t74 + d6*s245*t76 + d8*spa24*t76 + d14*t8 + t69*t83 + d19*spa34*t12*T(2) - d5*s245*t76*T(2) - s245*t103*t76*T(2) - d13*s25*t8*T(2) - d15*spa24*t8*T(2) + spa34*t25*t99*T(2) + d4*t76*T(4); 
complex<T> t6 = -(d1*spa23*t105) + d9*spa24*t40*t58 + t107*t40*t62 + d11*s245*t77 + d9*spa24*t77 + d14*t8 + t69*t83 + d3*spa23*t11*T(2) + spa23*t105*t30*T(2) - d10*s245*t77*T(2) - d13*s45*t8*T(2) - d16*spa24*t8*T(2) - s245*t77*t85*T(2) + d4*t77*T(4); 
complex<T> t7 = -(d6*s245*t76) - d8*spa24*t76 + t69*t76 - d11*s245*t77 - d9*spa24*t77 + t69*t77 + d8*spa24*t83 + d9*spa24*t83 + d5*s245*t76*T(2) + s245*t103*t76*T(2) + d10*s245*t77*T(2) - d14*t8*T(2) + d13*s25*t8*T(2) + d13*s45*t8*T(2) + d15*spa24*t8*T(2) + d16*spa24*t8*T(2) + s245*t77*t85*T(2) + d4*t83*T(8); 
complex<T> t18 = s25*(-(d14*t8) + t50*t57*t89 + t104*t50*t9 + t89*t91 + t127*t98); 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t7*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t5*(*CI_users[4]->get_value(mc,ind,mu)) + t4*(*CI_users[5]->get_value(mc,ind,mu)) + t18*(*CI_users[6]->get_value(mc,ind,mu)) + t1*(*CI_users[7]->get_value(mc,ind,mu)) + t1*(*CI_users[8]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phmmmm_G_wCI::\
C4g1ph_phmmmm_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phmmmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, m, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmmmm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:12:48 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s14 = S(1,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t2 = square(mH2); 
complex<T> t7 = cube(mH2); 
complex<T> d1 = spb23*spb25*spb34*spb45; d1 = T(1)/d1;
complex<T> d2 = spb23*spb25*spb34; d2 = T(1)/d2;
complex<T> d3 = spb23*spb34*spb45; d3 = T(1)/d3;
complex<T> d4 = spb23*spb25*spb34*T(2); d4 = T(1)/d4;
complex<T> d5 = spb23*spb25*spb34*spb45*T(2); d5 = T(1)/d5;
complex<T> d6 = spb23*spb34*spb45*T(2); d6 = T(1)/d6;
complex<T> d7 = spb25*spb34*spb45*T(2); d7 = T(1)/d7;
complex<T> d8 = spb23*spb25*spb45*T(2); d8 = T(1)/d8;
complex<T> d9 = spb23*spb34*T(2); d9 = T(1)/d9;
complex<T> d10 = spb34*spb45*T(2); d10 = T(1)/d10;
complex<T> d11 = spb25*spb45*T(2); d11 = T(1)/d11;
complex<T> d12 = spb23*spb25*T(2); d12 = T(1)/d12;
complex<T> t3 = -(d1*T(3)); 
complex<T> t4 = d2*spa45; 
complex<T> t5 = d5*s234; 
complex<T> t6 = d3*spa25; 
complex<T> t10 = d1*T(3); 
complex<T> t11 = d5*s245; 
complex<T> t13 = d1*t2; 
complex<T> t18 = s235*t2; 
complex<T> t20 = s24*t2; 
complex<T> t22 = s345*t2; 
complex<T> t1 = s15*t13*T(4) - d1*t7*T(4); 
complex<T> t12 = t10*t20 + s15*t2*t3; 
complex<T> t16 = t2*t4; 
complex<T> t21 = t18*t5 + d7*spa23*t7; 
complex<T> t23 = t11*t18 + d6*spa25*t7; 
complex<T> t24 = t2*t6; 
complex<T> t25 = t22*t5 + d8*spa34*t7; 
complex<T> t26 = t11*t22 + d4*spa45*t7; 
complex<T> t15 = s12*t13 + t16; 
complex<T> t17 = s13*t13 + t16; 
complex<T> t19 = s14*t13 + t24; 
complex<T> co1 = t20*t3; 
complex<T> co2 = -t24; 
complex<T> co3 = -(t16*T(2)); 
complex<T> co4 = d9*spa25*spa45*t2; 
complex<T> co5 = d10*spa23*spa25*t2; 
complex<T> co6 = d11*spa23*spa34*t2; 
complex<T> co7 = d12*spa34*spa45*t2; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + t17*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + t19*(*CI_users[4]->get_value(mc,ind,mu)) + t12*(*CI_users[5]->get_value(mc,ind,mu)) + t15*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + t26*(*CI_users[8]->get_value(mc,ind,mu)) + t23*(*CI_users[9]->get_value(mc,ind,mu)) + t21*(*CI_users[10]->get_value(mc,ind,mu)) + t25*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + co5*(*CI_users[13]->get_value(mc,ind,mu)) + co6*(*CI_users[14]->get_value(mc,ind,mu)) + co7*(*CI_users[15]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C4g1ph_phmmmm_nf_wCI::\
C4g1ph_phmmmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phmmmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phmmmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C4g1ph_phdpppp_G_wCI::\
C4g1ph_phdpppp_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpppp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpppp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:12:50 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s14 = S(1,4);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t2 = square(mH2); 
complex<T> t7 = cube(mH2); 
complex<T> d1 = spa23*spa25*spa34*spa45; d1 = T(1)/d1;
complex<T> d2 = spa23*spa25*spa34; d2 = T(1)/d2;
complex<T> d3 = spa23*spa34*spa45; d3 = T(1)/d3;
complex<T> d4 = spa23*spa25*spa34*spa45*T(2); d4 = T(1)/d4;
complex<T> d5 = spa23*spa25*spa34*T(2); d5 = T(1)/d5;
complex<T> d6 = spa23*spa34*spa45*T(2); d6 = T(1)/d6;
complex<T> d7 = spa25*spa34*spa45*T(2); d7 = T(1)/d7;
complex<T> d8 = spa23*spa25*spa45*T(2); d8 = T(1)/d8;
complex<T> d9 = spa23*spa34*T(2); d9 = T(1)/d9;
complex<T> d10 = spa34*spa45*T(2); d10 = T(1)/d10;
complex<T> d11 = spa25*spa45*T(2); d11 = T(1)/d11;
complex<T> d12 = spa23*spa25*T(2); d12 = T(1)/d12;
complex<T> t3 = -(d1*T(3)); 
complex<T> t4 = d2*spb45; 
complex<T> t5 = d4*s234; 
complex<T> t6 = d3*spb25; 
complex<T> t10 = d1*T(3); 
complex<T> t11 = d4*s245; 
complex<T> t13 = d1*t2; 
complex<T> t18 = s235*t2; 
complex<T> t20 = s24*t2; 
complex<T> t22 = s345*t2; 
complex<T> t1 = s15*t13*T(4) - d1*t7*T(4); 
complex<T> t12 = t10*t20 + s15*t2*t3; 
complex<T> t16 = t2*t4; 
complex<T> t21 = t18*t5 + d7*spb23*t7; 
complex<T> t23 = t11*t18 + d6*spb25*t7; 
complex<T> t24 = t2*t6; 
complex<T> t25 = t22*t5 + d8*spb34*t7; 
complex<T> t26 = t11*t22 + d5*spb45*t7; 
complex<T> t15 = s12*t13 + t16; 
complex<T> t17 = s13*t13 + t16; 
complex<T> t19 = s14*t13 + t24; 
complex<T> co1 = t20*t3; 
complex<T> co2 = -t24; 
complex<T> co3 = -(t16*T(2)); 
complex<T> co4 = d9*spb25*spb45*t2; 
complex<T> co5 = d10*spb23*spb25*t2; 
complex<T> co6 = d11*spb23*spb34*t2; 
complex<T> co7 = d12*spb34*spb45*t2; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + co1*(*CI_users[1]->get_value(mc,ind,mu)) + t17*(*CI_users[2]->get_value(mc,ind,mu)) + co2*(*CI_users[3]->get_value(mc,ind,mu)) + t19*(*CI_users[4]->get_value(mc,ind,mu)) + t12*(*CI_users[5]->get_value(mc,ind,mu)) + t15*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + t26*(*CI_users[8]->get_value(mc,ind,mu)) + t23*(*CI_users[9]->get_value(mc,ind,mu)) + t21*(*CI_users[10]->get_value(mc,ind,mu)) + t25*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + co5*(*CI_users[13]->get_value(mc,ind,mu)) + co6*(*CI_users[14]->get_value(mc,ind,mu)) + co7*(*CI_users[15]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C4g1ph_phdpppp_nf_wCI::\
C4g1ph_phdpppp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdpppp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpppp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C4g1ph_phdmppp_G_wCI::\
C4g1ph_phdmppp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c45));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c25, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdmppp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmppp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:16:12 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> s14 = S(1,4);
complex<T> s235 = SS(2,3,5);
complex<T> s35 = -(spa35*spb35);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s345 = SS(3,4,5);
complex<T> s15 = S(1,5);
complex<T> s12 = S(1,2);
complex<T> t18 = square(spa23*spb34 - spa25*spb45); 
complex<T> t23 = square(square(spa23*spb34 - spa25*spb45)); 
complex<T> t24 = spa23*spa25; 
complex<T> t27 = s235*T(2); 
complex<T> t42 = square(mH2); 
complex<T> t43 = square(spb45); 
complex<T> t44 = square(spb34); 
complex<T> t45 = square(spa23); 
complex<T> t46 = square(spa25); 
complex<T> t47 = cube(s345); 
complex<T> t48 = square(square(spb35)); 
complex<T> t49 = square(spa24); 
complex<T> t51 = spa24*spb23 + spa45*spb35; 
complex<T> t52 = spa24*spb25 + spa34*spb35; 
complex<T> t53 = spa35*spb23 + spa45*spb24; 
complex<T> t54 = spa34*spb24 + spa35*spb25; 
complex<T> t55 = spa25*spb24 + spa35*spb34; 
complex<T> t56 = spa23*spb24 + spa35*spb45; 
complex<T> t57 = -(spb24*T(2)); 
complex<T> t58 = s23*s25; 
complex<T> t65 = s235*T(3); 
complex<T> t72 = cube(spb45); 
complex<T> t73 = cube(spb34); 
complex<T> t78 = square(spb24); 
complex<T> t87 = spb23*spb25; 
complex<T> t99 = coeff3mass<T>(mc,ind,2,3,4,5); 
complex<T> t100 = coeff3mass<T>(mc,ind,2,5,4,3); 
complex<T> t111 = square(square(s345)); 
complex<T> t115 = spa23*spa34; 
complex<T> t125 = spb34*spb45; 
complex<T> t142 = s23*T(2); 
complex<T> t181 = spa23*spb45; 
complex<T> t191 = spa25*spb35; 
complex<T> t195 = spa23*spb35; 
complex<T> d1 = spb35*cube(spa35); d1 = T(1)/d1;
complex<T> d2 = spa34*square(s24 + s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = spa35*square(s25 + s35)*T(3); d3 = T(1)/d3;
complex<T> d4 = square(s25 + s35)*square(spa35); d4 = T(1)/d4;
complex<T> d5 = spa34*cube(s24 + s34)*T(3); d5 = T(1)/d5;
complex<T> d6 = (s25 + s35)*spb35*cube(spa35); d6 = T(1)/d6;
complex<T> d7 = spa35*cube(s25 + s35)*T(3); d7 = T(1)/d7;
complex<T> d8 = s234*(s24 + s34)*spa34*T(3); d8 = T(1)/d8;
complex<T> d9 = s235*square(spa35); d9 = T(1)/d9;
complex<T> d10 = s235*spb35*cube(spa35); d10 = T(1)/d10;
complex<T> d12 = square(s23 + s35)*square(spa35); d12 = T(1)/d12;
complex<T> d13 = (s23 + s35)*spb35*cube(spa35); d13 = T(1)/d13;
complex<T> d14 = spa35*cube(s23 + s35)*T(3); d14 = T(1)/d14;
complex<T> d15 = spa35*square(s23 + s35)*T(3); d15 = T(1)/d15;
complex<T> d17 = spa45*cube(s24 + s45)*T(3); d17 = T(1)/d17;
complex<T> d18 = spa45*square(s24 + s45)*T(3); d18 = T(1)/d18;
complex<T> d19 = s245*(s24 + s45)*spa45*T(3); d19 = T(1)/d19;
complex<T> d27 = s235*square(square(spa35)); d27 = T(1)/d27;
complex<T> d29 = s235*cube(spa35); d29 = T(1)/d29;
complex<T> t19 = square(spa24*spb34 + t191); 
complex<T> t20 = square(spa24*spb45 + t195); 
complex<T> t22 = cube(spa24*spb34 + t191); 
complex<T> t26 = cube(spa24*spb45 + t195); 
complex<T> t28 = t51*t52; 
complex<T> t29 = t55*t56; 
complex<T> t32 = t53*t54; 
complex<T> t62 = d29*s14; 
complex<T> t64 = d5*s234; 
complex<T> t67 = d17*s245; 
complex<T> t75 = -(t45*T(2)); 
complex<T> t83 = s245*t51; 
complex<T> t89 = t53*t55; 
complex<T> t90 = t54*t56; 
complex<T> t95 = d10*s25; 
complex<T> t107 = -t125; 
complex<T> t108 = t42*t48; 
complex<T> t109 = -(t46*T(2)); 
complex<T> t110 = d9*t18; 
complex<T> t112 = d29*spb24; 
complex<T> t130 = spa25*t45; 
complex<T> t133 = d27*t78; 
complex<T> t137 = d2*spa24; 
complex<T> t141 = t125*t24; 
complex<T> t143 = s45*t42; 
complex<T> t155 = d1*T(4); 
complex<T> t163 = t42*t72; 
complex<T> t180 = t42*t73; 
complex<T> t188 = spb34*t43; 
complex<T> t193 = spb45*t44; 
complex<T> d11 = (s25 + s35)*spa35*t65; d11 = T(1)/d11;
complex<T> d16 = (s23 + s35)*spa35*t65; d16 = T(1)/d16;
complex<T> d24 = s234*t115*t52; d24 = T(1)/d24;
complex<T> d33 = s234*spa34*t52*T(2); d33 = T(1)/d33;
complex<T> d40 = s234*t115*t52*T(2); d40 = T(1)/d40;
complex<T> d47 = s234*spa23*t52*T(2); d47 = T(1)/d47;
complex<T> d52 = spa25*spa45*t51*T(2); d52 = T(1)/d52;
complex<T> d57 = t115*t52*T(2); d57 = T(1)/d57;
complex<T> d58 = t27*square(square(spa35)); d58 = T(1)/d58;
complex<T> d59 = t27*cube(spa35); d59 = T(1)/d59;
complex<T> d61 = t27*square(spa35); d61 = T(1)/d61;
complex<T> d64 = s234*t52*T(2); d64 = T(1)/d64;
complex<T> t74 = -t108; 
complex<T> t98 = d8*t20; 
complex<T> t138 = d19*t19; 
complex<T> t140 = d59*t58; 
complex<T> t148 = s25*t112; 
complex<T> t157 = t133*t46; 
complex<T> t173 = -t180; 
complex<T> t206 = t18*t95; 
complex<T> d20 = spa34*spa45*t32; d20 = T(1)/d20;
complex<T> d21 = s234*spb23*t89; d21 = T(1)/d21;
complex<T> d22 = spa25*spa45*t83; d22 = T(1)/d22;
complex<T> d23 = s235*t28*t87; d23 = T(1)/d23;
complex<T> d25 = s245*spb25*t90; d25 = T(1)/d25;
complex<T> d26 = s235*t24*t29; d26 = T(1)/d26;
complex<T> d28 = spa34*spa45*t32*T(2); d28 = T(1)/d28;
complex<T> d30 = s234*t89*T(2); d30 = T(1)/d30;
complex<T> d31 = spa25*spa45*t83*T(2); d31 = T(1)/d31;
complex<T> d32 = spb25*t27*t28; d32 = T(1)/d32;
complex<T> d34 = s245*spb25*t90*T(2); d34 = T(1)/d34;
complex<T> d35 = spa25*t27*t29; d35 = T(1)/d35;
complex<T> d36 = spa25*t83; d36 = T(1)/d36;
complex<T> d37 = s234*spb23*t89*T(2); d37 = T(1)/d37;
complex<T> d38 = spa45*t83*T(2); d38 = T(1)/d38;
complex<T> d39 = spb23*t27*t28; d39 = T(1)/d39;
complex<T> d41 = s245*t90*T(2); d41 = T(1)/d41;
complex<T> d42 = spa23*t27*t29; d42 = T(1)/d42;
complex<T> d43 = s235*spb23*t28; d43 = T(1)/d43;
complex<T> d44 = s235*spa23*t29; d44 = T(1)/d44;
complex<T> d45 = spa45*t32*T(2); d45 = T(1)/d45;
complex<T> d46 = t27*t28*t87; d46 = T(1)/d46;
complex<T> d48 = t24*t27*t29; d48 = T(1)/d48;
complex<T> d49 = spa34*t32; d49 = T(1)/d49;
complex<T> d50 = spa34*t32*T(2); d50 = T(1)/d50;
complex<T> d51 = spa25*t83*T(2); d51 = T(1)/d51;
complex<T> d53 = spb25*t90*T(2); d53 = T(1)/d53;
complex<T> d54 = t28*t87*T(2); d54 = T(1)/d54;
complex<T> d55 = t24*t29*T(2); d55 = T(1)/d55;
complex<T> d56 = spb23*t89*T(2); d56 = T(1)/d56;
complex<T> d60 = s235*t28*T(4); d60 = T(1)/d60;
complex<T> d62 = s235*t29*T(4); d62 = T(1)/d62;
complex<T> d63 = t32*T(4); d63 = T(1)/d63;
complex<T> d65 = t83*T(2); d65 = T(1)/d65;
complex<T> t1 = d60*t108*t24 + spb24*t140*(spb34*t130 + t181*t46) + t110*t58 + d61*t141*t58 + d58*t45*t46*t58*t78 + d62*t23*t87; 
complex<T> t2 = d18*spa24*t107*t191 - d18*t193*t49 + t193*t49*t67*T(2) + spb45*t138*T(11); 
complex<T> t3 = s23*(d53*t163 + d52*t22); 
complex<T> t5 = t107*t137*t195 - d2*t188*t49 + t188*t49*t64*T(2) + spb34*t98*T(11); 
complex<T> t6 = s25*(d56*t180 + d57*t26); 
complex<T> t7 = -t110 + t141*t155 + t125*t137*t195 - d4*s235*t43*t46 + d6*t27*t43*t46 + d7*spb35*t27*t43*t46 + d3*spb35*(t141 - t43*t46) + d2*t188*t49 + t206*T(2) - t188*t49*t64*T(2) - d1*t43*t46*T(4) + d11*spb35*t18*T(11) - spb34*t98*T(11); 
complex<T> t8 = d15*spb35*t107*t24 + d3*spb35*t107*t24 + d6*s235*t109*t43 + d7*s235*spb35*t109*t43 + d12*s235*t44*t45 + d15*spb35*t44*t45 + t155*t44*t45 + d4*s235*t43*t46 + d3*spb35*t43*t46 + t155*t43*t46 + d13*s235*t44*t75 + d14*s235*spb35*t44*t75 + t110*T(2) - d10*s23*t18*T(2) - t206*T(2) - d1*t141*T(8) - d11*spb35*t18*T(11) - d16*spb35*t18*T(11); 
complex<T> t9 = -t110 + d15*spb35*t141 + t141*t155 + d10*t142*t18 + d18*spa24*t125*t191 - d12*s235*t44*t45 - d15*spb35*t44*t45 + d13*t27*t44*t45 + d14*spb35*t27*t44*t45 + d18*t193*t49 - t193*t49*t67*T(2) - d1*t44*t45*T(4) - spb45*t138*T(11) + d16*spb35*t18*T(11); 
complex<T> t12 = d34*s34*t163 + d37*s34*t180 - d31*s34*t22 + d48*s34*t23 + d47*spb34*t26 - d45*spb34*t47 + d46*s34*t74; 
complex<T> t13 = s34*(d54*t108 + d55*t23); 
complex<T> t14 = -(d51*spb45*t22) + d48*s45*t23 - d40*s45*t26 + d37*t143*t73 + d46*s45*t74 - d50*spb45*t47*T(3) + d34*t143*t72*T(3); 
complex<T> t17 = d23*s14*t108 + d43*spa25*t108 + d26*s14*t23 + d44*spb25*t23 + spb34*t130*t57*t62 + t181*t46*t57*t62 + s14*t157*t75 - d9*s14*t141*T(2) + d9*s25*t141*T(2) + spb34*t130*t148*T(2) + s25*t157*t45*T(2) + t148*t181*t46*T(2) - s14*t110*T(4) + s25*t110*T(4); 
complex<T> t35 = (d20*s12 + d49*spb45)*t47; 
complex<T> t37 = s34*(d53*t163 + d52*t22); 
complex<T> t38 = d65*spb25*spb45*t22 - d41*spa25*t143*t72; 
complex<T> t40 = s45*(d54*t108 + d55*t23); 
complex<T> t93 = d28*s23; 
complex<T> t101 = -(d25*t72); 
complex<T> t106 = d57*s45*t26 + d56*t143*t73; 
complex<T> t122 = d28*s25; 
complex<T> t190 = d30*spa23; 
complex<T> t4 = t101*t143 + d25*s13*t163 + d22*s13*t22 + d36*spb45*t22; 
complex<T> t10 = -(mH2*(d22*t22 + d26*t23 + d24*t26 + d20*t47)) + s15*(d23*t108 + d25*t163 + d21*t180 + d22*t22 + d26*t23 + d24*t26 + d20*t47) + (t101 - d23*t48 - d21*t73)*cube(mH2); 
complex<T> t11 = -(s24*(d22*t22 + d26*t23 - t101*t42 + d20*t47 - d23*t74)); 
complex<T> t15 = d32*spa23*t108 + spb34*t112*t130*t142 + d9*t141*t142 + t173*t190 - d35*spb23*t23 + d33*spb23*t26 + t142*t157*t45 + t112*t142*t181*t46 + t47*t93 + s23*(d34*t163 - d31*t22 + t110*T(4)); 
complex<T> t16 = -(d41*spa25*t163) + d37*s25*t180 + d38*spb25*t22 - d40*s25*t26 + t122*t47 + d39*spa25*t74 - d42*spb25*t23*T(3); 
complex<T> t39 = d23*s24*t108 + d25*s24*t163 - d22*s15*t22 + d22*s24*t22 - d26*s15*t23 + d26*s24*t23 + s15*t101*t42 - d20*s15*t47 + d20*s24*t47 + d23*s15*t74; 
complex<T> t41 = s34*t173*t190 + d64*spb23*spb34*t26; 
complex<T> co1 = -t99; 
complex<T> co2 = -t100; 
complex<T> co3 = t111*t93; 
complex<T> co4 = t111*t122; 
complex<T> co5 = d63*t125*t47; 
SeriesC<T> result = t7*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t8*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t9*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t10*(*CI_users[7]->get_value(mc,ind,mu)) + t15*(*CI_users[8]->get_value(mc,ind,mu)) + t11*(*CI_users[9]->get_value(mc,ind,mu)) + t4*(*CI_users[10]->get_value(mc,ind,mu)) + t16*(*CI_users[11]->get_value(mc,ind,mu)) + t17*(*CI_users[12]->get_value(mc,ind,mu)) + t39*(*CI_users[13]->get_value(mc,ind,mu)) + t12*(*CI_users[14]->get_value(mc,ind,mu)) + t35*(*CI_users[15]->get_value(mc,ind,mu)) + t14*(*CI_users[16]->get_value(mc,ind,mu)) + co3*(*CI_users[17]->get_value(mc,ind,mu)) + co4*(*CI_users[18]->get_value(mc,ind,mu)) + t3*(*CI_users[19]->get_value(mc,ind,mu)) + t37*(*CI_users[20]->get_value(mc,ind,mu)) + t13*(*CI_users[21]->get_value(mc,ind,mu)) + t40*(*CI_users[22]->get_value(mc,ind,mu)) + t6*(*CI_users[23]->get_value(mc,ind,mu)) + t106*(*CI_users[24]->get_value(mc,ind,mu)) + t1*(*CI_users[25]->get_value(mc,ind,mu)) + co5*(*CI_users[26]->get_value(mc,ind,mu)) + t41*(*CI_users[27]->get_value(mc,ind,mu)) + t38*(*CI_users[28]->get_value(mc,ind,mu)) + t1*(*CI_users[29]->get_value(mc,ind,mu)) + co5*(*CI_users[30]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdmppp_nf_wCI::\
C4g1ph_phdmppp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdmppp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmppp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:16:58 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> s24 = -(spa24*spb24);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s23 = S(2,3);
complex<T> s235 = SS(2,3,5);
complex<T> s14 = S(1,4);
complex<T> s25 = S(2,5);
complex<T> s35 = -(spa35*spb35);
complex<T> s45 = -(spa45*spb45);
complex<T> s245 = SS(2,4,5);
complex<T> t8 = square(spa23*spb34 - spa25*spb45); 
complex<T> t9 = -(spb34*T(2)); 
complex<T> t10 = -(s23*s25); 
complex<T> t11 = square(spa24*spb34 + spa25*spb35); 
complex<T> t12 = square(spa23*spb35 + spa24*spb45); 
complex<T> t13 = square(spa35); 
complex<T> t20 = square(spa23); 
complex<T> t21 = square(spa25); 
complex<T> t22 = square(spb34); 
complex<T> t23 = square(spb45); 
complex<T> t24 = square(spa24); 
complex<T> t25 = cube(spa35); 
complex<T> t27 = square(spb24); 
complex<T> t47 = spa35*T(9); 
complex<T> t57 = spa23*spb45; 
complex<T> t58 = spa25*spb34; 
complex<T> t106 = spb34*spb35; 
complex<T> t108 = spb35*spb45; 
complex<T> d2 = spa34*square(s24 + s34)*T(9); d2 = T(1)/d2;
complex<T> d5 = spa34*cube(s24 + s34)*T(9); d5 = T(1)/d5;
complex<T> d8 = s234*(s24 + s34)*spa34*T(9); d8 = T(1)/d8;
complex<T> d17 = spa45*cube(s24 + s45)*T(9); d17 = T(1)/d17;
complex<T> d18 = spa45*square(s24 + s45)*T(9); d18 = T(1)/d18;
complex<T> d19 = s245*(s24 + s45)*spa45*T(9); d19 = T(1)/d19;
complex<T> d20 = s235*square(square(spa35))*T(3); d20 = T(1)/d20;
complex<T> d22 = s235*square(square(spa35))*T(6); d22 = T(1)/d22;
complex<T> t30 = d5*s234; 
complex<T> t31 = d17*s245; 
complex<T> t40 = -t57; 
complex<T> t62 = d18*spa24; 
complex<T> t63 = d20*t27; 
complex<T> t67 = t20*t22; 
complex<T> t74 = d2*spa24; 
complex<T> t76 = t21*t23; 
complex<T> t82 = t57*t58; 
complex<T> t88 = -(t21*T(2)); 
complex<T> t89 = t22*t24; 
complex<T> t103 = spa25*t20; 
complex<T> t104 = t23*t24; 
complex<T> t123 = t57*t9; 
complex<T> d1 = spb35*t25*T(3); d1 = T(1)/d1;
complex<T> d3 = t47*square(s25 + s35); d3 = T(1)/d3;
complex<T> d4 = t13*square(s25 + s35)*T(3); d4 = T(1)/d4;
complex<T> d6 = (s25 + s35)*spb35*t25*T(3); d6 = T(1)/d6;
complex<T> d7 = t47*cube(s25 + s35); d7 = T(1)/d7;
complex<T> d9 = s235*t13*T(3); d9 = T(1)/d9;
complex<T> d10 = s235*spb35*t25*T(3); d10 = T(1)/d10;
complex<T> d11 = s235*(s25 + s35)*t47; d11 = T(1)/d11;
complex<T> d12 = t13*square(s23 + s35)*T(3); d12 = T(1)/d12;
complex<T> d13 = (s23 + s35)*spb35*t25*T(3); d13 = T(1)/d13;
complex<T> d14 = t47*cube(s23 + s35); d14 = T(1)/d14;
complex<T> d15 = t47*square(s23 + s35); d15 = T(1)/d15;
complex<T> d16 = s235*(s23 + s35)*t47; d16 = T(1)/d16;
complex<T> d21 = s235*t25*T(3); d21 = T(1)/d21;
complex<T> d23 = s235*t25*T(6); d23 = T(1)/d23;
complex<T> d24 = s235*t13*T(6); d24 = T(1)/d24;
complex<T> d25 = s235*t13*T(12); d25 = T(1)/d25;
complex<T> t2 = t108*t58*t62 - d19*spb45*t11*T(2) + spb45*t89*(d18 - t31*T(2)); 
complex<T> t3 = d2*spb34*t104 + t106*t57*t74 + d8*t12*t9 + t104*t30*t9; 
complex<T> t29 = d21*s14; 
complex<T> t39 = d23*t10; 
complex<T> t49 = d21*spb24; 
complex<T> t69 = -(d1*T(4)); 
complex<T> t80 = t20*t63; 
complex<T> t90 = d14*spb35; 
complex<T> t95 = d9*spa25; 
complex<T> t102 = d7*spb35; 
complex<T> t1 = d22*t10*t20*t21*t27 + spb24*t21*t39*t57 + spb24*t20*t39*t58 + d25*t10*t8 + d24*t10*t82; 
complex<T> t4 = s25*t49*t57*t88 + s25*t80*t88 + s25*t103*t49*t9 + s25*t123*t95 + spb24*t29*(t21*t57 + t20*t58)*T(2) + s14*t21*t80*T(2) + d9*(s14*t8 - s25*t8 + s14*t82*T(2)); 
complex<T> t5 = -(d2*spb34*t104) + d3*spb35*t40*t58 + t106*t40*t74 + d4*s235*t76 + d3*spb35*t76 + d9*t8 + t69*t82 + d8*spb34*t12*T(2) + spb34*t104*t30*T(2) - d6*s235*t76*T(2) - s235*t102*t76*T(2) - d10*s25*t8*T(2) - d11*spb35*t8*T(2) + d1*t76*T(4); 
complex<T> t6 = -(t108*t58*t62) + d12*s235*t67 + d15*spb35*(t40*t58 + t67) + d9*t8 + t69*t82 - d18*spb45*t89 + d19*spb45*t11*T(2) - d13*s235*t67*T(2) - d10*s23*t8*T(2) - d16*spb35*t8*T(2) + spb45*t31*t89*T(2) - s235*t67*t90*T(2) + d1*t67*T(4); 
complex<T> t7 = -(d12*s235*t67) - d15*spb35*t67 + t67*t69 - d4*s235*t76 - d3*spb35*t76 + t69*t76 + d15*spb35*t82 + d3*spb35*t82 + d13*s235*t67*T(2) + d6*s235*t76*T(2) + s235*t102*t76*T(2) - d9*t8*T(2) + d10*s23*t8*T(2) + d10*s25*t8*T(2) + d11*spb35*t8*T(2) + d16*spb35*t8*T(2) + s235*t67*t90*T(2) + d1*t82*T(8); 
complex<T> t18 = s23*(-(d9*t8) + t49*t57*t88 + t80*t88 + t103*t49*t9 + t123*t95); 
SeriesC<T> result = t5*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t7*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t6*(*CI_users[4]->get_value(mc,ind,mu)) + t18*(*CI_users[5]->get_value(mc,ind,mu)) + t4*(*CI_users[6]->get_value(mc,ind,mu)) + t1*(*CI_users[7]->get_value(mc,ind,mu)) + t1*(*CI_users[8]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdpmpp_G_wCI::\
C4g1ph_phdpmpp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c23, c45));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpmpp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:21:01 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> s14 = S(1,4);
complex<T> s235 = SS(2,3,5);
complex<T> s35 = -(spa35*spb35);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s345 = SS(3,4,5);
complex<T> s15 = S(1,5);
complex<T> s12 = S(1,2);
complex<T> t19 = square(spa23*spb25 - spa34*spb45); 
complex<T> t23 = s23*s34; 
complex<T> t25 = spa23*spa34; 
complex<T> t29 = square(square(spa23*spb25 - spa34*spb45)); 
complex<T> t30 = s234*T(2); 
complex<T> t45 = square(mH2); 
complex<T> t46 = square(spb25); 
complex<T> t47 = square(spb45); 
complex<T> t48 = square(spa23); 
complex<T> t49 = square(spa34); 
complex<T> t50 = cube(s245); 
complex<T> t51 = square(square(spb24)); 
complex<T> t53 = square(spa35); 
complex<T> t54 = spa35*spb23 + spa45*spb24; 
complex<T> t55 = spa25*spb24 + spa35*spb34; 
complex<T> t56 = -(spb35*T(2)); 
complex<T> t57 = spa24*spb34 + spa25*spb35; 
complex<T> t58 = spa24*spb25 + spa34*spb35; 
complex<T> t59 = spa24*spb23 + spa45*spb35; 
complex<T> t66 = s234*T(3); 
complex<T> t76 = cube(spb25); 
complex<T> t77 = cube(spb45); 
complex<T> t78 = square(spb35); 
complex<T> t92 = spb23*spb34; 
complex<T> t103 = coeff3mass<T>(mc,ind,3,2,5,4); 
complex<T> t104 = coeff3mass<T>(mc,ind,3,4,5,2); 
complex<T> t114 = spa23*spb35; 
complex<T> t115 = square(square(s245)); 
complex<T> t128 = spb25*spb45; 
complex<T> t168 = spa23*spb24; 
complex<T> t183 = spa34*spb24; 
complex<T> t188 = spb25*T(2); 
complex<T> d1 = spb24*cube(spa24); d1 = T(1)/d1;
complex<T> d2 = spa24*square(s24 + s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = spa25*square(s25 + s35)*T(3); d3 = T(1)/d3;
complex<T> d4 = square(s24 + s34)*square(spa24); d4 = T(1)/d4;
complex<T> d5 = (s24 + s34)*spb24*cube(spa24); d5 = T(1)/d5;
complex<T> d6 = spa24*cube(s24 + s34)*T(3); d6 = T(1)/d6;
complex<T> d7 = spa25*cube(s25 + s35)*T(3); d7 = T(1)/d7;
complex<T> d8 = s234*square(spa24); d8 = T(1)/d8;
complex<T> d9 = s234*spb24*cube(spa24); d9 = T(1)/d9;
complex<T> d11 = s235*(s25 + s35)*spa25*T(3); d11 = T(1)/d11;
complex<T> d12 = square(s23 + s24)*square(spa24); d12 = T(1)/d12;
complex<T> d13 = (s23 + s24)*spb24*cube(spa24); d13 = T(1)/d13;
complex<T> d14 = spa24*cube(s23 + s24)*T(3); d14 = T(1)/d14;
complex<T> d15 = spa24*square(s23 + s24)*T(3); d15 = T(1)/d15;
complex<T> d17 = spa45*square(s35 + s45)*T(3); d17 = T(1)/d17;
complex<T> d18 = spa45*cube(s35 + s45)*T(3); d18 = T(1)/d18;
complex<T> d19 = s345*(s35 + s45)*spa45*T(3); d19 = T(1)/d19;
complex<T> d29 = s234*cube(spa24); d29 = T(1)/d29;
complex<T> d30 = s234*square(square(spa24)); d30 = T(1)/d30;
complex<T> d36 = s234*spa24; d36 = T(1)/d36;
complex<T> t20 = square(spa35*spb25 + t183); 
complex<T> t21 = square(spa35*spb45 + t168); 
complex<T> t22 = d8*T(2); 
complex<T> t24 = d36*spb24; 
complex<T> t26 = cube(spa35*spb25 + t183); 
complex<T> t27 = cube(spa35*spb45 + t168); 
complex<T> t33 = s345*t54; 
complex<T> t34 = s235*t55; 
complex<T> t35 = t57*t59; 
complex<T> t60 = spa24*spb45 + t114; 
complex<T> t67 = d7*s235; 
complex<T> t72 = d18*s345; 
complex<T> t79 = -(t48*T(2)); 
complex<T> t89 = t54*t55; 
complex<T> t93 = t58*t59; 
complex<T> t96 = d29*s15; 
complex<T> t110 = -t128; 
complex<T> t112 = spb24*t49; 
complex<T> t130 = s45*t45; 
complex<T> t131 = spa34*t48; 
complex<T> t135 = d30*t78; 
complex<T> t139 = d3*spa35; 
complex<T> t143 = t128*t25; 
complex<T> t147 = d8*t19; 
complex<T> t154 = d29*t78; 
complex<T> t155 = t46*t48; 
complex<T> t157 = d1*T(4); 
complex<T> t166 = t47*t49; 
complex<T> t173 = t45*t76; 
complex<T> t181 = t45*t77; 
complex<T> t189 = spb45*t114; 
complex<T> t199 = t48*T(2); 
complex<T> t230 = s234*t46; 
complex<T> d10 = (s24 + s34)*spa24*t66; d10 = T(1)/d10;
complex<T> d16 = (s23 + s24)*spa24*t66; d16 = T(1)/d16;
complex<T> d52 = spa34*spa45*t54*T(2); d52 = T(1)/d52;
complex<T> d55 = spa23*spa25*t55*T(2); d55 = T(1)/d55;
complex<T> d59 = t30*cube(spa24); d59 = T(1)/d59;
complex<T> d60 = t30*square(square(spa24)); d60 = T(1)/d60;
complex<T> d61 = t30*square(spa24); d61 = T(1)/d61;
complex<T> t36 = t58*t60; 
complex<T> t94 = t57*t60; 
complex<T> t140 = d19*t20; 
complex<T> t159 = t135*t49; 
complex<T> t175 = spb35*t131; 
complex<T> t187 = t53*t67; 
complex<T> t198 = d6*t112; 
complex<T> t209 = t19*t24; 
complex<T> t213 = t189*t49; 
complex<T> d20 = spa34*spa45*t33; d20 = T(1)/d20;
complex<T> d21 = s234*t89*t92; d21 = T(1)/d21;
complex<T> d22 = spa25*spa45*t35; d22 = T(1)/d22;
complex<T> d23 = s235*spb23*t93; d23 = T(1)/d23;
complex<T> d26 = spa23*spa25*t34; d26 = T(1)/d26;
complex<T> d27 = spa34*spa45*t33*T(2); d27 = T(1)/d27;
complex<T> d28 = spb34*t30*t89; d28 = T(1)/d28;
complex<T> d31 = spa25*spa45*t35*T(2); d31 = T(1)/d31;
complex<T> d32 = s235*t93*T(2); d32 = T(1)/d32;
complex<T> d35 = spa25*t34*T(2); d35 = T(1)/d35;
complex<T> d37 = spa25*t35; d37 = T(1)/d37;
complex<T> d38 = t30*t89*t92; d38 = T(1)/d38;
complex<T> d39 = spa45*t35*T(2); d39 = T(1)/d39;
complex<T> d40 = s235*spb23*t93*T(2); d40 = T(1)/d40;
complex<T> d42 = spa23*t34*T(2); d42 = T(1)/d42;
complex<T> d43 = spa23*t34; d43 = T(1)/d43;
complex<T> d44 = spa45*t33*T(2); d44 = T(1)/d44;
complex<T> d45 = spb23*t30*t89; d45 = T(1)/d45;
complex<T> d48 = spa23*spa25*t34*T(2); d48 = T(1)/d48;
complex<T> d49 = spa34*t33; d49 = T(1)/d49;
complex<T> d50 = spa34*t33*T(2); d50 = T(1)/d50;
complex<T> d51 = spa25*t35*T(2); d51 = T(1)/d51;
complex<T> d54 = spb23*t93*T(2); d54 = T(1)/d54;
complex<T> d56 = t89*t92*T(2); d56 = T(1)/d56;
complex<T> d58 = s234*t89*T(4); d58 = T(1)/d58;
complex<T> d63 = t35*T(4); d63 = T(1)/d63;
complex<T> d64 = t34*T(2); d64 = T(1)/d64;
complex<T> d65 = t33*T(2); d65 = T(1)/d65;
complex<T> t2 = d17*spa35*t110*t183 + spb45*(-(t46*t53*(d17 - t72*T(2))) + t140*T(11)); 
complex<T> t5 = d12*s234*t155 + d15*spb24*t155 + t155*t157 + d4*s234*t166 + t157*t166 + t19*t22 + d15*spb24*t110*t25 + d2*spb24*t110*t25 + d2*t112*t47 + d13*t230*t79 + d14*spb24*t230*t79 - d5*s234*t166*T(2) - d9*s23*t19*T(2) - d9*s34*t19*T(2) - s234*t198*t47*T(2) - d1*t143*T(8) - d10*spb24*t19*T(11) - d16*spb24*t19*T(11); 
complex<T> t6 = -t147 + d15*spb24*(t143 - t155) - d12*s234*t155 + t143*t157 + d17*spa35*t128*t183 + d13*t155*t30 + d14*spb24*t155*t30 + d9*s23*t19*T(2) - d1*t155*T(4) + d16*spb24*t19*T(11) + spb45*(t46*t53*(d17 - t72*T(2)) - t140*T(11)); 
complex<T> t8 = d2*spb24*t143 - t147 + t143*t157 - d4*s234*t166 + t128*t139*t168 + d5*t166*t30 - d2*t112*t47 + t198*t30*t47 + d3*spb25*t47*t53 + d9*s34*t19*T(2) - spb25*t187*t47*T(2) - d1*t166*T(4) + d10*spb24*t19*T(11) - d11*spb25*t21*T(11); 
complex<T> t9 = t110*t139*t168 + t187*t188*t47 - d3*spb25*t47*t53 + d11*spb25*t21*T(11); 
complex<T> t14 = s34*(d54*t173 + d55*t27); 
complex<T> t18 = d23*(s14 - s25)*t173 + (d26*s14 + d43*spb25)*t27; 
complex<T> t39 = (d22*s13 + d37*spb45)*t50; 
complex<T> t43 = d55*s45*t27 + d54*t130*t76; 
complex<T> t44 = -(d32*s25*spa23*t173) + d64*spb23*spb25*t27; 
complex<T> t70 = d31*s34; 
complex<T> t97 = d31*s23; 
complex<T> d24 = s345*spb34*t94; d24 = T(1)/d24;
complex<T> d25 = s234*t25*t36; d25 = T(1)/d25;
complex<T> d33 = s345*spb34*t94*T(2); d33 = T(1)/d33;
complex<T> d34 = spa34*t30*t36; d34 = T(1)/d34;
complex<T> d41 = t25*t30*t36; d41 = T(1)/d41;
complex<T> d46 = s345*t94*T(2); d46 = T(1)/d46;
complex<T> d47 = spa23*t30*t36; d47 = T(1)/d47;
complex<T> d53 = spb34*t94*T(2); d53 = T(1)/d53;
complex<T> d57 = t25*t36*T(2); d57 = T(1)/d57;
complex<T> d62 = s234*t36*T(4); d62 = T(1)/d62;
complex<T> t1 = d61*t143*t23 + t147*t23 + d59*spb25*t175*t23 + d59*t213*t23 + d58*t25*t45*t51 + d60*t23*t48*t49*t78 + d62*t29*t92; 
complex<T> t3 = s23*(d53*t181 + d52*t26); 
complex<T> t4 = d24*s12*t181 + d20*s12*t26 + d49*spb45*t26 - d24*t130*t77; 
complex<T> t7 = s25*(d57*t29 + d56*t45*t51); 
complex<T> t10 = d23*s15*t173 + d24*s15*t181 - d20*mH2*t26 + d20*s15*t26 - d26*mH2*t27 + d26*s15*t27 - d25*mH2*t29 + d25*s15*t29 - d22*mH2*t50 + d22*s15*t50 + d21*s15*t45*t51 - d21*t51*cube(mH2) - d23*t76*cube(mH2) - d24*t77*cube(mH2); 
complex<T> t11 = -(d23*s24*t173) - d24*s24*t181 + t112*t154*t199 + t112*t189*t22 - d20*s24*t26 - d26*s24*t27 + spb25*spb35*t183*t22*t48 - d22*s24*t50 + t143*t24*T(2) + t209*T(4); 
complex<T> t12 = d23*s24*t173 + d24*s24*t181 + d20*s24*t26 + d26*s24*t27 + d22*s24*t50 + t112*t154*t79 + spb25*t131*t56*t96 + spa23*spb45*t49*t56*t96 - t143*t24*T(2) + d8*(spb25*t183*t48*t56 + spb45*t168*t49*t56 - s15*t143*T(2)) - t209*T(4) - s15*(d23*t173 + d24*t181 + d20*t26 + d26*t27 + d22*t50 - t159*t79 + t147*T(4)); 
complex<T> t13 = -(d46*spa34*t181) + d44*spb34*t26 - d47*spb34*t29 + d45*spa34*t45*t51 + t50*t70 + s34*(d40*t173 + d29*t175*t188 + t159*t199 + t143*t22 - d48*t27 + d29*t213*T(2) + t147*T(4)); 
complex<T> t15 = -(d50*spb45*t26) - d48*s45*t27 + d41*s45*t29 - d38*t130*t51 + d40*t130*t76 - d51*spb45*t50*T(3) + d33*t130*t77*T(3); 
complex<T> t16 = -(d32*spa23*t173) + d35*spb23*t27 - d34*spb23*t29 + d28*spa23*t45*t51 + t50*t97 + s23*(d33*t181 + d29*t175*t188 + t159*t199 + t143*t22 - d27*t26 + d29*t213*T(2) + t147*T(4)); 
complex<T> t17 = d33*s25*t181 - d27*s25*t26 - d42*spb25*t27 + d41*s25*t29 - d39*spb25*t50 - d38*s25*t45*t51 + d40*s25*t173*T(3); 
complex<T> t41 = s25*(d53*t181 + d52*t26); 
complex<T> t42 = d65*spb34*spb45*t26 - d46*spa34*t130*t77; 
complex<T> t75 = d57*s45*t29 + d56*t130*t51; 
complex<T> co1 = -t103; 
complex<T> co2 = -t104; 
complex<T> co3 = t115*t97; 
complex<T> co4 = t115*t70; 
complex<T> co5 = d63*t128*t50; 
SeriesC<T> result = t8*(*CI_users[0]->get_value(mc,ind,mu)) + t5*(*CI_users[1]->get_value(mc,ind,mu)) + t9*(*CI_users[2]->get_value(mc,ind,mu)) + t6*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t10*(*CI_users[7]->get_value(mc,ind,mu)) + t16*(*CI_users[8]->get_value(mc,ind,mu)) + t11*(*CI_users[9]->get_value(mc,ind,mu)) + t39*(*CI_users[10]->get_value(mc,ind,mu)) + t17*(*CI_users[11]->get_value(mc,ind,mu)) + t18*(*CI_users[12]->get_value(mc,ind,mu)) + t12*(*CI_users[13]->get_value(mc,ind,mu)) + t13*(*CI_users[14]->get_value(mc,ind,mu)) + t4*(*CI_users[15]->get_value(mc,ind,mu)) + t15*(*CI_users[16]->get_value(mc,ind,mu)) + t3*(*CI_users[17]->get_value(mc,ind,mu)) + t41*(*CI_users[18]->get_value(mc,ind,mu)) + co3*(*CI_users[19]->get_value(mc,ind,mu)) + co4*(*CI_users[20]->get_value(mc,ind,mu)) + t14*(*CI_users[21]->get_value(mc,ind,mu)) + t43*(*CI_users[22]->get_value(mc,ind,mu)) + t7*(*CI_users[23]->get_value(mc,ind,mu)) + t75*(*CI_users[24]->get_value(mc,ind,mu)) + t1*(*CI_users[25]->get_value(mc,ind,mu)) + co5*(*CI_users[26]->get_value(mc,ind,mu)) + t1*(*CI_users[27]->get_value(mc,ind,mu)) + co5*(*CI_users[28]->get_value(mc,ind,mu)) + t44*(*CI_users[29]->get_value(mc,ind,mu)) + t42*(*CI_users[30]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdpmpp_nf_wCI::\
C4g1ph_phdpmpp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpmpp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpmpp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:22:01 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> s23 = S(2,3);
complex<T> s234 = SS(2,3,4);
complex<T> s34 = S(3,4);
complex<T> s15 = S(1,5);
complex<T> s24 = -(spa24*spb24);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = -(spa35*spb35);
complex<T> s235 = SS(2,3,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> t9 = square(spa23*spb25 - spa34*spb45); 
complex<T> t10 = -(spb25*T(2)); 
complex<T> t11 = -(s23*s34); 
complex<T> t12 = square(spa34*spb24 + spa35*spb25); 
complex<T> t13 = square(spa23*spb24 + spa35*spb45); 
complex<T> t15 = square(spa24); 
complex<T> t22 = square(spa23); 
complex<T> t23 = square(spa34); 
complex<T> t24 = square(spb25); 
complex<T> t25 = square(spb45); 
complex<T> t26 = cube(spa24); 
complex<T> t27 = square(spa35); 
complex<T> t29 = square(spb35); 
complex<T> t49 = spa24*T(9); 
complex<T> t60 = spa23*spb45; 
complex<T> t61 = spa34*spb25; 
complex<T> t104 = spb24*spb25; 
complex<T> d3 = spa25*square(s25 + s35)*T(9); d3 = T(1)/d3;
complex<T> d7 = spa25*cube(s25 + s35)*T(9); d7 = T(1)/d7;
complex<T> d11 = s235*(s25 + s35)*spa25*T(9); d11 = T(1)/d11;
complex<T> d17 = spa45*square(s35 + s45)*T(9); d17 = T(1)/d17;
complex<T> d18 = spa45*cube(s35 + s45)*T(9); d18 = T(1)/d18;
complex<T> d19 = s345*(s35 + s45)*spa45*T(9); d19 = T(1)/d19;
complex<T> d21 = s234*square(square(spa24))*T(3); d21 = T(1)/d21;
complex<T> d22 = s234*spa24*T(3); d22 = T(1)/d22;
complex<T> d24 = s234*square(square(spa24))*T(6); d24 = T(1)/d24;
complex<T> t32 = d7*s235; 
complex<T> t35 = d18*s345; 
complex<T> t42 = -t60; 
complex<T> t64 = d17*spa35; 
complex<T> t66 = d21*t29; 
complex<T> t69 = t22*t24; 
complex<T> t76 = d3*spa35; 
complex<T> t78 = t23*t25; 
complex<T> t83 = t22*t29; 
complex<T> t85 = t60*t61; 
complex<T> t92 = spa34*t10; 
complex<T> t93 = spb45*t27; 
complex<T> t97 = -(t23*T(2)); 
complex<T> t102 = t23*t60; 
complex<T> d1 = spb24*t26*T(3); d1 = T(1)/d1;
complex<T> d2 = t49*square(s24 + s34); d2 = T(1)/d2;
complex<T> d4 = t15*square(s24 + s34)*T(3); d4 = T(1)/d4;
complex<T> d5 = (s24 + s34)*spb24*t26*T(3); d5 = T(1)/d5;
complex<T> d6 = t49*cube(s24 + s34); d6 = T(1)/d6;
complex<T> d8 = s234*t15*T(3); d8 = T(1)/d8;
complex<T> d9 = s234*spb24*t26*T(3); d9 = T(1)/d9;
complex<T> d10 = s234*(s24 + s34)*t49; d10 = T(1)/d10;
complex<T> d12 = t15*square(s23 + s24)*T(3); d12 = T(1)/d12;
complex<T> d13 = (s23 + s24)*spb24*t26*T(3); d13 = T(1)/d13;
complex<T> d14 = t49*cube(s23 + s24); d14 = T(1)/d14;
complex<T> d15 = t49*square(s23 + s24); d15 = T(1)/d15;
complex<T> d16 = s234*(s23 + s24)*t49; d16 = T(1)/d16;
complex<T> d20 = s234*t26*T(3); d20 = T(1)/d20;
complex<T> d23 = s234*t26*T(6); d23 = T(1)/d23;
complex<T> d25 = s234*t15*T(6); d25 = T(1)/d25;
complex<T> d26 = s234*t15*T(12); d26 = T(1)/d26;
complex<T> t7 = d11*t10*t13 + d3*spb25*t25*t27 + t10*t25*t27*t32 + t104*t60*t76; 
complex<T> t31 = d20*s15; 
complex<T> t41 = d23*t11; 
complex<T> t50 = -(d8*t9); 
complex<T> t52 = d20*spb35; 
complex<T> t53 = d9*s23; 
complex<T> t68 = d9*s34; 
complex<T> t71 = -(d1*T(4)); 
complex<T> t73 = d8*spb35; 
complex<T> t88 = d14*s234; 
complex<T> t94 = t22*t66; 
complex<T> t99 = d6*s234; 
complex<T> t111 = spb24*t64; 
complex<T> t112 = t23*t83; 
complex<T> t117 = t60*t92; 
complex<T> t1 = d24*t11*t112 + spb35*t41*(t102 + t22*t61) + t11*(d25*t85 + d26*t9); 
complex<T> t2 = spb45*t111*t61 - d19*spb45*t12*T(2) + t24*t93*(d17 - t35*T(2)); 
complex<T> t5 = spb24*(d22*(t117 - t9) + t22*t73*t92 + t60*t73*t97 + d20*t83*t97); 
complex<T> t20 = s23*(d8*t117 + t50 + t22*t52*t92 + t52*t60*t97 + t94*t97); 
complex<T> t59 = s34*(d8*t117 + t50 + t22*t52*t92 + t52*t60*t97 + t94*t97); 
complex<T> t90 = spb35*t31; 
complex<T> t91 = t53*t9; 
complex<T> t96 = t68*t9; 
complex<T> t3 = -(spb45*t111*t61) + d15*spb24*t42*t61 + d12*s234*t69 + d15*spb24*t69 + t71*t85 + d8*t9 - d17*t24*t93 + d19*spb45*t12*T(2) - d13*s234*t69*T(2) - spb24*t69*t88*T(2) - d16*spb24*t9*T(2) - t91*T(2) + t24*t35*t93*T(2) + d1*t69*T(4); 
complex<T> t4 = -(d12*s234*t69) - d15*spb24*t69 + t69*t71 - d4*s234*t78 - d2*spb24*t78 + t71*t78 + d15*spb24*t85 + d2*spb24*t85 + d13*s234*t69*T(2) + d5*s234*t78*T(2) + spb24*t69*t88*T(2) - d8*t9*T(2) + d10*spb24*t9*T(2) + d16*spb24*t9*T(2) + t91*T(2) + t96*T(2) + spb24*t78*t99*T(2) + d1*t85*T(8); 
complex<T> t6 = d8*s15*t9 + d20*spb24*t112*T(2) + d8*s15*t85*T(2) + t102*t90*T(2) + t22*t61*t90*T(2) + s15*t23*t94*T(2) + spb24*(d22*t9 + t102*t73*T(2) + t22*t61*t73*T(2) + d22*t85*T(2)); 
complex<T> t8 = -(d3*spb25*t25*t27) + d2*spb24*t42*t61 + t104*t42*t76 + d4*s234*t78 + d2*spb24*t78 + t71*t85 + d8*t9 + d11*spb25*t13*T(2) + spb25*t25*t27*t32*T(2) - d5*s234*t78*T(2) - d10*spb24*t9*T(2) - t96*T(2) - spb24*t78*t99*T(2) + d1*t78*T(4); 
SeriesC<T> result = t8*(*CI_users[0]->get_value(mc,ind,mu)) + t4*(*CI_users[1]->get_value(mc,ind,mu)) + t7*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t20*(*CI_users[5]->get_value(mc,ind,mu)) + t5*(*CI_users[6]->get_value(mc,ind,mu)) + t6*(*CI_users[7]->get_value(mc,ind,mu)) + t59*(*CI_users[8]->get_value(mc,ind,mu)) + t1*(*CI_users[9]->get_value(mc,ind,mu)) + t1*(*CI_users[10]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdppmp_G_wCI::\
C4g1ph_phdppmp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c34, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c45, c23));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdppmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdppmp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:25:30 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb23 = SPB(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> s235 = SS(2,3,5);
complex<T> s14 = S(1,4);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s345 = SS(3,4,5);
complex<T> s12 = S(1,2);
complex<T> s15 = S(1,5);
complex<T> s35 = -(spa35*spb35);
complex<T> t18 = square(spa34*spb23 - spa45*spb25); 
complex<T> t22 = square(square(spa34*spb23 - spa45*spb25)); 
complex<T> t24 = spa34*spa45; 
complex<T> t27 = s345*T(2); 
complex<T> t43 = square(mH2); 
complex<T> t44 = square(spb25); 
complex<T> t45 = square(spb23); 
complex<T> t46 = square(spa34); 
complex<T> t47 = square(spa45); 
complex<T> t48 = cube(s235); 
complex<T> t49 = square(square(spb35)); 
complex<T> t50 = square(spa24); 
complex<T> t52 = spa24*spb34 + spa25*spb35; 
complex<T> t53 = spa23*spb35 + spa24*spb45; 
complex<T> t55 = spa23*spb24 + spa35*spb45; 
complex<T> t57 = spa25*spb24 + spa35*spb34; 
complex<T> t58 = -(spb24*T(2)); 
complex<T> t59 = s34*s45; 
complex<T> t73 = cube(spb25); 
complex<T> t74 = cube(spb23); 
complex<T> t79 = square(spb24); 
complex<T> t82 = s24 + s25; 
complex<T> t83 = spa35*T(3); 
complex<T> t89 = spb34*spb45; 
 complex<T> t100 = coeff3mass<T>(mc,ind,4,3,2,5); 
 complex<T> t101 = coeff3mass<T>(mc,ind,4,5,2,3); 
complex<T> t113 = square(square(s235)); 
complex<T> t114 = spa34*spb24; 
complex<T> t126 = spb23*spb25; 
complex<T> t128 = spa45*spb24; 
complex<T> t158 = spa25*spa45; 
complex<T> t181 = spb23*T(2); 
complex<T> t183 = spa34*spb35; 
complex<T> t190 = spa45*spb35; 
complex<T> d1 = spa23*cube(s23 + s24)*T(3); d1 = T(1)/d1;
complex<T> d2 = spa23*square(s23 + s24)*T(3); d2 = T(1)/d2;
complex<T> d3 = s234*(s23 + s24)*spa23*T(3); d3 = T(1)/d3;
complex<T> d7 = square(s35 + s45)*square(spa35); d7 = T(1)/d7;
complex<T> d8 = s345*square(spa35); d8 = T(1)/d8;
complex<T> d9 = spb35*cube(spa35); d9 = T(1)/d9;
complex<T> d10 = (s35 + s45)*spb35*cube(spa35); d10 = T(1)/d10;
complex<T> d11 = s345*spb35*cube(spa35); d11 = T(1)/d11;
complex<T> d15 = square(s34 + s35)*square(spa35); d15 = T(1)/d15;
complex<T> d16 = (s34 + s35)*spb35*cube(spa35); d16 = T(1)/d16;
complex<T> d40 = s345*cube(spa35); d40 = T(1)/d40;
complex<T> d41 = s345*square(square(spa35)); d41 = T(1)/d41;
complex<T> t19 = square(spa24*spb23 + t190); 
complex<T> t20 = square(spa24*spb25 + t183); 
complex<T> t23 = cube(spa24*spb23 + t190); 
complex<T> t26 = cube(spa24*spb25 + t183); 
complex<T> t29 = t52*t53; 
complex<T> t32 = t55*t57; 
complex<T> t54 = spa35*spb25 + t114; 
complex<T> t56 = spa35*spb23 + t128; 
complex<T> t63 = d40*s12; 
complex<T> t65 = d1*s234; 
complex<T> t76 = -(t46*T(2)); 
complex<T> t87 = s234*t53; 
complex<T> t96 = d11*s45; 
complex<T> t109 = -t126; 
complex<T> t110 = t43*t49; 
complex<T> t111 = -(t47*T(2)); 
complex<T> t112 = d8*t18; 
complex<T> t144 = d9*T(4); 
complex<T> t146 = t47*t79; 
complex<T> t154 = t43*t73; 
complex<T> t166 = t126*t24; 
complex<T> t182 = spb25*t47; 
complex<T> t189 = spb25*t45; 
complex<T> t192 = d41*t46; 
complex<T> t196 = t128*t46; 
complex<T> d4 = spa25*cube(t82)*T(3); d4 = T(1)/d4;
complex<T> d5 = spa25*square(t82)*T(3); d5 = T(1)/d5;
complex<T> d6 = s245*spa25*t82*T(3); d6 = T(1)/d6;
complex<T> d12 = t83*square(s35 + s45); d12 = T(1)/d12;
complex<T> d13 = t83*cube(s35 + s45); d13 = T(1)/d13;
complex<T> d14 = s345*(s35 + s45)*t83; d14 = T(1)/d14;
complex<T> d17 = t83*cube(s34 + s35); d17 = T(1)/d17;
complex<T> d18 = t83*square(s34 + s35); d18 = T(1)/d18;
complex<T> d19 = s345*(s34 + s35)*t83; d19 = T(1)/d19;
complex<T> d22 = s245*t158*t52; d22 = T(1)/d22;
complex<T> d29 = s245*t158*t52*T(2); d29 = T(1)/d29;
complex<T> d34 = s245*spa25*t52; d34 = T(1)/d34;
complex<T> d36 = s245*spa45*t52*T(2); d36 = T(1)/d36;
complex<T> d50 = s245*spa25*t52*T(2); d50 = T(1)/d50;
complex<T> d55 = t158*t52*T(2); d55 = T(1)/d55;
complex<T> d58 = spa23*spa34*t53*T(2); d58 = T(1)/d58;
complex<T> d60 = s245*t52*T(2); d60 = T(1)/d60;
complex<T> d62 = t27*cube(spa35); d62 = T(1)/d62;
complex<T> d63 = t27*square(square(spa35)); d63 = T(1)/d63;
complex<T> d64 = t27*square(spa35); d64 = T(1)/d64;
complex<T> t28 = t54*t56; 
complex<T> t68 = d4*s245; 
complex<T> t75 = -t110; 
complex<T> t88 = t54*t55; 
complex<T> t91 = t56*t57; 
complex<T> t99 = d3*t20; 
complex<T> t139 = d5*spa24; 
complex<T> t186 = t50*t65; 
complex<T> t203 = t18*t96; 
complex<T> t205 = t114*t182; 
complex<T> d23 = spa23*spa34*t87; d23 = T(1)/d23;
complex<T> d24 = s345*t29*t89; d24 = T(1)/d24;
complex<T> d25 = spa23*spa25*t32; d25 = T(1)/d25;
complex<T> d30 = spa34*t87*T(2); d30 = T(1)/d30;
complex<T> d31 = t27*t29*t89; d31 = T(1)/d31;
complex<T> d32 = spa25*t32*T(2); d32 = T(1)/d32;
complex<T> d37 = spa23*spa34*t87*T(2); d37 = T(1)/d37;
complex<T> d38 = spa23*t32*T(2); d38 = T(1)/d38;
complex<T> d39 = spa23*t32; d39 = T(1)/d39;
complex<T> d44 = spa23*t87*T(2); d44 = T(1)/d44;
complex<T> d45 = spb45*t27*t29; d45 = T(1)/d45;
complex<T> d46 = spa23*spa25*t32*T(2); d46 = T(1)/d46;
complex<T> d48 = s345*spb34*t29; d48 = T(1)/d48;
complex<T> d51 = spb34*t27*t29; d51 = T(1)/d51;
complex<T> d54 = t29*t89*T(2); d54 = T(1)/d54;
complex<T> d59 = t87*T(2); d59 = T(1)/d59;
complex<T> d61 = t32*T(4); d61 = T(1)/d61;
complex<T> d66 = s345*t29*T(4); d66 = T(1)/d66;
complex<T> t2 = d12*spb35*t109*t24 + d18*spb35*t109*t24 + d10*s345*t111*t44 + d13*s345*spb35*t111*t44 + d15*s345*t45*t46 + d18*spb35*t45*t46 + t144*t45*t46 + d7*s345*t44*t47 + d12*spb35*t44*t47 + t144*t44*t47 + d16*s345*t45*t76 + d17*s345*spb35*t45*t76 + t112*T(2) - d11*s34*t18*T(2) - t203*T(2) - d9*t166*T(8) - d14*spb35*t18*T(11) - d19*spb35*t18*T(11); 
complex<T> t4 = -t112 + t144*t166 + d2*spa24*t126*t183 - d7*s345*t44*t47 + d10*t27*t44*t47 + d13*spb35*t27*t44*t47 + d12*spb35*(t166 - t44*t47) + d2*spb23*t44*t50 + t203*T(2) - spb23*t186*t44*T(2) - d9*t44*t47*T(4) + d14*spb35*t18*T(11) - spb23*t99*T(11); 
complex<T> t5 = d2*spa24*t109*t183 + t181*t186*t44 - d2*spb23*t44*t50 + spb23*t99*T(11); 
complex<T> t35 = (d25*s14 + d39*spb25)*t48; 
complex<T> t121 = d46*s34; 
complex<T> t138 = d46*s45; 
complex<T> t193 = t50*t68; 
complex<T> d20 = s345*t24*t28; d20 = T(1)/d20;
complex<T> d21 = s234*spb34*t91; d21 = T(1)/d21;
complex<T> d26 = s245*spb45*t88; d26 = T(1)/d26;
complex<T> d27 = t24*t27*t28; d27 = T(1)/d27;
complex<T> d28 = s234*spb34*t91*T(2); d28 = T(1)/d28;
complex<T> d33 = s245*spb45*t88*T(2); d33 = T(1)/d33;
complex<T> d35 = s245*t88; d35 = T(1)/d35;
complex<T> d42 = spa45*t27*t28; d42 = T(1)/d42;
complex<T> d43 = s234*t91*T(2); d43 = T(1)/d43;
complex<T> d47 = s345*spa34*t28; d47 = T(1)/d47;
complex<T> d49 = spa34*t27*t28; d49 = T(1)/d49;
complex<T> d52 = s245*t88*T(2); d52 = T(1)/d52;
complex<T> d53 = t24*t28*T(2); d53 = T(1)/d53;
complex<T> d56 = spb45*t88*T(2); d56 = T(1)/d56;
complex<T> d57 = spb34*t91*T(2); d57 = T(1)/d57;
complex<T> d65 = s345*t28*T(4); d65 = T(1)/d65;
complex<T> t1 = d66*t110*t24 + t112*t59 + d64*t166*t59 + d62*spb23*t196*t59 + d62*t205*t59 + d63*t146*t46*t59 + d65*t22*t89; 
complex<T> t3 = s23*(d54*t110 + d53*t22); 
complex<T> t6 = s25*(d58*t26 + d57*t43*t74); 
complex<T> t9 = d24*s15*t110 + d26*s15*t154 - d20*mH2*t22 + d20*s15*t22 - d22*mH2*t23 + d22*s15*t23 - d23*mH2*t26 + d23*s15*t26 - d25*mH2*t48 + d25*s15*t48 + d21*s15*t43*t74 - d24*t49*cube(mH2) - d26*t73*cube(mH2) - d21*t74*cube(mH2); 
complex<T> t10 = d33*s23*t154 + d27*s23*t22 - d29*s23*t23 + d30*spb23*t26 - d32*spb23*t48 + d28*s23*t43*t74 + d31*s23*t75; 
complex<T> t11 = -(s24*(d26*t154 + d20*t22 + d22*t23 + d25*t48 - d24*t75)); 
complex<T> t12 = d45*spa34*t110 - d42*spb34*t22 + d44*spb34*t26 + t121*t48 - d43*spa34*t43*t74 + s34*(d33*t154 + d40*t181*t196 - d29*t23 + d8*t166*T(2) + t146*t192*T(2) + d40*t205*T(2) + t112*T(4)); 
complex<T> t13 = s34*(d56*t154 + d55*t23); 
complex<T> t14 = d33*s25*t154 + d27*s25*t22 + d36*spb25*t23 - d37*s25*t26 + d28*s25*t43*t74 + d31*s25*t75 - d38*spb25*t48*T(3); 
complex<T> t15 = d24*s12*t110 + d48*spa45*t110 + d40*s45*t181*t196 + d20*s12*t22 + d47*spb45*t22 + spa34*t182*t58*t63 + spa45*spb23*t46*t58*t63 + d41*s12*t146*t76 - d8*s12*t166*T(2) + d8*s45*t166*T(2) + s45*t146*t192*T(2) + d40*s45*t205*T(2) - s12*t112*T(4) + s45*t112*T(4); 
complex<T> t16 = d26*s13*t154 + d35*spa45*t154 + d22*s13*t23 + d34*spb45*t23; 
complex<T> t17 = -(d50*spb45*t23) - d37*s45*t26 + t138*t48 + d28*s45*t43*t74 + d51*spa45*t75 - d52*spa45*t154*T(3) - d49*spb45*t22*T(3); 
complex<T> t37 = s25*(d54*t110 + d53*t22); 
complex<T> t38 = s23*(d56*t154 + d55*t23); 
complex<T> t39 = -((s15 - s24)*(d26*t154 + d20*t22 + d22*t23 + d25*t48)) + d24*(s24*t110 + s15*t75); 
complex<T> t40 = -(d52*s25*spa45*t154) + d60*spb25*spb45*t23; 
complex<T> t41 = s45*(d58*t26 + d57*t43*t74); 
complex<T> t42 = d59*spb23*spb34*t26 - d43*s23*spa34*t43*t74; 
complex<T> t220 = t189*t193; 
complex<T> t7 = -t112 + t144*t166 + t126*t139*t190 - d15*s345*t45*t46 + d16*t27*t45*t46 + d17*spb35*t27*t45*t46 + d18*spb35*(t166 - t45*t46) + d5*t189*t50 + d11*s34*t18*T(2) - t220*T(2) - d9*t45*t46*T(4) + d19*spb35*t18*T(11) - d6*spb25*t19*T(11); 
complex<T> t8 = t109*t139*t190 - d5*t189*t50 + t220*T(2) + d6*spb25*t19*T(11); 
complex<T> co1 = -t100; 
complex<T> co2 = -t101; 
complex<T> co3 = t113*t121; 
complex<T> co4 = t113*t138; 
complex<T> co5 = d61*t126*t48; 
SeriesC<T> result = t5*(*CI_users[0]->get_value(mc,ind,mu)) + t8*(*CI_users[1]->get_value(mc,ind,mu)) + t4*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t9*(*CI_users[7]->get_value(mc,ind,mu)) + t10*(*CI_users[8]->get_value(mc,ind,mu)) + t11*(*CI_users[9]->get_value(mc,ind,mu)) + t16*(*CI_users[10]->get_value(mc,ind,mu)) + t14*(*CI_users[11]->get_value(mc,ind,mu)) + t35*(*CI_users[12]->get_value(mc,ind,mu)) + t39*(*CI_users[13]->get_value(mc,ind,mu)) + t12*(*CI_users[14]->get_value(mc,ind,mu)) + t15*(*CI_users[15]->get_value(mc,ind,mu)) + t17*(*CI_users[16]->get_value(mc,ind,mu)) + t3*(*CI_users[17]->get_value(mc,ind,mu)) + t37*(*CI_users[18]->get_value(mc,ind,mu)) + t38*(*CI_users[19]->get_value(mc,ind,mu)) + t13*(*CI_users[20]->get_value(mc,ind,mu)) + co3*(*CI_users[21]->get_value(mc,ind,mu)) + co4*(*CI_users[22]->get_value(mc,ind,mu)) + t6*(*CI_users[23]->get_value(mc,ind,mu)) + t41*(*CI_users[24]->get_value(mc,ind,mu)) + t42*(*CI_users[25]->get_value(mc,ind,mu)) + t40*(*CI_users[26]->get_value(mc,ind,mu)) + co5*(*CI_users[27]->get_value(mc,ind,mu)) + t1*(*CI_users[28]->get_value(mc,ind,mu)) + co5*(*CI_users[29]->get_value(mc,ind,mu)) + t1*(*CI_users[30]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdppmp_nf_wCI::\
C4g1ph_phdppmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdppmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdppmp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:26:16 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb24 = SPB(2,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s234 = SS(2,3,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s245 = SS(2,4,5);
complex<T> s34 = S(3,4);
complex<T> s345 = SS(3,4,5);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> s35 = -(spa35*spb35);
complex<T> t8 = square(spa34*spb23 - spa45*spb25); 
complex<T> t9 = -(spb23*T(2)); 
complex<T> t10 = -(s34*s45); 
complex<T> t11 = square(spa24*spb25 + spa34*spb35); 
complex<T> t12 = square(spa24*spb23 + spa45*spb35); 
complex<T> t13 = square(spa35); 
complex<T> t20 = square(spa34); 
complex<T> t21 = square(spa45); 
complex<T> t22 = square(spb23); 
complex<T> t23 = square(spb25); 
complex<T> t24 = square(spa24); 
complex<T> t25 = cube(spa35); 
complex<T> t27 = square(spb24); 
complex<T> t47 = spa35*T(9); 
complex<T> t57 = spa34*spb25; 
complex<T> t58 = spa45*spb23; 
complex<T> t107 = spb23*spb35; 
complex<T> t109 = spb25*spb35; 
complex<T> d1 = spa23*cube(s23 + s24)*T(9); d1 = T(1)/d1;
complex<T> d2 = spa23*square(s23 + s24)*T(9); d2 = T(1)/d2;
complex<T> d3 = s234*(s23 + s24)*spa23*T(9); d3 = T(1)/d3;
complex<T> d4 = spa25*cube(s24 + s25)*T(9); d4 = T(1)/d4;
complex<T> d5 = spa25*square(s24 + s25)*T(9); d5 = T(1)/d5;
complex<T> d6 = s245*(s24 + s25)*spa25*T(9); d6 = T(1)/d6;
complex<T> d21 = s345*square(square(spa35))*T(3); d21 = T(1)/d21;
complex<T> d23 = s345*square(square(spa35))*T(6); d23 = T(1)/d23;
complex<T> t30 = d1*s234; 
complex<T> t32 = d4*s245; 
complex<T> t40 = -t57; 
complex<T> t62 = d2*spa24; 
complex<T> t63 = d21*t27; 
complex<T> t74 = d5*spa24; 
complex<T> t76 = t20*t22; 
complex<T> t77 = t21*t23; 
complex<T> t83 = t57*t58; 
complex<T> t89 = -(t21*T(2)); 
complex<T> t104 = spa45*t20; 
complex<T> t105 = t23*t24; 
complex<T> t126 = t57*t9; 
complex<T> d7 = t13*square(s35 + s45)*T(3); d7 = T(1)/d7;
complex<T> d8 = s345*t13*T(3); d8 = T(1)/d8;
complex<T> d9 = spb35*t25*T(3); d9 = T(1)/d9;
complex<T> d10 = (s35 + s45)*spb35*t25*T(3); d10 = T(1)/d10;
complex<T> d11 = s345*spb35*t25*T(3); d11 = T(1)/d11;
complex<T> d12 = t47*square(s35 + s45); d12 = T(1)/d12;
complex<T> d13 = t47*cube(s35 + s45); d13 = T(1)/d13;
complex<T> d14 = s345*(s35 + s45)*t47; d14 = T(1)/d14;
complex<T> d15 = t13*square(s34 + s35)*T(3); d15 = T(1)/d15;
complex<T> d16 = (s34 + s35)*spb35*t25*T(3); d16 = T(1)/d16;
complex<T> d17 = t47*cube(s34 + s35); d17 = T(1)/d17;
complex<T> d18 = t47*square(s34 + s35); d18 = T(1)/d18;
complex<T> d19 = s345*(s34 + s35)*t47; d19 = T(1)/d19;
complex<T> d20 = s345*t25*T(3); d20 = T(1)/d20;
complex<T> d22 = s345*t25*T(6); d22 = T(1)/d22;
complex<T> d24 = s345*t13*T(6); d24 = T(1)/d24;
complex<T> d25 = s345*t13*T(12); d25 = T(1)/d25;
complex<T> t4 = d2*spb23*t105 + t107*t57*t62 + d3*t11*t9 + t105*t30*t9; 
complex<T> t29 = d20*s12; 
complex<T> t39 = d22*t10; 
complex<T> t49 = d20*spb24; 
complex<T> t69 = -(d9*T(4)); 
complex<T> t85 = d13*spb35; 
complex<T> t91 = t20*t63; 
complex<T> t97 = d8*spa45; 
complex<T> t98 = d17*spb35; 
complex<T> t99 = t22*t32; 
complex<T> t1 = d23*t10*t20*t21*t27 + spb24*t21*t39*t57 + spb24*t20*t39*t58 + d25*t10*t8 + d24*t10*t83; 
complex<T> t2 = s45*t49*t57*t89 + s45*t104*t49*t9 + s45*t89*t91 + s45*t126*t97 + spb24*t29*(t21*t57 + t20*t58)*T(2) + s12*t21*t91*T(2) + d8*(s12*t8 - s45*t8 + s12*t83*T(2)); 
complex<T> t3 = -(d15*s345*t76) - d18*spb35*t76 + t69*t76 - d7*s345*t77 - d12*spb35*t77 + t69*t77 + d12*spb35*t83 + d18*spb35*t83 + d16*s345*t76*T(2) + d10*s345*t77*T(2) - d8*t8*T(2) + d11*s34*t8*T(2) + d11*s45*t8*T(2) + d14*spb35*t8*T(2) + d19*spb35*t8*T(2) + s345*t77*t85*T(2) + s345*t76*t98*T(2) + d9*t83*T(8); 
complex<T> t5 = -(d2*spb23*t105) + d12*spb35*t40*t58 + t107*t40*t62 + d7*s345*t77 + d12*spb35*t77 + d8*t8 + t69*t83 + d3*spb23*t11*T(2) + spb23*t105*t30*T(2) - d10*s345*t77*T(2) - d11*s45*t8*T(2) - d14*spb35*t8*T(2) - s345*t77*t85*T(2) + d9*t77*T(4); 
complex<T> t6 = d5*spb25*t22*t24 + t109*t58*t74 - d6*spb25*t12*T(2) - spb25*t24*t99*T(2); 
complex<T> t7 = -(d5*spb25*t22*t24) + d18*spb35*t40*t58 - t109*t58*t74 + d15*s345*t76 + d18*spb35*t76 + d8*t8 + t69*t83 + d6*spb25*t12*T(2) - d16*s345*t76*T(2) - d11*s34*t8*T(2) - d19*spb35*t8*T(2) - s345*t76*t98*T(2) + spb25*t24*t99*T(2) + d9*t76*T(4); 
complex<T> t18 = s34*(-(d8*t8) + t49*t57*t89 + t104*t49*t9 + t89*t91 + t126*t97); 
SeriesC<T> result = t4*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + t18*(*CI_users[5]->get_value(mc,ind,mu)) + t2*(*CI_users[6]->get_value(mc,ind,mu)) + t1*(*CI_users[7]->get_value(mc,ind,mu)) + t1*(*CI_users[8]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdpppm_G_wCI::\
C4g1ph_phdpppm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c25, c34));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c45, c23));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpppm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpppm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:29:50 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s14 = S(1,4);
complex<T> s35 = -(spa35*spb35);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s24 = -(spa24*spb24);
complex<T> s345 = SS(3,4,5);
complex<T> s12 = S(1,2);
complex<T> s15 = S(1,5);
complex<T> t18 = square(spa25*spb23 - spa45*spb34); 
complex<T> t24 = square(square(spa25*spb23 - spa45*spb34)); 
complex<T> t25 = spa25*spa45; 
complex<T> t27 = s245*T(2); 
complex<T> t42 = square(mH2); 
complex<T> t43 = square(spb23); 
complex<T> t44 = square(spb34); 
complex<T> t45 = square(spa25); 
complex<T> t46 = square(spa45); 
complex<T> t47 = square(square(spb24)); 
complex<T> t48 = cube(s234); 
complex<T> t50 = square(spa35); 
complex<T> t51 = spa34*spb24 + spa35*spb25; 
complex<T> t52 = spa23*spb24 + spa35*spb45; 
complex<T> t54 = spa24*spb25 + spa34*spb35; 
complex<T> t56 = spa23*spb35 + spa24*spb45; 
complex<T> t57 = -(spb35*T(2)); 
complex<T> t58 = s25*s45; 
complex<T> t66 = s245*T(3); 
complex<T> t73 = cube(spb23); 
complex<T> t74 = cube(spb34); 
complex<T> t79 = square(spb35); 
complex<T> t89 = square(square(s234)); 
complex<T> t91 = spb25*spb45; 
complex<T> t100 = coeff3mass<T>(mc,ind,5,2,3,4); 
complex<T> t101 = coeff3mass<T>(mc,ind,5,4,3,2); 
complex<T> t111 = spa25*spb35; 
complex<T> t124 = spb23*spb34; 
complex<T> t127 = spa45*spb35; 
complex<T> t130 = spa23*spa25; 
complex<T> t156 = spa34*spa45; 
complex<T> t187 = spb34*T(2); 
complex<T> t188 = spa25*spb24; 
complex<T> t196 = spa45*spb24; 
complex<T> d1 = spa23*square(s23 + s35)*T(3); d1 = T(1)/d1;
complex<T> d2 = spa23*cube(s23 + s35)*T(3); d2 = T(1)/d2;
complex<T> d3 = s235*(s23 + s35)*spa23*T(3); d3 = T(1)/d3;
complex<T> d4 = square(s24 + s25)*square(spa24); d4 = T(1)/d4;
complex<T> d5 = spb24*cube(spa24); d5 = T(1)/d5;
complex<T> d6 = (s24 + s25)*spb24*cube(spa24); d6 = T(1)/d6;
complex<T> d7 = spa24*cube(s24 + s25)*T(3); d7 = T(1)/d7;
complex<T> d8 = spa24*square(s24 + s25)*T(3); d8 = T(1)/d8;
complex<T> d9 = spa24*square(s24 + s45)*T(3); d9 = T(1)/d9;
complex<T> d10 = square(s24 + s45)*square(spa24); d10 = T(1)/d10;
complex<T> d11 = (s24 + s45)*spb24*cube(spa24); d11 = T(1)/d11;
complex<T> d12 = spa24*cube(s24 + s45)*T(3); d12 = T(1)/d12;
complex<T> d13 = s245*square(spa24); d13 = T(1)/d13;
complex<T> d14 = s245*spb24*cube(spa24); d14 = T(1)/d14;
complex<T> d17 = spa34*cube(s34 + s35)*T(3); d17 = T(1)/d17;
complex<T> d18 = spa34*square(s34 + s35)*T(3); d18 = T(1)/d18;
complex<T> d19 = s345*(s34 + s35)*spa34*T(3); d19 = T(1)/d19;
complex<T> d34 = s245*cube(spa24); d34 = T(1)/d34;
complex<T> d35 = s245*square(square(spa24)); d35 = T(1)/d35;
complex<T> t19 = square(spa35*spb23 + t196); 
complex<T> t20 = square(spa35*spb34 + t188); 
complex<T> t22 = cube(spa35*spb23 + t196); 
complex<T> t23 = cube(spa35*spb34 + t188); 
complex<T> t29 = t51*t52; 
complex<T> t32 = t54*t56; 
complex<T> t53 = spa24*spb34 + t111; 
complex<T> t55 = spa24*spb23 + t127; 
complex<T> t62 = d34*s13; 
complex<T> t64 = d2*s235; 
complex<T> t68 = d17*s345; 
complex<T> t76 = -(t45*T(2)); 
complex<T> t97 = d14*s45; 
complex<T> t107 = -t124; 
complex<T> t108 = t42*t47; 
complex<T> t109 = -(t46*T(2)); 
complex<T> t110 = d13*t18; 
complex<T> t137 = d18*spa35; 
complex<T> t141 = d5*T(4); 
complex<T> t143 = t46*t79; 
complex<T> t152 = t42*t74; 
complex<T> t155 = spb23*t45; 
complex<T> t165 = t42*t73; 
complex<T> t175 = t124*t25; 
complex<T> t202 = t111*t46; 
complex<T> d15 = (s24 + s25)*spa24*t66; d15 = T(1)/d15;
complex<T> d16 = (s24 + s45)*spa24*t66; d16 = T(1)/d16;
complex<T> d20 = s345*t156*t51; d20 = T(1)/d20;
complex<T> d25 = s235*t130*t52; d25 = T(1)/d25;
complex<T> d27 = s345*t156*t51*T(2); d27 = T(1)/d27;
complex<T> d32 = s235*spa25*t52*T(2); d32 = T(1)/d32;
complex<T> d41 = s235*spa23*t52*T(2); d41 = T(1)/d41;
complex<T> d44 = s235*spa23*t52; d44 = T(1)/d44;
complex<T> d45 = s345*spa45*t51*T(2); d45 = T(1)/d45;
complex<T> d47 = s235*t130*t52*T(2); d47 = T(1)/d47;
complex<T> d48 = s345*spa34*t51; d48 = T(1)/d48;
complex<T> d50 = s345*spa34*t51*T(2); d50 = T(1)/d50;
complex<T> d54 = t156*t51*T(2); d54 = T(1)/d54;
complex<T> d59 = t130*t52*T(2); d59 = T(1)/d59;
complex<T> d61 = t27*square(spa24); d61 = T(1)/d61;
complex<T> d62 = t27*cube(spa24); d62 = T(1)/d62;
complex<T> d63 = t27*square(square(spa24)); d63 = T(1)/d63;
complex<T> d66 = s235*t52*T(2); d66 = T(1)/d66;
complex<T> d67 = s345*t51*T(2); d67 = T(1)/d67;
complex<T> t28 = t53*t55; 
complex<T> t75 = -t108; 
complex<T> t87 = t54*t55; 
complex<T> t88 = t53*t56; 
complex<T> t161 = d3*t20; 
complex<T> t182 = d35*t143; 
complex<T> t186 = t50*t64; 
complex<T> t193 = t50*t68; 
complex<T> t200 = t127*t155; 
complex<T> t209 = t18*t97; 
complex<T> d23 = spa23*spa34*t32; d23 = T(1)/d23;
complex<T> d26 = s245*t29*t91; d26 = T(1)/d26;
complex<T> d30 = spa34*t32*T(2); d30 = T(1)/d30;
complex<T> d33 = t27*t29*t91; d33 = T(1)/d33;
complex<T> d37 = s245*spb25*t29; d37 = T(1)/d37;
complex<T> d40 = spa23*spa34*t32*T(2); d40 = T(1)/d40;
complex<T> d42 = spb45*t27*t29; d42 = T(1)/d42;
complex<T> d46 = spa23*t32*T(2); d46 = T(1)/d46;
complex<T> d53 = spb25*t27*t29; d53 = T(1)/d53;
complex<T> d57 = t29*t91*T(2); d57 = T(1)/d57;
complex<T> d60 = t32*T(4); d60 = T(1)/d60;
complex<T> d65 = s245*t29*T(4); d65 = T(1)/d65;
complex<T> t2 = t107*t137*t196 + t187*t193*t43 - d18*spb34*t43*t50 + d19*spb34*t19*T(11); 
complex<T> t4 = d1*spa35*t107*t188 + spb23*(-(d1*t44*t50) + t186*t44*T(2) + t161*T(11)); 
complex<T> t7 = -t110 + t141*t175 + t124*t137*t196 - d4*s245*t43*t45 + d6*t27*t43*t45 + d7*spb24*t27*t43*t45 + d8*spb24*(t175 - t43*t45) + d18*spb34*t43*t50 + d14*s25*t18*T(2) - spb34*t193*t43*T(2) - d5*t43*t45*T(4) + d15*spb24*t18*T(11) - d19*spb34*t19*T(11); 
complex<T> t8 = d8*spb24*t107*t25 + d9*spb24*t107*t25 + d11*s245*t109*t44 + d12*s245*spb24*t109*t44 + d4*s245*t43*t45 + d8*spb24*t43*t45 + t141*t43*t45 + d10*s245*t44*t46 + d9*spb24*t44*t46 + t141*t44*t46 + d6*s245*t43*t76 + d7*s245*spb24*t43*t76 + t110*T(2) - d14*s25*t18*T(2) - t209*T(2) - d5*t175*T(8) - d15*spb24*t18*T(11) - d16*spb24*t18*T(11); 
complex<T> t9 = -t110 + d9*spb24*t175 + t141*t175 + d1*spa35*t124*t188 - d10*s245*t44*t46 - d9*spb24*t44*t46 + d11*t27*t44*t46 + d12*spb24*t27*t44*t46 + t209*T(2) - d5*t44*t46*T(4) + d16*spb24*t18*T(11) + spb23*(d1*t44*t50 - t186*t44*T(2) - t161*T(11)); 
complex<T> t120 = d40*s25; 
complex<T> t136 = d40*s45; 
complex<T> d21 = s245*t25*t28; d21 = T(1)/d21;
complex<T> d22 = s235*spb25*t87; d22 = T(1)/d22;
complex<T> d24 = s345*spb45*t88; d24 = T(1)/d24;
complex<T> d28 = t25*t27*t28; d28 = T(1)/d28;
complex<T> d29 = s235*spb25*t87*T(2); d29 = T(1)/d29;
complex<T> d31 = s345*spb45*t88*T(2); d31 = T(1)/d31;
complex<T> d36 = s245*spa25*t28; d36 = T(1)/d36;
complex<T> d38 = spa45*t27*t28; d38 = T(1)/d38;
complex<T> d39 = s235*t87*T(2); d39 = T(1)/d39;
complex<T> d43 = s235*t87; d43 = T(1)/d43;
complex<T> d49 = s345*t88; d49 = T(1)/d49;
complex<T> d51 = spa25*t27*t28; d51 = T(1)/d51;
complex<T> d52 = s345*t88*T(2); d52 = T(1)/d52;
complex<T> d55 = spb45*t88*T(2); d55 = T(1)/d55;
complex<T> d56 = t25*t28*T(2); d56 = T(1)/d56;
complex<T> d58 = spb25*t87*T(2); d58 = T(1)/d58;
complex<T> d64 = s245*t28*T(4); d64 = T(1)/d64;
complex<T> t1 = d65*t108*t25 + t110*t58 + d61*t175*t58 + d62*t200*t58 + d62*spb34*t202*t58 + d63*t143*t45*t58 + d64*t24*t91; 
complex<T> t3 = s23*(d55*t152 + d54*t22); 
complex<T> t5 = s34*(d58*t165 + d59*t23); 
complex<T> t6 = d22*s14*t165 + d43*spa25*t165 + d25*s14*t23 + d44*spb25*t23; 
complex<T> t10 = d26*s15*t108 + d24*s15*t152 + d22*s15*t165 - d20*mH2*t22 + d20*s15*t22 - d25*mH2*t23 + d25*s15*t23 - d21*mH2*t24 + d21*s15*t24 - d23*mH2*t48 + d23*s15*t48 - d26*t47*cube(mH2) - d22*t73*cube(mH2) - d24*t74*cube(mH2); 
complex<T> t11 = d31*s23*t152 + d29*s23*t165 - d27*s23*t22 + d32*spb23*t23 + d28*s23*t24 - d30*spb23*t48 + d33*s23*t75; 
complex<T> t12 = -(s24*(d24*t152 + d22*t165 + d20*t22 + d25*t23 + d21*t24 - d26*t75)); 
complex<T> t13 = d31*s34*t152 + d29*s34*t165 + d45*spb34*t22 - d47*s34*t23 + d28*s34*t24 - d46*spb34*t48 + d33*s34*t75; 
complex<T> t14 = d42*spa25*t108 - d41*spb25*t23 - d38*spb25*t24 + t120*t48 - d39*spa25*t165*T(3) + s25*(d31*t152 + d34*t187*t202 - d27*t22 + d13*t175*T(2) + d34*t200*T(2) + t182*t45*T(2) + t110*T(4)); 
complex<T> t15 = d24*s12*t152 + d49*spa45*t152 + d20*s12*t22 + d48*spb45*t22; 
complex<T> t16 = d26*s13*t108 + d37*spa45*t108 + d34*s45*t187*t202 + d21*s13*t24 + d36*spb45*t24 + spa45*t155*t57*t62 + spa25*spb34*t46*t57*t62 + s13*t182*t76 - d13*s13*t175*T(2) + d13*s45*t175*T(2) + d34*s45*t200*T(2) + s45*t182*t45*T(2) - s13*t110*T(4) + s45*t110*T(4); 
complex<T> t17 = d29*s45*t165 - d50*spb45*t22 - d47*s45*t23 + t136*t48 + d53*spa45*t75 - d52*spa45*t152*T(3) - d51*spb45*t24*T(3); 
complex<T> t36 = s25*(d55*t152 + d54*t22); 
complex<T> t37 = -(d52*s34*spa45*t152) + d67*spb34*spb45*t22; 
complex<T> t38 = s45*(d58*t165 + d59*t23); 
complex<T> t39 = -(d39*s23*spa25*t165) + d66*spb23*spb25*t23; 
complex<T> t40 = s23*(d57*t108 + d56*t24); 
complex<T> t41 = -((s15 - s24)*(d24*t152 + d22*t165 + d20*t22 + d25*t23 + d21*t24)) + d26*(s24*t108 + s15*t75); 
complex<T> t72 = s34*(d57*t108 + d56*t24); 
complex<T> co1 = -t100; 
complex<T> co2 = -t101; 
complex<T> co3 = t120*t89; 
complex<T> co4 = t136*t89; 
complex<T> co5 = d60*t124*t48; 
SeriesC<T> result = t4*(*CI_users[0]->get_value(mc,ind,mu)) + t8*(*CI_users[1]->get_value(mc,ind,mu)) + t9*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t10*(*CI_users[7]->get_value(mc,ind,mu)) + t11*(*CI_users[8]->get_value(mc,ind,mu)) + t12*(*CI_users[9]->get_value(mc,ind,mu)) + t16*(*CI_users[10]->get_value(mc,ind,mu)) + t14*(*CI_users[11]->get_value(mc,ind,mu)) + t6*(*CI_users[12]->get_value(mc,ind,mu)) + t41*(*CI_users[13]->get_value(mc,ind,mu)) + t13*(*CI_users[14]->get_value(mc,ind,mu)) + t15*(*CI_users[15]->get_value(mc,ind,mu)) + t17*(*CI_users[16]->get_value(mc,ind,mu)) + t3*(*CI_users[17]->get_value(mc,ind,mu)) + t36*(*CI_users[18]->get_value(mc,ind,mu)) + t40*(*CI_users[19]->get_value(mc,ind,mu)) + t72*(*CI_users[20]->get_value(mc,ind,mu)) + t5*(*CI_users[21]->get_value(mc,ind,mu)) + t38*(*CI_users[22]->get_value(mc,ind,mu)) + co3*(*CI_users[23]->get_value(mc,ind,mu)) + co4*(*CI_users[24]->get_value(mc,ind,mu)) + co5*(*CI_users[25]->get_value(mc,ind,mu)) + t1*(*CI_users[26]->get_value(mc,ind,mu)) + t39*(*CI_users[27]->get_value(mc,ind,mu)) + t37*(*CI_users[28]->get_value(mc,ind,mu)) + co5*(*CI_users[29]->get_value(mc,ind,mu)) + t1*(*CI_users[30]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdpppm_nf_wCI::\
C4g1ph_phdpppm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpppm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpppm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:30:35 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s35 = -(spa35*spb35);
complex<T> s235 = SS(2,3,5);
complex<T> s25 = S(2,5);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s45 = S(4,5);
complex<T> s24 = -(spa24*spb24);
complex<T> s34 = -(spa34*spb34);
complex<T> s345 = SS(3,4,5);
complex<T> t8 = square(spa25*spb23 - spa45*spb34); 
complex<T> t9 = -(spb23*T(2)); 
complex<T> t10 = -(s25*s45); 
complex<T> t11 = square(spa35*spb23 + spa45*spb24); 
complex<T> t12 = square(spa25*spb24 + spa35*spb34); 
complex<T> t13 = square(spa24); 
complex<T> t20 = square(spa25); 
complex<T> t21 = square(spa45); 
complex<T> t22 = square(spb23); 
complex<T> t23 = square(spb34); 
complex<T> t24 = cube(spa24); 
complex<T> t25 = square(spa35); 
complex<T> t27 = square(spb35); 
complex<T> t47 = spa24*T(9); 
complex<T> t57 = spa25*spb34; 
complex<T> t58 = spa45*spb23; 
complex<T> t107 = spb23*spb24; 
complex<T> t109 = spb24*spb34; 
complex<T> d1 = spa23*square(s23 + s35)*T(9); d1 = T(1)/d1;
complex<T> d2 = spa23*cube(s23 + s35)*T(9); d2 = T(1)/d2;
complex<T> d3 = s235*(s23 + s35)*spa23*T(9); d3 = T(1)/d3;
complex<T> d17 = spa34*cube(s34 + s35)*T(9); d17 = T(1)/d17;
complex<T> d18 = spa34*square(s34 + s35)*T(9); d18 = T(1)/d18;
complex<T> d19 = s345*(s34 + s35)*spa34*T(9); d19 = T(1)/d19;
complex<T> d21 = s245*square(square(spa24))*T(3); d21 = T(1)/d21;
complex<T> d25 = s245*square(square(spa24))*T(6); d25 = T(1)/d25;
complex<T> t30 = d2*s235; 
complex<T> t32 = d17*s345; 
complex<T> t40 = -t57; 
complex<T> t62 = d1*spa35; 
complex<T> t63 = d21*t27; 
complex<T> t74 = d18*spa35; 
complex<T> t76 = t20*t22; 
complex<T> t77 = t21*t23; 
complex<T> t83 = t57*t58; 
complex<T> t89 = -(t21*T(2)); 
complex<T> t104 = spa45*t20; 
complex<T> t105 = t23*t25; 
complex<T> t127 = t57*t9; 
complex<T> d4 = t13*square(s24 + s25)*T(3); d4 = T(1)/d4;
complex<T> d5 = spb24*t24*T(3); d5 = T(1)/d5;
complex<T> d6 = (s24 + s25)*spb24*t24*T(3); d6 = T(1)/d6;
complex<T> d7 = t47*cube(s24 + s25); d7 = T(1)/d7;
complex<T> d8 = t47*square(s24 + s25); d8 = T(1)/d8;
complex<T> d9 = t47*square(s24 + s45); d9 = T(1)/d9;
complex<T> d10 = t13*square(s24 + s45)*T(3); d10 = T(1)/d10;
complex<T> d11 = (s24 + s45)*spb24*t24*T(3); d11 = T(1)/d11;
complex<T> d12 = t47*cube(s24 + s45); d12 = T(1)/d12;
complex<T> d13 = s245*t13*T(3); d13 = T(1)/d13;
complex<T> d14 = s245*spb24*t24*T(3); d14 = T(1)/d14;
complex<T> d15 = s245*(s24 + s25)*t47; d15 = T(1)/d15;
complex<T> d16 = s245*(s24 + s45)*t47; d16 = T(1)/d16;
complex<T> d20 = s245*t24*T(3); d20 = T(1)/d20;
complex<T> d22 = s245*t13*T(6); d22 = T(1)/d22;
complex<T> d23 = s245*t13*T(12); d23 = T(1)/d23;
complex<T> d24 = s245*t24*T(6); d24 = T(1)/d24;
complex<T> t3 = d1*spb23*t105 + t107*t57*t62 + d3*t12*t9 + t105*t30*t9; 
complex<T> t29 = d20*s13; 
complex<T> t39 = d24*t10; 
complex<T> t50 = d20*spb35; 
complex<T> t69 = -(d5*T(4)); 
complex<T> t85 = d12*spb24; 
complex<T> t91 = t20*t63; 
complex<T> t97 = d13*spa45; 
complex<T> t99 = t22*t32; 
complex<T> t103 = d7*spb24; 
complex<T> t1 = d25*t10*t20*t21*t27 + spb35*t21*t39*t57 + spb35*t20*t39*t58 + d23*t10*t8 + d22*t10*t83; 
complex<T> t2 = d18*spb34*t22*t25 + t109*t58*t74 - d19*spb34*t11*T(2) - spb34*t25*t99*T(2); 
complex<T> t4 = -(d18*spb34*t22*t25) + d8*spb24*t40*t58 - t109*t58*t74 + d4*s245*t76 + d8*spb24*t76 + d13*t8 + t69*t83 + d19*spb34*t11*T(2) - d6*s245*t76*T(2) - s245*t103*t76*T(2) - d14*s25*t8*T(2) - d15*spb24*t8*T(2) + spb34*t25*t99*T(2) + d5*t76*T(4); 
complex<T> t5 = -(d1*spb23*t105) + d9*spb24*t40*t58 + t107*t40*t62 + d10*s245*t77 + d9*spb24*t77 + d13*t8 + t69*t83 + d3*spb23*t12*T(2) + spb23*t105*t30*T(2) - d11*s245*t77*T(2) - d14*s45*t8*T(2) - d16*spb24*t8*T(2) - s245*t77*t85*T(2) + d5*t77*T(4); 
complex<T> t6 = -(d4*s245*t76) - d8*spb24*t76 + t69*t76 - d10*s245*t77 - d9*spb24*t77 + t69*t77 + d8*spb24*t83 + d9*spb24*t83 + d6*s245*t76*T(2) + s245*t103*t76*T(2) + d11*s245*t77*T(2) - d13*t8*T(2) + d14*s25*t8*T(2) + d14*s45*t8*T(2) + d15*spb24*t8*T(2) + d16*spb24*t8*T(2) + s245*t77*t85*T(2) + d5*t83*T(8); 
complex<T> t7 = s45*t50*t57*t89 + s45*t104*t50*t9 + s45*t89*t91 + s45*t127*t97 + spb35*t29*(t21*t57 + t20*t58)*T(2) + s13*t21*t91*T(2) + d13*(s13*t8 - s45*t8 + s13*t83*T(2)); 
complex<T> t18 = s25*(-(d13*t8) + t50*t57*t89 + t104*t50*t9 + t89*t91 + t127*t97); 
SeriesC<T> result = t3*(*CI_users[0]->get_value(mc,ind,mu)) + t6*(*CI_users[1]->get_value(mc,ind,mu)) + t5*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t4*(*CI_users[4]->get_value(mc,ind,mu)) + t7*(*CI_users[5]->get_value(mc,ind,mu)) + t18*(*CI_users[6]->get_value(mc,ind,mu)) + t1*(*CI_users[7]->get_value(mc,ind,mu)) + t1*(*CI_users[8]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdmmpp_G_wCI::\
C4g1ph_phdmmpp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdmmpp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmmpp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:30:43 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> s45 = S(4,5);
complex<T> s24 = S(2,4);
complex<T> s14 = S(1,4);
complex<T> s15 = S(1,5);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s35 = S(3,5);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spb45); 
complex<T> t11 = square(spa23); 
complex<T> t23 = cube(spb45); 
complex<T> t24 = cube(spa23); 
complex<T> t25 = spb25*spb34; 
complex<T> t26 = spb24*spb35; 
complex<T> d1 = spb23*cube(s23 + s24)*T(3); d1 = T(1)/d1;
complex<T> d2 = (s23 + s24)*spb23*T(3); d2 = T(1)/d2;
complex<T> d3 = spb23*cube(s23 + s35)*T(3); d3 = T(1)/d3;
complex<T> d4 = (s23 + s35)*spb23*T(3); d4 = T(1)/d4;
complex<T> d6 = spb23*spb34; d6 = T(1)/d6;
complex<T> d8 = spb23*spb34*T(2); d8 = T(1)/d8;
complex<T> d10 = spb23*spb25*T(2); d10 = T(1)/d10;
complex<T> d11 = spb34*T(2); d11 = T(1)/d11;
complex<T> d12 = spb25*T(2); d12 = T(1)/d12;
complex<T> t3 = square(s24*spb45 - spa23*t26); 
complex<T> t4 = square(s35*spb45 - spa23*t26); 
complex<T> t12 = -t25; 
complex<T> t19 = d6*spa25; 
complex<T> t20 = d10*spa34; 
complex<T> t30 = d3*s35; 
complex<T> t35 = d8*spa25; 
complex<T> t39 = t24*t26; 
complex<T> t48 = spb45*t11; 
complex<T> d5 = spb23*t25; d5 = T(1)/d5;
complex<T> d7 = spb23*t25*T(2); d7 = T(1)/d7;
complex<T> d9 = t25*T(2); d9 = T(1)/d9;
complex<T> t1 = d1*(-(spa23*t3) + t12*t39 + s24*t25*t48) - d2*spa23*t10*T(11); 
complex<T> t9 = d1*(spa23*t3 + t25*t39 + s24*t12*t48) + d2*spa23*t10*T(11); 
complex<T> t14 = -(d5*T(3)); 
complex<T> t16 = d7*s234; 
complex<T> t33 = -(d5*s45); 
complex<T> t34 = d7*s245; 
complex<T> t37 = d3*t12*t39 - d3*spa23*t4 + t25*t30*t48 - d4*spa23*t10*T(11); 
complex<T> t38 = d5*t23; 
complex<T> t43 = d3*t25*t39 + d3*spa23*t4 + t12*t30*t48 + d4*spa23*t10*T(11); 
complex<T> t63 = t19*t23; 
complex<T> t8 = (-mH2 + s15)*t38*T(4); 
complex<T> t42 = (d9*mH2*spa23 + s235*t16)*t23; 
complex<T> t47 = (s345*t16 + mH2*t20)*t23; 
complex<T> t49 = -(d7*mH2*s45*t23) + s345*t23*t34; 
complex<T> t52 = t23*(s235*t34 + mH2*t35); 
complex<T> t55 = t23*t33 + s12*t38; 
complex<T> t58 = t23*t33 + s13*t38; 
complex<T> t61 = s14*t38 + t63; 
complex<T> t62 = t14*t23; 
complex<T> t36 = s15*t62 + s24*t38*T(3); 
complex<T> co1 = s24*t62; 
complex<T> co2 = -t63; 
complex<T> co3 = s45*t38*T(2); 
complex<T> co4 = -(s45*t23*t35); 
complex<T> co5 = d11*spa23*spa25*t23; 
complex<T> co6 = d12*spa23*spa34*t23; 
complex<T> co7 = -(s45*t20*t23); 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t37*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t58*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t61*(*CI_users[8]->get_value(mc,ind,mu)) + t36*(*CI_users[9]->get_value(mc,ind,mu)) + t55*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t49*(*CI_users[12]->get_value(mc,ind,mu)) + t52*(*CI_users[13]->get_value(mc,ind,mu)) + t42*(*CI_users[14]->get_value(mc,ind,mu)) + t47*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdmmpp_nf_wCI::\
C4g1ph_phdmmpp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c234, c15));
CI_users.push_back(new Cached_Bubble_Integral_User(c235, c14));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdmmpp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmmpp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:30:46 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = S(2,4);
complex<T> s35 = S(3,5);
complex<T> t9 = square(spa23); 
complex<T> t11 = square(spb45); 
complex<T> t18 = cube(spa23); 
complex<T> t19 = spb24*spb35; 
complex<T> t20 = spb25*spb34; 
complex<T> d1 = spb23*cube(s23 + s24)*T(9); d1 = T(1)/d1;
complex<T> d2 = (s23 + s24)*spb23*T(9); d2 = T(1)/d2;
complex<T> d3 = spb23*cube(s23 + s35)*T(9); d3 = T(1)/d3;
complex<T> d4 = (s23 + s35)*spb23*T(9); d4 = T(1)/d4;
complex<T> t3 = square(s24*spb45 - spa23*t19); 
complex<T> t4 = square(s35*spb45 - spa23*t19); 
complex<T> t10 = -t20; 
complex<T> t24 = d3*s35; 
complex<T> t25 = -(t11*T(2)); 
complex<T> t27 = t18*t19; 
complex<T> t31 = spb45*t9; 
complex<T> t1 = d2*spa23*t25 + d1*(t10*t27 - spa23*t3 + s24*t20*t31); 
complex<T> t8 = d1*(t20*t27 + spa23*t3 + s24*t10*t31) + d2*spa23*t11*T(2); 
complex<T> t26 = d3*t20*t27 + t10*t24*t31 + d3*spa23*t4 + d4*spa23*t11*T(2); 
complex<T> t29 = -(d3*t4); 
complex<T> t30 = d4*spa23*t25 + d3*t10*t27 + spa23*t29 + t20*t24*t31; 
SeriesC<T> result = t8*(*CI_users[0]->get_value(mc,ind,mu)) + t26*(*CI_users[1]->get_value(mc,ind,mu)) + t30*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdpmmp_G_wCI::\
C4g1ph_phdpmmp_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpmmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpmmp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:30:53 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> s25 = S(2,5);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s14 = S(1,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = S(3,5);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spb25); 
complex<T> t11 = square(spa34); 
complex<T> t23 = cube(spb25); 
complex<T> t24 = cube(spa34); 
complex<T> t25 = spb23*spb45; 
complex<T> t26 = spb24*spb35; 
complex<T> d1 = (s24 + s34)*spb34*T(3); d1 = T(1)/d1;
complex<T> d2 = spb34*cube(s24 + s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s34 + s35)*spb34*T(3); d3 = T(1)/d3;
complex<T> d4 = spb34*cube(s34 + s35)*T(3); d4 = T(1)/d4;
complex<T> d6 = spb23*spb34; d6 = T(1)/d6;
complex<T> d7 = spb23*spb34*T(2); d7 = T(1)/d7;
complex<T> d9 = spb34*spb45*T(2); d9 = T(1)/d9;
complex<T> d11 = spb45*T(2); d11 = T(1)/d11;
complex<T> d12 = spb23*T(2); d12 = T(1)/d12;
complex<T> t3 = square(s24*spb25 - spa34*t26); 
complex<T> t4 = square(s35*spb25 - spa34*t26); 
complex<T> t12 = -t25; 
complex<T> t15 = d6*spa45; 
complex<T> t20 = d9*spa23; 
complex<T> t30 = d4*s35; 
complex<T> t35 = d7*spa45; 
complex<T> t39 = t24*t26; 
complex<T> t51 = spb25*t11; 
complex<T> t53 = -(mH2*t23); 
complex<T> t63 = s24*t23; 
complex<T> d5 = spb34*t25; d5 = T(1)/d5;
complex<T> d8 = spb34*t25*T(2); d8 = T(1)/d8;
complex<T> d10 = t25*T(2); d10 = T(1)/d10;
complex<T> t1 = d2*(-(spa34*t3) + t12*t39 + s24*t25*t51) - d1*spa34*t10*T(11); 
complex<T> t9 = d2*(spa34*t3 + t25*t39 + s24*t12*t51) + d1*spa34*t10*T(11); 
complex<T> t14 = -(d5*T(3)); 
complex<T> t17 = d8*s234; 
complex<T> t33 = d5*T(3); 
complex<T> t34 = d8*s245; 
complex<T> t37 = d4*t12*t39 - d4*spa34*t4 + t25*t30*t51 - d3*spa34*t10*T(11); 
complex<T> t43 = d4*t25*t39 + d4*spa34*t4 + t12*t30*t51 + d3*spa34*t10*T(11); 
complex<T> t44 = d5*t23; 
complex<T> t59 = t15*t23; 
complex<T> t7 = (-mH2 + s15)*t44*T(4); 
complex<T> t36 = (s14 - s25)*t44; 
complex<T> t42 = s12*t44 + t59; 
complex<T> t48 = (d10*mH2*spa34 + s345*t17)*t23; 
complex<T> t50 = (s235*t17 + mH2*t20)*t23; 
complex<T> t52 = s15*t14*t23 + t33*t63; 
complex<T> t55 = t23*(s345*t34 + mH2*t35); 
complex<T> t58 = s13*t44 + t59; 
complex<T> t61 = s235*t23*t34 + d8*s25*t53; 
complex<T> co1 = t14*t63; 
complex<T> co2 = s25*t44; 
complex<T> co3 = -(t59*T(2)); 
complex<T> co4 = -(s25*t23*t35); 
complex<T> co5 = -(s25*t20*t23); 
complex<T> co6 = d11*spa23*spa34*t23; 
complex<T> co7 = d12*spa34*spa45*t23; 
SeriesC<T> result = t9*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t37*(*CI_users[2]->get_value(mc,ind,mu)) + t43*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t58*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t36*(*CI_users[8]->get_value(mc,ind,mu)) + t52*(*CI_users[9]->get_value(mc,ind,mu)) + t42*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t55*(*CI_users[12]->get_value(mc,ind,mu)) + t61*(*CI_users[13]->get_value(mc,ind,mu)) + t50*(*CI_users[14]->get_value(mc,ind,mu)) + t48*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdpmmp_nf_wCI::\
C4g1ph_phdpmmp_nf_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpmmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpmmp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:30:56 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s24 = S(2,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = S(3,5);
complex<T> t9 = square(spa34); 
complex<T> t10 = square(spb25); 
complex<T> t18 = cube(spa34); 
complex<T> t19 = spb24*spb35; 
complex<T> t20 = spb23*spb45; 
complex<T> d1 = (s24 + s34)*spb34*T(9); d1 = T(1)/d1;
complex<T> d2 = spb34*cube(s24 + s34)*T(9); d2 = T(1)/d2;
complex<T> d3 = (s34 + s35)*spb34*T(9); d3 = T(1)/d3;
complex<T> d4 = spb34*cube(s34 + s35)*T(9); d4 = T(1)/d4;
complex<T> t3 = square(s24*spb25 - spa34*t19); 
complex<T> t4 = square(s35*spb25 - spa34*t19); 
complex<T> t11 = -t20; 
complex<T> t24 = d4*s35; 
complex<T> t25 = -(t10*T(2)); 
complex<T> t27 = t18*t19; 
complex<T> t31 = spb25*t9; 
complex<T> t1 = d1*spa34*t25 + d2*(t11*t27 - spa34*t3 + s24*t20*t31); 
complex<T> t8 = d2*(t20*t27 + spa34*t3 + s24*t11*t31) + d1*spa34*t10*T(2); 
complex<T> t26 = d4*t20*t27 + t11*t24*t31 + d4*spa34*t4 + d3*spa34*t10*T(2); 
complex<T> t29 = -(d4*t4); 
complex<T> t30 = d3*spa34*t25 + d4*t11*t27 + spa34*t29 + t20*t24*t31; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t8*(*CI_users[1]->get_value(mc,ind,mu)) + t26*(*CI_users[2]->get_value(mc,ind,mu)) + t30*(*CI_users[3]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdppmm_G_wCI::\
C4g1ph_phdppmm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdppmm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdppmm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:31:05 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s14 = S(1,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = S(3,5);
complex<T> s23 = S(2,3);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spb23); 
complex<T> t11 = square(spa45); 
complex<T> t12 = -(spb25*spb34); 
complex<T> t23 = cube(spb23); 
complex<T> t24 = cube(spa45); 
complex<T> t25 = spb24*spb35; 
complex<T> d1 = (s24 + s45)*spb45*T(3); d1 = T(1)/d1;
complex<T> d2 = spb45*cube(s24 + s45)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s35 + s45)*spb45*T(3); d3 = T(1)/d3;
complex<T> d4 = spb45*cube(s35 + s45)*T(3); d4 = T(1)/d4;
complex<T> d5 = spb25*spb34*spb45; d5 = T(1)/d5;
complex<T> d6 = spb25*spb34; d6 = T(1)/d6;
complex<T> d7 = spb34*spb45; d7 = T(1)/d7;
complex<T> d8 = spb25*spb34*T(2); d8 = T(1)/d8;
complex<T> d9 = spb25*spb34*spb45*T(2); d9 = T(1)/d9;
complex<T> d10 = spb34*spb45*T(2); d10 = T(1)/d10;
complex<T> d11 = spb25*spb45*T(2); d11 = T(1)/d11;
complex<T> d12 = spb34*T(2); d12 = T(1)/d12;
complex<T> d13 = spb25*T(2); d13 = T(1)/d13;
complex<T> t3 = square(s24*spb23 - spa45*t25); 
complex<T> t4 = square(s35*spb23 - spa45*t25); 
complex<T> t14 = -(d5*T(3)); 
complex<T> t16 = d9*s234; 
complex<T> t19 = d10*spa25; 
complex<T> t20 = d11*spa34; 
complex<T> t29 = d4*s35; 
complex<T> t32 = d5*T(3); 
complex<T> t33 = d9*s245; 
complex<T> t34 = d7*spa25; 
complex<T> t37 = spa45*t23; 
complex<T> t49 = spb23*t11; 
complex<T> t56 = spb25*t24; 
complex<T> t58 = d5*t23; 
complex<T> t67 = s24*t23; 
complex<T> t7 = (-mH2 + s15)*t58*T(4); 
complex<T> t9 = d2*(spa45*t3 + s24*t12*t49 + spb34*t25*t56) + d1*spa45*t10*T(11); 
complex<T> t42 = -(d9*mH2*s23*t23) + s235*t16*t23; 
complex<T> t43 = d4*spa45*t4 + t12*t29*t49 + d4*spb34*t25*t56 + d3*spa45*t10*T(11); 
complex<T> t48 = (s345*t16 + mH2*t20)*t23; 
complex<T> t51 = s15*t14*t23 + t32*t67; 
complex<T> t54 = t23*(mH2*t19 + s235*t33); 
complex<T> t60 = s345*t23*t33 + d8*mH2*t37; 
complex<T> t61 = d6*t37; 
complex<T> t69 = t23*t34; 
complex<T> t70 = spb25*t49; 
complex<T> t1 = d2*(t12*t24*t25 - spa45*t3 + s24*spb34*t70) - d1*spa45*t10*T(11); 
complex<T> t35 = s12*t58 + t61; 
complex<T> t36 = d4*t12*t24*t25 - d4*spa45*t4 + spb34*t29*t70 - d3*spa45*t10*T(11); 
complex<T> t57 = s14*t58 + t69; 
complex<T> t62 = s13*t58 + t61; 
complex<T> co1 = t14*t67; 
complex<T> co2 = -t69; 
complex<T> co3 = -(t61*T(2)); 
complex<T> co4 = d12*spa25*t37; 
complex<T> co5 = -(s23*t19*t23); 
complex<T> co6 = -(s23*t20*t23); 
complex<T> co7 = d13*spa34*t37; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t36*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t62*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t57*(*CI_users[8]->get_value(mc,ind,mu)) + t51*(*CI_users[9]->get_value(mc,ind,mu)) + t35*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t60*(*CI_users[12]->get_value(mc,ind,mu)) + t54*(*CI_users[13]->get_value(mc,ind,mu)) + t42*(*CI_users[14]->get_value(mc,ind,mu)) + t48*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdppmm_nf_wCI::\
C4g1ph_phdppmm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c345, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdppmm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdppmm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:31:08 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> s24 = S(2,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = S(3,5);
complex<T> t9 = square(spa45); 
complex<T> t10 = square(spb23); 
complex<T> t18 = cube(spa45); 
complex<T> t19 = spb24*spb35; 
complex<T> t20 = spb25*spb34; 
complex<T> d1 = (s24 + s45)*spb45*T(9); d1 = T(1)/d1;
complex<T> d2 = spb45*cube(s24 + s45)*T(9); d2 = T(1)/d2;
complex<T> d3 = (s35 + s45)*spb45*T(9); d3 = T(1)/d3;
complex<T> d4 = spb45*cube(s35 + s45)*T(9); d4 = T(1)/d4;
complex<T> t3 = square(s24*spb23 - spa45*t19); 
complex<T> t4 = square(s35*spb23 - spa45*t19); 
complex<T> t11 = -t20; 
complex<T> t24 = d4*s35; 
complex<T> t25 = -(t10*T(2)); 
complex<T> t27 = t18*t19; 
complex<T> t31 = spb23*t9; 
complex<T> t1 = d1*spa45*t25 + d2*(t11*t27 - spa45*t3 + s24*t20*t31); 
complex<T> t8 = d2*(t20*t27 + spa45*t3 + s24*t11*t31) + d1*spa45*t10*T(2); 
complex<T> t26 = d4*t20*t27 + t11*t24*t31 + d4*spa45*t4 + d3*spa45*t10*T(2); 
complex<T> t29 = -(d4*t4); 
complex<T> t30 = d3*spa45*t25 + d4*t11*t27 + spa45*t29 + t20*t24*t31; 
SeriesC<T> result = t8*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t30*(*CI_users[2]->get_value(mc,ind,mu)) + t26*(*CI_users[3]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdmppm_G_wCI::\
C4g1ph_phdmppm_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdmppm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmppm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:31:17 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> s34 = S(3,4);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s14 = S(1,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = S(3,5);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spb34); 
complex<T> t11 = square(spa25); 
complex<T> t12 = -(spb23*spb45); 
complex<T> t23 = cube(spb34); 
complex<T> t24 = cube(spa25); 
complex<T> t25 = spb24*spb35; 
complex<T> d1 = (s25 + s35)*spb25*T(3); d1 = T(1)/d1;
complex<T> d2 = spb25*cube(s25 + s35)*T(3); d2 = T(1)/d2;
complex<T> d3 = (s24 + s25)*spb25*T(3); d3 = T(1)/d3;
complex<T> d4 = spb25*cube(s24 + s25)*T(3); d4 = T(1)/d4;
complex<T> d5 = spb23*spb25*spb45; d5 = T(1)/d5;
complex<T> d6 = spb23*spb25; d6 = T(1)/d6;
complex<T> d7 = spb23*spb45; d7 = T(1)/d7;
complex<T> d8 = spb23*spb25*T(2); d8 = T(1)/d8;
complex<T> d9 = spb23*spb25*spb45*T(2); d9 = T(1)/d9;
complex<T> d10 = spb23*spb45*T(2); d10 = T(1)/d10;
complex<T> d11 = spb25*spb45*T(2); d11 = T(1)/d11;
complex<T> d12 = spb23*T(2); d12 = T(1)/d12;
complex<T> d13 = spb45*T(2); d13 = T(1)/d13;
complex<T> t3 = square(s24*spb34 - spa25*t25); 
complex<T> t4 = square(s35*spb34 - spa25*t25); 
complex<T> t14 = -(d5*T(3)); 
complex<T> t15 = d6*spa45; 
complex<T> t17 = d9*s234; 
complex<T> t20 = d11*spa23; 
complex<T> t28 = d4*s24; 
complex<T> t29 = d2*s35; 
complex<T> t31 = d1*T(11); 
complex<T> t32 = d5*T(3); 
complex<T> t33 = d9*s245; 
complex<T> t34 = d8*spa45; 
complex<T> t37 = d5*t23; 
complex<T> t38 = t24*t25; 
complex<T> t43 = spa25*t23; 
complex<T> t49 = spb34*t11; 
complex<T> t62 = s24*t23; 
complex<T> t8 = (-mH2 + s15)*t37*T(4); 
complex<T> t36 = d4*spa25*t3 + d4*spb23*spb45*t38 + t12*t28*t49 + d3*spa25*t10*T(11); 
complex<T> t41 = -(d9*mH2*s34*t23) + s345*t17*t23; 
complex<T> t42 = spa25*t10*t31 + d2*spb23*spb45*t38 + d2*spa25*t4 + t12*t29*t49; 
complex<T> t46 = (s235*t17 + mH2*t20)*t23; 
complex<T> t48 = s15*t14*t23 + t32*t62; 
complex<T> t50 = t23*(s345*t33 + mH2*t34); 
complex<T> t57 = t15*t23; 
complex<T> t59 = s235*t23*t33 + d10*mH2*t43; 
complex<T> t66 = d7*t43; 
complex<T> t67 = spb23*t49; 
complex<T> t1 = d2*t12*t38 - d2*spa25*t4 + spb45*t29*t67 - d1*spa25*t10*T(11); 
complex<T> t9 = -(d4*spa25*t3) + d4*t12*t38 + spb45*t28*t67 - d3*spa25*t10*T(11); 
complex<T> t35 = s14*t37 + t66; 
complex<T> t53 = s12*t37 + t57; 
complex<T> t56 = s13*t37 + t57; 
complex<T> co1 = t14*t62; 
complex<T> co2 = -t66; 
complex<T> co3 = -(t57*T(2)); 
complex<T> co4 = d12*spa45*t43; 
complex<T> co5 = d13*spa23*t43; 
complex<T> co6 = -(s34*t20*t23); 
complex<T> co7 = -(s34*t23*t34); 
SeriesC<T> result = t42*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t9*(*CI_users[2]->get_value(mc,ind,mu)) + t36*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t56*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t35*(*CI_users[8]->get_value(mc,ind,mu)) + t48*(*CI_users[9]->get_value(mc,ind,mu)) + t53*(*CI_users[10]->get_value(mc,ind,mu)) + co3*(*CI_users[11]->get_value(mc,ind,mu)) + t50*(*CI_users[12]->get_value(mc,ind,mu)) + t59*(*CI_users[13]->get_value(mc,ind,mu)) + t46*(*CI_users[14]->get_value(mc,ind,mu)) + t41*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdmppm_nf_wCI::\
C4g1ph_phdmppm_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdmppm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmppm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:31:20 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> s24 = S(2,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = S(3,5);
complex<T> t9 = square(spa25); 
complex<T> t10 = square(spb34); 
complex<T> t18 = cube(spa25); 
complex<T> t19 = spb24*spb35; 
complex<T> t20 = spb23*spb45; 
complex<T> d1 = (s25 + s35)*spb25*T(9); d1 = T(1)/d1;
complex<T> d2 = spb25*cube(s25 + s35)*T(9); d2 = T(1)/d2;
complex<T> d3 = (s24 + s25)*spb25*T(9); d3 = T(1)/d3;
complex<T> d4 = spb25*cube(s24 + s25)*T(9); d4 = T(1)/d4;
complex<T> t3 = square(s24*spb34 - spa25*t19); 
complex<T> t4 = square(s35*spb34 - spa25*t19); 
complex<T> t11 = -t20; 
complex<T> t23 = d4*s24; 
complex<T> t24 = d2*s35; 
complex<T> t27 = t18*t19; 
complex<T> t28 = t10*T(2); 
complex<T> t31 = spb34*t9; 
complex<T> t1 = d2*t11*t27 + t20*t24*t31 - d2*spa25*t4 - d1*spa25*t10*T(2); 
complex<T> t8 = d4*t11*t27 - d4*spa25*t3 + t20*t23*t31 - d3*spa25*t10*T(2); 
complex<T> t26 = d4*t20*t27 + d3*spa25*t28 + d4*spa25*t3 + t11*t23*t31; 
complex<T> t30 = d2*t20*t27 + d1*spa25*t28 + t11*t24*t31 + d2*spa25*t4; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t30*(*CI_users[1]->get_value(mc,ind,mu)) + t26*(*CI_users[2]->get_value(mc,ind,mu)) + t8*(*CI_users[3]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdmpmp_G_wCI::\
C4g1ph_phdmpmp_G_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c1, c5, c234));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdmpmp_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmpmp G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:31:58 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa24 = SPA(2,4);
complex<T> s12 = S(1,2);
complex<T> s14 = S(1,4);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s13 = S(1,3);
complex<T> s45 = -(spa45*spb45);
complex<T> s24 = -(spa24*spb24);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t18 = square(spb35); 
complex<T> t21 = square(spa24); 
complex<T> t34 = square(square(spb35)); 
complex<T> t35 = spb23*spb45; 
complex<T> t36 = spb25*spb34; 
complex<T> t37 = -(spa24*T(2)); 
complex<T> t39 = s23*s34; 
complex<T> t40 = s25*s45; 
complex<T> t45 = cube(spa24); 
complex<T> t46 = square(spb23); 
complex<T> t47 = square(spb25); 
complex<T> t48 = square(spb34); 
complex<T> t49 = square(spb45); 
complex<T> t54 = spa24*T(2); 
complex<T> d1 = -((s24 + s34)*spb24); d1 = T(1)/d1;
complex<T> d2 = -((s24 + s34)*cube(spb24)); d2 = T(1)/d2;
complex<T> d3 = square(s24 + s34)*square(spb24); d3 = T(1)/d3;
complex<T> d4 = -(spb24*cube(s24 + s34)*T(3)); d4 = T(1)/d4;
complex<T> d5 = -(spb24*cube(s23 + s24)*T(3)); d5 = T(1)/d5;
complex<T> d6 = -((s23 + s24)*spb24); d6 = T(1)/d6;
complex<T> d7 = -((s23 + s24)*cube(spb24)); d7 = T(1)/d7;
complex<T> d8 = square(s23 + s24)*square(spb24); d8 = T(1)/d8;
complex<T> d9 = -(spb24*cube(s24 + s45)*T(3)); d9 = T(1)/d9;
complex<T> d10 = -((s24 + s25)*spb24); d10 = T(1)/d10;
complex<T> d11 = -((s24 + s45)*spb24); d11 = T(1)/d11;
complex<T> d12 = -((s24 + s25)*cube(spb24)); d12 = T(1)/d12;
complex<T> d13 = -((s24 + s45)*cube(spb24)); d13 = T(1)/d13;
complex<T> d14 = square(s24 + s25)*square(spb24); d14 = T(1)/d14;
complex<T> d15 = square(s24 + s45)*square(spb24); d15 = T(1)/d15;
complex<T> d16 = -(spb24*cube(s24 + s25)*T(3)); d16 = T(1)/d16;
complex<T> d18 = square(spb24); d18 = T(1)/d18;
complex<T> d19 = square(square(spb24)); d19 = T(1)/d19;
complex<T> d20 = spb24; d20 = T(1)/d20;
complex<T> d21 = cube(spb24); d21 = T(1)/d21;
complex<T> d29 = square(square(spb24))*T(2); d29 = T(1)/d29;
complex<T> d30 = spb23*spb34*T(2); d30 = T(1)/d30;
complex<T> d31 = spb34*spb45*T(2); d31 = T(1)/d31;
complex<T> d32 = spb25*spb45*T(2); d32 = T(1)/d32;
complex<T> d33 = spb23*spb25*T(2); d33 = T(1)/d33;
complex<T> t19 = -t36; 
complex<T> t26 = -(d18*T(4)); 
complex<T> t27 = -t39; 
complex<T> t28 = -t40; 
complex<T> t32 = d32*spa23; 
complex<T> t56 = d19*t35; 
complex<T> t57 = t47*t48; 
complex<T> t62 = -(t18*T(4)); 
complex<T> t63 = t21*t35; 
complex<T> t64 = d29*t36; 
complex<T> t68 = -(t36*T(2)); 
complex<T> t69 = spa24*t18; 
complex<T> t70 = t46*t49; 
complex<T> t73 = t36*T(2); 
complex<T> t74 = mH2*t34; 
complex<T> t98 = spa25*t34; 
complex<T> d17 = t35*t36; d17 = T(1)/d17;
complex<T> d22 = spb23*t36; d22 = T(1)/d22;
complex<T> d23 = spb34*t35; d23 = T(1)/d23;
complex<T> d26 = spb34*t35*T(2); d26 = T(1)/d26;
complex<T> d28 = spb25*t35*T(2); d28 = T(1)/d28;
complex<T> t29 = -(d17*T(3)); 
complex<T> t30 = d22*spa45; 
complex<T> t41 = d23*spa25; 
complex<T> t43 = spa34*t32; 
complex<T> t55 = d17*t34; 
complex<T> t58 = d18*t27; 
complex<T> t61 = s23*(t18*t26 + t56*t73); 
complex<T> t66 = s34*(t18*t26 + t56*t73); 
complex<T> t75 = t45*t70; 
complex<T> t86 = t35*t64; 
complex<T> t90 = t45*t57; 
complex<T> d24 = spb23*t73; d24 = T(1)/d24;
complex<T> d25 = t35*t73; d25 = T(1)/d25;
complex<T> d27 = spb45*t73; d27 = T(1)/d27;
complex<T> t9 = (-mH2 + s15)*t55*T(4); 
complex<T> t10 = s25*t18*t26 - t34*t41 + s25*t56*t73; 
complex<T> t11 = t18*t58 + t39*t86; 
complex<T> t12 = t34*t43 + t18*t58 + t39*t86; 
complex<T> t15 = s15*t29*t34 + d21*t35*t36*t37 + s15*t56*t68 + s24*t55*T(3) + d18*s15*t18*T(4) + d20*t69*T(4); 
complex<T> t16 = s24*t29*t34 + d21*t35*t36*t54 + d20*spa24*t62; 
complex<T> t17 = d13*t35*t36*t37 + d15*t19*t63 + d9*t90*T(2) + d11*t69*T(4); 
complex<T> t31 = d25*s234; 
complex<T> t42 = d25*s245; 
complex<T> t44 = d12*t35*t36*t54 + d13*t35*t36*t54 + d10*spa24*t62 + d11*spa24*t62 + d14*t36*t63 + d15*t36*t63 - d16*t75*T(2) - d9*t90*T(2); 
complex<T> t52 = d18*t18*t28 + t40*t86; 
complex<T> t60 = d12*t35*t36*t37 + d14*t19*t63 + d16*t75*T(2) + d10*t69*T(4); 
complex<T> t67 = d2*t35*t36*t37 + d3*t19*t63 + d4*t75*T(2) + d1*t69*T(4); 
complex<T> t72 = d7*t35*t36*t37 + d8*t19*t63 + d5*t90*T(2) + d6*t69*T(4); 
complex<T> t76 = t34*t41 + s14*t55; 
complex<T> t77 = d2*t35*t36*t54 + d7*t35*t36*t54 + d1*spa24*t62 + d6*spa24*t62 + d3*t36*t63 + d8*t36*t63 - d4*t75*T(2) - d5*t90*T(2); 
complex<T> t83 = t30*t34; 
complex<T> t13 = s45*t18*t26 + s13*t55 + s13*t56*t68 + s45*t56*t73 + t83 + d18*s13*t18*T(4); 
complex<T> t14 = t52 + d30*spa45*t98; 
complex<T> t51 = s12*t55 + t83; 
complex<T> t59 = s235*t31*t34 + d27*spa23*t74; 
complex<T> t65 = s345*t31*t34 + d28*spa34*t74; 
complex<T> t71 = s235*t34*t42 + d26*spa25*t74; 
complex<T> t81 = s345*t34*t42 + d24*spa45*t74; 
complex<T> co1 = -(t83*T(2)); 
complex<T> co2 = d31*spa23*t98; 
complex<T> co3 = d33*spa34*spa45*t34; 
SeriesC<T> result = t67*(*CI_users[0]->get_value(mc,ind,mu)) + t77*(*CI_users[1]->get_value(mc,ind,mu)) + t44*(*CI_users[2]->get_value(mc,ind,mu)) + t17*(*CI_users[3]->get_value(mc,ind,mu)) + t72*(*CI_users[4]->get_value(mc,ind,mu)) + t60*(*CI_users[5]->get_value(mc,ind,mu)) + t9*(*CI_users[6]->get_value(mc,ind,mu)) + t61*(*CI_users[7]->get_value(mc,ind,mu)) + t16*(*CI_users[8]->get_value(mc,ind,mu)) + t13*(*CI_users[9]->get_value(mc,ind,mu)) + t10*(*CI_users[10]->get_value(mc,ind,mu)) + t76*(*CI_users[11]->get_value(mc,ind,mu)) + t15*(*CI_users[12]->get_value(mc,ind,mu)) + t66*(*CI_users[13]->get_value(mc,ind,mu)) + t51*(*CI_users[14]->get_value(mc,ind,mu)) + co1*(*CI_users[15]->get_value(mc,ind,mu)) + t81*(*CI_users[16]->get_value(mc,ind,mu)) + t71*(*CI_users[17]->get_value(mc,ind,mu)) + t59*(*CI_users[18]->get_value(mc,ind,mu)) + t65*(*CI_users[19]->get_value(mc,ind,mu)) + t11*(*CI_users[20]->get_value(mc,ind,mu)) + t14*(*CI_users[21]->get_value(mc,ind,mu)) + co2*(*CI_users[22]->get_value(mc,ind,mu)) + t12*(*CI_users[23]->get_value(mc,ind,mu)) + t52*(*CI_users[24]->get_value(mc,ind,mu)) + co3*(*CI_users[25]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdmpmp_nf_wCI::\
C4g1ph_phdmpmp_nf_wCI
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
CI_users.push_back(new Cached_Bubble_Integral_User(c245, c13));
CI_users.push_back(new Cached_Bubble_Integral_User(c25, c134));
CI_users.push_back(new Cached_Bubble_Integral_User(c34, c125));
CI_users.push_back(new Cached_Bubble_Integral_User(c45, c123));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdmpmp_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, m, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmpmp nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:32:16 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> s23 = S(2,3);
complex<T> s25 = S(2,5);
complex<T> s34 = S(3,4);
complex<T> s13 = S(1,3);
complex<T> s45 = S(4,5);
complex<T> s15 = S(1,5);
complex<T> s24 = -(spa24*spb24);
complex<T> t15 = square(spb35); 
complex<T> t16 = square(spa24); 
complex<T> t19 = -s24 - s25; 
complex<T> t20 = -s24 - s34; 
complex<T> t21 = -s24 - s45; 
complex<T> t22 = -s23 - s24; 
complex<T> t23 = s25*s45; 
complex<T> t28 = spb23*spb45; 
complex<T> t29 = spb25*spb34; 
complex<T> t30 = -(spa24*T(2)); 
complex<T> t36 = cube(spa24); 
complex<T> t37 = square(spb23); 
complex<T> t38 = square(spb25); 
complex<T> t39 = square(spb34); 
complex<T> t40 = square(spb45); 
complex<T> d17 = square(spb24)*T(3); d17 = T(1)/d17;
complex<T> d18 = square(square(spb24))*T(3); d18 = T(1)/d18;
complex<T> d19 = spb24*T(3); d19 = T(1)/d19;
complex<T> d20 = cube(spb24)*T(3); d20 = T(1)/d20;
complex<T> d21 = square(spb24)*T(12); d21 = T(1)/d21;
complex<T> d22 = square(square(spb24))*T(6); d22 = T(1)/d22;
complex<T> t17 = -t29; 
complex<T> t34 = d21*t23; 
complex<T> t44 = t28*t29; 
complex<T> t46 = -(d18*T(2)); 
complex<T> t51 = spa24*t15; 
complex<T> t55 = -(t36*T(2)); 
complex<T> t64 = t38*t39; 
complex<T> t65 = t37*t40; 
complex<T> d1 = spb24*t20*T(3); d1 = T(1)/d1;
complex<T> d2 = t20*cube(spb24)*T(3); d2 = T(1)/d2;
complex<T> d3 = square(spb24)*square(t20)*T(3); d3 = T(1)/d3;
complex<T> d4 = spb24*cube(t20)*T(9); d4 = T(1)/d4;
complex<T> d5 = spb24*cube(t22)*T(9); d5 = T(1)/d5;
complex<T> d6 = spb24*t22*T(3); d6 = T(1)/d6;
complex<T> d7 = t22*cube(spb24)*T(3); d7 = T(1)/d7;
complex<T> d8 = square(spb24)*square(t22)*T(3); d8 = T(1)/d8;
complex<T> d9 = spb24*cube(t21)*T(9); d9 = T(1)/d9;
complex<T> d10 = spb24*t19*T(3); d10 = T(1)/d10;
complex<T> d11 = spb24*t21*T(3); d11 = T(1)/d11;
complex<T> d12 = t19*cube(spb24)*T(3); d12 = T(1)/d12;
complex<T> d13 = t21*cube(spb24)*T(3); d13 = T(1)/d13;
complex<T> d14 = square(spb24)*square(t19)*T(3); d14 = T(1)/d14;
complex<T> d15 = square(spb24)*square(t21)*T(3); d15 = T(1)/d15;
complex<T> d16 = spb24*cube(t19)*T(9); d16 = T(1)/d16;
complex<T> t12 = s25*(d17*t15 + t44*t46); 
complex<T> t13 = d20*t30*t44 + d19*t51; 
complex<T> t25 = -(d16*T(2)); 
complex<T> t42 = d16*T(2); 
complex<T> t43 = s23*(d17*t15 + t44*t46); 
complex<T> t45 = -t51; 
complex<T> t49 = s34*(d17*t15 + t44*t46); 
complex<T> t50 = t44*T(2); 
complex<T> t54 = t17*t28; 
complex<T> t7 = s23*s34*(d21*t15 + d22*t54); 
complex<T> t27 = d17*(-s13 + s45)*t15 + s45*t44*t46 + d18*s13*t50; 
complex<T> t41 = t15*t34 + d22*t23*t54; 
complex<T> t48 = d12*t30*t44 + d13*t30*t44 + d10*t51 + d11*t51 + d14*t16*t54 + d15*t16*t54 + t36*t42*t65 + d9*t36*t64*T(2); 
complex<T> t58 = d2*t30*t44 + d7*t30*t44 + d1*t51 + d6*t51 + d3*t16*t54 + d8*t16*t54 + d5*t36*t64*T(2) + d4*t36*t65*T(2); 
complex<T> t61 = spa24*t50; 
complex<T> t14 = -(d17*s15*t15) + d19*t45 + d18*s15*t50 + d20*t61; 
complex<T> t35 = d14*t16*t44 + d10*t45 + d12*t61 + t25*t36*t65; 
complex<T> t53 = d15*t16*t44 + d11*t45 + d13*t61 + d9*t55*t64; 
complex<T> t60 = d8*t16*t44 + d6*t45 + d7*t61 + d5*t55*t64; 
complex<T> t62 = d3*t16*t44 + d1*t45 + d2*t61 + d4*t55*t65; 
SeriesC<T> result = t62*(*CI_users[0]->get_value(mc,ind,mu)) + t58*(*CI_users[1]->get_value(mc,ind,mu)) + t48*(*CI_users[2]->get_value(mc,ind,mu)) + t53*(*CI_users[3]->get_value(mc,ind,mu)) + t60*(*CI_users[4]->get_value(mc,ind,mu)) + t35*(*CI_users[5]->get_value(mc,ind,mu)) + t43*(*CI_users[6]->get_value(mc,ind,mu)) + t13*(*CI_users[7]->get_value(mc,ind,mu)) + t27*(*CI_users[8]->get_value(mc,ind,mu)) + t12*(*CI_users[9]->get_value(mc,ind,mu)) + t14*(*CI_users[10]->get_value(mc,ind,mu)) + t49*(*CI_users[11]->get_value(mc,ind,mu)) + t7*(*CI_users[12]->get_value(mc,ind,mu)) + t41*(*CI_users[13]->get_value(mc,ind,mu)) + t7*(*CI_users[14]->get_value(mc,ind,mu)) + t41*(*CI_users[15]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdpmpm_G_wCI::\
C4g1ph_phdpmpm_G_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpmpm_G_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpmpm G");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:32:53 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa35 = SPA(3,5);
complex<T> s13 = S(1,3);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s14 = S(1,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s34 = -(spa34*spb34);
complex<T> s12 = S(1,2);
complex<T> s45 = -(spa45*spb45);
complex<T> s35 = -(spa35*spb35);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t17 = square(spb24); 
complex<T> t20 = square(spa35); 
complex<T> t33 = square(square(spb24)); 
complex<T> t34 = spb23*spb45; 
complex<T> t35 = spb25*spb34; 
complex<T> t36 = -(spa35*T(2)); 
complex<T> t38 = s23*s25; 
complex<T> t39 = s34*s45; 
complex<T> t43 = cube(spa35); 
complex<T> t44 = square(spb23); 
complex<T> t45 = square(spb25); 
complex<T> t46 = square(spb34); 
complex<T> t47 = square(spb45); 
complex<T> t52 = spa35*T(2); 
complex<T> d1 = -((s25 + s35)*spb35); d1 = T(1)/d1;
complex<T> d2 = -((s25 + s35)*cube(spb35)); d2 = T(1)/d2;
complex<T> d3 = square(s25 + s35)*square(spb35); d3 = T(1)/d3;
complex<T> d4 = -(spb35*cube(s25 + s35)*T(3)); d4 = T(1)/d4;
complex<T> d5 = -((s23 + s35)*spb35); d5 = T(1)/d5;
complex<T> d6 = -(spb35*cube(s23 + s35)*T(3)); d6 = T(1)/d6;
complex<T> d7 = -((s23 + s35)*cube(spb35)); d7 = T(1)/d7;
complex<T> d8 = square(s23 + s35)*square(spb35); d8 = T(1)/d8;
complex<T> d9 = -((s35 + s45)*spb35); d9 = T(1)/d9;
complex<T> d10 = -(spb35*cube(s35 + s45)*T(3)); d10 = T(1)/d10;
complex<T> d11 = -((s35 + s45)*cube(spb35)); d11 = T(1)/d11;
complex<T> d12 = square(s35 + s45)*square(spb35); d12 = T(1)/d12;
complex<T> d13 = -((s34 + s35)*spb35); d13 = T(1)/d13;
complex<T> d14 = -((s34 + s35)*cube(spb35)); d14 = T(1)/d14;
complex<T> d15 = square(s34 + s35)*square(spb35); d15 = T(1)/d15;
complex<T> d16 = -(spb35*cube(s34 + s35)*T(3)); d16 = T(1)/d16;
complex<T> d18 = square(spb35); d18 = T(1)/d18;
complex<T> d19 = square(square(spb35)); d19 = T(1)/d19;
complex<T> d27 = spb23*spb34*T(2); d27 = T(1)/d27;
complex<T> d28 = spb34*spb45*T(2); d28 = T(1)/d28;
complex<T> d29 = square(square(spb35))*T(2); d29 = T(1)/d29;
complex<T> d30 = spb25*spb45*T(2); d30 = T(1)/d30;
complex<T> d31 = spb23*spb25*T(2); d31 = T(1)/d31;
complex<T> t18 = -t35; 
complex<T> t25 = -(d18*T(4)); 
complex<T> t26 = -t38; 
complex<T> t27 = -t39; 
complex<T> t31 = d28*spa23; 
complex<T> t54 = d19*t34; 
complex<T> t55 = t45*t46; 
complex<T> t60 = -(t17*T(4)); 
complex<T> t61 = t20*t34; 
complex<T> t62 = d29*t35; 
complex<T> t66 = -(t35*T(2)); 
complex<T> t67 = spa35*t17; 
complex<T> t68 = t44*t47; 
complex<T> t71 = t35*T(2); 
complex<T> t72 = mH2*t33; 
complex<T> t81 = d18*t17; 
complex<T> t97 = spa34*t33; 
complex<T> d17 = t34*t35; d17 = T(1)/d17;
complex<T> d20 = spb23*t35; d20 = T(1)/d20;
complex<T> d21 = spb34*t34; d21 = T(1)/d21;
complex<T> d24 = spb34*t34*T(2); d24 = T(1)/d24;
complex<T> d26 = spb25*t34*T(2); d26 = T(1)/d26;
complex<T> t28 = -(d17*T(3)); 
complex<T> t29 = d20*spa45; 
complex<T> t32 = d21*spa25; 
complex<T> t41 = spa25*t31; 
complex<T> t53 = d17*t33; 
complex<T> t56 = d18*t26; 
complex<T> t59 = s23*(t17*t25 + t54*t71); 
complex<T> t64 = s34*(t17*t25 + t54*t71); 
complex<T> t73 = t43*t68; 
complex<T> t77 = t43*t55; 
complex<T> t87 = t34*t62; 
complex<T> d22 = spb23*t71; d22 = T(1)/d22;
complex<T> d23 = t34*t71; d23 = T(1)/d23;
complex<T> d25 = spb45*t71; d25 = T(1)/d25;
complex<T> t10 = (-mH2 + s15)*t53*T(4); 
complex<T> t12 = t17*t56 + t38*t87; 
complex<T> t13 = t33*t41 + t17*t56 + t38*t87; 
complex<T> t15 = t27*t81 + t39*t87 + d31*spa45*t97; 
complex<T> t16 = d11*t34*t35*t36 + d12*t18*t61 + d10*t77*T(2) + d9*t67*T(4); 
complex<T> t30 = d23*s234; 
complex<T> t40 = d23*s245; 
complex<T> t42 = d11*t34*t35*t52 + d14*t34*t35*t52 + d13*spa35*t60 + d9*spa35*t60 + d12*t35*t61 + d15*t35*t61 - d16*t73*T(2) - d10*t77*T(2); 
complex<T> t50 = t27*t81 + t39*t87; 
complex<T> t58 = d14*t34*t35*t36 + d15*t18*t61 + d16*t73*T(2) + d13*t67*T(4); 
complex<T> t84 = t29*t33; 
complex<T> t99 = t28*t33; 
complex<T> t100 = d4*t73; 
complex<T> t102 = t32*t33; 
complex<T> t105 = d6*t77; 
complex<T> t11 = t102 + s25*(t17*t25 + t54*t71) + s14*(t53 + t54*t66 + t81*T(4)); 
complex<T> t14 = s45*(t17*t25 + t54*t71) + t84 + s12*(t53 + t54*t66 + t81*T(4)); 
complex<T> t49 = s15*t99 + s24*t53*T(3); 
complex<T> t57 = s235*t30*t33 + d25*spa23*t72; 
complex<T> t63 = s345*t30*t33 + d26*spa34*t72; 
complex<T> t65 = d2*t34*t35*t36 + d3*t18*t61 + t100*T(2) + d1*t67*T(4); 
complex<T> t69 = s235*t33*t40 + d24*spa25*t72; 
complex<T> t70 = d7*t34*t35*t36 + d8*t18*t61 + t105*T(2) + d5*t67*T(4); 
complex<T> t74 = s13*t53 + t84; 
complex<T> t75 = d2*t34*t35*t52 + d7*t34*t35*t52 + d1*spa35*t60 + d5*spa35*t60 + d3*t35*t61 + d8*t35*t61 - t100*T(2) - t105*T(2); 
complex<T> t79 = s345*t33*t40 + d22*spa45*t72; 
complex<T> co1 = s24*t99; 
complex<T> co2 = -t102; 
complex<T> co3 = -(t84*T(2)); 
complex<T> co4 = d27*spa25*spa45*t33; 
complex<T> co5 = d30*spa23*t97; 
SeriesC<T> result = t65*(*CI_users[0]->get_value(mc,ind,mu)) + t75*(*CI_users[1]->get_value(mc,ind,mu)) + t70*(*CI_users[2]->get_value(mc,ind,mu)) + t16*(*CI_users[3]->get_value(mc,ind,mu)) + t42*(*CI_users[4]->get_value(mc,ind,mu)) + t58*(*CI_users[5]->get_value(mc,ind,mu)) + t10*(*CI_users[6]->get_value(mc,ind,mu)) + t59*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t74*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t11*(*CI_users[11]->get_value(mc,ind,mu)) + t49*(*CI_users[12]->get_value(mc,ind,mu)) + t64*(*CI_users[13]->get_value(mc,ind,mu)) + t14*(*CI_users[14]->get_value(mc,ind,mu)) + co3*(*CI_users[15]->get_value(mc,ind,mu)) + t79*(*CI_users[16]->get_value(mc,ind,mu)) + t69*(*CI_users[17]->get_value(mc,ind,mu)) + t57*(*CI_users[18]->get_value(mc,ind,mu)) + t63*(*CI_users[19]->get_value(mc,ind,mu)) + co4*(*CI_users[20]->get_value(mc,ind,mu)) + t13*(*CI_users[21]->get_value(mc,ind,mu)) + t50*(*CI_users[22]->get_value(mc,ind,mu)) + co5*(*CI_users[23]->get_value(mc,ind,mu)) + t12*(*CI_users[24]->get_value(mc,ind,mu)) + t15*(*CI_users[25]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C4g1ph_phdpmpm_nf_wCI::\
C4g1ph_phdpmpm_nf_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C4g1ph_phdpmpm_nf_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, p, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpmpm nf");
#endif
 
//#define TimeStamp "Wed 8 Dec 2010 21:33:07 on n2173"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> s23 = S(2,3);
complex<T> s14 = S(1,4);
complex<T> s25 = S(2,5);
complex<T> s34 = S(3,4);
complex<T> s12 = S(1,2);
complex<T> s45 = S(4,5);
complex<T> s35 = -(spa35*spb35);
complex<T> t6 = cube(spb35); 
complex<T> t14 = square(spb24); 
complex<T> t15 = square(spa35); 
complex<T> t21 = -s23 - s35; 
complex<T> t22 = s34*s45; 
complex<T> t27 = spb23*spb45; 
complex<T> t28 = spb25*spb34; 
complex<T> t29 = cube(spa35); 
complex<T> t36 = square(spb23); 
complex<T> t37 = square(spb25); 
complex<T> t38 = square(spb34); 
complex<T> t39 = square(spb45); 
complex<T> t44 = -(spa35*T(2)); 
complex<T> t54 = spa35*T(2); 
complex<T> d1 = -((s25 + s35)*spb35*T(3)); d1 = T(1)/d1;
complex<T> d3 = square(s25 + s35)*square(spb35)*T(3); d3 = T(1)/d3;
complex<T> d4 = -(spb35*cube(s25 + s35)*T(9)); d4 = T(1)/d4;
complex<T> d9 = -((s35 + s45)*spb35*T(3)); d9 = T(1)/d9;
complex<T> d10 = -(spb35*cube(s35 + s45)*T(9)); d10 = T(1)/d10;
complex<T> d12 = square(s35 + s45)*square(spb35)*T(3); d12 = T(1)/d12;
complex<T> d13 = -((s34 + s35)*spb35*T(3)); d13 = T(1)/d13;
complex<T> d15 = square(s34 + s35)*square(spb35)*T(3); d15 = T(1)/d15;
complex<T> d16 = -(spb35*cube(s34 + s35)*T(9)); d16 = T(1)/d16;
complex<T> d17 = square(spb35)*T(3); d17 = T(1)/d17;
complex<T> d18 = square(square(spb35))*T(3); d18 = T(1)/d18;
complex<T> d19 = square(spb35)*T(12); d19 = T(1)/d19;
complex<T> d20 = square(square(spb35))*T(6); d20 = T(1)/d20;
complex<T> t16 = -t28; 
complex<T> t31 = -(d18*T(2)); 
complex<T> t34 = d19*t22; 
complex<T> t40 = -(t29*T(2)); 
complex<T> t43 = t27*t28; 
complex<T> t45 = d17*t14; 
complex<T> t51 = d18*T(2); 
complex<T> t55 = t37*t38; 
complex<T> t56 = t36*t39; 
complex<T> t58 = spa35*t14; 
complex<T> d2 = -((s25 + s35)*t6*T(3)); d2 = T(1)/d2;
complex<T> d5 = spb35*t21*T(3); d5 = T(1)/d5;
complex<T> d6 = spb35*cube(t21)*T(9); d6 = T(1)/d6;
complex<T> d7 = t21*t6*T(3); d7 = T(1)/d7;
complex<T> d8 = square(spb35)*square(t21)*T(3); d8 = T(1)/d8;
complex<T> d11 = -((s35 + s45)*t6*T(3)); d11 = T(1)/d11;
complex<T> d14 = -((s34 + s35)*t6*T(3)); d14 = T(1)/d14;
complex<T> t12 = s25*t31*t43 - s14*t45 + s25*t45 + s14*t43*t51; 
complex<T> t26 = s45*t31*t43 - s12*t45 + s45*t45 + s12*t43*t51; 
complex<T> t42 = s23*(t31*t43 + t45); 
complex<T> t48 = s34*(t31*t43 + t45); 
complex<T> t49 = -t58; 
complex<T> t50 = t16*t27; 
complex<T> t7 = s23*s25*(d19*t14 + d20*t50); 
complex<T> t13 = d12*t15*t43 + d9*t49 + d11*t43*t54 + d10*t40*t55; 
complex<T> t35 = d8*t15*t43 + d5*t49 + d7*t43*t54 + d6*t40*t55; 
complex<T> t41 = t14*t34 + d20*t22*t50; 
complex<T> t47 = d15*t15*t43 + d13*t49 + d14*t43*t54 + d16*t40*t56; 
complex<T> t53 = d3*t15*t43 + d1*t49 + d2*t43*t54 + d4*t40*t56; 
complex<T> t57 = d11*t43*t44 + d14*t43*t44 + d12*t15*t50 + d15*t15*t50 + d13*t58 + d9*t58 + d10*t29*t55*T(2) + d16*t29*t56*T(2); 
complex<T> t60 = d2*t43*t44 + d7*t43*t44 + d3*t15*t50 + d8*t15*t50 + d1*t58 + d5*t58 + d6*t29*t55*T(2) + d4*t29*t56*T(2); 
SeriesC<T> result = t53*(*CI_users[0]->get_value(mc,ind,mu)) + t60*(*CI_users[1]->get_value(mc,ind,mu)) + t35*(*CI_users[2]->get_value(mc,ind,mu)) + t13*(*CI_users[3]->get_value(mc,ind,mu)) + t57*(*CI_users[4]->get_value(mc,ind,mu)) + t47*(*CI_users[5]->get_value(mc,ind,mu)) + t42*(*CI_users[6]->get_value(mc,ind,mu)) + t12*(*CI_users[7]->get_value(mc,ind,mu)) + t48*(*CI_users[8]->get_value(mc,ind,mu)) + t26*(*CI_users[9]->get_value(mc,ind,mu)) + t7*(*CI_users[10]->get_value(mc,ind,mu)) + t41*(*CI_users[11]->get_value(mc,ind,mu)) + t7*(*CI_users[12]->get_value(mc,ind,mu)) + t41*(*CI_users[13]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C4g1ph_phdpmmm_G_wCI::\
C4g1ph_phdpmmm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdpmmm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpmmm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdpmmm_nf_wCI::\
C4g1ph_phdpmmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdpmmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdpmmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdmpmm_G_wCI::\
C4g1ph_phdmpmm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdmpmm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmpmm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdmpmm_nf_wCI::\
C4g1ph_phdmpmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdmpmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmpmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdmmpm_G_wCI::\
C4g1ph_phdmmpm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdmmpm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmmpm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdmmpm_nf_wCI::\
C4g1ph_phdmmpm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdmmpm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmmpm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdmmmp_G_wCI::\
C4g1ph_phdmmmp_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdmmmp_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmmmp G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdmmmp_nf_wCI::\
C4g1ph_phdmmmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdmmmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmmmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdmmmm_G_wCI::\
C4g1ph_phdmmmm_G_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdmmmm_G_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmmmm G");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C4g1ph_phdmmmm_nf_wCI::\
C4g1ph_phdmmmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C4g1ph_phdmmmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C4g1ph :  phdmmmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
 // *************** table of switch values ************* 
 
#define _C_phpppp_G C4g1ph_1021_G
#define _C_phpppp_nf C4g1ph_1021_nf
#define _C_phmppp_G C4g1ph_1009_G
#define _C_phmppp_nf C4g1ph_1009_nf
#define _C_phpmpp_G C4g1ph_973_G
#define _C_phpmpp_nf C4g1ph_973_nf
#define _C_phppmp_G C4g1ph_829_G
#define _C_phppmp_nf C4g1ph_829_nf
#define _C_phpppm_G C4g1ph_253_G
#define _C_phpppm_nf C4g1ph_253_nf
#define _C_phmmpp_G C4g1ph_961_G
#define _C_phmmpp_nf C4g1ph_961_nf
#define _C_phpmmp_G C4g1ph_781_G
#define _C_phpmmp_nf C4g1ph_781_nf
#define _C_phppmm_G C4g1ph_61_G
#define _C_phppmm_nf C4g1ph_61_nf
#define _C_phmppm_G C4g1ph_241_G
#define _C_phmppm_nf C4g1ph_241_nf
#define _C_phmpmp_G C4g1ph_817_G
#define _C_phmpmp_nf C4g1ph_817_nf
#define _C_phpmpm_G C4g1ph_205_G
#define _C_phpmpm_nf C4g1ph_205_nf
#define _C_phpmmm_G C4g1ph_13_G
#define _C_phpmmm_nf C4g1ph_13_nf
#define _C_phmpmm_G C4g1ph_49_G
#define _C_phmpmm_nf C4g1ph_49_nf
#define _C_phmmpm_G C4g1ph_193_G
#define _C_phmmpm_nf C4g1ph_193_nf
#define _C_phmmmp_G C4g1ph_769_G
#define _C_phmmmp_nf C4g1ph_769_nf
#define _C_phmmmm_G C4g1ph_1_G
#define _C_phmmmm_nf C4g1ph_1_nf
#define _C_phdpppp_G C4g1ph_1022_G
#define _C_phdpppp_nf C4g1ph_1022_nf
#define _C_phdmppp_G C4g1ph_1010_G
#define _C_phdmppp_nf C4g1ph_1010_nf
#define _C_phdpmpp_G C4g1ph_974_G
#define _C_phdpmpp_nf C4g1ph_974_nf
#define _C_phdppmp_G C4g1ph_830_G
#define _C_phdppmp_nf C4g1ph_830_nf
#define _C_phdpppm_G C4g1ph_254_G
#define _C_phdpppm_nf C4g1ph_254_nf
#define _C_phdmmpp_G C4g1ph_962_G
#define _C_phdmmpp_nf C4g1ph_962_nf
#define _C_phdpmmp_G C4g1ph_782_G
#define _C_phdpmmp_nf C4g1ph_782_nf
#define _C_phdppmm_G C4g1ph_62_G
#define _C_phdppmm_nf C4g1ph_62_nf
#define _C_phdmppm_G C4g1ph_242_G
#define _C_phdmppm_nf C4g1ph_242_nf
#define _C_phdmpmp_G C4g1ph_818_G
#define _C_phdmpmp_nf C4g1ph_818_nf
#define _C_phdpmpm_G C4g1ph_206_G
#define _C_phdpmpm_nf C4g1ph_206_nf
#define _C_phdpmmm_G C4g1ph_14_G
#define _C_phdpmmm_nf C4g1ph_14_nf
#define _C_phdmpmm_G C4g1ph_50_G
#define _C_phdmpmm_nf C4g1ph_50_nf
#define _C_phdmmpm_G C4g1ph_194_G
#define _C_phdmmpm_nf C4g1ph_194_nf
#define _C_phdmmmp_G C4g1ph_770_G
#define _C_phdmmmp_nf C4g1ph_770_nf
#define _C_phdmmmm_G C4g1ph_2_G
#define _C_phdmmmm_nf C4g1ph_2_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_phpppp_G case 1021
 
#define _CASE_phpppp_nf case 1021
 
#define _CASE_phmppp_G case 1009
 
#define _CASE_phmppp_nf case 1009
 
#define _CASE_phpmpp_G case 973
 
#define _CASE_phpmpp_nf case 973
 
#define _CASE_phppmp_G case 829
 
#define _CASE_phppmp_nf case 829
 
#define _CASE_phpppm_G case 253
 
#define _CASE_phpppm_nf case 253
 
#define _CASE_phmmpp_G case 961
 
#define _CASE_phmmpp_nf case 961
 
#define _CASE_phpmmp_G case 781
 
#define _CASE_phpmmp_nf case 781
 
#define _CASE_phppmm_G case 61
 
#define _CASE_phppmm_nf case 61
 
#define _CASE_phmppm_G case 241
 
#define _CASE_phmppm_nf case 241
 
#define _CASE_phmpmp_G case 817
 
#define _CASE_phmpmp_nf case 817
 
#define _CASE_phpmpm_G case 205
 
#define _CASE_phpmpm_nf case 205
 
#define _CASE_phpmmm_G case 13
 
#define _CASE_phpmmm_nf case 13
 
#define _CASE_phmpmm_G case 49
 
#define _CASE_phmpmm_nf case 49
 
#define _CASE_phmmpm_G case 193
 
#define _CASE_phmmpm_nf case 193
 
#define _CASE_phmmmp_G case 769
 
#define _CASE_phmmmp_nf case 769
 
#define _CASE_phmmmm_G case 1
 
#define _CASE_phmmmm_nf case 1
 
#define _CASE_phdpppp_G case 1022
 
#define _CASE_phdpppp_nf case 1022
 
#define _CASE_phdmppp_G case 1010
 
#define _CASE_phdmppp_nf case 1010
 
#define _CASE_phdpmpp_G case 974
 
#define _CASE_phdpmpp_nf case 974
 
#define _CASE_phdppmp_G case 830
 
#define _CASE_phdppmp_nf case 830
 
#define _CASE_phdpppm_G case 254
 
#define _CASE_phdpppm_nf case 254
 
#define _CASE_phdmmpp_G case 962
 
#define _CASE_phdmmpp_nf case 962
 
#define _CASE_phdpmmp_G case 782
 
#define _CASE_phdpmmp_nf case 782
 
#define _CASE_phdppmm_G case 62
 
#define _CASE_phdppmm_nf case 62
 
#define _CASE_phdmppm_G case 242
 
#define _CASE_phdmppm_nf case 242
 
#define _CASE_phdmpmp_G case 818
 
#define _CASE_phdmpmp_nf case 818
 
#define _CASE_phdpmpm_G case 206
 
#define _CASE_phdpmpm_nf case 206
 
#define _CASE_phdpmmm_G case 14
 
#define _CASE_phdpmmm_nf case 14
 
#define _CASE_phdmpmm_G case 50
 
#define _CASE_phdmpmm_nf case 50
 
#define _CASE_phdmmpm_G case 194
 
#define _CASE_phdmmpm_nf case 194
 
#define _CASE_phdmmmp_G case 770
 
#define _CASE_phdmmmp_nf case 770
 
#define _CASE_phdmmmm_G case 2
 
#define _CASE_phdmmmm_nf case 2
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_4g1ph_G( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_phpppp_G: return new 
                       C4g1ph_phpppp_G_wCI(ind);
    _CASE_phmppp_G: return new 
                       C4g1ph_phmppp_G_wCI(ind);
    _CASE_phpmpp_G: return new 
                       C4g1ph_phpmpp_G_wCI(ind);
    _CASE_phppmp_G: return new 
                       C4g1ph_phppmp_G_wCI(ind);
    _CASE_phpppm_G: return new 
                       C4g1ph_phpppm_G_wCI(ind);
    _CASE_phmmpp_G: return new 
                       C4g1ph_phmmpp_G_wCI(ind);
    _CASE_phpmmp_G: return new 
                       C4g1ph_phpmmp_G_wCI(ind);
    _CASE_phppmm_G: return new 
                       C4g1ph_phppmm_G_wCI(ind);
    _CASE_phmppm_G: return new 
                       C4g1ph_phmppm_G_wCI(ind);
    _CASE_phmpmp_G: return new 
                       C4g1ph_phmpmp_G_wCI(ind);
    _CASE_phpmpm_G: return new 
                       C4g1ph_phpmpm_G_wCI(ind);
    _CASE_phpmmm_G: return new 
                       C4g1ph_phpmmm_G_wCI(ind);
    _CASE_phmpmm_G: return new 
                       C4g1ph_phmpmm_G_wCI(ind);
    _CASE_phmmpm_G: return new 
                       C4g1ph_phmmpm_G_wCI(ind);
    _CASE_phmmmp_G: return new 
                       C4g1ph_phmmmp_G_wCI(ind);
    _CASE_phmmmm_G: return new 
                       C4g1ph_phmmmm_G_wCI(ind);
    _CASE_phdpppp_G: return new 
                       C4g1ph_phdpppp_G_wCI(ind);
    _CASE_phdmppp_G: return new 
                       C4g1ph_phdmppp_G_wCI(ind);
    _CASE_phdpmpp_G: return new 
                       C4g1ph_phdpmpp_G_wCI(ind);
    _CASE_phdppmp_G: return new 
                       C4g1ph_phdppmp_G_wCI(ind);
    _CASE_phdpppm_G: return new 
                       C4g1ph_phdpppm_G_wCI(ind);
    _CASE_phdmmpp_G: return new 
                       C4g1ph_phdmmpp_G_wCI(ind);
    _CASE_phdpmmp_G: return new 
                       C4g1ph_phdpmmp_G_wCI(ind);
    _CASE_phdppmm_G: return new 
                       C4g1ph_phdppmm_G_wCI(ind);
    _CASE_phdmppm_G: return new 
                       C4g1ph_phdmppm_G_wCI(ind);
    _CASE_phdmpmp_G: return new 
                       C4g1ph_phdmpmp_G_wCI(ind);
    _CASE_phdpmpm_G: return new 
                       C4g1ph_phdpmpm_G_wCI(ind);
    _CASE_phdpmmm_G: return new 
                       C4g1ph_phdpmmm_G_wCI(ind);
    _CASE_phdmpmm_G: return new 
                       C4g1ph_phdmpmm_G_wCI(ind);
    _CASE_phdmmpm_G: return new 
                       C4g1ph_phdmmpm_G_wCI(ind);
    _CASE_phdmmmp_G: return new 
                       C4g1ph_phdmmmp_G_wCI(ind);
    _CASE_phdmmmm_G: return new 
                       C4g1ph_phdmmmm_G_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_4g1ph_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_phpppp_nf: return new 
                       C4g1ph_phpppp_nf_wCI(ind);
    _CASE_phmppp_nf: return new 
                       C4g1ph_phmppp_nf_wCI(ind);
    _CASE_phpmpp_nf: return new 
                       C4g1ph_phpmpp_nf_wCI(ind);
    _CASE_phppmp_nf: return new 
                       C4g1ph_phppmp_nf_wCI(ind);
    _CASE_phpppm_nf: return new 
                       C4g1ph_phpppm_nf_wCI(ind);
    _CASE_phmmpp_nf: return new 
                       C4g1ph_phmmpp_nf_wCI(ind);
    _CASE_phpmmp_nf: return new 
                       C4g1ph_phpmmp_nf_wCI(ind);
    _CASE_phppmm_nf: return new 
                       C4g1ph_phppmm_nf_wCI(ind);
    _CASE_phmppm_nf: return new 
                       C4g1ph_phmppm_nf_wCI(ind);
    _CASE_phmpmp_nf: return new 
                       C4g1ph_phmpmp_nf_wCI(ind);
    _CASE_phpmpm_nf: return new 
                       C4g1ph_phpmpm_nf_wCI(ind);
    _CASE_phpmmm_nf: return new 
                       C4g1ph_phpmmm_nf_wCI(ind);
    _CASE_phmpmm_nf: return new 
                       C4g1ph_phmpmm_nf_wCI(ind);
    _CASE_phmmpm_nf: return new 
                       C4g1ph_phmmpm_nf_wCI(ind);
    _CASE_phmmmp_nf: return new 
                       C4g1ph_phmmmp_nf_wCI(ind);
    _CASE_phmmmm_nf: return new 
                       C4g1ph_phmmmm_nf_wCI(ind);
    _CASE_phdpppp_nf: return new 
                       C4g1ph_phdpppp_nf_wCI(ind);
    _CASE_phdmppp_nf: return new 
                       C4g1ph_phdmppp_nf_wCI(ind);
    _CASE_phdpmpp_nf: return new 
                       C4g1ph_phdpmpp_nf_wCI(ind);
    _CASE_phdppmp_nf: return new 
                       C4g1ph_phdppmp_nf_wCI(ind);
    _CASE_phdpppm_nf: return new 
                       C4g1ph_phdpppm_nf_wCI(ind);
    _CASE_phdmmpp_nf: return new 
                       C4g1ph_phdmmpp_nf_wCI(ind);
    _CASE_phdpmmp_nf: return new 
                       C4g1ph_phdpmmp_nf_wCI(ind);
    _CASE_phdppmm_nf: return new 
                       C4g1ph_phdppmm_nf_wCI(ind);
    _CASE_phdmppm_nf: return new 
                       C4g1ph_phdmppm_nf_wCI(ind);
    _CASE_phdmpmp_nf: return new 
                       C4g1ph_phdmpmp_nf_wCI(ind);
    _CASE_phdpmpm_nf: return new 
                       C4g1ph_phdpmpm_nf_wCI(ind);
    _CASE_phdpmmm_nf: return new 
                       C4g1ph_phdpmmm_nf_wCI(ind);
    _CASE_phdmpmm_nf: return new 
                       C4g1ph_phdmpmm_nf_wCI(ind);
    _CASE_phdmmpm_nf: return new 
                       C4g1ph_phdmmpm_nf_wCI(ind);
    _CASE_phdmmmp_nf: return new 
                       C4g1ph_phdmmmp_nf_wCI(ind);
    _CASE_phdmmmm_nf: return new 
                       C4g1ph_phdmmmm_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
