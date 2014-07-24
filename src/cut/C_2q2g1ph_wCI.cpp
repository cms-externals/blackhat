/*
*C_2q2g1ph_wCI.cpp
*
* Created on 12/11, 2010
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
 
#define mH2 mc.s(ind[2-1],ind[3-1],ind[4-1],ind[5-1])
#define SPA(i,j) mc.spa(ind[i-1],ind[j-1])
#define SPB(i,j) mc.spb(ind[i-1],ind[j-1])
#define S(i,j) mc.s(ind[i-1],ind[j-1])
#define SS(i,j,k) mc.s(ind[i-1],ind[j-1],ind[k-1])

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}
 

template<class T> complex<T> coeff3mass1234(momentum_configuration<T>& mc,
							const Index_Vector& ind, int i1, int i2, int i3, int i4){

  complex<T> K1sq,K2sq,K1dK2,gammaP,gammaM,resultP,resultM,spab1243,spab1234;
  spab1243 = SPA(i1,i2)*SPB(i2,i3) + SPA(i1,i4)*SPB(i4,i3);
  spab1234 = SPA(i1,i2)*SPB(i2,i4) + SPA(i1,i3)*SPB(i3,i4);
  K1sq = mH2;
  K2sq = S(i1,i2);
  K1dK2 = S(i1,i2) + (S(i1,i3)+S(i1,i4)+S(i2,i3)+S(i2,i4))/T(2);
  gammaP = K1dK2 + sqrt(pow(K1dK2,2)-K1sq*K2sq);
  gammaM = K1dK2 - sqrt(pow(K1dK2,2)-K1sq*K2sq);
  resultP = (
	     pow(mH2,2)*pow(SPA(i3,i4),3)*(mH2*SPA(i1,i2)*SPB(i2,i3)-gammaP*spab1243)*(mH2*SPA(i1,i2)*SPB(i2,i4)-gammaP*spab1234)/T(4)/gammaP/(gammaP-mH2)/SPA(i1,i2)
	     /(mH2*(S(i1,i3)+S(i2,i3))/T(2)+gammaP*(SS(i1,i2,i4)-mH2)/T(2))
	     /(mH2*(S(i1,i4)+S(i2,i4))/T(2)+gammaP*(SS(i1,i2,i3)-mH2)/T(2))
	     );

  resultM = (
	     pow(mH2,2)*pow(SPA(i3,i4),3)*(mH2*SPA(i1,i2)*SPB(i2,i3)-gammaM*spab1243)*(mH2*SPA(i1,i2)*SPB(i2,i4)-gammaM*spab1234)/T(4)/gammaM/(gammaM-mH2)/SPA(i1,i2)
	     /(mH2*(S(i1,i3)+S(i2,i3))/T(2)+gammaM*(SS(i1,i2,i4)-mH2)/T(2))
	     /(mH2*(S(i1,i4)+S(i2,i4))/T(2)+gammaM*(SS(i1,i2,i3)-mH2)/T(2))
	     );

  return (resultP+resultM);
}


  template complex<R> coeff3mass1234(momentum_configuration<R>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);
  template complex<RHP> coeff3mass1234(momentum_configuration<RHP>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);
  template complex<RVHP> coeff3mass1234(momentum_configuration<RVHP>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);

template<class T> complex<T> coeff3mass4123(momentum_configuration<T>& mc,
							const Index_Vector& ind, int i1, int i2, int i3, int i4){

  complex<T> K1sq,K2sq,K1dK2,gammaP,gammaM,resultP,resultM,spab3241,spab3142;
  spab3241 = SPA(i3,i2)*SPB(i2,i1) + SPA(i3,i4)*SPB(i4,i1);
  spab3142 = SPA(i3,i1)*SPB(i1,i2) + SPA(i3,i4)*SPB(i4,i2);
  K1sq = mH2;
  K2sq = S(i1,i4);
  K1dK2 = S(i1,i4) + (S(i1,i3)+S(i1,i2)+S(i3,i4)+S(i2,i4))/T(2);
  gammaP = K1dK2 + sqrt(pow(K1dK2,2)-K1sq*K2sq);
  gammaM = K1dK2 - sqrt(pow(K1dK2,2)-K1sq*K2sq);
  resultP = (
	     -pow(mH2,2)*pow(SPA(i1,i4),2)*(mH2*SPA(i3,i4)*SPB(i4,i1)-gammaP*spab3241)*(mH2-gammaP)*spab3142/T(8)/gammaP/(gammaP-mH2)
	     /(mH2*S(i1,i4)/T(2)+gammaP*(SS(i2,i3,i4)-mH2)/T(2))
	     /(mH2*(S(i1,i2)+S(i2,i4))/T(2)+gammaP*(SS(i1,i3,i4)-mH2)/T(2))
	     );

  resultM = (
	     -pow(mH2,2)*pow(SPA(i1,i4),2)*(mH2*SPA(i3,i4)*SPB(i4,i1)-gammaM*spab3241)*(mH2-gammaM)*spab3142/T(8)/gammaM/(gammaM-mH2)
	     /(mH2*S(i1,i4)/T(2)+gammaM*(SS(i2,i3,i4)-mH2)/T(2))
	     /(mH2*(S(i1,i2)+S(i2,i4))/T(2)+gammaM*(SS(i1,i3,i4)-mH2)/T(2))
	     );


  return (resultP+resultM);
}


  template complex<R> coeff3mass4123(momentum_configuration<R>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);
  template complex<RHP> coeff3mass4123(momentum_configuration<RHP>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);
  template complex<RVHP> coeff3mass4123(momentum_configuration<RVHP>& mc,const Index_Vector& ind,int i1, int i2, int i3, int i4);


 
class C2q2g1ph_phqmqppm_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqppm_LT_wCI
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

 
class C2q2g1ph_phqmqppm_RT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqppm_RT_wCI
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

 
class C2q2g1ph_phqmqppm_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqppm_nfLT_wCI
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

 
class C2q2g1ph_phqmqpmp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqpmp_LT_wCI
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

 
class C2q2g1ph_phqmqpmp_RT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqpmp_RT_wCI
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

 
class C2q2g1ph_phqmqpmp_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqpmp_nfLT_wCI
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

 
class C2q2g1ph_phqmqpmm_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqpmm_LT_wCI
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

 
class C2q2g1ph_phqmqpmm_RT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqpmm_RT_wCI
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

 
class C2q2g1ph_phqmqpmm_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqpmm_nfLT_wCI
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

 
class C2q2g1ph_phqmqppp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqppp_LT_wCI
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

 
class C2q2g1ph_phqmqppp_RT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqppp_RT_wCI
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

 
class C2q2g1ph_phqmqppp_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmqppp_nfLT_wCI
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

 
class C2q2g1ph_phqmpqpm_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmpqpm_LT_wCI
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

 
class C2q2g1ph_phqmpqpm_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmpqpm_nfLT_wCI
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

 
class C2q2g1ph_phqmmqpp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmmqpp_LT_wCI
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

 
class C2q2g1ph_phqmmqpp_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phqmmqpp_nfLT_wCI
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

 
class C2q2g1ph_phdqmqppm_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqppm_LT_wCI
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

 
class C2q2g1ph_phdqmqppm_RT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqppm_RT_wCI
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

 
class C2q2g1ph_phdqmqppm_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqppm_nfLT_wCI
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

 
class C2q2g1ph_phdqmqpmp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqpmp_LT_wCI
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

 
class C2q2g1ph_phdqmqpmp_RT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqpmp_RT_wCI
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

 
class C2q2g1ph_phdqmqpmp_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqpmp_nfLT_wCI
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

 
class C2q2g1ph_phdqmqpmm_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqpmm_LT_wCI
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

 
class C2q2g1ph_phdqmqpmm_RT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqpmm_RT_wCI
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

 
class C2q2g1ph_phdqmqpmm_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqpmm_nfLT_wCI
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

 
class C2q2g1ph_phdqmqppp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqppp_LT_wCI
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

 
class C2q2g1ph_phdqmqppp_RT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqppp_RT_wCI
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

 
class C2q2g1ph_phdqmqppp_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmqppp_nfLT_wCI
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

 
class C2q2g1ph_phdqmpqpm_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmpqpm_LT_wCI
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

 
class C2q2g1ph_phdqmpqpm_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmpqpm_nfLT_wCI
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

 
class C2q2g1ph_phdqmmqpp_LT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmmqpp_LT_wCI
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

 
class C2q2g1ph_phdqmmqpp_nfLT_wCI : public Cut_Part_wCI {
public:
       C2q2g1ph_phdqmmqpp_nfLT_wCI
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


C2q2g1ph_phqmqppm_LT_wCI::\
C2q2g1ph_phqmqppm_LT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqppm_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, p, m}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqppm LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:49:46 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> s15 = S(1,5);
complex<T> s14 = S(1,4);
complex<T> s25 = S(2,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s24 = -(spa24*spb24);
complex<T> s13 = S(1,3);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> t3 = square(spa25*spb24 + spa35*spb34); 
complex<T> t4 = square(s23 - s234); 
complex<T> t6 = square(spa24); 
complex<T> t12 = square(spa25); 
complex<T> t13 = square(spb34); 
complex<T> t14 = -(spa23*spa45); 
complex<T> t25 = cube(spa25); 
complex<T> t27 = square(spa23); 
complex<T> t28 = square(spa45); 
complex<T> t34 = cube(spb34); 
complex<T> t36 = spa35*spa45; 
complex<T> t39 = -(spb34*T(2)); 
complex<T> t48 = spb34*T(2); 
complex<T> t54 = spa25*spa35; 
complex<T> d1 = spa23*spa24*spa45*T(6); d1 = T(1)/d1;
complex<T> d2 = spa23*spa34*spa45*T(2); d2 = T(1)/d2;
complex<T> d3 = (s23 - s234)*spa23*spa45; d3 = T(1)/d3;
complex<T> d4 = (s23 - s234)*spa34; d4 = T(1)/d4;
complex<T> d6 = cube(s23 - s234); d6 = T(1)/d6;
complex<T> d7 = spa24*cube(s23 - s234)*T(3); d7 = T(1)/d7;
complex<T> d10 = spa23*spa45*cube(s23 - s234)*T(3); d10 = T(1)/d10;
complex<T> d11 = spa23*spa34*spa45*T(3); d11 = T(1)/d11;
complex<T> d12 = (-s345 + s45)*spa34; d12 = T(1)/d12;
complex<T> d13 = spa34*square(s345 - s45)*T(2); d13 = T(1)/d13;
complex<T> d14 = spa23*spa34*spa45; d14 = T(1)/d14;
complex<T> d15 = spa24*spa45; d15 = T(1)/d15;
complex<T> d16 = spa23*spa45; d16 = T(1)/d16;
complex<T> d17 = spa23*spa24*spa45; d17 = T(1)/d17;
complex<T> d18 = spa23*spa24; d18 = T(1)/d18;
complex<T> d19 = spa23*spa34; d19 = T(1)/d19;
complex<T> d20 = spa34*spa45*T(2); d20 = T(1)/d20;
complex<T> d21 = spa23*spa45*T(2); d21 = T(1)/d21;
complex<T> d22 = spa23*spa24*T(2); d22 = T(1)/d22;
complex<T> d23 = spa24*spa45*T(2); d23 = T(1)/d23;
complex<T> d24 = spa23*T(2); d24 = T(1)/d24;
complex<T> t16 = -(d14*spa35); 
complex<T> t20 = d2*s234; 
complex<T> t21 = d20*spb23; 
complex<T> t23 = d19*spb45; 
complex<T> t31 = d16*spb24; 
complex<T> t45 = d10*t6; 
complex<T> t47 = spa35*t12; 
complex<T> t55 = d13*t13; 
complex<T> t61 = spa25*t13; 
complex<T> t73 = d17*t25; 
complex<T> t77 = -(spb45*t25); 
complex<T> d5 = spa24*t4*T(2); d5 = T(1)/d5;
complex<T> d8 = spa23*spa45*t4; d8 = T(1)/d8;
complex<T> d9 = spa23*spa34*spa45*t4*T(2); d9 = T(1)/d9;
complex<T> t10 = d14*(-mH2 + s15)*t47*T(3); 
complex<T> t22 = -t31; 
complex<T> t38 = -t47; 
complex<T> t46 = s25*t12*t16 + d14*s14*t47; 
complex<T> t52 = s24*t12*t16 + t25*t31; 
complex<T> t53 = d12*t39*t54 + spa35*t14*t55; 
complex<T> t60 = d12*t48*t54 + spa23*t36*t55; 
complex<T> t65 = -t73; 
complex<T> t68 = spb34*t47; 
complex<T> t69 = (s235*t20 + mH2*t21)*t47; 
complex<T> t1 = -(d3*spb24*t25) + d8*spa24*spa25*t3 - d7*t27*t28*t34 + d6*spa23*t34*t36 + d2*t38 + d3*spb34*t38 + d4*t39*t54 + d9*spa35*t3*t6 + d5*spa23*spa45*t61 + d6*spa23*spa45*spb24*t61 - t45*cube(spa25*spb24 + spa35*spb34)*T(4) - d1*t25*T(5); 
complex<T> t2 = d3*spb24*t25 - d8*spa24*spa25*t3 + d6*spa35*t14*t34 + d7*t27*t28*t34 + d4*t48*t54 - d9*spa35*t3*t6 + d5*t14*t61 + d6*spb24*t14*t61 + d3*t68 + t45*cube(spa25*spb24 + spa35*spb34)*T(4) + d1*t25*T(5) - d11*t47*T(5); 
complex<T> t37 = d16*t68 + s34*t73; 
complex<T> t59 = t22*t25 + d14*s24*t47 + s15*(t12*t16 + t65); 
complex<T> t64 = d2*(mH2*s25*t38 + s235*s245*t47); 
complex<T> t67 = s345*t20*t47 + d21*mH2*t68; 
complex<T> t72 = d14*s13*t47 + t23*t47 + s13*t65 + d18*t77; 
complex<T> co1 = -(d15*spb23*t25); 
complex<T> co2 = s25*t73; 
complex<T> co3 = t23*t38; 
complex<T> co4 = d22*s25*t77; 
complex<T> co5 = s25*t21*t38; 
complex<T> co6 = -(d23*s34*spb23*t25); 
complex<T> co7 = d24*spb45*t68; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t60*(*CI_users[2]->get_value(mc,ind,mu)) + t53*(*CI_users[3]->get_value(mc,ind,mu)) + t10*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t52*(*CI_users[6]->get_value(mc,ind,mu)) + t72*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu)) + t46*(*CI_users[9]->get_value(mc,ind,mu)) + t59*(*CI_users[10]->get_value(mc,ind,mu)) + t37*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + t64*(*CI_users[13]->get_value(mc,ind,mu)) + t69*(*CI_users[14]->get_value(mc,ind,mu)) + t67*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu)) + co6*(*CI_users[18]->get_value(mc,ind,mu)) + co7*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmqppm_RT_wCI::\
C2q2g1ph_phqmqppm_RT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqppm_RT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, p, m}, RT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqppm RT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:49:50 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> s25 = S(2,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s234 = SS(2,3,4);
complex<T> s345 = SS(3,4,5);
complex<T> s245 = SS(2,4,5);
complex<T> t2 = spa23*spa45; 
complex<T> t3 = square(spa24*spb34 + spa25*spb35); 
complex<T> t7 = cube(spa35); 
complex<T> t11 = square(spa25); 
complex<T> t15 = square(spa23); 
complex<T> t16 = square(spa45); 
complex<T> t17 = square(spb34); 
complex<T> t23 = -(spb34*T(2)); 
complex<T> t31 = spb34*T(2); 
complex<T> d1 = spa24*spa34*T(2); d1 = T(1)/d1;
complex<T> d3 = (s23 - s234)*spa24*spa34; d3 = T(1)/d3;
complex<T> d4 = spa24*spa34*square(s23 - s234)*T(2); d4 = T(1)/d4;
complex<T> d5 = (-s345 + s45)*spa34; d5 = T(1)/d5;
complex<T> d8 = spa24*spa34; d8 = T(1)/d8;
complex<T> d9 = spa34*spa45; d9 = T(1)/d9;
complex<T> d10 = spa34; d10 = T(1)/d10;
complex<T> d11 = spa23*spa34; d11 = T(1)/d11;
complex<T> d12 = spa24; d12 = T(1)/d12;
complex<T> d13 = spa23*spa34*T(2); d13 = T(1)/d13;
complex<T> d14 = spa24*T(2); d14 = T(1)/d14;
complex<T> t10 = d4*t15*t16*t17 + d3*spa25*t2*t31 + d1*t11*T(3); 
complex<T> t18 = d11*spb45; 
complex<T> t25 = d10*spb24; 
complex<T> t30 = spa35*t11; 
complex<T> t39 = d8*t11; 
complex<T> d2 = spa34*t2*T(2); d2 = T(1)/d2;
complex<T> d6 = spa34*t2*square(s345 - s45)*T(2); d6 = T(1)/d6;
complex<T> d7 = spa34*t2; d7 = T(1)/d7;
complex<T> t1 = d2*t30 + d5*spa25*spa35*t31 - d6*t3*t7; 
complex<T> t9 = -(d4*t15*t16*t17) + d3*spa25*t2*t23 - d1*t11*T(3) + d2*t30*T(3); 
complex<T> t12 = -(d7*spa35); 
complex<T> t20 = -t25; 
complex<T> t22 = -t30; 
complex<T> t28 = d9*spb23*t30 + s23*t39; 
complex<T> t35 = -t39; 
complex<T> t38 = (d2*s245*s345 + d13*mH2*spb45)*t30; 
complex<T> t42 = d7*t30; 
complex<T> t29 = d2*t22 + d5*spa25*spa35*t23 + d6*t3*t7; 
complex<T> t41 = t18*t30 + s12*t42; 
complex<T> t43 = t18*t30 + s45*t39 + s13*(t35 + t42); 
complex<T> t44 = t11*t12; 
complex<T> t8 = s15*t42 + mH2*t44; 
complex<T> t34 = t11*t25 + s24*t44; 
complex<T> t45 = t11*t20 + s24*t42 + s15*(t35 + t44); 
complex<T> co1 = s25*t39; 
complex<T> co2 = -(d12*spb34*t11); 
complex<T> co3 = t18*t22; 
complex<T> co4 = d1*s25*s45*t11; 
complex<T> co5 = -(d14*s23*spb34*t11); 
SeriesC<T> result = t9*(*CI_users[0]->get_value(mc,ind,mu)) + t10*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t29*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + t28*(*CI_users[5]->get_value(mc,ind,mu)) + t34*(*CI_users[6]->get_value(mc,ind,mu)) + t43*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t45*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t41*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + t38*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[14]->get_value(mc,ind,mu)) + co5*(*CI_users[15]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmqppm_nfLT_wCI::\
C2q2g1ph_phqmqppm_nfLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqppm_nfLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, p, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqppm nfLT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:49:50 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> s23 = S(2,3);
complex<T> s234 = SS(2,3,4);
complex<T> t4 = cube(spa25*spb24 + spa35*spb34); 
complex<T> t5 = spa23*spa45; 
complex<T> t7 = square(spa25); 
complex<T> t8 = square(spa23); 
complex<T> t9 = square(spa45); 
complex<T> t10 = cube(spb34); 
complex<T> d1 = spa24*spa34*T(3); d1 = T(1)/d1;
complex<T> d3 = spa24*cube(s23 - s234)*T(3); d3 = T(1)/d3;
complex<T> d2 = spa34*t5*T(3); d2 = T(1)/d2;
complex<T> d4 = t5*cube(s23 - s234)*T(3); d4 = T(1)/d4;
complex<T> t12 = d2*spa35; 
complex<T> t1 = d1*t7 + t12*t7 - d3*t10*t8*t9 - d4*t4*square(spa24); 
complex<T> t2 = -(d1*t7) + t12*t7 + d3*t10*t8*t9 + d4*t4*square(spa24); 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmqpmp_LT_wCI::\
C2q2g1ph_phqmqpmp_LT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqpmp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, m, p}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqpmp LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:50:20 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s34 = S(3,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s14 = S(1,4);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s235 = SS(2,3,5);
complex<T> s234 = SS(2,3,4);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> t5 = -(spa23*spb35); 
complex<T> t10 = square(spa24*spb23 + spa45*spb35); 
complex<T> t11 = square(spa23*spb35 + spa24*spb45); 
complex<T> t14 = spa45*T(3); 
complex<T> t22 = square(spa24); 
complex<T> t23 = square(spa34); 
complex<T> t24 = square(spb35); 
complex<T> t26 = square(spa25); 
complex<T> t27 = s23 - s235; 
complex<T> t41 = cube(spa24); 
complex<T> t42 = cube(spa34); 
complex<T> t51 = square(spa23); 
complex<T> t52 = square(spa45); 
complex<T> t53 = cube(spb35); 
complex<T> t60 = -(spa45*spb35); 
complex<T> t68 = spa34*spb35; 
complex<T> d2 = spa23*spa35*spa45*T(2); d2 = T(1)/d2;
complex<T> d5 = spb25*square(spa35); d5 = T(1)/d5;
complex<T> d6 = spa23*spa45*spb25*square(spa35); d6 = T(1)/d6;
complex<T> d13 = spa23*spa25*spa45*T(6); d13 = T(1)/d13;
complex<T> d14 = (-s235 + s25)*spa35; d14 = T(1)/d14;
complex<T> d15 = spa25*spb23*square(spa35); d15 = T(1)/d15;
complex<T> d16 = spa35*square(s235 - s25)*T(2); d16 = T(1)/d16;
complex<T> d17 = spa25*spb23*square(s235 - s25); d17 = T(1)/d17;
complex<T> d18 = (s34 - s345)*spa35; d18 = T(1)/d18;
complex<T> d19 = spa35*square(s34 - s345)*T(2); d19 = T(1)/d19;
complex<T> d20 = spa25*spb45*square(spa35); d20 = T(1)/d20;
complex<T> d21 = spa25*spb45*square(s34 - s345); d21 = T(1)/d21;
complex<T> d22 = (-s345 + s45)*square(spa35); d22 = T(1)/d22;
complex<T> d23 = spa35*square(s345 - s45)*T(2); d23 = T(1)/d23;
complex<T> d24 = spa23*spa25*spa45; d24 = T(1)/d24;
complex<T> d25 = spa45*cube(spa35); d25 = T(1)/d25;
complex<T> d26 = spa23*spa45; d26 = T(1)/d26;
complex<T> d27 = spa23*spa45*cube(spa35); d27 = T(1)/d27;
complex<T> d28 = spa23*spa25; d28 = T(1)/d28;
complex<T> d29 = spa23*cube(spa35); d29 = T(1)/d29;
complex<T> d30 = spa23*spa25*spa45*T(2); d30 = T(1)/d30;
complex<T> d31 = spa23*spa45*T(2); d31 = T(1)/d31;
complex<T> d32 = spa25*spa45*T(2); d32 = T(1)/d32;
complex<T> d33 = spa23*T(2); d33 = T(1)/d33;
complex<T> d34 = spa45*cube(spa35)*T(2); d34 = T(1)/d34;
complex<T> d35 = spa23*cube(spa35)*T(2); d35 = T(1)/d35;
complex<T> t2 = -(d16*spa25*t23*t24) + d17*t10*t5 + d15*spa23*spb35*t52 + d14*spa24*t68*T(3); 
complex<T> t4 = square(spa24*spb25 + t68); 
complex<T> t29 = -(d24*T(2)); 
complex<T> t34 = d30*s234; 
complex<T> t35 = d32*spb23; 
complex<T> t36 = d26*spb25; 
complex<T> t39 = d28*spb45; 
complex<T> t64 = d21*t11; 
complex<T> t69 = d24*t41; 
complex<T> t72 = d22*spa45; 
complex<T> t78 = spb25*t41; 
complex<T> t88 = spb35*t23; 
complex<T> t97 = t51*t52; 
complex<T> t98 = t26*t42; 
complex<T> t102 = d6*spb35; 
complex<T> d1 = spa23*spa25*t14; d1 = T(1)/d1;
complex<T> d3 = spa23*spa45*t27; d3 = T(1)/d3;
complex<T> d4 = spa25*t27*T(3); d4 = T(1)/d4;
complex<T> d7 = spa25*square(t27)*T(3); d7 = T(1)/d7;
complex<T> d8 = spa25*cube(t27)*T(3); d8 = T(1)/d8;
complex<T> d9 = spa35*square(t27)*T(2); d9 = T(1)/d9;
complex<T> d10 = spa23*t14*square(t27); d10 = T(1)/d10;
complex<T> d11 = spa35*spb25*t27; d11 = T(1)/d11;
complex<T> d12 = spa23*spa45*spb25*t27; d12 = T(1)/d12;
complex<T> t1 = -(d2*spa34*t22) + d4*spb35*t22 - d11*t4 + d9*spa25*t4 - d10*spa24*spa25*t4 - d1*t41 - d5*t88 - d7*spa23*spa24*spa45*t24*T(2) + d12*spa24*t4*T(2) - spa25*t102*t42*T(2) + d8*t53*t97*T(2) - d3*t22*t68*T(3) - d3*t78*T(3); 
complex<T> t3 = d17*spa23*spb35*t10 + d2*spa34*t22 - d4*spb35*t22 + d16*spa25*t23*t24 + d11*t4 - d9*spa25*t4 + d10*spa24*spa25*t4 + d15*t5*t52 + d5*t88 + d7*spa23*spa24*spa45*t24*T(2) - d12*spa24*t4*T(2) + spa25*t102*t42*T(2) - d8*t53*t97*T(2) - d14*spa24*t68*T(3) + d3*t22*t68*T(3) + d3*t78*T(3) - d13*t41*T(11); 
complex<T> t18 = s34*(-t69 + d27*t98); 
complex<T> t20 = -(d23*spa23*spa34*spa45*t24) + spa23*t68*t72 + d22*spa25*t88*T(2); 
complex<T> t48 = -t64; 
complex<T> t57 = d30*s235*s245*t41 + d31*mH2*t78; 
complex<T> t58 = -t98; 
complex<T> t66 = (-mH2 + s15)*t69*T(3); 
complex<T> t67 = d23*spa23*spa34*spa45*t24 + d19*spa25*t23*t24 + d20*t51*t60 + spa45*spb35*t64 + spa34*t5*t72 - d22*spa25*t88*T(2) - d18*spa24*t68*T(3); 
complex<T> t75 = -(d30*mH2*s34*t41) + s345*t34*t41; 
complex<T> t109 = t29*t41; 
complex<T> t114 = t35*t41; 
complex<T> t116 = t36*t41; 
complex<T> t118 = t39*t41; 
complex<T> t17 = s15*t109 + s24*t69*T(2); 
complex<T> t19 = d27*s14*t58 + d27*s25*t98 + t116*T(2) + s14*t69*T(2); 
complex<T> t50 = -(d19*spa25*t23*t24) + spa45*spb35*(t48 + d20*t51) + d18*spa24*t68*T(3); 
complex<T> t82 = mH2*t114 + s235*t34*t41; 
complex<T> t110 = spb45*t58; 
complex<T> t21 = d29*t110 + t118 + s12*(d27*t58 + t69); 
complex<T> co1 = d25*spb23*t58; 
complex<T> co2 = s24*t109; 
complex<T> co3 = -t116; 
complex<T> co4 = -t118; 
complex<T> co5 = d33*spb45*t78; 
complex<T> co6 = d34*s25*spb23*t58; 
complex<T> co7 = -(s34*t114); 
complex<T> co8 = d35*s34*t110; 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t2*(*CI_users[2]->get_value(mc,ind,mu)) + t50*(*CI_users[3]->get_value(mc,ind,mu)) + t67*(*CI_users[4]->get_value(mc,ind,mu)) + t20*(*CI_users[5]->get_value(mc,ind,mu)) + t66*(*CI_users[6]->get_value(mc,ind,mu)) + co1*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu)) + co3*(*CI_users[9]->get_value(mc,ind,mu)) + t19*(*CI_users[10]->get_value(mc,ind,mu)) + t17*(*CI_users[11]->get_value(mc,ind,mu)) + t18*(*CI_users[12]->get_value(mc,ind,mu)) + t21*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[14]->get_value(mc,ind,mu)) + t57*(*CI_users[15]->get_value(mc,ind,mu)) + t82*(*CI_users[16]->get_value(mc,ind,mu)) + t75*(*CI_users[17]->get_value(mc,ind,mu)) + co5*(*CI_users[18]->get_value(mc,ind,mu)) + co6*(*CI_users[19]->get_value(mc,ind,mu)) + co7*(*CI_users[20]->get_value(mc,ind,mu)) + co8*(*CI_users[21]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmqpmp_RT_wCI::\
C2q2g1ph_phqmqpmp_RT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqpmp_RT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, m, p}, RT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqpmp RT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:50:27 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s45 = -(spa45*spb45);
complex<T> s23 = -(spa23*spb23);
complex<T> s14 = S(1,4);
complex<T> s13 = S(1,3);
complex<T> s12 = S(1,2);
complex<T> s235 = SS(2,3,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s345 = SS(3,4,5);
complex<T> s245 = SS(2,4,5);
complex<T> t4 = spa23*spa45; 
complex<T> t7 = square(spa24*spb34 + spa25*spb35); 
complex<T> t8 = square(spa23*spb35 + spa24*spb45); 
complex<T> t9 = spa25*T(2); 
complex<T> t21 = square(spa23); 
complex<T> t22 = square(spa45); 
complex<T> t23 = square(spb35); 
complex<T> t24 = square(spa34); 
complex<T> t34 = cube(spa24); 
complex<T> d2 = (s23 - s235)*spa25*square(spa35); d2 = T(1)/d2;
complex<T> d4 = (-s235 + s25)*square(spa35); d4 = T(1)/d4;
complex<T> d5 = spa35*square(s235 - s25)*T(2); d5 = T(1)/d5;
complex<T> d6 = spa35*square(s34 - s345)*T(2); d6 = T(1)/d6;
complex<T> d7 = spa25*spb45*square(spa35); d7 = T(1)/d7;
complex<T> d8 = spa25*spb45*square(s34 - s345); d8 = T(1)/d8;
complex<T> d10 = spb34*square(spa35); d10 = T(1)/d10;
complex<T> d11 = spb34*square(s345 - s45); d11 = T(1)/d11;
complex<T> d14 = spa25*cube(spa35); d14 = T(1)/d14;
complex<T> d15 = spa25*spa45; d15 = T(1)/d15;
complex<T> d16 = spa23*spa25; d16 = T(1)/d16;
complex<T> d17 = cube(spa35); d17 = T(1)/d17;
complex<T> d19 = cube(spa35)*T(2); d19 = T(1)/d19;
complex<T> t3 = -(d7*spa45*spb35*t21) + d6*spa25*t23*t24 + d8*spa45*spb35*t8; 
complex<T> t27 = d16*spb45; 
complex<T> t37 = -(spa25*t24); 
complex<T> t44 = t21*t22; 
complex<T> t47 = d4*spb35; 
complex<T> d1 = t4*t9; d1 = T(1)/d1;
complex<T> d3 = spa35*t9*square(s23 - s235); d3 = T(1)/d3;
complex<T> d9 = spa35*t4*T(2); d9 = T(1)/d9;
complex<T> d12 = spa35*t4*square(s345 - s45)*T(2); d12 = T(1)/d12;
complex<T> d13 = spa25*t4; d13 = T(1)/d13;
complex<T> d18 = spa23*t9; d18 = T(1)/d18;
complex<T> d20 = t9*cube(spa35); d20 = T(1)/d20;
complex<T> t1 = d10*spb35*t4 - d12*t7*cube(spa34) - d11*t4*cube(spb35) + d9*spa34*square(spa24); 
complex<T> t2 = d7*spa45*spb35*t21 + d6*t23*t37 - d10*spb35*t4 - d8*spa45*spb35*t8 + d12*t7*cube(spa34) + d11*t4*cube(spb35) - d9*spa34*square(spa24); 
complex<T> t19 = d5*spa25*t23*t24 + d2*spb35*t44 + d3*t23*t44 + t37*t47 - spa34*t4*t47*T(2); 
complex<T> t28 = -(d13*s15); 
complex<T> t35 = -t44; 
complex<T> t42 = (d1*s245*s345 + d18*mH2*spb45)*t34; 
complex<T> t43 = d5*t23*t37 + spa25*t24*t47 + spa34*t4*t47*T(2); 
complex<T> t52 = d13*t34; 
complex<T> t59 = t27*t34; 
complex<T> t61 = d14*t44; 
complex<T> t15 = d15*spb23*t34 + s23*t61; 
complex<T> t17 = d2*spb35*t35 + d3*t23*t35 + d1*t34*T(3); 
complex<T> t18 = d14*s12*t35 + s12*t52 + t59 + s45*t61; 
complex<T> t45 = -t52; 
complex<T> t50 = s13*t52 + t59; 
complex<T> t56 = t28*t34 + s24*t52; 
complex<T> t70 = spb25*t35; 
complex<T> t14 = mH2*t45 + s15*t52; 
complex<T> t16 = d14*s14*t35 + d17*t70; 
complex<T> co1 = s24*t45; 
complex<T> co2 = s34*t61; 
complex<T> co3 = -t59; 
complex<T> co4 = d19*s23*t70; 
complex<T> co5 = d20*s34*s45*t44; 
SeriesC<T> result = t17*(*CI_users[0]->get_value(mc,ind,mu)) + t19*(*CI_users[1]->get_value(mc,ind,mu)) + t43*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + t1*(*CI_users[5]->get_value(mc,ind,mu)) + t14*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + co1*(*CI_users[8]->get_value(mc,ind,mu)) + t50*(*CI_users[9]->get_value(mc,ind,mu)) + t16*(*CI_users[10]->get_value(mc,ind,mu)) + t56*(*CI_users[11]->get_value(mc,ind,mu)) + co2*(*CI_users[12]->get_value(mc,ind,mu)) + t18*(*CI_users[13]->get_value(mc,ind,mu)) + co3*(*CI_users[14]->get_value(mc,ind,mu)) + t42*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmqpmp_nfLT_wCI::\
C2q2g1ph_phqmqpmp_nfLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqpmp_nfLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, m, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqpmp nfLT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:50:28 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> s23 = S(2,3);
complex<T> s235 = SS(2,3,5);
complex<T> t4 = cube(spa24*spb25 + spa34*spb35); 
complex<T> t7 = square(spa23); 
complex<T> t8 = cube(spa24); 
complex<T> t9 = square(spa45); 
complex<T> t10 = cube(spb35); 
complex<T> d1 = spa23*spa25*spa45*T(3); d1 = T(1)/d1;
complex<T> d2 = spa25*cube(s23 - s235)*T(3); d2 = T(1)/d2;
complex<T> d3 = spa23*spa45*cube(s23 - s235)*T(3); d3 = T(1)/d3;
complex<T> t1 = d1*t8 + d2*t10*t7*t9 - d3*t4*square(spa25); 
complex<T> t2 = d1*t8 - d2*t10*t7*t9 + d3*t4*square(spa25); 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmqpmm_LT_wCI::\
C2q2g1ph_phqmqpmm_LT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqpmm_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, m, m}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqpmm LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:53:28 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb25 = SPB(2,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = -(spa34*spb34);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s14 = S(1,4);
complex<T> s35 = -(spa35*spb35);
complex<T> s245 = SS(2,4,5);
complex<T> s15 = S(1,5);
complex<T> s345 = SS(3,4,5);
complex<T> t17 = square(spa25*spb23 - spa45*spb34); 
complex<T> t18 = -(spa24*spb34); 
complex<T> t19 = square(spa24*spb34 + spa25*spb35); 
complex<T> t21 = -(spa45*T(3)); 
complex<T> t23 = cube(-(spa25*spb23) + spa45*spb34); 
complex<T> t28 = square(spb24); 
complex<T> t31 = spb23*spb34; 
complex<T> t33 = s235*spb25; 
complex<T> t34 = spa23*(spa23*spb35 + spa24*spb45); 
complex<T> t39 = square(mH2); 
complex<T> t40 = square(spa25); 
complex<T> t41 = square(spa24); 
complex<T> t42 = square(s245); 
complex<T> t43 = square(spa45); 
complex<T> t44 = square(spb23); 
complex<T> t45 = spa23*spb24 + spa35*spb45; 
complex<T> t47 = square(spb34); 
complex<T> t48 = spa25*spb24 + spa35*spb34; 
complex<T> t51 = s234*T(2); 
complex<T> t62 = cube(spa24); 
complex<T> t83 = coeff3mass1234<T>(mc,ind,2,3,4,5); 
complex<T> t84 = coeff3mass4123<T>(mc,ind,2,3,4,5); 
complex<T> t90 = spa24*spb23; 
complex<T> t107 = spa25*spa45; 
complex<T> t114 = cube(s245); 
complex<T> t119 = spa24*spa35; 
complex<T> t132 = spa34*spb35; 
complex<T> t162 = spa35*spb23; 
complex<T> d1 = cube(s24 + s34)*T(3); d1 = T(1)/d1;
complex<T> d2 = square(s24 + s34)*T(6); d2 = T(1)/d2;
complex<T> d4 = square(s25 + s35)*T(6); d4 = T(1)/d4;
complex<T> d6 = (s25 + s35)*spb25*T(2); d6 = T(1)/d6;
complex<T> d8 = square(s24 + s34)*T(3); d8 = T(1)/d8;
complex<T> d9 = spb24*square(s24 + s34)*T(2); d9 = T(1)/d9;
complex<T> d11 = s23*s234*spb34*T(6); d11 = T(1)/d11;
complex<T> d12 = s234*(s24 + s34)*spb34*T(3); d12 = T(1)/d12;
complex<T> d13 = cube(s25 + s35)*T(3); d13 = T(1)/d13;
complex<T> d14 = square(s25 + s35)*T(3); d14 = T(1)/d14;
complex<T> d15 = spb25*square(s25 + s35)*T(2); d15 = T(1)/d15;
complex<T> d17 = s23*s235*spb35*T(6); d17 = T(1)/d17;
complex<T> d19 = spa23*spb34*spb35*spb45*T(6); d19 = T(1)/d19;
complex<T> d20 = square(s23 + s24)*T(2); d20 = T(1)/d20;
complex<T> d22 = spb24*square(s23 + s24)*T(2); d22 = T(1)/d22;
complex<T> d25 = (s35 + s45)*spb45*T(2); d25 = T(1)/d25;
complex<T> d26 = spb45*square(s35 + s45)*T(2); d26 = T(1)/d26;
complex<T> d31 = spb34*spb45*(spa23*spb35 + spa24*spb45)*T(2); d31 = T(1)/d31;
complex<T> d42 = s23*s234*spb34; d42 = T(1)/d42;
complex<T> d44 = s235*spb35; d44 = T(1)/d44;
complex<T> d45 = s23*s235*spb35; d45 = T(1)/d45;
complex<T> d47 = spa23*spb34*spb35*spb45; d47 = T(1)/d47;
complex<T> d54 = s23*s234; d54 = T(1)/d54;
complex<T> d57 = spa23*spb35*spb45; d57 = T(1)/d57;
complex<T> d58 = spb45*(spa23*spb35 + spa24*spb45); d58 = T(1)/d58;
complex<T> d59 = spa23*spb34*spb35; d59 = T(1)/d59;
complex<T> d60 = spb34*(spa23*spb35 + spa24*spb45)*T(2); d60 = T(1)/d60;
complex<T> d67 = (spa23*spb35 + spa24*spb45)*T(2); d67 = T(1)/d67;
complex<T> t16 = square(spa45*spb35 + t90); 
complex<T> t20 = square(spa45*spb24 + t162); 
complex<T> t24 = cube(spa45*spb24 + t162); 
complex<T> t29 = s234*t48; 
complex<T> t30 = cube(spa45*spb35 + t90); 
complex<T> t46 = spa24*spb25 + t132; 
complex<T> t58 = d13*s235; 
complex<T> t63 = -t107; 
complex<T> t69 = -(d31*t19); 
complex<T> t76 = t45*t48; 
complex<T> t92 = spa35*t39; 
complex<T> t98 = d26*s345; 
complex<T> t106 = -t119; 
complex<T> t116 = d15*spb23; 
complex<T> t121 = spa25*t44; 
complex<T> t125 = -(spa34*t19); 
complex<T> t126 = -(d42*t17); 
complex<T> t145 = spa45*t31; 
complex<T> t200 = spa45*t18; 
complex<T> d3 = (s24 + s34)*t51; d3 = T(1)/d3;
complex<T> d5 = (s24 + s34)*spb24*t51; d5 = T(1)/d5;
complex<T> d7 = (s24 + s34)*spb34*t51; d7 = T(1)/d7;
complex<T> d10 = s234*(s24 + s34)*t28; d10 = T(1)/d10;
complex<T> d16 = (s25 + s35)*t33*T(3); d16 = T(1)/d16;
complex<T> d18 = spb35*t33*T(6); d18 = T(1)/d18;
complex<T> d21 = (s23 + s24)*t51; d21 = T(1)/d21;
complex<T> d23 = (s23 + s24)*spb24*t51; d23 = T(1)/d23;
complex<T> d24 = s234*(s23 + s24)*t28; d24 = T(1)/d24;
complex<T> d33 = spb25*spb45*t45*T(2); d33 = T(1)/d33;
complex<T> d50 = spb45*t45; d50 = T(1)/d50;
complex<T> d56 = spb35*t33; d56 = T(1)/d56;
complex<T> d61 = spb25*t45*T(2); d61 = T(1)/d61;
complex<T> d65 = t31*t48*T(2); d65 = T(1)/d65;
complex<T> d69 = t45*T(2); d69 = T(1)/d69;
complex<T> t6 = d14*spa25*spb35*t43 + t116*t132*t63 + d4*t63*t90 - spa25*spb35*t43*t58*T(2) + d6*(spb35*t43 + spa45*t90)*T(3) + d16*spa35*t16*T(11); 
complex<T> t66 = -t92; 
complex<T> t89 = d25*(spa24*spb34 + spa25*spb35)*t21 + spb35*t63*t98; 
complex<T> t108 = spa34*t16; 
complex<T> t118 = s23*t69; 
complex<T> t120 = spa45*t16; 
complex<T> t122 = t24*t47; 
complex<T> t129 = t40*t92; 
complex<T> t139 = d45*t16; 
complex<T> t173 = d21*spa24; 
complex<T> t196 = spa35*t145; 
complex<T> d27 = t29*t31; d27 = T(1)/d27;
complex<T> d28 = s235*spb23*t46; d28 = T(1)/d28;
complex<T> d29 = t33*t46; d29 = T(1)/d29;
complex<T> d30 = s234*t34*t46; d30 = T(1)/d30;
complex<T> d32 = s235*spa23*t76; d32 = T(1)/d32;
complex<T> d34 = t29*cube(spb24); d34 = T(1)/d34;
complex<T> d35 = spb34*t29*T(2); d35 = T(1)/d35;
complex<T> d36 = s235*t46*T(2); d36 = T(1)/d36;
complex<T> d37 = t33*t46*T(2); d37 = T(1)/d37;
complex<T> d38 = s235*spb35*t46; d38 = T(1)/d38;
complex<T> d39 = (spa23*spb35 + spa24*spb45)*t46*t51; d39 = T(1)/d39;
complex<T> d40 = s235*t76*T(2); d40 = T(1)/d40;
complex<T> d41 = spb23*t28*t29; d41 = T(1)/d41;
complex<T> d43 = t29*t31*T(2); d43 = T(1)/d43;
complex<T> d46 = s235*spb23*t46*T(2); d46 = T(1)/d46;
complex<T> d48 = t34*t46*t51; d48 = T(1)/d48;
complex<T> d49 = s235*spa23*t76*T(2); d49 = T(1)/d49;
complex<T> d51 = s235*t46; d51 = T(1)/d51;
complex<T> d52 = s235*spb23*spb35*t46; d52 = T(1)/d52;
complex<T> d53 = spb23*t29*cube(spb24); d53 = T(1)/d53;
complex<T> d55 = spb23*t29*T(2); d55 = T(1)/d55;
complex<T> d62 = spb23*t46*T(2); d62 = T(1)/d62;
complex<T> d63 = spb25*t46*T(2); d63 = T(1)/d63;
complex<T> d64 = spa23*t76*T(2); d64 = T(1)/d64;
complex<T> d66 = t34*t46*T(2); d66 = T(1)/d66;
complex<T> d68 = t29*cube(spb24)*T(2); d68 = T(1)/d68;
complex<T> d70 = s235*spb35*t46*T(2); d70 = T(1)/d70;
complex<T> t2 = d5*t106*t121 + d22*t119*t121 + d5*t119*t145 + d9*t119*t145 + d10*spa24*spb34*t20 + d24*spa24*spb34*t20 + d8*spa24*spb34*t43 - d7*spa24*t40*t44 + d2*t107*t90 + d20*t107*t90 + d3*t107*t90 + d21*spa25*t21*t90 - d1*t196*t41*T(2) + d1*s24*t107*t90*T(2) - d23*t119*t121*T(3) + d23*t119*t145*T(3) + spb34*t173*t43*T(3) + d12*spa24*t17*T(8); 
complex<T> t4 = d22*t106*t121 + d24*t18*t20 + d23*t119*t21*t31 + d20*t63*t90 + spb35*t107*t98 + d25*(spa24*spa45*spb34 + spb35*t107)*T(3) + d23*t119*t121*T(3) - spb34*t173*t43*T(3) + d21*t107*t90*T(3); 
complex<T> t5 = d5*t119*t121 + t107*t116*t132 + d10*t18*t20 + d5*t162*t200 + d9*t162*t200 - d14*spa25*spb35*t43 + d8*t18*t43 + d7*spa24*t40*t44 + d4*t107*t90 + d6*t21*t90 + d2*t63*t90 + d3*t63*t90 + d1*t196*t41*T(2) + spa25*spb35*t43*t58*T(2) - d1*s24*t107*t90*T(2) - d6*spb35*t43*T(3) - d12*spa24*t17*T(8) - d16*spa35*t16*T(11) - d18*t16*T(13) + d17*spa25*t16*T(13) - d11*spa24*t17*T(13) + d19*t19*T(13); 
complex<T> t50 = -(d32*spa35); 
complex<T> t65 = -t108; 
complex<T> t73 = -(d27*t23); 
complex<T> t91 = -t120; 
complex<T> t99 = d29*spa34; 
complex<T> t111 = d48*t39; 
complex<T> t128 = -(d39*t39); 
complex<T> t135 = d40*spb23; 
complex<T> t143 = d66*t39; 
complex<T> t3 = -(d65*s25*t23) + s25*t143*t62; 
complex<T> t10 = spa25*(d44*t16 + d50*t42) - d36*spa25*t108*T(3) + s25*(spa24*t126 + spa25*t139 + d47*t19 + d46*t16*t21 - d43*t23 + t111*t62 + d49*t129*T(3)); 
complex<T> t11 = d37*s34*t108 - d53*s34*t122 + d57*t125 + d49*s34*t129 + s34*spa25*t139 - d56*s34*t16 + d54*spa24*spa34*t17 + d58*spa34*t19 - d55*spa34*t23 + s34*t111*t62 + d46*s34*t91; 
complex<T> t12 = d37*s45*t108 + spa45*(-(d59*t19) + d60*t19 + d61*t42) + s45*(spa24*t126 + d49*t129 + spa25*t139 - d56*t16 - d43*t23 + t111*t62 + d46*t91); 
complex<T> t13 = s45*(d63*t108 + d64*t129 + d62*t91); 
complex<T> t36 = s34*(d63*t108 + d64*t129 + d62*t91); 
complex<T> t37 = -(d65*s45*t23) + s45*t143*t62; 
complex<T> t87 = t40*t50; 
complex<T> t184 = spb23*t128; 
complex<T> t201 = t135*t40; 
complex<T> t1 = s34*(d68*spa23*t122 + t184*t62); 
complex<T> t7 = d28*mH2*t120 + d32*s15*t129 + d31*mH2*t19 + d27*mH2*t23 + d33*mH2*t42 - d33*s15*t42 + d30*s15*t39*t62 + s15*t69 + s15*t73 + d28*s15*t91 - mH2*t16*t99 + s15*t16*t99 - d30*t62*cube(mH2) + t87*cube(mH2); 
complex<T> t8 = d28*s15*t120 + d53*s15*t122 + d41*spa24*t122 + d32*s24*t129 + d31*s15*t19 + d27*s24*t23 + d33*s15*t42 - d33*s24*t42 + s24*t69 + s15*t73 + s15*t39*t87 + d28*s24*t91 - s15*t16*t99 + s24*t16*t99; 
complex<T> t9 = d28*s24*t120 - d41*spa24*t122 + s24*(d31*t19 + d33*t42 + t73 + t39*t87 - t16*t99); 
complex<T> t14 = d32*s14*t129 + d52*s14*t30 - d52*s25*t30 + s25*t39*t87 + d51*spa25*t108*T(2) + d28*(-s14 + s25)*t120*T(2) + s14*t16*t99*T(2); 
complex<T> t15 = t118 + d34*spa23*t122 - d35*spa23*t23 + d38*spa23*t30 - d33*s23*t42 + t184*t62 + d37*s23*t65 + t201*t66 + d36*spa23*t91; 
complex<T> t38 = s25*(d70*spa23*t30 + t201*t66); 
complex<T> co1 = -t83; 
complex<T> co2 = -t84; 
complex<T> co3 = s345*t118; 
complex<T> co4 = -(d33*s23*t114); 
complex<T> co5 = d67*spa45*t125; 
complex<T> co6 = d69*t42*t63; 
SeriesC<T> result = t5*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu)) + t6*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + t89*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + t7*(*CI_users[7]->get_value(mc,ind,mu)) + t15*(*CI_users[8]->get_value(mc,ind,mu)) + t9*(*CI_users[9]->get_value(mc,ind,mu)) + t10*(*CI_users[10]->get_value(mc,ind,mu)) + t14*(*CI_users[11]->get_value(mc,ind,mu)) + t8*(*CI_users[12]->get_value(mc,ind,mu)) + t11*(*CI_users[13]->get_value(mc,ind,mu)) + t12*(*CI_users[14]->get_value(mc,ind,mu)) + co3*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + t36*(*CI_users[17]->get_value(mc,ind,mu)) + t13*(*CI_users[18]->get_value(mc,ind,mu)) + t3*(*CI_users[19]->get_value(mc,ind,mu)) + t37*(*CI_users[20]->get_value(mc,ind,mu)) + co5*(*CI_users[21]->get_value(mc,ind,mu)) + t1*(*CI_users[22]->get_value(mc,ind,mu)) + co6*(*CI_users[23]->get_value(mc,ind,mu)) + t38*(*CI_users[24]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmqpmm_RT_wCI::\
C2q2g1ph_phqmqpmm_RT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqpmm_RT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, m, m}, RT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqpmm RT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:54:04 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> s12 = S(1,2);
complex<T> s23 = -(spa23*spb23);
complex<T> s24 = -(spa24*spb24);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = -(spa35*spb35);
complex<T> s14 = S(1,4);
complex<T> s245 = SS(2,4,5);
complex<T> s13 = S(1,3);
complex<T> s15 = S(1,5);
complex<T> s345 = SS(3,4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> t15 = square(spb23); 
complex<T> t16 = square(spa24*spb34 + spa25*spb35); 
complex<T> t18 = square(spa24*spb23 + spa45*spb35); 
complex<T> t19 = square(spa25); 
complex<T> t20 = square(spa35*spb23 + spa45*spb24); 
complex<T> t21 = square(spa25*spb23 - spa45*spb34); 
complex<T> t22 = square(spa24*spb25 + spa34*spb35); 
complex<T> t23 = square(s23); 
complex<T> t24 = square(spa24); 
complex<T> t34 = -(spa24*T(2)); 
complex<T> t35 = square(s245); 
complex<T> t37 = square(spa35); 
complex<T> t38 = spa23*spb24 + spa35*spb45; 
complex<T> t39 = square(s345); 
complex<T> t40 = square(spa45); 
complex<T> t41 = square(spb35); 
complex<T> t54 = spa24*T(2); 
complex<T> t70 = coeff3mass4123<T>(mc,ind,2,3,4,5); 
complex<T> t71 = spa45*spb23; 
complex<T> t72 = spa25*spa35; 
complex<T> t84 = spa35*spb34; 
complex<T> t86 = cube(s245); 
complex<T> d1 = s234*s34; d1 = T(1)/d1;
complex<T> d2 = s34*(s24 + s34); d2 = T(1)/d2;
complex<T> d3 = s234*square(spb24); d3 = T(1)/d3;
complex<T> d4 = s234*spa34*square(spb24); d4 = T(1)/d4;
complex<T> d5 = s234*spb24; d5 = T(1)/d5;
complex<T> d6 = (s24 + s34)*spb24; d6 = T(1)/d6;
complex<T> d7 = s234*spa34*spb24; d7 = T(1)/d7;
complex<T> d8 = (s24 + s34)*spa34*spb24; d8 = T(1)/d8;
complex<T> d9 = s234*s34*spb24; d9 = T(1)/d9;
complex<T> d10 = (s25 + s35)*spb25; d10 = T(1)/d10;
complex<T> d11 = s234*spb24*spb34; d11 = T(1)/d11;
complex<T> d12 = (s24 + s34)*square(spb24); d12 = T(1)/d12;
complex<T> d13 = s234*spb24*spb34*square(s24 + s34)*T(2); d13 = T(1)/d13;
complex<T> d14 = s23*s234*spb34*T(2); d14 = T(1)/d14;
complex<T> d15 = s235*spb25*spb35*square(s25 + s35)*T(2); d15 = T(1)/d15;
complex<T> d16 = s23*s235*spb35*T(2); d16 = T(1)/d16;
complex<T> d17 = s235*spb25*spb35; d17 = T(1)/d17;
complex<T> d18 = spa23*spb34*spb35*spb45*T(2); d18 = T(1)/d18;
complex<T> d19 = (s23 + s24)*s34; d19 = T(1)/d19;
complex<T> d20 = (s23 + s24)*spa34*square(spb24); d20 = T(1)/d20;
complex<T> d21 = (s23 + s24)*spa34*spb24; d21 = T(1)/d21;
complex<T> d22 = (s23 + s24)*s34*spb24; d22 = T(1)/d22;
complex<T> d23 = s234*cube(spb24)*T(2); d23 = T(1)/d23;
complex<T> d24 = s234*cube(spb24)*square(s23 + s24)*T(2); d24 = T(1)/d24;
complex<T> d25 = s235*spb25*spb35*T(2); d25 = T(1)/d25;
complex<T> d26 = (s35 + s45)*spb45; d26 = T(1)/d26;
complex<T> d27 = spb34*spb45*(spa23*spb35 + spa24*spb45)*square(s35 + s45)*T(2); d27 = T(1)/d27;
complex<T> d28 = spb34*spb45*(spa23*spb35 + spa24*spb45)*T(2); d28 = T(1)/d28;
complex<T> d30 = s234*spb34*cube(spb24); d30 = T(1)/d30;
complex<T> d31 = s234*spb34; d31 = T(1)/d31;
complex<T> d32 = s235*spb35; d32 = T(1)/d32;
complex<T> d33 = spb34*spb35*spb45; d33 = T(1)/d33;
complex<T> d34 = s234*spb34*square(spb24); d34 = T(1)/d34;
complex<T> d38 = s234*cube(spb24); d38 = T(1)/d38;
complex<T> d39 = spb45*(spa23*spb35 + spa24*spb45)*T(2); d39 = T(1)/d39;
complex<T> d40 = spb34*(spa23*spb35 + spa24*spb45); d40 = T(1)/d40;
complex<T> d41 = spb34*spb45*(spa23*spb35 + spa24*spb45); d41 = T(1)/d41;
complex<T> d42 = s235*spb35*T(2); d42 = T(1)/d42;
complex<T> t9 = -((d41*s12 + d40*spa45)*t16); 
complex<T> t11 = -(d25*t18) + d15*t15*t19*t22 + d10*t54*t71 + d10*spb35*t40*T(2); 
complex<T> t17 = square(spa25*spb24 + t84); 
complex<T> t51 = d20*t15; 
complex<T> t53 = -t71; 
complex<T> t57 = -(d28*t16); 
complex<T> t83 = spa25*t71; 
complex<T> t91 = spa35*t71; 
complex<T> t92 = spa25*t18; 
complex<T> t98 = t15*t72; 
complex<T> t99 = d26*spa45; 
complex<T> t107 = d4*t37; 
complex<T> t117 = s23*t15; 
complex<T> d29 = spb25*spb45*t38*T(2); d29 = T(1)/d29;
complex<T> d35 = spb25*t38; d35 = T(1)/d35;
complex<T> d36 = spb25*spb45*t38; d36 = T(1)/d36;
complex<T> d37 = spb45*t38*T(2); d37 = T(1)/d37;
complex<T> t1 = -(d11*t15*t19) - d23*spb34*t20 + d24*spb34*t20*t23 + d13*t15*t17*t24 + d6*spa25*t53 + t37*t51*t54 + d5*t83 + d19*t34*t83 + d2*t54*t83 + d12*t53*t84 + d3*t71*t84 + d8*t34*t91 + d21*t54*t91 - d3*t98 + d22*t34*t98; 
complex<T> t2 = d27*t19*t39*t41 + t57 + spb34*t34*t99 - spa25*spb35*t99*T(2); 
complex<T> t3 = d28*t16 + d23*spb34*t20 - d24*spb34*t20*t23 - d27*t19*t39*t41 + t34*t37*t51 + t107*t15*t54 + d1*t34*t83 + d19*t54*t83 + d21*t34*t91 + d7*t54*t91 + d9*t34*t98 + d22*t54*t98 + spb34*t54*t99 + spa25*spb35*t99*T(2); 
complex<T> t4 = d28*mH2*t16 + d29*(mH2 - s15)*t35 + s15*t57; 
complex<T> t5 = d28*s15*t16 - d30*s15*t15*t17 - d34*spa24*t15*t17 + d29*s15*t35 - d29*s24*t35 + s24*t57; 
complex<T> t6 = d28*s24*t16 + d34*spa24*t15*t17 + d29*s24*t35; 
complex<T> t10 = spa45*(d40*t16 + d35*t35); 
complex<T> t12 = d11*t15*t19 - d15*t15*t19*t22 - d13*t15*t17*t24 + t107*t15*t34 + d5*spa25*t53 + d10*t34*t71 + d6*t83 + d2*t34*t83 + d1*t54*t83 + d3*t53*t84 + d12*t71*t84 + d7*t34*t91 + d8*t54*t91 + d3*t98 + d9*t54*t98 + d17*t18*T(2) - d10*spb35*t40*T(2) - d18*t16*T(3) + d14*spa24*t21*T(3) - d16*t92*T(3); 
complex<T> t13 = -(d17*s14*t18) - d32*t92; 
complex<T> t14 = -(d33*spb23*t16) + d30*t117*t17 - d31*spa24*t21 + d32*t92; 
complex<T> t32 = -((d36*s13 + d35*spa45)*t35); 
complex<T> t101 = -(spa34*t17); 
complex<T> t104 = s25*t57; 
complex<T> t7 = t104 - d37*spa25*t35; 
complex<T> t8 = d38*t101*t15 - d39*spa34*t16 - d29*s34*t35; 
complex<T> co1 = -t70; 
complex<T> co2 = s345*t104; 
complex<T> co3 = -(d29*s34*t86); 
complex<T> co4 = -(d42*s23*t92); 
complex<T> co5 = d23*t101*t117; 
SeriesC<T> result = t12*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t11*(*CI_users[2]->get_value(mc,ind,mu)) + t3*(*CI_users[3]->get_value(mc,ind,mu)) + t2*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t4*(*CI_users[6]->get_value(mc,ind,mu)) + t14*(*CI_users[7]->get_value(mc,ind,mu)) + t6*(*CI_users[8]->get_value(mc,ind,mu)) + t32*(*CI_users[9]->get_value(mc,ind,mu)) + t7*(*CI_users[10]->get_value(mc,ind,mu)) + t13*(*CI_users[11]->get_value(mc,ind,mu)) + t5*(*CI_users[12]->get_value(mc,ind,mu)) + t8*(*CI_users[13]->get_value(mc,ind,mu)) + t9*(*CI_users[14]->get_value(mc,ind,mu)) + t10*(*CI_users[15]->get_value(mc,ind,mu)) + co2*(*CI_users[16]->get_value(mc,ind,mu)) + co3*(*CI_users[17]->get_value(mc,ind,mu)) + co4*(*CI_users[18]->get_value(mc,ind,mu)) + co5*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmqpmm_nfLT_wCI::\
C2q2g1ph_phqmqpmm_nfLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmqpmm_nfLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, qp, m, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqpmm nfLT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:54:08 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> s24 = S(2,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = -(spa35*spb35);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> t4 = square(spa24*spb23 + spa45*spb35); 
complex<T> t6 = square(spa25*spb23 - spa45*spb34); 
complex<T> t10 = spa45*spb23; 
complex<T> t13 = s23*T(3); 
complex<T> t30 = square(spa45); 
complex<T> t35 = square(spa24); 
complex<T> t38 = square(spa25); 
complex<T> t43 = spa35*spb34; 
complex<T> d1 = cube(s24 + s34)*T(3); d1 = T(1)/d1;
complex<T> d2 = square(s24 + s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = cube(s25 + s35)*T(3); d3 = T(1)/d3;
complex<T> d4 = square(s25 + s35)*T(3); d4 = T(1)/d4;
complex<T> d7 = spa23*spb34*spb35*spb45*T(3); d7 = T(1)/d7;
complex<T> d10 = s235*spb25*spb35*T(3); d10 = T(1)/d10;
complex<T> t17 = d3*spa34; 
complex<T> t25 = d1*s24; 
complex<T> t26 = d3*s25; 
complex<T> t34 = spa25*t10; 
complex<T> t51 = t35*t43; 
complex<T> d5 = (s24 + s34)*spb34*t13; d5 = T(1)/d5;
complex<T> d6 = (s25 + s35)*spb35*t13; d6 = T(1)/d6;
complex<T> d8 = s234*spb34*t13; d8 = T(1)/d8;
complex<T> d9 = s235*spb35*t13; d9 = T(1)/d9;
complex<T> t1 = d2*spa24*(-(spb34*t30) + t34) + (d1*t10*t51 - spa24*(t25*t34 + d5*t6 - d8*t6))*T(2); 
complex<T> t28 = d6*t4; 
complex<T> t40 = t10*t17; 
complex<T> t2 = d2*spa24*spb34*t30 - d4*spa25*spb35*t30 - d2*spa24*t34 - d4*spa24*t34 - spa25*t28*T(2) + spa24*t25*t34*T(2) + spa24*t26*t34*T(2) - spb35*t38*t40*T(2) - d1*t10*t51*T(2) + d5*spa24*t6*T(2) - d7*square(spa24*spb34 + spa25*spb35)*T(2); 
complex<T> t3 = d4*spa24*t34 - spa24*t26*t34*T(2) + d10*t4*T(2) + spb35*t38*t40*T(2) + spa25*(d4*spb35*t30 + t28*T(2) - d9*t4*T(2)); 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C2q2g1ph_phqmqppp_LT_wCI::\
C2q2g1ph_phqmqppp_LT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1ph_phqmqppp_LT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqppp LT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1ph_phqmqppp_RT_wCI::\
C2q2g1ph_phqmqppp_RT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1ph_phqmqppp_RT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqppp RT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1ph_phqmqppp_nfLT_wCI::\
C2q2g1ph_phqmqppp_nfLT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1ph_phqmqppp_nfLT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmqppp nfLT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1ph_phqmpqpm_LT_wCI::\
C2q2g1ph_phqmpqpm_LT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c3, c2, c5, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmpqpm_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, p, qp, m}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmpqpm LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:54:09 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spa24 = SPA(2,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa45 = SPA(4,5);
complex<T> s15 = S(1,5);
complex<T> s25 = S(2,5);
complex<T> s45 = S(4,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s345 = SS(3,4,5);
complex<T> t3 = square(spa25*spb24 + spa35*spb34); 
complex<T> t9 = square(spa25); 
complex<T> t12 = square(spa45); 
complex<T> t13 = square(spb34); 
complex<T> t20 = spa45*spb34; 
complex<T> d1 = spa23*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = (s23 - s234)*spa34; d2 = T(1)/d2;
complex<T> d3 = spa23*spa34*square(s23 - s234)*T(2); d3 = T(1)/d3;
complex<T> d4 = spa23*spa34*T(3); d4 = T(1)/d4;
complex<T> d5 = (-s345 + s45)*spa34; d5 = T(1)/d5;
complex<T> d6 = spa34*square(s345 - s45)*T(2); d6 = T(1)/d6;
complex<T> d7 = spa23*spa34; d7 = T(1)/d7;
complex<T> d8 = spa34*T(2); d8 = T(1)/d8;
complex<T> d9 = spa23*T(2); d9 = T(1)/d9;
complex<T> t1 = -(d3*t3*square(spa24)) + d2*spa25*t20*T(2) - d4*t9*T(5); 
complex<T> t2 = -(d1*t9) + d3*t3*square(spa24) - d2*spa25*t20*T(2); 
complex<T> t7 = d7*(-mH2 + s15)*t9*T(2); 
complex<T> t15 = d1*s234; 
complex<T> t17 = d8*spb23; 
complex<T> t21 = d6*spa23; 
complex<T> t26 = d9*spb34; 
complex<T> t8 = -(t12*t13*t21) - d5*spa25*t20*T(2); 
complex<T> t23 = t12*t13*t21 + d5*spa25*t20*T(2); 
complex<T> t24 = (s345*t15 + mH2*t26)*t9; 
complex<T> t28 = (s235*t15 + mH2*t17)*t9; 
complex<T> co1 = -(s25*t17*t9); 
complex<T> co2 = -(s45*t26*t9); 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t23*(*CI_users[2]->get_value(mc,ind,mu)) + t8*(*CI_users[3]->get_value(mc,ind,mu)) + t7*(*CI_users[4]->get_value(mc,ind,mu)) + t28*(*CI_users[5]->get_value(mc,ind,mu)) + t24*(*CI_users[6]->get_value(mc,ind,mu)) + co1*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmpqpm_nfLT_wCI::\
C2q2g1ph_phqmpqpm_nfLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmpqpm_nfLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, p, qp, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmpqpm nfLT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:54:10 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> t1 = square(spa25); 
complex<T> d1 = spa23*spa34*T(3); d1 = T(1)/d1;
complex<T> co1 = d1*t1*T(2); 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phqmmqpp_LT_wCI::\
C2q2g1ph_phqmmqpp_LT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c23, c5));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phqmmqpp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{ph, qm, m, qp, p}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmmqpp LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:54:12 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa23 = SPA(2,3);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa24 = SPA(2,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s23 = S(2,3);
complex<T> s34 = S(3,4);
complex<T> s14 = S(1,4);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t2 = square(spa23*spb35 + spa24*spb45); 
complex<T> t5 = square(spa35); 
complex<T> t10 = square(spa23); 
complex<T> t13 = square(spa34); 
complex<T> t14 = square(spb45); 
complex<T> t21 = -(spa34*T(2)); 
complex<T> t29 = spa34*T(2); 
complex<T> d1 = spa25*spa45*T(2); d1 = T(1)/d1;
complex<T> d2 = (-s245 + s25)*spa45; d2 = T(1)/d2;
complex<T> d3 = spa45*square(s245 - s25)*T(2); d3 = T(1)/d3;
complex<T> d4 = (s34 - s345)*spa45; d4 = T(1)/d4;
complex<T> d5 = spa25*spa45*square(s34 - s345)*T(2); d5 = T(1)/d5;
complex<T> d6 = spa25*spa45; d6 = T(1)/d6;
complex<T> d7 = spa25; d7 = T(1)/d7;
complex<T> d8 = spa45; d8 = T(1)/d8;
complex<T> d9 = T(2); d9 = T(1)/d9;
complex<T> t1 = d1*t10 + d4*spa23*spb45*t29 - d5*t2*t5; 
complex<T> t7 = -(d3*spa25*t13*t14) + d2*spa23*spb45*t21 - d1*t10*T(3); 
complex<T> t8 = d3*spa25*t13*t14 + d2*spa23*spb45*t29; 
complex<T> t16 = d1*mH2; 
complex<T> t22 = -(d1*s234); 
complex<T> t24 = d8*spb25; 
complex<T> t28 = d6*t10; 
complex<T> t33 = -(d1*t10); 
complex<T> t6 = (mH2 - s15)*t28*T(2); 
complex<T> t9 = d4*spa23*spb45*t21 + t33 + d5*t2*t5; 
complex<T> t18 = -t24; 
complex<T> t20 = -t28; 
complex<T> t27 = d7*spb45*t10 + s13*t28; 
complex<T> t41 = t10*t16; 
complex<T> t32 = -(d7*spb45*t10) + s12*t20; 
complex<T> t36 = t10*t18 + s14*t20; 
complex<T> t39 = s235*t10*t22 + s23*t41; 
complex<T> t42 = s24*t20 + s15*t28; 
complex<T> t44 = s345*t10*t22 + s34*t41; 
complex<T> co1 = s24*t28; 
complex<T> co2 = t10*t24; 
complex<T> co3 = -(d9*spb25*spb45*t10); 
complex<T> co4 = s23*s34*t33; 
SeriesC<T> result = t7*(*CI_users[0]->get_value(mc,ind,mu)) + t8*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t9*(*CI_users[3]->get_value(mc,ind,mu)) + t6*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t27*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t36*(*CI_users[8]->get_value(mc,ind,mu)) + t42*(*CI_users[9]->get_value(mc,ind,mu)) + t32*(*CI_users[10]->get_value(mc,ind,mu)) + t39*(*CI_users[11]->get_value(mc,ind,mu)) + t44*(*CI_users[12]->get_value(mc,ind,mu)) + co3*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[14]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C2q2g1ph_phqmmqpp_nfLT_wCI::\
C2q2g1ph_phqmmqpp_nfLT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1ph_phqmmqpp_nfLT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phqmmqpp nfLT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1ph_phdqmqppm_LT_wCI::\
C2q2g1ph_phdqmqppm_LT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c23, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqppm_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, p, m}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqppm LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:54:31 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa35 = SPA(3,5);
complex<T> s34 = S(3,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s12 = S(1,2);
complex<T> s14 = S(1,4);
complex<T> s25 = -(spa25*spb25);
complex<T> s23 = -(spa23*spb23);
complex<T> s235 = SS(2,3,5);
complex<T> s234 = SS(2,3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t3 = square(spa25*spb24 + spa35*spb34); 
complex<T> t4 = square(s23 - s235); 
complex<T> t6 = square(spb35); 
complex<T> t9 = cube(spa25*spb24 + spa35*spb34); 
complex<T> t12 = square(spb34); 
complex<T> t13 = square(spa25); 
complex<T> t14 = -(spb23*spb45); 
complex<T> t16 = -(spb24*T(2)); 
complex<T> t23 = square(s34); 
complex<T> t27 = cube(spa25); 
complex<T> t28 = square(spb23); 
complex<T> t29 = square(spb45); 
complex<T> t38 = cube(spb34); 
complex<T> t58 = spa25*spb34; 
complex<T> t60 = s34*spa45; 
complex<T> d1 = (s23 - s235)*spb25; d1 = T(1)/d1;
complex<T> d2 = (s23 - s235)*spb23*spb45; d2 = T(1)/d2;
complex<T> d3 = spb23*spb25*spb45*T(2); d3 = T(1)/d3;
complex<T> d4 = spb23*spb35*spb45*T(6); d4 = T(1)/d4;
complex<T> d7 = spb23*spb45*cube(s23 - s235)*T(3); d7 = T(1)/d7;
complex<T> d8 = cube(s23 - s235); d8 = T(1)/d8;
complex<T> d10 = spb35*cube(s23 - s235)*T(3); d10 = T(1)/d10;
complex<T> d11 = spb23*spb25*spb45*T(3); d11 = T(1)/d11;
complex<T> d12 = (-s245 + s45)*spb25; d12 = T(1)/d12;
complex<T> d13 = spb25*square(s245 - s45)*T(2); d13 = T(1)/d13;
complex<T> d14 = spb23*spb25*spb45; d14 = T(1)/d14;
complex<T> d15 = spb35*spb45; d15 = T(1)/d15;
complex<T> d16 = spb23*spb45; d16 = T(1)/d16;
complex<T> d17 = spb23*spb35*spb45; d17 = T(1)/d17;
complex<T> d18 = spb23*spb25; d18 = T(1)/d18;
complex<T> d19 = spb23*spb35; d19 = T(1)/d19;
complex<T> d20 = spb23*T(2); d20 = T(1)/d20;
complex<T> d21 = spb23*spb45*T(2); d21 = T(1)/d21;
complex<T> d22 = spb25*T(2); d22 = T(1)/d22;
complex<T> d23 = spb25*spb45*T(2); d23 = T(1)/d23;
complex<T> d24 = spb23*spb25*T(2); d24 = T(1)/d24;
complex<T> d25 = spb23*spb35*T(2); d25 = T(1)/d25;
complex<T> d26 = spb35*spb45*T(2); d26 = T(1)/d26;
complex<T> t32 = d3*s235; 
complex<T> t33 = d23*spa23; 
complex<T> t35 = d18*spa45; 
complex<T> t42 = spb24*t12; 
complex<T> t43 = d21*spa25; 
complex<T> t44 = d14*T(2); 
complex<T> t45 = d8*spa35; 
complex<T> t46 = d20*spa45; 
complex<T> t51 = d13*spb24*t13*t14 + d12*t16*t58; 
complex<T> t54 = spb23*t13; 
complex<T> t69 = t12*t16; 
complex<T> t73 = d16*spa25; 
complex<T> t74 = d8*spb24; 
complex<T> d5 = spb23*spb45*t4; d5 = T(1)/d5;
complex<T> d6 = spb23*spb25*spb45*t4*T(2); d6 = T(1)/d6;
complex<T> d9 = spb35*t4*T(2); d9 = T(1)/d9;
complex<T> t2 = d9*spb34*t13*t14 + d10*t27*t28*t29 - d5*spb34*spb35*t3 + d2*spa35*t38 + d3*t42 + d2*spa25*t42 + spb34*t13*t14*t45 - d6*spb24*t3*t6 + t14*t27*t74 + d1*spb24*t58*T(2) + d7*t6*t9*T(4) + d4*t38*T(5); 
complex<T> t10 = d14*(mH2 - s15)*t42*T(3); 
complex<T> t18 = -t33; 
complex<T> t20 = -t32; 
complex<T> t25 = -t35; 
complex<T> t26 = -t42; 
complex<T> t41 = spb24*(d13*spb45*t54 + d12*t58*T(2)); 
complex<T> t52 = d17*(s14 - s25)*t38 + t69*(d14*s14 + t73); 
complex<T> t53 = s34*t42; 
complex<T> t63 = s15*t42*t44 + d14*s24*t69; 
complex<T> t84 = spa25*t46; 
complex<T> t86 = s24*t42; 
complex<T> t1 = d2*spa25*t26 - d10*t27*t28*t29 + d5*spb34*spb35*t3 - d2*spa35*t38 + d9*spb34*spb45*t54 + spb34*spb45*t45*t54 + d1*t16*t58 + d6*spb24*t3*t6 + spb23*spb45*t27*t74 - d7*t6*t9*T(4) - d4*t38*T(5) + d11*t42*T(5); 
complex<T> t37 = -(d17*s34*t38) + d14*t53; 
complex<T> t57 = d14*s12*t26 + d17*s12*t38 + d19*spa45*t38 + t25*t42; 
complex<T> t68 = s245*t20*t42 + s235*t26*t43 + s24*t26*t43 + s34*t26*t43 + t42*t84; 
complex<T> t71 = d22*spa23*spa45*t42 + s235*t18*t42 + s234*t20*t42 + t18*t53 + t18*t86; 
complex<T> t76 = t32*t53 + d3*(s234*s345*t26 + t23*t42 + s24*t53) + d24*t26*t60; 
complex<T> co1 = d15*spa23*t38; 
complex<T> co2 = t44*t86; 
complex<T> co3 = t42*t73; 
complex<T> co4 = t35*t42; 
complex<T> co5 = t33*t53; 
complex<T> co6 = d25*t38*t60; 
complex<T> co7 = t26*t84; 
complex<T> co8 = d26*s25*spa23*t38; 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t51*(*CI_users[2]->get_value(mc,ind,mu)) + t41*(*CI_users[3]->get_value(mc,ind,mu)) + t10*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + co2*(*CI_users[6]->get_value(mc,ind,mu)) + co3*(*CI_users[7]->get_value(mc,ind,mu)) + t52*(*CI_users[8]->get_value(mc,ind,mu)) + t63*(*CI_users[9]->get_value(mc,ind,mu)) + t37*(*CI_users[10]->get_value(mc,ind,mu)) + t57*(*CI_users[11]->get_value(mc,ind,mu)) + co4*(*CI_users[12]->get_value(mc,ind,mu)) + t68*(*CI_users[13]->get_value(mc,ind,mu)) + t71*(*CI_users[14]->get_value(mc,ind,mu)) + t76*(*CI_users[15]->get_value(mc,ind,mu)) + co5*(*CI_users[16]->get_value(mc,ind,mu)) + co6*(*CI_users[17]->get_value(mc,ind,mu)) + co7*(*CI_users[18]->get_value(mc,ind,mu)) + co8*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmqppm_RT_wCI::\
C2q2g1ph_phdqmqppm_RT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c3, c145));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c4, c135));
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c45, c13));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqppm_RT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, p, m}, RT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqppm RT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:54:36 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> s14 = S(1,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s34 = S(3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s15 = S(1,5);
complex<T> s24 = -(spa24*spb24);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t2 = spb23*spb45; 
complex<T> t3 = square(spa24*spb34 + spa25*spb35); 
complex<T> t7 = cube(spb24); 
complex<T> t12 = square(spb34); 
complex<T> t16 = square(spa25); 
complex<T> t18 = square(spb23); 
complex<T> t19 = square(spb45); 
complex<T> t23 = -(spa25*T(2)); 
complex<T> t31 = spa25*T(2); 
complex<T> d1 = spb25*spb35*T(2); d1 = T(1)/d1;
complex<T> d3 = (s23 - s235)*spb25*spb35; d3 = T(1)/d3;
complex<T> d4 = spb25*spb35*square(s23 - s235)*T(2); d4 = T(1)/d4;
complex<T> d5 = (-s245 + s45)*spb25; d5 = T(1)/d5;
complex<T> d8 = spb25*spb35; d8 = T(1)/d8;
complex<T> d9 = spb25*spb45; d9 = T(1)/d9;
complex<T> d10 = spb23*spb25; d10 = T(1)/d10;
complex<T> d11 = spb35; d11 = T(1)/d11;
complex<T> d12 = spb23*spb25*T(2); d12 = T(1)/d12;
complex<T> d13 = spb35*T(2); d13 = T(1)/d13;
complex<T> t9 = -(d4*t16*t18*t19) + d3*spb34*t2*t23 - d1*t12*T(3); 
complex<T> t17 = -(d12*spa45); 
complex<T> t29 = spb24*t12; 
complex<T> t30 = d10*spa45; 
complex<T> t35 = -(d8*t12); 
complex<T> t47 = spa25*t12; 
complex<T> d2 = spb25*t2*T(2); d2 = T(1)/d2;
complex<T> d6 = spb25*t2*square(s245 - s45)*T(2); d6 = T(1)/d6;
complex<T> d7 = spb25*t2; d7 = T(1)/d7;
complex<T> t1 = d2*t29 + d5*spb24*spb34*t31 - d6*t3*t7; 
complex<T> t8 = d8*s14*t12 + d11*t47; 
complex<T> t10 = d4*t16*t18*t19 + d3*spb34*t2*t31 + d1*t12*T(3) - d2*t29*T(3); 
complex<T> t13 = -(d7*spb24); 
complex<T> t22 = -t29; 
complex<T> t46 = d7*t29; 
complex<T> t11 = d2*t22 + d5*spb24*spb34*t23 + d6*t3*t7; 
complex<T> t34 = d9*spa23*t22 + s23*t35; 
complex<T> t39 = t12*t13; 
complex<T> t41 = d2*s245*s345*t22 + (s235 + s24 + s34 + s45)*t17*t29; 
complex<T> t28 = s15*t39 + mH2*t46; 
complex<T> t38 = s24*t39 + s15*t46; 
complex<T> t43 = d8*s12*t12 + t22*t30 + s45*t35 + s12*t39; 
complex<T> t45 = t22*t30 + s13*t39; 
complex<T> co1 = s24*t46; 
complex<T> co2 = s34*t35; 
complex<T> co3 = t29*t30; 
complex<T> co4 = -(d1*s34*s45*t12); 
complex<T> co5 = d13*s23*t47; 
SeriesC<T> result = t10*(*CI_users[0]->get_value(mc,ind,mu)) + t9*(*CI_users[1]->get_value(mc,ind,mu)) + t11*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t28*(*CI_users[4]->get_value(mc,ind,mu)) + t34*(*CI_users[5]->get_value(mc,ind,mu)) + co1*(*CI_users[6]->get_value(mc,ind,mu)) + t45*(*CI_users[7]->get_value(mc,ind,mu)) + t8*(*CI_users[8]->get_value(mc,ind,mu)) + t38*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t43*(*CI_users[11]->get_value(mc,ind,mu)) + co3*(*CI_users[12]->get_value(mc,ind,mu)) + t41*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[14]->get_value(mc,ind,mu)) + co5*(*CI_users[15]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmqppm_nfLT_wCI::\
C2q2g1ph_phdqmqppm_nfLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqppm_nfLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, p, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqppm nfLT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:54:36 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb45 = SPB(4,5);
complex<T> spa25 = SPA(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> s23 = S(2,3);
complex<T> s235 = SS(2,3,5);
complex<T> t3 = cube(s23 - s235); 
complex<T> t4 = cube(spa25*spb24 + spa35*spb34); 
complex<T> t7 = square(spb34); 
complex<T> t8 = cube(spa25); 
complex<T> t9 = square(spb23); 
complex<T> t10 = square(spb45); 
complex<T> d1 = spb25*spb35*T(3); d1 = T(1)/d1;
complex<T> d2 = spb23*spb25*spb45*T(3); d2 = T(1)/d2;
complex<T> t12 = -(d2*spb24); 
complex<T> d3 = spb23*spb45*t3*T(3); d3 = T(1)/d3;
complex<T> d4 = spb35*t3*T(3); d4 = T(1)/d4;
complex<T> t1 = d1*t7 + t12*t7 - d4*t10*t8*t9 - d3*t4*square(spb35); 
complex<T> t2 = -(d1*t7) + t12*t7 + d4*t10*t8*t9 + d3*t4*square(spb35); 
SeriesC<T> result = t1*(*CI_users[0]->get_value(mc,ind,mu)) + t2*(*CI_users[1]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmqpmp_LT_wCI::\
C2q2g1ph_phdqmqpmp_LT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c25, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c23, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqpmp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, m, p}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqpmp LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:55:15 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb24 = SPB(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb34 = SPB(3,4);
complex<T> spa34 = SPA(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> s25 = S(2,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s13 = S(1,3);
complex<T> s15 = S(1,5);
complex<T> s24 = -(spa24*spb24);
complex<T> s14 = S(1,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s245 = SS(2,4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> t10 = square(spa23*spb35 + spa24*spb45); 
complex<T> t14 = spb45*T(3); 
complex<T> t15 = square(spb24); 
complex<T> t23 = square(spb35); 
complex<T> t24 = square(spb25); 
complex<T> t25 = square(spa24); 
complex<T> t26 = square(spb34); 
complex<T> t28 = s23 - s234; 
complex<T> t29 = cube(spb24); 
complex<T> t43 = cube(spb35); 
complex<T> t44 = cube(spb25); 
complex<T> t45 = spa24*spb23; 
complex<T> t47 = square(spb45); 
complex<T> t58 = square(spb23); 
complex<T> t74 = cube(spa24); 
complex<T> t88 = spa24*spb25; 
complex<T> d7 = spb23*spb24*spb45*T(2); d7 = T(1)/d7;
complex<T> d13 = spb24*square(s234 - s34)*T(2); d13 = T(1)/d13;
complex<T> d14 = (-s234 + s34)*spb24; d14 = T(1)/d14;
complex<T> d15 = spb23*spb34*spb45*T(6); d15 = T(1)/d15;
complex<T> d17 = spa23*spb34*square(s234 - s34); d17 = T(1)/d17;
complex<T> d19 = spb24*square(s245 - s25)*T(2); d19 = T(1)/d19;
complex<T> d20 = (-s245 + s25)*spb24; d20 = T(1)/d20;
complex<T> d21 = spb24*square(s245 - s45)*T(2); d21 = T(1)/d21;
complex<T> d23 = spa45*spb34*square(s245 - s25); d23 = T(1)/d23;
complex<T> d24 = spb23*spb34*spb45; d24 = T(1)/d24;
complex<T> d28 = spb23*spb34; d28 = T(1)/d28;
complex<T> d30 = spb23*spb45; d30 = T(1)/d30;
complex<T> d31 = spb23*spb34*T(2); d31 = T(1)/d31;
complex<T> d32 = spb23*spb45*T(2); d32 = T(1)/d32;
complex<T> d33 = spb23*spb34*spb45*T(2); d33 = T(1)/d33;
complex<T> d34 = spb34*T(2); d34 = T(1)/d34;
complex<T> d35 = spb45*T(2); d35 = T(1)/d35;
complex<T> d36 = spb34*spb45*T(2); d36 = T(1)/d36;
complex<T> d37 = spb23*T(2); d37 = T(1)/d37;
complex<T> t4 = square(spa34*spb35 + t88); 
complex<T> t5 = -t45; 
complex<T> t9 = square(spa45*spb35 + t45); 
complex<T> t34 = -(d32*spa34); 
complex<T> t38 = d24*s15; 
complex<T> t49 = d33*s235; 
complex<T> t50 = d36*spa23; 
complex<T> t51 = d28*spa45; 
complex<T> t52 = d17*t10; 
complex<T> t68 = -(d14*T(3)); 
complex<T> t70 = d37*spa34; 
complex<T> t77 = t47*t58; 
complex<T> t78 = t26*t44; 
complex<T> t79 = s25*t43; 
complex<T> t97 = spb23*t25; 
complex<T> t106 = spa24*t24; 
complex<T> d1 = spa34*t15; d1 = T(1)/d1;
complex<T> d2 = spb34*t28*T(3); d2 = T(1)/d2;
complex<T> d3 = spa34*spb24*t28; d3 = T(1)/d3;
complex<T> d4 = spb24*square(t28)*T(2); d4 = T(1)/d4;
complex<T> d5 = spa34*spb23*spb45*t15; d5 = T(1)/d5;
complex<T> d6 = spb23*spb45*t28; d6 = T(1)/d6;
complex<T> d8 = spb23*spb34*t14; d8 = T(1)/d8;
complex<T> d9 = spa34*spb23*spb45*t28; d9 = T(1)/d9;
complex<T> d10 = spb23*t14*square(t28); d10 = T(1)/d10;
complex<T> d11 = spb34*square(t28)*T(3); d11 = T(1)/d11;
complex<T> d12 = spb34*cube(t28)*T(3); d12 = T(1)/d12;
complex<T> d16 = spa23*spb34*t15; d16 = T(1)/d16;
complex<T> d18 = (-s245 + s45)*t15; d18 = T(1)/d18;
complex<T> d22 = spa45*spb34*t15; d22 = T(1)/d22;
complex<T> d25 = spb45*t29; d25 = T(1)/d25;
complex<T> d26 = spb23*spb45*t15; d26 = T(1)/d26;
complex<T> d27 = spb23*t29; d27 = T(1)/d27;
complex<T> d29 = spb23*spb45*t29; d29 = T(1)/d29;
complex<T> d38 = spb45*t29*T(2); d38 = T(1)/d38;
complex<T> d39 = spb23*t29*T(2); d39 = T(1)/d39;
complex<T> t17 = -(d24*s24*t43) + t38*t43 + d29*s15*t78 + d26*spa24*t78; 
complex<T> t19 = d24*(-(s14*t43) + t79); 
complex<T> t35 = -t49; 
complex<T> t36 = -t50; 
complex<T> t57 = d19*spb34*t24*t25 - d22*spa24*spb45*t58 + d23*spa24*spb45*t9 - d20*spb35*t88*T(3); 
complex<T> t64 = -t78; 
complex<T> t65 = -t79; 
complex<T> t82 = d18*spb45; 
complex<T> t85 = spa45*t70; 
complex<T> t86 = (d24*mH2 - t38)*t43*T(3); 
complex<T> t87 = d13*spb34*t24*t25 + d16*t47*t5 + t45*t52 + spb35*t68*t88; 
complex<T> t109 = spb45*t97; 
complex<T> t120 = d5*spa24; 
complex<T> t129 = spa45*t78; 
complex<T> t1 = -(d21*spb25*t109) - d19*spb34*t24*t25 + d22*spa24*spb45*t58 + spb25*t45*t82 - d23*spa24*spb45*t9 + d18*spb34*t106*T(2) + d20*spb35*t88*T(3); 
complex<T> t2 = d1*t106 - d2*spa24*t23 + d7*spb25*t23 + d3*t4 - d4*spb34*t4 + d10*spb34*spb35*t4 + d8*t43 + d11*spb35*t109*T(2) - d9*spb35*t4*T(2) + spb34*t120*t44*T(2) - d12*t74*t77*T(2) + d6*spa34*t43*T(3) + d6*t23*t88*T(3); 
complex<T> t3 = -(d1*t106) + d2*spa24*t23 - d7*spb25*t23 - d13*spb34*t24*t25 - d3*t4 + d4*spb34*t4 - d10*spb34*spb35*t4 + d16*t45*t47 + t5*t52 - d11*spb35*t109*T(2) + d9*spb35*t4*T(2) - spb34*t120*t44*T(2) + d12*t74*t77*T(2) - d6*spa34*t43*T(3) + d14*spb35*t88*T(3) - d6*t23*t88*T(3) + d15*t43*T(11); 
complex<T> t18 = d24*s24*t43 + d26*spa24*t64; 
complex<T> t21 = d27*t129 - d24*s13*t43 - t43*t51 + d29*s13*t78; 
complex<T> t22 = d21*spb25*t109 + spb25*t5*t82 - d18*spb34*t106*T(2); 
complex<T> t76 = (d35*spa23*spa34 + d34*spa23*spa45 + s234*t35 + s235*t36 + s24*t36)*t43; 
complex<T> t94 = s245*t35*t43 + d31*spa45*t65 + (d33*s24 + t34 + t49)*t79; 
complex<T> t122 = d29*t64; 
complex<T> t134 = t43*t85; 
complex<T> t20 = s34*t122 - d30*spa34*t43; 
complex<T> t63 = t134 + (-(d33*s234*s345) + (s235 + s24 + s34)*t34)*t43; 
complex<T> co1 = d25*spa23*t78; 
complex<T> co2 = s25*t122; 
complex<T> co3 = t43*t51; 
complex<T> co4 = d38*s34*spa23*t78; 
complex<T> co5 = -t134; 
complex<T> co6 = d39*s25*t129; 
complex<T> co7 = t50*t79; 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t3*(*CI_users[1]->get_value(mc,ind,mu)) + t1*(*CI_users[2]->get_value(mc,ind,mu)) + t57*(*CI_users[3]->get_value(mc,ind,mu)) + t87*(*CI_users[4]->get_value(mc,ind,mu)) + t22*(*CI_users[5]->get_value(mc,ind,mu)) + t86*(*CI_users[6]->get_value(mc,ind,mu)) + co1*(*CI_users[7]->get_value(mc,ind,mu)) + t18*(*CI_users[8]->get_value(mc,ind,mu)) + t21*(*CI_users[9]->get_value(mc,ind,mu)) + co2*(*CI_users[10]->get_value(mc,ind,mu)) + t19*(*CI_users[11]->get_value(mc,ind,mu)) + t17*(*CI_users[12]->get_value(mc,ind,mu)) + t20*(*CI_users[13]->get_value(mc,ind,mu)) + co3*(*CI_users[14]->get_value(mc,ind,mu)) + t94*(*CI_users[15]->get_value(mc,ind,mu)) + t76*(*CI_users[16]->get_value(mc,ind,mu)) + t63*(*CI_users[17]->get_value(mc,ind,mu)) + co4*(*CI_users[18]->get_value(mc,ind,mu)) + co5*(*CI_users[19]->get_value(mc,ind,mu)) + co6*(*CI_users[20]->get_value(mc,ind,mu)) + co7*(*CI_users[21]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmqpmp_RT_wCI::\
C2q2g1ph_phdqmqpmp_RT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c45, c12));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c45, c3));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqpmp_RT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, m, p}, RT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqpmp RT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:55:24 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa34 = SPA(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spa23 = SPA(2,3);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa25 = SPA(2,5);
complex<T> s12 = S(1,2);
complex<T> s15 = S(1,5);
complex<T> s23 = -(spa23*spb23);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s13 = S(1,3);
complex<T> s24 = -(spa24*spb24);
complex<T> s234 = SS(2,3,4);
complex<T> s34 = -(spa34*spb34);
complex<T> s245 = SS(2,4,5);
complex<T> s235 = SS(2,3,5);
complex<T> s345 = SS(3,4,5);
complex<T> t4 = spb23*spb45; 
complex<T> t8 = square(spa24*spb34 + spa25*spb35); 
complex<T> t9 = square(spa24*spb23 + spa45*spb35); 
complex<T> t22 = square(spb23); 
complex<T> t23 = square(spb45); 
complex<T> t24 = square(spa24); 
complex<T> t25 = square(spb25); 
complex<T> t35 = cube(spb35); 
complex<T> d2 = (s23 - s234)*spb34*square(spb24); d2 = T(1)/d2;
complex<T> d3 = spb24*spb34*square(s23 - s234)*T(2); d3 = T(1)/d3;
complex<T> d4 = (-s234 + s34)*square(spb24); d4 = T(1)/d4;
complex<T> d5 = spb24*square(s234 - s34)*T(2); d5 = T(1)/d5;
complex<T> d6 = spb24*square(s245 - s25)*T(2); d6 = T(1)/d6;
complex<T> d9 = spa25*square(s245 - s45); d9 = T(1)/d9;
complex<T> d10 = spa25*square(spb24); d10 = T(1)/d10;
complex<T> d11 = spa45*spb34*square(spb24); d11 = T(1)/d11;
complex<T> d12 = spa45*spb34*square(s245 - s25); d12 = T(1)/d12;
complex<T> d14 = spb34*spb45; d14 = T(1)/d14;
complex<T> d15 = spb34*cube(spb24); d15 = T(1)/d15;
complex<T> d16 = spb34*square(spb24); d16 = T(1)/d16;
complex<T> d17 = spb23*spb34; d17 = T(1)/d17;
complex<T> d18 = cube(spb24); d18 = T(1)/d18;
complex<T> d19 = spb23*T(2); d19 = T(1)/d19;
complex<T> d20 = spb23*spb34*T(2); d20 = T(1)/d20;
complex<T> d21 = cube(spb24)*T(2); d21 = T(1)/d21;
complex<T> d22 = spb34*cube(spb24)*T(2); d22 = T(1)/d22;
complex<T> t41 = d17*spa45; 
complex<T> t47 = t22*t23; 
complex<T> t48 = d4*spa24; 
complex<T> t50 = spb34*t25; 
complex<T> t51 = -(d20*spa45); 
complex<T> d1 = spb34*t4*T(2); d1 = T(1)/d1;
complex<T> d7 = spb24*t4*T(2); d7 = T(1)/d7;
complex<T> d8 = spb24*t4*square(s245 - s45)*T(2); d8 = T(1)/d8;
complex<T> d13 = spb34*t4; d13 = T(1)/d13;
complex<T> t1 = -(d10*spa24*t4) + d9*t4*cube(spa24) + d8*t8*cube(spb25) - d7*spb25*square(spb35); 
complex<T> t3 = -(d11*spa24*spb45*t22) + d10*spa24*t4 + d6*t24*t50 + d12*spa24*spb45*t9 - d9*t4*cube(spa24) - d8*t8*cube(spb25) + d7*spb25*square(spb35); 
complex<T> t19 = d2*spa24*t47 + d3*t24*t47 - d1*t35*T(3); 
complex<T> t28 = -t41; 
complex<T> t36 = -t47; 
complex<T> t37 = -t50; 
complex<T> t45 = t35*(-(d1*s245*s345) + d19*spa34*spa45 + (s235 + s24 + s45)*t51); 
complex<T> t49 = -(d13*t35); 
complex<T> t2 = d11*spa24*spb45*t22 + d6*t24*t37 - d12*spa24*spb45*t9; 
complex<T> t14 = d13*mH2*t35 + s15*t49; 
complex<T> t17 = d13*s24*t35 + d16*spa24*t36; 
complex<T> t18 = d13*s15*t35 + d15*s15*t47 + d16*spa24*t47 + s24*t49; 
complex<T> t20 = t37*t48 + d5*t24*t50 - spb25*t4*t48*T(2); 
complex<T> t54 = t28*t35 + s12*t49; 
complex<T> t55 = d2*spa24*t36 + d3*t24*t36 + d5*t24*t37 + t48*t50 + spb25*t4*t48*T(2); 
complex<T> t65 = d15*t36; 
complex<T> t15 = -(d14*spa23*t35) + s23*t65; 
complex<T> t16 = t28*t35 + d15*s13*t47 + s13*t49 + s45*t65; 
complex<T> co1 = s25*t65; 
complex<T> co2 = d18*spa34*t47; 
complex<T> co3 = t35*t41; 
complex<T> co4 = d21*s23*spa34*t47; 
complex<T> co5 = d22*s25*s45*t36; 
SeriesC<T> result = t19*(*CI_users[0]->get_value(mc,ind,mu)) + t55*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu)) + t2*(*CI_users[3]->get_value(mc,ind,mu)) + t20*(*CI_users[4]->get_value(mc,ind,mu)) + t1*(*CI_users[5]->get_value(mc,ind,mu)) + t14*(*CI_users[6]->get_value(mc,ind,mu)) + t15*(*CI_users[7]->get_value(mc,ind,mu)) + t17*(*CI_users[8]->get_value(mc,ind,mu)) + t16*(*CI_users[9]->get_value(mc,ind,mu)) + co1*(*CI_users[10]->get_value(mc,ind,mu)) + t18*(*CI_users[11]->get_value(mc,ind,mu)) + co2*(*CI_users[12]->get_value(mc,ind,mu)) + t54*(*CI_users[13]->get_value(mc,ind,mu)) + co3*(*CI_users[14]->get_value(mc,ind,mu)) + t45*(*CI_users[15]->get_value(mc,ind,mu)) + co4*(*CI_users[16]->get_value(mc,ind,mu)) + co5*(*CI_users[17]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmqpmp_nfLT_wCI::\
C2q2g1ph_phdqmqpmp_nfLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqpmp_nfLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, m, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqpmp nfLT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:55:25 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb45 = SPB(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa34 = SPA(3,4);
complex<T> s23 = S(2,3);
complex<T> s234 = SS(2,3,4);
complex<T> t4 = cube(spa24*spb25 + spa34*spb35); 
complex<T> t7 = cube(spa24); 
complex<T> t8 = square(spb23); 
complex<T> t9 = cube(spb35); 
complex<T> t10 = square(spb45); 
complex<T> d1 = spb23*spb34*spb45*T(3); d1 = T(1)/d1;
complex<T> d2 = spb23*spb45*cube(s23 - s234)*T(3); d2 = T(1)/d2;
complex<T> d3 = spb34*cube(s23 - s234)*T(3); d3 = T(1)/d3;
complex<T> t1 = -(d3*t10*t7*t8) - d1*t9 + d2*t4*square(spb34); 
complex<T> t2 = d3*t10*t7*t8 - d1*t9 - d2*t4*square(spb34); 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C2q2g1ph_phdqmqpmm_LT_wCI::\
C2q2g1ph_phdqmqpmm_LT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1ph_phdqmqpmm_LT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqpmm LT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1ph_phdqmqpmm_RT_wCI::\
C2q2g1ph_phdqmqpmm_RT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1ph_phdqmqpmm_RT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqpmm RT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q2g1ph_phdqmqpmm_nfLT_wCI::\
C2q2g1ph_phdqmqpmm_nfLT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1ph_phdqmqpmm_nfLT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqpmm nfLT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1ph_phdqmqppp_LT_wCI::\
C2q2g1ph_phdqmqppp_LT_wCI
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
CI_users.push_back(new Cached_Triangle_Integral_User(c2, c5, c134));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c14, c25));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c15, c24));
CI_users.push_back(new Cached_Triangle_Integral_User(c3, c4, c125));
CI_users.push_back(new Cached_Triangle_Integral_User(c4, c5, c123));
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c3, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c2, c45));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c3, c25));
CI_users.push_back(new Cached_Box_Integral_User(c1, c4, c5, c23));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c2, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c4, c23));
CI_users.push_back(new Cached_Box_Integral_User(c2, c5, c4, c13));
CI_users.push_back(new Cached_Box_Integral_User(c4, c3, c2, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
CI_users.push_back(new Cached_Box_Integral_User(c5, c4, c3, c12));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqppp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, p, p}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqppp LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:58:27 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb25 = SPB(2,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spa25 = SPA(2,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spa45 = SPA(4,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spa34 = SPA(3,4);
complex<T> spa23 = SPA(2,3);
complex<T> s25 = -(spa25*spb25);
complex<T> s45 = -(spa45*spb45);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> s23 = -(spa23*spb23);
complex<T> s235 = SS(2,3,5);
complex<T> s35 = -(spa35*spb35);
complex<T> s14 = S(1,4);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> s15 = S(1,5);
complex<T> t18 = -(spa25*spb35); 
complex<T> t19 = square(spa24*spb34 + spa25*spb35); 
complex<T> t20 = square(spa23*spb24 + spa35*spb45); 
complex<T> t28 = cube(spa23*spb24 + spa35*spb45); 
complex<T> t29 = spa23*spa25; 
complex<T> t38 = square(mH2); 
complex<T> t39 = square(spb35); 
complex<T> t40 = square(spb34); 
complex<T> t41 = square(s345); 
complex<T> t42 = square(spb45); 
complex<T> t43 = square(spa23); 
complex<T> t44 = spa24*spb23 + spa45*spb35; 
complex<T> t45 = spa35*spb23 + spa45*spb24; 
complex<T> t46 = spa25*spb24 + spa35*spb34; 
complex<T> t47 = spa24*spb25 + spa34*spb35; 
complex<T> t48 = square(spa25); 
complex<T> t50 = s235*T(2); 
complex<T> t51 = square(spa35); 
complex<T> t62 = cube(spb35); 
complex<T> t85 = coeff3mass1234<T>(mc,ind,3,2,5,4); 
complex<T> t86 = coeff3mass4123<T>(mc,ind,3,2,5,4); 
complex<T> t91 = spb35*spb45; 
complex<T> t92 = spa23*spb34; 
complex<T> t117 = cube(s345); 
complex<T> t141 = spa23*spb35; 
complex<T> t142 = spb34*spb45; 
complex<T> t153 = spb24*spb45; 
complex<T> d3 = spa24*spa25*spa45*spb23*T(6); d3 = T(1)/d3;
complex<T> d4 = spa34*square(s24 + s34)*T(2); d4 = T(1)/d4;
complex<T> d5 = (s24 + s34)*spa34*T(2); d5 = T(1)/d5;
complex<T> d6 = spa35*square(s25 + s35)*T(2); d6 = T(1)/d6;
complex<T> d7 = square(s24 + s34)*T(6); d7 = T(1)/d7;
complex<T> d8 = cube(s25 + s35)*T(3); d8 = T(1)/d8;
complex<T> d9 = square(s25 + s35)*T(6); d9 = T(1)/d9;
complex<T> d11 = cube(s24 + s34)*T(3); d11 = T(1)/d11;
complex<T> d12 = square(s24 + s34)*T(3); d12 = T(1)/d12;
complex<T> d13 = square(s25 + s35)*T(3); d13 = T(1)/d13;
complex<T> d14 = s234*spa24*spa34*T(6); d14 = T(1)/d14;
complex<T> d15 = s234*(s24 + s34)*spa34*T(3); d15 = T(1)/d15;
complex<T> d16 = s23*s234*spa24*T(6); d16 = T(1)/d16;
complex<T> d17 = s23*s235*spa25*T(6); d17 = T(1)/d17;
complex<T> d18 = s235*(s25 + s35)*spa25*T(3); d18 = T(1)/d18;
complex<T> d20 = spa35*square(s23 + s35)*T(2); d20 = T(1)/d20;
complex<T> d22 = square(s23 + s35)*T(2); d22 = T(1)/d22;
complex<T> d25 = spa45*square(s24 + s45)*T(2); d25 = T(1)/d25;
complex<T> d26 = (s24 + s45)*spa45*T(2); d26 = T(1)/d26;
complex<T> d42 = spa24*spa45*spb23; d42 = T(1)/d42;
complex<T> d45 = s234*spa24*spa34; d45 = T(1)/d45;
complex<T> d46 = s23*s234*spa24; d46 = T(1)/d46;
complex<T> d48 = s23*s235; d48 = T(1)/d48;
complex<T> d54 = spa24*spa25*spa45*spb23; d54 = T(1)/d54;
complex<T> d55 = s234*spa24; d55 = T(1)/d55;
complex<T> d56 = s23*s235*spa25; d56 = T(1)/d56;
complex<T> d59 = spa24*spa25*spb23; d59 = T(1)/d59;
complex<T> t16 = square(spa24*spb45 + t141); 
complex<T> t17 = square(-(spa25*spb45) + t92); 
complex<T> t22 = cube(-(spa25*spb45) + t92); 
complex<T> t23 = cube(spa24*spb45 + t141); 
complex<T> t27 = s235*t46; 
complex<T> t57 = d11*s234; 
complex<T> t63 = -t142; 
complex<T> t68 = -(spb24*t38); 
complex<T> t77 = s235*t44; 
complex<T> t89 = spb24*t40; 
complex<T> t101 = d25*s245; 
complex<T> t103 = s234*t46; 
complex<T> t122 = spb35*t43; 
complex<T> t126 = spb25*t19; 
complex<T> t127 = s25*t38; 
complex<T> t128 = -(t28*t48); 
complex<T> t131 = t91*t92; 
complex<T> t172 = t29*t91; 
complex<T> t210 = t153*t29; 
complex<T> d1 = (s25 + s35)*spa35*t50; d1 = T(1)/d1;
complex<T> d2 = (s25 + s35)*spa25*t50; d2 = T(1)/d2;
complex<T> d10 = (s25 + s35)*t50; d10 = T(1)/d10;
complex<T> d19 = s235*(s25 + s35)*t51; d19 = T(1)/d19;
complex<T> d21 = (s23 + s35)*spa35*t50; d21 = T(1)/d21;
complex<T> d23 = (s23 + s35)*t50; d23 = T(1)/d23;
complex<T> d24 = s235*(s23 + s35)*t51; d24 = T(1)/d24;
complex<T> d27 = spa34*spa45*t45*T(2); d27 = T(1)/d27;
complex<T> d29 = spa25*spa45*t44*T(2); d29 = T(1)/d29;
complex<T> d31 = s234*spa34*t47; d31 = T(1)/d31;
complex<T> d32 = s234*spa23*t47; d32 = T(1)/d32;
complex<T> d35 = t44*t47*t50; d35 = T(1)/d35;
complex<T> d36 = s234*spa34*t47*T(2); d36 = T(1)/d36;
complex<T> d37 = s234*t47*T(2); d37 = T(1)/d37;
complex<T> d38 = s234*spa24*t47; d38 = T(1)/d38;
complex<T> d43 = spa45*t44; d43 = T(1)/d43;
complex<T> d44 = spb23*t44*t47*t50; d44 = T(1)/d44;
complex<T> d47 = s234*spa23*t47*T(2); d47 = T(1)/d47;
complex<T> d52 = s234*spa23*spa24*t47; d52 = T(1)/d52;
complex<T> d53 = spa45*t45; d53 = T(1)/d53;
complex<T> d58 = spa34*t45*T(2); d58 = T(1)/d58;
complex<T> d60 = spa25*t44*T(2); d60 = T(1)/d60;
complex<T> d61 = spb23*t44*t47*T(2); d61 = T(1)/d61;
complex<T> d62 = t29*t46*T(2); d62 = T(1)/d62;
complex<T> d63 = spb23*t45*t46*T(2); d63 = T(1)/d63;
complex<T> d64 = spa34*t47*T(2); d64 = T(1)/d64;
complex<T> d65 = spa23*t47*T(2); d65 = T(1)/d65;
complex<T> d66 = t44*T(2); d66 = T(1)/d66;
complex<T> d67 = s234*spa24*t47*T(2); d67 = T(1)/d67;
complex<T> d69 = t45*T(2); d69 = T(1)/d69;
complex<T> t9 = -(s45*(d62*t22 + d61*t38*t62)); 
complex<T> t36 = -(s34*(d62*t22 + d61*t38*t62)); 
complex<T> t65 = -(spb25*t16); 
complex<T> t90 = spa24*t101*t142 + d26*(spa24*t142 + spa25*t91)*T(3); 
complex<T> t95 = -(d32*t16); 
complex<T> t96 = d29*t19; 
complex<T> t100 = d31*s15; 
complex<T> t111 = t40*t68; 
complex<T> t114 = d32*s15; 
complex<T> t118 = d18*t17; 
complex<T> t121 = spb45*t16; 
complex<T> t130 = d56*t17; 
complex<T> t143 = -(d46*t16); 
complex<T> t160 = spb34*t122; 
complex<T> t183 = d23*spb35; 
complex<T> t191 = spb34*t57; 
complex<T> d28 = spb23*t103*t45; d28 = T(1)/d28;
complex<T> d30 = spb23*t47*t77; d30 = T(1)/d30;
complex<T> d33 = t27*t29; d33 = T(1)/d33;
complex<T> d34 = t103*t45*T(2); d34 = T(1)/d34;
complex<T> d39 = spa25*t27*T(2); d39 = T(1)/d39;
complex<T> d40 = t27*cube(spa35); d40 = T(1)/d40;
complex<T> d41 = spb23*t103*t45*T(2); d41 = T(1)/d41;
complex<T> d49 = spa23*t27*T(2); d49 = T(1)/d49;
complex<T> d50 = spa23*t27; d50 = T(1)/d50;
complex<T> d51 = spa23*t27*cube(spa35); d51 = T(1)/d51;
complex<T> d57 = t27*t29*T(2); d57 = T(1)/d57;
complex<T> d68 = t27*cube(spa35)*T(2); d68 = T(1)/d68;
complex<T> t1 = d7*t131 - d12*spa24*spb34*t42 + d4*spa24*spb25*spb45*t92 + spa24*t191*t42*T(2) - d5*spa24*t42*T(3) - d5*spa23*t91*T(3) - d15*spb24*t16*T(11); 
complex<T> t2 = s25*(d63*t111 + d65*t121 + d64*t65); 
complex<T> t4 = d28*s15*t111 + t114*t121 + d31*mH2*spb25*t16 + d33*mH2*t22 - d33*s15*t22 - d27*mH2*t41 + d27*s15*t41 - d30*s15*t38*t62 + t100*t65 + mH2*spb45*t95 - mH2*t96 + s15*t96 + d30*t62*cube(mH2) + d28*t89*cube(mH2); 
complex<T> t5 = t114*t121 + d31*s24*spb25*t16 + d33*s15*t22 - d33*s24*t22 - d52*s15*t23 - d32*spb24*t23 - d27*s15*t41 + d27*s24*t41 + d30*s15*t38*t62 - d30*s24*t38*t62 + t100*t65 + s24*spb45*t95 - s15*t96 + s24*t96; 
complex<T> t6 = d32*(s24*t121 + spb24*t23) + s24*(d33*t22 - d27*t41 + d30*t38*t62 + d31*t65 - t96); 
complex<T> t7 = d41*s34*t111 - d47*s34*t121 + s34*spb35*t130 + s34*spb34*t143 - d55*spb34*t16 - d54*s34*t19 - d57*s34*t22 + d52*s34*t23 - d53*spb34*t41 - d44*s34*t38*t62 + d37*spb34*t65; 
complex<T> t8 = d41*s45*t111 + d47*s45*t121 + s45*spb35*t130 + s45*spb34*t143 + d45*s45*t16 + d59*spb45*t19 - d60*spb45*t19 - d57*s45*t22 - d58*spb45*t41 - d44*s45*t38*t62 + d36*s45*t65; 
complex<T> t10 = d41*s25*t111 + d47*s25*t121 + d42*t126 - d43*t126 + s25*spb34*t143 + d45*s25*t16 - d48*spb25*spb35*t17 + d36*s25*t65 + d49*spb25*t22*T(3) - d44*t127*t62*T(3); 
complex<T> t11 = d10*t131 + d9*t131 - d1*spb24*t160 + d1*spb24*t172 + d6*spb24*t172 - d2*t122*t40 + d12*spa24*spb34*t42 + d4*spa23*spa24*spb25*t63 + d7*t141*t63 + d8*s35*t131*T(2) - d8*t210*t39*T(2) - spa24*t191*t42*T(2) + d5*spa24*t42*T(3) + d5*spa23*t91*T(3) + d15*spb24*t16*T(11) + d14*t16*T(13) - d16*spb34*t16*T(13) - d3*t19*T(13) + spb35*(d19*spa25*t20 + d13*spa25*t42 + t118*T(8) + d17*t17*T(13)); 
complex<T> t12 = d1*spb24*t160 - d20*spb24*t160 + d1*spa23*t153*t18 + d6*spa23*t153*t18 + d19*t18*t20 + d24*t18*t20 + d2*t122*t40 + d13*t18*t42 + d10*t141*t63 + d22*t141*t63 + d9*t141*t63 - d8*s35*t131*T(2) + d8*t210*t39*T(2) + d23*t131*T(3) + d21*spb24*t160*T(3) - d21*spb24*t172*T(3) - spa25*t183*t42*T(3) - spb35*t118*T(8); 
complex<T> t13 = d22*t131 + d20*spb24*t160 + d24*spa25*spb35*t20 + spa24*t101*t63 - d23*t131*T(3) - d26*spa24*t142*T(3) - d21*spb24*t160*T(3) + d21*spb24*t172*T(3) + spa25*t183*t42*T(3) - d26*spa25*t91*T(3); 
complex<T> t14 = d51*s14*t128 + d51*s25*t28*t48 + d30*t127*t62 - d30*s14*t38*t62 - d33*s14*t22*T(2) - d50*spb25*t22*T(2); 
complex<T> t35 = s45*(d63*t111 + d65*t121 + d64*t65); 
complex<T> t37 = d68*s25*spb23*t128 + d35*spa23*t127*t62; 
complex<T> t139 = d34*t38; 
complex<T> t140 = s23*t96; 
complex<T> t189 = spa23*t139; 
complex<T> t3 = -(d67*s34*spb23*t23) + s34*t189*t89; 
complex<T> t15 = d37*spb23*t121 + d40*spb23*t128 + t140 + d36*s23*spb25*t16 - d39*spb23*t22 - d38*spb23*t23 + d27*s23*t41 + d35*spa23*t38*t62 + t189*t89; 
complex<T> co1 = d27*s23*t117; 
complex<T> co2 = s245*t140; 
complex<T> co3 = d66*spb45*t126; 
complex<T> co4 = d69*t142*t41; 
SeriesC<T> result = t11*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t12*(*CI_users[2]->get_value(mc,ind,mu)) + t90*(*CI_users[3]->get_value(mc,ind,mu)) + t13*(*CI_users[4]->get_value(mc,ind,mu)) + t85*(*CI_users[5]->get_value(mc,ind,mu)) + t86*(*CI_users[6]->get_value(mc,ind,mu)) + t4*(*CI_users[7]->get_value(mc,ind,mu)) + t15*(*CI_users[8]->get_value(mc,ind,mu)) + t6*(*CI_users[9]->get_value(mc,ind,mu)) + t10*(*CI_users[10]->get_value(mc,ind,mu)) + t14*(*CI_users[11]->get_value(mc,ind,mu)) + t5*(*CI_users[12]->get_value(mc,ind,mu)) + t7*(*CI_users[13]->get_value(mc,ind,mu)) + t8*(*CI_users[14]->get_value(mc,ind,mu)) + co1*(*CI_users[15]->get_value(mc,ind,mu)) + co2*(*CI_users[16]->get_value(mc,ind,mu)) + t36*(*CI_users[17]->get_value(mc,ind,mu)) + t9*(*CI_users[18]->get_value(mc,ind,mu)) + t2*(*CI_users[19]->get_value(mc,ind,mu)) + t35*(*CI_users[20]->get_value(mc,ind,mu)) + co3*(*CI_users[21]->get_value(mc,ind,mu)) + t3*(*CI_users[22]->get_value(mc,ind,mu)) + t37*(*CI_users[23]->get_value(mc,ind,mu)) + co4*(*CI_users[24]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmqppp_RT_wCI::\
C2q2g1ph_phdqmqppp_RT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c2, c5, c34));
CI_users.push_back(new Cached_Box_Integral_User(c1, c3, c4, c25));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqppp_RT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, p, p}, RT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqppp RT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:59:02 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spa35 = SPA(3,5);
complex<T> spb24 = SPB(2,4);
complex<T> s13 = S(1,3);
complex<T> s24 = -(spa24*spb24);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s23 = -(spa23*spb23);
complex<T> s14 = S(1,4);
complex<T> s235 = SS(2,3,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s35 = -(spa35*spb35);
complex<T> s245 = SS(2,4,5);
complex<T> s45 = -(spa45*spb45);
complex<T> s345 = SS(3,4,5);
complex<T> s15 = S(1,5);
complex<T> s12 = S(1,2);
complex<T> t15 = square(spa23); 
complex<T> t16 = square(spa24*spb34 + spa25*spb35); 
complex<T> t17 = square(spa23*spb35 + spa24*spb45); 
complex<T> t18 = square(spa25*spb24 + spa35*spb34); 
complex<T> t19 = square(spb34); 
complex<T> t20 = square(spa23*spb24 + spa35*spb45); 
complex<T> t21 = square(spa24*spb25 + spa34*spb35); 
complex<T> t22 = square(spa23*spb34 - spa25*spb45); 
complex<T> t23 = square(s23); 
complex<T> t24 = square(spb35); 
complex<T> t25 = spa24*T(2); 
complex<T> t33 = -(spb35*T(2)); 
complex<T> t34 = square(s345); 
complex<T> t36 = square(spb24); 
complex<T> t37 = spa35*spb23 + spa45*spb24; 
complex<T> t38 = square(s245); 
complex<T> t39 = square(spa24); 
complex<T> t41 = square(spb45); 
complex<T> t53 = spb35*T(2); 
complex<T> t69 = coeff3mass4123<T>(mc,ind,3,2,5,4); 
complex<T> t70 = spa23*spb45; 
complex<T> t71 = spb24*spb34; 
complex<T> t85 = cube(s345); 
complex<T> d1 = s235*square(spa35); d1 = T(1)/d1;
complex<T> d2 = s235*spa25*spa35; d2 = T(1)/d2;
complex<T> d3 = s235*spb25*square(spa35); d3 = T(1)/d3;
complex<T> d4 = s235*s25*spa35; d4 = T(1)/d4;
complex<T> d5 = s235*spa25*spa35*square(s25 + s35)*T(2); d5 = T(1)/d5;
complex<T> d8 = (s25 + s35)*square(spa35); d8 = T(1)/d8;
complex<T> d9 = s235*spa35; d9 = T(1)/d9;
complex<T> d10 = (s25 + s35)*spa35; d10 = T(1)/d10;
complex<T> d11 = (s24 + s34)*spa34; d11 = T(1)/d11;
complex<T> d12 = s235*spa35*spb25; d12 = T(1)/d12;
complex<T> d13 = (s25 + s35)*spa35*spb25; d13 = T(1)/d13;
complex<T> d14 = s235*s25; d14 = T(1)/d14;
complex<T> d15 = s25*(s25 + s35); d15 = T(1)/d15;
complex<T> d16 = s234*spa24*spa34; d16 = T(1)/d16;
complex<T> d18 = s23*s235*spa25*T(2); d18 = T(1)/d18;
complex<T> d20 = (s23 + s35)*spb25*square(spa35); d20 = T(1)/d20;
complex<T> d21 = s25*(s23 + s35)*spa35; d21 = T(1)/d21;
complex<T> d22 = (s23 + s35)*spa35*spb25; d22 = T(1)/d22;
complex<T> d23 = s25*(s23 + s35); d23 = T(1)/d23;
complex<T> d24 = s235*cube(spa35)*T(2); d24 = T(1)/d24;
complex<T> d25 = s235*cube(spa35)*square(s23 + s35)*T(2); d25 = T(1)/d25;
complex<T> d26 = spa25*spa45*(spa24*spb23 + spa45*spb35)*square(s24 + s45)*T(2); d26 = T(1)/d26;
complex<T> d27 = spa25*spa45*(spa24*spb23 + spa45*spb35)*T(2); d27 = T(1)/d27;
complex<T> d28 = (s24 + s45)*spa45; d28 = T(1)/d28;
complex<T> d30 = s235*spa25*cube(spa35); d30 = T(1)/d30;
complex<T> d31 = spa24*spa25*spa45; d31 = T(1)/d31;
complex<T> d32 = s234*spa24; d32 = T(1)/d32;
complex<T> d33 = s235*spa25; d33 = T(1)/d33;
complex<T> d34 = s234*spa34; d34 = T(1)/d34;
complex<T> d35 = spa25*spa45*(spa24*spb23 + spa45*spb35); d35 = T(1)/d35;
complex<T> d36 = spa25*(spa24*spb23 + spa45*spb35); d36 = T(1)/d36;
complex<T> d37 = spa45*(spa24*spb23 + spa45*spb35)*T(2); d37 = T(1)/d37;
complex<T> d38 = s235*cube(spa35); d38 = T(1)/d38;
complex<T> t6 = (d35*s13 + d36*spb45)*t16; 
complex<T> t50 = d20*t15; 
complex<T> t52 = -t70; 
complex<T> t72 = d27*t16; 
complex<T> t82 = spb24*t70; 
complex<T> t84 = spb34*t17; 
complex<T> t89 = spb34*t70; 
complex<T> t95 = t15*t71; 
complex<T> t96 = d28*spb45; 
complex<T> t98 = spb25*t18; 
complex<T> t108 = d3*t36; 
complex<T> t114 = s23*t15; 
complex<T> d6 = spa25*spa45*spb23*t25; d6 = T(1)/d6;
complex<T> d7 = s234*spa34*t25*square(s24 + s34); d7 = T(1)/d7;
complex<T> d17 = s23*s234*t25; d17 = T(1)/d17;
complex<T> d19 = s234*spa34*t25; d19 = T(1)/d19;
complex<T> d29 = spa34*spa45*t37*T(2); d29 = T(1)/d29;
complex<T> d39 = spa45*t37*T(2); d39 = T(1)/d39;
complex<T> d40 = spa34*spa45*t37; d40 = T(1)/d40;
complex<T> d41 = spa34*t37; d41 = T(1)/d41;
complex<T> d42 = s234*t25; d42 = T(1)/d42;
complex<T> t1 = t15*(d30*s14*t18 + d38*t98); 
complex<T> t2 = -((mH2 - s15)*(d29*t34 + t72)); 
complex<T> t3 = d37*spb25*t16 + d29*s25*t34; 
complex<T> t4 = -(d26*t19*t38*t39) + t72 + spb34*t25*t96 + spa25*t53*t96; 
complex<T> t5 = -(spb45*(d36*t16 + d41*t34)); 
complex<T> t7 = d19*t17 - d7*t15*t19*t21 + d11*t33*t70 - d11*spa24*t41*T(2); 
complex<T> t8 = -(d34*spb24*t17) - d29*s24*t34 - s24*t72; 
complex<T> t9 = d16*s15*t17 + d34*spb24*t17 - (s15 - s24)*(d29*t34 + t72); 
complex<T> t11 = -(d2*t15*t19) + d7*t15*t19*t21 + d5*t15*t18*t24 + d11*t25*t41 + d8*spa25*spb24*t52 + d10*spb34*t52 + t108*t15*t53 + d11*t53*t70 + d1*spa25*t82 + d13*t33*t82 + d12*t53*t82 + d9*t89 + d14*t33*t89 + d15*t53*t89 - d1*t95 + d4*t33*t95 - d16*t17*T(2) + d6*t16*T(3) - d18*spb35*t22*T(3) + d17*t84*T(3); 
complex<T> t12 = d31*spa23*t16 - d30*t114*t18 + d33*spb35*t22 - d32*t84; 
complex<T> t13 = d2*t15*t19 + d24*spa25*t20 - d25*spa25*t20*t23 - d5*t15*t18*t24 + t33*t36*t50 + d1*spa25*spb24*t52 + d9*spb34*t52 + d8*spa25*t82 + d22*t33*t82 + d13*t53*t82 + d10*t89 + d15*t33*t89 + d23*t53*t89 + d1*t95 + d21*t53*t95; 
complex<T> t14 = -(d24*spa25*t20) + d25*spa25*t20*t23 + t108*t15*t33 + d26*t19*t38*t39 + t36*t50*t53 - t72 + d12*t33*t82 + d22*t53*t82 + d23*t33*t89 + d14*t53*t89 + d21*t33*t95 + d4*t53*t95 + spa25*t33*t96 - spa24*spb34*t96*T(2); 
complex<T> t32 = (d40*s12 + d41*spb45)*t34; 
complex<T> t101 = s34*t72; 
complex<T> t10 = t101 + d39*spb34*t34 + d32*t84; 
complex<T> co1 = d29*s25*t85; 
complex<T> co2 = s245*t101; 
complex<T> co3 = d42*s23*t84; 
complex<T> co4 = d24*t114*t98; 
SeriesC<T> result = t11*(*CI_users[0]->get_value(mc,ind,mu)) + t7*(*CI_users[1]->get_value(mc,ind,mu)) + t13*(*CI_users[2]->get_value(mc,ind,mu)) + t4*(*CI_users[3]->get_value(mc,ind,mu)) + t14*(*CI_users[4]->get_value(mc,ind,mu)) + t69*(*CI_users[5]->get_value(mc,ind,mu)) + t2*(*CI_users[6]->get_value(mc,ind,mu)) + t12*(*CI_users[7]->get_value(mc,ind,mu)) + t8*(*CI_users[8]->get_value(mc,ind,mu)) + t6*(*CI_users[9]->get_value(mc,ind,mu)) + t3*(*CI_users[10]->get_value(mc,ind,mu)) + t1*(*CI_users[11]->get_value(mc,ind,mu)) + t9*(*CI_users[12]->get_value(mc,ind,mu)) + t10*(*CI_users[13]->get_value(mc,ind,mu)) + t32*(*CI_users[14]->get_value(mc,ind,mu)) + t5*(*CI_users[15]->get_value(mc,ind,mu)) + co1*(*CI_users[16]->get_value(mc,ind,mu)) + co2*(*CI_users[17]->get_value(mc,ind,mu)) + co3*(*CI_users[18]->get_value(mc,ind,mu)) + co4*(*CI_users[19]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmqppp_nfLT_wCI::\
C2q2g1ph_phdqmqppp_nfLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmqppp_nfLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, qp, p, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmqppp nfLT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:59:06 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa24 = SPA(2,4);
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa34 = SPA(3,4);
complex<T> s24 = -(spa24*spb24);
complex<T> s34 = -(spa34*spb34);
complex<T> s35 = S(3,5);
complex<T> s25 = -(spa25*spb25);
complex<T> s23 = -(spa23*spb23);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> t4 = square(spa23*spb35 + spa24*spb45); 
complex<T> t5 = -(spb34*T(2)); 
complex<T> t6 = square(spa23*spb34 - spa25*spb45); 
complex<T> t10 = spa23*spb45; 
complex<T> t13 = s23*T(3); 
complex<T> t24 = spb34*T(2); 
complex<T> t30 = square(spb45); 
complex<T> t35 = square(spb34); 
complex<T> t38 = square(spb35); 
complex<T> t42 = spa24*spb25; 
complex<T> d1 = spa24*spa25*spa45*spb23*T(3); d1 = T(1)/d1;
complex<T> d2 = cube(s24 + s34)*T(3); d2 = T(1)/d2;
complex<T> d3 = square(s24 + s34)*T(3); d3 = T(1)/d3;
complex<T> d4 = cube(s25 + s35)*T(3); d4 = T(1)/d4;
complex<T> d5 = square(s25 + s35)*T(3); d5 = T(1)/d5;
complex<T> d8 = s234*spa24*spa34*T(3); d8 = T(1)/d8;
complex<T> t17 = d4*spb24; 
complex<T> t34 = spb35*t10; 
complex<T> t50 = t10*t35; 
complex<T> d6 = (s24 + s34)*spa24*t13; d6 = T(1)/d6;
complex<T> d7 = (s25 + s35)*spa25*t13; d7 = T(1)/d7;
complex<T> d9 = s234*spa24*t13; d9 = T(1)/d9;
complex<T> d10 = s235*spa25*t13; d10 = T(1)/d10;
complex<T> t1 = d2*s34*t24*t34 - d3*spb34*(spa24*t30 + t34) + d9*t24*t4 + d6*t4*t5 - d8*t4*T(2) - d2*t42*t50*T(2); 
complex<T> t29 = d7*t6; 
complex<T> t40 = t10*t17; 
complex<T> t2 = d3*spa24*spb34*t30 + d3*spb34*t34 + d5*spb34*t34 + d6*t24*t4 + d2*s34*t34*t5 + d4*s35*t34*t5 + spa25*t38*t40*T(2) + d2*t42*t50*T(2) + d1*square(spa24*spb34 + spa25*spb35)*T(2) - spb35*(d5*spa25*t30 + t29*T(2)); 
complex<T> t3 = -(d5*spb34*t34) + d4*s35*t24*t34 - spa25*t38*t40*T(2) + spb35*(d5*spa25*t30 + t29*T(2) - d10*t6*T(2)); 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t3*(*CI_users[2]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmpqpm_LT_wCI::\
C2q2g1ph_phdqmpqpm_LT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c23, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c2, c3, c4, c15));
CI_users.push_back(new Cached_Box_Integral_User(c4, c5, c2, c13));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmpqpm_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, p, qp, m}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmpqpm LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:59:09 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa25 = SPA(2,5);
complex<T> spa45 = SPA(4,5);
complex<T> spb34 = SPB(3,4);
complex<T> spb25 = SPB(2,5);
complex<T> spb45 = SPB(4,5);
complex<T> spb24 = SPB(2,4);
complex<T> spa35 = SPA(3,5);
complex<T> spb35 = SPB(3,5);
complex<T> spb23 = SPB(2,3);
complex<T> s12 = S(1,2);
complex<T> s13 = S(1,3);
complex<T> s14 = S(1,4);
complex<T> s15 = S(1,5);
complex<T> s24 = S(2,4);
complex<T> s23 = S(2,3);
complex<T> s34 = S(3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s234 = SS(2,3,4);
complex<T> s45 = -(spa45*spb45);
complex<T> s245 = SS(2,4,5);
complex<T> s345 = SS(3,4,5);
complex<T> t2 = square(spa25*spb24 + spa35*spb34); 
complex<T> t5 = square(spb35); 
complex<T> t10 = square(spb34); 
complex<T> t13 = square(spa25); 
complex<T> t14 = square(spb23); 
complex<T> t21 = -(spa25*T(2)); 
complex<T> t29 = spa25*T(2); 
complex<T> d1 = (s23 - s235)*spb25; d1 = T(1)/d1;
complex<T> d2 = spb25*spb45*T(2); d2 = T(1)/d2;
complex<T> d3 = spb25*spb45*square(s23 - s235)*T(2); d3 = T(1)/d3;
complex<T> d4 = (-s245 + s45)*spb25; d4 = T(1)/d4;
complex<T> d5 = spb25*square(s245 - s45)*T(2); d5 = T(1)/d5;
complex<T> d6 = spb25*spb45; d6 = T(1)/d6;
complex<T> d7 = spb25; d7 = T(1)/d7;
complex<T> d8 = spb45; d8 = T(1)/d8;
complex<T> d9 = T(2); d9 = T(1)/d9;
complex<T> t7 = -(d5*spb45*t13*t14) + d4*spb23*spb34*t21; 
complex<T> t16 = -(d2*mH2); 
complex<T> t22 = d2*s234; 
complex<T> t23 = d7*spa45; 
complex<T> t28 = d6*t10; 
complex<T> t35 = d8*spa25; 
complex<T> t37 = d2*t10; 
complex<T> t1 = d1*spb23*spb34*t29 + t37 - d3*t2*t5; 
complex<T> t6 = (-mH2 + s15)*t28*T(2); 
complex<T> t8 = d5*spb45*t13*t14 + d4*spb23*spb34*t29 + t37*T(3); 
complex<T> t9 = d1*spb23*spb34*t21 - t37 + d3*t2*t5; 
complex<T> t17 = -t23; 
complex<T> t20 = -t28; 
complex<T> t36 = t10*(s23*t16 + s235*t22); 
complex<T> t39 = t10*(s34*t16 + s345*t22); 
complex<T> t41 = t10*t23 + s12*t28; 
complex<T> t46 = t10*t35; 
complex<T> t27 = s14*t28 + t46; 
complex<T> t32 = t10*t17 + s13*t20; 
complex<T> t43 = s15*t20 + s24*t28; 
complex<T> co1 = s24*t20; 
complex<T> co2 = -t46; 
complex<T> co3 = s23*s34*t37; 
complex<T> co4 = d9*spa25*spa45*t10; 
SeriesC<T> result = t9*(*CI_users[0]->get_value(mc,ind,mu)) + t1*(*CI_users[1]->get_value(mc,ind,mu)) + t8*(*CI_users[2]->get_value(mc,ind,mu)) + t7*(*CI_users[3]->get_value(mc,ind,mu)) + t6*(*CI_users[4]->get_value(mc,ind,mu)) + co1*(*CI_users[5]->get_value(mc,ind,mu)) + t32*(*CI_users[6]->get_value(mc,ind,mu)) + co2*(*CI_users[7]->get_value(mc,ind,mu)) + t27*(*CI_users[8]->get_value(mc,ind,mu)) + t43*(*CI_users[9]->get_value(mc,ind,mu)) + t41*(*CI_users[10]->get_value(mc,ind,mu)) + t36*(*CI_users[11]->get_value(mc,ind,mu)) + t39*(*CI_users[12]->get_value(mc,ind,mu)) + co3*(*CI_users[13]->get_value(mc,ind,mu)) + co4*(*CI_users[14]->get_value(mc,ind,mu));  
 return(result);
} 
  
  


C2q2g1ph_phdqmpqpm_nfLT_wCI::\
C2q2g1ph_phdqmpqpm_nfLT_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q2g1ph_phdqmpqpm_nfLT_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmpqpm nfLT");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 

C2q2g1ph_phdqmmqpp_LT_wCI::\
C2q2g1ph_phdqmmqpp_LT_wCI
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
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c23, c4));
CI_users.push_back(new Cached_Box_Integral_User(c1, c5, c34, c2));
CI_users.push_back(new Cached_Box_Integral_User(c3, c4, c5, c12));
CI_users.push_back(new Cached_Box_Integral_User(c5, c2, c3, c14));
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmmqpp_LT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, m, qp, p}, LT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmmqpp LT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:59:10 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa34 = SPA(3,4);
complex<T> spb23 = SPB(2,3);
complex<T> spb45 = SPB(4,5);
complex<T> spa23 = SPA(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb24 = SPB(2,4);
complex<T> spb35 = SPB(3,5);
complex<T> spa24 = SPA(2,4);
complex<T> spb25 = SPB(2,5);
complex<T> s45 = S(4,5);
complex<T> s25 = S(2,5);
complex<T> s15 = S(1,5);
complex<T> s34 = -(spa34*spb34);
complex<T> s234 = SS(2,3,4);
complex<T> s235 = SS(2,3,5);
complex<T> s345 = SS(3,4,5);
complex<T> t3 = square(spa23*spb35 + spa24*spb45); 
complex<T> t9 = square(spb45); 
complex<T> t12 = square(spa23); 
complex<T> t13 = square(spb25); 
complex<T> t20 = spb25*spb45; 
complex<T> d1 = (-s234 + s34)*spb23; d1 = T(1)/d1;
complex<T> d2 = spb23*spb34*T(3); d2 = T(1)/d2;
complex<T> d3 = spb23*spb34*square(s234 - s34)*T(2); d3 = T(1)/d3;
complex<T> d4 = spb23*square(s235 - s25)*T(2); d4 = T(1)/d4;
complex<T> d5 = (-s235 + s25)*spb23; d5 = T(1)/d5;
complex<T> d6 = spb23*spb34*T(2); d6 = T(1)/d6;
complex<T> d7 = spb23*spb34; d7 = T(1)/d7;
complex<T> d8 = spb34*T(2); d8 = T(1)/d8;
complex<T> d9 = spb23*T(2); d9 = T(1)/d9;
complex<T> t1 = d6*t9 - d3*t3*square(spb24) + d1*spa23*t20*T(2); 
complex<T> t2 = d3*t3*square(spb24) - d1*spa23*t20*T(2) + d2*t9*T(5); 
complex<T> t8 = d7*(mH2 - s15)*t9*T(2); 
complex<T> t15 = -(d6*s234); 
complex<T> t21 = d9*spa34; 
complex<T> t22 = d4*spb34; 
complex<T> t32 = d8*spa23; 
complex<T> t7 = -(t12*t13*t22) - d5*spa23*t20*T(2); 
complex<T> t16 = -t21; 
complex<T> t24 = t12*t13*t22 + d5*spa23*t20*T(2); 
complex<T> t25 = (s235*t15 - mH2*t32)*t9; 
complex<T> t29 = (s345*t15 + mH2*t16)*t9; 
complex<T> co1 = s45*t21*t9; 
complex<T> co2 = s25*t32*t9; 
SeriesC<T> result = t2*(*CI_users[0]->get_value(mc,ind,mu)) + t7*(*CI_users[1]->get_value(mc,ind,mu)) + t24*(*CI_users[2]->get_value(mc,ind,mu)) + t1*(*CI_users[3]->get_value(mc,ind,mu)) + t8*(*CI_users[4]->get_value(mc,ind,mu)) + t25*(*CI_users[5]->get_value(mc,ind,mu)) + t29*(*CI_users[6]->get_value(mc,ind,mu)) + co1*(*CI_users[7]->get_value(mc,ind,mu)) + co2*(*CI_users[8]->get_value(mc,ind,mu));  
 return(result);
} 
  
  

C2q2g1ph_phdqmmqpp_nfLT_wCI::\
C2q2g1ph_phdqmmqpp_nfLT_wCI
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
} 
  
  
template <class T> SeriesC<T> 
     C2q2g1ph_phdqmmqpp_nfLT_wCI::eval_fn
        (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
//{{phd, qm, m, qp, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("C2q2g1ph :  phdqmmqpp nfLT");
#endif
 
//#define TimeStamp "Sat 11 Dec 2010 20:59:11 on n2175"
//#define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spb23 = SPB(2,3);
complex<T> spb34 = SPB(3,4);
complex<T> spb45 = SPB(4,5);
complex<T> t1 = square(spb45); 
complex<T> d1 = spb23*spb34*T(3); d1 = T(1)/d1;
complex<T> co1 = -(d1*t1*T(2)); 
SeriesC<T> result = co1*(*CI_users[0]->get_value(mc,ind,mu));  
 return(result);
} 
  
  
 
 
 // *************** table of switch values ************* 
 
#define _C_phqmqppm_LT C2q2g1ph_598800_LT
#define _C_phqmqppm_RT C2q2g1ph_598800_RT
#define _C_phqmqppm_nfLT C2q2g1ph_598800_nfLT
#define _C_phqmqpmp_LT C2q2g1ph_598785_LT
#define _C_phqmqpmp_RT C2q2g1ph_598785_RT
#define _C_phqmqpmp_nfLT C2q2g1ph_598785_nfLT
#define _C_phqmqpmm_LT C2q2g1ph_598784_LT
#define _C_phqmqpmm_RT C2q2g1ph_598784_RT
#define _C_phqmqpmm_nfLT C2q2g1ph_598784_nfLT
#define _C_phqmqppp_LT C2q2g1ph_598801_LT
#define _C_phqmqppp_RT C2q2g1ph_598801_RT
#define _C_phqmqppp_nfLT C2q2g1ph_598801_nfLT
#define _C_phqmpqpm_LT C2q2g1ph_598320_LT
#define _C_phqmpqpm_nfLT C2q2g1ph_598320_nfLT
#define _C_phqmmqpp_LT C2q2g1ph_598065_LT
#define _C_phqmmqpp_nfLT C2q2g1ph_598065_nfLT
#define _C_phdqmqppm_LT C2q2g1ph_664336_LT
#define _C_phdqmqppm_RT C2q2g1ph_664336_RT
#define _C_phdqmqppm_nfLT C2q2g1ph_664336_nfLT
#define _C_phdqmqpmp_LT C2q2g1ph_664321_LT
#define _C_phdqmqpmp_RT C2q2g1ph_664321_RT
#define _C_phdqmqpmp_nfLT C2q2g1ph_664321_nfLT
#define _C_phdqmqpmm_LT C2q2g1ph_664320_LT
#define _C_phdqmqpmm_RT C2q2g1ph_664320_RT
#define _C_phdqmqpmm_nfLT C2q2g1ph_664320_nfLT
#define _C_phdqmqppp_LT C2q2g1ph_664337_LT
#define _C_phdqmqppp_RT C2q2g1ph_664337_RT
#define _C_phdqmqppp_nfLT C2q2g1ph_664337_nfLT
#define _C_phdqmpqpm_LT C2q2g1ph_663856_LT
#define _C_phdqmpqpm_nfLT C2q2g1ph_663856_nfLT
#define _C_phdqmmqpp_LT C2q2g1ph_663601_LT
#define _C_phdqmmqpp_nfLT C2q2g1ph_663601_nfLT
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_phqmqppm_LT case 598800
 
#define _CASE_phqmqppm_RT case 598800
 
#define _CASE_phqmqppm_nfLT case 598800
 
#define _CASE_phqmqpmp_LT case 598785
 
#define _CASE_phqmqpmp_RT case 598785
 
#define _CASE_phqmqpmp_nfLT case 598785
 
#define _CASE_phqmqpmm_LT case 598784
 
#define _CASE_phqmqpmm_RT case 598784
 
#define _CASE_phqmqpmm_nfLT case 598784
 
#define _CASE_phqmqppp_LT case 598801
 
#define _CASE_phqmqppp_RT case 598801
 
#define _CASE_phqmqppp_nfLT case 598801
 
#define _CASE_phqmpqpm_LT case 598320
 
#define _CASE_phqmpqpm_nfLT case 598320
 
#define _CASE_phqmmqpp_LT case 598065
 
#define _CASE_phqmmqpp_nfLT case 598065
 
#define _CASE_phdqmqppm_LT case 664336
 
#define _CASE_phdqmqppm_RT case 664336
 
#define _CASE_phdqmqppm_nfLT case 664336
 
#define _CASE_phdqmqpmp_LT case 664321
 
#define _CASE_phdqmqpmp_RT case 664321
 
#define _CASE_phdqmqpmp_nfLT case 664321
 
#define _CASE_phdqmqpmm_LT case 664320
 
#define _CASE_phdqmqpmm_RT case 664320
 
#define _CASE_phdqmqpmm_nfLT case 664320
 
#define _CASE_phdqmqppp_LT case 664337
 
#define _CASE_phdqmqppp_RT case 664337
 
#define _CASE_phdqmqppp_nfLT case 664337
 
#define _CASE_phdqmpqpm_LT case 663856
 
#define _CASE_phdqmpqpm_nfLT case 663856
 
#define _CASE_phdqmmqpp_LT case 663601
 
#define _CASE_phdqmmqpp_nfLT case 663601
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q2g1ph_LT( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_phqmqppm_LT: return new 
                       C2q2g1ph_phqmqppm_LT_wCI(ind);
    _CASE_phqmqpmp_LT: return new 
                       C2q2g1ph_phqmqpmp_LT_wCI(ind);
    _CASE_phqmqpmm_LT: return new 
                       C2q2g1ph_phqmqpmm_LT_wCI(ind);
    _CASE_phqmqppp_LT: return new 
                       C2q2g1ph_phqmqppp_LT_wCI(ind);
    _CASE_phqmpqpm_LT: return new 
                       C2q2g1ph_phqmpqpm_LT_wCI(ind);
    _CASE_phqmmqpp_LT: return new 
                       C2q2g1ph_phqmmqpp_LT_wCI(ind);
    _CASE_phdqmqppm_LT: return new 
                       C2q2g1ph_phdqmqppm_LT_wCI(ind);
    _CASE_phdqmqpmp_LT: return new 
                       C2q2g1ph_phdqmqpmp_LT_wCI(ind);
    _CASE_phdqmqpmm_LT: return new 
                       C2q2g1ph_phdqmqpmm_LT_wCI(ind);
    _CASE_phdqmqppp_LT: return new 
                       C2q2g1ph_phdqmqppp_LT_wCI(ind);
    _CASE_phdqmpqpm_LT: return new 
                       C2q2g1ph_phdqmpqpm_LT_wCI(ind);
    _CASE_phdqmmqpp_LT: return new 
                       C2q2g1ph_phdqmmqpp_LT_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g1ph_nfLT( int hc,const std::vector<int>& ind) 
{ 
    switch (hc) {
    _CASE_phqmqppm_nfLT: return new 
                       C2q2g1ph_phqmqppm_nfLT_wCI(ind);
    _CASE_phqmqpmp_nfLT: return new 
                       C2q2g1ph_phqmqpmp_nfLT_wCI(ind);
    _CASE_phqmqpmm_nfLT: return new 
                       C2q2g1ph_phqmqpmm_nfLT_wCI(ind);
    _CASE_phqmqppp_nfLT: return new 
                       C2q2g1ph_phqmqppp_nfLT_wCI(ind);
    _CASE_phqmpqpm_nfLT: return new 
                       C2q2g1ph_phqmpqpm_nfLT_wCI(ind);
    _CASE_phqmmqpp_nfLT: return new 
                       C2q2g1ph_phqmmqpp_nfLT_wCI(ind);
    _CASE_phdqmqppm_nfLT: return new 
                       C2q2g1ph_phdqmqppm_nfLT_wCI(ind);
    _CASE_phdqmqpmp_nfLT: return new 
                       C2q2g1ph_phdqmqpmp_nfLT_wCI(ind);
    _CASE_phdqmqpmm_nfLT: return new 
                       C2q2g1ph_phdqmqpmm_nfLT_wCI(ind);
    _CASE_phdqmqppp_nfLT: return new 
                       C2q2g1ph_phdqmqppp_nfLT_wCI(ind);
    _CASE_phdqmpqpm_nfLT: return new 
                       C2q2g1ph_phdqmpqpm_nfLT_wCI(ind);
    _CASE_phdqmmqpp_nfLT: return new 
                       C2q2g1ph_phdqmmqpp_nfLT_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q2g1ph_RT( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_phqmqppm_RT: return new 
                       C2q2g1ph_phqmqppm_RT_wCI(ind);
    _CASE_phqmqpmp_RT: return new 
                       C2q2g1ph_phqmqpmp_RT_wCI(ind);
    _CASE_phqmqpmm_RT: return new 
                       C2q2g1ph_phqmqpmm_RT_wCI(ind);
    _CASE_phqmqppp_RT: return new 
                       C2q2g1ph_phqmqppp_RT_wCI(ind);
    _CASE_phdqmqppm_RT: return new 
                       C2q2g1ph_phdqmqppm_RT_wCI(ind);
    _CASE_phdqmqpmp_RT: return new 
                       C2q2g1ph_phdqmqpmp_RT_wCI(ind);
    _CASE_phdqmqpmm_RT: return new 
                       C2q2g1ph_phdqmqpmm_RT_wCI(ind);
    _CASE_phdqmqppp_RT: return new 
                       C2q2g1ph_phdqmqppp_RT_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
