/*
* R_4g1ph_eval.cpp
*
* Created on 12/8, 2010
*      Author: Zvi's script
*/
 
#include "R_4g1ph_eval.h"
#include "eval_param.h"
 
using namespace std;
 
namespace BH  {
 
 
#define _VERBOSE 1
#define mH2 ep.s(1,2,3,4)
 
#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}

 
 
template <class T> complex<T> R4g1ph_phpppp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpppp G");
#endif
 
 return( (complex<T>(2,0)*pow(mH2,2))/(SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpppp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R4g1ph_phmppp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmppp G");
#endif
 
 return(-((pow(SPA(2,4),2)*SPA(2,3)*SPB(4,3))/(complex<T>(3,0)*pow(SPA(3,4),2)*
SPA(2,5)*SPA(4,5)))-(pow(SPA(2,4),2)*SPA(2,5)*SPB(5,4))/
(complex<T>(3,0)*pow(SPA(4,5),2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(-2,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),3))/
(SPA(2,3)*SPA(3,4)*(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*
 SS(2,3,4))+(complex<T>(2,0)*pow(mH2,2)*pow(SPB(5,3),4))/
(SPB(3,2)*SPB(5,2)*(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*
 (SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3))*SS(2,3,5))+
 (complex<T>(2,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),3))/
(SPA(2,5)*SPA(4,5)*(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3))*
 SS(2,4,5))-(complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5))/
(complex<T>(3,0)*SPA(3,4)*SPA(4,5)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3))/
(complex<T>(3,0)*SPA(3,4)*SPA(4,5)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(4,3)*SPB(5,4))/
(complex<T>(3,0)*SPA(3,4)*SPA(4,5)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPB(5,3)*SPB(5,4))/
(complex<T>(3,0)*SPA(3,4)*SPA(4,5)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmppp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(2,3)*SPB(4,3))/
(complex<T>(9,0)*pow(SPA(3,4),2)*SPA(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(2,5)*SPB(5,4))/
(complex<T>(9,0)*pow(SPA(4,5),2)*SPA(2,3)*SPA(3,4))+
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5))/(complex<T>(9,0)*SPA(3,4)*SPA(4,5)*
 SS(3,4,5))+(SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3))/
(complex<T>(9,0)*SPA(3,4)*SPA(4,5)*SS(3,4,5))+
 (pow(SPA(2,4),2)*SPB(4,3)*SPB(5,4))/(complex<T>(9,0)*SPA(3,4)*SPA(4,5)*
 SS(3,4,5))+(SPA(2,4)*SPA(2,5)*SPB(5,3)*SPB(5,4))/
(complex<T>(9,0)*SPA(3,4)*SPA(4,5)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpmpp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpmpp G");
#endif
 
 return(-((pow(SPA(3,5),2)*SPA(2,3)*SPB(5,2))/(complex<T>(3,0)*pow(SPA(2,5),2)*
SPA(3,4)*SPA(4,5)))-(pow(SPA(3,5),2)*SPA(3,4)*SPB(5,4))/
(complex<T>(3,0)*pow(SPA(4,5),2)*SPA(2,3)*SPA(2,5))+
 (complex<T>(-2,0)*pow(mH2,2)*pow(SPB(4,2),4))/
(SPB(3,2)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(4,3)*
 (SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4),3))/
(SPA(2,3)*SPA(2,5)*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*
 SS(2,3,5))-(pow(SPB(4,2),2)*SPA(2,3)*SPA(3,4))/
(complex<T>(3,0)*SPA(2,5)*SPA(4,5)*SS(2,4,5))-
 (SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(5,2))/(complex<T>(3,0)*SPA(2,5)*SPA(4,5)*
 SS(2,4,5))-(SPA(3,4)*SPA(3,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(3,0)*SPA(2,5)*SPA(4,5)*SS(2,4,5))-
 (pow(SPA(3,5),2)*SPB(5,2)*SPB(5,4))/(complex<T>(3,0)*SPA(2,5)*SPA(4,5)*
 SS(2,4,5))+(complex<T>(-2,0)*pow(SPA(3,4)*SPB(4,2)+
 SPA(3,5)*SPB(5,2),3))/(SPA(3,4)*SPA(4,5)*
 (-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpmpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpmpp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(2,3)*SPB(5,2))/
(complex<T>(9,0)*pow(SPA(2,5),2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(3,4)*SPB(5,4))/
(complex<T>(9,0)*pow(SPA(4,5),2)*SPA(2,3)*SPA(2,5))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(3,4))/
(complex<T>(9,0)*SPA(2,5)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(9,0)*SPA(2,5)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(9,0)*SPA(2,5)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(5,2)*SPB(5,4))/
(complex<T>(9,0)*SPA(2,5)*SPA(4,5)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phppmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phppmp G");
#endif
 
 return(-((pow(SPA(2,4),2)*SPA(3,4)*SPB(3,2))/(complex<T>(3,0)*pow(SPA(2,3),2)*
SPA(2,5)*SPA(4,5)))-(pow(SPA(2,4),2)*SPA(4,5)*SPB(5,2))/
(complex<T>(3,0)*pow(SPA(2,5),2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(-2,0)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),3))/
(SPA(2,3)*SPA(3,4)*(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4))*
 SS(2,3,4))-(pow(SPB(5,3),2)*SPA(3,4)*SPA(4,5))/
(complex<T>(3,0)*SPA(2,3)*SPA(2,5)*SS(2,3,5))-
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(5,2))/(complex<T>(3,0)*SPA(2,3)*SPA(2,5)*
 SS(2,3,5))-(SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,3))/
(complex<T>(3,0)*SPA(2,3)*SPA(2,5)*SS(2,3,5))-
 (SPA(2,4)*SPA(4,5)*SPB(5,2)*SPB(5,3))/(complex<T>(3,0)*SPA(2,3)*SPA(2,5)*
 SS(2,3,5))+(complex<T>(2,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),
3))/(SPA(2,5)*SPA(4,5)*(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3))*
 SS(2,4,5))+(complex<T>(-2,0)*pow(mH2,2)*pow(SPB(5,3),4))/
(SPB(4,3)*(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3))*SPB(5,4)*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4))*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phppmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phppmp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(3,4)*SPB(3,2))/
(complex<T>(9,0)*pow(SPA(2,3),2)*SPA(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(4,5)*SPB(5,2))/
(complex<T>(9,0)*pow(SPA(2,5),2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(3,4)*SPA(4,5))/
(complex<T>(9,0)*SPA(2,3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(5,2))/
(complex<T>(9,0)*SPA(2,3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,3))/
(complex<T>(9,0)*SPA(2,3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(4,5)*SPB(5,2)*SPB(5,3))/
(complex<T>(9,0)*SPA(2,3)*SPA(2,5)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpppm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpppm G");
#endif
 
 return(-((pow(SPA(3,5),2)*SPA(2,5)*SPB(3,2))/(complex<T>(3,0)*pow(SPA(2,3),2)*
SPA(3,4)*SPA(4,5)))-(pow(SPA(3,5),2)*SPA(4,5)*SPB(4,3))/
(complex<T>(3,0)*pow(SPA(3,4),2)*SPA(2,3)*SPA(2,5))-
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,5)*SPA(4,5))/
(complex<T>(3,0)*SPA(2,3)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,5)*SPA(3,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(3,0)*SPA(2,3)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3))/
(complex<T>(3,0)*SPA(2,3)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(3,0)*SPA(2,3)*SPA(3,4)*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3))/
(SPA(2,3)*SPA(2,5)*(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4))*
 SS(2,3,5))+(complex<T>(2,0)*pow(mH2,2)*pow(SPB(4,2),4))/
(SPB(5,2)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*SPB(5,4)*
 (SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4))*SS(2,4,5))+
 (complex<T>(-2,0)*pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),3))/
(SPA(3,4)*SPA(4,5)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*
 SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpppm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(2,5)*SPB(3,2))/
(complex<T>(9,0)*pow(SPA(2,3),2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(4,5)*SPB(4,3))/
(complex<T>(9,0)*pow(SPA(3,4),2)*SPA(2,3)*SPA(2,5))+
 (pow(SPB(4,2),2)*SPA(2,5)*SPA(4,5))/(complex<T>(9,0)*SPA(2,3)*SPA(3,4)*
 SS(2,3,4))+(SPA(2,5)*SPA(3,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(9,0)*SPA(2,3)*SPA(3,4)*SS(2,3,4))+
 (pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3))/(complex<T>(9,0)*SPA(2,3)*SPA(3,4)*
 SS(2,3,4))+(SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(9,0)*SPA(2,3)*SPA(3,4)*SS(2,3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmmpp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmmpp G");
#endif
 
 return(-(pow(SPA(2,3),2)/(complex<T>(3,0)*pow(SPA(4,5),2)))+
 (complex<T>(2,0)*pow(SPA(2,3),3))/(SPA(2,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*S(2,4)*SPA(2,3)*SPA(2,5)*SPA(3,4))/
(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*S(3,5)*SPA(2,3)*SPA(2,5))/
(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*SPB(4,3))+
 (complex<T>(-2,0)*pow(SPB(5,4),3))/(SPB(3,2)*SPB(4,3)*SPB(5,2))-
 (pow(SPA(2,3),2)*SPB(5,4))/(complex<T>(6,0)*S(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(S(2,4)*SPA(2,3)+SPA(2,4)*SPA(3,5)*SPB(5,4),2)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(4,5))-
 (pow(SPA(2,3),2)*SPB(5,4))/(complex<T>(6,0)*S(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(4,5))-
 (SPA(2,3)*SPB(5,4))/(complex<T>(3,0)*SPA(4,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(S(3,5)*SPA(2,3)+SPA(2,4)*SPA(3,5)*SPB(5,4),2)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*SPA(3,4)*SPA(4,5)*
 SPB(4,3))+(complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,4)*SPA(2,5)*SPA(3,5)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*SPB(4,3))-
 (S(2,4)*SPA(2,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(4,5),2)*SPB(3,2)*
 SPB(4,3))-(S(3,4)*SPA(2,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(4,5),2)*
 SPB(3,2)*SPB(4,3))-(S(2,5)*SPA(3,4)*SPB(5,4))/
(complex<T>(3,0)*pow(SPA(4,5),2)*SPB(3,2)*SPB(5,2))-
 (S(3,5)*SPA(3,4)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(4,5),2)*SPB(3,2)*
 SPB(5,2))+(complex<T>(1,0)*pow(SPB(5,4),2)*S(2,4)*SPA(2,3)*SPA(2,5)*
 SPA(3,4))/(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*SPA(2,5))/
(complex<T>(3,0)*pow(SPA(4,5),2)*SPB(5,2)*SS(2,4,5))-
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,3)*SPA(3,4))/
(complex<T>(6,0)*SPA(4,5)*SPB(5,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(5,4))/(complex<T>(6,0)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(S(2,4)*SPA(2,3)+SPA(2,4)*SPA(3,5)*SPB(5,4),2)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*S(3,5)*SPA(2,3)*SPA(2,5)*SPA(3,4))/
(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*SS(3,4,5))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPA(3,4))/
(complex<T>(3,0)*pow(SPA(4,5),2)*SPB(4,3)*SS(3,4,5))-
 (pow(SPB(5,4),2)*SPA(2,3)*SPA(2,5))/(complex<T>(6,0)*SPA(4,5)*SPB(4,3)*
 SS(3,4,5))+(pow(SPA(2,3),2)*SPB(5,4))/(complex<T>(6,0)*SPA(4,5)*
 SS(3,4,5))+(complex<T>(1,0)*pow(S(3,5)*SPA(2,3)+SPA(2,4)*SPA(3,5)*
  SPB(5,4),2)*SPB(5,4))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*
 SPA(4,5)*SS(3,4,5))+(complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,4)*
 SPA(2,5)*SPA(3,4)*SPA(3,5)*SPB(5,4))/
(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmmpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmmpp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,3),2))/(complex<T>(9,0)*pow(SPA(4,5),2))-
 (pow(SPB(5,4),2)*S(2,4)*SPA(2,3)*SPA(2,5)*SPA(3,4))/
(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(4,5))-
 (pow(SPB(5,4),2)*S(3,5)*SPA(2,3)*SPA(2,5))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(5,4))/(complex<T>(18,0)*S(2,5)*SPA(4,5))-
 (pow(S(2,4)*SPA(2,3)+SPA(2,4)*SPA(3,5)*SPB(5,4),2)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(5,4))/(complex<T>(18,0)*S(3,4)*SPA(4,5))-
 (pow(SPB(5,4),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(4,5))+
 (complex<T>(1,0)*SPA(2,3)*SPB(5,4))/(complex<T>(9,0)*SPA(4,5)*SPB(3,2))-
 (pow(S(3,5)*SPA(2,3)+SPA(2,4)*SPA(3,5)*SPB(5,4),2)*SPB(5,4))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPA(3,4)*SPA(4,5)*SPB(4,3))-
 (pow(SPB(5,4),2)*SPA(2,4)*SPA(2,5)*SPA(3,5)*SPB(5,4))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*SPB(4,3))+
 (complex<T>(1,0)*S(2,4)*SPA(2,5)*SPB(5,4))/(complex<T>(9,0)*pow(SPA(4,5),2)*
 SPB(3,2)*SPB(4,3))+(complex<T>(1,0)*S(3,4)*SPA(2,5)*SPB(5,4))/
(complex<T>(9,0)*pow(SPA(4,5),2)*SPB(3,2)*SPB(4,3))+
 (complex<T>(1,0)*S(2,5)*SPA(3,4)*SPB(5,4))/(complex<T>(9,0)*pow(SPA(4,5),2)*
 SPB(3,2)*SPB(5,2))+(complex<T>(1,0)*S(3,5)*SPA(3,4)*SPB(5,4))/
(complex<T>(9,0)*pow(SPA(4,5),2)*SPB(3,2)*SPB(5,2))-
 (pow(SPB(5,4),2)*S(2,4)*SPA(2,3)*SPA(2,5)*SPA(3,4))/
(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPA(4,5)*SS(2,4,5))-
 (pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*SPA(2,5))/
(complex<T>(9,0)*pow(SPA(4,5),2)*SPB(5,2)*SS(2,4,5))+
 (pow(SPB(5,4),2)*SPA(2,3)*SPA(3,4))/(complex<T>(18,0)*SPA(4,5)*SPB(5,2)*
 SS(2,4,5))-(pow(SPA(2,3),2)*SPB(5,4))/(complex<T>(18,0)*SPA(4,5)*
 SS(2,4,5))-(pow(S(2,4)*SPA(2,3)+SPA(2,4)*SPA(3,5)*SPB(5,4),
2)*SPB(5,4))/(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPA(4,5)*
 SS(2,4,5))-(pow(SPB(5,4),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*
 SPA(3,5)*SPB(5,4))/(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPA(4,5)*
 SS(2,4,5))-(pow(SPB(5,4),2)*S(3,5)*SPA(2,3)*SPA(2,5)*
 SPA(3,4))/(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*
 SS(3,4,5))-(pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*
 SPA(3,4))/(complex<T>(9,0)*pow(SPA(4,5),2)*SPB(4,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,3)*SPA(2,5))/
(complex<T>(18,0)*SPA(4,5)*SPB(4,3)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(5,4))/(complex<T>(18,0)*SPA(4,5)*SS(3,4,5))-
 (pow(S(3,5)*SPA(2,3)+SPA(2,4)*SPA(3,5)*SPB(5,4),2)*SPB(5,4))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*SS(3,4,5))-
 (pow(SPB(5,4),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*SPB(5,4))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpmmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpmmp G");
#endif
 
 return(-(pow(SPA(3,4),2)/(complex<T>(3,0)*pow(SPA(2,5),2)))+
 (complex<T>(2,0)*pow(SPA(3,4),3))/(SPA(2,3)*SPA(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*S(3,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SPB(3,2))-
 (pow(SPA(3,4),2)*SPB(5,2))/(complex<T>(6,0)*S(2,3)*SPA(2,5))-
 (pow(SPA(3,4),2)*SPB(5,2))/(complex<T>(6,0)*S(4,5)*SPA(2,5))+
 (complex<T>(1,0)*pow(S(3,5)*SPA(3,4)+SPA(2,4)*SPA(3,5)*SPB(5,2),2)*
 SPB(5,2))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*SPA(2,3)*SPA(2,5)*
 SPB(3,2))+(complex<T>(1,0)*pow(SPB(5,2),2)*SPA(2,4)*SPA(3,5)*SPA(4,5)*
 SPB(5,2))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SPB(3,2))-
 (SPA(3,4)*SPB(5,2))/(complex<T>(3,0)*SPA(2,5)*SPB(4,3))-
 (S(2,3)*SPA(4,5)*SPB(5,2))/(complex<T>(3,0)*pow(SPA(2,5),2)*SPB(3,2)*
 SPB(4,3))-(S(2,4)*SPA(4,5)*SPB(5,2))/(complex<T>(3,0)*pow(SPA(2,5),2)*
 SPB(3,2)*SPB(4,3))+(complex<T>(1,0)*pow(SPB(5,2),2)*S(2,4)*SPA(2,3)*
 SPA(3,4))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SPB(5,4))+
 (complex<T>(-2,0)*pow(SPB(5,2),3))/(SPB(3,2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPB(5,2))/
(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SPB(5,4))+
 (complex<T>(1,0)*pow(S(2,4)*SPA(3,4)+SPA(2,4)*SPA(3,5)*SPB(5,2),2)*
 SPB(5,2))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SPA(4,5)*
 SPB(5,4))-(S(3,5)*SPA(2,3)*SPB(5,2))/(complex<T>(3,0)*pow(SPA(2,5),2)*
 SPB(4,3)*SPB(5,4))-(S(4,5)*SPA(2,3)*SPB(5,2))/
(complex<T>(3,0)*pow(SPA(2,5),2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*S(3,5)*SPA(2,3)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*SPA(2,3))/
(complex<T>(3,0)*pow(SPA(2,5),2)*SPB(3,2)*SS(2,3,5))-
 (pow(SPB(5,2),2)*SPA(3,4)*SPA(4,5))/(complex<T>(6,0)*SPA(2,5)*SPB(3,2)*
 SS(2,3,5))+(pow(SPA(3,4),2)*SPB(5,2))/(complex<T>(6,0)*SPA(2,5)*
 SS(2,3,5))+(complex<T>(1,0)*pow(S(3,5)*SPA(3,4)+SPA(2,4)*SPA(3,5)*
  SPB(5,2),2)*SPB(5,2))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*
 SPA(2,5)*SS(2,3,5))+(complex<T>(1,0)*pow(SPB(5,2),2)*SPA(2,3)*
 SPA(2,4)*SPA(3,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*S(2,4)*SPA(2,3)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SS(2,4,5))+
 (pow(SPA(3,4),2)*SPB(5,2))/(complex<T>(6,0)*SPA(2,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(S(2,4)*SPA(3,4)+SPA(2,4)*SPA(3,5)*SPB(5,2),2)*
 SPB(5,2))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPA(4,5)*
 SPB(5,2))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SS(2,4,5))-
 (pow(SPB(5,2),2)*SPA(2,3)*SPA(3,4))/(complex<T>(6,0)*SPA(2,5)*SPB(5,4)*
 SS(2,4,5))+(complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),
2)*SPA(4,5))/(complex<T>(3,0)*pow(SPA(2,5),2)*SPB(5,4)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpmmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpmmp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(3,4),2))/(complex<T>(9,0)*pow(SPA(2,5),2))-
 (pow(SPB(5,2),2)*S(3,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPB(5,2))/(complex<T>(18,0)*S(2,3)*SPA(2,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPB(5,2))/(complex<T>(18,0)*S(4,5)*SPA(2,5))-
 (pow(S(3,5)*SPA(3,4)+SPA(2,4)*SPA(3,5)*SPB(5,2),2)*SPB(5,2))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPA(2,3)*SPA(2,5)*SPB(3,2))-
 (pow(SPB(5,2),2)*SPA(2,4)*SPA(3,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SPB(3,2))+
 (complex<T>(1,0)*SPA(3,4)*SPB(5,2))/(complex<T>(9,0)*SPA(2,5)*SPB(4,3))+
 (complex<T>(1,0)*S(2,3)*SPA(4,5)*SPB(5,2))/(complex<T>(9,0)*pow(SPA(2,5),2)*
 SPB(3,2)*SPB(4,3))+(complex<T>(1,0)*S(2,4)*SPA(4,5)*SPB(5,2))/
(complex<T>(9,0)*pow(SPA(2,5),2)*SPB(3,2)*SPB(4,3))-
 (pow(SPB(5,2),2)*S(2,4)*SPA(2,3)*SPA(3,4))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SPB(5,4))-
 (pow(SPB(5,2),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPB(5,2))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SPB(5,4))-
 (pow(S(2,4)*SPA(3,4)+SPA(2,4)*SPA(3,5)*SPB(5,2),2)*SPB(5,2))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SPA(4,5)*SPB(5,4))+
 (complex<T>(1,0)*S(3,5)*SPA(2,3)*SPB(5,2))/(complex<T>(9,0)*pow(SPA(2,5),2)*
 SPB(4,3)*SPB(5,4))+(complex<T>(1,0)*S(4,5)*SPA(2,3)*SPB(5,2))/
(complex<T>(9,0)*pow(SPA(2,5),2)*SPB(4,3)*SPB(5,4))-
 (pow(SPB(5,2),2)*S(3,5)*SPA(2,3)*SPA(3,4)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SS(2,3,5))-
 (pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*SPA(2,3))/
(complex<T>(9,0)*pow(SPA(2,5),2)*SPB(3,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*SPA(3,4)*SPA(4,5))/
(complex<T>(18,0)*SPA(2,5)*SPB(3,2)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPB(5,2))/(complex<T>(18,0)*SPA(2,5)*SS(2,3,5))-
 (pow(S(3,5)*SPA(3,4)+SPA(2,4)*SPA(3,5)*SPB(5,2),2)*SPB(5,2))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SS(2,3,5))-
 (pow(SPB(5,2),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5)*SS(2,3,5))-
 (pow(SPB(5,2),2)*S(2,4)*SPA(2,3)*SPA(3,4)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SS(2,4,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPB(5,2))/(complex<T>(18,0)*SPA(2,5)*SS(2,4,5))-
 (pow(S(2,4)*SPA(3,4)+SPA(2,4)*SPA(3,5)*SPB(5,2),2)*SPB(5,2))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SS(2,4,5))-
 (pow(SPB(5,2),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*SPA(2,3)*SPA(3,4))/
(complex<T>(18,0)*SPA(2,5)*SPB(5,4)*SS(2,4,5))-
 (pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*SPA(4,5))/
(complex<T>(9,0)*pow(SPA(2,5),2)*SPB(5,4)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phppmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phppmm G");
#endif
 
 return(-(pow(SPA(4,5),2)/(complex<T>(3,0)*pow(SPA(2,3),2)))+
 (complex<T>(2,0)*pow(SPA(4,5),3))/(SPA(2,3)*SPA(2,5)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*S(3,5)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(2,3))-
 (pow(SPA(4,5),2)*SPB(3,2))/(complex<T>(6,0)*S(2,5)*SPA(2,3))+
 (complex<T>(1,0)*pow(S(3,5)*SPA(4,5)+SPA(2,4)*SPA(3,5)*SPB(3,2),2)*
 SPB(3,2))/(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(2,3))-
 (pow(SPA(4,5),2)*SPB(3,2))/(complex<T>(6,0)*S(3,4)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*
 SPB(3,2))/(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*S(2,4)*SPA(2,5)*SPA(4,5))/
(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SPB(4,3))+
 (complex<T>(1,0)*pow(S(2,4)*SPA(4,5)+SPA(2,4)*SPA(3,5)*SPB(3,2),2)*
 SPB(3,2))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SPA(3,4)*
 SPB(4,3))+(complex<T>(1,0)*pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,5)*
 SPB(3,2))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SPB(4,3))-
 (SPA(4,5)*SPB(3,2))/(complex<T>(3,0)*SPA(2,3)*SPB(5,4))-
 (S(3,4)*SPA(2,5)*SPB(3,2))/(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(4,3)*
 SPB(5,4))-(S(3,5)*SPA(2,5)*SPB(3,2))/(complex<T>(3,0)*pow(SPA(2,3),2)*
 SPB(4,3)*SPB(5,4))-(S(2,4)*SPA(3,4)*SPB(3,2))/
(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(5,2)*SPB(5,4))-
 (S(2,5)*SPA(3,4)*SPB(3,2))/(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(5,2)*
 SPB(5,4))+(complex<T>(-2,0)*pow(SPB(3,2),3))/(SPB(4,3)*SPB(5,2)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPB(3,2),2)*S(2,4)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SS(2,3,4))+
 (pow(SPA(4,5),2)*SPB(3,2))/(complex<T>(6,0)*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(S(2,4)*SPA(4,5)+SPA(2,4)*SPA(3,5)*SPB(3,2),2)*
 SPB(3,2))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*
 SPB(3,2))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPA(3,4))/
(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(4,3)*SS(2,3,4))-
 (pow(SPB(3,2),2)*SPA(2,5)*SPA(4,5))/(complex<T>(6,0)*SPA(2,3)*SPB(4,3)*
 SS(2,3,4))+(complex<T>(1,0)*pow(SPB(3,2),2)*S(3,5)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))/(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPA(2,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPB(3,2))/(complex<T>(6,0)*SPA(2,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(S(3,5)*SPA(4,5)+SPA(2,4)*SPA(3,5)*SPB(3,2),2)*
 SPB(3,2))/(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPA(2,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*
 SPB(3,2))/(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPA(2,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*SPA(2,5))/
(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(5,2)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*SPA(2,3)*SPB(5,2)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phppmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phppmm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(4,5),2))/(complex<T>(9,0)*pow(SPA(2,3),2))-
 (pow(SPB(3,2),2)*S(3,5)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPB(3,2))/(complex<T>(18,0)*S(2,5)*SPA(2,3))-
 (pow(S(3,5)*SPA(4,5)+SPA(2,4)*SPA(3,5)*SPB(3,2),2)*SPB(3,2))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPB(3,2))/(complex<T>(18,0)*S(3,4)*SPA(2,3))-
 (pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*SPB(3,2))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(2,3))-
 (pow(SPB(3,2),2)*S(2,4)*SPA(2,5)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SPB(4,3))-
 (pow(S(2,4)*SPA(4,5)+SPA(2,4)*SPA(3,5)*SPB(3,2),2)*SPB(3,2))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SPA(3,4)*SPB(4,3))-
 (pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,5)*SPB(3,2))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SPB(4,3))+
 (complex<T>(1,0)*SPA(4,5)*SPB(3,2))/(complex<T>(9,0)*SPA(2,3)*SPB(5,4))+
 (complex<T>(1,0)*S(3,4)*SPA(2,5)*SPB(3,2))/(complex<T>(9,0)*pow(SPA(2,3),2)*
 SPB(4,3)*SPB(5,4))+(complex<T>(1,0)*S(3,5)*SPA(2,5)*SPB(3,2))/
(complex<T>(9,0)*pow(SPA(2,3),2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(1,0)*S(2,4)*SPA(3,4)*SPB(3,2))/(complex<T>(9,0)*pow(SPA(2,3),2)*
 SPB(5,2)*SPB(5,4))+(complex<T>(1,0)*S(2,5)*SPA(3,4)*SPB(3,2))/
(complex<T>(9,0)*pow(SPA(2,3),2)*SPB(5,2)*SPB(5,4))-
 (pow(SPB(3,2),2)*S(2,4)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPB(3,2))/(complex<T>(18,0)*SPA(2,3)*SS(2,3,4))-
 (pow(S(2,4)*SPA(4,5)+SPA(2,4)*SPA(3,5)*SPB(3,2),2)*SPB(3,2))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SS(2,3,4))-
 (pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*SPB(3,2))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3)*SS(2,3,4))-
 (pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPA(3,4))/
(complex<T>(9,0)*pow(SPA(2,3),2)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(2,5)*SPA(4,5))/
(complex<T>(18,0)*SPA(2,3)*SPB(4,3)*SS(2,3,4))-
 (pow(SPB(3,2),2)*S(3,5)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPA(2,3)*SS(2,3,5))-
 (pow(SPA(4,5),2)*SPB(3,2))/(complex<T>(18,0)*SPA(2,3)*SS(2,3,5))-
 (pow(S(3,5)*SPA(4,5)+SPA(2,4)*SPA(3,5)*SPB(3,2),2)*SPB(3,2))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPA(2,3)*SS(2,3,5))-
 (pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,4)*SPA(3,5)*SPB(3,2))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPA(2,3)*SS(2,3,5))-
 (pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*SPA(2,5))/
(complex<T>(9,0)*pow(SPA(2,3),2)*SPB(5,2)*SS(2,3,5))+
 (pow(SPB(3,2),2)*SPA(3,4)*SPA(4,5))/(complex<T>(18,0)*SPA(2,3)*SPB(5,2)*
 SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmppm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmppm G");
#endif
 
 return(-(pow(SPA(2,5),2)/(complex<T>(3,0)*pow(SPA(3,4),2)))+
 (complex<T>(2,0)*pow(SPA(2,5),3))/(SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (pow(SPB(4,3),2)*S(2,4)*SPA(2,5)*SPA(4,5))/
(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SPB(3,2))-
 (pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(6,0)*S(2,3)*SPA(3,4))-
 (pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(6,0)*S(4,5)*SPA(3,4))+
 (complex<T>(1,0)*pow(-(S(2,4)*SPA(2,5))-SPA(2,4)*SPA(3,5)*SPB(4,3),2)*
 SPB(4,3))/(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPA(2,3)*SPA(3,4)*
 SPB(3,2))+(pow(SPB(4,3),2)*SPA(2,4)*SPA(3,5)*SPA(4,5)*
 SPB(4,3))/(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SPB(3,2))-
 (SPA(2,5)*SPB(4,3))/(complex<T>(3,0)*SPA(3,4)*SPB(5,2))-
 (S(2,3)*SPA(4,5)*SPB(4,3))/(complex<T>(3,0)*pow(SPA(3,4),2)*SPB(3,2)*
 SPB(5,2))-(S(3,5)*SPA(4,5)*SPB(4,3))/(complex<T>(3,0)*pow(SPA(3,4),2)*
 SPB(3,2)*SPB(5,2))+(pow(SPB(4,3),2)*S(3,5)*SPA(2,3)*SPA(2,5))/
(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SPB(5,4))+
 (pow(SPB(4,3),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPB(4,3))/
(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SPB(5,4))+
 (complex<T>(1,0)*pow(-(S(3,5)*SPA(2,5))-SPA(2,4)*SPA(3,5)*SPB(4,3),2)*
 SPB(4,3))/(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SPA(4,5)*
 SPB(5,4))+(complex<T>(-2,0)*pow(SPB(4,3),3))/(SPB(3,2)*SPB(5,2)*
 SPB(5,4))-(S(2,4)*SPA(2,3)*SPB(4,3))/(complex<T>(3,0)*pow(SPA(3,4),2)*
 SPB(5,2)*SPB(5,4))-(S(4,5)*SPA(2,3)*SPB(4,3))/
(complex<T>(3,0)*pow(SPA(3,4),2)*SPB(5,2)*SPB(5,4))+
 (pow(SPB(4,3),2)*S(2,4)*SPA(2,3)*SPA(2,5)*SPA(4,5))/
(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPA(2,3))/
(complex<T>(3,0)*pow(SPA(3,4),2)*SPB(3,2)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,5)*SPA(4,5))/
(complex<T>(6,0)*SPA(3,4)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(6,0)*SPA(3,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(-(S(2,4)*SPA(2,5))-SPA(2,4)*SPA(3,5)*SPB(4,3),2)*
 SPB(4,3))/(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SS(2,3,4))+
 (pow(SPB(4,3),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPA(4,5)*SPB(4,3))/
(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SS(2,3,4))+
 (pow(SPB(4,3),2)*S(3,5)*SPA(2,3)*SPA(2,5)*SPA(4,5))/
(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(6,0)*SPA(3,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(-(S(3,5)*SPA(2,5))-SPA(2,4)*SPA(3,5)*SPB(4,3),2)*
 SPB(4,3))/(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SS(3,4,5))+
 (pow(SPB(4,3),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPA(4,5)*SPB(4,3))/
(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,3)*SPA(2,5))/
(complex<T>(6,0)*SPA(3,4)*SPB(5,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPA(4,5))/
(complex<T>(3,0)*pow(SPA(3,4),2)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmppm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,5),2))/(complex<T>(9,0)*pow(SPA(3,4),2))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*S(2,4)*SPA(2,5)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(18,0)*S(2,3)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(18,0)*S(4,5)*SPA(3,4))-
 (pow(-(S(2,4)*SPA(2,5))-SPA(2,4)*SPA(3,5)*SPB(4,3),2)*SPB(4,3))/
(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*SPA(2,3)*SPA(3,4)*SPB(3,2))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,4)*SPA(3,5)*SPA(4,5)*SPB(4,3))/
(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SPB(3,2))+
 (complex<T>(1,0)*SPA(2,5)*SPB(4,3))/(complex<T>(9,0)*SPA(3,4)*SPB(5,2))+
 (complex<T>(1,0)*S(2,3)*SPA(4,5)*SPB(4,3))/(complex<T>(9,0)*pow(SPA(3,4),2)*
 SPB(3,2)*SPB(5,2))+(complex<T>(1,0)*S(3,5)*SPA(4,5)*SPB(4,3))/
(complex<T>(9,0)*pow(SPA(3,4),2)*SPB(3,2)*SPB(5,2))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*S(3,5)*SPA(2,3)*SPA(2,5))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SPB(5,4))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPB(4,3))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SPB(5,4))-
 (pow(-(S(3,5)*SPA(2,5))-SPA(2,4)*SPA(3,5)*SPB(4,3),2)*SPB(4,3))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SPA(4,5)*SPB(5,4))+
 (complex<T>(1,0)*S(2,4)*SPA(2,3)*SPB(4,3))/(complex<T>(9,0)*pow(SPA(3,4),2)*
 SPB(5,2)*SPB(5,4))+(complex<T>(1,0)*S(4,5)*SPA(2,3)*SPB(4,3))/
(complex<T>(9,0)*pow(SPA(3,4),2)*SPB(5,2)*SPB(5,4))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*S(2,4)*SPA(2,3)*SPA(2,5)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SS(2,3,4))-
 (pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPA(2,3))/
(complex<T>(9,0)*pow(SPA(3,4),2)*SPB(3,2)*SS(2,3,4))+
 (pow(SPB(4,3),2)*SPA(2,5)*SPA(4,5))/(complex<T>(18,0)*SPA(3,4)*SPB(3,2)*
 SS(2,3,4))-(pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(18,0)*SPA(3,4)*
 SS(2,3,4))-(pow(-(S(2,4)*SPA(2,5))-SPA(2,4)*SPA(3,5)*
  SPB(4,3),2)*SPB(4,3))/(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*
 SPA(3,4)*SS(2,3,4))-(complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,3)*
 SPA(2,4)*SPA(3,5)*SPA(4,5)*SPB(4,3))/
(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*S(3,5)*SPA(2,3)*SPA(2,5)*SPA(4,5))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SS(3,4,5))-
 (pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(18,0)*SPA(3,4)*SS(3,4,5))-
 (pow(-(S(3,5)*SPA(2,5))-SPA(2,4)*SPA(3,5)*SPB(4,3),2)*SPB(4,3))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,3)*SPA(2,4)*SPA(3,5)*SPA(4,5)*
 SPB(4,3))/(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4)*
 SS(3,4,5))+(pow(SPB(4,3),2)*SPA(2,3)*SPA(2,5))/
(complex<T>(18,0)*SPA(3,4)*SPB(5,4)*SS(3,4,5))-
 (pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPA(4,5))/
(complex<T>(9,0)*pow(SPA(3,4),2)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmpmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmpmp G");
#endif
 
 return(-((pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
 (complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPA(3,5)))-
 (pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),3))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(3,5))-
 (pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPA(3,5))+
 (complex<T>(0,-2)*pow(SPA(2,4),4))/(SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*S(2,5)*(S(2,3)+S(3,5)))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*S(2,3)*(S(2,5)+S(3,5)))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(3,4)+S(3,5))*S(4,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,3),3)*SPA(3,4))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*SPA(3,5)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(3,5)+S(4,5))*SPB(4,3))-
 (pow(SPB(5,3),2)*SPA(3,4)*SPA(4,5))/(complex<T>(2,0)*pow(SPA(3,5),2)*
 SPB(3,2)*SPB(5,2))-(pow(SPB(5,3),2)*S(2,5)*SPA(2,3))/
(complex<T>(2,0)*pow(SPA(3,5),2)*SPB(4,3)*SPB(5,2)*SPB(5,4))+
 (complex<T>(0,2)*pow(SPB(5,3),4))/(SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),3))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*SPA(3,5)*SS(2,3,5))-
 (pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPA(3,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(2,3)+S(3,5))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(2,5)+S(3,5))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,3),3)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*SPA(3,5)*SPB(3,2)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),3))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*SPA(3,5)*SS(3,4,5))-
 (pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*SPA(3,5)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(3,4)+S(3,5))*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(3,5)+S(4,5))*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),3)*S(2,5)*SPA(2,3))/(complex<T>(6,0)*SPA(3,5)*
 SPB(4,3)*SPB(5,2)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmpmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmpmp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPA(3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),3))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(3,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPA(3,5))-
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*S(2,5)*(S(2,3)+S(3,5)))-
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*S(2,3)*(S(2,5)+S(3,5)))-
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*(S(3,4)+S(3,5))*S(4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,3),3)*SPA(3,4))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*SPA(3,5)*SPB(4,3))-
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*(S(3,5)+S(4,5))*SPB(4,3))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*S(2,5)*SPA(2,3))/(complex<T>(6,0)*pow(SPA(3,5),2)*
 SPB(4,3)*SPB(5,2)*SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,5),2)*
 pow(SPA(3,4),2)*pow(SPB(5,3),3))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*SPA(3,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*SPA(3,5)*SS(2,3,5))-
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*(S(2,3)+S(3,5))*SS(2,3,5))-
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*(S(2,5)+S(3,5))*SS(2,3,5))-
 (pow(SPB(5,3),3)*SPA(3,4)*SPA(4,5))/(complex<T>(18,0)*SPA(3,5)*SPB(3,2)*
 SPB(5,2)*SS(2,3,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(3,4),2)*
 pow(SPB(5,3),3))/(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*SPA(3,5)*
 SS(3,4,5))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPA(4,5),2)*
 pow(SPB(5,3),3))/(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*SPA(3,5)*
 SS(3,4,5))-(pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))/(complex<T>(6,0)*pow(SPA(3,5),2)*(S(3,4)+S(3,5))*
 SS(3,4,5))-(pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))/(complex<T>(6,0)*pow(SPA(3,5),2)*(S(3,5)+S(4,5))*
 SS(3,4,5))-(pow(SPB(5,3),3)*S(2,5)*SPA(2,3))/
(complex<T>(18,0)*SPA(3,5)*SPB(4,3)*SPB(5,2)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpmpm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpmpm G");
#endif
 
 return(-((pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(4,2),3))/
 (complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPA(2,4)))+
 (complex<T>(0,-2)*pow(SPA(3,5),4))/(SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,3)+S(2,4))*S(3,4))-
 (pow(SPA(4,5),2)*pow(SPB(4,2),3)*SPA(2,3))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPA(2,4)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(3,4))*SPB(3,2))-
 (pow(SPA(3,4),2)*pow(SPB(4,2),3)*SPA(2,5))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*SPA(2,4)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(4,5))*SPB(5,2))-
 (pow(SPB(4,2),2)*S(2,5)*SPA(4,5))/(complex<T>(2,0)*pow(SPA(2,4),2)*
 SPB(3,2)*SPB(4,3)*SPB(5,2))+(complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*
 SPA(2,5)*SPA(3,4))/(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(2,5))*
 SPB(5,4))-(pow(SPA(2,3),2)*pow(SPB(4,2),3)*SPA(4,5))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*SPA(2,4)*SPB(5,4))-
 (pow(SPB(4,2),2)*SPA(2,3)*SPA(3,4))/(complex<T>(2,0)*pow(SPA(2,4),2)*
 SPB(5,2)*SPB(5,4))+(complex<T>(0,2)*pow(SPB(4,2),4))/
(SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(4,2),3))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*SPA(2,4)*SS(2,3,4))-
 (pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(4,2),3))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPA(2,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,3)+S(2,4))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(3,4))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(4,2),4)*S(2,5)*SPA(4,5))/(complex<T>(6,0)*S(2,4)*SPB(3,2)*
 SPB(4,3)*SPB(5,2)*SS(2,3,4))-
 (pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(4,2),3))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*SPA(2,4)*SS(2,4,5))-
 (pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(4,2),3))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*SPA(2,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(2,5))*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(4,5))*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPB(4,2),4)*SPA(2,3)*SPA(3,4))/
(complex<T>(6,0)*S(2,4)*SPB(5,2)*SPB(5,4)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpmpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpmpm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(4,2),3))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPA(2,4))-
 (pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,3)+S(2,4))*S(3,4))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,2),3)*SPA(2,3))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*SPA(2,4)*SPB(3,2))-
 (pow(SPB(4,2),2)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,4)+S(3,4))*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(4,2),3)*SPA(2,5))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*SPA(2,4)*SPB(5,2))-
 (pow(SPB(4,2),2)*SPA(2,3)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,4)+S(4,5))*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*S(2,5)*SPA(4,5))/(complex<T>(6,0)*pow(SPA(2,4),2)*
 SPB(3,2)*SPB(4,3)*SPB(5,2))-(pow(SPB(4,2),2)*SPA(2,3)*
 SPA(2,5)*SPA(3,4))/(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,4)+S(2,5))*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,2),3)*SPA(4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*SPA(2,4)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(3,4))/
(complex<T>(6,0)*pow(SPA(2,4),2)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(4,2),3))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*SPA(2,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(4,2),3))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*SPA(2,4)*SS(2,3,4))-
 (pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,3)+S(2,4))*SS(2,3,4))-
 (pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,4)+S(3,4))*SS(2,3,4))-
 (pow(SPB(4,2),4)*S(2,5)*SPA(4,5))/(complex<T>(18,0)*S(2,4)*SPB(3,2)*
 SPB(4,3)*SPB(5,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(4,2),3))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*SPA(2,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPA(4,5),2)*pow(SPB(4,2),3))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*SPA(2,4)*SS(2,4,5))-
 (pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,4)+S(2,5))*SS(2,4,5))-
 (pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,4)+S(4,5))*SS(2,4,5))-
 (pow(SPB(4,2),4)*SPA(2,3)*SPA(3,4))/(complex<T>(18,0)*S(2,4)*SPB(5,2)*
 SPB(5,4)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpmmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpmmm G");
#endif
 
 return(-((pow(SPA(4,5),2)*pow(SPB(5,2),2))/(complex<T>(2,0)*pow(SPB(5,3),2)*
S(2,3)))-(pow(SPA(3,4),2)*pow(SPB(3,2),2))/
(complex<T>(2,0)*pow(SPB(5,3),2)*S(2,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2))/(complex<T>(2,0)*pow(SPB(5,3),2)*
 (S(2,3)+S(3,5)))+(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(2,5)+S(3,5)))-
 (pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPA(3,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPB(4,3))+
 (pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPA(3,4))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,3))-
 (pow(SPA(3,5),2)*SPA(4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPB(4,3))+
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPB(5,3))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(6,0)*S(2,5)*(S(2,3)+S(3,5))*SPB(5,3))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5))*SPB(5,3))-
 (SPA(3,4)*SPA(3,5)*SPA(4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*
 SPB(5,3))+(complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*
 SPB(5,2))/(complex<T>(6,0)*S(2,5)*(S(2,3)+S(3,5))*SPB(5,3))+
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5))*SPB(5,3))-
 (pow(SPA(3,5),2)*SPA(3,4))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPB(5,4))-
 (pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5))/
(complex<T>(6,0)*S(2,5)*(S(2,4)+S(4,5))*SPB(5,4))-
 (pow(SPA(3,5),2)*S(3,4))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*
 SPB(5,4))-(pow(SPA(3,5),2)*S(4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*
 SPB(4,3)*SPB(5,4))+(complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*
 SPB(4,2)*SPB(5,2))/(complex<T>(6,0)*S(2,5)*(S(2,4)+S(4,5))*SPB(5,4))-
 (S(3,4)*S(4,5)*SPA(3,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*
 SPB(5,3)*SPB(5,4))+(pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPA(3,4))/
(complex<T>(6,0)*S(2,3)*SPB(4,3)*SS(2,3,4))+
 (pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPA(3,4))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPB(4,3)*SS(2,3,4))+
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(6,0)*S(2,3)*SPB(4,3)*SS(2,3,4))+
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPB(4,3)*SS(2,3,4))+
 (complex<T>(2,0)*pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),3))/
(SPB(3,2)*SPB(4,3)*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*
 SS(2,3,4))-(pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPA(3,4)*
 SS(2,3,4))/(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,3))-
 pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)/
(pow(SPB(5,3),2)*SS(2,3,5))-(pow(SPA(4,5),2)*pow(SPB(5,2),2)*
 SPA(3,5))/(complex<T>(6,0)*S(2,3)*SPB(5,3)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(6,0)*S(2,5)*SPB(5,3)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(6,0)*(S(2,3)+S(3,5))*SPB(5,3)*SS(2,3,5))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SPB(5,3)*SS(2,3,5))+
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(6,0)*S(2,3)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(6,0)*S(2,5)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(6,0)*(S(2,3)+S(3,5))*SPB(5,3)*SS(2,3,5))+
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SPB(5,3)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(mH2,2)*pow(SPA(3,5),4))/(SPA(2,3)*SPA(2,5)*
 (SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*(SPA(2,3)*SPB(4,2)+
SPA(3,5)*SPB(5,4))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SS(2,3,5))/
(complex<T>(2,0)*pow(SPB(5,3),2)*S(2,5)*(S(2,3)+S(3,5)))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SS(2,3,5))/
(complex<T>(2,0)*pow(SPB(5,3),2)*S(2,3)*(S(2,5)+S(3,5)))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5)*SS(2,3,5))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5)*SS(2,3,5))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5))/
(complex<T>(6,0)*S(2,5)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5))/
(complex<T>(6,0)*(S(2,4)+S(4,5))*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(6,0)*S(2,5)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(6,0)*(S(2,4)+S(4,5))*SPB(5,4)*SS(2,4,5))+
 (complex<T>(-2,0)*pow(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2),3))/
(SPB(5,2)*SPB(5,4)*(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4))*
 SS(2,4,5))-(pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5)*
 SS(2,4,5))/(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phpmmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phpmmm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2))/
(complex<T>(6,0)*pow(SPB(5,3),2)*S(2,3))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2))/(complex<T>(6,0)*pow(SPB(5,3),2)*
 S(2,5))-(pow(SPA(3,4),2)*pow(SPB(3,2),2))/
(complex<T>(6,0)*pow(SPB(5,3),2)*(S(2,3)+S(3,5)))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2))/(complex<T>(6,0)*pow(SPB(5,3),2)*
 (S(2,5)+S(3,5)))+(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,2),2)*
 SPA(3,4))/(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*SPB(4,3))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPA(3,4))/
(complex<T>(18,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(4,5))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPB(4,3))-(complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*
 SPB(4,2))/(complex<T>(18,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,3))-
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*SPB(5,3))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*SPB(5,3))+
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(18,0)*S(2,5)*(S(2,3)+S(3,5))*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(18,0)*S(2,3)*(S(2,5)+S(3,5))*SPB(5,3))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPB(5,3))-(SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(18,0)*S(2,5)*(S(2,3)+S(3,5))*SPB(5,3))-
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(18,0)*S(2,3)*(S(2,5)+S(3,5))*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(3,4))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*SPB(5,4))-
 (pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5))/
(complex<T>(18,0)*S(2,5)*(S(2,4)+S(4,5))*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*S(3,4))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPB(4,3)*SPB(5,4))+(complex<T>(1,0)*pow(SPA(3,5),2)*S(4,5))/
(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,4))-
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(18,0)*S(2,5)*(S(2,4)+S(4,5))*SPB(5,4))+
 (complex<T>(1,0)*S(3,4)*S(4,5)*SPA(3,5))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPB(4,3)*SPB(5,3)*SPB(5,4))-(complex<T>(1,0)*pow(SPA(4,5),2)*
 pow(SPB(4,2),2)*SPA(3,4))/(complex<T>(18,0)*S(2,3)*SPB(4,3)*
 SS(2,3,4))-(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPA(3,4))/
(complex<T>(18,0)*(S(2,4)+S(3,4))*SPB(4,3)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(18,0)*S(2,3)*SPB(4,3)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(18,0)*(S(2,4)+S(3,4))*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPA(3,4)*SS(2,3,4))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,3))+
 (complex<T>(1,0)*pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2))/
(complex<T>(3,0)*pow(SPB(5,3),2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(18,0)*S(2,3)*SPB(5,3)*SS(2,3,5))+
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(18,0)*S(2,5)*SPB(5,3)*SS(2,3,5))+
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(18,0)*(S(2,3)+S(3,5))*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(18,0)*(S(2,5)+S(3,5))*SPB(5,3)*SS(2,3,5))-
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(18,0)*S(2,3)*SPB(5,3)*SS(2,3,5))-
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(18,0)*S(2,5)*SPB(5,3)*SS(2,3,5))-
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(18,0)*(S(2,3)+S(3,5))*SPB(5,3)*SS(2,3,5))-
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(18,0)*(S(2,5)+S(3,5))*SPB(5,3)*SS(2,3,5))-
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SS(2,3,5))/
(complex<T>(6,0)*pow(SPB(5,3),2)*S(2,5)*(S(2,3)+S(3,5)))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SS(2,3,5))/
(complex<T>(6,0)*pow(SPB(5,3),2)*S(2,3)*(S(2,5)+S(3,5)))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5)*SS(2,3,5))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,3))-
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPA(3,5)*SS(2,3,5))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPB(5,3))-
 (pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5))/
(complex<T>(18,0)*S(2,5)*SPB(5,4)*SS(2,4,5))-
 (pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5))/
(complex<T>(18,0)*(S(2,4)+S(4,5))*SPB(5,4)*SS(2,4,5))-
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(18,0)*S(2,5)*SPB(5,4)*SS(2,4,5))-
 (SPA(3,4)*SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(18,0)*(S(2,4)+S(4,5))*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(4,2),2)*SPA(4,5)*SS(2,4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmpmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmpmm G");
#endif
 
 return(-((pow(SPA(4,5),2)*pow(SPB(4,3),2))/(complex<T>(2,0)*pow(SPB(4,2),2)*
S(2,3)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,3)+S(2,4)))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2))/(complex<T>(2,0)*pow(SPB(4,2),2)*
 S(3,4))+(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,4)+S(3,4)))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPB(4,2))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(4,2))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,2))+
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(4,2))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,2))-
 (pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5))*SPB(5,2))-
 (pow(SPA(2,4),2)*SPA(4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*SPB(5,2))-
 (S(2,5)*SPA(2,4)*SPA(4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*
 SPB(5,2))+(complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*
 SPB(5,3))/(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5))*SPB(5,2))-
 (pow(SPA(2,4),2)*SPA(2,5))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*SPB(5,4))-
 (pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5))/
(complex<T>(6,0)*S(3,4)*(S(3,5)+S(4,5))*SPB(5,4))-
 (S(4,5)*SPA(2,4)*SPA(2,5))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*
 SPB(5,4))-(pow(SPA(2,4),2)*S(2,5))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*
 SPB(5,2)*SPB(5,4))-(pow(SPA(2,4),2)*S(4,5))/
(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(4,3)*SPB(5,3))/
(complex<T>(6,0)*S(3,4)*(S(3,5)+S(4,5))*SPB(5,4))-
 pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)/
(pow(SPB(4,2),2)*SS(2,3,4))-(complex<T>(1,0)*pow(SPA(4,5),2)*
 pow(SPB(4,3),2)*SPA(2,4))/(complex<T>(6,0)*S(2,3)*SPB(4,2)*SS(2,3,4))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*SPB(4,2)*SS(2,3,4))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4))/
(complex<T>(6,0)*S(3,4)*SPB(4,2)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPB(4,2)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*S(2,3)*SPB(4,2)*SS(2,3,4))+
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*SPB(4,2)*SS(2,3,4))+
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*S(3,4)*SPB(4,2)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPB(4,2)*SS(2,3,4))+
 (complex<T>(2,0)*pow(mH2,2)*pow(SPA(2,4),4))/(SPA(2,3)*SPA(3,4)*
 (SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*(-(SPA(2,3)*SPB(5,3))-
SPA(2,4)*SPB(5,4))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SS(2,3,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,3)+S(2,4))*S(3,4))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SS(2,3,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*S(2,3)*(S(2,4)+S(3,4)))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4)*SS(2,3,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4)*SS(2,3,4))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(6,0)*S(2,3)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(6,0)*S(2,3)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SPB(5,2)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),3))/
(SPB(3,2)*SPB(5,2)*(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*
 SS(2,3,5))-(pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5)*
 SS(2,3,5))/(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5))/
(complex<T>(6,0)*S(3,4)*SPB(5,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5))/
(complex<T>(6,0)*(S(3,5)+S(4,5))*SPB(5,4)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(4,3)*SPB(5,3))/
(complex<T>(6,0)*S(3,4)*SPB(5,4)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(4,3)*SPB(5,3))/
(complex<T>(6,0)*(S(3,5)+S(4,5))*SPB(5,4)*SS(3,4,5))+
 (complex<T>(2,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),3))/
(SPB(4,3)*SPB(5,4)*(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4))*
 SS(3,4,5))-(pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5)*
 SS(3,4,5))/(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmpmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmpmm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2))/
(complex<T>(6,0)*pow(SPB(4,2),2)*S(2,3))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2))/(complex<T>(6,0)*pow(SPB(4,2),2)*
 (S(2,3)+S(2,4)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2))/
(complex<T>(6,0)*pow(SPB(4,2),2)*S(3,4))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2))/(complex<T>(6,0)*pow(SPB(4,2),2)*
 (S(2,4)+S(3,4)))-(pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*SPB(4,2))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(4,2))+
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4))/
(complex<T>(18,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,2))-
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(4,2))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(18,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*SPB(5,2))-
 (pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(18,0)*S(2,3)*(S(2,5)+S(3,5))*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(4,5))/(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*
 SPB(5,2))+(complex<T>(1,0)*S(2,5)*SPA(2,4)*SPA(4,5))/
(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(18,0)*S(2,3)*(S(2,5)+S(3,5))*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(2,5))/(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*SPB(5,4))-
 (pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5))/
(complex<T>(18,0)*S(3,4)*(S(3,5)+S(4,5))*SPB(5,4))+
 (complex<T>(1,0)*S(4,5)*SPA(2,4)*SPA(2,5))/(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*
 SPB(4,2)*SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,4),2)*S(2,5))/
(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*S(4,5))/(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*
 SPB(5,2)*SPB(5,4))-(SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(4,3)*
 SPB(5,3))/(complex<T>(18,0)*S(3,4)*(S(3,5)+S(4,5))*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2))/
(complex<T>(3,0)*pow(SPB(4,2),2)*SS(2,3,4))+
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4))/
(complex<T>(18,0)*S(2,3)*SPB(4,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*SPB(4,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4))/
(complex<T>(18,0)*S(3,4)*SPB(4,2)*SS(2,3,4))+
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4))/
(complex<T>(18,0)*(S(2,4)+S(3,4))*SPB(4,2)*SS(2,3,4))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(18,0)*S(2,3)*SPB(4,2)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*SPB(4,2)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(18,0)*S(3,4)*SPB(4,2)*SS(2,3,4))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,4)+S(3,4))*SPB(4,2)*SS(2,3,4))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SS(2,3,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*(S(2,3)+S(2,4))*S(3,4))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SS(2,3,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*S(2,3)*(S(2,4)+S(3,4)))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,4)*SS(2,3,4))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,2))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPA(2,4)*SS(2,3,4))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(4,2))-
 (pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(18,0)*S(2,3)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(18,0)*(S(2,5)+S(3,5))*SPB(5,2)*SS(2,3,5))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(18,0)*S(2,3)*SPB(5,2)*SS(2,3,5))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(18,0)*(S(2,5)+S(3,5))*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(2,5)*SS(2,3,5))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,2))-
 (pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5))/
(complex<T>(18,0)*S(3,4)*SPB(5,4)*SS(3,4,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5))/
(complex<T>(18,0)*(S(3,5)+S(4,5))*SPB(5,4)*SS(3,4,5))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(4,3)*SPB(5,3))/
(complex<T>(18,0)*S(3,4)*SPB(5,4)*SS(3,4,5))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(4,3)*SPB(5,3))/
(complex<T>(18,0)*(S(3,5)+S(4,5))*SPB(5,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,3),2)*SPA(4,5)*SS(3,4,5))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmmpm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmmpm G");
#endif
 
 return(-((pow(SPA(2,5),2)*pow(SPB(5,4),2))/(complex<T>(2,0)*pow(SPB(5,3),2)*
S(3,4)))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(3,4)+S(3,5)))-
 (pow(SPA(2,3),2)*pow(SPB(4,3),2))/(complex<T>(2,0)*pow(SPB(5,3),2)*
 S(4,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(3,5)+S(4,5)))-
 (pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(3,2))-
 (pow(SPA(3,5),2)*SPA(2,5))/(complex<T>(6,0)*SPA(3,4)*SPA(4,5)*SPB(3,2))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(3,2))-
 (pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*S(4,5)*SPB(5,2))-
 (pow(SPA(3,5),2)*SPA(2,3))/(complex<T>(6,0)*SPA(3,4)*SPA(4,5)*SPB(5,2))-
 (pow(SPA(3,5),2)*S(2,3))/(complex<T>(6,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*
 SPB(5,2))-(pow(SPA(3,5),2)*S(2,5))/(complex<T>(6,0)*SPA(3,4)*SPA(4,5)*
 SPB(3,2)*SPB(5,2))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*
 SPA(3,5))/(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*SPB(5,3))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*S(4,5)*SPB(5,3))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5))/
(complex<T>(6,0)*S(3,4)*(S(3,5)+S(4,5))*SPB(5,3))-
 (S(2,3)*SPA(2,5)*SPA(3,5))/(complex<T>(6,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*
 SPB(5,3))-(S(2,5)*SPA(2,3)*SPA(3,5))/(complex<T>(6,0)*SPA(3,4)*
 SPA(4,5)*SPB(5,2)*SPB(5,3))+(complex<T>(1,0)*SPA(2,3)*SPA(2,5)*
 SPA(3,5)*SPB(4,2)*SPB(5,4))/(complex<T>(6,0)*(S(2,4)+S(2,5))*S(4,5)*
 SPB(5,2))+(complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*
 SPB(5,4))/(complex<T>(6,0)*(S(3,4)+S(3,5))*S(4,5)*SPB(5,3))+
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*S(3,4)*(S(3,5)+S(4,5))*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3))/
(complex<T>(6,0)*S(3,4)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3))/
(SPB(3,2)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(4,3)*
 SS(2,3,4))+(complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*
 SPB(4,3))/(complex<T>(6,0)*(S(2,3)+S(2,4))*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*S(3,4)*SPB(3,2)*SS(2,3,4))-
 (pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3)*SS(2,3,4))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*SPB(5,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5))/
(complex<T>(6,0)*S(4,5)*SPB(5,2)*SS(2,4,5))+
 (complex<T>(-2,0)*pow(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4),3))/
(SPB(5,2)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*SPB(5,4)*
 SS(2,4,5))+(complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*
 SPB(5,4))/(complex<T>(6,0)*(S(2,4)+S(2,5))*SPB(5,2)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(6,0)*S(4,5)*SPB(5,2)*SS(2,4,5))-
 (pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5)*SS(2,4,5))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(5,2))-
 pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)/
(pow(SPB(5,3),2)*SS(3,4,5))+
 (complex<T>(2,0)*pow(mH2,2)*pow(SPA(3,5),4))/(SPA(3,4)*SPA(4,5)*
 (-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*(SPA(3,4)*SPB(4,2)+
SPA(3,5)*SPB(5,2))*SS(3,4,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5))/
(complex<T>(6,0)*S(3,4)*SPB(5,3)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*SPB(5,3)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5))/
(complex<T>(6,0)*S(4,5)*SPB(5,3)*SS(3,4,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5))/
(complex<T>(6,0)*(S(3,5)+S(4,5))*SPB(5,3)*SS(3,4,5))+
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*S(3,4)*SPB(5,3)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*SPB(5,3)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*S(4,5)*SPB(5,3)*SS(3,4,5))+
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*(S(3,5)+S(4,5))*SPB(5,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SS(3,4,5))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(3,4)+S(3,5))*S(4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SS(3,4,5))/
(complex<T>(2,0)*pow(SPB(5,3),2)*S(3,4)*(S(3,5)+S(4,5)))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5)*SS(3,4,5))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5)*SS(3,4,5))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(5,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmmpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmmpm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2))/
(complex<T>(6,0)*pow(SPB(5,3),2)*S(3,4))-
 (pow(SPA(2,3),2)*pow(SPB(4,3),2))/(complex<T>(6,0)*pow(SPB(5,3),2)*
 (S(3,4)+S(3,5)))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2))/
(complex<T>(6,0)*pow(SPB(5,3),2)*S(4,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2))/(complex<T>(6,0)*pow(SPB(5,3),2)*
 (S(3,5)+S(4,5)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(4,2),2)*
 SPA(2,3))/(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*SPB(3,2))-
 (pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(2,5))/(complex<T>(18,0)*SPA(3,4)*SPA(4,5)*
 SPB(3,2))-(SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*SPB(5,2))-
 (pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*S(4,5)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(2,3))/(complex<T>(18,0)*SPA(3,4)*SPA(4,5)*
 SPB(5,2))+(complex<T>(1,0)*pow(SPA(3,5),2)*S(2,3))/
(complex<T>(18,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*S(2,5))/(complex<T>(18,0)*SPA(3,4)*SPA(4,5)*
 SPB(3,2)*SPB(5,2))-(pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*SPB(5,3))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*SPB(5,3))+
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*S(4,5)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5))/
(complex<T>(18,0)*S(3,4)*(S(3,5)+S(4,5))*SPB(5,3))+
 (complex<T>(1,0)*S(2,3)*SPA(2,5)*SPA(3,5))/(complex<T>(18,0)*SPA(3,4)*SPA(4,5)*
 SPB(3,2)*SPB(5,3))+(complex<T>(1,0)*S(2,5)*SPA(2,3)*SPA(3,5))/
(complex<T>(18,0)*SPA(3,4)*SPA(4,5)*SPB(5,2)*SPB(5,3))-
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*S(4,5)*SPB(5,2))-
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*S(4,5)*SPB(5,3))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*S(3,4)*(S(3,5)+S(4,5))*SPB(5,3))-
 (pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*SPB(3,2)*SS(2,3,4))-
 (pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3))/
(complex<T>(18,0)*S(3,4)*SPB(3,2)*SS(2,3,4))-
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*SPB(3,2)*SS(2,3,4))-
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*S(3,4)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(2,3)*SS(2,3,4))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(3,2))-
 (pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*SPB(5,2)*SS(2,4,5))-
 (pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5))/
(complex<T>(18,0)*S(4,5)*SPB(5,2)*SS(2,4,5))-
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*SPB(5,2)*SS(2,4,5))-
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(18,0)*S(4,5)*SPB(5,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,2),2)*SPA(2,5)*SS(2,4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(5,2))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2))/
(complex<T>(3,0)*pow(SPB(5,3),2)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5))/
(complex<T>(18,0)*S(3,4)*SPB(5,3)*SS(3,4,5))+
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*SPB(5,3)*SS(3,4,5))+
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5))/
(complex<T>(18,0)*S(4,5)*SPB(5,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5))/
(complex<T>(18,0)*(S(3,5)+S(4,5))*SPB(5,3)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*S(3,4)*SPB(5,3)*SS(3,4,5))-
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*SPB(5,3)*SS(3,4,5))-
 (SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*S(4,5)*SPB(5,3)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*(S(3,5)+S(4,5))*SPB(5,3)*SS(3,4,5))-
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SS(3,4,5))/
(complex<T>(6,0)*pow(SPB(5,3),2)*(S(3,4)+S(3,5))*S(4,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SS(3,4,5))/
(complex<T>(6,0)*pow(SPB(5,3),2)*S(3,4)*(S(3,5)+S(4,5)))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5)*SS(3,4,5))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,3))-
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPA(3,5)*SS(3,4,5))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(5,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmmmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmmmp G");
#endif
 
 return(-((pow(SPA(3,4),2)*pow(SPB(5,4),2))/(complex<T>(2,0)*pow(SPB(4,2),2)*
S(2,5)))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,4)+S(2,5)))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2))/(complex<T>(2,0)*pow(SPB(4,2),2)*
 S(4,5))+(complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,4)+S(4,5)))-
 (pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(6,0)*S(2,5)*(S(2,3)+S(3,5))*SPB(3,2))-
 (pow(SPA(2,4),2)*SPA(3,4))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*SPB(4,2))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*S(4,5)*SPB(4,2))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4))/
(complex<T>(6,0)*S(2,5)*(S(2,4)+S(4,5))*SPB(4,2))-
 (SPA(2,3)*SPA(2,4)*SPA(3,4))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*
 SPB(4,2))-(pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*SPB(4,3))+
 (pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*S(4,5)*SPB(4,3))-
 (pow(SPA(2,4),2)*SPA(2,3))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*SPB(4,3))-
 (pow(SPA(2,4),2)*S(2,3))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*
 SPB(4,3))-(pow(SPA(2,4),2)*S(3,4))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*
 SPB(3,2)*SPB(4,3))-(S(2,3)*S(3,4)*SPA(2,4))/
(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*S(2,5)*(S(2,3)+S(3,5))*SPB(3,2))+
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*S(4,5)*SPB(4,2))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(2,5)*(S(2,4)+S(4,5))*SPB(4,2))+
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*S(4,5)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(6,0)*S(2,5)*SPB(3,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(6,0)*(S(2,3)+S(3,5))*SPB(3,2)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*S(2,5)*SPB(3,2)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*(S(2,3)+S(3,5))*SPB(3,2)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),3))/
(SPB(3,2)*SPB(5,2)*(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3))*
 SS(2,3,5))-(pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3)*
 SS(2,3,5))/(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPB(3,2))-
 pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)/
(pow(SPB(4,2),2)*SS(2,4,5))-(complex<T>(1,0)*pow(SPA(3,4),2)*
 pow(SPB(5,4),2)*SPA(2,4))/(complex<T>(6,0)*S(2,5)*SPB(4,2)*SS(2,4,5))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*SPB(4,2)*SS(2,4,5))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(6,0)*S(4,5)*SPB(4,2)*SS(2,4,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4))/
(complex<T>(6,0)*(S(2,4)+S(4,5))*SPB(4,2)*SS(2,4,5))+
 (complex<T>(-2,0)*pow(mH2,2)*pow(SPA(2,4),4))/(SPA(2,5)*SPA(4,5)*
 (SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3))*(SPA(2,4)*SPB(3,2)+
SPA(4,5)*SPB(5,3))*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(2,5)*SPB(4,2)*SS(2,4,5))+
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*SPB(4,2)*SS(2,4,5))+
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(4,5)*SPB(4,2)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*(S(2,4)+S(4,5))*SPB(4,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SS(2,4,5))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,4)+S(2,5))*S(4,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SS(2,4,5))/
(complex<T>(2,0)*pow(SPB(4,2),2)*S(2,5)*(S(2,4)+S(4,5)))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4)*SS(2,4,5))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4)*SS(2,4,5))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(4,2))+
 (pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*SPB(4,3)*SS(3,4,5))+
 (pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4))/
(complex<T>(6,0)*S(4,5)*SPB(4,3)*SS(3,4,5))+
 (complex<T>(2,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),3))/
(SPB(4,3)*(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3))*SPB(5,4)*
 SS(3,4,5))+(SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*SPB(4,3)*SS(3,4,5))+
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(4,5)*SPB(4,3)*SS(3,4,5))-
 (pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4)*SS(3,4,5))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmmmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmmmp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2))/
(complex<T>(6,0)*pow(SPB(4,2),2)*S(2,5))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2))/(complex<T>(6,0)*pow(SPB(4,2),2)*
 (S(2,4)+S(2,5)))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2))/
(complex<T>(6,0)*pow(SPB(4,2),2)*S(4,5))-
 (pow(SPA(3,4),2)*pow(SPB(5,4),2))/(complex<T>(6,0)*pow(SPB(4,2),2)*
 (S(2,4)+S(4,5)))+(complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*
 SPA(2,3))/(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*SPB(3,2))-
 (pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(18,0)*S(2,5)*(S(2,3)+S(3,5))*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(3,4))/(complex<T>(18,0)*SPA(2,5)*SPA(4,5)*
 SPB(3,2))-(pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*SPB(4,2))-
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*S(4,5)*SPB(4,2))+
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4))/
(complex<T>(18,0)*S(2,5)*(S(2,4)+S(4,5))*SPB(4,2))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4))/(complex<T>(18,0)*SPA(2,5)*SPA(4,5)*
 SPB(4,2))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*SPB(4,3))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*S(4,5)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(2,3))/(complex<T>(18,0)*SPA(2,5)*SPA(4,5)*
 SPB(4,3))+(complex<T>(1,0)*pow(SPA(2,4),2)*S(2,3))/
(complex<T>(18,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*S(3,4))/(complex<T>(18,0)*SPA(2,5)*SPA(4,5)*
 SPB(3,2)*SPB(4,3))+(complex<T>(1,0)*S(2,3)*S(3,4)*SPA(2,4))/
(complex<T>(18,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))-
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*S(2,5)*(S(2,3)+S(3,5))*SPB(3,2))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*S(4,5)*SPB(4,2))-
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(2,5)*(S(2,4)+S(4,5))*SPB(4,2))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*S(4,5)*SPB(4,3))-
 (pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(18,0)*S(2,5)*SPB(3,2)*SS(2,3,5))-
 (pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(18,0)*(S(2,3)+S(3,5))*SPB(3,2)*SS(2,3,5))-
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*S(2,5)*SPB(3,2)*SS(2,3,5))-
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*(S(2,3)+S(3,5))*SPB(3,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,3)*SS(2,3,5))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2))/
(complex<T>(3,0)*pow(SPB(4,2),2)*SS(2,4,5))+
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4))/
(complex<T>(18,0)*S(2,5)*SPB(4,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*SPB(4,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(18,0)*S(4,5)*SPB(4,2)*SS(2,4,5))+
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4))/
(complex<T>(18,0)*(S(2,4)+S(4,5))*SPB(4,2)*SS(2,4,5))-
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(2,5)*SPB(4,2)*SS(2,4,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*SPB(4,2)*SS(2,4,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(4,5)*SPB(4,2)*SS(2,4,5))-
 (SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(4,5))*SPB(4,2)*SS(2,4,5))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SS(2,4,5))/
(complex<T>(6,0)*pow(SPB(4,2),2)*(S(2,4)+S(2,5))*S(4,5))-
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SS(2,4,5))/
(complex<T>(6,0)*pow(SPB(4,2),2)*S(2,5)*(S(2,4)+S(4,5)))-
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPA(2,4)*SS(2,4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPB(4,2))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4)*SS(2,4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(4,2))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*SPB(4,3)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4))/
(complex<T>(18,0)*S(4,5)*SPB(4,3)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*SPB(4,3)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*S(4,5)*SPB(4,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(3,4)*SS(3,4,5))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmmmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmmmm G");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(4,5),2))/(complex<T>(3,0)*pow(SPB(3,2),2))+
 (complex<T>(1,0)*pow(SPA(2,5),2))/(complex<T>(3,0)*pow(SPB(4,3),2))+
 (complex<T>(1,0)*pow(SPA(3,4),2))/(complex<T>(3,0)*pow(SPB(5,2),2))+
 (complex<T>(1,0)*pow(SPA(2,3),2))/(complex<T>(3,0)*pow(SPB(5,4),2))+
 (complex<T>(-2,0)*pow(mH2,2))/(SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(-2,0)*SPA(2,5)*SPA(4,5))/(complex<T>(3,0)*SPB(3,2)*SPB(4,3))+
 (complex<T>(-2,0)*SPA(3,4)*SPA(4,5))/(complex<T>(3,0)*SPB(3,2)*SPB(5,2))-
 (SPA(2,3)*SPA(4,5))/(complex<T>(3,0)*SPB(4,3)*SPB(5,2))-
 (SPA(2,5)*SPA(3,4))/(complex<T>(3,0)*SPB(3,2)*SPB(5,4))+
 (complex<T>(-2,0)*SPA(2,3)*SPA(2,5))/(complex<T>(3,0)*SPB(4,3)*SPB(5,4))+
 (complex<T>(-2,0)*SPA(2,3)*SPA(3,4))/(complex<T>(3,0)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(3,2))/(complex<T>(6,0)*SPB(4,3)*SPB(5,2)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(3,4),2)*SPB(4,3))/
(complex<T>(6,0)*SPB(3,2)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(5,2))/(complex<T>(6,0)*SPB(3,2)*SPB(4,3)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(4,5),2)*SPB(5,4))/
(complex<T>(6,0)*SPB(3,2)*SPB(4,3)*SPB(5,2))-
 (pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*S(2,4))/
(complex<T>(3,0)*pow(SPB(3,2),2)*pow(SPB(4,3),2)*SS(2,3,4))-
 (pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*S(3,5))/
(complex<T>(3,0)*pow(SPB(3,2),2)*pow(SPB(5,2),2)*SS(2,3,5))-
 (SS(2,3,4)*SS(2,3,5))/(complex<T>(6,0)*SPB(3,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,4))-(pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*
 S(2,4))/(complex<T>(3,0)*pow(SPB(5,2),2)*pow(SPB(5,4),2)*SS(2,4,5))-
 (SS(2,3,5)*SS(2,4,5))/(complex<T>(6,0)*SPB(3,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,4))-(pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*
 S(3,5))/(complex<T>(3,0)*pow(SPB(4,3),2)*pow(SPB(5,4),2)*SS(3,4,5))-
 (SS(2,3,4)*SS(3,4,5))/(complex<T>(6,0)*SPB(3,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,4))-(SS(2,4,5)*SS(3,4,5))/(complex<T>(6,0)*SPB(3,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phmmmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phmmmm nf");
#endif
 
 return(-(pow(SPA(4,5),2)/(complex<T>(9,0)*pow(SPB(3,2),2)))-
 pow(SPA(2,5),2)/(complex<T>(9,0)*pow(SPB(4,3),2))-
 pow(SPA(3,4),2)/(complex<T>(9,0)*pow(SPB(5,2),2))-
 pow(SPA(2,3),2)/(complex<T>(9,0)*pow(SPB(5,4),2))+
 (complex<T>(2,0)*SPA(2,5)*SPA(4,5))/(complex<T>(9,0)*SPB(3,2)*SPB(4,3))+
 (complex<T>(2,0)*SPA(3,4)*SPA(4,5))/(complex<T>(9,0)*SPB(3,2)*SPB(5,2))+
 (complex<T>(1,0)*SPA(2,3)*SPA(4,5))/(complex<T>(9,0)*SPB(4,3)*SPB(5,2))+
 (complex<T>(1,0)*SPA(2,5)*SPA(3,4))/(complex<T>(9,0)*SPB(3,2)*SPB(5,4))+
 (complex<T>(2,0)*SPA(2,3)*SPA(2,5))/(complex<T>(9,0)*SPB(4,3)*SPB(5,4))+
 (complex<T>(2,0)*SPA(2,3)*SPA(3,4))/(complex<T>(9,0)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(2,3),2)*SPB(3,2))/(complex<T>(18,0)*SPB(4,3)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(3,4),2)*SPB(4,3))/(complex<T>(18,0)*SPB(3,2)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(2,5),2)*SPB(5,2))/(complex<T>(18,0)*SPB(3,2)*SPB(4,3)*SPB(5,4))-
 (pow(SPA(4,5),2)*SPB(5,4))/(complex<T>(18,0)*SPB(3,2)*SPB(4,3)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*S(2,4))/
(complex<T>(9,0)*pow(SPB(3,2),2)*pow(SPB(4,3),2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*S(3,5))/
(complex<T>(9,0)*pow(SPB(3,2),2)*pow(SPB(5,2),2)*SS(2,3,5))+
 (complex<T>(1,0)*SS(2,3,4)*SS(2,3,5))/(complex<T>(18,0)*SPB(3,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*S(2,4))/
(complex<T>(9,0)*pow(SPB(5,2),2)*pow(SPB(5,4),2)*SS(2,4,5))+
 (complex<T>(1,0)*SS(2,3,5)*SS(2,4,5))/(complex<T>(18,0)*SPB(3,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*S(3,5))/
(complex<T>(9,0)*pow(SPB(4,3),2)*pow(SPB(5,4),2)*SS(3,4,5))+
 (complex<T>(1,0)*SS(2,3,4)*SS(3,4,5))/(complex<T>(18,0)*SPB(3,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,4))+(complex<T>(1,0)*SS(2,4,5)*SS(3,4,5))/
(complex<T>(18,0)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpppp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpppp G");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(3,2),2))/(complex<T>(3,0)*pow(SPA(4,5),2))+
 (complex<T>(1,0)*pow(SPB(4,3),2))/(complex<T>(3,0)*pow(SPA(2,5),2))+
 (complex<T>(1,0)*pow(SPB(5,2),2))/(complex<T>(3,0)*pow(SPA(3,4),2))+
 (complex<T>(1,0)*pow(SPB(5,4),2))/(complex<T>(3,0)*pow(SPA(2,3),2))-
 (S(2,5)*S(3,4))/(complex<T>(3,0)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))-
 (S(2,3)*S(4,5))/(complex<T>(3,0)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(2,3))/(complex<T>(6,0)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))+(complex<T>(1,0)*pow(SPB(5,2),2)*SPA(2,5))/
(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(3,4))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*
 SPA(4,5))+(complex<T>(1,0)*pow(SPB(5,4),2)*SPA(4,5))/
(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPA(3,4))+
 (complex<T>(-2,0)*SPB(3,2)*SPB(4,3))/(complex<T>(3,0)*SPA(2,5)*SPA(4,5))+
 (complex<T>(-2,0)*SPB(3,2)*SPB(5,2))/(complex<T>(3,0)*SPA(3,4)*SPA(4,5))+
 (complex<T>(-2,0)*pow(mH2,2))/(SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))+
 (complex<T>(-2,0)*SPB(4,3)*SPB(5,4))/(complex<T>(3,0)*SPA(2,3)*SPA(2,5))+
 (complex<T>(-2,0)*SPB(5,2)*SPB(5,4))/(complex<T>(3,0)*SPA(2,3)*SPA(3,4))-
 (pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*S(2,4))/
(complex<T>(3,0)*pow(SPA(2,3),2)*pow(SPA(3,4),2)*SS(2,3,4))-
 (pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*S(3,5))/
(complex<T>(3,0)*pow(SPA(2,3),2)*pow(SPA(2,5),2)*SS(2,3,5))-
 (SS(2,3,4)*SS(2,3,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))-(pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*
 S(2,4))/(complex<T>(3,0)*pow(SPA(2,5),2)*pow(SPA(4,5),2)*SS(2,4,5))-
 (SS(2,3,5)*SS(2,4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))-(pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*
 S(3,5))/(complex<T>(3,0)*pow(SPA(3,4),2)*pow(SPA(4,5),2)*SS(3,4,5))-
 (SS(2,3,4)*SS(3,4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))-(SS(2,4,5)*SS(3,4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpppp nf");
#endif
 
 return(-(pow(SPB(3,2),2)/(complex<T>(9,0)*pow(SPA(4,5),2)))-
 pow(SPB(4,3),2)/(complex<T>(9,0)*pow(SPA(2,5),2))-
 pow(SPB(5,2),2)/(complex<T>(9,0)*pow(SPA(3,4),2))-
 pow(SPB(5,4),2)/(complex<T>(9,0)*pow(SPA(2,3),2))+
 (complex<T>(1,0)*S(2,5)*S(3,4))/(complex<T>(9,0)*SPA(2,3)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))+(complex<T>(1,0)*S(2,3)*S(4,5))/(complex<T>(9,0)*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5))-(pow(SPB(3,2),2)*SPA(2,3))/
(complex<T>(18,0)*SPA(2,5)*SPA(3,4)*SPA(4,5))-
 (pow(SPB(5,2),2)*SPA(2,5))/(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*SPA(4,5))-
 (pow(SPB(4,3),2)*SPA(3,4))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*SPA(4,5))-
 (pow(SPB(5,4),2)*SPA(4,5))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*SPA(3,4))+
 (complex<T>(2,0)*SPB(3,2)*SPB(4,3))/(complex<T>(9,0)*SPA(2,5)*SPA(4,5))+
 (complex<T>(2,0)*SPB(3,2)*SPB(5,2))/(complex<T>(9,0)*SPA(3,4)*SPA(4,5))+
 (complex<T>(2,0)*SPB(4,3)*SPB(5,4))/(complex<T>(9,0)*SPA(2,3)*SPA(2,5))+
 (complex<T>(2,0)*SPB(5,2)*SPB(5,4))/(complex<T>(9,0)*SPA(2,3)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*S(2,4))/
(complex<T>(9,0)*pow(SPA(2,3),2)*pow(SPA(3,4),2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*S(3,5))/
(complex<T>(9,0)*pow(SPA(2,3),2)*pow(SPA(2,5),2)*SS(2,3,5))+
 (complex<T>(1,0)*SS(2,3,4)*SS(2,3,5))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*S(2,4))/
(complex<T>(9,0)*pow(SPA(2,5),2)*pow(SPA(4,5),2)*SS(2,4,5))+
 (complex<T>(1,0)*SS(2,3,5)*SS(2,4,5))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*S(3,5))/
(complex<T>(9,0)*pow(SPA(3,4),2)*pow(SPA(4,5),2)*SS(3,4,5))+
 (complex<T>(1,0)*SS(2,3,4)*SS(3,4,5))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5))+(complex<T>(1,0)*SS(2,4,5)*SS(3,4,5))/
(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmppp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmppp G");
#endif
 
 return(-((pow(SPA(2,5),2)*pow(SPB(5,4),2))/(complex<T>(2,0)*pow(SPA(3,5),2)*
S(2,3)))-(pow(SPA(2,3),2)*pow(SPB(4,3),2))/
(complex<T>(2,0)*pow(SPA(3,5),2)*S(2,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2))/(complex<T>(2,0)*pow(SPA(3,5),2)*
 (S(2,3)+S(3,5)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(2,5)+S(3,5)))-
 (pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4))+
 (pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(3,4))-
 (pow(SPB(5,3),2)*S(3,4))/(complex<T>(6,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*
 SPB(5,2))-(pow(SPB(5,3),2)*S(4,5))/(complex<T>(6,0)*SPA(3,4)*SPA(4,5)*
 SPB(3,2)*SPB(5,2))-(pow(SPB(5,3),2)*SPB(4,3))/
(complex<T>(6,0)*SPA(4,5)*SPB(3,2)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*SPA(3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPA(3,5))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3))/
(complex<T>(6,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(3,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(3,5))-
 (S(3,4)*S(4,5)*SPB(5,3))/(complex<T>(6,0)*SPA(3,4)*SPA(3,5)*SPA(4,5)*
 SPB(3,2)*SPB(5,2))-(pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4))/
(complex<T>(6,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(4,5))-
 (pow(SPB(5,3),2)*SPB(5,4))/(complex<T>(6,0)*SPA(3,4)*SPB(3,2)*SPB(5,2))+
 (SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(3,4))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(3,5))+
 (SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(3,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(4,5))-
 (SPB(4,3)*SPB(5,3)*SPB(5,4))/(complex<T>(6,0)*SPA(3,5)*SPB(3,2)*
 SPB(5,2))+(pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3))/
(complex<T>(6,0)*S(2,3)*SPA(3,4)*SS(2,3,4))+
 (pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPA(3,4)*SS(2,3,4))+
 (complex<T>(2,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),3))/
(SPA(2,3)*SPA(3,4)*(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*
 SS(2,3,4))+(SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*SPA(3,4)*SS(2,3,4))+
 (SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPA(3,4)*SS(2,3,4))-
 (pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3)*SS(2,3,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPA(3,4))-
 pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)/
(pow(SPA(3,5),2)*SS(2,3,5))-(pow(SPA(2,5),2)*pow(SPB(5,4),2)*
 SPB(5,3))/(complex<T>(6,0)*S(2,3)*SPA(3,5)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3))/
(complex<T>(6,0)*S(2,5)*SPA(3,5)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3))/
(complex<T>(6,0)*(S(2,3)+S(3,5))*SPA(3,5)*SS(2,3,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SPA(3,5)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(mH2,2)*pow(SPB(5,3),4))/(SPB(3,2)*SPB(5,2)*
 (SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*(SPA(2,4)*SPB(3,2)+
SPA(4,5)*SPB(5,3))*SS(2,3,5))+
 (SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*SPA(3,5)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,5)*SPA(3,5)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*(S(2,3)+S(3,5))*SPA(3,5)*SS(2,3,5))+
 (SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SPA(3,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SS(2,3,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*S(2,5)*(S(2,3)+S(3,5)))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SS(2,3,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*S(2,3)*(S(2,5)+S(3,5)))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3)*SS(2,3,5))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPA(3,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3)*SS(2,3,5))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(3,5))+
 (complex<T>(-2,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),3))/
(SPA(2,5)*SPA(4,5)*(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3))*
 SS(2,4,5))+(complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4))/
(complex<T>(6,0)*S(2,5)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4))/
(complex<T>(6,0)*(S(2,4)+S(4,5))*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,5)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*(S(2,4)+S(4,5))*SPA(4,5)*SS(2,4,5))-
 (pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4)*SS(2,4,5))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmppp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2))/
(complex<T>(6,0)*pow(SPA(3,5),2)*S(2,3))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,3),2))/(complex<T>(6,0)*pow(SPA(3,5),2)*
 S(2,5))-(pow(SPA(2,3),2)*pow(SPB(4,3),2))/
(complex<T>(6,0)*pow(SPA(3,5),2)*(S(2,3)+S(3,5)))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2))/(complex<T>(6,0)*pow(SPA(3,5),2)*
 (S(2,5)+S(3,5)))+(complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*
 SPB(4,3))/(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*SPA(3,4))-
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3))/
(complex<T>(18,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(3,4))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*S(3,4))/(complex<T>(18,0)*SPA(3,4)*SPA(4,5)*
 SPB(3,2)*SPB(5,2))+(complex<T>(1,0)*pow(SPB(5,3),2)*S(4,5))/
(complex<T>(18,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPB(4,3))/(complex<T>(18,0)*SPA(4,5)*SPB(3,2)*
 SPB(5,2))-(pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*SPA(3,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*SPA(3,5))+
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3))/
(complex<T>(18,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3))/
(complex<T>(18,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(3,5))+
 (complex<T>(1,0)*S(3,4)*S(4,5)*SPB(5,3))/(complex<T>(18,0)*SPA(3,4)*SPA(3,5)*
 SPA(4,5)*SPB(3,2)*SPB(5,2))+(complex<T>(1,0)*pow(SPA(2,4),2)*
 pow(SPB(4,3),2)*SPB(5,4))/(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*
 SPA(4,5))-(pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4))/
(complex<T>(18,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPB(5,4))/(complex<T>(18,0)*SPA(3,4)*SPB(3,2)*
 SPB(5,2))-(complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3)*
 SPB(5,4))/(complex<T>(18,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(3,4))-
 (SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(3,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(3,5))-
 (SPA(2,4)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(4,5))+
 (complex<T>(1,0)*SPB(4,3)*SPB(5,3)*SPB(5,4))/(complex<T>(18,0)*SPA(3,5)*SPB(3,2)*
 SPB(5,2))-(complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3))/
(complex<T>(18,0)*S(2,3)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,4)+S(3,4))*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*S(2,3)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(3,4))*SPA(3,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(4,3)*SS(2,3,4))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPA(3,4))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2))/
(complex<T>(3,0)*pow(SPA(3,5),2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3))/
(complex<T>(18,0)*S(2,3)*SPA(3,5)*SS(2,3,5))+
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3))/
(complex<T>(18,0)*S(2,5)*SPA(3,5)*SS(2,3,5))+
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3))/
(complex<T>(18,0)*(S(2,3)+S(3,5))*SPA(3,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3))/
(complex<T>(18,0)*(S(2,5)+S(3,5))*SPA(3,5)*SS(2,3,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*S(2,3)*SPA(3,5)*SS(2,3,5))-
 (SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*S(2,5)*SPA(3,5)*SS(2,3,5))-
 (SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*(S(2,3)+S(3,5))*SPA(3,5)*SS(2,3,5))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*(S(2,5)+S(3,5))*SPA(3,5)*SS(2,3,5))-
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SS(2,3,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*S(2,5)*(S(2,3)+S(3,5)))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SS(2,3,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*S(2,3)*(S(2,5)+S(3,5)))-
 (pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPB(5,3)*SS(2,3,5))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPA(3,5))-
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*SPB(5,3)*SS(2,3,5))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(3,5))-
 (pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4))/
(complex<T>(18,0)*S(2,5)*SPA(4,5)*SS(2,4,5))-
 (pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(4,5))*SPA(4,5)*SS(2,4,5))-
 (SPA(2,4)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*S(2,5)*SPA(4,5)*SS(2,4,5))-
 (SPA(2,4)*SPA(2,5)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(4,5))*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(4,3),2)*SPB(5,4)*SS(2,4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpmpp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpmpp G");
#endif
 
 return(-((pow(SPA(3,4),2)*pow(SPB(5,4),2))/(complex<T>(2,0)*pow(SPA(2,4),2)*
S(2,3)))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,3)+S(2,4)))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2))/(complex<T>(2,0)*pow(SPA(2,4),2)*
 S(3,4))+(complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(3,4)))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*SPA(2,4))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPA(2,4))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*S(3,4)*SPA(2,4))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(2,4))-
 (pow(SPB(4,2),2)*S(2,5))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*
 SPB(4,3))-(pow(SPB(4,2),2)*S(4,5))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*
 SPB(3,2)*SPB(4,3))-(pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(2,5))-
 (pow(SPB(4,2),2)*SPB(5,2))/(complex<T>(6,0)*SPA(4,5)*SPB(3,2)*SPB(4,3))-
 (S(4,5)*SPB(4,2)*SPB(5,2))/(complex<T>(6,0)*SPA(2,4)*SPA(4,5)*SPB(3,2)*
 SPB(4,3))-(pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4))/
(complex<T>(6,0)*S(3,4)*(S(3,5)+S(4,5))*SPA(4,5))-
 (pow(SPB(4,2),2)*SPB(5,4))/(complex<T>(6,0)*SPA(2,5)*SPB(3,2)*SPB(4,3))-
 (S(2,5)*SPB(4,2)*SPB(5,4))/(complex<T>(6,0)*SPA(2,4)*SPA(2,5)*SPB(3,2)*
 SPB(4,3))+(SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*S(3,4)*SPA(2,4))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(2,4))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(2,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(3,4)*(S(3,5)+S(4,5))*SPA(4,5))-
 pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)/
(pow(SPA(2,4),2)*SS(2,3,4))-(complex<T>(1,0)*pow(SPA(3,4),2)*
 pow(SPB(5,4),2)*SPB(4,2))/(complex<T>(6,0)*S(2,3)*SPA(2,4)*SS(2,3,4))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*SPA(2,4)*SS(2,3,4))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(6,0)*S(3,4)*SPA(2,4)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPA(2,4)*SS(2,3,4))+
 (complex<T>(2,0)*pow(mH2,2)*pow(SPB(4,2),4))/
(SPB(3,2)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(4,3)*
 (SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*SPA(2,4)*SS(2,3,4))+
 (SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*SPA(2,4)*SS(2,3,4))+
 (SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(3,4)*SPA(2,4)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPA(2,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SS(2,3,4))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,3)+S(2,4))*S(3,4))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SS(2,3,4))/
(complex<T>(2,0)*pow(SPA(2,4),2)*S(2,3)*(S(2,4)+S(3,4)))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2)*SS(2,3,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPA(2,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2)*SS(2,3,4))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPA(2,4))+
 (complex<T>(-2,0)*pow(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4),3))/
(SPA(2,3)*SPA(2,5)*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*
 SS(2,3,5))+(complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2))/
(complex<T>(6,0)*S(2,3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SPA(2,5)*SS(2,3,5))-
 (pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2)*SS(2,3,5))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPA(2,5))+
 (complex<T>(2,0)*pow(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2),3))/
(SPA(3,4)*SPA(4,5)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*
 SS(3,4,5))+(complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4))/
(complex<T>(6,0)*S(3,4)*SPA(4,5)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4))/
(complex<T>(6,0)*(S(3,5)+S(4,5))*SPA(4,5)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*S(3,4)*SPA(4,5)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*(S(3,5)+S(4,5))*SPA(4,5)*SS(3,4,5))-
 (pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4)*SS(3,4,5))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpmpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpmpp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2))/
(complex<T>(6,0)*pow(SPA(2,4),2)*S(2,3))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2))/(complex<T>(6,0)*pow(SPA(2,4),2)*
 (S(2,3)+S(2,4)))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2))/
(complex<T>(6,0)*pow(SPA(2,4),2)*S(3,4))-
 (pow(SPA(3,4),2)*pow(SPB(5,4),2))/(complex<T>(6,0)*pow(SPA(2,4),2)*
 (S(2,4)+S(3,4)))-(pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*SPA(2,4))-
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*SPA(2,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*S(3,4)*SPA(2,4))+
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2))/
(complex<T>(18,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(2,4))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*S(2,5))/(complex<T>(18,0)*SPA(2,5)*SPA(4,5)*
 SPB(3,2)*SPB(4,3))+(complex<T>(1,0)*pow(SPB(4,2),2)*S(4,5))/
(complex<T>(18,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*SPA(2,5))-
 (pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2))/
(complex<T>(18,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(2,5))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPB(5,2))/(complex<T>(18,0)*SPA(4,5)*SPB(3,2)*
 SPB(4,3))+(complex<T>(1,0)*S(4,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(18,0)*SPA(2,4)*SPA(4,5)*SPB(3,2)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*SPA(4,5))-
 (pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4))/
(complex<T>(18,0)*S(3,4)*(S(3,5)+S(4,5))*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPB(5,4))/(complex<T>(18,0)*SPA(2,5)*SPB(3,2)*
 SPB(4,3))+(complex<T>(1,0)*S(2,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(18,0)*SPA(2,4)*SPA(2,5)*SPB(3,2)*SPB(4,3))-
 (complex<T>(1,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*S(3,4)*SPA(2,4))-
 (SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(2,4))-
 (SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(2,5))-
 (SPA(3,4)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(3,4)*(S(3,5)+S(4,5))*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2))/
(complex<T>(3,0)*pow(SPA(2,4),2)*SS(2,3,4))+
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2))/
(complex<T>(18,0)*S(2,3)*SPA(2,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*SPA(2,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(18,0)*S(3,4)*SPA(2,4)*SS(2,3,4))+
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2))/
(complex<T>(18,0)*(S(2,4)+S(3,4))*SPA(2,4)*SS(2,3,4))-
 (SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(2,3)*SPA(2,4)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*SPA(2,4)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(3,4)*SPA(2,4)*SS(2,3,4))-
 (SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,4)+S(3,4))*SPA(2,4)*SS(2,3,4))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SS(2,3,4))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,3)+S(2,4))*S(3,4))-
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SS(2,3,4))/
(complex<T>(6,0)*pow(SPA(2,4),2)*S(2,3)*(S(2,4)+S(3,4)))-
 (pow(SPA(3,4),2)*pow(SPB(5,4),2)*SPB(4,2)*SS(2,3,4))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPA(2,4))-
 (pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(4,2)*SS(2,3,4))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPA(2,4))-
 (pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2))/
(complex<T>(18,0)*S(2,3)*SPA(2,5)*SS(2,3,5))-
 (pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2))/
(complex<T>(18,0)*(S(2,5)+S(3,5))*SPA(2,5)*SS(2,3,5))-
 (SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(2,3)*SPA(2,5)*SS(2,3,5))-
 (SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*(S(2,5)+S(3,5))*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,4),2)*SPB(5,2)*SS(2,3,5))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPA(2,5))-
 (pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4))/
(complex<T>(18,0)*S(3,4)*SPA(4,5)*SS(3,4,5))-
 (pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4))/
(complex<T>(18,0)*(S(3,5)+S(4,5))*SPA(4,5)*SS(3,4,5))-
 (SPA(3,4)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*S(3,4)*SPA(4,5)*SS(3,4,5))-
 (SPA(3,4)*SPA(3,5)*SPB(4,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*(S(3,5)+S(4,5))*SPA(4,5)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(5,4)*SS(3,4,5))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdppmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdppmp G");
#endif
 
 return(-((pow(SPA(4,5),2)*pow(SPB(5,2),2))/(complex<T>(2,0)*pow(SPA(3,5),2)*
S(3,4)))+(complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(3,4)+S(3,5)))-
 (pow(SPA(3,4),2)*pow(SPB(3,2),2))/(complex<T>(2,0)*pow(SPA(3,5),2)*
 S(4,5))+(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(3,5)+S(4,5)))-
 (pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*S(3,4)*SPA(2,3))-
 (pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*S(4,5)*SPA(2,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*SPA(3,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*SPA(3,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*S(4,5)*SPA(3,5))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3))/
(complex<T>(6,0)*S(3,4)*(S(3,5)+S(4,5))*SPA(3,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*S(3,4)*SPA(2,3))+
 (complex<T>(1,0)*SPA(2,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*S(4,5)*SPA(2,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*S(4,5)*SPA(3,5))+
 (SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*S(3,4)*(S(3,5)+S(4,5))*SPA(3,5))-
 (pow(SPB(5,3),2)*S(2,3))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*
 SPB(5,4))-(pow(SPB(5,3),2)*S(2,5))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*
 SPB(4,3)*SPB(5,4))-(pow(SPB(5,3),2)*SPB(3,2))/
(complex<T>(6,0)*SPA(2,5)*SPB(4,3)*SPB(5,4))-(pow(SPB(5,3),2)*SPB(5,2))/
(complex<T>(6,0)*SPA(2,3)*SPB(4,3)*SPB(5,4))-(S(2,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(6,0)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))-
 (S(2,3)*SPB(5,2)*SPB(5,3))/(complex<T>(6,0)*SPA(2,3)*SPA(3,5)*SPB(4,3)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2))/
(complex<T>(6,0)*S(3,4)*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*(S(2,3)+S(2,4))*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*S(3,4)*SPA(2,3)*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),3))/
(SPA(2,3)*SPA(3,4)*(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4))*
 SS(2,3,4))-(pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2)*
 SS(2,3,4))/(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*SPA(2,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2))/
(complex<T>(6,0)*S(4,5)*SPA(2,5)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*SPA(2,5)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*S(4,5)*SPA(2,5)*SS(2,4,5))+
 (complex<T>(-2,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),3))/
(SPA(2,5)*SPA(4,5)*(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3))*
 SS(2,4,5))-(pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2)*
 SS(2,4,5))/(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPA(2,5))-
 pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)/
(pow(SPA(3,5),2)*SS(3,4,5))-(pow(SPA(4,5),2)*pow(SPB(5,2),2)*
 SPB(5,3))/(complex<T>(6,0)*S(3,4)*SPA(3,5)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*SPA(3,5)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3))/
(complex<T>(6,0)*S(4,5)*SPA(3,5)*SS(3,4,5))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3))/
(complex<T>(6,0)*(S(3,5)+S(4,5))*SPA(3,5)*SS(3,4,5))+
 (SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*S(3,4)*SPA(3,5)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*SPA(3,5)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*S(4,5)*SPA(3,5)*SS(3,4,5))+
 (SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*(S(3,5)+S(4,5))*SPA(3,5)*SS(3,4,5))+
 (complex<T>(2,0)*pow(mH2,2)*pow(SPB(5,3),4))/
(SPB(4,3)*(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3))*SPB(5,4)*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4))*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SS(3,4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*(S(3,4)+S(3,5))*S(4,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SS(3,4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*S(3,4)*(S(3,5)+S(4,5)))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3)*SS(3,4,5))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPA(3,5))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3)*SS(3,4,5))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPA(3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdppmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdppmp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2))/
(complex<T>(6,0)*pow(SPA(3,5),2)*S(3,4))-
 (pow(SPA(3,4),2)*pow(SPB(3,2),2))/(complex<T>(6,0)*pow(SPA(3,5),2)*
 (S(3,4)+S(3,5)))+(complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(3,2),2))/
(complex<T>(6,0)*pow(SPA(3,5),2)*S(4,5))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2))/(complex<T>(6,0)*pow(SPA(3,5),2)*
 (S(3,5)+S(4,5)))+(complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*
 SPB(3,2))/(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*SPA(2,3))-
 (pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*S(3,4)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*SPA(2,5))-
 (pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*S(4,5)*SPA(2,5))-
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*SPA(3,5))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*SPA(3,5))+
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*S(4,5)*SPA(3,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3))/
(complex<T>(18,0)*S(3,4)*(S(3,5)+S(4,5))*SPA(3,5))-
 (SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*S(3,4)*SPA(2,3))-
 (SPA(2,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*S(4,5)*SPA(2,5))-
 (SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*S(4,5)*SPA(3,5))-
 (complex<T>(1,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*S(3,4)*(S(3,5)+S(4,5))*SPA(3,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*S(2,3))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPB(4,3)*SPB(5,4))+(complex<T>(1,0)*pow(SPB(5,3),2)*S(2,5))/
(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPB(3,2))/(complex<T>(18,0)*SPA(2,5)*SPB(4,3)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPB(5,3),2)*SPB(5,2))/
(complex<T>(18,0)*SPA(2,3)*SPB(4,3)*SPB(5,4))+
 (complex<T>(1,0)*S(2,5)*SPB(3,2)*SPB(5,3))/(complex<T>(18,0)*SPA(2,5)*SPA(3,5)*
 SPB(4,3)*SPB(5,4))+(complex<T>(1,0)*S(2,3)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*SPA(2,3)*SPA(3,5)*SPB(4,3)*SPB(5,4))-
 (pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*SPA(2,3)*SS(2,3,4))-
 (pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2))/
(complex<T>(18,0)*S(3,4)*SPA(2,3)*SS(2,3,4))-
 (SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*(S(2,3)+S(2,4))*SPA(2,3)*SS(2,3,4))-
 (SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*S(3,4)*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(3,2)*SS(2,3,4))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPA(2,3))-
 (pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*SPA(2,5)*SS(2,4,5))-
 (pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2))/
(complex<T>(18,0)*S(4,5)*SPA(2,5)*SS(2,4,5))-
 (SPA(2,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*SPA(2,5)*SS(2,4,5))-
 (SPA(2,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*S(4,5)*SPA(2,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,2)*SS(2,4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPA(2,5))+
 (complex<T>(1,0)*pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2))/
(complex<T>(3,0)*pow(SPA(3,5),2)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3))/
(complex<T>(18,0)*S(3,4)*SPA(3,5)*SS(3,4,5))+
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*SPA(3,5)*SS(3,4,5))+
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3))/
(complex<T>(18,0)*S(4,5)*SPA(3,5)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3))/
(complex<T>(18,0)*(S(3,5)+S(4,5))*SPA(3,5)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*S(3,4)*SPA(3,5)*SS(3,4,5))-
 (SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*SPA(3,5)*SS(3,4,5))-
 (SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*S(4,5)*SPA(3,5)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*(S(3,5)+S(4,5))*SPA(3,5)*SS(3,4,5))-
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SS(3,4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*(S(3,4)+S(3,5))*S(4,5))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SS(3,4,5))/
(complex<T>(6,0)*pow(SPA(3,5),2)*S(3,4)*(S(3,5)+S(4,5)))-
 (pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3)*SS(3,4,5))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPA(3,5))-
 (pow(SPA(3,4),2)*pow(SPB(3,2),2)*SPB(5,3)*SS(3,4,5))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPA(3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpppm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpppm G");
#endif
 
 return(-((pow(SPA(4,5),2)*pow(SPB(4,3),2))/(complex<T>(2,0)*pow(SPA(2,4),2)*
S(2,5)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(2,5)))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2))/(complex<T>(2,0)*pow(SPA(2,4),2)*
 S(4,5))+(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(4,5)))-
 (pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2))/
(complex<T>(6,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*SPA(2,4))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*SPA(2,4))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*S(4,5)*SPA(2,4))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2))/
(complex<T>(6,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(2,4))-
 (pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4))+
 (pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*S(4,5)*SPA(3,4))+
 (complex<T>(1,0)*SPA(2,5)*SPA(3,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(2,3))+
 (SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*S(4,5)*SPA(2,4))+
 (complex<T>(1,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(2,4))+
 (SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*S(4,5)*SPA(3,4))-
 (pow(SPB(4,2),2)*S(2,3))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*SPB(5,2)*
 SPB(5,4))-(pow(SPB(4,2),2)*S(3,4))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*
 SPB(5,2)*SPB(5,4))-(pow(SPB(4,2),2)*SPB(3,2))/
(complex<T>(6,0)*SPA(3,4)*SPB(5,2)*SPB(5,4))-(S(2,3)*S(3,4)*SPB(4,2))/
(complex<T>(6,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))-
 (pow(SPB(4,2),2)*SPB(4,3))/(complex<T>(6,0)*SPA(2,3)*SPB(5,2)*SPB(5,4))-
 (SPB(3,2)*SPB(4,2)*SPB(4,3))/(complex<T>(6,0)*SPA(2,4)*SPB(5,2)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2))/
(complex<T>(6,0)*S(2,5)*SPA(2,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2))/
(complex<T>(6,0)*(S(2,3)+S(3,5))*SPA(2,3)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,5)*SPA(3,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*S(2,5)*SPA(2,3)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,5)*SPA(3,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,3)+S(3,5))*SPA(2,3)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3))/
(SPA(2,3)*SPA(2,5)*(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4))*
 SS(2,3,5))-(pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2)*
 SS(2,3,5))/(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(2,3))-
 pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)/
(pow(SPA(2,4),2)*SS(2,4,5))-(complex<T>(1,0)*pow(SPA(4,5),2)*
 pow(SPB(4,3),2)*SPB(4,2))/(complex<T>(6,0)*S(2,5)*SPA(2,4)*SS(2,4,5))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*SPA(2,4)*SS(2,4,5))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2))/
(complex<T>(6,0)*S(4,5)*SPA(2,4)*SS(2,4,5))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2))/
(complex<T>(6,0)*(S(2,4)+S(4,5))*SPA(2,4)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*S(2,5)*SPA(2,4)*SS(2,4,5))+
 (SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,4)+S(2,5))*SPA(2,4)*SS(2,4,5))+
 (SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*S(4,5)*SPA(2,4)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*(S(2,4)+S(4,5))*SPA(2,4)*SS(2,4,5))+
 (complex<T>(-2,0)*pow(mH2,2)*pow(SPB(4,2),4))/
(SPB(5,2)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*SPB(5,4)*
 (SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4))*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SS(2,4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*(S(2,4)+S(2,5))*S(4,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SS(2,4,5))/
(complex<T>(2,0)*pow(SPA(2,4),2)*S(2,5)*(S(2,4)+S(4,5)))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2)*SS(2,4,5))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(2,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2)*SS(2,4,5))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPA(2,4))+
 (pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*SPA(3,4)*SS(3,4,5))+
 (pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3))/
(complex<T>(6,0)*S(4,5)*SPA(3,4)*SS(3,4,5))+
 (SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*(S(3,4)+S(3,5))*SPA(3,4)*SS(3,4,5))+
 (SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(6,0)*S(4,5)*SPA(3,4)*SS(3,4,5))+
 (complex<T>(2,0)*pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),3))/
(SPA(3,4)*SPA(4,5)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*
 SS(3,4,5))-(pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3)*
 SS(3,4,5))/(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpppm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2))/
(complex<T>(6,0)*pow(SPA(2,4),2)*S(2,5))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2))/(complex<T>(6,0)*pow(SPA(2,4),2)*
 (S(2,4)+S(2,5)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2))/
(complex<T>(6,0)*pow(SPA(2,4),2)*S(4,5))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2))/(complex<T>(6,0)*pow(SPA(2,4),2)*
 (S(2,4)+S(4,5)))+(complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*
 SPB(3,2))/(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*SPA(2,3))-
 (pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2))/
(complex<T>(18,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(2,3))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*SPA(2,4))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*SPA(2,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*S(4,5)*SPA(2,4))+
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2))/
(complex<T>(18,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(2,4))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*SPA(3,4))-
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*S(4,5)*SPA(3,4))-
 (SPA(2,5)*SPA(3,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(2,3))-
 (complex<T>(1,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*S(4,5)*SPA(2,4))-
 (SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(2,4))-
 (complex<T>(1,0)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*S(4,5)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*S(2,3))/(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*
 SPB(5,2)*SPB(5,4))+(complex<T>(1,0)*pow(SPB(4,2),2)*S(3,4))/
(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPB(3,2))/(complex<T>(18,0)*SPA(3,4)*SPB(5,2)*
 SPB(5,4))+(complex<T>(1,0)*S(2,3)*S(3,4)*SPB(4,2))/
(complex<T>(18,0)*SPA(2,3)*SPA(2,4)*SPA(3,4)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPB(4,3))/(complex<T>(18,0)*SPA(2,3)*SPB(5,2)*
 SPB(5,4))+(complex<T>(1,0)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*SPA(2,4)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2))/
(complex<T>(18,0)*S(2,5)*SPA(2,3)*SS(2,3,5))-
 (pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2))/
(complex<T>(18,0)*(S(2,3)+S(3,5))*SPA(2,3)*SS(2,3,5))-
 (SPA(2,5)*SPA(3,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*S(2,5)*SPA(2,3)*SS(2,3,5))-
 (SPA(2,5)*SPA(3,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,3)+S(3,5))*SPA(2,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*SPB(3,2)*SS(2,3,5))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*S(2,5)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2))/
(complex<T>(3,0)*pow(SPA(2,4),2)*SS(2,4,5))+
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2))/
(complex<T>(18,0)*S(2,5)*SPA(2,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*SPA(2,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2))/
(complex<T>(18,0)*S(4,5)*SPA(2,4)*SS(2,4,5))+
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2))/
(complex<T>(18,0)*(S(2,4)+S(4,5))*SPA(2,4)*SS(2,4,5))-
 (SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*S(2,5)*SPA(2,4)*SS(2,4,5))-
 (complex<T>(1,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,4)+S(2,5))*SPA(2,4)*SS(2,4,5))-
 (complex<T>(1,0)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*S(4,5)*SPA(2,4)*SS(2,4,5))-
 (SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*(S(2,4)+S(4,5))*SPA(2,4)*SS(2,4,5))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SS(2,4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*(S(2,4)+S(2,5))*S(4,5))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SS(2,4,5))/
(complex<T>(6,0)*pow(SPA(2,4),2)*S(2,5)*(S(2,4)+S(4,5)))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2)*SS(2,4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPA(2,4))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*SPB(4,2)*SS(2,4,5))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPA(2,4))-
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*SPA(3,4)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3))/
(complex<T>(18,0)*S(4,5)*SPA(3,4)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*(S(3,4)+S(3,5))*SPA(3,4)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(18,0)*S(4,5)*SPA(3,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPB(4,3)*SS(3,4,5))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmmpp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmmpp G");
#endif
 
 return(-(pow(SPB(5,4),2)/(complex<T>(3,0)*pow(SPB(3,2),2)))+
 (complex<T>(-2,0)*pow(SPA(2,3),3))/(SPA(2,5)*SPA(3,4)*SPA(4,5))-
 (pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(6,0)*S(2,5)*SPB(3,2))-
 (pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(6,0)*S(3,4)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(4,2)*SPB(5,3)+S(2,4)*SPB(5,4),2)*
 SPA(2,3))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(3,2))-
 (S(2,4)*SPA(2,3)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(2,5)*
 SPA(4,5))-(S(2,5)*SPA(2,3)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(3,2),2)*
 SPA(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(4,2)*SPB(5,3)+S(3,5)*SPB(5,4),2)*
 SPA(2,3))/(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPA(2,5)*SPB(3,2)*
 SPB(5,2))+(complex<T>(2,0)*pow(SPB(5,4),3))/(SPB(3,2)*SPB(4,3)*
 SPB(5,2))-(S(3,4)*SPA(2,3)*SPB(5,2))/(complex<T>(3,0)*pow(SPB(3,2),2)*
 SPA(3,4)*SPA(4,5))-(S(3,5)*SPA(2,3)*SPB(5,2))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,3))/
(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPA(2,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,3))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(3,2))-
 (SPA(2,3)*SPB(5,4))/(complex<T>(3,0)*SPA(4,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*S(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPA(2,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*S(2,4)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(3,2))+
 (pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(6,0)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(4,2)*SPB(5,3)+S(2,4)*SPB(5,4),2)*
 SPA(2,3))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*SPB(4,3))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(3,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,3))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPB(3,2)*SS(2,3,4))-
 (pow(SPA(2,3),2)*SPB(5,2)*SPB(5,4))/(complex<T>(6,0)*SPA(3,4)*SPB(3,2)*
 SS(2,3,4))+(complex<T>(1,0)*pow(SPA(2,3),2)*S(2,4)*SPB(4,3)*SPB(5,2)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,3)+S(2,4),2)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(6,0)*SPB(3,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(4,2)*SPB(5,3)+S(3,5)*SPB(5,4),2)*
 SPA(2,3))/(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPB(3,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPB(5,2))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,3))/(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPB(3,2)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*SPA(2,5)*SPB(3,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*S(3,5)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(S(2,3)+S(3,5),2)*SPB(3,2)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmmpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmmpp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(5,4),2))/(complex<T>(9,0)*pow(SPB(3,2),2))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(18,0)*S(2,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(18,0)*S(3,4)*SPB(3,2))-
 (pow(SPA(2,3)*SPB(4,2)*SPB(5,3)+S(2,4)*SPB(5,4),2)*SPA(2,3))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(3,2))+
 (complex<T>(1,0)*S(2,4)*SPA(2,3)*SPB(4,3))/(complex<T>(9,0)*pow(SPB(3,2),2)*
 SPA(2,5)*SPA(4,5))+(complex<T>(1,0)*S(2,5)*SPA(2,3)*SPB(4,3))/
(complex<T>(9,0)*pow(SPB(3,2),2)*SPA(2,5)*SPA(4,5))-
 (pow(SPA(2,3)*SPB(4,2)*SPB(5,3)+S(3,5)*SPB(5,4),2)*SPA(2,3))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPA(2,5)*SPB(3,2)*SPB(5,2))+
 (complex<T>(1,0)*S(3,4)*SPA(2,3)*SPB(5,2))/(complex<T>(9,0)*pow(SPB(3,2),2)*
 SPA(3,4)*SPA(4,5))+(complex<T>(1,0)*S(3,5)*SPA(2,3)*SPB(5,2))/
(complex<T>(9,0)*pow(SPB(3,2),2)*SPA(3,4)*SPA(4,5))-
 (pow(SPA(2,3),2)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,3))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPA(2,5)*SPB(3,2))-
 (pow(SPA(2,3),2)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(3,2))+
 (complex<T>(1,0)*SPA(2,3)*SPB(5,4))/(complex<T>(9,0)*SPA(4,5)*SPB(3,2))-
 (pow(SPA(2,3),2)*S(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPA(2,5)*SPB(3,2))-
 (pow(SPA(2,3),2)*S(2,4)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*S(3,4)*SPB(3,2))-
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(18,0)*SPB(3,2)*SS(2,3,4))-
 (pow(SPA(2,3)*SPB(4,2)*SPB(5,3)+S(2,4)*SPB(5,4),2)*SPA(2,3))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPB(3,2)*SS(2,3,4))-
 (pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*SPB(4,3))/
(complex<T>(9,0)*pow(SPB(3,2),2)*SPA(3,4)*SS(2,3,4))-
 (pow(SPA(2,3),2)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*SPA(3,4)*SPB(3,2)*SS(2,3,4))-
 (pow(SPA(2,3),2)*S(2,4)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,3)+S(2,4),2)*SPB(3,2)*SS(2,3,4))-
 (pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(18,0)*SPB(3,2)*SS(2,3,5))-
 (pow(SPA(2,3)*SPB(4,2)*SPB(5,3)+S(3,5)*SPB(5,4),2)*SPA(2,3))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPB(3,2)*SS(2,3,5))-
 (pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPB(5,2))/
(complex<T>(9,0)*pow(SPB(3,2),2)*SPA(2,5)*SS(2,3,5))-
 (pow(SPA(2,3),2)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPB(3,2)*SS(2,3,5))+
 (pow(SPA(2,3),2)*SPB(4,3)*SPB(5,4))/(complex<T>(18,0)*SPA(2,5)*SPB(3,2)*
 SS(2,3,5))-(pow(SPA(2,3),2)*S(3,5)*SPB(4,3)*SPB(5,2)*
 SPB(5,4))/(complex<T>(18,0)*pow(S(2,3)+S(3,5),2)*SPB(3,2)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpmmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpmmp G");
#endif
 
 return(-(pow(SPB(5,2),2)/(complex<T>(3,0)*pow(SPB(4,3),2)))+
 (complex<T>(-2,0)*pow(SPA(3,4),3))/(SPA(2,3)*SPA(2,5)*SPA(4,5))-
 (S(2,4)*SPA(3,4)*SPB(3,2))/(complex<T>(3,0)*pow(SPB(4,3),2)*SPA(2,5)*
 SPA(4,5))-(S(4,5)*SPA(3,4)*SPB(3,2))/(complex<T>(3,0)*pow(SPB(4,3),2)*
 SPA(2,5)*SPA(4,5))-(pow(SPB(5,2),2)*SPA(3,4))/
(complex<T>(6,0)*S(2,3)*SPB(4,3))+
 (complex<T>(1,0)*pow(-(S(2,4)*SPB(5,2))-SPA(3,4)*SPB(4,2)*SPB(5,3),2)*
 SPA(3,4))/(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,3))-
 (pow(SPB(5,2),2)*SPA(3,4))/(complex<T>(6,0)*S(4,5)*SPB(4,3))+
 (complex<T>(1,0)*pow(-(S(3,5)*SPB(5,2))-SPA(3,4)*SPB(4,2)*SPB(5,3),2)*
 SPA(3,4))/(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(4,3))-
 (SPA(3,4)*SPB(5,2))/(complex<T>(3,0)*SPA(2,5)*SPB(4,3))+
 (complex<T>(2,0)*pow(SPB(5,2),3))/(SPB(3,2)*SPB(4,3)*SPB(5,4))-
 (S(2,3)*SPA(3,4)*SPB(5,4))/(complex<T>(3,0)*pow(SPB(4,3),2)*SPA(2,3)*
 SPA(2,5))-(S(3,5)*SPA(3,4)*SPB(5,4))/(complex<T>(3,0)*pow(SPB(4,3),2)*
 SPA(2,3)*SPA(2,5))+(pow(SPA(3,4),2)*S(2,4)*SPB(3,2)*SPB(5,2)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,3))+
 (pow(SPA(3,4),2)*S(3,5)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(4,3))+
 (pow(SPA(3,4),2)*SPA(3,4)*SPB(3,2)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,3))+
 (pow(SPA(3,4),2)*SPA(3,4)*SPB(3,2)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*SPB(3,2))/
(complex<T>(3,0)*pow(SPB(4,3),2)*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*SPA(3,4))/(complex<T>(6,0)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(-(S(2,4)*SPB(5,2))-SPA(3,4)*SPB(4,2)*SPB(5,3),2)*
 SPA(3,4))/(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPB(4,3)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*SPA(2,3)*SPB(4,3)*SS(2,3,4))+
 (pow(SPA(3,4),2)*S(2,4)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPB(4,3)*SS(2,3,4))+
 (pow(SPA(3,4),2)*SPA(3,4)*SPB(3,2)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*pow(S(2,4)+S(3,4),2)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*SPA(3,4))/(complex<T>(6,0)*SPB(4,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(-(S(3,5)*SPB(5,2))-SPA(3,4)*SPB(4,2)*SPB(5,3),2)*
 SPA(3,4))/(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPB(4,3)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPB(3,2)*SPB(5,2))/
(complex<T>(6,0)*SPA(4,5)*SPB(4,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*SPB(5,4))/
(complex<T>(3,0)*pow(SPB(4,3),2)*SPA(4,5)*SS(3,4,5))+
 (pow(SPA(3,4),2)*S(3,5)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPB(4,3)*SS(3,4,5))+
 (pow(SPA(3,4),2)*SPA(3,4)*SPB(3,2)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*pow(S(3,4)+S(3,5),2)*SPB(4,3)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpmmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpmmp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(5,2),2))/(complex<T>(9,0)*pow(SPB(4,3),2))+
 (complex<T>(1,0)*S(2,4)*SPA(3,4)*SPB(3,2))/(complex<T>(9,0)*pow(SPB(4,3),2)*
 SPA(2,5)*SPA(4,5))+(complex<T>(1,0)*S(4,5)*SPA(3,4)*SPB(3,2))/
(complex<T>(9,0)*pow(SPB(4,3),2)*SPA(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*SPA(3,4))/(complex<T>(18,0)*S(2,3)*SPB(4,3))-
 (pow(-(S(2,4)*SPB(5,2))-SPA(3,4)*SPB(4,2)*SPB(5,3),2)*SPA(3,4))/
(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*SPA(3,4))/(complex<T>(18,0)*S(4,5)*SPB(4,3))-
 (pow(-(S(3,5)*SPB(5,2))-SPA(3,4)*SPB(4,2)*SPB(5,3),2)*SPA(3,4))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(4,3))+
 (complex<T>(1,0)*SPA(3,4)*SPB(5,2))/(complex<T>(9,0)*SPA(2,5)*SPB(4,3))+
 (complex<T>(1,0)*S(2,3)*SPA(3,4)*SPB(5,4))/(complex<T>(9,0)*pow(SPB(4,3),2)*
 SPA(2,3)*SPA(2,5))+(complex<T>(1,0)*S(3,5)*SPA(3,4)*SPB(5,4))/
(complex<T>(9,0)*pow(SPB(4,3),2)*SPA(2,3)*SPA(2,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*S(2,4)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,3))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*S(3,5)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(4,3))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPA(3,4)*SPB(3,2)*SPB(4,2)*SPB(5,3)*
 SPB(5,4))/(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,3))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPA(3,4)*SPB(3,2)*SPB(4,2)*SPB(5,3)*
 SPB(5,4))/(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*S(4,5)*SPB(4,3))-
 (pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)*SPB(3,2))/
(complex<T>(9,0)*pow(SPB(4,3),2)*SPA(2,3)*SS(2,3,4))-
 (pow(SPB(5,2),2)*SPA(3,4))/(complex<T>(18,0)*SPB(4,3)*SS(2,3,4))-
 (pow(-(S(2,4)*SPB(5,2))-SPA(3,4)*SPB(4,2)*SPB(5,3),2)*SPA(3,4))/
(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*SPB(4,3)*SS(2,3,4))+
 (pow(SPA(3,4),2)*SPB(5,2)*SPB(5,4))/(complex<T>(18,0)*SPA(2,3)*SPB(4,3)*
 SS(2,3,4))-(complex<T>(1,0)*pow(SPA(3,4),2)*S(2,4)*SPB(3,2)*SPB(5,2)*
 SPB(5,4))/(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*SPB(4,3)*
 SS(2,3,4))-(complex<T>(1,0)*pow(SPA(3,4),2)*SPA(3,4)*SPB(3,2)*
 SPB(4,2)*SPB(5,3)*SPB(5,4))/(complex<T>(18,0)*pow(S(2,4)+S(3,4),2)*
 SPB(4,3)*SS(2,3,4))-(pow(SPB(5,2),2)*SPA(3,4))/
(complex<T>(18,0)*SPB(4,3)*SS(3,4,5))-
 (pow(-(S(3,5)*SPB(5,2))-SPA(3,4)*SPB(4,2)*SPB(5,3),2)*SPA(3,4))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPB(4,3)*SS(3,4,5))+
 (pow(SPA(3,4),2)*SPB(3,2)*SPB(5,2))/(complex<T>(18,0)*SPA(4,5)*SPB(4,3)*
 SS(3,4,5))-(pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*
 SPB(5,4))/(complex<T>(9,0)*pow(SPB(4,3),2)*SPA(4,5)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*S(3,5)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPB(4,3)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(3,4),2)*SPA(3,4)*SPB(3,2)*SPB(4,2)*SPB(5,3)*
 SPB(5,4))/(complex<T>(18,0)*pow(S(3,4)+S(3,5),2)*SPB(4,3)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdppmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdppmm G");
#endif
 
 return(-(pow(SPB(3,2),2)/(complex<T>(3,0)*pow(SPB(5,4),2)))+
 (complex<T>(-2,0)*pow(SPA(4,5),3))/(SPA(2,3)*SPA(2,5)*SPA(3,4))-
 (S(2,5)*SPA(4,5)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*SPA(2,3)*
 SPA(2,5))-(S(3,5)*SPA(4,5)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*
 SPA(2,3)*SPA(2,5))-(S(2,4)*SPA(4,5)*SPB(5,2))/
(complex<T>(3,0)*pow(SPB(5,4),2)*SPA(2,3)*SPA(3,4))-
 (S(3,4)*SPA(4,5)*SPB(5,2))/(complex<T>(3,0)*pow(SPB(5,4),2)*SPA(2,3)*
 SPA(3,4))-(pow(SPB(3,2),2)*SPA(4,5))/(complex<T>(6,0)*S(2,5)*
 SPB(5,4))-(pow(SPB(3,2),2)*SPA(4,5))/(complex<T>(6,0)*S(3,4)*
 SPB(5,4))+(complex<T>(1,0)*pow(S(3,5)*SPB(3,2)+SPA(4,5)*SPB(4,2)*
  SPB(5,3),2)*SPA(4,5))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*
 SPB(5,4))-(SPA(4,5)*SPB(3,2))/(complex<T>(3,0)*SPA(2,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*S(2,4)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPA(2,5)*SPB(5,4))+
 (complex<T>(1,0)*pow(S(2,4)*SPB(3,2)+SPA(4,5)*SPB(4,2)*SPB(5,3),2)*
 SPA(4,5))/(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPA(2,5)*SPB(5,2)*
 SPB(5,4))+(complex<T>(2,0)*pow(SPB(3,2),3))/(SPB(4,3)*SPB(5,2)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(4,5),2)*S(3,5)*SPB(3,2)*SPB(4,3)*
 SPB(5,2))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPA(4,5)*SPB(4,2)*SPB(4,3)*SPB(5,3))/
(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPA(2,5)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPA(4,5)*SPB(4,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,3))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPB(5,2))/
(complex<T>(3,0)*pow(SPB(5,4),2)*SPA(2,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(4,5))/(complex<T>(6,0)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(S(2,4)*SPB(3,2)+SPA(4,5)*SPB(4,2)*SPB(5,3),2)*
 SPA(4,5))/(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPB(5,4)*SS(2,4,5))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*SPA(2,5)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*S(2,4)*SPB(3,2)*SPB(4,3)*SPB(5,2))/
(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPA(4,5)*SPB(4,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,3))/(complex<T>(6,0)*pow(S(2,4)+S(4,5),2)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*SPB(4,3))/
(complex<T>(3,0)*pow(SPB(5,4),2)*SPA(3,4)*SS(3,4,5))+
 (pow(SPB(3,2),2)*SPA(4,5))/(complex<T>(6,0)*SPB(5,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(S(3,5)*SPB(3,2)+SPA(4,5)*SPB(4,2)*SPB(5,3),2)*
 SPA(4,5))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*SPB(5,4)*SS(3,4,5))-
 (pow(SPA(4,5),2)*SPB(3,2)*SPB(5,2))/(complex<T>(6,0)*SPA(3,4)*SPB(5,4)*
 SS(3,4,5))+(complex<T>(1,0)*pow(SPA(4,5),2)*S(3,5)*SPB(3,2)*SPB(4,3)*
 SPB(5,2))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*SPB(5,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPA(4,5)*SPB(4,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,3))/(complex<T>(6,0)*pow(S(3,5)+S(4,5),2)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdppmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdppmm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(3,2),2))/(complex<T>(9,0)*pow(SPB(5,4),2))+
 (complex<T>(1,0)*S(2,5)*SPA(4,5)*SPB(4,3))/(complex<T>(9,0)*pow(SPB(5,4),2)*
 SPA(2,3)*SPA(2,5))+(complex<T>(1,0)*S(3,5)*SPA(4,5)*SPB(4,3))/
(complex<T>(9,0)*pow(SPB(5,4),2)*SPA(2,3)*SPA(2,5))+
 (complex<T>(1,0)*S(2,4)*SPA(4,5)*SPB(5,2))/(complex<T>(9,0)*pow(SPB(5,4),2)*
 SPA(2,3)*SPA(3,4))+(complex<T>(1,0)*S(3,4)*SPA(4,5)*SPB(5,2))/
(complex<T>(9,0)*pow(SPB(5,4),2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(4,5))/(complex<T>(18,0)*S(2,5)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(4,5))/(complex<T>(18,0)*S(3,4)*SPB(5,4))-
 (pow(S(3,5)*SPB(3,2)+SPA(4,5)*SPB(4,2)*SPB(5,3),2)*SPA(4,5))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,4))+
 (complex<T>(1,0)*SPA(4,5)*SPB(3,2))/(complex<T>(9,0)*SPA(2,3)*SPB(5,4))-
 (pow(SPA(4,5),2)*S(2,4)*SPB(3,2)*SPB(4,3))/
(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPA(2,5)*SPB(5,4))-
 (pow(S(2,4)*SPB(3,2)+SPA(4,5)*SPB(4,2)*SPB(5,3),2)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPA(2,5)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(4,5),2)*S(3,5)*SPB(3,2)*SPB(4,3)*SPB(5,2))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,4))-
 (pow(SPA(4,5),2)*SPA(4,5)*SPB(4,2)*SPB(4,3)*SPB(5,3))/
(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPA(2,5)*SPB(5,4))-
 (pow(SPA(4,5),2)*SPA(4,5)*SPB(4,2)*SPB(4,3)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,4))-
 (pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPB(5,2))/
(complex<T>(9,0)*pow(SPB(5,4),2)*SPA(2,5)*SS(2,4,5))-
 (pow(SPB(3,2),2)*SPA(4,5))/(complex<T>(18,0)*SPB(5,4)*SS(2,4,5))-
 (pow(S(2,4)*SPB(3,2)+SPA(4,5)*SPB(4,2)*SPB(5,3),2)*SPA(4,5))/
(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPB(5,4)*SS(2,4,5))+
 (pow(SPA(4,5),2)*SPB(3,2)*SPB(4,3))/(complex<T>(18,0)*SPA(2,5)*SPB(5,4)*
 SS(2,4,5))-(pow(SPA(4,5),2)*S(2,4)*SPB(3,2)*SPB(4,3)*
 SPB(5,2))/(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPB(5,4)*
 SS(2,4,5))-(pow(SPA(4,5),2)*SPA(4,5)*SPB(4,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,3))/(complex<T>(18,0)*pow(S(2,4)+S(4,5),2)*SPB(5,4)*
 SS(2,4,5))-(pow(-(SPA(3,4)*SPB(3,2))+SPA(4,5)*SPB(5,2),2)*
 SPB(4,3))/(complex<T>(9,0)*pow(SPB(5,4),2)*SPA(3,4)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(4,5))/(complex<T>(18,0)*SPB(5,4)*SS(3,4,5))-
 (pow(S(3,5)*SPB(3,2)+SPA(4,5)*SPB(4,2)*SPB(5,3),2)*SPA(4,5))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPB(5,4)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPB(3,2)*SPB(5,2))/
(complex<T>(18,0)*SPA(3,4)*SPB(5,4)*SS(3,4,5))-
 (pow(SPA(4,5),2)*S(3,5)*SPB(3,2)*SPB(4,3)*SPB(5,2))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPB(5,4)*SS(3,4,5))-
 (pow(SPA(4,5),2)*SPA(4,5)*SPB(4,2)*SPB(4,3)*SPB(5,2)*SPB(5,3))/
(complex<T>(18,0)*pow(S(3,5)+S(4,5),2)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmppm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmppm G");
#endif
 
 return(-(pow(SPB(4,3),2)/(complex<T>(3,0)*pow(SPB(5,2),2)))+
 (complex<T>(-2,0)*pow(SPA(2,5),3))/(SPA(2,3)*SPA(3,4)*SPA(4,5))-
 (S(3,5)*SPA(2,5)*SPB(3,2))/(complex<T>(3,0)*pow(SPB(5,2),2)*SPA(3,4)*
 SPA(4,5))-(S(4,5)*SPA(2,5)*SPB(3,2))/(complex<T>(3,0)*pow(SPB(5,2),2)*
 SPA(3,4)*SPA(4,5))-(pow(SPB(4,3),2)*SPA(2,5))/
(complex<T>(6,0)*S(2,3)*SPB(5,2))+
 (complex<T>(1,0)*pow(S(3,5)*SPB(4,3)+SPA(2,5)*SPB(4,2)*SPB(5,3),2)*
 SPA(2,5))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,2))-
 (pow(SPB(4,3),2)*SPA(2,5))/(complex<T>(6,0)*S(4,5)*SPB(5,2))+
 (complex<T>(1,0)*pow(S(2,4)*SPB(4,3)+SPA(2,5)*SPB(4,2)*SPB(5,3),2)*
 SPA(2,5))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(5,2))-
 (SPA(2,5)*SPB(4,3))/(complex<T>(3,0)*SPA(3,4)*SPB(5,2))+
 (complex<T>(2,0)*pow(SPB(4,3),3))/(SPB(3,2)*SPB(5,2)*SPB(5,4))-
 (S(2,3)*SPA(2,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPB(5,2),2)*SPA(2,3)*
 SPA(3,4))-(S(2,4)*SPA(2,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPB(5,2),2)*
 SPA(2,3)*SPA(3,4))+(complex<T>(1,0)*pow(SPA(2,5),2)*S(3,5)*SPB(3,2)*
 SPB(4,3)*SPB(5,4))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*
 SPB(5,2))+(complex<T>(1,0)*pow(SPA(2,5),2)*S(2,4)*SPB(3,2)*SPB(4,3)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2)*SPB(5,3)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2)*SPB(5,3)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(5,2))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPB(3,2))/
(complex<T>(3,0)*pow(SPB(5,2),2)*SPA(2,3)*SS(2,3,5))+
 (pow(SPB(4,3),2)*SPA(2,5))/(complex<T>(6,0)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(S(3,5)*SPB(4,3)+SPA(2,5)*SPB(4,2)*SPB(5,3),2)*
 SPA(2,5))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(2,5),2)*SPB(4,3)*SPB(5,4))/(complex<T>(6,0)*SPA(2,3)*SPB(5,2)*
 SS(2,3,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*S(3,5)*SPB(3,2)*SPB(4,3)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2)*SPB(5,3)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,5)+S(3,5),2)*SPB(5,2)*SS(2,3,5))+
 (pow(SPB(4,3),2)*SPA(2,5))/(complex<T>(6,0)*SPB(5,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(S(2,4)*SPB(4,3)+SPA(2,5)*SPB(4,2)*SPB(5,3),2)*
 SPA(2,5))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPB(5,2)*SS(2,4,5))-
 (pow(SPA(2,5),2)*SPB(3,2)*SPB(4,3))/(complex<T>(6,0)*SPA(4,5)*SPB(5,2)*
 SS(2,4,5))+(complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),
2)*SPB(5,4))/(complex<T>(3,0)*pow(SPB(5,2),2)*SPA(4,5)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*S(2,4)*SPB(3,2)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPB(5,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2)*SPB(5,3)*
 SPB(5,4))/(complex<T>(6,0)*pow(S(2,4)+S(2,5),2)*SPB(5,2)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmppm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(4,3),2))/(complex<T>(9,0)*pow(SPB(5,2),2))+
 (complex<T>(1,0)*S(3,5)*SPA(2,5)*SPB(3,2))/(complex<T>(9,0)*pow(SPB(5,2),2)*
 SPA(3,4)*SPA(4,5))+(complex<T>(1,0)*S(4,5)*SPA(2,5)*SPB(3,2))/
(complex<T>(9,0)*pow(SPB(5,2),2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,5))/(complex<T>(18,0)*S(2,3)*SPB(5,2))-
 (pow(S(3,5)*SPB(4,3)+SPA(2,5)*SPB(4,2)*SPB(5,3),2)*SPA(2,5))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,5))/(complex<T>(18,0)*S(4,5)*SPB(5,2))-
 (pow(S(2,4)*SPB(4,3)+SPA(2,5)*SPB(4,2)*SPB(5,3),2)*SPA(2,5))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(5,2))+
 (complex<T>(1,0)*SPA(2,5)*SPB(4,3))/(complex<T>(9,0)*SPA(3,4)*SPB(5,2))+
 (complex<T>(1,0)*S(2,3)*SPA(2,5)*SPB(5,4))/(complex<T>(9,0)*pow(SPB(5,2),2)*
 SPA(2,3)*SPA(3,4))+(complex<T>(1,0)*S(2,4)*SPA(2,5)*SPB(5,4))/
(complex<T>(9,0)*pow(SPB(5,2),2)*SPA(2,3)*SPA(3,4))-
 (pow(SPA(2,5),2)*S(3,5)*SPB(3,2)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,2))-
 (pow(SPA(2,5),2)*S(2,4)*SPB(3,2)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(5,2))-
 (pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*S(2,3)*SPB(5,2))-
 (pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(5,2))-
 (pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPB(3,2))/
(complex<T>(9,0)*pow(SPB(5,2),2)*SPA(2,3)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,5))/(complex<T>(18,0)*SPB(5,2)*SS(2,3,5))-
 (pow(S(3,5)*SPB(4,3)+SPA(2,5)*SPB(4,2)*SPB(5,3),2)*SPA(2,5))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*SPA(2,3)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(2,5),2)*S(3,5)*SPB(3,2)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,5)+S(3,5),2)*SPB(5,2)*SS(2,3,5))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,5))/(complex<T>(18,0)*SPB(5,2)*SS(2,4,5))-
 (pow(S(2,4)*SPB(4,3)+SPA(2,5)*SPB(4,2)*SPB(5,3),2)*SPA(2,5))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPB(5,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(3,2)*SPB(4,3))/
(complex<T>(18,0)*SPA(4,5)*SPB(5,2)*SS(2,4,5))-
 (pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPB(5,4))/
(complex<T>(9,0)*pow(SPB(5,2),2)*SPA(4,5)*SS(2,4,5))-
 (pow(SPA(2,5),2)*S(2,4)*SPB(3,2)*SPB(4,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPB(5,2)*SS(2,4,5))-
 (pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(18,0)*pow(S(2,4)+S(2,5),2)*SPB(5,2)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmpmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmpmp G");
#endif
 
 return( (complex<T>(0,2)*pow(SPA(2,4),4))/(SPA(2,3)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))-(pow(SPA(2,4),2)*S(2,3)*S(3,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))-
 (pow(SPA(2,4),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,2))-
 (pow(SPA(2,4),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPB(4,2))-
 (pow(SPA(2,4),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(4,2))-
 (pow(SPA(2,4),3)*pow(SPB(5,2),2)*SPB(4,3))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*SPA(3,4)*SPB(4,2))-
 (pow(SPA(2,4),2)*S(4,5)*SPB(5,2))/(complex<T>(2,0)*pow(SPB(4,2),2)*
 SPA(2,3)*SPA(3,4)*SPA(4,5))+(complex<T>(0,-2)*pow(SPB(5,3),4))/
(SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,3)+S(2,4))*SPA(3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*S(2,3)*(S(2,4)+S(3,4)))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,4)+S(2,5))*S(4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*S(2,5)*(S(2,4)+S(4,5)))-
 (pow(SPA(2,4),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(3,0)*pow(S(2,3)+S(2,4),2)*SPB(4,2)*SS(2,3,4))-
 (pow(SPA(2,4),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPB(4,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*S(4,5)*SPB(5,2))/(complex<T>(6,0)*SPA(2,3)*
 SPA(3,4)*SPA(4,5)*SPB(4,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,3)+S(2,4))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,4)+S(3,4))*SS(2,3,4))-
 (pow(SPA(2,4),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(3,0)*pow(S(2,4)+S(4,5),2)*SPB(4,2)*SS(2,4,5))-
 (pow(SPA(2,4),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(3,0)*pow(S(2,4)+S(2,5),2)*SPB(4,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*S(2,3)*S(3,4))/(complex<T>(6,0)*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5)*SPB(4,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,4)+S(2,5))*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(4,2),2)*(S(2,4)+S(4,5))*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmpmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmpmp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,4),2)*S(2,3)*S(3,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*S(2,3)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*S(2,5)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*S(4,5)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(5,2),2)*SPB(4,3))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*SPA(3,4)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*S(4,5)*SPB(5,2))/(complex<T>(6,0)*pow(SPB(4,2),2)*
 SPA(2,3)*SPA(3,4)*SPA(4,5))-(pow(SPA(2,4),2)*SPB(3,2)*
 SPB(5,2)*SPB(5,4))/(complex<T>(6,0)*pow(SPB(4,2),2)*(S(2,3)+S(2,4))*
 SPA(3,4))-(pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*
 SPB(5,4))/(complex<T>(6,0)*pow(SPB(4,2),2)*S(2,3)*(S(2,4)+S(3,4)))-
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*(S(2,4)+S(2,5))*S(4,5))-
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*S(2,5)*(S(2,4)+S(4,5)))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(9,0)*pow(S(2,3)+S(2,4),2)*SPB(4,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(9,0)*pow(S(2,4)+S(3,4),2)*SPB(4,2)*SS(2,3,4))-
 (pow(SPA(2,4),3)*S(4,5)*SPB(5,2))/(complex<T>(18,0)*SPA(2,3)*SPA(3,4)*
 SPA(4,5)*SPB(4,2)*SS(2,3,4))-
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*(S(2,3)+S(2,4))*SS(2,3,4))-
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*(S(2,4)+S(3,4))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(9,0)*pow(S(2,4)+S(4,5),2)*SPB(4,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(9,0)*pow(S(2,4)+S(2,5),2)*SPB(4,2)*SS(2,4,5))-
 (pow(SPA(2,4),3)*S(2,3)*S(3,4))/(complex<T>(18,0)*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5)*SPB(4,2)*SS(2,4,5))-
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*(S(2,4)+S(2,5))*SS(2,4,5))-
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(4,2),2)*(S(2,4)+S(4,5))*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpmpm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpmpm G");
#endif
 
 return( (complex<T>(0,2)*pow(SPA(3,5),4))/(SPA(2,3)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))-(pow(SPA(3,5),2)*S(3,4)*S(4,5))/
(complex<T>(2,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))-
 (pow(SPA(3,5),2)*S(2,3)*SPB(5,2))/(complex<T>(2,0)*pow(SPB(5,3),2)*
 SPA(2,3)*SPA(3,4)*SPA(4,5))+(complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*
 SPB(4,3)*SPB(5,2))/(complex<T>(2,0)*pow(SPB(5,3),2)*(S(3,4)+S(3,5))*
 SPA(4,5))-(pow(SPA(3,5),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*SPB(5,3))-
 (pow(SPA(3,5),3)*pow(SPB(5,4),2)*SPB(3,2))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPA(2,3)*SPB(5,3))-
 (pow(SPA(3,5),3)*pow(SPB(4,3),2)*SPB(5,2))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*SPA(2,5)*SPB(5,3))+
 (complex<T>(0,-2)*pow(SPB(4,2),4))/(SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(2,3)+S(3,5))*SPA(2,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(2,5)+S(3,5))*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(5,3),2)*S(3,4)*(S(3,5)+S(4,5)))-
 (pow(SPA(3,5),3)*pow(SPB(3,2),2)*SPB(5,4))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*SPA(4,5)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(3,5),4)*S(3,4)*S(4,5))/(complex<T>(6,0)*S(3,5)*SPA(2,3)*
 SPA(2,5)*SPA(3,4)*SPA(4,5)*SS(2,3,5))-
 (pow(SPA(3,5),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(3,0)*pow(S(2,3)+S(3,5),2)*SPB(5,3)*SS(2,3,5))-
 (pow(SPA(3,5),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(2,3)+S(3,5))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(2,5)+S(3,5))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),4)*S(2,3)*SPB(5,2))/(complex<T>(6,0)*S(3,5)*SPA(2,3)*
 SPA(3,4)*SPA(4,5)*SS(3,4,5))-
 (pow(SPA(3,5),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(3,0)*pow(S(3,5)+S(4,5),2)*SPB(5,3)*SS(3,4,5))-
 (pow(SPA(3,5),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(3,0)*pow(S(3,4)+S(3,5),2)*SPB(5,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(3,4)+S(3,5))*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(5,3),2)*(S(3,5)+S(4,5))*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpmpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpmpm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(3,5),2)*S(3,4)*S(4,5))/
(complex<T>(6,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*S(2,3)*SPB(5,2))/(complex<T>(6,0)*pow(SPB(5,3),2)*
 SPA(2,3)*SPA(3,4)*SPA(4,5))-(pow(SPA(3,5),2)*SPB(3,2)*
 SPB(4,3)*SPB(5,2))/(complex<T>(6,0)*pow(SPB(5,3),2)*(S(3,4)+S(3,5))*
 SPA(4,5))+(complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPB(4,3),2)*
 pow(SPB(5,2),2))/(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*S(3,4)*
 SPB(5,3))+(complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPB(5,4),2)*SPB(3,2))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*SPA(2,3)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPB(4,3),2)*SPB(5,2))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*SPA(2,5)*SPB(5,3))-
 (pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(5,3),2)*(S(2,3)+S(3,5))*SPA(2,5))-
 (pow(SPA(3,5),2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(5,3),2)*(S(2,5)+S(3,5))*SPA(2,3))-
 (pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(5,3),2)*S(3,4)*(S(3,5)+S(4,5)))+
 (complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPB(3,2),2)*SPB(5,4))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*SPA(4,5)*SPB(5,3))-
 (pow(SPA(3,5),4)*S(3,4)*S(4,5))/(complex<T>(18,0)*S(3,5)*SPA(2,3)*
 SPA(2,5)*SPA(3,4)*SPA(4,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(9,0)*pow(S(2,3)+S(3,5),2)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(9,0)*pow(S(2,5)+S(3,5),2)*SPB(5,3)*SS(2,3,5))-
 (pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(5,3),2)*(S(2,3)+S(3,5))*SS(2,3,5))-
 (pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(5,3),2)*(S(2,5)+S(3,5))*SS(2,3,5))-
 (pow(SPA(3,5),4)*S(2,3)*SPB(5,2))/(complex<T>(18,0)*S(3,5)*SPA(2,3)*
 SPA(3,4)*SPA(4,5)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(9,0)*pow(S(3,5)+S(4,5),2)*SPB(5,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(9,0)*pow(S(3,4)+S(3,5),2)*SPB(5,3)*SS(3,4,5))-
 (pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(5,3),2)*(S(3,4)+S(3,5))*SS(3,4,5))-
 (pow(SPA(3,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPB(5,3),2)*(S(3,5)+S(4,5))*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpmmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpmmm G");
#endif
 
 return(-((pow(SPB(4,2),2)*SPA(4,5)*SPB(5,2))/(complex<T>(3,0)*pow(SPB(5,4),2)*
SPB(3,2)*SPB(4,3)))-(pow(SPB(4,2),2)*SPA(3,4)*SPB(3,2))/
(complex<T>(3,0)*pow(SPB(4,3),2)*SPB(5,2)*SPB(5,4))+
 (complex<T>(-2,0)*pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),3))/
(SPB(3,2)*SPB(4,3)*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*
 SS(2,3,4))+(complex<T>(2,0)*pow(mH2,2)*pow(SPA(3,5),4))/
(SPA(2,3)*SPA(2,5)*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*
 (SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4))*SS(2,3,5))+
 (complex<T>(2,0)*pow(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2),3))/
(SPB(5,2)*SPB(5,4)*(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4))*
 SS(2,4,5))-(complex<T>(1,0)*pow(SPB(4,2),2)*SPA(3,4)*SPA(4,5))/
(complex<T>(3,0)*SPB(4,3)*SPB(5,4)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(3,4)*SPA(3,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(3,0)*SPB(4,3)*SPB(5,4)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(3,2)*SPB(5,2))/
(complex<T>(3,0)*SPB(4,3)*SPB(5,4)*SS(3,4,5))-
 (complex<T>(1,0)*SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(3,0)*SPB(4,3)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdpmmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdpmmm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(4,5)*SPB(5,2))/
(complex<T>(9,0)*pow(SPB(5,4),2)*SPB(3,2)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(3,4)*SPB(3,2))/
(complex<T>(9,0)*pow(SPB(4,3),2)*SPB(5,2)*SPB(5,4))+
 (pow(SPB(4,2),2)*SPA(3,4)*SPA(4,5))/(complex<T>(9,0)*SPB(4,3)*SPB(5,4)*
 SS(3,4,5))+(SPA(3,4)*SPA(3,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(9,0)*SPB(4,3)*SPB(5,4)*SS(3,4,5))+
 (pow(SPA(3,5),2)*SPB(3,2)*SPB(5,2))/(complex<T>(9,0)*SPB(4,3)*SPB(5,4)*
 SS(3,4,5))+(SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,2))/
(complex<T>(9,0)*SPB(4,3)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmpmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmpmm G");
#endif
 
 return(-((pow(SPB(5,3),2)*SPA(4,5)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*
SPB(3,2)*SPB(5,2)))-(pow(SPB(5,3),2)*SPA(2,5)*SPB(3,2))/
(complex<T>(3,0)*pow(SPB(5,2),2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(-2,0)*pow(mH2,2)*pow(SPA(2,4),4))/(SPA(2,3)*SPA(3,4)*
 (SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*(-(SPA(2,3)*SPB(5,3))-
SPA(2,4)*SPB(5,4))*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),3))/
(SPB(3,2)*SPB(5,2)*(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*
 SS(2,3,5))-(pow(SPB(5,3),2)*SPA(2,5)*SPA(4,5))/
(complex<T>(3,0)*SPB(5,2)*SPB(5,4)*SS(2,4,5))-
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3))/(complex<T>(3,0)*SPB(5,2)*SPB(5,4)*
 SS(2,4,5))-(SPA(2,4)*SPA(2,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(3,0)*SPB(5,2)*SPB(5,4)*SS(2,4,5))-
 (SPA(2,4)*SPA(4,5)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*SPB(5,2)*SPB(5,4)*
 SS(2,4,5))+(complex<T>(-2,0)*pow(SPA(2,4)*SPB(4,3)+
 SPA(2,5)*SPB(5,3),3))/(SPB(4,3)*SPB(5,4)*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4))*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmpmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmpmm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(4,5)*SPB(4,3))/
(complex<T>(9,0)*pow(SPB(5,4),2)*SPB(3,2)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,5)*SPB(3,2))/
(complex<T>(9,0)*pow(SPB(5,2),2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,5)*SPA(4,5))/
(complex<T>(9,0)*SPB(5,2)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(4,3))/
(complex<T>(9,0)*SPB(5,2)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(9,0)*SPB(5,2)*SPB(5,4)*SS(2,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(4,5)*SPB(4,3)*SPB(5,3))/
(complex<T>(9,0)*SPB(5,2)*SPB(5,4)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmmpm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmmpm G");
#endif
 
 return(-((pow(SPB(4,2),2)*SPA(2,3)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(3,2),2)*
SPB(5,2)*SPB(5,4)))-(pow(SPB(4,2),2)*SPA(2,5)*SPB(5,4))/
(complex<T>(3,0)*pow(SPB(5,2),2)*SPB(3,2)*SPB(4,3))+
 (complex<T>(-2,0)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3))/
(SPB(3,2)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(4,3)*
 SS(2,3,4))-(pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5))/
(complex<T>(3,0)*SPB(3,2)*SPB(5,2)*SS(2,3,5))-
 (SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(4,3))/(complex<T>(3,0)*SPB(3,2)*SPB(5,2)*
 SS(2,3,5))-(SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(3,0)*SPB(3,2)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(3,5),2)*SPB(4,3)*SPB(5,4))/(complex<T>(3,0)*SPB(3,2)*SPB(5,2)*
 SS(2,3,5))+(complex<T>(2,0)*pow(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4),
3))/(SPB(5,2)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*SPB(5,4)*
 SS(2,4,5))+(complex<T>(-2,0)*pow(mH2,2)*pow(SPA(3,5),4))/
(SPA(3,4)*SPA(4,5)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*
 (SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmmpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmmpm nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPB(4,3))/
(complex<T>(9,0)*pow(SPB(3,2),2)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,5)*SPB(5,4))/
(complex<T>(9,0)*pow(SPB(5,2),2)*SPB(3,2)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(2,5))/
(complex<T>(9,0)*SPB(3,2)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(9,0)*SPB(3,2)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,5)*SPA(3,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(9,0)*SPB(3,2)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPB(4,3)*SPB(5,4))/
(complex<T>(9,0)*SPB(3,2)*SPB(5,2)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmmmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmmmp G");
#endif
 
 return(-((pow(SPB(5,3),2)*SPA(2,3)*SPB(5,2))/(complex<T>(3,0)*pow(SPB(3,2),2)*
SPB(4,3)*SPB(5,4)))-(pow(SPB(5,3),2)*SPA(3,4)*SPB(5,4))/
(complex<T>(3,0)*pow(SPB(4,3),2)*SPB(3,2)*SPB(5,2))-
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(3,4))/
(complex<T>(3,0)*SPB(3,2)*SPB(4,3)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(3,0)*SPB(3,2)*SPB(4,3)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(5,2)*SPB(5,4))/
(complex<T>(3,0)*SPB(3,2)*SPB(4,3)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(3,0)*SPB(3,2)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),3))/
(SPB(3,2)*SPB(5,2)*(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3))*
 SS(2,3,5))+(complex<T>(2,0)*pow(mH2,2)*pow(SPA(2,4),4))/
(SPA(2,5)*SPA(4,5)*(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3))*
 (SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3))*SS(2,4,5))+
 (complex<T>(-2,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),3))/
(SPB(4,3)*(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3))*SPB(5,4)*
 SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmmmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmmmp nf");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPB(5,2))/
(complex<T>(9,0)*pow(SPB(3,2),2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(3,4)*SPB(5,4))/
(complex<T>(9,0)*pow(SPB(4,3),2)*SPB(3,2)*SPB(5,2))+
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(3,4))/(complex<T>(9,0)*SPB(3,2)*SPB(4,3)*
 SS(2,3,4))+(SPA(2,3)*SPA(2,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(9,0)*SPB(3,2)*SPB(4,3)*SS(2,3,4))+
 (pow(SPA(2,4),2)*SPB(5,2)*SPB(5,4))/(complex<T>(9,0)*SPB(3,2)*SPB(4,3)*
 SS(2,3,4))+(SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(9,0)*SPB(3,2)*SPB(4,3)*SS(2,3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmmmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmmmm G");
#endif
 
 return( (complex<T>(2,0)*pow(mH2,2))/(SPA(2,3)*SPA(2,5)*SPA(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R4g1ph_phdmmmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R4g1ph :  phdmmmm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 // *************** table of switch values ************* 
 
#define _R_phpppp_G R4g1ph_1021_G
#define _R_phpppp_nf R4g1ph_1021_nf
#define _R_phmppp_G R4g1ph_1009_G
#define _R_phmppp_nf R4g1ph_1009_nf
#define _R_phpmpp_G R4g1ph_973_G
#define _R_phpmpp_nf R4g1ph_973_nf
#define _R_phppmp_G R4g1ph_829_G
#define _R_phppmp_nf R4g1ph_829_nf
#define _R_phpppm_G R4g1ph_253_G
#define _R_phpppm_nf R4g1ph_253_nf
#define _R_phmmpp_G R4g1ph_961_G
#define _R_phmmpp_nf R4g1ph_961_nf
#define _R_phpmmp_G R4g1ph_781_G
#define _R_phpmmp_nf R4g1ph_781_nf
#define _R_phppmm_G R4g1ph_61_G
#define _R_phppmm_nf R4g1ph_61_nf
#define _R_phmppm_G R4g1ph_241_G
#define _R_phmppm_nf R4g1ph_241_nf
#define _R_phmpmp_G R4g1ph_817_G
#define _R_phmpmp_nf R4g1ph_817_nf
#define _R_phpmpm_G R4g1ph_205_G
#define _R_phpmpm_nf R4g1ph_205_nf
#define _R_phpmmm_G R4g1ph_13_G
#define _R_phpmmm_nf R4g1ph_13_nf
#define _R_phmpmm_G R4g1ph_49_G
#define _R_phmpmm_nf R4g1ph_49_nf
#define _R_phmmpm_G R4g1ph_193_G
#define _R_phmmpm_nf R4g1ph_193_nf
#define _R_phmmmp_G R4g1ph_769_G
#define _R_phmmmp_nf R4g1ph_769_nf
#define _R_phmmmm_G R4g1ph_1_G
#define _R_phmmmm_nf R4g1ph_1_nf
#define _R_phdpppp_G R4g1ph_1022_G
#define _R_phdpppp_nf R4g1ph_1022_nf
#define _R_phdmppp_G R4g1ph_1010_G
#define _R_phdmppp_nf R4g1ph_1010_nf
#define _R_phdpmpp_G R4g1ph_974_G
#define _R_phdpmpp_nf R4g1ph_974_nf
#define _R_phdppmp_G R4g1ph_830_G
#define _R_phdppmp_nf R4g1ph_830_nf
#define _R_phdpppm_G R4g1ph_254_G
#define _R_phdpppm_nf R4g1ph_254_nf
#define _R_phdmmpp_G R4g1ph_962_G
#define _R_phdmmpp_nf R4g1ph_962_nf
#define _R_phdpmmp_G R4g1ph_782_G
#define _R_phdpmmp_nf R4g1ph_782_nf
#define _R_phdppmm_G R4g1ph_62_G
#define _R_phdppmm_nf R4g1ph_62_nf
#define _R_phdmppm_G R4g1ph_242_G
#define _R_phdmppm_nf R4g1ph_242_nf
#define _R_phdmpmp_G R4g1ph_818_G
#define _R_phdmpmp_nf R4g1ph_818_nf
#define _R_phdpmpm_G R4g1ph_206_G
#define _R_phdpmpm_nf R4g1ph_206_nf
#define _R_phdpmmm_G R4g1ph_14_G
#define _R_phdpmmm_nf R4g1ph_14_nf
#define _R_phdmpmm_G R4g1ph_50_G
#define _R_phdmpmm_nf R4g1ph_50_nf
#define _R_phdmmpm_G R4g1ph_194_G
#define _R_phdmmpm_nf R4g1ph_194_nf
#define _R_phdmmmp_G R4g1ph_770_G
#define _R_phdmmmp_nf R4g1ph_770_nf
#define _R_phdmmmm_G R4g1ph_2_G
#define _R_phdmmmm_nf R4g1ph_2_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_phpppp_G case 1021 : \
          return &R4g1ph_1021_G
#define _CASE_phpppp_nf case 1021 : \
          return &R4g1ph_1021_nf
#define _CASE_phmppp_G case 1009 : \
          return &R4g1ph_1009_G
#define _CASE_phmppp_nf case 1009 : \
          return &R4g1ph_1009_nf
#define _CASE_phpmpp_G case 973 : \
          return &R4g1ph_973_G
#define _CASE_phpmpp_nf case 973 : \
          return &R4g1ph_973_nf
#define _CASE_phppmp_G case 829 : \
          return &R4g1ph_829_G
#define _CASE_phppmp_nf case 829 : \
          return &R4g1ph_829_nf
#define _CASE_phpppm_G case 253 : \
          return &R4g1ph_253_G
#define _CASE_phpppm_nf case 253 : \
          return &R4g1ph_253_nf
#define _CASE_phmmpp_G case 961 : \
          return &R4g1ph_961_G
#define _CASE_phmmpp_nf case 961 : \
          return &R4g1ph_961_nf
#define _CASE_phpmmp_G case 781 : \
          return &R4g1ph_781_G
#define _CASE_phpmmp_nf case 781 : \
          return &R4g1ph_781_nf
#define _CASE_phppmm_G case 61 : \
          return &R4g1ph_61_G
#define _CASE_phppmm_nf case 61 : \
          return &R4g1ph_61_nf
#define _CASE_phmppm_G case 241 : \
          return &R4g1ph_241_G
#define _CASE_phmppm_nf case 241 : \
          return &R4g1ph_241_nf
#define _CASE_phmpmp_G case 817 : \
          return &R4g1ph_817_G
#define _CASE_phmpmp_nf case 817 : \
          return &R4g1ph_817_nf
#define _CASE_phpmpm_G case 205 : \
          return &R4g1ph_205_G
#define _CASE_phpmpm_nf case 205 : \
          return &R4g1ph_205_nf
#define _CASE_phpmmm_G case 13 : \
          return &R4g1ph_13_G
#define _CASE_phpmmm_nf case 13 : \
          return &R4g1ph_13_nf
#define _CASE_phmpmm_G case 49 : \
          return &R4g1ph_49_G
#define _CASE_phmpmm_nf case 49 : \
          return &R4g1ph_49_nf
#define _CASE_phmmpm_G case 193 : \
          return &R4g1ph_193_G
#define _CASE_phmmpm_nf case 193 : \
          return &R4g1ph_193_nf
#define _CASE_phmmmp_G case 769 : \
          return &R4g1ph_769_G
#define _CASE_phmmmp_nf case 769 : \
          return &R4g1ph_769_nf
#define _CASE_phmmmm_G case 1 : \
          return &R4g1ph_1_G
#define _CASE_phmmmm_nf case 1 : \
          return &R4g1ph_1_nf
#define _CASE_phdpppp_G case 1022 : \
          return &R4g1ph_1022_G
#define _CASE_phdpppp_nf case 1022 : \
          return &R4g1ph_1022_nf
#define _CASE_phdmppp_G case 1010 : \
          return &R4g1ph_1010_G
#define _CASE_phdmppp_nf case 1010 : \
          return &R4g1ph_1010_nf
#define _CASE_phdpmpp_G case 974 : \
          return &R4g1ph_974_G
#define _CASE_phdpmpp_nf case 974 : \
          return &R4g1ph_974_nf
#define _CASE_phdppmp_G case 830 : \
          return &R4g1ph_830_G
#define _CASE_phdppmp_nf case 830 : \
          return &R4g1ph_830_nf
#define _CASE_phdpppm_G case 254 : \
          return &R4g1ph_254_G
#define _CASE_phdpppm_nf case 254 : \
          return &R4g1ph_254_nf
#define _CASE_phdmmpp_G case 962 : \
          return &R4g1ph_962_G
#define _CASE_phdmmpp_nf case 962 : \
          return &R4g1ph_962_nf
#define _CASE_phdpmmp_G case 782 : \
          return &R4g1ph_782_G
#define _CASE_phdpmmp_nf case 782 : \
          return &R4g1ph_782_nf
#define _CASE_phdppmm_G case 62 : \
          return &R4g1ph_62_G
#define _CASE_phdppmm_nf case 62 : \
          return &R4g1ph_62_nf
#define _CASE_phdmppm_G case 242 : \
          return &R4g1ph_242_G
#define _CASE_phdmppm_nf case 242 : \
          return &R4g1ph_242_nf
#define _CASE_phdmpmp_G case 818 : \
          return &R4g1ph_818_G
#define _CASE_phdmpmp_nf case 818 : \
          return &R4g1ph_818_nf
#define _CASE_phdpmpm_G case 206 : \
          return &R4g1ph_206_G
#define _CASE_phdpmpm_nf case 206 : \
          return &R4g1ph_206_nf
#define _CASE_phdpmmm_G case 14 : \
          return &R4g1ph_14_G
#define _CASE_phdpmmm_nf case 14 : \
          return &R4g1ph_14_nf
#define _CASE_phdmpmm_G case 50 : \
          return &R4g1ph_50_G
#define _CASE_phdmpmm_nf case 50 : \
          return &R4g1ph_50_nf
#define _CASE_phdmmpm_G case 194 : \
          return &R4g1ph_194_G
#define _CASE_phdmmpm_nf case 194 : \
          return &R4g1ph_194_nf
#define _CASE_phdmmmp_G case 770 : \
          return &R4g1ph_770_G
#define _CASE_phdmmmp_nf case 770 : \
          return &R4g1ph_770_nf
#define _CASE_phdmmmm_G case 2 : \
          return &R4g1ph_2_G
#define _CASE_phdmmmm_nf case 2 : \
          return &R4g1ph_2_nf
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_phpppp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpppp_G(ep,mpc);}
 
template <class T> complex<T> _R_phpppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmppp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmppp_G(ep,mpc);}
 
template <class T> complex<T> _R_phmppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phpmpp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpmpp_G(ep,mpc);}
 
template <class T> complex<T> _R_phpmpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpmpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phppmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phppmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phppmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phppmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phpppm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpppm_G(ep,mpc);}
 
template <class T> complex<T> _R_phpppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmmpp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmmpp_G(ep,mpc);}
 
template <class T> complex<T> _R_phmmpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmmpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phpmmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpmmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phpmmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpmmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phppmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phppmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phppmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phppmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmppm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmppm_G(ep,mpc);}
 
template <class T> complex<T> _R_phmppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmpmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmpmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phmpmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmpmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phpmpm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpmpm_G(ep,mpc);}
 
template <class T> complex<T> _R_phpmpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpmpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phpmmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpmmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phpmmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phpmmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmpmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmpmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phmpmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmpmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmmpm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmmpm_G(ep,mpc);}
 
template <class T> complex<T> _R_phmmpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmmpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmmmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmmmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phmmmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmmmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmmmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmmmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phmmmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phmmmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdpppp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpppp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdpppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmppp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmppp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdpmpp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpmpp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdpmpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpmpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdppmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdppmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdppmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdppmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdpppm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpppm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdpppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmmpp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmmpp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmmpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmmpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdpmmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpmmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdpmmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpmmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdppmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdppmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdppmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdppmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmppm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmppm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmpmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmpmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmpmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmpmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdpmpm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpmpm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdpmpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpmpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdpmmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpmmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdpmmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdpmmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmpmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmpmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmpmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmpmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmmpm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmmpm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmmpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmmpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmmmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmmmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmmmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmmmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmmmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmmmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmmmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R4g1ph_phdmmmm_nf(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R4g1ph_G_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phpppp_G;
       _CASE_phmppp_G;
       _CASE_phpmpp_G;
       _CASE_phppmp_G;
       _CASE_phpppm_G;
       _CASE_phmmpp_G;
       _CASE_phpmmp_G;
       _CASE_phppmm_G;
       _CASE_phmppm_G;
       _CASE_phmpmp_G;
       _CASE_phpmpm_G;
       _CASE_phpmmm_G;
       _CASE_phmpmm_G;
       _CASE_phmmpm_G;
       _CASE_phmmmp_G;
       _CASE_phmmmm_G;
       _CASE_phdpppp_G;
       _CASE_phdmppp_G;
       _CASE_phdpmpp_G;
       _CASE_phdppmp_G;
       _CASE_phdpppm_G;
       _CASE_phdmmpp_G;
       _CASE_phdpmmp_G;
       _CASE_phdppmm_G;
       _CASE_phdmppm_G;
       _CASE_phdmpmp_G;
       _CASE_phdpmpm_G;
       _CASE_phdpmmm_G;
       _CASE_phdmpmm_G;
       _CASE_phdmmpm_G;
       _CASE_phdmmmp_G;
       _CASE_phdmmmm_G;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R4g1ph_nf_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phpppp_nf;
       _CASE_phmppp_nf;
       _CASE_phpmpp_nf;
       _CASE_phppmp_nf;
       _CASE_phpppm_nf;
       _CASE_phmmpp_nf;
       _CASE_phpmmp_nf;
       _CASE_phppmm_nf;
       _CASE_phmppm_nf;
       _CASE_phmpmp_nf;
       _CASE_phpmpm_nf;
       _CASE_phpmmm_nf;
       _CASE_phmpmm_nf;
       _CASE_phmmpm_nf;
       _CASE_phmmmp_nf;
       _CASE_phmmmm_nf;
       _CASE_phdpppp_nf;
       _CASE_phdmppp_nf;
       _CASE_phdpmpp_nf;
       _CASE_phdppmp_nf;
       _CASE_phdpppm_nf;
       _CASE_phdmmpp_nf;
       _CASE_phdpmmp_nf;
       _CASE_phdppmm_nf;
       _CASE_phdmppm_nf;
       _CASE_phdmpmp_nf;
       _CASE_phdpmpm_nf;
       _CASE_phdpmmm_nf;
       _CASE_phdmpmm_nf;
       _CASE_phdmmpm_nf;
       _CASE_phdmmmp_nf;
       _CASE_phdmmmm_nf;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R4g1ph_G_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R4g1ph_G_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R4g1ph_G_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R4g1ph_nf_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R4g1ph_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R4g1ph_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> ( *R4g1ph_G_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R4g1ph_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);

#endif


}
