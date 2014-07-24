/*
* R_3g1ph_eval.cpp
*
* Created on 12/9, 2010
*      Author: Zvi's script
*/
 
#include "R_3g1ph_eval.h"
#include "eval_param.h"
 
using namespace std;
 
namespace BH  {
 
 
#define _VERBOSE 1
#define mH2 ep.s(1,2,3)
 
#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}



 
 
template <class T> complex<T> R3g1ph_phppp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phppp G");
#endif
 
 return( (complex<T>(2,0)*pow(mH2,2))/(SPA(2,3)*SPA(2,4)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phppp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R3g1ph_phmpp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phmpp G");
#endif
 
 return( (complex<T>(2,0)*pow(SPB(4,3),3))/(SPB(3,2)*SPB(4,2))-
 (SPA(2,3)*SPA(2,4)*SPB(4,3))/(complex<T>(3,0)*pow(SPA(3,4),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phmpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phmpp nf");
#endif
 
 return( (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPB(4,3))/(complex<T>(9,0)*pow(SPA(3,4),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phpmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phpmp G");
#endif
 
 return(-((SPA(2,3)*SPA(3,4)*SPB(4,2))/(complex<T>(3,0)*pow(SPA(2,4),2)))+
 (complex<T>(2,0)*pow(SPB(4,2),3))/(SPB(3,2)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phpmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phpmp nf");
#endif
 
 return( (complex<T>(1,0)*SPA(2,3)*SPA(3,4)*SPB(4,2))/(complex<T>(9,0)*pow(SPA(2,4),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phppm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phppm G");
#endif
 
 return(-((SPA(2,4)*SPA(3,4)*SPB(3,2))/(complex<T>(3,0)*pow(SPA(2,3),2)))+
 (complex<T>(2,0)*pow(SPB(3,2),3))/(SPB(4,2)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phppm nf");
#endif
 
 return( (complex<T>(1,0)*SPA(2,4)*SPA(3,4)*SPB(3,2))/(complex<T>(9,0)*pow(SPA(2,3),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phmmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phmmp G");
#endif
 
 return( (complex<T>(-2,0)*pow(SPA(2,3),3))/(SPA(2,4)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phmmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phmmp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R3g1ph_phmpm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phmpm G");
#endif
 
 return( (complex<T>(-2,0)*pow(SPA(2,4),3))/(SPA(2,3)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phmpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phmpm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R3g1ph_phpmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phpmm G");
#endif
 
 return( (complex<T>(-2,0)*pow(SPA(3,4),3))/(SPA(2,3)*SPA(2,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phpmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phpmm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R3g1ph_phmmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phmmm G");
#endif
 
 return(-((SPA(2,4)*SPA(3,4))/(complex<T>(3,0)*SPB(3,2)))-
 (SPA(2,3)*SPA(3,4))/(complex<T>(3,0)*SPB(4,2))-(SPA(2,3)*SPA(2,4))/
(complex<T>(3,0)*SPB(4,3))+(complex<T>(-2,0)*pow(mH2,2))/(SPB(3,2)*SPB(4,2)*
 SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phmmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phmmm nf");
#endif
 
 return( (complex<T>(1,0)*SPA(2,4)*SPA(3,4))/(complex<T>(9,0)*SPB(3,2))+
 (complex<T>(1,0)*SPA(2,3)*SPA(3,4))/(complex<T>(9,0)*SPB(4,2))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4))/(complex<T>(9,0)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdppp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdppp G");
#endif
 
 return( (complex<T>(-2,0)*pow(mH2,2))/(SPA(2,3)*SPA(2,4)*SPA(3,4))-
 (S(2,3)*S(2,4))/(complex<T>(3,0)*SPA(2,3)*SPA(2,4)*SPA(3,4))-
 (S(2,3)*S(3,4))/(complex<T>(3,0)*SPA(2,3)*SPA(2,4)*SPA(3,4))-
 (S(2,4)*S(3,4))/(complex<T>(3,0)*SPA(2,3)*SPA(2,4)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdppp nf");
#endif
 
 return( (complex<T>(1,0)*S(2,3)*S(2,4))/(complex<T>(9,0)*SPA(2,3)*SPA(2,4)*SPA(3,4))+
 (complex<T>(1,0)*S(2,3)*S(3,4))/(complex<T>(9,0)*SPA(2,3)*SPA(2,4)*SPA(3,4))+
 (complex<T>(1,0)*S(2,4)*S(3,4))/(complex<T>(9,0)*SPA(2,3)*SPA(2,4)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdmpp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, p}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdmpp G");
#endif
 
 return( (complex<T>(-2,0)*pow(SPB(4,3),3))/(SPB(3,2)*SPB(4,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdmpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdmpp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R3g1ph_phdpmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdpmp G");
#endif
 
 return( (complex<T>(-2,0)*pow(SPB(4,2),3))/(SPB(3,2)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdpmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdpmp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R3g1ph_phdppm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdppm G");
#endif
 
 return( (complex<T>(-2,0)*pow(SPB(3,2),3))/(SPB(4,2)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdppm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R3g1ph_phdmmp_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, p}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdmmp G");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(2,3),3))/(SPA(2,4)*SPA(3,4))-
 (SPA(2,3)*SPB(4,2)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(3,2),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdmmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdmmp nf");
#endif
 
 return( (complex<T>(1,0)*SPA(2,3)*SPB(4,2)*SPB(4,3))/(complex<T>(9,0)*pow(SPB(3,2),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdmpm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, m}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdmpm G");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(2,4),3))/(SPA(2,3)*SPA(3,4))-
 (SPA(2,4)*SPB(3,2)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(4,2),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdmpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdmpm nf");
#endif
 
 return( (complex<T>(1,0)*SPA(2,4)*SPB(3,2)*SPB(4,3))/(complex<T>(9,0)*pow(SPB(4,2),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdpmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdpmm G");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(3,4),3))/(SPA(2,3)*SPA(2,4))-
 (SPA(3,4)*SPB(3,2)*SPB(4,2))/(complex<T>(3,0)*pow(SPB(4,3),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdpmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdpmm nf");
#endif
 
 return( (complex<T>(1,0)*SPA(3,4)*SPB(3,2)*SPB(4,2))/(complex<T>(9,0)*pow(SPB(4,3),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdmmm_G
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, m}, G}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdmmm G");
#endif
 
 return( (complex<T>(2,0)*pow(mH2,2))/(SPB(3,2)*SPB(4,2)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R3g1ph_phdmmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R3g1ph :  phdmmm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 // *************** table of switch values ************* 
 
#define _R_phppp_G R3g1ph_253_G
#define _R_phppp_nf R3g1ph_253_nf
#define _R_phmpp_G R3g1ph_241_G
#define _R_phmpp_nf R3g1ph_241_nf
#define _R_phpmp_G R3g1ph_205_G
#define _R_phpmp_nf R3g1ph_205_nf
#define _R_phppm_G R3g1ph_61_G
#define _R_phppm_nf R3g1ph_61_nf
#define _R_phmmp_G R3g1ph_193_G
#define _R_phmmp_nf R3g1ph_193_nf
#define _R_phmpm_G R3g1ph_49_G
#define _R_phmpm_nf R3g1ph_49_nf
#define _R_phpmm_G R3g1ph_13_G
#define _R_phpmm_nf R3g1ph_13_nf
#define _R_phmmm_G R3g1ph_1_G
#define _R_phmmm_nf R3g1ph_1_nf
#define _R_phdppp_G R3g1ph_254_G
#define _R_phdppp_nf R3g1ph_254_nf
#define _R_phdmpp_G R3g1ph_242_G
#define _R_phdmpp_nf R3g1ph_242_nf
#define _R_phdpmp_G R3g1ph_206_G
#define _R_phdpmp_nf R3g1ph_206_nf
#define _R_phdppm_G R3g1ph_62_G
#define _R_phdppm_nf R3g1ph_62_nf
#define _R_phdmmp_G R3g1ph_194_G
#define _R_phdmmp_nf R3g1ph_194_nf
#define _R_phdmpm_G R3g1ph_50_G
#define _R_phdmpm_nf R3g1ph_50_nf
#define _R_phdpmm_G R3g1ph_14_G
#define _R_phdpmm_nf R3g1ph_14_nf
#define _R_phdmmm_G R3g1ph_2_G
#define _R_phdmmm_nf R3g1ph_2_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_phppp_G case 253 : \
          return &R3g1ph_253_G
#define _CASE_phppp_nf case 253 : \
          return &R3g1ph_253_nf
#define _CASE_phmpp_G case 241 : \
          return &R3g1ph_241_G
#define _CASE_phmpp_nf case 241 : \
          return &R3g1ph_241_nf
#define _CASE_phpmp_G case 205 : \
          return &R3g1ph_205_G
#define _CASE_phpmp_nf case 205 : \
          return &R3g1ph_205_nf
#define _CASE_phppm_G case 61 : \
          return &R3g1ph_61_G
#define _CASE_phppm_nf case 61 : \
          return &R3g1ph_61_nf
#define _CASE_phmmp_G case 193 : \
          return &R3g1ph_193_G
#define _CASE_phmmp_nf case 193 : \
          return &R3g1ph_193_nf
#define _CASE_phmpm_G case 49 : \
          return &R3g1ph_49_G
#define _CASE_phmpm_nf case 49 : \
          return &R3g1ph_49_nf
#define _CASE_phpmm_G case 13 : \
          return &R3g1ph_13_G
#define _CASE_phpmm_nf case 13 : \
          return &R3g1ph_13_nf
#define _CASE_phmmm_G case 1 : \
          return &R3g1ph_1_G
#define _CASE_phmmm_nf case 1 : \
          return &R3g1ph_1_nf
#define _CASE_phdppp_G case 254 : \
          return &R3g1ph_254_G
#define _CASE_phdppp_nf case 254 : \
          return &R3g1ph_254_nf
#define _CASE_phdmpp_G case 242 : \
          return &R3g1ph_242_G
#define _CASE_phdmpp_nf case 242 : \
          return &R3g1ph_242_nf
#define _CASE_phdpmp_G case 206 : \
          return &R3g1ph_206_G
#define _CASE_phdpmp_nf case 206 : \
          return &R3g1ph_206_nf
#define _CASE_phdppm_G case 62 : \
          return &R3g1ph_62_G
#define _CASE_phdppm_nf case 62 : \
          return &R3g1ph_62_nf
#define _CASE_phdmmp_G case 194 : \
          return &R3g1ph_194_G
#define _CASE_phdmmp_nf case 194 : \
          return &R3g1ph_194_nf
#define _CASE_phdmpm_G case 50 : \
          return &R3g1ph_50_G
#define _CASE_phdmpm_nf case 50 : \
          return &R3g1ph_50_nf
#define _CASE_phdpmm_G case 14 : \
          return &R3g1ph_14_G
#define _CASE_phdpmm_nf case 14 : \
          return &R3g1ph_14_nf
#define _CASE_phdmmm_G case 2 : \
          return &R3g1ph_2_G
#define _CASE_phdmmm_nf case 2 : \
          return &R3g1ph_2_nf
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_phppp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phppp_G(ep,mpc);}
 
template <class T> complex<T> _R_phppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmpp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phmpp_G(ep,mpc);}
 
template <class T> complex<T> _R_phmpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phmpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phpmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phpmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phpmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phpmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phppm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phppm_G(ep,mpc);}
 
template <class T> complex<T> _R_phppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phmmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phmmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phmmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmpm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phmpm_G(ep,mpc);}
 
template <class T> complex<T> _R_phmpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phmpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phpmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phpmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phpmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phpmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phmmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phmmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phmmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phmmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdppp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdppp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmpp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdmpp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdmpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdpmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdpmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdpmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdpmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdppm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdppm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmmp_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdmmp_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdmmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmpm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdmpm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdmpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdpmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdpmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdpmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdpmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phdmmm_G(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdmmm_G(ep,mpc);}
 
template <class T> complex<T> _R_phdmmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R3g1ph_phdmmm_nf(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R3g1ph_G_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phppp_G;
       _CASE_phmpp_G;
       _CASE_phpmp_G;
       _CASE_phppm_G;
       _CASE_phmmp_G;
       _CASE_phmpm_G;
       _CASE_phpmm_G;
       _CASE_phmmm_G;
       _CASE_phdppp_G;
       _CASE_phdmpp_G;
       _CASE_phdpmp_G;
       _CASE_phdppm_G;
       _CASE_phdmmp_G;
       _CASE_phdmpm_G;
       _CASE_phdpmm_G;
       _CASE_phdmmm_G;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R3g1ph_nf_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phppp_nf;
       _CASE_phmpp_nf;
       _CASE_phpmp_nf;
       _CASE_phppm_nf;
       _CASE_phmmp_nf;
       _CASE_phmpm_nf;
       _CASE_phpmm_nf;
       _CASE_phmmm_nf;
       _CASE_phdppp_nf;
       _CASE_phdmpp_nf;
       _CASE_phdpmp_nf;
       _CASE_phdppm_nf;
       _CASE_phdmmp_nf;
       _CASE_phdmpm_nf;
       _CASE_phdpmm_nf;
       _CASE_phdmmm_nf;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R3g1ph_G_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R3g1ph_G_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R3g1ph_G_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R3g1ph_nf_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R3g1ph_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R3g1ph_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> ( *R3g1ph_G_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R3g1ph_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);

#endif


}
