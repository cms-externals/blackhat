/*
* R_2q2Q1ph_eval.cpp
*
* Created on 12/7, 2010
*      Author: Zvi's script
*/
 
#include "R_2q2Q1ph_eval.h"
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


 
 
template <class T> complex<T> R2q2Q1ph_phqmqpQpQm_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, Qp, Qm}, lc}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1ph :  phqmqpQpQm lc");
#endif
 
 return( (complex<T>(-20,0)*pow(SPA(2,5),2))/(complex<T>(9,0)*SPA(2,3)*SPA(4,5))-
 (pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1)*SPA(2,3))/
(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(4,5))-
 (pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*SPA(4,5))/
(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,3))+(complex<T>(2,0)*pow(SPB(4,3),2))/
(SPB(3,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1ph_phqmqpQpQm_slc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, Qp, Qm}, slc}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1ph :  phqmqpQpQm slc");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(2,5),2))/(SPA(2,3)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1)*
 SPA(2,3))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*
 SPA(4,5))/(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1ph_phqmqpQpQm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, Qp, Qm}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1ph :  phqmqpQpQm nf");
#endif
 
 return( (complex<T>(-4,0)*pow(SPA(2,5),2))/(complex<T>(9,0)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1ph_phqmqpQmQp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, Qm, Qp}, lc}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1ph :  phqmqpQmQp lc");
#endif
 
 return( (complex<T>(-20,0)*pow(SPA(2,4),2))/(complex<T>(9,0)*SPA(2,3)*SPA(4,5))-
 (pow(SPA(3,4),2)*pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),
-1)*SPA(2,3))/(complex<T>(2,0)*pow(SPA(3,5),2)*pow(SPB(5,4),2)*
 SPA(4,5))-(pow(SPA(2,5),2)*pow(SPB(5,2),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1)*SPA(4,5))/
(complex<T>(2,0)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*SPA(2,3))-
 (S(2,5)*SPA(4,5))/(complex<T>(2,0)*pow(SPA(3,5),2)*SPB(3,2))+
 (complex<T>(1,0)*SPA(4,5)*SPB(5,3))/(complex<T>(2,0)*SPA(3,5)*SPB(3,2))-
 (S(3,4)*SPA(2,3))/(complex<T>(2,0)*pow(SPA(3,5),2)*SPB(5,4))+
 (complex<T>(2,0)*pow(SPB(5,3),2))/(SPB(3,2)*SPB(5,4))+
 (complex<T>(1,0)*SPA(2,3)*SPB(5,3))/(complex<T>(2,0)*SPA(3,5)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1ph_phqmqpQmQp_slc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, Qm, Qp}, slc}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1ph :  phqmqpQmQp slc");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(2,4),2))/(SPA(2,3)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1)*
 SPA(2,3))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1)*
 SPA(4,5))/(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1ph_phqmqpQmQp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, Qm, Qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1ph :  phqmqpQmQp nf");
#endif
 
 return( (complex<T>(-4,0)*pow(SPA(2,4),2))/(complex<T>(9,0)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 // *************** table of switch values ************* 
 
#define _R_phqmqpQpQm_lc R2q2Q1ph_6356_lc
#define _R_phqmqpQpQm_slc R2q2Q1ph_6356_slc
#define _R_phqmqpQpQm_nf R2q2Q1ph_6356_nf
#define _R_phqmqpQmQp_lc R2q2Q1ph_7436_lc
#define _R_phqmqpQmQp_slc R2q2Q1ph_7436_slc
#define _R_phqmqpQmQp_nf R2q2Q1ph_7436_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_phqmqpQpQm_lc case 6356 : \
          return &R2q2Q1ph_6356_lc
#define _CASE_phqmqpQpQm_slc case 6356 : \
          return &R2q2Q1ph_6356_slc
#define _CASE_phqmqpQpQm_nf case 6356 : \
          return &R2q2Q1ph_6356_nf
#define _CASE_phqmqpQmQp_lc case 7436 : \
          return &R2q2Q1ph_7436_lc
#define _CASE_phqmqpQmQp_slc case 7436 : \
          return &R2q2Q1ph_7436_slc
#define _CASE_phqmqpQmQp_nf case 7436 : \
          return &R2q2Q1ph_7436_nf
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_phqmqpQpQm_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1ph_phqmqpQpQm_lc(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpQpQm_slc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1ph_phqmqpQpQm_slc(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpQpQm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1ph_phqmqpQpQm_nf(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpQmQp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1ph_phqmqpQmQp_lc(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpQmQp_slc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1ph_phqmqpQmQp_slc(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpQmQp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1ph_phqmqpQmQp_nf(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q2Q1ph_lc_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phqmqpQpQm_lc;
       _CASE_phqmqpQmQp_lc;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q1ph_nf_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phqmqpQpQm_nf;
       _CASE_phqmqpQmQp_nf;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q1ph_slc_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phqmqpQpQm_slc;
       _CASE_phqmqpQmQp_slc;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q2Q1ph_lc_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q1ph_lc_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1ph_lc_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2Q1ph_nf_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q1ph_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1ph_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2Q1ph_slc_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q1ph_slc_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1ph_slc_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q1ph_lc_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2Q1ph_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2Q1ph_slc_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);


#endif


}
