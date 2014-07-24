/*
* R_2q2g_eval.cpp
*
* Created on 11/12, 2010
*      Author: Zvi's script
*/
 
#include "R_2q2g_eval.h"
#include "eval_param.h"
 
using namespace std;
 
namespace BH  {
 
 
#define _VERBOSE 0
 
#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}


 
 
template <class T> complex<T> R2q2g_qmpmqp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, qp}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmpmqp LT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,3),3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*
 SPA(2,3))+(complex<T>(0,-1)*pow(SPA(1,3),3)*S(1,4))/
(complex<T>(2,0)*S(2,4)*SPA(1,2)*SPA(1,4)*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g_qmmpqp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, qp}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmmpqp LT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,2),2)*SPA(2,4))/(complex<T>(2,0)*SPA(1,4)*SPA(2,3)*
SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g_qmpqpm_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, m}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmpqpm LT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g_qmmqpp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, p}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmmqpp LT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g_qmqppm_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, m}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmqppm LT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,4),2)*SPA(2,4))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*
SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g_qmqpmp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, p}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmqpmp LT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,3),3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*
 SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPB(2,1))/
(complex<T>(2,0)*SPA(1,4)*SPA(3,4)*SPB(3,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g_qmpqpm_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmpqpm nfLT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g_qmmqpp_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmmqpp nfLT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g_qmqppm_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmqppm nfLT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g_qmqpmp_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g :  qmqpmp nfLT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 // *************** table of switch values ************* 
 
#define _R_qmpmqp_LT R2q2g_141_LT
#define _R_qmmpqp_LT R2q2g_177_LT
#define _R_qmpqpm_LT R2q2g_45_LT
#define _R_qmmqpp_LT R2q2g_225_LT
#define _R_qmqppm_LT R2q2g_57_LT
#define _R_qmqpmp_LT R2q2g_201_LT
#define _R_qmpqpm_nfLT R2q2g_45_nfLT
#define _R_qmmqpp_nfLT R2q2g_225_nfLT
#define _R_qmqppm_nfLT R2q2g_57_nfLT
#define _R_qmqpmp_nfLT R2q2g_201_nfLT
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmpmqp_LT case 141 : \
          return &R2q2g_141_LT
#define _CASE_qmmpqp_LT case 177 : \
          return &R2q2g_177_LT
#define _CASE_qmpqpm_LT case 45 : \
          return &R2q2g_45_LT
#define _CASE_qmmqpp_LT case 225 : \
          return &R2q2g_225_LT
#define _CASE_qmqppm_LT case 57 : \
          return &R2q2g_57_LT
#define _CASE_qmqpmp_LT case 201 : \
          return &R2q2g_201_LT
#define _CASE_qmpqpm_nfLT case 45 : \
          return &R2q2g_45_nfLT
#define _CASE_qmmqpp_nfLT case 225 : \
          return &R2q2g_225_nfLT
#define _CASE_qmqppm_nfLT case 57 : \
          return &R2q2g_57_nfLT
#define _CASE_qmqpmp_nfLT case 201 : \
          return &R2q2g_201_nfLT
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_qmpmqp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmpmqp_LT(ep,mpc);}
 
template <class T> complex<T> _R_qmmpqp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmmpqp_LT(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpm_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmpqpm_LT(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmmqpp_LT(ep,mpc);}
 
template <class T> complex<T> _R_qmqppm_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmqppm_LT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmqpmp_LT(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpm_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmpqpm_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpp_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmmqpp_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqppm_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmqppm_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmp_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g_qmqpmp_nfLT(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q2g_LT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmpmqp_LT;
       _CASE_qmmpqp_LT;
       _CASE_qmpqpm_LT;
       _CASE_qmmqpp_LT;
       _CASE_qmqppm_LT;
       _CASE_qmqpmp_LT;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g_nfLT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmpqpm_nfLT;
       _CASE_qmmqpp_nfLT;
       _CASE_qmqppm_nfLT;
       _CASE_qmqpmp_nfLT;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q2g_LT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2g_LT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2g_LT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2g_nfLT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2g_nfLT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2g_nfLT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


#if BH_USE_GMP
template complex<RGMP> ( *R2q2g_LT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2g_nfLT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
#endif

}
