/*
* R_2q2Q_eval.cpp
*
* Created on 11/12, 2010
*      Author: Zvi's script
* convention: first fermion is always the antifermion
*/
 
#include "R_2q2Q_eval.h"
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



 
 
template <class T> complex<T> R2q2Q_qmQpQmqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmQpQmqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,3),2))/(complex<T>(2,0)*SPA(1,4)*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmQmQpqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmQmQpqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,2),2))/(complex<T>(2,0)*SPA(1,4)*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmqpQpQm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmqpQpQm LLT");
#endif
 
 return( (complex<T>(0,-7)*pow(SPA(1,4),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmqpQmQp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmqpQmQp LLT");
#endif
 
 return( (complex<T>(0,-7)*pow(SPA(1,3),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmQpQmqp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmQpQmqp LRT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2Q_qmQmQpqp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmQmQpqp LRT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2Q_qmqpQpQm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmqpQpQm LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,4),2))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmqpQmQp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmqpQmQp LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,3),2))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmqpQpQm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmqpQpQm nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,4),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmqpQmQp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmqpQmQp nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,3),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmqpQpQm_nfLRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm}, nfLRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmqpQpQm nfLRT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,4),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q_qmqpQmQp_nfLRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp}, nfLRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q :  qmqpQmQp nfLRT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,3),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 // *************** table of switch values ************* 
 
#define _R_qmQpQmqp_LLT R2q2Q_607_LLT
#define _R_qmQmQpqp_LLT R2q2Q_637_LLT
#define _R_qmqpQpQm_LLT R2q2Q_1057_LLT
#define _R_qmqpQmQp_LLT R2q2Q_1237_LLT
#define _R_qmQpQmqp_LRT R2q2Q_607_LRT
#define _R_qmQmQpqp_LRT R2q2Q_637_LRT
#define _R_qmqpQpQm_LRT R2q2Q_1057_LRT
#define _R_qmqpQmQp_LRT R2q2Q_1237_LRT
#define _R_qmqpQpQm_nfLLT R2q2Q_1057_nfLLT
#define _R_qmqpQmQp_nfLLT R2q2Q_1237_nfLLT
#define _R_qmqpQpQm_nfLRT R2q2Q_1057_nfLRT
#define _R_qmqpQmQp_nfLRT R2q2Q_1237_nfLRT
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmQpQmqp_LLT case 607 : \
          return &R2q2Q_607_LLT
#define _CASE_qmQmQpqp_LLT case 637 : \
          return &R2q2Q_637_LLT
#define _CASE_qmqpQpQm_LLT case 1057 : \
          return &R2q2Q_1057_LLT
#define _CASE_qmqpQmQp_LLT case 1237 : \
          return &R2q2Q_1237_LLT
#define _CASE_qmQpQmqp_LRT case 607 : \
          return &R2q2Q_607_LRT
#define _CASE_qmQmQpqp_LRT case 637 : \
          return &R2q2Q_637_LRT
#define _CASE_qmqpQpQm_LRT case 1057 : \
          return &R2q2Q_1057_LRT
#define _CASE_qmqpQmQp_LRT case 1237 : \
          return &R2q2Q_1237_LRT
#define _CASE_qmqpQpQm_nfLLT case 1057 : \
          return &R2q2Q_1057_nfLLT
#define _CASE_qmqpQmQp_nfLLT case 1237 : \
          return &R2q2Q_1237_nfLLT
#define _CASE_qmqpQpQm_nfLRT case 1057 : \
          return &R2q2Q_1057_nfLRT
#define _CASE_qmqpQmQp_nfLRT case 1237 : \
          return &R2q2Q_1237_nfLRT
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_qmQpQmqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmQpQmqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmQmQpqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmqpQpQm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmqpQmQp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQpQmqp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmQpQmqp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpqp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmQmQpqp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmqpQpQm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmqpQmQp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmqpQpQm_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmqpQmQp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQm_nfLRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmqpQpQm_nfLRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQp_nfLRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q_qmqpQmQp_nfLRT(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q2Q_LLT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmQpQmqp_LLT;
       _CASE_qmQmQpqp_LLT;
       _CASE_qmqpQpQm_LLT;
       _CASE_qmqpQmQp_LLT;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q_LRT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmQpQmqp_LRT;
       _CASE_qmQmQpqp_LRT;
       _CASE_qmqpQpQm_LRT;
       _CASE_qmqpQmQp_LRT;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q_nfLLT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmqpQpQm_nfLLT;
       _CASE_qmqpQmQp_nfLLT;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q_nfLRT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmqpQpQm_nfLRT;
       _CASE_qmqpQmQp_nfLRT;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q2Q_LLT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q_LLT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q_LLT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2Q_LRT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q_LRT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q_LRT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2Q_nfLLT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q_nfLLT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q_nfLLT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2Q_nfLRT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q_nfLRT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q_nfLRT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> ( *R2q2Q_LLT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2Q_LRT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2Q_nfLLT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2Q_nfLRT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);

#endif



}
