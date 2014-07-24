/*
* R_2q1g2l.cpp
*
* Created on 12/13, 2008
*      Author: Zvi's script
*/

#include "R_2q1g2l_eval.h"
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



template <class T> complex<T> R2q1g2l_qppqmemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, qm, em, ep}, L}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep) :  qppqmemep L");
#endif

 return( (complex<T>(0,1)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*pow(SPA(1,3),2)*
pow(SPB(5,1),2))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(1,2)*SPA(2,3)*
SPA(4,5))
        );

}

template <class T> complex<T> R2q1g2l_qmpqpemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, em, ep}, L}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep) :  qmpqpemep L");
#endif
// added Harald
 return( (complex<T>(0,1)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*pow(SPB(1,3),2)*
pow(SPA(4,1),2))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(1,2)*SPB(2,3)*
SPB(5,4))
        );

}




template <class T> complex<T> R2q1g2l_qpmqmemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, qm, em, ep}, L}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep) :  qpmqmemep L");
#endif
 
  return( (complex<T>(0,1)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(3,4),2)*
pow(SPB(3,1),2))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(2,1)*SPB(3,2)*
SPB(5,4))
        );

}


template <class T> complex<T> R2q1g2l_qmmqpemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, em, ep}, L}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep) :  qmmqpemep L");
#endif
// added Harald
 return( (complex<T>(0,1)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPB(3,5),2)*
pow(SPA(3,1),2))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(2,1)*SPA(3,2)*
SPA(4,5))
        );

}



template <class T> complex<T> R2q1g2l_qpqmmemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, em, ep}, SLC}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep) :  qpqmmemep SLC");
#endif

 return( complex<T>(0,-1)*(-(((SPA(1,4)*SPA(2,3)+SPA(1,3)*SPA(2,4))*SPA(3,4))/
 (complex<T>(2,0)*SPA(1,2)*SPA(1,3)*SPA(4,5)*SPB(3,2)))+
 (pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*SPA(2,3)*SPA(3,4)*SPB(5,3))/
(pow(SPA(1,2),2)*SPB(2,1)*SPB(3,2))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,3),-1)*pow(SPA(2,3),2)*
 pow(SPB(5,2),2))/(complex<T>(2,0)*pow(SPA(1,3),2)*SPB(3,1)*SPB(3,2)*
 SPB(5,4)))+(complex<T>(0,1)*pow(SPB(5,1),2))/(complex<T>(2,0)*SPB(3,1)*
 SPB(3,2)*SPB(5,4))
        );

}




template <class T> complex<T> R2q1g2l_qpqmpemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, em, ep}, SLC}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep) :  qpqmpemep SLC");
#endif

 return( (complex<T>(0,1)*pow(SPA(2,4),2))/(complex<T>(2,0)*SPA(1,3)*SPA(2,3)*
 SPA(4,5))+complex<T>(0,1)*(-((pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*
  pow(SPA(1,4),2)*pow(SPB(3,1),2))/(complex<T>(2,0)*pow(SPB(3,2),2)*
  SPA(1,3)*SPA(2,3)*SPA(4,5)))-
 (pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*SPA(3,4)*SPB(3,1)*SPB(5,3))/
(pow(SPB(2,1),2)*SPA(1,2)*SPA(1,3))-
 ((SPB(3,2)*SPB(5,1)+SPB(3,1)*SPB(5,2))*SPB(5,3))/
(complex<T>(2,0)*SPA(1,3)*SPB(2,1)*SPB(3,2)*SPB(5,4)))
        );

}




template <class T> complex<T> R2q1g2l_qmqppemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, em, ep}, AX}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep):  qmqppemep AX");
#endif

 complex<T> mt = pow(complex <T> (10.,0), 20);
 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( (complex<T>(1,0)*SPA(1,4)*SPB(3,2)*SPB(5,3))/(mtsq*complex<T>(12,0)*S(4,5))+
 (SPA(1,4)*SPB(3,2)*SPB(5,3))/((S(1,2)-S(4,5))*S(4,5))
        );

}




template <class T> complex<T> R2q1g2l_qmqpmemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, em, ep}, AX}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep):  qmqpmemep AX");
#endif

 complex<T> mt = pow(complex <T> (10.,0), 20);
 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( (complex<T>(1,0)*SPA(1,3)*SPA(3,4)*SPB(5,2))/(mtsq*complex<T>(12,0)*S(4,5))+
 (SPA(1,3)*SPA(3,4)*SPB(5,2))/((S(1,2)-S(4,5))*S(4,5))
        );

}




template <class T> complex<T> R2q1g2l_qpqmpemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, em, ep}, AX}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep) :  qpqmpemep AX");
#endif

 complex<T> mt = pow(complex <T> (10.,0), 20);
 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( (complex<T>(1,0)*SPA(2,4)*SPB(3,1)*SPB(5,3))/(mtsq*complex<T>(12,0)*S(4,5))+
 (SPA(2,4)*SPB(3,1)*SPB(5,3))/((S(1,2)-S(4,5))*S(4,5))
        );

}




template <class T> complex<T> R2q1g2l_qpqmmemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, em, ep}, AX}

#if _VERBOSE
  _MESSAGE("R2q1g2l (ep) :  qpqmmemep AX");
#endif

 complex<T> mt = pow(complex <T> (10.,0), 20);
 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( (complex<T>(1,0)*SPA(2,3)*SPA(3,4)*SPB(5,1))/(mtsq*complex<T>(12,0)*S(4,5))+
 (SPA(2,3)*SPA(3,4)*SPB(5,1))/((S(1,2)-S(4,5))*S(4,5))
        );

}



 // *************** table of switch values *************

#define _R_qppqmemep_L R2q1g2l_7400_L
#define _R_qpmqmemep_L R2q1g2l_7382_L

#define _R_qmmqpemep_L R2q1g2l_7435_L
#define _R_qmpqpemep_L R2q1g2l_7417_L

#define _R_qpqmmemep_SLC R2q1g2l_7352_SLC
#define _R_qpqmpemep_SLC R2q1g2l_7460_SLC
#define _R_qmqppemep_AX R2q1g2l_7465_AX
#define _R_qmqpmemep_AX R2q1g2l_7357_AX
#define _R_qpqmpemep_AX R2q1g2l_7460_AX
#define _R_qpqmmemep_AX R2q1g2l_7352_AX


 // *************** more macro definitions *************

#define _CASE_qppqmemep_L case 7400 : \
          return &R2q1g2l_7400_L
#define _CASE_qpmqmemep_L case 7382 : \
          return &R2q1g2l_7382_L


#define _CASE_qmmqpemep_L case 7435 : \
          return &R2q1g2l_7435_L
#define _CASE_qmpqpemep_L case 7417 : \
          return &R2q1g2l_7417_L


#define _CASE_qpqmmemep_SLC case 7352 : \
          return &R2q1g2l_7352_SLC
#define _CASE_qpqmpemep_SLC case 7460 : \
          return &R2q1g2l_7460_SLC
#define _CASE_qmqppemep_AX case 7465 : \
          return &R2q1g2l_7465_AX
#define _CASE_qmqpmemep_AX case 7357 : \
          return &R2q1g2l_7357_AX
#define _CASE_qpqmpemep_AX case 7460 : \
          return &R2q1g2l_7460_AX
#define _CASE_qpqmmemep_AX case 7352 : \
          return &R2q1g2l_7352_AX


 // *************** function definitions using macros *************

template <class T> complex<T> _R_qppqmemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qppqmemep_L(ep,mpc);}

template <class T> complex<T> _R_qpmqmemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qpmqmemep_L(ep,mpc);}


template <class T> complex<T> _R_qmmqpemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qmmqpemep_L(ep,mpc);}

template <class T> complex<T> _R_qmpqpemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qmpqpemep_L(ep,mpc);}



template <class T> complex<T> _R_qpqmmemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qpqmmemep_SLC(ep,mpc);}

template <class T> complex<T> _R_qpqmpemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qpqmpemep_SLC(ep,mpc);}

template <class T> complex<T> _R_qmqppemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qmqppemep_AX(ep,mpc);}

template <class T> complex<T> _R_qmqpmemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qmqpmemep_AX(ep,mpc);}

template <class T> complex<T> _R_qpqmpemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qpqmpemep_AX(ep,mpc);}

template <class T> complex<T> _R_qpqmmemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g2l_qpqmmemep_AX(ep,mpc);}




 // *************** define pointers *************

template <class T> complex<T> ( *R2q1g2l_AX_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qmqppemep_AX;
       _CASE_qmqpmemep_AX;
       _CASE_qpqmpemep_AX;
       _CASE_qpqmmemep_AX;

       default: return 0;
        }
 }

template <class T> complex<T> ( *R2q1g2l_L_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qppqmemep_L;
       _CASE_qpmqmemep_L;
       _CASE_qmmqpemep_L;
       _CASE_qmpqpemep_L;

       default: return 0;
        }
 }

template <class T> complex<T> ( *R2q1g2l_SLC_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpqmmemep_SLC;
       _CASE_qpqmpemep_SLC;

       default: return 0;
        }
 }


 // *************** definitions for template *************

template complex<R> ( *R2q1g2l_AX_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q1g2l_AX_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q1g2l_AX_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q1g2l_AX_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q1g2l_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q1g2l_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q1g2l_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q1g2l_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q1g2l_SLC_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q1g2l_SLC_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q1g2l_SLC_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q1g2l_SLC_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif




}

