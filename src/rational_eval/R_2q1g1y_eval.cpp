/*
* R_2q1g1y.cpp
*
* Created on 12/13, 2008
*      Author: Zvi's script
*/

#include "R_2q1g1y_eval.h"
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





template <class T> complex<T> R2q1g1y_qmpqpgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, gap}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmpqpgap L");
#endif

 return( (complex<T>(0,-1)*SPA(1,3)*SPB(4,2))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4))
        );

}




template <class T> complex<T> R2q1g1y_qmpqpgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, gam}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmpqpgam L");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmmqpgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, gap}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmmqpgap L");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmmqpgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, gam}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmmqpgam L");
#endif

 return( (complex<T>(0,1)*SPA(2,4)*SPB(3,1))/(complex<T>(2,0)*SPB(2,1)*SPB(4,1))
        );

}




template <class T> complex<T> R2q1g1y_qmgapqpp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gap, qp, p}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmgapqpp L");
#endif

 return( (complex<T>(0,1)*SPA(1,3)*SPB(4,2))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4))
        );

}




template <class T> complex<T> R2q1g1y_qmgapqpm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gap, qp, m}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmgapqpm L");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmgamqpp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, p}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmgamqpp L");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmgamqpm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, m}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmgamqpm L");
#endif

 return( (complex<T>(0,-1)*SPA(2,4)*SPB(3,1))/(complex<T>(2,0)*SPB(2,1)*SPB(4,1))
        );

}




template <class T> complex<T> R2q1g1y_qmpqpgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, gap}, nf}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmpqpgap nf");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmpqpgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, gam}, nf}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmpqpgam nf");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmmqpgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, gap}, nf}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmmqpgap nf");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmmqpgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, gam}, nf}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmmqpgam nf");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmgapqpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gap, qp, p}, nf}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmgapqpp nf");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmgapqpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gap, qp, m}, nf}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmgapqpm nf");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmgamqpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, p}, nf}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmgamqpp nf");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmgamqpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, m}, nf}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmgamqpm nf");
#endif

   return( complex<T>(0,0));

}


template <class T> complex<T> R2q1g1y_qmqppgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, gap}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmqppgap L");
#endif

 return( (complex<T>(0,1)*SPA(1,2)*SPB(3,2))/(complex<T>(2,0)*SPA(2,4)*SPA(3,4))+
 (complex<T>(0,-1)*SPA(1,2)*SPB(4,2))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4))
        );

}




template <class T> complex<T> R2q1g1y_qmqppgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, gam}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmqppgam L");
#endif

 return( (complex<T>(0,-1)*pow(SPA(1,4),3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,3)*
 SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,4),2)*SPA(2,4))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*SPB(2,1))/(complex<T>(2,0)*SPA(1,3)*SPA(3,4)*
 SPB(4,1))
        );

}




template <class T> complex<T> R2q1g1y_qmqpmgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, gap}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmqpmgap L");
#endif

 return( (complex<T>(0,1)*pow(SPA(1,3),3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*
 SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPA(2,3))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,4)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*SPB(2,1))/(complex<T>(2,0)*SPA(1,4)*SPA(3,4)*
 SPB(3,1))
        );

}




template <class T> complex<T> R2q1g1y_qmqpmgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, gam}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmqpmgam L");
#endif

 return( (complex<T>(0,-1)*SPA(1,4)*SPB(2,1))/(complex<T>(2,0)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,1)*SPA(1,3)*SPB(2,1))/(complex<T>(2,0)*SPB(4,1)*SPB(4,3))
        );

}




template <class T> complex<T> R2q1g1y_qmqpgapp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, gap, p}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmqpgapp L");
#endif

 return( (complex<T>(0,1)*SPA(1,2)*SPB(3,2))/(complex<T>(2,0)*SPA(2,4)*SPA(3,4))+
 (complex<T>(0,-1)*SPA(1,2)*SPB(4,2))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4))
        );

}




template <class T> complex<T> R2q1g1y_qmqpgapm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, gap, m}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmqpgapm L");
#endif

 return( (complex<T>(0,1)*pow(SPA(1,3),3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*
 SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPA(2,3))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,4)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*SPB(2,1))/(complex<T>(2,0)*SPA(1,4)*SPA(3,4)*
 SPB(3,1))
        );

}




template <class T> complex<T> R2q1g1y_qmqpgamp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, gam, p}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmqpgamp L");
#endif

 return( (complex<T>(0,1)*pow(SPA(1,3),3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*
 SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPA(2,3))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,4)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*SPB(2,1))/(complex<T>(2,0)*SPA(1,4)*SPA(3,4)*
 SPB(3,1))
        );

}




template <class T> complex<T> R2q1g1y_qmqpgamm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, gam, m}, L}

#if _VERBOSE
  _MESSAGE("R2q1g1y :  qmqpgamm L");
#endif

 return( (complex<T>(0,-1)*SPA(1,4)*SPB(2,1))/(complex<T>(2,0)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,1)*SPA(1,3)*SPB(2,1))/(complex<T>(2,0)*SPB(4,1)*SPB(4,3))
        );

}



 // *************** table of switch values *************

#define _R_qmpqpgap_L R2q1g1y_1171_L
#define _R_qmpqpgam_L R2q1g1y_955_L
#define _R_qmmqpgap_L R2q1g1y_1153_L
#define _R_qmmqpgam_L R2q1g1y_937_L
#define _R_qmgapqpp_L R2q1g1y_751_L
#define _R_qmgapqpm_L R2q1g1y_103_L
#define _R_qmgamqpp_L R2q1g1y_745_L
#define _R_qmgamqpm_L R2q1g1y_97_L
#define _R_qmpqpgap_nf R2q1g1y_1171_nf
#define _R_qmpqpgam_nf R2q1g1y_955_nf
#define _R_qmmqpgap_nf R2q1g1y_1153_nf
#define _R_qmmqpgam_nf R2q1g1y_937_nf
#define _R_qmgapqpp_nf R2q1g1y_751_nf
#define _R_qmgapqpm_nf R2q1g1y_103_nf
#define _R_qmgamqpp_nf R2q1g1y_745_nf
#define _R_qmgamqpm_nf R2q1g1y_97_nf
#define _R_qmqppgap_L R2q1g1y_1201_L
#define _R_qmqppgam_L R2q1g1y_985_L
#define _R_qmqpmgap_L R2q1g1y_1093_L
#define _R_qmqpmgam_L R2q1g1y_877_L
#define _R_qmqpgapp_L R2q1g1y_841_L
#define _R_qmqpgapm_L R2q1g1y_193_L
#define _R_qmqpgamp_L R2q1g1y_805_L
#define _R_qmqpgamm_L R2q1g1y_157_L


 // *************** more macro definitions *************

#define _CASE_qmpqpgap_L case 1171 : \
          return &R2q1g1y_1171_L
#define _CASE_qmpqpgam_L case 955 : \
          return &R2q1g1y_955_L
#define _CASE_qmmqpgap_L case 1153 : \
          return &R2q1g1y_1153_L
#define _CASE_qmmqpgam_L case 937 : \
          return &R2q1g1y_937_L
#define _CASE_qmgapqpp_L case 751 : \
          return &R2q1g1y_751_L
#define _CASE_qmgapqpm_L case 103 : \
          return &R2q1g1y_103_L
#define _CASE_qmgamqpp_L case 745 : \
          return &R2q1g1y_745_L
#define _CASE_qmgamqpm_L case 97 : \
          return &R2q1g1y_97_L
#define _CASE_qmpqpgap_nf case 1171 : \
          return &R2q1g1y_1171_nf
#define _CASE_qmpqpgam_nf case 955 : \
          return &R2q1g1y_955_nf
#define _CASE_qmmqpgap_nf case 1153 : \
          return &R2q1g1y_1153_nf
#define _CASE_qmmqpgam_nf case 937 : \
          return &R2q1g1y_937_nf
#define _CASE_qmgapqpp_nf case 751 : \
          return &R2q1g1y_751_nf
#define _CASE_qmgapqpm_nf case 103 : \
          return &R2q1g1y_103_nf
#define _CASE_qmgamqpp_nf case 745 : \
          return &R2q1g1y_745_nf
#define _CASE_qmgamqpm_nf case 97 : \
          return &R2q1g1y_97_nf
#define _CASE_qmqppgap_L case 1201 : \
          return &R2q1g1y_1201_L
#define _CASE_qmqppgam_L case 985 : \
          return &R2q1g1y_985_L
#define _CASE_qmqpmgap_L case 1093 : \
          return &R2q1g1y_1093_L
#define _CASE_qmqpmgam_L case 877 : \
          return &R2q1g1y_877_L
#define _CASE_qmqpgapp_L case 841 : \
          return &R2q1g1y_841_L
#define _CASE_qmqpgapm_L case 193 : \
          return &R2q1g1y_193_L
#define _CASE_qmqpgamp_L case 805 : \
          return &R2q1g1y_805_L
#define _CASE_qmqpgamm_L case 157 : \
          return &R2q1g1y_157_L


 // *************** function definitions using macros *************

template <class T> complex<T> _R_qmpqpgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmpqpgap_L(ep,mpc);}

template <class T> complex<T> _R_qmpqpgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmpqpgam_L(ep,mpc);}

template <class T> complex<T> _R_qmmqpgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmmqpgap_L(ep,mpc);}

template <class T> complex<T> _R_qmmqpgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmmqpgam_L(ep,mpc);}

template <class T> complex<T> _R_qmgapqpp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmgapqpp_L(ep,mpc);}

template <class T> complex<T> _R_qmgapqpm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmgapqpm_L(ep,mpc);}

template <class T> complex<T> _R_qmgamqpp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmgamqpp_L(ep,mpc);}

template <class T> complex<T> _R_qmgamqpm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmgamqpm_L(ep,mpc);}

template <class T> complex<T> _R_qmpqpgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmpqpgap_nf(ep,mpc);}

template <class T> complex<T> _R_qmpqpgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmpqpgam_nf(ep,mpc);}

template <class T> complex<T> _R_qmmqpgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmmqpgap_nf(ep,mpc);}

template <class T> complex<T> _R_qmmqpgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmmqpgam_nf(ep,mpc);}

template <class T> complex<T> _R_qmgapqpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmgapqpp_nf(ep,mpc);}

template <class T> complex<T> _R_qmgapqpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmgapqpm_nf(ep,mpc);}

template <class T> complex<T> _R_qmgamqpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmgamqpp_nf(ep,mpc);}

template <class T> complex<T> _R_qmgamqpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmgamqpm_nf(ep,mpc);}

template <class T> complex<T> _R_qmqppgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmqppgap_L(ep,mpc);}

template <class T> complex<T> _R_qmqppgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmqppgam_L(ep,mpc);}

template <class T> complex<T> _R_qmqpmgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmqpmgap_L(ep,mpc);}

template <class T> complex<T> _R_qmqpmgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmqpmgam_L(ep,mpc);}

template <class T> complex<T> _R_qmqpgapp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmqpgapp_L(ep,mpc);}

template <class T> complex<T> _R_qmqpgapm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmqpgapm_L(ep,mpc);}

template <class T> complex<T> _R_qmqpgamp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmqpgamp_L(ep,mpc);}

template <class T> complex<T> _R_qmqpgamm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q1g1y_qmqpgamm_L(ep,mpc);}




 // *************** define pointers *************

template <class T> complex<T> ( *R2q1g1y_SLC_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qmpqpgap_L;
       _CASE_qmpqpgam_L;
       _CASE_qmmqpgap_L;
       _CASE_qmmqpgam_L;
       _CASE_qmgapqpp_L;
       _CASE_qmgapqpm_L;
       _CASE_qmgamqpp_L;
       _CASE_qmgamqpm_L;
       _CASE_qmqppgap_L;
       _CASE_qmqppgam_L;
       _CASE_qmqpmgap_L;
       _CASE_qmqpmgam_L;
       _CASE_qmqpgapp_L;
       _CASE_qmqpgapm_L;
       _CASE_qmqpgamp_L;
       _CASE_qmqpgamm_L;

       default: return 0;
        }
 }

template <class T> complex<T> ( *R2q1g1y_L_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qmpqpgap_L;
       _CASE_qmpqpgam_L;
       _CASE_qmmqpgap_L;
       _CASE_qmmqpgam_L;
       _CASE_qmgapqpp_L;
       _CASE_qmgapqpm_L;
       _CASE_qmgamqpp_L;
       _CASE_qmgamqpm_L;
       _CASE_qmqppgap_L;
       _CASE_qmqppgam_L;
       _CASE_qmqpmgap_L;
       _CASE_qmqpmgam_L;
       _CASE_qmqpgapp_L;
       _CASE_qmqpgapm_L;
       _CASE_qmqpgamp_L;
       _CASE_qmqpgamm_L;

       default: return 0;
        }
 }

template <class T> complex<T> ( *R2q1g1y_nf_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qmpqpgap_nf;
       _CASE_qmpqpgam_nf;
       _CASE_qmmqpgap_nf;
       _CASE_qmmqpgam_nf;
       _CASE_qmgapqpp_nf;
       _CASE_qmgapqpm_nf;
       _CASE_qmgamqpp_nf;
       _CASE_qmgamqpm_nf;

       default: return 0;
        }
 }


 // *************** definitions for template *************

template complex<R> ( *R2q1g1y_SLC_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q1g1y_SLC_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q1g1y_SLC_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q1g1y_SLC_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif

template complex<R> ( *R2q1g1y_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q1g1y_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q1g1y_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q1g1y_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q1g1y_nf_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q1g1y_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q1g1y_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q1g1y_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif




}

