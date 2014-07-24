/*
* R_2q2Q1y_eval.cpp
*
* Created on 2/6, 2011
*      Author: Zvi's script
*/
 
#include "R_2q2Q1y_eval.h"
#include "eval_param.h"
 
using namespace std;
 
namespace BH  {
 
 
#define _VERBOSE 0

 
#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)


 
 
template <class T> complex<T> R2q2Q1y_qmQpQmqpgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmQpQmqpgam nf");
#endif
 
 return( (complex<T>(0,2)*pow(SPB(4,2),2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmQmQpqpgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmQmQpqpgam nf");
#endif
 
 return( (complex<T>(0,2)*pow(SPB(4,3),2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpQmQpqmgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qm, Qp, qm, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpQmQpqmgam nf");
#endif
 
 return( (complex<T>(0,2)*pow(SPB(3,1),2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpQpQmqmgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qp, Qm, qm, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpQpQmqmgam nf");
#endif
 
 return( (complex<T>(0,2)*pow(SPB(2,1),2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmQpQmqpgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmQpQmqpgam L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPA(4,5))/
 (complex<T>(2,0)*S(1,2)*S(3,4)*SPB(5,1))+
(complex<T>(1,0)*S(4,5)*SPA(1,5)*SPA(3,5)*SPB(2,1))/
 (complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPA(4,5)*SPB(5,1))+
(complex<T>(1,0)*SPA(1,3)*SPA(1,5)*SPB(4,1))/(complex<T>(2,0)*S(4,5)*SPA(1,2)*
SPB(5,1))+(SPA(3,5)*SPB(4,2))/(S(2,3)*SPB(5,1))+
(complex<T>(1,0)*SPA(3,5)*SPB(3,2)*SPB(4,1)*SPB(4,2))/
 (complex<T>(2,0)*S(4,5)*SPB(2,1)*SPB(4,3)*SPB(5,1))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(3,5),2)*
SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*SPB(4,3)*SPB(5,1))+
(complex<T>(1,0)*(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)+
 pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))*pow(SPA(3,5),2)*
pow(SPB(3,2),2)*pow(SPB(5,4),2)*SPB(3,1))/
 (complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(4,3)*SPB(5,1)*
SPB(5,3))+(((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
pow(SPA(3,5),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
 (pow(SPA(1,2),3)*pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-
 S(4,5)/S(1,2))*SPB(5,4))+(complex<T>(1,0)*pow(SPB(4,2),2)*SPA(3,5))/
 (complex<T>(2,0)*pow(SPB(2,1),2)*SPA(1,2)*SPB(5,4))-
(pow(SPB(4,2),2)*SPA(3,5))/(complex<T>(2,0)*S(3,4)*SPB(2,1)*SPB(5,4))-
(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*SPA(3,5)*SPB(4,2)*
(S(3,5)*SPB(4,2)+complex<T>(-2,0)*SPA(1,3)*SPB(2,1)*SPB(4,3)))/
 (complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(5,4))+
(complex<T>(7,0)*pow(SPB(4,2),2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))+
(complex<T>(1,0)*pow(SPA(1,3),2)*SPB(4,1))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*
SPB(5,1)*SPB(5,4))-(S(4,5)*SPA(1,3)*SPB(4,1)*SPB(4,2))/
 (complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPB(5,1)*SPB(5,4))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(3,5),2)*
pow(SPB(3,2),2)*SPB(4,1))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(2,1)*
SPB(4,3)*SPB(5,1)*SPB(5,4))+(complex<T>(1,0)*pow(SPA(4,5),2)*
pow(SPB(4,2),2)*SPB(4,1)*SPB(5,4))/(complex<T>(2,0)*pow(SPB(2,1),2)*
pow(SPB(4,3),2)*SPA(1,2)*SPA(3,4)*SPB(5,1))-
(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*SPA(1,5)*SPA(3,5)*SPB(3,1)*
SPB(4,2)*SPB(5,4))/(pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPB(4,3)*
SPB(5,1))+(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
(-((pow(SPA(3,5),2)*SPB(3,2)*SPB(4,2))/(complex<T>(2,0)*SPB(2,1)))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,2),2)*SPB(4,1)*SPB(5,4))/
  (complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,1))))/(pow(SPA(3,4),2)*
pow(SPB(4,3),2)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmQmQpqpgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmQmQpqpgam L");
#endif
 
 return( complex<T>(0,1)*(-((S(4,5)*SPA(1,5)*SPA(2,5)*SPB(3,1))/
(complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPA(4,5)*SPB(5,1)))-
(SPA(2,5)*SPB(4,3))/(S(2,3)*SPB(5,1))+(complex<T>(7,0)*pow(SPB(4,3),2))/
 (complex<T>(9,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))+
(complex<T>(1,0)*S(4,5)*SPA(1,2)*SPB(4,1)*SPB(4,3))/
 (complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpQmQpqmgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qm, Qp, qm, gam}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpQmQpqmgam L");
#endif
 
 return( complex<T>(0,-1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(2,5),3)*pow(SPB(2,1),2)*pow(SPB(5,3),2))/
(pow(SPA(3,4),3)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(1,5)/S(3,4))*SPB(5,1)))-(pow(SPB(3,1),2)*SPA(2,5))/
 (complex<T>(2,0)*pow(SPB(4,3),2)*SPA(3,4)*SPB(5,1))+
(complex<T>(1,0)*pow(SPB(3,1),2)*SPA(2,5))/(complex<T>(2,0)*S(1,2)*SPB(4,3)*
SPB(5,1))-(pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*SPA(2,5)*SPB(3,1)*
(-(S(2,5)*SPB(3,1))+complex<T>(2,0)*SPA(2,4)*SPB(2,1)*SPB(4,3)))/
 (complex<T>(2,0)*pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPB(5,1))+
(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
((complex<T>(1,0)*pow(SPA(2,5),2)*SPB(3,1)*SPB(3,2))/
  (complex<T>(2,0)*SPB(4,3))-(pow(SPA(1,5),2)*pow(SPB(3,1),2)*
   SPB(4,1)*SPB(5,1))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,4))))/
 (pow(SPA(1,2),2)*pow(SPB(2,1),2))-
(pow(SPB(3,1),2)*SPA(1,5)*SPA(2,3))/(complex<T>(2,0)*S(1,2)*S(3,4)*
SPB(5,4))-(SPA(2,5)*SPB(3,1))/(S(2,3)*SPB(5,4))-
(SPA(2,4)*SPA(4,5)*SPB(4,1))/(complex<T>(2,0)*S(1,5)*SPA(3,4)*SPB(5,4))-
(SPA(2,5)*SPB(3,1)*SPB(3,2)*SPB(4,1))/(complex<T>(2,0)*S(1,5)*SPB(2,1)*
SPB(4,3)*SPB(5,4))-(S(1,5)*SPA(2,5)*SPA(4,5)*SPB(4,3))/
 (complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPA(1,5)*SPB(5,4))+
(complex<T>(-7,0)*pow(SPB(3,1),2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))-
(pow(SPA(2,4),2)*SPB(4,1))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*SPB(5,1)*
SPB(5,4))+(complex<T>(1,0)*S(1,5)*SPA(2,4)*SPB(3,1)*SPB(4,1))/
 (complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPB(5,1)*SPB(5,4))-
(pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*
SPB(4,1))/(complex<T>(2,0)*pow(SPA(1,5),2)*SPB(2,1)*SPB(4,3)*SPB(5,1)*
SPB(5,4))-(pow(SPA(1,5),2)*pow(SPB(3,1),2)*SPB(4,1)*
SPB(5,1))/(complex<T>(2,0)*pow(SPB(2,1),2)*pow(SPB(4,3),2)*SPA(1,2)*
SPA(3,4)*SPB(5,4))+(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*SPA(2,5)*
SPA(4,5)*SPB(3,1)*SPB(4,2)*SPB(5,1))/(pow(SPA(3,4),2)*
pow(SPB(4,3),2)*SPB(2,1)*SPB(5,4))-
((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)+pow(complex<T>(1,0)-S(1,5)/S(3,4),
 -1))*pow(SPA(2,5),2)*pow(SPB(3,2),2)*pow(SPB(5,1),2)*
SPB(4,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPB(2,1)*
SPB(5,2)*SPB(5,4))-(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
pow(SPA(2,5),2)*SPB(5,3))/(complex<T>(2,0)*pow(SPA(1,2),2)*SPB(2,1)*
SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpQpQmqmgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qp, Qm, qm, gam}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpQpQmqmgam L");
#endif
 
 return( complex<T>(0,-1)*((SPA(3,5)*SPB(2,1))/(S(2,3)*SPB(5,4))+
(complex<T>(1,0)*S(1,5)*SPA(3,5)*SPA(4,5)*SPB(4,2))/
 (complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPA(1,5)*SPB(5,4))+
(complex<T>(-7,0)*pow(SPB(2,1),2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))-
(S(1,5)*SPA(3,4)*SPB(2,1)*SPB(4,1))/(complex<T>(2,0)*S(2,3)*
(-S(1,5)+S(2,3))*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmqpQmQpgam_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, gam}, sl}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmqpQmQpgam sl");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(4,2),2))/(complex<T>(2,0)*SPB(4,3)*SPB(5,1)*
 SPB(5,2))+complex<T>(0,1)*
((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(1,5),3)*pow(SPB(2,1),2)*pow(SPB(5,4),2))/
(pow(SPA(3,4),3)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(2,5)/S(3,4))*SPB(5,2))+(complex<T>(1,0)*pow(SPB(4,2),2)*SPA(1,5))/
(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(3,4)*SPB(5,2))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1)*SPA(1,5)*SPB(4,2)*
 (-(S(1,5)*SPB(4,2))+complex<T>(2,0)*SPA(1,3)*SPB(2,1)*SPB(4,3)))/
(complex<T>(2,0)*pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPB(5,2))-
 (pow(SPB(4,2),2)*SPB(3,2))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,2)*
 SPB(5,3)))+complex<T>(0,-1)*
(-((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(3,5),2)*SPB(5,2))/
 (complex<T>(2,0)*pow(SPA(3,4),2)*SPB(4,3)*SPB(5,1)))-
 (pow(SPB(4,2),2)*SPB(4,1))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,1)*
 SPB(5,4)))+(complex<T>(0,-1)*pow(SPB(4,2),2))/(complex<T>(2,0)*SPB(2,1)*
 SPB(5,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmqpQpQmgam_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, gam}, sl}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmqpQpQmgam sl");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(3,2),2))/(complex<T>(2,0)*SPB(4,3)*SPB(5,1)*
 SPB(5,2))+complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
 pow(SPA(4,5),2)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*SPB(4,3)*
 SPB(5,1))+(complex<T>(1,0)*pow(SPB(3,2),2)*SPB(3,1))/
(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,1)*SPB(5,3)))+
 complex<T>(0,-1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
   (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1))/complex<T>(2,0))*
  pow(SPA(1,5),3)*pow(SPB(2,1),2)*pow(SPB(5,3),2))/
 (pow(SPA(3,4),3)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
   S(2,5)/S(3,4))*SPB(5,2)))-(pow(SPB(3,2),2)*SPA(1,5))/
(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(3,4)*SPB(5,2))-
 (pow(complex<T>(1,0)-S(2,5)/S(3,4),-1)*SPA(1,5)*SPB(3,2)*
 (-(S(1,5)*SPB(3,2))+complex<T>(-2,0)*SPA(1,4)*SPB(2,1)*SPB(4,3)))/
(complex<T>(2,0)*pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPB(4,2))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*
 SPB(5,2)*SPB(5,4)))+(complex<T>(0,1)*pow(SPB(3,2),2))/
(complex<T>(2,0)*SPB(2,1)*SPB(5,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpqmQmQpgam_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, Qm, Qp, gam}, sl}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpqmQmQpgam sl");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(4,1),2))/(complex<T>(2,0)*SPB(4,3)*SPB(5,1)*
 SPB(5,2))+complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
 pow(SPA(3,5),2)*SPB(5,1))/(complex<T>(2,0)*S(3,4)*SPA(3,4)*SPB(5,2))-
 (pow(SPB(4,1),2)*SPB(3,2))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,2)*
 SPB(5,3)))+complex<T>(0,-1)*
(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
   (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
  pow(SPA(2,5),3)*pow(SPB(2,1),2)*pow(SPB(5,4),2))/
 (pow(SPA(3,4),3)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
   S(1,5)/S(3,4))*SPB(5,1)))-(pow(SPB(4,1),2)*SPA(2,5))/
(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(3,4)*SPB(5,1))-
 (pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*SPA(2,5)*SPB(4,1)*
 (-(S(2,5)*SPB(4,1))+complex<T>(-2,0)*SPA(2,3)*SPB(2,1)*SPB(4,3)))/
(complex<T>(2,0)*pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPB(5,1))-
 pow(SPB(4,1),3)/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,1)*SPB(5,4)))+
 (complex<T>(0,-1)*pow(SPB(4,1),2))/(complex<T>(2,0)*SPB(2,1)*SPB(5,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpqmQpQmgam_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, Qp, Qm, gam}, sl}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpqmQpQmgam sl");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(3,1),2))/(complex<T>(2,0)*SPB(4,3)*SPB(5,1)*
 SPB(5,2))+complex<T>(0,1)*
((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(2,5),3)*pow(SPB(2,1),2)*pow(SPB(5,3),2))/
(pow(SPA(3,4),3)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(1,5)/S(3,4))*SPB(5,1))+(complex<T>(1,0)*pow(SPB(3,1),2)*SPA(2,5))/
(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(3,4)*SPB(5,1))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*SPA(2,5)*SPB(3,1)*
 (-(S(2,5)*SPB(3,1))+complex<T>(2,0)*SPA(2,4)*SPB(2,1)*SPB(4,3)))/
(complex<T>(2,0)*pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPB(5,1))+
 (complex<T>(1,0)*pow(SPB(3,1),3))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,1)*
 SPB(5,3)))+complex<T>(0,-1)*
(-((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(4,5),2)*SPB(5,1))/
 (complex<T>(2,0)*pow(SPA(3,4),2)*SPB(4,3)*SPB(5,2)))+
 (complex<T>(1,0)*pow(SPB(3,1),2)*SPB(4,2))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*
 SPB(5,2)*SPB(5,4)))+(complex<T>(0,1)*pow(SPB(3,1),2))/
(complex<T>(2,0)*SPB(2,1)*SPB(5,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpQmQpqmgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qm, Qp, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpQmQpqmgap nf");
#endif
 
 return( (complex<T>(0,2)*pow(SPA(2,4),2))/(complex<T>(9,0)*SPA(1,5)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpQpQmqmgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qp, Qm, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpQpQmqmgap nf");
#endif
 
 return( (complex<T>(0,2)*pow(SPA(3,4),2))/(complex<T>(9,0)*SPA(1,5)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmQpQmqpgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmQpQmqpgap nf");
#endif
 
 return( (complex<T>(0,2)*pow(SPA(1,3),2))/(complex<T>(9,0)*SPA(1,5)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmQmQpqpgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmQmQpqpgap nf");
#endif
 
 return( (complex<T>(0,2)*pow(SPA(1,2),2))/(complex<T>(9,0)*SPA(1,5)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpQmQpqmgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qm, Qp, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpQmQpqmgap L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
pow(SPB(5,3),2)*SPA(2,5))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(1,5)*
SPA(3,4))+(complex<T>(1,0)*(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)+
 pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))*pow(SPA(2,3),2)*
pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(1,3))/
 (complex<T>(2,0)*pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPA(1,5)*SPA(3,4)*
SPA(3,5))+(((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),3))/
 (pow(SPA(1,2),4)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
 S(4,5)/S(1,2))*SPA(4,5))+(complex<T>(7,0)*pow(SPA(2,4),2))/
 (complex<T>(9,0)*SPA(1,5)*SPA(2,3)*SPA(4,5))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(2,3),2)*
pow(SPB(5,3),2)*SPA(1,4))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(1,2)*
SPA(1,5)*SPA(3,4)*SPA(4,5))+(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
(-((pow(SPB(5,3),2)*SPA(2,3)*SPA(2,4))/(complex<T>(2,0)*SPA(1,2)))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPA(1,4)*SPA(4,5))/
  (complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(3,4))))/(pow(SPA(3,4),2)*
pow(SPB(4,3),2))+(complex<T>(1,0)*pow(SPB(3,1),2)*SPA(1,4))/
 (complex<T>(2,0)*SPA(1,5)*SPA(4,5)*SPB(2,1)*SPB(4,3))+
(complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPA(1,4)*SPA(4,5))/
 (complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPA(3,4),2)*SPA(1,5)*SPB(2,1)*
SPB(4,3))+(complex<T>(1,0)*SPA(1,4)*SPB(3,1)*SPB(5,1))/
 (complex<T>(2,0)*S(4,5)*SPA(1,5)*SPB(2,1))+(SPA(2,4)*SPB(5,3))/
 (S(2,3)*SPA(1,5))+(complex<T>(1,0)*SPA(1,4)*SPA(2,3)*SPA(2,4)*
SPB(5,3))/(complex<T>(2,0)*S(4,5)*SPA(1,2)*SPA(1,5)*SPA(3,4))-
(pow(SPA(2,4),2)*SPB(5,3))/(complex<T>(2,0)*S(3,4)*SPA(1,2)*SPA(4,5))+
(complex<T>(1,0)*pow(SPA(2,4),2)*SPB(5,3))/(complex<T>(2,0)*pow(SPA(1,2),2)*
SPA(4,5)*SPB(2,1))-(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*SPA(2,4)*
(S(3,5)*SPA(2,4)+complex<T>(-2,0)*SPA(1,2)*SPA(3,4)*SPB(3,1))*
SPB(5,3))/(complex<T>(2,0)*pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPA(4,5))+
(complex<T>(1,0)*SPA(1,2)*SPA(4,5)*SPB(5,1)*SPB(5,3))/
 (complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPA(1,5))-
(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*SPA(1,3)*SPA(2,4)*SPA(4,5)*
SPB(5,1)*SPB(5,3))/(pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPA(1,5)*
SPA(3,4))-(SPA(1,4)*SPA(2,4)*SPB(3,1)*SPB(5,4))/
 (complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPA(1,5))+
(complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*SPB(5,4))/
 (complex<T>(2,0)*S(1,2)*S(3,4)*SPA(1,5)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpQpQmqmgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qp, Qm, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpQpQmqmgap L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(7,0)*pow(SPA(3,4),2))/(complex<T>(9,0)*SPA(1,5)*SPA(2,3)*
SPA(4,5))-(SPA(3,4)*SPB(5,2))/(S(2,3)*SPA(1,5))-
(SPA(1,3)*SPA(4,5)*SPB(5,1)*SPB(5,2))/(complex<T>(2,0)*S(2,3)*
(S(2,3)-S(4,5))*SPA(1,5))+(complex<T>(1,0)*SPA(1,4)*SPA(3,4)*
SPB(2,1)*SPB(5,4))/(complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPA(1,5)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmQpQmqpgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, gap}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmQpQmqpgap L");
#endif
 
 return( complex<T>(0,-1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(1,2),2)*pow(SPA(3,5),2)*pow(SPB(5,2),3))/
(pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(1,5)/S(3,4))*SPA(1,5)))+(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
((complex<T>(1,0)*pow(SPB(5,2),2)*SPA(1,3)*SPA(2,3))/
  (complex<T>(2,0)*SPA(3,4))-(pow(SPA(1,3),2)*pow(SPB(5,1),2)*
   SPA(1,4)*SPA(1,5))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*SPA(4,5))))/
 (pow(SPA(1,2),2)*pow(SPB(2,1),2))+(complex<T>(-7,0)*pow(SPA(1,3),2))/
 (complex<T>(9,0)*SPA(1,5)*SPA(2,3)*SPA(4,5))-
((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)+pow(complex<T>(1,0)-S(1,5)/S(3,4),
 -1))*pow(SPA(1,5),2)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*
SPA(2,4))/(complex<T>(2,0)*pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(1,2)*
SPA(2,5)*SPA(4,5))-(pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*
pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(1,4))/
 (complex<T>(2,0)*pow(SPB(5,1),2)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(4,5))-
(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*pow(SPB(5,2),2)*SPA(3,5))/
 (complex<T>(2,0)*pow(SPB(2,1),2)*SPA(1,2)*SPA(4,5))-
(pow(SPB(4,2),2)*SPA(1,4))/(complex<T>(2,0)*SPA(1,5)*SPA(4,5)*SPB(2,1)*
SPB(4,3))-(pow(SPA(1,3),2)*pow(SPB(5,1),2)*SPA(1,4)*
SPA(1,5))/(complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPA(3,4),2)*SPA(4,5)*
SPB(2,1)*SPB(4,3))-(pow(SPA(1,3),2)*SPB(3,2)*SPB(5,1))/
 (complex<T>(2,0)*S(1,2)*S(3,4)*SPA(4,5))+
(complex<T>(1,0)*SPA(1,3)*SPA(1,4)*SPB(4,2)*SPB(5,1))/
 (complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPA(4,5))+
(complex<T>(1,0)*pow(SPA(1,3),2)*SPB(5,2))/(complex<T>(2,0)*S(1,2)*SPA(1,5)*
SPA(3,4))-(SPA(1,3)*SPB(5,2))/(S(2,3)*SPA(4,5))-
(SPA(1,3)*SPA(1,4)*SPA(2,3)*SPB(5,2))/(complex<T>(2,0)*S(1,5)*SPA(1,2)*
SPA(3,4)*SPA(4,5))-(pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*SPA(1,3)*
(-(S(2,5)*SPA(1,3))+complex<T>(2,0)*SPA(1,2)*SPA(3,4)*SPB(4,2))*
SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(1,5))-
(pow(SPA(1,3),2)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*SPA(1,5)*
SPB(4,3))-(SPA(1,4)*SPB(4,2)*SPB(5,4))/
 (complex<T>(2,0)*S(1,5)*SPA(4,5)*SPB(4,3))+
(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*SPA(1,3)*SPA(1,5)*SPA(2,4)*
SPB(5,2)*SPB(5,4))/(pow(SPA(3,4),2)*pow(SPB(4,3),2)*SPA(1,2)*
SPA(4,5))-(SPA(1,5)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
 (complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPA(4,5)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmQmQpqpgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, gap}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmQmQpqpgap L");
#endif
 
 return( complex<T>(0,-1)*((complex<T>(-7,0)*pow(SPA(1,2),2))/(complex<T>(9,0)*SPA(1,5)*SPA(2,3)*
SPA(4,5))-(SPA(1,2)*SPA(1,4)*SPB(4,3)*SPB(5,1))/
 (complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPA(4,5))+
(SPA(1,2)*SPB(5,3))/(S(2,3)*SPA(4,5))+
(complex<T>(1,0)*SPA(1,5)*SPA(2,4)*SPB(5,3)*SPB(5,4))/
 (complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPA(4,5)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpqmQpQmgap_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, Qp, Qm, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpqmQpQmgap sl");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(2,4),2))/(complex<T>(2,0)*SPA(1,5)*SPA(2,5)*
 SPA(3,4))+complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
  pow(SPB(5,3),2)*SPA(2,5))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(1,5)*
  SPA(3,4)))-(pow(SPA(2,4),2)*SPA(1,4))/(complex<T>(2,0)*SPA(1,2)*
 SPA(1,5)*SPA(3,4)*SPA(4,5)))+(complex<T>(0,-1)*pow(SPA(2,4),2))/
(complex<T>(2,0)*SPA(1,2)*SPA(3,5)*SPA(4,5))+
 complex<T>(0,1)*((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(1,2),2)*pow(SPA(4,5),2)*pow(SPB(5,1),3))/
(pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(2,5)/S(3,4))*SPA(2,5))-(pow(SPA(2,4),2)*SPA(2,3))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,5)*SPA(3,4)*SPA(3,5))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1)*SPA(2,4)*
 (-(S(1,5)*SPA(2,4))+complex<T>(2,0)*SPA(1,2)*SPA(3,4)*SPB(3,1))*
 SPB(5,1))/(complex<T>(2,0)*pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(2,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(5,1))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 SPA(2,5)*SPB(4,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qpqmQmQpgap_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, Qm, Qp, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qpqmQmQpgap sl");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(2,3),2))/(complex<T>(2,0)*SPA(1,5)*SPA(2,5)*
 SPA(3,4))+complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
 pow(SPB(5,4),2)*SPA(2,5))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(1,5)*
 SPA(3,4))+(complex<T>(1,0)*pow(SPA(2,3),2)*SPA(1,3))/
(complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(3,5)))+
 (complex<T>(0,1)*pow(SPA(2,3),2))/(complex<T>(2,0)*SPA(1,2)*SPA(3,5)*SPA(4,5))+
 complex<T>(0,-1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
   (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1))/complex<T>(2,0))*
  pow(SPA(1,2),2)*pow(SPA(3,5),2)*pow(SPB(5,1),3))/
 (pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
   S(2,5)/S(3,4))*SPA(2,5)))+(complex<T>(1,0)*pow(SPA(2,3),2)*
 SPA(2,4))/(complex<T>(2,0)*SPA(1,2)*SPA(2,5)*SPA(3,4)*SPA(4,5))-
 (pow(complex<T>(1,0)-S(2,5)/S(3,4),-1)*SPA(2,3)*(-(S(1,5)*SPA(2,3))+
  complex<T>(-2,0)*SPA(1,2)*SPA(3,4)*SPB(4,1))*SPB(5,1))/
(complex<T>(2,0)*pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(2,5))-
 (pow(SPA(2,3),2)*SPB(5,1))/(complex<T>(2,0)*pow(SPA(3,4),2)*SPA(2,5)*
 SPB(4,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmqpQpQmgap_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmqpQpQmgap sl");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,4),2))/(complex<T>(2,0)*SPA(1,5)*SPA(2,5)*
 SPA(3,4))+complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
 pow(SPB(5,3),2)*SPA(1,5))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(2,5)*
 SPA(3,4))-(pow(SPA(1,4),2)*SPA(2,3))/(complex<T>(2,0)*SPA(1,2)*
 SPA(2,5)*SPA(3,4)*SPA(3,5)))+(complex<T>(0,-1)*pow(SPA(1,4),2))/
(complex<T>(2,0)*SPA(1,2)*SPA(3,5)*SPA(4,5))+
 complex<T>(0,-1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
   (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
  pow(SPA(1,2),2)*pow(SPA(4,5),2)*pow(SPB(5,2),3))/
 (pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
   S(1,5)/S(3,4))*SPA(1,5)))-pow(SPA(1,4),3)/
(complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(4,5))-
 (pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*SPA(1,4)*(-(S(2,5)*SPA(1,4))+
  complex<T>(-2,0)*SPA(1,2)*SPA(3,4)*SPB(3,2))*SPB(5,2))/
(complex<T>(2,0)*pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(1,5))-
 (pow(SPA(1,4),2)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*SPA(1,5)*
 SPB(4,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1y_qmqpQmQpgap_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, gap}, sl}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1y :  qmqpQmQpgap sl");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,3),2))/(complex<T>(2,0)*SPA(1,5)*SPA(2,5)*
 SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,3),2))/(complex<T>(2,0)*SPA(1,2)*SPA(3,5)*
 SPA(4,5))+complex<T>(0,-1)*((complex<T>(1,0)*pow(SPA(1,3),2)*SPA(2,4))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,5)*SPA(3,4)*SPA(4,5))-
 (pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPB(5,4),2)*SPA(1,5))/
(complex<T>(2,0)*S(3,4)*SPA(2,5)*SPB(4,3)))+
 complex<T>(0,1)*((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(1,2),2)*pow(SPA(3,5),2)*pow(SPB(5,2),3))/
(pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(1,5)/S(3,4))*SPA(1,5))+(complex<T>(1,0)*pow(SPA(1,3),3))/
(complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(3,5))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*SPA(1,3)*
 (-(S(2,5)*SPA(1,3))+complex<T>(2,0)*SPA(1,2)*SPA(3,4)*SPB(4,2))*
 SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(1,5))+
 (complex<T>(1,0)*pow(SPA(1,3),2)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 SPA(1,5)*SPB(4,3)))
        ); 
  
} 
  
  
 
 // *************** table of switch values ************* 
 
#define _R_qmQpQmqpgam_nf R2q2Q1y_607_nf
#define _R_qmQmQpqpgam_nf R2q2Q1y_637_nf
#define _R_qpQmQpqmgam_nf R2q2Q1y_422_nf
#define _R_qpQpQmqmgam_nf R2q2Q1y_392_nf
#define _R_qmQpQmqpgam_L R2q2Q1y_607_L
#define _R_qmQmQpqpgam_L R2q2Q1y_637_L
#define _R_qpQmQpqmgam_L R2q2Q1y_422_L
#define _R_qpQpQmqmgam_L R2q2Q1y_392_L
#define _R_qmqpQmQpgam_sl R2q2Q1y_1237_sl
#define _R_qmqpQpQmgam_sl R2q2Q1y_1057_sl
#define _R_qpqmQmQpgam_sl R2q2Q1y_1232_sl
#define _R_qpqmQpQmgam_sl R2q2Q1y_1052_sl
#define _R_qpQmQpqmgap_nf R2q2Q1y_4310_nf
#define _R_qpQpQmqmgap_nf R2q2Q1y_4280_nf
#define _R_qmQpQmqpgap_nf R2q2Q1y_4495_nf
#define _R_qmQmQpqpgap_nf R2q2Q1y_4525_nf
#define _R_qpQmQpqmgap_L R2q2Q1y_4310_L
#define _R_qpQpQmqmgap_L R2q2Q1y_4280_L
#define _R_qmQpQmqpgap_L R2q2Q1y_4495_L
#define _R_qmQmQpqpgap_L R2q2Q1y_4525_L
#define _R_qpqmQpQmgap_sl R2q2Q1y_4940_sl
#define _R_qpqmQmQpgap_sl R2q2Q1y_5120_sl
#define _R_qmqpQpQmgap_sl R2q2Q1y_4945_sl
#define _R_qmqpQmQpgap_sl R2q2Q1y_5125_sl
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmQpQmqpgam_nf case 607 : \
          return &R2q2Q1y_607_nf
#define _CASE_qmQmQpqpgam_nf case 637 : \
          return &R2q2Q1y_637_nf
#define _CASE_qpQmQpqmgam_nf case 422 : \
          return &R2q2Q1y_422_nf
#define _CASE_qpQpQmqmgam_nf case 392 : \
          return &R2q2Q1y_392_nf
#define _CASE_qmQpQmqpgam_L case 607 : \
          return &R2q2Q1y_607_L
#define _CASE_qmQmQpqpgam_L case 637 : \
          return &R2q2Q1y_637_L
#define _CASE_qpQmQpqmgam_L case 422 : \
          return &R2q2Q1y_422_L
#define _CASE_qpQpQmqmgam_L case 392 : \
          return &R2q2Q1y_392_L
#define _CASE_qmqpQmQpgam_sl case 1237 : \
          return &R2q2Q1y_1237_sl
#define _CASE_qmqpQpQmgam_sl case 1057 : \
          return &R2q2Q1y_1057_sl
#define _CASE_qpqmQmQpgam_sl case 1232 : \
          return &R2q2Q1y_1232_sl
#define _CASE_qpqmQpQmgam_sl case 1052 : \
          return &R2q2Q1y_1052_sl
#define _CASE_qpQmQpqmgap_nf case 4310 : \
          return &R2q2Q1y_4310_nf
#define _CASE_qpQpQmqmgap_nf case 4280 : \
          return &R2q2Q1y_4280_nf
#define _CASE_qmQpQmqpgap_nf case 4495 : \
          return &R2q2Q1y_4495_nf
#define _CASE_qmQmQpqpgap_nf case 4525 : \
          return &R2q2Q1y_4525_nf
#define _CASE_qpQmQpqmgap_L case 4310 : \
          return &R2q2Q1y_4310_L
#define _CASE_qpQpQmqmgap_L case 4280 : \
          return &R2q2Q1y_4280_L
#define _CASE_qmQpQmqpgap_L case 4495 : \
          return &R2q2Q1y_4495_L
#define _CASE_qmQmQpqpgap_L case 4525 : \
          return &R2q2Q1y_4525_L
#define _CASE_qpqmQpQmgap_sl case 4940 : \
          return &R2q2Q1y_4940_sl
#define _CASE_qpqmQmQpgap_sl case 5120 : \
          return &R2q2Q1y_5120_sl
#define _CASE_qmqpQpQmgap_sl case 4945 : \
          return &R2q2Q1y_4945_sl
#define _CASE_qmqpQmQpgap_sl case 5125 : \
          return &R2q2Q1y_5125_sl
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_qmQpQmqpgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmQpQmqpgam_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpqpgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmQmQpqpgam_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpQmQpqmgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpQmQpqmgam_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpQpQmqmgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpQpQmqmgam_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmQpQmqpgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmQpQmqpgam_L(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpqpgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmQmQpqpgam_L(ep,mpc);}
 
template <class T> complex<T> _R_qpQmQpqmgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpQmQpqmgam_L(ep,mpc);}
 
template <class T> complex<T> _R_qpQpQmqmgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpQpQmqmgam_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQpgam_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmqpQmQpgam_sl(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQmgam_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmqpQpQmgam_sl(ep,mpc);}
 
template <class T> complex<T> _R_qpqmQmQpgam_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpqmQmQpgam_sl(ep,mpc);}
 
template <class T> complex<T> _R_qpqmQpQmgam_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpqmQpQmgam_sl(ep,mpc);}
 
template <class T> complex<T> _R_qpQmQpqmgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpQmQpqmgap_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpQpQmqmgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpQpQmqmgap_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmQpQmqpgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmQpQmqpgap_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpqpgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmQmQpqpgap_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpQmQpqmgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpQmQpqmgap_L(ep,mpc);}
 
template <class T> complex<T> _R_qpQpQmqmgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpQpQmqmgap_L(ep,mpc);}
 
template <class T> complex<T> _R_qmQpQmqpgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmQpQmqpgap_L(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpqpgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmQmQpqpgap_L(ep,mpc);}
 
template <class T> complex<T> _R_qpqmQpQmgap_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpqmQpQmgap_sl(ep,mpc);}
 
template <class T> complex<T> _R_qpqmQmQpgap_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qpqmQmQpgap_sl(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQmgap_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmqpQpQmgap_sl(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQpgap_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1y_qmqpQmQpgap_sl(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q2Q1y_L_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmQpQmqpgam_L;
       _CASE_qmQmQpqpgam_L;
       _CASE_qpQmQpqmgam_L;
       _CASE_qpQpQmqmgam_L;
       _CASE_qpQmQpqmgap_L;
       _CASE_qpQpQmqmgap_L;
       _CASE_qmQpQmqpgap_L;
       _CASE_qmQmQpqpgap_L;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q1y_nf_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmQpQmqpgam_nf;
       _CASE_qmQmQpqpgam_nf;
       _CASE_qpQmQpqmgam_nf;
       _CASE_qpQpQmqmgam_nf;
       _CASE_qpQmQpqmgap_nf;
       _CASE_qpQpQmqmgap_nf;
       _CASE_qmQpQmqpgap_nf;
       _CASE_qmQmQpqpgap_nf;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q1y_sl_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmqpQmQpgam_sl;
       _CASE_qmqpQpQmgam_sl;
       _CASE_qpqmQmQpgam_sl;
       _CASE_qpqmQpQmgam_sl;
       _CASE_qpqmQpQmgap_sl;
       _CASE_qpqmQmQpgap_sl;
       _CASE_qmqpQpQmgap_sl;
       _CASE_qmqpQmQpgap_sl;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q2Q1y_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q1y_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1y_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q1y_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2Q1y_nf_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q1y_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1y_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q1y_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2Q1y_sl_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q1y_sl_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1y_sl_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q1y_sl_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif




}
