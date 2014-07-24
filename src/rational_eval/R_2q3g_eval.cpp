/*
* R_2q3g.cpp
*
* Created on 12/12, 2008
*      Author: Zvi's script
*/
 
#include "R_2q3g_eval.h"
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

 
 
template <class T> complex<T> R2q3g_qmqpppp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpppp L");
#endif
 
 return( (complex<T>(0,-1)*(SPA(1,2)*SPA(1,3)*SPB(3,2)+SPA(1,4)*SPA(1,5)*
 SPB(5,4)))/(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*((pow(SPA(1,4),2)*SPA(1,3)*SPB(4,3))/(pow(SPA(3,4),2)*
  SPA(1,2)*SPA(1,5)*SPA(4,5))-(SPB(3,2)*SPB(5,2))/
 (SPA(3,4)*SPA(4,5)*SPB(2,1))+(SPA(1,4)*SPA(1,5)*SPA(2,4)*
  SPB(5,4))/(pow(SPA(4,5),2)*SPA(1,2)*SPA(2,3)*SPA(3,4))))/
complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpmpp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpmpp L");
#endif
 
 return( complex<T>(0,1)*((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
pow(SPA(3,5),2)*pow(SPB(5,2),3)*SPA(1,2)*SPA(2,3))/
 (pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
 S(1,5)/S(3,4))*SPA(4,5))-(pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*
pow(SPB(4,2),2)*SPA(1,3)*SPA(1,4)*SPA(2,3))/
 (complex<T>(2,0)*pow(SPA(1,5),3)*pow(SPB(5,1),2)*SPA(4,5))-
(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPB(5,4),2)*SPA(1,3)*
SPA(1,5))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(1,2)*SPA(3,4)*
SPA(4,5))+(((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))*
pow(SPB(4,2),3)*SPA(1,2)*SPA(1,4)*SPA(2,3)*SPA(3,4))/
 (pow(SPA(1,5),4)*pow(SPB(5,1),3)*(complex<T>(1,0)-S(2,3)/S(1,5)-
 S(3,4)/S(1,5))*SPA(4,5))-(pow(SPA(1,3),2)*SPA(2,3)*SPB(2,1)*
SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*S(1,5)*SPA(4,5)*SPB(4,3))-
(pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*
((pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(1,3)*SPA(1,5))/
  (SPA(1,2)*SPA(3,4)*SPA(4,5))+(complex<T>(1,0)*pow(SPA(1,3),2)*
   SPA(2,3)*SPB(2,1)*SPB(5,2))/(complex<T>(2,0)*SPA(3,4)*SPA(4,5))))/
 (pow(SPA(1,5),2)*pow(SPB(5,1),2))-
(pow(SPB(4,2),2)*SPA(1,4)*SPB(5,4))/(complex<T>(2,0)*S(1,5)*SPA(4,5)*
SPB(3,2)*SPB(4,3))-(SPA(1,3)*SPA(2,3)*SPB(5,2)*SPB(5,4))/
 (S(3,4)*SPA(1,2)*SPA(4,5)*SPB(5,1))+
(complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPB(5,4),2)*
SPA(1,3)*SPA(1,5))/(pow(SPB(4,3),2)*SPA(1,2)*SPA(3,4)*
SPA(4,5)))-(pow(SPB(5,4),2)*(-(S(1,2)/S(3,4))+
S(3,4)/S(1,2))*SPA(1,5)*SPA(3,5)*SPB(5,2))/
  (pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
   pow(SPB(4,3),3)*SPA(4,5))+(pow(SPB(4,2),2)*SPB(5,2))/
  (SPA(4,5)*SPB(2,1)*SPB(3,2)*SPB(4,3))-
 (SPA(1,3)*SPB(4,2)*SPB(5,4))/(S(1,2)*SPA(4,5)*SPB(4,3))))/
 complex<T>(3,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqppmp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqppmp L");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,4),3)*SPA(2,4))/(SPA(1,2)*SPA(1,5)*
 SPA(2,3)*SPA(3,4)*SPA(4,5))+
 complex<T>(0,1)*(-((pow(SPB(5,3),3)*(-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*
  SPA(1,3)*SPA(1,4))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*
  pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPA(1,2)))-
 (pow(SPB(5,3),3)*(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPA(1,4)*
 SPA(1,5))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPA(1,2))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(2,4),2)*pow(SPB(5,2),3)*SPA(1,2)*SPA(4,5))/
(pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(1,5)/S(3,4))*SPA(2,3))+(complex<T>(2,0)*pow(SPB(5,3),4)*
 (-(S(1,2)/(complex<T>(6,0)*S(3,4)))+(complex<T>(1,0)*(S(1,2)/S(3,4)-
 S(3,4)/S(1,2)))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPA(1,3)*SPA(1,5)*
 SPA(3,4)*SPA(4,5))/(pow(SPA(1,2),5)*pow(SPB(2,1),4)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2)))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPB(5,3),3)*SPA(1,3)*SPA(3,4)*SPA(4,5))/
(pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
  S(4,5)/S(1,2))*SPA(2,3))+
 (complex<T>(1,0)*((pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*pow(SPA(2,4),2)*
pow(SPB(5,2),2)*SPA(1,4))/(pow(SPA(1,5),2)*pow(SPB(5,1),2)*
SPA(2,3)*SPA(3,4))+(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
pow(SPB(5,3),2)*SPA(1,4))/(pow(SPB(4,3),2)*SPA(2,3)*
SPA(3,4))+(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
pow(SPB(5,2),2)*SPA(1,4)*SPA(2,4)*SPA(4,5))/
   (pow(SPB(2,1),2)*SPA(1,2)*SPA(1,5)*SPA(2,3)*SPA(3,4))+
  (pow(SPA(1,4),2)*SPA(2,4)*SPB(3,2))/(S(1,2)*SPA(1,5)*
SPA(2,3)*SPA(4,5))+(pow(SPA(2,4),2)*pow(SPB(5,2),2)*
SPA(1,4))/(pow(SPA(3,4),2)*S(1,5)*SPA(2,3)*SPB(4,3))+
  (pow(SPB(5,3),2)*SPA(1,3)*SPA(1,4)*SPB(5,1))/
   (S(1,2)*S(3,4)*SPA(2,3)*SPB(5,4))-
  (pow(SPB(5,2),2)*SPA(2,4)*SPB(5,3))/(S(3,4)*SPA(2,3)*
SPB(2,1)*SPB(5,4))))/complex<T>(2,0)-
 (pow(SPB(5,3),4)*SPA(1,3)*SPA(1,5))/(complex<T>(3,0)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(4,3)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpppm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpppm L");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,5),2)*SPA(2,5))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*
 SPA(3,4)*SPA(4,5))+complex<T>(0,1)*
((complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPB(4,3),2)*
 SPA(1,3)*SPA(3,5))/(pow(SPB(5,4),2)*SPA(2,3)*SPA(3,4)*
 SPA(4,5)))+(SPB(4,2)*SPB(4,3))/(SPA(3,4)*SPB(5,1)*
SPB(5,4))))/complex<T>(2,0)+
 (complex<T>(1,0)*((pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(4,3),3)*
(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),4)*
pow(SPB(2,1),3)*SPA(3,4))+
  (complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*pow(SPB(4,3),2)*
SPA(1,3)*SPA(1,5)*SPA(4,5))/(pow(SPA(1,2),3)*
pow(SPB(2,1),2)*SPA(3,4))+(complex<T>(2,0)*pow(SPA(1,5),2)*
SPB(4,3))/(pow(SPA(1,2),2)*SPA(3,4)*SPB(2,1))+
  (SPA(1,5)*SPB(3,1)*SPB(4,2))/(S(1,2)*SPA(3,4)*SPB(5,1))+
  (S(1,4)*SPB(4,2)*SPB(4,3))/(S(1,2)*SPA(3,4)*SPB(5,1)*
SPB(5,4))))/complex<T>(3,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpqppp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpqppp L");
#endif
 
 return( (complex<T>(0,-1)*SPA(1,4)*SPA(1,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(4,5),2)*
 SPA(1,2)*SPA(2,3))+
 (complex<T>(0,1)*(-((pow(SPA(1,3),2)*SPB(3,2))/(SPA(1,5)*SPA(2,3)*SPA(3,4)*
   SPA(4,5)))-(SPA(1,3)*SPA(1,4)*SPB(5,4))/
 (SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmqppp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmqppp L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*
pow(SPA(2,4),2)*pow(SPB(5,4),2))/(complex<T>(2,0)*pow(SPB(5,1),2)*
SPA(1,5)*SPA(3,4)*SPA(4,5))+(complex<T>(1,0)*SPB(4,3)*SPB(5,3))/
 (complex<T>(3,0)*SPA(4,5)*SPB(2,1)*SPB(3,2))+
(complex<T>(1,0)*SPB(3,1)*SPB(5,3)*SPB(5,4))/(complex<T>(2,0)*SPA(4,5)*SPB(2,1)*
SPB(3,2)*SPB(5,1)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpqpmp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpqpmp L");
#endif
 
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*pow(SPB(5,2),2)*
 SPA(1,4)*SPA(2,4))/(complex<T>(2,0)*pow(SPA(3,4),2)*pow(SPB(4,3),2)*
 SPA(2,3)))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
pow(SPB(5,3),2)*SPA(1,3)*SPA(1,4))/(complex<T>(2,0)*pow(SPB(4,3),2)*
SPA(1,2)*SPA(2,3)*SPA(3,4))+
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
pow(SPB(5,2),3)*SPA(1,2)*SPA(2,4)*SPA(4,5))/
 (pow(SPA(3,4),3)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
 S(1,5)/S(3,4))*SPA(2,3))+
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
pow(SPA(1,3),2)*pow(SPB(5,3),3)*SPA(3,4)*SPA(4,5))/
 (pow(SPA(1,2),4)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
 S(4,5)/S(1,2))*SPA(2,3))+(complex<T>(1,0)*SPA(1,4)*SPA(2,4)*
SPB(5,2))/(complex<T>(2,0)*S(3,4)*SPA(2,3)*SPA(2,5))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*SPA(1,2)*SPA(1,4)*
SPA(4,5)*SPB(5,1)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*
pow(SPB(4,3),2)*SPA(2,3)*SPA(2,5))-
(pow(SPA(1,3),2)*pow(SPB(5,3),3))/(complex<T>(2,0)*pow(SPA(1,2),2)*
SPA(2,3)*SPB(2,1)*SPB(4,3)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpqppm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpqppm L");
#endif
 
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(1,3),2)*
 pow(SPB(4,3),2)*SPA(3,5))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(1,2)*
 SPA(2,3)*SPA(3,4)*SPA(4,5)))-pow(SPB(4,2),2)/
 (complex<T>(2,0)*SPA(3,4)*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmppqpp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmppqpp L");
#endif
 
 return( (complex<T>(0,1)*(-((SPA(1,3)*SPA(1,4)*SPB(3,2))/(SPA(1,5)*SPA(2,3)*
  SPA(3,4)*SPA(4,5)))-(pow(SPA(1,4),2)*SPB(5,4))/
(SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmpqpp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmpqpp L");
#endif
 
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-S(1,5)/S(2,3),-1)*pow(SPA(1,4),2)*
 pow(SPB(4,3),2)*SPA(2,4))/(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(1,5)*
 SPA(2,3)*SPA(3,4)*SPA(4,5)))-pow(SPB(5,3),2)/
 (complex<T>(2,0)*SPA(3,4)*SPB(2,1)*SPB(3,2)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpmqpp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpmqpp L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*
pow(SPB(4,2),2)*SPA(1,3)*SPA(1,4))/(complex<T>(2,0)*pow(SPB(4,3),2)*
SPA(1,5)*SPA(3,4)*SPA(4,5))+
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))*
pow(SPA(1,4),2)*pow(SPB(4,2),3)*SPA(2,3)*SPA(3,4))/
 (pow(SPA(1,5),4)*pow(SPB(5,1),3)*(complex<T>(1,0)-S(2,3)/S(1,5)-
 S(3,4)/S(1,5))*SPA(4,5))-(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
pow(SPB(5,2),2)*SPA(1,3)*SPA(3,5))/(complex<T>(2,0)*pow(SPA(3,4),2)*
pow(SPB(4,3),2)*SPA(4,5))+
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
pow(SPB(5,2),3)*SPA(1,5)*SPA(2,3)*SPA(3,5))/
 (pow(SPA(3,4),3)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
 S(1,5)/S(3,4))*SPA(4,5))-(pow(SPA(1,4),2)*pow(SPB(4,2),3))/
 (complex<T>(2,0)*pow(SPA(1,5),2)*SPA(4,5)*SPB(3,2)*SPB(4,3)*SPB(5,1))+
(complex<T>(1,0)*SPA(1,3)*SPA(3,5)*SPB(5,2))/(complex<T>(2,0)*S(3,4)*SPA(2,5)*
SPA(4,5))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*SPA(1,3)*
SPA(1,5)*SPA(2,3)*SPB(2,1)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*
pow(SPB(4,3),2)*SPA(2,5)*SPA(4,5)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmppqpm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmppqpm L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*
 pow(SPA(3,5),2)*pow(SPB(3,2),2))/(complex<T>(2,0)*pow(SPB(2,1),2)*
 SPA(1,2)*SPA(2,3)*SPA(3,4))+(complex<T>(1,0)*SPB(3,2)*SPB(4,1)*
 SPB(4,2))/(complex<T>(2,0)*SPA(2,3)*SPB(2,1)*SPB(5,1)*SPB(5,4))+
 (complex<T>(1,0)*SPB(4,2)*SPB(4,3))/(complex<T>(3,0)*SPA(2,3)*SPB(5,1)*
 SPB(5,4)))+(complex<T>(0,-1)*SPB(4,2)*SPB(4,3))/
(complex<T>(3,0)*SPA(2,3)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpppqp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, p, qp}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpppqp L");
#endif
 
 return( (complex<T>(0,1)*(-(SPA(1,2)*SPA(1,3)*SPB(3,2))-
 SPA(1,4)*SPA(1,5)*SPB(5,4)))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*
SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmppqp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, p, qp}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmppqp L");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,2),2)*SPA(2,5))/(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*
 SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*((pow(SPA(1,4),2)*pow(SPA(2,3),2)*pow(SPB(4,3),3)*
  (S(1,5)/S(2,3)-S(2,3)/S(1,5)))/
 (pow(complex<T>(1,0)-S(2,3)/S(1,5),3)*pow(SPA(1,5),4)*
  pow(SPB(5,1),3)*SPA(3,4))+
(complex<T>(3,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*pow(SPB(4,3),2)*
  SPA(1,2)*SPA(1,4)*SPA(2,3))/(pow(SPA(1,5),3)*pow(SPB(5,1),2)*
  SPA(3,4))+(complex<T>(2,0)*pow(SPA(1,2),2)*SPB(4,3))/
 (pow(SPA(1,5),2)*SPA(3,4)*SPB(5,1))+
(SPA(1,2)*SPB(4,1)*SPB(5,3))/(S(1,5)*SPA(3,4)*SPB(2,1))+
(S(1,3)*SPB(4,3)*SPB(5,3))/(S(1,5)*SPA(3,4)*SPB(2,1)*
  SPB(3,2))))/complex<T>(3,0)+
 complex<T>(0,1)*((complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(1,5)/S(2,3),-1)*
 pow(SPB(4,3),2)*SPA(1,4)*SPA(2,4))/(pow(SPB(3,2),2)*
 SPA(2,3)*SPA(3,4)*SPA(4,5)))+(SPB(4,3)*SPB(5,3))/
   (SPA(3,4)*SPB(2,1)*SPB(3,2))))/complex<T>(2,0)+
 (complex<T>(1,0)*((pow(SPA(1,4),2)*pow(SPA(2,3),2)*pow(SPB(4,3),3)*
(S(1,5)/S(2,3)-S(2,3)/S(1,5)))/
   (pow(complex<T>(1,0)-S(2,3)/S(1,5),3)*pow(SPA(1,5),4)*
pow(SPB(5,1),3)*SPA(3,4))+
  (complex<T>(3,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*pow(SPB(4,3),2)*
SPA(1,2)*SPA(1,4)*SPA(2,3))/(pow(SPA(1,5),3)*
pow(SPB(5,1),2)*SPA(3,4))+(complex<T>(2,0)*pow(SPA(1,2),2)*
SPB(4,3))/(pow(SPA(1,5),2)*SPA(3,4)*SPB(5,1))+
  (SPA(1,2)*SPB(4,1)*SPB(5,3))/(S(1,5)*SPA(3,4)*SPB(2,1))+
  (S(1,3)*SPB(4,3)*SPB(5,3))/(S(1,5)*SPA(3,4)*SPB(2,1)*
SPB(3,2))))/complex<T>(3,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpmpqp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, p, qp}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpmpqp L");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,3),3)*SPA(3,5))/(SPA(1,2)*SPA(1,5)*
 SPA(2,3)*SPA(3,4)*SPA(4,5))+
 complex<T>(0,-1)*(-((pow(SPB(4,2),3)*(-(S(1,5)/S(3,4))+S(3,4)/S(1,5))*
  SPA(1,2)*SPA(1,3))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),3)*
  pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPA(1,5)))-
 (pow(SPB(4,2),3)*(-(S(1,5)/S(2,3))+S(2,3)/S(1,5))*SPA(1,3)*
 SPA(1,4))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),3)*
 pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPA(1,5))+
 (complex<T>(2,0)*pow(SPB(4,2),4)*(-(S(1,5)/(complex<T>(6,0)*S(2,3)))+
  (complex<T>(1,0)*(S(1,5)/S(2,3)-S(2,3)/S(1,5)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),3))-
  S(1,5)/(complex<T>(6,0)*S(3,4))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5))+
  (complex<T>(1,0)*(S(1,5)/S(3,4)-S(3,4)/S(1,5)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),3)))*SPA(1,2)*SPA(1,4)*
 SPA(2,3)*SPA(3,4))/(pow(SPA(1,5),5)*pow(SPB(5,1),4)*
 (complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5)))-
 (pow(SPB(4,2),4)*SPA(1,2)*SPA(1,4))/(complex<T>(3,0)*pow(SPA(1,5),3)*
 pow(SPB(5,1),2)*SPB(3,2)*SPB(4,3)))+
 complex<T>(0,1)*(-((pow(SPB(4,2),3)*(-(S(1,5)/S(3,4))+S(3,4)/S(1,5))*
  SPA(1,2)*SPA(1,3))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),3)*
  pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPA(1,5)))-
 (pow(SPB(4,2),3)*(-(S(1,5)/S(2,3))+S(2,3)/S(1,5))*SPA(1,3)*
 SPA(1,4))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),3)*
 pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPA(1,5))+
 (complex<T>(2,0)*pow(SPB(4,2),4)*(-(S(1,5)/(complex<T>(6,0)*S(2,3)))+
  (complex<T>(1,0)*(S(1,5)/S(2,3)-S(2,3)/S(1,5)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),3))-
  S(1,5)/(complex<T>(6,0)*S(3,4))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5))+
  (complex<T>(1,0)*(S(1,5)/S(3,4)-S(3,4)/S(1,5)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),3)))*SPA(1,2)*SPA(1,4)*
 SPA(2,3)*SPA(3,4))/(pow(SPA(1,5),5)*pow(SPB(5,1),4)*
 (complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5)))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(3,5),2)*pow(SPB(5,2),3)*SPA(1,5)*SPA(2,3))/
(pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(1,5)/S(3,4))*SPA(4,5))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))*
 pow(SPB(4,2),3)*SPA(1,4)*SPA(2,3)*SPA(3,4))/
(pow(SPA(1,5),3)*pow(SPB(5,1),3)*(complex<T>(1,0)-S(2,3)/S(1,5)-
  S(3,4)/S(1,5))*SPA(4,5))-(pow(SPB(4,2),4)*SPA(1,2)*
 SPA(1,4))/(complex<T>(3,0)*pow(SPA(1,5),3)*pow(SPB(5,1),2)*SPB(3,2)*
 SPB(4,3))+(complex<T>(1,0)*((pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*
pow(SPB(4,2),2)*SPA(1,3))/(pow(SPB(4,3),2)*SPA(3,4)*
SPA(4,5))+(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPA(1,3))/
   (pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPA(3,4)*SPA(4,5))+
  (pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPA(1,3))/(S(1,2)*S(3,4)*
SPA(3,4)*SPA(4,5))+(pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*
pow(SPB(5,2),2)*SPA(1,3)*SPA(2,3)*SPA(3,5))/
   (pow(SPB(5,1),2)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(4,5))+
  (pow(SPB(4,2),2)*SPA(1,3)*SPA(1,4)*SPB(2,1))/
   (S(1,5)*S(3,4)*SPA(4,5)*SPB(3,2))-
  (pow(SPB(5,2),2)*SPA(3,5)*SPB(4,2))/(S(3,4)*SPA(4,5)*
SPB(3,2)*SPB(5,1))+(pow(SPA(1,3),2)*SPA(3,5)*SPB(5,4))/
   (S(1,5)*SPA(1,2)*SPA(2,3)*SPA(4,5))))/complex<T>(2,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmppmqp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, m, qp}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmppmqp L");
#endif
 
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*pow(SPB(3,2),2)*
  SPA(1,2)*SPA(1,4))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(1,5)*SPA(2,3)*
  SPA(3,4)))-(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*pow(SPB(5,3),2)*
 SPA(1,3)*SPA(1,4)*SPA(4,5))/(complex<T>(2,0)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPA(2,3))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(2,4),2)*pow(SPB(5,2),3)*SPA(1,5)*SPA(4,5))/
(pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(1,5)/S(3,4))*SPA(2,3))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPB(5,3),3)*SPA(1,3)*SPA(1,5)*SPA(3,4)*SPA(4,5))/
(pow(SPA(1,2),4)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
  S(4,5)/S(1,2))*SPA(2,3))-(SPA(1,4)*SPA(4,5)*SPB(3,2)*
 SPB(5,2))/(S(3,4)*SPA(1,5)*SPA(2,3)*SPB(2,1))-
 (pow(SPA(1,4),2)*SPA(4,5)*SPB(5,1)*SPB(5,2))/
(complex<T>(2,0)*pow(SPA(3,4),2)*S(1,2)*SPA(2,3)*SPB(4,3))-
 (pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
 ((pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(1,2)*SPA(1,4))/
   (SPA(1,5)*SPA(2,3)*SPA(3,4))+(complex<T>(1,0)*pow(SPA(1,4),2)*
SPA(4,5)*SPB(5,1)*SPB(5,2))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4))))/
(pow(SPA(1,2),2)*pow(SPB(2,1),2))+
 (complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*pow(SPB(3,2),2)*
 SPA(1,2)*SPA(1,4))/(pow(SPB(4,3),2)*SPA(1,5)*SPA(2,3)*
 SPA(3,4)))-(pow(SPB(3,2),2)*(-(S(1,5)/S(3,4))+
 S(3,4)/S(1,5))*SPA(1,2)*SPA(2,4)*SPB(5,2))/
   (pow(complex<T>(1,0)-S(1,5)/S(3,4),3)*pow(SPA(3,4),2)*
pow(SPB(4,3),3)*SPA(2,3))-(SPA(1,4)*SPB(3,2)*SPB(5,3))/
   (S(1,5)*SPA(2,3)*SPB(4,3))+(pow(SPB(5,3),2)*SPB(5,2))/
   (SPA(2,3)*SPB(4,3)*SPB(5,1)*SPB(5,4))))/complex<T>(3,0)-
 (pow(SPB(5,3),2)*SPA(1,3)*SPB(3,2))/(complex<T>(2,0)*S(1,2)*SPA(2,3)*
 SPB(4,3)*SPB(5,4)))+
 (complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*pow(SPB(3,2),2)*
   SPA(1,2)*SPA(1,4))/(pow(SPB(4,3),2)*SPA(1,5)*SPA(2,3)*
   SPA(3,4)))-(pow(SPB(3,2),2)*(-(S(1,5)/S(3,4))+
   S(3,4)/S(1,5))*SPA(1,2)*SPA(2,4)*SPB(5,2))/
 (pow(complex<T>(1,0)-S(1,5)/S(3,4),3)*pow(SPA(3,4),2)*
  pow(SPB(4,3),3)*SPA(2,3))-(SPA(1,4)*SPB(3,2)*SPB(5,3))/
 (S(1,5)*SPA(2,3)*SPB(4,3))+(pow(SPB(5,3),2)*SPB(5,2))/
 (SPA(2,3)*SPB(4,3)*SPB(5,1)*SPB(5,4))))/complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpmmm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpmmm L");
#endif
 
 return( (complex<T>(0,-1)*(-((pow(SPB(4,2),2)*SPA(4,5)*SPB(5,2))/
  (pow(SPB(5,4),2)*SPB(2,1)*SPB(3,2)*SPB(4,3)))+
(SPA(1,3)*SPA(1,5))/(SPA(1,2)*SPB(4,3)*SPB(5,4))-
(SPA(3,4)*SPB(3,2)*SPB(4,1)*SPB(4,2))/(pow(SPB(4,3),2)*
  SPB(2,1)*SPB(5,1)*SPB(5,4))))/complex<T>(3,0)+
 (complex<T>(0,-1)*(-(SPA(3,4)*SPB(3,2)*SPB(4,2))-SPA(1,5)*SPB(2,1)*
 SPB(5,2)))/(complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqppmm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqppmm L");
#endif
 
 return(
complex<T>(0,1)*((complex<T>(1,0)*(-((pow(SPA(4,5),3)*pow(SPB(4,3),2)*
 pow(SPB(5,2),2)*(S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
(pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4)*SPB(5,4)))+(complex<T>(-2,0)*pow(SPB(3,2),2)*
SPA(4,5))/(pow(SPB(2,1),2)*SPA(1,2)*SPB(5,4))-
  (S(2,4)*SPA(1,4)*SPA(4,5))/(S(1,2)*SPA(2,3)*SPA(3,4)*
SPB(5,4))-(SPA(1,4)*SPA(2,5)*SPB(3,2))/(S(1,2)*SPA(2,3)*
SPB(5,4))+(complex<T>(-3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
pow(SPA(4,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2))/
   (pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(5,4))))/complex<T>(3,0)+
 (complex<T>(1,0)*(-((SPA(1,4)*SPA(4,5))/(SPA(2,3)*SPA(3,4)*SPB(5,4)))+
  (pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(4,5),2)*SPB(5,2)*
SPB(5,3))/(pow(SPA(3,4),2)*SPB(4,3)*SPB(5,1)*SPB(5,4))))/
complex<T>(2,0))+(complex<T>(0,1)*pow(SPB(3,2),2)*SPB(3,1))/
(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpmpm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpmpm L");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(4,2),3)*SPB(4,1))/(SPB(2,1)*SPB(3,2)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+complex<T>(0,1)*
((complex<T>(1,0)*pow(SPA(3,5),3)*(-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*
 SPB(3,2)*SPB(4,2))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*
 pow(SPA(4,5),3)*pow(SPB(5,4),2)*SPB(2,1))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1))/complex<T>(2,0))*
 pow(SPA(1,3),3)*pow(SPB(4,1),2)*SPB(2,1)*SPB(4,3))/
(pow(SPA(4,5),3)*pow(SPB(5,4),4)*(complex<T>(1,0)-S(1,2)/S(4,5)-
  S(2,3)/S(4,5))*SPB(5,1))+(complex<T>(1,0)*pow(SPA(3,5),4)*SPB(3,2)*
 SPB(5,2))/(complex<T>(3,0)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPA(3,4)*
 SPA(4,5))+(complex<T>(1,0)*pow(SPA(3,5),3)*(-(S(1,2)/S(3,4))+
  S(3,4)/S(1,2))*SPB(4,2)*SPB(5,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(2,1))+
 (complex<T>(1,0)*((pow(SPA(1,3),2)*SPA(3,5)*SPB(4,1))/(S(4,5)*SPA(1,2)*
SPA(3,4)*SPB(5,1))-(pow(SPA(1,3),2)*pow(SPB(4,1),2)*
SPB(4,2))/(pow(SPB(5,4),2)*S(2,3)*SPA(4,5)*SPB(5,1))-
  (pow(SPB(4,2),2)*SPA(1,5)*SPB(4,1))/(S(1,2)*SPB(3,2)*
SPB(4,3)*SPB(5,1))-(pow(SPA(3,5),2)*SPA(2,3)*SPB(4,2)*
SPB(5,2))/(S(1,2)*S(4,5)*SPA(3,4)*SPB(5,1))-
  (pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(3,5),2)*SPB(4,2))/
   (pow(SPA(4,5),2)*SPB(5,1)*SPB(5,4))-
  (pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*pow(SPA(1,3),2)*
pow(SPB(4,1),2)*SPB(4,2))/(pow(SPA(2,3),2)*pow(SPB(3,2),2)*
SPB(5,1)*SPB(5,4))-(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*
pow(SPA(1,3),2)*SPB(4,1)*SPB(4,2)*SPB(4,3))/
   (pow(SPA(1,2),2)*SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,4))))/
complex<T>(2,0)+(complex<T>(-2,0)*pow(SPA(3,5),4)*(-(S(1,2)/(complex<T>(6,0)*S(3,4)))+
  (complex<T>(1,0)*(S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPB(3,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,4))/(pow(SPA(1,2),4)*pow(SPB(2,1),5)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2)))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPA(3,5),3)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
  S(4,5)/S(1,2))*SPB(5,1)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpmmp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpmmp L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(1,4),2)*SPA(3,4)*SPB(4,2))/
 (complex<T>(2,0)*S(2,3)*SPA(1,5)*SPA(4,5)*SPB(4,3))+
(complex<T>(1,0)*pow(SPB(5,2),2)*SPA(1,2)*SPA(1,3)*SPB(5,1))/
 (complex<T>(2,0)*pow(SPB(5,4),2)*S(2,3)*SPA(4,5)*SPB(4,3))-
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1))/complex<T>(2,0))*
pow(SPA(1,3),3)*pow(SPB(5,3),2)*SPB(2,1)*SPB(5,1))/
 (pow(SPA(4,5),3)*pow(SPB(5,4),4)*(complex<T>(1,0)-S(1,2)/S(4,5)-
 S(2,3)/S(4,5))*SPB(4,3))+(SPA(1,3)*SPA(3,4)*SPB(5,1)*
SPB(5,2))/(S(4,5)*SPA(2,3)*SPB(2,1)*SPB(4,3))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),-1)*pow(SPA(1,4),2)*SPB(4,2)*
SPB(5,1)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(2,3),2)*pow(SPB(3,2),3)*
SPB(4,3))+(complex<T>(1,0)*(-((pow(SPA(1,4),2)*SPA(1,3))/
   (SPA(1,2)*SPA(1,5)*SPA(4,5)*SPB(4,3)))+
 (SPA(1,4)*SPA(3,4)*SPB(5,2))/(S(1,2)*SPA(4,5)*SPB(4,3))+
 (pow(SPA(3,4),2)*(-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*SPA(1,3)*
   SPB(3,2)*SPB(5,3))/(pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*
   pow(SPA(4,5),3)*pow(SPB(5,4),2)*SPB(4,3))+
 (pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(3,4),2)*SPB(3,2)*
   SPB(5,2))/(pow(SPA(4,5),2)*SPB(2,1)*SPB(4,3)*SPB(5,4))))/
 complex<T>(3,0)-(pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*
(-((pow(SPB(5,2),2)*SPA(1,2)*SPA(1,3)*SPB(5,1))/
   (complex<T>(2,0)*SPB(4,3)*SPB(5,4)))-(pow(SPA(1,3),2)*
   pow(SPB(5,1),2)*SPB(3,2)*SPB(5,2))/(SPB(2,1)*SPB(4,3)*
   SPB(5,4))))/(pow(SPA(2,3),2)*pow(SPB(3,2),2))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(3,4),2)*SPB(3,2)*
SPB(5,2))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(2,1)*SPB(4,3)*
SPB(5,4))-(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
pow(SPA(1,4),3)*SPB(2,1)*SPB(4,2)*SPB(5,1)*SPB(5,4))/
 (pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-
 S(4,5)/S(2,3))*SPB(4,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmqpmm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmqpmm L");
#endif
 
 return( (complex<T>(0,1)*SPA(4,5)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*
 SPB(2,1)*SPB(3,2))+
 (complex<T>(0,1)*((pow(SPB(3,1),2)*SPA(1,2))/(SPB(2,1)*SPB(4,3)*SPB(5,1)*
  SPB(5,4))+(SPA(4,5)*SPB(3,1)*SPB(5,3))/(SPB(2,1)*SPB(3,2)*
  SPB(5,1)*SPB(5,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpqpmm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpqpmm L");
#endif
 
 return( complex<T>(0,1)*(-((SPA(1,4)*SPA(1,5))/(complex<T>(3,0)*SPA(1,2)*SPA(2,3)*
 SPB(5,4)))-(SPA(1,3)*SPA(1,4)*SPA(4,5))/
 (complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(5,4))-
(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(4,5),2)*
pow(SPB(5,2),2))/(complex<T>(2,0)*pow(SPA(3,4),2)*SPB(4,3)*SPB(5,1)*
SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmqppm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmqppm L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(2,5),2))/(complex<T>(2,0)*SPA(3,4)*SPA(4,5)*
SPB(5,1))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*
pow(SPA(1,5),2)*pow(SPB(3,1),2)*SPB(4,1))/
 (complex<T>(2,0)*pow(SPA(4,5),2)*SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmqpmp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmqpmp L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(1,4),3)*pow(SPB(3,1),2))/
 (complex<T>(2,0)*pow(SPB(3,2),2)*SPA(1,5)*SPA(2,3)*SPA(4,5)*SPB(2,1))-
(pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*pow(SPA(1,4),2)*SPB(3,1)*
SPB(5,3))/(complex<T>(2,0)*pow(SPA(1,5),2)*SPB(2,1)*SPB(3,2)*
SPB(5,1))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*
pow(SPA(2,4),2)*SPB(5,2)*SPB(5,3))/(complex<T>(2,0)*pow(SPA(1,5),2)*
pow(SPB(5,1),2)*SPB(2,1))-(SPA(2,4)*SPB(5,2)*SPB(5,3))/
 (complex<T>(2,0)*S(1,5)*SPB(2,1)*SPB(4,2))-
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
pow(SPA(1,4),3)*pow(SPB(3,1),2)*SPB(5,1)*SPB(5,4))/
 (pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-
 S(4,5)/S(2,3))*SPB(2,1))-
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))*
pow(SPA(2,4),3)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
 (pow(SPA(1,5),3)*pow(SPB(5,1),3)*(complex<T>(1,0)-S(2,3)/S(1,5)-
 S(3,4)/S(1,5))*SPB(2,1))-(pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*
SPA(2,4)*SPA(3,4)*SPB(3,2)*SPB(5,3)*SPB(5,4))/
 (complex<T>(2,0)*pow(SPA(1,5),2)*pow(SPB(5,1),2)*SPB(2,1)*SPB(4,2)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmmqpm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmmqpm L");
#endif
 
 return( (complex<T>(0,1)*((pow(SPB(4,1),2)*SPA(1,5))/(SPB(2,1)*SPB(3,2)*
 SPB(4,3)*SPB(5,1))+(SPA(2,3)*SPB(4,1)*SPB(4,2))/
(SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpmqpm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpmqpm L");
#endif
 
 return( complex<T>(0,1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1))/complex<T>(2,0))*
 pow(SPA(1,3),3)*pow(SPB(4,1),2)*SPB(2,1)*SPB(3,2))/
(pow(SPA(4,5),3)*pow(SPB(5,4),4)*(complex<T>(1,0)-S(1,2)/S(4,5)-
  S(2,3)/S(4,5))*SPB(5,1)))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*pow(SPA(3,5),2)*SPB(4,2)*
SPB(5,2))/(complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPB(5,1))-
(SPA(3,5)*SPB(4,2)*SPB(5,2))/(complex<T>(2,0)*S(1,2)*SPB(5,1)*SPB(5,3))+
(complex<T>(1,0)*pow(SPA(1,3),3)*pow(SPB(4,1),2))/(complex<T>(2,0)*S(4,5)*
SPA(1,2)*SPA(2,3)*SPB(5,1)*SPB(5,4))-
(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*pow(SPA(1,3),2)*SPB(4,1)*
SPB(4,2))/(complex<T>(2,0)*pow(SPA(1,2),2)*SPB(2,1)*SPB(5,1)*
SPB(5,4))-(((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
pow(SPA(3,5),3)*SPB(3,2)*SPB(5,2)*SPB(5,4))/
 (pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
 S(4,5)/S(1,2))*SPB(5,1))-(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*
SPA(3,4)*SPA(3,5)*SPB(3,2)*SPB(4,2)*SPB(5,4))/
 (complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPB(5,1)*SPB(5,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmpqpm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, qp, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmpqpm L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(2,5),2))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*
SPB(2,1))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*
pow(SPA(1,2),2)*pow(SPB(4,1),2)*SPB(3,1))/
 (complex<T>(2,0)*pow(SPA(2,3),2)*SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmmqpp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, qp, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmmqpp L");
#endif
 
 return( (complex<T>(0,1)*SPA(1,2)*SPA(1,3))/(complex<T>(3,0)*SPA(1,5)*SPA(4,5)*
 SPB(3,2))+complex<T>(0,1)*(-((SPA(1,2)*SPA(1,3))/(complex<T>(3,0)*SPA(1,5)*
  SPA(4,5)*SPB(3,2)))-(SPA(1,3)*SPA(1,4)*SPA(2,3))/
(complex<T>(2,0)*SPA(1,5)*SPA(3,4)*SPA(4,5)*SPB(3,2))-
 (pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*pow(SPA(2,3),2)*
 pow(SPB(5,2),2))/(complex<T>(2,0)*pow(SPA(3,4),2)*SPB(2,1)*SPB(3,2)*
 SPB(4,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmmmqp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, m, qp}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmmmqp L");
#endif
 
 return( (complex<T>(0,1)*(SPA(1,2)*SPB(5,1)*SPB(5,2)+SPA(3,4)*SPB(5,3)*
SPB(5,4)))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpmmqp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, m, qp}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpmmqp L");
#endif
 
 return( (complex<T>(0,-1)*(-((pow(SPA(1,3),2)*SPA(1,4))/(SPA(1,2)*SPA(1,5)*
   SPA(2,3)*SPB(4,3)))+(SPA(1,3)*SPA(3,4)*SPB(5,2))/
 (S(1,5)*SPA(2,3)*SPB(4,3))+(pow(SPA(3,4),2)*
  (-(S(1,5)/S(2,3))+S(2,3)/S(1,5))*SPA(1,4)*SPB(4,2)*
  SPB(5,4))/(pow(complex<T>(1,0)-S(1,5)/S(2,3),3)*pow(SPA(2,3),3)*
  pow(SPB(3,2),2)*SPB(4,3))+(pow(complex<T>(1,0)-S(1,5)/S(2,3),-1)*
  pow(SPA(3,4),2)*SPB(5,2)*SPB(5,4))/(pow(SPA(2,3),2)*SPB(3,2)*
  SPB(4,3)*SPB(5,1))))/complex<T>(3,0)+
 complex<T>(0,1)*((complex<T>(1,0)*pow(SPB(5,2),2)*SPA(1,4)*SPA(1,5)*SPB(2,1))/
(complex<T>(2,0)*pow(SPB(3,2),2)*S(4,5)*SPA(2,3)*SPB(4,3))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
 pow(SPA(1,4),3)*pow(SPB(4,2),2)*SPB(2,1)*SPB(5,1))/
(pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-
  S(4,5)/S(2,3))*SPB(4,3))+(SPA(1,4)*SPA(3,4)*SPB(2,1)*
 SPB(5,2))/(S(2,3)*SPA(4,5)*SPB(4,3)*SPB(5,1))+
 (complex<T>(1,0)*pow(SPA(1,3),2)*SPA(3,4)*SPB(5,3))/
(complex<T>(2,0)*S(4,5)*SPA(1,2)*SPA(2,3)*SPB(4,3))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1))/complex<T>(2,0))*
 pow(SPA(1,3),3)*SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,3))/
(pow(SPA(4,5),3)*pow(SPB(5,4),4)*(complex<T>(1,0)-S(1,2)/S(4,5)-
  S(2,3)/S(4,5))*SPB(4,3))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(1,3),2)*
 SPB(2,1)*SPB(5,2)*SPB(5,3))/(complex<T>(2,0)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPB(4,3))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),-1)*pow(SPA(3,4),2)*
 SPB(5,2)*SPB(5,4))/(complex<T>(2,0)*pow(SPA(2,3),2)*SPB(3,2)*SPB(4,3)*
 SPB(5,1))+(complex<T>(1,0)*(-((pow(SPA(1,3),2)*SPA(1,4))/
(SPA(1,2)*SPA(1,5)*SPA(2,3)*SPB(4,3)))+
  (SPA(1,3)*SPA(3,4)*SPB(5,2))/(S(1,5)*SPA(2,3)*SPB(4,3))+
  (pow(SPA(3,4),2)*(-(S(1,5)/S(2,3))+S(2,3)/S(1,5))*SPA(1,4)*
SPB(4,2)*SPB(5,4))/(pow(complex<T>(1,0)-S(1,5)/S(2,3),3)*
pow(SPA(2,3),3)*pow(SPB(3,2),2)*SPB(4,3))+
  (pow(complex<T>(1,0)-S(1,5)/S(2,3),-1)*pow(SPA(3,4),2)*SPB(5,2)*
SPB(5,4))/(pow(SPA(2,3),2)*SPB(3,2)*SPB(4,3)*SPB(5,1))))/
complex<T>(3,0)-(pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*
 (-((pow(SPB(5,2),2)*SPA(1,4)*SPA(1,5)*SPB(2,1))/
(complex<T>(2,0)*SPB(3,2)*SPB(4,3)))-(pow(SPA(1,4),2)*
pow(SPB(2,1),2)*SPB(5,2)*SPB(5,4))/(SPB(3,2)*SPB(4,3)*
SPB(5,1))))/(pow(SPA(4,5),2)*pow(SPB(5,4),2)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmpmqp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, m, qp}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmpmqp L");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(5,3),3)*SPB(3,1))/(SPB(2,1)*SPB(3,2)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+complex<T>(0,-1)*
((complex<T>(1,0)*pow(SPA(2,4),3)*(-(S(1,5)/S(3,4))+S(3,4)/S(1,5))*
 SPB(5,2)*SPB(5,3))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),3)*
 pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPB(5,1))+
 (complex<T>(1,0)*pow(SPA(2,4),4)*SPB(5,2)*SPB(5,4))/
(complex<T>(3,0)*pow(SPA(1,5),2)*pow(SPB(5,1),3)*SPA(2,3)*SPA(3,4))+
 (complex<T>(-2,0)*pow(SPA(2,4),4)*(-(S(1,5)/(complex<T>(6,0)*S(2,3)))+
  (complex<T>(1,0)*(S(1,5)/S(2,3)-S(2,3)/S(1,5)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),3))-
  S(1,5)/(complex<T>(6,0)*S(3,4))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5))+
  (complex<T>(1,0)*(S(1,5)/S(3,4)-S(3,4)/S(1,5)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),3)))*SPB(3,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,4))/(pow(SPA(1,5),4)*pow(SPB(5,1),5)*
 (complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5)))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*(-(S(1,5)/S(2,3))+S(2,3)/S(1,5))*
 SPB(5,3)*SPB(5,4))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),3)*
 pow(SPA(2,3),3)*pow(SPB(3,2),2)*SPB(5,1)))+
 complex<T>(0,1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),-1))/complex<T>(2,0)+
   (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
  pow(SPA(1,4),3)*pow(SPB(3,1),2)*SPB(4,3)*SPB(5,1))/
 (pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-
   S(4,5)/S(2,3))*SPB(2,1)))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))*
 pow(SPA(2,4),3)*SPB(3,2)*SPB(4,3)*SPB(5,2))/
(pow(SPA(1,5),3)*pow(SPB(5,1),3)*(complex<T>(1,0)-S(2,3)/S(1,5)-
  S(3,4)/S(1,5))*SPB(2,1))+(complex<T>(1,0)*pow(SPA(2,4),3)*
 (-(S(1,5)/S(3,4))+S(3,4)/S(1,5))*SPB(5,2)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),3)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(5,1))+
 (complex<T>(1,0)*((pow(SPA(1,4),2)*SPA(2,4)*SPB(3,1))/(S(2,3)*SPA(1,5)*
SPA(3,4)*SPB(2,1))-(pow(complex<T>(1,0)-S(1,5)/S(2,3),-1)*
pow(SPA(2,4),2)*SPB(5,3))/(pow(SPA(2,3),2)*SPB(2,1)*
SPB(3,2))-(pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*
pow(SPA(1,4),2)*pow(SPB(3,1),2)*SPB(5,3))/
   (pow(SPA(4,5),2)*pow(SPB(5,4),2)*SPB(2,1)*SPB(3,2))-
  (pow(SPA(1,4),2)*pow(SPB(3,1),2)*SPB(5,3))/(S(2,3)*S(4,5)*
SPB(2,1)*SPB(3,2))-(pow(SPA(2,4),2)*SPA(4,5)*SPB(5,2)*
SPB(5,3))/(S(1,5)*S(2,3)*SPA(3,4)*SPB(2,1))-
  (pow(SPB(5,3),2)*SPA(1,2)*SPB(3,1))/(S(1,5)*SPB(2,1)*
SPB(4,3)*SPB(5,4))-(pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*
pow(SPA(1,4),2)*SPB(3,1)*SPB(4,3)*SPB(5,3))/
   (pow(SPA(1,5),2)*SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,4))))/
complex<T>(2,0)+(complex<T>(1,0)*pow(SPA(2,4),4)*SPB(5,2)*SPB(5,4))/
(complex<T>(3,0)*pow(SPA(1,5),2)*pow(SPB(5,1),3)*SPA(2,3)*SPA(3,4))+
 (complex<T>(-2,0)*pow(SPA(2,4),4)*(-(S(1,5)/(complex<T>(6,0)*S(2,3)))+
  (complex<T>(1,0)*(S(1,5)/S(2,3)-S(2,3)/S(1,5)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),3))-
  S(1,5)/(complex<T>(6,0)*S(3,4))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5))+
  (complex<T>(1,0)*(S(1,5)/S(3,4)-S(3,4)/S(1,5)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),3)))*SPB(3,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,4))/(pow(SPA(1,5),4)*pow(SPB(5,1),5)*
 (complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5)))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*(-(S(1,5)/S(2,3))+S(2,3)/S(1,5))*
 SPB(5,3)*SPB(5,4))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),3)*
 pow(SPA(2,3),3)*pow(SPB(3,2),2)*SPB(5,1)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmmpqp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, p, qp}, L}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmmpqp L");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(5,4),2)*SPB(4,1))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*
 SPB(4,3)*SPB(5,1))+
 (complex<T>(0,-1)*(-((pow(SPA(2,3),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2)*
   (S(1,5)/S(3,4)-S(3,4)/S(1,5)))/
  (pow(complex<T>(1,0)-S(3,4)/S(1,5),3)*pow(SPA(1,5),3)*
   pow(SPB(5,1),4)*SPB(3,2)))+(complex<T>(-2,0)*pow(SPB(5,4),2)*
  SPA(2,3))/(pow(SPB(5,1),2)*SPA(1,5)*SPB(3,2))-
(S(3,5)*SPA(1,3)*SPA(2,3))/(S(1,5)*SPA(3,4)*SPA(4,5)*
  SPB(3,2))-(SPA(1,3)*SPA(2,5)*SPB(5,4))/(S(1,5)*SPA(4,5)*
  SPB(3,2))+(complex<T>(-3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*
  pow(SPA(2,3),2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
 (pow(SPA(1,5),2)*pow(SPB(5,1),3)*SPB(3,2))))/complex<T>(3,0)+
 complex<T>(0,1)*((complex<T>(1,0)*(-((SPA(1,3)*SPA(2,3))/(SPA(3,4)*SPA(4,5)*
 SPB(3,2)))+(pow(complex<T>(1,0)-S(1,5)/S(3,4),-1)*
pow(SPA(2,3),2)*SPB(4,2)*SPB(5,2))/(pow(SPA(3,4),2)*
SPB(2,1)*SPB(3,2)*SPB(4,3))))/complex<T>(2,0)+
 (complex<T>(1,0)*(-((pow(SPA(2,3),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2)*
 (S(1,5)/S(3,4)-S(3,4)/S(1,5)))/
(pow(complex<T>(1,0)-S(3,4)/S(1,5),3)*pow(SPA(1,5),3)*
 pow(SPB(5,1),4)*SPB(3,2)))+(complex<T>(-2,0)*pow(SPB(5,4),2)*
SPA(2,3))/(pow(SPB(5,1),2)*SPA(1,5)*SPB(3,2))-
  (S(3,5)*SPA(1,3)*SPA(2,3))/(S(1,5)*SPA(3,4)*SPA(4,5)*
SPB(3,2))-(SPA(1,3)*SPA(2,5)*SPB(5,4))/(S(1,5)*SPA(4,5)*
SPB(3,2))+(complex<T>(-3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*
pow(SPA(2,3),2)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
   (pow(SPA(1,5),2)*pow(SPB(5,1),3)*SPB(3,2))))/complex<T>(3,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpppp nf");
#endif
 
 return( (complex<T>(0,1)*((pow(SPA(1,4),2)*SPA(1,3)*SPB(4,3))/
(pow(SPA(3,4),2)*SPA(1,2)*SPA(1,5)*SPA(4,5))-
 (SPB(3,2)*SPB(5,2))/(SPA(3,4)*SPA(4,5)*SPB(2,1))+
 (SPA(1,4)*SPA(1,5)*SPA(2,4)*SPB(5,4))/(pow(SPA(4,5),2)*SPA(1,2)*
 SPA(2,3)*SPA(3,4))))/complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpmpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpmpp nf");
#endif
 
 return( (complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPB(5,4),2)*
  SPA(1,3)*SPA(1,5))/(pow(SPB(4,3),2)*SPA(1,2)*SPA(3,4)*
  SPA(4,5)))-(pow(SPB(5,4),2)*(-(S(1,2)/S(3,4))+
  S(3,4)/S(1,2))*SPA(1,5)*SPA(3,5)*SPB(5,2))/
(pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
 pow(SPB(4,3),3)*SPA(4,5))+(pow(SPB(4,2),2)*SPB(5,2))/
(SPA(4,5)*SPB(2,1)*SPB(3,2)*SPB(4,3))-
 (SPA(1,3)*SPB(4,2)*SPB(5,4))/(S(1,2)*SPA(4,5)*SPB(4,3))))/complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqppmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqppmp nf");
#endif
 
 return( complex<T>(0,-1)*(-((pow(SPB(5,3),3)*(-(S(1,2)/S(4,5))+
  S(4,5)/S(1,2))*SPA(1,3)*SPA(1,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(1,2)))-
(pow(SPB(5,3),3)*(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPA(1,4)*
SPA(1,5))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPA(1,2))+
(complex<T>(2,0)*pow(SPB(5,3),4)*(-(S(1,2)/(complex<T>(6,0)*S(3,4)))+
 (complex<T>(1,0)*(S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
  (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
 S(1,2)/(complex<T>(6,0)*S(4,5))+
 ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
   (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
  (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))+
 (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
  (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPA(1,3)*SPA(1,5)*
SPA(3,4)*SPA(4,5))/(pow(SPA(1,2),5)*pow(SPB(2,1),4)*
(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2)))-
(pow(SPB(5,3),4)*SPA(1,3)*SPA(1,5))/(complex<T>(3,0)*pow(SPA(1,2),3)*
pow(SPB(2,1),2)*SPB(4,3)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpppm nf");
#endif
 
 return( (complex<T>(0,-1)*((pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(4,3),3)*
 (S(1,2)/S(4,5)-S(4,5)/S(1,2)))/(pow(complex<T>(1,0)-S(4,5)/S(1,2),
  3)*pow(SPA(1,2),4)*pow(SPB(2,1),3)*SPA(3,4))+
 (complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*pow(SPB(4,3),2)*
 SPA(1,3)*SPA(1,5)*SPA(4,5))/(pow(SPA(1,2),3)*pow(SPB(2,1),2)*
 SPA(3,4))+(complex<T>(2,0)*pow(SPA(1,5),2)*SPB(4,3))/
(pow(SPA(1,2),2)*SPA(3,4)*SPB(2,1))+
 (SPA(1,5)*SPB(3,1)*SPB(4,2))/(S(1,2)*SPA(3,4)*SPB(5,1))+
 (S(1,4)*SPB(4,2)*SPB(4,3))/(S(1,2)*SPA(3,4)*SPB(5,1)*
 SPB(5,4))))/complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpqppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpqppp nf");
#endif
 
 return( (complex<T>(0,1)*SPA(1,4)*SPA(1,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(4,5),2)*
SPA(1,2)*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmqppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmqppp nf");
#endif
 
 return( (complex<T>(0,-1)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*SPA(4,5)*SPB(2,1)*
SPB(3,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpqpmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpqpmp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmpqppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpqppm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmppqpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmppqpp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmpqpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmpqpp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmpmqpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpmqpp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmppqpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmppqpm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmpppqp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpppqp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmppqp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmppqp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmpmpqp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpmpqp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmppmqp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, m, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmppmqp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmqpmmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpmmm nf");
#endif
 
 return( (complex<T>(0,1)*(-((pow(SPB(4,2),2)*SPA(4,5)*SPB(5,2))/
 (pow(SPB(5,4),2)*SPB(2,1)*SPB(3,2)*SPB(4,3)))+
 (SPA(1,3)*SPA(1,5))/(SPA(1,2)*SPB(4,3)*SPB(5,4))-
 (SPA(3,4)*SPB(3,2)*SPB(4,1)*SPB(4,2))/(pow(SPB(4,3),2)*SPB(2,1)*
 SPB(5,1)*SPB(5,4))))/complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqppmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqppmm nf");
#endif
 
 return( (complex<T>(0,-1)*(-((pow(SPA(4,5),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2)*
  (S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
 (pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPA(1,2),3)*
  pow(SPB(2,1),4)*SPB(5,4)))+(complex<T>(-2,0)*pow(SPB(3,2),2)*
 SPA(4,5))/(pow(SPB(2,1),2)*SPA(1,2)*SPB(5,4))-
 (S(2,4)*SPA(1,4)*SPA(4,5))/(S(1,2)*SPA(2,3)*SPA(3,4)*
 SPB(5,4))-(SPA(1,4)*SPA(2,5)*SPB(3,2))/(S(1,2)*SPA(2,3)*
 SPB(5,4))+(complex<T>(-3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
 pow(SPA(4,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2))/
(pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(5,4))))/complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpmpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpmpm nf");
#endif
 
 return( complex<T>(0,-1)*((complex<T>(1,0)*pow(SPA(3,5),3)*(-(S(1,2)/S(4,5))+
 S(4,5)/S(1,2))*SPB(3,2)*SPB(4,2))/
 (complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),3)*
pow(SPB(5,4),2)*SPB(2,1))+(complex<T>(1,0)*pow(SPA(3,5),4)*SPB(3,2)*
SPB(5,2))/(complex<T>(3,0)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPA(3,4)*
SPA(4,5))+(complex<T>(1,0)*pow(SPA(3,5),3)*(-(S(1,2)/S(3,4))+
 S(3,4)/S(1,2))*SPB(4,2)*SPB(5,2))/
 (complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),3)*
pow(SPB(4,3),2)*SPB(2,1))+(complex<T>(-2,0)*pow(SPA(3,5),4)*
(-(S(1,2)/(complex<T>(6,0)*S(3,4)))+(complex<T>(1,0)*(S(1,2)/S(3,4)-
S(3,4)/S(1,2)))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
 S(1,2)/(complex<T>(6,0)*S(4,5))+
 ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
   (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
  (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))+
 (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
  (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPB(3,2)*SPB(4,3)*
SPB(5,2)*SPB(5,4))/(pow(SPA(1,2),4)*pow(SPB(2,1),5)*
(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmqpmmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmqpmmp nf");
#endif
 
 return( (complex<T>(0,-1)*(-((pow(SPA(1,4),2)*SPA(1,3))/(SPA(1,2)*SPA(1,5)*
  SPA(4,5)*SPB(4,3)))+(SPA(1,4)*SPA(3,4)*SPB(5,2))/
(S(1,2)*SPA(4,5)*SPB(4,3))+
 (pow(SPA(3,4),2)*(-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*SPA(1,3)*
 SPB(3,2)*SPB(5,3))/(pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*
 pow(SPA(4,5),3)*pow(SPB(5,4),2)*SPB(4,3))+
 (pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(3,4),2)*SPB(3,2)*
 SPB(5,2))/(pow(SPA(4,5),2)*SPB(2,1)*SPB(4,3)*SPB(5,4))))/
 complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmqpmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmqpmm nf");
#endif
 
 return( (complex<T>(0,-1)*SPA(4,5)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*
SPB(2,1)*SPB(3,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmpqpmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpqpmm nf");
#endif
 
 return( (complex<T>(0,1)*SPA(1,4)*SPA(1,5))/(complex<T>(3,0)*SPA(1,2)*SPA(2,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q3g_qmmqppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmqppm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmqpmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmqpmp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmmqpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmmqpm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmpmqpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpmqpm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmpqpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, qp, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmpqpm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmmqpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, qp, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmmqpp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmmmqp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, m, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmmmqp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmpmmqp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, m, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmpmmqp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmpmqp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, m, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmpmqp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q3g_qmmmpqp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, p, qp}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q3g :  qmmmpqp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 // *************** table of switch values ************* 
 
#define _R_qmqpppp_L R2q3g_1017_L
#define _R_qmqpmpp_L R2q3g_969_L
#define _R_qmqppmp_L R2q3g_825_L
#define _R_qmqpppm_L R2q3g_249_L
#define _R_qmpqppp_L R2q3g_1005_L
#define _R_qmmqppp_L R2q3g_993_L
#define _R_qmpqpmp_L R2q3g_813_L
#define _R_qmpqppm_L R2q3g_237_L
#define _R_qmppqpp_L R2q3g_957_L
#define _R_qmmpqpp_L R2q3g_945_L
#define _R_qmpmqpp_L R2q3g_909_L
#define _R_qmppqpm_L R2q3g_189_L
#define _R_qmpppqp_L R2q3g_765_L
#define _R_qmmppqp_L R2q3g_753_L
#define _R_qmpmpqp_L R2q3g_717_L
#define _R_qmppmqp_L R2q3g_573_L
#define _R_qmqpmmm_L R2q3g_9_L
#define _R_qmqppmm_L R2q3g_57_L
#define _R_qmqpmpm_L R2q3g_201_L
#define _R_qmqpmmp_L R2q3g_777_L
#define _R_qmmqpmm_L R2q3g_33_L
#define _R_qmpqpmm_L R2q3g_45_L
#define _R_qmmqppm_L R2q3g_225_L
#define _R_qmmqpmp_L R2q3g_801_L
#define _R_qmmmqpm_L R2q3g_129_L
#define _R_qmpmqpm_L R2q3g_141_L
#define _R_qmmpqpm_L R2q3g_177_L
#define _R_qmmmqpp_L R2q3g_897_L
#define _R_qmmmmqp_L R2q3g_513_L
#define _R_qmpmmqp_L R2q3g_525_L
#define _R_qmmpmqp_L R2q3g_561_L
#define _R_qmmmpqp_L R2q3g_705_L
#define _R_qmqpppp_nf R2q3g_1017_nf
#define _R_qmqpmpp_nf R2q3g_969_nf
#define _R_qmqppmp_nf R2q3g_825_nf
#define _R_qmqpppm_nf R2q3g_249_nf
#define _R_qmpqppp_nf R2q3g_1005_nf
#define _R_qmmqppp_nf R2q3g_993_nf
#define _R_qmpqpmp_nf R2q3g_813_nf
#define _R_qmpqppm_nf R2q3g_237_nf
#define _R_qmppqpp_nf R2q3g_957_nf
#define _R_qmmpqpp_nf R2q3g_945_nf
#define _R_qmpmqpp_nf R2q3g_909_nf
#define _R_qmppqpm_nf R2q3g_189_nf
#define _R_qmpppqp_nf R2q3g_765_nf
#define _R_qmmppqp_nf R2q3g_753_nf
#define _R_qmpmpqp_nf R2q3g_717_nf
#define _R_qmppmqp_nf R2q3g_573_nf
#define _R_qmqpmmm_nf R2q3g_9_nf
#define _R_qmqppmm_nf R2q3g_57_nf
#define _R_qmqpmpm_nf R2q3g_201_nf
#define _R_qmqpmmp_nf R2q3g_777_nf
#define _R_qmmqpmm_nf R2q3g_33_nf
#define _R_qmpqpmm_nf R2q3g_45_nf
#define _R_qmmqppm_nf R2q3g_225_nf
#define _R_qmmqpmp_nf R2q3g_801_nf
#define _R_qmmmqpm_nf R2q3g_129_nf
#define _R_qmpmqpm_nf R2q3g_141_nf
#define _R_qmmpqpm_nf R2q3g_177_nf
#define _R_qmmmqpp_nf R2q3g_897_nf
#define _R_qmmmmqp_nf R2q3g_513_nf
#define _R_qmpmmqp_nf R2q3g_525_nf
#define _R_qmmpmqp_nf R2q3g_561_nf
#define _R_qmmmpqp_nf R2q3g_705_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmqpppp_L case 1017 : \
          return &R2q3g_1017_L
#define _CASE_qmqpmpp_L case 969 : \
          return &R2q3g_969_L
#define _CASE_qmqppmp_L case 825 : \
          return &R2q3g_825_L
#define _CASE_qmqpppm_L case 249 : \
          return &R2q3g_249_L
#define _CASE_qmpqppp_L case 1005 : \
          return &R2q3g_1005_L
#define _CASE_qmmqppp_L case 993 : \
          return &R2q3g_993_L
#define _CASE_qmpqpmp_L case 813 : \
          return &R2q3g_813_L
#define _CASE_qmpqppm_L case 237 : \
          return &R2q3g_237_L
#define _CASE_qmppqpp_L case 957 : \
          return &R2q3g_957_L
#define _CASE_qmmpqpp_L case 945 : \
          return &R2q3g_945_L
#define _CASE_qmpmqpp_L case 909 : \
          return &R2q3g_909_L
#define _CASE_qmppqpm_L case 189 : \
          return &R2q3g_189_L
#define _CASE_qmpppqp_L case 765 : \
          return &R2q3g_765_L
#define _CASE_qmmppqp_L case 753 : \
          return &R2q3g_753_L
#define _CASE_qmpmpqp_L case 717 : \
          return &R2q3g_717_L
#define _CASE_qmppmqp_L case 573 : \
          return &R2q3g_573_L
#define _CASE_qmqpmmm_L case 9 : \
          return &R2q3g_9_L
#define _CASE_qmqppmm_L case 57 : \
          return &R2q3g_57_L
#define _CASE_qmqpmpm_L case 201 : \
          return &R2q3g_201_L
#define _CASE_qmqpmmp_L case 777 : \
          return &R2q3g_777_L
#define _CASE_qmmqpmm_L case 33 : \
          return &R2q3g_33_L
#define _CASE_qmpqpmm_L case 45 : \
          return &R2q3g_45_L
#define _CASE_qmmqppm_L case 225 : \
          return &R2q3g_225_L
#define _CASE_qmmqpmp_L case 801 : \
          return &R2q3g_801_L
#define _CASE_qmmmqpm_L case 129 : \
          return &R2q3g_129_L
#define _CASE_qmpmqpm_L case 141 : \
          return &R2q3g_141_L
#define _CASE_qmmpqpm_L case 177 : \
          return &R2q3g_177_L
#define _CASE_qmmmqpp_L case 897 : \
          return &R2q3g_897_L
#define _CASE_qmmmmqp_L case 513 : \
          return &R2q3g_513_L
#define _CASE_qmpmmqp_L case 525 : \
          return &R2q3g_525_L
#define _CASE_qmmpmqp_L case 561 : \
          return &R2q3g_561_L
#define _CASE_qmmmpqp_L case 705 : \
          return &R2q3g_705_L
#define _CASE_qmqpppp_nf case 1017 : \
          return &R2q3g_1017_nf
#define _CASE_qmqpmpp_nf case 969 : \
          return &R2q3g_969_nf
#define _CASE_qmqppmp_nf case 825 : \
          return &R2q3g_825_nf
#define _CASE_qmqpppm_nf case 249 : \
          return &R2q3g_249_nf
#define _CASE_qmpqppp_nf case 1005 : \
          return &R2q3g_1005_nf
#define _CASE_qmmqppp_nf case 993 : \
          return &R2q3g_993_nf
#define _CASE_qmpqpmp_nf case 813 : \
          return &R2q3g_813_nf
#define _CASE_qmpqppm_nf case 237 : \
          return &R2q3g_237_nf
#define _CASE_qmppqpp_nf case 957 : \
          return &R2q3g_957_nf
#define _CASE_qmmpqpp_nf case 945 : \
          return &R2q3g_945_nf
#define _CASE_qmpmqpp_nf case 909 : \
          return &R2q3g_909_nf
#define _CASE_qmppqpm_nf case 189 : \
          return &R2q3g_189_nf
#define _CASE_qmpppqp_nf case 765 : \
          return &R2q3g_765_nf
#define _CASE_qmmppqp_nf case 753 : \
          return &R2q3g_753_nf
#define _CASE_qmpmpqp_nf case 717 : \
          return &R2q3g_717_nf
#define _CASE_qmppmqp_nf case 573 : \
          return &R2q3g_573_nf
#define _CASE_qmqpmmm_nf case 9 : \
          return &R2q3g_9_nf
#define _CASE_qmqppmm_nf case 57 : \
          return &R2q3g_57_nf
#define _CASE_qmqpmpm_nf case 201 : \
          return &R2q3g_201_nf
#define _CASE_qmqpmmp_nf case 777 : \
          return &R2q3g_777_nf
#define _CASE_qmmqpmm_nf case 33 : \
          return &R2q3g_33_nf
#define _CASE_qmpqpmm_nf case 45 : \
          return &R2q3g_45_nf
#define _CASE_qmmqppm_nf case 225 : \
          return &R2q3g_225_nf
#define _CASE_qmmqpmp_nf case 801 : \
          return &R2q3g_801_nf
#define _CASE_qmmmqpm_nf case 129 : \
          return &R2q3g_129_nf
#define _CASE_qmpmqpm_nf case 141 : \
          return &R2q3g_141_nf
#define _CASE_qmmpqpm_nf case 177 : \
          return &R2q3g_177_nf
#define _CASE_qmmmqpp_nf case 897 : \
          return &R2q3g_897_nf
#define _CASE_qmmmmqp_nf case 513 : \
          return &R2q3g_513_nf
#define _CASE_qmpmmqp_nf case 525 : \
          return &R2q3g_525_nf
#define _CASE_qmmpmqp_nf case 561 : \
          return &R2q3g_561_nf
#define _CASE_qmmmpqp_nf case 705 : \
          return &R2q3g_705_nf
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_qmqpppp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpppp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmpp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpmpp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqppmp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqppmp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqpppm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpppm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpqppp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpqppp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmqppp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmqppp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpmp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpqpmp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpqppm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpqppm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmppqpp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmppqpp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmpqpp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmpqpp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpmqpp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpmqpp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmppqpm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmppqpm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpppqp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpppqp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmppqp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmppqp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpmpqp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpmpqp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmppmqp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmppmqp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmmm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpmmm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqppmm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqppmm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmpm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpmpm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmmp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpmmp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpmm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmqpmm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpmm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpqpmm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmqppm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmqppm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpmp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmqpmp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmmqpm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmmqpm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpmqpm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpmqpm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmpqpm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmpqpm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmmqpp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmmqpp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmmmqp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmmmqp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpmmqp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpmmqp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmpmqp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmpmqp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmmpqp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmmpqp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmqpppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpmpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmqppmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqppmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmqpppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpqppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpqppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmqppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmqppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpqpmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpqppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpqppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmppqpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmppqpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmpqpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmpqpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpmqpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpmqpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmppqpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmppqpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpppqp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpppqp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmppqp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmppqp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpmpqp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpmpqp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmppmqp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmppmqp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpmmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmqppmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqppmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpmpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmqpmmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmqpmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpqpmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmqppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmqppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmqpmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmmqpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmmqpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpmqpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpmqpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmpqpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmpqpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmmqpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmmqpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmmmqp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmmmqp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpmmqp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmpmmqp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmpmqp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmpmqp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmmpqp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g_qmmmpqp_nf(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q3g_L_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qmqpppp_L;
       _CASE_qmqpmpp_L;
       _CASE_qmqppmp_L;
       _CASE_qmqpppm_L;
       _CASE_qmpqppp_L;
       _CASE_qmmqppp_L;
       _CASE_qmpqpmp_L;
       _CASE_qmpqppm_L;
       _CASE_qmppqpp_L;
       _CASE_qmmpqpp_L;
       _CASE_qmpmqpp_L;
       _CASE_qmppqpm_L;
       _CASE_qmpppqp_L;
       _CASE_qmmppqp_L;
       _CASE_qmpmpqp_L;
       _CASE_qmppmqp_L;
       _CASE_qmqpmmm_L;
       _CASE_qmqppmm_L;
       _CASE_qmqpmpm_L;
       _CASE_qmqpmmp_L;
       _CASE_qmmqpmm_L;
       _CASE_qmpqpmm_L;
       _CASE_qmmqppm_L;
       _CASE_qmmqpmp_L;
       _CASE_qmmmqpm_L;
       _CASE_qmpmqpm_L;
       _CASE_qmmpqpm_L;
       _CASE_qmmmqpp_L;
       _CASE_qmmmmqp_L;
       _CASE_qmpmmqp_L;
       _CASE_qmmpmqp_L;
       _CASE_qmmmpqp_L;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q3g_nf_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qmqpppp_nf;
       _CASE_qmqpmpp_nf;
       _CASE_qmqppmp_nf;
       _CASE_qmqpppm_nf;
       _CASE_qmpqppp_nf;
       _CASE_qmmqppp_nf;
       _CASE_qmpqpmp_nf;
       _CASE_qmpqppm_nf;
       _CASE_qmppqpp_nf;
       _CASE_qmmpqpp_nf;
       _CASE_qmpmqpp_nf;
       _CASE_qmppqpm_nf;
       _CASE_qmpppqp_nf;
       _CASE_qmmppqp_nf;
       _CASE_qmpmpqp_nf;
       _CASE_qmppmqp_nf;
       _CASE_qmqpmmm_nf;
       _CASE_qmqppmm_nf;
       _CASE_qmqpmpm_nf;
       _CASE_qmqpmmp_nf;
       _CASE_qmmqpmm_nf;
       _CASE_qmpqpmm_nf;
       _CASE_qmmqppm_nf;
       _CASE_qmmqpmp_nf;
       _CASE_qmmmqpm_nf;
       _CASE_qmpmqpm_nf;
       _CASE_qmmpqpm_nf;
       _CASE_qmmmqpp_nf;
       _CASE_qmmmmqp_nf;
       _CASE_qmpmmqp_nf;
       _CASE_qmmpmqp_nf;
       _CASE_qmmmpqp_nf;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q3g_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q3g_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q3g_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q3g_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q3g_nf_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q3g_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q3g_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q3g_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif




}

