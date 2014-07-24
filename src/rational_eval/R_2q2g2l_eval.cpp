/*
* R_2q2g2l.cpp
*
* Created on 12/2, 2008
*      Author: Zvi's script
*/
 
#include "R_2q2g2l_eval.h"
#include "eval_param.h"
#include "constants.h"
 
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


 
 
template <class T> complex<T> R2q2g2l_qpppqmemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, p, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpppqmemep L");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(4,5),2))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*
 SPA(5,6))+complex<T>(0,1)*
(-((pow(SPA(4,5),2)*(-((pow(SPA(2,5),2)*pow(SPA(3,4),2)*
pow(SPB(3,2),2)*pow(complex<T>(1,0)-S(3,4)/SS(2,3,4),-1))/
 (complex<T>(2,0)*pow(SPA(4,5),2)*pow(SS(2,3,4),2)))+
   (complex<T>(1,0)*pow(-(SPA(1,5)*SPA(2,4)*SPB(2,1))-SPA(1,5)*
  SPA(3,4)*SPB(3,1),2)*pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),
-1))/(complex<T>(2,0)*pow(SPA(4,5),2)*pow(SS(2,3,4),2))))/
 (SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(5,6)))+
 (complex<T>(1,0)*(-((SPA(4,5)*(S(1,3)*SPB(6,1)+SPA(2,3)*SPB(3,1)*SPB(6,
  2)))/SS(1,2,3))-((-(S(3,4)*SPA(4,5))+
 SPA(2,5)*SPA(3,4)*SPB(3,2))*SPB(6,1))/SS(2,3,4)))/
(complex<T>(3,0)*pow(SPA(2,3),2)*S(5,6))+
 (complex<T>(1,0)*((SPA(2,5)*SPA(4,5)*SPB(3,2))/(SPA(1,2)*SPA(2,3)*
SPA(5,6)*SS(2,3,4))-pow(SPA(2,4)*SPB(6,2)+
 SPA(3,4)*SPB(6,3),2)/(SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(6,5)*
SS(2,3,4))+(SPA(4,5)*SPB(3,2)*(SPA(2,5)*SPB(3,2)-
 SPA(4,5)*SPB(4,3)))/(SPA(1,2)*SPA(5,6)*SS(1,2,3)*
SS(2,3,4))-(SPA(2,4)*SPB(3,2)*SPB(6,1)*
(-(SPA(1,2)*SPB(6,2))-SPA(1,3)*SPB(6,3)))/
   (SPA(1,2)*SPA(2,3)*SPB(6,5)*SS(1,2,3)*SS(2,3,4))))/complex<T>(2,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qppmqmemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, m, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qppmqmemep L");
#endif
 
 return( complex<T>(0,1)*(-(pow(SPA(3,5),2)/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(5,6)*
  SPB(4,3)))-pow(SPB(6,2),2)/(complex<T>(2,0)*SPA(1,2)*SPB(3,2)*
 SPB(4,3)*SPB(6,5))+(complex<T>(1,0)*(pow(SPA(3,4),2)*pow(SPB(6,4),2)*
   pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)+
  (pow(SPB(2,1),2)*pow(-(SPA(1,2)*SPB(6,2))-SPA(1,3)*SPB(6,
  3),2)*pow(complex<T>(1,0)-SS(1,2,3)/S(2,3),-1))/
   pow(SPB(3,2),2))*SPA(1,3))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SPB(6,5)*
 SS(1,2,3))+(complex<T>(1,0)*(pow(SPA(1,5),2)*pow(SPB(2,1),2)*
   pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)+
  (pow(SPA(3,4),2)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),
 2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/pow(SPA(2,3),2))*
 SPB(4,2))/(complex<T>(2,0)*SPA(5,6)*SPB(3,2)*SPB(4,3)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SS(2,3,4)))+
 (complex<T>(0,-1)*(((-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*
  (SPA(1,3)*SPB(6,1)+SPA(2,3)*SPB(6,2)))/(S(2,3)*S(5,6)*
  SPA(1,2)*SPB(4,3))+(SPA(1,3)*SPA(4,5)*SPB(2,1)*
  (SPA(1,3)*SPB(6,1)+SPA(2,3)*SPB(6,2)))/(S(2,3)*S(5,6)*
  SPA(1,2)*SS(1,2,3))-(SPA(3,4)*SPB(4,2)*
  (-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(6,1))/
 (S(2,3)*S(5,6)*SPB(4,3)*SS(2,3,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpmpqmemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, p, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpmpqmemep L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
   pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
   complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(1,5)*
 (SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*(SPA(2,3)*SPB(3,1)+
  SPA(2,4)*SPB(4,1))*(-((S(1,2)-S(3,4)-S(5,6))*SPA(1,5)*
SPB(2,1))+(-S(1,2)-S(3,4)+S(5,6))*SPA(5,6)*SPB(6,2)))/
(complex<T>(2,0)*SPA(3,4)*SPA(5,6)*SPB(2,1)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2))*(-(SPA(1,2)*SPB(4,2))-
  SPA(1,3)*SPB(4,3)))-(SPA(1,5)*SPB(3,1)*
 (-(SPA(1,5)*(-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1)))+
  SPA(2,4)*SPA(5,6)*SPB(6,2)))/(complex<T>(2,0)*SPA(3,4)*SPA(5,6)*
 SPB(2,1)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3)))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2))*(SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*
 SPB(6,4)*((-S(1,2)+S(3,4)-S(5,6))*SPA(3,4)*SPB(6,4)-
  (-S(1,2)-S(3,4)+S(5,6))*SPA(3,5)*SPB(6,5)))/
(complex<T>(2,0)*SPA(3,4)*SPB(2,1)*(SPA(1,3)*SPB(4,1)+
  SPA(2,3)*SPB(4,2))*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*
 SPB(6,5))+(complex<T>(1,0)*SPA(2,4)*SPB(6,4)*
 ((-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SPB(6,4)-
  SPA(3,5)*SPB(3,1)*SPB(6,5)))/(complex<T>(2,0)*SPA(3,4)*SPB(2,1)*
 (SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SPB(6,5))-
 (pow(SPB(3,1),2)*pow(SPA(1,3)*SPB(6,1)+SPA(2,3)*SPB(6,2),2)*
 pow(complex<T>(1,0)-SS(1,2,3)/S(1,2),-1))/(complex<T>(2,0)*pow(SPB(2,1),2)*
 SPA(1,3)*(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*SPB(6,5)*
 SS(1,2,3))+(complex<T>(1,0)*pow(SPB(3,1),2)*
 pow(-(SPA(1,2)*SPB(6,2))-SPA(1,3)*SPB(6,3),2)*
 pow(complex<T>(1,0)-SS(1,2,3)/S(2,3),-1)*SPA(1,2))/
(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(1,3)*SPA(2,3)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SPB(6,5)*
 SS(1,2,3))+(complex<T>(1,0)*pow(SPB(3,1),2)*
 pow(-(SPA(1,2)*SPB(6,2))-SPA(1,3)*SPB(6,3),2))/
(complex<T>(2,0)*SPA(1,3)*SPB(2,1)*SPB(3,2)*(-(SPA(1,2)*SPB(4,2))-
  SPA(1,3)*SPB(4,3))*SPB(6,5)*SS(1,2,3))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(6,4),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)*(SPA(1,2)*SPB(4,1)-
  SPA(2,3)*SPB(4,3)))/(complex<T>(2,0)*SPA(2,3)*(SPA(1,3)*SPB(4,1)+
  SPA(2,3)*SPB(4,2))*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*
 SPB(6,5)*SS(1,2,3))-(pow(SPA(2,4),2)*
 pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(3,4),-1))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 SPA(5,6)*SPB(4,2)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
 SS(2,3,4))+(complex<T>(1,0)*pow(SPA(2,4),2)*
 pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2))/
(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*SPA(5,6)*SPB(4,2)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+
   SPA(3,5)*SPB(4,3),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*
 SPB(4,3))/(complex<T>(2,0)*pow(SPA(2,3),2)*SPA(5,6)*SPB(3,2)*SPB(4,2)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(1,5),2)*pow(SPB(3,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)*(-(SPA(1,2)*SPB(3,2))+
  SPA(1,4)*SPB(4,3)))/(complex<T>(2,0)*SPA(5,6)*SPB(3,2)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SS(2,3,4)))+
 (complex<T>(0,-1)*(-((SPA(2,4)*SPA(4,5)*SPB(3,1)*SPB(6,1))/
  (S(2,3)*S(5,6)*SPA(3,4)*SPB(2,1)))+
(pow(SPB(3,1),2)*SPA(4,5)*(SPA(1,2)*SPB(6,1)-
   SPA(2,3)*SPB(6,3)))/(S(2,3)*S(5,6)*SPB(2,1)*SS(1,2,3))-
(pow(SPA(2,4),2)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*
  SPB(6,1))/(S(2,3)*S(5,6)*SPA(3,4)*SS(2,3,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpmmqmemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, m, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpmmqmemep L");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(6,1),2))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*SPB(4,3)*
 SPB(6,5))+complex<T>(0,1)*
(-((pow(SPB(6,1),2)*(-((pow(SPA(2,3),2)*pow(SPB(2,1),2)*
pow(SPB(6,3),2)*pow(complex<T>(1,0)-S(1,2)/SS(1,2,3),-1))/
 (complex<T>(2,0)*pow(SPB(6,1),2)*pow(SS(1,2,3),2)))+
   (complex<T>(1,0)*pow(-(SPA(2,4)*SPB(2,1)*SPB(6,4))-SPA(3,4)*
  SPB(3,1)*SPB(6,4),2)*pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),
-1))/(complex<T>(2,0)*pow(SPB(6,1),2)*pow(SS(1,2,3),2))))/
 (SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(6,5)))+
 (complex<T>(1,0)*(-((SPA(4,5)*(-(S(1,2)*SPB(6,1))+SPA(2,3)*SPB(2,1)*SPB(
  6,3)))/SS(1,2,3))-((S(2,4)*SPA(4,5)+
 SPA(2,4)*SPA(3,5)*SPB(3,2))*SPB(6,1))/SS(2,3,4)))/
(complex<T>(3,0)*pow(SPB(3,2),2)*S(5,6))+
 (complex<T>(1,0)*(-(pow(-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,1),2)/
(SPA(5,6)*SPB(2,1)*SPB(3,2)*SPB(4,3)*SS(1,2,3)))+
  (SPA(2,3)*SPB(6,1)*SPB(6,3))/(SPB(3,2)*SPB(4,3)*SPB(6,5)*
SS(1,2,3))+(SPA(2,3)*SPA(4,5)*SPB(3,1)*
(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3)))/(SPA(5,6)*SPB(3,2)*
SPB(4,3)*SS(1,2,3)*SS(2,3,4))-
  (SPA(2,3)*SPB(6,1)*(SPA(1,2)*SPB(6,1)-SPA(2,3)*SPB(6,3)))/
   (SPB(4,3)*SPB(6,5)*SS(1,2,3)*SS(2,3,4))))/complex<T>(2,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qmppqpemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qmppqpemep L");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,5),2))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*
 SPA(5,6))+complex<T>(0,-1)*
((pow(SPA(1,5),2)*(-((pow(SPA(1,2),2)*pow(SPA(3,5),2)*
 pow(SPB(3,2),2)*pow(complex<T>(1,0)-S(1,2)/SS(1,2,3),-1))/
(complex<T>(2,0)*pow(SPA(1,5),2)*pow(SS(1,2,3),2)))+
  (complex<T>(1,0)*pow(-(SPA(1,2)*SPA(4,5)*SPB(4,2))-SPA(1,3)*SPA(4,
  5)*SPB(4,3),2)*pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1))/
   (complex<T>(2,0)*pow(SPA(1,5),2)*pow(SS(1,2,3),2))))/
(SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(5,6))+
 (complex<T>(1,0)*(-(((-(S(1,2)*SPA(1,5))+SPA(1,2)*SPA(3,5)*SPB(3,2))*
 SPB(6,4))/SS(1,2,3))-
  (SPA(1,5)*(SPA(2,3)*SPB(4,2)*SPB(6,3)+S(2,4)*SPB(6,4)))/
   SS(2,3,4)))/(complex<T>(3,0)*pow(SPA(2,3),2)*S(5,6))+
 (complex<T>(1,0)*(-((SPA(1,5)*SPA(3,5)*SPB(3,2))/(SPA(2,3)*SPA(3,4)*
 SPA(5,6)*SS(1,2,3)))+pow(-(SPA(1,2)*SPB(6,2))-
 SPA(1,3)*SPB(6,3),2)/(SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(6,5)*
SS(1,2,3))+(SPA(1,5)*SPB(3,2)*(SPA(1,5)*SPB(2,1)-
 SPA(3,5)*SPB(3,2)))/(SPA(3,4)*SPA(5,6)*SS(1,2,3)*
SS(2,3,4))-(SPA(1,3)*SPB(3,2)*(SPA(2,4)*SPB(6,2)+
 SPA(3,4)*SPB(6,3))*SPB(6,4))/(SPA(2,3)*SPA(3,4)*SPB(6,5)*
SS(1,2,3)*SS(2,3,4))))/complex<T>(2,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qmpmqpemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qmpmqpemep L");
#endif
 
 return( (complex<T>(0,1)*(-((SPA(1,3)*SPA(1,5)*SPB(4,2)*SPB(6,4))/
  (S(2,3)*S(5,6)*SPA(1,2)*SPB(4,3)))+
(pow(SPA(1,3),2)*(SPA(1,5)*SPB(2,1)-SPA(3,5)*SPB(3,2))*
  SPB(6,4))/(S(2,3)*S(5,6)*SPA(1,2)*SS(1,2,3))-
(pow(SPB(4,2),2)*SPA(1,5)*(SPA(2,3)*SPB(6,2)-
   SPA(3,4)*SPB(6,4)))/(S(2,3)*S(5,6)*SPB(4,3)*SS(2,3,4))))/
complex<T>(2,0)+complex<T>(0,-1)*((complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
   pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
   complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(4,5)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*(SPA(1,3)*SPB(4,1)+
  SPA(2,3)*SPB(4,2))*((-S(1,2)+S(3,4)-S(5,6))*SPA(4,5)*
   SPB(4,3)+(-S(1,2)-S(3,4)+S(5,6))*SPA(5,6)*SPB(6,3)))/
(complex<T>(2,0)*SPA(1,2)*SPA(5,6)*(-(SPA(2,4)*SPB(2,1))-
  SPA(3,4)*SPB(3,1))*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*
 SPB(4,3))-(SPA(4,5)*SPB(4,2)*
 (SPA(4,5)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))+
  SPA(1,3)*SPA(5,6)*SPB(6,3)))/(complex<T>(2,0)*SPA(1,2)*SPA(5,6)*
 (-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*
 (SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SPB(4,3))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2))*(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*
 SPB(6,1)*(-((S(1,2)-S(3,4)-S(5,6))*SPA(1,2)*SPB(6,1))-
  (-S(1,2)-S(3,4)+S(5,6))*SPA(2,5)*SPB(6,5)))/
(complex<T>(2,0)*SPA(1,2)*(-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*
 (SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(4,3)*SPB(6,5))+
 (complex<T>(1,0)*SPA(1,3)*SPB(6,1)*
 (-((-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SPB(6,1))-
  SPA(2,5)*SPB(4,2)*SPB(6,5)))/(complex<T>(2,0)*SPA(1,2)*
 (-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*
 (SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(4,3)*SPB(6,5))-
 (pow(SPA(1,3),2)*pow(-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,1),
  2))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(5,6)*SPB(3,1)*
 (-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SS(1,2,3))-
 (pow(SPA(1,3),2)*pow(-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,1),
  2)*pow(complex<T>(1,0)-SS(1,2,3)/S(2,3),-1)*SPB(2,1))/
(complex<T>(2,0)*pow(SPA(2,3),2)*SPA(5,6)*SPB(3,1)*
 (-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SPB(3,2)*
 SS(1,2,3))+(complex<T>(1,0)*pow(SPA(1,3),2)*
 pow(SPA(1,5)*SPB(3,1)+SPA(2,5)*SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(1,2,3)/S(1,2),-1))/(complex<T>(2,0)*pow(SPA(1,2),2)*
 SPA(5,6)*SPB(3,1)*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*
 SS(1,2,3))-(pow(SPA(4,5),2)*pow(SPB(4,2),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)*(SPA(1,4)*SPB(2,1)-
  SPA(3,4)*SPB(3,2)))/(complex<T>(2,0)*SPA(5,6)*(-(SPA(2,4)*SPB(2,1))-
  SPA(3,4)*SPB(3,1))*SPB(3,2)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2))*SS(1,2,3))-
 (pow(SPB(4,2),2)*pow(SPA(2,4)*SPB(6,2)+SPA(3,4)*SPB(6,3),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*SPA(3,4))/
(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,3)*SPA(2,4)*
 (-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SPB(6,5)*
 SS(2,3,4))+(complex<T>(1,0)*pow(SPB(4,2),2)*
 pow(-(SPA(2,3)*SPB(6,3))-SPA(2,4)*SPB(6,4),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(3,4),-1))/(complex<T>(2,0)*pow(SPB(4,3),2)*
 SPA(2,4)*(SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(6,5)*
 SS(2,3,4))-(pow(SPA(1,3),2)*pow(SPB(6,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)*(-(SPA(2,3)*SPB(2,1))+
  SPA(3,4)*SPB(4,1)))/(complex<T>(2,0)*SPA(2,3)*(-(SPA(2,4)*SPB(2,1))-
  SPA(3,4)*SPB(3,1))*(SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*
 SPB(6,5)*SS(2,3,4))-(pow(SPB(4,2),2)*
 pow(SPA(2,4)*SPB(6,2)+SPA(3,4)*SPB(6,3),2))/
(complex<T>(2,0)*SPA(2,4)*(-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*
 SPB(3,2)*SPB(4,3)*SPB(6,5)*SS(2,3,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qmmpqpemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qmmpqpemep L");
#endif
 
 return( (complex<T>(0,1)*(((SPA(1,5)*SPB(3,1)+SPA(2,5)*SPB(3,2))*
  (-(SPA(2,3)*SPB(6,3))-SPA(2,4)*SPB(6,4)))/
 (S(2,3)*S(5,6)*SPA(3,4)*SPB(2,1))+
(SPA(1,2)*SPB(3,1)*(SPA(1,5)*SPB(3,1)+SPA(2,5)*SPB(3,2))*
  SPB(6,4))/(S(2,3)*S(5,6)*SPB(2,1)*SS(1,2,3))-
(SPA(1,5)*SPA(2,4)*SPB(4,3)*(-(SPA(2,3)*SPB(6,3))-
   SPA(2,4)*SPB(6,4)))/(S(2,3)*S(5,6)*SPA(3,4)*SS(2,3,4))))/
complex<T>(2,0)+complex<T>(0,-1)*((complex<T>(1,0)*pow(SPA(2,5),2))/(complex<T>(2,0)*SPA(2,3)*
 SPA(3,4)*SPA(5,6)*SPB(2,1))+(complex<T>(1,0)*pow(SPB(6,3),2))/
(complex<T>(2,0)*SPA(3,4)*SPB(2,1)*SPB(3,2)*SPB(6,5))-
 ((pow(SPA(4,5),2)*pow(SPB(4,3),2)*
   pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)+
  (pow(SPA(1,2),2)*pow(-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,
  1),2)*pow(complex<T>(1,0)-SS(1,2,3)/S(2,3),-1))/
   pow(SPA(2,3),2))*SPB(3,1))/(complex<T>(2,0)*SPA(5,6)*SPB(2,1)*
 (-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SPB(3,2)*
 SS(1,2,3))-((pow(SPA(1,2),2)*pow(SPB(6,1),2)*
   pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)+
  (pow(SPB(4,3),2)*pow(SPA(2,4)*SPB(6,2)+SPA(3,4)*SPB(6,3),
 2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/pow(SPB(3,2),2))*
 SPA(2,4))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*(-(SPA(2,4)*SPB(2,1))-
  SPA(3,4)*SPB(3,1))*SPB(6,5)*SS(2,3,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qmmmqpemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, qp, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qmmmqpemep L");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(6,4),2))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*SPB(4,3)*
 SPB(6,5))+complex<T>(0,-1)*
((pow(SPB(6,4),2)*(-((pow(SPA(2,3),2)*pow(SPB(4,3),2)*
 pow(SPB(6,2),2)*pow(complex<T>(1,0)-S(3,4)/SS(2,3,4),-1))/
(complex<T>(2,0)*pow(SPB(6,4),2)*pow(SS(2,3,4),2)))+
  (complex<T>(1,0)*pow(-(SPA(1,2)*SPB(4,2)*SPB(6,1))-SPA(1,3)*SPB(4,
  3)*SPB(6,1),2)*pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1))/
   (complex<T>(2,0)*pow(SPB(6,4),2)*pow(SS(2,3,4),2))))/
(SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(6,5))+
 (complex<T>(1,0)*(-(((S(1,3)*SPA(1,5)+SPA(1,3)*SPA(2,5)*SPB(3,2))*
 SPB(6,4))/SS(1,2,3))-
  (SPA(1,5)*(SPA(2,3)*SPB(4,3)*SPB(6,2)-S(3,4)*SPB(6,4)))/
   SS(2,3,4)))/(complex<T>(3,0)*pow(SPB(3,2),2)*S(5,6))+
 (complex<T>(1,0)*(pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)/
   (SPA(5,6)*SPB(2,1)*SPB(3,2)*SPB(4,3)*SS(2,3,4))-
  (SPA(2,3)*SPB(6,2)*SPB(6,4))/(SPB(2,1)*SPB(3,2)*SPB(6,5)*
SS(2,3,4))+(SPA(1,5)*SPA(2,3)*(-(SPA(2,5)*SPB(2,1))-
 SPA(3,5)*SPB(3,1))*SPB(4,2))/(SPA(5,6)*SPB(2,1)*SPB(3,2)*
SS(1,2,3)*SS(2,3,4))-(SPA(2,3)*SPB(6,4)*
(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4)))/(SPB(2,1)*SPB(6,5)*
SS(1,2,3)*SS(2,3,4))))/complex<T>(2,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qppqmpemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, qm, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qppqmpemep SLC");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(3,5),2)*SPA(1,3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*
 SPA(2,3)*SPA(3,4)*SPA(5,6))+
 complex<T>(0,1)*(-((pow(SPA(4,5),2)*pow(SPB(4,2),2)*
  pow(complex<T>(1,0)-S(2,3)/SS(2,3,4),-1)*SPA(2,3))/
 (complex<T>(2,0)*pow(SS(2,3,4),2)*SPA(1,4)*SPA(2,4)*SPA(5,6)))+
 (complex<T>(1,0)*pow(SPA(1,3),3)*pow(SPB(6,1),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(5,6),-1))/(complex<T>(2,0)*pow(SPB(6,5),2)*
 SPA(1,2)*SPA(1,4)*SPA(2,3)*SPA(3,4)*SPA(5,6))-
 (pow(SPA(3,5),2)*SPA(1,3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*SPA(2,3)*
 SPA(3,4)*SPA(5,6))+(pow(SPB(6,4),2)*
 pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)*SPA(1,3)*SPA(3,4))/
(pow(SPB(6,5),2)*SPA(1,2)*SPA(1,4)*SPA(2,3)*SPA(5,6))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(4,2),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(2,3,4),-1)*SPA(3,4))/
(complex<T>(2,0)*pow(SS(2,3,4),2)*SPA(1,2)*SPA(2,4)*SPA(5,6))-
 (pow(SPB(4,2),2)*SPA(1,5)*SPA(2,4)*SPA(3,5))/
(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*SPA(5,6)*SS(1,2,4)*SS(2,3,4))-
 (SPB(4,2)*(-(SPA(1,2)*SPB(6,2))-SPA(1,4)*SPB(6,4))*
 (SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4)))/(complex<T>(2,0)*SPA(1,2)*
 SPA(1,4)*SPB(6,5)*SS(1,2,4)*SS(2,3,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qppqmmemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, qm, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qppqmmemep SLC");
#endif
 
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*SPA(3,4)*
  SPA(4,5)*SPB(4,2)*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3)))/
 (complex<T>(2,0)*pow(SPA(2,3),2)*SPA(5,6)*SPB(3,2)*SPB(4,3)*
  (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))))-
 (pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)*SPA(1,5)*SPB(2,1)*
 (-(SPA(1,3)*SPB(3,2)*SPB(6,1))-SPA(1,4)*SPB(4,2)*SPB(6,1)))/
(complex<T>(2,0)*pow(SS(2,3,4),2)*SPB(4,3)*(-(SPA(1,2)*SPB(4,2))-
  SPA(1,3)*SPB(4,3)))+
 (SPA(1,3)*(-((pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+
 pow(SPA(1,4),2)*pow(SPB(4,1),2)+pow(SPA(5,6),2)*
  pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*
  S(1,4)*S(5,6)+complex<T>(-2,0)*S(2,3)*S(5,6),-1)*SPA(4,5)*
 (-((-S(1,4)-S(2,3)+S(5,6))*SPA(1,3))+complex<T>(2,0)*SPA(1,
  4)*SPA(2,3)*SPB(4,2))*SPB(6,1))/(SPA(1,2)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))))-
  (pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+pow(SPA(1,4),2)*pow(
  SPB(4,1),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*S(1,4)*S(5,6)+
complex<T>(-2,0)*S(2,3)*S(5,6),-1)*(-((S(1,4)-S(2,3)-S(5,6))*
 SPA(1,3)*SPB(4,1))+(-S(1,4)+S(2,3)-S(5,6))*
SPA(2,3)*SPB(4,2))*((pow(SPA(4,5),2)*SPB(4,1))/SPA(5,6)+
 (pow(SPB(6,1),2)*SPA(1,4))/SPB(6,5)))/(complex<T>(2,0)*SPA(1,2)*
SPB(4,1)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3)))+
  (pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)*SPA(3,4)*SPA(4,5)*
SPB(6,4))/(SPA(1,2)*(-(SPA(1,2)*SPB(4,2))-
 SPA(1,3)*SPB(4,3))*SS(1,2,3))))/SPA(2,3)+
 (SPA(1,4)*(-((pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+
 pow(SPA(1,4),2)*pow(SPB(4,1),2)+pow(SPA(5,6),2)*
  pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*
  S(1,4)*S(5,6)+complex<T>(-2,0)*S(2,3)*S(5,6),-1)*SPA(2,4)*
 ((SPA(4,5)*(-((SPA(3,5)*(-((S(1,4)-S(2,3)-S(5,6))*
     SPA(1,4)*SPB(3,1))+(-S(1,4)+S(2,3)-S(5,6))*
    SPA(2,4)*SPB(3,2)))/SPA(5,6))-
   ((-S(1,4)-S(2,3)+S(5,6))*SPA(2,4)+complex<T>(-2,0)*
  SPA(1,4)*SPA(2,3)*SPB(3,1))*SPB(6,2)))/SPA(1,4)-
(SPB(6,1)*(-(SPA(3,5)*((-S(1,4)-S(2,3)+S(5,6))*
   SPB(3,1)+complex<T>(-2,0)*SPA(2,4)*SPB(3,2)*SPB(4,1)))-
   (((-S(1,4)+S(2,3)-S(5,6))*SPA(2,3)*SPB(3,1)-
  (S(1,4)-S(2,3)-S(5,6))*SPA(2,4)*SPB(4,1))*
 SPB(6,2))/SPB(6,5)))/SPB(4,1)))/(complex<T>(2,0)*SPA(1,2)*
 (SPA(1,2)*SPB(3,1)+SPA(2,4)*SPB(4,3))))+
  (complex<T>(1,0)*pow(SPA(2,4),2)*pow(complex<T>(1,0)-S(1,4)/SS(1,2,4),-1)*
SPB(2,1)*SPB(6,2)*(SPA(1,2)*SPB(6,1)-SPA(2,4)*SPB(6,4)))/
   (complex<T>(2,0)*pow(SS(1,2,4),2)*SPA(1,2)*(SPA(1,2)*SPB(3,1)+
 SPA(2,4)*SPB(4,3))*SPB(6,5))-(pow(SPA(3,4),2)*
pow(SPB(6,3),2)*pow(complex<T>(1,0)-S(5,6)/SS(1,2,4),-1)*
SPA(2,4))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*(SPA(1,2)*SPB(3,1)+
 SPA(2,4)*SPB(4,3))*SPB(6,5)*SS(1,2,4))))/SPA(2,4)+
 (complex<T>(1,0)*(-((SPA(1,3)*SPA(3,5)*SPA(4,5))/(SPA(1,2)*SPA(2,3)*
 SPA(5,6)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))))-
  (SPA(1,5)*SPA(3,5)*SPB(2,1))/(SPA(2,3)*SPA(5,6)*SPB(4,3)*
(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3)))-
  (pow(SPA(3,5),2)*pow(SPB(3,1),2))/(SPA(2,3)*SPA(5,6)*
SPB(4,1)*SPB(4,3)*(SPA(1,2)*SPB(3,1)+SPA(2,4)*
SPB(4,3)))-(S(3,4)*SPA(1,5)*SPA(4,5))/
   (SPA(1,2)*SPA(5,6)*(-(SPA(1,2)*SPB(4,2))-
 SPA(1,3)*SPB(4,3))*(SPA(1,2)*SPB(3,1)+
 SPA(2,4)*SPB(4,3)))-(SPA(1,5)*SPA(3,4)*SPA(3,5)*
SPB(3,1))/(SPA(2,3)*SPA(5,6)*(-(SPA(1,2)*SPB(4,2))-
 SPA(1,3)*SPB(4,3))*(SPA(1,2)*SPB(3,1)+
 SPA(2,4)*SPB(4,3)))+(pow(SPB(6,1),2)*SPA(1,3))/
   (SPA(1,2)*SPA(2,3)*SPB(4,1)*SPB(4,3)*SPB(6,5))-
  (SPA(1,3)*SPB(2,1)*SPB(6,1)*SPB(6,4))/(SPA(2,3)*SPB(4,1)*
SPB(4,3)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*
SPB(6,5))-(SPA(1,4)*SPA(3,4)*SPB(6,3)*SPB(6,4))/
   (SPA(1,2)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*
(SPA(1,2)*SPB(3,1)+SPA(2,4)*SPB(4,3))*SPB(6,5))-
  (S(1,3)*SPA(3,4)*SPB(6,3)*SPB(6,4))/(SPA(2,3)*SPB(4,3)*
(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*
(SPA(1,2)*SPB(3,1)+SPA(2,4)*SPB(4,3))*SPB(6,5))-
  (SPA(3,4)*(-(SPA(1,2)*SPB(6,2))-SPA(1,3)*SPB(6,3))*
SPB(6,4))/(SPA(1,2)*SPA(2,3)*SPB(4,3)*
(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SPB(6,5))-
  (SPA(2,4)*SPB(6,2)*(SPA(1,4)*SPB(6,1)+SPA(2,4)*SPB(6,2)))/
   (SPA(1,2)*(SPA(1,2)*SPB(3,1)+SPA(2,4)*SPB(4,3))*SPB(6,5)*
SS(1,2,4))-(SPA(1,5)*SPB(2,1)*(-(SPA(3,5)*SPB(3,2))-
 SPA(4,5)*SPB(4,2)))/(SPA(5,6)*SPB(4,3)*
(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SS(2,3,4))))/
complex<T>(2,0))+(complex<T>(0,-1)*(((-(SPA(3,5)*SPB(3,1))-SPA(4,5)*SPB(4,1))*
  (SPA(1,3)*SPB(6,1)+SPA(2,3)*SPB(6,2)))/(S(5,6)*SPA(1,2)*
  SPA(2,3)*SPB(4,1)*SPB(4,3))+(SPA(3,5)*SPB(2,1)*
  (SPA(1,4)*SPB(6,1)+SPA(2,4)*SPB(6,2)))/(S(5,6)*SPA(1,2)*
  SPB(4,1)*SS(1,2,4))-(SPA(3,4)*(-(SPA(3,5)*SPB(3,2))-
   SPA(4,5)*SPB(4,2))*SPB(6,1))/(S(5,6)*SPA(2,3)*SPB(4,3)*
  SS(2,3,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpmqmpemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, qm, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpmqmpemep SLC");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(5,6)/SS(1,2,4),-1)*SPA(2,3)*
 (-(SPA(1,2)*SPA(3,5)*SPB(3,1))-SPA(2,4)*SPA(3,5)*SPB(4,3))*
 SPB(6,3))/(complex<T>(2,0)*pow(SS(1,2,4),2)*SPA(1,4)*
 (SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2)))-
 (pow(complex<T>(1,0)-SS(1,2,4)/S(1,2),-1)*SPA(2,4)*SPB(4,1)*
 (SPA(1,4)*SPB(6,1)+SPA(2,4)*SPB(6,2))*SPB(6,4))/
(complex<T>(2,0)*pow(SPB(2,1),2)*SPA(1,2)*SPA(1,4)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2))*SPB(6,5))+
 (SPB(3,1)*((pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(
  SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(
  5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(3,5)*
((-S(1,2)-S(3,4)+S(5,6))*SPB(3,1)+complex<T>(-2,0)*SPA(2,4)*
SPB(2,1)*SPB(4,3))*SPB(6,4))/(SPB(3,2)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2)))+
  (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),
  2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-((S(1,2)-S(3,4)-S(5,6))*
 SPA(2,4)*SPB(2,1))+(-S(1,2)+S(3,4)-S(5,6))*
SPA(3,4)*SPB(3,1))*(-((pow(SPA(3,5),2)*SPB(4,3))/SPA(5,
  6))-(pow(SPB(6,4),2)*SPA(3,4))/SPB(6,5)))/
   (complex<T>(2,0)*SPA(3,4)*SPB(3,2)*(SPA(1,4)*SPB(3,1)+
 SPA(2,4)*SPB(3,2)))-(pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)*
SPA(4,5)*SPB(4,1)*SPB(6,4))/(SPB(3,2)*(SPA(1,4)*SPB(3,1)+
 SPA(2,4)*SPB(3,2))*SS(1,2,3))))/SPB(2,1)+
 (SPB(4,3)*(-((pow(SPB(4,2),2)*pow(complex<T>(1,0)-S(3,4)/SS(2,3,4),
-1)*SPA(2,3)*SPA(2,5)*(-(SPA(3,5)*SPB(3,2))-
SPA(4,5)*SPB(4,2)))/(complex<T>(2,0)*pow(SS(2,3,4),2)*SPA(5,6)*
 SPB(3,2)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))))+
  (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),
  2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPB(4,2)*
(-((SPA(3,5)*((SPA(2,5)*((S(1,2)-S(3,4)-S(5,6))*
   SPA(1,3)*SPB(2,1)-(-S(1,2)+S(3,4)-S(5,6))*
   SPA(3,4)*SPB(4,2)))/SPA(5,6)-
   (-((-S(1,2)-S(3,4)+S(5,6))*SPA(1,3))+complex<T>(2,0)*
  SPA(1,2)*SPA(3,4)*SPB(4,2))*SPB(6,1)))/SPA(3,4))+
 (SPB(6,4)*(-(SPA(2,5)*((-S(1,2)-S(3,4)+S(5,6))*
  SPB(4,2)+complex<T>(-2,0)*SPA(1,3)*SPB(2,1)*SPB(4,3)))+
  ((-((S(1,2)-S(3,4)-S(5,6))*SPA(1,2)*SPB(4,2))+
 (-S(1,2)+S(3,4)-S(5,6))*SPA(1,3)*SPB(4,3))*
SPB(6,1))/SPB(6,5)))/SPB(4,3)))/(complex<T>(2,0)*SPB(3,2)*
(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2)))-
  (pow(SPA(1,5),2)*pow(SPB(4,1),2)*pow(complex<T>(1,0)-S(5,6)/SS(2,3,
  4),-1)*SPB(4,2))/(complex<T>(2,0)*SPA(5,6)*SPB(3,2)*
(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(4,3)*
SS(2,3,4))))/SPB(4,2)+
 (complex<T>(1,0)*((pow(SPA(3,5),2)*SPB(3,1))/(SPA(1,4)*SPA(3,4)*SPA(5,6)*
SPB(2,1)*SPB(3,2))-(SPA(2,3)*SPA(3,5)*SPA(4,5)*SPB(3,1))/
   (SPA(1,4)*SPA(3,4)*SPA(5,6)*SPB(2,1)*(SPA(1,4)*SPB(3,1)+
 SPA(2,4)*SPB(3,2)))+(SPA(4,5)*(SPA(1,5)*SPB(3,1)+
 SPA(2,5)*SPB(3,2))*SPB(4,1))/(SPA(1,4)*SPA(5,6)*SPB(2,1)*
SPB(3,2)*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2)))-
  (S(1,3)*SPA(1,5)*SPA(4,5)*SPB(4,1))/(SPA(1,4)*SPA(5,6)*
SPB(2,1)*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*
(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2)))-
  (SPA(1,5)*SPA(4,5)*SPB(4,1)*SPB(4,3))/(SPA(5,6)*SPB(3,2)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*
(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2)))-
  (pow(SPA(1,3),2)*pow(SPB(6,1),2))/(SPA(1,4)*SPA(3,4)*
SPB(2,1)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
SPB(6,5))+(SPA(2,3)*SPB(6,1)*SPB(6,3))/(SPA(1,4)*SPB(2,1)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SPB(6,5))+
  (SPA(1,3)*SPB(4,1)*SPB(6,1)*SPB(6,3))/(SPB(2,1)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*
(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(6,5))-
  (SPB(3,1)*SPB(6,1)*SPB(6,4))/(SPB(2,1)*SPB(3,2)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SPB(6,5))-
  (S(1,4)*SPB(6,3)*SPB(6,4))/(SPB(3,2)*(SPA(1,4)*SPB(3,1)+
 SPA(2,4)*SPB(3,2))*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
SPB(6,5))-(SPA(2,3)*SPB(6,3)*(SPA(1,2)*SPB(6,1)-
 SPA(2,4)*SPB(6,4)))/(SPA(1,4)*(SPA(1,4)*SPB(3,1)+
 SPA(2,4)*SPB(3,2))*SPB(6,5)*SS(1,2,4))-
  (SPA(2,5)*SPB(4,2)*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3)))/
   (SPA(5,6)*SPB(3,2)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
SS(2,3,4))))/complex<T>(2,0))+
 (complex<T>(0,-1)*(((-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,1))*
  (SPA(1,3)*SPB(6,1)-SPA(3,4)*SPB(6,4)))/(S(5,6)*SPA(1,4)*
  SPA(3,4)*SPB(2,1)*SPB(3,2))+(SPA(3,5)*SPB(4,1)*
  (SPA(1,2)*SPB(6,1)-SPA(2,4)*SPB(6,4)))/(S(5,6)*SPA(1,4)*
  SPB(2,1)*SS(1,2,4))+(SPA(2,3)*(SPA(2,5)*SPB(4,2)+
   SPA(3,5)*SPB(4,3))*SPB(6,1))/(S(5,6)*SPA(3,4)*SPB(3,2)*
  SS(2,3,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpmqmmemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, qm, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpmqmmemep SLC");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(6,1),2)*SPB(3,1))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*
 SPB(4,1)*SPB(4,3)*SPB(6,5))+
 complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(6,2),2)*
 pow(complex<T>(1,0)-S(1,4)/SS(1,2,4),-1)*SPB(4,1))/
(complex<T>(2,0)*pow(SS(1,2,4),2)*SPB(3,2)*SPB(4,2)*SPB(6,5))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SPB(3,1),3)*
 pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1))/(complex<T>(2,0)*pow(SPA(5,6),2)*
 SPB(2,1)*SPB(3,2)*SPB(4,1)*SPB(4,3)*SPB(6,5))-
 (pow(SPB(6,1),2)*SPB(3,1))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*SPB(4,1)*
 SPB(4,3)*SPB(6,5))+(pow(SPA(4,5),2)*
 pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)*SPB(3,1)*SPB(4,1))/
(pow(SPA(5,6),2)*SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(6,5))-
 (pow(SPA(2,4),2)*pow(SPB(6,4),2)*
 pow(complex<T>(1,0)-S(1,2)/SS(1,2,4),-1)*SPB(2,1))/
(complex<T>(2,0)*pow(SS(1,2,4),2)*SPB(4,2)*SPB(4,3)*SPB(6,5))-
 (SPA(2,4)*(-(SPA(2,5)*SPB(2,1))-SPA(4,5)*SPB(4,1))*
 (SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3)))/(complex<T>(2,0)*SPA(5,6)*
 SPB(3,2)*SPB(4,3)*SS(1,2,4)*SS(2,3,4))-
 (pow(SPA(2,4),2)*SPB(4,2)*SPB(6,1)*SPB(6,3))/
(complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(6,5)*SS(1,2,4)*SS(2,3,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmppemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmppemep SLC");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(2,5),2))/(complex<T>(2,0)*SPA(1,4)*SPA(2,3)*SPA(3,4)*
 SPA(5,6))+complex<T>(0,1)*((pow(SPA(2,4),2)*pow(SPB(6,4),2)*
 pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1))/(pow(SPB(6,5),2)*SPA(1,4)*
 SPA(2,3)*SPA(3,4)*SPA(5,6))+(complex<T>(1,0)*pow(SPA(1,2),2)*
 pow(SPB(6,1),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(5,6),-1))/
(complex<T>(2,0)*pow(SPB(6,5),2)*SPA(1,4)*SPA(2,3)*SPA(3,4)*SPA(5,6))+
 (pow(SPB(6,3),2)*pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1)*SPA(2,3))/
(pow(SPB(6,5),2)*SPA(1,4)*SPA(3,4)*SPA(5,6))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*
 pow(complex<T>(1,0)-S(2,3)/SS(2,3,4),-1)*SPA(2,3))/
(complex<T>(2,0)*pow(SS(2,3,4),2)*SPA(1,4)*SPA(3,4)*SPA(5,6))+
 (complex<T>(1,0)*(-(pow(SPA(2,5),2)/(SPA(1,4)*SPA(2,3)*SPA(3,4)*
 SPA(5,6)))-(SPA(2,5)*SPA(4,5)*SPB(4,3))/
   (SPA(1,4)*SPA(3,4)*SPA(5,6)*SS(2,3,4))-
  (SPA(1,5)*SPA(2,5)*SPB(3,1)*SPB(4,3))/(SPA(1,4)*SPA(5,6)*
SS(1,3,4)*SS(2,3,4))-(SPB(4,3)*(SPA(1,4)*SPB(6,1)+
 SPA(3,4)*SPB(6,3))*(-(SPA(2,3)*SPB(6,3))-
 SPA(2,4)*SPB(6,4)))/(SPA(1,4)*SPA(3,4)*SPB(6,5)*SS(1,3,4)*
SS(2,3,4))))/complex<T>(2,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmpmemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmpmemep SLC");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*SPA(1,5)*SPA(4,5)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2)))/(complex<T>(2,0)*SPA(3,4)*SPA(5,6)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3)))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(3,5)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2))*(-((SPA(3,5)*SPB(2,1)*
 ((-S(1,2)-S(3,4)+S(5,6))*SPA(2,4)+complex<T>(-2,0)*SPA(1,
  2)*SPA(3,4)*SPB(3,1)))/SPA(3,4))+
  ((SPA(1,5)*SPB(4,1)+SPA(2,5)*SPB(4,2))*
(-((-S(1,2)-S(3,4)+S(5,6))*SPB(3,1))+complex<T>(2,0)*SPA(2,4)*
SPB(2,1)*SPB(4,3)))/SPB(4,3)))/(complex<T>(2,0)*SPA(5,6)*
 (SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*(SPA(1,3)*SPB(2,1)+
  SPA(3,4)*SPB(4,2)))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(1,5)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2))*(SPB(3,2)*((-S(1,2)-S(3,4)+S(5,6))*
 SPA(2,5)+complex<T>(2,0)*SPA(1,2)*(-(SPA(3,5)*SPB(3,1))-
SPA(4,5)*SPB(4,1)))+
  (SPA(1,4)*(-((S(1,2)-S(3,4)-S(5,6))*SPA(2,5)*SPB(2,1))-
 (-S(1,2)-S(3,4)+S(5,6))*SPA(5,6)*SPB(6,1)))/SPA(3,4)))/
(complex<T>(2,0)*SPA(5,6)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3)))+
 (complex<T>(1,0)*pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+pow(SPA(1,4),2)*
pow(SPB(4,1),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*S(1,4)*S(5,6)+
   complex<T>(-2,0)*S(2,3)*S(5,6),-1)*(SPA(2,3)*SPB(3,1)+
  SPA(2,4)*SPB(4,1))*(((S(1,4)-S(2,3)-S(5,6))*SPA(2,5)*
SPA(4,5))/(SPA(2,3)*SPA(5,6))+
  (complex<T>(1,0)*(-S(1,4)-S(2,3)+S(5,6))*SPA(2,5)*SPB(6,1))/
   (complex<T>(2,0)*SPA(2,3)*SPB(4,1))-SPA(4,5)*SPB(6,3)))/
(complex<T>(2,0)*(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2)))+
 (complex<T>(1,0)*(-((SPA(2,5)*SPA(4,5)*(SPA(1,2)*SPB(4,1)-
SPA(2,3)*SPB(4,3)))/(SPA(2,3)*SPA(5,6)))-
  (SPA(2,4)*SPA(4,5)*SPB(6,4))/SPA(3,4)+
  (pow(SPB(3,1),2)*pow(SPB(6,4),2)*SPA(1,2))/
   (SPB(4,1)*SPB(4,3)*SPB(6,5))))/(complex<T>(2,0)*(SPA(1,3)*SPB(4,1)+
  SPA(2,3)*SPB(4,2))*(-(SPA(1,2)*SPB(4,2))-
  SPA(1,3)*SPB(4,3)))+
 (complex<T>(1,0)*pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+pow(SPA(1,4),2)*
pow(SPB(4,1),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*S(1,4)*S(5,6)+
   complex<T>(-2,0)*S(2,3)*S(5,6),-1)*(SPA(2,3)*SPB(3,1)+
  SPA(2,4)*SPB(4,1))*((complex<T>(1,0)*(-S(1,4)-S(2,3)+S(5,6))*
SPA(2,5)*SPB(6,1))/(complex<T>(2,0)*SPA(2,3)*SPB(4,1))-
  SPA(4,5)*SPB(6,3)-((-S(1,4)+S(2,3)-S(5,6))*SPB(6,1)*
SPB(6,3))/(SPB(4,1)*SPB(6,5))))/
(complex<T>(2,0)*(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2)))+
 (complex<T>(1,0)*(-((pow(SPA(2,4),2)*pow(SPA(3,5),2)*SPB(2,1))/
(SPA(2,3)*SPA(3,4)*SPA(5,6)))+(SPA(3,5)*SPB(3,1)*
SPB(6,3))/SPB(4,3)+((-(SPA(2,3)*SPB(2,1))+
 SPA(3,4)*SPB(4,1))*SPB(6,1)*SPB(6,3))/
   (SPB(4,1)*SPB(6,5))))/(complex<T>(2,0)*(SPA(1,3)*SPB(4,1)+
  SPA(2,3)*SPB(4,2))*(SPA(1,3)*SPB(2,1)+SPA(3,4)*SPB(4,2)))+
 (complex<T>(1,0)*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SPB(6,2)*
 SPB(6,3))/(complex<T>(2,0)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
 (SPA(1,3)*SPB(2,1)+SPA(3,4)*SPB(4,2))*SPB(4,3)*SPB(6,5))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2))*SPB(6,4)*
 ((((-S(1,2)-S(3,4)+S(5,6))*SPA(2,4)+complex<T>(-2,0)*SPA(1,2)*
SPA(3,4)*SPB(3,1))*(SPA(1,3)*SPB(6,1)+
 SPA(2,3)*SPB(6,2)))/SPA(3,4)+
  (SPA(1,2)*(-((-S(1,2)-S(3,4)+S(5,6))*SPB(3,1))+
 complex<T>(2,0)*SPA(2,4)*SPB(2,1)*SPB(4,3))*SPB(6,4))/SPB(4,3)))/
(complex<T>(2,0)*(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SPB(6,5))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2))*SPB(6,2)*
 (-(SPA(1,4)*(-((-S(1,2)-S(3,4)+S(5,6))*SPB(6,1))+
 complex<T>(2,0)*SPB(2,1)*(-(SPA(2,3)*SPB(6,3))-SPA(2,4)*
  SPB(6,4))))-(SPB(3,2)*((S(1,2)-S(3,4)-S(5,6))*
SPA(1,2)*SPB(6,1)+(-S(1,2)-S(3,4)+S(5,6))*SPA(2,5)*
SPB(6,5)))/SPB(4,3)))/(complex<T>(2,0)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2))*(SPA(1,3)*SPB(2,1)+SPA(3,4)*SPB(4,2))*
 SPB(6,5))+(SPA(1,3)*(SPA(1,2)*SPB(4,1)-SPA(2,3)*SPB(4,3))*
 (-((pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+pow(SPA(1,4),2)*
  pow(SPB(4,1),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*S(1,4)*S(5,6)+complex<T>(-2,0)*
  S(2,3)*S(5,6),-1)*SPA(4,5)*(-((-S(1,4)-S(2,3)+
   S(5,6))*SPA(1,2))+complex<T>(-2,0)*SPA(1,4)*SPA(2,3)*SPB(4,
  3))*SPB(6,1))/(SPA(1,3)*(-(SPA(1,2)*SPB(4,2))-
SPA(1,3)*SPB(4,3))))-
  (pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+pow(SPA(1,4),2)*pow(
  SPB(4,1),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*S(1,4)*S(5,6)+
complex<T>(-2,0)*S(2,3)*S(5,6),-1)*(-((S(1,4)-S(2,3)-S(5,6))*
 SPA(1,2)*SPB(4,1))-(-S(1,4)+S(2,3)-S(5,6))*
SPA(2,3)*SPB(4,3))*((pow(SPA(4,5),2)*SPB(4,1))/SPA(5,6)+
 (pow(SPB(6,1),2)*SPA(1,4))/SPB(6,5)))/(complex<T>(2,0)*SPA(1,3)*
SPB(4,1)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3)))+
  (pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)*SPA(2,4)*SPA(4,5)*
SPB(6,4))/(SPA(1,3)*(-(SPA(1,2)*SPB(4,2))-
 SPA(1,3)*SPB(4,3))*SS(1,2,3))))/(SPA(2,3)*
 (SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2)))+
 ((-(SPA(2,3)*SPB(2,1))+SPA(3,4)*SPB(4,1))*SPB(4,2)*
 ((pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+pow(SPA(1,4),2)*pow(
  SPB(4,1),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*S(1,4)*S(5,6)+
complex<T>(-2,0)*S(2,3)*S(5,6),-1)*SPA(2,5)*
(-((-S(1,4)-S(2,3)+S(5,6))*SPB(2,1))+complex<T>(-2,0)*SPA(3,4)*
SPB(3,2)*SPB(4,1))*SPB(6,3))/(SPB(4,2)*
(SPA(1,3)*SPB(2,1)+SPA(3,4)*SPB(4,2)))-
  (pow(pow(SPA(2,3),2)*pow(SPB(3,2),2)+pow(SPA(1,4),2)*pow(
  SPB(4,1),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,4)*S(2,3)+complex<T>(-2,0)*S(1,4)*S(5,6)+
complex<T>(-2,0)*S(2,3)*S(5,6),-1)*((-S(1,4)+S(2,3)-S(5,6))*
SPA(2,3)*SPB(2,1)+(S(1,4)-S(2,3)-S(5,6))*SPA(3,4)*
SPB(4,1))*(-((pow(SPA(2,5),2)*SPB(3,2))/SPA(5,6))-
 (pow(SPB(6,3),2)*SPA(2,3))/SPB(6,5)))/(complex<T>(2,0)*SPA(2,3)*
SPB(4,2)*(SPA(1,3)*SPB(2,1)+SPA(3,4)*SPB(4,2)))+
  (pow(complex<T>(1,0)-S(5,6)/SS(1,2,4),-1)*SPA(3,5)*SPB(3,1)*
SPB(6,3))/(SPB(4,2)*(SPA(1,3)*SPB(2,1)+SPA(3,4)*SPB(4,2))*
SS(1,2,4))))/(SPB(4,1)*(SPA(1,3)*SPB(4,1)+
  SPA(2,3)*SPB(4,2)))-(pow(SPA(2,4),2)*pow(SPB(6,2),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(1,3,4),-1)*(SPA(1,4)*SPB(2,1)-
  SPA(3,4)*SPB(3,2)))/(complex<T>(2,0)*SPA(3,4)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2))*(SPA(1,3)*SPB(2,1)+SPA(3,4)*SPB(4,2))*
 SPB(6,5)*SS(1,3,4))+(complex<T>(1,0)*pow(SPB(3,1),2)*SPA(1,4)*
 ((pow(-(SPA(1,3)*SPB(6,3))-SPA(1,4)*SPB(6,4),2)*
pow(complex<T>(1,0)-SS(1,3,4)/S(3,4),-1))/(pow(SPB(4,3),2)*
SPA(3,4)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2)))+
  (pow(SPA(1,3)*SPB(6,1)-SPA(3,4)*SPB(6,4),2)*
pow(complex<T>(1,0)-SS(1,3,4)/S(1,4),-1))/(pow(SPB(4,1),2)*
SPA(1,4)*(SPA(1,3)*SPB(2,1)+SPA(3,4)*SPB(4,2)))))/
(complex<T>(2,0)*SPA(1,3)*SPB(6,5)*SS(1,3,4))+
 (complex<T>(1,0)*pow(SPB(3,1),2)*(-(SPA(1,3)*SPB(6,3))-
  SPA(1,4)*SPB(6,4))*(SPA(1,3)*SPB(6,1)-SPA(3,4)*SPB(6,4)))/
(complex<T>(2,0)*SPA(1,3)*SPB(4,1)*(SPA(1,3)*SPB(2,1)+
  SPA(3,4)*SPB(4,2))*SPB(4,3)*SPB(6,5)*SS(1,3,4))-
 (pow(SPA(1,5),2)*pow(SPB(3,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)*(-(SPA(1,2)*SPB(3,2))+
  SPA(1,4)*SPB(4,3)))/(complex<T>(2,0)*SPA(5,6)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2))*SPB(4,3)*(-(SPA(1,2)*SPB(4,2))-
  SPA(1,3)*SPB(4,3))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*(-(SPA(3,5)*SPB(3,2))-
  SPA(4,5)*SPB(4,2))*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3)))/
(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*SPA(5,6)*SPB(4,2)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(3,2)*
 ((pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),2)*
pow(complex<T>(1,0)-SS(2,3,4)/S(3,4),-1))/(pow(SPA(3,4),2)*
(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(4,3))-
  (pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*
pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/(pow(SPA(2,3),2)*
SPB(3,2)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3)))))/
(complex<T>(2,0)*SPA(5,6)*SPB(4,2)*SS(2,3,4)))+
 (complex<T>(0,1)*((SPA(2,4)*SPA(2,5)*SPB(3,1)*SPB(6,1))/
 (S(3,4)*S(5,6)*SPA(2,3)*SPB(4,1))-
(pow(SPB(3,1),2)*SPA(2,5)*(SPA(1,4)*SPB(6,1)+
   SPA(3,4)*SPB(6,3)))/(S(3,4)*S(5,6)*SPB(4,1)*SS(1,3,4))-
(pow(SPA(2,4),2)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*
  SPB(6,1))/(S(3,4)*S(5,6)*SPA(2,3)*SS(2,3,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmmpemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, p, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmmpemep SLC");
#endif
 
 return( complex<T>(0,1)*(-((SPA(3,5)*(-(SPA(2,5)*SPA(3,4))+
   SPA(2,3)*SPA(4,5)))/(complex<T>(2,0)*SPA(1,4)*SPA(3,4)*SPA(5,6)*
  (SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))))-
 (pow(SPB(4,2),2)*pow(complex<T>(1,0)-S(3,4)/SS(2,3,4),-1)*SPA(2,3)*
 SPA(2,5)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2)))/
(complex<T>(2,0)*pow(SS(2,3,4),2)*SPA(5,6)*SPB(3,2)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2)))+
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(3,5)*
 ((-S(1,2)-S(3,4)+S(5,6))*SPA(2,4)+complex<T>(-2,0)*SPA(1,2)*
   SPA(3,4)*SPB(3,1))*SPB(6,4))/(SPA(1,4)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2)))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(3,5)*
 (-((-S(1,2)-S(3,4)+S(5,6))*SPB(3,1))+complex<T>(2,0)*SPA(2,4)*
   SPB(2,1)*SPB(4,3))*SPB(6,4))/(SPB(3,2)*(SPA(1,4)*SPB(3,1)+
  SPA(2,4)*SPB(3,2)))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(1,3)*
 (-((SPA(3,5)*(-((SPA(2,5)*((S(1,2)-S(3,4)-S(5,6))*SPA(1,3)*
 SPB(2,1)-(-S(1,2)+S(3,4)-S(5,6))*SPA(3,4)*
 SPB(4,2)))/SPA(5,6))-((-S(1,2)-S(3,4)+S(5,6))*
   SPA(1,3)+complex<T>(-2,0)*SPA(1,2)*SPA(3,4)*SPB(4,2))*SPB(6,
  1)))/SPA(3,4))+(SPB(6,4)*
(-(SPA(2,5)*(-((-S(1,2)-S(3,4)+S(5,6))*SPB(4,2))+
  complex<T>(2,0)*SPA(1,3)*SPB(2,1)*SPB(4,3)))-
 ((-((S(1,2)-S(3,4)-S(5,6))*SPA(1,2)*SPB(4,2))+
  (-S(1,2)+S(3,4)-S(5,6))*SPA(1,3)*SPB(4,3))*SPB(6,
  1))/SPB(6,5)))/SPB(4,3)))/(complex<T>(2,0)*SPA(1,4)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2)))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPB(4,2)*
 (-((SPA(3,5)*((SPA(2,5)*((S(1,2)-S(3,4)-S(5,6))*SPA(1,3)*
SPB(2,1)-(-S(1,2)+S(3,4)-S(5,6))*SPA(3,4)*
SPB(4,2)))/SPA(5,6)+((-S(1,2)-S(3,4)+S(5,6))*
   SPA(1,3)+complex<T>(-2,0)*SPA(1,2)*SPA(3,4)*SPB(4,2))*SPB(6,
  1)))/SPA(3,4))+(SPB(6,4)*
(SPA(2,5)*(-((-S(1,2)-S(3,4)+S(5,6))*SPB(4,2))+complex<T>(2,0)*
  SPA(1,3)*SPB(2,1)*SPB(4,3))+
 ((-((S(1,2)-S(3,4)-S(5,6))*SPA(1,2)*SPB(4,2))+
  (-S(1,2)+S(3,4)-S(5,6))*SPA(1,3)*SPB(4,3))*SPB(6,
  1))/SPB(6,5)))/SPB(4,3)))/(complex<T>(2,0)*SPB(3,2)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2)))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-((S(1,2)-S(3,4)-S(5,6))*
SPA(2,4)*SPB(2,1))+(-S(1,2)+S(3,4)-S(5,6))*SPA(3,4)*
   SPB(3,1))*(-((pow(SPA(3,5),2)*SPB(4,3))/SPA(5,6))-
  (pow(SPB(6,4),2)*SPA(3,4))/SPB(6,5)))/(complex<T>(2,0)*SPA(3,4)*
 SPB(3,2)*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2)))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*((S(1,2)-S(3,4)-S(5,6))*
   SPA(1,2)*SPB(3,1)-(-S(1,2)+S(3,4)-S(5,6))*SPA(2,4)*
   SPB(4,3))*(-((pow(SPA(3,5),2)*SPB(4,3))/SPA(5,6))-
  (pow(SPB(6,4),2)*SPA(3,4))/SPB(6,5)))/(complex<T>(2,0)*SPA(1,4)*
 (SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SPB(4,3))-
 ((SPB(4,3)*SPB(6,1)+SPB(4,1)*SPB(6,3))*SPB(6,4))/
(complex<T>(2,0)*SPB(3,2)*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*
 SPB(4,3)*SPB(6,5))+(complex<T>(1,0)*pow(SPA(1,3),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(1,3,4),-1)*SPB(4,1)*SPB(6,1)*
 (-(SPA(1,3)*SPB(6,3))-SPA(1,4)*SPB(6,4)))/
(complex<T>(2,0)*pow(SS(1,3,4),2)*SPA(1,4)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2))*SPB(6,5))-
 (pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)*SPA(4,5)*SPB(4,1)*SPB(6,4))/
(SPB(3,2)*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SS(1,2,3))-
 (pow(complex<T>(1,0)-S(5,6)/SS(1,2,4),-1)*SPA(2,3)*SPA(3,5)*SPB(6,3))/
(SPA(1,4)*(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SS(1,2,4))-
 (SPA(1,3)*SPB(6,4)*(-(SPB(6,4)/SPB(4,3))+
  (SPA(1,3)*SPB(6,1)-SPA(3,4)*SPB(6,4))/SS(1,3,4)))/
(complex<T>(2,0)*SPA(1,4)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*
 SPB(6,5))-(pow(SPA(2,3),2)*pow(SPB(6,2),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(1,3,4),-1)*SPA(1,3))/
(complex<T>(2,0)*SPA(1,4)*SPA(3,4)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2))*SPB(6,5)*SS(1,3,4))+
 (complex<T>(1,0)*SPA(3,5)*SPB(4,2)*(SPA(3,5)/SPA(3,4)+
  (SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))/SS(2,3,4)))/
(complex<T>(2,0)*SPA(5,6)*SPB(3,2)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2)))-(pow(SPA(1,5),2)*pow(SPB(4,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)*SPB(4,2))/
(complex<T>(2,0)*SPA(5,6)*SPB(3,2)*(SPA(1,3)*SPB(3,2)+
  SPA(1,4)*SPB(4,2))*SPB(4,3)*SS(2,3,4)))+
 (complex<T>(0,-1)*(-(((SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*
   (SPA(1,3)*SPB(6,1)-SPA(3,4)*SPB(6,4)))/(S(3,4)*S(5,6)*
   SPA(1,4)*SPB(3,2)))+(SPA(1,3)*SPA(2,5)*SPB(4,1)*
  (SPA(1,3)*SPB(6,1)-SPA(3,4)*SPB(6,4)))/(S(3,4)*S(5,6)*
  SPA(1,4)*SS(1,3,4))+(SPA(2,3)*SPB(4,2)*
  (SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*SPB(6,1))/
 (S(3,4)*S(5,6)*SPB(3,2)*SS(2,3,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmmmemep_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, m, em, ep}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmmmemep SLC");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(6,1),2))/(complex<T>(2,0)*SPB(3,2)*SPB(4,1)*SPB(4,3)*
 SPB(6,5))+complex<T>(0,1)*((pow(SPA(3,5),2)*pow(SPB(3,1),2)*
 pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1))/(pow(SPA(5,6),2)*SPB(3,2)*
 SPB(4,1)*SPB(4,3)*SPB(6,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*
 pow(SPB(2,1),2)*pow(complex<T>(1,0)-SS(1,3,4)/S(5,6),-1))/
(complex<T>(2,0)*pow(SPA(5,6),2)*SPB(3,2)*SPB(4,1)*SPB(4,3)*SPB(6,5))+
 (pow(SPA(4,5),2)*pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)*SPB(4,1))/
(pow(SPA(5,6),2)*SPB(3,2)*SPB(4,3)*SPB(6,5))-
 (pow(SPA(3,4),2)*pow(SPB(6,3),2)*
 pow(complex<T>(1,0)-S(1,4)/SS(1,3,4),-1)*SPB(4,1))/
(complex<T>(2,0)*pow(SS(1,3,4),2)*SPB(3,2)*SPB(4,3)*SPB(6,5))+
 (complex<T>(1,0)*(-(pow(SPB(6,1),2)/(SPB(3,2)*SPB(4,1)*SPB(4,3)*
 SPB(6,5)))+(SPA(3,4)*SPB(6,1)*SPB(6,3))/
   (SPB(3,2)*SPB(4,3)*SPB(6,5)*SS(1,3,4))+
  (SPA(3,4)*(-(SPA(3,5)*SPB(3,1))-SPA(4,5)*SPB(4,1))*
(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3)))/(SPA(5,6)*SPB(3,2)*
SPB(4,3)*SS(1,3,4)*SS(2,3,4))-
  (SPA(2,4)*SPA(3,4)*SPB(6,1)*SPB(6,2))/(SPB(3,2)*SPB(6,5)*
SS(1,3,4)*SS(2,3,4))))/complex<T>(2,0))
        ); 
  
} 
  
  
 
template <class T> complex<T> R2q2g2l_qpppqmemep_nf_top
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, p, qm, em, ep}, nf_top}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpppqmemep nf_top");
#endif

complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);
complex<T> mtsq = pow(mt, 2); 
complex<T> top = (complex<T>(1,0)/complex<T>(20,0)*S(2,3)/mtsq);

 return( top*(complex<T>(0,-1)*(-((SPA(4,5)*SPB(3,1)*(SPA(1,3)*SPB(6,1)+
   SPA(2,3)*SPB(6,2)))/SS(1,2,3))-
 (SPA(3,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*SPB(6,1))/
SS(2,3,4)))/(complex<T>(3,0)*pow(SPA(2,3),2)*S(5,6))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qppmqmemep_nf_top
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, m, qm, em, ep}, nf_top}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qppmqmemep nf_top");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g2l_qpmpqmemep_nf_top
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, p, qm, em, ep}, nf_top}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpmpqmemep nf_top");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g2l_qpmmqmemep_nf_top
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, m, qm, em, ep}, nf_top}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpmmqmemep nf_top");
#endif

complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);
complex<T> mtsq = pow(mt, 2); 
complex<T> top = (complex<T>(1,0)/complex<T>(20,0)*S(2,3)/mtsq);

 return( top*(complex<T>(0,-1)*((SPA(4,5)*SPB(2,1)*(SPA(1,2)*SPB(6,1)-
  SPA(2,3)*SPB(6,3)))/SS(1,2,3)+
 (SPA(2,4)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(6,1))/
SS(2,3,4)))/(complex<T>(3,0)*pow(SPB(3,2),2)*S(5,6))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpppqmemep_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, p, qm, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpppqmemep nf");
#endif
 
 return( (complex<T>(0,-1)*(-((SPA(4,5)*SPB(3,1)*(SPA(1,3)*SPB(6,1)+
   SPA(2,3)*SPB(6,2)))/SS(1,2,3))-
 (SPA(3,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*SPB(6,1))/
SS(2,3,4)))/(complex<T>(3,0)*pow(SPA(2,3),2)*S(5,6))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qppmqmemep_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, m, qm, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qppmqmemep nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g2l_qpmpqmemep_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, p, qm, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpmpqmemep nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g2l_qpmmqmemep_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, m, qm, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpmmqmemep nf");
#endif
 
 return( (complex<T>(0,-1)*((SPA(4,5)*SPB(2,1)*(SPA(1,2)*SPB(6,1)-
  SPA(2,3)*SPB(6,3)))/SS(1,2,3)+
 (SPA(2,4)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(6,1))/
SS(2,3,4)))/(complex<T>(3,0)*pow(SPB(3,2),2)*S(5,6))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmppemep_VECT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, p, em, ep}, VECT}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmppemep VECT");
#endif
 
 return( complex<T>(0,-1)*(-((pow(SPA(3,5),2)*pow(SPB(3,1),2)*
 pow(complex<T>(1,0)-SS(1,2,3)/S(1,2),-1))/(pow(SPA(3,4),2)*
 pow(SPB(2,1),2)*SPA(1,2)*SPA(5,6)))-
(pow(SPA(4,5),2)*pow(SPB(4,1),2)*pow(complex<T>(1,0)-SS(1,2,4)/S(1,2),
-1))/(pow(SPA(3,4),2)*pow(SPB(2,1),2)*SPA(1,2)*SPA(5,6))-
(pow(SPA(2,3),2)*pow(SPB(6,3),2)*pow(complex<T>(1,0)-SS(3,5,6)/S(5,6),
-1))/(pow(SPA(3,4),2)*pow(SPB(6,5),2)*SPA(1,2)*SPA(5,6))-
(pow(SPA(2,4),2)*pow(SPB(6,4),2)*pow(complex<T>(1,0)-SS(4,5,6)/S(5,6),
-1))/(pow(SPA(3,4),2)*pow(SPB(6,5),2)*SPA(1,2)*SPA(5,6))+
(pow(SPA(2,5),2)/(SPA(1,2)*SPA(5,6))-pow(SPB(6,1),2)/
 (SPB(2,1)*SPB(6,5)))/pow(SPA(3,4),2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmpmemep_VECT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, m, em, ep}, VECT}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmpmemep VECT");
#endif
 
 return( complex<T>(0,-1)*((pow(SPA(4,5),2)*pow(SPB(4,1),2)*
pow(complex<T>(1,0)-SS(1,2,4)/S(1,2),-1))/
 (pow(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2),2)*SPA(5,6)*
SPB(2,1))+(pow(SPA(3,5),2)*pow(SPB(3,1),2)*
pow(complex<T>(1,0)-SS(3,5,6)/S(5,6),-1))/
 (pow(SPA(3,5)*SPB(5,4)+SPA(3,6)*SPB(6,4),2)*SPA(5,6)*
SPB(2,1))+(pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
  pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
   pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
  complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(2,5)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*
((-S(1,2)+S(3,4)-S(5,6))*SPA(2,5)+complex<T>(2,0)*SPA(1,2)*
  SPA(5,6)*SPB(6,1)))/(SPA(1,2)*SPA(5,6)*(SPA(1,3)*SPB(4,1)+
 SPA(2,3)*SPB(4,2)))-
(pow(-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,1),2)-
SPA(2,5)*SPA(5,6)*SPB(2,1)*SPB(6,1))/
 (pow(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2),2)*SPA(5,6)*
SPB(2,1))+(pow(SPA(2,3),2)*pow(SPB(6,3),2)*
pow(complex<T>(1,0)-SS(1,2,3)/S(1,2),-1))/
 (pow(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2),2)*SPA(1,2)*
SPB(6,5))+(pow(SPA(2,4),2)*pow(SPB(6,4),2)*
pow(complex<T>(1,0)-SS(4,5,6)/S(5,6),-1))/
 (pow(SPA(3,5)*SPB(5,4)+SPA(3,6)*SPB(6,4),2)*SPA(1,2)*
SPB(6,5))-(pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
  pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
   pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
  complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SPB(6,1)*
(-((-S(1,2)+S(3,4)-S(5,6))*SPB(6,1))+complex<T>(-2,0)*SPA(2,5)*
  SPB(2,1)*SPB(6,5)))/(SPB(2,1)*(SPA(1,3)*SPB(4,1)+
 SPA(2,3)*SPB(4,2))*SPB(6,5))-
(pow(SPA(1,2)*SPB(6,1)-SPA(2,4)*SPB(6,4),2)-
SPA(1,2)*SPA(2,5)*SPB(6,1)*SPB(6,5))/
 (pow(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2),2)*SPA(1,2)*
SPB(6,5)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmmpemep_VECT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, p, em, ep}, VECT}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmmpemep VECT");
#endif
 
 return( complex<T>(0,-1)*((pow(SPA(3,5),2)*pow(SPB(3,1),2)*
pow(complex<T>(1,0)-SS(1,2,3)/S(1,2),-1))/
 (pow(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2),2)*SPA(5,6)*
SPB(2,1))+(pow(SPA(4,5),2)*pow(SPB(4,1),2)*
pow(complex<T>(1,0)-SS(4,5,6)/S(5,6),-1))/
 (pow(SPA(4,5)*SPB(5,3)+SPA(4,6)*SPB(6,3),2)*SPA(5,6)*
SPB(2,1))+(pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
  pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
   pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
  complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(2,5)*
(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*
((-S(1,2)+S(3,4)-S(5,6))*SPA(2,5)+complex<T>(2,0)*SPA(1,2)*
  SPA(5,6)*SPB(6,1)))/(SPA(1,2)*SPA(5,6)*(SPA(1,4)*SPB(3,1)+
 SPA(2,4)*SPB(3,2)))-
(pow(-(SPA(2,5)*SPB(2,1))-SPA(4,5)*SPB(4,1),2)-
SPA(2,5)*SPA(5,6)*SPB(2,1)*SPB(6,1))/
 (pow(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2),2)*SPA(5,6)*
SPB(2,1))+(pow(SPA(2,4),2)*pow(SPB(6,4),2)*
pow(complex<T>(1,0)-SS(1,2,4)/S(1,2),-1))/
 (pow(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2),2)*SPA(1,2)*
SPB(6,5))+(pow(SPA(2,3),2)*pow(SPB(6,3),2)*
pow(complex<T>(1,0)-SS(3,5,6)/S(5,6),-1))/
 (pow(SPA(4,5)*SPB(5,3)+SPA(4,6)*SPB(6,3),2)*SPA(1,2)*
SPB(6,5))-(pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
  pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
   pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
  complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*SPB(6,1)*
(-((-S(1,2)+S(3,4)-S(5,6))*SPB(6,1))+complex<T>(-2,0)*SPA(2,5)*
  SPB(2,1)*SPB(6,5)))/(SPB(2,1)*(SPA(1,4)*SPB(3,1)+
 SPA(2,4)*SPB(3,2))*SPB(6,5))-
(pow(SPA(1,2)*SPB(6,1)-SPA(2,3)*SPB(6,3),2)-
SPA(1,2)*SPA(2,5)*SPB(6,1)*SPB(6,5))/
 (pow(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2),2)*SPA(1,2)*
SPB(6,5)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmmmemep_VECT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, m, em, ep}, VECT}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmmmemep VECT");
#endif
 
 return( complex<T>(0,-1)*((-(pow(SPA(2,5),2)/(SPA(1,2)*SPA(5,6)))+
pow(SPB(6,1),2)/(SPB(2,1)*SPB(6,5)))/pow(SPB(4,3),2)-
(pow(SPA(2,3),2)*pow(SPB(6,3),2)*pow(complex<T>(1,0)-SS(1,2,3)/S(1,2),
-1))/(pow(SPA(1,2),2)*pow(SPB(4,3),2)*SPB(2,1)*SPB(6,5))-
(pow(SPA(2,4),2)*pow(SPB(6,4),2)*pow(complex<T>(1,0)-SS(1,2,4)/S(1,2),
-1))/(pow(SPA(1,2),2)*pow(SPB(4,3),2)*SPB(2,1)*SPB(6,5))-
(pow(SPA(3,5),2)*pow(SPB(3,1),2)*pow(complex<T>(1,0)-SS(3,5,6)/S(5,6),
-1))/(pow(SPA(5,6),2)*pow(SPB(4,3),2)*SPB(2,1)*SPB(6,5))-
(pow(SPA(4,5),2)*pow(SPB(4,1),2)*pow(complex<T>(1,0)-SS(4,5,6)/S(5,6),
-1))/(pow(SPA(5,6),2)*pow(SPB(4,3),2)*SPB(2,1)*SPB(6,5)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmppemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, p, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmppemep AX");
#endif

 complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);


 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1)*SPA(2,3)*
 SPA(2,5)*SPB(3,1)*SPB(6,3))/(pow(SPA(5,6),2)*pow(SPB(6,5),2)*
 SPA(2,4)*SPA(3,4)))-(pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1)*
SPA(2,3)*SPA(2,5)*SPB(4,3)*SPB(6,3))/(pow(SPA(5,6),2)*
pow(SPB(6,5),2)*SPA(1,2)*SPA(3,4))-
(pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)*(S(1,4)+S(3,4))*SPA(2,5)*
SPB(6,4))/(pow(SPA(5,6),2)*pow(SPB(6,5),2)*SPA(1,3)*
SPA(3,4))+(pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)*SPA(2,4)*
SPA(2,5)*SPB(4,3)*SPB(6,4))/(pow(SPA(5,6),2)*pow(SPB(6,5),2)*
SPA(1,2)*SPA(3,4))-
(SPA(2,5)*(-(((SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(6,3))/
   SPA(2,4))-(SPB(6,4)*SS(1,3,4))/SPA(1,3)))/
 (mtsq*complex<T>(12,0)*S(5,6)*SPA(3,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmpmemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, m, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmpmemep AX");
#endif

 complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);

 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( complex<T>(0,1)*((complex<T>(1,0)*SPB(3,1)*(SPA(2,5)*SPB(3,2)-
 SPA(4,5)*SPB(4,3))*SPB(6,3))/(mtsq*complex<T>(12,0)*S(5,6)*SPB(4,2)*
SPB(4,3))-(SPA(2,4)*SPA(4,5)*(SPA(1,4)*SPB(6,1)+
 SPA(3,4)*SPB(6,3)))/(mtsq*complex<T>(12,0)*S(5,6)*SPA(1,3)*SPA(3,4))+
(SPA(2,4)*SPA(3,5)*(SPA(1,4)*SPB(6,1)+SPA(3,4)*SPB(6,3)))/
 (S(5,6)*SPA(1,3)*SPA(3,4)*(SPA(1,3)*SPB(4,1)+
 SPA(2,3)*SPB(4,2)))-(SPB(3,1)*(SPA(2,5)*SPB(3,2)-
 SPA(4,5)*SPB(4,3))*SPB(6,4))/(S(5,6)*SPB(4,2)*
(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*SPB(4,3))-
(SPA(4,5)*SPB(3,1)*(-((pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(
  SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(
  5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,2)-S(3,4)-
 S(5,6))*SPA(2,4)*SPA(3,5))/(SPA(3,4)*SPA(5,6)))+
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
 complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)+S(3,4)-S(5,6))*
   SPA(3,5)*SPB(3,1))/(SPA(5,6)*SPB(2,1))+
 complex<T>(2,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(2,4)*SPB(6,4)+
 (SPB(3,1)*SPB(6,4))/(S(5,6)*SPB(2,1)*SPB(4,3))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
 complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)-S(3,4)+S(5,6))*
   SPB(3,1)*SPB(6,4))/(SPB(2,1)*SPB(4,3))))/
 (SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))-
(SPA(2,4)*SPB(6,3)*((SPA(2,4)*SPA(3,5))/(S(5,6)*SPA(1,2)*
   SPA(3,4))-(pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
 pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
 complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
   (-S(1,2)-S(3,4)+S(5,6))*SPA(2,4)*SPA(3,5))/
  (SPA(1,2)*SPA(3,4))+complex<T>(2,0)*
  pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(3,5)*SPB(3,1)+
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
 complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)+S(3,4)-S(5,6))*
   SPA(2,4)*SPB(6,4))/(SPA(1,2)*SPB(6,5))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
 complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,2)-S(3,4)-S(5,6))*
   SPB(3,1)*SPB(6,4))/(SPB(4,3)*SPB(6,5))))/
 (SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))-
(pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)*SPA(1,4)*SPA(2,4)*
(SPA(1,2)*SPB(6,1)-SPA(2,3)*SPB(6,3))*SPB(6,4))/
 (SPA(1,2)*SPA(1,3)*(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*
SPB(6,5)*SS(1,2,3))-(pow(complex<T>(1,0)-S(5,6)/SS(1,2,4),-1)*
SPA(3,5)*SPB(3,1)*SPB(3,2)*(-(SPA(2,5)*SPB(2,1))-
 SPA(4,5)*SPB(4,1)))/(SPA(5,6)*SPB(2,1)*SPB(4,2)*
(SPA(1,3)*SPB(4,1)+SPA(2,3)*SPB(4,2))*SS(1,2,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmmpemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, p, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmmpemep AX");
#endif

 complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);


 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( complex<T>(0,1)*(-((pow(SPA(2,3),2)*SPA(3,5)*SPB(6,1))/
(mtsq*complex<T>(12,0)*S(5,6)*SPA(2,4)*SPA(3,4)))+
(pow(SPA(2,3),2)*SPA(4,5)*SPB(6,1))/(S(5,6)*SPA(2,4)*SPA(3,4)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2)))-
(pow(SPB(4,1),2)*SPA(2,5)*SPB(6,3))/(S(5,6)*SPB(3,1)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SPB(4,3))+
(SPA(3,5)*SPB(4,1)*((pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
 pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
 complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
   (S(1,2)-S(3,4)-S(5,6))*SPA(2,3)*SPA(4,5))/
  (SPA(3,4)*SPA(5,6))+(pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
 pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
 complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
   (-S(1,2)+S(3,4)-S(5,6))*SPA(4,5)*SPB(4,1))/
  (SPA(5,6)*SPB(2,1))+complex<T>(2,0)*
  pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(2,3)*SPB(6,3)-
 (SPB(4,1)*SPB(6,3))/(S(5,6)*SPB(2,1)*SPB(4,3))+
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
 complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)-S(3,4)+S(5,6))*
   SPB(4,1)*SPB(6,3))/(SPB(2,1)*SPB(4,3))))/
 (SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))+
(complex<T>(1,0)*pow(SPB(4,1),2)*SPA(2,5)*SPB(6,4))/(mtsq*complex<T>(12,0)*S(5,6)*
SPB(3,1)*SPB(4,3))+(SPA(2,3)*SPB(6,4)*
(-((SPA(2,3)*SPA(4,5))/(S(5,6)*SPA(1,2)*SPA(3,4)))+
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
 complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)-S(3,4)+S(5,6))*
   SPA(2,3)*SPA(4,5))/(SPA(1,2)*SPA(3,4))+
 complex<T>(2,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*SPA(4,5)*SPB(4,1)+
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
 complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)+S(3,4)-S(5,6))*
   SPA(2,3)*SPB(6,3))/(SPA(1,2)*SPB(6,5))+
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
 complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
 complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,2)-S(3,4)-S(5,6))*
   SPB(4,1)*SPB(6,3))/(SPB(4,3)*SPB(6,5))))/
 (SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))+
(pow(SPB(4,1),2)*pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)*SPA(4,5)*
(-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,1)))/
 (SPA(5,6)*SPB(2,1)*SPB(3,1)*(SPA(1,4)*SPB(3,1)+
 SPA(2,4)*SPB(3,2))*SS(1,2,3))+
(pow(SPA(2,3),2)*pow(complex<T>(1,0)-S(5,6)/SS(1,2,4),-1)*SPB(6,3)*
(SPA(1,2)*SPB(6,1)-SPA(2,4)*SPB(6,4)))/(SPA(1,2)*SPA(2,4)*
(SPA(1,4)*SPB(3,1)+SPA(2,4)*SPB(3,2))*SPB(6,5)*SS(1,2,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmmmemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, m, em, ep}, AX}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmmmemep AX");
#endif
 complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);

 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1)*SPA(3,4)*
 SPA(3,5)*SPB(3,1)*SPB(6,1))/(pow(SPA(5,6),2)*pow(SPB(6,5),2)*
 SPB(2,1)*SPB(4,3)))+(pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)*
SPA(3,4)*SPA(4,5)*SPB(4,1)*SPB(6,1))/(pow(SPA(5,6),2)*
pow(SPB(6,5),2)*SPB(2,1)*SPB(4,3))+
(pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)*SPA(2,4)*SPA(4,5)*SPB(4,1)*
SPB(6,1))/(pow(SPA(5,6),2)*pow(SPB(6,5),2)*SPB(3,1)*
SPB(4,3))+(pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1)*
(S(2,3)+S(3,4))*SPA(3,5)*SPB(6,1))/(pow(SPA(5,6),2)*
pow(SPB(6,5),2)*SPB(4,2)*SPB(4,3))+
(complex<T>(1,0)*SPB(6,1)*(-((SPA(4,5)*(SPA(2,3)*SPB(3,1)+
 SPA(2,4)*SPB(4,1)))/SPB(3,1))-(SPA(3,5)*SS(2,3,4))/
  SPB(4,2)))/(mtsq*complex<T>(12,0)*S(5,6)*SPB(4,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmppemep_AXSL
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, p, em, ep}, AXSL}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmppemep AXSL");
#endif
 complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);

 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( complex<T>(0,1)*((complex<T>(-2,0)*(complex<T>(1,0)/(mtsq*complex<T>(24,0))-
 pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1)/(complex<T>(2,0)*S(5,6)))*SPA(2,5)*
(SPA(1,2)*SPB(3,1)+SPA(2,4)*SPB(4,3))*SPB(6,3))/
 (S(5,6)*SPA(1,4)*SPA(2,4))+
(complex<T>(-2,0)*(complex<T>(1,0)/(mtsq*complex<T>(24,0))-pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),
  -1)/(complex<T>(2,0)*S(5,6)))*SPA(2,5)*(SPA(1,2)*SPB(4,1)-
 SPA(2,3)*SPB(4,3))*SPB(6,4))/(S(5,6)*SPA(1,3)*SPA(2,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmpmemep_AXSL
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, m, em, ep}, AXSL}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmpmemep AXSL");
#endif
 complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);

 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( complex<T>(0,1)*((complex<T>(2,0)*(complex<T>(1,0)/(mtsq*complex<T>(24,0))-
 pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),-1)/(complex<T>(2,0)*S(5,6)))*SPB(3,1)*
(-(SPA(2,5)*SPB(2,1))-SPA(4,5)*SPB(4,1))*SPB(6,3))/
 (S(5,6)*SPB(4,1)*SPB(4,2))+
(complex<T>(2,0)*(complex<T>(1,0)/(mtsq*complex<T>(24,0))-pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),
  -1)/(complex<T>(2,0)*S(5,6)))*SPA(2,4)*SPA(4,5)*(SPA(1,2)*SPB(6,1)-
 SPA(2,3)*SPB(6,3)))/(S(5,6)*SPA(1,3)*SPA(2,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmmpemep_AXSL
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, p, em, ep}, AXSL}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmmpemep AXSL");
#endif
 complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);

 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( complex<T>(0,1)*((complex<T>(2,0)*(complex<T>(1,0)/(mtsq*complex<T>(24,0))-
 pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)/(complex<T>(2,0)*S(5,6)))*
(-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,1))*SPB(4,1)*SPB(6,4))/
 (S(5,6)*SPB(3,1)*SPB(3,2))+
(complex<T>(2,0)*(complex<T>(1,0)/(mtsq*complex<T>(24,0))-pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),
  -1)/(complex<T>(2,0)*S(5,6)))*SPA(2,3)*SPA(3,5)*(SPA(1,2)*SPB(6,1)-
 SPA(2,4)*SPB(6,4)))/(S(5,6)*SPA(1,4)*SPA(2,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g2l_qpqmmmemep_AXSL
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, m, em, ep}, AXSL}
 
#if _VERBOSE
  _MESSAGE("R2q2g2l :  qpqmmmemep AXSL");
#endif
 complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);

 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( complex<T>(0,1)*((complex<T>(-2,0)*(complex<T>(1,0)/(mtsq*complex<T>(24,0))-
 pow(complex<T>(1,0)-SS(1,2,3)/S(5,6),-1)/(complex<T>(2,0)*S(5,6)))*SPA(4,5)*
(-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SPB(6,1))/
 (S(5,6)*SPB(3,1)*SPB(3,2))+
(complex<T>(-2,0)*(complex<T>(1,0)/(mtsq*complex<T>(24,0))-pow(complex<T>(1,0)-SS(1,2,4)/S(5,6),
  -1)/(complex<T>(2,0)*S(5,6)))*SPA(3,5)*(-(SPA(2,3)*SPB(2,1))+
 SPA(3,4)*SPB(4,1))*SPB(6,1))/(S(5,6)*SPB(4,1)*SPB(4,2)))
        ); 
  
} 
  
  
 
 // *************** table of switch values ************* 
 
#define _R_qpppqmemep_L R2q2g2l_44408_L
#define _R_qppmqmemep_L R2q2g2l_44300_L
#define _R_qpmpqmemep_L R2q2g2l_44390_L
#define _R_qpmmqmemep_L R2q2g2l_44282_L
#define _R_qmppqpemep_L R2q2g2l_44623_L
#define _R_qmpmqpemep_L R2q2g2l_44515_L
#define _R_qmmpqpemep_L R2q2g2l_44605_L
#define _R_qmmmqpemep_L R2q2g2l_44497_L
#define _R_qppqmpemep_SLC R2q2g2l_44768_SLC
#define _R_qppqmmemep_SLC R2q2g2l_44120_SLC
#define _R_qpmqmpemep_SLC R2q2g2l_44750_SLC
#define _R_qpmqmmemep_SLC R2q2g2l_44102_SLC
#define _R_qpqmppemep_SLC R2q2g2l_44828_SLC
#define _R_qpqmpmemep_SLC R2q2g2l_44180_SLC
#define _R_qpqmmpemep_SLC R2q2g2l_44720_SLC
#define _R_qpqmmmemep_SLC R2q2g2l_44072_SLC
#define _R_qpppqmemep_nf_top R2q2g2l_44408_nf_top
#define _R_qppmqmemep_nf_top R2q2g2l_44300_nf_top
#define _R_qpmpqmemep_nf_top R2q2g2l_44390_nf_top
#define _R_qpmmqmemep_nf_top R2q2g2l_44282_nf_top
#define _R_qpppqmemep_nf R2q2g2l_44408_nf
#define _R_qppmqmemep_nf R2q2g2l_44300_nf
#define _R_qpmpqmemep_nf R2q2g2l_44390_nf
#define _R_qpmmqmemep_nf R2q2g2l_44282_nf

#define _R_qpqmppemep_VECT R2q2g2l_44828_VECT
#define _R_qpqmpmemep_VECT R2q2g2l_44180_VECT
#define _R_qpqmmpemep_VECT R2q2g2l_44720_VECT
#define _R_qpqmmmemep_VECT R2q2g2l_44072_VECT
#define _R_qpqmppemep_AX R2q2g2l_44828_AX
#define _R_qpqmpmemep_AX R2q2g2l_44180_AX
#define _R_qpqmmpemep_AX R2q2g2l_44720_AX
#define _R_qpqmmmemep_AX R2q2g2l_44072_AX
#define _R_qpqmppemep_AXSL R2q2g2l_44828_AXSL
#define _R_qpqmpmemep_AXSL R2q2g2l_44180_AXSL
#define _R_qpqmmpemep_AXSL R2q2g2l_44720_AXSL
#define _R_qpqmmmemep_AXSL R2q2g2l_44072_AXSL
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qpppqmemep_L case 44408 : \
          return &R2q2g2l_44408_L
#define _CASE_qppmqmemep_L case 44300 : \
          return &R2q2g2l_44300_L
#define _CASE_qpmpqmemep_L case 44390 : \
          return &R2q2g2l_44390_L
#define _CASE_qpmmqmemep_L case 44282 : \
          return &R2q2g2l_44282_L
#define _CASE_qmppqpemep_L case 44623 : \
          return &R2q2g2l_44623_L
#define _CASE_qmpmqpemep_L case 44515 : \
          return &R2q2g2l_44515_L
#define _CASE_qmmpqpemep_L case 44605 : \
          return &R2q2g2l_44605_L
#define _CASE_qmmmqpemep_L case 44497 : \
          return &R2q2g2l_44497_L
#define _CASE_qppqmpemep_SLC case 44768 : \
          return &R2q2g2l_44768_SLC
#define _CASE_qppqmmemep_SLC case 44120 : \
          return &R2q2g2l_44120_SLC
#define _CASE_qpmqmpemep_SLC case 44750 : \
          return &R2q2g2l_44750_SLC
#define _CASE_qpmqmmemep_SLC case 44102 : \
          return &R2q2g2l_44102_SLC
#define _CASE_qpqmppemep_SLC case 44828 : \
          return &R2q2g2l_44828_SLC
#define _CASE_qpqmpmemep_SLC case 44180 : \
          return &R2q2g2l_44180_SLC
#define _CASE_qpqmmpemep_SLC case 44720 : \
          return &R2q2g2l_44720_SLC
#define _CASE_qpqmmmemep_SLC case 44072 : \
          return &R2q2g2l_44072_SLC
#define _CASE_qpppqmemep_nf_top case 44408 : \
          return &R2q2g2l_44408_nf_top
#define _CASE_qppmqmemep_nf_top case 44300 : \
          return &R2q2g2l_44300_nf_top
#define _CASE_qpmpqmemep_nf_top case 44390 : \
          return &R2q2g2l_44390_nf_top
#define _CASE_qpmmqmemep_nf_top case 44282 : \
          return &R2q2g2l_44282_nf_top
#define _CASE_qpppqmemep_nf case 44408 : \
          return &R2q2g2l_44408_nf
#define _CASE_qppmqmemep_nf case 44300 : \
          return &R2q2g2l_44300_nf
#define _CASE_qpmpqmemep_nf case 44390 : \
          return &R2q2g2l_44390_nf
#define _CASE_qpmmqmemep_nf case 44282 : \
          return &R2q2g2l_44282_nf
#define _CASE_qpqmppemep_VECT case 44828 : \
          return &R2q2g2l_44828_VECT
#define _CASE_qpqmpmemep_VECT case 44180 : \
          return &R2q2g2l_44180_VECT
#define _CASE_qpqmmpemep_VECT case 44720 : \
          return &R2q2g2l_44720_VECT
#define _CASE_qpqmmmemep_VECT case 44072 : \
          return &R2q2g2l_44072_VECT
#define _CASE_qpqmppemep_AX case 44828 : \
          return &R2q2g2l_44828_AX
#define _CASE_qpqmpmemep_AX case 44180 : \
          return &R2q2g2l_44180_AX
#define _CASE_qpqmmpemep_AX case 44720 : \
          return &R2q2g2l_44720_AX
#define _CASE_qpqmmmemep_AX case 44072 : \
          return &R2q2g2l_44072_AX
#define _CASE_qpqmppemep_AXSL case 44828 : \
          return &R2q2g2l_44828_AXSL
#define _CASE_qpqmpmemep_AXSL case 44180 : \
          return &R2q2g2l_44180_AXSL
#define _CASE_qpqmmpemep_AXSL case 44720 : \
          return &R2q2g2l_44720_AXSL
#define _CASE_qpqmmmemep_AXSL case 44072 : \
          return &R2q2g2l_44072_AXSL
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_qpppqmemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpppqmemep_L(ep,mpc);}
 
template <class T> complex<T> _R_qppmqmemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qppmqmemep_L(ep,mpc);}
 
template <class T> complex<T> _R_qpmpqmemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpmpqmemep_L(ep,mpc);}
 
template <class T> complex<T> _R_qpmmqmemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpmmqmemep_L(ep,mpc);}
 
template <class T> complex<T> _R_qmppqpemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qmppqpemep_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpmqpemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qmpmqpemep_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmpqpemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qmmpqpemep_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmmqpemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qmmmqpemep_L(ep,mpc);}
 
template <class T> complex<T> _R_qppqmpemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qppqmpemep_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qppqmmemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qppqmmemep_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpmqmpemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpmqmpemep_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpmqmmemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpmqmmemep_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpqmppemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmppemep_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpqmpmemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmpmemep_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmpemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmmpemep_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmmemep_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmmmemep_SLC(ep,mpc);}

template <class T> complex<T> _R_qpppqmemep_nf_top(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpppqmemep_nf_top(ep,mpc);}
 
template <class T> complex<T> _R_qppmqmemep_nf_top(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qppmqmemep_nf_top(ep,mpc);}
 
template <class T> complex<T> _R_qpmpqmemep_nf_top(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpmpqmemep_nf_top(ep,mpc);}
 
template <class T> complex<T> _R_qpmmqmemep_nf_top(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpmmqmemep_nf_top(ep,mpc);}

template <class T> complex<T> _R_qpppqmemep_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpppqmemep_nf(ep,mpc);}
 
template <class T> complex<T> _R_qppmqmemep_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qppmqmemep_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpmpqmemep_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpmpqmemep_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpmmqmemep_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpmmqmemep_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpqmppemep_VECT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmppemep_VECT(ep,mpc);}
 
template <class T> complex<T> _R_qpqmpmemep_VECT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmpmemep_VECT(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmpemep_VECT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmmpemep_VECT(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmmemep_VECT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmmmemep_VECT(ep,mpc);}
 
template <class T> complex<T> _R_qpqmppemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmppemep_AX(ep,mpc);}
 
template <class T> complex<T> _R_qpqmpmemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmpmemep_AX(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmpemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmmpemep_AX(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmmemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmmmemep_AX(ep,mpc);}
 
template <class T> complex<T> _R_qpqmppemep_AXSL(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmppemep_AXSL(ep,mpc);}
 
template <class T> complex<T> _R_qpqmpmemep_AXSL(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmpmemep_AXSL(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmpemep_AXSL(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmmpemep_AXSL(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmmemep_AXSL(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g2l_qpqmmmemep_AXSL(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q2g2l_AX_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpqmppemep_AX;
       _CASE_qpqmpmemep_AX;
       _CASE_qpqmmpemep_AX;
       _CASE_qpqmmmemep_AX;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g2l_AXSL_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpqmppemep_AXSL;
       _CASE_qpqmpmemep_AXSL;
       _CASE_qpqmmpemep_AXSL;
       _CASE_qpqmmmemep_AXSL;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g2l_L_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpppqmemep_L;
       _CASE_qppmqmemep_L;
       _CASE_qpmpqmemep_L;
       _CASE_qpmmqmemep_L;
       _CASE_qmppqpemep_L;
       _CASE_qmpmqpemep_L;
       _CASE_qmmpqpemep_L;
       _CASE_qmmmqpemep_L;
 
       default: return 0;
        }
 }

template <class T> complex<T> ( *R2q2g2l_nf_top_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpppqmemep_nf_top;
       _CASE_qppmqmemep_nf_top;
       _CASE_qpmpqmemep_nf_top;
       _CASE_qpmmqmemep_nf_top;
 
       default: return 0;
        }
 }


template <class T> complex<T> ( *R2q2g2l_nf_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpppqmemep_nf;
       _CASE_qppmqmemep_nf;
       _CASE_qpmpqmemep_nf;
       _CASE_qpmmqmemep_nf;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g2l_SLC_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qppqmpemep_SLC;
       _CASE_qppqmmemep_SLC;
       _CASE_qpmqmpemep_SLC;
       _CASE_qpmqmmemep_SLC;
       _CASE_qpqmppemep_SLC;
       _CASE_qpqmpmemep_SLC;
       _CASE_qpqmmpemep_SLC;
       _CASE_qpqmmmemep_SLC;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g2l_VECT_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpqmppemep_VECT;
       _CASE_qpqmpmemep_VECT;
       _CASE_qpqmmpemep_VECT;
       _CASE_qpqmmmemep_VECT;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q2g2l_AX_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g2l_AX_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g2l_AX_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g2l_AX_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2g2l_AXSL_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g2l_AXSL_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g2l_AXSL_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g2l_AXSL_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2g2l_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g2l_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g2l_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g2l_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2g2l_nf_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g2l_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g2l_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g2l_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif

template complex<R> ( *R2q2g2l_nf_top_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g2l_nf_top_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g2l_nf_top_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g2l_nf_top_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2g2l_SLC_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g2l_SLC_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g2l_SLC_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g2l_SLC_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2g2l_VECT_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g2l_VECT_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g2l_VECT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g2l_VECT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif




}

