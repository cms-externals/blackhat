/*
* R_2q2g1y_wCI.cpp
*
* Created on 2/5, 2011
*      Author: Zvi's script
*/
 
#include "R_2q2g1y_eval.h"
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


 
 
template <class T> complex<T> R2q2g1y_qmmmqpgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmmmqpgam L");
#endif
 
 return( (complex<T>(0,-1)*((pow(SPB(4,1),2)*SPA(1,5))/(SPB(2,1)*SPB(3,2)*
  SPB(4,3)*SPB(5,1))+(SPA(2,3)*SPB(4,1)*SPB(4,2))/
 (SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,4))))/complex<T>(2,0)+
 (complex<T>(0,-1)*SPA(2,3)*SPB(4,2)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(3,2),2)*
 SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmmpqpgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmmpqpgam L");
#endif
 
 return( complex<T>(0,-1)*((complex<T>(1,0)*pow(SPA(2,5),2))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*
SPB(2,1))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*
pow(SPA(1,2),2)*pow(SPB(4,1),2)*SPB(3,1))/
 (complex<T>(2,0)*pow(SPA(2,3),2)*SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmpmqpgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmpmqpgam L");
#endif
 
 return( complex<T>(0,-1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
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
  
  
 
 
template <class T> complex<T> R2q2g1y_qmppqpgam_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, qp, gam}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmppqpgam L");
#endif
 
 return( complex<T>(0,-1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*
pow(SPA(3,5),2)*pow(SPB(3,2),2))/(complex<T>(2,0)*pow(SPB(2,1),2)*
SPA(1,2)*SPA(2,3)*SPA(3,4))+(complex<T>(1,0)*SPB(3,2)*SPB(4,1)*
SPB(4,2))/(complex<T>(2,0)*SPA(2,3)*SPB(2,1)*SPB(5,1)*SPB(5,4))+
(complex<T>(1,0)*SPB(4,2)*SPB(4,3))/(complex<T>(3,0)*SPA(2,3)*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmmqpmgam_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, m, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmmqpmgam SLC");
#endif
 
 return( (complex<T>(0,-1)*((SPA(4,5)*SPB(3,1)*SPB(4,3))/(SPB(2,1)*SPB(3,2)*
  SPB(4,1)*SPB(5,4))-(pow(SPB(3,1),2)*SPA(1,2))/
 (SPB(2,1)*SPB(4,1)*SPB(5,3)*SPB(5,4))))/complex<T>(2,0)+
 (complex<T>(0,-1)*((pow(SPB(3,1),2)*SPA(1,2))/(SPB(2,1)*SPB(4,3)*SPB(5,1)*
  SPB(5,4))+(SPA(4,5)*SPB(3,1)*SPB(5,3))/(SPB(2,1)*SPB(3,2)*
  SPB(5,1)*SPB(5,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmmqppgam_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, p, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmmqppgam SLC");
#endif
 
 return( complex<T>(0,-1)*((complex<T>(1,0)*pow(SPA(2,5),2))/(complex<T>(2,0)*SPA(3,4)*SPA(4,5)*
 SPB(5,1))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*
 pow(SPA(1,5),2)*pow(SPB(3,1),2)*SPB(4,1))/
(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(2,1)*SPB(3,2)*SPB(5,1)*
 SPB(5,4)))+complex<T>(0,-1)*(-((pow(SPA(1,5),3)*pow(SPB(3,1),2))/
 (complex<T>(2,0)*pow(SPB(3,2),2)*SPA(1,4)*SPA(2,3)*SPA(4,5)*
  SPB(2,1)))-(pow(complex<T>(1,0)-S(2,3)/S(1,4),-1)*pow(SPA(1,5),2)*
 SPB(3,1)*SPB(4,3))/(complex<T>(2,0)*pow(SPA(1,4),2)*SPB(2,1)*SPB(3,2)*
 SPB(4,1))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,4),-1)*
 pow(SPA(2,5),2)*SPB(4,2)*SPB(4,3))/(complex<T>(2,0)*pow(SPA(1,4),2)*
 pow(SPB(4,1),2)*SPB(2,1))-(SPA(2,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(2,0)*S(1,4)*SPB(2,1)*SPB(5,2))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(2,3),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
 pow(SPA(1,5),3)*pow(SPB(3,1),2)*SPB(4,1)*SPB(5,4))/
(pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,4)/S(2,3)-
  S(4,5)/S(2,3))*SPB(2,1))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,4),-1))/complex<T>(2,0))*
 pow(SPA(2,5),3)*SPB(3,2)*SPB(4,2)*SPB(5,4))/
(pow(SPA(1,4),3)*pow(SPB(4,1),3)*(complex<T>(1,0)-S(2,3)/S(1,4)-
  S(3,5)/S(1,4))*SPB(2,1))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,4),-1)*SPA(2,5)*SPA(3,5)*
 SPB(3,2)*SPB(4,3)*SPB(5,4))/(complex<T>(2,0)*pow(SPA(1,4),2)*
 pow(SPB(4,1),2)*SPB(2,1)*SPB(5,2)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmpqpmgam_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, m, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmpqpmgam SLC");
#endif
 
 return( complex<T>(0,-1)*(-((SPA(1,4)*SPA(1,5))/(complex<T>(3,0)*SPA(1,2)*SPA(2,3)*
  SPB(5,4)))-(SPA(1,3)*SPA(1,4)*SPA(4,5))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(5,4))-
 (pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(4,5),2)*
 pow(SPB(5,2),2))/(complex<T>(2,0)*pow(SPA(3,4),2)*SPB(4,3)*SPB(5,1)*
 SPB(5,4)))+complex<T>(0,-1)*((complex<T>(1,0)*SPA(1,4)*SPA(1,5))/
(complex<T>(3,0)*SPA(1,2)*SPA(2,3)*SPB(5,4))-
 (SPA(1,3)*SPA(1,5)*SPA(4,5))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,5)*
 SPB(5,4))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
 pow(SPA(4,5),2)*pow(SPB(4,2),2))/(complex<T>(2,0)*pow(SPA(3,5),2)*
 SPB(4,1)*SPB(5,3)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmpqppgam_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, p, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmpqppgam SLC");
#endif
 
 return( complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(1,3),2)*
  pow(SPB(4,3),2)*SPA(3,5))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(1,2)*
  SPA(2,3)*SPA(3,4)*SPA(4,5)))-pow(SPB(4,2),2)/
(complex<T>(2,0)*SPA(3,4)*SPB(5,1)*SPB(5,4)))+
 complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(1,4)/S(3,5),-1)*pow(SPB(4,2),2)*
  SPA(1,5)*SPA(2,5))/(complex<T>(2,0)*pow(SPA(3,5),2)*pow(SPB(5,3),2)*
  SPA(2,3)))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
 pow(SPB(4,3),2)*SPA(1,3)*SPA(1,5))/(complex<T>(2,0)*pow(SPB(5,3),2)*
 SPA(1,2)*SPA(2,3)*SPA(3,5))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(3,5),-1))/complex<T>(2,0))*
 pow(SPB(4,2),3)*SPA(1,2)*SPA(2,5)*SPA(4,5))/
(pow(SPA(3,5),3)*pow(SPB(5,3),3)*(complex<T>(1,0)-S(1,2)/S(3,5)-
  S(1,4)/S(3,5))*SPA(2,3))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPA(1,3),2)*pow(SPB(4,3),3)*SPA(3,5)*SPA(4,5))/
(pow(SPA(1,2),4)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,5)/S(1,2)-
  S(4,5)/S(1,2))*SPA(2,3))+(complex<T>(1,0)*SPA(1,5)*SPA(2,5)*
 SPB(4,2))/(complex<T>(2,0)*S(3,5)*SPA(2,3)*SPA(2,4))-
 (pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*SPA(1,2)*SPA(1,5)*SPA(4,5)*
 SPB(4,1)*SPB(4,2))/(complex<T>(2,0)*pow(SPA(3,5),2)*pow(SPB(5,3),2)*
 SPA(2,3)*SPA(2,4))+(complex<T>(1,0)*pow(SPA(1,3),2)*pow(SPB(4,3),3))/
(complex<T>(2,0)*pow(SPA(1,2),2)*SPA(2,3)*SPB(2,1)*SPB(5,3)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmqpmmgam_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, m, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmqpmmgam SLC");
#endif
 
 return( (complex<T>(0,-1)*(-(SPA(1,4)*SPB(2,1)*SPB(4,2))+SPA(3,5)*SPB(3,2)*
 SPB(5,2)))/(complex<T>(2,0)*SPB(4,1)*SPB(4,3)*SPB(5,2)*SPB(5,3))+
 (complex<T>(0,1)*(-(SPA(3,4)*SPB(3,2)*SPB(4,2))-SPA(1,5)*SPB(2,1)*
 SPB(5,2)))/(complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*(-(SPA(1,4)*SPB(2,1)*SPB(4,2))-SPA(3,5)*SPB(3,2)*
 SPB(5,2)))/(complex<T>(2,0)*SPB(3,2)*SPB(4,1)*SPB(5,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmqpmpgam_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, p, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmqpmpgam SLC");
#endif
 
 return( (complex<T>(0,1)*((pow(SPA(1,5),2)*SPA(1,3))/(SPA(1,2)*SPA(1,4)*
  SPA(4,5)*SPB(5,3))-(SPA(1,5)*SPA(3,5)*SPB(4,2))/
 (S(1,2)*SPA(4,5)*SPB(5,3))-(pow(SPA(3,5),2)*
  (-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*SPA(1,3)*SPB(3,2)*
  SPB(4,3))/(pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),3)*
  pow(SPB(5,4),2)*SPB(5,3))-(pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*
  pow(SPA(3,5),2)*SPB(3,2)*SPB(4,2))/(pow(SPA(4,5),2)*SPB(2,1)*
  SPB(5,3)*SPB(5,4))))/complex<T>(3,0)+
 (complex<T>(0,-1)*pow(SPB(4,2),3)*SPB(4,1))/(SPB(2,1)*SPB(3,2)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+complex<T>(0,1)*
((complex<T>(1,0)*pow(SPA(3,5),3)*(-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*
 SPB(3,2)*SPB(4,2))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*
 pow(SPA(4,5),3)*pow(SPB(5,4),2)*SPB(2,1))+
 (complex<T>(1,0)*pow(SPA(3,5),4)*SPB(3,2)*SPB(5,2))/
(complex<T>(3,0)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(3,5),3)*(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*
 SPB(4,2)*SPB(5,2))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPB(2,1))+
 (complex<T>(-2,0)*pow(SPA(3,5),4)*(-(S(1,2)/(complex<T>(6,0)*S(3,4)))+
  (complex<T>(1,0)*(S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPB(3,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,4))/(pow(SPA(1,2),4)*pow(SPB(2,1),5)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))))+
 complex<T>(0,-1)*((complex<T>(1,0)*pow(SPA(3,5),3)*(-(S(1,2)/S(4,5))+
  S(4,5)/S(1,2))*SPB(3,2)*SPB(4,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),3)*
 pow(SPB(5,4),2)*SPB(2,1))-
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
  S(4,5)/S(1,2))*SPB(5,1)))+
 (complex<T>(0,1)*((pow(SPA(1,3),2)*SPA(1,5))/(SPA(1,2)*SPA(1,4)*SPA(3,4)*
  SPB(5,3))+(SPA(1,3)*SPA(3,5)*SPB(4,2))/(S(1,2)*SPA(3,4)*
  SPB(5,3))-(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(3,5),2)*
  SPB(4,2)*SPB(5,2))/(pow(SPA(3,4),2)*SPB(2,1)*SPB(4,3)*
  SPB(5,3))+(pow(SPA(3,5),2)*(-(S(1,2)/S(3,4))+
   S(3,4)/S(1,2))*SPA(1,5)*SPB(5,2)*SPB(5,4))/
 (pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),3)*
  pow(SPB(4,3),2)*SPB(5,3))))/complex<T>(3,0)+
 complex<T>(0,-1)*(-((pow(SPB(4,2),2)*SPA(1,2)*SPA(1,3)*SPB(4,1))/
 (complex<T>(2,0)*pow(SPB(5,4),2)*S(2,3)*SPA(4,5)*SPB(5,3)))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1))/complex<T>(2,0))*
 pow(SPA(1,3),3)*pow(SPB(4,3),2)*SPB(2,1)*SPB(4,1))/
(pow(SPA(4,5),3)*pow(SPB(5,4),4)*(complex<T>(1,0)-S(1,2)/S(4,5)-
  S(2,3)/S(4,5))*SPB(5,3))+(SPA(1,3)*SPA(3,5)*SPB(4,1)*
 SPB(4,2))/(S(4,5)*SPA(2,3)*SPB(2,1)*SPB(5,3))-
 (pow(SPA(1,5),2)*SPA(3,5)*SPB(5,2))/(complex<T>(2,0)*S(2,3)*SPA(1,4)*
 SPA(4,5)*SPB(5,3))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(2,3),-1)*
 pow(SPA(1,5),2)*SPB(4,1)*SPB(4,2)*SPB(5,2))/
(complex<T>(2,0)*pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPB(5,3))+
 (complex<T>(1,0)*((pow(SPA(1,5),2)*SPA(1,3))/(SPA(1,2)*SPA(1,4)*SPA(4,5)*
SPB(5,3))-(SPA(1,5)*SPA(3,5)*SPB(4,2))/(S(1,2)*SPA(4,5)*
SPB(5,3))-(pow(SPA(3,5),2)*(-(S(1,2)/S(4,5))+
 S(4,5)/S(1,2))*SPA(1,3)*SPB(3,2)*SPB(4,3))/
   (pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),3)*
pow(SPB(5,4),2)*SPB(5,3))-(pow(complex<T>(1,0)-S(1,2)/S(4,5),
-1)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,2))/(pow(SPA(4,5),2)*
SPB(2,1)*SPB(5,3)*SPB(5,4))))/complex<T>(3,0)-
 (pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*
 ((complex<T>(1,0)*pow(SPB(4,2),2)*SPA(1,2)*SPA(1,3)*SPB(4,1))/
   (complex<T>(2,0)*SPB(5,3)*SPB(5,4))+(pow(SPA(1,3),2)*
pow(SPB(4,1),2)*SPB(3,2)*SPB(4,2))/(SPB(2,1)*SPB(5,3)*
SPB(5,4))))/(pow(SPA(2,3),2)*pow(SPB(3,2),2))-
 (pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(3,5),2)*SPB(3,2)*
 SPB(4,2))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(2,1)*SPB(5,3)*
 SPB(5,4))+(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(2,3),-1))/
   complex<T>(2,0)+(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
 pow(SPA(1,5),3)*SPB(2,1)*SPB(4,1)*SPB(5,2)*SPB(5,4))/
(pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,4)/S(2,3)-
  S(4,5)/S(2,3))*SPB(5,3)))+
 complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(3,4)/S(2,5),-1)*
  ((complex<T>(1,0)*pow(SPB(4,2),2)*SPA(1,2)*SPA(1,5)*SPB(4,1))/
(complex<T>(2,0)*SPB(4,3)*SPB(5,3))+(pow(SPA(1,5),2)*
 pow(SPB(4,1),2)*SPB(4,2)*SPB(5,2))/(SPB(2,1)*SPB(4,3)*
 SPB(5,3))))/(pow(SPA(2,5),2)*pow(SPB(5,2),2)))+
 (complex<T>(1,0)*pow(SPA(1,3),2)*SPA(3,5)*SPB(3,2))/
(complex<T>(2,0)*S(2,5)*SPA(1,4)*SPA(3,4)*SPB(5,3))-
 (pow(SPB(4,2),2)*SPA(1,2)*SPA(1,5)*SPB(4,1))/
(complex<T>(2,0)*pow(SPB(4,3),2)*S(2,5)*SPA(3,4)*SPB(5,3))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(1,5),3)*pow(SPB(5,4),2)*SPB(2,1)*SPB(4,1))/
(pow(SPA(3,4),3)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(2,5)/S(3,4))*SPB(5,3))+(SPA(1,5)*SPA(3,5)*SPB(4,1)*
 SPB(4,2))/(S(3,4)*SPA(2,5)*SPB(2,1)*SPB(5,3))-
 (pow(complex<T>(1,0)-S(1,4)/S(2,5),-1)*pow(SPA(1,3),2)*SPB(3,2)*
 SPB(4,1)*SPB(4,2))/(complex<T>(2,0)*pow(SPA(2,5),2)*pow(SPB(5,2),3)*
 SPB(5,3))+(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(2,5),-1))/
   complex<T>(2,0)+(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(2,5),-1))/complex<T>(2,0))*
 pow(SPA(1,3),3)*SPB(2,1)*SPB(3,2)*SPB(4,1)*SPB(4,3))/
(pow(SPA(2,5),3)*pow(SPB(5,2),4)*(complex<T>(1,0)-S(1,4)/S(2,5)-
  S(3,4)/S(2,5))*SPB(5,3))-(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
 pow(SPA(3,5),2)*SPB(4,2)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 SPB(2,1)*SPB(4,3)*SPB(5,3))+
 (complex<T>(1,0)*((pow(SPA(1,3),2)*SPA(1,5))/(SPA(1,2)*SPA(1,4)*SPA(3,4)*
SPB(5,3))+(SPA(1,3)*SPA(3,5)*SPB(4,2))/(S(1,2)*SPA(3,4)*
SPB(5,3))-(pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
pow(SPA(3,5),2)*SPB(4,2)*SPB(5,2))/(pow(SPA(3,4),2)*
SPB(2,1)*SPB(4,3)*SPB(5,3))+(pow(SPA(3,5),2)*
(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPA(1,5)*SPB(5,2)*
SPB(5,4))/(pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),3)*
pow(SPB(4,3),2)*SPB(5,3))))/complex<T>(3,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmqppmgam_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, m, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmqppmgam SLC");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(3,2),3)*SPB(3,1))/(SPB(2,1)*SPB(4,1)*SPB(4,3)*
 SPB(5,2)*SPB(5,3))+complex<T>(0,1)*
((complex<T>(1,0)*pow(SPA(4,5),3)*(-(S(1,2)/S(3,5))+S(3,5)/S(1,2))*
 SPB(3,2)*SPB(4,2))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*
 pow(SPA(3,5),3)*pow(SPB(5,3),2)*SPB(2,1))-
 (pow(SPA(4,5),3)*(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPB(3,2)*
 SPB(5,2))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPB(2,1))-
 (pow(SPA(4,5),4)*SPB(4,2)*SPB(5,2))/(complex<T>(3,0)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(3,4)*SPA(3,5))+
 (complex<T>(2,0)*pow(SPA(4,5),4)*(-(S(1,2)/(complex<T>(6,0)*S(3,4)))+
  (complex<T>(1,0)*(S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(3,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(3,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),3)))*SPB(4,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,3))/(pow(SPA(1,2),4)*pow(SPB(2,1),5)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(3,5)/S(1,2))))+
 complex<T>(0,-1)*((complex<T>(1,0)*pow(SPA(4,5),3)*(-(S(1,2)/S(3,5))+
  S(3,5)/S(1,2))*SPB(3,2)*SPB(4,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*pow(SPA(3,5),3)*
 pow(SPB(5,3),2)*SPB(2,1))-
 (pow(SPA(4,5),3)*(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPB(3,2)*
 SPB(5,2))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPB(2,1))-
 (pow(SPA(4,5),4)*SPB(4,2)*SPB(5,2))/(complex<T>(3,0)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(3,4)*SPA(3,5))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(1,5),3)*pow(SPB(3,1),2)*SPB(2,1)*SPB(5,3))/
(pow(SPA(3,4),3)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(2,5)/S(3,4))*SPB(4,1))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPA(4,5),3)*SPB(4,2)*SPB(4,3)*SPB(5,3))/
(pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
  S(3,5)/S(1,2))*SPB(4,1))+(complex<T>(2,0)*pow(SPA(4,5),4)*
 (-(S(1,2)/(complex<T>(6,0)*S(3,4)))+(complex<T>(1,0)*(S(1,2)/S(3,4)-
 S(3,4)/S(1,2)))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(3,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(3,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),3)))*SPB(4,2)*SPB(4,3)*
 SPB(5,2)*SPB(5,3))/(pow(SPA(1,2),4)*pow(SPB(2,1),5)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(3,5)/S(1,2)))+
 (complex<T>(1,0)*((pow(SPA(1,5),2)*SPA(4,5)*SPB(3,1))/(S(3,4)*SPA(1,2)*
SPA(3,5)*SPB(4,1))-(pow(SPA(1,5),2)*pow(SPB(3,1),2)*
SPB(3,2))/(pow(SPB(4,3),2)*S(2,5)*SPA(3,4)*SPB(4,1))+
  (pow(SPA(4,5),2)*SPA(2,5)*SPB(3,2)*SPB(4,2))/
   (S(1,2)*S(3,4)*SPA(3,5)*SPB(4,1))-
  (pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(4,5),2)*SPB(3,2))/
   (pow(SPA(3,4),2)*SPB(4,1)*SPB(4,3))-
  (pow(complex<T>(1,0)-S(3,4)/S(2,5),-1)*pow(SPA(1,5),2)*
pow(SPB(3,1),2)*SPB(3,2))/(pow(SPA(2,5),2)*pow(SPB(5,2),2)*
SPB(4,1)*SPB(4,3))+(pow(SPB(3,2),2)*SPA(1,4)*SPB(3,1))/
   (S(1,2)*SPB(4,1)*SPB(5,2)*SPB(5,3))+
  (pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*pow(SPA(1,5),2)*SPB(3,1)*
SPB(3,2)*SPB(5,3))/(pow(SPA(1,2),2)*SPB(2,1)*SPB(4,1)*
SPB(4,3)*SPB(5,2))))/complex<T>(2,0))+
 complex<T>(0,-1)*((complex<T>(1,0)*(-((SPA(1,5)*SPA(4,5))/(SPA(2,3)*SPA(3,5)*
 SPB(5,4)))-(pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
pow(SPA(4,5),2)*SPB(4,2)*SPB(4,3))/(pow(SPA(3,5),2)*
SPB(4,1)*SPB(5,3)*SPB(5,4))))/complex<T>(2,0)+
 (complex<T>(1,0)*(-((pow(SPA(4,5),3)*pow(SPB(4,2),2)*pow(SPB(5,3),2)*
 (S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
(pow(complex<T>(1,0)-S(3,5)/S(1,2),3)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4)*SPB(5,4)))+(complex<T>(-2,0)*pow(SPB(3,2),2)*
SPA(4,5))/(pow(SPB(2,1),2)*SPA(1,2)*SPB(5,4))-
  (S(2,5)*SPA(1,5)*SPA(4,5))/(S(1,2)*SPA(2,3)*SPA(3,5)*
SPB(5,4))+(SPA(1,5)*SPA(2,4)*SPB(3,2))/(S(1,2)*SPA(2,3)*
SPB(5,4))+(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1)*
pow(SPA(4,5),2)*SPB(3,2)*SPB(4,2)*SPB(5,3))/
   (pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(5,4))))/complex<T>(3,0))+
 complex<T>(0,-1)*((complex<T>(1,0)*(-((pow(SPA(4,5),3)*pow(SPB(4,3),2)*
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
complex<T>(2,0))+(complex<T>(0,1)*(-((pow(SPA(4,5),3)*pow(SPB(4,3),2)*
   pow(SPB(5,2),2)*(S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
  (pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPA(1,2),3)*
   pow(SPB(2,1),4)*SPB(5,4)))+(complex<T>(-2,0)*pow(SPB(3,2),2)*
  SPA(4,5))/(pow(SPB(2,1),2)*SPA(1,2)*SPB(5,4))-
(S(2,4)*SPA(1,4)*SPA(4,5))/(S(1,2)*SPA(2,3)*SPA(3,4)*
  SPB(5,4))-(SPA(1,4)*SPA(2,5)*SPB(3,2))/(S(1,2)*SPA(2,3)*
  SPB(5,4))+(complex<T>(-3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
  pow(SPA(4,5),2)*SPB(3,2)*SPB(4,3)*SPB(5,2))/
 (pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(5,4))))/complex<T>(3,0)+
 (complex<T>(0,1)*(-((pow(SPA(4,5),3)*pow(SPB(4,2),2)*pow(SPB(5,3),2)*
   (S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
  (pow(complex<T>(1,0)-S(3,5)/S(1,2),3)*pow(SPA(1,2),3)*
   pow(SPB(2,1),4)*SPB(5,4)))+(complex<T>(-2,0)*pow(SPB(3,2),2)*
  SPA(4,5))/(pow(SPB(2,1),2)*SPA(1,2)*SPB(5,4))-
(S(2,5)*SPA(1,5)*SPA(4,5))/(S(1,2)*SPA(2,3)*SPA(3,5)*
  SPB(5,4))+(SPA(1,5)*SPA(2,4)*SPB(3,2))/(S(1,2)*SPA(2,3)*
  SPB(5,4))+(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1)*
  pow(SPA(4,5),2)*SPB(3,2)*SPB(4,2)*SPB(5,3))/
 (pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(5,4))))/complex<T>(3,0)+
 (complex<T>(0,-1)*pow(SPB(3,2),2)*SPB(3,1))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*pow(SPB(3,2),2)*SPB(3,1))/
(complex<T>(2,0)*SPB(2,1)*SPB(4,1)*SPB(5,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmqpppgam_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, p, gam}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmqpppgam SLC");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,5),2)*SPA(2,5))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*
 SPA(3,4)*SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,5),3)*SPA(2,5))/
(SPA(1,2)*SPA(1,4)*SPA(2,3)*SPA(3,5)*SPA(4,5))+
 (complex<T>(0,1)*((pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*pow(SPB(4,3),2)*
  SPA(1,4)*SPA(1,5))/(pow(SPB(5,3),2)*SPA(1,2)*SPA(3,4)*
  SPA(3,5))-(pow(SPB(4,3),2)*(-(S(1,2)/S(3,5))+
   S(3,5)/S(1,2))*SPA(1,4)*SPA(4,5)*SPB(4,2))/
 (pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*pow(SPA(3,5),2)*
  pow(SPB(5,3),3)*SPA(3,4))+(SPA(1,5)*SPB(3,2)*SPB(4,3))/
 (S(1,2)*SPA(3,4)*SPB(5,3))-(pow(SPB(3,2),2)*SPB(4,2))/
 (SPA(3,4)*SPB(2,1)*SPB(5,2)*SPB(5,3))))/complex<T>(3,0)+
 complex<T>(0,-1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1))/complex<T>(2,0)+
   (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(3,5),-1))/complex<T>(2,0))*
  pow(SPA(4,5),2)*pow(SPB(4,2),3)*SPA(1,2)*SPA(2,5))/
 (pow(SPA(3,5),4)*pow(SPB(5,3),3)*(complex<T>(1,0)-S(1,2)/S(3,5)-
   S(1,4)/S(3,5))*SPA(3,4)))-(pow(complex<T>(1,0)-S(2,5)/S(1,4),-1)*
 pow(SPB(3,2),2)*SPA(1,3)*SPA(1,5)*SPA(2,5))/
(complex<T>(2,0)*pow(SPA(1,4),3)*pow(SPB(4,1),2)*SPA(3,4))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*pow(SPB(4,3),2)*
 SPA(1,4)*SPA(1,5))/(complex<T>(2,0)*pow(SPB(5,3),2)*SPA(1,2)*SPA(3,4)*
 SPA(3,5))-(((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(1,4),-1))/
   complex<T>(2,0)+(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,4),-1))/complex<T>(2,0))*
 pow(SPB(3,2),3)*SPA(1,2)*SPA(1,3)*SPA(2,5)*SPA(3,5))/
(pow(SPA(1,4),4)*pow(SPB(4,1),3)*(complex<T>(1,0)-S(2,5)/S(1,4)-
  S(3,5)/S(1,4))*SPA(3,4))-(pow(complex<T>(1,0)-S(3,5)/S(1,4),-1)*
 (-((pow(SPA(2,5),2)*pow(SPB(4,2),2)*SPA(1,4)*SPA(1,5))/
(SPA(1,2)*SPA(3,4)*SPA(3,5)))-(pow(SPA(1,5),2)*SPA(2,5)*
SPB(2,1)*SPB(4,2))/(complex<T>(2,0)*SPA(3,4)*SPA(3,5))))/
(pow(SPA(1,4),2)*pow(SPB(4,1),2))-
 (SPA(1,5)*SPA(2,5)*SPB(4,2)*SPB(4,3))/(S(3,5)*SPA(1,2)*SPA(3,4)*
 SPB(4,1))+(complex<T>(1,0)*((pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
pow(SPB(4,3),2)*SPA(1,4)*SPA(1,5))/(pow(SPB(5,3),2)*
SPA(1,2)*SPA(3,4)*SPA(3,5))-(pow(SPB(4,3),2)*
(-(S(1,2)/S(3,5))+S(3,5)/S(1,2))*SPA(1,4)*SPA(4,5)*
SPB(4,2))/(pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*pow(SPA(3,5),2)*
pow(SPB(5,3),3)*SPA(3,4))+(SPA(1,5)*SPB(3,2)*SPB(4,3))/
   (S(1,2)*SPA(3,4)*SPB(5,3))-(pow(SPB(3,2),2)*SPB(4,2))/
   (SPA(3,4)*SPB(2,1)*SPB(5,2)*SPB(5,3))))/complex<T>(3,0)+
 (complex<T>(1,0)*pow(SPA(1,5),2)*SPA(2,5)*SPB(2,1)*SPB(4,2))/
(complex<T>(2,0)*pow(SPA(3,5),2)*S(1,4)*SPA(3,4)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(1,3)*SPB(4,3))/
(complex<T>(2,0)*S(1,4)*SPA(3,4)*SPB(5,2)*SPB(5,3)))+
 complex<T>(0,-1)*((complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*
 pow(SPB(4,3),2)*SPA(1,3)*SPA(3,5))/(pow(SPB(5,4),2)*
 SPA(2,3)*SPA(3,4)*SPA(4,5)))+(SPB(4,2)*SPB(4,3))/
   (SPA(3,4)*SPB(5,1)*SPB(5,4))))/complex<T>(2,0)+
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
SPB(5,4))))/complex<T>(3,0))+
 (complex<T>(0,1)*((pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(4,3),3)*
  (S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
 (pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),4)*
  pow(SPB(2,1),3)*SPA(3,4))+
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*pow(SPB(4,3),2)*
  SPA(1,3)*SPA(1,5)*SPA(4,5))/(pow(SPA(1,2),3)*pow(SPB(2,1),2)*
  SPA(3,4))+(complex<T>(2,0)*pow(SPA(1,5),2)*SPB(4,3))/
 (pow(SPA(1,2),2)*SPA(3,4)*SPB(2,1))+
(SPA(1,5)*SPB(3,1)*SPB(4,2))/(S(1,2)*SPA(3,4)*SPB(5,1))+
(S(1,4)*SPB(4,2)*SPB(4,3))/(S(1,2)*SPA(3,4)*SPB(5,1)*
  SPB(5,4))))/complex<T>(3,0)+
 complex<T>(0,1)*((complex<T>(1,0)*pow(SPB(4,3),3)*(-(S(1,2)/S(4,5))+
  S(4,5)/S(1,2))*SPA(1,3)*SPA(1,5))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(1,2))-
 (pow(SPB(4,3),3)*(-(S(1,2)/S(3,5))+S(3,5)/S(1,2))*SPA(1,4)*
 SPA(1,5))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*
 pow(SPA(3,5),2)*pow(SPB(5,3),3)*SPA(1,2))+
 (complex<T>(-2,0)*pow(SPB(4,3),4)*(-(S(1,2)/(complex<T>(6,0)*S(3,5)))+
  (complex<T>(1,0)*(S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,5)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPA(1,3)*SPA(1,4)*
 SPA(3,5)*SPA(4,5))/(pow(SPA(1,2),5)*pow(SPB(2,1),4)*
 (complex<T>(1,0)-S(3,5)/S(1,2)-S(4,5)/S(1,2)))+
 (complex<T>(1,0)*pow(SPB(4,3),4)*SPA(1,3)*SPA(1,4))/
(complex<T>(3,0)*pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPB(5,3)*SPB(5,4)))+
 complex<T>(0,-1)*((complex<T>(1,0)*pow(SPB(4,3),3)*(-(S(1,2)/S(4,5))+
  S(4,5)/S(1,2))*SPA(1,3)*SPA(1,5))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(1,2))-
 (pow(SPB(4,3),3)*(-(S(1,2)/S(3,5))+S(3,5)/S(1,2))*SPA(1,4)*
 SPA(1,5))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*
 pow(SPA(3,5),2)*pow(SPB(5,3),3)*SPA(1,2))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(3,5),-1))/complex<T>(2,0))*
 pow(SPA(2,5),2)*pow(SPB(4,2),3)*SPA(1,2)*SPA(4,5))/
(pow(SPA(3,5),4)*pow(SPB(5,3),3)*(complex<T>(1,0)-S(1,2)/S(3,5)-
  S(1,4)/S(3,5))*SPA(2,3))+(complex<T>(-2,0)*pow(SPB(4,3),4)*
 (-(S(1,2)/(complex<T>(6,0)*S(3,5)))+(complex<T>(1,0)*(S(1,2)/S(3,5)-
 S(3,5)/S(1,2)))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,5)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPA(1,3)*SPA(1,4)*
 SPA(3,5)*SPA(4,5))/(pow(SPA(1,2),5)*pow(SPB(2,1),4)*
 (complex<T>(1,0)-S(3,5)/S(1,2)-S(4,5)/S(1,2)))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPB(4,3),3)*SPA(1,3)*SPA(3,5)*SPA(4,5))/
(pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,5)/S(1,2)-
  S(4,5)/S(1,2))*SPA(2,3))+
 (complex<T>(1,0)*((pow(complex<T>(1,0)-S(3,5)/S(1,4),-1)*pow(SPA(2,5),2)*
pow(SPB(4,2),2)*SPA(1,5))/(pow(SPA(1,4),2)*pow(SPB(4,1),2)*
SPA(2,3)*SPA(3,5))+(pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
pow(SPB(4,3),2)*SPA(1,5))/(pow(SPB(5,3),2)*SPA(2,3)*
SPA(3,5))-(pow(complex<T>(1,0)-S(3,5)/S(1,2),-1)*
pow(SPB(4,2),2)*SPA(1,5)*SPA(2,5)*SPA(4,5))/
   (pow(SPB(2,1),2)*SPA(1,2)*SPA(1,4)*SPA(2,3)*SPA(3,5))-
  (pow(SPA(1,5),2)*SPA(2,5)*SPB(3,2))/(S(1,2)*SPA(1,4)*
SPA(2,3)*SPA(4,5))+(pow(SPA(2,5),2)*pow(SPB(4,2),2)*
SPA(1,5))/(pow(SPA(3,5),2)*S(1,4)*SPA(2,3)*SPB(5,3))-
  (pow(SPB(4,3),2)*SPA(1,3)*SPA(1,5)*SPB(4,1))/
   (S(1,2)*S(3,5)*SPA(2,3)*SPB(5,4))+
  (pow(SPB(4,2),2)*SPA(2,5)*SPB(4,3))/(S(3,5)*SPA(2,3)*
SPB(2,1)*SPB(5,4))))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(SPB(4,3),4)*SPA(1,3)*SPA(1,4))/
(complex<T>(3,0)*pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPB(5,3)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmmmqpgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, m, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmmmqpgam nf");
#endif
 
 return( (complex<T>(0,1)*SPA(2,3)*SPB(4,2)*SPB(4,3))/(complex<T>(3,0)*pow(SPB(3,2),2)*
SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmmpqpgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, p, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmmpqpgam nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1y_qmpmqpgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, m, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmpmqpgam nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1y_qmppqpgam_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, p, qp, gam}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmppqpgam nf");
#endif
 
 return( (complex<T>(0,1)*SPB(4,2)*SPB(4,3))/(complex<T>(3,0)*SPA(2,3)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmgamqpmm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmgamqpmm L");
#endif
 
 return( (complex<T>(0,1)*SPA(4,5)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*
 SPB(2,1)*SPB(3,2))+
 (complex<T>(0,1)*((pow(SPB(3,1),2)*SPA(1,2))/(SPB(2,1)*SPB(4,3)*SPB(5,1)*
  SPB(5,4))+(SPA(4,5)*SPB(3,1)*SPB(5,3))/(SPB(2,1)*SPB(3,2)*
  SPB(5,1)*SPB(5,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmgamqpmp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmgamqpmp L");
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
  
  
 
 
template <class T> complex<T> R2q2g1y_qmgamqppm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmgamqppm L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(2,5),2))/(complex<T>(2,0)*SPA(3,4)*SPA(4,5)*
SPB(5,1))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*
pow(SPA(1,5),2)*pow(SPB(3,1),2)*SPB(4,1))/
 (complex<T>(2,0)*pow(SPA(4,5),2)*SPB(2,1)*SPB(3,2)*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmgamqppp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmgamqppp L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*
pow(SPA(2,4),2)*pow(SPB(5,4),2))/(complex<T>(2,0)*pow(SPB(5,1),2)*
SPA(1,5)*SPA(3,4)*SPA(4,5))+(complex<T>(1,0)*SPB(4,3)*SPB(5,3))/
 (complex<T>(3,0)*SPA(4,5)*SPB(2,1)*SPB(3,2))+
(complex<T>(1,0)*SPB(3,1)*SPB(5,3)*SPB(5,4))/(complex<T>(2,0)*SPA(4,5)*SPB(2,1)*
SPB(3,2)*SPB(5,1)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmgamqpmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmgamqpmm nf");
#endif
 
 return( (complex<T>(0,-1)*SPA(4,5)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*
SPB(2,1)*SPB(3,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qmgamqpmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmgamqpmp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1y_qmgamqppm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmgamqppm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1y_qmgamqppp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, gam, qp, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qmgamqppp nf");
#endif
 
 return( (complex<T>(0,-1)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*SPA(4,5)*SPB(2,1)*
SPB(3,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpppqmgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, p, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpppqmgap L");
#endif
 
 return( (complex<T>(0,1)*SPA(2,4)*SPA(3,4)*SPB(3,2))/(complex<T>(3,0)*pow(SPA(2,3),2)*
 SPA(1,5)*SPA(4,5))+
 (complex<T>(0,1)*((SPA(1,4)*SPA(2,4)*SPB(3,2))/(SPA(1,2)*SPA(1,5)*SPA(2,3)*
  SPA(4,5))+(pow(SPA(1,4),2)*SPB(5,1))/(SPA(1,2)*SPA(1,5)*
  SPA(2,3)*SPA(3,4))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qppmqmgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, m, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qppmqmgap L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*
pow(SPA(1,4),2)*pow(SPB(2,1),2)*SPA(1,3))/
 (complex<T>(2,0)*pow(SPB(3,2),2)*SPA(1,2)*SPA(1,5)*SPA(2,3)*SPA(4,5))+
(complex<T>(1,0)*pow(SPB(5,2),2))/(complex<T>(2,0)*SPA(1,2)*SPB(3,2)*SPB(4,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpmpqmgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, p, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpmpqmgap L");
#endif
 
 return( complex<T>(0,1)*(-((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1))/complex<T>(2,0))*
 pow(SPA(1,4),2)*pow(SPB(3,1),3)*SPA(1,2)*SPA(2,3))/
(pow(SPA(4,5),4)*pow(SPB(5,4),3)*(complex<T>(1,0)-S(1,2)/S(4,5)-
  S(2,3)/S(4,5))*SPA(1,5)))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*pow(SPB(5,3),2)*SPA(2,4)*
SPA(2,5))/(complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPA(1,5))-
(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*pow(SPB(3,1),2)*SPA(1,4)*
SPA(2,4))/(complex<T>(2,0)*pow(SPB(2,1),2)*SPA(1,2)*SPA(1,5)*
SPA(4,5))-(((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
pow(SPB(5,3),3)*SPA(2,3)*SPA(2,5)*SPA(4,5))/
 (pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
 S(4,5)/S(1,2))*SPA(1,5))-(SPA(2,4)*SPA(2,5)*SPB(5,3))/
 (complex<T>(2,0)*S(1,2)*SPA(1,5)*SPA(3,5))-
(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*SPA(2,3)*SPA(2,4)*SPA(4,5)*
SPB(4,3)*SPB(5,3))/(complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*
SPA(1,5)*SPA(3,5))+(complex<T>(1,0)*pow(SPA(1,4),2)*pow(SPB(3,1),3))/
 (complex<T>(2,0)*pow(SPA(4,5),2)*SPA(1,5)*SPB(2,1)*SPB(3,2)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpmmqmgap_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, m, qm, gap}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpmmqmgap L");
#endif
 
 return( complex<T>(0,1)*((complex<T>(1,0)*SPA(1,4)*SPA(2,3)*SPA(2,4))/
 (complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(4,5)*SPB(3,2))+
(complex<T>(1,0)*SPA(2,4)*SPA(3,4))/(complex<T>(3,0)*SPA(1,5)*SPA(4,5)*SPB(3,2))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*pow(SPA(2,3),2)*
pow(SPB(5,3),2))/(complex<T>(2,0)*pow(SPA(1,2),2)*SPB(2,1)*SPB(3,2)*
SPB(4,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qppqmpgap_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, qm, p, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qppqmpgap SLC");
#endif
 
 return( (complex<T>(0,1)*(-((pow(SPA(1,3),2)*SPB(2,1))/(SPA(1,2)*SPA(1,4)*
   SPA(3,5)*SPA(4,5)))+(SPA(1,3)*SPA(3,4)*SPB(5,4))/
 (SPA(1,2)*SPA(1,4)*SPA(2,3)*SPA(4,5))))/complex<T>(2,0)+
 (complex<T>(0,1)*((pow(SPA(1,3),2)*SPB(2,1))/(SPA(1,2)*SPA(1,5)*SPA(3,4)*
  SPA(4,5))+(SPA(1,3)*SPA(3,5)*SPB(5,4))/(SPA(1,2)*SPA(1,5)*
  SPA(2,3)*SPA(4,5))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qppqmmgap_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, qm, m, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qppqmmgap SLC");
#endif
 
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-S(2,3)/S(1,4),-1)*pow(SPB(5,1),2)*
  SPA(1,3)*SPA(3,4))/(complex<T>(2,0)*pow(SPB(4,1),2)*SPA(1,2)*SPA(1,4)*
  SPA(2,3)))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,4),-1)*
 pow(SPB(5,2),2)*SPA(2,4)*SPA(3,4))/(complex<T>(2,0)*pow(SPA(1,4),2)*
 pow(SPB(4,1),2)*SPA(1,2))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(2,3),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
 pow(SPA(1,3),2)*pow(SPB(5,1),3)*SPA(1,4)*SPA(4,5))/
(pow(SPA(2,3),4)*pow(SPB(3,2),3)*(complex<T>(1,0)-S(1,4)/S(2,3)-
  S(4,5)/S(2,3))*SPA(1,2))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,4),-1))/complex<T>(2,0))*
 pow(SPB(5,2),3)*SPA(2,3)*SPA(2,4)*SPA(4,5))/
(pow(SPA(1,4),3)*pow(SPB(4,1),3)*(complex<T>(1,0)-S(2,3)/S(1,4)-
  S(3,5)/S(1,4))*SPA(1,2))-(SPA(2,4)*SPA(3,4)*SPB(5,2))/
(complex<T>(2,0)*S(1,4)*SPA(1,2)*SPA(2,5))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,4),-1)*SPA(2,3)*SPA(3,4)*
 SPA(4,5)*SPB(5,2)*SPB(5,3))/(complex<T>(2,0)*pow(SPA(1,4),2)*
 pow(SPB(4,1),2)*SPA(1,2)*SPA(2,5))-
 (pow(SPA(1,3),2)*pow(SPB(5,1),3))/(complex<T>(2,0)*S(2,3)*SPA(1,2)*
 SPA(2,3)*SPB(4,1)*SPB(5,4)))+
 complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*pow(SPA(1,3),2)*
 pow(SPB(5,1),2)*SPA(1,4))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(1,2)*
 SPA(1,5)*SPA(2,3)*SPA(4,5))+(complex<T>(1,0)*pow(SPB(5,2),2))/
(complex<T>(2,0)*SPA(1,5)*SPB(4,3)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpmqmpgap_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, qm, p, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpmqmpgap SLC");
#endif
 
 return( complex<T>(0,1)*(-((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPA(2,5),2)*
  pow(SPB(5,4),2))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(1,5)*SPA(3,4)*
  SPA(4,5)))-(SPB(4,1)*SPB(5,1))/(complex<T>(3,0)*SPA(4,5)*SPB(2,1)*
 SPB(3,2))-(SPB(3,1)*SPB(4,1)*SPB(5,4))/(complex<T>(2,0)*SPA(4,5)*
 SPB(2,1)*SPB(3,2)*SPB(4,3)))+
 complex<T>(0,1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*pow(SPA(2,4),2)*
 pow(SPB(5,4),2))/(complex<T>(2,0)*pow(SPB(5,3),2)*SPA(1,4)*SPA(3,5)*
 SPA(4,5))+(complex<T>(1,0)*SPB(4,1)*SPB(5,1))/(complex<T>(3,0)*SPA(4,5)*
 SPB(2,1)*SPB(3,2))-(SPB(3,1)*SPB(5,1)*SPB(5,4))/
(complex<T>(2,0)*SPA(4,5)*SPB(2,1)*SPB(3,2)*SPB(5,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpmqmmgap_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, qm, m, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpmqmmgap SLC");
#endif
 
 return( complex<T>(0,1)*(-(pow(SPA(2,4),2)/(complex<T>(2,0)*SPA(1,5)*SPA(4,5)*
  SPB(4,3)))-(pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPA(3,4),2)*
 pow(SPB(3,1),2)*SPB(5,3))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(2,1)*
 SPB(3,2)*SPB(4,3)*SPB(5,4)))+
 complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(3,4),3)*pow(SPB(3,1),2))/
(complex<T>(2,0)*S(1,2)*SPA(3,5)*SPA(4,5)*SPB(2,1)*SPB(3,2))-
 (pow(complex<T>(1,0)-S(1,4)/S(3,5),-1)*pow(SPA(2,4),2)*SPB(5,1)*
 SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,5),2)*pow(SPB(5,3),2)*SPB(3,2))+
 (complex<T>(1,0)*SPA(2,4)*SPB(5,1)*SPB(5,2))/(complex<T>(2,0)*S(3,5)*SPB(3,2)*
 SPB(4,2))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
 pow(SPA(3,4),2)*SPB(3,1)*SPB(5,1))/(complex<T>(2,0)*pow(SPA(3,5),2)*
 SPB(2,1)*SPB(3,2)*SPB(5,3))-(pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
 SPA(1,4)*SPA(2,4)*SPB(2,1)*SPB(5,1)*SPB(5,4))/
(complex<T>(2,0)*pow(SPA(3,5),2)*pow(SPB(5,3),2)*SPB(3,2)*SPB(4,2))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(3,5),-1))/complex<T>(2,0))*
 pow(SPA(2,4),3)*SPB(2,1)*SPB(5,2)*SPB(5,4))/
(pow(SPA(3,5),3)*pow(SPB(5,3),3)*(complex<T>(1,0)-S(1,2)/S(3,5)-
  S(1,4)/S(3,5))*SPB(3,2))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPA(3,4),3)*pow(SPB(3,1),2)*SPB(5,3)*SPB(5,4))/
(pow(SPA(1,2),3)*pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,5)/S(1,2)-
  S(4,5)/S(1,2))*SPB(3,2)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpqmppgap_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, p, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpqmppgap SLC");
#endif
 
 return( (complex<T>(0,-1)*(-(SPA(2,3)*SPA(2,4)*SPB(4,3))-SPA(1,2)*SPA(2,5)*
 SPB(5,1)))/(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*(-(SPA(1,2)*SPA(2,4)*SPB(4,1))-SPA(2,3)*SPA(2,5)*
 SPB(5,3)))/(complex<T>(2,0)*SPA(1,4)*SPA(2,3)*SPA(3,5)*SPA(4,5))+
 (complex<T>(0,1)*(-(SPA(1,2)*SPA(2,4)*SPB(4,1))+SPA(2,3)*SPA(2,5)*
 SPB(5,3)))/(complex<T>(2,0)*SPA(1,4)*SPA(2,5)*SPA(3,4)*SPA(3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpqmpmgap_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, p, m, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpqmpmgap SLC");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(2,4),3)*SPA(1,4))/(SPA(1,2)*SPA(1,5)*SPA(2,3)*
 SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPB(5,3),2)*
   SPA(2,4)*SPA(2,5))/(pow(SPB(4,3),2)*SPA(1,2)*SPA(3,4)*
   SPA(3,5)))+(pow(SPB(5,3),2)*(-(S(1,2)/S(3,4))+
   S(3,4)/S(1,2))*SPA(2,5)*SPA(4,5)*SPB(5,1))/
 (pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
  pow(SPB(4,3),3)*SPA(3,5))+(pow(SPB(3,1),2)*SPB(5,1))/
 (SPA(3,5)*SPB(2,1)*SPB(4,1)*SPB(4,3))+
(SPA(2,4)*SPB(3,1)*SPB(5,3))/(S(1,2)*SPA(3,5)*SPB(4,3))))/
complex<T>(3,0)+complex<T>(0,1)*((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/
   complex<T>(2,0)+(complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(4,5),2)*pow(SPB(5,1),3)*SPA(1,2)*SPA(1,4))/
(pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(2,5)/S(3,4))*SPA(3,5))-(pow(complex<T>(1,0)-S(1,4)/S(2,5),-1)*
 pow(SPB(3,1),2)*SPA(1,4)*SPA(2,3)*SPA(2,4))/
(complex<T>(2,0)*pow(SPA(2,5),3)*pow(SPB(5,2),2)*SPA(3,5))-
 (pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPB(5,3),2)*SPA(2,4)*
 SPA(2,5))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(1,2)*SPA(3,4)*
 SPA(3,5))+(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(2,5),-1))/
   complex<T>(2,0)+(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(2,5),-1))/complex<T>(2,0))*
 pow(SPB(3,1),3)*SPA(1,2)*SPA(1,4)*SPA(2,3)*SPA(3,4))/
(pow(SPA(2,5),4)*pow(SPB(5,2),3)*(complex<T>(1,0)-S(1,4)/S(2,5)-
  S(3,4)/S(2,5))*SPA(3,5))-(pow(SPA(2,4),2)*SPA(1,4)*
 SPB(2,1)*SPB(5,1))/(complex<T>(2,0)*pow(SPA(3,4),2)*S(2,5)*SPA(3,5)*
 SPB(4,3))-(pow(complex<T>(1,0)-S(3,4)/S(2,5),-1)*
 ((pow(SPA(1,4),2)*pow(SPB(5,1),2)*SPA(2,4)*SPA(2,5))/
   (SPA(1,2)*SPA(3,4)*SPA(3,5))+(complex<T>(1,0)*pow(SPA(2,4),2)*
SPA(1,4)*SPB(2,1)*SPB(5,1))/(complex<T>(2,0)*SPA(3,4)*SPA(3,5))))/
(pow(SPA(2,5),2)*pow(SPB(5,2),2))+
 (complex<T>(1,0)*pow(SPB(3,1),2)*SPA(2,3)*SPB(5,3))/
(complex<T>(2,0)*S(2,5)*SPA(3,5)*SPB(4,1)*SPB(4,3))+
 (SPA(1,4)*SPA(2,4)*SPB(5,1)*SPB(5,3))/(S(3,4)*SPA(1,2)*SPA(3,5)*
 SPB(5,2))+(complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*
 pow(SPB(5,3),2)*SPA(2,4)*SPA(2,5))/(pow(SPB(4,3),2)*
 SPA(1,2)*SPA(3,4)*SPA(3,5)))+(pow(SPB(5,3),2)*
(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPA(2,5)*SPA(4,5)*
SPB(5,1))/(pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
pow(SPB(4,3),3)*SPA(3,5))+(pow(SPB(3,1),2)*SPB(5,1))/
   (SPA(3,5)*SPB(2,1)*SPB(4,1)*SPB(4,3))+
  (SPA(2,4)*SPB(3,1)*SPB(5,3))/(S(1,2)*SPA(3,5)*SPB(4,3))))/
complex<T>(3,0))+complex<T>(0,-1)*((complex<T>(1,0)*pow(SPB(5,3),3)*(-(S(1,2)/S(4,5))+
  S(4,5)/S(1,2))*SPA(2,3)*SPA(2,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(1,2))+(complex<T>(1,0)*pow(SPB(5,3),3)*
 (-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPA(2,4)*SPA(2,5))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
 pow(SPB(4,3),3)*SPA(1,2))+(complex<T>(-2,0)*pow(SPB(5,3),4)*
 (-(S(1,2)/(complex<T>(6,0)*S(3,4)))+(complex<T>(1,0)*(S(1,2)/S(3,4)-
 S(3,4)/S(1,2)))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5))/(pow(SPA(1,2),5)*pow(SPB(2,1),4)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2)))+
 (complex<T>(1,0)*pow(SPB(5,3),4)*SPA(2,3)*SPA(2,5))/
(complex<T>(3,0)*pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPB(4,3)*SPB(5,4)))+
 complex<T>(0,1)*((complex<T>(1,0)*pow(SPB(5,3),3)*(-(S(1,2)/S(4,5))+
  S(4,5)/S(1,2))*SPA(2,3)*SPA(2,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(1,2))+(complex<T>(1,0)*pow(SPB(5,3),3)*
 (-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPA(2,4)*SPA(2,5))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
 pow(SPB(4,3),3)*SPA(1,2))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1))/complex<T>(2,0))*
 pow(SPA(1,4),2)*pow(SPB(3,1),3)*SPA(1,2)*SPA(3,4))/
(pow(SPA(4,5),4)*pow(SPB(5,4),3)*(complex<T>(1,0)-S(1,2)/S(4,5)-
  S(2,3)/S(4,5))*SPA(1,5))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPB(5,3),3)*SPA(2,5)*SPA(3,4)*SPA(4,5))/
(pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
  S(4,5)/S(1,2))*SPA(1,5))+(complex<T>(-2,0)*pow(SPB(5,3),4)*
 (-(S(1,2)/(complex<T>(6,0)*S(3,4)))+(complex<T>(1,0)*(S(1,2)/S(3,4)-
 S(3,4)/S(1,2)))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPA(2,3)*SPA(2,5)*
 SPA(3,4)*SPA(4,5))/(pow(SPA(1,2),5)*pow(SPB(2,1),4)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2)))+
 (complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*pow(SPA(1,4),2)*
 pow(SPB(3,1),2)*SPA(2,4))/(pow(SPA(2,3),2)*
 pow(SPB(3,2),2)*SPA(1,5)*SPA(4,5)))-
  (pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPB(5,3),2)*SPA(2,4))/
   (pow(SPB(5,4),2)*SPA(1,5)*SPA(4,5))-
  (pow(SPA(1,4),2)*pow(SPB(3,1),2)*SPA(2,4))/(S(2,3)*S(4,5)*
SPA(1,5)*SPA(4,5))-(pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*
pow(SPB(3,1),2)*SPA(1,4)*SPA(2,4)*SPA(3,4))/
   (pow(SPB(2,1),2)*SPA(1,2)*SPA(1,5)*SPA(2,3)*SPA(4,5))-
  (pow(SPB(5,3),2)*SPA(2,4)*SPA(2,5)*SPB(3,2))/
   (S(1,2)*S(4,5)*SPA(1,5)*SPB(4,3))-
  (pow(SPA(2,4),2)*SPA(1,4)*SPB(5,1))/(S(1,2)*SPA(1,5)*
SPA(2,3)*SPA(3,4))+(pow(SPB(3,1),2)*SPA(1,4)*SPB(5,3))/
   (S(4,5)*SPA(1,5)*SPB(2,1)*SPB(4,3))))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(SPB(5,3),4)*SPA(2,3)*SPA(2,5))/
(complex<T>(3,0)*pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPB(4,3)*SPB(5,4)))+
 complex<T>(0,1)*((((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1))/complex<T>(2,0))*
 pow(SPA(3,4),2)*pow(SPB(3,1),3)*SPA(1,2)*SPA(1,4))/
(pow(SPA(4,5),4)*pow(SPB(5,4),3)*(complex<T>(1,0)-S(1,2)/S(4,5)-
  S(2,3)/S(4,5))*SPA(3,5))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(2,3),-1)*pow(SPB(5,1),2)*
 SPA(1,4)*SPA(2,4)*SPA(2,5))/(complex<T>(2,0)*pow(SPA(2,3),3)*
 pow(SPB(3,2),2)*SPA(3,5))-(pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*
 pow(SPB(5,3),2)*SPA(2,3)*SPA(2,4))/(complex<T>(2,0)*pow(SPB(5,4),2)*
 SPA(1,2)*SPA(3,5)*SPA(4,5))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(2,3),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
 pow(SPB(5,1),3)*SPA(1,2)*SPA(1,4)*SPA(2,5)*SPA(4,5))/
(pow(SPA(2,3),4)*pow(SPB(3,2),3)*(complex<T>(1,0)-S(1,4)/S(2,3)-
  S(4,5)/S(2,3))*SPA(3,5))-(pow(complex<T>(1,0)-S(4,5)/S(2,3),-1)*
 ((pow(SPA(1,4),2)*pow(SPB(3,1),2)*SPA(2,3)*SPA(2,4))/
   (SPA(1,2)*SPA(3,5)*SPA(4,5))+(complex<T>(1,0)*pow(SPA(2,4),2)*
SPA(1,4)*SPB(2,1)*SPB(3,1))/(complex<T>(2,0)*SPA(3,5)*SPA(4,5))))/
(pow(SPA(2,3),2)*pow(SPB(3,2),2))+
 (SPA(1,4)*SPA(2,4)*SPB(3,1)*SPB(5,3))/(S(4,5)*SPA(1,2)*SPA(3,5)*
 SPB(3,2))+(complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*
 pow(SPB(5,3),2)*SPA(2,3)*SPA(2,4))/(pow(SPB(5,4),2)*
 SPA(1,2)*SPA(3,5)*SPA(4,5)))-(pow(SPB(5,3),2)*
(-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*SPA(2,3)*SPA(3,4)*
SPB(3,1))/(pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),2)*
pow(SPB(5,4),3)*SPA(3,5))+(pow(SPB(5,1),2)*SPB(3,1))/
   (SPA(3,5)*SPB(2,1)*SPB(4,1)*SPB(5,4))-
  (SPA(2,4)*SPB(5,1)*SPB(5,3))/(S(1,2)*SPA(3,5)*SPB(5,4))))/
complex<T>(3,0)-(pow(SPA(2,4),2)*SPA(1,4)*SPB(2,1)*SPB(3,1))/
(complex<T>(2,0)*pow(SPA(4,5),2)*S(2,3)*SPA(3,5)*SPB(5,4))-
 (pow(SPB(5,1),2)*SPA(2,5)*SPB(5,3))/(complex<T>(2,0)*S(2,3)*SPA(3,5)*
 SPB(4,1)*SPB(5,4)))+
 (complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*pow(SPB(5,3),2)*
   SPA(2,3)*SPA(2,4))/(pow(SPB(5,4),2)*SPA(1,2)*SPA(3,5)*
   SPA(4,5)))-(pow(SPB(5,3),2)*(-(S(1,2)/S(4,5))+
   S(4,5)/S(1,2))*SPA(2,3)*SPA(3,4)*SPB(3,1))/
 (pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),2)*
  pow(SPB(5,4),3)*SPA(3,5))+(pow(SPB(5,1),2)*SPB(3,1))/
 (SPA(3,5)*SPB(2,1)*SPB(4,1)*SPB(5,4))-
(SPA(2,4)*SPB(5,1)*SPB(5,3))/(S(1,2)*SPA(3,5)*SPB(5,4))))/
complex<T>(3,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpqmmpgap_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, p, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpqmmpgap SLC");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(2,3),3)*SPA(1,3))/(SPA(1,2)*SPA(1,4)*
 SPA(2,5)*SPA(3,4)*SPA(3,5))+(complex<T>(0,1)*pow(SPA(2,3),2)*SPA(1,3))/
(complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(2,3),2)*SPA(1,3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,4)*
 SPA(3,5)*SPA(4,5))+complex<T>(0,-1)*
((complex<T>(1,0)*pow(SPB(5,4),3)*(-(S(1,2)/S(3,5))+S(3,5)/S(1,2))*
 SPA(2,3)*SPA(2,4))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*
 pow(SPA(3,5),2)*pow(SPB(5,3),3)*SPA(1,2))-
 (pow(SPB(5,4),3)*(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPA(2,3)*
 SPA(2,5))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPA(1,2))+
 (complex<T>(2,0)*pow(SPB(5,4),4)*(-(S(1,2)/(complex<T>(6,0)*S(3,4)))+
  (complex<T>(1,0)*(S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(3,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(3,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),3)))*SPA(2,4)*SPA(2,5)*
 SPA(3,4)*SPA(3,5))/(pow(SPA(1,2),5)*pow(SPB(2,1),4)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(3,5)/S(1,2)))-
 (pow(SPB(5,4),4)*SPA(2,4)*SPA(2,5))/(complex<T>(3,0)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(4,3)*SPB(5,3)))+
 (complex<T>(0,-1)*(-((pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(5,4),3)*
   (S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
  (pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPA(1,2),4)*
   pow(SPB(2,1),3)*SPA(4,5)))+
(complex<T>(-3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*pow(SPB(5,4),2)*
  SPA(2,3)*SPA(2,5)*SPA(3,4))/(pow(SPA(1,2),3)*pow(SPB(2,1),2)*
  SPA(4,5))-(SPA(2,3)*SPB(4,1)*SPB(5,2))/(S(1,2)*SPA(4,5)*
  SPB(3,2))+(complex<T>(-2,0)*pow(SPA(2,3),2)*SPB(5,4))/
 (pow(SPA(1,2),2)*SPA(4,5)*SPB(2,1))-
(S(2,4)*SPB(4,1)*SPB(5,4))/(S(1,2)*SPA(4,5)*SPB(3,2)*
  SPB(4,3))))/complex<T>(3,0)+
 (complex<T>(0,-1)*(-((pow(SPA(2,4),2)*pow(SPA(3,5),2)*pow(SPB(5,4),3)*
   (S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
  (pow(complex<T>(1,0)-S(3,5)/S(1,2),3)*pow(SPA(1,2),4)*
   pow(SPB(2,1),3)*SPA(4,5)))+
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1)*pow(SPB(5,4),2)*
  SPA(2,3)*SPA(2,4)*SPA(3,5))/(pow(SPA(1,2),3)*pow(SPB(2,1),2)*
  SPA(4,5))+(SPA(2,3)*SPB(4,2)*SPB(5,1))/(S(1,2)*SPA(4,5)*
  SPB(3,2))+(complex<T>(-2,0)*pow(SPA(2,3),2)*SPB(5,4))/
 (pow(SPA(1,2),2)*SPA(4,5)*SPB(2,1))-
(S(2,5)*SPB(5,1)*SPB(5,4))/(S(1,2)*SPA(4,5)*SPB(3,2)*
  SPB(5,3))))/complex<T>(3,0)+
 complex<T>(0,1)*((complex<T>(1,0)*((pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPB(5,4),2)*
SPA(2,5)*SPA(3,5))/(pow(SPB(4,3),2)*SPA(1,5)*SPA(3,4)*
SPA(4,5))-(SPB(4,1)*SPB(5,4))/(SPA(4,5)*SPB(3,2)*
SPB(4,3))))/complex<T>(2,0)+
 (complex<T>(1,0)*(-((pow(SPA(2,5),2)*pow(SPA(3,4),2)*pow(SPB(5,4),3)*
 (S(1,2)/S(3,4)-S(3,4)/S(1,2)))/
(pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3)*SPA(4,5)))+
  (complex<T>(-3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*pow(SPB(5,4),2)*
SPA(2,3)*SPA(2,5)*SPA(3,4))/(pow(SPA(1,2),3)*
pow(SPB(2,1),2)*SPA(4,5))-(SPA(2,3)*SPB(4,1)*SPB(5,2))/
   (S(1,2)*SPA(4,5)*SPB(3,2))+(complex<T>(-2,0)*pow(SPA(2,3),2)*
SPB(5,4))/(pow(SPA(1,2),2)*SPA(4,5)*SPB(2,1))-
  (S(2,4)*SPB(4,1)*SPB(5,4))/(S(1,2)*SPA(4,5)*SPB(3,2)*
SPB(4,3))))/complex<T>(3,0))+
 complex<T>(0,1)*((complex<T>(1,0)*pow(SPB(5,4),3)*(-(S(1,2)/S(3,5))+
  S(3,5)/S(1,2))*SPA(2,3)*SPA(2,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*pow(SPA(3,5),2)*
 pow(SPB(5,3),3)*SPA(1,2))-
 (pow(SPB(5,4),3)*(-(S(1,2)/S(3,4))+S(3,4)/S(1,2))*SPA(2,3)*
 SPA(2,5))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPA(1,2))+
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(3,4),-1))/complex<T>(2,0))*
 pow(SPA(1,3),2)*pow(SPB(5,1),3)*SPA(1,2)*SPA(3,5))/
(pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
  S(2,5)/S(3,4))*SPA(1,4))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPB(5,4),3)*SPA(2,4)*SPA(3,4)*SPA(3,5))/
(pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
  S(3,5)/S(1,2))*SPA(1,4))+(complex<T>(2,0)*pow(SPB(5,4),4)*
 (-(S(1,2)/(complex<T>(6,0)*S(3,4)))+(complex<T>(1,0)*(S(1,2)/S(3,4)-
 S(3,4)/S(1,2)))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(3,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,4)/S(1,2)-S(3,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),3)))*SPA(2,4)*SPA(2,5)*
 SPA(3,4)*SPA(3,5))/(pow(SPA(1,2),5)*pow(SPB(2,1),4)*
 (complex<T>(1,0)-S(3,4)/S(1,2)-S(3,5)/S(1,2)))-
 (pow(SPB(5,4),4)*SPA(2,4)*SPA(2,5))/(complex<T>(3,0)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(4,3)*SPB(5,3))+
 (complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(3,4)/S(2,5),-1)*pow(SPA(1,3),2)*
 pow(SPB(5,1),2)*SPA(2,3))/(pow(SPA(2,5),2)*
 pow(SPB(5,2),2)*SPA(1,4)*SPA(3,4)))-
  (pow(complex<T>(1,0)-S(1,2)/S(3,4),-1)*pow(SPB(5,4),2)*SPA(2,3))/
   (pow(SPB(4,3),2)*SPA(1,4)*SPA(3,4))-
  (pow(SPA(1,3),2)*pow(SPB(5,1),2)*SPA(2,3))/(S(2,5)*S(3,4)*
SPA(1,4)*SPA(3,4))+(pow(complex<T>(1,0)-S(3,4)/S(1,2),-1)*
pow(SPB(5,1),2)*SPA(1,3)*SPA(2,3)*SPA(3,5))/
   (pow(SPB(2,1),2)*SPA(1,2)*SPA(1,4)*SPA(2,5)*SPA(3,4))+
  (pow(SPA(2,3),2)*SPA(1,3)*SPB(4,1))/(S(1,2)*SPA(1,4)*
SPA(2,5)*SPA(3,5))+(pow(SPB(5,4),2)*SPA(2,3)*SPA(2,4)*
SPB(5,2))/(S(1,2)*S(3,4)*SPA(1,4)*SPB(5,3))+
  (pow(SPB(5,1),2)*SPA(1,3)*SPB(5,4))/(S(3,4)*SPA(1,4)*
SPB(2,1)*SPB(5,3))))/complex<T>(2,0))+
 complex<T>(0,1)*((complex<T>(1,0)*(-((pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
 pow(SPB(5,4),2)*SPA(2,4)*SPA(3,4))/(pow(SPB(5,3),2)*
 SPA(1,4)*SPA(3,5)*SPA(4,5)))-(SPB(5,1)*SPB(5,4))/
   (SPA(4,5)*SPB(3,2)*SPB(5,3))))/complex<T>(2,0)+
 (complex<T>(1,0)*(-((pow(SPA(2,4),2)*pow(SPA(3,5),2)*pow(SPB(5,4),3)*
 (S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
(pow(complex<T>(1,0)-S(3,5)/S(1,2),3)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3)*SPA(4,5)))+
  (complex<T>(3,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1)*pow(SPB(5,4),2)*
SPA(2,3)*SPA(2,4)*SPA(3,5))/(pow(SPA(1,2),3)*
pow(SPB(2,1),2)*SPA(4,5))+(SPA(2,3)*SPB(4,2)*SPB(5,1))/
   (S(1,2)*SPA(4,5)*SPB(3,2))+(complex<T>(-2,0)*pow(SPA(2,3),2)*
SPB(5,4))/(pow(SPA(1,2),2)*SPA(4,5)*SPB(2,1))-
  (S(2,5)*SPB(5,1)*SPB(5,4))/(S(1,2)*SPA(4,5)*SPB(3,2)*
SPB(5,3))))/complex<T>(3,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpqmmmgap_SLC
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, m, m, gap}, SLC}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpqmmmgap SLC");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,1),2)*SPB(5,2))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*
 SPB(4,3)*SPB(5,4))+(complex<T>(0,1)*pow(SPB(5,1),3)*SPB(5,2))/
(SPB(2,1)*SPB(3,2)*SPB(4,1)*SPB(5,3)*SPB(5,4))+
 (complex<T>(0,-1)*(-((pow(SPA(2,3),2)*SPA(2,4))/(SPA(1,2)*SPA(2,5)*
   SPA(3,5)*SPB(4,3)))+(SPA(2,3)*SPA(3,4)*SPB(5,1))/
 (S(1,2)*SPA(3,5)*SPB(4,3))+(pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
  pow(SPA(3,4),2)*SPB(4,1)*SPB(5,1))/(pow(SPA(3,5),2)*SPB(2,1)*
  SPB(4,3)*SPB(5,3))-(pow(SPA(3,4),2)*(-(S(1,2)/S(3,5))+
   S(3,5)/S(1,2))*SPA(2,4)*SPB(4,1)*SPB(5,4))/
 (pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*pow(SPA(3,5),3)*
  pow(SPB(5,3),2)*SPB(4,3))))/complex<T>(3,0)+
 (complex<T>(0,-1)*((pow(SPA(3,4),3)*pow(SPB(3,1),2)*pow(SPB(5,4),2)*
  (S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
 (pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),3)*
  pow(SPB(2,1),4)*SPB(4,3))+(complex<T>(2,0)*pow(SPB(5,1),2)*
  SPA(3,4))/(pow(SPB(2,1),2)*SPA(1,2)*SPB(4,3))+
(S(1,4)*SPA(2,4)*SPA(3,4))/(S(1,2)*SPA(1,5)*SPA(4,5)*
  SPB(4,3))+(SPA(1,3)*SPA(2,4)*SPB(5,1))/(S(1,2)*SPA(1,5)*
  SPB(4,3))+(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*
  pow(SPA(3,4),2)*SPB(3,1)*SPB(5,1)*SPB(5,4))/
 (pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(4,3))))/complex<T>(3,0)+
 complex<T>(0,-1)*((complex<T>(1,0)*pow(SPA(3,4),4)*SPB(3,1)*SPB(4,1))/
(complex<T>(3,0)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPA(3,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(3,4),3)*(-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*
 SPB(3,1)*SPB(5,1))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*
 pow(SPA(4,5),3)*pow(SPB(5,4),2)*SPB(2,1))-
 (pow(SPA(3,4),3)*(-(S(1,2)/S(3,5))+S(3,5)/S(1,2))*SPB(4,1)*
 SPB(5,1))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*
 pow(SPA(3,5),3)*pow(SPB(5,3),2)*SPB(2,1))+
 (complex<T>(-2,0)*pow(SPA(3,4),4)*(-(S(1,2)/(complex<T>(6,0)*S(3,5)))+
  (complex<T>(1,0)*(S(1,2)/S(3,5)-S(3,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,5)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPB(3,1)*SPB(4,1)*
 SPB(5,3)*SPB(5,4))/(pow(SPA(1,2),4)*pow(SPB(2,1),5)*
 (complex<T>(1,0)-S(3,5)/S(1,2)-S(4,5)/S(1,2))))+
 complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(2,3),2)*SPA(3,4)*SPB(3,1))/
(complex<T>(2,0)*S(1,4)*SPA(2,5)*SPA(3,5)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPB(5,1),2)*SPA(1,2)*SPA(2,4)*SPB(5,2))/
(complex<T>(2,0)*pow(SPB(5,3),2)*S(1,4)*SPA(3,5)*SPB(4,3))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(3,5),-1))/complex<T>(2,0))*
 pow(SPA(2,4),3)*pow(SPB(5,4),2)*SPB(2,1)*SPB(5,2))/
(pow(SPA(3,5),3)*pow(SPB(5,3),4)*(complex<T>(1,0)-S(1,2)/S(3,5)-
  S(1,4)/S(3,5))*SPB(4,3))-(SPA(2,4)*SPA(3,4)*SPB(5,1)*
 SPB(5,2))/(S(3,5)*SPA(1,4)*SPB(2,1)*SPB(4,3))-
 (pow(complex<T>(1,0)-S(2,5)/S(1,4),-1)*pow(SPA(2,3),2)*SPB(3,1)*
 SPB(5,1)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(1,4),2)*pow(SPB(4,1),3)*
 SPB(4,3))-(pow(complex<T>(1,0)-S(3,5)/S(1,4),-1)*
 (-((pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(4,1)*SPB(5,1))/
(SPB(2,1)*SPB(4,3)*SPB(5,3)))-(pow(SPB(5,1),2)*SPA(1,2)*
SPA(2,4)*SPB(5,2))/(complex<T>(2,0)*SPB(4,3)*SPB(5,3))))/
(pow(SPA(1,4),2)*pow(SPB(4,1),2))+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*pow(SPA(3,4),2)*
 SPB(4,1)*SPB(5,1))/(complex<T>(2,0)*pow(SPA(3,5),2)*SPB(2,1)*SPB(4,3)*
 SPB(5,3))-(((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,5)/S(1,4),-1))/
   complex<T>(2,0)+(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,4),-1))/complex<T>(2,0))*
 pow(SPA(2,3),3)*SPB(2,1)*SPB(3,1)*SPB(5,2)*SPB(5,3))/
(pow(SPA(1,4),3)*pow(SPB(4,1),4)*(complex<T>(1,0)-S(2,5)/S(1,4)-
  S(3,5)/S(1,4))*SPB(4,3))+
 (complex<T>(1,0)*(-((pow(SPA(2,3),2)*SPA(2,4))/(SPA(1,2)*SPA(2,5)*
 SPA(3,5)*SPB(4,3)))+(SPA(2,3)*SPA(3,4)*SPB(5,1))/
   (S(1,2)*SPA(3,5)*SPB(4,3))+(pow(complex<T>(1,0)-S(1,2)/S(3,5),
-1)*pow(SPA(3,4),2)*SPB(4,1)*SPB(5,1))/(pow(SPA(3,5),2)*
SPB(2,1)*SPB(4,3)*SPB(5,3))-(pow(SPA(3,4),2)*
(-(S(1,2)/S(3,5))+S(3,5)/S(1,2))*SPA(2,4)*SPB(4,1)*
SPB(5,4))/(pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*pow(SPA(3,5),3)*
pow(SPB(5,3),2)*SPB(4,3))))/complex<T>(3,0))+
 complex<T>(0,1)*((complex<T>(1,0)*((SPA(2,4)*SPA(3,4))/(SPA(1,5)*SPA(4,5)*
SPB(4,3))-(pow(complex<T>(1,0)-S(1,2)/S(4,5),-1)*
pow(SPA(3,4),2)*SPB(3,1)*SPB(5,3))/(pow(SPA(4,5),2)*
SPB(3,2)*SPB(4,3)*SPB(5,4))))/complex<T>(2,0)+
 (complex<T>(1,0)*((pow(SPA(3,4),3)*pow(SPB(3,1),2)*pow(SPB(5,4),2)*
(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),3)*
pow(SPB(2,1),4)*SPB(4,3))+(complex<T>(2,0)*pow(SPB(5,1),2)*
SPA(3,4))/(pow(SPB(2,1),2)*SPA(1,2)*SPB(4,3))+
  (S(1,4)*SPA(2,4)*SPA(3,4))/(S(1,2)*SPA(1,5)*SPA(4,5)*
SPB(4,3))+(SPA(1,3)*SPA(2,4)*SPB(5,1))/(S(1,2)*SPA(1,5)*
SPB(4,3))+(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1)*
pow(SPA(3,4),2)*SPB(3,1)*SPB(5,1)*SPB(5,4))/
   (pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(4,3))))/complex<T>(3,0))+
 complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(3,4),4)*SPB(3,1)*SPB(4,1))/
(complex<T>(3,0)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPA(3,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(3,4),3)*(-(S(1,2)/S(4,5))+S(4,5)/S(1,2))*
 SPB(3,1)*SPB(5,1))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*
 pow(SPA(4,5),3)*pow(SPB(5,4),2)*SPB(2,1))-
 (pow(SPA(3,4),3)*(-(S(1,2)/S(3,5))+S(3,5)/S(1,2))*SPB(4,1)*
 SPB(5,1))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),3)*
 pow(SPA(3,5),3)*pow(SPB(5,3),2)*SPB(2,1))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,2)/S(3,5),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(1,4)/S(3,5),-1))/complex<T>(2,0))*
 pow(SPA(2,4),3)*pow(SPB(5,2),2)*SPB(2,1)*SPB(5,4))/
(pow(SPA(3,5),3)*pow(SPB(5,3),4)*(complex<T>(1,0)-S(1,2)/S(3,5)-
  S(1,4)/S(3,5))*SPB(3,2))-
 (((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0)+
  (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))*
 pow(SPA(3,4),3)*SPB(3,1)*SPB(5,3)*SPB(5,4))/
(pow(SPA(1,2),3)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,5)/S(1,2)-
  S(4,5)/S(1,2))*SPB(3,2))+(complex<T>(-2,0)*pow(SPA(3,4),4)*
 (-(S(1,2)/(complex<T>(6,0)*S(3,5)))+(complex<T>(1,0)*(S(1,2)/S(3,5)-
 S(3,5)/S(1,2)))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),3))-
  S(1,2)/(complex<T>(6,0)*S(4,5))+
  ((complex<T>(1,0)*pow(complex<T>(1,0)-S(3,5)/S(1,2),-1))/complex<T>(2,0)+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),-1))/complex<T>(2,0))/
   (complex<T>(1,0)-S(3,5)/S(1,2)-S(4,5)/S(1,2))+
  (complex<T>(1,0)*(S(1,2)/S(4,5)-S(4,5)/S(1,2)))/
   (complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)))*SPB(3,1)*SPB(4,1)*
 SPB(5,3)*SPB(5,4))/(pow(SPA(1,2),4)*pow(SPB(2,1),5)*
 (complex<T>(1,0)-S(3,5)/S(1,2)-S(4,5)/S(1,2)))+
 (complex<T>(1,0)*(-((pow(SPA(3,4),2)*SPA(1,4)*SPB(3,1)*SPB(5,1))/
(S(1,2)*S(3,5)*SPA(4,5)*SPB(3,2)))+
  (pow(SPA(2,4),2)*SPA(3,4)*SPB(5,2))/(S(3,5)*SPA(1,2)*
SPA(4,5)*SPB(3,2))+(pow(complex<T>(1,0)-S(1,2)/S(3,5),-1)*
pow(SPA(3,4),2)*SPB(5,1))/(pow(SPA(3,5),2)*SPB(3,2)*
SPB(5,3))+(pow(complex<T>(1,0)-S(3,5)/S(1,4),-1)*
pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(5,1))/
   (pow(SPA(1,4),2)*pow(SPB(4,1),2)*SPB(3,2)*SPB(5,3))+
  (pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPB(5,1))/(S(1,4)*S(3,5)*
SPB(3,2)*SPB(5,3))-(pow(SPB(5,1),2)*SPA(2,3)*SPB(5,2))/
   (S(1,2)*SPB(3,2)*SPB(4,1)*SPB(5,4))-
  (pow(complex<T>(1,0)-S(3,5)/S(1,2),-1)*pow(SPA(2,4),2)*SPB(5,1)*
SPB(5,2)*SPB(5,4))/(pow(SPA(1,2),2)*SPB(2,1)*SPB(3,2)*
SPB(4,1)*SPB(5,3))))/complex<T>(2,0))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpppqmgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, p, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpppqmgap nf");
#endif
 
 return( (complex<T>(0,-1)*SPA(2,4)*SPA(3,4)*SPB(3,2))/(complex<T>(3,0)*pow(SPA(2,3),2)*
SPA(1,5)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qppmqmgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, p, m, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qppmqmgap nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1y_qpmpqmgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, p, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpmpqmgap nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1y_qpmmqmgap_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, m, m, qm, gap}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpmmqmgap nf");
#endif
 
 return( (complex<T>(0,-1)*SPA(2,4)*SPA(3,4))/(complex<T>(3,0)*SPA(1,5)*SPA(4,5)*
SPB(3,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpgapqmpp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, gap, qm, p, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpgapqmpp L");
#endif
 
 return( (complex<T>(0,-1)*SPA(3,4)*SPA(3,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(4,5),2)*
 SPA(1,2)*SPA(2,3))+
 (complex<T>(0,-1)*((pow(SPA(1,3),2)*SPB(2,1))/(SPA(1,2)*SPA(1,5)*SPA(3,4)*
  SPA(4,5))+(SPA(1,3)*SPA(3,5)*SPB(5,4))/(SPA(1,2)*SPA(1,5)*
  SPA(2,3)*SPA(4,5))))/complex<T>(2,0)
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpgapqmpm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, gap, qm, p, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpgapqmpm L");
#endif
 
 return( complex<T>(0,-1)*(-((pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*pow(SPB(4,1),2)*
 SPA(1,3)*SPA(3,5))/(complex<T>(2,0)*pow(SPB(5,1),2)*SPA(1,2)*SPA(1,5)*
 SPA(2,3)))+(complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1)*
pow(SPB(4,2),2)*SPA(2,5)*SPA(3,5))/(complex<T>(2,0)*pow(SPA(1,5),2)*
pow(SPB(5,1),2)*SPA(1,2))-
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),-1))/complex<T>(2,0))*
pow(SPA(1,3),2)*pow(SPB(4,1),3)*SPA(1,5)*SPA(4,5))/
 (pow(SPA(2,3),4)*pow(SPB(3,2),3)*(complex<T>(1,0)-S(1,5)/S(2,3)-
 S(4,5)/S(2,3))*SPA(1,2))-
(((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1))/complex<T>(2,0)+
 (complex<T>(1,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),-1))/complex<T>(2,0))*
pow(SPB(4,2),3)*SPA(2,3)*SPA(2,5)*SPA(4,5))/
 (pow(SPA(1,5),3)*pow(SPB(5,1),3)*(complex<T>(1,0)-S(2,3)/S(1,5)-
 S(3,4)/S(1,5))*SPA(1,2))-(SPA(2,5)*SPA(3,5)*SPB(4,2))/
 (complex<T>(2,0)*S(1,5)*SPA(1,2)*SPA(2,4))-
(pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*SPA(2,3)*SPA(3,5)*SPA(4,5)*
SPB(4,2)*SPB(4,3))/(complex<T>(2,0)*pow(SPA(1,5),2)*pow(SPB(5,1),2)*
SPA(1,2)*SPA(2,4))+(complex<T>(1,0)*pow(SPA(1,3),2)*pow(SPB(4,1),3))/
 (complex<T>(2,0)*S(2,3)*SPA(1,2)*SPA(2,3)*SPB(5,1)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpgapqmmp_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, gap, qm, m, p}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpgapqmmp L");
#endif
 
 return( complex<T>(0,-1)*((complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),-1)*
pow(SPA(1,3),2)*pow(SPB(5,1),2)*SPA(1,4))/
 (complex<T>(2,0)*pow(SPB(5,4),2)*SPA(1,2)*SPA(1,5)*SPA(2,3)*SPA(4,5))+
(complex<T>(1,0)*pow(SPB(5,2),2))/(complex<T>(2,0)*SPA(1,5)*SPB(4,3)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpgapqmmm_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, gap, qm, m, m}, L}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpgapqmmm L");
#endif
 
 return( complex<T>(0,-1)*((complex<T>(1,0)*SPA(3,4)*SPA(3,5))/(complex<T>(3,0)*SPA(1,2)*SPA(2,3)*
SPB(5,4))+(complex<T>(1,0)*SPA(1,3)*SPA(3,5)*SPA(4,5))/
 (complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(2,3)*SPB(5,4))+
(complex<T>(1,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),-1)*pow(SPA(4,5),2)*
pow(SPB(4,2),2))/(complex<T>(2,0)*pow(SPA(1,5),2)*SPB(4,3)*SPB(5,1)*
SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpgapqmpp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, gap, qm, p, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpgapqmpp nf");
#endif
 
 return( (complex<T>(0,1)*SPA(3,4)*SPA(3,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(4,5),2)*
SPA(1,2)*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1y_qpgapqmpm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, gap, qm, p, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpgapqmpm nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1y_qpgapqmmp_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, gap, qm, m, p}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpgapqmmp nf");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1y_qpgapqmmm_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, gap, qm, m, m}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2g1y :  qpgapqmmm nf");
#endif
 
 return( (complex<T>(0,1)*SPA(3,4)*SPA(3,5))/(complex<T>(3,0)*SPA(1,2)*SPA(2,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 // *************** table of switch values ************* 
 
#define _R_qmmmqpgam_L R2q2g1y_5617_L
#define _R_qmmpqpgam_L R2q2g1y_5725_L
#define _R_qmpmqpgam_L R2q2g1y_5635_L
#define _R_qmppqpgam_L R2q2g1y_5743_L
#define _R_qmmqpmgam_SLC R2q2g1y_5257_SLC
#define _R_qmmqppgam_SLC R2q2g1y_5905_SLC
#define _R_qmpqpmgam_SLC R2q2g1y_5275_SLC
#define _R_qmpqppgam_SLC R2q2g1y_5923_SLC
#define _R_qmqpmmgam_SLC R2q2g1y_5197_SLC
#define _R_qmqpmpgam_SLC R2q2g1y_5845_SLC
#define _R_qmqppmgam_SLC R2q2g1y_5305_SLC
#define _R_qmqpppgam_SLC R2q2g1y_5953_SLC
#define _R_qmmmqpgam_nf R2q2g1y_5617_nf
#define _R_qmmpqpgam_nf R2q2g1y_5725_nf
#define _R_qmpmqpgam_nf R2q2g1y_5635_nf
#define _R_qmppqpgam_nf R2q2g1y_5743_nf
#define _R_qmgamqpmm_L R2q2g1y_97_L
#define _R_qmgamqpmp_L R2q2g1y_3985_L
#define _R_qmgamqppm_L R2q2g1y_745_L
#define _R_qmgamqppp_L R2q2g1y_4633_L
#define _R_qmgamqpmm_nf R2q2g1y_97_nf
#define _R_qmgamqpmp_nf R2q2g1y_3985_nf
#define _R_qmgamqppm_nf R2q2g1y_745_nf
#define _R_qmgamqppp_nf R2q2g1y_4633_nf
#define _R_qpppqmgap_L R2q2g1y_6824_L
#define _R_qppmqmgap_L R2q2g1y_6716_L
#define _R_qpmpqmgap_L R2q2g1y_6806_L
#define _R_qpmmqmgap_L R2q2g1y_6698_L
#define _R_qppqmpgap_SLC R2q2g1y_7184_SLC
#define _R_qppqmmgap_SLC R2q2g1y_6536_SLC
#define _R_qpmqmpgap_SLC R2q2g1y_7166_SLC
#define _R_qpmqmmgap_SLC R2q2g1y_6518_SLC
#define _R_qpqmppgap_SLC R2q2g1y_7244_SLC
#define _R_qpqmpmgap_SLC R2q2g1y_6596_SLC
#define _R_qpqmmpgap_SLC R2q2g1y_7136_SLC
#define _R_qpqmmmgap_SLC R2q2g1y_6488_SLC
#define _R_qpppqmgap_nf R2q2g1y_6824_nf
#define _R_qppmqmgap_nf R2q2g1y_6716_nf
#define _R_qpmpqmgap_nf R2q2g1y_6806_nf
#define _R_qpmmqmgap_nf R2q2g1y_6698_nf
#define _R_qpgapqmpp_L R2q2g1y_4604_L
#define _R_qpgapqmpm_L R2q2g1y_716_L
#define _R_qpgapqmmp_L R2q2g1y_3956_L
#define _R_qpgapqmmm_L R2q2g1y_68_L
#define _R_qpgapqmpp_nf R2q2g1y_4604_nf
#define _R_qpgapqmpm_nf R2q2g1y_716_nf
#define _R_qpgapqmmp_nf R2q2g1y_3956_nf
#define _R_qpgapqmmm_nf R2q2g1y_68_nf
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmmmqpgam_L case 5617 : \
          return &R2q2g1y_5617_L
#define _CASE_qmmpqpgam_L case 5725 : \
          return &R2q2g1y_5725_L
#define _CASE_qmpmqpgam_L case 5635 : \
          return &R2q2g1y_5635_L
#define _CASE_qmppqpgam_L case 5743 : \
          return &R2q2g1y_5743_L
#define _CASE_qmmqpmgam_SLC case 5257 : \
          return &R2q2g1y_5257_SLC
#define _CASE_qmmqppgam_SLC case 5905 : \
          return &R2q2g1y_5905_SLC
#define _CASE_qmpqpmgam_SLC case 5275 : \
          return &R2q2g1y_5275_SLC
#define _CASE_qmpqppgam_SLC case 5923 : \
          return &R2q2g1y_5923_SLC
#define _CASE_qmqpmmgam_SLC case 5197 : \
          return &R2q2g1y_5197_SLC
#define _CASE_qmqpmpgam_SLC case 5845 : \
          return &R2q2g1y_5845_SLC
#define _CASE_qmqppmgam_SLC case 5305 : \
          return &R2q2g1y_5305_SLC
#define _CASE_qmqpppgam_SLC case 5953 : \
          return &R2q2g1y_5953_SLC
#define _CASE_qmmmqpgam_nf case 5617 : \
          return &R2q2g1y_5617_nf
#define _CASE_qmmpqpgam_nf case 5725 : \
          return &R2q2g1y_5725_nf
#define _CASE_qmpmqpgam_nf case 5635 : \
          return &R2q2g1y_5635_nf
#define _CASE_qmppqpgam_nf case 5743 : \
          return &R2q2g1y_5743_nf
#define _CASE_qmgamqpmm_L case 97 : \
          return &R2q2g1y_97_L
#define _CASE_qmgamqpmp_L case 3985 : \
          return &R2q2g1y_3985_L
#define _CASE_qmgamqppm_L case 745 : \
          return &R2q2g1y_745_L
#define _CASE_qmgamqppp_L case 4633 : \
          return &R2q2g1y_4633_L
#define _CASE_qmgamqpmm_nf case 97 : \
          return &R2q2g1y_97_nf
#define _CASE_qmgamqpmp_nf case 3985 : \
          return &R2q2g1y_3985_nf
#define _CASE_qmgamqppm_nf case 745 : \
          return &R2q2g1y_745_nf
#define _CASE_qmgamqppp_nf case 4633 : \
          return &R2q2g1y_4633_nf
#define _CASE_qpppqmgap_L case 6824 : \
          return &R2q2g1y_6824_L
#define _CASE_qppmqmgap_L case 6716 : \
          return &R2q2g1y_6716_L
#define _CASE_qpmpqmgap_L case 6806 : \
          return &R2q2g1y_6806_L
#define _CASE_qpmmqmgap_L case 6698 : \
          return &R2q2g1y_6698_L
#define _CASE_qppqmpgap_SLC case 7184 : \
          return &R2q2g1y_7184_SLC
#define _CASE_qppqmmgap_SLC case 6536 : \
          return &R2q2g1y_6536_SLC
#define _CASE_qpmqmpgap_SLC case 7166 : \
          return &R2q2g1y_7166_SLC
#define _CASE_qpmqmmgap_SLC case 6518 : \
          return &R2q2g1y_6518_SLC
#define _CASE_qpqmppgap_SLC case 7244 : \
          return &R2q2g1y_7244_SLC
#define _CASE_qpqmpmgap_SLC case 6596 : \
          return &R2q2g1y_6596_SLC
#define _CASE_qpqmmpgap_SLC case 7136 : \
          return &R2q2g1y_7136_SLC
#define _CASE_qpqmmmgap_SLC case 6488 : \
          return &R2q2g1y_6488_SLC
#define _CASE_qpppqmgap_nf case 6824 : \
          return &R2q2g1y_6824_nf
#define _CASE_qppmqmgap_nf case 6716 : \
          return &R2q2g1y_6716_nf
#define _CASE_qpmpqmgap_nf case 6806 : \
          return &R2q2g1y_6806_nf
#define _CASE_qpmmqmgap_nf case 6698 : \
          return &R2q2g1y_6698_nf
#define _CASE_qpgapqmpp_L case 4604 : \
          return &R2q2g1y_4604_L
#define _CASE_qpgapqmpm_L case 716 : \
          return &R2q2g1y_716_L
#define _CASE_qpgapqmmp_L case 3956 : \
          return &R2q2g1y_3956_L
#define _CASE_qpgapqmmm_L case 68 : \
          return &R2q2g1y_68_L
#define _CASE_qpgapqmpp_nf case 4604 : \
          return &R2q2g1y_4604_nf
#define _CASE_qpgapqmpm_nf case 716 : \
          return &R2q2g1y_716_nf
#define _CASE_qpgapqmmp_nf case 3956 : \
          return &R2q2g1y_3956_nf
#define _CASE_qpgapqmmm_nf case 68 : \
          return &R2q2g1y_68_nf
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_qmmmqpgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmmmqpgam_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmpqpgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmmpqpgam_L(ep,mpc);}
 
template <class T> complex<T> _R_qmpmqpgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmpmqpgam_L(ep,mpc);}
 
template <class T> complex<T> _R_qmppqpgam_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmppqpgam_L(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpmgam_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmmqpmgam_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qmmqppgam_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmmqppgam_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpmgam_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmpqpmgam_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qmpqppgam_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmpqppgam_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmmgam_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmqpmmgam_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmpgam_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmqpmpgam_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qmqppmgam_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmqppmgam_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qmqpppgam_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmqpppgam_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qmmmqpgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmmmqpgam_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmmpqpgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmmpqpgam_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmpmqpgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmpmqpgam_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmppqpgam_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmppqpgam_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmgamqpmm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmgamqpmm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmgamqpmp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmgamqpmp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmgamqppm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmgamqppm_L(ep,mpc);}
 
template <class T> complex<T> _R_qmgamqppp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmgamqppp_L(ep,mpc);}
 
template <class T> complex<T> _R_qmgamqpmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmgamqpmm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmgamqpmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmgamqpmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmgamqppm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmgamqppm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qmgamqppp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qmgamqppp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpppqmgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpppqmgap_L(ep,mpc);}
 
template <class T> complex<T> _R_qppmqmgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qppmqmgap_L(ep,mpc);}
 
template <class T> complex<T> _R_qpmpqmgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpmpqmgap_L(ep,mpc);}
 
template <class T> complex<T> _R_qpmmqmgap_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpmmqmgap_L(ep,mpc);}
 
template <class T> complex<T> _R_qppqmpgap_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qppqmpgap_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qppqmmgap_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qppqmmgap_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpmqmpgap_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpmqmpgap_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpmqmmgap_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpmqmmgap_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpqmppgap_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpqmppgap_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpqmpmgap_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpqmpmgap_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmpgap_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpqmmpgap_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpqmmmgap_SLC(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpqmmmgap_SLC(ep,mpc);}
 
template <class T> complex<T> _R_qpppqmgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpppqmgap_nf(ep,mpc);}
 
template <class T> complex<T> _R_qppmqmgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qppmqmgap_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpmpqmgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpmpqmgap_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpmmqmgap_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpmmqmgap_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpgapqmpp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpgapqmpp_L(ep,mpc);}
 
template <class T> complex<T> _R_qpgapqmpm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpgapqmpm_L(ep,mpc);}
 
template <class T> complex<T> _R_qpgapqmmp_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpgapqmmp_L(ep,mpc);}
 
template <class T> complex<T> _R_qpgapqmmm_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpgapqmmm_L(ep,mpc);}
 
template <class T> complex<T> _R_qpgapqmpp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpgapqmpp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpgapqmpm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpgapqmpm_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpgapqmmp_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpgapqmmp_nf(ep,mpc);}
 
template <class T> complex<T> _R_qpgapqmmm_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1y_qpgapqmmm_nf(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q2g1y_L_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmmmqpgam_L;
       _CASE_qmmpqpgam_L;
       _CASE_qmpmqpgam_L;
       _CASE_qmppqpgam_L;
       _CASE_qmgamqpmm_L;
       _CASE_qmgamqpmp_L;
       _CASE_qmgamqppm_L;
       _CASE_qmgamqppp_L;
       _CASE_qpppqmgap_L;
       _CASE_qppmqmgap_L;
       _CASE_qpmpqmgap_L;
       _CASE_qpmmqmgap_L;
       _CASE_qpgapqmpp_L;
       _CASE_qpgapqmpm_L;
       _CASE_qpgapqmmp_L;
       _CASE_qpgapqmmm_L;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g1y_nf_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmmmqpgam_nf;
       _CASE_qmmpqpgam_nf;
       _CASE_qmpmqpgam_nf;
       _CASE_qmppqpgam_nf;
       _CASE_qmgamqpmm_nf;
       _CASE_qmgamqpmp_nf;
       _CASE_qmgamqppm_nf;
       _CASE_qmgamqppp_nf;
       _CASE_qpppqmgap_nf;
       _CASE_qppmqmgap_nf;
       _CASE_qpmpqmgap_nf;
       _CASE_qpmmqmgap_nf;
       _CASE_qpgapqmpp_nf;
       _CASE_qpgapqmpm_nf;
       _CASE_qpgapqmmp_nf;
       _CASE_qpgapqmmm_nf;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g1y_SLC_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmmqpmgam_SLC;
       _CASE_qmmqppgam_SLC;
       _CASE_qmpqpmgam_SLC;
       _CASE_qmpqppgam_SLC;
       _CASE_qmqpmmgam_SLC;
       _CASE_qmqpmpgam_SLC;
       _CASE_qmqppmgam_SLC;
       _CASE_qmqpppgam_SLC;
       _CASE_qppqmpgap_SLC;
       _CASE_qppqmmgap_SLC;
       _CASE_qpmqmpgap_SLC;
       _CASE_qpmqmmgap_SLC;
       _CASE_qpqmppgap_SLC;
       _CASE_qpqmpmgap_SLC;
       _CASE_qpqmmpgap_SLC;
       _CASE_qpqmmmgap_SLC;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q2g1y_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g1y_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g1y_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g1y_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2g1y_nf_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g1y_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g1y_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g1y_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2g1y_SLC_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2g1y_SLC_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2g1y_SLC_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2g1y_SLC_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif




}
