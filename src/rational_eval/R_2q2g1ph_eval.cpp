/*
* R_2q2g1ph_eval.cpp
*
* Created on 12/11, 2010
*      Author: Zvi's script
*/
 
#include "R_2q2g1ph_eval.h"
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



 
 
template <class T> complex<T> R2q2g1ph_phqmqppm_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, p, m}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqppm LT");
#endif
 
 return( (pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*SPA(2,4)*SPA(2,5))/
(pow(SPA(2,3),3)*pow(SPB(3,2),2)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+
 SPA(3,5)*SPB(4,3),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*
 SPA(3,5))/(complex<T>(2,0)*pow(SPA(2,3),3)*pow(SPB(3,2),2)*SPA(3,4)*
 SPA(4,5))+(complex<T>(-19,0)*pow(SPA(2,5),2)*SPA(3,5))/
(complex<T>(9,0)*SPA(2,3)*SPA(3,4)*SPA(4,5))-
 (pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1)*SPA(2,3)*
 SPA(3,5))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*
 SPA(2,5)*SPA(4,5))/(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,3)*SPA(2,4))+
 (complex<T>(-2,0)*pow(SPA(2,5),2)*SPA(3,5)*SPB(3,2))/
(complex<T>(3,0)*S(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,2))/(complex<T>(2,0)*S(2,3)*SPA(3,4))+
 (complex<T>(5,0)*pow(SPA(2,5),2)*SPA(2,5)*SPB(4,2))/
(complex<T>(6,0)*pow(SPA(2,3),2)*SPA(4,5)*SPB(3,2))-
 (complex<T>(-2,0)*SPA(2,5)*SPA(3,5)*SPB(4,3))/(complex<T>(3,0)*S(2,3)*SPA(3,4))+
 (complex<T>(5,0)*pow(SPA(2,5),2)*SPA(3,5)*SPB(4,3))/
(complex<T>(6,0)*pow(SPA(2,3),2)*SPA(4,5)*SPB(3,2))+
 (complex<T>(2,0)*SPA(2,5)*SPA(4,5)*SPB(4,3))/(complex<T>(3,0)*SPA(2,4)*SPA(3,4)*
 SPB(3,2))+(complex<T>(1,0)*pow(SPA(2,5),2)*SPA(4,5)*SPB(4,3))/
(complex<T>(6,0)*SPA(2,4)*SPA(3,4)*(SPA(2,5)*SPB(3,2)-
SPA(4,5)*SPB(4,3)))-(SPA(3,5)*SPB(4,3))/
(complex<T>(2,0)*SPA(3,4)*SPB(5,2))+(complex<T>(-2,0)*pow(SPB(4,3),2)*SPB(4,2))/
(SPB(3,2)*SPB(5,2)*SPB(5,4))-(SPA(2,3)*SPB(4,2)*SPB(4,3))/
(complex<T>(2,0)*SPA(3,4)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),3))/(complex<T>(6,0)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPA(2,4)*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+
 SPA(3,5)*SPB(4,3),3))/(complex<T>(3,0)*pow(SPA(2,3),3)*
 pow(SPB(3,2),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPA(4,5)*
 SS(2,3,4))+(pow(SPB(4,3),2)*SPA(3,5)*SPA(4,5))/
(complex<T>(2,0)*SPA(3,4)*SPB(3,2)*SS(2,3,4))-
 (pow(SPB(4,3),2)*SPA(3,5)*SPA(4,5)*SPB(3,2))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(3,4)*SS(2,3,4))-
 (pow(SPA(2,5),2)*SPB(4,2))/(complex<T>(2,0)*SPA(3,4)*SS(2,3,4))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPB(4,2))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(3,4)*SS(2,3,4))-
 (pow(SPB(4,3),2)*SPA(2,5)*SPA(4,5)*SPB(4,2))/
(complex<T>(2,0)*pow(SPB(3,2),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*
 SPA(2,3)*SS(2,3,4))-(SPA(2,5)*SPA(3,5)*SPB(4,3))/
(complex<T>(2,0)*SPA(3,4)*SS(2,3,4))+(SPA(2,5)*SPA(3,5)*SPB(4,3))/
(complex<T>(3,0)*SPA(3,4)*SS(2,3,4))-(pow(SPB(4,3),2)*SPA(3,5)*
 SPA(4,5)*SPB(4,3))/(complex<T>(2,0)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPA(2,3)*SS(2,3,4))+
 (SPA(2,5)*SPA(4,5)*SPB(4,2)*SPB(4,3))/(complex<T>(2,0)*SPA(3,4)*SPB(3,2)*
 SS(2,3,4))+(SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2)*SPB(4,3))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(-2,0)*pow(SPA(2,5),3)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(3,0)*SPA(3,4)*SPA(4,5)*(SPA(2,5)*SPB(3,2)-
SPA(4,5)*SPB(4,3))*SS(2,3,4))-
 (complex<T>(-3,0)*pow(SPB(4,3),2)*SPA(2,5)*SPA(3,5)*SPA(4,5))/
(complex<T>(2,0)*SPA(3,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*
 SS(2,3,4))-(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),3)*SPA(3,5)*
 SPB(3,2))/(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(3,4)*
 (SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),3)*SPA(4,5)*SPB(4,2))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(3,4)*(SPA(2,5)*SPB(3,2)-
SPA(4,5)*SPB(4,3))*SS(2,3,4))-
 (complex<T>(-3,0)*pow(SPA(4,5),2)*pow(SPB(4,3),2)*SPA(2,5)*SPB(4,2))/
(complex<T>(2,0)*SPA(3,4)*SPB(3,2)*(SPA(2,5)*SPB(3,2)-
SPA(4,5)*SPB(4,3))*SS(2,3,4))-
 (complex<T>(-2,0)*pow(SPA(2,5),3)*SPB(3,2)*SPB(4,2))/
(complex<T>(3,0)*SPA(3,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*
 SS(2,3,4))-(complex<T>(2,0)*pow(SPA(2,5),2)*SPA(3,5)*SPB(3,2)*
 SPB(4,3))/(SPA(3,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*
 SS(2,3,4))-(complex<T>(2,0)*pow(SPA(2,5),2)*SPA(4,5)*SPB(4,2)*
 SPB(4,3))/(SPA(3,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*
 SS(2,3,4))-(pow(SPA(4,5),2)*pow(SPB(4,3),3)*SS(2,3,4))/
(complex<T>(6,0)*pow(SPA(2,3),2)*pow(SPB(3,2),4)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPA(2,4))+
 (complex<T>(-2,0)*pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+
 SPA(3,5)*SPB(4,3),3)*SS(2,3,4))/(complex<T>(3,0)*pow(SPA(2,3),5)*
 pow(SPB(3,2),4)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,5)*SPA(4,5)*SPB(4,2)*SS(2,3,4))/
(complex<T>(2,0)*pow(SPA(2,3),3)*pow(SPB(3,2),4)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(3,5)*SPA(4,5)*SPB(4,3)*SS(2,3,4))/
(complex<T>(2,0)*pow(SPA(2,3),3)*pow(SPB(3,2),4)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqppm_RT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, p, m}, RT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqppm RT");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(3,5),3)*pow(SPA(2,4)*SPB(4,3)+
 SPA(2,5)*SPB(5,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1))/
(complex<T>(2,0)*pow(SPA(4,5),3)*pow(SPB(5,4),2)*SPA(2,3)*SPA(3,4))-
 (pow(SPA(4,5),2)*pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),
-1))/(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,4)*SPA(3,4))+
 (pow(SPA(2,5),2)*SPA(3,5))/(SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(1,0)*SPA(2,5)*SPA(4,5)*SPB(4,3))/(complex<T>(2,0)*SPA(2,4)*SPA(3,4)*
 SPB(3,2))+(complex<T>(1,0)*SPA(3,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(2,0)*SPB(5,2)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2)))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*SPA(2,3)*SPB(4,3))/
(complex<T>(2,0)*SPB(5,2)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(3,5),2)*SPB(3,2))/
(complex<T>(2,0)*pow(SPA(4,5),2)*SPA(3,4)*(SPA(2,3)*SPB(4,2)+
SPA(3,5)*SPB(5,4)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(3,5),2)*
 SPB(4,3))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPA(2,3)*
 (SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4)))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(3,5),3)*SPB(5,3))/
(complex<T>(2,0)*pow(SPA(4,5),2)*SPA(2,3)*SPA(3,4)*(SPA(2,3)*SPB(4,2)+
SPA(3,5)*SPB(5,4)))+(complex<T>(1,0)*pow(SPA(3,5),3)*
 pow(SS(2,4,5),2)*SPB(4,3))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPA(3,4)*
 (SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*SPB(5,4)*
 (SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4)))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(3,5)*SPA(4,5))/
(complex<T>(2,0)*SPA(3,4)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,2))/(complex<T>(2,0)*SPA(3,4)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,5)*SPA(3,5)*SPB(4,3))/(complex<T>(2,0)*SPA(3,4)*
 SS(2,3,4))-(complex<T>(1,0)*SPA(2,5)*SPA(4,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(2,0)*SPA(3,4)*SPB(3,2)*SS(2,3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqppm_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, p, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqppm nfLT");
#endif
 
 return( (complex<T>(-2,0)*pow(SPA(2,5),2)*SPA(3,5))/(complex<T>(9,0)*SPA(2,3)*SPA(3,4)*
 SPA(4,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,3))/
(complex<T>(2,0)*S(2,3)*SPA(2,4))-(pow(SPA(4,5),2)*pow(SPB(4,3),3))/
(complex<T>(6,0)*pow(SPB(3,2),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*
 SPA(2,4)*SS(2,3,4))-(pow(SPA(2,4),2)*
 pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3))/
(complex<T>(6,0)*pow(SPA(2,3),3)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPA(4,5)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,5)*SPA(4,5))/
(complex<T>(6,0)*SPA(2,4)*SPB(3,2)*SS(2,3,4))-
 (pow(SPA(2,5),2)*SPA(2,5)*SPB(3,2))/(complex<T>(6,0)*SPA(2,4)*SPA(4,5)*
 SS(2,3,4))+(pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(6,0)*SPA(2,4)*
 SS(2,3,4))+(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(4,3),3)*
 SS(2,3,4))/(complex<T>(6,0)*pow(SPA(2,3),2)*pow(SPB(3,2),4)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPA(2,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+
 SPA(3,5)*SPB(4,3),3)*SS(2,3,4))/(complex<T>(6,0)*pow(SPA(2,3),5)*
 pow(SPB(3,2),4)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,5),3)*SS(2,3,4))/(complex<T>(6,0)*pow(SPA(2,3),2)*
 SPA(2,4)*SPA(4,5)*SPB(3,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqpmp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, m, p}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqpmp LT");
#endif
 
 return(-((pow(SPA(3,4),2)*pow(SPB(5,3),2)*
pow(complex<T>(1,0)-SS(2,3,5)/S(2,5),-1))/(complex<T>(2,0)*pow(SPB(5,2),2)*
SPA(2,5)*SPA(3,5)))+
 (complex<T>(1,0)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1)*SPA(2,5))/
(complex<T>(2,0)*pow(SPA(2,3),2)*pow(SPB(3,2),2)*SPA(3,5))-
 (pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(3,4),-1)*SPA(2,5))/
(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(3,5))+(complex<T>(-19,0)*pow(SPA(2,4),3))/
(complex<T>(9,0)*SPA(2,3)*SPA(2,5)*SPA(4,5))-
 (pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1)*SPA(2,4)*SPA(2,5))/
(complex<T>(3,0)*pow(SPA(2,3),3)*pow(SPB(3,2),2)*SPA(4,5))-
 (pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1)*SPA(2,3)*
 SPA(3,4))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(3,5)*SPA(4,5))+
 (complex<T>(-2,0)*pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1)*
 SPA(2,4)*SPA(4,5))/(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(2,3)*SPA(2,5))+
 (pow(SPB(5,3),2)*SPA(3,4)*SPA(4,5))/(S(2,5)*SPA(3,5)*SPB(3,2))-
 (pow(SPB(5,3),2)*SPA(2,4)*SPA(4,5))/(pow(SPA(2,5),2)*SPB(3,2)*
 SPB(5,2))+(complex<T>(1,0)*pow(SPA(2,4),2)*SPB(5,2))/
(complex<T>(2,0)*S(2,3)*SPA(3,5))-(pow(SPA(2,4),2)*SPA(2,4)*SPB(5,2))/
(complex<T>(3,0)*pow(SPA(2,3),2)*SPA(4,5)*SPB(3,2))+
 (complex<T>(5,0)*pow(SPA(2,4),2)*SPB(5,3))/(complex<T>(6,0)*S(2,3)*SPA(2,5))-
 (pow(SPA(2,4),2)*SPB(5,3))/(S(3,4)*SPA(2,5))-
 (pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,5),-1)*SPA(2,3)*SPB(5,3))/
(pow(SPA(2,5),3)*pow(SPB(5,2),2)*SPB(3,2))-
 (pow(SPA(2,4),2)*SPA(3,4)*SPB(5,3))/(complex<T>(3,0)*pow(SPA(2,3),2)*
 SPA(4,5)*SPB(3,2))-(SPA(2,4)*SPA(4,5)*SPB(5,3))/
(complex<T>(2,0)*SPA(2,5)*SPA(3,5)*SPB(3,2))-(pow(SPA(2,4),2)*SPB(5,3))/
(pow(SPA(2,5),2)*SPB(5,2))-(pow(SPB(5,3),2)*SPA(2,3)*SPA(2,4))/
(S(3,4)*SPA(2,5)*SPB(5,4))+(complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(2,0)*SPA(3,5)*SPB(4,3)*SPB(5,4))+(complex<T>(-2,0)*pow(SPB(5,3),3))/
(SPB(3,2)*SPB(4,3)*SPB(5,4))-
 (pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*
 pow(complex<T>(1,0)-SS(3,4,5)/S(3,4),-1)*SPA(4,5)*SPB(5,3))/
(pow(SPA(3,4),2)*pow(SPB(4,3),2)*SPA(2,5)*SPB(5,4))-
 (pow(SPA(4,5),2)*pow(SPB(5,3),3))/(complex<T>(3,0)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPB(5,3),2)*SPA(2,4)*SPA(4,5))/
(complex<T>(3,0)*SPA(2,5)*SPB(3,2)*SS(2,3,5))-
 (pow(SPB(5,3),2)*SPA(3,4)*SPA(4,5))/(complex<T>(2,0)*SPA(3,5)*SPB(3,2)*
 SS(2,3,5))+(pow(SPA(3,4),2)*pow(SPB(5,3),2))/
(complex<T>(2,0)*SPA(3,5)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(2,4),2)*SPB(5,3))/(complex<T>(3,0)*SPA(2,5)*SS(2,3,5))+
 (SPA(2,4)*SPA(3,4)*SPB(5,3))/(complex<T>(2,0)*SPA(3,5)*SS(2,3,5))-
 (SPA(2,4)*SPA(4,5)*SPB(5,2)*SPB(5,3))/(complex<T>(2,0)*SPA(3,5)*SPB(3,2)*
 SS(2,3,5))+(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,3),3)*
 SS(2,3,5))/(complex<T>(3,0)*pow(SPA(2,3),2)*pow(SPB(3,2),4)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPA(2,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqpmp_RT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, m, p}, RT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqpmp RT");
#endif
 
 return(-((pow(SPA(3,4),3)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),
 2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1))/
 (complex<T>(2,0)*pow(SPA(4,5),3)*pow(SPB(5,4),2)*SPA(2,3)*SPA(3,5)))-
 (pow(SPA(4,5),2)*pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),
-1))/(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,5)*SPA(3,5))-
 (pow(SPA(3,4),2)*pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,5),
-1))/(complex<T>(2,0)*pow(SPB(5,2),2)*SPA(2,5)*SPA(3,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(3,4),-1)*
 SPA(2,5))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(3,5))+
 pow(SPA(2,4),3)/(SPA(2,3)*SPA(2,5)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(3,4)*SPA(4,5))/
(complex<T>(2,0)*S(2,5)*SPA(3,5)*SPB(3,2))-
 (pow(SPB(5,3),3)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1)*SPA(2,3))/
(pow(SPB(5,4),2)*SPA(4,5)*SPB(4,3))+
 (SPA(2,3)*SPA(2,4)*SPA(4,5)*SPB(5,3))/(S(3,4)*SPA(2,5)*
 SPA(3,5))+(SPA(2,4)*SPB(5,3))/(SPA(3,5)*SPB(4,3))+
 (pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPA(4,5))/
(S(3,4)*SPA(2,5)*SPA(3,5)*SPB(5,4))-
 (pow(SPA(2,4),2)*pow(SPA(3,4),2)*SPB(4,3))/
(complex<T>(2,0)*pow(SPA(4,5),2)*SPA(2,3)*SPA(3,5)*SPB(5,4))-
 (pow(SPA(3,4),2)*SPA(2,4)*SPA(2,5)*SPB(5,3))/
(complex<T>(2,0)*pow(SPA(4,5),2)*SPA(2,3)*SPA(3,5)*SPB(5,4))+
 (pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*
 pow(complex<T>(1,0)-SS(3,4,5)/S(3,4),-1)*SPA(4,5)*SPB(5,3))/
(pow(SPA(3,4),2)*pow(SPB(4,3),2)*SPA(2,5)*SPB(5,4))-
 (pow(SPB(5,3),2)*SPA(2,4)*SPA(3,4))/(complex<T>(2,0)*S(2,5)*SS(2,3,5))-
 (pow(SPB(5,3),2)*SPA(2,4)*SPA(4,5)*SPB(5,2))/
(complex<T>(2,0)*S(2,5)*SPB(3,2)*SS(2,3,5))-
 (pow(SPA(2,4),2)*SPB(5,2)*SPB(5,3))/(complex<T>(2,0)*S(2,5)*SS(2,3,5))-
	(SPA(3,4)*SPA(4,5)*pow(SPB(5,3),3))/(complex<T>(2,0)*S(2,5)*SPB(3,2)*
 SS(2,3,5))+(complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(2,0)*S(4,5)*SPA(3,5)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,4))/(complex<T>(2,0)*SPB(4,3)*SS(3,4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,5)*SPA(3,4))/
(complex<T>(2,0)*SPA(3,5)*SPB(4,3)*SS(3,4,5))+
 (complex<T>(1,0)*SPA(2,4)*SPA(3,4)*SPB(5,3))/(complex<T>(2,0)*SPA(3,5)*
 SS(3,4,5))+(complex<T>(1,0)*pow(SPA(3,4),2)*SPA(2,4)*SPB(4,3)*
 SPB(5,3))/(complex<T>(2,0)*S(4,5)*SPA(3,5)*SS(3,4,5))-
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPB(5,3))/
(complex<T>(2,0)*SPB(4,3)*SPB(5,4)*SS(3,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqpmp_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, m, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqpmp nfLT");
#endif
 
 return( (complex<T>(-2,0)*pow(SPA(2,4),3))/(complex<T>(9,0)*SPA(2,3)*SPA(2,5)*
 SPA(4,5))-(pow(SPA(2,4),2)*SPB(5,3))/(complex<T>(2,0)*S(2,3)*
 SPA(2,5))+(complex<T>(1,0)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
(complex<T>(6,0)*pow(SPB(3,2),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*
 SPA(2,5)*SS(2,3,5))-(pow(SPA(2,5),2)*
 pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),3))/
(complex<T>(6,0)*pow(SPA(2,3),3)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPA(4,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,4)*SPA(4,5))/
(complex<T>(6,0)*SPA(2,5)*SPB(3,2)*SS(2,3,5))-
 (pow(SPA(2,4),2)*SPA(2,4)*SPB(3,2))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*
 SS(2,3,5))-(pow(SPA(2,4),2)*SPB(5,3))/(complex<T>(6,0)*SPA(2,5)*
 SS(2,3,5))-(pow(SPA(4,5),2)*pow(SPB(5,3),3)*SS(2,3,5))/
(complex<T>(6,0)*pow(SPA(2,3),2)*pow(SPB(3,2),4)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPA(2,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(2,4)*SPB(5,2)+
 SPA(3,4)*SPB(5,3),3)*SS(2,3,5))/(complex<T>(6,0)*pow(SPA(2,3),5)*
 pow(SPB(3,2),4)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPA(4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*SS(2,3,5))/(complex<T>(6,0)*pow(SPA(2,3),2)*
 SPA(2,5)*SPA(4,5)*SPB(3,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqpmm_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, m, m}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqpmm LT");
#endif
 
 return( (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5))/(complex<T>(4,0)*(S(2,5)+S(3,5))*
 SPA(2,3))-(complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/
(complex<T>(4,0)*(S(2,3)+S(2,4))*S(3,4))+(complex<T>(5,0)*pow(SPA(4,5),2))/
(complex<T>(12,0)*SPA(3,4)*SPB(4,2))-(complex<T>(1,0)*pow(SPB(3,2),2)*SPA(2,4)*
 SPA(2,5)*SPA(3,5))/(complex<T>(4,0)*(S(2,3)+S(2,4))*S(3,4)*SPB(4,2))-
 (SPA(2,4)*SPA(2,5)*SPA(3,5)*SPB(3,2))/(complex<T>(3,0)*S(3,4)*SPA(2,3)*
 SPB(4,2))+(complex<T>(5,0)*SPA(2,5)*SPA(4,5)*SPB(3,2))/
(complex<T>(6,0)*S(3,4)*SPB(4,2))-(S(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/
(complex<T>(12,0)*S(2,3)*(S(2,4)+S(3,4))*SPB(4,2))-
 (pow(SPA(2,4),2)*SPA(2,5)*SPA(4,5)*SPB(4,2))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2))/
(complex<T>(6,0)*S(2,3)*SPB(4,2)*SPB(4,3))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPA(2,4)*SPB(4,3))/
(complex<T>(4,0)*S(2,3)*(S(2,4)+S(3,4)))-
 (pow(SPA(2,4),2)*SPA(3,5)*SPA(4,5)*SPB(4,3))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPA(2,3))-
 (complex<T>(5,0)*pow(SPA(4,5),2)*SPB(4,3))/(complex<T>(6,0)*S(3,4)*SPB(4,2))+
 (pow(SPA(4,5),2)*S(2,4)*SPB(4,3))/(complex<T>(12,0)*S(2,3)*
 (S(2,4)+S(3,4))*SPB(4,2))-(complex<T>(1,0)*SPA(2,4)*SPA(3,5)*
 SPA(4,5)*SPB(3,2)*SPB(4,3))/(complex<T>(4,0)*S(2,3)*(S(2,4)+S(3,4))*
 SPB(4,2))-(SPA(2,4)*SPA(4,5))/(complex<T>(12,0)*SPA(2,3)*SPB(5,2))-
 (S(2,5)*SPA(2,4)*SPA(4,5)*SPB(3,2))/(complex<T>(6,0)*S(2,3)*
 (S(2,5)+S(3,5))*SPB(5,2))+(complex<T>(1,0)*pow(SPA(4,5),2)*SPA(2,5)*
 SPB(5,3))/(complex<T>(3,0)*pow(S(2,5)+S(3,5),2))-
 (pow(SPA(4,5),2)*SPB(5,3))/(complex<T>(6,0)*S(2,3)*SPB(5,2))-
 (pow(SPA(4,5),2)*S(2,5)*SPB(5,3))/(complex<T>(6,0)*S(2,3)*
 (S(2,5)+S(3,5))*SPB(5,2))+(complex<T>(1,0)*SPA(2,5)*SPA(3,4)*
 SPA(4,5)*SPB(5,3))/(complex<T>(4,0)*(S(2,5)+S(3,5))*SPA(2,3)*
 SPB(5,2))+(complex<T>(1,0)*S(2,5)*SPA(2,4)*SPA(3,5))/
(complex<T>(4,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SPB(5,4))+
 (complex<T>(1,0)*S(4,5)*SPA(2,4)*SPA(3,5))/(complex<T>(4,0)*SPA(2,3)*SPA(3,4)*
 SPB(4,2)*SPB(5,4))-(S(2,5)*SPA(2,4))/(complex<T>(3,0)*SPA(2,3)*
 SPB(5,2)*SPB(5,4))-(S(4,5)*SPA(2,4))/(complex<T>(3,0)*SPA(2,3)*
 SPB(5,2)*SPB(5,4))+(complex<T>(-5,0)*pow(SPA(2,5),2)*SPB(5,2))/
(complex<T>(12,0)*SPA(2,3)*SPB(4,2)*SPB(5,4))+
 (complex<T>(-19,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),2))/
(complex<T>(9,0)*SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))-
 (SPA(2,5)*SPA(4,5)*SPB(5,3))/(complex<T>(4,0)*(S(3,5)+S(4,5))*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,5),2)*SPB(5,3))/
(complex<T>(6,0)*SPA(2,3)*SPB(4,3)*SPB(5,4))+
 (pow(SPA(2,4),2)*SPA(3,5)*SPA(4,5))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*
 SS(2,3,4))-(complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/
(complex<T>(4,0)*(S(2,3)+S(2,4))*SS(2,3,4))-
 (pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPA(3,4))/
(complex<T>(6,0)*S(2,3)*SPB(4,2)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPB(3,2),2)*SPA(2,4)*SPA(2,5)*SPA(3,5))/
(complex<T>(4,0)*(S(2,3)+S(2,4))*SPB(4,2)*SS(2,3,4))+
 (pow(SPA(2,4),2)*pow(SPA(3,5),2)*SPB(3,2))/
(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*SPB(4,2)*SS(2,3,4))-
 (S(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/(complex<T>(12,0)*(S(2,4)+S(3,4))*
 SPB(4,2)*SS(2,3,4))-(pow(SPA(2,4),2)*SPA(2,5)*SPA(4,5)*
 SPB(3,2)*SPB(4,2))/(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SS(2,3,4))+
 (complex<T>(19,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPA(2,4))/
(complex<T>(9,0)*S(2,3)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2))/
(complex<T>(4,0)*SPB(4,2)*SPB(4,3)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPA(2,4)*SPB(4,3))/
(complex<T>(4,0)*(S(2,4)+S(3,4))*SS(2,3,4))-
 (pow(SPA(2,4),2)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SS(2,3,4))+
 (pow(SPA(4,5),2)*S(2,4)*SPB(4,3))/(complex<T>(12,0)*(S(2,4)+S(3,4))*
 SPB(4,2)*SS(2,3,4))-(complex<T>(1,0)*SPA(2,4)*SPA(3,5)*SPA(4,5)*
 SPB(3,2)*SPB(4,3))/(complex<T>(4,0)*(S(2,4)+S(3,4))*SPB(4,2)*
 SS(2,3,4))+(complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/
(complex<T>(4,0)*(S(2,5)+S(3,5))*SS(2,3,5))-
 (SPA(2,4)*SPA(4,5)*SPB(3,2))/(complex<T>(12,0)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*S(3,5)*SPA(2,4)*SPA(4,5)*SPB(3,2))/
(complex<T>(6,0)*S(2,3)*SPB(5,2)*SS(2,3,5))-
 (S(2,5)*SPA(2,4)*SPA(4,5)*SPB(3,2))/(complex<T>(6,0)*(S(2,5)+S(3,5))*
 SPB(5,2)*SS(2,3,5))+
 (complex<T>(-19,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*SPA(2,5))/
(complex<T>(9,0)*S(2,3)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(-19,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*SPA(2,5))/
(complex<T>(9,0)*S(2,5)*SPB(5,3)*SS(2,3,5))-(pow(SPA(4,5),2)*SPB(5,3))/
(complex<T>(12,0)*SPB(5,2)*SS(2,3,5))+(complex<T>(1,0)*pow(SPA(4,5),2)*S(3,5)*
 SPB(5,3))/(complex<T>(6,0)*S(2,3)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(4,5),2)*S(2,5)*SPB(5,3))/(complex<T>(6,0)*(S(2,5)+S(3,5))*
 SPB(5,2)*SS(2,3,5))+(complex<T>(1,0)*SPA(2,5)*SPA(3,4)*SPA(4,5)*
 SPB(3,2)*SPB(5,3))/(complex<T>(4,0)*(S(2,5)+S(3,5))*SPB(5,2)*
 SS(2,3,5))+(complex<T>(1,0)*pow(SPA(4,5),2)*SPA(2,5)*SPB(5,3)*
 SS(2,3,5))/(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*S(2,3))-
 (SPA(2,5)*SPA(4,5)*SPB(5,3)*SS(3,4,5))/(complex<T>(4,0)*S(3,4)*
 (S(3,5)+S(4,5))*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqpmm_RT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, m, m}, RT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqpmm RT");
#endif
 
 return(-((pow(SPA(3,5),2)*pow(SPB(3,2),2))/(complex<T>(2,0)*pow(SPB(4,2),3)*
SPA(3,4)))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(3,2),2)*
 pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),2))/
(complex<T>(4,0)*pow(SPB(4,2),3)*pow(SS(2,3,4),2)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*
 pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2))/
(complex<T>(4,0)*pow(SS(2,3,4),2)*(S(2,4)+S(3,4))*SPB(4,2)*SPB(4,3))-
 (pow(SPA(2,5),2)*SPB(3,2))/(complex<T>(2,0)*SPA(2,3)*SPB(4,2)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+
 SPA(3,5)*SPB(4,3),2)*SPB(3,2))/(complex<T>(4,0)*pow(SS(2,3,4),2)*
 SPA(2,3)*SPB(4,2)*SPB(4,3))+(complex<T>(1,0)*pow(SPA(2,3),2)*
 pow(SPB(3,2),2)*pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),2)*
 SPB(4,3))/(complex<T>(4,0)*pow(SPB(4,2),3)*pow(SS(2,3,4),2)*
 (S(2,3)+S(2,4)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(3,2),2)*
 pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),2))/
(complex<T>(4,0)*pow(SS(2,3,5),2)*(S(2,5)+S(3,5))*SPB(5,2)*SPB(5,3))-
 (pow(SPA(2,4),2)*SPB(3,2))/(complex<T>(2,0)*SPA(2,3)*SPB(5,2)*SPB(5,3))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(2,4)*SPB(5,2)+
 SPA(3,4)*SPB(5,3),2)*SPB(3,2))/(complex<T>(4,0)*pow(SS(2,3,5),2)*
 SPA(2,3)*SPB(5,2)*SPB(5,3))+
 pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),2)/
(SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,3),2)*pow(SS(3,4,5),2))/
(complex<T>(4,0)*pow(SPB(4,3),2)*(S(3,5)+S(4,5))*SPA(3,4)*SPB(5,4)*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))-
 (pow(SPA(2,5),2)*pow(SPB(5,3),2))/(complex<T>(2,0)*SPB(4,3)*SPB(5,4)*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))-
 (SPA(2,4)*SPA(2,5)*SPB(5,3))/(complex<T>(2,0)*SPB(5,4)*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(3,2),2)*
 pow(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2),2))/
(complex<T>(4,0)*pow(SPB(4,2),3)*(S(2,3)+S(2,4))*SPA(3,4)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(3,2),2)*pow(SPA(2,5)*SPB(4,2)+
 SPA(3,5)*SPB(4,3),2))/(complex<T>(2,0)*pow(SPB(4,2),3)*SPB(4,3)*
 SS(2,3,4))-(pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*
 SPA(2,4))/(S(2,3)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+
 SPA(3,5)*SPB(4,3),2)*SPB(3,2))/(complex<T>(4,0)*(S(2,4)+S(3,4))*
 SPA(2,3)*SPB(4,2)*SPB(4,3)*SS(2,3,4))+
 (pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*SPA(2,5))/
(S(2,3)*SPB(5,3)*SS(2,3,5))+
 (pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*SPA(2,5))/
(S(2,5)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2))/
(complex<T>(2,0)*SPB(5,2)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPA(2,4)*SPB(5,2)+
 SPA(3,4)*SPB(5,3),2)*SPB(3,2))/(complex<T>(4,0)*(S(2,5)+S(3,5))*
 SPA(2,3)*SPB(5,2)*SPB(5,3)*SS(2,3,5))-
 (pow(SPA(2,5),2)*pow(SPB(5,3),2)*SS(3,4,5))/
(complex<T>(4,0)*pow(SPB(4,3),2)*SPA(3,4)*SPB(5,4)*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,3),2)*SS(3,4,5))/
(complex<T>(4,0)*(S(3,5)+S(4,5))*SPB(4,3)*SPB(5,4)*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqpmm_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, m, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqpmm nfLT");
#endif
 
 return(-((SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/
 (complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))))-
 (SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/(complex<T>(6,0)*S(2,3)*
 (S(2,5)+S(3,5)))+(complex<T>(1,0)*pow(SPA(2,4),2)*SPA(2,5)*SPA(4,5)*
 SPB(4,2))/(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPA(2,3))+
 (pow(SPA(4,5),2)*SPA(2,4)*SPB(4,3))/(complex<T>(6,0)*S(2,3)*
 (S(2,4)+S(3,4)))+(complex<T>(1,0)*pow(SPA(2,4),2)*SPA(3,5)*SPA(4,5)*
 SPB(4,3))/(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPA(2,3))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPA(2,4)*SPA(4,5)*SPB(5,2))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPA(2,3))-
 (pow(SPA(4,5),2)*SPA(2,5)*SPB(5,3))/(complex<T>(6,0)*S(2,3)*
 (S(2,5)+S(3,5)))+(complex<T>(1,0)*pow(SPA(2,5),2)*SPA(3,4)*SPA(4,5)*
 SPB(5,3))/(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPA(2,3))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5))/(complex<T>(3,0)*SPA(2,3)*SPB(5,4))+
 (complex<T>(-2,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),2))/
(complex<T>(9,0)*SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/
(complex<T>(6,0)*S(2,3)*SS(2,3,4))-(SPA(2,4)*SPA(2,5)*SPA(4,5)*
 SPB(3,2))/(complex<T>(6,0)*(S(2,4)+S(3,4))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(2,5)*SPA(4,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)*SPA(2,4))/
(complex<T>(9,0)*S(2,3)*SPB(4,3)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPA(2,4)*SPB(4,3))/
(complex<T>(6,0)*S(2,3)*SS(2,3,4))+(pow(SPA(4,5),2)*SPA(2,4)*SPB(4,3))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))/
(complex<T>(6,0)*S(2,3)*SS(2,3,5))-(SPA(2,4)*SPA(2,5)*SPA(4,5)*
 SPB(3,2))/(complex<T>(6,0)*(S(2,5)+S(3,5))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPA(2,4)*SPA(4,5)*SPB(3,2)*SPB(5,2))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*SPA(2,5))/
(complex<T>(9,0)*S(2,3)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*SPA(2,5))/
(complex<T>(9,0)*S(2,5)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(4,5),2)*SPA(2,5)*SPB(5,3))/
(complex<T>(6,0)*S(2,3)*SS(2,3,5))-(pow(SPA(4,5),2)*SPA(2,5)*SPB(5,3))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqppp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, p, p}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqppp LT");
#endif
 
 return( (complex<T>(-2,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),2))/
(SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))-
 (SPA(2,3)*SPA(2,4)*SPB(4,3))/(complex<T>(2,0)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))-(complex<T>(1,0)*SPA(2,3)*SPB(5,3))/(complex<T>(2,0)*SPA(3,4)*
 SPA(4,5))-(SPA(2,4)*SPA(2,5)*SPB(5,4))/
(complex<T>(3,0)*pow(SPA(4,5),2)*SPA(2,3))-(complex<T>(1,0)*SPA(2,4)*SPB(5,4))/
(complex<T>(2,0)*SPA(3,4)*SPA(4,5))+
 (complex<T>(-2,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2))/
(SPA(2,4)*SPA(3,4)*SS(2,3,4))+
 (complex<T>(-2,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*
 SPB(4,3))/(S(2,3)*SPA(2,4)*SS(2,3,4))+
 (complex<T>(2,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPB(5,3))/
(S(2,3)*SPA(2,5)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqppp_RT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, p, p}, RT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqppp RT");
#endif
 
 return( (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPB(4,3))/(complex<T>(2,0)*SPA(2,5)*SPA(3,4)*
 SPA(4,5))+(SPA(2,3)*SPB(5,3))/(complex<T>(2,0)*SPA(3,4)*SPA(4,5))+
 (SPA(2,4)*SPB(5,4))/(complex<T>(2,0)*SPA(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmqppp_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, qp, p, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmqppp nfLT");
#endif
 
 return( (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(4,5),2)*
SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmpqpm_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, p, qp, m}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmpqpm LT");
#endif
 
 return( (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+
 SPA(3,5)*SPB(4,3),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/
(complex<T>(2,0)*pow(SPA(2,3),3)*pow(SPB(3,2),2)*SPA(3,4))+
 (complex<T>(-19,0)*pow(SPA(2,5),2))/(complex<T>(9,0)*SPA(2,3)*SPA(3,4))-
 (pow(SPB(4,3),2)*pow(complex<T>(1,0)-SS(3,4,5)/S(4,5),-1)*SPA(2,3))/
(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*S(2,4))/(complex<T>(2,0)*pow(SPA(2,3),2)*SPA(3,4)*
 SPB(3,2))+(SPA(2,5)*SPA(4,5)*SPB(4,3))/(S(2,3)*SPA(3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,3))/(complex<T>(2,0)*pow(SPA(2,3),2)*
 SPB(3,2))-(SPA(3,5)*SPA(4,5)*SPB(4,3))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2)))-
 (SPA(2,3)*SPA(4,5)*SPB(4,2)*SPB(4,3))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2))/
(complex<T>(2,0)*SPA(3,4)*SPB(3,2)*SS(2,3,4))-
 (SPB(4,2)*SPB(4,3)*SS(2,4,5))/(complex<T>(2,0)*SPB(5,2)*
 (SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2))*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmpqpm_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, p, qp, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmpqpm nfLT");
#endif
 
 return( (complex<T>(-2,0)*pow(SPA(2,5),2))/(complex<T>(9,0)*SPA(2,3)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmmqpp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, m, qp, p}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmmqpp LT");
#endif
 
 return(-(pow(SPA(2,3),2)/(SPA(2,5)*SPA(4,5)))+
 (complex<T>(1,0)*pow(SPA(3,4),2)*pow(SPB(5,4),2)*
 pow(complex<T>(1,0)-SS(2,4,5)/S(2,5),-1))/(complex<T>(2,0)*pow(SPB(5,2),2)*
 SPA(2,5)*SPA(4,5))-(pow(SPA(3,5),2)*
 pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*
 pow(complex<T>(1,0)-SS(3,4,5)/S(3,4),-1))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPA(2,5)*SPA(4,5))+(complex<T>(-2,0)*pow(SPB(5,4),2))/
(SPB(3,2)*SPB(4,3))-(pow(SPA(2,3),2)*SPA(3,5)*SPB(5,2))/
(complex<T>(2,0)*SPA(3,4)*SPA(4,5)*(SPA(2,5)*SPB(4,2)+
SPA(3,5)*SPB(4,3)))-(pow(SPA(2,3),2)*pow(SPA(3,5),2)*
 SPB(5,3))/(complex<T>(2,0)*SPA(2,5)*SPA(3,4)*SPA(4,5)*
 (SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3)))-
 (pow(SPA(2,3),2)*SPA(3,5)*SPB(5,4))/(complex<T>(2,0)*SPA(2,5)*SPA(3,4)*
 (SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3)))+
 (complex<T>(1,0)*pow(SPA(3,5),2)*pow(SS(2,3,4),2)*SPB(5,4))/
(complex<T>(2,0)*S(3,4)*SPA(4,5)*(-(SPA(3,5)*SPB(3,2))-
SPA(4,5)*SPB(4,2))*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3)))+
 (complex<T>(1,0)*SPB(4,2)*SPB(5,4)*SS(2,3,4))/(complex<T>(2,0)*SPB(3,2)*
 (-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2))/
(complex<T>(2,0)*SPA(4,5)*SPB(5,2)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phqmmqpp_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{ph, qm, m, qp, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phqmmqpp nfLT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1ph_phdqmqppm_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, p, m}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqppm LT");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(2,5),2)*SPA(3,5))/(SPA(2,3)*SPA(3,4)*
 SPA(4,5))-(pow(SPB(4,3),2)*SPA(3,5))/(complex<T>(2,0)*S(2,3)*
 SPB(5,2))+(SPA(2,5)*SPA(3,5)*SPB(3,2))/(complex<T>(2,0)*SPA(3,4)*
 SPA(4,5)*SPB(5,2))+(SPA(2,5)*SPB(4,2))/(complex<T>(2,0)*SPA(3,4)*
 SPB(5,2))+(complex<T>(-2,0)*SPA(2,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(3,0)*S(2,3)*SPB(5,2))+(complex<T>(-5,0)*pow(SPB(4,3),2)*SPA(2,5)*
 SPB(4,2))/(complex<T>(6,0)*pow(SPB(3,2),2)*SPA(2,3)*SPB(5,4))+
 (complex<T>(-5,0)*pow(SPB(4,3),2)*SPA(3,5)*SPB(4,3))/
(complex<T>(6,0)*pow(SPB(3,2),2)*SPA(2,3)*SPB(5,4))-
 (pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*pow(SPB(5,3),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1)*SPB(4,2))/
(complex<T>(2,0)*pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPB(5,2)*SPB(5,4))-
 (complex<T>(-2,0)*pow(SPB(4,3),2)*SPA(2,3)*SPB(4,2))/
(complex<T>(3,0)*S(2,3)*SPB(5,2)*SPB(5,4))+
 (complex<T>(19,0)*pow(SPB(4,3),2)*SPB(4,2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,2)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,5),2)*
 pow(complex<T>(1,0)-SS(2,4,5)/S(4,5),-1)*SPB(3,2)*SPB(4,2))/
(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1)*SPB(4,3)*SPB(5,3))/
(pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPB(5,4))-
 (pow(SPA(2,5),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1)*SPB(4,3)*
 SPB(5,4))/(complex<T>(2,0)*pow(SPA(2,3),2)*SPB(3,2)*SPB(5,3))+
 (complex<T>(-2,0)*SPA(2,5)*SPB(4,3)*SPB(5,4))/(complex<T>(3,0)*SPA(2,3)*SPB(5,2)*
 SPB(5,3))+(complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,5)*SPB(5,4))/
(complex<T>(6,0)*SPB(5,2)*SPB(5,3)*(-(SPA(2,3)*SPB(4,3))+
SPA(2,5)*SPB(5,4)))+(pow(SPB(4,3),2)*SPA(3,5))/
(complex<T>(2,0)*SPB(5,2)*SS(2,3,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*
 pow(SPB(5,4),2)*SPA(3,5))/(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(5,2)*
 SS(2,3,5))+(SPA(2,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(2,0)*SPB(5,2)*SS(2,3,5))-(complex<T>(1,0)*SPA(2,5)*SPB(4,2)*
 SPB(4,3))/(complex<T>(3,0)*SPB(5,2)*SS(2,3,5))-
 (pow(SPA(2,5),3)*pow(SPB(5,4),2))/(complex<T>(6,0)*pow(SPA(2,3),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3)*
 pow(SPB(5,3),2))/(complex<T>(3,0)*pow(SPA(2,3),2)*pow(SPB(3,2),3)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPB(5,4)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*SPA(2,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(2,0)*pow(SPA(2,3),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*
 SPB(3,2)*SS(2,3,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*SPA(3,5)*
 SPB(4,3)*SPB(5,4))/(complex<T>(2,0)*pow(SPA(2,3),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPB(3,2)*SS(2,3,5))-
 (pow(SPA(2,5),2)*SPB(4,2)*SPB(5,4))/(complex<T>(2,0)*SPA(2,3)*SPB(5,2)*
 SS(2,3,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*SPA(2,3)*SPB(4,2)*
 SPB(5,4))/(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(5,2)*SS(2,3,5))-
 (SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/(complex<T>(2,0)*SPA(2,3)*SPB(5,2)*
 SS(2,3,5))-(complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPA(3,5)*SPB(4,3)*
 SPB(5,4))/(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(2,0)*pow(SPB(4,3),3)*SPA(2,3)*SPA(3,5))/
(complex<T>(3,0)*SPB(5,2)*(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4))*
 SS(2,3,5))-(pow(SPA(2,5),3)*pow(SPB(5,4),2)*SPA(2,3)*
 SPB(4,2))/(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(5,2)*
 (-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4))*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPB(4,3),2)*SPA(2,3)*SPA(2,5)*SPB(4,2))/
(SPB(5,2)*(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4))*
 SS(2,3,5))+(complex<T>(3,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*SPA(3,5)*
 SPB(4,3))/(complex<T>(2,0)*SPA(2,3)*SPB(5,2)*(-(SPA(2,3)*SPB(4,3))+
SPA(2,5)*SPB(5,4))*SS(2,3,5))+
 (complex<T>(2,0)*pow(SPA(2,3),2)*pow(SPB(4,3),3)*SPB(4,2))/
(complex<T>(3,0)*SPB(5,2)*SPB(5,4)*(-(SPA(2,3)*SPB(4,3))+
SPA(2,5)*SPB(5,4))*SS(2,3,5))-
 (pow(SPA(2,5),3)*pow(SPB(5,4),2)*SPA(3,5)*SPB(5,4))/
(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(5,2)*(-(SPA(2,3)*SPB(4,3))+
SPA(2,5)*SPB(5,4))*SS(2,3,5))+
 (complex<T>(-2,0)*pow(SPB(4,3),2)*SPA(2,5)*SPA(3,5)*SPB(5,4))/
(SPB(5,2)*(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4))*
 SS(2,3,5))+(complex<T>(3,0)*pow(SPA(2,5),2)*SPB(4,2)*SPB(4,3)*
 SPB(5,4))/(complex<T>(2,0)*SPB(5,2)*(-(SPA(2,3)*SPB(4,3))+
SPA(2,5)*SPB(5,4))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5),3)*pow(SPB(5,4),2)*SS(2,3,5))/
(complex<T>(6,0)*pow(SPA(2,3),4)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPB(5,3))+
 (complex<T>(2,0)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3)*
 pow(SPB(5,3),2)*SS(2,3,5))/(complex<T>(3,0)*pow(SPA(2,3),4)*
 pow(SPB(3,2),5)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPB(5,4))-
 (pow(SPA(2,5),2)*SPA(2,5)*SPB(4,2)*SPB(5,4)*SS(2,3,5))/
(complex<T>(2,0)*pow(SPA(2,3),4)*pow(SPB(3,2),3)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3))-
 (pow(SPA(2,5),2)*SPA(3,5)*SPB(4,3)*SPB(5,4)*SS(2,3,5))/
(complex<T>(2,0)*pow(SPA(2,3),4)*pow(SPB(3,2),3)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqppm_RT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, p, m}, RT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqppm RT");
#endif
 
 return(-((pow(SPB(4,2),3)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),
 2)*pow(complex<T>(1,0)-SS(2,4,5)/S(4,5),-1))/
 (complex<T>(2,0)*pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPB(3,2)*SPB(5,2)))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*pow(SPB(4,3),2)*SPA(2,3))/
(complex<T>(2,0)*pow(SPB(5,4),2)*(-(SPA(3,5)*SPB(3,2))-
SPA(4,5)*SPB(4,2))*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*pow(SPB(4,3),2)*S(2,4))/
(complex<T>(2,0)*pow(SPB(5,4),2)*SPB(3,2)*(-(SPA(3,5)*SPB(3,2))-
SPA(4,5)*SPB(4,2))*SPB(5,2))+
 (complex<T>(1,0)*pow(SPB(4,2),2)*pow(SPB(4,3),2)*S(2,5))/
(complex<T>(2,0)*pow(SPB(5,4),2)*SPB(3,2)*(-(SPA(3,5)*SPB(3,2))-
SPA(4,5)*SPB(4,2))*SPB(5,2))-
 (complex<T>(1,0)*pow(SPA(3,5),2)*SPA(2,5)*SPB(3,2))/
(complex<T>(2,0)*SPA(3,4)*SPA(4,5)*(SPA(3,4)*SPB(4,2)+
SPA(3,5)*SPB(5,2)))-(complex<T>(1,0)*SPA(2,5)*SPA(3,5)*SPB(4,2))/
(complex<T>(2,0)*SPA(3,4)*(SPA(3,4)*SPB(4,2)+SPA(3,5)*SPB(5,2)))+
 (complex<T>(1,0)*pow(SPB(4,2),3)*pow(SS(3,4,5),2)*SPA(2,5))/
(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(4,5)*(-(SPA(3,5)*SPB(3,2))-
SPA(4,5)*SPB(4,2))*SPB(5,2)*(SPA(3,4)*SPB(4,2)+
SPA(3,5)*SPB(5,2)))+(complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(5,4),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1))/(complex<T>(2,0)*pow(SPA(2,3),2)*
 SPB(5,2)*SPB(5,3))-(pow(SPB(4,3),2)*SPB(4,2))/
(SPB(3,2)*SPB(5,2)*SPB(5,4))-(SPA(2,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(2,0)*SPA(2,3)*SPB(5,2)*SPB(5,3))-
 (complex<T>(1,0)*pow(SPB(4,3),2)*SPA(3,5))/(complex<T>(2,0)*SPB(5,2)*SS(2,3,5))-
 (complex<T>(1,0)*SPA(2,5)*SPB(4,2)*SPB(4,3))/(complex<T>(2,0)*SPB(5,2)*
 SS(2,3,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*SPB(4,2)*SPB(5,4))/
(complex<T>(2,0)*SPA(2,3)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,5)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(2,0)*SPA(2,3)*SPB(5,2)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqppm_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, p, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqppm nfLT");
#endif
 
 return(-((pow(SPB(4,3),2)*SPA(2,5))/(complex<T>(2,0)*S(2,3)*SPB(5,3)))+
 (complex<T>(2,0)*pow(SPB(4,3),2)*SPB(4,2))/(complex<T>(9,0)*SPB(3,2)*SPB(5,2)*
 SPB(5,4))+(complex<T>(1,0)*pow(SPA(2,5),3)*pow(SPB(5,4),2))/
(complex<T>(6,0)*pow(SPA(2,3),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*
 SPB(5,3)*SS(2,3,5))-(pow(SPB(4,3),2)*SPA(2,5))/
(complex<T>(6,0)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3)*
 pow(SPB(5,3),2))/(complex<T>(6,0)*pow(SPA(2,3),2)*pow(SPB(3,2),3)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPB(5,4)*SS(2,3,5))+
 (pow(SPB(4,3),2)*SPA(2,3)*SPB(4,3))/(complex<T>(6,0)*SPB(5,3)*SPB(5,4)*
 SS(2,3,5))-(pow(SPA(2,5),2)*SPB(4,3)*SPB(5,4))/
(complex<T>(6,0)*SPA(2,3)*SPB(5,3)*SS(2,3,5))-
 (pow(SPA(2,5),3)*pow(SPB(5,4),2)*SS(2,3,5))/
(complex<T>(6,0)*pow(SPA(2,3),4)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPB(5,3))-
 (pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),3)*pow(SPB(5,3),2)*
 SS(2,3,5))/(complex<T>(6,0)*pow(SPA(2,3),4)*pow(SPB(3,2),5)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),3)*SPB(5,4))-
 (pow(SPB(4,3),3)*SS(2,3,5))/(complex<T>(6,0)*pow(SPB(3,2),2)*SPA(2,3)*
 SPB(5,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqpmp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, m, p}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqpmp LT");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(2,4),3))/(SPA(2,3)*SPA(2,5)*SPA(4,5))+
 (pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(3,4),-1)*SPA(2,4)*SPB(3,2))/
(pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPA(2,3))-
 (pow(SPB(5,3),2)*SPA(3,4))/(complex<T>(2,0)*S(2,3)*SPB(4,2))-
 (pow(SPA(2,4),2)*SPB(3,2))/(complex<T>(2,0)*SPA(2,5)*SPA(4,5)*SPB(4,2))+
 (complex<T>(-5,0)*pow(SPB(5,3),2)*SPA(2,4))/(complex<T>(6,0)*S(2,3)*SPB(4,3))+
 (pow(SPB(5,3),2)*SPA(2,4))/(S(2,5)*SPB(4,3))+
 (pow(SPB(5,3),2)*SPA(2,4))/(S(3,4)*SPB(4,3))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(3,4),-1))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 SPB(4,2)*SPB(4,3))-(pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),
2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*SPB(4,3))/
(complex<T>(2,0)*pow(SPA(2,3),2)*pow(SPB(3,2),2)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(complex<T>(1,0)-SS(2,4,5)/S(2,5),-1)*
 SPB(4,3))/(complex<T>(2,0)*pow(SPA(2,5),2)*SPB(4,2))+
 (pow(SPA(2,4),2)*SPB(3,2)*SPB(5,3))/(S(2,5)*SPA(4,5)*SPB(4,3))+
 (complex<T>(19,0)*pow(SPB(5,3),3))/(complex<T>(9,0)*SPB(3,2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,4)*SPB(5,2))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(2,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(complex<T>(1,0)-SS(2,4,5)/S(4,5),-1)*
 SPB(3,2)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(4,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(3,4)*SPB(5,3))/
(complex<T>(3,0)*pow(SPB(3,2),2)*SPA(2,3)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*SPB(4,3)*SPB(5,3))/
(complex<T>(3,0)*pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPB(5,4))+
 (pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*
 pow(complex<T>(1,0)-SS(2,4,5)/S(2,5),-1)*SPA(2,4)*SPB(5,4))/
(pow(SPA(2,5),2)*pow(SPB(5,2),2)*SPA(4,5)*SPB(4,3))-
 (pow(SPA(2,4),2)*SPB(5,2)*SPB(5,4))/(S(3,4)*SPA(2,3)*SPB(4,2))+
 (pow(SPA(2,4),2)*SPB(5,3)*SPB(5,4))/(S(3,4)*SPA(2,3)*SPB(4,3))+
 (complex<T>(2,0)*pow(SPA(2,4),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1)*
 SPB(5,3)*SPB(5,4))/(complex<T>(3,0)*pow(SPA(2,3),2)*SPB(3,2)*SPB(4,3))+
 (complex<T>(1,0)*SPA(2,4)*SPB(5,3)*SPB(5,4))/(complex<T>(2,0)*SPA(2,3)*SPB(4,2)*
 SPB(4,3))-(pow(SPA(2,4),2)*pow(SPB(5,2),2))/
(complex<T>(2,0)*SPA(3,4)*SPB(4,2)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(5,4),2))/(complex<T>(3,0)*pow(SPA(2,3),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,4))/(complex<T>(3,0)*SPB(4,3)*SS(2,3,4))-
 (SPA(2,4)*SPB(5,2)*SPB(5,3))/(complex<T>(2,0)*SPB(4,2)*SS(2,3,4))+
 (pow(SPA(2,4),2)*SPB(5,2)*SPB(5,4))/(complex<T>(2,0)*SPA(2,3)*SPB(4,2)*
 SS(2,3,4))+(SPA(2,4)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(2,0)*SPA(2,3)*SPB(4,2)*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,4),2)*SPB(5,3)*SPB(5,4))/
(complex<T>(3,0)*SPA(2,3)*SPB(4,3)*SS(2,3,4))-
 (pow(SPA(2,4),3)*pow(SPB(5,4),2)*SS(2,3,4))/
(complex<T>(3,0)*pow(SPA(2,3),4)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqpmp_RT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, m, p}, RT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqpmp RT");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(5,2),3)*pow(SPA(2,4)*SPB(4,3)+
 SPA(2,5)*SPB(5,3),2)*pow(complex<T>(1,0)-SS(2,4,5)/S(4,5),-1))/
(complex<T>(2,0)*pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPB(3,2)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPB(5,2),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(4,5)*SPB(3,2)*SPB(4,2))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/(complex<T>(2,0)*pow(SPA(2,3),2)*
 SPB(4,2)*SPB(4,3))+(complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(3,4),-1))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 SPB(4,2)*SPB(4,3))-(pow(SPA(2,4),2)*
 pow(complex<T>(1,0)-SS(2,4,5)/S(2,5),-1)*SPB(4,3))/
(complex<T>(2,0)*pow(SPA(2,5),2)*SPB(4,2))-(SPA(2,4)*SPB(5,3))/
(SPA(2,5)*SPB(4,2))+(complex<T>(1,0)*pow(SPB(5,2),2)*SPA(2,4)*SPB(4,3)*
 SPB(5,3))/(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(4,5)*SPB(3,2)*SPB(4,2))+
 (pow(SPA(2,4),3)*pow(complex<T>(1,0)-SS(2,4,5)/S(4,5),-1)*SPB(3,2))/
(pow(SPA(4,5),2)*SPA(2,5)*SPB(5,4))-pow(SPB(5,3),3)/
(SPB(3,2)*SPB(4,3)*SPB(5,4))-
 (pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*
 pow(complex<T>(1,0)-SS(2,4,5)/S(2,5),-1)*SPA(2,4)*SPB(5,4))/
(pow(SPA(2,5),2)*pow(SPB(5,2),2)*SPA(4,5)*SPB(4,3))-
 (pow(SPA(2,4),2)*pow(SPB(3,2),2)*SPB(5,4))/
(S(2,5)*SPA(4,5)*SPB(4,2)*SPB(4,3))-
 (pow(SPA(2,4),2)*SPB(5,2)*SPB(5,4))/(complex<T>(2,0)*S(3,4)*SPA(2,3)*
 SPB(4,2))-(SPA(2,4)*SPB(3,2)*SPB(5,3)*SPB(5,4))/
(S(2,5)*SPB(4,2)*SPB(4,3))+(pow(SPB(5,3),2)*SPA(2,4)*SPA(3,4))/
(complex<T>(2,0)*S(3,4)*SS(2,3,4))+(pow(SPA(2,4),2)*SPB(5,2)*SPB(5,3))/
	 (complex<T>(2,0)*S(3,4)*SS(2,3,4))+(pow(SPA(2,4),3)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*S(3,4)*SPA(2,3)*SS(2,3,4))+
 (pow(SPA(2,4),2)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(2,0)*S(3,4)*SPA(2,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPA(2,4)*SPB(3,2))/
(complex<T>(2,0)*SPA(2,5)*SPA(4,5)*SS(2,4,5))-
 (pow(SPA(2,4),2)*SPB(4,3)*SPB(5,2))/(complex<T>(2,0)*SPA(2,5)*SPB(4,2)*
 SS(2,4,5))-(pow(SPA(2,4),2)*S(2,5)*SPB(4,3)*SPB(5,2))/
(complex<T>(2,0)*S(4,5)*SPA(2,5)*SPB(4,2)*SS(2,4,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*SPB(5,3))/(complex<T>(2,0)*SPA(2,5)*SS(2,4,5))-
 (SPA(2,4)*SPB(5,2)*SPB(5,3))/(complex<T>(2,0)*SPB(4,2)*SS(2,4,5))-
 (S(2,5)*SPA(2,4)*SPB(5,2)*SPB(5,3))/(complex<T>(2,0)*S(4,5)*SPB(4,2)*
 SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqpmp_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, m, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqpmp nfLT");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,4))/(complex<T>(2,0)*S(2,3)*SPB(4,3))+
 (complex<T>(2,0)*pow(SPB(5,3),3))/(complex<T>(9,0)*SPB(3,2)*SPB(4,3)*SPB(5,4))-
 (pow(SPA(2,4),3)*pow(SPB(5,4),2))/(complex<T>(6,0)*pow(SPA(2,3),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPB(4,3)*SS(2,3,4))+
 (pow(SPB(5,3),2)*SPA(2,4))/(complex<T>(6,0)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPB(4,3),2)*pow(SPA(2,4)*SPB(5,2)+
 SPA(3,4)*SPB(5,3),3))/(complex<T>(6,0)*pow(SPA(2,3),2)*
 pow(SPB(3,2),3)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPB(5,4)*
 SS(2,3,4))+(pow(SPB(5,3),2)*SPA(2,3)*SPB(5,3))/
(complex<T>(6,0)*SPB(4,3)*SPB(5,4)*SS(2,3,4))-
 (pow(SPA(2,4),2)*SPB(5,3)*SPB(5,4))/(complex<T>(6,0)*SPA(2,3)*SPB(4,3)*
 SS(2,3,4))+(complex<T>(1,0)*pow(SPA(2,4),3)*pow(SPB(5,4),2)*
 SS(2,3,4))/(complex<T>(6,0)*pow(SPA(2,3),4)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPB(4,3))-
 (pow(SPB(4,3),2)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),3)*
 SS(2,3,4))/(complex<T>(6,0)*pow(SPA(2,3),4)*pow(SPB(3,2),5)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),3)*SPB(5,4))-
 (pow(SPB(5,3),3)*SS(2,3,4))/(complex<T>(6,0)*pow(SPB(3,2),2)*SPA(2,3)*
 SPB(4,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqpmm_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, m, m}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqpmm LT");
#endif
 
 return( (complex<T>(1,0)*SPA(4,5)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*
 SPB(3,2))+(complex<T>(1,0)*SPA(2,4)*SPB(3,2))/(complex<T>(2,0)*SPB(5,2)*
 SPB(5,4))+(complex<T>(2,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),
2))/(SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))+
 (complex<T>(1,0)*SPA(4,5)*SPB(5,3))/(complex<T>(2,0)*SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*SPA(2,5)*SPB(3,2)*SPB(5,3))/(complex<T>(2,0)*SPB(4,3)*SPB(5,2)*
 SPB(5,4))+(complex<T>(-2,0)*pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),
2)*SPA(2,4))/(S(2,3)*SPB(4,3)*SS(2,3,4))+
 (complex<T>(2,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*SPA(2,5))/
(S(2,3)*SPB(5,3)*SS(2,3,5))+
 (complex<T>(2,0)*pow(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3),2)*SPA(2,5))/
(S(2,5)*SPB(5,3)*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqpmm_RT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, m, m}, RT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqpmm RT");
#endif
 
 return(-((SPA(2,4)*SPB(3,2))/(complex<T>(2,0)*SPB(5,2)*SPB(5,4)))-
 (SPA(4,5)*SPB(5,3))/(complex<T>(2,0)*SPB(5,2)*SPB(5,4))-
 (SPA(2,5)*SPB(3,2)*SPB(5,3))/(complex<T>(2,0)*SPB(4,3)*SPB(5,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqpmm_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, m, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqpmm nfLT");
#endif
 
 return(-((SPA(4,5)*SPB(4,3)*SPB(5,3))/(complex<T>(3,0)*pow(SPB(5,4),2)*
 SPB(3,2)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqppp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, p, p}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqppp LT");
#endif
 
 return( (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,4))/(complex<T>(6,0)*S(2,3)*SPA(3,4))+
 (pow(SPB(5,4),2)*S(3,4)*SPA(2,4))/(complex<T>(6,0)*S(2,3)*
 (S(2,4)+S(3,4))*SPA(3,4))-
 pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)/
(complex<T>(6,0)*S(2,3)*SPA(2,5)*SPA(3,5))+
 (complex<T>(5,0)*pow(SPB(5,4),2)*SPA(2,5))/(complex<T>(6,0)*S(2,5)*SPA(3,5))-
 (pow(SPB(5,4),2)*S(3,5)*SPA(2,5))/(complex<T>(12,0)*S(2,3)*
 (S(2,5)+S(3,5))*SPA(3,5))+
 (complex<T>(19,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),2))/
(complex<T>(9,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))-
 (pow(SPB(4,3),2)*SPA(2,4))/(complex<T>(6,0)*SPA(2,5)*SPA(4,5)*SPB(3,2))+
 (complex<T>(5,0)*pow(SPB(4,3),2)*SPA(3,4))/(complex<T>(12,0)*SPA(3,5)*SPA(4,5)*
 SPB(3,2))-(pow(SPB(5,4),2)*SPA(2,4)*SPB(4,3))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2))+(complex<T>(-5,0)*pow(SPB(5,4),2))/
(complex<T>(12,0)*SPA(3,5)*SPB(5,2))+(complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,5)*
 SPB(5,3))/(complex<T>(4,0)*S(2,3)*(S(2,5)+S(3,5)))+
 (complex<T>(1,0)*S(3,4)*SPB(5,3))/(complex<T>(3,0)*SPA(3,4)*SPA(4,5)*SPB(3,2))+
 (complex<T>(1,0)*S(4,5)*SPB(5,3))/(complex<T>(3,0)*SPA(3,4)*SPA(4,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(4,2)*SPB(4,3)*SPB(5,3))/
(complex<T>(4,0)*S(2,5)*(S(2,3)+S(3,5))*SPA(3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,3))/
(complex<T>(3,0)*S(2,5)*SPA(3,5)*SPB(3,2))-(S(3,4)*SPB(4,2)*SPB(5,3))/
(complex<T>(4,0)*SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,2))-
 (S(4,5)*SPB(4,2)*SPB(5,3))/(complex<T>(4,0)*SPA(3,5)*SPA(4,5)*SPB(3,2)*
 SPB(5,2))+(complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPB(3,2))-
 (complex<T>(5,0)*SPA(2,3)*SPB(4,3)*SPB(5,4))/(complex<T>(6,0)*S(2,5)*SPA(3,5))+
 (S(3,5)*SPA(2,3)*SPB(4,3)*SPB(5,4))/(complex<T>(12,0)*S(2,3)*
 (S(2,5)+S(3,5))*SPA(3,5))+(complex<T>(1,0)*SPA(2,4)*SPB(4,3)*
 SPB(5,4))/(complex<T>(4,0)*(S(2,4)+S(4,5))*SPA(4,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPB(3,2))-
 (SPA(2,4)*SPB(4,3)*SPB(5,2)*SPB(5,4))/(complex<T>(4,0)*(S(2,4)+S(3,4))*
 SPA(3,4)*SPB(3,2))+(S(3,4)*SPA(2,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*(S(2,4)+S(3,4))*SPA(3,4))+
 (complex<T>(1,0)*SPB(5,3)*SPB(5,4))/(complex<T>(12,0)*SPA(3,4)*SPB(3,2))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(4,0)*S(2,3)*(S(2,5)+S(3,5))*SPA(3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(4,0)*S(2,5)*(S(2,3)+S(3,5)))-(SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(4,0)*(S(2,4)+S(3,4))*SPB(3,2))+
 (complex<T>(19,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2))/
(complex<T>(9,0)*SPA(2,4)*SPA(3,4)*SS(2,3,4))+
 (pow(SPB(5,4),2)*SPA(2,4))/(complex<T>(12,0)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPB(5,4),2)*S(2,4)*SPA(2,4))/(complex<T>(6,0)*S(2,3)*SPA(3,4)*
 SS(2,3,4))+(pow(SPB(5,4),2)*S(3,4)*SPA(2,4))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SPA(3,4)*SS(2,3,4))+
 (complex<T>(19,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*
 SPB(4,3))/(complex<T>(9,0)*S(2,3)*SPA(2,4)*SS(2,3,4))-
 (SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,2)*SPB(5,4))/
(complex<T>(4,0)*(S(2,4)+S(3,4))*SPA(3,4)*SS(2,3,4))+
 (SPA(2,3)*SPB(5,3)*SPB(5,4))/(complex<T>(12,0)*SPA(3,4)*SS(2,3,4))-
 (complex<T>(1,0)*S(2,4)*SPA(2,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*SPA(3,4)*SS(2,3,4))+
 (S(3,4)*SPA(2,3)*SPB(5,3)*SPB(5,4))/(complex<T>(6,0)*(S(2,4)+S(3,4))*
 SPA(3,4)*SS(2,3,4))-(SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(4,0)*(S(2,4)+S(3,4))*SS(2,3,4))-
 (pow(SPB(5,4),2)*SPA(2,4)*SPB(4,3)*SS(2,3,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*S(2,3))-
 pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)/
(complex<T>(4,0)*SPA(2,5)*SPA(3,5)*SS(2,3,5))-
 (pow(SPB(5,4),2)*S(3,5)*SPA(2,5))/(complex<T>(12,0)*(S(2,5)+S(3,5))*
 SPA(3,5)*SS(2,3,5))-(pow(SPB(4,2),2)*pow(SPB(5,3),2)*
 SPA(2,3))/(complex<T>(6,0)*SPA(3,5)*SPB(3,2)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*SPB(5,2))/
(complex<T>(6,0)*S(2,3)*SPA(3,5)*SS(2,3,5))+
 (complex<T>(-19,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*
 SPB(5,3))/(complex<T>(9,0)*S(2,3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,5)*SPB(5,3))/
(complex<T>(4,0)*(S(2,5)+S(3,5))*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*SPB(4,2)*SPB(4,3)*SPB(5,3))/
(complex<T>(4,0)*(S(2,3)+S(3,5))*SPA(3,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SS(2,3,5))+
 (S(3,5)*SPA(2,3)*SPB(4,3)*SPB(5,4))/(complex<T>(12,0)*(S(2,5)+S(3,5))*
 SPA(3,5)*SS(2,3,5))+(complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3)*
 SPA(3,5)*SPB(4,3)*SPB(5,4))/(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*
 SS(2,3,5))-(pow(SPB(5,3),2)*SPB(4,2)*SPB(5,4))/
(complex<T>(6,0)*SPB(3,2)*SPB(5,2)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,5)*SPB(4,2)*SPB(5,3)*SPB(5,4))/
(complex<T>(4,0)*(S(2,5)+S(3,5))*SPA(3,5)*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(4,0)*(S(2,3)+S(3,5))*SS(2,3,5))+
 (complex<T>(1,0)*SPA(2,4)*SPB(4,3)*SPB(5,4)*SS(2,4,5))/
(complex<T>(4,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqppp_RT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, p, p}, RT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqppp RT");
#endif
 
 return(-((pow(SPA(2,3),2)*pow(SPB(3,2),2)*
pow(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4),2)*SPA(2,5))/
 (complex<T>(4,0)*pow(SPA(3,5),3)*pow(SS(2,3,5),2)*(S(2,3)+S(3,5))))-
 (pow(SPA(2,3),2)*pow(SPB(4,3),2)*pow(SPA(2,4)*SPB(5,2)+
 SPA(3,4)*SPB(5,3),2))/(complex<T>(4,0)*pow(SS(2,3,4),2)*
 (S(2,4)+S(3,4))*SPA(2,4)*SPA(3,4))-
 (pow(SPA(2,3),2)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*
 pow(SPB(5,3),2))/(complex<T>(4,0)*pow(SS(2,3,5),2)*(S(2,5)+S(3,5))*
 SPA(2,5)*SPA(3,5))+(complex<T>(1,0)*pow(SPB(5,3),2)*SPA(2,3))/
(complex<T>(2,0)*SPA(2,4)*SPA(3,4)*SPB(3,2))-
 (pow(SPB(4,3),2)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),2)*
 SPA(2,3))/(complex<T>(4,0)*pow(SS(2,3,4),2)*SPA(2,4)*SPA(3,4)*
 SPB(3,2))+(complex<T>(1,0)*pow(SPB(4,3),2)*SPA(2,3))/
(complex<T>(2,0)*SPA(2,5)*SPA(3,5)*SPB(3,2))-
 (pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*pow(SPB(5,3),2)*
 SPA(2,3))/(complex<T>(4,0)*pow(SS(2,3,5),2)*SPA(2,5)*SPA(3,5)*
 SPB(3,2))-pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),2)/
(SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))+
 (complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(4,2),2))/(complex<T>(2,0)*pow(SPA(3,5),3)*
 SPB(5,2))-(pow(SPA(2,3),2)*pow(SPB(3,2),2)*
 pow(SPA(2,3)*SPB(4,2)+SPA(3,5)*SPB(5,4),2))/
(complex<T>(4,0)*pow(SPA(3,5),3)*pow(SS(2,3,5),2)*SPB(5,2))-
 (pow(SPA(2,4),2)*pow(SPB(4,3),2))/(complex<T>(2,0)*SPA(2,5)*SPA(4,5)*
 (SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3)))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(4,3),2)*pow(SS(2,4,5),2))/
(complex<T>(4,0)*S(2,5)*(S(2,4)+S(4,5))*SPA(2,5)*SPA(4,5)*
 (SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3)))-
 (SPA(2,4)*SPB(4,3)*SPB(5,3))/(complex<T>(2,0)*SPA(4,5)*
 (SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3)))+
 (complex<T>(-3,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2))/
(complex<T>(2,0)*SPA(2,4)*SPA(3,4)*SS(2,3,4))-
 (pow(SPB(4,3),2)*pow(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3),2)*
 SPA(2,3))/(complex<T>(4,0)*(S(2,4)+S(3,4))*SPA(2,4)*SPA(3,4)*SPB(3,2)*
 SS(2,3,4))-(pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*
 SPB(4,3))/(S(2,3)*SPA(2,4)*SS(2,3,4))-
 (pow(SPA(2,3),2)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2))/
(complex<T>(2,0)*pow(SPA(3,5),3)*SPA(2,5)*SS(2,3,5))-
 (pow(SPA(2,3),2)*pow(SPB(3,2),2)*pow(SPA(2,3)*SPB(4,2)+
 SPA(3,5)*SPB(5,4),2)*SPA(2,5))/(complex<T>(4,0)*pow(SPA(3,5),3)*
 S(2,5)*(S(2,3)+S(3,5))*SS(2,3,5))-
 (pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*pow(SPB(5,3),2)*
 SPA(2,3))/(complex<T>(4,0)*(S(2,5)+S(3,5))*SPA(2,5)*SPA(3,5)*SPB(3,2)*
 SS(2,3,5))+(pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*
 SPB(5,3))/(S(2,3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(4,3),2)*SS(2,4,5))/
(complex<T>(4,0)*(S(2,4)+S(4,5))*SPA(2,5)*SPA(4,5)*
 (SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3)))-
 (pow(SPA(2,4),2)*pow(SPB(4,3),2)*SS(2,4,5))/
(complex<T>(4,0)*pow(SPA(2,5),2)*SPA(4,5)*SPB(5,2)*(SPA(2,4)*SPB(3,2)+
SPA(4,5)*SPB(5,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmqppp_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, qp, p, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmqppp nfLT");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(2,4)*SPB(4,3)+SPA(2,5)*SPB(5,3),2))/
(complex<T>(9,0)*SPA(2,4)*SPA(2,5)*SPA(4,5)*SPB(3,2))+
 (pow(SPB(5,4),2)*SPA(2,4)*SPB(4,3))/(complex<T>(6,0)*S(2,3)*
 (S(2,4)+S(3,4)))-(pow(SPB(5,4),2)*SPA(2,5)*SPB(5,3))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5)))-(SPB(4,3)*SPB(5,3))/
(complex<T>(3,0)*SPA(4,5)*SPB(3,2))-(pow(SPB(5,3),2)*SPA(2,5)*SPB(4,2)*
 SPB(5,4))/(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPB(3,2))-
 (pow(SPB(5,3),2)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SPB(3,2))-
 (pow(SPB(4,3),2)*SPA(2,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPB(3,2))-
 (pow(SPB(4,3),2)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SPB(3,2))+
 (SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))/(complex<T>(6,0)*S(2,3)*
 (S(2,4)+S(3,4)))+(SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*(S(2,5)+S(3,5)))+
 (complex<T>(2,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2))/
(complex<T>(9,0)*SPA(2,4)*SPA(3,4)*SS(2,3,4))+
 (complex<T>(2,0)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*SPB(4,3))/
(complex<T>(9,0)*S(2,3)*SPA(2,4)*SS(2,3,4))-
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,4)*SPB(4,3))/
(complex<T>(6,0)*S(2,3)*SS(2,3,4))+(pow(SPB(5,4),2)*SPA(2,4)*SPB(4,3))/
(complex<T>(6,0)*(S(2,4)+S(3,4))*SS(2,3,4))-
 (pow(SPB(4,3),2)*SPA(2,3)*SPA(2,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SS(2,3,4))-
 (pow(SPB(4,3),2)*SPA(2,3)*SPA(3,4)*SPB(5,3)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,4)+S(3,4),2)*SS(2,3,4))-
 (complex<T>(1,0)*SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*SS(2,3,4))+(SPA(2,3)*SPB(4,3)*SPB(5,3)*
 SPB(5,4))/(complex<T>(6,0)*(S(2,4)+S(3,4))*SS(2,3,4))+
 (complex<T>(-2,0)*pow(-(SPA(2,3)*SPB(4,3))+SPA(2,5)*SPB(5,4),2)*
 SPB(5,3))/(complex<T>(9,0)*S(2,3)*SPA(2,5)*SS(2,3,5))+
 (complex<T>(1,0)*pow(SPB(5,4),2)*SPA(2,5)*SPB(5,3))/
(complex<T>(6,0)*S(2,3)*SS(2,3,5))-(pow(SPB(5,4),2)*SPA(2,5)*SPB(5,3))/
(complex<T>(6,0)*(S(2,5)+S(3,5))*SS(2,3,5))-
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(2,5)*SPB(4,2)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SS(2,3,5))-
 (pow(SPB(5,3),2)*SPA(2,3)*SPA(3,5)*SPB(4,3)*SPB(5,4))/
(complex<T>(3,0)*pow(S(2,5)+S(3,5),2)*SS(2,3,5))-
 (complex<T>(1,0)*SPA(2,3)*SPB(4,3)*SPB(5,3)*SPB(5,4))/
(complex<T>(6,0)*S(2,3)*SS(2,3,5))+(SPA(2,3)*SPB(4,3)*SPB(5,3)*
 SPB(5,4))/(complex<T>(6,0)*(S(2,5)+S(3,5))*SS(2,3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmpqpm_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, p, qp, m}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmpqpm LT");
#endif
 
 return( (complex<T>(2,0)*pow(SPA(2,5),2))/(SPA(2,3)*SPA(3,4))+
 pow(SPB(4,3),2)/(SPB(5,2)*SPB(5,4))+
 (complex<T>(1,0)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*
 pow(SPB(5,3),2)*pow(complex<T>(1,0)-SS(2,3,5)/S(2,3),-1))/
(complex<T>(2,0)*pow(SPA(2,3),2)*pow(SPB(3,2),2)*SPB(5,2)*SPB(5,4))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*pow(complex<T>(1,0)-SS(2,4,5)/S(4,5),
-1))/(complex<T>(2,0)*pow(SPA(4,5),2)*SPB(5,2)*SPB(5,4))-
 (pow(SPB(5,3),2)*pow(SS(2,3,4),2)*SPA(2,5))/
(complex<T>(2,0)*S(2,3)*SPB(5,2)*(SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3))*
 (-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))-
 (pow(SPB(4,3),2)*S(2,5)*SPB(5,3))/(complex<T>(2,0)*SPB(3,2)*SPB(5,2)*
 SPB(5,4)*(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))-
 (pow(SPB(4,3),2)*S(3,5)*SPB(5,3))/(complex<T>(2,0)*SPB(3,2)*SPB(5,2)*
 SPB(5,4)*(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))-
 (pow(SPB(4,3),2)*S(4,5)*SPB(5,3))/(complex<T>(2,0)*SPB(3,2)*SPB(5,2)*
 SPB(5,4)*(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4)))+
 (complex<T>(1,0)*SPA(2,4)*SPA(2,5)*SS(2,3,4))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*
 (SPA(2,4)*SPB(5,2)+SPA(3,4)*SPB(5,3)))-
 pow(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3),2)/
(complex<T>(2,0)*SPA(4,5)*SPB(5,2)*SS(2,4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmpqpm_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, p, qp, m}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmpqpm nfLT");
#endif
 
   return( complex<T>(0,0));
  
} 
 
 
template <class T> complex<T> R2q2g1ph_phdqmmqpp_LT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, m, qp, p}, LT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmmqpp LT");
#endif
 
 return(-((pow(SPB(5,4),2)*SPA(2,3))/(complex<T>(2,0)*pow(SPB(4,3),2)*
SPA(3,4)))-(pow(SPB(4,2),2)*pow(-(SPA(2,3)*SPB(5,3))-
 SPA(2,4)*SPB(5,4),2)*pow(complex<T>(1,0)-SS(2,3,4)/S(3,4),-1))/
(complex<T>(2,0)*pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPB(3,2))-
 (pow(SPB(5,4),2)*S(2,4))/(complex<T>(2,0)*pow(SPB(4,3),2)*SPA(3,4)*
 SPB(3,2))+(complex<T>(19,0)*pow(SPB(5,4),2))/(complex<T>(9,0)*SPB(3,2)*
 SPB(4,3))+(complex<T>(1,0)*pow(SPA(2,3),2)*
 pow(complex<T>(1,0)-SS(2,3,5)/S(2,5),-1)*SPB(4,3))/
(complex<T>(2,0)*pow(SPA(2,5),2)*SPB(3,2))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,2))/
(complex<T>(2,0)*SPA(2,5)*SPB(3,2)*(SPA(2,4)*SPB(3,2)+
SPA(4,5)*SPB(5,3)))+(complex<T>(1,0)*SPA(2,3)*SPB(5,2)*SPB(5,3))/
(complex<T>(2,0)*SPB(3,2)*(SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3)))-
 (SPA(2,3)*SPB(5,2)*SPB(5,4))/(S(3,4)*SPB(3,2))-
 pow(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4),2)/
(complex<T>(2,0)*SPA(3,4)*SPB(3,2)*SS(2,3,4))+
 (complex<T>(1,0)*SPA(2,3)*SPA(2,4)*SS(2,4,5))/(complex<T>(2,0)*SPA(2,5)*SPA(4,5)*
 (SPA(2,4)*SPB(3,2)+SPA(4,5)*SPB(5,3)))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2g1ph_phdqmmqpp_nfLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{phd, qm, m, qp, p}, nfLT}
 
#if _VERBOSE
  _MESSAGE("R2q2g1ph :  phdqmmqpp nfLT");
#endif
 
 return( (complex<T>(2,0)*pow(SPB(5,4),2))/(complex<T>(9,0)*SPB(3,2)*SPB(4,3))
        ); 
  
} 
  
  
 
 // *************** table of switch values ************* 
 
#define _R_phqmqppm_LT R2q2g1ph_598800_LT
#define _R_phqmqppm_RT R2q2g1ph_598800_RT
#define _R_phqmqppm_nfLT R2q2g1ph_598800_nfLT
#define _R_phqmqpmp_LT R2q2g1ph_598785_LT
#define _R_phqmqpmp_RT R2q2g1ph_598785_RT
#define _R_phqmqpmp_nfLT R2q2g1ph_598785_nfLT
#define _R_phqmqpmm_LT R2q2g1ph_598784_LT
#define _R_phqmqpmm_RT R2q2g1ph_598784_RT
#define _R_phqmqpmm_nfLT R2q2g1ph_598784_nfLT
#define _R_phqmqppp_LT R2q2g1ph_598801_LT
#define _R_phqmqppp_RT R2q2g1ph_598801_RT
#define _R_phqmqppp_nfLT R2q2g1ph_598801_nfLT
#define _R_phqmpqpm_LT R2q2g1ph_598320_LT
#define _R_phqmpqpm_nfLT R2q2g1ph_598320_nfLT
#define _R_phqmmqpp_LT R2q2g1ph_598065_LT
#define _R_phqmmqpp_nfLT R2q2g1ph_598065_nfLT
#define _R_phdqmqppm_LT R2q2g1ph_664336_LT
#define _R_phdqmqppm_RT R2q2g1ph_664336_RT
#define _R_phdqmqppm_nfLT R2q2g1ph_664336_nfLT
#define _R_phdqmqpmp_LT R2q2g1ph_664321_LT
#define _R_phdqmqpmp_RT R2q2g1ph_664321_RT
#define _R_phdqmqpmp_nfLT R2q2g1ph_664321_nfLT
#define _R_phdqmqpmm_LT R2q2g1ph_664320_LT
#define _R_phdqmqpmm_RT R2q2g1ph_664320_RT
#define _R_phdqmqpmm_nfLT R2q2g1ph_664320_nfLT
#define _R_phdqmqppp_LT R2q2g1ph_664337_LT
#define _R_phdqmqppp_RT R2q2g1ph_664337_RT
#define _R_phdqmqppp_nfLT R2q2g1ph_664337_nfLT
#define _R_phdqmpqpm_LT R2q2g1ph_663856_LT
#define _R_phdqmpqpm_nfLT R2q2g1ph_663856_nfLT
#define _R_phdqmmqpp_LT R2q2g1ph_663601_LT
#define _R_phdqmmqpp_nfLT R2q2g1ph_663601_nfLT
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_phqmqppm_LT case 598800 : \
          return &R2q2g1ph_598800_LT
#define _CASE_phqmqppm_RT case 598800 : \
          return &R2q2g1ph_598800_RT
#define _CASE_phqmqppm_nfLT case 598800 : \
          return &R2q2g1ph_598800_nfLT
#define _CASE_phqmqpmp_LT case 598785 : \
          return &R2q2g1ph_598785_LT
#define _CASE_phqmqpmp_RT case 598785 : \
          return &R2q2g1ph_598785_RT
#define _CASE_phqmqpmp_nfLT case 598785 : \
          return &R2q2g1ph_598785_nfLT
#define _CASE_phqmqpmm_LT case 598784 : \
          return &R2q2g1ph_598784_LT
#define _CASE_phqmqpmm_RT case 598784 : \
          return &R2q2g1ph_598784_RT
#define _CASE_phqmqpmm_nfLT case 598784 : \
          return &R2q2g1ph_598784_nfLT
#define _CASE_phqmqppp_LT case 598801 : \
          return &R2q2g1ph_598801_LT
#define _CASE_phqmqppp_RT case 598801 : \
          return &R2q2g1ph_598801_RT
#define _CASE_phqmqppp_nfLT case 598801 : \
          return &R2q2g1ph_598801_nfLT
#define _CASE_phqmpqpm_LT case 598320 : \
          return &R2q2g1ph_598320_LT
#define _CASE_phqmpqpm_nfLT case 598320 : \
          return &R2q2g1ph_598320_nfLT
#define _CASE_phqmmqpp_LT case 598065 : \
          return &R2q2g1ph_598065_LT
#define _CASE_phqmmqpp_nfLT case 598065 : \
          return &R2q2g1ph_598065_nfLT
#define _CASE_phdqmqppm_LT case 664336 : \
          return &R2q2g1ph_664336_LT
#define _CASE_phdqmqppm_RT case 664336 : \
          return &R2q2g1ph_664336_RT
#define _CASE_phdqmqppm_nfLT case 664336 : \
          return &R2q2g1ph_664336_nfLT
#define _CASE_phdqmqpmp_LT case 664321 : \
          return &R2q2g1ph_664321_LT
#define _CASE_phdqmqpmp_RT case 664321 : \
          return &R2q2g1ph_664321_RT
#define _CASE_phdqmqpmp_nfLT case 664321 : \
          return &R2q2g1ph_664321_nfLT
#define _CASE_phdqmqpmm_LT case 664320 : \
          return &R2q2g1ph_664320_LT
#define _CASE_phdqmqpmm_RT case 664320 : \
          return &R2q2g1ph_664320_RT
#define _CASE_phdqmqpmm_nfLT case 664320 : \
          return &R2q2g1ph_664320_nfLT
#define _CASE_phdqmqppp_LT case 664337 : \
          return &R2q2g1ph_664337_LT
#define _CASE_phdqmqppp_RT case 664337 : \
          return &R2q2g1ph_664337_RT
#define _CASE_phdqmqppp_nfLT case 664337 : \
          return &R2q2g1ph_664337_nfLT
#define _CASE_phdqmpqpm_LT case 663856 : \
          return &R2q2g1ph_663856_LT
#define _CASE_phdqmpqpm_nfLT case 663856 : \
          return &R2q2g1ph_663856_nfLT
#define _CASE_phdqmmqpp_LT case 663601 : \
          return &R2q2g1ph_663601_LT
#define _CASE_phdqmmqpp_nfLT case 663601 : \
          return &R2q2g1ph_663601_nfLT
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_phqmqppm_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqppm_LT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqppm_RT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqppm_RT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqppm_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqppm_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpmp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqpmp_LT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpmp_RT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqpmp_RT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpmp_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqpmp_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpmm_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqpmm_LT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpmm_RT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqpmm_RT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqpmm_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqpmm_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqppp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqppp_LT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqppp_RT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqppp_RT(ep,mpc);}
 
template <class T> complex<T> _R_phqmqppp_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmqppp_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phqmpqpm_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmpqpm_LT(ep,mpc);}
 
template <class T> complex<T> _R_phqmpqpm_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmpqpm_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phqmmqpp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmmqpp_LT(ep,mpc);}
 
template <class T> complex<T> _R_phqmmqpp_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phqmmqpp_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqppm_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqppm_LT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqppm_RT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqppm_RT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqppm_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqppm_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqpmp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqpmp_LT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqpmp_RT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqpmp_RT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqpmp_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqpmp_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqpmm_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqpmm_LT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqpmm_RT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqpmm_RT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqpmm_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqpmm_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqppp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqppp_LT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqppp_RT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqppp_RT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmqppp_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmqppp_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmpqpm_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmpqpm_LT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmpqpm_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmpqpm_nfLT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmmqpp_LT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmmqpp_LT(ep,mpc);}
 
template <class T> complex<T> _R_phdqmmqpp_nfLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2g1ph_phdqmmqpp_nfLT(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q2g1ph_LT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phqmqppm_LT;
       _CASE_phqmqpmp_LT;
       _CASE_phqmqpmm_LT;
       _CASE_phqmqppp_LT;
       _CASE_phqmpqpm_LT;
       _CASE_phqmmqpp_LT;
       _CASE_phdqmqppm_LT;
       _CASE_phdqmqpmp_LT;
       _CASE_phdqmqpmm_LT;
       _CASE_phdqmqppp_LT;
       _CASE_phdqmpqpm_LT;
       _CASE_phdqmmqpp_LT;
 
       default: cout << "Missing rational term: " << hc << endl; return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g1ph_nfLT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phqmqppm_nfLT;
       _CASE_phqmqpmp_nfLT;
       _CASE_phqmqpmm_nfLT;
       _CASE_phqmqppp_nfLT;
       _CASE_phqmpqpm_nfLT;
       _CASE_phqmmqpp_nfLT;
       _CASE_phdqmqppm_nfLT;
       _CASE_phdqmqpmp_nfLT;
       _CASE_phdqmqpmm_nfLT;
       _CASE_phdqmqppp_nfLT;
       _CASE_phdqmpqpm_nfLT;
       _CASE_phdqmmqpp_nfLT;
 
       default: cout << "Missing rational term: " << hc << endl; return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2g1ph_RT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_phqmqppm_RT;
       _CASE_phqmqpmp_RT;
       _CASE_phqmqpmm_RT;
       _CASE_phqmqppp_RT;
       _CASE_phdqmqppm_RT;
       _CASE_phdqmqpmp_RT;
       _CASE_phdqmqpmm_RT;
       _CASE_phdqmqppp_RT;
 
       default: cout << "Missing rational term: " << hc << endl; return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q2g1ph_LT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2g1ph_LT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2g1ph_LT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2g1ph_nfLT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2g1ph_nfLT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2g1ph_nfLT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2g1ph_RT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2g1ph_RT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2g1ph_RT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> ( *R2q2g1ph_LT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2g1ph_nfLT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2g1ph_RT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
#endif


}
