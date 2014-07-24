/*
* R_2q2Q1g_eval.cpp
*
* Created on 12/22, 2010
*      Author: Zvi's script
*/
 
#include "R_2q2Q1g_eval.h"
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


 
 
template <class T> complex<T> R2q2Q1g_qmqpmQmQp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpmQmQp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(5,2),2)*SPA(1,5)*SPA(2,3))/
(complex<T>(2,0)*S(1,2)*S(4,5)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(3,2),2)*pow(SPB(5,1),2)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPB(2,1)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(3,2),2)*pow(SPB(5,1),2)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),2)*
 pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPB(2,1)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(3,2),2)*pow(SPB(5,1),2)*S(1,2)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),3)*pow(SPB(5,4),4)*SPB(2,1)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(3,2),2)*pow(SPB(5,1),2)*S(2,3)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),2)*
 pow(SPA(4,5),3)*pow(SPB(5,4),4)*SPB(2,1)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,1)*SPA(1,4)*SPA(3,4)*SPB(4,2))/(complex<T>(2,0)*S(2,3)*SPA(4,5)*
 SPB(4,3))+(complex<T>(0,1)*pow(SPA(1,4),2)*SPB(4,2))/
(complex<T>(2,0)*SPA(1,2)*SPA(4,5)*SPB(3,2)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(3,2)*SPB(4,2))/
(complex<T>(2,0)*pow(SPB(2,1),2)*pow(SPB(5,4),2)*SPA(1,2)*SPA(4,5)*
 SPB(4,3))+(complex<T>(0,1)*SPA(1,3)*SPA(1,4)*SPB(5,2))/
(complex<T>(3,0)*S(4,5)*SPA(1,2)*SPB(3,2))+
 (complex<T>(0,-1)*SPA(1,3)*SPA(3,4)*SPB(3,2)*SPB(4,1)*SPB(5,2))/
(pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*pow(SPB(5,4),2)*
 SPB(2,1)*SPB(4,3))+(complex<T>(0,1)*S(1,2)*SPA(1,3)*SPA(3,4)*SPB(3,2)*
 SPB(4,1)*SPB(5,2))/(pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),3)*pow(SPB(5,4),3)*SPB(2,1)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*S(4,5)*SPA(3,4)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3))+(complex<T>(0,1)*pow(SPA(1,3),2)*SPA(3,4)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),2)*S(4,5)*
 SPB(2,1))+(complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(1,3))/
(complex<T>(2,0)*S(1,2)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(5,1),2)*S(4,5)*SPB(4,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),2)*SPB(2,1)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(5,1),2)*SPB(4,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),2)*SPB(2,1)*
 SPB(3,2)*SPB(4,3)*SPB(5,4))+(complex<T>(0,-7)*pow(SPB(5,2),2)*
 SPB(4,2))/(complex<T>(9,0)*SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPB(3,2)*SPB(4,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*S(4,5)*SPB(3,2)*
 SPB(4,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),4)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*SPB(5,1)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPB(5,4))+(complex<T>(0,1)*pow(SPA(1,3),2)*S(4,5)*
 SPB(5,1)*SPB(5,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),3)*SPB(5,4))+
 (complex<T>(0,1)*SPA(1,3)*SPB(4,2)*SPB(5,1)*SPB(5,2))/
(complex<T>(2,0)*S(2,3)*SPB(2,1)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*S(4,5)*SPB(5,2)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*SPB(5,2)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*SPB(2,1)*
 SPB(3,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqppQmQp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqppQmQp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,3),2)*SPA(1,4)*SPA(1,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),2)*SPA(1,2))+(complex<T>(0,1)*pow(SPB(5,3),2)*S(1,2)*
 SPA(1,4)*SPA(1,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),3)*pow(SPB(5,4),3)*SPA(1,2))+
 (complex<T>(0,1)*pow(SPB(5,3),2)*S(1,2)*SPA(1,3)*SPA(1,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(1,2)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(4,3),2)*SPA(2,4)*SPA(3,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),3)*
 pow(SPB(5,4),2)*SPA(1,2)*SPA(2,3))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(4,3),2)*S(1,2)*SPA(2,4)*
 SPA(3,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),4)*pow(SPB(5,4),3)*SPA(1,2)*SPA(2,3))+
 (complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPB(5,3),2)*S(1,2)*SPA(2,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*
 pow(SPB(4,3),3)*SPA(1,2)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPB(5,3),2)*SPA(1,3)*SPA(1,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPB(5,4),2)*SPA(1,2)*
 SPA(3,4)*SPA(4,5))+(complex<T>(0,-7)*pow(SPA(1,4),2)*SPA(2,4))/
(complex<T>(9,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPB(5,3),2)*SPA(2,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPB(4,3),2)*SPA(1,2)*
 SPA(2,3)*SPA(3,4)*SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,5),2)*
 pow(SPA(3,4),2)*pow(SPB(5,3),2)*SPA(2,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPA(2,3)*SPA(3,5)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*
 SPA(2,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPA(2,3)*SPA(3,5)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*S(3,4)*
 SPA(2,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3)*SPA(2,3)*SPA(3,5)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),2)*S(4,5)*
 SPA(2,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3)*SPA(2,3)*SPA(3,5)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPB(5,3),2)*S(1,2)*SPA(1,3)*SPB(3,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),3)*
 pow(SPB(5,4),4))+(complex<T>(0,1)*pow(SPB(5,3),2)*SPA(1,3)*SPB(3,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPB(5,4),2)*S(1,2)*
 SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,4),2)*SPB(4,3)*SPB(5,1))/
(complex<T>(2,0)*S(1,2)*S(4,5)*SPA(2,3))+
 (complex<T>(0,1)*SPA(2,4)*SPB(3,2)*SPB(5,2))/(complex<T>(2,0)*S(3,4)*SPA(2,3)*
 SPB(2,1))+(complex<T>(0,-1)*pow(SPA(1,4),2)*SPB(5,3))/
(complex<T>(2,0)*S(4,5)*SPA(1,2)*SPA(3,4))+
 (complex<T>(0,1)*SPA(1,4)*SPA(1,5)*SPA(2,4)*SPB(5,3))/
(complex<T>(2,0)*S(3,4)*SPA(1,2)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,-1)*SPA(1,4)*SPA(2,5)*SPA(3,4)*SPB(3,2)*SPB(5,3))/
(pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*
 SPA(2,3)*SPA(4,5))+(complex<T>(0,1)*S(4,5)*SPA(1,4)*SPA(2,5)*SPA(3,4)*
 SPB(3,2)*SPB(5,3))/(pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),3)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPB(5,2),2)*SPA(2,4))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*
 SPB(2,1)*SPB(5,4))+(complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(4,3),2)*
 SPA(2,4)*SPA(3,4))/(complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPA(4,5),2)*
 SPA(2,3)*SPB(2,1)*SPB(5,4))+(complex<T>(0,1)*SPA(1,4)*SPB(5,2)*
 SPB(5,3))/(complex<T>(3,0)*S(1,2)*SPA(3,4)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpmQpQm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpmQpQm LLT");
#endif
 
 return( (complex<T>(0,5)*SPA(1,3)*SPA(3,5)*SPB(4,2))/(complex<T>(6,0)*S(1,2)*S(4,5))+
 (complex<T>(0,1)*pow(SPB(4,2),2)*S(1,3)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),2)*SPB(2,1)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPB(4,2),2)*S(1,2)*S(1,3)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),3)*
 pow(SPB(5,4),3)*SPB(2,1)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(3,2),2)*pow(SPB(4,1),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPB(2,1)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(3,2),2)*pow(SPB(4,1),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPB(2,1)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(3,2),2)*pow(SPB(4,1),3)*S(1,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),3)*
 pow(SPB(5,4),4)*SPB(2,1)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(3,2),2)*pow(SPB(4,1),3)*S(2,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),2)*pow(SPA(4,5),3)*
 pow(SPB(5,4),4)*SPB(2,1)*SPB(3,1)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPB(4,2),2)*SPA(1,3)*SPB(4,1))/
(complex<T>(2,0)*pow(SPB(5,4),2)*SPA(4,5)*SPB(2,1)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(4,1),2)*SPB(3,2)*SPB(4,2))/
(complex<T>(2,0)*pow(SPB(2,1),2)*pow(SPB(5,4),2)*SPA(1,2)*SPA(4,5)*
 SPB(4,3))+(complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(4,1),2)*SPB(3,2)*
 SPB(4,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),2)*
 pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPB(2,1)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(4,1),2)*S(2,3)*SPB(3,2)*
 SPB(4,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),2)*
 pow(SPA(4,5),3)*pow(SPB(5,4),4)*SPB(2,1)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*S(4,5)*SPA(3,5)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPA(3,5)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),2)*S(4,5)*
 SPB(2,1))+(complex<T>(0,1)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(4,2)*
 SPB(4,3))/(complex<T>(2,0)*pow(SPB(2,1),2)*pow(SPB(5,4),2)*SPA(1,2)*
 SPA(4,5)*SPB(3,2))+(complex<T>(0,1)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*
 SPB(4,2)*SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPB(2,1)*SPB(3,2))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*S(1,2)*SPB(4,2)*
 SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),3)*pow(SPB(5,4),4)*SPB(2,1)*SPB(3,2))+
 (complex<T>(0,1)*pow(SPB(4,2),2)*S(3,5)*SPA(1,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPB(4,2),2)*S(3,5)*S(4,5)*SPA(1,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),3)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*pow(SPB(5,2),2)*SPA(3,4)*
 SPB(4,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),4)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,-7)*pow(SPB(4,2),3))/(complex<T>(9,0)*SPB(2,1)*SPB(3,2)*SPB(4,3)*
 SPB(5,4))+(complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(4,1),2)*SPB(3,2)*
 SPB(4,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(4,1),2)*S(4,5)*SPB(3,2)*
 SPB(4,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),4)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(5,2),2)*SPB(4,2)*SPB(4,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPB(4,2),2)*SPA(3,5)*SPB(5,2))/
(complex<T>(2,0)*pow(SPB(2,1),2)*SPA(1,2)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*pow(SPB(5,2),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPB(3,2)*SPB(5,3)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*pow(SPB(5,2),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPB(3,2)*SPB(5,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(4,3),2)*pow(SPB(5,2),3)*S(4,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4)*SPB(3,2)*SPB(5,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(4,3),3)*pow(SPB(5,2),3)*
 SPA(3,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),4)*SPB(3,2)*SPB(5,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*SPA(1,3)*SPB(3,2)*SPB(5,4))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*S(4,5))+(complex<T>(0,1)*pow(SPA(3,5),2)*S(4,5)*
 SPA(1,3)*SPB(3,2)*SPB(5,4))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),
3)*pow(SPA(1,2),4)*pow(SPB(2,1),4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqppQpQm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqppQpQm LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(3,2),2)*S(4,5)*SPA(3,5))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(3,4))+(complex<T>(0,1)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(3,4))+(complex<T>(0,-23)*pow(SPA(1,5),2)*SPA(2,4))/
(complex<T>(18,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*S(3,4)*S(4,5)*SPA(3,5)*SPB(3,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,-1)*S(3,4)*SPA(3,5)*SPB(3,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(2,3)*SPA(3,4))+(complex<T>(0,1)*SPA(1,5)*SPA(2,5)*SPB(3,2))/
(complex<T>(2,0)*S(1,2)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,1)*SPA(1,4)*SPA(1,5)*SPB(4,3))/(complex<T>(2,0)*S(4,5)*SPA(1,2)*
 SPA(3,4))+(complex<T>(0,-1)*S(1,2)*S(3,4)*SPA(1,3)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,1)*S(3,4)*SPA(1,3)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPB(5,4),2)*SPA(2,3)*
 SPA(3,4)*SPA(4,5))+(complex<T>(0,1)*SPA(1,3)*SPA(2,5)*SPB(3,2)*
 SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPA(2,3))+
 (complex<T>(0,-1)*S(4,5)*SPA(1,3)*SPA(2,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),3)*SPA(2,3))+(complex<T>(0,1)*S(1,2)*SPA(1,3)*SPB(3,2)*
 SPB(4,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPA(3,4))+
 (complex<T>(0,1)*SPA(1,4)*SPA(3,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),2)*SPA(3,4))+(complex<T>(0,-1)*S(1,2)*SPA(1,4)*SPA(3,5)*
 SPB(3,2)*SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),3)*pow(SPB(5,4),3)*SPA(3,4))+
 (complex<T>(0,1)*SPA(1,3)*SPA(2,4)*SPA(3,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,-1)*pow(SPA(4,5),2)*pow(SPB(5,4),2)*SPA(1,3)*SPA(2,4)*
 SPA(3,5)*SPB(3,2)*SPB(4,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),
3)*pow(SPA(1,2),4)*pow(SPB(2,1),4)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,1)*SPA(1,3)*SPA(2,4)*SPA(3,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*S(1,2)*S(4,5)*SPA(2,3)*
 SPA(3,4))+(complex<T>(0,-1)*S(4,5)*SPA(1,3)*SPA(2,4)*SPA(3,5)*SPB(3,2)*
 SPB(4,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*
 pow(SPA(1,2),3)*pow(SPB(2,1),3)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,-1)*SPA(1,3)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPB(5,4),2)*SPA(3,4)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*SPB(2,1))/
(complex<T>(6,0)*pow(SPA(4,5),2)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(5,1))+
 (complex<T>(0,1)*pow(SPA(2,4),2)*SPA(1,5)*SPB(4,2))/
(complex<T>(6,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(1,5),2)*SPA(2,4)*SPB(3,1)*SPB(5,3))/
(complex<T>(6,0)*pow(SPA(1,2),2)*pow(SPA(4,5),2)*pow(SPB(5,1),2))+
 (complex<T>(0,-1)*SPA(1,5)*SPB(3,1)*SPB(4,2)*SPB(5,3))/
(complex<T>(6,0)*S(1,2)*S(4,5)*SPB(5,1))+
 (complex<T>(0,1)*pow(SPB(4,2),2)*SPA(2,4))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*
 SPB(2,1)*SPB(5,4))+(complex<T>(0,-1)*SPA(1,3)*SPA(3,5)*SPB(3,2)*
 SPB(4,3)*SPB(5,1))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*
 pow(SPB(2,1),2)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,4),3)*SPB(2,1)*SPB(5,4))/
(complex<T>(6,0)*pow(SPA(1,2),2)*pow(SPA(4,5),2)*pow(SPB(5,1),2)*
 SPA(2,3)*SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*
 SPB(5,4))/(complex<T>(6,0)*pow(SPA(1,2),2)*SPA(2,3)*SPA(3,4)*SPA(4,5)*
 SPB(5,1))+(complex<T>(0,1)*pow(SPA(4,5),2)*SPA(1,3)*SPA(3,5)*SPB(3,2)*
 SPB(4,3)*SPB(5,1)*SPB(5,4))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),
3)*pow(SPA(1,2),3)*pow(SPB(2,1),4)*SPA(2,3)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmmQp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, m, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmmQp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,4),3)*pow(SPB(4,2),2)*pow(SPB(5,1),2)*
 SPA(4,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),4)*pow(SPB(3,2),5)*(complex<T>(1,0)-S(1,5)/S(2,3)-
S(4,5)/S(2,3)))+(complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(5,2),2)*
 SPA(4,5)*SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),4))+
 (complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(1,2)*SPA(4,5))/
(complex<T>(2,0)*S(1,5)*S(2,3)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(2,1),2)*pow(SPB(5,4),2)*SPA(1,5)*
 SPB(3,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),4)*SPB(4,1)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*S(2,3)*SPB(4,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),3)*
 pow(SPB(5,1),2)*SPB(4,3))+(complex<T>(0,-1)*pow(SPA(1,4),2)*
 pow(SPB(2,1),2)*pow(SPB(5,4),2)*SPB(3,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),2)*
 pow(SPB(3,2),3)*SPB(4,1)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(2,1),2)*pow(SPB(5,4),2)*
 SPB(3,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPB(4,1)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(2,1),2)*pow(SPB(5,4),3)*SPA(4,5)*
 SPB(3,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),4)*SPB(4,1)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*SPB(4,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),2)*SPB(4,3)*
 SPB(5,1))+(complex<T>(0,-1)*pow(SPA(1,4),2)*SPA(2,3)*SPB(2,1)*
 SPB(5,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*
 pow(SPA(1,5),3)*pow(SPB(5,1),3))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*SPB(2,1)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),2)*
 pow(SPB(5,1),2)*SPB(3,2))+(complex<T>(0,-1)*SPA(1,4)*SPB(5,2))/
(S(1,2)*SPB(4,3))+(complex<T>(0,1)*SPA(1,3)*SPA(1,4)*SPA(4,5)*SPB(5,1)*
 SPB(5,2))/(pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),3))+(complex<T>(0,-1)*SPA(1,3)*SPA(3,4)*SPB(5,3))/
(complex<T>(2,0)*S(4,5)*SPA(2,3)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(2,1),2)*SPA(2,3)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),2)*pow(SPA(4,5),3)*
 pow(SPB(5,4),2)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,1)*SPA(1,3)*SPA(4,5)*SPB(5,2)*SPB(5,3))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(4,5))*SPB(4,3))+
 (complex<T>(0,-1)*SPA(1,4)*SPB(2,1)*SPB(5,2)*SPB(5,3))/
(complex<T>(2,0)*S(4,5)*SPB(3,2)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(1,4),3)*pow(SPB(4,2),2)*pow(SPB(5,1),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-S(4,5)/S(2,3))*
 SPB(5,4))+(complex<T>(0,-1)*pow(SPA(1,4),3)*pow(SPB(4,2),2)*
 pow(SPB(5,1),2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-
S(4,5)/S(2,3))*SPB(5,4))+(complex<T>(0,1)*pow(SPA(1,4),3)*
 pow(SPB(4,2),2)*pow(SPB(5,1),3)*SPA(1,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),4)*
 pow(SPB(3,2),5)*(complex<T>(1,0)-S(1,5)/S(2,3)-S(4,5)/S(2,3))*
 SPB(5,4))+(complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(1,4))/
(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,3)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPB(5,2),2)*SPA(1,4))/(complex<T>(2,0)*S(1,5)*SPB(3,2)*
 SPB(5,4))+(complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(5,2),2)*SPB(4,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),2)*
 pow(SPB(3,2),3)*SPB(5,4))+(complex<T>(0,-7)*pow(SPB(5,2),2))/
(complex<T>(9,0)*SPB(2,1)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,-1)*SPA(1,3)*SPA(1,4)*SPB(5,1)*SPB(5,2))/
(pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),2)*pow(SPB(3,2),2)*
 SPB(5,4))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPB(5,3))/
(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(2,1),2)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(4,5),2)*pow(SPA(4,5),2)*SPB(3,2)*
 SPB(4,3)*SPB(5,1)*SPB(5,4))+(complex<T>(0,-1)*SPA(1,4)*SPA(3,4)*
 SPB(3,2)*SPB(5,4))/(complex<T>(2,0)*S(1,2)*(S(1,2)-S(4,5))*SPB(4,3))+
 (complex<T>(0,-1)*SPA(1,4)*SPA(1,5)*SPA(3,4)*SPB(3,1)*SPB(5,2)*SPB(5,4))/
(pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),3)*pow(SPB(3,2),3)*
 SPB(4,3))+(complex<T>(0,1)*SPA(1,4)*SPA(3,4)*SPB(3,1)*SPB(5,2)*
 SPB(5,4))/(pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),2)*
 pow(SPB(3,2),2)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3)*SPB(5,4))/
(complex<T>(2,0)*pow(SPB(3,2),2)*pow(SPB(5,1),2)*SPA(1,5)*SPA(2,3)*
 SPB(4,3))+(complex<T>(0,1)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(2,3)*
 SPB(5,3)*SPB(5,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*
 pow(SPA(1,5),3)*pow(SPB(5,1),4)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPB(5,3)*SPB(5,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),2)*
 pow(SPB(5,1),3)*SPB(3,2)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmpQp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, p, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmpQp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(4,2),2)*SPA(1,2)*SPA(1,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),2)*
 pow(SPB(3,2),2)*SPA(1,5))+(complex<T>(0,-1)*pow(SPB(4,2),2)*S(1,5)*
 SPA(1,2)*SPA(1,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),3)*SPA(1,5))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPA(2,3),2)*pow(SPB(4,2),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),4)*
 pow(SPB(5,1),3)*(complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5))*
 SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPA(2,3),2)*
 pow(SPB(4,2),3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),4)*pow(SPB(5,1),3)*(complex<T>(1,0)-S(2,3)/S(1,5)-
S(3,4)/S(1,5))*SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,4),2)*
 pow(SPA(2,3),2)*pow(SPB(4,2),3)*S(2,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),5)*
 pow(SPB(5,1),4)*(complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5))*
 SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(4,2),2)*SPA(2,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),3)*
 pow(SPB(5,1),2)*SPA(3,4))+(complex<T>(0,1)*pow(SPB(4,2),2)*S(1,5)*
 SPA(1,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*
 pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPB(4,2),2)*SPA(1,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPB(3,2),2)*SPA(2,3)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,2),2)*pow(SPA(3,4),2)*
 pow(SPB(4,2),2)*SPA(2,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),
2)*pow(SPA(1,5),3)*pow(SPB(5,1),2)*SPA(2,3)*SPA(2,4)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,2),2)*pow(SPA(3,4),2)*
 pow(SPB(4,2),2)*SPA(2,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),
2)*pow(SPA(1,5),3)*pow(SPB(5,1),2)*SPA(2,3)*SPA(2,4)*
 SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,2),2)*pow(SPA(3,4),2)*
 pow(SPB(4,2),2)*S(2,3)*SPA(2,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),4)*
 pow(SPB(5,1),3)*SPA(2,3)*SPA(2,4)*SPA(4,5))+
 (complex<T>(0,-7)*pow(SPA(1,3),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,2),2)*pow(SPB(4,2),2)*S(1,5)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),2)*
 pow(SPB(4,3),3)*SPA(1,5)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,2),2)*pow(SPB(4,2),2)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPB(4,3),2)*SPA(1,5)*
 SPA(2,3)*SPA(3,4)*SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,3),2)*
 pow(SPB(4,3),2)*SPA(3,4)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),2)*SPA(1,5)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(4,3),2)*S(1,5)*SPA(3,4)*
 SPA(3,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*
 pow(SPA(2,3),4)*pow(SPB(3,2),3)*SPA(1,5)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*SPB(4,2))/(complex<T>(2,0)*S(2,3)*SPA(1,5)*
 SPA(3,4))+(complex<T>(0,-1)*SPA(1,3)*SPB(4,2))/(S(1,2)*SPA(4,5))+
 (complex<T>(0,-1)*SPA(1,2)*SPA(1,3)*SPA(3,5)*SPB(4,2))/
(complex<T>(2,0)*S(3,4)*SPA(1,5)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPA(2,3),2)*pow(SPB(4,2),3)*
 SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),5)*pow(SPB(5,1),4)*(complex<T>(1,0)-S(2,3)/S(1,5)-
S(3,4)/S(1,5)))+(complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(4,2),2)*
 SPA(2,4)*SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),4)*pow(SPB(5,1),3))+
 (complex<T>(0,1)*pow(SPA(1,2),2)*pow(SPA(3,4),3)*pow(SPB(4,2),2)*SPA(2,5)*
 SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),4)*pow(SPB(5,1),3)*SPA(2,3)*SPA(2,4)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*SPB(2,1)*SPB(4,3))/
(complex<T>(2,0)*S(1,5)*S(2,3)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(3,5))/(complex<T>(2,0)*SPA(3,4)*SPA(4,5)*
 SPB(3,2)*SPB(5,1))+(complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(4,3),2)*
 SPA(3,4)*SPA(3,5))/(complex<T>(2,0)*pow(SPA(1,5),2)*pow(SPA(2,3),2)*
 SPA(4,5)*SPB(3,2)*SPB(5,1))+(complex<T>(0,-1)*pow(SPA(1,3),2)*
 SPB(4,2))/(complex<T>(2,0)*pow(SPA(1,5),2)*SPA(3,4)*SPB(5,1))+
 (complex<T>(0,-1)*SPA(1,3)*SPA(2,3)*SPB(4,2)*SPB(5,2))/
(pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),2)*pow(SPB(5,1),2)*
 SPA(3,4))+(complex<T>(0,1)*SPA(1,3)*SPA(3,5)*SPB(4,3)*SPB(5,2))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(3,4))*SPA(4,5))+
 (complex<T>(0,1)*SPA(1,3)*SPA(2,3)*SPB(4,2)*SPB(4,3)*SPB(5,2))/
(pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),3)*
 pow(SPB(5,1),3))+(complex<T>(0,-1)*SPA(1,5)*SPA(3,4)*SPB(4,2)*
 SPB(5,4))/(complex<T>(2,0)*S(1,2)*(S(1,2)-S(3,4))*SPA(4,5))+
 (complex<T>(0,1)*SPA(1,3)*SPA(2,5)*SPA(3,4)*SPB(4,2)*SPB(5,4))/
(pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),2)*pow(SPB(5,1),2)*
 SPA(2,3)*SPA(4,5))+(complex<T>(0,-1)*S(2,3)*SPA(1,3)*SPA(2,5)*SPA(3,4)*
 SPB(4,2)*SPB(5,4))/(pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*
 pow(SPA(1,5),3)*pow(SPB(5,1),3)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,-1)*SPA(3,5)*SPB(5,2)*SPB(5,4))/(complex<T>(2,0)*S(3,4)*SPA(4,5)*
 SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpmQm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, m, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpmQm LLT");
#endif
 
 return( (complex<T>(0,1)*SPA(1,4)*SPB(3,2))/(S(1,2)*SPB(5,4))+
 (complex<T>(0,-7)*pow(SPB(3,2),2))/(complex<T>(9,0)*SPB(2,1)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,1)*S(3,4)*SPA(1,4)*SPA(4,5)*SPB(5,2))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(3,4))*SPA(3,4)*SPB(5,4))+
 (complex<T>(0,-1)*S(3,4)*SPA(1,5)*SPB(3,2)*SPB(5,3))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(3,4))*SPB(4,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQppQm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, p, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQppQm LLT");
#endif
 
 return( (complex<T>(0,-7)*pow(SPA(1,5),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4)*
 SPA(4,5))+(complex<T>(0,-1)*S(4,5)*SPA(1,5)*SPA(3,5)*SPB(3,2))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(4,5))*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*SPA(1,5)*SPB(4,2))/(S(1,2)*SPA(3,4))+
 (complex<T>(0,1)*S(4,5)*SPA(1,3)*SPB(4,2)*SPB(4,3))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(4,5))*SPA(3,4)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmQpm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, m}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmQpm LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(3,5),2)*SPB(3,2)*SPB(4,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPB(2,1))+(complex<T>(0,1)*pow(SPA(3,5),2)*S(1,2)*
 SPB(3,2)*SPB(4,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),3)*SPB(2,1))+
 (complex<T>(0,1)*pow(SPB(4,2),2)*SPA(2,3)*SPA(4,5))/(complex<T>(2,0)*S(1,2)*S(3,4)*
 SPB(5,1))+(complex<T>(0,1)*SPA(1,3)*SPA(1,5)*SPB(4,1))/
(complex<T>(2,0)*S(4,5)*SPA(1,2)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*S(1,2)*SPB(4,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),3)*
 pow(SPB(5,4),2)*SPB(2,1)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,1)*SPA(3,5)*SPB(3,2)*SPB(4,1)*SPB(4,2))/
(complex<T>(2,0)*S(4,5)*SPB(2,1)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*S(1,2)*SPA(1,5)*SPB(5,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3))+(complex<T>(0,1)*pow(SPA(3,5),2)*SPA(1,5)*SPB(5,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*S(1,2)*
 SPB(4,3))+(complex<T>(0,1)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*
 pow(SPB(5,4),2)*SPB(3,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),
2)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(4,3)*SPB(5,1)*
 SPB(5,3))+(complex<T>(0,1)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*
 pow(SPB(5,4),2)*SPB(3,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),
2)*pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(4,3)*SPB(5,1)*
 SPB(5,3))+(complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*
 pow(SPB(5,4),2)*S(3,4)*SPB(3,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4)*SPB(4,3)*SPB(5,1)*SPB(5,3))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*pow(SPB(5,4),2)*S(4,5)*
 SPB(3,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),4)*SPB(4,3)*SPB(5,1)*SPB(5,3))+
 (complex<T>(0,-1)*pow(SPB(4,2),2)*SPA(3,5))/(complex<T>(2,0)*S(3,4)*SPB(2,1)*
 SPB(5,4))+(complex<T>(0,1)*SPA(1,3)*SPA(3,5)*SPB(4,2))/
(complex<T>(3,0)*S(1,2)*SPA(3,4)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*SPB(4,1))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*pow(SPA(3,5),2)*pow(SPB(3,2),2)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),2)*SPB(2,1)*SPB(4,3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-7)*pow(SPB(4,2),2)*SPB(4,1))/(complex<T>(9,0)*SPB(2,1)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*pow(SPA(3,5),2)*S(1,2)*SPB(4,2)*
 SPB(5,2))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPB(2,1)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*SPB(4,2)*SPB(5,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*SPB(2,1)*
 SPB(4,3)*SPB(5,4))+(complex<T>(0,1)*pow(SPA(4,5),2)*pow(SPB(4,2),2)*
 SPB(4,1)*SPB(5,4))/(complex<T>(2,0)*pow(SPB(2,1),2)*pow(SPB(4,3),2)*
 SPA(1,2)*SPA(3,4)*SPB(5,1))+(complex<T>(0,1)*pow(SPA(4,5),2)*
 pow(SPB(4,2),2)*SPB(4,1)*SPB(5,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*
 pow(SPB(4,3),3)*SPB(2,1)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(4,5),2)*pow(SPB(4,2),2)*S(1,2)*SPB(4,1)*
 SPB(5,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),4)*SPB(2,1)*SPB(5,1))+
 (complex<T>(0,-1)*SPA(1,5)*SPA(3,5)*SPB(3,1)*SPB(4,2)*SPB(5,4))/
(pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*
 SPB(4,3)*SPB(5,1))+(complex<T>(0,1)*S(3,4)*SPA(1,5)*SPA(3,5)*SPB(3,1)*
 SPB(4,2)*SPB(5,4))/(pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),3)*SPB(4,3)*SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmQpp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, p}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmQpp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(1,3)*SPA(2,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPA(3,4))+(complex<T>(0,1)*pow(SPB(5,2),2)*S(3,4)*
 SPA(1,3)*SPA(2,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),3)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPB(5,2),2)*S(3,4)*SPA(1,3)*SPA(3,5))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(1,5)*SPA(3,4))+
 (complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(1,3)*SPA(3,5))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(1,5)*SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,3),2)*
 pow(SPB(5,2),2)*SPA(2,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),
2)*pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(1,2)*SPA(2,5)*
 SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,3),2)*
 pow(SPB(5,2),2)*SPA(2,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),
2)*pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(1,2)*SPA(2,5)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,3),2)*
 pow(SPB(5,2),2)*S(1,2)*SPA(2,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3)*SPA(1,2)*SPA(2,5)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*S(3,4)*SPA(1,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),2)*
 pow(SPB(5,1),3)*SPA(1,2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-7)*pow(SPA(1,3),2)*SPA(1,4))/(complex<T>(9,0)*SPA(1,2)*SPA(1,5)*
 SPA(3,4)*SPA(4,5))+(complex<T>(0,1)*pow(SPA(2,3),2)*pow(SPB(5,2),2)*
 SPA(1,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPB(5,1),2)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(5,1),2)*SPA(1,4)*SPA(1,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(5,1),2)*S(3,4)*SPA(1,4)*
 SPA(1,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPB(4,2),2)*SPA(1,4))/(complex<T>(2,0)*SPA(1,5)*SPA(4,5)*
 SPB(2,1)*SPB(4,3))+(complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(5,1),2)*
 SPA(1,4)*SPA(1,5))/(complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPA(3,4),2)*
 SPA(4,5)*SPB(2,1)*SPB(4,3))+(complex<T>(0,-1)*pow(SPA(1,5),3)*
 pow(SPA(2,3),2)*pow(SPB(5,2),2)*SPA(2,4)*SPB(5,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3)*SPA(1,2)*SPA(2,5)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*SPB(3,2)*SPB(5,1))/(complex<T>(2,0)*S(1,2)*S(3,4)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPB(5,2))/
(complex<T>(2,0)*S(1,2)*SPA(1,5)*SPA(3,4))+
 (complex<T>(0,1)*SPA(1,3)*SPA(1,4)*SPA(2,3)*SPB(5,2))/
(complex<T>(2,0)*S(1,5)*SPA(1,2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*SPA(1,3)*SPB(4,2)*SPB(5,2))/(complex<T>(3,0)*S(3,4)*SPA(1,5)*
 SPB(2,1))+(complex<T>(0,-1)*pow(SPB(5,2),2)*S(3,4)*SPA(3,5)*SPB(5,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4))+(complex<T>(0,1)*pow(SPB(5,2),2)*SPA(3,5)*SPB(5,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPB(2,1),2)*S(3,4)*
 SPA(1,2))+(complex<T>(0,1)*SPA(1,4)*SPB(4,2)*SPB(5,4))/
(complex<T>(2,0)*S(1,5)*SPA(4,5)*SPB(4,3))+
 (complex<T>(0,-1)*SPA(1,3)*SPA(1,5)*SPA(2,4)*SPB(5,2)*SPB(5,4))/
(pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*pow(SPB(4,3),2)*
 SPA(1,2)*SPA(4,5))+(complex<T>(0,1)*S(1,2)*SPA(1,3)*SPA(1,5)*SPA(2,4)*
 SPB(5,2)*SPB(5,4))/(pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),3)*SPA(1,2)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpQmm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, m}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpQmm LLT");
#endif
 
 return( (complex<T>(0,-1)*SPA(1,4)*SPA(2,5)*SPA(3,5)*SPB(3,2))/
(complex<T>(6,0)*S(1,2)*S(3,4)*SPA(2,3))+
 (complex<T>(0,-1)*pow(SPB(3,2),2)*SPA(2,5)*SPA(3,5)*SPB(4,1))/
(complex<T>(6,0)*pow(SPA(2,3),2)*pow(SPB(2,1),2)*pow(SPB(4,3),2))+
 (complex<T>(0,1)*SPA(1,5)*SPB(3,1)*SPB(3,2))/(complex<T>(2,0)*S(1,2)*SPB(4,3)*
 SPB(5,1))+(complex<T>(0,1)*SPA(1,5)*SPA(4,5)*SPB(3,1)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPB(5,1))+(complex<T>(0,-1)*S(3,4)*SPA(1,5)*SPA(4,5)*
 SPB(3,1)*SPB(5,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),3)*SPB(5,1))+
 (complex<T>(0,1)*SPA(4,5)*SPB(3,2)*SPB(4,2))/(complex<T>(2,0)*S(3,4)*SPB(2,1)*
 SPB(5,4))+(complex<T>(0,1)*pow(SPB(3,2),2)*pow(SPB(4,1),3)*SPA(1,2)*
 SPA(3,4))/(complex<T>(6,0)*pow(SPA(2,3),2)*pow(SPB(2,1),2)*
 pow(SPB(4,3),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPB(3,2),2)*pow(SPB(4,1),2)*SPA(1,2))/
(complex<T>(6,0)*pow(SPB(4,3),2)*SPA(2,3)*SPB(2,1)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*SPB(4,1))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,-1)*pow(SPB(3,2),2)*pow(SPB(4,1),2)*
 SPA(3,4))/(complex<T>(6,0)*pow(SPB(2,1),2)*SPA(2,3)*SPB(4,3)*SPB(5,1)*
 SPB(5,4))+(complex<T>(0,1)*pow(SPB(4,1),2)*SPA(1,4)*SPB(3,2))/
(complex<T>(6,0)*SPA(2,3)*SPB(2,1)*SPB(4,3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-23)*pow(SPB(3,2),2)*SPB(4,1))/(complex<T>(18,0)*SPB(2,1)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*S(1,2)*S(1,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*S(1,2)*S(4,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*S(1,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*S(4,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*SPA(1,5)*SPA(4,5)*SPB(4,2)*
 SPB(5,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),2)*pow(SPB(4,3),2)*SPB(5,4))+
 (complex<T>(0,-1)*S(1,2)*SPA(1,5)*SPA(4,5)*SPB(4,2)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),3)*SPB(5,4))+(complex<T>(0,-1)*S(1,5)*S(3,4)*SPA(1,5)*
 SPB(5,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*S(3,4)*S(4,5)*SPA(1,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*S(1,5)*SPA(1,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*SPB(2,1)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,-1)*S(4,5)*SPA(1,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*SPB(2,1)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*S(1,2)*SPA(1,5)*SPA(2,3)*SPA(4,5)*
 SPB(2,1)*SPB(5,2)*SPB(5,3))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),
3)*pow(SPA(3,4),4)*pow(SPB(4,3),3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*SPA(1,5)*SPA(4,5)*SPB(4,1)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPA(1,5)*SPA(4,5)*
 SPB(4,1)*SPB(5,2)*SPB(5,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),
3)*pow(SPA(3,4),4)*pow(SPB(4,3),4)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*S(1,2)*SPA(1,5)*SPA(4,5)*SPB(4,1)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),3)*
 pow(SPB(4,3),3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*SPA(1,5)*SPA(4,5)*SPB(4,1)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*S(1,2)*S(3,4)*SPB(5,1)*
 SPB(5,4))+(complex<T>(0,-1)*SPA(1,5)*SPA(2,3)*SPA(4,5)*SPB(2,1)*
 SPB(5,2)*SPB(5,3))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),2)*S(1,2)*SPB(4,3)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpQmp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, p}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpQmp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,3),3)*pow(SPA(4,5),2)*pow(SPB(5,3),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPA(1,5)*SPA(3,4)*SPA(3,5))+
 (complex<T>(0,1)*pow(SPA(1,3),3)*pow(SPA(4,5),2)*pow(SPB(5,3),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPA(1,5)*SPA(3,4)*SPA(3,5))+
 (complex<T>(0,-1)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPA(1,4)*SPA(1,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPA(1,2)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,4),3)*pow(SPB(5,2),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPA(1,2)*SPA(2,5)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,4),3)*pow(SPB(5,2),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPA(1,2)*SPA(2,5)*SPA(4,5))+
 (complex<T>(0,-7)*pow(SPA(1,4),3))/(complex<T>(9,0)*SPA(1,2)*SPA(1,5)*SPA(3,4)*
 SPA(4,5))+(complex<T>(0,1)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPA(1,4)*
 SPA(1,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(5,3),2)*SPA(1,4)*SPA(4,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPA(1,2)*SPA(1,5))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(5,3),2)*SPA(1,4)*SPA(4,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPA(1,5)*SPA(3,4))+
 (complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,4),3)*pow(SPB(5,2),2)*
 SPB(2,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3)*SPA(2,5)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPB(5,3),2)*SPA(1,4)*SPA(4,5)*
 SPB(2,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3)*SPA(1,5))+
 (complex<T>(0,1)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPA(1,4)*SPA(1,5))/
(complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPA(3,4),2)*SPA(4,5)*SPB(2,1)*
 SPB(4,3))+(complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPB(5,3),2)*SPA(1,4)*
 SPA(4,5))/(complex<T>(2,0)*pow(SPA(1,2),2)*pow(SPA(3,4),2)*SPA(1,5)*
 SPB(2,1)*SPB(4,3))+(complex<T>(0,-1)*pow(SPA(1,3),3)*pow(SPA(4,5),2)*
 pow(SPB(5,3),2)*SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),
2)*pow(SPA(1,2),4)*pow(SPB(2,1),3)*SPA(1,5)*SPA(3,5))+
 (complex<T>(0,-1)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPA(1,4)*SPA(1,5)*
 SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*pow(SPB(5,2),2)*SPA(1,4)*
 SPB(5,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3)*SPA(1,2)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,5),3)*pow(SPA(2,4),3)*pow(SPB(5,2),2)*
 SPB(5,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3)*SPA(1,2)*SPA(2,5)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPB(5,3),2)*S(1,2)*SPA(1,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),3)*
 pow(SPB(4,3),4))+(complex<T>(0,1)*pow(SPA(1,4),2)*S(3,5)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPA(1,5)*SPA(3,4))+
 (complex<T>(0,-1)*pow(SPB(5,3),2)*SPA(1,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPB(4,3),2)*S(1,2)*
 SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,4),2)*SPA(2,4)*SPB(5,2))/
(complex<T>(2,0)*pow(SPA(3,4),2)*SPA(1,2)*SPA(4,5)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*S(3,5)*SPB(4,3)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),3)*SPA(1,5))+(complex<T>(0,1)*pow(SPA(1,4),2)*S(2,5)*
 SPB(5,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),2)*pow(SPB(4,3),2)*SPA(1,2)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*SPA(1,3)*SPB(5,3))/
(complex<T>(2,0)*pow(SPA(1,2),2)*SPA(1,5)*SPA(3,4)*SPB(2,1))+
 (complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(4,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPB(2,1))+(complex<T>(0,-1)*pow(SPA(1,4),2)*S(2,5)*
 SPB(2,1)*SPB(5,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),3)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,2),2)*pow(SPB(5,2),2)*SPA(4,5)*SPB(2,1)*
 SPB(5,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),4)*pow(SPB(4,3),4))+
 (complex<T>(0,5)*SPA(1,4)*SPB(5,2)*SPB(5,3))/(complex<T>(6,0)*S(1,2)*S(3,4))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),2)*SPA(1,4)*
 SPB(5,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3)*SPA(1,5)*SPA(3,4))+
 (complex<T>(0,-1)*pow(SPA(1,3),3)*pow(SPA(4,5),3)*pow(SPB(5,3),2)*
 SPB(5,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3)*SPA(1,5)*SPA(3,4)*SPA(3,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmmqpQmQp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmmqpQmQp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(2,4),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2)*
 SPA(2,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*
 pow(SPA(1,5),4)*pow(SPB(5,1),5)*(complex<T>(1,0)-S(2,3)/S(1,5)-
S(3,4)/S(1,5)))+(complex<T>(0,-1)*pow(SPB(5,3),2)*SPA(2,3)*SPA(4,5))/
(complex<T>(2,0)*S(1,5)*S(3,4)*SPB(2,1))+
 (complex<T>(0,-1)*SPA(1,2)*SPA(1,4)*SPB(3,1))/(complex<T>(2,0)*S(2,3)*SPA(1,5)*
 SPB(2,1))+(complex<T>(0,-1)*pow(SPA(2,4),3)*pow(SPB(4,3),2)*
 pow(SPB(5,2),2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*
 pow(SPA(1,5),3)*pow(SPB(5,1),4)*(complex<T>(1,0)-S(2,3)/S(1,5)-
S(3,4)/S(1,5))*SPB(3,2))+(complex<T>(0,-1)*pow(SPA(2,4),3)*
 pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),3)*
 pow(SPB(5,1),4)*(complex<T>(1,0)-S(2,3)/S(1,5)-S(3,4)/S(1,5))*
 SPB(3,2))+(complex<T>(0,1)*pow(SPA(2,4),3)*pow(SPB(4,3),2)*
 pow(SPB(5,2),2)*S(3,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),4)*pow(SPB(5,1),5)*(complex<T>(1,0)-S(2,3)/S(1,5)-
S(3,4)/S(1,5))*SPB(3,2))+(complex<T>(0,-1)*pow(SPB(5,3),2)*SPA(2,4))/
(complex<T>(2,0)*pow(SPB(5,1),2)*SPA(1,5)*SPB(3,2))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*SPB(3,1))/(complex<T>(2,0)*SPA(1,5)*SPA(3,4)*
 SPB(2,1)*SPB(3,2))+(complex<T>(0,-1)*pow(SPA(2,3),2)*pow(SPB(5,3),2)*
 SPB(3,1)*SPB(3,2))/(complex<T>(2,0)*pow(SPB(4,3),2)*pow(SPB(5,1),2)*
 SPA(1,5)*SPA(3,4)*SPB(2,1))+(complex<T>(0,-1)*pow(SPA(2,4),2)*
 pow(SPB(5,3),2)*SPA(2,3)*SPB(4,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),3)*
 pow(SPB(5,1),4))+(complex<T>(0,1)*pow(SPA(2,4),2)*pow(SPB(5,3),2)*
 SPB(4,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*
 pow(SPA(1,5),2)*pow(SPB(5,1),3)*SPB(3,2))+
 (complex<T>(0,-1)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*pow(SPB(5,4),2)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*
 pow(SPA(1,5),2)*pow(SPB(5,1),3)*SPB(2,1)*SPB(4,2)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*pow(SPB(5,4),2)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),2)*pow(SPB(5,1),3)*SPB(2,1)*SPB(4,2)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(2,4),2)*pow(SPB(3,2),2)*pow(SPB(5,4),2)*S(3,4)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),3)*pow(SPB(5,1),4)*SPB(2,1)*SPB(4,2)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPA(2,4),2)*pow(SPB(3,2),3)*pow(SPB(5,4),2)*SPA(2,3)*
 SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*
 pow(SPA(1,5),3)*pow(SPB(5,1),4)*SPB(2,1)*SPB(4,2)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPB(5,3),2)*SPA(2,4))/(complex<T>(2,0)*S(3,4)*SPB(3,2)*
 SPB(5,1))+(complex<T>(0,-1)*pow(SPA(2,3),2)*pow(SPB(5,3),2)*SPB(3,1)*
 SPB(3,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),2)*pow(SPB(4,3),3)*SPB(2,1)*SPB(5,1))+
 (complex<T>(0,1)*pow(SPA(2,3),2)*pow(SPB(5,3),2)*S(1,5)*SPB(3,1)*
 SPB(3,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),4)*SPB(2,1)*SPB(5,1))+
 (complex<T>(0,1)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*S(1,5)*SPB(3,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),2)*SPB(2,1)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(2,4),2)*pow(SPB(5,4),2)*SPB(3,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),2)*SPB(2,1)*
 SPB(3,2)*SPB(4,3)*SPB(5,1))+(complex<T>(0,-1)*SPA(1,2)*SPA(2,4)*
 SPB(3,2)*SPB(5,1))/(complex<T>(2,0)*S(4,5)*(-S(2,3)+S(4,5))*
 SPB(2,1))+(complex<T>(0,1)*pow(SPA(2,4),2)*S(1,5)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(2,1))+(complex<T>(0,-1)*pow(SPA(2,4),2)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),2)*SPB(2,1)*
 SPB(4,3))+(complex<T>(0,-1)*SPA(2,4)*SPB(5,3))/(S(4,5)*SPB(2,1))+
 (complex<T>(0,1)*SPA(1,4)*SPA(2,3)*SPB(3,1)*SPB(5,3))/
(complex<T>(2,0)*S(4,5)*(-S(2,3)+S(4,5))*SPB(2,1))+
 (complex<T>(0,1)*SPA(1,2)*SPA(2,4)*SPB(3,2)*SPB(4,1)*SPB(5,3))/
(pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),2)*pow(SPB(5,1),2)*
 SPB(2,1)*SPB(4,3))+(complex<T>(0,-1)*S(3,4)*SPA(1,2)*SPA(2,4)*SPB(3,2)*
 SPB(4,1)*SPB(5,3))/(pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),3)*pow(SPB(5,1),3)*SPB(2,1)*SPB(4,3))+
 (complex<T>(0,1)*SPA(1,4)*SPA(2,3)*SPA(2,4)*SPB(4,3)*SPB(5,3))/
(pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),3)*
 pow(SPB(5,1),3))+(complex<T>(0,-1)*SPA(1,4)*SPA(2,4)*SPB(4,3)*
 SPB(5,3))/(pow(complex<T>(1,0)-S(2,3)/S(1,5),2)*pow(SPA(1,5),2)*
 pow(SPB(5,1),2)*SPB(3,2))+(complex<T>(0,-7)*pow(SPB(5,3),2))/
(complex<T>(9,0)*SPB(2,1)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(2,4),2)*SPB(5,3)*SPB(5,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPB(5,1))+(complex<T>(0,-1)*pow(SPA(2,4),2)*S(1,5)*
 SPB(5,3)*SPB(5,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),3)*SPB(5,1))+
 (complex<T>(0,-1)*SPA(2,4)*SPB(3,1)*SPB(5,3)*SPB(5,4))/
(complex<T>(2,0)*S(2,3)*SPB(2,1)*SPB(4,3)*SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmpqpQmQp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, Qm, Qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmpqpQmQp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*pow(SPB(5,2),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-S(1,5)/S(3,4))*
 SPA(1,2))+(complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*
 pow(SPB(5,2),3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
S(1,5)/S(3,4))*SPA(1,2))+(complex<T>(0,1)*pow(SPB(5,2),2)*S(3,4)*
 SPA(2,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),2)*pow(SPB(5,1),3)*SPA(2,3))+
 (complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(2,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPB(5,1),2)*SPA(1,5)*
 SPA(2,3))+(complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(5,2),2)*SPA(2,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPA(1,2))+(complex<T>(0,-1)*pow(SPA(1,4),2)*
 pow(SPB(2,1),2)*SPA(1,2)*SPA(1,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),3)*
 pow(SPB(5,1),2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,-1)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(1,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(1,5)*SPA(2,3)*SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,2),2)*
 pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPA(1,5)*SPA(2,3)*SPA(2,5))+
 (complex<T>(0,-1)*pow(SPA(1,2),2)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*
 SPA(3,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPA(1,5)*SPA(2,3)*SPA(2,5))+
 (complex<T>(0,-7)*pow(SPA(1,4),2))/(complex<T>(9,0)*SPA(1,2)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPB(5,2),2)*SPA(1,4)*SPA(4,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),2)*
 pow(SPB(5,1),2)*SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,5),2)*
 pow(SPA(2,4),2)*pow(SPB(5,2),3)*SPB(2,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),5)*
 pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-S(1,5)/S(3,4)))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(5,2),2)*SPA(2,5)*SPB(2,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3))+(complex<T>(0,1)*pow(SPA(1,2),3)*pow(SPA(4,5),2)*
 pow(SPB(5,2),2)*SPA(3,5)*SPB(2,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3)*SPA(1,5)*SPA(2,3)*SPA(2,5))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(2,1),2)*SPA(1,2)*SPA(1,3)*
 SPB(4,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*
 pow(SPA(1,5),4)*pow(SPB(5,1),3)*SPA(2,3))+
 (complex<T>(0,1)*pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(1,3)*SPB(4,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(1,5)*SPA(2,3))+
 (complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(1,4)*SPA(4,5)*SPB(4,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,5),2)*pow(SPA(1,5),3)*
 pow(SPB(5,1),3))+(complex<T>(0,-1)*pow(SPB(5,3),2)*SPA(1,3))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPB(4,3)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(2,1),2)*SPA(1,2)*SPA(1,3))/
(complex<T>(2,0)*pow(SPA(1,5),2)*pow(SPA(3,4),2)*SPA(2,3)*SPB(4,3)*
 SPB(5,1))+(complex<T>(0,1)*pow(SPA(1,5),3)*pow(SPA(2,4),2)*
 pow(SPB(5,2),3)*SPB(5,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),
2)*pow(SPA(3,4),5)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
S(1,5)/S(3,4))*SPA(1,2))+(complex<T>(0,1)*pow(SPA(1,2),2)*
 pow(SPA(4,5),2)*pow(SPB(5,2),2)*SPA(3,5)*SPB(5,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3)*SPA(2,3)*SPA(2,5))+
 (complex<T>(0,-1)*SPA(1,4)*SPB(5,2))/(S(4,5)*SPA(2,3))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*SPB(5,2))/(complex<T>(2,0)*S(1,5)*SPA(1,2)*
 SPA(3,4))+(complex<T>(0,-1)*SPA(1,3)*SPA(1,4)*SPA(4,5)*SPB(5,2))/
(complex<T>(2,0)*S(1,2)*SPA(1,5)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,-1)*SPA(1,2)*SPA(3,4)*SPB(3,2)*SPB(5,2))/
(complex<T>(2,0)*S(4,5)*(-S(1,2)+S(4,5))*SPA(2,3))+
 (complex<T>(0,1)*SPA(1,2)*SPA(1,4)*SPA(3,5)*SPB(3,2)*SPB(5,2))/
(pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),2)*pow(SPB(4,3),2)*
 SPA(1,5)*SPA(2,3))+(complex<T>(0,-1)*pow(SPA(1,4),2)*SPB(5,2))/
(complex<T>(2,0)*pow(SPA(3,4),2)*SPA(1,2)*SPB(4,3))+
 (complex<T>(0,-1)*SPA(1,2)*SPA(1,4)*SPA(3,5)*SPB(3,2)*SPB(5,1)*SPB(5,2))/
(pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),3)*pow(SPB(4,3),3)*
 SPA(2,3))+(complex<T>(0,1)*SPA(1,3)*SPA(1,4)*SPB(2,1)*SPB(5,3))/
(complex<T>(2,0)*S(4,5)*(-S(1,2)+S(4,5))*SPA(2,3))+
 (complex<T>(0,-1)*SPA(1,3)*SPB(3,2)*SPB(5,3))/(complex<T>(2,0)*S(1,2)*SPA(2,3)*
 SPB(4,3))+(complex<T>(0,-1)*SPA(1,4)*SPA(1,5)*SPB(5,2)*SPB(5,3))/
(pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*pow(SPB(4,3),2)*
 SPA(1,2))+(complex<T>(0,1)*SPA(1,4)*SPA(1,5)*SPB(2,1)*SPB(5,2)*
 SPB(5,3))/(pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),3))+(complex<T>(0,-1)*pow(SPA(1,4),2)*SPB(2,1)*SPB(5,4))/
(complex<T>(2,0)*S(1,5)*S(3,4)*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmmqpQpQm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmmqpQpQm LLT");
#endif
 
 return( (complex<T>(0,1)*SPA(1,2)*SPA(2,5)*SPB(3,2)*SPB(4,1))/
(complex<T>(2,0)*S(4,5)*(-S(2,3)+S(4,5))*SPB(2,1))+
 (complex<T>(0,1)*SPA(2,5)*SPB(4,3))/(S(4,5)*SPB(2,1))+
 (complex<T>(0,-1)*SPA(1,5)*SPA(2,3)*SPB(3,1)*SPB(4,3))/
(complex<T>(2,0)*S(4,5)*(-S(2,3)+S(4,5))*SPB(2,1))+
 (complex<T>(0,-7)*pow(SPB(4,3),2))/(complex<T>(9,0)*SPB(2,1)*SPB(3,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmpqpQpQm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, Qp, Qm}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmpqpQpQm LLT");
#endif
 
 return( (complex<T>(0,-7)*pow(SPA(1,5),2))/(complex<T>(9,0)*SPA(1,2)*SPA(2,3)*
 SPA(4,5))+(complex<T>(0,1)*SPA(1,5)*SPB(4,2))/(S(4,5)*SPA(2,3))+
 (complex<T>(0,1)*SPA(1,2)*SPA(3,5)*SPB(3,2)*SPB(4,2))/
(complex<T>(2,0)*S(4,5)*(-S(1,2)+S(4,5))*SPA(2,3))+
 (complex<T>(0,-1)*SPA(1,3)*SPA(1,5)*SPB(2,1)*SPB(4,3))/
(complex<T>(2,0)*S(4,5)*(-S(1,2)+S(4,5))*SPA(2,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpmQmQp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpmQmQp LRT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,3),2)*S(4,5)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(4,3))+(complex<T>(0,1)*pow(SPA(1,3),2)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*SPB(2,1)*
 SPB(4,3))+(complex<T>(0,1)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqppQmQp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqppQmQp LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPA(4,5),2)*pow(SPB(5,3),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))*
 SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,3),2)*pow(SPA(4,5),2)*
 pow(SPB(5,3),3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2))*SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,3),2)*
 pow(SPA(4,5),2)*pow(SPB(5,3),3)*S(3,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),5)*
 pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))*
 SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,3),2)*pow(SPA(4,5),2)*
 pow(SPB(5,3),3)*S(4,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),5)*pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2))*SPA(3,4))+(complex<T>(0,-1)*pow(SPA(1,4),2)*
 pow(SPB(5,3),2)*SPA(3,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),
2)*pow(SPA(1,2),3)*pow(SPB(2,1),2)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(5,3),2)*S(3,4)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3)*SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,4),2)*SPA(2,4))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*SPB(5,3))/(complex<T>(2,0)*pow(SPA(1,2),2)*
 SPA(3,4)*SPB(2,1))+(complex<T>(0,1)*SPA(1,4)*SPA(4,5)*SPB(5,2)*
 SPB(5,3))/(pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPA(3,4))+(complex<T>(0,-1)*S(3,4)*SPA(1,4)*SPA(4,5)*
 SPB(5,2)*SPB(5,3))/(pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),3)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpmQpQm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpmQpQm LRT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(3,5),3)*pow(SPB(3,2),2)*pow(SPB(5,4),2)*
 SPA(3,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),5)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2)))+(complex<T>(0,1)*pow(SPA(3,5),3)*pow(SPB(3,2),2)*
 pow(SPB(5,4),2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2))*SPB(4,3))+(complex<T>(0,1)*pow(SPA(3,5),3)*
 pow(SPB(3,2),2)*pow(SPB(5,4),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))*
 SPB(4,3))+(complex<T>(0,1)*pow(SPB(4,2),2)*SPA(3,5))/
(complex<T>(2,0)*pow(SPB(2,1),2)*SPA(1,2)*SPB(4,3))+
 (complex<T>(0,-1)*pow(SPA(3,5),3)*pow(SPB(3,2),2)*pow(SPB(5,4),3)*
 SPA(4,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),5)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2))*SPB(4,3))+(complex<T>(0,1)*pow(SPA(3,5),2)*
 pow(SPB(4,2),2)*SPA(3,4)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4))+(complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(4,2),2)*
 SPB(5,3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*
 pow(SPA(1,2),2)*pow(SPB(2,1),3)*SPB(4,3))+
 (complex<T>(0,1)*pow(SPB(4,2),3))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*SPB(4,3)*
 SPB(5,4))+(complex<T>(0,1)*SPA(1,5)*SPA(3,4)*SPA(3,5)*SPB(4,2)*
 SPB(5,4))/(pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),3))+(complex<T>(0,-1)*SPA(1,5)*SPA(3,5)*SPB(4,2)*
 SPB(5,4))/(pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqppQpQm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqppQpQm LRT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(3,2),2)*S(4,5)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,5),2)*SPA(2,4))/
(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*S(1,2)*SPA(3,4)*SPB(2,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmmQp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, m, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmmQp LRT");
#endif
 
 return( (complex<T>(0,1)*S(4,5)*SPA(1,4)*SPA(3,4)*SPB(3,2))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(4,5))*SPA(4,5)*SPB(4,3))+
 (complex<T>(0,1)*SPA(1,4)*SPB(5,2))/((S(1,2)-S(4,5))*SPB(4,3))+
 (complex<T>(0,-1)*S(4,5)*SPA(1,4)*SPB(5,2))/(complex<T>(2,0)*S(1,2)*
 (S(1,2)-S(4,5))*SPB(4,3))+(complex<T>(0,1)*pow(SPB(5,2),2))/
(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPB(5,2),2)*S(4,5))/(complex<T>(2,0)*(S(1,2)-S(4,5))*
 SPB(2,1)*SPB(4,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmpQp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, p, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmpQp LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,3),2))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*
 SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,3),2)*S(3,4))/
(complex<T>(2,0)*(S(1,2)-S(3,4))*SPA(1,2)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*SPA(1,3)*SPB(4,2))/((S(1,2)-S(3,4))*SPA(4,5))+
 (complex<T>(0,-1)*S(3,4)*SPA(1,3)*SPB(4,2))/(complex<T>(2,0)*S(1,2)*
 (S(1,2)-S(3,4))*SPA(4,5))+(complex<T>(0,1)*S(3,4)*SPA(1,5)*SPB(4,2)*
 SPB(5,4))/(complex<T>(2,0)*S(1,2)*(S(1,2)-S(3,4))*SPA(4,5)*SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpmQm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, m, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpmQm LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(3,2),2)*SPA(3,4))/(complex<T>(2,0)*(S(1,2)-S(3,4))*
 SPB(2,1)*SPB(5,4))+(complex<T>(0,-1)*SPA(1,4)*SPB(3,2))/
((S(1,2)-S(3,4))*SPB(5,4))+(complex<T>(0,1)*S(3,4)*SPA(1,4)*SPB(3,2))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(3,4))*SPB(5,4))+
 (complex<T>(0,1)*pow(SPB(3,2),2))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,-1)*SPA(1,4)*SPA(4,5)*SPB(4,3)*SPB(5,2))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(3,4))*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQppQm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, p, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQppQm LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,5),2))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*
 SPA(4,5))+(complex<T>(0,-1)*SPA(1,5)*SPB(4,2))/((S(1,2)-S(4,5))*
 SPA(3,4))+(complex<T>(0,1)*S(4,5)*SPA(1,5)*SPB(4,2))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(4,5))*SPA(3,4))+
 (complex<T>(0,-1)*SPA(1,3)*SPA(4,5)*SPB(4,2)*SPB(4,3))/
(complex<T>(2,0)*S(1,2)*(S(1,2)-S(4,5))*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*SPB(5,4))/(complex<T>(2,0)*(S(1,2)-S(4,5))*
 SPA(1,2)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmQpm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, m}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmQpm LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(3,5),3)*pow(SPB(4,3),2)*pow(SPB(5,2),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))*
 SPB(5,4))+(complex<T>(0,1)*pow(SPA(3,5),3)*pow(SPB(4,3),2)*
 pow(SPB(5,2),2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),3)*pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2))*SPB(5,4))+(complex<T>(0,-1)*pow(SPA(3,5),3)*
 pow(SPB(4,3),2)*pow(SPB(5,2),2)*S(3,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),4)*
 pow(SPB(2,1),5)*(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))*
 SPB(5,4))+(complex<T>(0,-1)*pow(SPA(3,5),3)*pow(SPB(4,3),2)*
 pow(SPB(5,2),2)*S(4,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),5)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2))*SPB(5,4))+(complex<T>(0,1)*pow(SPB(4,2),2)*SPA(3,5))/
(complex<T>(2,0)*pow(SPB(2,1),2)*SPA(1,2)*SPB(5,4))+
 (complex<T>(0,1)*SPA(1,3)*SPA(3,5)*SPB(4,2)*SPB(4,3))/
(pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*
 SPB(5,4))+(complex<T>(0,-1)*S(4,5)*SPA(1,3)*SPA(3,5)*SPB(4,2)*
 SPB(4,3))/(pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),3)*SPB(5,4))+(complex<T>(0,1)*pow(SPB(4,2),2)*SPB(4,1))/
(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*pow(SPB(4,2),2)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPB(5,4))+(complex<T>(0,1)*pow(SPA(3,5),2)*
 pow(SPB(4,2),2)*S(4,5)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmQpp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, p}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmQpp LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,3),2)*SPA(1,4))/(complex<T>(2,0)*SPA(1,2)*SPA(1,5)*
 SPA(3,4)*SPA(4,5))+(complex<T>(0,-1)*pow(SPB(5,2),2)*S(3,4)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(4,5))+(complex<T>(0,1)*pow(SPB(5,2),2)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpQmm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, m}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpQmm LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(3,2),2)*SPB(4,1))/(complex<T>(2,0)*SPB(2,1)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,-1)*pow(SPA(1,5),2)*S(3,4)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(5,4))+(complex<T>(0,1)*pow(SPA(1,5),2)*SPB(5,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*S(1,2)*SPA(1,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpQmp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, p}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpQmp LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(3,4),2)*pow(SPB(5,3),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))*
 SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(3,4),2)*
 pow(SPB(5,3),3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2))*SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,4),3))/
(complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(5,3),2)*SPA(3,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,5),2)*
 pow(SPA(3,4),3)*pow(SPB(5,3),3)*SPB(4,3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),5)*
 pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-S(4,5)/S(1,2))*
 SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,4),2)*SPB(5,3))/
(complex<T>(2,0)*pow(SPA(1,2),2)*SPA(4,5)*SPB(2,1))+
 (complex<T>(0,-1)*SPA(1,4)*SPA(3,4)*SPB(3,2)*SPB(5,3))/
(pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(3,4),2)*
 pow(SPB(5,3),3)*SPB(5,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),
2)*pow(SPA(1,2),5)*pow(SPB(2,1),4)*(complex<T>(1,0)-S(3,4)/S(1,2)-
S(4,5)/S(1,2)))+(complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(5,3),2)*
 SPA(3,5)*SPB(5,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*
 pow(SPA(1,2),4)*pow(SPB(2,1),3))+
 (complex<T>(0,1)*SPA(1,4)*SPA(3,4)*SPB(3,2)*SPB(5,3)*SPB(5,4))/
(pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*pow(SPB(2,1),3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmmqpQmQp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmmqpQmQp LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(5,3),2))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmpqpQmQp_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, Qm, Qp}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmpqpQmQp LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,4),2))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmmqpQpQm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmmqpQpQm LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(4,3),2))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmpqpQpQm_LRT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, Qp, Qm}, LRT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmpqpQpQm LRT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,5),2))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmmQmQpqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, Qm, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmmQmQpqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,4),2)*SPB(3,1))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*
 SPB(4,3)*SPB(5,1))+(complex<T>(0,1)*pow(SPA(2,3),2)*S(1,5)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(2,1))+(complex<T>(0,-1)*pow(SPA(2,3),2)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*S(3,4)*SPA(3,4)*SPB(2,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmmQpQmqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, Qp, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmmQpQmqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,3),2)*SPB(3,1))/(complex<T>(2,0)*SPB(2,1)*SPB(3,2)*
 SPB(4,3)*SPB(5,1))+(complex<T>(0,1)*pow(SPA(2,4),2)*S(1,5)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(2,1))+(complex<T>(0,-1)*pow(SPA(2,4),2)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),2)*SPB(2,1)*
 SPB(4,3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmpQmQpqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, Qm, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmpQmQpqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,3),2)*pow(SPB(5,2),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-S(1,5)/S(3,4))*
 SPA(1,2))+(complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,3),2)*
 pow(SPB(5,2),3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
S(1,5)/S(3,4))*SPA(1,2))+(complex<T>(0,1)*pow(SPA(1,5),2)*
 pow(SPA(2,3),2)*pow(SPB(5,2),3)*S(1,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),5)*
 pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-S(1,5)/S(3,4))*
 SPA(1,2))+(complex<T>(0,-1)*pow(SPA(1,3),3))/(complex<T>(2,0)*SPA(1,2)*SPA(1,5)*
 SPA(2,3)*SPA(3,4))+(complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,3),2)*
 pow(SPB(5,2),3)*SPB(2,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),
2)*pow(SPA(3,4),5)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
S(1,5)/S(3,4)))+(complex<T>(0,1)*pow(SPA(1,3),2)*S(2,5)*SPB(5,2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPA(1,2))+(complex<T>(0,-1)*pow(SPA(1,3),2)*S(2,5)*
 SPB(2,1)*SPB(5,2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 SPA(1,2)*SPB(4,3))+(complex<T>(0,1)*SPA(1,3)*SPA(1,5)*SPB(5,2)*
 SPB(5,4))/(pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPA(1,2))+(complex<T>(0,-1)*SPA(1,3)*SPA(1,5)*SPB(2,1)*
 SPB(5,2)*SPB(5,4))/(pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmpQpQmqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, Qp, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmpQpQmqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*pow(SPB(5,2),3))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-S(1,5)/S(3,4))*
 SPA(1,2))+(complex<T>(0,-1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*
 pow(SPB(5,2),3))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3)*(complex<T>(1,0)-S(1,2)/S(3,4)-
S(1,5)/S(3,4))*SPA(1,2))+(complex<T>(0,1)*pow(SPA(1,5),2)*
 pow(SPA(2,4),2)*pow(SPB(5,2),3)*S(1,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(3,4),2)*pow(SPA(3,4),5)*
 pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-S(1,5)/S(3,4))*
 SPA(1,2))+(complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(5,2),2)*SPA(2,5))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPA(1,2))+(complex<T>(0,-1)*pow(SPA(1,4),2)*SPA(1,3))/
(complex<T>(2,0)*SPA(1,2)*SPA(1,5)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*pow(SPB(5,2),3)*
 SPB(2,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),5)*pow(SPB(4,3),4)*(complex<T>(1,0)-S(1,2)/S(3,4)-
S(1,5)/S(3,4)))+(complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(5,2),2)*
 SPA(2,5)*SPB(2,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),4)*pow(SPB(4,3),3))+
 (complex<T>(0,-1)*pow(SPA(1,4),2)*SPB(5,2))/(complex<T>(2,0)*pow(SPA(3,4),2)*
 SPA(1,2)*SPB(4,3))+(complex<T>(0,-1)*SPA(1,4)*SPA(1,5)*SPB(5,2)*
 SPB(5,3))/(pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPA(1,2))+(complex<T>(0,1)*SPA(1,4)*SPA(1,5)*SPB(2,1)*
 SPB(5,2)*SPB(5,3))/(pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),3))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQmmQpqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, m, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQmmQpqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,4),2))/(complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQpmQmqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, m, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQpmQmqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,2),2))/(complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQmpQpqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, p, Qp, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQmpQpqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,2),2))/(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQppQmqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, p, Qm, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQppQmqp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,4),2))/(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQmQpqpm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, m}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQmQpqpm LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(4,3),2)*SPA(4,5))/(complex<T>(2,0)*(S(2,3)-S(4,5))*
 SPB(3,2)*SPB(5,1))+(complex<T>(0,1)*SPA(2,5)*SPB(4,3))/
((S(2,3)-S(4,5))*SPB(5,1))+(complex<T>(0,-1)*S(4,5)*SPA(2,5)*
 SPB(4,3))/(complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPB(4,3),2))/(complex<T>(2,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*SPA(1,5)*SPA(2,5)*SPB(3,1)*SPB(5,4))/
(complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQmQpqpp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, p}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQmQpqpp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,2),2))/(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,2),2)*S(1,5))/
(complex<T>(2,0)*(-S(1,5)+S(2,3))*SPA(1,5)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,1)*SPA(1,2)*SPB(5,3))/((-S(1,5)+S(2,3))*SPA(4,5))+
 (complex<T>(0,-1)*S(1,5)*SPA(1,2)*SPB(5,3))/(complex<T>(2,0)*S(2,3)*
 (-S(1,5)+S(2,3))*SPA(4,5))+(complex<T>(0,1)*S(1,5)*SPA(2,4)*SPB(5,3)*
 SPB(5,4))/(complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPA(4,5)*SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQpQmqpm_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, m}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQpQmqpm LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(4,2),2)*SPA(4,5))/(complex<T>(2,0)*(S(2,3)-S(4,5))*
 SPB(3,2)*SPB(5,1))+(complex<T>(0,-1)*SPA(3,5)*SPB(4,2))/
((S(2,3)-S(4,5))*SPB(5,1))+(complex<T>(0,1)*S(4,5)*SPA(3,5)*SPB(4,2))/
(complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPB(4,2),2))/(complex<T>(2,0)*SPB(3,2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*SPA(1,5)*SPA(3,5)*SPB(2,1)*SPB(5,4))/
(complex<T>(2,0)*S(2,3)*(S(2,3)-S(4,5))*SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQpQmqpp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, p}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQpQmqpp LLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPA(1,3),2))/(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,3),2)*S(1,5))/
(complex<T>(2,0)*(-S(1,5)+S(2,3))*SPA(1,5)*SPA(2,3)*SPA(4,5))+
 (complex<T>(0,-1)*SPA(1,3)*SPB(5,2))/((-S(1,5)+S(2,3))*SPA(4,5))+
 (complex<T>(0,1)*S(1,5)*SPA(1,3)*SPB(5,2))/(complex<T>(2,0)*S(2,3)*
 (-S(1,5)+S(2,3))*SPA(4,5))+
 (complex<T>(0,-1)*S(1,5)*SPA(3,4)*SPB(5,2)*SPB(5,4))/
(complex<T>(2,0)*S(2,3)*(-S(1,5)+S(2,3))*SPA(4,5)*SPB(5,1))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQmQpmqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, m, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQmQpmqp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,4),3)*pow(SPB(4,3),2)*pow(SPB(5,1),2)*
 SPA(4,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),4)*pow(SPB(3,2),5)*(complex<T>(1,0)-S(1,5)/S(2,3)-
S(4,5)/S(2,3)))+(complex<T>(0,-1)*pow(SPB(5,3),2)*S(1,4)*SPA(1,4)*
 SPA(4,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),4))+
 (complex<T>(0,-1)*SPA(1,2)*SPA(1,4)*SPA(4,5)*SPB(5,1)*SPB(5,3))/
(pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),3))+(complex<T>(0,-1)*pow(SPA(1,4),3)*pow(SPB(4,3),2)*
 pow(SPB(5,1),2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-
S(4,5)/S(2,3))*SPB(5,4))+(complex<T>(0,-1)*pow(SPA(1,4),3)*
 pow(SPB(4,3),2)*pow(SPB(5,1),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-S(4,5)/S(2,3))*
 SPB(5,4))+(complex<T>(0,1)*pow(SPA(1,4),3)*pow(SPB(4,3),2)*
 pow(SPB(5,1),2)*S(1,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*
 pow(SPA(2,3),4)*pow(SPB(3,2),5)*(complex<T>(1,0)-S(1,5)/S(2,3)-
S(4,5)/S(2,3))*SPB(5,4))+(complex<T>(0,1)*pow(SPB(5,3),2)*S(1,4)*
 SPA(1,4))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),2)*pow(SPB(3,2),3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPB(5,3),2)*SPA(1,4))/(complex<T>(2,0)*pow(SPB(3,2),2)*
 SPA(2,3)*SPB(5,4))+(complex<T>(0,-1)*pow(SPB(5,3),3))/
(complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*SPA(1,2)*SPA(1,4)*SPB(5,1)*SPB(5,3))/
(pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),2)*pow(SPB(3,2),2)*
 SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQmQppqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, p, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQmQppqp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(4,3),2)*S(1,5)*SPA(1,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),2)*
 pow(SPB(3,2),3)*SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,2),2)*SPA(3,5))/
(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*pow(SPB(4,3),2)*SPA(1,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*S(2,3)*SPA(4,5)*SPB(3,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQpQmmqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, m, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQpQmmqp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(1,4),3)*pow(SPB(4,2),2)*pow(SPB(5,1),2)*
 SPA(4,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),4)*pow(SPB(3,2),5)*(complex<T>(1,0)-S(1,5)/S(2,3)-
S(4,5)/S(2,3)))+(complex<T>(0,-1)*pow(SPA(1,4),2)*pow(SPB(5,2),2)*
 SPA(4,5)*SPB(4,1))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),4))+
 (complex<T>(0,1)*SPA(1,3)*SPA(1,4)*SPA(4,5)*SPB(5,1)*SPB(5,2))/
(pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),3))+(complex<T>(0,-1)*pow(SPA(1,4),3)*pow(SPB(4,2),2)*
 pow(SPB(5,1),2))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*
 pow(SPA(2,3),3)*pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-
S(4,5)/S(2,3))*SPB(5,4))+(complex<T>(0,-1)*pow(SPA(1,4),3)*
 pow(SPB(4,2),2)*pow(SPB(5,1),2))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),3)*
 pow(SPB(3,2),4)*(complex<T>(1,0)-S(1,5)/S(2,3)-S(4,5)/S(2,3))*
 SPB(5,4))+(complex<T>(0,1)*pow(SPA(1,4),3)*pow(SPB(4,2),2)*
 pow(SPB(5,1),2)*S(1,5))/(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*
 pow(SPA(2,3),4)*pow(SPB(3,2),5)*(complex<T>(1,0)-S(1,5)/S(2,3)-
S(4,5)/S(2,3))*SPB(5,4))+(complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(1,4))/
(complex<T>(2,0)*pow(SPB(3,2),2)*SPA(2,3)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(1,4),2)*pow(SPB(5,2),2)*SPB(4,1))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),2)*
 pow(SPB(3,2),3)*SPB(5,4))+(complex<T>(0,-1)*SPA(1,3)*SPA(1,4)*SPB(5,1)*
 SPB(5,2))/(pow(complex<T>(1,0)-S(4,5)/S(2,3),2)*pow(SPA(2,3),2)*
 pow(SPB(3,2),2)*SPB(5,4))+(complex<T>(0,-1)*pow(SPB(5,2),2)*SPB(5,3))/
(complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmQpQmpqp_LLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, p, qp}, LLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmQpQmpqp LLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(4,2),2)*S(1,5)*SPA(1,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPA(2,3),2)*
 pow(SPB(3,2),3)*SPA(4,5))+(complex<T>(0,-1)*pow(SPB(4,2),2)*SPA(1,4))/
(complex<T>(2,0)*pow(complex<T>(1,0)-S(1,5)/S(2,3),2)*pow(SPB(3,2),2)*SPA(2,3)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPA(3,5))/
(complex<T>(2,0)*SPA(1,5)*SPA(2,3)*SPA(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpmQmQp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpmQmQp nfLLT");
#endif
 
 return( (complex<T>(0,-1)*SPA(1,3)*SPA(1,4)*SPB(5,2))/(complex<T>(3,0)*S(4,5)*SPA(1,2)*
 SPB(3,2))+(complex<T>(0,1)*pow(SPA(1,3),2)*S(4,5)*SPA(3,4)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3))+(complex<T>(0,-1)*pow(SPA(1,3),2)*SPA(3,4)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),2)*S(4,5)*
 SPB(2,1))+(complex<T>(0,-2)*pow(SPB(5,2),2)*SPB(4,2))/
(complex<T>(9,0)*SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*S(4,5)*SPB(5,2)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(3,2)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(1,3),2)*SPB(5,2)*SPB(5,3))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*SPB(2,1)*
 SPB(3,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqppQmQp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqppQmQp nfLLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,3),2)*S(1,2)*SPA(1,3)*SPA(1,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(1,2)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPB(5,3),2)*SPA(1,3)*SPA(1,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPB(5,4),2)*SPA(1,2)*
 SPA(3,4)*SPA(4,5))+(complex<T>(0,-2)*pow(SPA(1,4),2)*SPA(2,4))/
(complex<T>(9,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,1)*pow(SPB(5,3),2)*S(1,2)*SPA(1,3)*SPB(3,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPA(4,5),3)*
 pow(SPB(5,4),4))+(complex<T>(0,-1)*pow(SPB(5,3),2)*SPA(1,3)*SPB(3,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),3)*pow(SPB(5,4),2)*S(1,2)*
 SPA(4,5))+(complex<T>(0,-1)*SPA(1,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(3,0)*S(1,2)*SPA(3,4)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpmQpQm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, m, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpmQpQm nfLLT");
#endif
 
 return( (complex<T>(0,1)*SPA(1,3)*SPA(3,5)*SPB(4,2))/(complex<T>(6,0)*S(1,2)*S(4,5))+
 (complex<T>(0,-1)*pow(SPA(1,3),2)*S(4,5)*SPA(3,5)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),4)*
 pow(SPB(2,1),3))+(complex<T>(0,1)*pow(SPA(1,3),2)*SPA(3,5)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),2)*S(4,5)*
 SPB(2,1))+(complex<T>(0,-2)*pow(SPB(4,2),3))/(complex<T>(9,0)*SPB(2,1)*SPB(3,2)*
 SPB(4,3)*SPB(5,4))+(complex<T>(0,1)*pow(SPA(3,5),2)*SPA(1,3)*SPB(3,2)*
 SPB(5,4))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*
 pow(SPA(1,2),2)*pow(SPB(2,1),2)*S(4,5))+
 (complex<T>(0,-1)*pow(SPA(3,5),2)*S(4,5)*SPA(1,3)*SPB(3,2)*SPB(5,4))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),4)*
 pow(SPB(2,1),4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqppQpQm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, p, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqppQpQm nfLLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPB(3,2),2)*S(4,5)*SPA(3,5))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(3,4))+(complex<T>(0,-1)*pow(SPB(3,2),2)*SPA(3,5))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(3,4))+(complex<T>(0,-2)*pow(SPA(1,5),2)*SPA(2,4))/
(complex<T>(9,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*S(3,4)*S(4,5)*SPA(3,5)*SPB(3,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,1)*S(3,4)*SPA(3,5)*SPB(3,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),2)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(2,3)*SPA(3,4))+(complex<T>(0,1)*S(1,2)*S(3,4)*SPA(1,3)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPA(4,5),2)*
 pow(SPB(5,4),3)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,-1)*S(3,4)*SPA(1,3)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPB(5,4),2)*SPA(2,3)*
 SPA(3,4)*SPA(4,5))+(complex<T>(0,-1)*S(1,2)*SPA(1,3)*SPB(3,2)*
 SPB(4,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*
 pow(SPA(4,5),2)*pow(SPB(5,4),3)*SPA(3,4))+
 (complex<T>(0,-1)*SPA(1,3)*SPA(2,4)*SPA(3,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPA(1,2),2)*
 pow(SPB(2,1),2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(4,5),2)*pow(SPB(5,4),2)*SPA(1,3)*SPA(2,4)*
 SPA(3,5)*SPB(3,2)*SPB(4,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),
3)*pow(SPA(1,2),4)*pow(SPB(2,1),4)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,-1)*SPA(1,3)*SPA(2,4)*SPA(3,5)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*S(1,2)*S(4,5)*SPA(2,3)*
 SPA(3,4))+(complex<T>(0,1)*S(4,5)*SPA(1,3)*SPA(2,4)*SPA(3,5)*SPB(3,2)*
 SPB(4,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*
 pow(SPA(1,2),3)*pow(SPB(2,1),3)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,1)*SPA(1,3)*SPB(3,2)*SPB(4,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(4,5),2)*pow(SPB(5,4),2)*SPA(3,4)*
 SPA(4,5))+(complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*SPB(2,1))/
(complex<T>(6,0)*pow(SPA(4,5),2)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(2,4),2)*SPA(1,5)*SPB(4,2))/
(complex<T>(6,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPB(5,1))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*SPA(2,4)*SPB(3,1)*SPB(5,3))/
(complex<T>(6,0)*pow(SPA(1,2),2)*pow(SPA(4,5),2)*pow(SPB(5,1),2))+
 (complex<T>(0,1)*SPA(1,5)*SPB(3,1)*SPB(4,2)*SPB(5,3))/
(complex<T>(6,0)*S(1,2)*S(4,5)*SPB(5,1))+
 (complex<T>(0,1)*SPA(1,3)*SPA(3,5)*SPB(3,2)*SPB(4,3)*SPB(5,1))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(2,3)*SPA(3,4)*SPB(5,4))+(complex<T>(0,-1)*pow(SPA(1,5),2)*
 pow(SPA(2,4),3)*SPB(2,1)*SPB(5,4))/(complex<T>(6,0)*pow(SPA(1,2),2)*
 pow(SPA(4,5),2)*pow(SPB(5,1),2)*SPA(2,3)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPA(1,5),2)*pow(SPA(2,4),2)*SPB(5,4))/
(complex<T>(6,0)*pow(SPA(1,2),2)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPB(5,1))+
 (complex<T>(0,-1)*pow(SPA(4,5),2)*SPA(1,3)*SPA(3,5)*SPB(3,2)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(4,5)/S(1,2),3)*
 pow(SPA(1,2),3)*pow(SPB(2,1),4)*SPA(2,3)*SPA(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmmQp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, m, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmmQp nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPB(5,2),2))/(complex<T>(9,0)*SPB(2,1)*SPB(4,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmpQp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, p, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmpQp nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,3),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpmQm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, m, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpmQm nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPB(3,2),2))/(complex<T>(9,0)*SPB(2,1)*SPB(4,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQppQm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, p, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQppQm nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,5),2))/(complex<T>(9,0)*SPA(1,2)*SPA(3,4)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmQpm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, m}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmQpm nfLLT");
#endif
 
 return( (complex<T>(0,1)*pow(SPA(3,5),2)*S(1,2)*SPA(1,5)*SPB(5,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),4)*
 pow(SPB(4,3),3))+(complex<T>(0,-1)*pow(SPA(3,5),2)*SPA(1,5)*SPB(5,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*S(1,2)*
 SPB(4,3))+(complex<T>(0,-1)*SPA(1,3)*SPA(3,5)*SPB(4,2))/
(complex<T>(3,0)*S(1,2)*SPA(3,4)*SPB(5,4))+
 (complex<T>(0,-2)*pow(SPB(4,2),2)*SPB(4,1))/(complex<T>(9,0)*SPB(2,1)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,-1)*pow(SPA(3,5),2)*S(1,2)*SPB(4,2)*
 SPB(5,2))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*
 pow(SPA(3,4),3)*pow(SPB(4,3),2)*SPB(2,1)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(3,5),2)*SPB(4,2)*SPB(5,2))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*SPB(2,1)*
 SPB(4,3)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQmQpp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, p}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQmQpp nfLLT");
#endif
 
 return( (complex<T>(0,-1)*pow(SPB(5,2),2)*S(3,4)*SPA(1,3)*SPA(3,5))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*
 pow(SPB(2,1),3)*SPA(1,5)*SPA(3,4))+
 (complex<T>(0,1)*pow(SPB(5,2),2)*SPA(1,3)*SPA(3,5))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPB(2,1),2)*SPA(1,2)*
 SPA(1,5)*SPA(3,4))+(complex<T>(0,-2)*pow(SPA(1,3),2)*SPA(1,4))/
(complex<T>(9,0)*SPA(1,2)*SPA(1,5)*SPA(3,4)*SPA(4,5))+
 (complex<T>(0,-1)*SPA(1,3)*SPB(4,2)*SPB(5,2))/(complex<T>(3,0)*S(3,4)*SPA(1,5)*
 SPB(2,1))+(complex<T>(0,1)*pow(SPB(5,2),2)*S(3,4)*SPA(3,5)*SPB(5,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPA(1,2),3)*
 pow(SPB(2,1),4))+(complex<T>(0,-1)*pow(SPB(5,2),2)*SPA(3,5)*SPB(5,4))/
(complex<T>(3,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),3)*pow(SPB(2,1),2)*S(3,4)*
 SPA(1,2))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpQmm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, m}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpQmm nfLLT");
#endif
 
 return( (complex<T>(0,1)*SPA(1,4)*SPA(2,5)*SPA(3,5)*SPB(3,2))/
(complex<T>(6,0)*S(1,2)*S(3,4)*SPA(2,3))+
 (complex<T>(0,1)*pow(SPB(3,2),2)*SPA(2,5)*SPA(3,5)*SPB(4,1))/
(complex<T>(6,0)*pow(SPA(2,3),2)*pow(SPB(2,1),2)*pow(SPB(4,3),2))+
 (complex<T>(0,-1)*pow(SPB(3,2),2)*pow(SPB(4,1),3)*SPA(1,2)*SPA(3,4))/
(complex<T>(6,0)*pow(SPA(2,3),2)*pow(SPB(2,1),2)*pow(SPB(4,3),2)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*pow(SPB(3,2),2)*pow(SPB(4,1),2)*
 SPA(1,2))/(complex<T>(6,0)*pow(SPB(4,3),2)*SPA(2,3)*SPB(2,1)*SPB(5,1)*
 SPB(5,4))+(complex<T>(0,1)*pow(SPB(3,2),2)*pow(SPB(4,1),2)*SPA(3,4))/
(complex<T>(6,0)*pow(SPB(2,1),2)*SPA(2,3)*SPB(4,3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*pow(SPB(4,1),2)*SPA(1,4)*SPB(3,2))/
(complex<T>(6,0)*SPA(2,3)*SPB(2,1)*SPB(4,3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-2)*pow(SPB(3,2),2)*SPB(4,1))/(complex<T>(9,0)*SPB(2,1)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,-1)*S(1,2)*S(1,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*S(1,2)*S(4,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),3)*
 pow(SPB(4,3),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*S(1,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,-1)*S(4,5)*SPA(4,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),2)*pow(SPA(3,4),2)*SPB(4,3)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*S(1,5)*S(3,4)*SPA(1,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*S(3,4)*S(4,5)*SPA(1,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),3)*
 pow(SPB(2,1),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*S(1,5)*SPA(1,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*SPB(2,1)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,1)*S(4,5)*SPA(1,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(3,4)/S(1,2),2)*pow(SPA(1,2),2)*SPB(2,1)*
 SPB(5,1)*SPB(5,4))+(complex<T>(0,-1)*S(1,2)*SPA(1,5)*SPA(2,3)*SPA(4,5)*
 SPB(2,1)*SPB(5,2)*SPB(5,3))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),
3)*pow(SPA(3,4),4)*pow(SPB(4,3),3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*SPA(1,5)*SPA(4,5)*SPB(4,1)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*pow(SPA(1,2),2)*pow(SPB(2,1),2)*SPA(1,5)*SPA(4,5)*
 SPB(4,1)*SPB(5,2)*SPB(5,3))/(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),
3)*pow(SPA(3,4),4)*pow(SPB(4,3),4)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,1)*S(1,2)*SPA(1,5)*SPA(4,5)*SPB(4,1)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),3)*
 pow(SPB(4,3),3)*SPB(5,1)*SPB(5,4))+
 (complex<T>(0,-1)*SPA(1,5)*SPA(4,5)*SPB(4,1)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*S(1,2)*S(3,4)*SPB(5,1)*
 SPB(5,4))+(complex<T>(0,1)*SPA(1,5)*SPA(2,3)*SPA(4,5)*SPB(2,1)*
 SPB(5,2)*SPB(5,3))/(complex<T>(3,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*
 pow(SPA(3,4),2)*S(1,2)*SPB(4,3)*SPB(5,1)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmqpQpQmp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, p}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmqpQpQmp nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,4),3))/(complex<T>(9,0)*SPA(1,2)*SPA(1,5)*SPA(3,4)*
 SPA(4,5))+(complex<T>(0,-1)*pow(SPB(5,3),2)*S(1,2)*SPA(1,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),3)*
 pow(SPB(4,3),4))+(complex<T>(0,1)*pow(SPB(5,3),2)*SPA(1,5)*SPB(5,2))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPB(4,3),2)*S(1,2)*
 SPA(3,4))+(complex<T>(0,1)*pow(SPB(5,2),2)*SPA(4,5)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),2)*
 pow(SPB(4,3),2)*SPB(2,1))+(complex<T>(0,-1)*pow(SPA(1,2),2)*
 pow(SPB(5,2),2)*SPA(4,5)*SPB(2,1)*SPB(5,3))/
(complex<T>(6,0)*pow(complex<T>(1,0)-S(1,2)/S(3,4),3)*pow(SPA(3,4),4)*
 pow(SPB(4,3),4))+(complex<T>(0,1)*SPA(1,4)*SPB(5,2)*SPB(5,3))/
(complex<T>(6,0)*S(1,2)*S(3,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmmqpQmQp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmmqpQmQp nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPB(5,3),2))/(complex<T>(9,0)*SPB(2,1)*SPB(3,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmpqpQmQp_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, Qm, Qp}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmpqpQmQp nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,4),2))/(complex<T>(9,0)*SPA(1,2)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmmqpQpQm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, m, qp, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmmqpQpQm nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPB(4,3),2))/(complex<T>(9,0)*SPB(2,1)*SPB(3,2)*SPB(5,4))
        ); 
  
} 
  
  
 
 
template <class T> complex<T> R2q2Q1g_qmpqpQpQm_nfLLT
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, p, qp, Qp, Qm}, nfLLT}
 
#if _VERBOSE
  _MESSAGE("R2q2Q1g :  qmpqpQpQm nfLLT");
#endif
 
 return( (complex<T>(0,-2)*pow(SPA(1,5),2))/(complex<T>(9,0)*SPA(1,2)*SPA(2,3)*SPA(4,5))
        ); 
  
} 
  
  
 
 // *************** table of switch values ************* 
 
#define _R_qmqpmQmQp_LLT R2q2Q1g_7357_LLT
#define _R_qmqppQmQp_LLT R2q2Q1g_7465_LLT
#define _R_qmqpmQpQm_LLT R2q2Q1g_6277_LLT
#define _R_qmqppQpQm_LLT R2q2Q1g_6385_LLT
#define _R_qmqpQmmQp_LLT R2q2Q1g_6637_LLT
#define _R_qmqpQmpQp_LLT R2q2Q1g_7285_LLT
#define _R_qmqpQpmQm_LLT R2q2Q1g_5377_LLT
#define _R_qmqpQppQm_LLT R2q2Q1g_6025_LLT
#define _R_qmqpQmQpm_LLT R2q2Q1g_1237_LLT
#define _R_qmqpQmQpp_LLT R2q2Q1g_5125_LLT
#define _R_qmqpQpQmm_LLT R2q2Q1g_1057_LLT
#define _R_qmqpQpQmp_LLT R2q2Q1g_4945_LLT
#define _R_qmmqpQmQp_LLT R2q2Q1g_7417_LLT
#define _R_qmpqpQmQp_LLT R2q2Q1g_7435_LLT
#define _R_qmmqpQpQm_LLT R2q2Q1g_6337_LLT
#define _R_qmpqpQpQm_LLT R2q2Q1g_6355_LLT
#define _R_qmqpmQmQp_LRT R2q2Q1g_7357_LRT
#define _R_qmqppQmQp_LRT R2q2Q1g_7465_LRT
#define _R_qmqpmQpQm_LRT R2q2Q1g_6277_LRT
#define _R_qmqppQpQm_LRT R2q2Q1g_6385_LRT
#define _R_qmqpQmmQp_LRT R2q2Q1g_6637_LRT
#define _R_qmqpQmpQp_LRT R2q2Q1g_7285_LRT
#define _R_qmqpQpmQm_LRT R2q2Q1g_5377_LRT
#define _R_qmqpQppQm_LRT R2q2Q1g_6025_LRT
#define _R_qmqpQmQpm_LRT R2q2Q1g_1237_LRT
#define _R_qmqpQmQpp_LRT R2q2Q1g_5125_LRT
#define _R_qmqpQpQmm_LRT R2q2Q1g_1057_LRT
#define _R_qmqpQpQmp_LRT R2q2Q1g_4945_LRT
#define _R_qmmqpQmQp_LRT R2q2Q1g_7417_LRT
#define _R_qmpqpQmQp_LRT R2q2Q1g_7435_LRT
#define _R_qmmqpQpQm_LRT R2q2Q1g_6337_LRT
#define _R_qmpqpQpQm_LRT R2q2Q1g_6355_LRT
#define _R_qmmQmQpqp_LLT R2q2Q1g_3817_LLT
#define _R_qmmQpQmqp_LLT R2q2Q1g_3637_LLT
#define _R_qmpQmQpqp_LLT R2q2Q1g_3835_LLT
#define _R_qmpQpQmqp_LLT R2q2Q1g_3655_LLT
#define _R_qmQmmQpqp_LLT R2q2Q1g_3697_LLT
#define _R_qmQpmQmqp_LLT R2q2Q1g_3487_LLT
#define _R_qmQmpQpqp_LLT R2q2Q1g_3805_LLT
#define _R_qmQppQmqp_LLT R2q2Q1g_3595_LLT
#define _R_qmQmQpqpm_LLT R2q2Q1g_637_LLT
#define _R_qmQmQpqpp_LLT R2q2Q1g_4525_LLT
#define _R_qmQpQmqpm_LLT R2q2Q1g_607_LLT
#define _R_qmQpQmqpp_LLT R2q2Q1g_4495_LLT
#define _R_qmQmQpmqp_LLT R2q2Q1g_2797_LLT
#define _R_qmQmQppqp_LLT R2q2Q1g_3445_LLT
#define _R_qmQpQmmqp_LLT R2q2Q1g_2767_LLT
#define _R_qmQpQmpqp_LLT R2q2Q1g_3415_LLT
#define _R_qmqpmQmQp_nfLLT R2q2Q1g_7357_nfLLT
#define _R_qmqppQmQp_nfLLT R2q2Q1g_7465_nfLLT
#define _R_qmqpmQpQm_nfLLT R2q2Q1g_6277_nfLLT
#define _R_qmqppQpQm_nfLLT R2q2Q1g_6385_nfLLT
#define _R_qmqpQmmQp_nfLLT R2q2Q1g_6637_nfLLT
#define _R_qmqpQmpQp_nfLLT R2q2Q1g_7285_nfLLT
#define _R_qmqpQpmQm_nfLLT R2q2Q1g_5377_nfLLT
#define _R_qmqpQppQm_nfLLT R2q2Q1g_6025_nfLLT
#define _R_qmqpQmQpm_nfLLT R2q2Q1g_1237_nfLLT
#define _R_qmqpQmQpp_nfLLT R2q2Q1g_5125_nfLLT
#define _R_qmqpQpQmm_nfLLT R2q2Q1g_1057_nfLLT
#define _R_qmqpQpQmp_nfLLT R2q2Q1g_4945_nfLLT
#define _R_qmmqpQmQp_nfLLT R2q2Q1g_7417_nfLLT
#define _R_qmpqpQmQp_nfLLT R2q2Q1g_7435_nfLLT
#define _R_qmmqpQpQm_nfLLT R2q2Q1g_6337_nfLLT
#define _R_qmpqpQpQm_nfLLT R2q2Q1g_6355_nfLLT
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qmqpmQmQp_LLT case 7357 : \
          return &R2q2Q1g_7357_LLT
#define _CASE_qmqppQmQp_LLT case 7465 : \
          return &R2q2Q1g_7465_LLT
#define _CASE_qmqpmQpQm_LLT case 6277 : \
          return &R2q2Q1g_6277_LLT
#define _CASE_qmqppQpQm_LLT case 6385 : \
          return &R2q2Q1g_6385_LLT
#define _CASE_qmqpQmmQp_LLT case 6637 : \
          return &R2q2Q1g_6637_LLT
#define _CASE_qmqpQmpQp_LLT case 7285 : \
          return &R2q2Q1g_7285_LLT
#define _CASE_qmqpQpmQm_LLT case 5377 : \
          return &R2q2Q1g_5377_LLT
#define _CASE_qmqpQppQm_LLT case 6025 : \
          return &R2q2Q1g_6025_LLT
#define _CASE_qmqpQmQpm_LLT case 1237 : \
          return &R2q2Q1g_1237_LLT
#define _CASE_qmqpQmQpp_LLT case 5125 : \
          return &R2q2Q1g_5125_LLT
#define _CASE_qmqpQpQmm_LLT case 1057 : \
          return &R2q2Q1g_1057_LLT
#define _CASE_qmqpQpQmp_LLT case 4945 : \
          return &R2q2Q1g_4945_LLT
#define _CASE_qmmqpQmQp_LLT case 7417 : \
          return &R2q2Q1g_7417_LLT
#define _CASE_qmpqpQmQp_LLT case 7435 : \
          return &R2q2Q1g_7435_LLT
#define _CASE_qmmqpQpQm_LLT case 6337 : \
          return &R2q2Q1g_6337_LLT
#define _CASE_qmpqpQpQm_LLT case 6355 : \
          return &R2q2Q1g_6355_LLT
#define _CASE_qmqpmQmQp_LRT case 7357 : \
          return &R2q2Q1g_7357_LRT
#define _CASE_qmqppQmQp_LRT case 7465 : \
          return &R2q2Q1g_7465_LRT
#define _CASE_qmqpmQpQm_LRT case 6277 : \
          return &R2q2Q1g_6277_LRT
#define _CASE_qmqppQpQm_LRT case 6385 : \
          return &R2q2Q1g_6385_LRT
#define _CASE_qmqpQmmQp_LRT case 6637 : \
          return &R2q2Q1g_6637_LRT
#define _CASE_qmqpQmpQp_LRT case 7285 : \
          return &R2q2Q1g_7285_LRT
#define _CASE_qmqpQpmQm_LRT case 5377 : \
          return &R2q2Q1g_5377_LRT
#define _CASE_qmqpQppQm_LRT case 6025 : \
          return &R2q2Q1g_6025_LRT
#define _CASE_qmqpQmQpm_LRT case 1237 : \
          return &R2q2Q1g_1237_LRT
#define _CASE_qmqpQmQpp_LRT case 5125 : \
          return &R2q2Q1g_5125_LRT
#define _CASE_qmqpQpQmm_LRT case 1057 : \
          return &R2q2Q1g_1057_LRT
#define _CASE_qmqpQpQmp_LRT case 4945 : \
          return &R2q2Q1g_4945_LRT
#define _CASE_qmmqpQmQp_LRT case 7417 : \
          return &R2q2Q1g_7417_LRT
#define _CASE_qmpqpQmQp_LRT case 7435 : \
          return &R2q2Q1g_7435_LRT
#define _CASE_qmmqpQpQm_LRT case 6337 : \
          return &R2q2Q1g_6337_LRT
#define _CASE_qmpqpQpQm_LRT case 6355 : \
          return &R2q2Q1g_6355_LRT
#define _CASE_qmmQmQpqp_LLT case 3817 : \
          return &R2q2Q1g_3817_LLT
#define _CASE_qmmQpQmqp_LLT case 3637 : \
          return &R2q2Q1g_3637_LLT
#define _CASE_qmpQmQpqp_LLT case 3835 : \
          return &R2q2Q1g_3835_LLT
#define _CASE_qmpQpQmqp_LLT case 3655 : \
          return &R2q2Q1g_3655_LLT
#define _CASE_qmQmmQpqp_LLT case 3697 : \
          return &R2q2Q1g_3697_LLT
#define _CASE_qmQpmQmqp_LLT case 3487 : \
          return &R2q2Q1g_3487_LLT
#define _CASE_qmQmpQpqp_LLT case 3805 : \
          return &R2q2Q1g_3805_LLT
#define _CASE_qmQppQmqp_LLT case 3595 : \
          return &R2q2Q1g_3595_LLT
#define _CASE_qmQmQpqpm_LLT case 637 : \
          return &R2q2Q1g_637_LLT
#define _CASE_qmQmQpqpp_LLT case 4525 : \
          return &R2q2Q1g_4525_LLT
#define _CASE_qmQpQmqpm_LLT case 607 : \
          return &R2q2Q1g_607_LLT
#define _CASE_qmQpQmqpp_LLT case 4495 : \
          return &R2q2Q1g_4495_LLT
#define _CASE_qmQmQpmqp_LLT case 2797 : \
          return &R2q2Q1g_2797_LLT
#define _CASE_qmQmQppqp_LLT case 3445 : \
          return &R2q2Q1g_3445_LLT
#define _CASE_qmQpQmmqp_LLT case 2767 : \
          return &R2q2Q1g_2767_LLT
#define _CASE_qmQpQmpqp_LLT case 3415 : \
          return &R2q2Q1g_3415_LLT
#define _CASE_qmqpmQmQp_nfLLT case 7357 : \
          return &R2q2Q1g_7357_nfLLT
#define _CASE_qmqppQmQp_nfLLT case 7465 : \
          return &R2q2Q1g_7465_nfLLT
#define _CASE_qmqpmQpQm_nfLLT case 6277 : \
          return &R2q2Q1g_6277_nfLLT
#define _CASE_qmqppQpQm_nfLLT case 6385 : \
          return &R2q2Q1g_6385_nfLLT
#define _CASE_qmqpQmmQp_nfLLT case 6637 : \
          return &R2q2Q1g_6637_nfLLT
#define _CASE_qmqpQmpQp_nfLLT case 7285 : \
          return &R2q2Q1g_7285_nfLLT
#define _CASE_qmqpQpmQm_nfLLT case 5377 : \
          return &R2q2Q1g_5377_nfLLT
#define _CASE_qmqpQppQm_nfLLT case 6025 : \
          return &R2q2Q1g_6025_nfLLT
#define _CASE_qmqpQmQpm_nfLLT case 1237 : \
          return &R2q2Q1g_1237_nfLLT
#define _CASE_qmqpQmQpp_nfLLT case 5125 : \
          return &R2q2Q1g_5125_nfLLT
#define _CASE_qmqpQpQmm_nfLLT case 1057 : \
          return &R2q2Q1g_1057_nfLLT
#define _CASE_qmqpQpQmp_nfLLT case 4945 : \
          return &R2q2Q1g_4945_nfLLT
#define _CASE_qmmqpQmQp_nfLLT case 7417 : \
          return &R2q2Q1g_7417_nfLLT
#define _CASE_qmpqpQmQp_nfLLT case 7435 : \
          return &R2q2Q1g_7435_nfLLT
#define _CASE_qmmqpQpQm_nfLLT case 6337 : \
          return &R2q2Q1g_6337_nfLLT
#define _CASE_qmpqpQpQm_nfLLT case 6355 : \
          return &R2q2Q1g_6355_nfLLT
 
 
 // *************** function definitions using macros ************* 
 
template <class T> complex<T> _R_qmqpmQmQp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpmQmQp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqppQmQp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqppQmQp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmQpQm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpmQpQm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqppQpQm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqppQpQm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmmQp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmmQp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmpQp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmpQp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpmQm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpmQm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQppQm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQppQm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQpm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmQpm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQpp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmQpp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQmm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpQmm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQmp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpQmp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpQmQp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmmqpQmQp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpQmQp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmpqpQmQp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpQpQm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmmqpQpQm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpQpQm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmpqpQpQm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmQmQp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpmQmQp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqppQmQp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqppQmQp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmQpQm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpmQpQm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqppQpQm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqppQpQm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmmQp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmmQp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmpQp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmpQp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpmQm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpmQm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQppQm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQppQm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQpm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmQpm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQpp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmQpp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQmm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpQmm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQmp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpQmp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpQmQp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmmqpQmQp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpQmQp_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmpqpQmQp_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpQpQm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmmqpQpQm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpQpQm_LRT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmpqpQpQm_LRT(ep,mpc);}
 
template <class T> complex<T> _R_qmmQmQpqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmmQmQpqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmmQpQmqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmmQpQmqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmpQmQpqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmpQmQpqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmpQpQmqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmpQpQmqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQmmQpqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQmmQpqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQpmQmqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQpmQmqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQmpQpqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQmpQpqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQppQmqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQppQmqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpqpm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQmQpqpm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpqpp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQmQpqpp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQpQmqpm_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQpQmqpm_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQpQmqpp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQpQmqpp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQpmqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQmQpmqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQmQppqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQmQppqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQpQmmqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQpQmmqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmQpQmpqp_LLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmQpQmpqp_LLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmQmQp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpmQmQp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqppQmQp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqppQmQp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpmQpQm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpmQpQm_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqppQpQm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqppQpQm_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmmQp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmmQp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmpQp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmpQp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpmQm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpmQm_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQppQm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQppQm_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQpm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmQpm_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQmQpp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQmQpp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQmm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpQmm_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmqpQpQmp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmqpQpQmp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpQmQp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmmqpQmQp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpQmQp_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmpqpQmQp_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmmqpQpQm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmmqpQpQm_nfLLT(ep,mpc);}
 
template <class T> complex<T> _R_qmpqpQpQm_nfLLT(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g_qmpqpQpQm_nfLLT(ep,mpc);}
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> complex<T> ( *R2q2Q1g_LLT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmqpmQmQp_LLT;
       _CASE_qmqppQmQp_LLT;
       _CASE_qmqpmQpQm_LLT;
       _CASE_qmqppQpQm_LLT;
       _CASE_qmqpQmmQp_LLT;
       _CASE_qmqpQmpQp_LLT;
       _CASE_qmqpQpmQm_LLT;
       _CASE_qmqpQppQm_LLT;
       _CASE_qmqpQmQpm_LLT;
       _CASE_qmqpQmQpp_LLT;
       _CASE_qmqpQpQmm_LLT;
       _CASE_qmqpQpQmp_LLT;
       _CASE_qmmqpQmQp_LLT;
       _CASE_qmpqpQmQp_LLT;
       _CASE_qmmqpQpQm_LLT;
       _CASE_qmpqpQpQm_LLT;
       _CASE_qmmQmQpqp_LLT;
       _CASE_qmmQpQmqp_LLT;
       _CASE_qmpQmQpqp_LLT;
       _CASE_qmpQpQmqp_LLT;
       _CASE_qmQmmQpqp_LLT;
       _CASE_qmQpmQmqp_LLT;
       _CASE_qmQmpQpqp_LLT;
       _CASE_qmQppQmqp_LLT;
       _CASE_qmQmQpqpm_LLT;
       _CASE_qmQmQpqpp_LLT;
       _CASE_qmQpQmqpm_LLT;
       _CASE_qmQpQmqpp_LLT;
       _CASE_qmQmQpmqp_LLT;
       _CASE_qmQmQppqp_LLT;
       _CASE_qmQpQmmqp_LLT;
       _CASE_qmQpQmpqp_LLT;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q1g_LRT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmqpmQmQp_LRT;
       _CASE_qmqppQmQp_LRT;
       _CASE_qmqpmQpQm_LRT;
       _CASE_qmqppQpQm_LRT;
       _CASE_qmqpQmmQp_LRT;
       _CASE_qmqpQmpQp_LRT;
       _CASE_qmqpQpmQm_LRT;
       _CASE_qmqpQppQm_LRT;
       _CASE_qmqpQmQpm_LRT;
       _CASE_qmqpQmQpp_LRT;
       _CASE_qmqpQpQmm_LRT;
       _CASE_qmqpQpQmp_LRT;
       _CASE_qmmqpQmQp_LRT;
       _CASE_qmpqpQmQp_LRT;
       _CASE_qmmqpQpQm_LRT;
       _CASE_qmpqpQpQm_LRT;
 
       default: return 0;
        }
 }
 
template <class T> complex<T> ( *R2q2Q1g_nfLLT_Ptr_eval( int hc))
      (  const eval_param<T>& ep,  const mass_param_coll& mpc){
       switch (hc) {
       _CASE_qmqpmQmQp_nfLLT;
       _CASE_qmqppQmQp_nfLLT;
       _CASE_qmqpmQpQm_nfLLT;
       _CASE_qmqppQpQm_nfLLT;
       _CASE_qmqpQmmQp_nfLLT;
       _CASE_qmqpQmpQp_nfLLT;
       _CASE_qmqpQpmQm_nfLLT;
       _CASE_qmqpQppQm_nfLLT;
       _CASE_qmqpQmQpm_nfLLT;
       _CASE_qmqpQmQpp_nfLLT;
       _CASE_qmqpQpQmm_nfLLT;
       _CASE_qmqpQpQmp_nfLLT;
       _CASE_qmmqpQmQp_nfLLT;
       _CASE_qmpqpQmQp_nfLLT;
       _CASE_qmmqpQpQm_nfLLT;
       _CASE_qmpqpQpQm_nfLLT;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template complex<R> ( *R2q2Q1g_LLT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q1g_LLT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1g_LLT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2Q1g_LRT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q1g_LRT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1g_LRT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);


template complex<R> ( *R2q2Q1g_nfLLT_Ptr_eval(int hc))
             (const eval_param<R>&,  const mass_param_coll&);
template complex<RHP> ( *R2q2Q1g_nfLLT_Ptr_eval(int hc))
             (const eval_param<RHP>&,  const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1g_nfLLT_Ptr_eval(int hc))
             (const eval_param<RVHP>&,  const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q1g_LLT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2Q1g_LRT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);
template complex<RGMP> ( *R2q2Q1g_nfLLT_Ptr_eval(int hc))
             (const eval_param<RGMP>&,  const mass_param_coll&);

#endif


}
