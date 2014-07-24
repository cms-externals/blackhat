/*
* R_2q2Q2l.cpp
*
* Created on 11/16, 2008
*      Author: Zvi's script
*/

#include "R_2q2Q2l_eval.h"
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



template <class T> complex<T> R2q2Q2l_qpQmQpqmemep_nf_top
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qm, Qp, qm, em, ep}, nf_top}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpQmQpqmemep nf_top");
#endif

complex<T> mt = complex<T>(constants::s_GeV,0)*complex<T>(constants::Mtop,0);
complex<T> mtsq = pow(mt, 2);
complex<T> top = (complex<T>(9,0)/complex<T>(-2,0)*complex<T>(-2,0)/complex<T>(15,0)*S(2,3)/mtsq);

 return( top*(complex<T>(0,-2)*((SPA(4,5)*SPB(3,1)*(SPA(1,2)*SPB(6,1)-
  SPA(2,3)*SPB(6,3)))/(S(2,3)*S(5,6)*SS(1,2,3))+
 (SPA(2,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*SPB(6,1))/
(S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}




template <class T> complex<T> R2q2Q2l_qpQpQmqmemep_nf_top
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qp, Qm, qm, em, ep}, nf_top}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpQpQmqmemep nf_top");
#endif

complex<T> mt = complex<T>(BH::constants::s_GeV,0)*complex<T>(BH::constants::Mtop,0);
complex<T> mtsq = pow(mt, 2);
complex<T> top = (complex<T>(9,0)/complex<T>(-2,0)*complex<T>(-2,0)/complex<T>(15,0)*S(2,3)/mtsq);

 return( top*(complex<T>(0,2)*((SPA(4,5)*SPB(2,1)*(SPA(1,3)*SPB(6,1)+
  SPA(2,3)*SPB(6,2)))/(S(2,3)*S(5,6)*SS(1,2,3))+
 (SPA(3,4)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(6,1))/
(S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}




template <class T> complex<T> R2q2Q2l_qmQpQmqpemep_nf_top
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, em, ep}, nf_top}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmQpQmqpemep nf_top");
#endif

complex<T> mt = complex<T>(BH::constants::s_GeV,0)*complex<T>(BH::constants::Mtop,0);
complex<T> mtsq = pow(mt, 2);
complex<T> top = (complex<T>(9,0)/complex<T>(-2,0)*complex<T>(-2,0)/complex<T>(15,0)*S(2,3)/mtsq);

 return( top*(complex<T>(0,-2)*((SPA(1,3)*(SPA(1,5)*SPB(2,1)-SPA(3,5)*SPB(3,2))*
 SPB(6,4))/(S(2,3)*S(5,6)*SS(1,2,3))+
 (SPA(1,5)*SPB(4,2)*(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4)))/
(S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}


template <class T> complex<T> R2q2Q2l_qmQmQpqpemep_nf_top
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, em, ep}, nf_top}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmQmQpqpemep nf_top");
#endif

complex<T> mt = complex<T>(BH::constants::s_GeV,0)*complex<T>(BH::constants::Mtop,0);
complex<T> mtsq = pow(mt, 2);
complex<T> top = (complex<T>(9,0)/complex<T>(-2,0)*complex<T>(-2,0)/complex<T>(15,0)*S(2,3)/mtsq);

 return( top*(complex<T>(0,2)*((SPA(1,2)*(SPA(1,5)*SPB(3,1)+SPA(2,5)*SPB(3,2))*
 SPB(6,4))/(S(2,3)*S(5,6)*SS(1,2,3))+
 (SPA(1,5)*SPB(4,3)*(-(SPA(2,3)*SPB(6,3))-SPA(2,4)*SPB(6,4)))/
(S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}


template <class T> complex<T> R2q2Q2l_qpQmQpqmemep_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qm, Qp, qm, em, ep}, nf}
 
#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpQmQpqmemep nf");
#endif
 
 return( (complex<T>(0,-2)*((SPA(4,5)*SPB(3,1)*(SPA(1,2)*SPB(6,1)-
  SPA(2,3)*SPB(6,3)))/(S(2,3)*S(5,6)*SS(1,2,3))+
 (SPA(2,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*SPB(6,1))/
(S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        ); 
  
} 


template <class T> complex<T> R2q2Q2l_qpQpQmqmemep_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qp, Qm, qm, em, ep}, nf}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpQpQmqmemep nf");
#endif

 return( (complex<T>(0,2)*((SPA(4,5)*SPB(2,1)*(SPA(1,3)*SPB(6,1)+
  SPA(2,3)*SPB(6,2)))/(S(2,3)*S(5,6)*SS(1,2,3))+
 (SPA(3,4)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(6,1))/
(S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}




template <class T> complex<T> R2q2Q2l_qmQpQmqpemep_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, em, ep}, nf}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmQpQmqpemep nf");
#endif

 return( (complex<T>(0,-2)*((SPA(1,3)*(SPA(1,5)*SPB(2,1)-SPA(3,5)*SPB(3,2))*
 SPB(6,4))/(S(2,3)*S(5,6)*SS(1,2,3))+
 (SPA(1,5)*SPB(4,2)*(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4)))/
(S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}




template <class T> complex<T> R2q2Q2l_qmQmQpqpemep_nf
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, em, ep}, nf}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmQmQpqpemep nf");
#endif

 return( (complex<T>(0,2)*((SPA(1,2)*(SPA(1,5)*SPB(3,1)+SPA(2,5)*SPB(3,2))*
 SPB(6,4))/(S(2,3)*S(5,6)*SS(1,2,3))+
 (SPA(1,5)*SPB(4,3)*(-(SPA(2,3)*SPB(6,3))-SPA(2,4)*SPB(6,4)))/
(S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}




template <class T> complex<T> R2q2Q2l_qpQmQpqmemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qm, Qp, qm, em, ep}, L}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpQmQpqmemep L");
#endif

 return( complex<T>(0,1)*((complex<T>(1,0)*pow(SPA(2,4),2)*pow(SPB(6,4),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1))/(complex<T>(2,0)*SPA(2,3)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SPB(6,5)*
 SS(1,2,3))+(complex<T>(1,0)*pow(SPB(3,1),2)*
 pow(-(SPA(1,2)*SPB(6,2))-SPA(1,3)*SPB(6,3),2)*
 pow(complex<T>(1,0)-SS(1,2,3)/S(2,3),-1))/(complex<T>(2,0)*pow(SPB(3,2),2)*
 SPA(2,3)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SPB(6,5)*
 SS(1,2,3))-(pow(SPA(1,5),2)*pow(SPB(3,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1))/(complex<T>(2,0)*SPA(5,6)*SPB(3,2)*
 (-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*SS(2,3,4))-
 (pow(SPA(2,4),2)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/(complex<T>(2,0)*pow(SPA(2,3),2)*
 SPA(5,6)*SPB(3,2)*(-(SPA(1,2)*SPB(4,2))-SPA(1,3)*SPB(4,3))*
 SS(2,3,4)))+
 (complex<T>(0,2)*((SPA(4,5)*SPB(3,1)*(SPA(1,2)*SPB(6,1)-
   SPA(2,3)*SPB(6,3)))/(S(2,3)*S(5,6)*SS(1,2,3))+
(SPA(2,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*SPB(6,1))/
 (S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}




template <class T> complex<T> R2q2Q2l_qpQpQmqmemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, Qp, Qm, qm, em, ep}, L}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpQpQmqmemep L");
#endif

 return( complex<T>(0,1)*((complex<T>(1,0)*(pow(SPA(3,4),2)*pow(SPB(6,4),2)*
   pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)+
  (pow(SPB(2,1),2)*pow(-(SPA(1,2)*SPB(6,2))-SPA(1,3)*SPB(6,
  3),2)*pow(complex<T>(1,0)-SS(1,2,3)/S(2,3),-1))/
   pow(SPB(3,2),2)))/(complex<T>(2,0)*SPA(2,3)*(-(SPA(1,2)*SPB(4,2))-
  SPA(1,3)*SPB(4,3))*SPB(6,5)*SS(1,2,3))-
 (pow(SPA(1,5),2)*pow(SPB(2,1),2)*
  pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)+
 (pow(SPA(3,4),2)*pow(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3),2)*
   pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/pow(SPA(2,3),2))/
(complex<T>(2,0)*SPA(5,6)*SPB(3,2)*(-(SPA(1,2)*SPB(4,2))-
  SPA(1,3)*SPB(4,3))*SS(2,3,4)))+
 (complex<T>(0,-2)*((SPA(4,5)*SPB(2,1)*(SPA(1,3)*SPB(6,1)+
   SPA(2,3)*SPB(6,2)))/(S(2,3)*S(5,6)*SS(1,2,3))+
(SPA(3,4)*(-(SPA(3,5)*SPB(3,2))-SPA(4,5)*SPB(4,2))*SPB(6,1))/
 (S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)
        );

}




template <class T> complex<T> R2q2Q2l_qmQpQmqpemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qp, Qm, qp, em, ep}, L}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmQpQmqpemep L");
#endif

 return( (complex<T>(0,2)*((SPA(1,3)*(SPA(1,5)*SPB(2,1)-SPA(3,5)*SPB(3,2))*
  SPB(6,4))/(S(2,3)*S(5,6)*SS(1,2,3))+
(SPA(1,5)*SPB(4,2)*(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4)))/
 (S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)+
 complex<T>(0,1)*(-((pow(SPA(4,5),2)*pow(SPB(4,2),2)*
  pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1))/(complex<T>(2,0)*SPA(5,6)*
  (-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SPB(3,2)*
  SS(1,2,3)))-(pow(SPA(1,3),2)*pow(-(SPA(2,5)*SPB(2,1))-
   SPA(3,5)*SPB(3,1),2)*pow(complex<T>(1,0)-SS(1,2,3)/S(2,3),-1))/
(complex<T>(2,0)*pow(SPA(2,3),2)*SPA(5,6)*(-(SPA(2,4)*SPB(2,1))-
  SPA(3,4)*SPB(3,1))*SPB(3,2)*SS(1,2,3))+
 (complex<T>(1,0)*pow(SPA(1,3),2)*pow(SPB(6,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1))/(complex<T>(2,0)*SPA(2,3)*
 (-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SPB(6,5)*
 SS(2,3,4))+(complex<T>(1,0)*pow(SPB(4,2),2)*
 pow(SPA(2,4)*SPB(6,2)+SPA(3,4)*SPB(6,3),2)*
 pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/(complex<T>(2,0)*pow(SPB(3,2),2)*
 SPA(2,3)*(-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*SPB(6,5)*
 SS(2,3,4)))
        );

}




template <class T> complex<T> R2q2Q2l_qmQmQpqpemep_L
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, Qm, Qp, qp, em, ep}, L}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmQmQpqpemep L");
#endif

 return( (complex<T>(0,-2)*((SPA(1,2)*(SPA(1,5)*SPB(3,1)+SPA(2,5)*SPB(3,2))*
  SPB(6,4))/(S(2,3)*S(5,6)*SS(1,2,3))+
(SPA(1,5)*SPB(4,3)*(-(SPA(2,3)*SPB(6,3))-SPA(2,4)*SPB(6,4)))/
 (S(2,3)*S(5,6)*SS(2,3,4))))/complex<T>(9,0)+
 complex<T>(0,1)*(-((pow(SPA(4,5),2)*pow(SPB(4,3),2)*
   pow(complex<T>(1,0)-S(5,6)/SS(1,2,3),-1)+
  (pow(SPA(1,2),2)*pow(-(SPA(2,5)*SPB(2,1))-SPA(3,5)*SPB(3,
  1),2)*pow(complex<T>(1,0)-SS(1,2,3)/S(2,3),-1))/
   pow(SPA(2,3),2))/(complex<T>(2,0)*SPA(5,6)*(-(SPA(2,4)*SPB(2,1))-
   SPA(3,4)*SPB(3,1))*SPB(3,2)*SS(1,2,3)))+
 (complex<T>(1,0)*(pow(SPA(1,2),2)*pow(SPB(6,1),2)*
   pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1)+
  (pow(SPB(4,3),2)*pow(SPA(2,4)*SPB(6,2)+SPA(3,4)*SPB(6,3),
 2)*pow(complex<T>(1,0)-SS(2,3,4)/S(2,3),-1))/pow(SPB(3,2),2)))/
(complex<T>(2,0)*SPA(2,3)*(-(SPA(2,4)*SPB(2,1))-SPA(3,4)*SPB(3,1))*
 SPB(6,5)*SS(2,3,4)))
        );

}




template <class T> complex<T> R2q2Q2l_qpqmQpQmemep_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, Qp, Qm, em, ep}, sl}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpqmQpQmemep sl");
#endif

 return( (complex<T>(0,-3)*(-((SPA(2,5)*SPB(3,1)*(SPA(1,4)*SPB(6,1)+
SPA(3,4)*SPB(6,3)))/(S(3,4)*S(5,6)*SS(1,3,4)))+
(SPA(2,4)*(SPA(2,5)*SPB(3,2)-SPA(4,5)*SPB(4,3))*SPB(6,1))/
 (S(3,4)*S(5,6)*SS(2,3,4))))/complex<T>(2,0)+
 complex<T>(0,1)*(-((pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,5)-S(1,6)+S(2,5)+
   S(2,6))*SPA(4,5)*(-((S(1,2)-S(3,4)-S(5,6))*SPA(4,5))+
   complex<T>(-2,0)*SPA(3,4)*SPA(5,6)*SPB(6,3)))/(complex<T>(4,0)*SPA(3,4)*
  SPA(5,6)*(SPA(1,5)*SPB(5,2)+SPA(1,6)*SPB(6,2))))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,3)-S(1,4)+S(2,3)+
  S(2,4))*SPA(4,5)*((S(1,2)-S(3,4)-S(5,6))*SPA(4,5)+
  complex<T>(2,0)*SPA(3,4)*SPA(5,6)*SPB(6,3)))/(complex<T>(4,0)*SPA(3,4)*
 SPA(5,6)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2)))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,3)+S(1,4)-S(2,3)-
  S(2,4))*SPB(6,3)*(-((S(1,2)-S(3,4)-S(5,6))*SPB(6,3))+
  complex<T>(-2,0)*SPA(4,5)*SPB(4,3)*SPB(6,5)))/
(complex<T>(4,0)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(4,3)*
 SPB(6,5))-(pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
   pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
   complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
 (S(1,5)+S(1,6)-S(2,5)-S(2,6))*SPB(6,3)*
 ((S(1,2)-S(3,4)-S(5,6))*SPB(6,3)+complex<T>(2,0)*SPA(4,5)*
   SPB(4,3)*SPB(6,5)))/(complex<T>(4,0)*SPB(4,3)*(SPA(1,5)*SPB(5,2)+
  SPA(1,6)*SPB(6,2))*SPB(6,5))-
 (pow(SPA(2,4),2)*pow(SPB(6,2),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(1,3,4),-1))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(6,5)*SS(1,3,4))-
 (pow(SPA(2,5),2)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(1,5,6),-1))/(complex<T>(2,0)*SPA(5,6)*SPB(4,3)*
 (SPA(1,5)*SPB(5,2)+SPA(1,6)*SPB(6,2))*SS(1,5,6))+
 (complex<T>(1,0)*pow(SPA(1,5),2)*pow(SPB(3,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1))/(complex<T>(2,0)*SPA(5,6)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(4,3)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(1,4),2)*pow(SPB(6,1),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(2,5,6),-1))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(1,5)*SPB(5,2)+SPA(1,6)*SPB(6,2))*SPB(6,5)*SS(2,5,6)))
        );

}




template <class T> complex<T> R2q2Q2l_qpqmQmQpemep_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, Qm, Qp, em, ep}, sl}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpqmQmQpemep sl");
#endif

 return( (complex<T>(0,3)*(-((SPA(2,5)*SPB(4,1)*(SPA(1,3)*SPB(6,1)-
SPA(3,4)*SPB(6,4)))/(S(3,4)*S(5,6)*SS(1,3,4)))+
(SPA(2,3)*(SPA(2,5)*SPB(4,2)+SPA(3,5)*SPB(4,3))*SPB(6,1))/
 (S(3,4)*S(5,6)*SS(2,3,4))))/complex<T>(2,0)+
 complex<T>(0,-1)*(-((pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
 pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
  (-S(1,3)-S(1,4)+S(2,3)+S(2,4))*SPA(3,5)*
  ((S(1,2)-S(3,4)-S(5,6))*SPA(3,5)+complex<T>(-2,0)*SPA(3,4)*
SPA(5,6)*SPB(6,4)))/(complex<T>(4,0)*SPA(3,4)*SPA(5,6)*
  (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,5)-S(1,6)+S(2,5)+
  S(2,6))*SPA(3,5)*(-((S(1,2)-S(3,4)-S(5,6))*SPA(3,5))+
  complex<T>(2,0)*SPA(3,4)*SPA(5,6)*SPB(6,4)))/(complex<T>(4,0)*SPA(3,4)*
 SPA(5,6)*(SPA(1,5)*SPB(5,2)+SPA(1,6)*SPB(6,2)))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,5)+S(1,6)-S(2,5)-
  S(2,6))*SPB(6,4)*((S(1,2)-S(3,4)-S(5,6))*SPB(6,4)+
  complex<T>(-2,0)*SPA(3,5)*SPB(4,3)*SPB(6,5)))/(complex<T>(4,0)*SPB(4,3)*
 (SPA(1,5)*SPB(5,2)+SPA(1,6)*SPB(6,2))*SPB(6,5))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,3)+S(1,4)-S(2,3)-
  S(2,4))*SPB(6,4)*(-((S(1,2)-S(3,4)-S(5,6))*SPB(6,4))+
  complex<T>(2,0)*SPA(3,5)*SPB(4,3)*SPB(6,5)))/
(complex<T>(4,0)*(SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(4,3)*
 SPB(6,5))+(complex<T>(1,0)*pow(SPA(2,3),2)*pow(SPB(6,2),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(1,3,4),-1))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(6,5)*SS(1,3,4))+
 (complex<T>(1,0)*pow(SPA(2,5),2)*pow(SPB(4,2),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(1,5,6),-1))/(complex<T>(2,0)*SPA(5,6)*SPB(4,3)*
 (SPA(1,5)*SPB(5,2)+SPA(1,6)*SPB(6,2))*SS(1,5,6))-
 (pow(SPA(1,5),2)*pow(SPB(4,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1))/(complex<T>(2,0)*SPA(5,6)*
 (SPA(1,3)*SPB(3,2)+SPA(1,4)*SPB(4,2))*SPB(4,3)*SS(2,3,4))-
 (pow(SPA(1,3),2)*pow(SPB(6,1),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(2,5,6),-1))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(1,5)*SPB(5,2)+SPA(1,6)*SPB(6,2))*SPB(6,5)*SS(2,5,6)))
        );

}




template <class T> complex<T> R2q2Q2l_qmqpQpQmemep_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, em, ep}, sl}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmqpQpQmemep sl");
#endif

 return( (complex<T>(0,3)*(-((SPA(1,4)*(SPA(1,5)*SPB(3,1)-SPA(4,5)*SPB(4,3))*
   SPB(6,2))/(S(3,4)*S(5,6)*SS(1,3,4)))+
(SPA(1,5)*SPB(3,2)*(SPA(2,4)*SPB(6,2)+SPA(3,4)*SPB(6,3)))/
 (S(3,4)*S(5,6)*SS(2,3,4))))/complex<T>(2,0)+
 complex<T>(0,-1)*((complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
   pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
   complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
 (S(1,3)+S(1,4)-S(2,3)-S(2,4))*SPA(4,5)*
 (-((S(1,2)-S(3,4)-S(5,6))*SPA(4,5))+complex<T>(-2,0)*SPA(3,4)*
   SPA(5,6)*SPB(6,3)))/(complex<T>(4,0)*SPA(3,4)*SPA(5,6)*
 (SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1)))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,5)+S(1,6)-S(2,5)-
  S(2,6))*SPA(4,5)*((S(1,2)-S(3,4)-S(5,6))*SPA(4,5)+
  complex<T>(2,0)*SPA(3,4)*SPA(5,6)*SPB(6,3)))/(complex<T>(4,0)*SPA(3,4)*
 SPA(5,6)*(SPA(2,5)*SPB(5,1)+SPA(2,6)*SPB(6,1)))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,5)-S(1,6)+S(2,5)+
  S(2,6))*SPB(6,3)*(-((S(1,2)-S(3,4)-S(5,6))*SPB(6,3))+
  complex<T>(-2,0)*SPA(4,5)*SPB(4,3)*SPB(6,5)))/(complex<T>(4,0)*SPB(4,3)*
 (SPA(2,5)*SPB(5,1)+SPA(2,6)*SPB(6,1))*SPB(6,5))+
 (complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,3)-S(1,4)+S(2,3)+
  S(2,4))*SPB(6,3)*((S(1,2)-S(3,4)-S(5,6))*SPB(6,3)+
  complex<T>(2,0)*SPA(4,5)*SPB(4,3)*SPB(6,5)))/
(complex<T>(4,0)*(SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(4,3)*
 SPB(6,5))-(pow(SPA(2,5),2)*pow(SPB(3,2),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(1,3,4),-1))/(complex<T>(2,0)*SPA(5,6)*
 (SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(4,3)*SS(1,3,4))-
 (pow(SPA(2,4),2)*pow(SPB(6,2),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(1,5,6),-1))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(2,5)*SPB(5,1)+SPA(2,6)*SPB(6,1))*SPB(6,5)*SS(1,5,6))+
 (complex<T>(1,0)*pow(SPA(1,4),2)*pow(SPB(6,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(6,5)*SS(2,3,4))+
 (complex<T>(1,0)*pow(SPA(1,5),2)*pow(SPB(3,1),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(2,5,6),-1))/(complex<T>(2,0)*SPA(5,6)*SPB(4,3)*
 (SPA(2,5)*SPB(5,1)+SPA(2,6)*SPB(6,1))*SS(2,5,6)))
        );

}




template <class T> complex<T> R2q2Q2l_qmqpQmQpemep_sl
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, em, ep}, sl}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmqpQmQpemep sl");
#endif

 return( (complex<T>(0,-3)*(-((SPA(1,3)*(SPA(1,5)*SPB(4,1)+SPA(3,5)*SPB(4,3))*
   SPB(6,2))/(S(3,4)*S(5,6)*SS(1,3,4)))+
(SPA(1,5)*SPB(4,2)*(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4)))/
 (S(3,4)*S(5,6)*SS(2,3,4))))/complex<T>(2,0)+
 complex<T>(0,1)*((complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
   pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
   complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
 (S(1,5)+S(1,6)-S(2,5)-S(2,6))*SPA(3,5)*
 ((S(1,2)-S(3,4)-S(5,6))*SPA(3,5)+complex<T>(-2,0)*SPA(3,4)*
   SPA(5,6)*SPB(6,4)))/(complex<T>(4,0)*SPA(3,4)*SPA(5,6)*
 (SPA(2,5)*SPB(5,1)+SPA(2,6)*SPB(6,1)))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,3)+S(1,4)-S(2,3)-
  S(2,4))*SPA(3,5)*(-((S(1,2)-S(3,4)-S(5,6))*SPA(3,5))+
  complex<T>(2,0)*SPA(3,4)*SPA(5,6)*SPB(6,4)))/(complex<T>(4,0)*SPA(3,4)*
 SPA(5,6)*(SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1)))-
 (pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
   complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
   complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,3)-S(1,4)+S(2,3)+
  S(2,4))*SPB(6,4)*((S(1,2)-S(3,4)-S(5,6))*SPB(6,4)+
  complex<T>(-2,0)*SPA(3,5)*SPB(4,3)*SPB(6,5)))/
(complex<T>(4,0)*(SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(4,3)*
 SPB(6,5))+(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+
   pow(SPA(3,4),2)*pow(SPB(4,3),2)+pow(SPA(5,6),2)*
pow(SPB(6,5),2)+complex<T>(-2,0)*S(1,2)*S(3,4)+
   complex<T>(-2,0)*S(1,2)*S(5,6)+complex<T>(-2,0)*S(3,4)*S(5,6),-1)*
 (-S(1,5)-S(1,6)+S(2,5)+S(2,6))*SPB(6,4)*
 (-((S(1,2)-S(3,4)-S(5,6))*SPB(6,4))+complex<T>(2,0)*SPA(3,5)*
   SPB(4,3)*SPB(6,5)))/(complex<T>(4,0)*SPB(4,3)*(SPA(2,5)*SPB(5,1)+
  SPA(2,6)*SPB(6,1))*SPB(6,5))+(complex<T>(1,0)*pow(SPA(2,5),2)*
 pow(SPB(4,2),2)*pow(complex<T>(1,0)-S(5,6)/SS(1,3,4),-1))/
(complex<T>(2,0)*SPA(5,6)*(SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*
 SPB(4,3)*SS(1,3,4))+(complex<T>(1,0)*pow(SPA(2,3),2)*
 pow(SPB(6,2),2)*pow(complex<T>(1,0)-S(3,4)/SS(1,5,6),-1))/
(complex<T>(2,0)*SPA(3,4)*(SPA(2,5)*SPB(5,1)+SPA(2,6)*SPB(6,1))*
 SPB(6,5)*SS(1,5,6))-(pow(SPA(1,3),2)*pow(SPB(6,1),2)*
 pow(complex<T>(1,0)-S(5,6)/SS(2,3,4),-1))/(complex<T>(2,0)*SPA(3,4)*
 (SPA(2,3)*SPB(3,1)+SPA(2,4)*SPB(4,1))*SPB(6,5)*SS(2,3,4))-
 (pow(SPA(1,5),2)*pow(SPB(4,1),2)*
 pow(complex<T>(1,0)-S(3,4)/SS(2,5,6),-1))/(complex<T>(2,0)*SPA(5,6)*SPB(4,3)*
 (SPA(2,5)*SPB(5,1)+SPA(2,6)*SPB(6,1))*SS(2,5,6)))
        );

}




template <class T> complex<T> R2q2Q2l_qpqmQpQmemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, Qp, Qm, em, ep}, AX}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpqmQpQmemep AX");
#endif

 complex<T> mt = complex<T>(BH::constants::s_GeV,0)*complex<T>(BH::constants::Mtop,0);
 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( (complex<T>(0,-2)*(complex<T>(1,0)/(complex<T>(24,0)*pow(mt,2))+
(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)+S(3,4)-S(5,6)))/
 complex<T>(2,0)+(complex<T>(1,0)*(S(1,2)+complex<T>(2,0)*S(3,4)+S(5,6)))/
 (complex<T>(360,0)*pow(mt,4)))*(-((SPA(4,5)*SPB(3,1)*SPB(6,1))/
  SPB(2,1))-(SPA(2,4)*SPA(2,5)*SPB(6,3))/SPA(1,2)))/S(5,6)+
 (complex<T>(0,-2)*(complex<T>(1,0)/(complex<T>(24,0)*pow(mt,2))+
(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,2)-S(3,4)-S(5,6)))/
 complex<T>(2,0)+(complex<T>(1,0)*(complex<T>(2,0)*S(1,2)+S(3,4)+S(5,6)))/
 (complex<T>(360,0)*pow(mt,4)))*((SPA(2,4)*SPA(4,5)*SPB(6,1))/SPA(3,4)+
(SPA(2,5)*SPB(3,1)*SPB(6,3))/SPB(4,3)))/S(5,6)
        );

}




template <class T> complex<T> R2q2Q2l_qmqpQpQmemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qp, Qm, em, ep}, AX}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmqpQpQmemep AX");
#endif

 complex<T> mt = complex<T>(BH::constants::s_GeV,0)*complex<T>(BH::constants::Mtop,0);
 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( (complex<T>(0,2)*(complex<T>(1,0)/(complex<T>(24,0)*pow(mt,2))+
(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)+S(3,4)-S(5,6)))/
 complex<T>(2,0)+(complex<T>(1,0)*(S(1,2)+complex<T>(2,0)*S(3,4)+S(5,6)))/
 (complex<T>(360,0)*pow(mt,4)))*((SPA(4,5)*SPB(3,2)*SPB(6,2))/SPB(2,1)+
(SPA(1,4)*SPA(1,5)*SPB(6,3))/SPA(1,2)))/S(5,6)+
 (complex<T>(0,2)*(complex<T>(1,0)/(complex<T>(24,0)*pow(mt,2))+
(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,2)-S(3,4)-S(5,6)))/
 complex<T>(2,0)+(complex<T>(1,0)*(complex<T>(2,0)*S(1,2)+S(3,4)+S(5,6)))/
 (complex<T>(360,0)*pow(mt,4)))*((SPA(1,4)*SPA(4,5)*SPB(6,2))/SPA(3,4)+
(SPA(1,5)*SPB(3,2)*SPB(6,3))/SPB(4,3)))/S(5,6)
        );

}




template <class T> complex<T> R2q2Q2l_qpqmQmQpemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qp, qm, Qm, Qp, em, ep}, AX}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qpqmQmQpemep AX");
#endif

 complex<T> mt = complex<T>(BH::constants::s_GeV,0)*complex<T>(BH::constants::Mtop,0);
 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( (complex<T>(0,2)*(complex<T>(1,0)/(complex<T>(24,0)*pow(mt,2))+
(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)+S(3,4)-S(5,6)))/
 complex<T>(2,0)+(complex<T>(1,0)*(S(1,2)+complex<T>(2,0)*S(3,4)+S(5,6)))/
 (complex<T>(360,0)*pow(mt,4)))*(-((SPA(3,5)*SPB(4,1)*SPB(6,1))/
  SPB(2,1))-(SPA(2,3)*SPA(2,5)*SPB(6,4))/SPA(1,2)))/S(5,6)+
 (complex<T>(0,2)*(complex<T>(1,0)/(complex<T>(24,0)*pow(mt,2))+
(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,2)-S(3,4)-S(5,6)))/
 complex<T>(2,0)+(complex<T>(1,0)*(complex<T>(2,0)*S(1,2)+S(3,4)+S(5,6)))/
 (complex<T>(360,0)*pow(mt,4)))*(-((SPA(2,3)*SPA(3,5)*SPB(6,1))/
  SPA(3,4))-(SPA(2,5)*SPB(4,1)*SPB(6,4))/SPB(4,3)))/S(5,6)
        );

}




template <class T> complex<T> R2q2Q2l_qmqpQmQpemep_AX
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
      // {{qm, qp, Qm, Qp, em, ep}, AX}

#if _VERBOSE
  _MESSAGE("R2q2Q2l :  qmqpQmQpemep AX");
#endif

 complex<T> mt = complex<T>(BH::constants::s_GeV,0)*complex<T>(BH::constants::Mtop,0);
 complex<T> mtsq = pow(mt, 2);     // top mass temporary setting
 return( (complex<T>(0,-2)*(complex<T>(1,0)/(complex<T>(24,0)*pow(mt,2))+
(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(-S(1,2)+S(3,4)-S(5,6)))/
 complex<T>(2,0)+(complex<T>(1,0)*(S(1,2)+complex<T>(2,0)*S(3,4)+S(5,6)))/
 (complex<T>(360,0)*pow(mt,4)))*((SPA(3,5)*SPB(4,2)*SPB(6,2))/SPB(2,1)+
(SPA(1,3)*SPA(1,5)*SPB(6,4))/SPA(1,2)))/S(5,6)+
 (complex<T>(0,-2)*(complex<T>(1,0)/(complex<T>(24,0)*pow(mt,2))+
(complex<T>(1,0)*pow(pow(SPA(1,2),2)*pow(SPB(2,1),2)+pow(SPA(3,4),2)*
 pow(SPB(4,3),2)+pow(SPA(5,6),2)*pow(SPB(6,5),2)+
complex<T>(-2,0)*S(1,2)*S(3,4)+complex<T>(-2,0)*S(1,2)*S(5,6)+
complex<T>(-2,0)*S(3,4)*S(5,6),-1)*(S(1,2)-S(3,4)-S(5,6)))/
 complex<T>(2,0)+(complex<T>(1,0)*(complex<T>(2,0)*S(1,2)+S(3,4)+S(5,6)))/
 (complex<T>(360,0)*pow(mt,4)))*(-((SPA(1,3)*SPA(3,5)*SPB(6,2))/
  SPA(3,4))-(SPA(1,5)*SPB(4,2)*SPB(6,4))/SPB(4,3)))/S(5,6)
        );

}



 // *************** table of switch values *************

#define _R_qpQmQpqmemep_nf R2q2Q2l_254161_nf
#define _R_qpQpQmqmemep_nf R2q2Q2l_254105_nf
#define _R_qmQpQmqpemep_nf R2q2Q2l_254616_nf
#define _R_qmQmQpqpemep_nf R2q2Q2l_254672_nf

#define _R_qpQmQpqmemep_nf_top R2q2Q2l_254161_nf_top
#define _R_qpQpQmqmemep_nf_top R2q2Q2l_254105_nf_top
#define _R_qmQpQmqpemep_nf_top R2q2Q2l_254616_nf_top
#define _R_qmQmQpqpemep_nf_top R2q2Q2l_254672_nf_top

#define _R_qpQmQpqmemep_L R2q2Q2l_254161_L
#define _R_qpQpQmqmemep_L R2q2Q2l_254105_L
#define _R_qmQpQmqpemep_L R2q2Q2l_254616_L
#define _R_qmQmQpqpemep_L R2q2Q2l_254672_L
#define _R_qpqmQpQmemep_sl R2q2Q2l_255169_sl
#define _R_qpqmQmQpemep_sl R2q2Q2l_255617_sl
#define _R_qmqpQpQmemep_sl R2q2Q2l_255176_sl
#define _R_qmqpQmQpemep_sl R2q2Q2l_255624_sl
#define _R_qpqmQpQmemep_AX R2q2Q2l_255169_AX
#define _R_qmqpQpQmemep_AX R2q2Q2l_255176_AX
#define _R_qpqmQmQpemep_AX R2q2Q2l_255617_AX
#define _R_qmqpQmQpemep_AX R2q2Q2l_255624_AX


 // *************** more macro definitions *************

#define _CASE_qpQmQpqmemep_nf_top case 254161 : \
          return &R2q2Q2l_254161_nf_top
#define _CASE_qpQpQmqmemep_nf_top case 254105 : \
          return &R2q2Q2l_254105_nf_top
#define _CASE_qmQpQmqpemep_nf_top case 254616 : \
          return &R2q2Q2l_254616_nf_top
#define _CASE_qmQmQpqpemep_nf_top case 254672 : \
          return &R2q2Q2l_254672_nf_top


#define _CASE_qpQmQpqmemep_nf case 254161 : \
          return &R2q2Q2l_254161_nf
#define _CASE_qpQpQmqmemep_nf case 254105 : \
          return &R2q2Q2l_254105_nf
#define _CASE_qmQpQmqpemep_nf case 254616 : \
          return &R2q2Q2l_254616_nf
#define _CASE_qmQmQpqpemep_nf case 254672 : \
          return &R2q2Q2l_254672_nf




#define _CASE_qpQmQpqmemep_L case 254161 : \
          return &R2q2Q2l_254161_L
#define _CASE_qpQpQmqmemep_L case 254105 : \
          return &R2q2Q2l_254105_L
#define _CASE_qmQpQmqpemep_L case 254616 : \
          return &R2q2Q2l_254616_L
#define _CASE_qmQmQpqpemep_L case 254672 : \
          return &R2q2Q2l_254672_L
#define _CASE_qpqmQpQmemep_sl case 255169 : \
          return &R2q2Q2l_255169_sl
#define _CASE_qpqmQmQpemep_sl case 255617 : \
          return &R2q2Q2l_255617_sl
#define _CASE_qmqpQpQmemep_sl case 255176 : \
          return &R2q2Q2l_255176_sl
#define _CASE_qmqpQmQpemep_sl case 255624 : \
          return &R2q2Q2l_255624_sl
#define _CASE_qpqmQpQmemep_AX case 255169 : \
          return &R2q2Q2l_255169_AX
#define _CASE_qmqpQpQmemep_AX case 255176 : \
          return &R2q2Q2l_255176_AX
#define _CASE_qpqmQmQpemep_AX case 255617 : \
          return &R2q2Q2l_255617_AX
#define _CASE_qmqpQmQpemep_AX case 255624 : \
          return &R2q2Q2l_255624_AX


 // *************** function definitions using macros *************


template <class T> complex<T> _R_qpQmQpqmemep_nf_top(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpQmQpqmemep_nf_top(ep,mpc);}

template <class T> complex<T> _R_qpQpQmqmemep_nf_top(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpQpQmqmemep_nf_top(ep,mpc);}

template <class T> complex<T> _R_qmQpQmqpemep_nf_top(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmQpQmqpemep_nf_top(ep,mpc);}

template <class T> complex<T> _R_qmQmQpqpemep_nf_top(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmQmQpqpemep_nf_top(ep,mpc);}



template <class T> complex<T> _R_qpQmQpqmemep_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpQmQpqmemep_nf(ep,mpc);}

template <class T> complex<T> _R_qpQpQmqmemep_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpQpQmqmemep_nf(ep,mpc);}

template <class T> complex<T> _R_qmQpQmqpemep_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmQpQmqpemep_nf(ep,mpc);}

template <class T> complex<T> _R_qmQmQpqpemep_nf(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmQmQpqpemep_nf(ep,mpc);}

template <class T> complex<T> _R_qpQmQpqmemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpQmQpqmemep_L(ep,mpc);}

template <class T> complex<T> _R_qpQpQmqmemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpQpQmqmemep_L(ep,mpc);}

template <class T> complex<T> _R_qmQpQmqpemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmQpQmqpemep_L(ep,mpc);}

template <class T> complex<T> _R_qmQmQpqpemep_L(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmQmQpqpemep_L(ep,mpc);}

template <class T> complex<T> _R_qpqmQpQmemep_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpqmQpQmemep_sl(ep,mpc);}

template <class T> complex<T> _R_qpqmQmQpemep_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpqmQmQpemep_sl(ep,mpc);}

template <class T> complex<T> _R_qmqpQpQmemep_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmqpQpQmemep_sl(ep,mpc);}

template <class T> complex<T> _R_qmqpQmQpemep_sl(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmqpQmQpemep_sl(ep,mpc);}

template <class T> complex<T> _R_qpqmQpQmemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpqmQpQmemep_AX(ep,mpc);}

template <class T> complex<T> _R_qmqpQpQmemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmqpQpQmemep_AX(ep,mpc);}

template <class T> complex<T> _R_qpqmQmQpemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qpqmQmQpemep_AX(ep,mpc);}

template <class T> complex<T> _R_qmqpQmQpemep_AX(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q2l_qmqpQmQpemep_AX(ep,mpc);}




 // *************** define pointers *************

template <class T> complex<T> ( *R2q2Q2l_AX_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpqmQpQmemep_AX;
       _CASE_qmqpQpQmemep_AX;
       _CASE_qpqmQmQpemep_AX;
       _CASE_qmqpQmQpemep_AX;

       default: return 0;
        }
 }

template <class T> complex<T> ( *R2q2Q2l_L_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpQmQpqmemep_L;
       _CASE_qpQpQmqmemep_L;
       _CASE_qmQpQmqpemep_L;
       _CASE_qmQmQpqpemep_L;

       default: return 0;
        }
 }

template <class T> complex<T> ( *R2q2Q2l_nf_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpQmQpqmemep_nf;
       _CASE_qpQpQmqmemep_nf;
       _CASE_qmQpQmqpemep_nf;
       _CASE_qmQmQpqpemep_nf;

       default: return 0;
        }
 }

template <class T> complex<T> ( *R2q2Q2l_nf_top_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpQmQpqmemep_nf_top;
       _CASE_qpQpQmqmemep_nf_top;
       _CASE_qmQpQmqpemep_nf_top;
       _CASE_qmQmQpqpemep_nf_top;

       default: return 0;
        }
 }



template <class T> complex<T> ( *R2q2Q2l_sl_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpqmQpQmemep_sl;
       _CASE_qpqmQmQpemep_sl;
       _CASE_qmqpQpQmemep_sl;
       _CASE_qmqpQmQpemep_sl;

       default: return 0;
        }
 }


 // *************** definitions for template *************

template complex<R> ( *R2q2Q2l_AX_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q2l_AX_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q2l_AX_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q2l_AX_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2Q2l_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q2l_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q2l_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q2l_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2Q2l_nf_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q2l_nf_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q2l_nf_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q2l_nf_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif



template complex<R> ( *R2q2Q2l_nf_top_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q2l_nf_top_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q2l_nf_top_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q2l_nf_top_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif


template complex<R> ( *R2q2Q2l_sl_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q2l_sl_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q2l_sl_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q2l_sl_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif




}

