/*
*C_2q2l.cpp
*
* Created on 9/11, 2009
*      Author: Zvi's script
*/
 
#include "C_2q2l_eval.h"
#include "integrals_ep.h"
 
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


 
 
template <class T> SeriesC<T> C2q2l_qpqmemep_L
      (const eval_param<T>& ep,
                 const T& mu){
//{{qp, qm, em, ep}, L}
 
#if _VERBOSE
  _MESSAGE("C2q2l :  qpqmemep L");
#endif
 
	//define corner vectors
	 vector<int> c1;  c1.push_back(1-1);
	 vector<int> c2;  c2.push_back(2-1);
	 vector<int> c3;  c3.push_back(3-1);
	 vector<int> c4;  c4.push_back(4-1);

	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
 	 vector<int> c13;  c13.push_back(1-1); c13.push_back(3-1);
         vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
         vector<int> c24;  c24.push_back(2-1); c24.push_back(4-1);
	 vector<int> c41;  c41.push_back(4-1); c41.push_back(1-1);
         vector<int> c14 = c41;
 // #define TimeStamp "Fri 11 Sep 2009 15:46:54 on login2"
 // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
#define Complex complex<T>
complex<T> spa12 = SPA(1,2);
complex<T> spa23 = SPA(2,3);
complex<T> spa34 = SPA(3,4);
complex<T> spb34 = SPB(3,4);
complex<T> t1 = square(spa23); 
complex<T> d1 = spa12*spa34*T(2); d1 = T(1)/d1;
complex<T> d2 = spa12; d2 = T(1)/d2;
complex<T> co1 = d1*t1*T(3); 
complex<T> co2 = d2*spb34*t1; 
complex<T> co3 = Complex(0,1); 
SeriesC<T> result = co3*(co1*Int(ep,mu,c34,c12) + co2*Int(ep,mu,c3,c4,c12));  
 return(result);
} 
  
  
 
 
  
 
 
 // *************** table of switch values ************* 
 
#define _C_qpqmemep_L C2q2l_1232_L
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qpqmemep_L case 1232 : \
          return &C2q2l_1232_L
 
 
 // *************** function definitions using macros ************* 
 
template <class T> SeriesC<T> _C_qpqmemep_L(
        const eval_param<T>& ep,  const T& mu){
          return C2q2l_qpqmemep_L(ep,mu);}
 
 
 
 
 
 // *************** define pointers ************* 
 
template <class T> SeriesC<T> ( *C2q2l_L_Ptr_eval( int hc))
     (const eval_param<T>& ,const T& mu) {
       switch (hc) {
       _CASE_qpqmemep_L;
 
       default: return 0;
        }
 }
 

 // *************** definitions for template ************* 

template SeriesC<R> ( *C2q2l_L_Ptr_eval(int hc))
             (const eval_param<R>& ,const R& mu);
template SeriesC<RHP> ( *C2q2l_L_Ptr_eval(int hc))
             (const eval_param<RHP>& ,const RHP& mu);
template SeriesC<RVHP> ( *C2q2l_L_Ptr_eval(int hc))
             (const eval_param<RVHP>& ,const RVHP& mu);


}
 
