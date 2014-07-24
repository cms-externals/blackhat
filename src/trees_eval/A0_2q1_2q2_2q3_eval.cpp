/* The 6 quark amplitudes with 3 distinguished flavors
  See bottom of file for table of switch numbers and values
*/

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;
namespace BH {


template <class T> complex<T> A2q1_2q2_2q3zero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q1_2q2_2q31865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
 2)*(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0)))/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
(ep.spa(1,3)*ep.spb(4,1)+
 ep.spa(2,3)*ep.spb(4,2)))/
 (ep.s(1,2,3)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q31870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2)*
 (-(ep.spa(3,5)*ep.spb(3,0))-
  ep.spa(4,5)*ep.spb(4,0)))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))-(pow(ep.spa(0,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*ep.spb(3,2)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))+
pow(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4),3)/
 (ep.s(0,4,5)*ep.spa(2,3)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q31895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q31905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(1,0),2))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q31930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q31935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
 2))/(ep.s(3,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0)))-(pow(ep.spa(1,3),2)*
pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
ep.spa(1,2)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0),3)/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2)*(ep.spa(0,1)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(5,2)))/(ep.s(0,1,5)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4)))/
 (ep.s(0,4,5)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q32050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spb(3,2)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
pow(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1),3)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(0,2),2)*
pow(ep.spb(4,3),2)*
(-(ep.spa(0,1)*ep.spb(5,1))-
 ep.spa(0,2)*ep.spb(5,2)))/
 (ep.s(0,1,2)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
 ep.spa(4,5)*ep.spb(4,1),2)/
 (ep.s(3,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0)))-
pow(ep.spa(0,2)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2),2)/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))+
pow(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(1,0),2))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
 2))/(ep.s(3,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0)))-(pow(ep.spa(1,3),2)*
pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
ep.spa(1,2)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
 2)*(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0)))/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
(ep.spa(1,3)*ep.spb(4,1)+
 ep.spa(2,3)*ep.spb(4,2)))/
 (ep.s(1,2,3)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2)*
 (-(ep.spa(3,5)*ep.spb(3,0))-
  ep.spa(4,5)*ep.spb(4,0)))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))-(pow(ep.spa(0,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*ep.spb(3,2)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))+
pow(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4),3)/
 (ep.s(0,4,5)*ep.spa(2,3)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q32470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
 ep.spa(4,5)*ep.spb(4,1),2)/
 (ep.s(3,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0)))-
pow(ep.spa(0,2)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2),2)/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))+
pow(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0),3)/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2)*(ep.spa(0,1)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(5,2)))/(ep.s(0,1,5)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4)))/
 (ep.s(0,4,5)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q32540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spb(3,2)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
pow(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1),3)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(0,2),2)*
pow(ep.spb(4,3),2)*
(-(ep.spa(0,1)*ep.spb(5,1))-
 ep.spa(0,2)*ep.spb(5,2)))/
 (ep.s(0,1,2)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q32945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q32985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
 ep.spa(3,4)*ep.spb(3,0),2)/
 (ep.s(2,3,4)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,3)+
 ep.spa(1,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q33310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1),2)/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))+
pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q33395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(4,5),2)*pow(ep.spb(3,1),
 2)*(ep.spa(0,2)*ep.spb(2,1)+
 ep.spa(0,3)*ep.spb(3,1)))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spb(2,1)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))+
pow(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3),3)/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(2,4),2)*
pow(ep.spb(1,0),2)*(ep.spa(2,4)*
  ep.spb(5,2)+ep.spa(3,4)*
  ep.spb(5,3)))/(ep.s(2,3,4)*
ep.spa(3,4)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*ep.spb(5,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q33430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spb(2,1)*
 (ep.spa(0,5)*ep.spb(2,0)+
  ep.spa(1,5)*ep.spb(2,1))*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(1,5),2)*pow(ep.spb(3,2),2)*
(ep.spa(0,5)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(4,1)))/
 (ep.s(0,1,5)*ep.spa(0,5)*
(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1))*
ep.spb(4,3)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,0),3)/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spb(5,0)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q33515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1),2)/
 (ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q33730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q33835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
  ep.spa(1,5)*ep.spb(3,1),2)/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4))))+
pow(ep.spa(2,4)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q34205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,2),2))/(ep.s(2,3,4)*
ep.spa(0,1)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q34390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))+
(pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q34475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
  ep.spa(3,5)*ep.spb(4,3),3)/
(ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)*(-(ep.spa(1,2)*
ep.spb(4,2))-ep.spa(1,3)*
   ep.spb(4,3))))+(pow(ep.spa(2,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,4)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,0)))/(ep.s(0,4,5)*
ep.spa(1,2)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(3,5),2)*pow(ep.spb(1,0),2)*
(ep.spa(3,4)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(5,2)))/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(2,1)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))))
); }

template <class T> complex<T> A2q1_2q2_2q34510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*ep.spb(2,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2)*
(ep.spa(2,5)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(4,3)))/
 (ep.s(2,3,4)*ep.spa(0,5)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))+(pow(ep.spa(1,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,2)*
  ep.spb(2,0)+ep.spa(1,3)*
  ep.spb(3,0)))/(ep.s(1,2,3)*
ep.spa(1,2)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q34805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q34945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))+
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q35020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q35495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1),2)/
 (ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q35805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q35935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
  ep.spa(1,5)*ep.spb(3,1),2)/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4))))+
pow(ep.spa(2,4)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q36020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
 ep.spa(3,4)*ep.spb(3,0),2)/
 (ep.s(2,3,4)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,3)+
 ep.spa(1,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q36320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1),2)/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))+
pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q36345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(4,5),2)*pow(ep.spb(3,1),
 2)*(ep.spa(0,2)*ep.spb(2,1)+
 ep.spa(0,3)*ep.spb(3,1)))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spb(2,1)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))+
pow(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3),3)/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(2,4),2)*
pow(ep.spb(1,0),2)*(ep.spa(2,4)*
  ep.spb(5,2)+ep.spa(3,4)*
  ep.spb(5,3)))/(ep.s(2,3,4)*
ep.spa(3,4)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*ep.spb(5,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q36380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spb(2,1)*
 (ep.spa(0,5)*ep.spb(2,0)+
  ep.spa(1,5)*ep.spb(2,1))*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(1,5),2)*pow(ep.spb(3,2),2)*
(ep.spa(0,5)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(4,1)))/
 (ep.s(0,1,5)*ep.spa(0,5)*
(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1))*
ep.spb(4,3)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,0),3)/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spb(5,0)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q36790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q36975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))+
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q37065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q37280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,2),2))/(ep.s(2,3,4)*
ep.spa(0,1)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q37400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))+
(pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q37425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
  ep.spa(3,5)*ep.spb(4,3),3)/
(ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)*(-(ep.spa(1,2)*
ep.spb(4,2))-ep.spa(1,3)*
   ep.spb(4,3))))+(pow(ep.spa(2,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,4)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,0)))/(ep.s(0,4,5)*
ep.spa(1,2)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(3,5),2)*pow(ep.spb(1,0),2)*
(ep.spa(3,4)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(5,2)))/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(2,1)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))))
); }

template <class T> complex<T> A2q1_2q2_2q37460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q37465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*ep.spb(2,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2)*
(ep.spa(2,5)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(4,3)))/
 (ep.s(2,3,4)*ep.spa(0,5)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))+(pow(ep.spa(1,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,2)*
  ep.spb(2,0)+ep.spa(1,3)*
  ep.spb(3,0)))/(ep.s(1,2,3)*
ep.spa(1,2)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q38345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2)*
 (-(ep.spa(1,5)*ep.spb(1,0))-
  ep.spa(2,5)*ep.spb(2,0)))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(ep.spa(1,3)*ep.spb(3,2)+
 ep.spa(1,4)*ep.spb(4,2),3)/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spb(3,2)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))-(pow(ep.spa(1,3),2)*
pow(ep.spb(5,0),2)*(ep.spa(1,3)*
  ep.spb(4,1)+ep.spa(2,3)*
  ep.spb(4,2)))/(ep.s(1,2,3)*
ep.spa(2,3)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(2,1),2)*(ep.spa(0,4)*
   ep.spb(3,0)+ep.spa(4,5)*
   ep.spb(5,3)))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+
pow(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4),3)/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q38375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q38385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1),2)/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))))+
pow(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q38415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
 ep.spa(4,5)*ep.spb(5,2),2)/
 (ep.s(0,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))-
pow(-(ep.spa(1,3)*ep.spb(5,3))-
 ep.spa(1,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(1,4)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(3,2),3)/
 (ep.s(1,2,3)*ep.spa(4,5)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(3,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,0),2)*(ep.spa(2,3)*
  ep.spb(3,1)+ep.spa(2,4)*
  ep.spb(4,1)))/(ep.s(2,3,4)*
ep.spa(2,3)*ep.spb(1,0)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q38530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(2,3)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(4,1)))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spb(1,0)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(0,2),2)*
pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q38585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q38600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),
 2))/(ep.s(0,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))-(pow(ep.spa(0,2),2)*
pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
ep.spa(1,2)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q38630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,2),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))+
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1),2)/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))))+
pow(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q38765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
 ep.spa(4,5)*ep.spb(5,2),2)/
 (ep.s(0,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))-
pow(-(ep.spa(1,3)*ep.spb(5,3))-
 ep.spa(1,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q38835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2)*
 (-(ep.spa(1,5)*ep.spb(1,0))-
  ep.spa(2,5)*ep.spb(2,0)))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(ep.spa(1,3)*ep.spb(3,2)+
 ep.spa(1,4)*ep.spb(4,2),3)/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spb(3,2)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))-(pow(ep.spa(1,3),2)*
pow(ep.spb(5,0),2)*(ep.spa(1,3)*
  ep.spb(4,1)+ep.spa(2,3)*
  ep.spb(4,2)))/(ep.s(1,2,3)*
ep.spa(2,3)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(2,1),2)*(ep.spa(0,4)*
   ep.spb(3,0)+ep.spa(4,5)*
   ep.spb(5,3)))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+
pow(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4),3)/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q38950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),
 2))/(ep.s(0,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))-(pow(ep.spa(0,2),2)*
pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
ep.spa(1,2)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q38980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,2),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))+
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q38990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q39015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(1,4)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(3,2),3)/
 (ep.s(1,2,3)*ep.spa(4,5)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(3,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,0),2)*(ep.spa(2,3)*
  ep.spb(3,1)+ep.spa(2,4)*
  ep.spb(4,1)))/(ep.s(2,3,4)*
ep.spa(2,3)*ep.spb(1,0)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q39020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(2,3)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(4,1)))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spb(1,0)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(0,2),2)*
pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q310505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q310510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q310535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q310545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q310570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q310575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(2,4),2)*
 pow(ep.spb(5,0),2))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))+
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q311050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q311165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2)*(ep.spa(0,2)*
   ep.spb(3,0)+ep.spa(1,2)*
   ep.spb(3,1)))/(ep.s(0,1,2)*
 ep.spa(1,2)*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))*ep.spb(4,3)*
 (ep.spa(0,2)*ep.spb(5,0)+
  ep.spa(1,2)*ep.spb(5,1))))+
(pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
(ep.spa(0,2)*ep.spb(2,1)+
 ep.spa(0,3)*ep.spb(3,1)))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spb(2,1)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2)*(ep.spa(0,4)*
  ep.spb(5,0)+ep.spa(1,4)*
  ep.spb(5,1)))/(ep.s(0,1,5)*
ep.spa(3,4)*ep.spb(5,0)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q311200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1),3)/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
pow(ep.spb(3,2),2)*(ep.spa(0,4)*
  ep.spb(4,1)+ep.spa(0,5)*
  ep.spb(5,1)))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3)))/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))))
); }

template <class T> complex<T> A2q1_2q2_2q311255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(2,1),2))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1)))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q311470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q311765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q311870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0),2)/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q312130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
 ep.spa(3,5)*ep.spb(5,1),2)/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2)))-
pow(-(ep.spa(0,2)*ep.spb(4,2))-
 ep.spa(0,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q312245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
  ep.spa(0,3)*ep.spb(3,1),3)/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spb(2,1)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))))-(pow(ep.spa(0,2),2)*
pow(ep.spb(5,4),2)*(ep.spa(0,2)*
  ep.spb(3,0)+ep.spa(1,2)*
  ep.spb(3,1)))/(ep.s(0,1,2)*
ep.spa(1,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*ep.spb(4,3)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1)))+
(pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,4)*ep.spb(5,0)+
 ep.spa(1,4)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(3,4)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(5,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q312280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(2,1),3)/
 (ep.s(0,1,2)*ep.spa(3,4)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(2,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))+(pow(ep.spa(0,1),2)*
pow(ep.spb(4,2),2)*(ep.spa(2,5)*
  ep.spb(4,2)+ep.spa(3,5)*
  ep.spb(4,3)))/(ep.s(2,3,4)*
ep.spa(0,5)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(4,3)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))+
(pow(ep.spa(1,3),2)*pow(ep.spb(5,4),2)*
(ep.spa(1,2)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(3,0)))/
 (ep.s(1,2,3)*ep.spa(1,2)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q312545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
  ep.spa(0,5)*ep.spb(5,2),2)/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
pow(ep.spa(0,3)*ep.spb(5,0)+
 ep.spa(1,3)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q312760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q312900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-
pow(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q313055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(2,1),2))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1)))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q313575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q313790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q313800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(2,4),2)*
 pow(ep.spb(5,0),2))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))+
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q314060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q314115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2)*(ep.spa(0,2)*
   ep.spb(3,0)+ep.spa(1,2)*
   ep.spb(3,1)))/(ep.s(0,1,2)*
 ep.spa(1,2)*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))*ep.spb(4,3)*
 (ep.spa(0,2)*ep.spb(5,0)+
  ep.spa(1,2)*ep.spb(5,1))))+
(pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
(ep.spa(0,2)*ep.spb(2,1)+
 ep.spa(0,3)*ep.spb(3,1)))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spb(2,1)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2)*(ep.spa(0,4)*
  ep.spb(5,0)+ep.spa(1,4)*
  ep.spb(5,1)))/(ep.s(0,1,5)*
ep.spa(3,4)*ep.spb(5,0)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q314150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1),3)/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
pow(ep.spb(3,2),2)*(ep.spa(0,4)*
  ep.spb(4,1)+ep.spa(0,5)*
  ep.spb(5,1)))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3)))/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))))
); }

template <class T> complex<T> A2q1_2q2_2q314350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
  ep.spa(0,5)*ep.spb(5,2),2)/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
pow(ep.spa(0,3)*ep.spb(5,0)+
 ep.spa(1,3)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q314835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q314980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q315000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-
pow(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q315050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q315060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q315135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0),2)/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q315140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
 ep.spa(3,5)*ep.spb(5,1),2)/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2)))-
pow(-(ep.spa(0,2)*ep.spb(4,2))-
 ep.spa(0,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q315195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q315210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
  ep.spa(0,3)*ep.spb(3,1),3)/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spb(2,1)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))))-(pow(ep.spa(0,2),2)*
pow(ep.spb(5,4),2)*(ep.spa(0,2)*
  ep.spb(3,0)+ep.spa(1,2)*
  ep.spb(3,1)))/(ep.s(0,1,2)*
ep.spa(1,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*ep.spb(4,3)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1)))+
(pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,4)*ep.spb(5,0)+
 ep.spa(1,4)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(3,4)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(5,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q315230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q315240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(2,1),3)/
 (ep.s(0,1,2)*ep.spa(3,4)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(2,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))+(pow(ep.spa(0,1),2)*
pow(ep.spb(4,2),2)*(ep.spa(2,5)*
  ep.spb(4,2)+ep.spa(3,5)*
  ep.spb(4,3)))/(ep.s(2,3,4)*
ep.spa(0,5)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(4,3)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))+
(pow(ep.spa(1,3),2)*pow(ep.spb(5,4),2)*
(ep.spa(1,2)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(3,0)))/
 (ep.s(1,2,3)*ep.spa(1,2)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q315905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
 ep.spa(3,4)*ep.spb(3,0),2)/
 (ep.s(2,3,4)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,3)+
 ep.spa(1,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q315910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1),2)/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))+
pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q315935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q315945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(4,5),2)*
pow(ep.spb(3,1),2)*(ep.spa(0,2)*
  ep.spb(2,1)+ep.spa(0,3)*
  ep.spb(3,1)))/(ep.s(1,2,3)*
ep.spa(0,5)*ep.spb(2,1)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))+
pow(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3),3)/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(2,4),2)*
pow(ep.spb(1,0),2)*(ep.spa(2,4)*
  ep.spb(5,2)+ep.spa(3,4)*
  ep.spb(5,3)))/(ep.s(2,3,4)*
ep.spa(3,4)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*ep.spb(5,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q315970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q315975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spb(2,1)*
 (ep.spa(0,5)*ep.spb(2,0)+
  ep.spa(1,5)*ep.spb(2,1))*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(1,5),2)*pow(ep.spb(3,2),2)*
(ep.spa(0,5)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(4,1)))/
 (ep.s(0,1,5)*ep.spa(0,5)*
(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1))*
ep.spb(4,3)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,0),3)/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spb(5,0)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q316265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1),2)/
 (ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q316615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
  ep.spa(1,5)*ep.spb(3,1),2)/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4))))+
pow(ep.spa(2,4)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q316795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q316985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2)))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q316990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))+
(pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q317015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
  ep.spa(3,5)*ep.spb(4,3),3)/
(ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)*(-(ep.spa(1,2)*
ep.spb(4,2))-ep.spa(1,3)*
   ep.spb(4,3))))+(pow(ep.spa(2,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,4)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,0)))/(ep.s(0,4,5)*
ep.spa(1,2)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(3,5),2)*pow(ep.spb(1,0),2)*
(ep.spa(3,4)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(5,2)))/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(2,1)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))))
); }

template <class T> complex<T> A2q1_2q2_2q317050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*ep.spb(2,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2)*
(ep.spa(2,5)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(4,3)))/
 (ep.s(2,3,4)*ep.spa(0,5)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))+(pow(ep.spa(1,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,2)*
  ep.spb(2,0)+ep.spa(1,3)*
  ep.spb(3,0)))/(ep.s(1,2,3)*
ep.spa(1,2)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q317525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))+
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q317910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q317955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q318040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q318060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q318075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q318090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q319505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2)*
(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0)))/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
(ep.spa(1,3)*ep.spb(4,1)+
 ep.spa(2,3)*ep.spb(4,2)))/
 (ep.s(1,2,3)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q319510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2)*
 (-(ep.spa(3,5)*ep.spb(3,0))-
  ep.spa(4,5)*ep.spb(4,0)))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))-(pow(ep.spa(0,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*ep.spb(3,2)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))+
pow(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4),3)/
 (ep.s(0,4,5)*ep.spa(2,3)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q319595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q319615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(1,0),2))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q319630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q319645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
ep.spa(4,5)*ep.spb(2,1)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q319685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0),3)/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2)*(ep.spa(0,1)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(5,2)))/(ep.s(0,1,5)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4)))/
 (ep.s(0,4,5)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q319690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spb(3,2)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
pow(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1),3)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(0,2),2)*
pow(ep.spb(4,3),2)*
(-(ep.spa(0,1)*ep.spb(5,1))-
 ep.spa(0,2)*ep.spb(5,2)))/
 (ep.s(0,1,2)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q319805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q319830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
 ep.spa(4,5)*ep.spb(4,1),2)/
 (ep.s(3,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0)))-
pow(ep.spa(0,2)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q319840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q319860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2),2)/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))+
pow(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q320315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(1,0),2))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q320335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
ep.spa(4,5)*ep.spb(2,1)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q320370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2)*
(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0)))/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
(ep.spa(1,3)*ep.spb(4,1)+
 ep.spa(2,3)*ep.spb(4,2)))/
 (ep.s(1,2,3)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q320490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2)*
 (-(ep.spa(3,5)*ep.spb(3,0))-
  ep.spa(4,5)*ep.spb(4,0)))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))-(pow(ep.spa(0,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*ep.spb(3,2)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))+
pow(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4),3)/
 (ep.s(0,4,5)*ep.spa(2,3)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q320530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
 ep.spa(4,5)*ep.spb(4,1),2)/
 (ep.s(3,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0)))-
pow(ep.spa(0,2)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q320545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2),2)/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))+
pow(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q320580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0),3)/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2)*(ep.spa(0,1)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(5,2)))/(ep.s(0,1,5)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4)))/
 (ep.s(0,4,5)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q320670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spb(3,2)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
pow(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1),3)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(0,2),2)*
pow(ep.spb(4,3),2)*
(-(ep.spa(0,1)*ep.spb(5,1))-
 ep.spa(0,2)*ep.spb(5,2)))/
 (ep.s(0,1,2)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q320795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1),2)/
 (ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q320855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q320985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
  ep.spa(1,5)*ep.spb(3,1),2)/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4))))+
pow(ep.spa(2,4)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q321065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(4,5),2)*
pow(ep.spb(3,1),2)*(ep.spa(0,2)*
  ep.spb(2,1)+ep.spa(0,3)*
  ep.spb(3,1)))/(ep.s(1,2,3)*
ep.spa(0,5)*ep.spb(2,1)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))+
pow(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3),3)/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(2,4),2)*
pow(ep.spb(1,0),2)*(ep.spa(2,4)*
  ep.spb(5,2)+ep.spa(3,4)*
  ep.spb(5,3)))/(ep.s(2,3,4)*
ep.spa(3,4)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*ep.spb(5,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q321835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spb(2,1)*
 (ep.spa(0,5)*ep.spb(2,0)+
  ep.spa(1,5)*ep.spb(2,1))*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(1,5),2)*pow(ep.spb(3,2),2)*
(ep.spa(0,5)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(4,1)))/
 (ep.s(0,1,5)*ep.spa(0,5)*
(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1))*
ep.spb(4,3)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,0),3)/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spb(5,0)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q321870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q321925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
 ep.spa(3,4)*ep.spb(3,0),2)/
 (ep.s(2,3,4)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,3)+
 ep.spa(1,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q321930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1),2)/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))+
pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q322090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))+
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q322150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q322360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
  ep.spa(3,5)*ep.spb(4,3),3)/
(ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)*(-(ep.spa(1,2)*
ep.spb(4,2))-ep.spa(1,3)*
   ep.spb(4,3))))+(pow(ep.spa(2,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,4)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,0)))/(ep.s(0,4,5)*
ep.spa(1,2)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(3,5),2)*pow(ep.spb(1,0),2)*
(ep.spa(3,4)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(5,2)))/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(2,1)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))))
); }

template <class T> complex<T> A2q1_2q2_2q322915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q322935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*ep.spb(2,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2)*
(ep.spa(2,5)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(4,3)))/
 (ep.s(2,3,4)*ep.spa(0,5)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))+(pow(ep.spa(1,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,2)*
  ep.spb(2,0)+ep.spa(1,3)*
  ep.spb(3,0)))/(ep.s(1,2,3)*
ep.spa(1,2)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q322950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q323005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2)))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q323010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))+
(pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q323645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(2,4),2)*
 pow(ep.spb(5,0),2))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))+
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q323650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q323705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q323720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2)*(ep.spa(0,2)*
   ep.spb(3,0)+ep.spa(1,2)*
   ep.spb(3,1)))/(ep.s(0,1,2)*
 ep.spa(1,2)*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))*ep.spb(4,3)*
 (ep.spa(0,2)*ep.spb(5,0)+
  ep.spa(1,2)*ep.spb(5,1))))+
(pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
(ep.spa(0,2)*ep.spb(2,1)+
 ep.spa(0,3)*ep.spb(3,1)))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spb(2,1)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2)*(ep.spa(0,4)*
  ep.spb(5,0)+ep.spa(1,4)*
  ep.spb(5,1)))/(ep.s(0,1,5)*
ep.spa(3,4)*ep.spb(5,0)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q323740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q323750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1),3)/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
pow(ep.spb(3,2),2)*(ep.spa(0,4)*
  ep.spb(4,1)+ep.spa(0,5)*
  ep.spb(5,1)))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3)))/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))))
); }

template <class T> complex<T> A2q1_2q2_2q323825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q323830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q323915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q323935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q323950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q323965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(2,1),2))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1)))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q324385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q324565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0),2)/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q324730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
 ep.spa(3,5)*ep.spb(5,1),2)/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2)))-
pow(-(ep.spa(0,2)*ep.spb(4,2))-
 ep.spa(0,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q324785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
  ep.spa(0,3)*ep.spb(3,1),3)/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spb(2,1)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))))-(pow(ep.spa(0,2),2)*
pow(ep.spb(5,4),2)*(ep.spa(0,2)*
  ep.spb(3,0)+ep.spa(1,2)*
  ep.spb(3,1)))/(ep.s(0,1,2)*
ep.spa(1,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*ep.spb(4,3)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1)))+
(pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,4)*ep.spb(5,0)+
 ep.spa(1,4)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(3,4)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(5,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q324820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q324830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(2,1),3)/
 (ep.s(0,1,2)*ep.spa(3,4)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(2,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))+(pow(ep.spa(0,1),2)*
pow(ep.spb(4,2),2)*(ep.spa(2,5)*
  ep.spb(4,2)+ep.spa(3,5)*
  ep.spb(4,3)))/(ep.s(2,3,4)*
ep.spa(0,5)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(4,3)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))+
(pow(ep.spa(1,3),2)*pow(ep.spb(5,4),2)*
(ep.spa(1,2)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(3,0)))/
 (ep.s(1,2,3)*ep.spa(1,2)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q325085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
  ep.spa(0,5)*ep.spb(5,2),2)/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
pow(ep.spa(0,3)*ep.spb(5,0)+
 ep.spa(1,3)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q325680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-
pow(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q325860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q325985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2)*
 (-(ep.spa(1,5)*ep.spb(1,0))-
  ep.spa(2,5)*ep.spb(2,0)))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(ep.spa(1,3)*ep.spb(3,2)+
 ep.spa(1,4)*ep.spb(4,2),3)/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spb(3,2)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))-(pow(ep.spa(1,3),2)*
pow(ep.spb(5,0),2)*(ep.spa(1,3)*
  ep.spb(4,1)+ep.spa(2,3)*
  ep.spb(4,2)))/(ep.s(1,2,3)*
ep.spa(2,3)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q325990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(2,1),2)*(ep.spa(0,4)*
   ep.spb(3,0)+ep.spa(4,5)*
   ep.spb(5,3)))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+
pow(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4),3)/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q326075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q326095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1),2)/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))))+
pow(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q326110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q326125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
 ep.spa(4,5)*ep.spb(5,2),2)/
 (ep.s(0,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))-
pow(-(ep.spa(1,3)*ep.spb(5,3))-
 ep.spa(1,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q326165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(1,4)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(3,2),3)/
 (ep.s(1,2,3)*ep.spa(4,5)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(3,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,0),2)*(ep.spa(2,3)*
  ep.spb(3,1)+ep.spa(2,4)*
  ep.spb(4,1)))/(ep.s(2,3,4)*
ep.spa(2,3)*ep.spb(1,0)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q326170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(2,3)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(4,1)))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spb(1,0)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(0,2),2)*
pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q326285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q326310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
ep.spa(4,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q326320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q326340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,2),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))+
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q326795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1),2)/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))))+
pow(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q326815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q326825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
 ep.spa(4,5)*ep.spb(5,2),2)/
 (ep.s(0,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))-
pow(-(ep.spa(1,3)*ep.spb(5,3))-
 ep.spa(1,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q326850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q326965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2)*
 (-(ep.spa(1,5)*ep.spb(1,0))-
  ep.spa(2,5)*ep.spb(2,0)))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(ep.spa(1,3)*ep.spb(3,2)+
 ep.spa(1,4)*ep.spb(4,2),3)/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spb(3,2)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))-(pow(ep.spa(1,3),2)*
pow(ep.spb(5,0),2)*(ep.spa(1,3)*
  ep.spb(4,1)+ep.spa(2,3)*
  ep.spb(4,2)))/(ep.s(1,2,3)*
ep.spa(2,3)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q326970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(2,1),2)*(ep.spa(0,4)*
   ep.spb(3,0)+ep.spa(4,5)*
   ep.spb(5,3)))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+
pow(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4),3)/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q327010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
ep.spa(4,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q327025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q327040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,2),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))+
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q327060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q327145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(1,4)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(3,2),3)/
 (ep.s(1,2,3)*ep.spa(4,5)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(3,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,0),2)*(ep.spa(2,3)*
  ep.spb(3,1)+ep.spa(2,4)*
  ep.spb(4,1)))/(ep.s(2,3,4)*
ep.spa(2,3)*ep.spb(1,0)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q327150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(2,3)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(4,1)))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spb(1,0)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(0,2),2)*
pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q328565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(2,1),2))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1)))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q328595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q328805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q328985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2)*(ep.spa(0,2)*
   ep.spb(3,0)+ep.spa(1,2)*
   ep.spb(3,1)))/(ep.s(0,1,2)*
 ep.spa(1,2)*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))*ep.spb(4,3)*
 (ep.spa(0,2)*ep.spb(5,0)+
  ep.spa(1,2)*ep.spb(5,1))))+
(pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
(ep.spa(0,2)*ep.spb(2,1)+
 ep.spa(0,3)*ep.spb(3,1)))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spb(2,1)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2)*(ep.spa(0,4)*
  ep.spb(5,0)+ep.spa(1,4)*
  ep.spb(5,1)))/(ep.s(0,1,5)*
ep.spa(3,4)*ep.spb(5,0)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q329605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1),3)/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
pow(ep.spb(3,2),2)*(ep.spa(0,4)*
  ep.spb(4,1)+ep.spa(0,5)*
  ep.spb(5,1)))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3)))/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))))
); }

template <class T> complex<T> A2q1_2q2_2q329640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(2,4),2)*
 pow(ep.spb(5,0),2))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))+
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q329670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q329860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
  ep.spa(0,5)*ep.spb(5,2),2)/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
pow(ep.spa(0,3)*ep.spb(5,0)+
 ep.spa(1,3)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q329890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q329965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-
pow(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q330100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
  ep.spa(0,3)*ep.spb(3,1),3)/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spb(2,1)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))))-(pow(ep.spa(0,2),2)*
pow(ep.spb(5,4),2)*(ep.spa(0,2)*
  ep.spb(3,0)+ep.spa(1,2)*
  ep.spb(3,1)))/(ep.s(0,1,2)*
ep.spa(1,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*ep.spb(4,3)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1)))+
(pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,4)*ep.spb(5,0)+
 ep.spa(1,4)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(3,4)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(5,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q330685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(2,1),3)/
 (ep.s(0,1,2)*ep.spa(3,4)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(2,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))+(pow(ep.spa(0,1),2)*
pow(ep.spb(4,2),2)*(ep.spa(2,5)*
  ep.spb(4,2)+ep.spa(3,5)*
  ep.spb(4,3)))/(ep.s(2,3,4)*
ep.spa(0,5)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(4,3)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))+
(pow(ep.spa(1,3),2)*pow(ep.spb(5,4),2)*
(ep.spa(1,2)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(3,0)))/
 (ep.s(1,2,3)*ep.spa(1,2)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q330720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q330745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0),2)/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q330750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
 ep.spa(3,5)*ep.spb(5,1),2)/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2)))-
pow(-(ep.spa(0,2)*ep.spb(4,2))-
 ep.spa(0,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q331415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(4,5),2)*
pow(ep.spb(3,1),2)*(ep.spa(0,2)*
  ep.spb(2,1)+ep.spa(0,3)*
  ep.spb(3,1)))/(ep.s(1,2,3)*
ep.spa(0,5)*ep.spb(2,1)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))+
pow(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3),3)/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(2,4),2)*
pow(ep.spb(1,0),2)*(ep.spa(2,4)*
  ep.spb(5,2)+ep.spa(3,4)*
  ep.spb(5,3)))/(ep.s(2,3,4)*
ep.spa(3,4)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*ep.spb(5,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q331425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spb(2,1)*
 (ep.spa(0,5)*ep.spb(2,0)+
  ep.spa(1,5)*ep.spb(2,1))*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(1,5),2)*pow(ep.spb(3,2),2)*
(ep.spa(0,5)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(4,1)))/
 (ep.s(0,1,5)*ep.spa(0,5)*
(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1))*
ep.spb(4,3)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,0),3)/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spb(5,0)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q331460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
 ep.spa(3,4)*ep.spb(3,0),2)/
 (ep.s(2,3,4)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,3)+
 ep.spa(1,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q331520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1),2)/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))+
pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q331595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1),2)/
 (ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q331675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
  ep.spa(1,5)*ep.spb(3,1),2)/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4))))+
pow(ep.spa(2,4)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q331855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q331945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
  ep.spa(3,5)*ep.spb(4,3),3)/
(ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)*(-(ep.spa(1,2)*
ep.spb(4,2))-ep.spa(1,3)*
   ep.spb(4,3))))+(pow(ep.spa(2,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,4)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,0)))/(ep.s(0,4,5)*
ep.spa(1,2)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(3,5),2)*pow(ep.spb(1,0),2)*
(ep.spa(3,4)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(5,2)))/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(2,1)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))))
); }

template <class T> complex<T> A2q1_2q2_2q332505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*ep.spb(2,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2)*
(ep.spa(2,5)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(4,3)))/
 (ep.s(2,3,4)*ep.spa(0,5)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))+(pow(ep.spa(1,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,2)*
  ep.spb(2,0)+ep.spa(1,3)*
  ep.spb(3,0)))/(ep.s(1,2,3)*
ep.spa(1,2)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q332540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2)))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q332600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))+
(pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q332855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q332945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))+
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q332970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q333150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1),2)/
 (ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q333765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q333935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
  ep.spa(1,5)*ep.spb(3,1),2)/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4))))+
pow(ep.spa(2,4)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,0),2)/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q333945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(4,5),2)*
pow(ep.spb(3,1),2)*(ep.spa(0,2)*
  ep.spb(2,1)+ep.spa(0,3)*
  ep.spb(3,1)))/(ep.s(1,2,3)*
ep.spa(0,5)*ep.spb(2,1)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))+
pow(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3),3)/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(2,4),2)*
pow(ep.spb(1,0),2)*(ep.spa(2,4)*
  ep.spb(5,2)+ep.spa(3,4)*
  ep.spb(5,3)))/(ep.s(2,3,4)*
ep.spa(3,4)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*ep.spb(5,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q334375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spb(2,1)*
 (ep.spa(0,5)*ep.spb(2,0)+
  ep.spa(1,5)*ep.spb(2,1))*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(1,5),2)*pow(ep.spb(3,2),2)*
(ep.spa(0,5)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(4,1)))/
 (ep.s(0,1,5)*ep.spa(0,5)*
(ep.spa(0,5)*ep.spb(2,0)+
 ep.spa(1,5)*ep.spb(2,1))*
ep.spb(4,3)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,0),3)/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spb(5,0)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q334410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
 ep.spa(3,4)*ep.spb(3,0),2)/
 (ep.s(2,3,4)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))-
pow(ep.spa(1,4)*ep.spb(4,3)+
 ep.spa(1,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q334530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1),2)/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))+
pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q334785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q334890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))+
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q335060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q335240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
  ep.spa(3,5)*ep.spb(4,3),3)/
(ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)*(-(ep.spa(1,2)*
ep.spb(4,2))-ep.spa(1,3)*
   ep.spb(4,3))))+(pow(ep.spa(2,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,4)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,0)))/(ep.s(0,4,5)*
ep.spa(1,2)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(3,5),2)*pow(ep.spb(1,0),2)*
(ep.spa(3,4)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(5,2)))/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(2,1)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0))))
); }

template <class T> complex<T> A2q1_2q2_2q335455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2)*(ep.spa(0,3)*
   ep.spb(2,0)+ep.spa(1,3)*
   ep.spb(2,1)))/(ep.s(0,1,2)*
 ep.spa(3,4)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*ep.spb(2,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2)*
(ep.spa(2,5)*ep.spb(4,2)+
 ep.spa(3,5)*ep.spb(4,3)))/
 (ep.s(2,3,4)*ep.spa(0,5)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))+(pow(ep.spa(1,3),2)*
pow(ep.spb(4,0),2)*(ep.spa(1,2)*
  ep.spb(2,0)+ep.spa(1,3)*
  ep.spb(3,0)))/(ep.s(1,2,3)*
ep.spa(1,2)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q335490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q335605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2)))-
(pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q335610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))))+
(pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q336080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q336085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q336110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q336120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q336145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q336150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q337635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2)*
(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0)))/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
(ep.spa(1,3)*ep.spb(4,1)+
 ep.spa(2,3)*ep.spb(4,2)))/
 (ep.s(1,2,3)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q337640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2)*
 (-(ep.spa(3,5)*ep.spb(3,0))-
  ep.spa(4,5)*ep.spb(4,0)))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))-(pow(ep.spa(0,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*ep.spb(3,2)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))+
pow(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4),3)/
 (ep.s(0,4,5)*ep.spa(2,3)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q337665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q337675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(1,0),2))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q337700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q337705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
ep.spa(4,5)*ep.spb(2,1)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q337815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0),3)/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2)*(ep.spa(0,1)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(5,2)))/(ep.s(0,1,5)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4)))/
 (ep.s(0,4,5)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q337820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spb(3,2)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
pow(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1),3)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(0,2),2)*
pow(ep.spb(4,3),2)*
(-(ep.spa(0,1)*ep.spb(5,1))-
 ep.spa(0,2)*ep.spb(5,2)))/
 (ep.s(0,1,2)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q337875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q337890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
 ep.spa(4,5)*ep.spb(4,1),2)/
 (ep.s(3,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0)))-
pow(ep.spa(0,2)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q337910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q337920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2),2)/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))+
pow(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q338025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(1,0),2))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0))))+
(pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q338035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q338055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
ep.spa(4,5)*ep.spb(2,1)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q338070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q338125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
pow(ep.spb(2,0),2)*
(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0)))/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))-
(pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
(ep.spa(1,3)*ep.spb(4,1)+
 ep.spa(2,3)*ep.spb(4,2)))/
 (ep.s(1,2,3)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q338130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,1),2)*
 (-(ep.spa(3,5)*ep.spb(3,0))-
  ep.spa(4,5)*ep.spb(4,0)))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))-(pow(ep.spa(0,5),2)*
pow(ep.spb(4,2),2)*(ep.spa(1,3)*
  ep.spb(3,2)+ep.spa(1,4)*
  ep.spb(4,2)))/(ep.s(2,3,4)*
ep.spa(0,1)*ep.spb(3,2)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))+
pow(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4),3)/
 (ep.s(0,4,5)*ep.spa(2,3)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q338240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
 ep.spa(4,5)*ep.spb(4,1),2)/
 (ep.s(3,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(3,4)*
  ep.spb(4,0)+ep.spa(3,5)*
  ep.spb(5,0)))-
pow(ep.spa(0,2)*ep.spb(4,0)+
 ep.spa(2,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q338245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q338270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2),2)/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))))+
pow(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4),2)/
 (ep.s(0,4,5)*ep.spa(1,2)*
(ep.spa(3,4)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,0))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q338280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q338305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(-(ep.spa(1,5)*ep.spb(1,0))-
 ep.spa(2,5)*ep.spb(2,0),3)/
 (ep.s(0,1,2)*ep.spa(4,5)*
ep.spb(1,0)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))-(pow(ep.spa(1,5),2)*
pow(ep.spb(4,3),2)*(ep.spa(0,1)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(5,2)))/(ep.s(0,1,5)*
ep.spa(0,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(3,2)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
(ep.spa(0,3)*ep.spb(4,0)+
 ep.spa(3,5)*ep.spb(5,4)))/
 (ep.s(0,4,5)*ep.spa(2,3)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(5,4)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q338310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spb(3,2)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
pow(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1),3)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(0,2),2)*
pow(ep.spb(4,3),2)*
(-(ep.spa(0,1)*ep.spb(5,1))-
 ep.spa(0,2)*ep.spb(5,2)))/
 (ep.s(0,1,2)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q339190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2)*(ep.spa(0,2)*
   ep.spb(3,0)+ep.spa(1,2)*
   ep.spb(3,1)))/(ep.s(0,1,2)*
 ep.spa(1,2)*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))*ep.spb(4,3)*
 (ep.spa(0,2)*ep.spb(5,0)+
  ep.spa(1,2)*ep.spb(5,1))))+
(pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
(ep.spa(0,2)*ep.spb(2,1)+
 ep.spa(0,3)*ep.spb(3,1)))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spb(2,1)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2)*(ep.spa(0,4)*
  ep.spb(5,0)+ep.spa(1,4)*
  ep.spb(5,1)))/(ep.s(0,1,5)*
ep.spa(3,4)*ep.spb(5,0)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q339195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1),3)/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
pow(ep.spb(3,2),2)*(ep.spa(0,4)*
  ep.spb(4,1)+ep.spa(0,5)*
  ep.spb(5,1)))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3)))/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))))
); }

template <class T> complex<T> A2q1_2q2_2q339230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(2,4),2)*
 pow(ep.spb(5,0),2))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))+
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q339260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q339370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(2,1),2))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1)))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q339445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q339625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q339865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
  ep.spa(0,3)*ep.spb(3,1),3)/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spb(2,1)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))))-(pow(ep.spa(0,2),2)*
pow(ep.spb(5,4),2)*(ep.spa(0,2)*
  ep.spb(3,0)+ep.spa(1,2)*
  ep.spb(3,1)))/(ep.s(0,1,2)*
ep.spa(1,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*ep.spb(4,3)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1)))+
(pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,4)*ep.spb(5,0)+
 ep.spa(1,4)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(3,4)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(5,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q340275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(2,1),3)/
 (ep.s(0,1,2)*ep.spa(3,4)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(2,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))+(pow(ep.spa(0,1),2)*
pow(ep.spb(4,2),2)*(ep.spa(2,5)*
  ep.spb(4,2)+ep.spa(3,5)*
  ep.spb(4,3)))/(ep.s(2,3,4)*
ep.spa(0,5)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(4,3)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))+
(pow(ep.spa(1,3),2)*pow(ep.spb(5,4),2)*
(ep.spa(1,2)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(3,0)))/
 (ep.s(1,2,3)*ep.spa(1,2)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q340310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0),2)/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q340340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
 ep.spa(3,5)*ep.spb(5,1),2)/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2)))-
pow(-(ep.spa(0,2)*ep.spb(4,2))-
 ep.spa(0,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q340630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
  ep.spa(0,5)*ep.spb(5,2),2)/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
pow(ep.spa(0,3)*ep.spb(5,0)+
 ep.spa(1,3)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q340740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-
pow(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q340920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q340980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(2,1),2))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(3,2)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1)))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q341535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q341715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q341850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2)*(ep.spa(0,2)*
   ep.spb(3,0)+ep.spa(1,2)*
   ep.spb(3,1)))/(ep.s(0,1,2)*
 ep.spa(1,2)*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))*ep.spb(4,3)*
 (ep.spa(0,2)*ep.spb(5,0)+
  ep.spa(1,2)*ep.spb(5,1))))+
(pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
(ep.spa(0,2)*ep.spb(2,1)+
 ep.spa(0,3)*ep.spb(3,1)))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spb(2,1)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2)*(ep.spa(0,4)*
  ep.spb(5,0)+ep.spa(1,4)*
  ep.spb(5,1)))/(ep.s(0,1,5)*
ep.spa(3,4)*ep.spb(5,0)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q342145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1),3)/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
pow(ep.spb(3,2),2)*(ep.spa(0,4)*
  ep.spb(4,1)+ep.spa(0,5)*
  ep.spb(5,1)))/(ep.s(0,4,5)*
ep.spa(0,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(1,0)+
 ep.spa(4,5)*ep.spb(5,1))*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
(ep.spa(2,4)*ep.spb(4,3)+
 ep.spa(2,5)*ep.spb(5,3)))/
 (ep.s(3,4,5)*ep.spa(1,2)*
ep.spb(4,3)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3))))
); }

template <class T> complex<T> A2q1_2q2_2q342180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(2,4),2)*
 pow(ep.spb(5,0),2))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))+
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q342270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(2,4),2)*
pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
ep.spa(3,4)*ep.spb(1,0)*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spb(4,3)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q342345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
  ep.spa(0,5)*ep.spb(5,2),2)/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
pow(ep.spa(0,3)*ep.spb(5,0)+
 ep.spa(1,3)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q342830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q342925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3),2)/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1)))-
pow(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1),2)/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(5,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q343010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
  ep.spa(0,3)*ep.spb(3,1),3)/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spb(2,1)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*(-(ep.spa(0,1)*
ep.spb(3,1))-ep.spa(0,2)*
   ep.spb(3,2))))-(pow(ep.spa(0,2),2)*
pow(ep.spb(5,4),2)*(ep.spa(0,2)*
  ep.spb(3,0)+ep.spa(1,2)*
  ep.spb(3,1)))/(ep.s(0,1,2)*
ep.spa(1,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2))*ep.spb(4,3)*
(ep.spa(0,2)*ep.spb(5,0)+
 ep.spa(1,2)*ep.spb(5,1)))+
(pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,4)*ep.spb(5,0)+
 ep.spa(1,4)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(3,4)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(5,0)*(ep.spa(0,2)*
  ep.spb(5,0)+ep.spa(1,2)*
  ep.spb(5,1))))
); }

template <class T> complex<T> A2q1_2q2_2q343225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(2,1),3)/
 (ep.s(0,1,2)*ep.spa(3,4)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
ep.spb(2,1)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1)))+(pow(ep.spa(0,1),2)*
pow(ep.spb(4,2),2)*(ep.spa(2,5)*
  ep.spb(4,2)+ep.spa(3,5)*
  ep.spb(4,3)))/(ep.s(2,3,4)*
ep.spa(0,5)*(ep.spa(0,5)*
  ep.spb(2,0)+ep.spa(1,5)*
  ep.spb(2,1))*ep.spb(4,3)*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3)))+
(pow(ep.spa(1,3),2)*pow(ep.spb(5,4),2)*
(ep.spa(1,2)*ep.spb(2,0)+
 ep.spa(1,3)*ep.spb(3,0)))/
 (ep.s(1,2,3)*ep.spa(1,2)*
(-(ep.spa(1,3)*ep.spb(1,0))-
 ep.spa(2,3)*ep.spb(2,0))*
(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3))*
ep.spb(5,0)))
); }

template <class T> complex<T> A2q1_2q2_2q343260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
  ep.spa(3,5)*ep.spb(5,0),2)/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(-(ep.spa(1,2)*ep.spb(4,2))-
 ep.spa(1,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q343350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
 ep.spa(3,5)*ep.spb(5,1),2)/
 (ep.s(3,4,5)*ep.spa(3,4)*
ep.spb(1,0)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2)))-
pow(-(ep.spa(0,2)*ep.spb(4,2))-
 ep.spa(0,3)*ep.spb(4,3),2)/
 (ep.s(2,3,4)*ep.spa(0,1)*
(-(ep.spa(3,5)*ep.spb(3,2))-
 ep.spa(4,5)*ep.spb(4,2))*
ep.spb(4,3)))
); }

template <class T> complex<T> A2q1_2q2_2q343640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q343710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2)*
 (-(ep.spa(1,5)*ep.spb(1,0))-
  ep.spa(2,5)*ep.spb(2,0)))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(ep.spa(1,3)*ep.spb(3,2)+
 ep.spa(1,4)*ep.spb(4,2),3)/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spb(3,2)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))-(pow(ep.spa(1,3),2)*
pow(ep.spb(5,0),2)*(ep.spa(1,3)*
  ep.spb(4,1)+ep.spa(2,3)*
  ep.spb(4,2)))/(ep.s(1,2,3)*
ep.spa(2,3)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(2,1),2)*(ep.spa(0,4)*
   ep.spb(3,0)+ep.spa(4,5)*
   ep.spb(5,3)))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+
pow(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4),3)/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q344145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1),2)/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))))+
pow(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
 ep.spa(4,5)*ep.spb(5,2),2)/
 (ep.s(0,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))-
pow(-(ep.spa(1,3)*ep.spb(5,3))-
 ep.spa(1,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(1,4)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(3,2),3)/
 (ep.s(1,2,3)*ep.spa(4,5)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(3,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,0),2)*(ep.spa(2,3)*
  ep.spb(3,1)+ep.spa(2,4)*
  ep.spb(4,1)))/(ep.s(2,3,4)*
ep.spa(2,3)*ep.spb(1,0)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q344300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(2,3)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(4,1)))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spb(1,0)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(0,2),2)*
pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q344355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
ep.spa(4,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,2),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))+
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1),2)/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))))+
pow(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
 ep.spa(4,5)*ep.spb(5,2),2)/
 (ep.s(0,4,5)*ep.spa(4,5)*
ep.spb(2,1)*(ep.spa(0,4)*
  ep.spb(4,3)+ep.spa(0,5)*
  ep.spb(5,3)))-
pow(-(ep.spa(1,3)*ep.spb(5,3))-
 ep.spa(1,4)*ep.spb(5,4),2)/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2)*
 (-(ep.spa(1,5)*ep.spb(1,0))-
  ep.spa(2,5)*ep.spb(2,0)))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spb(1,0)*(-(ep.spa(1,3)*
ep.spb(1,0))-ep.spa(2,3)*
   ep.spb(2,0))*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))))+
pow(ep.spa(1,3)*ep.spb(3,2)+
 ep.spa(1,4)*ep.spb(4,2),3)/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spb(3,2)*(-(ep.spa(3,5)*
   ep.spb(3,2))-ep.spa(4,5)*
  ep.spb(4,2))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3)))-(pow(ep.spa(1,3),2)*
pow(ep.spb(5,0),2)*(ep.spa(1,3)*
  ep.spb(4,1)+ep.spa(2,3)*
  ep.spb(4,2)))/(ep.s(1,2,3)*
ep.spa(2,3)*(-(ep.spa(1,3)*
   ep.spb(1,0))-ep.spa(2,3)*
  ep.spb(2,0))*(-(ep.spa(1,2)*
   ep.spb(4,2))-ep.spa(1,3)*
  ep.spb(4,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(2,1),2)*(ep.spa(0,4)*
   ep.spb(3,0)+ep.spa(4,5)*
   ep.spb(5,3)))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(0,2)*ep.spb(1,0)+
 ep.spa(2,5)*ep.spb(5,1)))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spb(1,0)*(ep.spa(0,4)*
  ep.spb(1,0)+ep.spa(4,5)*
  ep.spb(5,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+
pow(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4),3)/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q344720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
ep.spa(4,5)*ep.spb(2,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3)))-
(pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,2),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))+
(pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2))/
 (ep.s(3,4,5)*ep.spa(1,2)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)))
); }

template <class T> complex<T> A2q1_2q2_2q344760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1_2q2_2q344785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(1,4)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(3,2),3)/
 (ep.s(1,2,3)*ep.spa(4,5)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
ep.spb(3,2)*(-(ep.spa(0,1)*
   ep.spb(3,1))-ep.spa(0,2)*
  ep.spb(3,2)))-(pow(ep.spa(2,4),2)*
pow(ep.spb(5,0),2)*(ep.spa(2,3)*
  ep.spb(3,1)+ep.spa(2,4)*
  ep.spb(4,1)))/(ep.s(2,3,4)*
ep.spa(2,3)*ep.spb(1,0)*
(-(ep.spa(2,4)*ep.spb(2,1))-
 ep.spa(3,4)*ep.spb(3,1))*
(-(ep.spa(2,3)*ep.spb(5,3))-
 ep.spa(2,4)*ep.spb(5,4)))-
(pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(-(ep.spa(0,1)*ep.spb(3,1))-
 ep.spa(0,2)*ep.spb(3,2))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))
); }

template <class T> complex<T> A2q1_2q2_2q344790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(1,4)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(3,2)))/(ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.spa(2,4)*
ep.spb(2,1))-ep.spa(3,4)*
   ep.spb(3,1))*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))))-
(pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
(ep.spa(2,3)*ep.spb(3,1)+
 ep.spa(2,4)*ep.spb(4,1)))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spb(1,0)*(-(ep.spa(2,4)*
   ep.spb(2,1))-ep.spa(3,4)*
  ep.spb(3,1))*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4)))+(pow(ep.spa(0,2),2)*
pow(ep.spb(5,3),2)*
(-(ep.spa(0,3)*ep.spb(5,3))-
 ep.spa(0,4)*ep.spb(5,4)))/
 (ep.s(3,4,5)*ep.spa(0,1)*
(ep.spa(0,4)*ep.spb(4,3)+
 ep.spa(0,5)*ep.spb(5,3))*
ep.spb(5,4)*(-(ep.spa(2,3)*
   ep.spb(5,3))-ep.spa(2,4)*
  ep.spb(5,4))))

 ); }


template <class T> complex<T>  (*A2q1_2q2_2q3_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 1865 : return &A2q1_2q2_2q31865_eval;
case 1870 : return &A2q1_2q2_2q31870_eval;
case 1895 : return &A2q1_2q2_2q31895_eval;
case 1905 : return &A2q1_2q2_2q31905_eval;
case 1930 : return &A2q1_2q2_2q31930_eval;
case 1935 : return &A2q1_2q2_2q31935_eval;
case 2045 : return &A2q1_2q2_2q32045_eval;
case 2050 : return &A2q1_2q2_2q32050_eval;
case 2105 : return &A2q1_2q2_2q32105_eval;
case 2120 : return &A2q1_2q2_2q32120_eval;
case 2140 : return &A2q1_2q2_2q32140_eval;
case 2150 : return &A2q1_2q2_2q32150_eval;
case 2255 : return &A2q1_2q2_2q32255_eval;
case 2265 : return &A2q1_2q2_2q32265_eval;
case 2285 : return &A2q1_2q2_2q32285_eval;
case 2300 : return &A2q1_2q2_2q32300_eval;
case 2355 : return &A2q1_2q2_2q32355_eval;
case 2360 : return &A2q1_2q2_2q32360_eval;
case 2470 : return &A2q1_2q2_2q32470_eval;
case 2475 : return &A2q1_2q2_2q32475_eval;
case 2500 : return &A2q1_2q2_2q32500_eval;
case 2510 : return &A2q1_2q2_2q32510_eval;
case 2535 : return &A2q1_2q2_2q32535_eval;
case 2540 : return &A2q1_2q2_2q32540_eval;
case 2945 : return &A2q1_2q2_2q32945_eval;
case 2950 : return &A2q1_2q2_2q32950_eval;
case 2975 : return &A2q1_2q2_2q32975_eval;
case 2985 : return &A2q1_2q2_2q32985_eval;
case 3010 : return &A2q1_2q2_2q33010_eval;
case 3015 : return &A2q1_2q2_2q33015_eval;
case 3305 : return &A2q1_2q2_2q33305_eval;
case 3310 : return &A2q1_2q2_2q33310_eval;
case 3395 : return &A2q1_2q2_2q33395_eval;
case 3415 : return &A2q1_2q2_2q33415_eval;
case 3430 : return &A2q1_2q2_2q33430_eval;
case 3445 : return &A2q1_2q2_2q33445_eval;
case 3515 : return &A2q1_2q2_2q33515_eval;
case 3525 : return &A2q1_2q2_2q33525_eval;
case 3575 : return &A2q1_2q2_2q33575_eval;
case 3595 : return &A2q1_2q2_2q33595_eval;
case 3645 : return &A2q1_2q2_2q33645_eval;
case 3655 : return &A2q1_2q2_2q33655_eval;
case 3730 : return &A2q1_2q2_2q33730_eval;
case 3735 : return &A2q1_2q2_2q33735_eval;
case 3790 : return &A2q1_2q2_2q33790_eval;
case 3805 : return &A2q1_2q2_2q33805_eval;
case 3825 : return &A2q1_2q2_2q33825_eval;
case 3835 : return &A2q1_2q2_2q33835_eval;
case 4205 : return &A2q1_2q2_2q34205_eval;
case 4210 : return &A2q1_2q2_2q34210_eval;
case 4265 : return &A2q1_2q2_2q34265_eval;
case 4280 : return &A2q1_2q2_2q34280_eval;
case 4300 : return &A2q1_2q2_2q34300_eval;
case 4310 : return &A2q1_2q2_2q34310_eval;
case 4385 : return &A2q1_2q2_2q34385_eval;
case 4390 : return &A2q1_2q2_2q34390_eval;
case 4475 : return &A2q1_2q2_2q34475_eval;
case 4495 : return &A2q1_2q2_2q34495_eval;
case 4510 : return &A2q1_2q2_2q34510_eval;
case 4525 : return &A2q1_2q2_2q34525_eval;
case 4805 : return &A2q1_2q2_2q34805_eval;
case 4820 : return &A2q1_2q2_2q34820_eval;
case 4835 : return &A2q1_2q2_2q34835_eval;
case 4855 : return &A2q1_2q2_2q34855_eval;
case 4940 : return &A2q1_2q2_2q34940_eval;
case 4945 : return &A2q1_2q2_2q34945_eval;
case 5020 : return &A2q1_2q2_2q35020_eval;
case 5030 : return &A2q1_2q2_2q35030_eval;
case 5050 : return &A2q1_2q2_2q35050_eval;
case 5065 : return &A2q1_2q2_2q35065_eval;
case 5120 : return &A2q1_2q2_2q35120_eval;
case 5125 : return &A2q1_2q2_2q35125_eval;
case 5495 : return &A2q1_2q2_2q35495_eval;
case 5505 : return &A2q1_2q2_2q35505_eval;
case 5525 : return &A2q1_2q2_2q35525_eval;
case 5540 : return &A2q1_2q2_2q35540_eval;
case 5595 : return &A2q1_2q2_2q35595_eval;
case 5600 : return &A2q1_2q2_2q35600_eval;
case 5675 : return &A2q1_2q2_2q35675_eval;
case 5685 : return &A2q1_2q2_2q35685_eval;
case 5735 : return &A2q1_2q2_2q35735_eval;
case 5755 : return &A2q1_2q2_2q35755_eval;
case 5805 : return &A2q1_2q2_2q35805_eval;
case 5815 : return &A2q1_2q2_2q35815_eval;
case 5885 : return &A2q1_2q2_2q35885_eval;
case 5900 : return &A2q1_2q2_2q35900_eval;
case 5915 : return &A2q1_2q2_2q35915_eval;
case 5935 : return &A2q1_2q2_2q35935_eval;
case 6020 : return &A2q1_2q2_2q36020_eval;
case 6025 : return &A2q1_2q2_2q36025_eval;
case 6315 : return &A2q1_2q2_2q36315_eval;
case 6320 : return &A2q1_2q2_2q36320_eval;
case 6345 : return &A2q1_2q2_2q36345_eval;
case 6355 : return &A2q1_2q2_2q36355_eval;
case 6380 : return &A2q1_2q2_2q36380_eval;
case 6385 : return &A2q1_2q2_2q36385_eval;
case 6790 : return &A2q1_2q2_2q36790_eval;
case 6795 : return &A2q1_2q2_2q36795_eval;
case 6820 : return &A2q1_2q2_2q36820_eval;
case 6830 : return &A2q1_2q2_2q36830_eval;
case 6855 : return &A2q1_2q2_2q36855_eval;
case 6860 : return &A2q1_2q2_2q36860_eval;
case 6970 : return &A2q1_2q2_2q36970_eval;
case 6975 : return &A2q1_2q2_2q36975_eval;
case 7030 : return &A2q1_2q2_2q37030_eval;
case 7045 : return &A2q1_2q2_2q37045_eval;
case 7065 : return &A2q1_2q2_2q37065_eval;
case 7075 : return &A2q1_2q2_2q37075_eval;
case 7180 : return &A2q1_2q2_2q37180_eval;
case 7190 : return &A2q1_2q2_2q37190_eval;
case 7210 : return &A2q1_2q2_2q37210_eval;
case 7225 : return &A2q1_2q2_2q37225_eval;
case 7280 : return &A2q1_2q2_2q37280_eval;
case 7285 : return &A2q1_2q2_2q37285_eval;
case 7395 : return &A2q1_2q2_2q37395_eval;
case 7400 : return &A2q1_2q2_2q37400_eval;
case 7425 : return &A2q1_2q2_2q37425_eval;
case 7435 : return &A2q1_2q2_2q37435_eval;
case 7460 : return &A2q1_2q2_2q37460_eval;
case 7465 : return &A2q1_2q2_2q37465_eval;
case 8345 : return &A2q1_2q2_2q38345_eval;
case 8350 : return &A2q1_2q2_2q38350_eval;
case 8375 : return &A2q1_2q2_2q38375_eval;
case 8385 : return &A2q1_2q2_2q38385_eval;
case 8410 : return &A2q1_2q2_2q38410_eval;
case 8415 : return &A2q1_2q2_2q38415_eval;
case 8525 : return &A2q1_2q2_2q38525_eval;
case 8530 : return &A2q1_2q2_2q38530_eval;
case 8585 : return &A2q1_2q2_2q38585_eval;
case 8600 : return &A2q1_2q2_2q38600_eval;
case 8620 : return &A2q1_2q2_2q38620_eval;
case 8630 : return &A2q1_2q2_2q38630_eval;
case 8735 : return &A2q1_2q2_2q38735_eval;
case 8745 : return &A2q1_2q2_2q38745_eval;
case 8765 : return &A2q1_2q2_2q38765_eval;
case 8780 : return &A2q1_2q2_2q38780_eval;
case 8835 : return &A2q1_2q2_2q38835_eval;
case 8840 : return &A2q1_2q2_2q38840_eval;
case 8950 : return &A2q1_2q2_2q38950_eval;
case 8955 : return &A2q1_2q2_2q38955_eval;
case 8980 : return &A2q1_2q2_2q38980_eval;
case 8990 : return &A2q1_2q2_2q38990_eval;
case 9015 : return &A2q1_2q2_2q39015_eval;
case 9020 : return &A2q1_2q2_2q39020_eval;
case 10505 : return &A2q1_2q2_2q310505_eval;
case 10510 : return &A2q1_2q2_2q310510_eval;
case 10535 : return &A2q1_2q2_2q310535_eval;
case 10545 : return &A2q1_2q2_2q310545_eval;
case 10570 : return &A2q1_2q2_2q310570_eval;
case 10575 : return &A2q1_2q2_2q310575_eval;
case 11045 : return &A2q1_2q2_2q311045_eval;
case 11050 : return &A2q1_2q2_2q311050_eval;
case 11165 : return &A2q1_2q2_2q311165_eval;
case 11190 : return &A2q1_2q2_2q311190_eval;
case 11200 : return &A2q1_2q2_2q311200_eval;
case 11220 : return &A2q1_2q2_2q311220_eval;
case 11255 : return &A2q1_2q2_2q311255_eval;
case 11265 : return &A2q1_2q2_2q311265_eval;
case 11345 : return &A2q1_2q2_2q311345_eval;
case 11370 : return &A2q1_2q2_2q311370_eval;
case 11415 : return &A2q1_2q2_2q311415_eval;
case 11430 : return &A2q1_2q2_2q311430_eval;
case 11470 : return &A2q1_2q2_2q311470_eval;
case 11475 : return &A2q1_2q2_2q311475_eval;
case 11560 : return &A2q1_2q2_2q311560_eval;
case 11580 : return &A2q1_2q2_2q311580_eval;
case 11595 : return &A2q1_2q2_2q311595_eval;
case 11610 : return &A2q1_2q2_2q311610_eval;
case 11765 : return &A2q1_2q2_2q311765_eval;
case 11770 : return &A2q1_2q2_2q311770_eval;
case 11825 : return &A2q1_2q2_2q311825_eval;
case 11840 : return &A2q1_2q2_2q311840_eval;
case 11860 : return &A2q1_2q2_2q311860_eval;
case 11870 : return &A2q1_2q2_2q311870_eval;
case 12125 : return &A2q1_2q2_2q312125_eval;
case 12130 : return &A2q1_2q2_2q312130_eval;
case 12245 : return &A2q1_2q2_2q312245_eval;
case 12270 : return &A2q1_2q2_2q312270_eval;
case 12280 : return &A2q1_2q2_2q312280_eval;
case 12300 : return &A2q1_2q2_2q312300_eval;
case 12545 : return &A2q1_2q2_2q312545_eval;
case 12560 : return &A2q1_2q2_2q312560_eval;
case 12605 : return &A2q1_2q2_2q312605_eval;
case 12630 : return &A2q1_2q2_2q312630_eval;
case 12710 : return &A2q1_2q2_2q312710_eval;
case 12720 : return &A2q1_2q2_2q312720_eval;
case 12760 : return &A2q1_2q2_2q312760_eval;
case 12770 : return &A2q1_2q2_2q312770_eval;
case 12820 : return &A2q1_2q2_2q312820_eval;
case 12840 : return &A2q1_2q2_2q312840_eval;
case 12890 : return &A2q1_2q2_2q312890_eval;
case 12900 : return &A2q1_2q2_2q312900_eval;
case 13055 : return &A2q1_2q2_2q313055_eval;
case 13065 : return &A2q1_2q2_2q313065_eval;
case 13085 : return &A2q1_2q2_2q313085_eval;
case 13100 : return &A2q1_2q2_2q313100_eval;
case 13155 : return &A2q1_2q2_2q313155_eval;
case 13160 : return &A2q1_2q2_2q313160_eval;
case 13415 : return &A2q1_2q2_2q313415_eval;
case 13425 : return &A2q1_2q2_2q313425_eval;
case 13505 : return &A2q1_2q2_2q313505_eval;
case 13530 : return &A2q1_2q2_2q313530_eval;
case 13575 : return &A2q1_2q2_2q313575_eval;
case 13590 : return &A2q1_2q2_2q313590_eval;
case 13625 : return &A2q1_2q2_2q313625_eval;
case 13640 : return &A2q1_2q2_2q313640_eval;
case 13685 : return &A2q1_2q2_2q313685_eval;
case 13710 : return &A2q1_2q2_2q313710_eval;
case 13790 : return &A2q1_2q2_2q313790_eval;
case 13800 : return &A2q1_2q2_2q313800_eval;
case 14055 : return &A2q1_2q2_2q314055_eval;
case 14060 : return &A2q1_2q2_2q314060_eval;
case 14115 : return &A2q1_2q2_2q314115_eval;
case 14130 : return &A2q1_2q2_2q314130_eval;
case 14150 : return &A2q1_2q2_2q314150_eval;
case 14160 : return &A2q1_2q2_2q314160_eval;
case 14350 : return &A2q1_2q2_2q314350_eval;
case 14355 : return &A2q1_2q2_2q314355_eval;
case 14380 : return &A2q1_2q2_2q314380_eval;
case 14390 : return &A2q1_2q2_2q314390_eval;
case 14415 : return &A2q1_2q2_2q314415_eval;
case 14420 : return &A2q1_2q2_2q314420_eval;
case 14710 : return &A2q1_2q2_2q314710_eval;
case 14715 : return &A2q1_2q2_2q314715_eval;
case 14800 : return &A2q1_2q2_2q314800_eval;
case 14820 : return &A2q1_2q2_2q314820_eval;
case 14835 : return &A2q1_2q2_2q314835_eval;
case 14850 : return &A2q1_2q2_2q314850_eval;
case 14920 : return &A2q1_2q2_2q314920_eval;
case 14930 : return &A2q1_2q2_2q314930_eval;
case 14980 : return &A2q1_2q2_2q314980_eval;
case 15000 : return &A2q1_2q2_2q315000_eval;
case 15050 : return &A2q1_2q2_2q315050_eval;
case 15060 : return &A2q1_2q2_2q315060_eval;
case 15135 : return &A2q1_2q2_2q315135_eval;
case 15140 : return &A2q1_2q2_2q315140_eval;
case 15195 : return &A2q1_2q2_2q315195_eval;
case 15210 : return &A2q1_2q2_2q315210_eval;
case 15230 : return &A2q1_2q2_2q315230_eval;
case 15240 : return &A2q1_2q2_2q315240_eval;
case 15905 : return &A2q1_2q2_2q315905_eval;
case 15910 : return &A2q1_2q2_2q315910_eval;
case 15935 : return &A2q1_2q2_2q315935_eval;
case 15945 : return &A2q1_2q2_2q315945_eval;
case 15970 : return &A2q1_2q2_2q315970_eval;
case 15975 : return &A2q1_2q2_2q315975_eval;
case 16265 : return &A2q1_2q2_2q316265_eval;
case 16270 : return &A2q1_2q2_2q316270_eval;
case 16355 : return &A2q1_2q2_2q316355_eval;
case 16375 : return &A2q1_2q2_2q316375_eval;
case 16390 : return &A2q1_2q2_2q316390_eval;
case 16405 : return &A2q1_2q2_2q316405_eval;
case 16475 : return &A2q1_2q2_2q316475_eval;
case 16485 : return &A2q1_2q2_2q316485_eval;
case 16535 : return &A2q1_2q2_2q316535_eval;
case 16555 : return &A2q1_2q2_2q316555_eval;
case 16605 : return &A2q1_2q2_2q316605_eval;
case 16615 : return &A2q1_2q2_2q316615_eval;
case 16690 : return &A2q1_2q2_2q316690_eval;
case 16695 : return &A2q1_2q2_2q316695_eval;
case 16750 : return &A2q1_2q2_2q316750_eval;
case 16765 : return &A2q1_2q2_2q316765_eval;
case 16785 : return &A2q1_2q2_2q316785_eval;
case 16795 : return &A2q1_2q2_2q316795_eval;
case 16985 : return &A2q1_2q2_2q316985_eval;
case 16990 : return &A2q1_2q2_2q316990_eval;
case 17015 : return &A2q1_2q2_2q317015_eval;
case 17025 : return &A2q1_2q2_2q317025_eval;
case 17050 : return &A2q1_2q2_2q317050_eval;
case 17055 : return &A2q1_2q2_2q317055_eval;
case 17525 : return &A2q1_2q2_2q317525_eval;
case 17530 : return &A2q1_2q2_2q317530_eval;
case 17645 : return &A2q1_2q2_2q317645_eval;
case 17670 : return &A2q1_2q2_2q317670_eval;
case 17680 : return &A2q1_2q2_2q317680_eval;
case 17700 : return &A2q1_2q2_2q317700_eval;
case 17735 : return &A2q1_2q2_2q317735_eval;
case 17745 : return &A2q1_2q2_2q317745_eval;
case 17825 : return &A2q1_2q2_2q317825_eval;
case 17850 : return &A2q1_2q2_2q317850_eval;
case 17895 : return &A2q1_2q2_2q317895_eval;
case 17910 : return &A2q1_2q2_2q317910_eval;
case 17950 : return &A2q1_2q2_2q317950_eval;
case 17955 : return &A2q1_2q2_2q317955_eval;
case 18040 : return &A2q1_2q2_2q318040_eval;
case 18060 : return &A2q1_2q2_2q318060_eval;
case 18075 : return &A2q1_2q2_2q318075_eval;
case 18090 : return &A2q1_2q2_2q318090_eval;
case 19505 : return &A2q1_2q2_2q319505_eval;
case 19510 : return &A2q1_2q2_2q319510_eval;
case 19595 : return &A2q1_2q2_2q319595_eval;
case 19615 : return &A2q1_2q2_2q319615_eval;
case 19630 : return &A2q1_2q2_2q319630_eval;
case 19645 : return &A2q1_2q2_2q319645_eval;
case 19685 : return &A2q1_2q2_2q319685_eval;
case 19690 : return &A2q1_2q2_2q319690_eval;
case 19805 : return &A2q1_2q2_2q319805_eval;
case 19830 : return &A2q1_2q2_2q319830_eval;
case 19840 : return &A2q1_2q2_2q319840_eval;
case 19860 : return &A2q1_2q2_2q319860_eval;
case 20315 : return &A2q1_2q2_2q320315_eval;
case 20335 : return &A2q1_2q2_2q320335_eval;
case 20345 : return &A2q1_2q2_2q320345_eval;
case 20370 : return &A2q1_2q2_2q320370_eval;
case 20485 : return &A2q1_2q2_2q320485_eval;
case 20490 : return &A2q1_2q2_2q320490_eval;
case 20530 : return &A2q1_2q2_2q320530_eval;
case 20545 : return &A2q1_2q2_2q320545_eval;
case 20560 : return &A2q1_2q2_2q320560_eval;
case 20580 : return &A2q1_2q2_2q320580_eval;
case 20665 : return &A2q1_2q2_2q320665_eval;
case 20670 : return &A2q1_2q2_2q320670_eval;
case 20795 : return &A2q1_2q2_2q320795_eval;
case 20805 : return &A2q1_2q2_2q320805_eval;
case 20855 : return &A2q1_2q2_2q320855_eval;
case 20875 : return &A2q1_2q2_2q320875_eval;
case 20925 : return &A2q1_2q2_2q320925_eval;
case 20935 : return &A2q1_2q2_2q320935_eval;
case 20975 : return &A2q1_2q2_2q320975_eval;
case 20985 : return &A2q1_2q2_2q320985_eval;
case 21065 : return &A2q1_2q2_2q321065_eval;
case 21090 : return &A2q1_2q2_2q321090_eval;
case 21135 : return &A2q1_2q2_2q321135_eval;
case 21150 : return &A2q1_2q2_2q321150_eval;
case 21395 : return &A2q1_2q2_2q321395_eval;
case 21415 : return &A2q1_2q2_2q321415_eval;
case 21425 : return &A2q1_2q2_2q321425_eval;
case 21450 : return &A2q1_2q2_2q321450_eval;
case 21565 : return &A2q1_2q2_2q321565_eval;
case 21570 : return &A2q1_2q2_2q321570_eval;
case 21825 : return &A2q1_2q2_2q321825_eval;
case 21835 : return &A2q1_2q2_2q321835_eval;
case 21855 : return &A2q1_2q2_2q321855_eval;
case 21870 : return &A2q1_2q2_2q321870_eval;
case 21925 : return &A2q1_2q2_2q321925_eval;
case 21930 : return &A2q1_2q2_2q321930_eval;
case 22090 : return &A2q1_2q2_2q322090_eval;
case 22095 : return &A2q1_2q2_2q322095_eval;
case 22150 : return &A2q1_2q2_2q322150_eval;
case 22165 : return &A2q1_2q2_2q322165_eval;
case 22185 : return &A2q1_2q2_2q322185_eval;
case 22195 : return &A2q1_2q2_2q322195_eval;
case 22270 : return &A2q1_2q2_2q322270_eval;
case 22275 : return &A2q1_2q2_2q322275_eval;
case 22360 : return &A2q1_2q2_2q322360_eval;
case 22380 : return &A2q1_2q2_2q322380_eval;
case 22395 : return &A2q1_2q2_2q322395_eval;
case 22410 : return &A2q1_2q2_2q322410_eval;
case 22690 : return &A2q1_2q2_2q322690_eval;
case 22705 : return &A2q1_2q2_2q322705_eval;
case 22720 : return &A2q1_2q2_2q322720_eval;
case 22740 : return &A2q1_2q2_2q322740_eval;
case 22825 : return &A2q1_2q2_2q322825_eval;
case 22830 : return &A2q1_2q2_2q322830_eval;
case 22905 : return &A2q1_2q2_2q322905_eval;
case 22915 : return &A2q1_2q2_2q322915_eval;
case 22935 : return &A2q1_2q2_2q322935_eval;
case 22950 : return &A2q1_2q2_2q322950_eval;
case 23005 : return &A2q1_2q2_2q323005_eval;
case 23010 : return &A2q1_2q2_2q323010_eval;
case 23645 : return &A2q1_2q2_2q323645_eval;
case 23650 : return &A2q1_2q2_2q323650_eval;
case 23705 : return &A2q1_2q2_2q323705_eval;
case 23720 : return &A2q1_2q2_2q323720_eval;
case 23740 : return &A2q1_2q2_2q323740_eval;
case 23750 : return &A2q1_2q2_2q323750_eval;
case 23825 : return &A2q1_2q2_2q323825_eval;
case 23830 : return &A2q1_2q2_2q323830_eval;
case 23915 : return &A2q1_2q2_2q323915_eval;
case 23935 : return &A2q1_2q2_2q323935_eval;
case 23950 : return &A2q1_2q2_2q323950_eval;
case 23965 : return &A2q1_2q2_2q323965_eval;
case 24245 : return &A2q1_2q2_2q324245_eval;
case 24260 : return &A2q1_2q2_2q324260_eval;
case 24275 : return &A2q1_2q2_2q324275_eval;
case 24295 : return &A2q1_2q2_2q324295_eval;
case 24380 : return &A2q1_2q2_2q324380_eval;
case 24385 : return &A2q1_2q2_2q324385_eval;
case 24460 : return &A2q1_2q2_2q324460_eval;
case 24470 : return &A2q1_2q2_2q324470_eval;
case 24490 : return &A2q1_2q2_2q324490_eval;
case 24505 : return &A2q1_2q2_2q324505_eval;
case 24560 : return &A2q1_2q2_2q324560_eval;
case 24565 : return &A2q1_2q2_2q324565_eval;
case 24725 : return &A2q1_2q2_2q324725_eval;
case 24730 : return &A2q1_2q2_2q324730_eval;
case 24785 : return &A2q1_2q2_2q324785_eval;
case 24800 : return &A2q1_2q2_2q324800_eval;
case 24820 : return &A2q1_2q2_2q324820_eval;
case 24830 : return &A2q1_2q2_2q324830_eval;
case 25085 : return &A2q1_2q2_2q325085_eval;
case 25090 : return &A2q1_2q2_2q325090_eval;
case 25205 : return &A2q1_2q2_2q325205_eval;
case 25230 : return &A2q1_2q2_2q325230_eval;
case 25240 : return &A2q1_2q2_2q325240_eval;
case 25260 : return &A2q1_2q2_2q325260_eval;
case 25505 : return &A2q1_2q2_2q325505_eval;
case 25520 : return &A2q1_2q2_2q325520_eval;
case 25565 : return &A2q1_2q2_2q325565_eval;
case 25590 : return &A2q1_2q2_2q325590_eval;
case 25670 : return &A2q1_2q2_2q325670_eval;
case 25680 : return &A2q1_2q2_2q325680_eval;
case 25720 : return &A2q1_2q2_2q325720_eval;
case 25730 : return &A2q1_2q2_2q325730_eval;
case 25780 : return &A2q1_2q2_2q325780_eval;
case 25800 : return &A2q1_2q2_2q325800_eval;
case 25850 : return &A2q1_2q2_2q325850_eval;
case 25860 : return &A2q1_2q2_2q325860_eval;
case 25985 : return &A2q1_2q2_2q325985_eval;
case 25990 : return &A2q1_2q2_2q325990_eval;
case 26075 : return &A2q1_2q2_2q326075_eval;
case 26095 : return &A2q1_2q2_2q326095_eval;
case 26110 : return &A2q1_2q2_2q326110_eval;
case 26125 : return &A2q1_2q2_2q326125_eval;
case 26165 : return &A2q1_2q2_2q326165_eval;
case 26170 : return &A2q1_2q2_2q326170_eval;
case 26285 : return &A2q1_2q2_2q326285_eval;
case 26310 : return &A2q1_2q2_2q326310_eval;
case 26320 : return &A2q1_2q2_2q326320_eval;
case 26340 : return &A2q1_2q2_2q326340_eval;
case 26795 : return &A2q1_2q2_2q326795_eval;
case 26815 : return &A2q1_2q2_2q326815_eval;
case 26825 : return &A2q1_2q2_2q326825_eval;
case 26850 : return &A2q1_2q2_2q326850_eval;
case 26965 : return &A2q1_2q2_2q326965_eval;
case 26970 : return &A2q1_2q2_2q326970_eval;
case 27010 : return &A2q1_2q2_2q327010_eval;
case 27025 : return &A2q1_2q2_2q327025_eval;
case 27040 : return &A2q1_2q2_2q327040_eval;
case 27060 : return &A2q1_2q2_2q327060_eval;
case 27145 : return &A2q1_2q2_2q327145_eval;
case 27150 : return &A2q1_2q2_2q327150_eval;
case 28565 : return &A2q1_2q2_2q328565_eval;
case 28580 : return &A2q1_2q2_2q328580_eval;
case 28595 : return &A2q1_2q2_2q328595_eval;
case 28615 : return &A2q1_2q2_2q328615_eval;
case 28700 : return &A2q1_2q2_2q328700_eval;
case 28705 : return &A2q1_2q2_2q328705_eval;
case 28745 : return &A2q1_2q2_2q328745_eval;
case 28760 : return &A2q1_2q2_2q328760_eval;
case 28805 : return &A2q1_2q2_2q328805_eval;
case 28830 : return &A2q1_2q2_2q328830_eval;
case 28910 : return &A2q1_2q2_2q328910_eval;
case 28920 : return &A2q1_2q2_2q328920_eval;
case 28955 : return &A2q1_2q2_2q328955_eval;
case 28975 : return &A2q1_2q2_2q328975_eval;
case 28985 : return &A2q1_2q2_2q328985_eval;
case 29010 : return &A2q1_2q2_2q329010_eval;
case 29125 : return &A2q1_2q2_2q329125_eval;
case 29130 : return &A2q1_2q2_2q329130_eval;
case 29600 : return &A2q1_2q2_2q329600_eval;
case 29605 : return &A2q1_2q2_2q329605_eval;
case 29630 : return &A2q1_2q2_2q329630_eval;
case 29640 : return &A2q1_2q2_2q329640_eval;
case 29665 : return &A2q1_2q2_2q329665_eval;
case 29670 : return &A2q1_2q2_2q329670_eval;
case 29860 : return &A2q1_2q2_2q329860_eval;
case 29870 : return &A2q1_2q2_2q329870_eval;
case 29890 : return &A2q1_2q2_2q329890_eval;
case 29905 : return &A2q1_2q2_2q329905_eval;
case 29960 : return &A2q1_2q2_2q329960_eval;
case 29965 : return &A2q1_2q2_2q329965_eval;
case 30040 : return &A2q1_2q2_2q330040_eval;
case 30050 : return &A2q1_2q2_2q330050_eval;
case 30100 : return &A2q1_2q2_2q330100_eval;
case 30120 : return &A2q1_2q2_2q330120_eval;
case 30170 : return &A2q1_2q2_2q330170_eval;
case 30180 : return &A2q1_2q2_2q330180_eval;
case 30250 : return &A2q1_2q2_2q330250_eval;
case 30265 : return &A2q1_2q2_2q330265_eval;
case 30280 : return &A2q1_2q2_2q330280_eval;
case 30300 : return &A2q1_2q2_2q330300_eval;
case 30385 : return &A2q1_2q2_2q330385_eval;
case 30390 : return &A2q1_2q2_2q330390_eval;
case 30680 : return &A2q1_2q2_2q330680_eval;
case 30685 : return &A2q1_2q2_2q330685_eval;
case 30710 : return &A2q1_2q2_2q330710_eval;
case 30720 : return &A2q1_2q2_2q330720_eval;
case 30745 : return &A2q1_2q2_2q330745_eval;
case 30750 : return &A2q1_2q2_2q330750_eval;
case 31415 : return &A2q1_2q2_2q331415_eval;
case 31425 : return &A2q1_2q2_2q331425_eval;
case 31445 : return &A2q1_2q2_2q331445_eval;
case 31460 : return &A2q1_2q2_2q331460_eval;
case 31515 : return &A2q1_2q2_2q331515_eval;
case 31520 : return &A2q1_2q2_2q331520_eval;
case 31595 : return &A2q1_2q2_2q331595_eval;
case 31605 : return &A2q1_2q2_2q331605_eval;
case 31655 : return &A2q1_2q2_2q331655_eval;
case 31675 : return &A2q1_2q2_2q331675_eval;
case 31725 : return &A2q1_2q2_2q331725_eval;
case 31735 : return &A2q1_2q2_2q331735_eval;
case 31805 : return &A2q1_2q2_2q331805_eval;
case 31820 : return &A2q1_2q2_2q331820_eval;
case 31835 : return &A2q1_2q2_2q331835_eval;
case 31855 : return &A2q1_2q2_2q331855_eval;
case 31940 : return &A2q1_2q2_2q331940_eval;
case 31945 : return &A2q1_2q2_2q331945_eval;
case 32235 : return &A2q1_2q2_2q332235_eval;
case 32240 : return &A2q1_2q2_2q332240_eval;
case 32265 : return &A2q1_2q2_2q332265_eval;
case 32275 : return &A2q1_2q2_2q332275_eval;
case 32300 : return &A2q1_2q2_2q332300_eval;
case 32305 : return &A2q1_2q2_2q332305_eval;
case 32495 : return &A2q1_2q2_2q332495_eval;
case 32505 : return &A2q1_2q2_2q332505_eval;
case 32525 : return &A2q1_2q2_2q332525_eval;
case 32540 : return &A2q1_2q2_2q332540_eval;
case 32595 : return &A2q1_2q2_2q332595_eval;
case 32600 : return &A2q1_2q2_2q332600_eval;
case 32855 : return &A2q1_2q2_2q332855_eval;
case 32865 : return &A2q1_2q2_2q332865_eval;
case 32945 : return &A2q1_2q2_2q332945_eval;
case 32970 : return &A2q1_2q2_2q332970_eval;
case 33015 : return &A2q1_2q2_2q333015_eval;
case 33030 : return &A2q1_2q2_2q333030_eval;
case 33065 : return &A2q1_2q2_2q333065_eval;
case 33080 : return &A2q1_2q2_2q333080_eval;
case 33125 : return &A2q1_2q2_2q333125_eval;
case 33150 : return &A2q1_2q2_2q333150_eval;
case 33230 : return &A2q1_2q2_2q333230_eval;
case 33240 : return &A2q1_2q2_2q333240_eval;
case 33495 : return &A2q1_2q2_2q333495_eval;
case 33500 : return &A2q1_2q2_2q333500_eval;
case 33555 : return &A2q1_2q2_2q333555_eval;
case 33570 : return &A2q1_2q2_2q333570_eval;
case 33590 : return &A2q1_2q2_2q333590_eval;
case 33600 : return &A2q1_2q2_2q333600_eval;
case 33755 : return &A2q1_2q2_2q333755_eval;
case 33765 : return &A2q1_2q2_2q333765_eval;
case 33815 : return &A2q1_2q2_2q333815_eval;
case 33835 : return &A2q1_2q2_2q333835_eval;
case 33885 : return &A2q1_2q2_2q333885_eval;
case 33895 : return &A2q1_2q2_2q333895_eval;
case 33935 : return &A2q1_2q2_2q333935_eval;
case 33945 : return &A2q1_2q2_2q333945_eval;
case 34025 : return &A2q1_2q2_2q334025_eval;
case 34050 : return &A2q1_2q2_2q334050_eval;
case 34095 : return &A2q1_2q2_2q334095_eval;
case 34110 : return &A2q1_2q2_2q334110_eval;
case 34355 : return &A2q1_2q2_2q334355_eval;
case 34375 : return &A2q1_2q2_2q334375_eval;
case 34385 : return &A2q1_2q2_2q334385_eval;
case 34410 : return &A2q1_2q2_2q334410_eval;
case 34525 : return &A2q1_2q2_2q334525_eval;
case 34530 : return &A2q1_2q2_2q334530_eval;
case 34785 : return &A2q1_2q2_2q334785_eval;
case 34795 : return &A2q1_2q2_2q334795_eval;
case 34815 : return &A2q1_2q2_2q334815_eval;
case 34830 : return &A2q1_2q2_2q334830_eval;
case 34885 : return &A2q1_2q2_2q334885_eval;
case 34890 : return &A2q1_2q2_2q334890_eval;
case 35045 : return &A2q1_2q2_2q335045_eval;
case 35060 : return &A2q1_2q2_2q335060_eval;
case 35075 : return &A2q1_2q2_2q335075_eval;
case 35095 : return &A2q1_2q2_2q335095_eval;
case 35180 : return &A2q1_2q2_2q335180_eval;
case 35185 : return &A2q1_2q2_2q335185_eval;
case 35225 : return &A2q1_2q2_2q335225_eval;
case 35240 : return &A2q1_2q2_2q335240_eval;
case 35285 : return &A2q1_2q2_2q335285_eval;
case 35310 : return &A2q1_2q2_2q335310_eval;
case 35390 : return &A2q1_2q2_2q335390_eval;
case 35400 : return &A2q1_2q2_2q335400_eval;
case 35435 : return &A2q1_2q2_2q335435_eval;
case 35455 : return &A2q1_2q2_2q335455_eval;
case 35465 : return &A2q1_2q2_2q335465_eval;
case 35490 : return &A2q1_2q2_2q335490_eval;
case 35605 : return &A2q1_2q2_2q335605_eval;
case 35610 : return &A2q1_2q2_2q335610_eval;
case 36080 : return &A2q1_2q2_2q336080_eval;
case 36085 : return &A2q1_2q2_2q336085_eval;
case 36110 : return &A2q1_2q2_2q336110_eval;
case 36120 : return &A2q1_2q2_2q336120_eval;
case 36145 : return &A2q1_2q2_2q336145_eval;
case 36150 : return &A2q1_2q2_2q336150_eval;
case 37635 : return &A2q1_2q2_2q337635_eval;
case 37640 : return &A2q1_2q2_2q337640_eval;
case 37665 : return &A2q1_2q2_2q337665_eval;
case 37675 : return &A2q1_2q2_2q337675_eval;
case 37700 : return &A2q1_2q2_2q337700_eval;
case 37705 : return &A2q1_2q2_2q337705_eval;
case 37815 : return &A2q1_2q2_2q337815_eval;
case 37820 : return &A2q1_2q2_2q337820_eval;
case 37875 : return &A2q1_2q2_2q337875_eval;
case 37890 : return &A2q1_2q2_2q337890_eval;
case 37910 : return &A2q1_2q2_2q337910_eval;
case 37920 : return &A2q1_2q2_2q337920_eval;
case 38025 : return &A2q1_2q2_2q338025_eval;
case 38035 : return &A2q1_2q2_2q338035_eval;
case 38055 : return &A2q1_2q2_2q338055_eval;
case 38070 : return &A2q1_2q2_2q338070_eval;
case 38125 : return &A2q1_2q2_2q338125_eval;
case 38130 : return &A2q1_2q2_2q338130_eval;
case 38240 : return &A2q1_2q2_2q338240_eval;
case 38245 : return &A2q1_2q2_2q338245_eval;
case 38270 : return &A2q1_2q2_2q338270_eval;
case 38280 : return &A2q1_2q2_2q338280_eval;
case 38305 : return &A2q1_2q2_2q338305_eval;
case 38310 : return &A2q1_2q2_2q338310_eval;
case 39190 : return &A2q1_2q2_2q339190_eval;
case 39195 : return &A2q1_2q2_2q339195_eval;
case 39220 : return &A2q1_2q2_2q339220_eval;
case 39230 : return &A2q1_2q2_2q339230_eval;
case 39255 : return &A2q1_2q2_2q339255_eval;
case 39260 : return &A2q1_2q2_2q339260_eval;
case 39370 : return &A2q1_2q2_2q339370_eval;
case 39375 : return &A2q1_2q2_2q339375_eval;
case 39430 : return &A2q1_2q2_2q339430_eval;
case 39445 : return &A2q1_2q2_2q339445_eval;
case 39465 : return &A2q1_2q2_2q339465_eval;
case 39475 : return &A2q1_2q2_2q339475_eval;
case 39580 : return &A2q1_2q2_2q339580_eval;
case 39590 : return &A2q1_2q2_2q339590_eval;
case 39610 : return &A2q1_2q2_2q339610_eval;
case 39625 : return &A2q1_2q2_2q339625_eval;
case 39680 : return &A2q1_2q2_2q339680_eval;
case 39685 : return &A2q1_2q2_2q339685_eval;
case 39795 : return &A2q1_2q2_2q339795_eval;
case 39800 : return &A2q1_2q2_2q339800_eval;
case 39825 : return &A2q1_2q2_2q339825_eval;
case 39835 : return &A2q1_2q2_2q339835_eval;
case 39860 : return &A2q1_2q2_2q339860_eval;
case 39865 : return &A2q1_2q2_2q339865_eval;
case 40270 : return &A2q1_2q2_2q340270_eval;
case 40275 : return &A2q1_2q2_2q340275_eval;
case 40300 : return &A2q1_2q2_2q340300_eval;
case 40310 : return &A2q1_2q2_2q340310_eval;
case 40335 : return &A2q1_2q2_2q340335_eval;
case 40340 : return &A2q1_2q2_2q340340_eval;
case 40630 : return &A2q1_2q2_2q340630_eval;
case 40635 : return &A2q1_2q2_2q340635_eval;
case 40720 : return &A2q1_2q2_2q340720_eval;
case 40740 : return &A2q1_2q2_2q340740_eval;
case 40755 : return &A2q1_2q2_2q340755_eval;
case 40770 : return &A2q1_2q2_2q340770_eval;
case 40840 : return &A2q1_2q2_2q340840_eval;
case 40850 : return &A2q1_2q2_2q340850_eval;
case 40900 : return &A2q1_2q2_2q340900_eval;
case 40920 : return &A2q1_2q2_2q340920_eval;
case 40970 : return &A2q1_2q2_2q340970_eval;
case 40980 : return &A2q1_2q2_2q340980_eval;
case 41055 : return &A2q1_2q2_2q341055_eval;
case 41060 : return &A2q1_2q2_2q341060_eval;
case 41115 : return &A2q1_2q2_2q341115_eval;
case 41130 : return &A2q1_2q2_2q341130_eval;
case 41150 : return &A2q1_2q2_2q341150_eval;
case 41160 : return &A2q1_2q2_2q341160_eval;
case 41530 : return &A2q1_2q2_2q341530_eval;
case 41535 : return &A2q1_2q2_2q341535_eval;
case 41590 : return &A2q1_2q2_2q341590_eval;
case 41605 : return &A2q1_2q2_2q341605_eval;
case 41625 : return &A2q1_2q2_2q341625_eval;
case 41635 : return &A2q1_2q2_2q341635_eval;
case 41710 : return &A2q1_2q2_2q341710_eval;
case 41715 : return &A2q1_2q2_2q341715_eval;
case 41800 : return &A2q1_2q2_2q341800_eval;
case 41820 : return &A2q1_2q2_2q341820_eval;
case 41835 : return &A2q1_2q2_2q341835_eval;
case 41850 : return &A2q1_2q2_2q341850_eval;
case 42130 : return &A2q1_2q2_2q342130_eval;
case 42145 : return &A2q1_2q2_2q342145_eval;
case 42160 : return &A2q1_2q2_2q342160_eval;
case 42180 : return &A2q1_2q2_2q342180_eval;
case 42265 : return &A2q1_2q2_2q342265_eval;
case 42270 : return &A2q1_2q2_2q342270_eval;
case 42345 : return &A2q1_2q2_2q342345_eval;
case 42355 : return &A2q1_2q2_2q342355_eval;
case 42375 : return &A2q1_2q2_2q342375_eval;
case 42390 : return &A2q1_2q2_2q342390_eval;
case 42445 : return &A2q1_2q2_2q342445_eval;
case 42450 : return &A2q1_2q2_2q342450_eval;
case 42820 : return &A2q1_2q2_2q342820_eval;
case 42830 : return &A2q1_2q2_2q342830_eval;
case 42850 : return &A2q1_2q2_2q342850_eval;
case 42865 : return &A2q1_2q2_2q342865_eval;
case 42920 : return &A2q1_2q2_2q342920_eval;
case 42925 : return &A2q1_2q2_2q342925_eval;
case 43000 : return &A2q1_2q2_2q343000_eval;
case 43010 : return &A2q1_2q2_2q343010_eval;
case 43060 : return &A2q1_2q2_2q343060_eval;
case 43080 : return &A2q1_2q2_2q343080_eval;
case 43130 : return &A2q1_2q2_2q343130_eval;
case 43140 : return &A2q1_2q2_2q343140_eval;
case 43210 : return &A2q1_2q2_2q343210_eval;
case 43225 : return &A2q1_2q2_2q343225_eval;
case 43240 : return &A2q1_2q2_2q343240_eval;
case 43260 : return &A2q1_2q2_2q343260_eval;
case 43345 : return &A2q1_2q2_2q343345_eval;
case 43350 : return &A2q1_2q2_2q343350_eval;
case 43640 : return &A2q1_2q2_2q343640_eval;
case 43645 : return &A2q1_2q2_2q343645_eval;
case 43670 : return &A2q1_2q2_2q343670_eval;
case 43680 : return &A2q1_2q2_2q343680_eval;
case 43705 : return &A2q1_2q2_2q343705_eval;
case 43710 : return &A2q1_2q2_2q343710_eval;
case 44115 : return &A2q1_2q2_2q344115_eval;
case 44120 : return &A2q1_2q2_2q344120_eval;
case 44145 : return &A2q1_2q2_2q344145_eval;
case 44155 : return &A2q1_2q2_2q344155_eval;
case 44180 : return &A2q1_2q2_2q344180_eval;
case 44185 : return &A2q1_2q2_2q344185_eval;
case 44295 : return &A2q1_2q2_2q344295_eval;
case 44300 : return &A2q1_2q2_2q344300_eval;
case 44355 : return &A2q1_2q2_2q344355_eval;
case 44370 : return &A2q1_2q2_2q344370_eval;
case 44390 : return &A2q1_2q2_2q344390_eval;
case 44400 : return &A2q1_2q2_2q344400_eval;
case 44505 : return &A2q1_2q2_2q344505_eval;
case 44515 : return &A2q1_2q2_2q344515_eval;
case 44535 : return &A2q1_2q2_2q344535_eval;
case 44550 : return &A2q1_2q2_2q344550_eval;
case 44605 : return &A2q1_2q2_2q344605_eval;
case 44610 : return &A2q1_2q2_2q344610_eval;
case 44720 : return &A2q1_2q2_2q344720_eval;
case 44725 : return &A2q1_2q2_2q344725_eval;
case 44750 : return &A2q1_2q2_2q344750_eval;
case 44760 : return &A2q1_2q2_2q344760_eval;
case 44785 : return &A2q1_2q2_2q344785_eval;
case 44790 : return &A2q1_2q2_2q344790_eval;

default: _WARNING3("Unknown pointer amplitude (*A2q1_2q2_2q3_Tree_Ptr(int hc)) - case:",hc , " - throw BH error.");
		throw BHerror("case missing for tree amplitude!");

}
}


template complex<R>  (*A2q1_2q2_2q3_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q1_2q2_2q3_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q1_2q2_2q3_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q1_2q2_2q3_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
44790: Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 3]
38310: Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 3]
43710: Q[m, 1] Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 3]
30750: Q[m, 1] Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 2]
36150: Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 3]
29670: Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 2]
44610: Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 3]
38130: Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 3]
42450: Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 3]
23010: Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 2]
34890: Q[m, 1] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 3]
21930: Q[m, 1] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 2]
43350: Q[m, 1] Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 3]
30390: Q[m, 1] Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 2]
42270: Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 3]
22830: Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 2]
27150: Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 2]
20670: Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 2]
35610: Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 3]
29130: Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 2]
34530: Q[m, 1] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 3]
21570: Q[m, 1] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 2]
26970: Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 2]
20490: Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 2]
44760: Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 3] Q[p, 3]
38280: Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 3]
43680: Q[m, 1] Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 2] Q[p, 3]
30720: Q[m, 1] Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 3] Q[p, 2]
36120: Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 3]
29640: Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 3] Q[p, 2]
44400: Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 3] Q[p, 3]
37920: Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 3]
41160: Q[m, 1] Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 1] Q[p, 3]
15240: Q[m, 1] Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 3] Q[p, 1]
33600: Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 3]
14160: Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 3] Q[p, 1]
43140: Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 2] Q[p, 3]
30180: Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 3] Q[p, 2]
40980: Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 1] Q[p, 3]
15060: Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 3] Q[p, 1]
25860: Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 1] Q[p, 2]
12900: Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 2] Q[p, 1]
35400: Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 3]
28920: Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 3] Q[p, 2]
33240: Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 3]
13800: Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 3] Q[p, 1]
25680: Q[m, 1] Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 1] Q[p, 2]
12720: Q[m, 1] Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 2] Q[p, 1]
44550: Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 3]
38070: Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 3]
42390: Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 3]
22950: Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 2]
34830: Q[m, 1] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 3]
21870: Q[m, 1] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 2]
44370: Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 3]
37890: Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 3]
41130: Q[m, 1] Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 3]
15210: Q[m, 1] Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 1]
33570: Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 3]
14130: Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 1]
41850: Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 3]
22410: Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 2]
40770: Q[m, 1] Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 3]
14850: Q[m, 1] Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 1]
18090: Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 2]
11610: Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 1]
34110: Q[m, 1] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 3]
21150: Q[m, 1] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 2]
33030: Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 3]
13590: Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 1]
17910: Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 2]
11430: Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 1]
43260: Q[m, 1] Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 2] Q[p, 3]
30300: Q[m, 1] Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 3] Q[p, 2]
42180: Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 2] Q[p, 3]
22740: Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 2]
27060: Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 2] Q[p, 2]
20580: Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 2]
43080: Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 2] Q[p, 3]
30120: Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 3] Q[p, 2]
40920: Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 1] Q[p, 3]
15000: Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 3] Q[p, 1]
25800: Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 1] Q[p, 2]
12840: Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 2] Q[p, 1]
41820: Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 2] Q[p, 3]
22380: Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 2]
40740: Q[m, 1] Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 1] Q[p, 3]
14820: Q[m, 1] Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 3] Q[p, 1]
18060: Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 2]
11580: Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 2] Q[p, 1]
26340: Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 2] Q[p, 2]
19860: Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 2]
25260: Q[m, 1] Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 1] Q[p, 2]
12300: Q[m, 1] Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 2] Q[p, 1]
17700: Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 2]
11220: Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 2] Q[p, 1]
35490: Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 3]
29010: Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 2]
34410: Q[m, 1] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 3]
21450: Q[m, 1] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 2]
26850: Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 2]
20370: Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 2]
35310: Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 3]
28830: Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 2]
33150: Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 3]
13710: Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 1]
25590: Q[m, 1] Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 2]
12630: Q[m, 1] Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 1]
34050: Q[m, 1] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 3]
21090: Q[m, 1] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 2]
32970: Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 3]
13530: Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 1]
17850: Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 2]
11370: Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 1]
26310: Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 2]
19830: Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 2]
25230: Q[m, 1] Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 2]
12270: Q[m, 1] Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 1]
17670: Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 2]
11190: Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 1]
44785: Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 3]
38305: Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 3]
43705: Q[p, 1] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 3]
30745: Q[p, 1] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 2]
36145: Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 3]
29665: Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 2]
44605: Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 3]
38125: Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 3]
42445: Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 3]
23005: Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 2]
34885: Q[p, 1] Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 3]
21925: Q[p, 1] Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 2]
43345: Q[p, 1] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 3]
30385: Q[p, 1] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 2]
42265: Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 3]
22825: Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 2]
27145: Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 2]
20665: Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 2]
35605: Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 3]
29125: Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 2]
34525: Q[p, 1] Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 3]
21565: Q[p, 1] Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 2]
26965: Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 2]
20485: Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 2]
44725: Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 3]
38245: Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 3]
43645: Q[p, 1] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 3]
30685: Q[p, 1] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 2]
36085: Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 3]
29605: Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 2]
44185: Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 3]
37705: Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 3]
39865: Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 3]
7465: Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 1]
32305: Q[p, 1] Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 3]
6385: Q[p, 1] Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 1]
42925: Q[p, 1] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 3]
29965: Q[p, 1] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 2]
39685: Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 3]
7285: Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 1]
24565: Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 2]
5125: Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 1]
35185: Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 3]
28705: Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 2]
31945: Q[p, 1] Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 3]
6025: Q[p, 1] Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 1]
24385: Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 2]
4945: Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 1]
44515: Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 3]
38035: Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 3] Q[m, 3]
42355: Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 3]
22915: Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 3] Q[m, 2]
34795: Q[p, 1] Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 2] Q[m, 3]
21835: Q[p, 1] Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 3] Q[m, 2]
44155: Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 3]
37675: Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 3] Q[m, 3]
39835: Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 3]
7435: Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 3] Q[m, 1]
32275: Q[p, 1] Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 1] Q[m, 3]
6355: Q[p, 1] Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 3] Q[m, 1]
41635: Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 3]
22195: Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 3] Q[m, 2]
39475: Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 3]
7075: Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 3] Q[m, 1]
16795: Q[p, 1] Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 1] Q[m, 2]
3835: Q[p, 1] Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 2] Q[m, 1]
33895: Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 2] Q[m, 3]
20935: Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 3] Q[m, 2]
31735: Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 1] Q[m, 3]
5815: Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 3] Q[m, 1]
16615: Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 1] Q[m, 2]
3655: Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 2] Q[m, 1]
43225: Q[p, 1] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 3]
30265: Q[p, 1] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 2]
42145: Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 3]
22705: Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 2]
27025: Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 2]
20545: Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 2]
42865: Q[p, 1] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 3]
29905: Q[p, 1] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 2]
39625: Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 3]
7225: Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 1]
24505: Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 2]
5065: Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 1]
41605: Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 3]
22165: Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 2]
39445: Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 3]
7045: Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 1]
16765: Q[p, 1] Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 2]
3805: Q[p, 1] Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 1]
26125: Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 2]
19645: Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 2]
23965: Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 2]
4525: Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 1]
16405: Q[p, 1] Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 2]
3445: Q[p, 1] Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 1]
35455: Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 2] Q[m, 3]
28975: Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 2]
34375: Q[p, 1] Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 2] Q[m, 3]
21415: Q[p, 1] Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 3] Q[m, 2]
26815: Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 2]
20335: Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 2] Q[m, 2]
35095: Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 2] Q[m, 3]
28615: Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 2]
31855: Q[p, 1] Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 1] Q[m, 3]
5935: Q[p, 1] Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 3] Q[m, 1]
24295: Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 2]
4855: Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 2] Q[m, 1]
33835: Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 2] Q[m, 3]
20875: Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 3] Q[m, 2]
31675: Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 1] Q[m, 3]
5755: Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 3] Q[m, 1]
16555: Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 1] Q[m, 2]
3595: Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 2] Q[m, 1]
26095: Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 2]
19615: Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 2] Q[m, 2]
23935: Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 2]
4495: Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 2] Q[m, 1]
16375: Q[p, 1] Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 1] Q[m, 2]
3415: Q[p, 1] Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 2] Q[m, 1]
44750: Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 3] Q[p, 3]
38270: Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 3]
43670: Q[m, 2] Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 2] Q[p, 3]
30710: Q[m, 2] Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 3] Q[p, 2]
36110: Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 3]
29630: Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 3] Q[p, 2]
44390: Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 3] Q[p, 3]
37910: Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 3]
41150: Q[m, 2] Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 1] Q[p, 3]
15230: Q[m, 2] Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 3] Q[p, 1]
33590: Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 3]
14150: Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 3] Q[p, 1]
43130: Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 2] Q[p, 3]
30170: Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 3] Q[p, 2]
40970: Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 1] Q[p, 3]
15050: Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 3] Q[p, 1]
25850: Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 1] Q[p, 2]
12890: Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 2] Q[p, 1]
35390: Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 3]
28910: Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 3] Q[p, 2]
33230: Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 3]
13790: Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 3] Q[p, 1]
25670: Q[m, 2] Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 1] Q[p, 2]
12710: Q[m, 2] Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 2] Q[p, 1]
44720: Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 3]
38240: Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 3]
43640: Q[m, 2] Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 3]
30680: Q[m, 2] Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 2]
36080: Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 3]
29600: Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 2]
44180: Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 3]
37700: Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 3]
39860: Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 3]
7460: Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 1]
32300: Q[m, 2] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 3]
6380: Q[m, 2] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 1]
42920: Q[m, 2] Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 3]
29960: Q[m, 2] Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 2]
39680: Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 3]
7280: Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 1]
24560: Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 2]
5120: Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 1]
35180: Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 3]
28700: Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 2]
31940: Q[m, 2] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 3]
6020: Q[m, 2] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 1]
24380: Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 2]
4940: Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 1]
44300: Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 3]
37820: Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 3]
41060: Q[m, 2] Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 3]
15140: Q[m, 2] Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 1]
33500: Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 3]
14060: Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 1]
44120: Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 3]
37640: Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 3]
39800: Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 3]
7400: Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 1]
32240: Q[m, 2] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 3]
6320: Q[m, 2] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 1]
40340: Q[m, 2] Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 3]
14420: Q[m, 2] Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 1]
39260: Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 3]
6860: Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 1]
9020: Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 1]
2540: Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 1]
32600: Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 3]
13160: Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 1]
31520: Q[m, 2] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 3]
5600: Q[m, 2] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 1]
8840: Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 1]
2360: Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 1]
43010: Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 2] Q[p, 3]
30050: Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 3] Q[p, 2]
40850: Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 1] Q[p, 3]
14930: Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 3] Q[p, 1]
25730: Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 1] Q[p, 2]
12770: Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 2] Q[p, 1]
42830: Q[m, 2] Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 2] Q[p, 3]
29870: Q[m, 2] Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 3] Q[p, 2]
39590: Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 1] Q[p, 3]
7190: Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 1]
24470: Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 1] Q[p, 2]
5030: Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 1]
40310: Q[m, 2] Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 1] Q[p, 3]
14390: Q[m, 2] Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 3] Q[p, 1]
39230: Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 1] Q[p, 3]
6830: Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 1]
8990: Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 1] Q[p, 1]
2510: Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 1]
24830: Q[m, 2] Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 1] Q[p, 2]
11870: Q[m, 2] Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 2] Q[p, 1]
23750: Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 1] Q[p, 2]
4310: Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 1]
8630: Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 1] Q[p, 1]
2150: Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 1]
35240: Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 3]
28760: Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 2]
33080: Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 3]
13640: Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 1]
25520: Q[m, 2] Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 2]
12560: Q[m, 2] Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 1]
35060: Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 3]
28580: Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 2]
31820: Q[m, 2] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 3]
5900: Q[m, 2] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 1]
24260: Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 2]
4820: Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 1]
32540: Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 3]
13100: Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 1]
31460: Q[m, 2] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 3]
5540: Q[m, 2] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 1]
8780: Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 1]
2300: Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 1]
24800: Q[m, 2] Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 2]
11840: Q[m, 2] Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 1]
23720: Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 2]
4280: Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 1]
8600: Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 1]
2120: Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 1]
44535: Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 3]
38055: Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 3]
42375: Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 3]
22935: Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 2]
34815: Q[p, 2] Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 3]
21855: Q[p, 2] Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 2]
44355: Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 3]
37875: Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 3]
41115: Q[p, 2] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 3]
15195: Q[p, 2] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 1]
33555: Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 3]
14115: Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 1]
41835: Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 3]
22395: Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 2]
40755: Q[p, 2] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 3]
14835: Q[p, 2] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 1]
18075: Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 2]
11595: Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 1]
34095: Q[p, 2] Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 3]
21135: Q[p, 2] Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 2]
33015: Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 3]
13575: Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 1]
17895: Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 2]
11415: Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 1]
44505: Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 3]
38025: Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 3] Q[m, 3]
42345: Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 3]
22905: Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 3] Q[m, 2]
34785: Q[p, 2] Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 2] Q[m, 3]
21825: Q[p, 2] Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 3] Q[m, 2]
44145: Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 3]
37665: Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 3] Q[m, 3]
39825: Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 3]
7425: Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 3] Q[m, 1]
32265: Q[p, 2] Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 1] Q[m, 3]
6345: Q[p, 2] Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 3] Q[m, 1]
41625: Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 3]
22185: Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 3] Q[m, 2]
39465: Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 3]
7065: Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 3] Q[m, 1]
16785: Q[p, 2] Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 1] Q[m, 2]
3825: Q[p, 2] Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 2] Q[m, 1]
33885: Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 2] Q[m, 3]
20925: Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 3] Q[m, 2]
31725: Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 1] Q[m, 3]
5805: Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 3] Q[m, 1]
16605: Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 1] Q[m, 2]
3645: Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 2] Q[m, 1]
44295: Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 3]
37815: Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 3]
41055: Q[p, 2] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 3]
15135: Q[p, 2] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 3] Q[p, 1]
33495: Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 3]
14055: Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 3] Q[m, 3] Q[p, 1]
44115: Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 3]
37635: Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 3]
39795: Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 3]
7395: Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 3] Q[m, 1]
32235: Q[p, 2] Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 3]
6315: Q[p, 2] Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 3] Q[m, 1]
40335: Q[p, 2] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 3]
14415: Q[p, 2] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 1]
39255: Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 3]
6855: Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 1]
9015: Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 1]
2535: Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 1]
32595: Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 3]
13155: Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 1]
31515: Q[p, 2] Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 3]
5595: Q[p, 2] Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 1]
8835: Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 1]
2355: Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 1]
41715: Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 3]
22275: Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 2]
40635: Q[p, 2] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 3]
14715: Q[p, 2] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 1]
17955: Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 2]
11475: Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 1]
41535: Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 3]
22095: Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 2]
39375: Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 3]
6975: Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 1]
16695: Q[p, 2] Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 2]
3735: Q[p, 2] Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 1]
40275: Q[p, 2] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 3]
14355: Q[p, 2] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 1]
39195: Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 3]
6795: Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 1]
8955: Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 1]
2475: Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 1]
17055: Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 2]
10575: Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 1]
15975: Q[p, 2] Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 2]
3015: Q[p, 2] Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 1]
8415: Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 1]
1935: Q[p, 2] Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 1]
33945: Q[p, 2] Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 2] Q[m, 3]
20985: Q[p, 2] Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 3] Q[m, 2]
32865: Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 1] Q[m, 3]
13425: Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 1]
17745: Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 1] Q[m, 2]
11265: Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 1]
33765: Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 2] Q[m, 3]
20805: Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 3] Q[m, 2]
31605: Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 1] Q[m, 3]
5685: Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 3] Q[m, 1]
16485: Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 1] Q[m, 2]
3525: Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 2] Q[m, 1]
32505: Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 1] Q[m, 3]
13065: Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 1]
31425: Q[p, 2] Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 1] Q[m, 3]
5505: Q[p, 2] Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 3] Q[m, 1]
8745: Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 1]
2265: Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 1] Q[m, 1]
17025: Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 1] Q[m, 2]
10545: Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 1]
15945: Q[p, 2] Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 1] Q[m, 2]
2985: Q[p, 2] Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 2] Q[m, 1]
8385: Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 1]
1905: Q[p, 2] Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 1] Q[m, 1]
43240: Q[m, 3] Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 2] Q[p, 3]
30280: Q[m, 3] Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 3] Q[p, 2]
42160: Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 2] Q[p, 3]
22720: Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 2]
27040: Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 2] Q[p, 2]
20560: Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 2]
43060: Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 2] Q[p, 3]
30100: Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 3] Q[p, 2]
40900: Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 1] Q[p, 3]
14980: Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 3] Q[p, 1]
25780: Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 1] Q[p, 2]
12820: Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 2] Q[p, 1]
41800: Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 2] Q[p, 3]
22360: Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 2]
40720: Q[m, 3] Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 1] Q[p, 3]
14800: Q[m, 3] Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 3] Q[p, 1]
18040: Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 2]
11560: Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 2] Q[p, 1]
26320: Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 2] Q[p, 2]
19840: Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 2]
25240: Q[m, 3] Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 1] Q[p, 2]
12280: Q[m, 3] Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 2] Q[p, 1]
17680: Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 2]
11200: Q[m, 3] Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 2] Q[p, 1]
43210: Q[m, 3] Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 3]
30250: Q[m, 3] Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 2]
42130: Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 3]
22690: Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 2] Q[p, 3] Q[m, 2]
27010: Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 2]
20530: Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 3] Q[p, 2] Q[m, 2]
42850: Q[m, 3] Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 3]
29890: Q[m, 3] Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 2]
39610: Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 3]
7210: Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 1]
24490: Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 2]
5050: Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 1]
41590: Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 3]
22150: Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 2]
39430: Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 3]
7030: Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 1]
16750: Q[m, 3] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 2]
3790: Q[m, 3] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 1]
26110: Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 2]
19630: Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 2]
23950: Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 2]
4510: Q[m, 3] Q[p, 1] Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 1]
16390: Q[m, 3] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 2]
3430: Q[m, 3] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 1]
43000: Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 2] Q[p, 3]
30040: Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 3] Q[p, 2]
40840: Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 1] Q[p, 3]
14920: Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 3] Q[p, 1]
25720: Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 1] Q[p, 2]
12760: Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 2] Q[p, 1]
42820: Q[m, 3] Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 2] Q[p, 3]
29860: Q[m, 3] Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 3] Q[p, 2]
39580: Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 1] Q[p, 3]
7180: Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 2] Q[p, 3] Q[m, 1]
24460: Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 1] Q[p, 2]
5020: Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 3] Q[p, 2] Q[m, 1]
40300: Q[m, 3] Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 1] Q[p, 3]
14380: Q[m, 3] Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 3] Q[p, 1]
39220: Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 1] Q[p, 3]
6820: Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 1]
8980: Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 3] Q[m, 1] Q[p, 1]
2500: Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 1]
24820: Q[m, 3] Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 1] Q[p, 2]
11860: Q[m, 3] Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 2] Q[p, 1]
23740: Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 1] Q[p, 2]
4300: Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 1]
8620: Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 2] Q[m, 1] Q[p, 1]
2140: Q[m, 3] Q[m, 2] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 1]
41710: Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 3]
22270: Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 1] Q[p, 3] Q[m, 2]
40630: Q[m, 3] Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 3]
14710: Q[m, 3] Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 3] Q[p, 1]
17950: Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 3] Q[p, 1] Q[m, 2]
11470: Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 3] Q[m, 2] Q[p, 1]
41530: Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 3]
22090: Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 1] Q[p, 3] Q[m, 2]
39370: Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 3]
6970: Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 2] Q[p, 3] Q[m, 1]
16690: Q[m, 3] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 1] Q[m, 2]
3730: Q[m, 3] Q[p, 2] Q[p, 1] Q[p, 3] Q[m, 2] Q[m, 1]
40270: Q[m, 3] Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 3]
14350: Q[m, 3] Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 3] Q[p, 1]
39190: Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 3]
6790: Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 1] Q[p, 3] Q[m, 1]
8950: Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 3] Q[m, 1] Q[p, 1]
2470: Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 3] Q[p, 1] Q[m, 1]
17050: Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 2]
10570: Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 1]
15970: Q[m, 3] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 2]
3010: Q[m, 3] Q[p, 2] Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 1]
8410: Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 1]
1930: Q[m, 3] Q[p, 2] Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 1]
26170: Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 2]
19690: Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 2]
25090: Q[m, 3] Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 2]
12130: Q[m, 3] Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 1]
17530: Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 2]
11050: Q[m, 3] Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 1]
25990: Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 2]
19510: Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 2]
23830: Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 2]
4390: Q[m, 3] Q[p, 3] Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 1]
16270: Q[m, 3] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 2]
3310: Q[m, 3] Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 1]
24730: Q[m, 3] Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 2]
11770: Q[m, 3] Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 1]
23650: Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 2]
4210: Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 1]
8530: Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 1]
2050: Q[m, 3] Q[p, 3] Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 1]
16990: Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 2]
10510: Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 1]
15910: Q[m, 3] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 2]
2950: Q[m, 3] Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 1]
8350: Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 1]
1870: Q[m, 3] Q[p, 3] Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 1]
35465: Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 3]
28985: Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 2]
34385: Q[p, 3] Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 3]
21425: Q[p, 3] Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 2]
26825: Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 2]
20345: Q[p, 3] Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 2]
35285: Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 3]
28805: Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 2]
33125: Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 3]
13685: Q[p, 3] Q[m, 1] Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 1]
25565: Q[p, 3] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 2]
12605: Q[p, 3] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 1]
34025: Q[p, 3] Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 3]
21065: Q[p, 3] Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 2]
32945: Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 3]
13505: Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 1]
17825: Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 2]
11345: Q[p, 3] Q[m, 1] Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 1]
26285: Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 2]
19805: Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 2]
25205: Q[p, 3] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 2]
12245: Q[p, 3] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 1]
17645: Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 2]
11165: Q[p, 3] Q[m, 1] Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 1]
35435: Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 2] Q[m, 3]
28955: Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 2]
34355: Q[p, 3] Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 2] Q[m, 3]
21395: Q[p, 3] Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 3] Q[m, 2]
26795: Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 2]
20315: Q[p, 3] Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 2] Q[m, 2]
35075: Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 2] Q[m, 3]
28595: Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 2]
31835: Q[p, 3] Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 1] Q[m, 3]
5915: Q[p, 3] Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 3] Q[m, 1]
24275: Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 2]
4835: Q[p, 3] Q[p, 1] Q[m, 2] Q[m, 3] Q[p, 2] Q[m, 1]
33815: Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 2] Q[m, 3]
20855: Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 3] Q[m, 2]
31655: Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 1] Q[m, 3]
5735: Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 3] Q[m, 1]
16535: Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 1] Q[m, 2]
3575: Q[p, 3] Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 2] Q[m, 1]
26075: Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 2]
19595: Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 2] Q[m, 2]
23915: Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 2]
4475: Q[p, 3] Q[p, 1] Q[m, 3] Q[m, 2] Q[p, 2] Q[m, 1]
16355: Q[p, 3] Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 1] Q[m, 2]
3395: Q[p, 3] Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 2] Q[m, 1]
35225: Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 3]
28745: Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 1] Q[m, 3] Q[p, 2]
33065: Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 3]
13625: Q[p, 3] Q[m, 2] Q[m, 1] Q[p, 2] Q[m, 3] Q[p, 1]
25505: Q[p, 3] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 1] Q[p, 2]
12545: Q[p, 3] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 2] Q[p, 1]
35045: Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 3]
28565: Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 1] Q[m, 3] Q[p, 2]
31805: Q[p, 3] Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 3]
5885: Q[p, 3] Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 3] Q[m, 1]
24245: Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 3] Q[m, 1] Q[p, 2]
4805: Q[p, 3] Q[m, 2] Q[p, 1] Q[m, 3] Q[p, 2] Q[m, 1]
32525: Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 3]
13085: Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 1]
31445: Q[p, 3] Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 3]
5525: Q[p, 3] Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 1]
8765: Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 1]
2285: Q[p, 3] Q[m, 2] Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 1]
24785: Q[p, 3] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 2]
11825: Q[p, 3] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 1]
23705: Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 2]
4265: Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 1]
8585: Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 1]
2105: Q[p, 3] Q[m, 2] Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 1]
33935: Q[p, 3] Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Q[m, 3]
20975: Q[p, 3] Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 3] Q[m, 2]
32855: Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Q[m, 3]
13415: Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 2] Q[m, 3] Q[p, 1]
17735: Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 3] Q[p, 1] Q[m, 2]
11255: Q[p, 3] Q[p, 2] Q[m, 1] Q[m, 3] Q[m, 2] Q[p, 1]
33755: Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Q[m, 3]
20795: Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 3] Q[m, 2]
31595: Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Q[m, 3]
5675: Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 3] Q[m, 1]
16475: Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 1] Q[m, 2]
3515: Q[p, 3] Q[p, 2] Q[p, 1] Q[m, 3] Q[m, 2] Q[m, 1]
32495: Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Q[m, 3]
13055: Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 1] Q[m, 3] Q[p, 1]
31415: Q[p, 3] Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Q[m, 3]
5495: Q[p, 3] Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 3] Q[m, 1]
8735: Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 3] Q[m, 1] Q[p, 1]
2255: Q[p, 3] Q[p, 2] Q[m, 2] Q[m, 3] Q[p, 1] Q[m, 1]
17015: Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 1] Q[p, 1] Q[m, 2]
10535: Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 1]
15935: Q[p, 3] Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 1] Q[m, 2]
2975: Q[p, 3] Q[p, 2] Q[m, 3] Q[p, 1] Q[m, 2] Q[m, 1]
8375: Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 1]
1895: Q[p, 3] Q[p, 2] Q[m, 3] Q[m, 2] Q[p, 1] Q[m, 1]
26165: Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 1] Q[m, 2] Q[p, 2]
19685: Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 1] Q[p, 2] Q[m, 2]
25085: Q[p, 3] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 1] Q[p, 2]
12125: Q[p, 3] Q[m, 3] Q[m, 1] Q[m, 2] Q[p, 2] Q[p, 1]
17525: Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 2] Q[p, 1] Q[m, 2]
11045: Q[p, 3] Q[m, 3] Q[m, 1] Q[p, 2] Q[m, 2] Q[p, 1]
25985: Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 1] Q[m, 2] Q[p, 2]
19505: Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 1] Q[p, 2] Q[m, 2]
23825: Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 2] Q[m, 1] Q[p, 2]
4385: Q[p, 3] Q[m, 3] Q[p, 1] Q[m, 2] Q[p, 2] Q[m, 1]
16265: Q[p, 3] Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 1] Q[m, 2]
3305: Q[p, 3] Q[m, 3] Q[p, 1] Q[p, 2] Q[m, 2] Q[m, 1]
24725: Q[p, 3] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 1] Q[p, 2]
11765: Q[p, 3] Q[m, 3] Q[m, 2] Q[m, 1] Q[p, 2] Q[p, 1]
23645: Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 1] Q[m, 1] Q[p, 2]
4205: Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 1] Q[p, 2] Q[m, 1]
8525: Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 2] Q[m, 1] Q[p, 1]
2045: Q[p, 3] Q[m, 3] Q[m, 2] Q[p, 2] Q[p, 1] Q[m, 1]
16985: Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 1] Q[p, 1] Q[m, 2]
10505: Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 1] Q[m, 2] Q[p, 1]
15905: Q[p, 3] Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 1] Q[m, 2]
2945: Q[p, 3] Q[m, 3] Q[p, 2] Q[p, 1] Q[m, 2] Q[m, 1]
8345: Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 2] Q[m, 1] Q[p, 1]
1865: Q[p, 3] Q[m, 3] Q[p, 2] Q[m, 2] Q[p, 1] Q[m, 1]
*/

