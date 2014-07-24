/* The 6 quark amplitudes with 2 distinguished flavors
  See bottom of file for table of switch numbers and values
*/

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;

#define _FAST_COMPILE_eval 0

namespace BH {
#if _FAST_COMPILE_eval==1

template <class T> complex<T> A4q1_2q2zero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A4q1_2q2281_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2286_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2358_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2371_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2383_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2391_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2393_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2406_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2418_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2421_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2423_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2713_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2718_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
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

template <class T> complex<T> A4q1_2q2790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
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

template <class T> complex<T> A4q1_2q2803_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q2838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q2911_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2931_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2933_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2983_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q2995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21003_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21051_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21053_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21063_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q21126_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21138_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21141_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21143_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21198_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21213_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21231_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21233_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q21361_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21366_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21438_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21451_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21463_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21471_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21473_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21486_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21498_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21501_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21541_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21546_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21618_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21661_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21673_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21686_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21688_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21696_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21716_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21718_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21798_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
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

template <class T> complex<T> A4q1_2q21870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
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

template <class T> complex<T> A4q1_2q21883_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21903_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q21918_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21933_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q21973_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q21978_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22093_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22128_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22148_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22171_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22183_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22191_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22193_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22201_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22213_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22226_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22263_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22298_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22341_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22343_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22346_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22348_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q22358_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22386_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22398_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22401_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22416_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22428_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22436_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22438_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22458_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22473_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22488_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22508_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22521_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22523_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22526_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22528_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q22538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22873_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22878_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q22950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q22963_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22983_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q22985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q22998_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q23305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23382_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23407_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23442_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23447_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23523_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23587_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23597_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23643_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q23655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23718_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23733_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23807_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*pow(ep.spb(4,3),
 2))/(ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-(pow(ep.spa(1,2),2)*
pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
ep.spa(2,3)*ep.spb(5,0)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q23835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23953_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q23958_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2)))-
 (pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)))+
 complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2)*
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

template <class T> complex<T> A4q1_2q24030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,1),2))/
 (ep.s(3,4,5)*ep.spa(3,4)*
  ep.spb(1,0)*(-(ep.spa(3,5)*
 ep.spb(3,2))-ep.spa(4,5)*
ep.spb(4,2))))+(pow(ep.spa(0,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*ep.spb(4,3)))+
 complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,1),2)*
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

template <class T> complex<T> A4q1_2q24043_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24063_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
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
   ep.spb(5,0))))-
 complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(1,0),2))/
 (ep.s(3,4,5)*ep.spa(4,5)*
  ep.spb(2,1)*(ep.spa(3,4)*
ep.spb(4,0)+ep.spa(3,5)*
ep.spb(5,0))))+(pow(ep.spa(2,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(1,2)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q24078_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24093_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
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
   ep.spb(4,3))*ep.spb(5,0)))+
 complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0)))-(pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(1,2)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q24133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24138_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24253_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24288_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24308_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24462_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24487_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24522_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24527_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24642_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24697_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24732_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24763_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24775_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24783_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24818_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24847_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24857_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24877_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24892_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24933_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(5,4)))-
 complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(0,1,5)*ep.spa(0,5)*
  ep.spb(3,2)*(ep.spa(0,1)*
ep.spb(4,0)+ep.spa(1,5)*
ep.spb(5,4))))+(pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q24938_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24947_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24978_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q24993_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q24995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25008_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q25028_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q25062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25067_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q25092_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25102_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25113_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*pow(ep.spb(4,3),
  2))/(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spb(3,2)*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4)))-(pow(ep.spa(1,2),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))+
 complex<T>(0,1)*(pow(-(ep.spa(1,5)*ep.spb(1,0))-
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

template <class T> complex<T> A4q1_2q25118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q25125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q25127_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q25132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25231_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25251_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25253_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25303_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25323_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q25371_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25383_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25411_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25423_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25431_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25441_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25453_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25466_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25468_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q25513_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25581_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25583_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25586_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25588_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25593_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25598_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25663_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25683_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q25735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25747_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25757_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25803_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25863_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q25873_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25898_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25927_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25957_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q25972_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26018_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26027_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26032_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26271_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26283_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q26301_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26303_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26306_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26308_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26313_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q26318_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26343_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q26355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(4,5),2)*
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
   ep.spb(5,4))))+
 complex<T>(0,1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
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

template <class T> complex<T> A4q1_2q26378_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q26385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q26387_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q26392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26526_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26541_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26543_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26598_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q26631_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26643_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26706_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26718_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26721_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26736_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26748_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26756_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26758_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26778_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q26808_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26828_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26841_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26846_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26848_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26858_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26958_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26973_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q26975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*pow(ep.spb(4,3),
 2))/(ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-(pow(ep.spa(1,2),2)*
pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
ep.spa(2,3)*ep.spb(5,0)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q27030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27042_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27047_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27063_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27077_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27138_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*pow(ep.spb(4,3),
 2))/(ep.s(0,1,5)*ep.spa(0,5)*
ep.spb(3,2)*(ep.spa(0,1)*
  ep.spb(4,0)+ep.spa(1,5)*
  ep.spb(5,4)))-(pow(ep.spa(1,2),2)*
pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
ep.spa(2,3)*ep.spb(5,0)*
(ep.spa(0,1)*ep.spb(4,0)+
 ep.spa(1,5)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q27168_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27188_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27222_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27227_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27287_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27290_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27292_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27351_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27363_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q27381_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27383_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27386_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27388_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27393_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q27398_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27423_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q27435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27437_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27453_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
  2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2)))-(pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*ep.spb(4,3)))-
 complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2)*
  (ep.spa(0,3)*ep.spb(2,0)+
   ep.spa(1,3)*ep.spb(2,1)))/
 (ep.s(0,1,2)*ep.spa(3,4)*
  (-(ep.spa(1,3)*ep.spb(1,0))-
   ep.spa(2,3)*ep.spb(2,0))*
  ep.spb(2,1)*(-(ep.spa(3,5)*
 ep.spb(3,2))-ep.spa(4,5)*
ep.spb(4,2))))-(pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2)*(ep.spa(2,5)*
   ep.spb(4,2)+ep.spa(3,5)*
   ep.spb(4,3)))/(ep.s(2,3,4)*
 ep.spa(0,5)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*ep.spb(4,3)*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3)))+
 (pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
 (ep.spa(1,2)*ep.spb(2,0)+
  ep.spa(1,3)*ep.spb(3,0)))/
(ep.s(1,2,3)*ep.spa(1,2)*
 (-(ep.spa(1,3)*ep.spb(1,0))-
  ep.spa(2,3)*ep.spb(2,0))*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3))*
 ep.spb(5,0)))
); }

template <class T> complex<T> A4q1_2q27458_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q27465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q27467_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q27472_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27841_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27846_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27913_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27918_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27931_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27943_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27951_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27953_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27966_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27978_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27981_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q27983_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28021_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28026_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28093_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28098_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28141_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28166_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28168_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28176_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28188_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28196_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28198_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28363_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28383_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28398_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28413_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28453_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28458_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
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

template <class T> complex<T> A4q1_2q28530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
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

template <class T> complex<T> A4q1_2q28573_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28598_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q28608_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28628_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q28651_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28663_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28671_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28673_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28681_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28693_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28706_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28743_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28753_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28778_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28821_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28826_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28828_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q28866_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28878_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28881_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28883_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28896_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28908_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28916_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28918_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28938_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28953_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28968_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28988_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q28990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29001_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29003_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29006_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29008_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29018_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q29101_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29106_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29173_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29178_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29221_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29233_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29246_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29248_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29256_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29268_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29276_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
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

template <class T> complex<T> A4q1_2q29610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(-(ep.spa(2,4)*ep.spb(2,1))-
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

template <class T> complex<T> A4q1_2q29653_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29678_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q29688_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q29941_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29953_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29966_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q29968_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210038_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210116_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210128_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q210156_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210168_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210176_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210178_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210248_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210296_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210298_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210308_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q210433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210438_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q210510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q210523_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210543_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q210558_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210573_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q210613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210618_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2)))-
 (pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)))+
 complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2)*
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

template <class T> complex<T> A4q1_2q210690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,1),2))/
 (ep.s(3,4,5)*ep.spa(3,4)*
  ep.spb(1,0)*(-(ep.spa(3,5)*
 ep.spb(3,2))-ep.spa(4,5)*
ep.spb(4,2))))+(pow(ep.spa(0,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*ep.spb(4,3)))+
 complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,1),2)*
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

template <class T> complex<T> A4q1_2q210733_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210758_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
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
   ep.spb(5,0))))-
 complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(1,0),2))/
 (ep.s(3,4,5)*ep.spa(4,5)*
  ep.spb(2,1)*(ep.spa(3,4)*
ep.spb(4,0)+ep.spa(3,5)*
ep.spb(5,0))))+(pow(ep.spa(2,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(1,2)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q210768_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210788_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2)*
  (ep.spa(0,3)*ep.spb(2,0)+
   ep.spa(1,3)*ep.spb(2,1)))/
 (ep.s(0,1,2)*ep.spa(3,4)*
  (-(ep.spa(1,3)*ep.spb(1,0))-
   ep.spa(2,3)*ep.spb(2,0))*
  ep.spb(2,1)*(-(ep.spa(3,5)*
 ep.spb(3,2))-ep.spa(4,5)*
ep.spb(4,2))))-(pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2)*(ep.spa(2,5)*
   ep.spb(4,2)+ep.spa(3,5)*
   ep.spb(4,3)))/(ep.s(2,3,4)*
 ep.spa(0,5)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*ep.spb(4,3)*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3)))+
 (pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
 (ep.spa(1,2)*ep.spb(2,0)+
  ep.spa(1,3)*ep.spb(3,0)))/
(ep.s(1,2,3)*ep.spa(1,2)*
 (-(ep.spa(1,3)*ep.spb(1,0))-
  ep.spa(2,3)*ep.spb(2,0))*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3))*
 ep.spb(5,0)))+
 complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0)))-(pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(1,2)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q210865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q210870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q210937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210942_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q210977_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q210990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211002_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211007_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211117_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211122_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211212_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211222_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211263_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211298_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211327_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211413_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211418_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(5,4)))-
 complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(0,1,5)*ep.spa(0,5)*
  ep.spb(3,2)*(ep.spa(0,1)*
ep.spb(4,0)+ep.spa(1,5)*
ep.spb(5,4))))+(pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q211425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211427_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211432_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211458_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211473_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211488_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211508_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211542_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211547_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211572_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211582_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211593_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211598_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*
 pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4)))-
 (pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spb(5,0)*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4))))+
 complex<T>(0,1)*(pow(-(ep.spa(1,5)*ep.spb(1,0))-
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

template <class T> complex<T> A4q1_2q211605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q211612_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211693_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211698_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q211770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q211813_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q211848_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211868_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q211870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q212125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212197_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212202_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212257_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212292_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212302_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212558_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212632_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q212720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212722_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212748_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212768_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212842_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212888_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q212900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212971_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212983_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212991_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q212993_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213001_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213026_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213028_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213043_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213063_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213073_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213098_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q213141_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213143_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213146_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213148_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213158_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213181_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213193_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213206_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213208_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213253_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q213356_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213358_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213368_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213423_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213458_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q213475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213487_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213517_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213573_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213578_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213587_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213592_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213638_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q213685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213697_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213788_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q213802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214041_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214043_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214046_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214048_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214053_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214058_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q214076_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214078_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214088_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q214113_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q214118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(4,5),2)*pow(ep.spb(2,0),2)*
  (ep.spa(0,3)*ep.spb(2,0)+
   ep.spa(1,3)*ep.spb(2,1)))/
 (ep.s(0,1,2)*ep.spa(3,4)*
  ep.spb(2,1)*(ep.spa(0,5)*
ep.spb(2,0)+ep.spa(1,5)*
ep.spb(2,1))*(ep.spa(3,4)*
ep.spb(4,0)+ep.spa(3,5)*
ep.spb(5,0))))+(pow(ep.spa(1,5),2)*
 pow(ep.spb(3,2),2)*(ep.spa(0,5)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(4,1)))/(ep.s(0,1,5)*
 ep.spa(0,5)*(ep.spa(0,5)*
   ep.spb(2,0)+ep.spa(1,5)*
   ep.spb(2,1))*ep.spb(4,3)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4)))-
 pow(ep.spa(1,4)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,0),3)/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spb(5,0)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))*(ep.spa(0,1)*
   ep.spb(4,0)+ep.spa(1,5)*
   ep.spb(5,4))))+
 complex<T>(0,1)*(pow(-(ep.spa(2,4)*ep.spb(2,0))-
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

template <class T> complex<T> A4q1_2q214125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q214127_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q214132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214148_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(4,5),2)*
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

template <class T> complex<T> A4q1_2q214160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214266_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214281_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214283_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214296_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214308_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214316_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214318_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214338_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214368_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214388_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q214401_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214406_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214408_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214413_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214418_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214476_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214488_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214496_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214498_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214548_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214568_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q214616_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214618_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214628_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214698_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214713_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214728_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214748_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q214770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214782_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214787_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214812_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214847_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214908_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214928_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q214980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q214992_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215002_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215048_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215121_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215123_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215126_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215128_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215138_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q215156_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215158_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215168_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(2,5)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q215193_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q215198_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2)))-
 (pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 (-(ep.spa(3,5)*ep.spb(3,2))-
  ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3)))-
 complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2)*
  (ep.spa(0,3)*ep.spb(2,0)+
   ep.spa(1,3)*ep.spb(2,1)))/
 (ep.s(0,1,2)*ep.spa(3,4)*
  (-(ep.spa(1,3)*ep.spb(1,0))-
   ep.spa(2,3)*ep.spb(2,0))*
  ep.spb(2,1)*(-(ep.spa(3,5)*
 ep.spb(3,2))-ep.spa(4,5)*
ep.spb(4,2))))-(pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2)*(ep.spa(2,5)*
   ep.spb(4,2)+ep.spa(3,5)*
   ep.spb(4,3)))/(ep.s(2,3,4)*
 ep.spa(0,5)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*ep.spb(4,3)*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3)))+
 (pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
 (ep.spa(1,2)*ep.spb(2,0)+
  ep.spa(1,3)*ep.spb(3,0)))/
(ep.s(1,2,3)*ep.spa(1,2)*
 (-(ep.spa(1,3)*ep.spb(1,0))-
  ep.spa(2,3)*ep.spb(2,0))*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3))*
 ep.spb(5,0)))
); }

template <class T> complex<T> A4q1_2q215205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q215207_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q215212_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q215240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215242_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215943_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215958_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215973_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q215975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(2,4),2)*
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

template <class T> complex<T> A4q1_2q216270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(2,4),2)*
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

template <class T> complex<T> A4q1_2q216337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216367_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,2),2)*
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

template <class T> complex<T> A4q1_2q216377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
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

template <class T> complex<T> A4q1_2q216407_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216463_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216547_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216557_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216603_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q216617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216678_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216693_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216767_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216783_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q216797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216913_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216918_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q216985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q216990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217003_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217023_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217038_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217053_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217093_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217098_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217213_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217238_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217248_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217268_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,1)*(-((pow(ep.spa(2,4),2)*pow(ep.spb(5,0),2))/
 (ep.s(2,3,4)*ep.spa(3,4)*
  ep.spb(1,0)*(-(ep.spa(2,3)*
 ep.spb(5,3))-ep.spa(2,4)*
ep.spb(5,4))))+(pow(ep.spa(1,2),2)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spb(4,3)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))+
 complex<T>(0,1)*(pow(ep.spa(1,4)*ep.spb(3,1)+
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

template <class T> complex<T> A4q1_2q217350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))-
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spb(4,3)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))+
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
  (ep.spa(1,4)*ep.spb(3,1)+
   ep.spa(2,4)*ep.spb(3,2)))/
 (ep.s(1,2,3)*ep.spa(4,5)*
  (-(ep.spa(2,4)*ep.spb(2,1))-
   ep.spa(3,4)*ep.spb(3,1))*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))-(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2)*(ep.spa(2,3)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(4,1)))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spb(1,0)*
 (-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1))*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))+
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
 (-(ep.spa(0,3)*ep.spb(5,3))-
  ep.spa(0,4)*ep.spb(5,4)))/
(ep.s(3,4,5)*ep.spa(0,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q217417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217447_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
  (ep.spa(0,2)*ep.spb(3,0)+
   ep.spa(1,2)*ep.spb(3,1)))/
 (ep.s(0,1,2)*ep.spa(1,2)*
  (-(ep.spa(0,1)*ep.spb(3,1))-
   ep.spa(0,2)*ep.spb(3,2))*
  ep.spb(4,3)*(ep.spa(0,2)*
ep.spb(5,0)+ep.spa(1,2)*
ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(0,2)*
   ep.spb(2,1)+ep.spa(0,3)*
   ep.spb(3,1)))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spb(2,1)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1)))-
 (pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
 (ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1)))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
 complex<T>(0,1)*((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3)))-(pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(1,2)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q217457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217482_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
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
   ep.spb(5,3))))-
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,2),2))/
 (ep.s(0,4,5)*ep.spa(4,5)*
  ep.spb(2,1)*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))+(pow(ep.spa(0,1),2)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(1,2)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q217487_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217597_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217602_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217672_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217692_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217702_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217743_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217753_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217778_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217807_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217893_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217898_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
 pow(ep.spb(2,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1)))-
 (pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spb(5,0)*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(2,1),2)*
  (ep.spa(0,4)*ep.spb(3,0)+
   ep.spa(4,5)*ep.spb(5,3)))/
 (ep.s(0,4,5)*ep.spa(4,5)*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(1,0)+ep.spa(4,5)*
ep.spb(5,1))*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))-(pow(ep.spa(3,4),2)*
 pow(ep.spb(5,1),2)*(ep.spa(0,2)*
   ep.spb(1,0)+ep.spa(2,5)*
   ep.spb(5,1)))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spb(1,0)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))+
 pow(-(ep.spa(0,3)*ep.spb(5,3))-
  ep.spa(0,4)*ep.spb(5,4),3)/
(ep.s(3,4,5)*ep.spa(0,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q217907_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217912_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217938_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217953_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217968_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q217988_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q217990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q218022_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218027_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q218052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218073_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q218078_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q218085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(1,0)+ep.spa(4,5)*
ep.spb(5,1))))+(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
  (ep.spa(1,4)*ep.spb(3,1)+
   ep.spa(2,4)*ep.spb(3,2)))/
 (ep.s(1,2,3)*ep.spa(4,5)*
  (-(ep.spa(2,4)*ep.spb(2,1))-
   ep.spa(3,4)*ep.spb(3,1))*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))-(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2)*(ep.spa(2,3)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(4,1)))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spb(1,0)*
 (-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1))*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))+
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
 (-(ep.spa(0,3)*ep.spb(5,3))-
  ep.spa(0,4)*ep.spb(5,4)))/
(ep.s(3,4,5)*ep.spa(0,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q218087_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q218092_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
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

template <class T> complex<T> A4q1_2q218430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
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

template <class T> complex<T> A4q1_2q218497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218527_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q218537_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218562_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q218567_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218857_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218862_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218929_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218934_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218947_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218959_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218969_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218982_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218994_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218997_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q218999_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219067_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219077_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219127_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219139_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219147_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219149_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q219197_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219207_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219209_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219282_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219287_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219354_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219359_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q219377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219387_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219389_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
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

template <class T> complex<T> A4q1_2q219510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
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

template <class T> complex<T> A4q1_2q219577_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219582_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q219617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219642_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q219647_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219757_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219862_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q219942_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220009_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220014_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220027_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220039_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220047_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220049_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220074_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220077_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220079_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220117_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220122_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220189_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220194_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220249_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220264_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220284_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220292_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220294_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220327_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220387_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220399_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220407_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220409_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220429_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220442_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220444_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q220487_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220492_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220499_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220504_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220542_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220547_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220572_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220582_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220602_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220614_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220619_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220632_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220644_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220652_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220654_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q220667_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220672_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220677_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220679_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220684_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220783_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220803_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220867_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q220877_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220963_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220983_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q220993_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221018_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221047_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q221057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221077_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221092_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221138_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221147_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221152_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221227_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q221237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221287_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221299_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221307_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221309_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221367_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221369_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221407_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q221417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221437_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221452_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221467_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221479_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221487_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221489_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221509_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221522_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221524_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221567_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221572_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221577_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221579_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221582_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221584_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,2),2)*
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

template <class T> complex<T> A4q1_2q221837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q221858_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q221865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
  (ep.spa(0,2)*ep.spb(3,0)+
   ep.spa(1,2)*ep.spb(3,1)))/
 (ep.s(0,1,2)*ep.spa(1,2)*
  (-(ep.spa(0,1)*ep.spb(3,1))-
   ep.spa(0,2)*ep.spb(3,2))*
  ep.spb(4,3)*(ep.spa(0,2)*
ep.spb(5,0)+ep.spa(1,2)*
ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(0,2)*
   ep.spb(2,1)+ep.spa(0,3)*
   ep.spb(3,1)))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spb(2,1)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1)))-
 (pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
 (ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1)))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
 complex<T>(0,1)*((pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4)))-(pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spb(4,3)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q221867_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q221872_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
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

template <class T> complex<T> A4q1_2q221897_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221907_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221909_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
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

template <class T> complex<T> A4q1_2q221927_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221939_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221942_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q221944_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222078_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222093_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q222167_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222183_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222197_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222258_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222288_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222308_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q222347_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222382_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222393_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222398_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222407_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222522_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q222527_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222582_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222594_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222597_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222599_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222627_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222629_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222702_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q222707_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222732_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222774_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222777_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222779_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222792_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222804_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222812_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222814_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222827_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222839_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222842_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222844_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222903_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q222917_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222933_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q222938_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q222945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(4,3)))-
 complex<T>(0,1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q222947_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q222952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q222977_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222987_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q222989_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q223007_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223012_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223017_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223019_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223022_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223024_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223393_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223398_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223518_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223573_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223578_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223693_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223718_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223728_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223748_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223897_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223927_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223962_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q223965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q223967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,1)*(-((pow(ep.spa(2,4),2)*pow(ep.spb(5,0),2))/
 (ep.s(2,3,4)*ep.spa(3,4)*
  ep.spb(1,0)*(-(ep.spa(2,3)*
 ep.spb(5,3))-ep.spa(2,4)*
ep.spb(5,4))))+(pow(ep.spa(1,2),2)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spb(4,3)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))+
 complex<T>(0,1)*(pow(ep.spa(1,4)*ep.spb(3,1)+
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

template <class T> complex<T> A4q1_2q224010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spb(1,0)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))-
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spb(4,3)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))+
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
  (ep.spa(1,4)*ep.spb(3,1)+
   ep.spa(2,4)*ep.spb(3,2)))/
 (ep.s(1,2,3)*ep.spa(4,5)*
  (-(ep.spa(2,4)*ep.spb(2,1))-
   ep.spa(3,4)*ep.spb(3,1))*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))-(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2)*(ep.spa(2,3)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(4,1)))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spb(1,0)*
 (-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1))*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))+
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
 (-(ep.spa(0,3)*ep.spb(5,3))-
  ep.spa(0,4)*ep.spb(5,4)))/
(ep.s(3,4,5)*ep.spa(0,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q224077_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224082_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224137_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
  (ep.spa(0,2)*ep.spb(3,0)+
   ep.spa(1,2)*ep.spb(3,1)))/
 (ep.s(0,1,2)*ep.spa(1,2)*
  (-(ep.spa(0,1)*ep.spb(3,1))-
   ep.spa(0,2)*ep.spb(3,2))*
  ep.spb(4,3)*(ep.spa(0,2)*
ep.spb(5,0)+ep.spa(1,2)*
ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(0,2)*
   ep.spb(2,1)+ep.spa(0,3)*
   ep.spb(3,1)))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spb(2,1)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1)))-
 (pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
 (ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1)))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
 complex<T>(0,1)*((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3)))-(pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(1,2)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q224152_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224172_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
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
   ep.spb(5,3))))-
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,2),2))/
 (ep.s(0,4,5)*ep.spa(4,5)*
  ep.spb(2,1)*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))+(pow(ep.spa(0,1),2)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(1,2)*(ep.spa(0,4)*
   ep.spb(4,3)+ep.spa(0,5)*
   ep.spb(5,3))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q224182_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224203_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224223_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224233_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224258_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224287_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224297_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224317_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224332_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224378_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224387_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
 pow(ep.spb(2,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1)))-
 (pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spb(5,0)*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(2,1),2)*
  (ep.spa(0,4)*ep.spb(3,0)+
   ep.spa(4,5)*ep.spb(5,3)))/
 (ep.s(0,4,5)*ep.spa(4,5)*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(1,0)+ep.spa(4,5)*
ep.spb(5,1))*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))-(pow(ep.spa(3,4),2)*
 pow(ep.spb(5,1),2)*(ep.spa(0,2)*
   ep.spb(1,0)+ep.spa(2,5)*
   ep.spb(5,1)))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spb(1,0)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))+
 pow(-(ep.spa(0,3)*ep.spb(5,3))-
  ep.spa(0,4)*ep.spb(5,4),3)/
(ep.s(3,4,5)*ep.spa(0,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q224392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224418_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224448_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224468_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224507_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224542_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224553_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224558_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q224567_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(1,0)+ep.spa(4,5)*
ep.spb(5,1))))+(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
  (ep.spa(1,4)*ep.spb(3,1)+
   ep.spa(2,4)*ep.spb(3,2)))/
 (ep.s(1,2,3)*ep.spa(4,5)*
  (-(ep.spa(2,4)*ep.spb(2,1))-
   ep.spa(3,4)*ep.spb(3,1))*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))-(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2)*(ep.spa(2,3)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(4,1)))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spb(1,0)*
 (-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1))*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))+
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
 (-(ep.spa(0,3)*ep.spb(5,3))-
  ep.spa(0,4)*ep.spb(5,4)))/
(ep.s(3,4,5)*ep.spa(0,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q224572_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224653_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224658_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224773_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224798_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224808_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224828_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q224830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(2,4),2)*
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

template <class T> complex<T> A4q1_2q225090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(2,4),2)*
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

template <class T> complex<T> A4q1_2q225157_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225217_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,2),2)*
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

template <class T> complex<T> A4q1_2q225232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
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

template <class T> complex<T> A4q1_2q225262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225493_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225518_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225577_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225592_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225668_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q225682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225728_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225792_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225848_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q225862_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q225990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226087_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226097_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226122_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226127_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
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

template <class T> complex<T> A4q1_2q226170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
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

template <class T> complex<T> A4q1_2q226237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226242_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226297_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q226312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226332_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q226342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226489_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226494_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226507_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226519_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226527_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226529_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226542_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226554_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226557_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226559_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226597_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226602_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226669_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226674_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226717_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226729_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226744_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226752_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226764_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226772_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226774_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226807_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226867_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226879_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226887_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226889_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226897_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226909_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226922_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226924_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q226972_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226977_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226979_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226982_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q226984_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227022_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227027_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227082_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227094_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227097_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227099_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227112_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227124_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227134_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227147_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q227152_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227157_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227159_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227164_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(3,4)*ep.spb(4,0)+
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

template <class T> complex<T> A4q1_2q227250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(3,4)*ep.spb(4,1)+
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

template <class T> complex<T> A4q1_2q227317_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q227392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q227422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227677_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227749_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227754_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227809_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227824_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227844_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q227854_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228097_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228112_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228157_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228169_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228182_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228184_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q228262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228274_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228384_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228394_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q228442_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228452_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228454_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228523_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228543_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228553_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228578_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q228652_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228693_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228698_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228707_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228733_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228758_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q228832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228908_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228922_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228977_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q228997_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q229012_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229027_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229039_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229047_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229049_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229069_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229082_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229084_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229127_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229137_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229139_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229144_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q229192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229249_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229264_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229354_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229593_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q229598_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q229605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q229607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
  (ep.spa(0,2)*ep.spb(3,0)+
   ep.spa(1,2)*ep.spb(3,1)))/
 (ep.s(0,1,2)*ep.spa(1,2)*
  (-(ep.spa(0,1)*ep.spb(3,1))-
   ep.spa(0,2)*ep.spb(3,2))*
  ep.spb(4,3)*(ep.spa(0,2)*
ep.spb(5,0)+ep.spa(1,2)*
ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(0,2)*
   ep.spb(2,1)+ep.spa(0,3)*
   ep.spb(3,1)))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spb(2,1)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1)))-
 (pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
 (ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1)))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))+
 complex<T>(0,1)*((pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4)))-(pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spb(4,3)*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q229612_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229628_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,2),2)*
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

template <class T> complex<T> A4q1_2q229642_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229667_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
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

template <class T> complex<T> A4q1_2q229672_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229677_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229679_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229684_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(5,0)+
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

template <class T> complex<T> A4q1_2q229702_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229714_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229818_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229848_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229868_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229907_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q229942_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229953_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229958_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q229972_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230028_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230048_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230112_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q230122_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230168_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230182_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230267_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230292_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q230302_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230334_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230339_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230364_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230374_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230387_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230397_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230399_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230404_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230472_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230480_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q230482_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230544_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230552_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230554_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230602_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230612_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230614_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230673_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q230678_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q230685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q230687_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(4,3)))-
 complex<T>(0,1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q230692_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q230722_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230747_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q230752_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230757_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230759_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230764_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,3)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q230782_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230792_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q230794_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231151_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231163_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231171_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231173_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231223_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q231243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231291_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231293_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231303_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231331_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231343_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231351_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231361_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231386_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231388_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q231423_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231458_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231501_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231506_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231508_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231513_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231518_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231583_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q231603_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231667_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231677_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231737_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231763_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231775_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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
   ep.spb(5,4))))+
 complex<T>(0,1)*((pow(ep.spa(4,5),2)*pow(ep.spb(3,1),2)*
 (ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q231783_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q231818_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q231847_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231857_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q231877_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231892_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231933_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q231938_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q231945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q231947_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q231950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q231952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232191_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232193_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232203_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232221_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232223_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232226_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232233_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232238_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232263_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232277_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232293_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232298_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232307_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232411_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232423_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232431_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232441_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232453_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232466_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232468_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232513_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q232538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232581_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232583_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232586_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232588_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232593_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232598_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232621_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232646_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232648_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232693_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,5)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q232718_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232796_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232798_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232808_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232810_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q232863_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232873_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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
   ep.spb(5,4))))+
 complex<T>(0,1)*((pow(ep.spa(4,5),2)*pow(ep.spb(3,1),2)*
 (ep.spa(0,2)*ep.spb(2,1)+
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

template <class T> complex<T> A4q1_2q232898_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q232927_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q232957_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q232972_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q233018_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q233025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q233027_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q233032_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233053_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,5)*ep.spb(3,0)+
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

template <class T> complex<T> A4q1_2q233078_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233137_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233152_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233242_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233481_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233486_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233488_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233493_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233498_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233516_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233518_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233528_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233553_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233558_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233567_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233572_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233588_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233602_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233743_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233763_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q233827_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233883_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233897_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q233943_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233953_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q233978_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q233995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2)*
  (ep.spa(0,3)*ep.spb(2,0)+
   ep.spa(1,3)*ep.spb(2,1)))/
 (ep.s(0,1,2)*ep.spa(3,4)*
  (-(ep.spa(1,3)*ep.spb(1,0))-
   ep.spa(2,3)*ep.spb(2,0))*
  ep.spb(2,1)*(-(ep.spa(3,5)*
 ep.spb(3,2))-ep.spa(4,5)*
ep.spb(4,2))))-(pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2)*(ep.spa(2,5)*
   ep.spb(4,2)+ep.spa(3,5)*
   ep.spb(4,3)))/(ep.s(2,3,4)*
 ep.spa(0,5)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*ep.spb(4,3)*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3)))+
 (pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
 (ep.spa(1,2)*ep.spb(2,0)+
  ep.spa(1,3)*ep.spb(3,0)))/
(ep.s(1,2,3)*ep.spa(1,2)*
 (-(ep.spa(1,3)*ep.spb(1,0))-
  ep.spa(2,3)*ep.spb(2,0))*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3))*
 ep.spb(5,0)))-
 complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(0,1,5)*ep.spa(0,5)*
  ep.spb(3,2)*(ep.spa(0,1)*
ep.spb(4,0)+ep.spa(1,5)*
ep.spb(5,4))))+(pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q234007_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234017_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q234037_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234093_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q234098_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q234105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q234107_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q234112_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234175_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q234187_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234197_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234247_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234259_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234267_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234269_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234317_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234327_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234329_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q234367_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234397_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234427_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234439_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234447_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234449_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234469_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234482_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234484_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234527_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234537_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234539_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234542_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234544_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234783_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234813_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234818_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234827_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234857_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234867_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234869_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234887_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234892_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234897_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234899_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q234904_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235003_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q235023_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235033_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q235058_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q235087_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235097_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2)*
  (ep.spa(0,3)*ep.spb(2,0)+
   ep.spa(1,3)*ep.spb(2,1)))/
 (ep.s(0,1,2)*ep.spa(3,4)*
  (-(ep.spa(1,3)*ep.spb(1,0))-
   ep.spa(2,3)*ep.spb(2,0))*
  ep.spb(2,1)*(-(ep.spa(3,5)*
 ep.spb(3,2))-ep.spa(4,5)*
ep.spb(4,2))))-(pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),2)*(ep.spa(2,5)*
   ep.spb(4,2)+ep.spa(3,5)*
   ep.spb(4,3)))/(ep.s(2,3,4)*
 ep.spa(0,5)*(-(ep.spa(3,5)*
ep.spb(3,2))-ep.spa(4,5)*
   ep.spb(4,2))*ep.spb(4,3)*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3)))+
 (pow(ep.spa(1,3),2)*pow(ep.spb(4,0),2)*
 (ep.spa(1,2)*ep.spb(2,0)+
  ep.spa(1,3)*ep.spb(3,0)))/
(ep.s(1,2,3)*ep.spa(1,2)*
 (-(ep.spa(1,3)*ep.spb(1,0))-
  ep.spa(2,3)*ep.spb(2,0))*
 (-(ep.spa(1,2)*ep.spb(4,2))-
  ep.spa(1,3)*ep.spb(4,3))*
 ep.spb(5,0)))-
 complex<T>(0,1)*(-((pow(ep.spa(1,5),2)*pow(ep.spb(4,2),2))/
 (ep.s(0,1,5)*ep.spa(0,5)*
  ep.spb(3,2)*(ep.spa(0,1)*
ep.spb(4,0)+ep.spa(1,5)*
ep.spb(5,4))))+(pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,1)*ep.spb(4,0)+
  ep.spa(1,5)*ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q235117_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235173_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235175_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q235178_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q235185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q235187_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q235192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235213_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235238_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q235297_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235388_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235447_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q235477_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235492_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235507_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235519_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235527_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235529_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235537_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235549_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235562_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235564_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235612_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235619_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235624_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(1,5),2)*
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

template <class T> complex<T> A4q1_2q235657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235672_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235717_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235729_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235744_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q235834_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236073_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236078_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236087_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236092_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236108_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236122_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236147_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236152_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236157_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236159_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236164_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236182_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q236194_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237591_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237593_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237603_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q237621_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237623_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237626_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237628_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q237638_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237663_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q237675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237677_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237693_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(5,4)))-
 complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0)))-(pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(1,2)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q237698_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q237705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q237707_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q237712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237801_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237803_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237806_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237808_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237813_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237818_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q237836_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237848_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q237873_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q237878_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237880_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(5,4)))-
 complex<T>(0,1)*((pow(ep.spa(3,5),2)*pow(ep.spb(2,0),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spb(2,1)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0)))-(pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),2))/(ep.s(0,4,5)*
 ep.spa(1,2)*(ep.spa(3,4)*
   ep.spb(4,0)+ep.spa(3,5)*
   ep.spb(5,0))*ep.spb(5,4)))
); }

template <class T> complex<T> A4q1_2q237885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q237887_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q237892_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237908_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(3,5),2)*
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

template <class T> complex<T> A4q1_2q237920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q237922_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238023_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
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

template <class T> complex<T> A4q1_2q238037_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238053_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q238058_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q238065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
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
 ep.spb(5,4)))+
 complex<T>(0,1)*(-((pow(ep.spa(0,5),2)*pow(ep.spb(3,1),2)*
  (ep.spa(1,4)*ep.spb(3,1)+
   ep.spa(2,4)*ep.spb(3,2)))/
 (ep.s(1,2,3)*ep.spa(4,5)*
  ep.spb(3,2)*(-(ep.spa(0,1)*
 ep.spb(3,1))-ep.spa(0,2)*
ep.spb(3,2))*(ep.spa(0,4)*
ep.spb(1,0)+ep.spa(4,5)*
ep.spb(5,1))))+
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

template <class T> complex<T> A4q1_2q238067_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q238072_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
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

template <class T> complex<T> A4q1_2q238097_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238107_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238109_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
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

template <class T> complex<T> A4q1_2q238127_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238137_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238139_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238144_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238233_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q238238_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q238245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q238247_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
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
 ep.spb(5,4)))+
 complex<T>(0,1)*(-((pow(ep.spa(0,5),2)*pow(ep.spb(3,1),2)*
  (ep.spa(1,4)*ep.spb(3,1)+
   ep.spa(2,4)*ep.spb(3,2)))/
 (ep.s(1,2,3)*ep.spa(4,5)*
  ep.spb(3,2)*(-(ep.spa(0,1)*
 ep.spb(3,1))-ep.spa(0,2)*
ep.spb(3,2))*(ep.spa(0,4)*
ep.spb(1,0)+ep.spa(4,5)*
ep.spb(5,1))))+
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

template <class T> complex<T> A4q1_2q238252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238268_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(-(ep.spa(3,5)*ep.spb(3,1))-
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

template <class T> complex<T> A4q1_2q238282_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238307_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
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

template <class T> complex<T> A4q1_2q238312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238317_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238319_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238324_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(-(ep.spa(3,5)*ep.spb(3,2))-
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

template <class T> complex<T> A4q1_2q238342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238354_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238926_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238938_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238941_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238943_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q238998_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q239013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239031_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239033_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239043_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239106_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239121_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239123_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239136_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239148_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239156_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239158_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239178_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q239193_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239208_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239241_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239246_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239248_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239253_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239258_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239358_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q239373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239442_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239447_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239463_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239477_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(1,0)+ep.spa(4,5)*
ep.spb(5,1))))+(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
 complex<T>(0,1)*(-((pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
  (ep.spa(0,2)*ep.spb(3,0)+
   ep.spa(1,2)*ep.spb(3,1)))/
 (ep.s(0,1,2)*ep.spa(1,2)*
  (-(ep.spa(0,1)*ep.spb(3,1))-
   ep.spa(0,2)*ep.spb(3,2))*
  ep.spb(4,3)*(ep.spa(0,2)*
ep.spb(5,0)+ep.spa(1,2)*
ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(0,2)*
   ep.spb(2,1)+ep.spa(0,3)*
   ep.spb(3,1)))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spb(2,1)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1)))-
 (pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
 (ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1)))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))
); }

template <class T> complex<T> A4q1_2q239553_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239568_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q239588_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q239622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239627_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q239652_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239662_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239673_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q239678_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q239685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q239687_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q239692_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239751_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239753_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239763_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239781_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239783_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239786_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239788_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239798_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239858_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239867_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q239872_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240186_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240198_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240201_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240203_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240216_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240236_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240238_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240258_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240288_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q240308_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240321_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240323_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240326_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240328_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240333_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240338_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240396_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240408_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240416_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240418_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240468_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240480_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q240488_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240536_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240548_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240618_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q240633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240648_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(1,0)+ep.spa(4,5)*
ep.spb(5,1))))+(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spb(5,0)*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1))))+
 complex<T>(0,1)*(-((pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
  (ep.spa(0,2)*ep.spb(3,0)+
   ep.spa(1,2)*ep.spb(3,1)))/
 (ep.s(0,1,2)*ep.spa(1,2)*
  (-(ep.spa(0,1)*ep.spb(3,1))-
   ep.spa(0,2)*ep.spb(3,2))*
  ep.spb(4,3)*(ep.spa(0,2)*
ep.spb(5,0)+ep.spa(1,2)*
ep.spb(5,1))))+(pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2)*(ep.spa(0,2)*
   ep.spb(2,1)+ep.spa(0,3)*
   ep.spb(3,1)))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spb(2,1)*
 (-(ep.spa(0,1)*ep.spb(3,1))-
  ep.spa(0,2)*ep.spb(3,2))*
 (ep.spa(0,4)*ep.spb(1,0)+
  ep.spa(4,5)*ep.spb(5,1)))-
 (pow(ep.spa(2,4),2)*pow(ep.spb(5,1),2)*
 (ep.spa(0,4)*ep.spb(5,0)+
  ep.spa(1,4)*ep.spb(5,1)))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spb(5,0)*(ep.spa(0,2)*
   ep.spb(5,0)+ep.spa(1,2)*
   ep.spb(5,1))*(ep.spa(0,4)*
   ep.spb(1,0)+ep.spa(4,5)*
   ep.spb(5,1))))
); }

template <class T> complex<T> A4q1_2q240668_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q240702_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240707_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q240732_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240753_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q240758_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q240765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q240767_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q240772_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240828_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q240848_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240912_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240922_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240968_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q240982_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241041_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241043_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241046_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241048_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241053_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241058_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241076_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241078_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241088_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241113_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241127_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241148_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241518_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q241602_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241623_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241698_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q241713_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241728_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q241748_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(5,0)))-
 complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q241782_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241787_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q241812_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q241838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q241845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q241847_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q241852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q241962_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q241967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242022_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242034_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242037_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242039_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242067_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242069_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q242142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242147_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242172_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242182_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242202_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242214_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242217_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242219_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242244_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242254_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242267_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242277_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242279_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242282_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242284_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242343_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242378_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242387_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242427_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242429_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242447_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242452_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242459_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242462_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242464_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242778_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q242793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242808_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q242828_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q242862_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242867_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242880_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(5,0)))-
 complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q242892_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242913_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q242918_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q242925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q242927_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q242932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q242988_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243008_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(-(pow(ep.spa(0,4)*ep.spb(4,2)+
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

template <class T> complex<T> A4q1_2q243072_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243082_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243128_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243222_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243227_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q243252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243282_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243294_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243297_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243299_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243324_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243332_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243334_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243347_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243359_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243362_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243364_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(4,3)+
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

template <class T> complex<T> A4q1_2q243432_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243442_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243492_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243504_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243512_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243514_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243562_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243572_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243574_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243638_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243647_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243652_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243668_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243707_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243717_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243719_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243722_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243724_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243752_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q243754_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244071_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244073_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244083_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
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

template <class T> complex<T> A4q1_2q244101_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244103_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244106_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244108_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244113_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
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

template <class T> complex<T> A4q1_2q244118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244143_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q244155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244157_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244173_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244175_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,1)*(-((pow(ep.spa(3,4),2)*pow(ep.spb(2,0),2)*
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
   ep.spb(4,3))*ep.spb(5,4)))-
 complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q244178_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244187_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244281_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244283_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244286_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244288_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244293_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244298_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
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

template <class T> complex<T> A4q1_2q244316_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244318_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244328_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-(pow(ep.spa(0,4)*ep.spb(1,0)+
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

template <class T> complex<T> A4q1_2q244353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244358_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(
complex<T>(0,1)*(-((pow(ep.spa(3,4),2)*pow(ep.spb(2,0),2)*
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
   ep.spb(4,3))*ep.spb(5,4)))-
 complex<T>(0,1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q244365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244367_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244388_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*(pow(ep.spa(0,4)*ep.spb(2,0)+
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

template <class T> complex<T> A4q1_2q244400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q244517_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3)))-
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
(ep.s(3,4,5)*ep.spa(1,2)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)))+
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
  (ep.spa(1,4)*ep.spb(3,1)+
   ep.spa(2,4)*ep.spb(3,2)))/
 (ep.s(1,2,3)*ep.spa(4,5)*
  (-(ep.spa(2,4)*ep.spb(2,1))-
   ep.spa(3,4)*ep.spb(3,1))*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))-(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2)*(ep.spa(2,3)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(4,1)))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spb(1,0)*
 (-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1))*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))+
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
 (-(ep.spa(0,3)*ep.spb(5,3))-
  ep.spa(0,4)*ep.spb(5,4)))/
(ep.s(3,4,5)*ep.spa(0,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q244547_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244552_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q244577_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244587_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244589_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q244607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244612_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244619_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244624_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244713_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244718_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q1_2q244727_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3)))-
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2))/
(ep.s(3,4,5)*ep.spa(1,2)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)))+
 complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*pow(ep.spb(3,1),2)*
  (ep.spa(1,4)*ep.spb(3,1)+
   ep.spa(2,4)*ep.spb(3,2)))/
 (ep.s(1,2,3)*ep.spa(4,5)*
  (-(ep.spa(2,4)*ep.spb(2,1))-
   ep.spa(3,4)*ep.spb(3,1))*
  ep.spb(3,2)*(ep.spa(0,4)*
ep.spb(4,3)+ep.spa(0,5)*
ep.spb(5,3))))-(pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),2)*(ep.spa(2,3)*
   ep.spb(3,1)+ep.spa(2,4)*
   ep.spb(4,1)))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spb(1,0)*
 (-(ep.spa(2,4)*ep.spb(2,1))-
  ep.spa(3,4)*ep.spb(3,1))*
 (-(ep.spa(2,3)*ep.spb(5,3))-
  ep.spa(2,4)*ep.spb(5,4)))+
 (pow(ep.spa(0,2),2)*pow(ep.spb(5,3),2)*
 (-(ep.spa(0,3)*ep.spb(5,3))-
  ep.spa(0,4)*ep.spb(5,4)))/
(ep.s(3,4,5)*ep.spa(0,1)*
 (ep.spa(0,4)*ep.spb(4,3)+
  ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4)*(-(ep.spa(2,3)*
ep.spb(5,3))-ep.spa(2,4)*
   ep.spb(5,4))))
); }

template <class T> complex<T> A4q1_2q244732_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244748_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,-1)*((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q244762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244787_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q244792_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244799_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244804_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,1)*(-((pow(ep.spa(0,4),2)*
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

template <class T> complex<T> A4q1_2q244822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1_2q244834_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)

 ); }
#endif

template <class T> complex<T>  (*A4q1_2q2_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
#if _FAST_COMPILE_eval==1
case 281 : return &A4q1_2q2281_eval;
case 286 : return &A4q1_2q2286_eval;
case 353 : return &A4q1_2q2353_eval;
case 358 : return &A4q1_2q2358_eval;
case 371 : return &A4q1_2q2371_eval;
case 383 : return &A4q1_2q2383_eval;
case 391 : return &A4q1_2q2391_eval;
case 393 : return &A4q1_2q2393_eval;
case 406 : return &A4q1_2q2406_eval;
case 418 : return &A4q1_2q2418_eval;
case 421 : return &A4q1_2q2421_eval;
case 423 : return &A4q1_2q2423_eval;
case 713 : return &A4q1_2q2713_eval;
case 718 : return &A4q1_2q2718_eval;
case 785 : return &A4q1_2q2785_eval;
case 790 : return &A4q1_2q2790_eval;
case 803 : return &A4q1_2q2803_eval;
case 815 : return &A4q1_2q2815_eval;
case 823 : return &A4q1_2q2823_eval;
case 825 : return &A4q1_2q2825_eval;
case 838 : return &A4q1_2q2838_eval;
case 850 : return &A4q1_2q2850_eval;
case 853 : return &A4q1_2q2853_eval;
case 855 : return &A4q1_2q2855_eval;
case 911 : return &A4q1_2q2911_eval;
case 923 : return &A4q1_2q2923_eval;
case 931 : return &A4q1_2q2931_eval;
case 933 : return &A4q1_2q2933_eval;
case 983 : return &A4q1_2q2983_eval;
case 995 : return &A4q1_2q2995_eval;
case 1003 : return &A4q1_2q21003_eval;
case 1005 : return &A4q1_2q21005_eval;
case 1051 : return &A4q1_2q21051_eval;
case 1053 : return &A4q1_2q21053_eval;
case 1063 : return &A4q1_2q21063_eval;
case 1065 : return &A4q1_2q21065_eval;
case 1126 : return &A4q1_2q21126_eval;
case 1138 : return &A4q1_2q21138_eval;
case 1141 : return &A4q1_2q21141_eval;
case 1143 : return &A4q1_2q21143_eval;
case 1198 : return &A4q1_2q21198_eval;
case 1210 : return &A4q1_2q21210_eval;
case 1213 : return &A4q1_2q21213_eval;
case 1215 : return &A4q1_2q21215_eval;
case 1231 : return &A4q1_2q21231_eval;
case 1233 : return &A4q1_2q21233_eval;
case 1243 : return &A4q1_2q21243_eval;
case 1245 : return &A4q1_2q21245_eval;
case 1361 : return &A4q1_2q21361_eval;
case 1366 : return &A4q1_2q21366_eval;
case 1433 : return &A4q1_2q21433_eval;
case 1438 : return &A4q1_2q21438_eval;
case 1451 : return &A4q1_2q21451_eval;
case 1463 : return &A4q1_2q21463_eval;
case 1471 : return &A4q1_2q21471_eval;
case 1473 : return &A4q1_2q21473_eval;
case 1486 : return &A4q1_2q21486_eval;
case 1498 : return &A4q1_2q21498_eval;
case 1501 : return &A4q1_2q21501_eval;
case 1503 : return &A4q1_2q21503_eval;
case 1541 : return &A4q1_2q21541_eval;
case 1546 : return &A4q1_2q21546_eval;
case 1613 : return &A4q1_2q21613_eval;
case 1618 : return &A4q1_2q21618_eval;
case 1661 : return &A4q1_2q21661_eval;
case 1673 : return &A4q1_2q21673_eval;
case 1686 : return &A4q1_2q21686_eval;
case 1688 : return &A4q1_2q21688_eval;
case 1696 : return &A4q1_2q21696_eval;
case 1708 : return &A4q1_2q21708_eval;
case 1716 : return &A4q1_2q21716_eval;
case 1718 : return &A4q1_2q21718_eval;
case 1793 : return &A4q1_2q21793_eval;
case 1798 : return &A4q1_2q21798_eval;
case 1865 : return &A4q1_2q21865_eval;
case 1870 : return &A4q1_2q21870_eval;
case 1883 : return &A4q1_2q21883_eval;
case 1895 : return &A4q1_2q21895_eval;
case 1903 : return &A4q1_2q21903_eval;
case 1905 : return &A4q1_2q21905_eval;
case 1918 : return &A4q1_2q21918_eval;
case 1930 : return &A4q1_2q21930_eval;
case 1933 : return &A4q1_2q21933_eval;
case 1935 : return &A4q1_2q21935_eval;
case 1973 : return &A4q1_2q21973_eval;
case 1978 : return &A4q1_2q21978_eval;
case 2045 : return &A4q1_2q22045_eval;
case 2050 : return &A4q1_2q22050_eval;
case 2093 : return &A4q1_2q22093_eval;
case 2105 : return &A4q1_2q22105_eval;
case 2118 : return &A4q1_2q22118_eval;
case 2120 : return &A4q1_2q22120_eval;
case 2128 : return &A4q1_2q22128_eval;
case 2140 : return &A4q1_2q22140_eval;
case 2148 : return &A4q1_2q22148_eval;
case 2150 : return &A4q1_2q22150_eval;
case 2171 : return &A4q1_2q22171_eval;
case 2183 : return &A4q1_2q22183_eval;
case 2191 : return &A4q1_2q22191_eval;
case 2193 : return &A4q1_2q22193_eval;
case 2201 : return &A4q1_2q22201_eval;
case 2213 : return &A4q1_2q22213_eval;
case 2226 : return &A4q1_2q22226_eval;
case 2228 : return &A4q1_2q22228_eval;
case 2243 : return &A4q1_2q22243_eval;
case 2255 : return &A4q1_2q22255_eval;
case 2263 : return &A4q1_2q22263_eval;
case 2265 : return &A4q1_2q22265_eval;
case 2273 : return &A4q1_2q22273_eval;
case 2285 : return &A4q1_2q22285_eval;
case 2298 : return &A4q1_2q22298_eval;
case 2300 : return &A4q1_2q22300_eval;
case 2341 : return &A4q1_2q22341_eval;
case 2343 : return &A4q1_2q22343_eval;
case 2346 : return &A4q1_2q22346_eval;
case 2348 : return &A4q1_2q22348_eval;
case 2353 : return &A4q1_2q22353_eval;
case 2355 : return &A4q1_2q22355_eval;
case 2358 : return &A4q1_2q22358_eval;
case 2360 : return &A4q1_2q22360_eval;
case 2386 : return &A4q1_2q22386_eval;
case 2398 : return &A4q1_2q22398_eval;
case 2401 : return &A4q1_2q22401_eval;
case 2403 : return &A4q1_2q22403_eval;
case 2416 : return &A4q1_2q22416_eval;
case 2428 : return &A4q1_2q22428_eval;
case 2436 : return &A4q1_2q22436_eval;
case 2438 : return &A4q1_2q22438_eval;
case 2458 : return &A4q1_2q22458_eval;
case 2470 : return &A4q1_2q22470_eval;
case 2473 : return &A4q1_2q22473_eval;
case 2475 : return &A4q1_2q22475_eval;
case 2488 : return &A4q1_2q22488_eval;
case 2500 : return &A4q1_2q22500_eval;
case 2508 : return &A4q1_2q22508_eval;
case 2510 : return &A4q1_2q22510_eval;
case 2521 : return &A4q1_2q22521_eval;
case 2523 : return &A4q1_2q22523_eval;
case 2526 : return &A4q1_2q22526_eval;
case 2528 : return &A4q1_2q22528_eval;
case 2533 : return &A4q1_2q22533_eval;
case 2535 : return &A4q1_2q22535_eval;
case 2538 : return &A4q1_2q22538_eval;
case 2540 : return &A4q1_2q22540_eval;
case 2873 : return &A4q1_2q22873_eval;
case 2878 : return &A4q1_2q22878_eval;
case 2945 : return &A4q1_2q22945_eval;
case 2950 : return &A4q1_2q22950_eval;
case 2963 : return &A4q1_2q22963_eval;
case 2975 : return &A4q1_2q22975_eval;
case 2983 : return &A4q1_2q22983_eval;
case 2985 : return &A4q1_2q22985_eval;
case 2998 : return &A4q1_2q22998_eval;
case 3010 : return &A4q1_2q23010_eval;
case 3013 : return &A4q1_2q23013_eval;
case 3015 : return &A4q1_2q23015_eval;
case 3305 : return &A4q1_2q23305_eval;
case 3310 : return &A4q1_2q23310_eval;
case 3377 : return &A4q1_2q23377_eval;
case 3382 : return &A4q1_2q23382_eval;
case 3395 : return &A4q1_2q23395_eval;
case 3407 : return &A4q1_2q23407_eval;
case 3415 : return &A4q1_2q23415_eval;
case 3417 : return &A4q1_2q23417_eval;
case 3430 : return &A4q1_2q23430_eval;
case 3442 : return &A4q1_2q23442_eval;
case 3445 : return &A4q1_2q23445_eval;
case 3447 : return &A4q1_2q23447_eval;
case 3503 : return &A4q1_2q23503_eval;
case 3515 : return &A4q1_2q23515_eval;
case 3523 : return &A4q1_2q23523_eval;
case 3525 : return &A4q1_2q23525_eval;
case 3575 : return &A4q1_2q23575_eval;
case 3587 : return &A4q1_2q23587_eval;
case 3595 : return &A4q1_2q23595_eval;
case 3597 : return &A4q1_2q23597_eval;
case 3643 : return &A4q1_2q23643_eval;
case 3645 : return &A4q1_2q23645_eval;
case 3655 : return &A4q1_2q23655_eval;
case 3657 : return &A4q1_2q23657_eval;
case 3718 : return &A4q1_2q23718_eval;
case 3730 : return &A4q1_2q23730_eval;
case 3733 : return &A4q1_2q23733_eval;
case 3735 : return &A4q1_2q23735_eval;
case 3790 : return &A4q1_2q23790_eval;
case 3802 : return &A4q1_2q23802_eval;
case 3805 : return &A4q1_2q23805_eval;
case 3807 : return &A4q1_2q23807_eval;
case 3823 : return &A4q1_2q23823_eval;
case 3825 : return &A4q1_2q23825_eval;
case 3835 : return &A4q1_2q23835_eval;
case 3837 : return &A4q1_2q23837_eval;
case 3953 : return &A4q1_2q23953_eval;
case 3958 : return &A4q1_2q23958_eval;
case 4025 : return &A4q1_2q24025_eval;
case 4030 : return &A4q1_2q24030_eval;
case 4043 : return &A4q1_2q24043_eval;
case 4055 : return &A4q1_2q24055_eval;
case 4063 : return &A4q1_2q24063_eval;
case 4065 : return &A4q1_2q24065_eval;
case 4078 : return &A4q1_2q24078_eval;
case 4090 : return &A4q1_2q24090_eval;
case 4093 : return &A4q1_2q24093_eval;
case 4095 : return &A4q1_2q24095_eval;
case 4133 : return &A4q1_2q24133_eval;
case 4138 : return &A4q1_2q24138_eval;
case 4205 : return &A4q1_2q24205_eval;
case 4210 : return &A4q1_2q24210_eval;
case 4253 : return &A4q1_2q24253_eval;
case 4265 : return &A4q1_2q24265_eval;
case 4278 : return &A4q1_2q24278_eval;
case 4280 : return &A4q1_2q24280_eval;
case 4288 : return &A4q1_2q24288_eval;
case 4300 : return &A4q1_2q24300_eval;
case 4308 : return &A4q1_2q24308_eval;
case 4310 : return &A4q1_2q24310_eval;
case 4385 : return &A4q1_2q24385_eval;
case 4390 : return &A4q1_2q24390_eval;
case 4457 : return &A4q1_2q24457_eval;
case 4462 : return &A4q1_2q24462_eval;
case 4475 : return &A4q1_2q24475_eval;
case 4487 : return &A4q1_2q24487_eval;
case 4495 : return &A4q1_2q24495_eval;
case 4497 : return &A4q1_2q24497_eval;
case 4510 : return &A4q1_2q24510_eval;
case 4522 : return &A4q1_2q24522_eval;
case 4525 : return &A4q1_2q24525_eval;
case 4527 : return &A4q1_2q24527_eval;
case 4565 : return &A4q1_2q24565_eval;
case 4570 : return &A4q1_2q24570_eval;
case 4637 : return &A4q1_2q24637_eval;
case 4642 : return &A4q1_2q24642_eval;
case 4685 : return &A4q1_2q24685_eval;
case 4697 : return &A4q1_2q24697_eval;
case 4710 : return &A4q1_2q24710_eval;
case 4712 : return &A4q1_2q24712_eval;
case 4720 : return &A4q1_2q24720_eval;
case 4732 : return &A4q1_2q24732_eval;
case 4740 : return &A4q1_2q24740_eval;
case 4742 : return &A4q1_2q24742_eval;
case 4763 : return &A4q1_2q24763_eval;
case 4775 : return &A4q1_2q24775_eval;
case 4783 : return &A4q1_2q24783_eval;
case 4785 : return &A4q1_2q24785_eval;
case 4793 : return &A4q1_2q24793_eval;
case 4805 : return &A4q1_2q24805_eval;
case 4818 : return &A4q1_2q24818_eval;
case 4820 : return &A4q1_2q24820_eval;
case 4835 : return &A4q1_2q24835_eval;
case 4847 : return &A4q1_2q24847_eval;
case 4855 : return &A4q1_2q24855_eval;
case 4857 : return &A4q1_2q24857_eval;
case 4865 : return &A4q1_2q24865_eval;
case 4877 : return &A4q1_2q24877_eval;
case 4890 : return &A4q1_2q24890_eval;
case 4892 : return &A4q1_2q24892_eval;
case 4933 : return &A4q1_2q24933_eval;
case 4935 : return &A4q1_2q24935_eval;
case 4938 : return &A4q1_2q24938_eval;
case 4940 : return &A4q1_2q24940_eval;
case 4945 : return &A4q1_2q24945_eval;
case 4947 : return &A4q1_2q24947_eval;
case 4950 : return &A4q1_2q24950_eval;
case 4952 : return &A4q1_2q24952_eval;
case 4978 : return &A4q1_2q24978_eval;
case 4990 : return &A4q1_2q24990_eval;
case 4993 : return &A4q1_2q24993_eval;
case 4995 : return &A4q1_2q24995_eval;
case 5008 : return &A4q1_2q25008_eval;
case 5020 : return &A4q1_2q25020_eval;
case 5028 : return &A4q1_2q25028_eval;
case 5030 : return &A4q1_2q25030_eval;
case 5050 : return &A4q1_2q25050_eval;
case 5062 : return &A4q1_2q25062_eval;
case 5065 : return &A4q1_2q25065_eval;
case 5067 : return &A4q1_2q25067_eval;
case 5080 : return &A4q1_2q25080_eval;
case 5092 : return &A4q1_2q25092_eval;
case 5100 : return &A4q1_2q25100_eval;
case 5102 : return &A4q1_2q25102_eval;
case 5113 : return &A4q1_2q25113_eval;
case 5115 : return &A4q1_2q25115_eval;
case 5118 : return &A4q1_2q25118_eval;
case 5120 : return &A4q1_2q25120_eval;
case 5125 : return &A4q1_2q25125_eval;
case 5127 : return &A4q1_2q25127_eval;
case 5130 : return &A4q1_2q25130_eval;
case 5132 : return &A4q1_2q25132_eval;
case 5231 : return &A4q1_2q25231_eval;
case 5243 : return &A4q1_2q25243_eval;
case 5251 : return &A4q1_2q25251_eval;
case 5253 : return &A4q1_2q25253_eval;
case 5303 : return &A4q1_2q25303_eval;
case 5315 : return &A4q1_2q25315_eval;
case 5323 : return &A4q1_2q25323_eval;
case 5325 : return &A4q1_2q25325_eval;
case 5371 : return &A4q1_2q25371_eval;
case 5373 : return &A4q1_2q25373_eval;
case 5383 : return &A4q1_2q25383_eval;
case 5385 : return &A4q1_2q25385_eval;
case 5411 : return &A4q1_2q25411_eval;
case 5423 : return &A4q1_2q25423_eval;
case 5431 : return &A4q1_2q25431_eval;
case 5433 : return &A4q1_2q25433_eval;
case 5441 : return &A4q1_2q25441_eval;
case 5453 : return &A4q1_2q25453_eval;
case 5466 : return &A4q1_2q25466_eval;
case 5468 : return &A4q1_2q25468_eval;
case 5483 : return &A4q1_2q25483_eval;
case 5495 : return &A4q1_2q25495_eval;
case 5503 : return &A4q1_2q25503_eval;
case 5505 : return &A4q1_2q25505_eval;
case 5513 : return &A4q1_2q25513_eval;
case 5525 : return &A4q1_2q25525_eval;
case 5538 : return &A4q1_2q25538_eval;
case 5540 : return &A4q1_2q25540_eval;
case 5581 : return &A4q1_2q25581_eval;
case 5583 : return &A4q1_2q25583_eval;
case 5586 : return &A4q1_2q25586_eval;
case 5588 : return &A4q1_2q25588_eval;
case 5593 : return &A4q1_2q25593_eval;
case 5595 : return &A4q1_2q25595_eval;
case 5598 : return &A4q1_2q25598_eval;
case 5600 : return &A4q1_2q25600_eval;
case 5663 : return &A4q1_2q25663_eval;
case 5675 : return &A4q1_2q25675_eval;
case 5683 : return &A4q1_2q25683_eval;
case 5685 : return &A4q1_2q25685_eval;
case 5735 : return &A4q1_2q25735_eval;
case 5747 : return &A4q1_2q25747_eval;
case 5755 : return &A4q1_2q25755_eval;
case 5757 : return &A4q1_2q25757_eval;
case 5803 : return &A4q1_2q25803_eval;
case 5805 : return &A4q1_2q25805_eval;
case 5815 : return &A4q1_2q25815_eval;
case 5817 : return &A4q1_2q25817_eval;
case 5843 : return &A4q1_2q25843_eval;
case 5855 : return &A4q1_2q25855_eval;
case 5863 : return &A4q1_2q25863_eval;
case 5865 : return &A4q1_2q25865_eval;
case 5873 : return &A4q1_2q25873_eval;
case 5885 : return &A4q1_2q25885_eval;
case 5898 : return &A4q1_2q25898_eval;
case 5900 : return &A4q1_2q25900_eval;
case 5915 : return &A4q1_2q25915_eval;
case 5927 : return &A4q1_2q25927_eval;
case 5935 : return &A4q1_2q25935_eval;
case 5937 : return &A4q1_2q25937_eval;
case 5945 : return &A4q1_2q25945_eval;
case 5957 : return &A4q1_2q25957_eval;
case 5970 : return &A4q1_2q25970_eval;
case 5972 : return &A4q1_2q25972_eval;
case 6013 : return &A4q1_2q26013_eval;
case 6015 : return &A4q1_2q26015_eval;
case 6018 : return &A4q1_2q26018_eval;
case 6020 : return &A4q1_2q26020_eval;
case 6025 : return &A4q1_2q26025_eval;
case 6027 : return &A4q1_2q26027_eval;
case 6030 : return &A4q1_2q26030_eval;
case 6032 : return &A4q1_2q26032_eval;
case 6271 : return &A4q1_2q26271_eval;
case 6273 : return &A4q1_2q26273_eval;
case 6283 : return &A4q1_2q26283_eval;
case 6285 : return &A4q1_2q26285_eval;
case 6301 : return &A4q1_2q26301_eval;
case 6303 : return &A4q1_2q26303_eval;
case 6306 : return &A4q1_2q26306_eval;
case 6308 : return &A4q1_2q26308_eval;
case 6313 : return &A4q1_2q26313_eval;
case 6315 : return &A4q1_2q26315_eval;
case 6318 : return &A4q1_2q26318_eval;
case 6320 : return &A4q1_2q26320_eval;
case 6343 : return &A4q1_2q26343_eval;
case 6345 : return &A4q1_2q26345_eval;
case 6355 : return &A4q1_2q26355_eval;
case 6357 : return &A4q1_2q26357_eval;
case 6373 : return &A4q1_2q26373_eval;
case 6375 : return &A4q1_2q26375_eval;
case 6378 : return &A4q1_2q26378_eval;
case 6380 : return &A4q1_2q26380_eval;
case 6385 : return &A4q1_2q26385_eval;
case 6387 : return &A4q1_2q26387_eval;
case 6390 : return &A4q1_2q26390_eval;
case 6392 : return &A4q1_2q26392_eval;
case 6526 : return &A4q1_2q26526_eval;
case 6538 : return &A4q1_2q26538_eval;
case 6541 : return &A4q1_2q26541_eval;
case 6543 : return &A4q1_2q26543_eval;
case 6598 : return &A4q1_2q26598_eval;
case 6610 : return &A4q1_2q26610_eval;
case 6613 : return &A4q1_2q26613_eval;
case 6615 : return &A4q1_2q26615_eval;
case 6631 : return &A4q1_2q26631_eval;
case 6633 : return &A4q1_2q26633_eval;
case 6643 : return &A4q1_2q26643_eval;
case 6645 : return &A4q1_2q26645_eval;
case 6706 : return &A4q1_2q26706_eval;
case 6718 : return &A4q1_2q26718_eval;
case 6721 : return &A4q1_2q26721_eval;
case 6723 : return &A4q1_2q26723_eval;
case 6736 : return &A4q1_2q26736_eval;
case 6748 : return &A4q1_2q26748_eval;
case 6756 : return &A4q1_2q26756_eval;
case 6758 : return &A4q1_2q26758_eval;
case 6778 : return &A4q1_2q26778_eval;
case 6790 : return &A4q1_2q26790_eval;
case 6793 : return &A4q1_2q26793_eval;
case 6795 : return &A4q1_2q26795_eval;
case 6808 : return &A4q1_2q26808_eval;
case 6820 : return &A4q1_2q26820_eval;
case 6828 : return &A4q1_2q26828_eval;
case 6830 : return &A4q1_2q26830_eval;
case 6841 : return &A4q1_2q26841_eval;
case 6843 : return &A4q1_2q26843_eval;
case 6846 : return &A4q1_2q26846_eval;
case 6848 : return &A4q1_2q26848_eval;
case 6853 : return &A4q1_2q26853_eval;
case 6855 : return &A4q1_2q26855_eval;
case 6858 : return &A4q1_2q26858_eval;
case 6860 : return &A4q1_2q26860_eval;
case 6958 : return &A4q1_2q26958_eval;
case 6970 : return &A4q1_2q26970_eval;
case 6973 : return &A4q1_2q26973_eval;
case 6975 : return &A4q1_2q26975_eval;
case 7030 : return &A4q1_2q27030_eval;
case 7042 : return &A4q1_2q27042_eval;
case 7045 : return &A4q1_2q27045_eval;
case 7047 : return &A4q1_2q27047_eval;
case 7063 : return &A4q1_2q27063_eval;
case 7065 : return &A4q1_2q27065_eval;
case 7075 : return &A4q1_2q27075_eval;
case 7077 : return &A4q1_2q27077_eval;
case 7138 : return &A4q1_2q27138_eval;
case 7150 : return &A4q1_2q27150_eval;
case 7153 : return &A4q1_2q27153_eval;
case 7155 : return &A4q1_2q27155_eval;
case 7168 : return &A4q1_2q27168_eval;
case 7180 : return &A4q1_2q27180_eval;
case 7188 : return &A4q1_2q27188_eval;
case 7190 : return &A4q1_2q27190_eval;
case 7210 : return &A4q1_2q27210_eval;
case 7222 : return &A4q1_2q27222_eval;
case 7225 : return &A4q1_2q27225_eval;
case 7227 : return &A4q1_2q27227_eval;
case 7240 : return &A4q1_2q27240_eval;
case 7252 : return &A4q1_2q27252_eval;
case 7260 : return &A4q1_2q27260_eval;
case 7262 : return &A4q1_2q27262_eval;
case 7273 : return &A4q1_2q27273_eval;
case 7275 : return &A4q1_2q27275_eval;
case 7278 : return &A4q1_2q27278_eval;
case 7280 : return &A4q1_2q27280_eval;
case 7285 : return &A4q1_2q27285_eval;
case 7287 : return &A4q1_2q27287_eval;
case 7290 : return &A4q1_2q27290_eval;
case 7292 : return &A4q1_2q27292_eval;
case 7351 : return &A4q1_2q27351_eval;
case 7353 : return &A4q1_2q27353_eval;
case 7363 : return &A4q1_2q27363_eval;
case 7365 : return &A4q1_2q27365_eval;
case 7381 : return &A4q1_2q27381_eval;
case 7383 : return &A4q1_2q27383_eval;
case 7386 : return &A4q1_2q27386_eval;
case 7388 : return &A4q1_2q27388_eval;
case 7393 : return &A4q1_2q27393_eval;
case 7395 : return &A4q1_2q27395_eval;
case 7398 : return &A4q1_2q27398_eval;
case 7400 : return &A4q1_2q27400_eval;
case 7423 : return &A4q1_2q27423_eval;
case 7425 : return &A4q1_2q27425_eval;
case 7435 : return &A4q1_2q27435_eval;
case 7437 : return &A4q1_2q27437_eval;
case 7453 : return &A4q1_2q27453_eval;
case 7455 : return &A4q1_2q27455_eval;
case 7458 : return &A4q1_2q27458_eval;
case 7460 : return &A4q1_2q27460_eval;
case 7465 : return &A4q1_2q27465_eval;
case 7467 : return &A4q1_2q27467_eval;
case 7470 : return &A4q1_2q27470_eval;
case 7472 : return &A4q1_2q27472_eval;
case 7841 : return &A4q1_2q27841_eval;
case 7846 : return &A4q1_2q27846_eval;
case 7913 : return &A4q1_2q27913_eval;
case 7918 : return &A4q1_2q27918_eval;
case 7931 : return &A4q1_2q27931_eval;
case 7943 : return &A4q1_2q27943_eval;
case 7951 : return &A4q1_2q27951_eval;
case 7953 : return &A4q1_2q27953_eval;
case 7966 : return &A4q1_2q27966_eval;
case 7978 : return &A4q1_2q27978_eval;
case 7981 : return &A4q1_2q27981_eval;
case 7983 : return &A4q1_2q27983_eval;
case 8021 : return &A4q1_2q28021_eval;
case 8026 : return &A4q1_2q28026_eval;
case 8093 : return &A4q1_2q28093_eval;
case 8098 : return &A4q1_2q28098_eval;
case 8141 : return &A4q1_2q28141_eval;
case 8153 : return &A4q1_2q28153_eval;
case 8166 : return &A4q1_2q28166_eval;
case 8168 : return &A4q1_2q28168_eval;
case 8176 : return &A4q1_2q28176_eval;
case 8188 : return &A4q1_2q28188_eval;
case 8196 : return &A4q1_2q28196_eval;
case 8198 : return &A4q1_2q28198_eval;
case 8273 : return &A4q1_2q28273_eval;
case 8278 : return &A4q1_2q28278_eval;
case 8345 : return &A4q1_2q28345_eval;
case 8350 : return &A4q1_2q28350_eval;
case 8363 : return &A4q1_2q28363_eval;
case 8375 : return &A4q1_2q28375_eval;
case 8383 : return &A4q1_2q28383_eval;
case 8385 : return &A4q1_2q28385_eval;
case 8398 : return &A4q1_2q28398_eval;
case 8410 : return &A4q1_2q28410_eval;
case 8413 : return &A4q1_2q28413_eval;
case 8415 : return &A4q1_2q28415_eval;
case 8453 : return &A4q1_2q28453_eval;
case 8458 : return &A4q1_2q28458_eval;
case 8525 : return &A4q1_2q28525_eval;
case 8530 : return &A4q1_2q28530_eval;
case 8573 : return &A4q1_2q28573_eval;
case 8585 : return &A4q1_2q28585_eval;
case 8598 : return &A4q1_2q28598_eval;
case 8600 : return &A4q1_2q28600_eval;
case 8608 : return &A4q1_2q28608_eval;
case 8620 : return &A4q1_2q28620_eval;
case 8628 : return &A4q1_2q28628_eval;
case 8630 : return &A4q1_2q28630_eval;
case 8651 : return &A4q1_2q28651_eval;
case 8663 : return &A4q1_2q28663_eval;
case 8671 : return &A4q1_2q28671_eval;
case 8673 : return &A4q1_2q28673_eval;
case 8681 : return &A4q1_2q28681_eval;
case 8693 : return &A4q1_2q28693_eval;
case 8706 : return &A4q1_2q28706_eval;
case 8708 : return &A4q1_2q28708_eval;
case 8723 : return &A4q1_2q28723_eval;
case 8735 : return &A4q1_2q28735_eval;
case 8743 : return &A4q1_2q28743_eval;
case 8745 : return &A4q1_2q28745_eval;
case 8753 : return &A4q1_2q28753_eval;
case 8765 : return &A4q1_2q28765_eval;
case 8778 : return &A4q1_2q28778_eval;
case 8780 : return &A4q1_2q28780_eval;
case 8821 : return &A4q1_2q28821_eval;
case 8823 : return &A4q1_2q28823_eval;
case 8826 : return &A4q1_2q28826_eval;
case 8828 : return &A4q1_2q28828_eval;
case 8833 : return &A4q1_2q28833_eval;
case 8835 : return &A4q1_2q28835_eval;
case 8838 : return &A4q1_2q28838_eval;
case 8840 : return &A4q1_2q28840_eval;
case 8866 : return &A4q1_2q28866_eval;
case 8878 : return &A4q1_2q28878_eval;
case 8881 : return &A4q1_2q28881_eval;
case 8883 : return &A4q1_2q28883_eval;
case 8896 : return &A4q1_2q28896_eval;
case 8908 : return &A4q1_2q28908_eval;
case 8916 : return &A4q1_2q28916_eval;
case 8918 : return &A4q1_2q28918_eval;
case 8938 : return &A4q1_2q28938_eval;
case 8950 : return &A4q1_2q28950_eval;
case 8953 : return &A4q1_2q28953_eval;
case 8955 : return &A4q1_2q28955_eval;
case 8968 : return &A4q1_2q28968_eval;
case 8980 : return &A4q1_2q28980_eval;
case 8988 : return &A4q1_2q28988_eval;
case 8990 : return &A4q1_2q28990_eval;
case 9001 : return &A4q1_2q29001_eval;
case 9003 : return &A4q1_2q29003_eval;
case 9006 : return &A4q1_2q29006_eval;
case 9008 : return &A4q1_2q29008_eval;
case 9013 : return &A4q1_2q29013_eval;
case 9015 : return &A4q1_2q29015_eval;
case 9018 : return &A4q1_2q29018_eval;
case 9020 : return &A4q1_2q29020_eval;
case 9101 : return &A4q1_2q29101_eval;
case 9106 : return &A4q1_2q29106_eval;
case 9173 : return &A4q1_2q29173_eval;
case 9178 : return &A4q1_2q29178_eval;
case 9221 : return &A4q1_2q29221_eval;
case 9233 : return &A4q1_2q29233_eval;
case 9246 : return &A4q1_2q29246_eval;
case 9248 : return &A4q1_2q29248_eval;
case 9256 : return &A4q1_2q29256_eval;
case 9268 : return &A4q1_2q29268_eval;
case 9276 : return &A4q1_2q29276_eval;
case 9278 : return &A4q1_2q29278_eval;
case 9533 : return &A4q1_2q29533_eval;
case 9538 : return &A4q1_2q29538_eval;
case 9605 : return &A4q1_2q29605_eval;
case 9610 : return &A4q1_2q29610_eval;
case 9653 : return &A4q1_2q29653_eval;
case 9665 : return &A4q1_2q29665_eval;
case 9678 : return &A4q1_2q29678_eval;
case 9680 : return &A4q1_2q29680_eval;
case 9688 : return &A4q1_2q29688_eval;
case 9700 : return &A4q1_2q29700_eval;
case 9708 : return &A4q1_2q29708_eval;
case 9710 : return &A4q1_2q29710_eval;
case 9941 : return &A4q1_2q29941_eval;
case 9953 : return &A4q1_2q29953_eval;
case 9966 : return &A4q1_2q29966_eval;
case 9968 : return &A4q1_2q29968_eval;
case 10013 : return &A4q1_2q210013_eval;
case 10025 : return &A4q1_2q210025_eval;
case 10038 : return &A4q1_2q210038_eval;
case 10040 : return &A4q1_2q210040_eval;
case 10116 : return &A4q1_2q210116_eval;
case 10118 : return &A4q1_2q210118_eval;
case 10128 : return &A4q1_2q210128_eval;
case 10130 : return &A4q1_2q210130_eval;
case 10156 : return &A4q1_2q210156_eval;
case 10168 : return &A4q1_2q210168_eval;
case 10176 : return &A4q1_2q210176_eval;
case 10178 : return &A4q1_2q210178_eval;
case 10228 : return &A4q1_2q210228_eval;
case 10240 : return &A4q1_2q210240_eval;
case 10248 : return &A4q1_2q210248_eval;
case 10250 : return &A4q1_2q210250_eval;
case 10296 : return &A4q1_2q210296_eval;
case 10298 : return &A4q1_2q210298_eval;
case 10308 : return &A4q1_2q210308_eval;
case 10310 : return &A4q1_2q210310_eval;
case 10433 : return &A4q1_2q210433_eval;
case 10438 : return &A4q1_2q210438_eval;
case 10505 : return &A4q1_2q210505_eval;
case 10510 : return &A4q1_2q210510_eval;
case 10523 : return &A4q1_2q210523_eval;
case 10535 : return &A4q1_2q210535_eval;
case 10543 : return &A4q1_2q210543_eval;
case 10545 : return &A4q1_2q210545_eval;
case 10558 : return &A4q1_2q210558_eval;
case 10570 : return &A4q1_2q210570_eval;
case 10573 : return &A4q1_2q210573_eval;
case 10575 : return &A4q1_2q210575_eval;
case 10613 : return &A4q1_2q210613_eval;
case 10618 : return &A4q1_2q210618_eval;
case 10685 : return &A4q1_2q210685_eval;
case 10690 : return &A4q1_2q210690_eval;
case 10733 : return &A4q1_2q210733_eval;
case 10745 : return &A4q1_2q210745_eval;
case 10758 : return &A4q1_2q210758_eval;
case 10760 : return &A4q1_2q210760_eval;
case 10768 : return &A4q1_2q210768_eval;
case 10780 : return &A4q1_2q210780_eval;
case 10788 : return &A4q1_2q210788_eval;
case 10790 : return &A4q1_2q210790_eval;
case 10865 : return &A4q1_2q210865_eval;
case 10870 : return &A4q1_2q210870_eval;
case 10937 : return &A4q1_2q210937_eval;
case 10942 : return &A4q1_2q210942_eval;
case 10955 : return &A4q1_2q210955_eval;
case 10967 : return &A4q1_2q210967_eval;
case 10975 : return &A4q1_2q210975_eval;
case 10977 : return &A4q1_2q210977_eval;
case 10990 : return &A4q1_2q210990_eval;
case 11002 : return &A4q1_2q211002_eval;
case 11005 : return &A4q1_2q211005_eval;
case 11007 : return &A4q1_2q211007_eval;
case 11045 : return &A4q1_2q211045_eval;
case 11050 : return &A4q1_2q211050_eval;
case 11117 : return &A4q1_2q211117_eval;
case 11122 : return &A4q1_2q211122_eval;
case 11165 : return &A4q1_2q211165_eval;
case 11177 : return &A4q1_2q211177_eval;
case 11190 : return &A4q1_2q211190_eval;
case 11192 : return &A4q1_2q211192_eval;
case 11200 : return &A4q1_2q211200_eval;
case 11212 : return &A4q1_2q211212_eval;
case 11220 : return &A4q1_2q211220_eval;
case 11222 : return &A4q1_2q211222_eval;
case 11243 : return &A4q1_2q211243_eval;
case 11255 : return &A4q1_2q211255_eval;
case 11263 : return &A4q1_2q211263_eval;
case 11265 : return &A4q1_2q211265_eval;
case 11273 : return &A4q1_2q211273_eval;
case 11285 : return &A4q1_2q211285_eval;
case 11298 : return &A4q1_2q211298_eval;
case 11300 : return &A4q1_2q211300_eval;
case 11315 : return &A4q1_2q211315_eval;
case 11327 : return &A4q1_2q211327_eval;
case 11335 : return &A4q1_2q211335_eval;
case 11337 : return &A4q1_2q211337_eval;
case 11345 : return &A4q1_2q211345_eval;
case 11357 : return &A4q1_2q211357_eval;
case 11370 : return &A4q1_2q211370_eval;
case 11372 : return &A4q1_2q211372_eval;
case 11413 : return &A4q1_2q211413_eval;
case 11415 : return &A4q1_2q211415_eval;
case 11418 : return &A4q1_2q211418_eval;
case 11420 : return &A4q1_2q211420_eval;
case 11425 : return &A4q1_2q211425_eval;
case 11427 : return &A4q1_2q211427_eval;
case 11430 : return &A4q1_2q211430_eval;
case 11432 : return &A4q1_2q211432_eval;
case 11458 : return &A4q1_2q211458_eval;
case 11470 : return &A4q1_2q211470_eval;
case 11473 : return &A4q1_2q211473_eval;
case 11475 : return &A4q1_2q211475_eval;
case 11488 : return &A4q1_2q211488_eval;
case 11500 : return &A4q1_2q211500_eval;
case 11508 : return &A4q1_2q211508_eval;
case 11510 : return &A4q1_2q211510_eval;
case 11530 : return &A4q1_2q211530_eval;
case 11542 : return &A4q1_2q211542_eval;
case 11545 : return &A4q1_2q211545_eval;
case 11547 : return &A4q1_2q211547_eval;
case 11560 : return &A4q1_2q211560_eval;
case 11572 : return &A4q1_2q211572_eval;
case 11580 : return &A4q1_2q211580_eval;
case 11582 : return &A4q1_2q211582_eval;
case 11593 : return &A4q1_2q211593_eval;
case 11595 : return &A4q1_2q211595_eval;
case 11598 : return &A4q1_2q211598_eval;
case 11600 : return &A4q1_2q211600_eval;
case 11605 : return &A4q1_2q211605_eval;
case 11607 : return &A4q1_2q211607_eval;
case 11610 : return &A4q1_2q211610_eval;
case 11612 : return &A4q1_2q211612_eval;
case 11693 : return &A4q1_2q211693_eval;
case 11698 : return &A4q1_2q211698_eval;
case 11765 : return &A4q1_2q211765_eval;
case 11770 : return &A4q1_2q211770_eval;
case 11813 : return &A4q1_2q211813_eval;
case 11825 : return &A4q1_2q211825_eval;
case 11838 : return &A4q1_2q211838_eval;
case 11840 : return &A4q1_2q211840_eval;
case 11848 : return &A4q1_2q211848_eval;
case 11860 : return &A4q1_2q211860_eval;
case 11868 : return &A4q1_2q211868_eval;
case 11870 : return &A4q1_2q211870_eval;
case 12125 : return &A4q1_2q212125_eval;
case 12130 : return &A4q1_2q212130_eval;
case 12197 : return &A4q1_2q212197_eval;
case 12202 : return &A4q1_2q212202_eval;
case 12245 : return &A4q1_2q212245_eval;
case 12257 : return &A4q1_2q212257_eval;
case 12270 : return &A4q1_2q212270_eval;
case 12272 : return &A4q1_2q212272_eval;
case 12280 : return &A4q1_2q212280_eval;
case 12292 : return &A4q1_2q212292_eval;
case 12300 : return &A4q1_2q212300_eval;
case 12302 : return &A4q1_2q212302_eval;
case 12533 : return &A4q1_2q212533_eval;
case 12545 : return &A4q1_2q212545_eval;
case 12558 : return &A4q1_2q212558_eval;
case 12560 : return &A4q1_2q212560_eval;
case 12605 : return &A4q1_2q212605_eval;
case 12617 : return &A4q1_2q212617_eval;
case 12630 : return &A4q1_2q212630_eval;
case 12632 : return &A4q1_2q212632_eval;
case 12708 : return &A4q1_2q212708_eval;
case 12710 : return &A4q1_2q212710_eval;
case 12720 : return &A4q1_2q212720_eval;
case 12722 : return &A4q1_2q212722_eval;
case 12748 : return &A4q1_2q212748_eval;
case 12760 : return &A4q1_2q212760_eval;
case 12768 : return &A4q1_2q212768_eval;
case 12770 : return &A4q1_2q212770_eval;
case 12820 : return &A4q1_2q212820_eval;
case 12832 : return &A4q1_2q212832_eval;
case 12840 : return &A4q1_2q212840_eval;
case 12842 : return &A4q1_2q212842_eval;
case 12888 : return &A4q1_2q212888_eval;
case 12890 : return &A4q1_2q212890_eval;
case 12900 : return &A4q1_2q212900_eval;
case 12902 : return &A4q1_2q212902_eval;
case 12971 : return &A4q1_2q212971_eval;
case 12983 : return &A4q1_2q212983_eval;
case 12991 : return &A4q1_2q212991_eval;
case 12993 : return &A4q1_2q212993_eval;
case 13001 : return &A4q1_2q213001_eval;
case 13013 : return &A4q1_2q213013_eval;
case 13026 : return &A4q1_2q213026_eval;
case 13028 : return &A4q1_2q213028_eval;
case 13043 : return &A4q1_2q213043_eval;
case 13055 : return &A4q1_2q213055_eval;
case 13063 : return &A4q1_2q213063_eval;
case 13065 : return &A4q1_2q213065_eval;
case 13073 : return &A4q1_2q213073_eval;
case 13085 : return &A4q1_2q213085_eval;
case 13098 : return &A4q1_2q213098_eval;
case 13100 : return &A4q1_2q213100_eval;
case 13141 : return &A4q1_2q213141_eval;
case 13143 : return &A4q1_2q213143_eval;
case 13146 : return &A4q1_2q213146_eval;
case 13148 : return &A4q1_2q213148_eval;
case 13153 : return &A4q1_2q213153_eval;
case 13155 : return &A4q1_2q213155_eval;
case 13158 : return &A4q1_2q213158_eval;
case 13160 : return &A4q1_2q213160_eval;
case 13181 : return &A4q1_2q213181_eval;
case 13193 : return &A4q1_2q213193_eval;
case 13206 : return &A4q1_2q213206_eval;
case 13208 : return &A4q1_2q213208_eval;
case 13253 : return &A4q1_2q213253_eval;
case 13265 : return &A4q1_2q213265_eval;
case 13278 : return &A4q1_2q213278_eval;
case 13280 : return &A4q1_2q213280_eval;
case 13356 : return &A4q1_2q213356_eval;
case 13358 : return &A4q1_2q213358_eval;
case 13368 : return &A4q1_2q213368_eval;
case 13370 : return &A4q1_2q213370_eval;
case 13403 : return &A4q1_2q213403_eval;
case 13415 : return &A4q1_2q213415_eval;
case 13423 : return &A4q1_2q213423_eval;
case 13425 : return &A4q1_2q213425_eval;
case 13433 : return &A4q1_2q213433_eval;
case 13445 : return &A4q1_2q213445_eval;
case 13458 : return &A4q1_2q213458_eval;
case 13460 : return &A4q1_2q213460_eval;
case 13475 : return &A4q1_2q213475_eval;
case 13487 : return &A4q1_2q213487_eval;
case 13495 : return &A4q1_2q213495_eval;
case 13497 : return &A4q1_2q213497_eval;
case 13505 : return &A4q1_2q213505_eval;
case 13517 : return &A4q1_2q213517_eval;
case 13530 : return &A4q1_2q213530_eval;
case 13532 : return &A4q1_2q213532_eval;
case 13573 : return &A4q1_2q213573_eval;
case 13575 : return &A4q1_2q213575_eval;
case 13578 : return &A4q1_2q213578_eval;
case 13580 : return &A4q1_2q213580_eval;
case 13585 : return &A4q1_2q213585_eval;
case 13587 : return &A4q1_2q213587_eval;
case 13590 : return &A4q1_2q213590_eval;
case 13592 : return &A4q1_2q213592_eval;
case 13613 : return &A4q1_2q213613_eval;
case 13625 : return &A4q1_2q213625_eval;
case 13638 : return &A4q1_2q213638_eval;
case 13640 : return &A4q1_2q213640_eval;
case 13685 : return &A4q1_2q213685_eval;
case 13697 : return &A4q1_2q213697_eval;
case 13710 : return &A4q1_2q213710_eval;
case 13712 : return &A4q1_2q213712_eval;
case 13788 : return &A4q1_2q213788_eval;
case 13790 : return &A4q1_2q213790_eval;
case 13800 : return &A4q1_2q213800_eval;
case 13802 : return &A4q1_2q213802_eval;
case 14041 : return &A4q1_2q214041_eval;
case 14043 : return &A4q1_2q214043_eval;
case 14046 : return &A4q1_2q214046_eval;
case 14048 : return &A4q1_2q214048_eval;
case 14053 : return &A4q1_2q214053_eval;
case 14055 : return &A4q1_2q214055_eval;
case 14058 : return &A4q1_2q214058_eval;
case 14060 : return &A4q1_2q214060_eval;
case 14076 : return &A4q1_2q214076_eval;
case 14078 : return &A4q1_2q214078_eval;
case 14088 : return &A4q1_2q214088_eval;
case 14090 : return &A4q1_2q214090_eval;
case 14113 : return &A4q1_2q214113_eval;
case 14115 : return &A4q1_2q214115_eval;
case 14118 : return &A4q1_2q214118_eval;
case 14120 : return &A4q1_2q214120_eval;
case 14125 : return &A4q1_2q214125_eval;
case 14127 : return &A4q1_2q214127_eval;
case 14130 : return &A4q1_2q214130_eval;
case 14132 : return &A4q1_2q214132_eval;
case 14148 : return &A4q1_2q214148_eval;
case 14150 : return &A4q1_2q214150_eval;
case 14160 : return &A4q1_2q214160_eval;
case 14162 : return &A4q1_2q214162_eval;
case 14266 : return &A4q1_2q214266_eval;
case 14278 : return &A4q1_2q214278_eval;
case 14281 : return &A4q1_2q214281_eval;
case 14283 : return &A4q1_2q214283_eval;
case 14296 : return &A4q1_2q214296_eval;
case 14308 : return &A4q1_2q214308_eval;
case 14316 : return &A4q1_2q214316_eval;
case 14318 : return &A4q1_2q214318_eval;
case 14338 : return &A4q1_2q214338_eval;
case 14350 : return &A4q1_2q214350_eval;
case 14353 : return &A4q1_2q214353_eval;
case 14355 : return &A4q1_2q214355_eval;
case 14368 : return &A4q1_2q214368_eval;
case 14380 : return &A4q1_2q214380_eval;
case 14388 : return &A4q1_2q214388_eval;
case 14390 : return &A4q1_2q214390_eval;
case 14401 : return &A4q1_2q214401_eval;
case 14403 : return &A4q1_2q214403_eval;
case 14406 : return &A4q1_2q214406_eval;
case 14408 : return &A4q1_2q214408_eval;
case 14413 : return &A4q1_2q214413_eval;
case 14415 : return &A4q1_2q214415_eval;
case 14418 : return &A4q1_2q214418_eval;
case 14420 : return &A4q1_2q214420_eval;
case 14476 : return &A4q1_2q214476_eval;
case 14488 : return &A4q1_2q214488_eval;
case 14496 : return &A4q1_2q214496_eval;
case 14498 : return &A4q1_2q214498_eval;
case 14548 : return &A4q1_2q214548_eval;
case 14560 : return &A4q1_2q214560_eval;
case 14568 : return &A4q1_2q214568_eval;
case 14570 : return &A4q1_2q214570_eval;
case 14616 : return &A4q1_2q214616_eval;
case 14618 : return &A4q1_2q214618_eval;
case 14628 : return &A4q1_2q214628_eval;
case 14630 : return &A4q1_2q214630_eval;
case 14698 : return &A4q1_2q214698_eval;
case 14710 : return &A4q1_2q214710_eval;
case 14713 : return &A4q1_2q214713_eval;
case 14715 : return &A4q1_2q214715_eval;
case 14728 : return &A4q1_2q214728_eval;
case 14740 : return &A4q1_2q214740_eval;
case 14748 : return &A4q1_2q214748_eval;
case 14750 : return &A4q1_2q214750_eval;
case 14770 : return &A4q1_2q214770_eval;
case 14782 : return &A4q1_2q214782_eval;
case 14785 : return &A4q1_2q214785_eval;
case 14787 : return &A4q1_2q214787_eval;
case 14800 : return &A4q1_2q214800_eval;
case 14812 : return &A4q1_2q214812_eval;
case 14820 : return &A4q1_2q214820_eval;
case 14822 : return &A4q1_2q214822_eval;
case 14833 : return &A4q1_2q214833_eval;
case 14835 : return &A4q1_2q214835_eval;
case 14838 : return &A4q1_2q214838_eval;
case 14840 : return &A4q1_2q214840_eval;
case 14845 : return &A4q1_2q214845_eval;
case 14847 : return &A4q1_2q214847_eval;
case 14850 : return &A4q1_2q214850_eval;
case 14852 : return &A4q1_2q214852_eval;
case 14908 : return &A4q1_2q214908_eval;
case 14920 : return &A4q1_2q214920_eval;
case 14928 : return &A4q1_2q214928_eval;
case 14930 : return &A4q1_2q214930_eval;
case 14980 : return &A4q1_2q214980_eval;
case 14992 : return &A4q1_2q214992_eval;
case 15000 : return &A4q1_2q215000_eval;
case 15002 : return &A4q1_2q215002_eval;
case 15048 : return &A4q1_2q215048_eval;
case 15050 : return &A4q1_2q215050_eval;
case 15060 : return &A4q1_2q215060_eval;
case 15062 : return &A4q1_2q215062_eval;
case 15121 : return &A4q1_2q215121_eval;
case 15123 : return &A4q1_2q215123_eval;
case 15126 : return &A4q1_2q215126_eval;
case 15128 : return &A4q1_2q215128_eval;
case 15133 : return &A4q1_2q215133_eval;
case 15135 : return &A4q1_2q215135_eval;
case 15138 : return &A4q1_2q215138_eval;
case 15140 : return &A4q1_2q215140_eval;
case 15156 : return &A4q1_2q215156_eval;
case 15158 : return &A4q1_2q215158_eval;
case 15168 : return &A4q1_2q215168_eval;
case 15170 : return &A4q1_2q215170_eval;
case 15193 : return &A4q1_2q215193_eval;
case 15195 : return &A4q1_2q215195_eval;
case 15198 : return &A4q1_2q215198_eval;
case 15200 : return &A4q1_2q215200_eval;
case 15205 : return &A4q1_2q215205_eval;
case 15207 : return &A4q1_2q215207_eval;
case 15210 : return &A4q1_2q215210_eval;
case 15212 : return &A4q1_2q215212_eval;
case 15228 : return &A4q1_2q215228_eval;
case 15230 : return &A4q1_2q215230_eval;
case 15240 : return &A4q1_2q215240_eval;
case 15242 : return &A4q1_2q215242_eval;
case 15833 : return &A4q1_2q215833_eval;
case 15838 : return &A4q1_2q215838_eval;
case 15905 : return &A4q1_2q215905_eval;
case 15910 : return &A4q1_2q215910_eval;
case 15923 : return &A4q1_2q215923_eval;
case 15935 : return &A4q1_2q215935_eval;
case 15943 : return &A4q1_2q215943_eval;
case 15945 : return &A4q1_2q215945_eval;
case 15958 : return &A4q1_2q215958_eval;
case 15970 : return &A4q1_2q215970_eval;
case 15973 : return &A4q1_2q215973_eval;
case 15975 : return &A4q1_2q215975_eval;
case 16265 : return &A4q1_2q216265_eval;
case 16270 : return &A4q1_2q216270_eval;
case 16337 : return &A4q1_2q216337_eval;
case 16342 : return &A4q1_2q216342_eval;
case 16355 : return &A4q1_2q216355_eval;
case 16367 : return &A4q1_2q216367_eval;
case 16375 : return &A4q1_2q216375_eval;
case 16377 : return &A4q1_2q216377_eval;
case 16390 : return &A4q1_2q216390_eval;
case 16402 : return &A4q1_2q216402_eval;
case 16405 : return &A4q1_2q216405_eval;
case 16407 : return &A4q1_2q216407_eval;
case 16463 : return &A4q1_2q216463_eval;
case 16475 : return &A4q1_2q216475_eval;
case 16483 : return &A4q1_2q216483_eval;
case 16485 : return &A4q1_2q216485_eval;
case 16535 : return &A4q1_2q216535_eval;
case 16547 : return &A4q1_2q216547_eval;
case 16555 : return &A4q1_2q216555_eval;
case 16557 : return &A4q1_2q216557_eval;
case 16603 : return &A4q1_2q216603_eval;
case 16605 : return &A4q1_2q216605_eval;
case 16615 : return &A4q1_2q216615_eval;
case 16617 : return &A4q1_2q216617_eval;
case 16678 : return &A4q1_2q216678_eval;
case 16690 : return &A4q1_2q216690_eval;
case 16693 : return &A4q1_2q216693_eval;
case 16695 : return &A4q1_2q216695_eval;
case 16750 : return &A4q1_2q216750_eval;
case 16762 : return &A4q1_2q216762_eval;
case 16765 : return &A4q1_2q216765_eval;
case 16767 : return &A4q1_2q216767_eval;
case 16783 : return &A4q1_2q216783_eval;
case 16785 : return &A4q1_2q216785_eval;
case 16795 : return &A4q1_2q216795_eval;
case 16797 : return &A4q1_2q216797_eval;
case 16913 : return &A4q1_2q216913_eval;
case 16918 : return &A4q1_2q216918_eval;
case 16985 : return &A4q1_2q216985_eval;
case 16990 : return &A4q1_2q216990_eval;
case 17003 : return &A4q1_2q217003_eval;
case 17015 : return &A4q1_2q217015_eval;
case 17023 : return &A4q1_2q217023_eval;
case 17025 : return &A4q1_2q217025_eval;
case 17038 : return &A4q1_2q217038_eval;
case 17050 : return &A4q1_2q217050_eval;
case 17053 : return &A4q1_2q217053_eval;
case 17055 : return &A4q1_2q217055_eval;
case 17093 : return &A4q1_2q217093_eval;
case 17098 : return &A4q1_2q217098_eval;
case 17165 : return &A4q1_2q217165_eval;
case 17170 : return &A4q1_2q217170_eval;
case 17213 : return &A4q1_2q217213_eval;
case 17225 : return &A4q1_2q217225_eval;
case 17238 : return &A4q1_2q217238_eval;
case 17240 : return &A4q1_2q217240_eval;
case 17248 : return &A4q1_2q217248_eval;
case 17260 : return &A4q1_2q217260_eval;
case 17268 : return &A4q1_2q217268_eval;
case 17270 : return &A4q1_2q217270_eval;
case 17345 : return &A4q1_2q217345_eval;
case 17350 : return &A4q1_2q217350_eval;
case 17417 : return &A4q1_2q217417_eval;
case 17422 : return &A4q1_2q217422_eval;
case 17435 : return &A4q1_2q217435_eval;
case 17447 : return &A4q1_2q217447_eval;
case 17455 : return &A4q1_2q217455_eval;
case 17457 : return &A4q1_2q217457_eval;
case 17470 : return &A4q1_2q217470_eval;
case 17482 : return &A4q1_2q217482_eval;
case 17485 : return &A4q1_2q217485_eval;
case 17487 : return &A4q1_2q217487_eval;
case 17525 : return &A4q1_2q217525_eval;
case 17530 : return &A4q1_2q217530_eval;
case 17597 : return &A4q1_2q217597_eval;
case 17602 : return &A4q1_2q217602_eval;
case 17645 : return &A4q1_2q217645_eval;
case 17657 : return &A4q1_2q217657_eval;
case 17670 : return &A4q1_2q217670_eval;
case 17672 : return &A4q1_2q217672_eval;
case 17680 : return &A4q1_2q217680_eval;
case 17692 : return &A4q1_2q217692_eval;
case 17700 : return &A4q1_2q217700_eval;
case 17702 : return &A4q1_2q217702_eval;
case 17723 : return &A4q1_2q217723_eval;
case 17735 : return &A4q1_2q217735_eval;
case 17743 : return &A4q1_2q217743_eval;
case 17745 : return &A4q1_2q217745_eval;
case 17753 : return &A4q1_2q217753_eval;
case 17765 : return &A4q1_2q217765_eval;
case 17778 : return &A4q1_2q217778_eval;
case 17780 : return &A4q1_2q217780_eval;
case 17795 : return &A4q1_2q217795_eval;
case 17807 : return &A4q1_2q217807_eval;
case 17815 : return &A4q1_2q217815_eval;
case 17817 : return &A4q1_2q217817_eval;
case 17825 : return &A4q1_2q217825_eval;
case 17837 : return &A4q1_2q217837_eval;
case 17850 : return &A4q1_2q217850_eval;
case 17852 : return &A4q1_2q217852_eval;
case 17893 : return &A4q1_2q217893_eval;
case 17895 : return &A4q1_2q217895_eval;
case 17898 : return &A4q1_2q217898_eval;
case 17900 : return &A4q1_2q217900_eval;
case 17905 : return &A4q1_2q217905_eval;
case 17907 : return &A4q1_2q217907_eval;
case 17910 : return &A4q1_2q217910_eval;
case 17912 : return &A4q1_2q217912_eval;
case 17938 : return &A4q1_2q217938_eval;
case 17950 : return &A4q1_2q217950_eval;
case 17953 : return &A4q1_2q217953_eval;
case 17955 : return &A4q1_2q217955_eval;
case 17968 : return &A4q1_2q217968_eval;
case 17980 : return &A4q1_2q217980_eval;
case 17988 : return &A4q1_2q217988_eval;
case 17990 : return &A4q1_2q217990_eval;
case 18010 : return &A4q1_2q218010_eval;
case 18022 : return &A4q1_2q218022_eval;
case 18025 : return &A4q1_2q218025_eval;
case 18027 : return &A4q1_2q218027_eval;
case 18040 : return &A4q1_2q218040_eval;
case 18052 : return &A4q1_2q218052_eval;
case 18060 : return &A4q1_2q218060_eval;
case 18062 : return &A4q1_2q218062_eval;
case 18073 : return &A4q1_2q218073_eval;
case 18075 : return &A4q1_2q218075_eval;
case 18078 : return &A4q1_2q218078_eval;
case 18080 : return &A4q1_2q218080_eval;
case 18085 : return &A4q1_2q218085_eval;
case 18087 : return &A4q1_2q218087_eval;
case 18090 : return &A4q1_2q218090_eval;
case 18092 : return &A4q1_2q218092_eval;
case 18425 : return &A4q1_2q218425_eval;
case 18430 : return &A4q1_2q218430_eval;
case 18497 : return &A4q1_2q218497_eval;
case 18502 : return &A4q1_2q218502_eval;
case 18515 : return &A4q1_2q218515_eval;
case 18527 : return &A4q1_2q218527_eval;
case 18535 : return &A4q1_2q218535_eval;
case 18537 : return &A4q1_2q218537_eval;
case 18550 : return &A4q1_2q218550_eval;
case 18562 : return &A4q1_2q218562_eval;
case 18565 : return &A4q1_2q218565_eval;
case 18567 : return &A4q1_2q218567_eval;
case 18857 : return &A4q1_2q218857_eval;
case 18862 : return &A4q1_2q218862_eval;
case 18929 : return &A4q1_2q218929_eval;
case 18934 : return &A4q1_2q218934_eval;
case 18947 : return &A4q1_2q218947_eval;
case 18959 : return &A4q1_2q218959_eval;
case 18967 : return &A4q1_2q218967_eval;
case 18969 : return &A4q1_2q218969_eval;
case 18982 : return &A4q1_2q218982_eval;
case 18994 : return &A4q1_2q218994_eval;
case 18997 : return &A4q1_2q218997_eval;
case 18999 : return &A4q1_2q218999_eval;
case 19055 : return &A4q1_2q219055_eval;
case 19067 : return &A4q1_2q219067_eval;
case 19075 : return &A4q1_2q219075_eval;
case 19077 : return &A4q1_2q219077_eval;
case 19127 : return &A4q1_2q219127_eval;
case 19139 : return &A4q1_2q219139_eval;
case 19147 : return &A4q1_2q219147_eval;
case 19149 : return &A4q1_2q219149_eval;
case 19195 : return &A4q1_2q219195_eval;
case 19197 : return &A4q1_2q219197_eval;
case 19207 : return &A4q1_2q219207_eval;
case 19209 : return &A4q1_2q219209_eval;
case 19270 : return &A4q1_2q219270_eval;
case 19282 : return &A4q1_2q219282_eval;
case 19285 : return &A4q1_2q219285_eval;
case 19287 : return &A4q1_2q219287_eval;
case 19342 : return &A4q1_2q219342_eval;
case 19354 : return &A4q1_2q219354_eval;
case 19357 : return &A4q1_2q219357_eval;
case 19359 : return &A4q1_2q219359_eval;
case 19375 : return &A4q1_2q219375_eval;
case 19377 : return &A4q1_2q219377_eval;
case 19387 : return &A4q1_2q219387_eval;
case 19389 : return &A4q1_2q219389_eval;
case 19505 : return &A4q1_2q219505_eval;
case 19510 : return &A4q1_2q219510_eval;
case 19577 : return &A4q1_2q219577_eval;
case 19582 : return &A4q1_2q219582_eval;
case 19595 : return &A4q1_2q219595_eval;
case 19607 : return &A4q1_2q219607_eval;
case 19615 : return &A4q1_2q219615_eval;
case 19617 : return &A4q1_2q219617_eval;
case 19630 : return &A4q1_2q219630_eval;
case 19642 : return &A4q1_2q219642_eval;
case 19645 : return &A4q1_2q219645_eval;
case 19647 : return &A4q1_2q219647_eval;
case 19685 : return &A4q1_2q219685_eval;
case 19690 : return &A4q1_2q219690_eval;
case 19757 : return &A4q1_2q219757_eval;
case 19762 : return &A4q1_2q219762_eval;
case 19805 : return &A4q1_2q219805_eval;
case 19817 : return &A4q1_2q219817_eval;
case 19830 : return &A4q1_2q219830_eval;
case 19832 : return &A4q1_2q219832_eval;
case 19840 : return &A4q1_2q219840_eval;
case 19852 : return &A4q1_2q219852_eval;
case 19860 : return &A4q1_2q219860_eval;
case 19862 : return &A4q1_2q219862_eval;
case 19937 : return &A4q1_2q219937_eval;
case 19942 : return &A4q1_2q219942_eval;
case 20009 : return &A4q1_2q220009_eval;
case 20014 : return &A4q1_2q220014_eval;
case 20027 : return &A4q1_2q220027_eval;
case 20039 : return &A4q1_2q220039_eval;
case 20047 : return &A4q1_2q220047_eval;
case 20049 : return &A4q1_2q220049_eval;
case 20062 : return &A4q1_2q220062_eval;
case 20074 : return &A4q1_2q220074_eval;
case 20077 : return &A4q1_2q220077_eval;
case 20079 : return &A4q1_2q220079_eval;
case 20117 : return &A4q1_2q220117_eval;
case 20122 : return &A4q1_2q220122_eval;
case 20189 : return &A4q1_2q220189_eval;
case 20194 : return &A4q1_2q220194_eval;
case 20237 : return &A4q1_2q220237_eval;
case 20249 : return &A4q1_2q220249_eval;
case 20262 : return &A4q1_2q220262_eval;
case 20264 : return &A4q1_2q220264_eval;
case 20272 : return &A4q1_2q220272_eval;
case 20284 : return &A4q1_2q220284_eval;
case 20292 : return &A4q1_2q220292_eval;
case 20294 : return &A4q1_2q220294_eval;
case 20315 : return &A4q1_2q220315_eval;
case 20327 : return &A4q1_2q220327_eval;
case 20335 : return &A4q1_2q220335_eval;
case 20337 : return &A4q1_2q220337_eval;
case 20345 : return &A4q1_2q220345_eval;
case 20357 : return &A4q1_2q220357_eval;
case 20370 : return &A4q1_2q220370_eval;
case 20372 : return &A4q1_2q220372_eval;
case 20387 : return &A4q1_2q220387_eval;
case 20399 : return &A4q1_2q220399_eval;
case 20407 : return &A4q1_2q220407_eval;
case 20409 : return &A4q1_2q220409_eval;
case 20417 : return &A4q1_2q220417_eval;
case 20429 : return &A4q1_2q220429_eval;
case 20442 : return &A4q1_2q220442_eval;
case 20444 : return &A4q1_2q220444_eval;
case 20485 : return &A4q1_2q220485_eval;
case 20487 : return &A4q1_2q220487_eval;
case 20490 : return &A4q1_2q220490_eval;
case 20492 : return &A4q1_2q220492_eval;
case 20497 : return &A4q1_2q220497_eval;
case 20499 : return &A4q1_2q220499_eval;
case 20502 : return &A4q1_2q220502_eval;
case 20504 : return &A4q1_2q220504_eval;
case 20530 : return &A4q1_2q220530_eval;
case 20542 : return &A4q1_2q220542_eval;
case 20545 : return &A4q1_2q220545_eval;
case 20547 : return &A4q1_2q220547_eval;
case 20560 : return &A4q1_2q220560_eval;
case 20572 : return &A4q1_2q220572_eval;
case 20580 : return &A4q1_2q220580_eval;
case 20582 : return &A4q1_2q220582_eval;
case 20602 : return &A4q1_2q220602_eval;
case 20614 : return &A4q1_2q220614_eval;
case 20617 : return &A4q1_2q220617_eval;
case 20619 : return &A4q1_2q220619_eval;
case 20632 : return &A4q1_2q220632_eval;
case 20644 : return &A4q1_2q220644_eval;
case 20652 : return &A4q1_2q220652_eval;
case 20654 : return &A4q1_2q220654_eval;
case 20665 : return &A4q1_2q220665_eval;
case 20667 : return &A4q1_2q220667_eval;
case 20670 : return &A4q1_2q220670_eval;
case 20672 : return &A4q1_2q220672_eval;
case 20677 : return &A4q1_2q220677_eval;
case 20679 : return &A4q1_2q220679_eval;
case 20682 : return &A4q1_2q220682_eval;
case 20684 : return &A4q1_2q220684_eval;
case 20783 : return &A4q1_2q220783_eval;
case 20795 : return &A4q1_2q220795_eval;
case 20803 : return &A4q1_2q220803_eval;
case 20805 : return &A4q1_2q220805_eval;
case 20855 : return &A4q1_2q220855_eval;
case 20867 : return &A4q1_2q220867_eval;
case 20875 : return &A4q1_2q220875_eval;
case 20877 : return &A4q1_2q220877_eval;
case 20923 : return &A4q1_2q220923_eval;
case 20925 : return &A4q1_2q220925_eval;
case 20935 : return &A4q1_2q220935_eval;
case 20937 : return &A4q1_2q220937_eval;
case 20963 : return &A4q1_2q220963_eval;
case 20975 : return &A4q1_2q220975_eval;
case 20983 : return &A4q1_2q220983_eval;
case 20985 : return &A4q1_2q220985_eval;
case 20993 : return &A4q1_2q220993_eval;
case 21005 : return &A4q1_2q221005_eval;
case 21018 : return &A4q1_2q221018_eval;
case 21020 : return &A4q1_2q221020_eval;
case 21035 : return &A4q1_2q221035_eval;
case 21047 : return &A4q1_2q221047_eval;
case 21055 : return &A4q1_2q221055_eval;
case 21057 : return &A4q1_2q221057_eval;
case 21065 : return &A4q1_2q221065_eval;
case 21077 : return &A4q1_2q221077_eval;
case 21090 : return &A4q1_2q221090_eval;
case 21092 : return &A4q1_2q221092_eval;
case 21133 : return &A4q1_2q221133_eval;
case 21135 : return &A4q1_2q221135_eval;
case 21138 : return &A4q1_2q221138_eval;
case 21140 : return &A4q1_2q221140_eval;
case 21145 : return &A4q1_2q221145_eval;
case 21147 : return &A4q1_2q221147_eval;
case 21150 : return &A4q1_2q221150_eval;
case 21152 : return &A4q1_2q221152_eval;
case 21215 : return &A4q1_2q221215_eval;
case 21227 : return &A4q1_2q221227_eval;
case 21235 : return &A4q1_2q221235_eval;
case 21237 : return &A4q1_2q221237_eval;
case 21287 : return &A4q1_2q221287_eval;
case 21299 : return &A4q1_2q221299_eval;
case 21307 : return &A4q1_2q221307_eval;
case 21309 : return &A4q1_2q221309_eval;
case 21355 : return &A4q1_2q221355_eval;
case 21357 : return &A4q1_2q221357_eval;
case 21367 : return &A4q1_2q221367_eval;
case 21369 : return &A4q1_2q221369_eval;
case 21395 : return &A4q1_2q221395_eval;
case 21407 : return &A4q1_2q221407_eval;
case 21415 : return &A4q1_2q221415_eval;
case 21417 : return &A4q1_2q221417_eval;
case 21425 : return &A4q1_2q221425_eval;
case 21437 : return &A4q1_2q221437_eval;
case 21450 : return &A4q1_2q221450_eval;
case 21452 : return &A4q1_2q221452_eval;
case 21467 : return &A4q1_2q221467_eval;
case 21479 : return &A4q1_2q221479_eval;
case 21487 : return &A4q1_2q221487_eval;
case 21489 : return &A4q1_2q221489_eval;
case 21497 : return &A4q1_2q221497_eval;
case 21509 : return &A4q1_2q221509_eval;
case 21522 : return &A4q1_2q221522_eval;
case 21524 : return &A4q1_2q221524_eval;
case 21565 : return &A4q1_2q221565_eval;
case 21567 : return &A4q1_2q221567_eval;
case 21570 : return &A4q1_2q221570_eval;
case 21572 : return &A4q1_2q221572_eval;
case 21577 : return &A4q1_2q221577_eval;
case 21579 : return &A4q1_2q221579_eval;
case 21582 : return &A4q1_2q221582_eval;
case 21584 : return &A4q1_2q221584_eval;
case 21823 : return &A4q1_2q221823_eval;
case 21825 : return &A4q1_2q221825_eval;
case 21835 : return &A4q1_2q221835_eval;
case 21837 : return &A4q1_2q221837_eval;
case 21853 : return &A4q1_2q221853_eval;
case 21855 : return &A4q1_2q221855_eval;
case 21858 : return &A4q1_2q221858_eval;
case 21860 : return &A4q1_2q221860_eval;
case 21865 : return &A4q1_2q221865_eval;
case 21867 : return &A4q1_2q221867_eval;
case 21870 : return &A4q1_2q221870_eval;
case 21872 : return &A4q1_2q221872_eval;
case 21895 : return &A4q1_2q221895_eval;
case 21897 : return &A4q1_2q221897_eval;
case 21907 : return &A4q1_2q221907_eval;
case 21909 : return &A4q1_2q221909_eval;
case 21925 : return &A4q1_2q221925_eval;
case 21927 : return &A4q1_2q221927_eval;
case 21930 : return &A4q1_2q221930_eval;
case 21932 : return &A4q1_2q221932_eval;
case 21937 : return &A4q1_2q221937_eval;
case 21939 : return &A4q1_2q221939_eval;
case 21942 : return &A4q1_2q221942_eval;
case 21944 : return &A4q1_2q221944_eval;
case 22078 : return &A4q1_2q222078_eval;
case 22090 : return &A4q1_2q222090_eval;
case 22093 : return &A4q1_2q222093_eval;
case 22095 : return &A4q1_2q222095_eval;
case 22150 : return &A4q1_2q222150_eval;
case 22162 : return &A4q1_2q222162_eval;
case 22165 : return &A4q1_2q222165_eval;
case 22167 : return &A4q1_2q222167_eval;
case 22183 : return &A4q1_2q222183_eval;
case 22185 : return &A4q1_2q222185_eval;
case 22195 : return &A4q1_2q222195_eval;
case 22197 : return &A4q1_2q222197_eval;
case 22258 : return &A4q1_2q222258_eval;
case 22270 : return &A4q1_2q222270_eval;
case 22273 : return &A4q1_2q222273_eval;
case 22275 : return &A4q1_2q222275_eval;
case 22288 : return &A4q1_2q222288_eval;
case 22300 : return &A4q1_2q222300_eval;
case 22308 : return &A4q1_2q222308_eval;
case 22310 : return &A4q1_2q222310_eval;
case 22330 : return &A4q1_2q222330_eval;
case 22342 : return &A4q1_2q222342_eval;
case 22345 : return &A4q1_2q222345_eval;
case 22347 : return &A4q1_2q222347_eval;
case 22360 : return &A4q1_2q222360_eval;
case 22372 : return &A4q1_2q222372_eval;
case 22380 : return &A4q1_2q222380_eval;
case 22382 : return &A4q1_2q222382_eval;
case 22393 : return &A4q1_2q222393_eval;
case 22395 : return &A4q1_2q222395_eval;
case 22398 : return &A4q1_2q222398_eval;
case 22400 : return &A4q1_2q222400_eval;
case 22405 : return &A4q1_2q222405_eval;
case 22407 : return &A4q1_2q222407_eval;
case 22410 : return &A4q1_2q222410_eval;
case 22412 : return &A4q1_2q222412_eval;
case 22510 : return &A4q1_2q222510_eval;
case 22522 : return &A4q1_2q222522_eval;
case 22525 : return &A4q1_2q222525_eval;
case 22527 : return &A4q1_2q222527_eval;
case 22582 : return &A4q1_2q222582_eval;
case 22594 : return &A4q1_2q222594_eval;
case 22597 : return &A4q1_2q222597_eval;
case 22599 : return &A4q1_2q222599_eval;
case 22615 : return &A4q1_2q222615_eval;
case 22617 : return &A4q1_2q222617_eval;
case 22627 : return &A4q1_2q222627_eval;
case 22629 : return &A4q1_2q222629_eval;
case 22690 : return &A4q1_2q222690_eval;
case 22702 : return &A4q1_2q222702_eval;
case 22705 : return &A4q1_2q222705_eval;
case 22707 : return &A4q1_2q222707_eval;
case 22720 : return &A4q1_2q222720_eval;
case 22732 : return &A4q1_2q222732_eval;
case 22740 : return &A4q1_2q222740_eval;
case 22742 : return &A4q1_2q222742_eval;
case 22762 : return &A4q1_2q222762_eval;
case 22774 : return &A4q1_2q222774_eval;
case 22777 : return &A4q1_2q222777_eval;
case 22779 : return &A4q1_2q222779_eval;
case 22792 : return &A4q1_2q222792_eval;
case 22804 : return &A4q1_2q222804_eval;
case 22812 : return &A4q1_2q222812_eval;
case 22814 : return &A4q1_2q222814_eval;
case 22825 : return &A4q1_2q222825_eval;
case 22827 : return &A4q1_2q222827_eval;
case 22830 : return &A4q1_2q222830_eval;
case 22832 : return &A4q1_2q222832_eval;
case 22837 : return &A4q1_2q222837_eval;
case 22839 : return &A4q1_2q222839_eval;
case 22842 : return &A4q1_2q222842_eval;
case 22844 : return &A4q1_2q222844_eval;
case 22903 : return &A4q1_2q222903_eval;
case 22905 : return &A4q1_2q222905_eval;
case 22915 : return &A4q1_2q222915_eval;
case 22917 : return &A4q1_2q222917_eval;
case 22933 : return &A4q1_2q222933_eval;
case 22935 : return &A4q1_2q222935_eval;
case 22938 : return &A4q1_2q222938_eval;
case 22940 : return &A4q1_2q222940_eval;
case 22945 : return &A4q1_2q222945_eval;
case 22947 : return &A4q1_2q222947_eval;
case 22950 : return &A4q1_2q222950_eval;
case 22952 : return &A4q1_2q222952_eval;
case 22975 : return &A4q1_2q222975_eval;
case 22977 : return &A4q1_2q222977_eval;
case 22987 : return &A4q1_2q222987_eval;
case 22989 : return &A4q1_2q222989_eval;
case 23005 : return &A4q1_2q223005_eval;
case 23007 : return &A4q1_2q223007_eval;
case 23010 : return &A4q1_2q223010_eval;
case 23012 : return &A4q1_2q223012_eval;
case 23017 : return &A4q1_2q223017_eval;
case 23019 : return &A4q1_2q223019_eval;
case 23022 : return &A4q1_2q223022_eval;
case 23024 : return &A4q1_2q223024_eval;
case 23393 : return &A4q1_2q223393_eval;
case 23398 : return &A4q1_2q223398_eval;
case 23465 : return &A4q1_2q223465_eval;
case 23470 : return &A4q1_2q223470_eval;
case 23483 : return &A4q1_2q223483_eval;
case 23495 : return &A4q1_2q223495_eval;
case 23503 : return &A4q1_2q223503_eval;
case 23505 : return &A4q1_2q223505_eval;
case 23518 : return &A4q1_2q223518_eval;
case 23530 : return &A4q1_2q223530_eval;
case 23533 : return &A4q1_2q223533_eval;
case 23535 : return &A4q1_2q223535_eval;
case 23573 : return &A4q1_2q223573_eval;
case 23578 : return &A4q1_2q223578_eval;
case 23645 : return &A4q1_2q223645_eval;
case 23650 : return &A4q1_2q223650_eval;
case 23693 : return &A4q1_2q223693_eval;
case 23705 : return &A4q1_2q223705_eval;
case 23718 : return &A4q1_2q223718_eval;
case 23720 : return &A4q1_2q223720_eval;
case 23728 : return &A4q1_2q223728_eval;
case 23740 : return &A4q1_2q223740_eval;
case 23748 : return &A4q1_2q223748_eval;
case 23750 : return &A4q1_2q223750_eval;
case 23825 : return &A4q1_2q223825_eval;
case 23830 : return &A4q1_2q223830_eval;
case 23897 : return &A4q1_2q223897_eval;
case 23902 : return &A4q1_2q223902_eval;
case 23915 : return &A4q1_2q223915_eval;
case 23927 : return &A4q1_2q223927_eval;
case 23935 : return &A4q1_2q223935_eval;
case 23937 : return &A4q1_2q223937_eval;
case 23950 : return &A4q1_2q223950_eval;
case 23962 : return &A4q1_2q223962_eval;
case 23965 : return &A4q1_2q223965_eval;
case 23967 : return &A4q1_2q223967_eval;
case 24005 : return &A4q1_2q224005_eval;
case 24010 : return &A4q1_2q224010_eval;
case 24077 : return &A4q1_2q224077_eval;
case 24082 : return &A4q1_2q224082_eval;
case 24125 : return &A4q1_2q224125_eval;
case 24137 : return &A4q1_2q224137_eval;
case 24150 : return &A4q1_2q224150_eval;
case 24152 : return &A4q1_2q224152_eval;
case 24160 : return &A4q1_2q224160_eval;
case 24172 : return &A4q1_2q224172_eval;
case 24180 : return &A4q1_2q224180_eval;
case 24182 : return &A4q1_2q224182_eval;
case 24203 : return &A4q1_2q224203_eval;
case 24215 : return &A4q1_2q224215_eval;
case 24223 : return &A4q1_2q224223_eval;
case 24225 : return &A4q1_2q224225_eval;
case 24233 : return &A4q1_2q224233_eval;
case 24245 : return &A4q1_2q224245_eval;
case 24258 : return &A4q1_2q224258_eval;
case 24260 : return &A4q1_2q224260_eval;
case 24275 : return &A4q1_2q224275_eval;
case 24287 : return &A4q1_2q224287_eval;
case 24295 : return &A4q1_2q224295_eval;
case 24297 : return &A4q1_2q224297_eval;
case 24305 : return &A4q1_2q224305_eval;
case 24317 : return &A4q1_2q224317_eval;
case 24330 : return &A4q1_2q224330_eval;
case 24332 : return &A4q1_2q224332_eval;
case 24373 : return &A4q1_2q224373_eval;
case 24375 : return &A4q1_2q224375_eval;
case 24378 : return &A4q1_2q224378_eval;
case 24380 : return &A4q1_2q224380_eval;
case 24385 : return &A4q1_2q224385_eval;
case 24387 : return &A4q1_2q224387_eval;
case 24390 : return &A4q1_2q224390_eval;
case 24392 : return &A4q1_2q224392_eval;
case 24418 : return &A4q1_2q224418_eval;
case 24430 : return &A4q1_2q224430_eval;
case 24433 : return &A4q1_2q224433_eval;
case 24435 : return &A4q1_2q224435_eval;
case 24448 : return &A4q1_2q224448_eval;
case 24460 : return &A4q1_2q224460_eval;
case 24468 : return &A4q1_2q224468_eval;
case 24470 : return &A4q1_2q224470_eval;
case 24490 : return &A4q1_2q224490_eval;
case 24502 : return &A4q1_2q224502_eval;
case 24505 : return &A4q1_2q224505_eval;
case 24507 : return &A4q1_2q224507_eval;
case 24520 : return &A4q1_2q224520_eval;
case 24532 : return &A4q1_2q224532_eval;
case 24540 : return &A4q1_2q224540_eval;
case 24542 : return &A4q1_2q224542_eval;
case 24553 : return &A4q1_2q224553_eval;
case 24555 : return &A4q1_2q224555_eval;
case 24558 : return &A4q1_2q224558_eval;
case 24560 : return &A4q1_2q224560_eval;
case 24565 : return &A4q1_2q224565_eval;
case 24567 : return &A4q1_2q224567_eval;
case 24570 : return &A4q1_2q224570_eval;
case 24572 : return &A4q1_2q224572_eval;
case 24653 : return &A4q1_2q224653_eval;
case 24658 : return &A4q1_2q224658_eval;
case 24725 : return &A4q1_2q224725_eval;
case 24730 : return &A4q1_2q224730_eval;
case 24773 : return &A4q1_2q224773_eval;
case 24785 : return &A4q1_2q224785_eval;
case 24798 : return &A4q1_2q224798_eval;
case 24800 : return &A4q1_2q224800_eval;
case 24808 : return &A4q1_2q224808_eval;
case 24820 : return &A4q1_2q224820_eval;
case 24828 : return &A4q1_2q224828_eval;
case 24830 : return &A4q1_2q224830_eval;
case 25085 : return &A4q1_2q225085_eval;
case 25090 : return &A4q1_2q225090_eval;
case 25157 : return &A4q1_2q225157_eval;
case 25162 : return &A4q1_2q225162_eval;
case 25205 : return &A4q1_2q225205_eval;
case 25217 : return &A4q1_2q225217_eval;
case 25230 : return &A4q1_2q225230_eval;
case 25232 : return &A4q1_2q225232_eval;
case 25240 : return &A4q1_2q225240_eval;
case 25252 : return &A4q1_2q225252_eval;
case 25260 : return &A4q1_2q225260_eval;
case 25262 : return &A4q1_2q225262_eval;
case 25493 : return &A4q1_2q225493_eval;
case 25505 : return &A4q1_2q225505_eval;
case 25518 : return &A4q1_2q225518_eval;
case 25520 : return &A4q1_2q225520_eval;
case 25565 : return &A4q1_2q225565_eval;
case 25577 : return &A4q1_2q225577_eval;
case 25590 : return &A4q1_2q225590_eval;
case 25592 : return &A4q1_2q225592_eval;
case 25668 : return &A4q1_2q225668_eval;
case 25670 : return &A4q1_2q225670_eval;
case 25680 : return &A4q1_2q225680_eval;
case 25682 : return &A4q1_2q225682_eval;
case 25708 : return &A4q1_2q225708_eval;
case 25720 : return &A4q1_2q225720_eval;
case 25728 : return &A4q1_2q225728_eval;
case 25730 : return &A4q1_2q225730_eval;
case 25780 : return &A4q1_2q225780_eval;
case 25792 : return &A4q1_2q225792_eval;
case 25800 : return &A4q1_2q225800_eval;
case 25802 : return &A4q1_2q225802_eval;
case 25848 : return &A4q1_2q225848_eval;
case 25850 : return &A4q1_2q225850_eval;
case 25860 : return &A4q1_2q225860_eval;
case 25862 : return &A4q1_2q225862_eval;
case 25985 : return &A4q1_2q225985_eval;
case 25990 : return &A4q1_2q225990_eval;
case 26057 : return &A4q1_2q226057_eval;
case 26062 : return &A4q1_2q226062_eval;
case 26075 : return &A4q1_2q226075_eval;
case 26087 : return &A4q1_2q226087_eval;
case 26095 : return &A4q1_2q226095_eval;
case 26097 : return &A4q1_2q226097_eval;
case 26110 : return &A4q1_2q226110_eval;
case 26122 : return &A4q1_2q226122_eval;
case 26125 : return &A4q1_2q226125_eval;
case 26127 : return &A4q1_2q226127_eval;
case 26165 : return &A4q1_2q226165_eval;
case 26170 : return &A4q1_2q226170_eval;
case 26237 : return &A4q1_2q226237_eval;
case 26242 : return &A4q1_2q226242_eval;
case 26285 : return &A4q1_2q226285_eval;
case 26297 : return &A4q1_2q226297_eval;
case 26310 : return &A4q1_2q226310_eval;
case 26312 : return &A4q1_2q226312_eval;
case 26320 : return &A4q1_2q226320_eval;
case 26332 : return &A4q1_2q226332_eval;
case 26340 : return &A4q1_2q226340_eval;
case 26342 : return &A4q1_2q226342_eval;
case 26417 : return &A4q1_2q226417_eval;
case 26422 : return &A4q1_2q226422_eval;
case 26489 : return &A4q1_2q226489_eval;
case 26494 : return &A4q1_2q226494_eval;
case 26507 : return &A4q1_2q226507_eval;
case 26519 : return &A4q1_2q226519_eval;
case 26527 : return &A4q1_2q226527_eval;
case 26529 : return &A4q1_2q226529_eval;
case 26542 : return &A4q1_2q226542_eval;
case 26554 : return &A4q1_2q226554_eval;
case 26557 : return &A4q1_2q226557_eval;
case 26559 : return &A4q1_2q226559_eval;
case 26597 : return &A4q1_2q226597_eval;
case 26602 : return &A4q1_2q226602_eval;
case 26669 : return &A4q1_2q226669_eval;
case 26674 : return &A4q1_2q226674_eval;
case 26717 : return &A4q1_2q226717_eval;
case 26729 : return &A4q1_2q226729_eval;
case 26742 : return &A4q1_2q226742_eval;
case 26744 : return &A4q1_2q226744_eval;
case 26752 : return &A4q1_2q226752_eval;
case 26764 : return &A4q1_2q226764_eval;
case 26772 : return &A4q1_2q226772_eval;
case 26774 : return &A4q1_2q226774_eval;
case 26795 : return &A4q1_2q226795_eval;
case 26807 : return &A4q1_2q226807_eval;
case 26815 : return &A4q1_2q226815_eval;
case 26817 : return &A4q1_2q226817_eval;
case 26825 : return &A4q1_2q226825_eval;
case 26837 : return &A4q1_2q226837_eval;
case 26850 : return &A4q1_2q226850_eval;
case 26852 : return &A4q1_2q226852_eval;
case 26867 : return &A4q1_2q226867_eval;
case 26879 : return &A4q1_2q226879_eval;
case 26887 : return &A4q1_2q226887_eval;
case 26889 : return &A4q1_2q226889_eval;
case 26897 : return &A4q1_2q226897_eval;
case 26909 : return &A4q1_2q226909_eval;
case 26922 : return &A4q1_2q226922_eval;
case 26924 : return &A4q1_2q226924_eval;
case 26965 : return &A4q1_2q226965_eval;
case 26967 : return &A4q1_2q226967_eval;
case 26970 : return &A4q1_2q226970_eval;
case 26972 : return &A4q1_2q226972_eval;
case 26977 : return &A4q1_2q226977_eval;
case 26979 : return &A4q1_2q226979_eval;
case 26982 : return &A4q1_2q226982_eval;
case 26984 : return &A4q1_2q226984_eval;
case 27010 : return &A4q1_2q227010_eval;
case 27022 : return &A4q1_2q227022_eval;
case 27025 : return &A4q1_2q227025_eval;
case 27027 : return &A4q1_2q227027_eval;
case 27040 : return &A4q1_2q227040_eval;
case 27052 : return &A4q1_2q227052_eval;
case 27060 : return &A4q1_2q227060_eval;
case 27062 : return &A4q1_2q227062_eval;
case 27082 : return &A4q1_2q227082_eval;
case 27094 : return &A4q1_2q227094_eval;
case 27097 : return &A4q1_2q227097_eval;
case 27099 : return &A4q1_2q227099_eval;
case 27112 : return &A4q1_2q227112_eval;
case 27124 : return &A4q1_2q227124_eval;
case 27132 : return &A4q1_2q227132_eval;
case 27134 : return &A4q1_2q227134_eval;
case 27145 : return &A4q1_2q227145_eval;
case 27147 : return &A4q1_2q227147_eval;
case 27150 : return &A4q1_2q227150_eval;
case 27152 : return &A4q1_2q227152_eval;
case 27157 : return &A4q1_2q227157_eval;
case 27159 : return &A4q1_2q227159_eval;
case 27162 : return &A4q1_2q227162_eval;
case 27164 : return &A4q1_2q227164_eval;
case 27245 : return &A4q1_2q227245_eval;
case 27250 : return &A4q1_2q227250_eval;
case 27317 : return &A4q1_2q227317_eval;
case 27322 : return &A4q1_2q227322_eval;
case 27365 : return &A4q1_2q227365_eval;
case 27377 : return &A4q1_2q227377_eval;
case 27390 : return &A4q1_2q227390_eval;
case 27392 : return &A4q1_2q227392_eval;
case 27400 : return &A4q1_2q227400_eval;
case 27412 : return &A4q1_2q227412_eval;
case 27420 : return &A4q1_2q227420_eval;
case 27422 : return &A4q1_2q227422_eval;
case 27677 : return &A4q1_2q227677_eval;
case 27682 : return &A4q1_2q227682_eval;
case 27749 : return &A4q1_2q227749_eval;
case 27754 : return &A4q1_2q227754_eval;
case 27797 : return &A4q1_2q227797_eval;
case 27809 : return &A4q1_2q227809_eval;
case 27822 : return &A4q1_2q227822_eval;
case 27824 : return &A4q1_2q227824_eval;
case 27832 : return &A4q1_2q227832_eval;
case 27844 : return &A4q1_2q227844_eval;
case 27852 : return &A4q1_2q227852_eval;
case 27854 : return &A4q1_2q227854_eval;
case 28085 : return &A4q1_2q228085_eval;
case 28097 : return &A4q1_2q228097_eval;
case 28110 : return &A4q1_2q228110_eval;
case 28112 : return &A4q1_2q228112_eval;
case 28157 : return &A4q1_2q228157_eval;
case 28169 : return &A4q1_2q228169_eval;
case 28182 : return &A4q1_2q228182_eval;
case 28184 : return &A4q1_2q228184_eval;
case 28260 : return &A4q1_2q228260_eval;
case 28262 : return &A4q1_2q228262_eval;
case 28272 : return &A4q1_2q228272_eval;
case 28274 : return &A4q1_2q228274_eval;
case 28300 : return &A4q1_2q228300_eval;
case 28312 : return &A4q1_2q228312_eval;
case 28320 : return &A4q1_2q228320_eval;
case 28322 : return &A4q1_2q228322_eval;
case 28372 : return &A4q1_2q228372_eval;
case 28384 : return &A4q1_2q228384_eval;
case 28392 : return &A4q1_2q228392_eval;
case 28394 : return &A4q1_2q228394_eval;
case 28440 : return &A4q1_2q228440_eval;
case 28442 : return &A4q1_2q228442_eval;
case 28452 : return &A4q1_2q228452_eval;
case 28454 : return &A4q1_2q228454_eval;
case 28523 : return &A4q1_2q228523_eval;
case 28535 : return &A4q1_2q228535_eval;
case 28543 : return &A4q1_2q228543_eval;
case 28545 : return &A4q1_2q228545_eval;
case 28553 : return &A4q1_2q228553_eval;
case 28565 : return &A4q1_2q228565_eval;
case 28578 : return &A4q1_2q228578_eval;
case 28580 : return &A4q1_2q228580_eval;
case 28595 : return &A4q1_2q228595_eval;
case 28607 : return &A4q1_2q228607_eval;
case 28615 : return &A4q1_2q228615_eval;
case 28617 : return &A4q1_2q228617_eval;
case 28625 : return &A4q1_2q228625_eval;
case 28637 : return &A4q1_2q228637_eval;
case 28650 : return &A4q1_2q228650_eval;
case 28652 : return &A4q1_2q228652_eval;
case 28693 : return &A4q1_2q228693_eval;
case 28695 : return &A4q1_2q228695_eval;
case 28698 : return &A4q1_2q228698_eval;
case 28700 : return &A4q1_2q228700_eval;
case 28705 : return &A4q1_2q228705_eval;
case 28707 : return &A4q1_2q228707_eval;
case 28710 : return &A4q1_2q228710_eval;
case 28712 : return &A4q1_2q228712_eval;
case 28733 : return &A4q1_2q228733_eval;
case 28745 : return &A4q1_2q228745_eval;
case 28758 : return &A4q1_2q228758_eval;
case 28760 : return &A4q1_2q228760_eval;
case 28805 : return &A4q1_2q228805_eval;
case 28817 : return &A4q1_2q228817_eval;
case 28830 : return &A4q1_2q228830_eval;
case 28832 : return &A4q1_2q228832_eval;
case 28908 : return &A4q1_2q228908_eval;
case 28910 : return &A4q1_2q228910_eval;
case 28920 : return &A4q1_2q228920_eval;
case 28922 : return &A4q1_2q228922_eval;
case 28955 : return &A4q1_2q228955_eval;
case 28967 : return &A4q1_2q228967_eval;
case 28975 : return &A4q1_2q228975_eval;
case 28977 : return &A4q1_2q228977_eval;
case 28985 : return &A4q1_2q228985_eval;
case 28997 : return &A4q1_2q228997_eval;
case 29010 : return &A4q1_2q229010_eval;
case 29012 : return &A4q1_2q229012_eval;
case 29027 : return &A4q1_2q229027_eval;
case 29039 : return &A4q1_2q229039_eval;
case 29047 : return &A4q1_2q229047_eval;
case 29049 : return &A4q1_2q229049_eval;
case 29057 : return &A4q1_2q229057_eval;
case 29069 : return &A4q1_2q229069_eval;
case 29082 : return &A4q1_2q229082_eval;
case 29084 : return &A4q1_2q229084_eval;
case 29125 : return &A4q1_2q229125_eval;
case 29127 : return &A4q1_2q229127_eval;
case 29130 : return &A4q1_2q229130_eval;
case 29132 : return &A4q1_2q229132_eval;
case 29137 : return &A4q1_2q229137_eval;
case 29139 : return &A4q1_2q229139_eval;
case 29142 : return &A4q1_2q229142_eval;
case 29144 : return &A4q1_2q229144_eval;
case 29165 : return &A4q1_2q229165_eval;
case 29177 : return &A4q1_2q229177_eval;
case 29190 : return &A4q1_2q229190_eval;
case 29192 : return &A4q1_2q229192_eval;
case 29237 : return &A4q1_2q229237_eval;
case 29249 : return &A4q1_2q229249_eval;
case 29262 : return &A4q1_2q229262_eval;
case 29264 : return &A4q1_2q229264_eval;
case 29340 : return &A4q1_2q229340_eval;
case 29342 : return &A4q1_2q229342_eval;
case 29352 : return &A4q1_2q229352_eval;
case 29354 : return &A4q1_2q229354_eval;
case 29593 : return &A4q1_2q229593_eval;
case 29595 : return &A4q1_2q229595_eval;
case 29598 : return &A4q1_2q229598_eval;
case 29600 : return &A4q1_2q229600_eval;
case 29605 : return &A4q1_2q229605_eval;
case 29607 : return &A4q1_2q229607_eval;
case 29610 : return &A4q1_2q229610_eval;
case 29612 : return &A4q1_2q229612_eval;
case 29628 : return &A4q1_2q229628_eval;
case 29630 : return &A4q1_2q229630_eval;
case 29640 : return &A4q1_2q229640_eval;
case 29642 : return &A4q1_2q229642_eval;
case 29665 : return &A4q1_2q229665_eval;
case 29667 : return &A4q1_2q229667_eval;
case 29670 : return &A4q1_2q229670_eval;
case 29672 : return &A4q1_2q229672_eval;
case 29677 : return &A4q1_2q229677_eval;
case 29679 : return &A4q1_2q229679_eval;
case 29682 : return &A4q1_2q229682_eval;
case 29684 : return &A4q1_2q229684_eval;
case 29700 : return &A4q1_2q229700_eval;
case 29702 : return &A4q1_2q229702_eval;
case 29712 : return &A4q1_2q229712_eval;
case 29714 : return &A4q1_2q229714_eval;
case 29818 : return &A4q1_2q229818_eval;
case 29830 : return &A4q1_2q229830_eval;
case 29833 : return &A4q1_2q229833_eval;
case 29835 : return &A4q1_2q229835_eval;
case 29848 : return &A4q1_2q229848_eval;
case 29860 : return &A4q1_2q229860_eval;
case 29868 : return &A4q1_2q229868_eval;
case 29870 : return &A4q1_2q229870_eval;
case 29890 : return &A4q1_2q229890_eval;
case 29902 : return &A4q1_2q229902_eval;
case 29905 : return &A4q1_2q229905_eval;
case 29907 : return &A4q1_2q229907_eval;
case 29920 : return &A4q1_2q229920_eval;
case 29932 : return &A4q1_2q229932_eval;
case 29940 : return &A4q1_2q229940_eval;
case 29942 : return &A4q1_2q229942_eval;
case 29953 : return &A4q1_2q229953_eval;
case 29955 : return &A4q1_2q229955_eval;
case 29958 : return &A4q1_2q229958_eval;
case 29960 : return &A4q1_2q229960_eval;
case 29965 : return &A4q1_2q229965_eval;
case 29967 : return &A4q1_2q229967_eval;
case 29970 : return &A4q1_2q229970_eval;
case 29972 : return &A4q1_2q229972_eval;
case 30028 : return &A4q1_2q230028_eval;
case 30040 : return &A4q1_2q230040_eval;
case 30048 : return &A4q1_2q230048_eval;
case 30050 : return &A4q1_2q230050_eval;
case 30100 : return &A4q1_2q230100_eval;
case 30112 : return &A4q1_2q230112_eval;
case 30120 : return &A4q1_2q230120_eval;
case 30122 : return &A4q1_2q230122_eval;
case 30168 : return &A4q1_2q230168_eval;
case 30170 : return &A4q1_2q230170_eval;
case 30180 : return &A4q1_2q230180_eval;
case 30182 : return &A4q1_2q230182_eval;
case 30250 : return &A4q1_2q230250_eval;
case 30262 : return &A4q1_2q230262_eval;
case 30265 : return &A4q1_2q230265_eval;
case 30267 : return &A4q1_2q230267_eval;
case 30280 : return &A4q1_2q230280_eval;
case 30292 : return &A4q1_2q230292_eval;
case 30300 : return &A4q1_2q230300_eval;
case 30302 : return &A4q1_2q230302_eval;
case 30322 : return &A4q1_2q230322_eval;
case 30334 : return &A4q1_2q230334_eval;
case 30337 : return &A4q1_2q230337_eval;
case 30339 : return &A4q1_2q230339_eval;
case 30352 : return &A4q1_2q230352_eval;
case 30364 : return &A4q1_2q230364_eval;
case 30372 : return &A4q1_2q230372_eval;
case 30374 : return &A4q1_2q230374_eval;
case 30385 : return &A4q1_2q230385_eval;
case 30387 : return &A4q1_2q230387_eval;
case 30390 : return &A4q1_2q230390_eval;
case 30392 : return &A4q1_2q230392_eval;
case 30397 : return &A4q1_2q230397_eval;
case 30399 : return &A4q1_2q230399_eval;
case 30402 : return &A4q1_2q230402_eval;
case 30404 : return &A4q1_2q230404_eval;
case 30460 : return &A4q1_2q230460_eval;
case 30472 : return &A4q1_2q230472_eval;
case 30480 : return &A4q1_2q230480_eval;
case 30482 : return &A4q1_2q230482_eval;
case 30532 : return &A4q1_2q230532_eval;
case 30544 : return &A4q1_2q230544_eval;
case 30552 : return &A4q1_2q230552_eval;
case 30554 : return &A4q1_2q230554_eval;
case 30600 : return &A4q1_2q230600_eval;
case 30602 : return &A4q1_2q230602_eval;
case 30612 : return &A4q1_2q230612_eval;
case 30614 : return &A4q1_2q230614_eval;
case 30673 : return &A4q1_2q230673_eval;
case 30675 : return &A4q1_2q230675_eval;
case 30678 : return &A4q1_2q230678_eval;
case 30680 : return &A4q1_2q230680_eval;
case 30685 : return &A4q1_2q230685_eval;
case 30687 : return &A4q1_2q230687_eval;
case 30690 : return &A4q1_2q230690_eval;
case 30692 : return &A4q1_2q230692_eval;
case 30708 : return &A4q1_2q230708_eval;
case 30710 : return &A4q1_2q230710_eval;
case 30720 : return &A4q1_2q230720_eval;
case 30722 : return &A4q1_2q230722_eval;
case 30745 : return &A4q1_2q230745_eval;
case 30747 : return &A4q1_2q230747_eval;
case 30750 : return &A4q1_2q230750_eval;
case 30752 : return &A4q1_2q230752_eval;
case 30757 : return &A4q1_2q230757_eval;
case 30759 : return &A4q1_2q230759_eval;
case 30762 : return &A4q1_2q230762_eval;
case 30764 : return &A4q1_2q230764_eval;
case 30780 : return &A4q1_2q230780_eval;
case 30782 : return &A4q1_2q230782_eval;
case 30792 : return &A4q1_2q230792_eval;
case 30794 : return &A4q1_2q230794_eval;
case 31151 : return &A4q1_2q231151_eval;
case 31163 : return &A4q1_2q231163_eval;
case 31171 : return &A4q1_2q231171_eval;
case 31173 : return &A4q1_2q231173_eval;
case 31223 : return &A4q1_2q231223_eval;
case 31235 : return &A4q1_2q231235_eval;
case 31243 : return &A4q1_2q231243_eval;
case 31245 : return &A4q1_2q231245_eval;
case 31291 : return &A4q1_2q231291_eval;
case 31293 : return &A4q1_2q231293_eval;
case 31303 : return &A4q1_2q231303_eval;
case 31305 : return &A4q1_2q231305_eval;
case 31331 : return &A4q1_2q231331_eval;
case 31343 : return &A4q1_2q231343_eval;
case 31351 : return &A4q1_2q231351_eval;
case 31353 : return &A4q1_2q231353_eval;
case 31361 : return &A4q1_2q231361_eval;
case 31373 : return &A4q1_2q231373_eval;
case 31386 : return &A4q1_2q231386_eval;
case 31388 : return &A4q1_2q231388_eval;
case 31403 : return &A4q1_2q231403_eval;
case 31415 : return &A4q1_2q231415_eval;
case 31423 : return &A4q1_2q231423_eval;
case 31425 : return &A4q1_2q231425_eval;
case 31433 : return &A4q1_2q231433_eval;
case 31445 : return &A4q1_2q231445_eval;
case 31458 : return &A4q1_2q231458_eval;
case 31460 : return &A4q1_2q231460_eval;
case 31501 : return &A4q1_2q231501_eval;
case 31503 : return &A4q1_2q231503_eval;
case 31506 : return &A4q1_2q231506_eval;
case 31508 : return &A4q1_2q231508_eval;
case 31513 : return &A4q1_2q231513_eval;
case 31515 : return &A4q1_2q231515_eval;
case 31518 : return &A4q1_2q231518_eval;
case 31520 : return &A4q1_2q231520_eval;
case 31583 : return &A4q1_2q231583_eval;
case 31595 : return &A4q1_2q231595_eval;
case 31603 : return &A4q1_2q231603_eval;
case 31605 : return &A4q1_2q231605_eval;
case 31655 : return &A4q1_2q231655_eval;
case 31667 : return &A4q1_2q231667_eval;
case 31675 : return &A4q1_2q231675_eval;
case 31677 : return &A4q1_2q231677_eval;
case 31723 : return &A4q1_2q231723_eval;
case 31725 : return &A4q1_2q231725_eval;
case 31735 : return &A4q1_2q231735_eval;
case 31737 : return &A4q1_2q231737_eval;
case 31763 : return &A4q1_2q231763_eval;
case 31775 : return &A4q1_2q231775_eval;
case 31783 : return &A4q1_2q231783_eval;
case 31785 : return &A4q1_2q231785_eval;
case 31793 : return &A4q1_2q231793_eval;
case 31805 : return &A4q1_2q231805_eval;
case 31818 : return &A4q1_2q231818_eval;
case 31820 : return &A4q1_2q231820_eval;
case 31835 : return &A4q1_2q231835_eval;
case 31847 : return &A4q1_2q231847_eval;
case 31855 : return &A4q1_2q231855_eval;
case 31857 : return &A4q1_2q231857_eval;
case 31865 : return &A4q1_2q231865_eval;
case 31877 : return &A4q1_2q231877_eval;
case 31890 : return &A4q1_2q231890_eval;
case 31892 : return &A4q1_2q231892_eval;
case 31933 : return &A4q1_2q231933_eval;
case 31935 : return &A4q1_2q231935_eval;
case 31938 : return &A4q1_2q231938_eval;
case 31940 : return &A4q1_2q231940_eval;
case 31945 : return &A4q1_2q231945_eval;
case 31947 : return &A4q1_2q231947_eval;
case 31950 : return &A4q1_2q231950_eval;
case 31952 : return &A4q1_2q231952_eval;
case 32191 : return &A4q1_2q232191_eval;
case 32193 : return &A4q1_2q232193_eval;
case 32203 : return &A4q1_2q232203_eval;
case 32205 : return &A4q1_2q232205_eval;
case 32221 : return &A4q1_2q232221_eval;
case 32223 : return &A4q1_2q232223_eval;
case 32226 : return &A4q1_2q232226_eval;
case 32228 : return &A4q1_2q232228_eval;
case 32233 : return &A4q1_2q232233_eval;
case 32235 : return &A4q1_2q232235_eval;
case 32238 : return &A4q1_2q232238_eval;
case 32240 : return &A4q1_2q232240_eval;
case 32263 : return &A4q1_2q232263_eval;
case 32265 : return &A4q1_2q232265_eval;
case 32275 : return &A4q1_2q232275_eval;
case 32277 : return &A4q1_2q232277_eval;
case 32293 : return &A4q1_2q232293_eval;
case 32295 : return &A4q1_2q232295_eval;
case 32298 : return &A4q1_2q232298_eval;
case 32300 : return &A4q1_2q232300_eval;
case 32305 : return &A4q1_2q232305_eval;
case 32307 : return &A4q1_2q232307_eval;
case 32310 : return &A4q1_2q232310_eval;
case 32312 : return &A4q1_2q232312_eval;
case 32411 : return &A4q1_2q232411_eval;
case 32423 : return &A4q1_2q232423_eval;
case 32431 : return &A4q1_2q232431_eval;
case 32433 : return &A4q1_2q232433_eval;
case 32441 : return &A4q1_2q232441_eval;
case 32453 : return &A4q1_2q232453_eval;
case 32466 : return &A4q1_2q232466_eval;
case 32468 : return &A4q1_2q232468_eval;
case 32483 : return &A4q1_2q232483_eval;
case 32495 : return &A4q1_2q232495_eval;
case 32503 : return &A4q1_2q232503_eval;
case 32505 : return &A4q1_2q232505_eval;
case 32513 : return &A4q1_2q232513_eval;
case 32525 : return &A4q1_2q232525_eval;
case 32538 : return &A4q1_2q232538_eval;
case 32540 : return &A4q1_2q232540_eval;
case 32581 : return &A4q1_2q232581_eval;
case 32583 : return &A4q1_2q232583_eval;
case 32586 : return &A4q1_2q232586_eval;
case 32588 : return &A4q1_2q232588_eval;
case 32593 : return &A4q1_2q232593_eval;
case 32595 : return &A4q1_2q232595_eval;
case 32598 : return &A4q1_2q232598_eval;
case 32600 : return &A4q1_2q232600_eval;
case 32621 : return &A4q1_2q232621_eval;
case 32633 : return &A4q1_2q232633_eval;
case 32646 : return &A4q1_2q232646_eval;
case 32648 : return &A4q1_2q232648_eval;
case 32693 : return &A4q1_2q232693_eval;
case 32705 : return &A4q1_2q232705_eval;
case 32718 : return &A4q1_2q232718_eval;
case 32720 : return &A4q1_2q232720_eval;
case 32796 : return &A4q1_2q232796_eval;
case 32798 : return &A4q1_2q232798_eval;
case 32808 : return &A4q1_2q232808_eval;
case 32810 : return &A4q1_2q232810_eval;
case 32843 : return &A4q1_2q232843_eval;
case 32855 : return &A4q1_2q232855_eval;
case 32863 : return &A4q1_2q232863_eval;
case 32865 : return &A4q1_2q232865_eval;
case 32873 : return &A4q1_2q232873_eval;
case 32885 : return &A4q1_2q232885_eval;
case 32898 : return &A4q1_2q232898_eval;
case 32900 : return &A4q1_2q232900_eval;
case 32915 : return &A4q1_2q232915_eval;
case 32927 : return &A4q1_2q232927_eval;
case 32935 : return &A4q1_2q232935_eval;
case 32937 : return &A4q1_2q232937_eval;
case 32945 : return &A4q1_2q232945_eval;
case 32957 : return &A4q1_2q232957_eval;
case 32970 : return &A4q1_2q232970_eval;
case 32972 : return &A4q1_2q232972_eval;
case 33013 : return &A4q1_2q233013_eval;
case 33015 : return &A4q1_2q233015_eval;
case 33018 : return &A4q1_2q233018_eval;
case 33020 : return &A4q1_2q233020_eval;
case 33025 : return &A4q1_2q233025_eval;
case 33027 : return &A4q1_2q233027_eval;
case 33030 : return &A4q1_2q233030_eval;
case 33032 : return &A4q1_2q233032_eval;
case 33053 : return &A4q1_2q233053_eval;
case 33065 : return &A4q1_2q233065_eval;
case 33078 : return &A4q1_2q233078_eval;
case 33080 : return &A4q1_2q233080_eval;
case 33125 : return &A4q1_2q233125_eval;
case 33137 : return &A4q1_2q233137_eval;
case 33150 : return &A4q1_2q233150_eval;
case 33152 : return &A4q1_2q233152_eval;
case 33228 : return &A4q1_2q233228_eval;
case 33230 : return &A4q1_2q233230_eval;
case 33240 : return &A4q1_2q233240_eval;
case 33242 : return &A4q1_2q233242_eval;
case 33481 : return &A4q1_2q233481_eval;
case 33483 : return &A4q1_2q233483_eval;
case 33486 : return &A4q1_2q233486_eval;
case 33488 : return &A4q1_2q233488_eval;
case 33493 : return &A4q1_2q233493_eval;
case 33495 : return &A4q1_2q233495_eval;
case 33498 : return &A4q1_2q233498_eval;
case 33500 : return &A4q1_2q233500_eval;
case 33516 : return &A4q1_2q233516_eval;
case 33518 : return &A4q1_2q233518_eval;
case 33528 : return &A4q1_2q233528_eval;
case 33530 : return &A4q1_2q233530_eval;
case 33553 : return &A4q1_2q233553_eval;
case 33555 : return &A4q1_2q233555_eval;
case 33558 : return &A4q1_2q233558_eval;
case 33560 : return &A4q1_2q233560_eval;
case 33565 : return &A4q1_2q233565_eval;
case 33567 : return &A4q1_2q233567_eval;
case 33570 : return &A4q1_2q233570_eval;
case 33572 : return &A4q1_2q233572_eval;
case 33588 : return &A4q1_2q233588_eval;
case 33590 : return &A4q1_2q233590_eval;
case 33600 : return &A4q1_2q233600_eval;
case 33602 : return &A4q1_2q233602_eval;
case 33743 : return &A4q1_2q233743_eval;
case 33755 : return &A4q1_2q233755_eval;
case 33763 : return &A4q1_2q233763_eval;
case 33765 : return &A4q1_2q233765_eval;
case 33815 : return &A4q1_2q233815_eval;
case 33827 : return &A4q1_2q233827_eval;
case 33835 : return &A4q1_2q233835_eval;
case 33837 : return &A4q1_2q233837_eval;
case 33883 : return &A4q1_2q233883_eval;
case 33885 : return &A4q1_2q233885_eval;
case 33895 : return &A4q1_2q233895_eval;
case 33897 : return &A4q1_2q233897_eval;
case 33923 : return &A4q1_2q233923_eval;
case 33935 : return &A4q1_2q233935_eval;
case 33943 : return &A4q1_2q233943_eval;
case 33945 : return &A4q1_2q233945_eval;
case 33953 : return &A4q1_2q233953_eval;
case 33965 : return &A4q1_2q233965_eval;
case 33978 : return &A4q1_2q233978_eval;
case 33980 : return &A4q1_2q233980_eval;
case 33995 : return &A4q1_2q233995_eval;
case 34007 : return &A4q1_2q234007_eval;
case 34015 : return &A4q1_2q234015_eval;
case 34017 : return &A4q1_2q234017_eval;
case 34025 : return &A4q1_2q234025_eval;
case 34037 : return &A4q1_2q234037_eval;
case 34050 : return &A4q1_2q234050_eval;
case 34052 : return &A4q1_2q234052_eval;
case 34093 : return &A4q1_2q234093_eval;
case 34095 : return &A4q1_2q234095_eval;
case 34098 : return &A4q1_2q234098_eval;
case 34100 : return &A4q1_2q234100_eval;
case 34105 : return &A4q1_2q234105_eval;
case 34107 : return &A4q1_2q234107_eval;
case 34110 : return &A4q1_2q234110_eval;
case 34112 : return &A4q1_2q234112_eval;
case 34175 : return &A4q1_2q234175_eval;
case 34187 : return &A4q1_2q234187_eval;
case 34195 : return &A4q1_2q234195_eval;
case 34197 : return &A4q1_2q234197_eval;
case 34247 : return &A4q1_2q234247_eval;
case 34259 : return &A4q1_2q234259_eval;
case 34267 : return &A4q1_2q234267_eval;
case 34269 : return &A4q1_2q234269_eval;
case 34315 : return &A4q1_2q234315_eval;
case 34317 : return &A4q1_2q234317_eval;
case 34327 : return &A4q1_2q234327_eval;
case 34329 : return &A4q1_2q234329_eval;
case 34355 : return &A4q1_2q234355_eval;
case 34367 : return &A4q1_2q234367_eval;
case 34375 : return &A4q1_2q234375_eval;
case 34377 : return &A4q1_2q234377_eval;
case 34385 : return &A4q1_2q234385_eval;
case 34397 : return &A4q1_2q234397_eval;
case 34410 : return &A4q1_2q234410_eval;
case 34412 : return &A4q1_2q234412_eval;
case 34427 : return &A4q1_2q234427_eval;
case 34439 : return &A4q1_2q234439_eval;
case 34447 : return &A4q1_2q234447_eval;
case 34449 : return &A4q1_2q234449_eval;
case 34457 : return &A4q1_2q234457_eval;
case 34469 : return &A4q1_2q234469_eval;
case 34482 : return &A4q1_2q234482_eval;
case 34484 : return &A4q1_2q234484_eval;
case 34525 : return &A4q1_2q234525_eval;
case 34527 : return &A4q1_2q234527_eval;
case 34530 : return &A4q1_2q234530_eval;
case 34532 : return &A4q1_2q234532_eval;
case 34537 : return &A4q1_2q234537_eval;
case 34539 : return &A4q1_2q234539_eval;
case 34542 : return &A4q1_2q234542_eval;
case 34544 : return &A4q1_2q234544_eval;
case 34783 : return &A4q1_2q234783_eval;
case 34785 : return &A4q1_2q234785_eval;
case 34795 : return &A4q1_2q234795_eval;
case 34797 : return &A4q1_2q234797_eval;
case 34813 : return &A4q1_2q234813_eval;
case 34815 : return &A4q1_2q234815_eval;
case 34818 : return &A4q1_2q234818_eval;
case 34820 : return &A4q1_2q234820_eval;
case 34825 : return &A4q1_2q234825_eval;
case 34827 : return &A4q1_2q234827_eval;
case 34830 : return &A4q1_2q234830_eval;
case 34832 : return &A4q1_2q234832_eval;
case 34855 : return &A4q1_2q234855_eval;
case 34857 : return &A4q1_2q234857_eval;
case 34867 : return &A4q1_2q234867_eval;
case 34869 : return &A4q1_2q234869_eval;
case 34885 : return &A4q1_2q234885_eval;
case 34887 : return &A4q1_2q234887_eval;
case 34890 : return &A4q1_2q234890_eval;
case 34892 : return &A4q1_2q234892_eval;
case 34897 : return &A4q1_2q234897_eval;
case 34899 : return &A4q1_2q234899_eval;
case 34902 : return &A4q1_2q234902_eval;
case 34904 : return &A4q1_2q234904_eval;
case 35003 : return &A4q1_2q235003_eval;
case 35015 : return &A4q1_2q235015_eval;
case 35023 : return &A4q1_2q235023_eval;
case 35025 : return &A4q1_2q235025_eval;
case 35033 : return &A4q1_2q235033_eval;
case 35045 : return &A4q1_2q235045_eval;
case 35058 : return &A4q1_2q235058_eval;
case 35060 : return &A4q1_2q235060_eval;
case 35075 : return &A4q1_2q235075_eval;
case 35087 : return &A4q1_2q235087_eval;
case 35095 : return &A4q1_2q235095_eval;
case 35097 : return &A4q1_2q235097_eval;
case 35105 : return &A4q1_2q235105_eval;
case 35117 : return &A4q1_2q235117_eval;
case 35130 : return &A4q1_2q235130_eval;
case 35132 : return &A4q1_2q235132_eval;
case 35173 : return &A4q1_2q235173_eval;
case 35175 : return &A4q1_2q235175_eval;
case 35178 : return &A4q1_2q235178_eval;
case 35180 : return &A4q1_2q235180_eval;
case 35185 : return &A4q1_2q235185_eval;
case 35187 : return &A4q1_2q235187_eval;
case 35190 : return &A4q1_2q235190_eval;
case 35192 : return &A4q1_2q235192_eval;
case 35213 : return &A4q1_2q235213_eval;
case 35225 : return &A4q1_2q235225_eval;
case 35238 : return &A4q1_2q235238_eval;
case 35240 : return &A4q1_2q235240_eval;
case 35285 : return &A4q1_2q235285_eval;
case 35297 : return &A4q1_2q235297_eval;
case 35310 : return &A4q1_2q235310_eval;
case 35312 : return &A4q1_2q235312_eval;
case 35388 : return &A4q1_2q235388_eval;
case 35390 : return &A4q1_2q235390_eval;
case 35400 : return &A4q1_2q235400_eval;
case 35402 : return &A4q1_2q235402_eval;
case 35435 : return &A4q1_2q235435_eval;
case 35447 : return &A4q1_2q235447_eval;
case 35455 : return &A4q1_2q235455_eval;
case 35457 : return &A4q1_2q235457_eval;
case 35465 : return &A4q1_2q235465_eval;
case 35477 : return &A4q1_2q235477_eval;
case 35490 : return &A4q1_2q235490_eval;
case 35492 : return &A4q1_2q235492_eval;
case 35507 : return &A4q1_2q235507_eval;
case 35519 : return &A4q1_2q235519_eval;
case 35527 : return &A4q1_2q235527_eval;
case 35529 : return &A4q1_2q235529_eval;
case 35537 : return &A4q1_2q235537_eval;
case 35549 : return &A4q1_2q235549_eval;
case 35562 : return &A4q1_2q235562_eval;
case 35564 : return &A4q1_2q235564_eval;
case 35605 : return &A4q1_2q235605_eval;
case 35607 : return &A4q1_2q235607_eval;
case 35610 : return &A4q1_2q235610_eval;
case 35612 : return &A4q1_2q235612_eval;
case 35617 : return &A4q1_2q235617_eval;
case 35619 : return &A4q1_2q235619_eval;
case 35622 : return &A4q1_2q235622_eval;
case 35624 : return &A4q1_2q235624_eval;
case 35645 : return &A4q1_2q235645_eval;
case 35657 : return &A4q1_2q235657_eval;
case 35670 : return &A4q1_2q235670_eval;
case 35672 : return &A4q1_2q235672_eval;
case 35717 : return &A4q1_2q235717_eval;
case 35729 : return &A4q1_2q235729_eval;
case 35742 : return &A4q1_2q235742_eval;
case 35744 : return &A4q1_2q235744_eval;
case 35820 : return &A4q1_2q235820_eval;
case 35822 : return &A4q1_2q235822_eval;
case 35832 : return &A4q1_2q235832_eval;
case 35834 : return &A4q1_2q235834_eval;
case 36073 : return &A4q1_2q236073_eval;
case 36075 : return &A4q1_2q236075_eval;
case 36078 : return &A4q1_2q236078_eval;
case 36080 : return &A4q1_2q236080_eval;
case 36085 : return &A4q1_2q236085_eval;
case 36087 : return &A4q1_2q236087_eval;
case 36090 : return &A4q1_2q236090_eval;
case 36092 : return &A4q1_2q236092_eval;
case 36108 : return &A4q1_2q236108_eval;
case 36110 : return &A4q1_2q236110_eval;
case 36120 : return &A4q1_2q236120_eval;
case 36122 : return &A4q1_2q236122_eval;
case 36145 : return &A4q1_2q236145_eval;
case 36147 : return &A4q1_2q236147_eval;
case 36150 : return &A4q1_2q236150_eval;
case 36152 : return &A4q1_2q236152_eval;
case 36157 : return &A4q1_2q236157_eval;
case 36159 : return &A4q1_2q236159_eval;
case 36162 : return &A4q1_2q236162_eval;
case 36164 : return &A4q1_2q236164_eval;
case 36180 : return &A4q1_2q236180_eval;
case 36182 : return &A4q1_2q236182_eval;
case 36192 : return &A4q1_2q236192_eval;
case 36194 : return &A4q1_2q236194_eval;
case 37591 : return &A4q1_2q237591_eval;
case 37593 : return &A4q1_2q237593_eval;
case 37603 : return &A4q1_2q237603_eval;
case 37605 : return &A4q1_2q237605_eval;
case 37621 : return &A4q1_2q237621_eval;
case 37623 : return &A4q1_2q237623_eval;
case 37626 : return &A4q1_2q237626_eval;
case 37628 : return &A4q1_2q237628_eval;
case 37633 : return &A4q1_2q237633_eval;
case 37635 : return &A4q1_2q237635_eval;
case 37638 : return &A4q1_2q237638_eval;
case 37640 : return &A4q1_2q237640_eval;
case 37663 : return &A4q1_2q237663_eval;
case 37665 : return &A4q1_2q237665_eval;
case 37675 : return &A4q1_2q237675_eval;
case 37677 : return &A4q1_2q237677_eval;
case 37693 : return &A4q1_2q237693_eval;
case 37695 : return &A4q1_2q237695_eval;
case 37698 : return &A4q1_2q237698_eval;
case 37700 : return &A4q1_2q237700_eval;
case 37705 : return &A4q1_2q237705_eval;
case 37707 : return &A4q1_2q237707_eval;
case 37710 : return &A4q1_2q237710_eval;
case 37712 : return &A4q1_2q237712_eval;
case 37801 : return &A4q1_2q237801_eval;
case 37803 : return &A4q1_2q237803_eval;
case 37806 : return &A4q1_2q237806_eval;
case 37808 : return &A4q1_2q237808_eval;
case 37813 : return &A4q1_2q237813_eval;
case 37815 : return &A4q1_2q237815_eval;
case 37818 : return &A4q1_2q237818_eval;
case 37820 : return &A4q1_2q237820_eval;
case 37836 : return &A4q1_2q237836_eval;
case 37838 : return &A4q1_2q237838_eval;
case 37848 : return &A4q1_2q237848_eval;
case 37850 : return &A4q1_2q237850_eval;
case 37873 : return &A4q1_2q237873_eval;
case 37875 : return &A4q1_2q237875_eval;
case 37878 : return &A4q1_2q237878_eval;
case 37880 : return &A4q1_2q237880_eval;
case 37885 : return &A4q1_2q237885_eval;
case 37887 : return &A4q1_2q237887_eval;
case 37890 : return &A4q1_2q237890_eval;
case 37892 : return &A4q1_2q237892_eval;
case 37908 : return &A4q1_2q237908_eval;
case 37910 : return &A4q1_2q237910_eval;
case 37920 : return &A4q1_2q237920_eval;
case 37922 : return &A4q1_2q237922_eval;
case 38023 : return &A4q1_2q238023_eval;
case 38025 : return &A4q1_2q238025_eval;
case 38035 : return &A4q1_2q238035_eval;
case 38037 : return &A4q1_2q238037_eval;
case 38053 : return &A4q1_2q238053_eval;
case 38055 : return &A4q1_2q238055_eval;
case 38058 : return &A4q1_2q238058_eval;
case 38060 : return &A4q1_2q238060_eval;
case 38065 : return &A4q1_2q238065_eval;
case 38067 : return &A4q1_2q238067_eval;
case 38070 : return &A4q1_2q238070_eval;
case 38072 : return &A4q1_2q238072_eval;
case 38095 : return &A4q1_2q238095_eval;
case 38097 : return &A4q1_2q238097_eval;
case 38107 : return &A4q1_2q238107_eval;
case 38109 : return &A4q1_2q238109_eval;
case 38125 : return &A4q1_2q238125_eval;
case 38127 : return &A4q1_2q238127_eval;
case 38130 : return &A4q1_2q238130_eval;
case 38132 : return &A4q1_2q238132_eval;
case 38137 : return &A4q1_2q238137_eval;
case 38139 : return &A4q1_2q238139_eval;
case 38142 : return &A4q1_2q238142_eval;
case 38144 : return &A4q1_2q238144_eval;
case 38233 : return &A4q1_2q238233_eval;
case 38235 : return &A4q1_2q238235_eval;
case 38238 : return &A4q1_2q238238_eval;
case 38240 : return &A4q1_2q238240_eval;
case 38245 : return &A4q1_2q238245_eval;
case 38247 : return &A4q1_2q238247_eval;
case 38250 : return &A4q1_2q238250_eval;
case 38252 : return &A4q1_2q238252_eval;
case 38268 : return &A4q1_2q238268_eval;
case 38270 : return &A4q1_2q238270_eval;
case 38280 : return &A4q1_2q238280_eval;
case 38282 : return &A4q1_2q238282_eval;
case 38305 : return &A4q1_2q238305_eval;
case 38307 : return &A4q1_2q238307_eval;
case 38310 : return &A4q1_2q238310_eval;
case 38312 : return &A4q1_2q238312_eval;
case 38317 : return &A4q1_2q238317_eval;
case 38319 : return &A4q1_2q238319_eval;
case 38322 : return &A4q1_2q238322_eval;
case 38324 : return &A4q1_2q238324_eval;
case 38340 : return &A4q1_2q238340_eval;
case 38342 : return &A4q1_2q238342_eval;
case 38352 : return &A4q1_2q238352_eval;
case 38354 : return &A4q1_2q238354_eval;
case 38926 : return &A4q1_2q238926_eval;
case 38938 : return &A4q1_2q238938_eval;
case 38941 : return &A4q1_2q238941_eval;
case 38943 : return &A4q1_2q238943_eval;
case 38998 : return &A4q1_2q238998_eval;
case 39010 : return &A4q1_2q239010_eval;
case 39013 : return &A4q1_2q239013_eval;
case 39015 : return &A4q1_2q239015_eval;
case 39031 : return &A4q1_2q239031_eval;
case 39033 : return &A4q1_2q239033_eval;
case 39043 : return &A4q1_2q239043_eval;
case 39045 : return &A4q1_2q239045_eval;
case 39106 : return &A4q1_2q239106_eval;
case 39118 : return &A4q1_2q239118_eval;
case 39121 : return &A4q1_2q239121_eval;
case 39123 : return &A4q1_2q239123_eval;
case 39136 : return &A4q1_2q239136_eval;
case 39148 : return &A4q1_2q239148_eval;
case 39156 : return &A4q1_2q239156_eval;
case 39158 : return &A4q1_2q239158_eval;
case 39178 : return &A4q1_2q239178_eval;
case 39190 : return &A4q1_2q239190_eval;
case 39193 : return &A4q1_2q239193_eval;
case 39195 : return &A4q1_2q239195_eval;
case 39208 : return &A4q1_2q239208_eval;
case 39220 : return &A4q1_2q239220_eval;
case 39228 : return &A4q1_2q239228_eval;
case 39230 : return &A4q1_2q239230_eval;
case 39241 : return &A4q1_2q239241_eval;
case 39243 : return &A4q1_2q239243_eval;
case 39246 : return &A4q1_2q239246_eval;
case 39248 : return &A4q1_2q239248_eval;
case 39253 : return &A4q1_2q239253_eval;
case 39255 : return &A4q1_2q239255_eval;
case 39258 : return &A4q1_2q239258_eval;
case 39260 : return &A4q1_2q239260_eval;
case 39358 : return &A4q1_2q239358_eval;
case 39370 : return &A4q1_2q239370_eval;
case 39373 : return &A4q1_2q239373_eval;
case 39375 : return &A4q1_2q239375_eval;
case 39430 : return &A4q1_2q239430_eval;
case 39442 : return &A4q1_2q239442_eval;
case 39445 : return &A4q1_2q239445_eval;
case 39447 : return &A4q1_2q239447_eval;
case 39463 : return &A4q1_2q239463_eval;
case 39465 : return &A4q1_2q239465_eval;
case 39475 : return &A4q1_2q239475_eval;
case 39477 : return &A4q1_2q239477_eval;
case 39538 : return &A4q1_2q239538_eval;
case 39550 : return &A4q1_2q239550_eval;
case 39553 : return &A4q1_2q239553_eval;
case 39555 : return &A4q1_2q239555_eval;
case 39568 : return &A4q1_2q239568_eval;
case 39580 : return &A4q1_2q239580_eval;
case 39588 : return &A4q1_2q239588_eval;
case 39590 : return &A4q1_2q239590_eval;
case 39610 : return &A4q1_2q239610_eval;
case 39622 : return &A4q1_2q239622_eval;
case 39625 : return &A4q1_2q239625_eval;
case 39627 : return &A4q1_2q239627_eval;
case 39640 : return &A4q1_2q239640_eval;
case 39652 : return &A4q1_2q239652_eval;
case 39660 : return &A4q1_2q239660_eval;
case 39662 : return &A4q1_2q239662_eval;
case 39673 : return &A4q1_2q239673_eval;
case 39675 : return &A4q1_2q239675_eval;
case 39678 : return &A4q1_2q239678_eval;
case 39680 : return &A4q1_2q239680_eval;
case 39685 : return &A4q1_2q239685_eval;
case 39687 : return &A4q1_2q239687_eval;
case 39690 : return &A4q1_2q239690_eval;
case 39692 : return &A4q1_2q239692_eval;
case 39751 : return &A4q1_2q239751_eval;
case 39753 : return &A4q1_2q239753_eval;
case 39763 : return &A4q1_2q239763_eval;
case 39765 : return &A4q1_2q239765_eval;
case 39781 : return &A4q1_2q239781_eval;
case 39783 : return &A4q1_2q239783_eval;
case 39786 : return &A4q1_2q239786_eval;
case 39788 : return &A4q1_2q239788_eval;
case 39793 : return &A4q1_2q239793_eval;
case 39795 : return &A4q1_2q239795_eval;
case 39798 : return &A4q1_2q239798_eval;
case 39800 : return &A4q1_2q239800_eval;
case 39823 : return &A4q1_2q239823_eval;
case 39825 : return &A4q1_2q239825_eval;
case 39835 : return &A4q1_2q239835_eval;
case 39837 : return &A4q1_2q239837_eval;
case 39853 : return &A4q1_2q239853_eval;
case 39855 : return &A4q1_2q239855_eval;
case 39858 : return &A4q1_2q239858_eval;
case 39860 : return &A4q1_2q239860_eval;
case 39865 : return &A4q1_2q239865_eval;
case 39867 : return &A4q1_2q239867_eval;
case 39870 : return &A4q1_2q239870_eval;
case 39872 : return &A4q1_2q239872_eval;
case 40186 : return &A4q1_2q240186_eval;
case 40198 : return &A4q1_2q240198_eval;
case 40201 : return &A4q1_2q240201_eval;
case 40203 : return &A4q1_2q240203_eval;
case 40216 : return &A4q1_2q240216_eval;
case 40228 : return &A4q1_2q240228_eval;
case 40236 : return &A4q1_2q240236_eval;
case 40238 : return &A4q1_2q240238_eval;
case 40258 : return &A4q1_2q240258_eval;
case 40270 : return &A4q1_2q240270_eval;
case 40273 : return &A4q1_2q240273_eval;
case 40275 : return &A4q1_2q240275_eval;
case 40288 : return &A4q1_2q240288_eval;
case 40300 : return &A4q1_2q240300_eval;
case 40308 : return &A4q1_2q240308_eval;
case 40310 : return &A4q1_2q240310_eval;
case 40321 : return &A4q1_2q240321_eval;
case 40323 : return &A4q1_2q240323_eval;
case 40326 : return &A4q1_2q240326_eval;
case 40328 : return &A4q1_2q240328_eval;
case 40333 : return &A4q1_2q240333_eval;
case 40335 : return &A4q1_2q240335_eval;
case 40338 : return &A4q1_2q240338_eval;
case 40340 : return &A4q1_2q240340_eval;
case 40396 : return &A4q1_2q240396_eval;
case 40408 : return &A4q1_2q240408_eval;
case 40416 : return &A4q1_2q240416_eval;
case 40418 : return &A4q1_2q240418_eval;
case 40468 : return &A4q1_2q240468_eval;
case 40480 : return &A4q1_2q240480_eval;
case 40488 : return &A4q1_2q240488_eval;
case 40490 : return &A4q1_2q240490_eval;
case 40536 : return &A4q1_2q240536_eval;
case 40538 : return &A4q1_2q240538_eval;
case 40548 : return &A4q1_2q240548_eval;
case 40550 : return &A4q1_2q240550_eval;
case 40618 : return &A4q1_2q240618_eval;
case 40630 : return &A4q1_2q240630_eval;
case 40633 : return &A4q1_2q240633_eval;
case 40635 : return &A4q1_2q240635_eval;
case 40648 : return &A4q1_2q240648_eval;
case 40660 : return &A4q1_2q240660_eval;
case 40668 : return &A4q1_2q240668_eval;
case 40670 : return &A4q1_2q240670_eval;
case 40690 : return &A4q1_2q240690_eval;
case 40702 : return &A4q1_2q240702_eval;
case 40705 : return &A4q1_2q240705_eval;
case 40707 : return &A4q1_2q240707_eval;
case 40720 : return &A4q1_2q240720_eval;
case 40732 : return &A4q1_2q240732_eval;
case 40740 : return &A4q1_2q240740_eval;
case 40742 : return &A4q1_2q240742_eval;
case 40753 : return &A4q1_2q240753_eval;
case 40755 : return &A4q1_2q240755_eval;
case 40758 : return &A4q1_2q240758_eval;
case 40760 : return &A4q1_2q240760_eval;
case 40765 : return &A4q1_2q240765_eval;
case 40767 : return &A4q1_2q240767_eval;
case 40770 : return &A4q1_2q240770_eval;
case 40772 : return &A4q1_2q240772_eval;
case 40828 : return &A4q1_2q240828_eval;
case 40840 : return &A4q1_2q240840_eval;
case 40848 : return &A4q1_2q240848_eval;
case 40850 : return &A4q1_2q240850_eval;
case 40900 : return &A4q1_2q240900_eval;
case 40912 : return &A4q1_2q240912_eval;
case 40920 : return &A4q1_2q240920_eval;
case 40922 : return &A4q1_2q240922_eval;
case 40968 : return &A4q1_2q240968_eval;
case 40970 : return &A4q1_2q240970_eval;
case 40980 : return &A4q1_2q240980_eval;
case 40982 : return &A4q1_2q240982_eval;
case 41041 : return &A4q1_2q241041_eval;
case 41043 : return &A4q1_2q241043_eval;
case 41046 : return &A4q1_2q241046_eval;
case 41048 : return &A4q1_2q241048_eval;
case 41053 : return &A4q1_2q241053_eval;
case 41055 : return &A4q1_2q241055_eval;
case 41058 : return &A4q1_2q241058_eval;
case 41060 : return &A4q1_2q241060_eval;
case 41076 : return &A4q1_2q241076_eval;
case 41078 : return &A4q1_2q241078_eval;
case 41088 : return &A4q1_2q241088_eval;
case 41090 : return &A4q1_2q241090_eval;
case 41113 : return &A4q1_2q241113_eval;
case 41115 : return &A4q1_2q241115_eval;
case 41118 : return &A4q1_2q241118_eval;
case 41120 : return &A4q1_2q241120_eval;
case 41125 : return &A4q1_2q241125_eval;
case 41127 : return &A4q1_2q241127_eval;
case 41130 : return &A4q1_2q241130_eval;
case 41132 : return &A4q1_2q241132_eval;
case 41148 : return &A4q1_2q241148_eval;
case 41150 : return &A4q1_2q241150_eval;
case 41160 : return &A4q1_2q241160_eval;
case 41162 : return &A4q1_2q241162_eval;
case 41518 : return &A4q1_2q241518_eval;
case 41530 : return &A4q1_2q241530_eval;
case 41533 : return &A4q1_2q241533_eval;
case 41535 : return &A4q1_2q241535_eval;
case 41590 : return &A4q1_2q241590_eval;
case 41602 : return &A4q1_2q241602_eval;
case 41605 : return &A4q1_2q241605_eval;
case 41607 : return &A4q1_2q241607_eval;
case 41623 : return &A4q1_2q241623_eval;
case 41625 : return &A4q1_2q241625_eval;
case 41635 : return &A4q1_2q241635_eval;
case 41637 : return &A4q1_2q241637_eval;
case 41698 : return &A4q1_2q241698_eval;
case 41710 : return &A4q1_2q241710_eval;
case 41713 : return &A4q1_2q241713_eval;
case 41715 : return &A4q1_2q241715_eval;
case 41728 : return &A4q1_2q241728_eval;
case 41740 : return &A4q1_2q241740_eval;
case 41748 : return &A4q1_2q241748_eval;
case 41750 : return &A4q1_2q241750_eval;
case 41770 : return &A4q1_2q241770_eval;
case 41782 : return &A4q1_2q241782_eval;
case 41785 : return &A4q1_2q241785_eval;
case 41787 : return &A4q1_2q241787_eval;
case 41800 : return &A4q1_2q241800_eval;
case 41812 : return &A4q1_2q241812_eval;
case 41820 : return &A4q1_2q241820_eval;
case 41822 : return &A4q1_2q241822_eval;
case 41833 : return &A4q1_2q241833_eval;
case 41835 : return &A4q1_2q241835_eval;
case 41838 : return &A4q1_2q241838_eval;
case 41840 : return &A4q1_2q241840_eval;
case 41845 : return &A4q1_2q241845_eval;
case 41847 : return &A4q1_2q241847_eval;
case 41850 : return &A4q1_2q241850_eval;
case 41852 : return &A4q1_2q241852_eval;
case 41950 : return &A4q1_2q241950_eval;
case 41962 : return &A4q1_2q241962_eval;
case 41965 : return &A4q1_2q241965_eval;
case 41967 : return &A4q1_2q241967_eval;
case 42022 : return &A4q1_2q242022_eval;
case 42034 : return &A4q1_2q242034_eval;
case 42037 : return &A4q1_2q242037_eval;
case 42039 : return &A4q1_2q242039_eval;
case 42055 : return &A4q1_2q242055_eval;
case 42057 : return &A4q1_2q242057_eval;
case 42067 : return &A4q1_2q242067_eval;
case 42069 : return &A4q1_2q242069_eval;
case 42130 : return &A4q1_2q242130_eval;
case 42142 : return &A4q1_2q242142_eval;
case 42145 : return &A4q1_2q242145_eval;
case 42147 : return &A4q1_2q242147_eval;
case 42160 : return &A4q1_2q242160_eval;
case 42172 : return &A4q1_2q242172_eval;
case 42180 : return &A4q1_2q242180_eval;
case 42182 : return &A4q1_2q242182_eval;
case 42202 : return &A4q1_2q242202_eval;
case 42214 : return &A4q1_2q242214_eval;
case 42217 : return &A4q1_2q242217_eval;
case 42219 : return &A4q1_2q242219_eval;
case 42232 : return &A4q1_2q242232_eval;
case 42244 : return &A4q1_2q242244_eval;
case 42252 : return &A4q1_2q242252_eval;
case 42254 : return &A4q1_2q242254_eval;
case 42265 : return &A4q1_2q242265_eval;
case 42267 : return &A4q1_2q242267_eval;
case 42270 : return &A4q1_2q242270_eval;
case 42272 : return &A4q1_2q242272_eval;
case 42277 : return &A4q1_2q242277_eval;
case 42279 : return &A4q1_2q242279_eval;
case 42282 : return &A4q1_2q242282_eval;
case 42284 : return &A4q1_2q242284_eval;
case 42343 : return &A4q1_2q242343_eval;
case 42345 : return &A4q1_2q242345_eval;
case 42355 : return &A4q1_2q242355_eval;
case 42357 : return &A4q1_2q242357_eval;
case 42373 : return &A4q1_2q242373_eval;
case 42375 : return &A4q1_2q242375_eval;
case 42378 : return &A4q1_2q242378_eval;
case 42380 : return &A4q1_2q242380_eval;
case 42385 : return &A4q1_2q242385_eval;
case 42387 : return &A4q1_2q242387_eval;
case 42390 : return &A4q1_2q242390_eval;
case 42392 : return &A4q1_2q242392_eval;
case 42415 : return &A4q1_2q242415_eval;
case 42417 : return &A4q1_2q242417_eval;
case 42427 : return &A4q1_2q242427_eval;
case 42429 : return &A4q1_2q242429_eval;
case 42445 : return &A4q1_2q242445_eval;
case 42447 : return &A4q1_2q242447_eval;
case 42450 : return &A4q1_2q242450_eval;
case 42452 : return &A4q1_2q242452_eval;
case 42457 : return &A4q1_2q242457_eval;
case 42459 : return &A4q1_2q242459_eval;
case 42462 : return &A4q1_2q242462_eval;
case 42464 : return &A4q1_2q242464_eval;
case 42778 : return &A4q1_2q242778_eval;
case 42790 : return &A4q1_2q242790_eval;
case 42793 : return &A4q1_2q242793_eval;
case 42795 : return &A4q1_2q242795_eval;
case 42808 : return &A4q1_2q242808_eval;
case 42820 : return &A4q1_2q242820_eval;
case 42828 : return &A4q1_2q242828_eval;
case 42830 : return &A4q1_2q242830_eval;
case 42850 : return &A4q1_2q242850_eval;
case 42862 : return &A4q1_2q242862_eval;
case 42865 : return &A4q1_2q242865_eval;
case 42867 : return &A4q1_2q242867_eval;
case 42880 : return &A4q1_2q242880_eval;
case 42892 : return &A4q1_2q242892_eval;
case 42900 : return &A4q1_2q242900_eval;
case 42902 : return &A4q1_2q242902_eval;
case 42913 : return &A4q1_2q242913_eval;
case 42915 : return &A4q1_2q242915_eval;
case 42918 : return &A4q1_2q242918_eval;
case 42920 : return &A4q1_2q242920_eval;
case 42925 : return &A4q1_2q242925_eval;
case 42927 : return &A4q1_2q242927_eval;
case 42930 : return &A4q1_2q242930_eval;
case 42932 : return &A4q1_2q242932_eval;
case 42988 : return &A4q1_2q242988_eval;
case 43000 : return &A4q1_2q243000_eval;
case 43008 : return &A4q1_2q243008_eval;
case 43010 : return &A4q1_2q243010_eval;
case 43060 : return &A4q1_2q243060_eval;
case 43072 : return &A4q1_2q243072_eval;
case 43080 : return &A4q1_2q243080_eval;
case 43082 : return &A4q1_2q243082_eval;
case 43128 : return &A4q1_2q243128_eval;
case 43130 : return &A4q1_2q243130_eval;
case 43140 : return &A4q1_2q243140_eval;
case 43142 : return &A4q1_2q243142_eval;
case 43210 : return &A4q1_2q243210_eval;
case 43222 : return &A4q1_2q243222_eval;
case 43225 : return &A4q1_2q243225_eval;
case 43227 : return &A4q1_2q243227_eval;
case 43240 : return &A4q1_2q243240_eval;
case 43252 : return &A4q1_2q243252_eval;
case 43260 : return &A4q1_2q243260_eval;
case 43262 : return &A4q1_2q243262_eval;
case 43282 : return &A4q1_2q243282_eval;
case 43294 : return &A4q1_2q243294_eval;
case 43297 : return &A4q1_2q243297_eval;
case 43299 : return &A4q1_2q243299_eval;
case 43312 : return &A4q1_2q243312_eval;
case 43324 : return &A4q1_2q243324_eval;
case 43332 : return &A4q1_2q243332_eval;
case 43334 : return &A4q1_2q243334_eval;
case 43345 : return &A4q1_2q243345_eval;
case 43347 : return &A4q1_2q243347_eval;
case 43350 : return &A4q1_2q243350_eval;
case 43352 : return &A4q1_2q243352_eval;
case 43357 : return &A4q1_2q243357_eval;
case 43359 : return &A4q1_2q243359_eval;
case 43362 : return &A4q1_2q243362_eval;
case 43364 : return &A4q1_2q243364_eval;
case 43420 : return &A4q1_2q243420_eval;
case 43432 : return &A4q1_2q243432_eval;
case 43440 : return &A4q1_2q243440_eval;
case 43442 : return &A4q1_2q243442_eval;
case 43492 : return &A4q1_2q243492_eval;
case 43504 : return &A4q1_2q243504_eval;
case 43512 : return &A4q1_2q243512_eval;
case 43514 : return &A4q1_2q243514_eval;
case 43560 : return &A4q1_2q243560_eval;
case 43562 : return &A4q1_2q243562_eval;
case 43572 : return &A4q1_2q243572_eval;
case 43574 : return &A4q1_2q243574_eval;
case 43633 : return &A4q1_2q243633_eval;
case 43635 : return &A4q1_2q243635_eval;
case 43638 : return &A4q1_2q243638_eval;
case 43640 : return &A4q1_2q243640_eval;
case 43645 : return &A4q1_2q243645_eval;
case 43647 : return &A4q1_2q243647_eval;
case 43650 : return &A4q1_2q243650_eval;
case 43652 : return &A4q1_2q243652_eval;
case 43668 : return &A4q1_2q243668_eval;
case 43670 : return &A4q1_2q243670_eval;
case 43680 : return &A4q1_2q243680_eval;
case 43682 : return &A4q1_2q243682_eval;
case 43705 : return &A4q1_2q243705_eval;
case 43707 : return &A4q1_2q243707_eval;
case 43710 : return &A4q1_2q243710_eval;
case 43712 : return &A4q1_2q243712_eval;
case 43717 : return &A4q1_2q243717_eval;
case 43719 : return &A4q1_2q243719_eval;
case 43722 : return &A4q1_2q243722_eval;
case 43724 : return &A4q1_2q243724_eval;
case 43740 : return &A4q1_2q243740_eval;
case 43742 : return &A4q1_2q243742_eval;
case 43752 : return &A4q1_2q243752_eval;
case 43754 : return &A4q1_2q243754_eval;
case 44071 : return &A4q1_2q244071_eval;
case 44073 : return &A4q1_2q244073_eval;
case 44083 : return &A4q1_2q244083_eval;
case 44085 : return &A4q1_2q244085_eval;
case 44101 : return &A4q1_2q244101_eval;
case 44103 : return &A4q1_2q244103_eval;
case 44106 : return &A4q1_2q244106_eval;
case 44108 : return &A4q1_2q244108_eval;
case 44113 : return &A4q1_2q244113_eval;
case 44115 : return &A4q1_2q244115_eval;
case 44118 : return &A4q1_2q244118_eval;
case 44120 : return &A4q1_2q244120_eval;
case 44143 : return &A4q1_2q244143_eval;
case 44145 : return &A4q1_2q244145_eval;
case 44155 : return &A4q1_2q244155_eval;
case 44157 : return &A4q1_2q244157_eval;
case 44173 : return &A4q1_2q244173_eval;
case 44175 : return &A4q1_2q244175_eval;
case 44178 : return &A4q1_2q244178_eval;
case 44180 : return &A4q1_2q244180_eval;
case 44185 : return &A4q1_2q244185_eval;
case 44187 : return &A4q1_2q244187_eval;
case 44190 : return &A4q1_2q244190_eval;
case 44192 : return &A4q1_2q244192_eval;
case 44281 : return &A4q1_2q244281_eval;
case 44283 : return &A4q1_2q244283_eval;
case 44286 : return &A4q1_2q244286_eval;
case 44288 : return &A4q1_2q244288_eval;
case 44293 : return &A4q1_2q244293_eval;
case 44295 : return &A4q1_2q244295_eval;
case 44298 : return &A4q1_2q244298_eval;
case 44300 : return &A4q1_2q244300_eval;
case 44316 : return &A4q1_2q244316_eval;
case 44318 : return &A4q1_2q244318_eval;
case 44328 : return &A4q1_2q244328_eval;
case 44330 : return &A4q1_2q244330_eval;
case 44353 : return &A4q1_2q244353_eval;
case 44355 : return &A4q1_2q244355_eval;
case 44358 : return &A4q1_2q244358_eval;
case 44360 : return &A4q1_2q244360_eval;
case 44365 : return &A4q1_2q244365_eval;
case 44367 : return &A4q1_2q244367_eval;
case 44370 : return &A4q1_2q244370_eval;
case 44372 : return &A4q1_2q244372_eval;
case 44388 : return &A4q1_2q244388_eval;
case 44390 : return &A4q1_2q244390_eval;
case 44400 : return &A4q1_2q244400_eval;
case 44402 : return &A4q1_2q244402_eval;
case 44503 : return &A4q1_2q244503_eval;
case 44505 : return &A4q1_2q244505_eval;
case 44515 : return &A4q1_2q244515_eval;
case 44517 : return &A4q1_2q244517_eval;
case 44533 : return &A4q1_2q244533_eval;
case 44535 : return &A4q1_2q244535_eval;
case 44538 : return &A4q1_2q244538_eval;
case 44540 : return &A4q1_2q244540_eval;
case 44545 : return &A4q1_2q244545_eval;
case 44547 : return &A4q1_2q244547_eval;
case 44550 : return &A4q1_2q244550_eval;
case 44552 : return &A4q1_2q244552_eval;
case 44575 : return &A4q1_2q244575_eval;
case 44577 : return &A4q1_2q244577_eval;
case 44587 : return &A4q1_2q244587_eval;
case 44589 : return &A4q1_2q244589_eval;
case 44605 : return &A4q1_2q244605_eval;
case 44607 : return &A4q1_2q244607_eval;
case 44610 : return &A4q1_2q244610_eval;
case 44612 : return &A4q1_2q244612_eval;
case 44617 : return &A4q1_2q244617_eval;
case 44619 : return &A4q1_2q244619_eval;
case 44622 : return &A4q1_2q244622_eval;
case 44624 : return &A4q1_2q244624_eval;
case 44713 : return &A4q1_2q244713_eval;
case 44715 : return &A4q1_2q244715_eval;
case 44718 : return &A4q1_2q244718_eval;
case 44720 : return &A4q1_2q244720_eval;
case 44725 : return &A4q1_2q244725_eval;
case 44727 : return &A4q1_2q244727_eval;
case 44730 : return &A4q1_2q244730_eval;
case 44732 : return &A4q1_2q244732_eval;
case 44748 : return &A4q1_2q244748_eval;
case 44750 : return &A4q1_2q244750_eval;
case 44760 : return &A4q1_2q244760_eval;
case 44762 : return &A4q1_2q244762_eval;
case 44785 : return &A4q1_2q244785_eval;
case 44787 : return &A4q1_2q244787_eval;
case 44790 : return &A4q1_2q244790_eval;
case 44792 : return &A4q1_2q244792_eval;
case 44797 : return &A4q1_2q244797_eval;
case 44799 : return &A4q1_2q244799_eval;
case 44802 : return &A4q1_2q244802_eval;
case 44804 : return &A4q1_2q244804_eval;
case 44820 : return &A4q1_2q244820_eval;
case 44822 : return &A4q1_2q244822_eval;
case 44832 : return &A4q1_2q244832_eval;
case 44834 : return &A4q1_2q244834_eval;
#endif
default: _WARNING3("Unknown pointer amplitude (*A4q1_2q2_Tree_Ptr(int hc)) - case:",hc , " - throw BH error.");
		throw BHerror("case missing for tree amplitude!");

}
}



template complex<R>  (*A4q1_2q2_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A4q1_2q2_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A4q1_2q2_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A4q1_2q2_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
44316: Q[m, 1] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37836: Q[m, 1] Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44748: Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38268: Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41076: Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43668: Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15156: Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30708: Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33516: Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36108: Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14076: Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29628: Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44388: Q[m, 1] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37908: Q[m, 1] Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44820: Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38340: Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41148: Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43740: Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15228: Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30780: Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33588: Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36180: Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14148: Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29700: Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
40536: Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43128: Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14616: Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30168: Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40968: Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43560: Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2]
15048: Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30600: Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10296: Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25848: Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12888: Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28440: Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1]
32796: Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2]
35388: Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2]
13356: Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1]
28908: Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1]
33228: Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2]
35820: Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2]
13788: Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1]
29340: Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1]
10116: Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1]
25668: Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1]
12708: Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1]
28260: Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1]
44286: Q[m, 1] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37806: Q[m, 1] Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44718: Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38238: Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41046: Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43638: Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15126: Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30678: Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33486: Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36078: Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14046: Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29598: Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44106: Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37626: Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44538: Q[m, 1] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38058: Q[m, 1] Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39786: Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42378: Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7386: Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22938: Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32226: Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34818: Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6306: Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21858: Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44358: Q[m, 1] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37878: Q[m, 1] Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44790: Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38310: Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41118: Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43710: Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15198: Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30750: Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33558: Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36150: Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14118: Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29670: Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44178: Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37698: Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44610: Q[m, 1] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38130: Q[m, 1] Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39858: Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42450: Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7458: Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
23010: Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32298: Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34890: Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6378: Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21930: Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
40326: Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2]
42918: Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2]
14406: Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1]
29958: Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1]
39246: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2]
41838: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2]
6846: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1]
22398: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1]
40758: Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2]
43350: Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2]
14838: Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1]
30390: Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1]
39678: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2]
42270: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2]
7278: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1]
22830: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1]
9006: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1]
24558: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1]
2526: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1]
18078: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1]
11598: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1]
27150: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1]
5118: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1]
20670: Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1]
32586: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2]
35178: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2]
13146: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1]
28698: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1]
31506: Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2]
34098: Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2]
5586: Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1]
21138: Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1]
33018: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2]
35610: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2]
13578: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1]
29130: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1]
31938: Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2]
34530: Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2]
6018: Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1]
21570: Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1]
8826: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1]
24378: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1]
2346: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1]
17898: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1]
11418: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1]
26970: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1]
4938: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1]
20490: Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1]
44328: Q[m, 1] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37848: Q[m, 1] Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44760: Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38280: Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41088: Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43680: Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15168: Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30720: Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33528: Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36120: Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14088: Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29640: Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44400: Q[m, 1] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37920: Q[m, 1] Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44832: Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38352: Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41160: Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43752: Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15240: Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30792: Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33600: Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36192: Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14160: Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29712: Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
40548: Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43140: Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14628: Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30180: Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40980: Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43572: Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2]
15060: Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30612: Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10308: Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25860: Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12900: Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28452: Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1]
32808: Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2]
35400: Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2]
13368: Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1]
28920: Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1]
33240: Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2]
35832: Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2]
13800: Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1]
29352: Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1]
10128: Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1]
25680: Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1]
12720: Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1]
28272: Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1]
44298: Q[m, 1] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37818: Q[m, 1] Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44730: Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38250: Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41058: Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43650: Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15138: Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30690: Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33498: Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36090: Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14058: Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29610: Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44118: Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37638: Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44550: Q[m, 1] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38070: Q[m, 1] Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39798: Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42390: Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7398: Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22950: Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32238: Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34830: Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6318: Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21870: Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44370: Q[m, 1] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37890: Q[m, 1] Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44802: Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38322: Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41130: Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43722: Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15210: Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30762: Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33570: Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36162: Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14130: Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29682: Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44190: Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37710: Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44622: Q[m, 1] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38142: Q[m, 1] Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39870: Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42462: Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7470: Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
23022: Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32310: Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34902: Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6390: Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21942: Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
40338: Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2]
42930: Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2]
14418: Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1]
29970: Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1]
39258: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2]
41850: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2]
6858: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1]
22410: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1]
40770: Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2]
43362: Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2]
14850: Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1]
30402: Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1]
39690: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2]
42282: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2]
7290: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1]
22842: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1]
9018: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1]
24570: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1]
2538: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1]
18090: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1]
11610: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1]
27162: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1]
5130: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1]
20682: Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1]
32598: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2]
35190: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2]
13158: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1]
28710: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1]
31518: Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2]
34110: Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2]
5598: Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1]
21150: Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1]
33030: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2]
35622: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2]
13590: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1]
29142: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1]
31950: Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2]
34542: Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2]
6030: Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1]
21582: Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1]
8838: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1]
24390: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1]
2358: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1]
17910: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1]
11430: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1]
26982: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1]
4950: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1]
20502: Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1]
40416: Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43008: Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14496: Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30048: Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40848: Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43440: Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2]
14928: Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30480: Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10176: Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25728: Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12768: Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28320: Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1]
40236: Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42828: Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14316: Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29868: Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39156: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41748: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6756: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22308: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40668: Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43260: Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14748: Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30300: Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39588: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42180: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7188: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22740: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8916: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24468: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2436: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
17988: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11508: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27060: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5028: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20580: Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
40488: Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43080: Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14568: Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30120: Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40920: Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43512: Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2]
15000: Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30552: Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10248: Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25800: Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12840: Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28392: Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1]
40308: Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42900: Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14388: Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29940: Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39228: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41820: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6828: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22380: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40740: Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43332: Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14820: Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30372: Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39660: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42252: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7260: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22812: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8988: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24540: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2508: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
18060: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11580: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27132: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5100: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20652: Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
9276: Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1]
24828: Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1]
11868: Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1]
27420: Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1]
8196: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1]
23748: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1]
1716: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1]
17268: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1]
10788: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1]
26340: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1]
4308: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1]
19860: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1]
9708: Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1]
25260: Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1]
12300: Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1]
27852: Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1]
8628: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1]
24180: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1]
2148: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1]
17700: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1]
11220: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1]
26772: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1]
4740: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1]
20292: Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1]
32646: Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2]
35238: Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2]
13206: Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1]
28758: Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1]
33078: Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2]
35670: Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2]
13638: Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1]
29190: Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1]
9966: Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1]
25518: Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1]
12558: Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1]
28110: Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1]
32466: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35058: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13026: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28578: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31386: Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33978: Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5466: Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21018: Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32898: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35490: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13458: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29010: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31818: Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34410: Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5898: Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21450: Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8706: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24258: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2226: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17778: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11298: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26850: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4818: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20370: Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
32718: Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2]
35310: Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2]
13278: Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1]
28830: Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1]
33150: Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2]
35742: Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2]
13710: Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1]
29262: Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1]
10038: Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1]
25590: Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1]
12630: Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1]
28182: Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1]
32538: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35130: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13098: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28650: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31458: Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
34050: Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5538: Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21090: Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32970: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35562: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13530: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29082: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31890: Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34482: Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5970: Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21522: Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8778: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24330: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2298: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17850: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11370: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26922: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4890: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20442: Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
9246: Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1]
24798: Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1]
11838: Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1]
27390: Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1]
8166: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1]
23718: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1]
1686: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1]
17238: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1]
10758: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1]
26310: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1]
4278: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1]
19830: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1]
9678: Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1]
25230: Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1]
12270: Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1]
27822: Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1]
8598: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1]
24150: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1]
2118: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1]
17670: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1]
11190: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1]
26742: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1]
4710: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1]
20262: Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1]
44281: Qb[m, 1] Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37801: Qb[m, 1] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44713: Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38233: Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41041: Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43633: Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15121: Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30673: Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33481: Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36073: Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14041: Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29593: Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44101: Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37621: Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44533: Qb[m, 1] Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38053: Qb[m, 1] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39781: Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42373: Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7381: Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22933: Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32221: Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34813: Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6301: Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21853: Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44353: Qb[m, 1] Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37873: Qb[m, 1] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44785: Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38305: Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41113: Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43705: Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15193: Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30745: Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33553: Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36145: Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14113: Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29665: Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44173: Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37693: Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44605: Qb[m, 1] Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38125: Qb[m, 1] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39853: Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42445: Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7453: Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
23005: Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32293: Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34885: Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6373: Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21925: Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
40321: Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2]
42913: Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2]
14401: Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1]
29953: Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1]
39241: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2]
41833: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2]
6841: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1]
22393: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1]
40753: Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2]
43345: Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2]
14833: Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1]
30385: Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1]
39673: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2]
42265: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2]
7273: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1]
22825: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1]
9001: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1]
24553: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1]
2521: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1]
18073: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1]
11593: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1]
27145: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1]
5113: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1]
20665: Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1]
32581: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2]
35173: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2]
13141: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1]
28693: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1]
31501: Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2]
34093: Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2]
5581: Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1]
21133: Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1]
33013: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2]
35605: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2]
13573: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1]
29125: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1]
31933: Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2]
34525: Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2]
6013: Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1]
21565: Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1]
8821: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1]
24373: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1]
2341: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1]
17893: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1]
11413: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1]
26965: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1]
4933: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1]
20485: Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1]
44071: Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37591: Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44503: Qb[m, 1] Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38023: Qb[m, 1] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39751: Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42343: Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7351: Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22903: Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32191: Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34783: Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6271: Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21823: Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44143: Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37663: Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44575: Qb[m, 1] Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38095: Qb[m, 1] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39823: Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42415: Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7423: Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22975: Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32263: Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34855: Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6343: Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21895: Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
39031: Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2]
41623: Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2]
6631: Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1]
22183: Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1]
39463: Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2]
42055: Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2]
7063: Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1]
22615: Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1]
1231: Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1]
16783: Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1]
3823: Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1]
19375: Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1]
31291: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2]
33883: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2]
5371: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1]
20923: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1]
31723: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2]
34315: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2]
5803: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1]
21355: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1]
1051: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1]
16603: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1]
3643: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1]
19195: Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1]
44293: Qb[m, 1] Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37813: Qb[m, 1] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44725: Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38245: Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41053: Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43645: Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15133: Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30685: Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33493: Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36085: Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14053: Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29605: Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44113: Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37633: Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44545: Qb[m, 1] Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38065: Qb[m, 1] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39793: Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42385: Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7393: Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22945: Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32233: Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34825: Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6313: Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21865: Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44365: Qb[m, 1] Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37885: Qb[m, 1] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44797: Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38317: Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41125: Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43717: Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15205: Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30757: Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33565: Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36157: Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14125: Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29677: Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44185: Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37705: Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44617: Qb[m, 1] Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38137: Qb[m, 1] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39865: Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42457: Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7465: Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
23017: Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32305: Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34897: Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6385: Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21937: Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
40333: Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2]
42925: Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2]
14413: Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1]
29965: Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1]
39253: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2]
41845: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2]
6853: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1]
22405: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1]
40765: Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2]
43357: Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2]
14845: Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1]
30397: Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1]
39685: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2]
42277: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2]
7285: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1]
22837: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1]
9013: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1]
24565: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1]
2533: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1]
18085: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1]
11605: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1]
27157: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1]
5125: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1]
20677: Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1]
32593: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2]
35185: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2]
13153: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1]
28705: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1]
31513: Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2]
34105: Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2]
5593: Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1]
21145: Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1]
33025: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2]
35617: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2]
13585: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1]
29137: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1]
31945: Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2]
34537: Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2]
6025: Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1]
21577: Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1]
8833: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1]
24385: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1]
2353: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1]
17905: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1]
11425: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1]
26977: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1]
4945: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1]
20497: Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1]
44083: Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37603: Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44515: Qb[m, 1] Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38035: Qb[m, 1] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39763: Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42355: Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7363: Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22915: Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32203: Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34795: Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6283: Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21835: Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44155: Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37675: Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44587: Qb[m, 1] Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38107: Qb[m, 1] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39835: Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42427: Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7435: Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22987: Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32275: Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34867: Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6355: Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21907: Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
39043: Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2]
41635: Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2]
6643: Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1]
22195: Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1]
39475: Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2]
42067: Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2]
7075: Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1]
22627: Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1]
1243: Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1]
16795: Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1]
3835: Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1]
19387: Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1]
31303: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2]
33895: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2]
5383: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1]
20935: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1]
31735: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2]
34327: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2]
5815: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1]
21367: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1]
1063: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1]
16615: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1]
3655: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1]
19207: Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1]
40201: Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42793: Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14281: Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29833: Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39121: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41713: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6721: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22273: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40633: Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43225: Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14713: Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30265: Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39553: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42145: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7153: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22705: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8881: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24433: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2401: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
17953: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11473: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27025: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
4993: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20545: Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
38941: Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2]
41533: Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2]
6541: Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1]
22093: Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1]
39373: Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2]
41965: Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2]
6973: Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1]
22525: Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1]
1141: Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1]
16693: Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1]
3733: Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1]
19285: Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1]
40273: Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42865: Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14353: Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29905: Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39193: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41785: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6793: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22345: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40705: Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43297: Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14785: Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30337: Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39625: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42217: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7225: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22777: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8953: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24505: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2473: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
18025: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11545: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27097: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5065: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20617: Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
39013: Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2]
41605: Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2]
6613: Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1]
22165: Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1]
39445: Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2]
42037: Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2]
7045: Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1]
22597: Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1]
1213: Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1]
16765: Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1]
3805: Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1]
19357: Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1]
7981: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1]
23533: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1]
1501: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1]
17053: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1]
10573: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1]
26125: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1]
4093: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1]
19645: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1]
421: Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1]
15973: Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1]
3013: Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1]
18565: Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1]
8413: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1]
23965: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1]
1933: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1]
17485: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1]
11005: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1]
26557: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1]
4525: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1]
20077: Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1]
853: Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1]
16405: Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1]
3445: Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1]
18997: Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1]
32431: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35023: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
12991: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28543: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31351: Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33943: Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5431: Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
20983: Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32863: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35455: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13423: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
28975: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31783: Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34375: Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5863: Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21415: Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8671: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24223: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2191: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17743: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11263: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26815: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4783: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20335: Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
31171: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2]
33763: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2]
5251: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1]
20803: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1]
31603: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2]
34195: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2]
5683: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1]
21235: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1]
931: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1]
16483: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1]
3523: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1]
19075: Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1]
32503: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35095: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13063: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28615: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31423: Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
34015: Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5503: Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21055: Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32935: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35527: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13495: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29047: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31855: Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34447: Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5935: Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21487: Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8743: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24295: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2263: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17815: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11335: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26887: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4855: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20407: Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
31243: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2]
33835: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2]
5323: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1]
20875: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1]
31675: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2]
34267: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2]
5755: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1]
21307: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1]
1003: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1]
16555: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1]
3595: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1]
19147: Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1]
7951: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1]
23503: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1]
1471: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1]
17023: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1]
10543: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1]
26095: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1]
4063: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1]
19615: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1]
391: Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1]
15943: Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1]
2983: Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1]
18535: Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1]
8383: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1]
23935: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1]
1903: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1]
17455: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1]
10975: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1]
26527: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1]
4495: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1]
20047: Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1]
823: Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1]
16375: Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1]
3415: Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1]
18967: Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1]
44318: Q[p, 1] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37838: Q[p, 1] Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44750: Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38270: Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41078: Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43670: Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15158: Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30710: Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33518: Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36110: Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14078: Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29630: Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44390: Q[p, 1] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37910: Q[p, 1] Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44822: Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38342: Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41150: Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43742: Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15230: Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30782: Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33590: Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36182: Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14150: Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29702: Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
40538: Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43130: Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14618: Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30170: Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40970: Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43562: Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2]
15050: Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30602: Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10298: Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25850: Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12890: Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28442: Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1]
32798: Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2]
35390: Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2]
13358: Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1]
28910: Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1]
33230: Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2]
35822: Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2]
13790: Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1]
29342: Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1]
10118: Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1]
25670: Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1]
12710: Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1]
28262: Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1]
44288: Q[p, 1] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37808: Q[p, 1] Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44720: Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38240: Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41048: Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43640: Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15128: Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30680: Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33488: Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36080: Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14048: Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29600: Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44108: Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37628: Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44540: Q[p, 1] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38060: Q[p, 1] Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39788: Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42380: Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7388: Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22940: Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32228: Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34820: Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6308: Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21860: Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44360: Q[p, 1] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37880: Q[p, 1] Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44792: Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38312: Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41120: Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43712: Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15200: Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30752: Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33560: Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36152: Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14120: Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29672: Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44180: Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37700: Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44612: Q[p, 1] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38132: Q[p, 1] Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39860: Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42452: Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7460: Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
23012: Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32300: Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34892: Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6380: Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21932: Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
40328: Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2]
42920: Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2]
14408: Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1]
29960: Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1]
39248: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2]
41840: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2]
6848: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1]
22400: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1]
40760: Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2]
43352: Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2]
14840: Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1]
30392: Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1]
39680: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2]
42272: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2]
7280: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1]
22832: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1]
9008: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1]
24560: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1]
2528: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1]
18080: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1]
11600: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1]
27152: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1]
5120: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1]
20672: Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1]
32588: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2]
35180: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2]
13148: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1]
28700: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1]
31508: Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2]
34100: Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2]
5588: Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1]
21140: Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1]
33020: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2]
35612: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2]
13580: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1]
29132: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1]
31940: Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2]
34532: Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2]
6020: Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1]
21572: Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1]
8828: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1]
24380: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1]
2348: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1]
17900: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1]
11420: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1]
26972: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1]
4940: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1]
20492: Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1]
44330: Q[p, 1] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37850: Q[p, 1] Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44762: Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38282: Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41090: Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43682: Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15170: Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30722: Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33530: Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36122: Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14090: Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29642: Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44402: Q[p, 1] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37922: Q[p, 1] Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44834: Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38354: Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41162: Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43754: Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15242: Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30794: Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33602: Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36194: Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14162: Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29714: Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
40550: Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43142: Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14630: Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30182: Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40982: Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43574: Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2]
15062: Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30614: Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10310: Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25862: Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12902: Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28454: Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1]
32810: Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2]
35402: Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2]
13370: Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1]
28922: Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1]
33242: Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2]
35834: Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2]
13802: Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1]
29354: Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1]
10130: Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1]
25682: Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1]
12722: Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1]
28274: Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1]
44300: Q[p, 1] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37820: Q[p, 1] Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44732: Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38252: Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41060: Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43652: Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15140: Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30692: Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33500: Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36092: Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14060: Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29612: Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44120: Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37640: Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44552: Q[p, 1] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38072: Q[p, 1] Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39800: Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42392: Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7400: Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22952: Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32240: Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34832: Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6320: Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21872: Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44372: Q[p, 1] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37892: Q[p, 1] Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44804: Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38324: Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41132: Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43724: Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15212: Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30764: Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33572: Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36164: Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14132: Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29684: Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44192: Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37712: Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44624: Q[p, 1] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38144: Q[p, 1] Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39872: Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42464: Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7472: Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
23024: Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32312: Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34904: Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6392: Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21944: Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
40340: Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2]
42932: Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2]
14420: Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1]
29972: Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1]
39260: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2]
41852: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2]
6860: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1]
22412: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1]
40772: Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2]
43364: Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2]
14852: Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1]
30404: Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1]
39692: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2]
42284: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2]
7292: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1]
22844: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1]
9020: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1]
24572: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1]
2540: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1]
18092: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1]
11612: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1]
27164: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1]
5132: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1]
20684: Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1]
32600: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2]
35192: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2]
13160: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1]
28712: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1]
31520: Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2]
34112: Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2]
5600: Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1]
21152: Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1]
33032: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2]
35624: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2]
13592: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1]
29144: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1]
31952: Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2]
34544: Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2]
6032: Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1]
21584: Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1]
8840: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1]
24392: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1]
2360: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1]
17912: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1]
11432: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1]
26984: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1]
4952: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1]
20504: Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1]
40418: Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43010: Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14498: Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30050: Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40850: Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43442: Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2]
14930: Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30482: Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10178: Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25730: Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12770: Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28322: Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1]
40238: Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42830: Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14318: Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29870: Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39158: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41750: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6758: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22310: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40670: Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43262: Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14750: Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30302: Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39590: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42182: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7190: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22742: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8918: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24470: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2438: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
17990: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11510: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27062: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5030: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20582: Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
40490: Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43082: Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14570: Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30122: Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40922: Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43514: Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2]
15002: Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30554: Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10250: Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25802: Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12842: Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28394: Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1]
40310: Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42902: Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14390: Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29942: Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39230: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41822: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6830: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22382: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40742: Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43334: Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14822: Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30374: Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39662: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42254: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7262: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22814: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8990: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24542: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2510: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
18062: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11582: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27134: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5102: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20654: Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
9278: Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1]
24830: Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1]
11870: Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1]
27422: Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1]
8198: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1]
23750: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1]
1718: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1]
17270: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1]
10790: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1]
26342: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1]
4310: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1]
19862: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1]
9710: Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1]
25262: Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1]
12302: Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1]
27854: Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1]
8630: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1]
24182: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1]
2150: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1]
17702: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1]
11222: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1]
26774: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1]
4742: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1]
20294: Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1]
32648: Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2]
35240: Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2]
13208: Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1]
28760: Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1]
33080: Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2]
35672: Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2]
13640: Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1]
29192: Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1]
9968: Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1]
25520: Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1]
12560: Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1]
28112: Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1]
32468: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35060: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13028: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28580: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31388: Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33980: Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5468: Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21020: Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32900: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35492: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13460: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29012: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31820: Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34412: Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5900: Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21452: Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8708: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24260: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2228: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17780: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11300: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26852: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4820: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20372: Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
32720: Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2]
35312: Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2]
13280: Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1]
28832: Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1]
33152: Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2]
35744: Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2]
13712: Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1]
29264: Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1]
10040: Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1]
25592: Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1]
12632: Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1]
28184: Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1]
32540: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35132: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13100: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28652: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31460: Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
34052: Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5540: Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21092: Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32972: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35564: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13532: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29084: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31892: Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34484: Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5972: Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21524: Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8780: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24332: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2300: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17852: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11372: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26924: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4892: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20444: Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
9248: Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1]
24800: Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1]
11840: Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1]
27392: Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1]
8168: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1]
23720: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1]
1688: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1]
17240: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1]
10760: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1]
26312: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1]
4280: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1]
19832: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1]
9680: Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1]
25232: Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1]
12272: Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1]
27824: Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1]
8600: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1]
24152: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1]
2120: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1]
17672: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1]
11192: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1]
26744: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1]
4712: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1]
20264: Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1]
44283: Qb[p, 1] Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37803: Qb[p, 1] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44715: Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38235: Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41043: Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43635: Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15123: Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30675: Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33483: Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36075: Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14043: Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29595: Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44103: Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37623: Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44535: Qb[p, 1] Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38055: Qb[p, 1] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39783: Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42375: Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7383: Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22935: Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32223: Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34815: Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6303: Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21855: Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44355: Qb[p, 1] Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37875: Qb[p, 1] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44787: Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38307: Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41115: Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43707: Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15195: Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30747: Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33555: Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36147: Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14115: Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29667: Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44175: Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37695: Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44607: Qb[p, 1] Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38127: Qb[p, 1] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39855: Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42447: Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7455: Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
23007: Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32295: Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34887: Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6375: Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21927: Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
40323: Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2]
42915: Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2]
14403: Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1]
29955: Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1]
39243: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2]
41835: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2]
6843: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1]
22395: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1]
40755: Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2]
43347: Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2]
14835: Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1]
30387: Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1]
39675: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2]
42267: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2]
7275: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1]
22827: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1]
9003: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1]
24555: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1]
2523: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1]
18075: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1]
11595: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1]
27147: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1]
5115: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1]
20667: Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1]
32583: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2]
35175: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2]
13143: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1]
28695: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1]
31503: Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2]
34095: Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2]
5583: Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1]
21135: Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1]
33015: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2]
35607: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2]
13575: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1]
29127: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1]
31935: Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2]
34527: Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2]
6015: Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1]
21567: Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1]
8823: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1]
24375: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1]
2343: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1]
17895: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1]
11415: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1]
26967: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1]
4935: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1]
20487: Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1]
44073: Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37593: Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44505: Qb[p, 1] Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38025: Qb[p, 1] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39753: Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42345: Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7353: Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22905: Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32193: Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34785: Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6273: Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21825: Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44145: Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37665: Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44577: Qb[p, 1] Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38097: Qb[p, 1] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39825: Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42417: Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7425: Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22977: Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32265: Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34857: Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6345: Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21897: Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
39033: Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2]
41625: Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2]
6633: Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1]
22185: Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1]
39465: Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2]
42057: Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2]
7065: Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1]
22617: Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1]
1233: Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1]
16785: Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1]
3825: Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1]
19377: Qb[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1]
31293: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2]
33885: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2]
5373: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1]
20925: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1]
31725: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2]
34317: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2]
5805: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1]
21357: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1]
1053: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1]
16605: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1]
3645: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1]
19197: Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1]
44295: Qb[p, 1] Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37815: Qb[p, 1] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44727: Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38247: Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41055: Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43647: Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15135: Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30687: Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33495: Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36087: Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14055: Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29607: Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44115: Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37635: Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44547: Qb[p, 1] Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38067: Qb[p, 1] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39795: Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42387: Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7395: Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22947: Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32235: Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34827: Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6315: Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21867: Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44367: Qb[p, 1] Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 2]
37887: Qb[p, 1] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 2]
44799: Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2]
38319: Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2]
41127: Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2]
43719: Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2]
15207: Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1]
30759: Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1]
33567: Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2]
36159: Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2]
14127: Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1]
29679: Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1]
44187: Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37707: Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44619: Qb[p, 1] Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38139: Qb[p, 1] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39867: Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42459: Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7467: Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
23019: Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32307: Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34899: Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6387: Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21939: Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
40335: Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2]
42927: Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2]
14415: Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1]
29967: Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1]
39255: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2]
41847: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2]
6855: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1]
22407: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1]
40767: Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2]
43359: Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2]
14847: Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1]
30399: Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1]
39687: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2]
42279: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2]
7287: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1]
22839: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1]
9015: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1]
24567: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1]
2535: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1]
18087: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1]
11607: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1]
27159: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1]
5127: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1]
20679: Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1]
32595: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2]
35187: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2]
13155: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1]
28707: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1]
31515: Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2]
34107: Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2]
5595: Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1]
21147: Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1]
33027: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2]
35619: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2]
13587: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1]
29139: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1]
31947: Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2]
34539: Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2]
6027: Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1]
21579: Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1]
8835: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1]
24387: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1]
2355: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1]
17907: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1]
11427: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1]
26979: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1]
4947: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1]
20499: Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1]
44085: Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37605: Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44517: Qb[p, 1] Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38037: Qb[p, 1] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39765: Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42357: Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7365: Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22917: Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32205: Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34797: Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6285: Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21837: Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
44157: Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 2]
37677: Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Q[m, 2]
44589: Qb[p, 1] Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2]
38109: Qb[p, 1] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2]
39837: Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2]
42429: Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2]
7437: Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1]
22989: Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1]
32277: Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2]
34869: Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2]
6357: Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1]
21909: Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1]
39045: Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2]
41637: Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2]
6645: Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1]
22197: Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1]
39477: Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2]
42069: Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2]
7077: Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1]
22629: Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1]
1245: Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1]
16797: Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1]
3837: Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1]
19389: Qb[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1]
31305: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2]
33897: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2]
5385: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1]
20937: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1]
31737: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2]
34329: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2]
5817: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1]
21369: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1]
1065: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1]
16617: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1]
3657: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1]
19209: Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1]
40203: Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42795: Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14283: Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29835: Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39123: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41715: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6723: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22275: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40635: Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43227: Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14715: Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30267: Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39555: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42147: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7155: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22707: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8883: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24435: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2403: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
17955: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11475: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27027: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
4995: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20547: Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
38943: Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2]
41535: Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2]
6543: Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1]
22095: Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1]
39375: Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2]
41967: Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2]
6975: Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1]
22527: Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1]
1143: Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1]
16695: Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1]
3735: Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1]
19287: Qb[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1]
40275: Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42867: Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14355: Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29907: Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39195: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41787: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6795: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22347: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40707: Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43299: Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14787: Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30339: Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39627: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42219: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7227: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22779: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8955: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24507: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2475: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
18027: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11547: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27099: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5067: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20619: Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
39015: Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2]
41607: Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2]
6615: Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1]
22167: Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1]
39447: Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2]
42039: Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2]
7047: Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1]
22599: Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1]
1215: Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1]
16767: Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1]
3807: Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1]
19359: Qb[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1]
7983: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1]
23535: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1]
1503: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1]
17055: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1]
10575: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1]
26127: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1]
4095: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1]
19647: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1]
423: Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1]
15975: Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1]
3015: Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1]
18567: Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1]
8415: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1]
23967: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1]
1935: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1]
17487: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1]
11007: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1]
26559: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1]
4527: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1]
20079: Qb[p, 1] Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1]
855: Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1]
16407: Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1]
3447: Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1]
18999: Qb[p, 1] Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1]
32433: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35025: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
12993: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28545: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31353: Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33945: Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5433: Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
20985: Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32865: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35457: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13425: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
28977: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31785: Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34377: Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5865: Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21417: Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8673: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24225: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2193: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17745: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11265: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26817: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4785: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20337: Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
31173: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2]
33765: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2]
5253: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1]
20805: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1]
31605: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2]
34197: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2]
5685: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1]
21237: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1]
933: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1]
16485: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1]
3525: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1]
19077: Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1]
32505: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35097: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13065: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28617: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31425: Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
34017: Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5505: Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21057: Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32937: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35529: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13497: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29049: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31857: Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34449: Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5937: Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21489: Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8745: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24297: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2265: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17817: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11337: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26889: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4857: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20409: Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
31245: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2]
33837: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2]
5325: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1]
20877: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1]
31677: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2]
34269: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2]
5757: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1]
21309: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1]
1005: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1]
16557: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1]
3597: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1]
19149: Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1]
7953: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1]
23505: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1]
1473: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1]
17025: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1]
10545: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1]
26097: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1]
4065: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1]
19617: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1]
393: Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1]
15945: Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1]
2985: Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1]
18537: Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1]
8385: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1]
23937: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1]
1905: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1]
17457: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1]
10977: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1]
26529: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1]
4497: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1]
20049: Qb[p, 1] Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1]
825: Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1]
16377: Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1]
3417: Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1]
18969: Qb[p, 1] Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1]
40396: Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2]
42988: Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14476: Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30028: Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40828: Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43420: Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2]
14908: Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30460: Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10156: Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25708: Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12748: Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28300: Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1]
40216: Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42808: Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14296: Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29848: Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39136: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41728: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6736: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22288: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40648: Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43240: Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14728: Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30280: Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39568: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42160: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7168: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22720: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8896: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24448: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2416: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
17968: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11488: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27040: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5008: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20560: Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
40468: Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43060: Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14548: Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30100: Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40900: Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43492: Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2]
14980: Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30532: Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10228: Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25780: Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12820: Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28372: Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1]
40288: Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42880: Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14368: Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29920: Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39208: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41800: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6808: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22360: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40720: Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43312: Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14800: Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30352: Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39640: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42232: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7240: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22792: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8968: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24520: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2488: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
18040: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11560: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27112: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5080: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20632: Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
9256: Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1]
24808: Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1]
11848: Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1]
27400: Q[m, 2] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1]
8176: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1]
23728: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1]
1696: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1]
17248: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1]
10768: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1]
26320: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1]
4288: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1]
19840: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1]
9688: Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1]
25240: Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1]
12280: Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1]
27832: Q[m, 2] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1]
8608: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1]
24160: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1]
2128: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1]
17680: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1]
11200: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1]
26752: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1]
4720: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1]
20272: Q[m, 2] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1]
40186: Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42778: Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14266: Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29818: Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39106: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41698: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6706: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22258: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40618: Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43210: Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14698: Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30250: Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39538: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42130: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7138: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22690: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8866: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24418: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2386: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
17938: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11458: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27010: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
4978: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20530: Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
38926: Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2]
41518: Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2]
6526: Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1]
22078: Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1]
39358: Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2]
41950: Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2]
6958: Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1]
22510: Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1]
1126: Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1]
16678: Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1]
3718: Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1]
19270: Q[m, 2] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1]
40258: Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42850: Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14338: Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29890: Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39178: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41770: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6778: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22330: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40690: Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43282: Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14770: Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30322: Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39610: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42202: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7210: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22762: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8938: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24490: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2458: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
18010: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11530: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27082: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5050: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20602: Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
38998: Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2]
41590: Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2]
6598: Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1]
22150: Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1]
39430: Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2]
42022: Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2]
7030: Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1]
22582: Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1]
1198: Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1]
16750: Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1]
3790: Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1]
19342: Q[m, 2] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1]
7966: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1]
23518: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1]
1486: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1]
17038: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1]
10558: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1]
26110: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1]
4078: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1]
19630: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1]
406: Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1]
15958: Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1]
2998: Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1]
18550: Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1]
8398: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1]
23950: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1]
1918: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1]
17470: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1]
10990: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1]
26542: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1]
4510: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1]
20062: Q[m, 2] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1]
838: Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1]
16390: Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1]
3430: Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1]
18982: Q[m, 2] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1]
40408: Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43000: Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14488: Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30040: Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40840: Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43432: Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2]
14920: Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30472: Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10168: Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25720: Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12760: Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28312: Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1]
40228: Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42820: Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14308: Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29860: Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39148: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41740: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6748: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22300: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40660: Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43252: Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14740: Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30292: Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39580: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42172: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7180: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22732: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8908: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24460: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2428: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
17980: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11500: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27052: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5020: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20572: Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
40480: Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[m, 1] Qb[p, 2]
43072: Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 1] Qb[p, 2]
14560: Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[m, 1]
30112: Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Qb[p, 1]
40912: Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2]
43504: Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2]
14992: Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1]
30544: Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1]
10240: Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1]
25792: Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1]
12832: Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1]
28384: Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1]
40300: Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42892: Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14380: Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29932: Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39220: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41812: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6820: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22372: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40732: Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43324: Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14812: Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30364: Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39652: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42244: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7252: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22804: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8980: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24532: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2500: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
18052: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11572: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27124: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5092: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20644: Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
9268: Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1]
24820: Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1]
11860: Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1]
27412: Q[m, 2] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1]
8188: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1]
23740: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1]
1708: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1]
17260: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1]
10780: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1]
26332: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1]
4300: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1]
19852: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1]
9700: Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1]
25252: Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1]
12292: Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1]
27844: Q[m, 2] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1]
8620: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1]
24172: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1]
2140: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1]
17692: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1]
11212: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1]
26764: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1]
4732: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1]
20284: Q[m, 2] Q[p, 1] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1]
40198: Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42790: Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14278: Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29830: Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39118: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41710: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6718: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22270: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40630: Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43222: Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14710: Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30262: Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39550: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42142: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7150: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22702: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8878: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24430: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2398: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
17950: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11470: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27022: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
4990: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20542: Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
38938: Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 2]
41530: Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 2]
6538: Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[m, 1]
22090: Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2] Q[p, 1]
39370: Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 2]
41962: Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 2]
6970: Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[m, 1]
22522: Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2] Q[p, 1]
1138: Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[m, 1]
16690: Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1] Q[p, 1]
3730: Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[m, 1]
19282: Q[m, 2] Qb[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1] Q[p, 1]
40270: Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 2]
42862: Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 2]
14350: Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[m, 1]
29902: Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2] Qb[p, 1]
39190: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 2]
41782: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 2]
6790: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[m, 1]
22342: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2] Q[p, 1]
40702: Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 2]
43294: Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2]
14782: Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1]
30334: Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1]
39622: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2]
42214: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2]
7222: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1]
22774: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1]
8950: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1]
24502: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1]
2470: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1]
18022: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1]
11542: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1]
27094: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1]
5062: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1]
20614: Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1]
39010: Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 2]
41602: Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 2]
6610: Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[m, 1]
22162: Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 2] Q[p, 1]
39442: Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 2]
42034: Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 2]
7042: Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[m, 1]
22594: Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 2] Q[p, 1]
1210: Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1]
16762: Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1]
3802: Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1]
19354: Q[m, 2] Qb[p, 1] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1]
7978: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1]
23530: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1]
1498: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1]
17050: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1]
10570: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1]
26122: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1]
4090: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1]
19642: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1]
418: Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1]
15970: Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1]
3010: Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1]
18562: Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1]
8410: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1]
23962: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1]
1930: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1]
17482: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1]
11002: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1]
26554: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1]
4522: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1]
20074: Q[m, 2] Qb[p, 1] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1]
850: Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1]
16402: Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1]
3442: Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1]
18994: Q[m, 2] Qb[p, 1] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1]
9106: Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[m, 1]
24658: Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 1]
11698: Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[m, 1]
27250: Q[m, 2] Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 1]
8026: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[m, 1]
23578: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 1]
1546: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 1]
17098: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[p, 1]
10618: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[m, 1]
26170: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 1]
4138: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 1]
19690: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[p, 1]
9538: Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[m, 1]
25090: Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 1]
12130: Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[m, 1]
27682: Q[m, 2] Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 1]
8458: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[m, 1]
24010: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 1]
1978: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 1]
17530: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[p, 1]
11050: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[m, 1]
26602: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 1]
4570: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 1]
20122: Q[m, 2] Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[p, 1]
7846: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[m, 1]
23398: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 1]
1366: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 1]
16918: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[p, 1]
10438: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[m, 1]
25990: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 1]
3958: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 1]
19510: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[p, 1]
286: Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 1]
15838: Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[p, 1]
2878: Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 1]
18430: Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[p, 1]
8278: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[m, 1]
23830: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 1]
1798: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 1]
17350: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[p, 1]
10870: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[m, 1]
26422: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 1]
4390: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 1]
19942: Q[m, 2] Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[p, 1]
718: Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 1]
16270: Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[p, 1]
3310: Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 1]
18862: Q[m, 2] Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[p, 1]
9178: Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[m, 1]
24730: Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 1]
11770: Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[m, 1]
27322: Q[m, 2] Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 1]
8098: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[m, 1]
23650: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 1]
1618: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 1]
17170: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[p, 1]
10690: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[m, 1]
26242: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 1]
4210: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 1]
19762: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[p, 1]
9610: Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[m, 1]
25162: Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 1]
12202: Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[m, 1]
27754: Q[m, 2] Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 1]
8530: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[m, 1]
24082: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 1]
2050: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 1]
17602: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[p, 1]
11122: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[m, 1]
26674: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 1]
4642: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 1]
20194: Q[m, 2] Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[p, 1]
7918: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[m, 1]
23470: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 1]
1438: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 1]
16990: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[p, 1]
10510: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[m, 1]
26062: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 1]
4030: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 1]
19582: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[p, 1]
358: Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 1]
15910: Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[p, 1]
2950: Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 1]
18502: Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[p, 1]
8350: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[m, 1]
23902: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 1]
1870: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 1]
17422: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[p, 1]
10942: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[m, 1]
26494: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 1]
4462: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 1]
20014: Q[m, 2] Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[p, 1]
790: Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 1]
16342: Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[p, 1]
3382: Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 1]
18934: Q[m, 2] Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[p, 1]
32621: Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2]
35213: Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2]
13181: Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1]
28733: Qb[p, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1]
33053: Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2]
35645: Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2]
13613: Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1]
29165: Qb[p, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1]
9941: Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1]
25493: Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1]
12533: Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1]
28085: Qb[p, 2] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1]
32441: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35033: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13001: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28553: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31361: Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33953: Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5441: Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
20993: Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32873: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35465: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13433: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
28985: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31793: Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34385: Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5873: Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21425: Qb[p, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8681: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24233: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2201: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17753: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11273: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26825: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4793: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20345: Qb[p, 2] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
32693: Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2]
35285: Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2]
13253: Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1]
28805: Qb[p, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1]
33125: Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2]
35717: Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2]
13685: Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1]
29237: Qb[p, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1]
10013: Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1]
25565: Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1]
12605: Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1]
28157: Qb[p, 2] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1]
32513: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35105: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13073: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28625: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31433: Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
34025: Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5513: Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21065: Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32945: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35537: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13505: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29057: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31865: Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34457: Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5945: Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21497: Qb[p, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8753: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24305: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2273: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17825: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11345: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26897: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4865: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20417: Qb[p, 2] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
9221: Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1]
24773: Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1]
11813: Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1]
27365: Qb[p, 2] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1]
8141: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1]
23693: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1]
1661: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1]
17213: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1]
10733: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1]
26285: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1]
4253: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1]
19805: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1]
9653: Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1]
25205: Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1]
12245: Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1]
27797: Qb[p, 2] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1]
8573: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1]
24125: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1]
2093: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1]
17645: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1]
11165: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1]
26717: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1]
4685: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1]
20237: Qb[p, 2] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1]
32411: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35003: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
12971: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28523: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31331: Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33923: Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5411: Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
20963: Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32843: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35435: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13403: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
28955: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31763: Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34355: Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5843: Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21395: Qb[p, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8651: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24203: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2171: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17723: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11243: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26795: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4763: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20315: Qb[p, 2] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
31151: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2]
33743: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2]
5231: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1]
20783: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1]
31583: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2]
34175: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2]
5663: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1]
21215: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1]
911: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1]
16463: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1]
3503: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1]
19055: Qb[p, 2] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1]
32483: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35075: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13043: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28595: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31403: Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33995: Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5483: Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21035: Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32915: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35507: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13475: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29027: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31835: Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34427: Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5915: Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21467: Qb[p, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8723: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24275: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2243: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17795: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11315: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26867: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4835: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20387: Qb[p, 2] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
31223: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2]
33815: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2]
5303: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1]
20855: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1]
31655: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2]
34247: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2]
5735: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1]
21287: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1]
983: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1]
16535: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1]
3575: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1]
19127: Qb[p, 2] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1]
7931: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1]
23483: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1]
1451: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1]
17003: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1]
10523: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1]
26075: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1]
4043: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1]
19595: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1]
371: Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1]
15923: Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1]
2963: Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1]
18515: Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1]
8363: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1]
23915: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1]
1883: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1]
17435: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1]
10955: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1]
26507: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1]
4475: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1]
20027: Qb[p, 2] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1]
803: Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1]
16355: Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1]
3395: Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1]
18947: Qb[p, 2] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1]
32633: Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 2]
35225: Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 2]
13193: Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[m, 1]
28745: Qb[p, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Qb[p, 1]
33065: Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 2]
35657: Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 2]
13625: Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[m, 1]
29177: Qb[p, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Qb[p, 1]
9953: Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[m, 1]
25505: Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Qb[p, 1]
12545: Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[m, 1]
28097: Qb[p, 2] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Qb[p, 1]
32453: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35045: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13013: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28565: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31373: Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33965: Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5453: Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21005: Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32885: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35477: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13445: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
28997: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31805: Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34397: Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5885: Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21437: Qb[p, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8693: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24245: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2213: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17765: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11285: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26837: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4805: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20357: Qb[p, 2] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
32705: Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 2]
35297: Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 2]
13265: Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[m, 1]
28817: Qb[p, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Qb[p, 1]
33137: Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2]
35729: Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2]
13697: Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1]
29249: Qb[p, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1]
10025: Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1]
25577: Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1]
12617: Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1]
28169: Qb[p, 2] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1]
32525: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35117: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13085: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28637: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31445: Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
34037: Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5525: Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21077: Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32957: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35549: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13517: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29069: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31877: Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34469: Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5957: Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21509: Qb[p, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8765: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24317: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2285: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17837: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11357: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26909: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4877: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20429: Qb[p, 2] Q[p, 1] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
9233: Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1]
24785: Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1]
11825: Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1]
27377: Qb[p, 2] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1]
8153: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1]
23705: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1]
1673: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1]
17225: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1]
10745: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1]
26297: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1]
4265: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1]
19817: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1]
9665: Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1]
25217: Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1]
12257: Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1]
27809: Qb[p, 2] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1]
8585: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1]
24137: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1]
2105: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1]
17657: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1]
11177: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1]
26729: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1]
4697: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1]
20249: Qb[p, 2] Q[p, 1] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1]
32423: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35015: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
12983: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28535: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31343: Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
33935: Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5423: Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
20975: Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32855: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35447: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13415: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
28967: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31775: Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34367: Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5855: Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21407: Qb[p, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8663: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24215: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2183: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17735: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11255: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26807: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4775: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20327: Qb[p, 2] Qb[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
31163: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 1] Q[m, 2]
33755: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[p, 1] Q[m, 2]
5243: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Q[m, 1]
20795: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2] Q[p, 1]
31595: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 1] Q[m, 2]
34187: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[p, 1] Q[m, 2]
5675: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Q[m, 1]
21227: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2] Q[p, 1]
923: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Q[m, 1]
16475: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1] Q[p, 1]
3515: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Q[m, 1]
19067: Qb[p, 2] Qb[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1] Q[p, 1]
32495: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[m, 1] Q[m, 2]
35087: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 1] Q[m, 2]
13055: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[m, 1]
28607: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2] Qb[p, 1]
31415: Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 1] Q[m, 2]
34007: Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[p, 1] Q[m, 2]
5495: Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Q[m, 1]
21047: Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2] Q[p, 1]
32927: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[m, 1] Q[m, 2]
35519: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2]
13487: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1]
29039: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1]
31847: Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2]
34439: Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2]
5927: Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1]
21479: Qb[p, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1]
8735: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1]
24287: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1]
2255: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1]
17807: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1]
11327: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1]
26879: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1]
4847: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1]
20399: Qb[p, 2] Qb[p, 1] Q[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1]
31235: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 1] Q[m, 2]
33827: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[p, 1] Q[m, 2]
5315: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Q[m, 1]
20867: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 2] Q[p, 1]
31667: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 1] Q[m, 2]
34259: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[p, 1] Q[m, 2]
5747: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Q[m, 1]
21299: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 2] Q[p, 1]
995: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1]
16547: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1]
3587: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1]
19139: Qb[p, 2] Qb[p, 1] Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1]
7943: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1]
23495: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1]
1463: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1]
17015: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1]
10535: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1]
26087: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1]
4055: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1]
19607: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1]
383: Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1]
15935: Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1]
2975: Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1]
18527: Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1]
8375: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1]
23927: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1]
1895: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1]
17447: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1]
10967: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1]
26519: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1]
4487: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1]
20039: Qb[p, 2] Qb[p, 1] Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1]
815: Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1]
16367: Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1]
3407: Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1]
18959: Qb[p, 2] Qb[p, 1] Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1]
9101: Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[m, 1]
24653: Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1] Qb[m, 1] Qb[p, 1]
11693: Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[m, 1]
27245: Qb[p, 2] Q[m, 2] Q[m, 1] Q[m, 1] Qb[p, 1] Qb[p, 1]
8021: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[m, 1]
23573: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Q[m, 1] Qb[p, 1]
1541: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[m, 1]
17093: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[m, 1] Q[p, 1]
10613: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[m, 1]
26165: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Q[p, 1] Qb[p, 1]
4133: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[m, 1]
19685: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[m, 1] Qb[p, 1] Q[p, 1]
9533: Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[m, 1]
25085: Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Qb[m, 1] Qb[p, 1]
12125: Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[m, 1]
27677: Qb[p, 2] Q[m, 2] Q[m, 1] Q[p, 1] Qb[p, 1] Qb[p, 1]
8453: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[m, 1]
24005: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Q[m, 1] Qb[p, 1]
1973: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[m, 1]
17525: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[m, 1] Q[p, 1]
11045: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[m, 1]
26597: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Q[p, 1] Qb[p, 1]
4565: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[m, 1]
20117: Qb[p, 2] Q[m, 2] Q[m, 1] Qb[p, 1] Qb[p, 1] Q[p, 1]
7841: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[m, 1]
23393: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Q[m, 1] Qb[p, 1]
1361: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[m, 1]
16913: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[m, 1] Q[p, 1]
10433: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[m, 1]
25985: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Q[p, 1] Qb[p, 1]
3953: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[m, 1]
19505: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[m, 1] Qb[p, 1] Q[p, 1]
281: Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[m, 1]
15833: Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[m, 1] Q[p, 1]
2873: Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[m, 1]
18425: Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[m, 1] Q[p, 1] Q[p, 1]
8273: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[m, 1]
23825: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Q[m, 1] Qb[p, 1]
1793: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[m, 1]
17345: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[m, 1] Q[p, 1]
10865: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[m, 1]
26417: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Q[p, 1] Qb[p, 1]
4385: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[m, 1]
19937: Qb[p, 2] Q[m, 2] Qb[m, 1] Q[p, 1] Qb[p, 1] Q[p, 1]
713: Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[m, 1]
16265: Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[m, 1] Q[p, 1]
3305: Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[m, 1]
18857: Qb[p, 2] Q[m, 2] Qb[m, 1] Qb[p, 1] Q[p, 1] Q[p, 1]
9173: Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[m, 1]
24725: Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Qb[m, 1] Qb[p, 1]
11765: Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[m, 1]
27317: Qb[p, 2] Q[m, 2] Q[p, 1] Q[m, 1] Qb[p, 1] Qb[p, 1]
8093: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[m, 1]
23645: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Q[m, 1] Qb[p, 1]
1613: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[m, 1]
17165: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[m, 1] Q[p, 1]
10685: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[m, 1]
26237: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Q[p, 1] Qb[p, 1]
4205: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[m, 1]
19757: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[m, 1] Qb[p, 1] Q[p, 1]
9605: Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[m, 1]
25157: Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1] Qb[m, 1] Qb[p, 1]
12197: Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[m, 1]
27749: Qb[p, 2] Q[m, 2] Q[p, 1] Q[p, 1] Qb[p, 1] Qb[p, 1]
8525: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[m, 1]
24077: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Q[m, 1] Qb[p, 1]
2045: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[m, 1]
17597: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[m, 1] Q[p, 1]
11117: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[m, 1]
26669: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Q[p, 1] Qb[p, 1]
4637: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[m, 1]
20189: Qb[p, 2] Q[m, 2] Q[p, 1] Qb[p, 1] Qb[p, 1] Q[p, 1]
7913: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[m, 1]
23465: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Q[m, 1] Qb[p, 1]
1433: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[m, 1]
16985: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[m, 1] Q[p, 1]
10505: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[m, 1]
26057: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Q[p, 1] Qb[p, 1]
4025: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[m, 1]
19577: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[m, 1] Qb[p, 1] Q[p, 1]
353: Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[m, 1]
15905: Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[m, 1] Q[p, 1]
2945: Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[m, 1]
18497: Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[m, 1] Q[p, 1] Q[p, 1]
8345: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[m, 1]
23897: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Q[m, 1] Qb[p, 1]
1865: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[m, 1]
17417: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[m, 1] Q[p, 1]
10937: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[m, 1]
26489: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Q[p, 1] Qb[p, 1]
4457: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[m, 1]
20009: Qb[p, 2] Q[m, 2] Qb[p, 1] Q[p, 1] Qb[p, 1] Q[p, 1]
785: Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[m, 1]
16337: Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[m, 1] Q[p, 1]
3377: Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[m, 1]
18929: Qb[p, 2] Q[m, 2] Qb[p, 1] Qb[p, 1] Q[p, 1] Q[p, 1]
*/

