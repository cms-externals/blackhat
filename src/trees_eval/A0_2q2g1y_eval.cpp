/* The 2 quark Ngluon 1 photon amplitudes use a base 6 code to enumerate them
   m,qm,qp,p,gam,gap corresponds to 0,1,2,3,5.
   For example, (qp,qm,m,p,gam,gap) = 2+ 1*6 + 0*36 + 3*216 + 4*1296 + 5*7776
  See bottom of file for table of switch numbers and values
*/

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"
#include <iostream>

using namespace std;

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i1,i2,i3) ep.s(ind[i1-1],ind[i2-1], ind[i3-1])
#define SSS(i1,i2,i3,i4) ep.s(ind[i1-1],ind[i2-1], ind[i3-1], ind[i4-1])
#define SSSS(i1,i2,i3,i4,i5) ep.s(ind[i1-1],ind[i2-1], ind[i3-1], ind[i4-1], ind[i5-1])
#define SSSSS(i1,i2,i3,i4,i5,i6) ep.s(ind[i1-1],ind[i2-1], ind[i3-1], ind[i4-1], ind[i5-1], ind[i6-1], )

namespace BH {


template <class T> complex<T> A2q2g1yzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q2g1y62_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y68_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(1,0),2))/(ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y97_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y103_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y152_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y188_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),3)*ep.spb(2,1))/
(ep.spb(1,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(2,0),3))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(2,1))/
(ep.spb(1,0)*ep.spb(3,1)*
 ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y242_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y248_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),3))/(ep.spb(2,0)*
 ep.spb(2,1)*ep.spb(4,0)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(1,0),2)*
 ep.spb(3,1))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,0)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spb(2,0),2)*
 ep.spb(3,2))/(ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y356_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,2)*
 ep.spa(1,2)*ep.spa(1,3))-
 (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g1y362_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),3))/(ep.spb(2,0)*
 ep.spb(2,1)*ep.spb(4,0)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(1,0),2)*
 ep.spb(3,1))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y398_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,0)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spb(2,0),2)*
 ep.spb(3,2))/(ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y416_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,2)*
 ep.spa(1,2)*ep.spa(1,3))-
 (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g1y457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y463_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(1,0))/
(ep.spb(2,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(3,2),3))/(ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y571_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,2)*
 ep.spa(1,2)*ep.spa(1,3))-
 (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g1y577_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(1,0))/
(ep.spb(2,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(3,2),3))/(ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y631_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,2)*
 ep.spa(1,2)*ep.spa(1,3))-
 (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g1y710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y716_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y751_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(3,0),3)*
 ep.spb(3,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
 ep.spb(2,0)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y836_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),3))/(ep.spa(0,1)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(2,4))-(complex<T>(0,1)*pow(ep.spa(1,4),3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(1,4),3)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,3)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y872_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(1,0),2))/(ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(3,0)*
 ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(3,0)*
 ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),3)*ep.spb(2,1))/
(ep.spb(1,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(2,0),3))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
 ep.spb(3,0)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y1088_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(3,0),3)*
 ep.spb(3,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y1118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,0)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(3,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y1136_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,3)*
 ep.spa(1,2)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y1153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(3,2),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y1171_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3))/(ep.spa(0,1)*
 ep.spa(0,3)*ep.spa(1,2)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 ep.spa(2,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y1196_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),3))/(ep.spa(0,1)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(2,4))-(complex<T>(0,1)*pow(ep.spa(1,4),3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(1,4),3)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,3)*
 ep.spa(1,2)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y2617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y2623_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(1,0))/
(ep.spb(2,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,1),3)*ep.spb(1,0))/
(ep.spb(2,0)*ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(4,1),3))/
(ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y2725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,1))+(complex<T>(0,1)*pow(ep.spb(4,2),3))/
(ep.spb(2,1)*ep.spb(3,1)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(4,2),3)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y2731_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/(ep.spa(0,2)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(1,3))+(complex<T>(0,1)*pow(ep.spa(0,3),3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(0,3),3)*ep.spa(3,4))/
(ep.spa(0,2)*ep.spa(0,4)*
 ep.spa(1,3)*ep.spa(1,4)*
 ep.spa(2,3))
); }

template <class T> complex<T> A2q2g1y2737_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y2755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),3))/(ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(4,1),3)*
 ep.spb(1,0))/(ep.spb(2,0)*
 ep.spb(2,1)*ep.spb(3,1)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(4,1),3))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y2773_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(4,2),3)*
 ep.spb(2,0))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y2791_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/(ep.spa(0,2)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(1,3))+(complex<T>(0,1)*pow(ep.spa(0,3),3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(0,3),3)*ep.spa(3,4))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(2,4))
); }

template <class T> complex<T> A2q2g1y3265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(3,0))/
(ep.spb(2,0)*ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,3),3)*ep.spb(3,0))/
(ep.spb(2,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,1))
); }

template <class T> complex<T> A2q2g1y3271_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
(ep.spa(0,4)*ep.spa(1,3)*
 ep.spa(1,4)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
(ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(1,3)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(2,4))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y3373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),3))/(ep.spa(0,2)*
 ep.spa(0,4)*ep.spa(1,3)*
 ep.spa(2,3))-(complex<T>(0,1)*pow(ep.spa(0,1),3)*
 ep.spa(1,4))/(ep.spa(0,2)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(1,3)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(1,4))/
(ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y3379_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y3385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(3,0))/
(ep.spb(2,0)*ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,3),3)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,2))
); }

template <class T> complex<T> A2q2g1y3403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,3)*
 ep.spa(2,3))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(2,4))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 ep.spa(2,4))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y3421_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*ep.spa(1,4))/
(ep.spa(0,4)*ep.spa(1,3)*
 ep.spa(2,3)*ep.spa(2,4))-
 (complex<T>(0,1)*pow(ep.spa(0,1),3)*ep.spa(1,4))/
(ep.spa(0,2)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 ep.spa(1,4))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y3439_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y3457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y3475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),3)*ep.spb(1,0))/
(ep.spb(2,1)*ep.spb(3,0)*
 ep.spb(3,1)*ep.spb(4,0)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(4,1),3))/
(ep.spb(3,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,2))-
 (complex<T>(0,1)*pow(ep.spb(4,1),3))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y3565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
(ep.spb(2,1)*ep.spb(3,0)*
 ep.spb(3,1)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(4,2),3)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y3583_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(2,3))-(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(3,4))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(2,4))-(complex<T>(0,1)*pow(ep.spa(0,3),3)*
 ep.spa(3,4))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,3)*
 ep.spa(2,3)*ep.spa(2,4))
); }

template <class T> complex<T> A2q2g1y3673_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,3),3))/(ep.spb(2,1)*
 ep.spb(3,1)*ep.spb(4,0)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(4,3),3)*
 ep.spb(3,0))/(ep.spb(1,0)*
 ep.spb(3,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,2))
); }

template <class T> complex<T> A2q2g1y3691_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3))/(ep.spa(0,3)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(1,3))-(complex<T>(0,1)*pow(ep.spa(0,2),3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,3)*ep.spa(2,3))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(2,4))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y3781_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),3)*ep.spa(1,4))/
(ep.spa(0,3)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(2,4))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 ep.spa(1,4))/(ep.spa(0,4)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(2,4))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 ep.spa(1,4))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y3799_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y3950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y3956_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y3985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y3991_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,2))+
 (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(4,0),3)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(2,0)*
 ep.spb(3,1)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y4076_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(2,4))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 ep.spa(0,3))/(ep.spa(0,1)*
 ep.spa(0,2)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/(ep.spb(2,0)*
 ep.spb(2,1)*ep.spb(3,1))+
 (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g1y4136_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(2,3),3))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4238_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,1))/
(ep.spa(0,2)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4244_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y4250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/(ep.spb(2,0)*
 ep.spb(2,1)*ep.spb(3,1))+
 (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g1y4268_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(2,3),3))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4286_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,1))/
(ep.spa(0,2)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4304_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y4345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/(ep.spb(2,0)*
 ep.spb(2,1)*ep.spb(3,1))+
 (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g1y4351_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(2,3))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4453_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),3))/(ep.spa(0,2)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 ep.spa(1,3))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4459_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y4465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/(ep.spb(2,0)*
 ep.spb(2,1)*ep.spb(3,1))+
 (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g1y4483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(2,3))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4501_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),3))/(ep.spa(0,2)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 ep.spa(1,3))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4519_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y4598_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,4)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4604_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y4633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4639_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y4688_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),3)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(2,4))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 ep.spa(0,2))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,2),3))/
(ep.spa(0,1)*ep.spa(1,3)*
 ep.spa(2,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4724_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y4760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(3,1)*
 ep.spb(3,2)*ep.spb(4,2))+
 (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(4,0),3)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,0)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y4790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,0)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spb(4,0),2)*
 ep.spb(4,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y4808_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(2,3),3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(4,2),2)*
 ep.spb(4,0))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y4843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(2,3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4868_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(2,4))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,3),3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y4976_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),3)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(2,4))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 ep.spa(0,2))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 ep.spa(0,2))/(ep.spa(0,1)*
 ep.spa(0,3)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y5006_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,3)*
 ep.spa(2,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y5024_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y5041_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
 ep.spa(2,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y5059_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y5084_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y5192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y5222_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y5240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(1,0),2))/(ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(3,0)*
 ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5257_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y5275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(3,0)*
 ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),3)*ep.spb(2,1))/
(ep.spb(1,0)*ep.spb(3,0)*
 ep.spb(3,2)*ep.spb(4,1)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(2,0),3))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
 ep.spb(3,0)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y5420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(1,0),2)*ep.spb(3,1))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5528_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g1y5617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y5635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5743_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g1y5840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,1)*ep.spb(4,2))+
 (complex<T>(0,1)*pow(ep.spb(3,0),3)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(3,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,0)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(3,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5888_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,3)*
 ep.spa(1,2)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y5905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(3,2),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y5923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3))/(ep.spa(0,1)*
 ep.spa(0,3)*ep.spa(1,2)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 ep.spa(2,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y5948_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,3)*
 ep.spa(2,3)*ep.spa(2,4))-
 (complex<T>(0,1)*pow(ep.spa(1,4),3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,4),3)*
 ep.spa(0,4))/(ep.spa(0,1)*
 ep.spa(0,3)*ep.spa(1,2)*
 ep.spa(2,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y6488_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/(ep.spb(1,0)*
 ep.spb(3,0)*ep.spb(3,2)*
 ep.spb(4,2))+(complex<T>(0,1)*pow(ep.spb(4,0),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(4,0),3)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,0)*
 ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y6518_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,0)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spb(4,0),2)*
 ep.spb(4,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y6536_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(2,3),3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y6553_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spb(4,2),2)*
 ep.spb(4,0))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,0)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g1y6571_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(2,3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y6596_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
 ep.spa(1,4)*ep.spa(2,3)*
 ep.spa(2,4))-(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 ep.spa(0,3))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y6698_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g1y6716_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y6806_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y6824_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y6913_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g1y6931_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y7021_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(1,3))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y7039_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y7136_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),3)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,3)*
 ep.spa(1,4)*ep.spa(2,3)*
 ep.spa(2,4))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 ep.spa(0,2))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 ep.spa(0,2))/(ep.spa(0,1)*
 ep.spa(0,3)*ep.spa(2,4)*
 ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y7166_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,3)*
 ep.spa(2,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y7184_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y7201_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
 ep.spa(2,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g1y7219_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g1y7244_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
 ); }


template <class T> complex<T>  (*A2q2g1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 62 : return &A2q2g1y62_eval;
case 68 : return &A2q2g1y68_eval;
case 97 : return &A2q2g1y97_eval;
case 103 : return &A2q2g1y103_eval;
case 152 : return &A2q2g1y152_eval;
case 188 : return &A2q2g1y188_eval;
case 242 : return &A2q2g1y242_eval;
case 248 : return &A2q2g1y248_eval;
case 350 : return &A2q2g1y350_eval;
case 356 : return &A2q2g1y356_eval;
case 362 : return &A2q2g1y362_eval;
case 380 : return &A2q2g1y380_eval;
case 398 : return &A2q2g1y398_eval;
case 416 : return &A2q2g1y416_eval;
case 457 : return &A2q2g1y457_eval;
case 463 : return &A2q2g1y463_eval;
case 565 : return &A2q2g1y565_eval;
case 571 : return &A2q2g1y571_eval;
case 577 : return &A2q2g1y577_eval;
case 595 : return &A2q2g1y595_eval;
case 613 : return &A2q2g1y613_eval;
case 631 : return &A2q2g1y631_eval;
case 710 : return &A2q2g1y710_eval;
case 716 : return &A2q2g1y716_eval;
case 745 : return &A2q2g1y745_eval;
case 751 : return &A2q2g1y751_eval;
case 800 : return &A2q2g1y800_eval;
case 836 : return &A2q2g1y836_eval;
case 872 : return &A2q2g1y872_eval;
case 902 : return &A2q2g1y902_eval;
case 920 : return &A2q2g1y920_eval;
case 937 : return &A2q2g1y937_eval;
case 955 : return &A2q2g1y955_eval;
case 980 : return &A2q2g1y980_eval;
case 1088 : return &A2q2g1y1088_eval;
case 1118 : return &A2q2g1y1118_eval;
case 1136 : return &A2q2g1y1136_eval;
case 1153 : return &A2q2g1y1153_eval;
case 1171 : return &A2q2g1y1171_eval;
case 1196 : return &A2q2g1y1196_eval;
case 2617 : return &A2q2g1y2617_eval;
case 2623 : return &A2q2g1y2623_eval;
case 2725 : return &A2q2g1y2725_eval;
case 2731 : return &A2q2g1y2731_eval;
case 2737 : return &A2q2g1y2737_eval;
case 2755 : return &A2q2g1y2755_eval;
case 2773 : return &A2q2g1y2773_eval;
case 2791 : return &A2q2g1y2791_eval;
case 3265 : return &A2q2g1y3265_eval;
case 3271 : return &A2q2g1y3271_eval;
case 3373 : return &A2q2g1y3373_eval;
case 3379 : return &A2q2g1y3379_eval;
case 3385 : return &A2q2g1y3385_eval;
case 3403 : return &A2q2g1y3403_eval;
case 3421 : return &A2q2g1y3421_eval;
case 3439 : return &A2q2g1y3439_eval;
case 3457 : return &A2q2g1y3457_eval;
case 3475 : return &A2q2g1y3475_eval;
case 3565 : return &A2q2g1y3565_eval;
case 3583 : return &A2q2g1y3583_eval;
case 3673 : return &A2q2g1y3673_eval;
case 3691 : return &A2q2g1y3691_eval;
case 3781 : return &A2q2g1y3781_eval;
case 3799 : return &A2q2g1y3799_eval;
case 3950 : return &A2q2g1y3950_eval;
case 3956 : return &A2q2g1y3956_eval;
case 3985 : return &A2q2g1y3985_eval;
case 3991 : return &A2q2g1y3991_eval;
case 4040 : return &A2q2g1y4040_eval;
case 4076 : return &A2q2g1y4076_eval;
case 4130 : return &A2q2g1y4130_eval;
case 4136 : return &A2q2g1y4136_eval;
case 4238 : return &A2q2g1y4238_eval;
case 4244 : return &A2q2g1y4244_eval;
case 4250 : return &A2q2g1y4250_eval;
case 4268 : return &A2q2g1y4268_eval;
case 4286 : return &A2q2g1y4286_eval;
case 4304 : return &A2q2g1y4304_eval;
case 4345 : return &A2q2g1y4345_eval;
case 4351 : return &A2q2g1y4351_eval;
case 4453 : return &A2q2g1y4453_eval;
case 4459 : return &A2q2g1y4459_eval;
case 4465 : return &A2q2g1y4465_eval;
case 4483 : return &A2q2g1y4483_eval;
case 4501 : return &A2q2g1y4501_eval;
case 4519 : return &A2q2g1y4519_eval;
case 4598 : return &A2q2g1y4598_eval;
case 4604 : return &A2q2g1y4604_eval;
case 4633 : return &A2q2g1y4633_eval;
case 4639 : return &A2q2g1y4639_eval;
case 4688 : return &A2q2g1y4688_eval;
case 4724 : return &A2q2g1y4724_eval;
case 4760 : return &A2q2g1y4760_eval;
case 4790 : return &A2q2g1y4790_eval;
case 4808 : return &A2q2g1y4808_eval;
case 4825 : return &A2q2g1y4825_eval;
case 4843 : return &A2q2g1y4843_eval;
case 4868 : return &A2q2g1y4868_eval;
case 4976 : return &A2q2g1y4976_eval;
case 5006 : return &A2q2g1y5006_eval;
case 5024 : return &A2q2g1y5024_eval;
case 5041 : return &A2q2g1y5041_eval;
case 5059 : return &A2q2g1y5059_eval;
case 5084 : return &A2q2g1y5084_eval;
case 5192 : return &A2q2g1y5192_eval;
case 5222 : return &A2q2g1y5222_eval;
case 5240 : return &A2q2g1y5240_eval;
case 5257 : return &A2q2g1y5257_eval;
case 5275 : return &A2q2g1y5275_eval;
case 5300 : return &A2q2g1y5300_eval;
case 5402 : return &A2q2g1y5402_eval;
case 5420 : return &A2q2g1y5420_eval;
case 5510 : return &A2q2g1y5510_eval;
case 5528 : return &A2q2g1y5528_eval;
case 5617 : return &A2q2g1y5617_eval;
case 5635 : return &A2q2g1y5635_eval;
case 5725 : return &A2q2g1y5725_eval;
case 5743 : return &A2q2g1y5743_eval;
case 5840 : return &A2q2g1y5840_eval;
case 5870 : return &A2q2g1y5870_eval;
case 5888 : return &A2q2g1y5888_eval;
case 5905 : return &A2q2g1y5905_eval;
case 5923 : return &A2q2g1y5923_eval;
case 5948 : return &A2q2g1y5948_eval;
case 6488 : return &A2q2g1y6488_eval;
case 6518 : return &A2q2g1y6518_eval;
case 6536 : return &A2q2g1y6536_eval;
case 6553 : return &A2q2g1y6553_eval;
case 6571 : return &A2q2g1y6571_eval;
case 6596 : return &A2q2g1y6596_eval;
case 6698 : return &A2q2g1y6698_eval;
case 6716 : return &A2q2g1y6716_eval;
case 6806 : return &A2q2g1y6806_eval;
case 6824 : return &A2q2g1y6824_eval;
case 6913 : return &A2q2g1y6913_eval;
case 6931 : return &A2q2g1y6931_eval;
case 7021 : return &A2q2g1y7021_eval;
case 7039 : return &A2q2g1y7039_eval;
case 7136 : return &A2q2g1y7136_eval;
case 7166 : return &A2q2g1y7166_eval;
case 7184 : return &A2q2g1y7184_eval;
case 7201 : return &A2q2g1y7201_eval;
case 7219 : return &A2q2g1y7219_eval;
case 7244 : return &A2q2g1y7244_eval;

default: return 0;

}
}



template complex<R>  (*A2q2g1y_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q2g1y_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q2g1y_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q2g1y_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
5617: qm m m qp gam
6913: qm m m qp gap
3457: qm m m gam qp
3673: qm m m gap qp
5257: qm m qp m gam
6553: qm m qp m gap
5905: qm m qp p gam
7201: qm m qp p gap
937: qm m qp gam m
4825: qm m qp gam p
1153: qm m qp gap m
5041: qm m qp gap p
5725: qm m p qp gam
7021: qm m p qp gap
3565: qm m p gam qp
3781: qm m p gap qp
2737: qm m gam m qp
577: qm m gam qp m
4465: qm m gam qp p
3385: qm m gam p qp
2773: qm m gap m qp
613: qm m gap qp m
4501: qm m gap qp p
3421: qm m gap p qp
5635: qm p m qp gam
6931: qm p m qp gap
3475: qm p m gam qp
3691: qm p m gap qp
5275: qm p qp m gam
6571: qm p qp m gap
5923: qm p qp p gam
7219: qm p qp p gap
955: qm p qp gam m
4843: qm p qp gam p
1171: qm p qp gap m
5059: qm p qp gap p
5743: qm p p qp gam
7039: qm p p qp gap
3583: qm p p gam qp
3799: qm p p gap qp
2755: qm p gam m qp
595: qm p gam qp m
4483: qm p gam qp p
3403: qm p gam p qp
2791: qm p gap m qp
631: qm p gap qp m
4519: qm p gap qp p
3439: qm p gap p qp
2617: qm gam m m qp
457: qm gam m qp m
4345: qm gam m qp p
3265: qm gam m p qp
97: qm gam qp m m
3985: qm gam qp m p
745: qm gam qp p m
4633: qm gam qp p p
2725: qm gam p m qp
565: qm gam p qp m
4453: qm gam p qp p
3373: qm gam p p qp
2623: qm gap m m qp
463: qm gap m qp m
4351: qm gap m qp p
3271: qm gap m p qp
103: qm gap qp m m
3991: qm gap qp m p
751: qm gap qp p m
4639: qm gap qp p p
2731: qm gap p m qp
571: qm gap p qp m
4459: qm gap p qp p
3379: qm gap p p qp
5402: qp m m qm gam
6698: qp m m qm gap
5222: qp m qm m gam
6518: qp m qm m gap
5870: qp m qm p gam
7166: qp m qm p gap
902: qp m qm gam m
4790: qp m qm gam p
1118: qp m qm gap m
5006: qp m qm gap p
5510: qp m p qm gam
6806: qp m p qm gap
362: qp m gam qm m
4250: qp m gam qm p
398: qp m gap qm m
4286: qp m gap qm p
5192: qp qm m m gam
6488: qp qm m m gap
5840: qp qm m p gam
7136: qp qm m p gap
872: qp qm m gam m
4760: qp qm m gam p
1088: qp qm m gap m
4976: qp qm m gap p
5300: qp qm p m gam
6596: qp qm p m gap
5948: qp qm p p gam
7244: qp qm p p gap
980: qp qm p gam m
4868: qp qm p gam p
1196: qp qm p gap m
5084: qp qm p gap p
152: qp qm gam m m
4040: qp qm gam m p
800: qp qm gam p m
4688: qp qm gam p p
188: qp qm gap m m
4076: qp qm gap m p
836: qp qm gap p m
4724: qp qm gap p p
5420: qp p m qm gam
6716: qp p m qm gap
5240: qp p qm m gam
6536: qp p qm m gap
5888: qp p qm p gam
7184: qp p qm p gap
920: qp p qm gam m
4808: qp p qm gam p
1136: qp p qm gap m
5024: qp p qm gap p
5528: qp p p qm gam
6824: qp p p qm gap
380: qp p gam qm m
4268: qp p gam qm p
416: qp p gap qm m
4304: qp p gap qm p
242: qp gam m qm m
4130: qp gam m qm p
62: qp gam qm m m
3950: qp gam qm m p
710: qp gam qm p m
4598: qp gam qm p p
350: qp gam p qm m
4238: qp gam p qm p
248: qp gap m qm m
4136: qp gap m qm p
68: qp gap qm m m
3956: qp gap qm m p
716: qp gap qm p m
4604: qp gap qm p p
356: qp gap p qm m
4244: qp gap p qm p
*/

