/*
*R_2q3g2l_num.cpp
*
* Created on 11/13, 2008
*/

#include "R_2q3g2l_eval.h"
#include "eval_param.h"

using namespace std;

namespace BH  {

#define _USE_OPTIMIZED 1
#define _VERBOSE 0

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)
#define SSS(i,j,k,l) ep.s(i-1,j-1,k-1,l-1)
#define SSSS(i,j,k,l,m) ep.s(i-1,j-1,k-1,l-1,m-1)


template <class T> complex<T> R2q3g2l_qpmmmqbmlmlbp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
//{{qp, m, m, m, qbm, lm, lbp}, lc}

#if _VERBOSE
  _MESSAGE("R2q3g2l_eval :  qpmmmqbmlmlbp lc");
#endif

complex<T> recursive;
recursive=complex<T>(1,0);
//new
#if _USE_OPTIMIZED

{
  complex<T> t1;
  complex<T> t11;
  complex<T> t115;
  complex<T> t116;
  complex<T> t119;
  complex<T> t12;
  complex<T> t120;
  complex<T> t125;
  complex<T> t128;
  complex<T> t131;
  complex<T> t133;
  complex<T> t134;
  complex<T> t135;
  complex<T> t136;
  complex<T> t137;
  complex<T> t14;
  complex<T> t144;
  complex<T> t147;
  complex<T> t15;
  complex<T> t150;
  complex<T> t151;
  complex<T> t153;
  complex<T> t155;
  complex<T> t156;
  complex<T> t160;
  complex<T> t162;
  complex<T> t166;
  complex<T> t17;
  complex<T> t174;
  complex<T> t176;
  complex<T> t177;
  complex<T> t178;
  complex<T> t180;
  complex<T> t181;
  complex<T> t182;
  complex<T> t189;
  complex<T> t19;
  complex<T> t194;
  complex<T> t197;
  complex<T> t198;
  complex<T> t199;
  complex<T> t2;
  complex<T> t200;
  complex<T> t208;
  complex<T> t209;
  complex<T> t21;
  complex<T> t210;
  complex<T> t212;
  complex<T> t214;
  complex<T> t22;
  complex<T> t222;
  complex<T> t23;
  complex<T> t233;
  complex<T> t24;
  complex<T> t242;
  complex<T> t25;
  complex<T> t260;
  complex<T> t265;
  complex<T> t270;
  complex<T> t271;
  complex<T> t275;
  complex<T> t278;
  complex<T> t279;
  complex<T> t285;
  complex<T> t291;
  complex<T> t296;
  complex<T> t299;
  complex<T> t3;
  complex<T> t309;
  complex<T> t31;
  complex<T> t311;
  complex<T> t315;
  complex<T> t33;
  complex<T> t333;
  complex<T> t334;
  complex<T> t336;
  complex<T> t338;
  complex<T> t34;
  complex<T> t357;
  complex<T> t359;
  complex<T> t36;
  complex<T> t364;
  complex<T> t371;
  complex<T> t372;
  complex<T> t377;
  complex<T> t379;
  complex<T> t38;
  complex<T> t385;
  complex<T> t4;
  complex<T> t40;
  complex<T> t413;
  complex<T> t414;
  complex<T> t420;
  complex<T> t421;
  complex<T> t434;
  complex<T> t435;
  complex<T> t436;
  complex<T> t437;
  complex<T> t439;
  complex<T> t44;
  complex<T> t440;
  complex<T> t446;
  complex<T> t447;
  complex<T> t45;
  complex<T> t474;
  complex<T> t475;
  complex<T> t486;
  complex<T> t494;
  complex<T> t5;
  complex<T> t51;
  complex<T> t53;
  complex<T> t54;
  complex<T> t56;
  complex<T> t57;
  complex<T> t6;
  complex<T> t61;
  complex<T> t62;
  complex<T> t64;
  complex<T> t67;
  complex<T> t69;
  complex<T> t70;
  complex<T> t71;
  complex<T> t72;
  complex<T> t74;
  complex<T> t75;
  complex<T> t8;
  complex<T> t83;
  complex<T> t84;
  complex<T> t85;
  complex<T> t87;
  complex<T> t89;
  complex<T> t9;
  complex<T> t92;
  complex<T> t97;
  complex<T> t99;
  {
    t1 = complex<T>(0,1);
    t2 = complex<T>(-2,0);
    t3 = SPA(2,3);
    t4 = SPB(4,3);
    t5 = t3*t4;
    t6 = SPB(5,3);
    t8 = SPA(1,6);
    t9 = SPB(2,1);
    t11 = SPA(6,7);
    t12 = SPB(7,2);
    t14 = t8*t9+t11*t12;
    t15 = SPB(3,1);
    t17 = SPB(7,3);
    t19 = t8*t15+t11*t17;
    t21 = SPB(3,2);
    t22 = pow(t21,2);
    t23 = complex<T>(1,0)/t22;
    t24 = SPB(5,4);
    t25 = complex<T>(1,0)/t24;
    t31 = SPA(2,4);
    t33 = SPA(3,4);
    t34 = t33*t19;
    t36 = SPB(4,1);
    t38 = SPB(7,4);
    t40 = t8*t36+t11*t38;
    t44 = complex<T>(1,0)/t21;
    t45 = complex<T>(1,0)/t14;
    t51 = pow(-t3*t19-t31*t40,2);
    t53 = SPA(1,7);
    t54 = SPB(7,1);
    t56 = SPA(5,6);
    t57 = SPB(7,5);
    t61 = complex<T>(-1,0);
    t62 = SPA(1,5);
    t64 = SPB(5,1);
    t67 = t8*t64+t11*t57;
    t69 = t61*t8;
    t70 = SPB(6,1);
    t71 = t8*t70;
    t72 = SPB(6,5);
    t74 = SPB(7,6);
    t75 = t11*t74;
    t83 = t53*t54;
    t84 = t71+t83+t75;
    t85 = complex<T>(1,0)/t84;
    t87 = complex<T>(-3,0);
    t89 = SPA(3,6);
    t92 = SPA(4,6);
    t97 = SPA(3,5);
    t99 = SPA(4,5);
    t115 = complex<T>(6,0);
    t116 = complex<T>(1,0)/t115;
    t119 = pow(t4,2);
    t120 = complex<T>(1,0)/t119;
    t125 = complex<T>(1,0)/(t8*t72+t53*t57);
    t128 = SPA(1,2);
    t131 = t128*t54-t3*t17;
    t133 = t128*t21;
    t134 = SPA(1,3);
    t135 = t134*t54;
    t136 = t3*t12;
    t137 = t135+t136;
    t144 = complex<T>(1,0)/t128;
    t147 = SPB(4,2);
    t150 = -t128*t147-t134*t4;
    t151 = complex<T>(1,0)/t150;
    t153 = complex<T>(1,0)/t74;
    t155 = SS(1,2,3);
    t156 = complex<T>(1,0)/t155;
    t160 = pow(t54,2);
    t162 = t87*t128;
    t166 = complex<T>(2,0);
    t174 = t97*t57;
    t176 = -t174-t89*t74;
    t177 = t21*t176;
    t178 = t99*t57;
    t180 = -t178-t92*t74;
    t181 = t147*t180;
    t182 = -t177-t181;
    t189 = SS(3,4,5);
    t194 = complex<T>(1,0)/t189;
    t197 = pow(t33,2);
    t198 = t197*t22;
    t199 = pow(t38,2);
    t200 = t21*t137;
    t208 = pow(-t97*t21-t99*t147,2);
    t209 = pow(t57,2);
    t210 = t208*t209;
    t212 = t61*t21;
    t214 = t61*t147;
    t222 = complex<T>(3,0);
    t233 = t33*t21;
    t242 = pow(-t89*t21-t92*t147,2);
    t260 = complex<T>(1,0)/t9;
    t265 = complex<T>(0,-1);
    t270 = -t62*t57-t8*t74;
    t271 = SPA(1,4);
    t275 = t61*t31;
    t278 = t3*t21;
    t279 = t33*t4;
    t285 = (-t133+t271*t4)*t54-t275*t4*t12-t61*(t278+t279)*t17-t275*t21*t38;
    t291 = t31*t147;
    t296 = t147*t17;
    t299 = (t134*t21+t271*t147)*t54-t61*(t278+t291)*t12-t61*t33*t296-t233*t38;
    t309 = t150*t54-t5*t12-t61*t3*t296-t61*(t291+t279)*t38;
    t311 = SPB(5,2);
    t315 = -t128*t311-t134*t6-t271*t24;
    t333 = t62*(t315*t54-t61*(-t3*t6-t31*t24)*t12-t61*(t3*t311-t33*t24)*t17-t61
*(t31*t311+t33*t6)*t38);
    t334 = SPB(6,2);
    t336 = SPB(6,3);
    t338 = SPB(6,4);
    t357 = (-t128*t334-t134*t336-t271*t338)*t54-t61*(-t3*t336-t31*t338)*t12-t61
*(t3*t334-t33*t338)*t17-t61*(t31*t334+t33*t336)*t38;
    t359 = t83*t270;
    t364 = complex<T>(1,0)/t299;
    t371 = (t62*t24+t8*t338)*t270+t53*t38*t270;
    t372 = complex<T>(1,0)/t371;
    t377 = SPA(2,5);
    t379 = SPA(2,6);
    t385 = pow(t299,2);
    t413 = SSS(1,2,3,4);
    t414 = complex<T>(1,0)/t413;
    t420 = pow(t3,2);
    t421 = pow(t17,2);
    t434 = t271*t54;
    t435 = t31*t12;
    t436 = t33*t17;
    t437 = t434+t435+t436;
    t439 = t33*t38;
    t440 = t135+t136-t439;
    t446 = t147*t437;
    t447 = t21*t440;
    t474 = SS(2,3,4);
    t475 = pow(t474,2);
    t486 = pow(-t377*t9-t97*t15-t99*t36,2);
    t494 = pow(-t31*t9-t33*t15,2);

    return(recursive*(-t1*(t2*(-t5*t6*t14*t19*t23*t25+t19*(-t3*t14*t19+(-t31*
t14-t34)*t40)*t44*t45+t8*t51/(t53*(t8*t54+t56*t57)-t61*t62*t67-t69*(t71+t56*t72
+t75)))*t85-t87*t4*t67*(t89*t19*t84+t40*(t34-t61*t92*t84)+t67*(t97*t19-t61*t99*
t40-t61*t56*t84))*t44*t25*t45*t85)*t116/t11*t120*t125-t1*t131*(t133*t137+t87*
t134*t21*t131)*t116*t144*t23*t151*t25*t153*t156-t1*t160*(t162*t4*t44*t25+(t166*
(t128*(t97*t92*t4-t61*t97*t56*t6)*t182-t56*t21*(t128*t137-t128*t33*t38)*t189)*
t45*t194+(t162*t4*(t198*t199/(-t200-t177-t181)+t210/(t14*t74-t212*t176-t214*
t180))+t222*t128*t4*(-t33*t56*t21*t147*(t89*t6+t92*t24)*t74-t233*t4*t14*t176+(-
t233*t14*t38-t61*t242*t74)*t189)*t45*t194)*t44*t25)/t182)*t116*t144*t260*t120*
t153-t265*(t222*t21*t4*t270*(-t3*t285*t299*t270+t309*t180*(-t333-t8*t357+t359))
*t364*t372+t166*(-t21*t4*(-t377*t57-t379*t74)*t180+t285*t270*(-t5*t385*((t62*t6
+t8*t336)*t270+t53*t17*t270)-t233*t285*t309*t371)*t364/(t333-t69*t357-t359)*
t372))*t116*t23*t120*t125*t153*t414)-t265*(-t420*t421*(t61*t144*t25*t156-t9*t25
*t156/(-t128*t9+t155))-(-t160*(t198*t199/(-t200-t214*t437-t212*t440)+t210/(-
t446-t447-t61*t311*(t62*t54+t377*t12+t97*t17+t99*t38)-t212*(t135+t136-t439-t174
)-t214*(t434+t435+t436-t178)))*t260/(-t446-t447)-(-t315*t437*t199-t56*t150*t209
*t74)*t475*t151/t315*t364*t414+(t486*t209/(-t75+t413)+t494*t199/(t155-t413))*
t260*t414)*t25)/t166*t44/t4*t153);
  }
}




#else

return
(

recursive*(-((complex<T>(0,1)*((complex<T>(-2,0)*((-SPA(2,3)*SPB(4,3)*SPB(5,3)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))/(pow(SPB(3,2),2)*SPB(5,4))-
(-(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3))*(-(SPA(2,3)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))+
(-(SPA(2,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2)))-
SPA(3,4)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))))/(SPB(3,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2)))-
(-SPA(1,6)*pow(-(SPA(2,3)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))-
SPA(2,4)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4)),2))/(-(-SPA(1,7)*(SPA(1,6)*SPB(7,1)+
SPA(5,6)*SPB(7,5)))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,6)*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(6,7)*SPB(7,6)))))/(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))-
(complex<T>(-3,0)*SPB(4,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(-(-SPA(3,6)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3))*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6)))+
(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))*(-(-SPA(3,4)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))-
complex<T>(-1,0)*SPA(4,6)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6)))+
(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(-(-SPA(3,5)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))-
complex<T>(-1,0)*SPA(4,5)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))-
complex<T>(-1,0)*SPA(5,6)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6)))))/(SPB(3,2)*SPB(5,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6)))))/(complex<T>(6,0)*SPA(6,7)*pow(SPB(4,3),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))))-
(complex<T>(0,1)*(SPA(1,2)*SPB(7,1)-SPA(2,3)*SPB(7,3))*(SPA(1,2)*SPB(3,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
complex<T>(-3,0)*SPA(1,3)*SPB(3,2)*(SPA(1,2)*SPB(7,1)-SPA(2,3)*SPB(7,3))))/(complex<T>(6,0)*SPA(1,2)*pow(SPB(3,2),2)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(5,4)*SPB(7,6)*SS(1,2,3))-
(complex<T>(0,1)*pow(SPB(7,1),2)*((complex<T>(-3,0)*SPA(1,2)*SPB(4,3))/(SPB(3,2)*SPB(5,4))+
((complex<T>(2,0)*(-(-SPA(1,2)*(-(-SPA(3,5)*SPA(4,6)*SPB(4,3))-
complex<T>(-1,0)*SPA(3,5)*SPA(5,6)*SPB(5,3))*(-(SPB(3,2)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6)))-
SPB(4,2)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))))-SPA(5,6)*SPB(3,2)*(-(-SPA(1,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)))-
SPA(1,2)*SPA(3,4)*SPB(7,4))*SS(3,4,5)))/((SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SS(3,4,5))+
(complex<T>(-3,0)*SPA(1,2)*SPB(4,3)*((pow(SPA(3,4),2)*pow(SPB(3,2),2)*pow(SPB(7,4),2))/(-(SPB(3,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)))-
SPB(3,2)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))-
SPB(4,2)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)))+
(pow(-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2),2)*pow(SPB(7,5),2))/(-(-(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,6))-
complex<T>(-1,0)*SPB(3,2)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))-
complex<T>(-1,0)*SPB(4,2)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))))+
(complex<T>(3,0)*SPA(1,2)*SPB(4,3)*(-(SPA(3,4)*SPA(5,6)*SPB(3,2)*SPB(4,2)*(SPA(3,6)*SPB(5,3)+
SPA(4,6)*SPB(5,4))*SPB(7,6))-
SPA(3,4)*SPB(3,2)*SPB(4,3)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))+
(-(SPA(3,4)*SPB(3,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,4))-
complex<T>(-1,0)*pow(-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2),2)*SPB(7,6))*SS(3,4,5)))/((SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SS(3,4,5)))/(SPB(3,2)*SPB(5,4)))/(-(SPB(3,2)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6)))-
SPB(4,2)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)))))/(complex<T>(6,0)*SPA(1,2)*SPB(2,1)*pow(SPB(4,3),2)*SPB(7,6))-
(complex<T>(0,-1)*((complex<T>(3,0)*SPB(3,2)*SPB(4,3)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))*(-(SPA(2,3)*(-(-(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*SPB(7,1))-
complex<T>(-1,0)*SPA(2,4)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))*SPB(7,3)-
complex<T>(-1,0)*SPA(2,4)*SPB(3,2)*SPB(7,4))*(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))+
(-(-(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(7,1))-
SPA(2,3)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*SPA(2,3)*SPB(4,2)*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(7,4))*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))*(-(SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
SPA(1,7)*SPB(7,1)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))))/((-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4))*((SPA(1,5)*SPB(5,4)+
SPA(1,6)*SPB(6,4))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))+
SPA(1,7)*SPB(7,4)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))))+
complex<T>(2,0)*(-SPB(3,2)*SPB(4,3)*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))+
((-(-(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*SPB(7,1))-
complex<T>(-1,0)*SPA(2,4)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))*SPB(7,3)-
complex<T>(-1,0)*SPA(2,4)*SPB(3,2)*SPB(7,4))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))*(-SPA(2,3)*SPB(4,3)*pow(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2)*((SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))+
SPA(1,7)*SPB(7,3)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))-SPA(3,4)*SPB(3,2)*(-(-(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*SPB(7,1))-
complex<T>(-1,0)*SPA(2,4)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))*SPB(7,3)-
complex<T>(-1,0)*SPA(2,4)*SPB(3,2)*SPB(7,4))*(-(-(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(7,1))-
SPA(2,3)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*SPA(2,3)*SPB(4,2)*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(7,4))*((SPA(1,5)*SPB(5,4)+
SPA(1,6)*SPB(6,4))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))+
SPA(1,7)*SPB(7,4)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))))/((-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4))*(-(-SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))-SPA(1,7)*SPB(7,1)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))*((SPA(1,5)*SPB(5,4)+
SPA(1,6)*SPB(6,4))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))+
SPA(1,7)*SPB(7,4)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))))))/(complex<T>(6,0)*pow(SPB(3,2),2)*pow(SPB(4,3),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*SPB(7,6)*SSS(1,2,3,4)))-
(complex<T>(0,-1)*(-pow(SPA(2,3),2)*pow(SPB(7,3),2)*(complex<T>(-1,0)/(SPA(1,2)*SPB(5,4)*SS(1,2,3))+
(-SPB(2,1))/(SPB(5,4)*SS(1,2,3)*(-SPA(1,2)*SPB(2,1)+
SS(1,2,3))))-
(-(-((-pow(SPB(7,1),2)*((pow(SPA(3,4),2)*pow(SPB(3,2),2)*pow(SPB(7,4),2))/(-(SPB(3,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)))-
complex<T>(-1,0)*SPB(4,2)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))-
complex<T>(-1,0)*SPB(3,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)-SPA(3,4)*SPB(7,4)))+
(pow(-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2),2)*pow(SPB(7,5),2))/(-(SPB(4,2)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3)))-
SPB(3,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)-SPA(3,4)*SPB(7,4))-
complex<T>(-1,0)*SPB(5,2)*(SPA(1,5)*SPB(7,1)+
SPA(2,5)*SPB(7,2)+
SPA(3,5)*SPB(7,3)+
SPA(4,5)*SPB(7,4))-
complex<T>(-1,0)*SPB(3,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)-SPA(3,4)*SPB(7,4)-SPA(3,5)*SPB(7,5))-
complex<T>(-1,0)*SPB(4,2)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3)-SPA(4,5)*SPB(7,5)))))/(SPB(2,1)*(-(SPB(4,2)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3)))-
SPB(3,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)-SPA(3,4)*SPB(7,4)))))+
((-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))*pow(SPB(7,4),2)-SPA(5,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*pow(SPB(7,5),2)*SPB(7,6))*pow(SS(2,3,4),2))/((-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4))*SSS(1,2,3,4))-
((pow(-SPA(2,5)*SPB(2,1)-SPA(3,5)*SPB(3,1)-SPA(4,5)*SPB(4,1),2)*pow(SPB(7,5),2))/(-SPA(6,7)*SPB(7,6)+
SSS(1,2,3,4))+
(pow(-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1),2)*pow(SPB(7,4),2))/(SS(1,2,3)-SSS(1,2,3,4)))/(SPB(2,1)*SSS(1,2,3,4))))/SPB(5,4)))/(complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(7,6))
)

	;

#endif
}


template <class T> complex<T> R2q3g2l_qppmmqbmlmlbp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
//{{qp, p, m, m, qbm, lm, lbp}, lc}

#if _VERBOSE
  _MESSAGE("R2q3g2l_eval :  qppmmqbmlmlbp lc");
#endif


complex<T> recursive;
recursive=complex<T>(1,0);
//new
#if _USE_OPTIMIZED

{
  complex<T> t1;
  complex<T> t100;
  complex<T> t101;
  complex<T> t113;
  complex<T> t115;
  complex<T> t117;
  complex<T> t12;
  complex<T> t121;
  complex<T> t125;
  complex<T> t126;
  complex<T> t127;
  complex<T> t128;
  complex<T> t130;
  complex<T> t131;
  complex<T> t135;
  complex<T> t139;
  complex<T> t141;
  complex<T> t142;
  complex<T> t145;
  complex<T> t148;
  complex<T> t151;
  complex<T> t152;
  complex<T> t153;
  complex<T> t155;
  complex<T> t165;
  complex<T> t171;
  complex<T> t172;
  complex<T> t177;
  complex<T> t179;
  complex<T> t18;
  complex<T> t180;
  complex<T> t184;
  complex<T> t19;
  complex<T> t195;
  complex<T> t196;
  complex<T> t198;
  complex<T> t199;
  complex<T> t2;
  complex<T> t20;
  complex<T> t200;
  complex<T> t204;
  complex<T> t206;
  complex<T> t209;
  complex<T> t21;
  complex<T> t211;
  complex<T> t213;
  complex<T> t214;
  complex<T> t217;
  complex<T> t218;
  complex<T> t219;
  complex<T> t22;
  complex<T> t220;
  complex<T> t229;
  complex<T> t233;
  complex<T> t24;
  complex<T> t241;
  complex<T> t244;
  complex<T> t25;
  complex<T> t256;
  complex<T> t257;
  complex<T> t261;
  complex<T> t263;
  complex<T> t266;
  complex<T> t27;
  complex<T> t274;
  complex<T> t280;
  complex<T> t286;
  complex<T> t289;
  complex<T> t29;
  complex<T> t3;
  complex<T> t30;
  complex<T> t303;
  complex<T> t307;
  complex<T> t314;
  complex<T> t32;
  complex<T> t322;
  complex<T> t323;
  complex<T> t325;
  complex<T> t341;
  complex<T> t342;
  complex<T> t343;
  complex<T> t344;
  complex<T> t345;
  complex<T> t346;
  complex<T> t350;
  complex<T> t356;
  complex<T> t36;
  complex<T> t361;
  complex<T> t363;
  complex<T> t37;
  complex<T> t370;
  complex<T> t371;
  complex<T> t373;
  complex<T> t378;
  complex<T> t379;
  complex<T> t387;
  complex<T> t388;
  complex<T> t39;
  complex<T> t390;
  complex<T> t397;
  complex<T> t399;
  complex<T> t4;
  complex<T> t40;
  complex<T> t417;
  complex<T> t419;
  complex<T> t421;
  complex<T> t422;
  complex<T> t424;
  complex<T> t431;
  complex<T> t433;
  complex<T> t434;
  complex<T> t444;
  complex<T> t451;
  complex<T> t452;
  complex<T> t454;
  complex<T> t459;
  complex<T> t46;
  complex<T> t460;
  complex<T> t470;
  complex<T> t471;
  complex<T> t475;
  complex<T> t477;
  complex<T> t478;
  complex<T> t480;
  complex<T> t483;
  complex<T> t485;
  complex<T> t486;
  complex<T> t49;
  complex<T> t495;
  complex<T> t498;
  complex<T> t499;
  complex<T> t5;
  complex<T> t50;
  complex<T> t505;
  complex<T> t506;
  complex<T> t507;
  complex<T> t508;
  complex<T> t51;
  complex<T> t518;
  complex<T> t52;
  complex<T> t521;
  complex<T> t522;
  complex<T> t528;
  complex<T> t531;
  complex<T> t536;
  complex<T> t538;
  complex<T> t553;
  complex<T> t554;
  complex<T> t558;
  complex<T> t560;
  complex<T> t567;
  complex<T> t57;
  complex<T> t570;
  complex<T> t571;
  complex<T> t575;
  complex<T> t576;
  complex<T> t58;
  complex<T> t581;
  complex<T> t584;
  complex<T> t587;
  complex<T> t595;
  complex<T> t597;
  complex<T> t598;
  complex<T> t6;
  complex<T> t600;
  complex<T> t604;
  complex<T> t61;
  complex<T> t611;
  complex<T> t616;
  complex<T> t617;
  complex<T> t62;
  complex<T> t625;
  complex<T> t626;
  complex<T> t64;
  complex<T> t643;
  complex<T> t647;
  complex<T> t65;
  complex<T> t657;
  complex<T> t659;
  complex<T> t662;
  complex<T> t665;
  complex<T> t68;
  complex<T> t681;
  complex<T> t682;
  complex<T> t69;
  complex<T> t696;
  complex<T> t7;
  complex<T> t70;
  complex<T> t700;
  complex<T> t701;
  complex<T> t72;
  complex<T> t73;
  complex<T> t732;
  complex<T> t734;
  complex<T> t752;
  complex<T> t755;
  complex<T> t76;
  complex<T> t765;
  complex<T> t78;
  complex<T> t8;
  complex<T> t80;
  complex<T> t82;
  complex<T> t83;
  complex<T> t86;
  complex<T> t87;
  complex<T> t9;
  complex<T> t90;
  complex<T> t93;
  complex<T> t94;
  complex<T> t96;
  {
    t1 = complex<T>(0,-1);
    t2 = SPA(1,3);
    t3 = SPB(3,2);
    t4 = t2*t3;
    t5 = SPB(7,2);
    t6 = pow(t5,2);
    t7 = SPA(1,2);
    t8 = pow(t7,2);
    t9 = SPB(4,2);
    t12 = SPB(4,3);
    t18 = t7*t5;
    t19 = SPB(7,3);
    t20 = t2*t19;
    t21 = -t18-t20;
    t22 = pow(t21,2);
    t24 = pow(t3,2);
    t25 = SPB(7,1);
    t27 = SPA(2,3);
    t29 = t2*t25+t27*t5;
    t30 = pow(t29,2);
    t32 = complex<T>(1,0)/t27;
    t36 = SPB(2,1);
    t37 = pow(t36,2);
    t39 = t27*t3;
    t40 = SS(1,2,3);
    t46 = complex<T>(1,0)/t7;
    t49 = -t7*t9-t2*t12;
    t50 = complex<T>(1,0)/t49;
    t51 = t46*t50;
    t52 = complex<T>(1,0)/t40;
    t57 = complex<T>(2,0);
    t58 = complex<T>(1,0)/t57;
    t61 = SPB(5,4);
    t62 = complex<T>(1,0)/t61;
    t64 = SPB(7,6);
    t65 = complex<T>(1,0)/t64;
    t68 = t1*t2;
    t69 = SPA(4,5);
    t70 = pow(t69,2);
    t72 = SPB(7,5);
    t73 = pow(t72,2);
    t76 = SPA(1,4);
    t78 = SPA(2,4);
    t80 = SPA(3,4);
    t82 = t76*t25+t78*t5+t80*t19;
    t83 = complex<T>(1,0)/t82;
    t86 = SPA(1,5);
    t87 = SPA(4,6);
    t90 = SPA(5,6);
    t93 = SPA(1,6);
    t94 = SPB(6,4);
    t96 = SPB(6,5);
    t100 = SPA(1,7);
    t101 = SPB(7,4);
    t113 = pow(t80,2);
    t115 = SPA(2,5);
    t117 = SPA(3,5);
    t121 = pow(t86*t25+t115*t5+t117*t19+t69*t101,2);
    t125 = SPB(5,2);
    t126 = t7*t125;
    t127 = SPB(5,3);
    t128 = t2*t127;
    t130 = -t126-t128-t76*t61;
    t131 = complex<T>(1,0)/t130;
    t135 = pow(t101,2);
    t139 = SPB(6,1);
    t141 = SPB(6,2);
    t142 = t27*t141;
    t145 = SPA(5,7);
    t148 = pow(t90*(t2*t139+t142)+t145*t29,2);
    t151 = SPA(6,7);
    t152 = t151*t64;
    t153 = SSS(1,2,3,4);
    t155 = complex<T>(1,0)/(t152-t153);
    t165 = complex<T>(1,0)/t153;
    t171 = t50*t62;
    t172 = pow(t64,2);
    t177 = -t49;
    t179 = t76*t9;
    t180 = t4+t179;
    t184 = complex<T>(1,0)/t12;
    t195 = t93*t36+t151*t5;
    t196 = pow(t195,2);
    t198 = complex<T>(1,0)/t151;
    t199 = t198*t62;
    t200 = SS(1,6,7);
    t204 = pow(t9,2);
    t206 = t78*t9;
    t209 = t80*t9;
    t211 = t80*t3;
    t213 = t180*t25+(t39+t206)*t5+t209*t19-t211*t101;
    t214 = pow(t213,2);
    t217 = t117*t3;
    t218 = t69*t9;
    t219 = -t217-t218;
    t220 = pow(t219,2);
    t229 = SSS(1,5,6,7);
    t233 = pow(t93,2);
    t241 = pow(-t18-t20-t76*t101,2);
    t244 = t65*t165;
    t256 = complex<T>(1,0)/t3;
    t257 = t58*t256;
    t261 = pow(t90,2);
    t263 = pow(-t126-t128,2);
    t266 = t87*t9;
    t274 = t7*(-t266-t90*t125)+t2*(-t87*t12-t90*t127);
    t280 = t274*t64+t49*(-t69*t72-t87*t64);
    t286 = complex<T>(1,0)/t274;
    t289 = t7*t3;
    t303 = complex<T>(0,1);
    t307 = -t86*t72-t93*t64;
    t314 = t27*t125;
    t322 = t130*t25+(-t27*t127-t78*t61)*t5+(t314-t80*t61)*t19+(t78*t125+t80*
t127)*t101;
    t323 = t86*t322;
    t325 = SPB(6,3);
    t341 = (-t7*t141-t2*t325-t76*t94)*t25+(-t27*t325-t78*t94)*t5+(t142-t80*t94)
*t19+(t78*t141+t80*t325)*t101;
    t342 = t93*t341;
    t343 = t100*t25;
    t344 = t343*t307;
    t345 = -t323-t342+t344;
    t346 = complex<T>(1,0)/t345;
    t350 = pow(t307,2);
    t356 = complex<T>(1,0)/t307;
    t361 = SPA(3,6);
    t363 = t117*t72+t361*t64;
    t370 = t80*t12;
    t371 = t206+t370;
    t373 = t49*t25-t27*t12*t5+t27*t9*t19+t371*t101;
    t378 = -t323-t342+(-t39+t343)*t307;
    t379 = complex<T>(1,0)/t378;
    t387 = complex<T>(-2,0);
    t388 = t387*t80;
    t390 = complex<T>(-3,0);
    t397 = t2*t80;
    t399 = -t289+t76*t12;
    t417 = pow(t80,3);
    t419 = pow(t373,2);
    t421 = pow(t86,2);
    t422 = pow(t322,2);
    t424 = pow(t341,2);
    t431 = pow(t27,2);
    t433 = pow(t100,2);
    t434 = pow(t25,2);
    t444 = pow(t378,3);
    t451 = complex<T>(3,0);
    t452 = complex<T>(1,0)/t451;
    t454 = pow(t345,2);
    t459 = complex<T>(1,0);
    t460 = complex<T>(1,0)/t459;
    t470 = t93*t96+t100*t72;
    t471 = complex<T>(1,0)/t470;
    t475 = t57*t2;
    t477 = t3*t61;
    t478 = t93*t27;
    t480 = SPB(4,1);
    t483 = t93*t480+t151*t101;
    t485 = t478*t9+t2*t483;
    t486 = pow(t485,2);
    t495 = SPB(5,1);
    t498 = t93*t495+t151*t72;
    t499 = t2*t498;
    t505 = complex<T>(6,0);
    t506 = complex<T>(1,0)/t505;
    t507 = pow(t49,2);
    t508 = complex<T>(1,0)/t507;
    t518 = t93*t139;
    t521 = -t478*t141-t2*t100*t25-t100*t27*t5-t2*(t518+t152);
    t522 = pow(t2,2);
    t528 = complex<T>(1,0)/t233;
    t531 = complex<T>(-1,0)/t521;
    t536 = complex<T>(1,0)/t93;
    t538 = t314+t499*t536;
    t553 = -t361*t3-t266;
    t554 = pow(t553,2);
    t558 = t57*t80;
    t560 = pow(t9,3);
    t567 = t506*t528;
    t570 = complex<T>(1,0)/t219;
    t571 = pow(t12,2);
    t575 = SS(2,3,4);
    t576 = complex<T>(1,0)/t575;
    t581 = t390*t125;
    t584 = t581*t256*t62*t471;
    t587 = t93*t94+t100*t101;
    t595 = t93*t141+t100*t5;
    t597 = t80*t587;
    t598 = t57*t27;
    t600 = t598*t3+t206+t370;
    t604 = t518+t343+t152;
    t611 = complex<T>(1,0)/(t371*t595+t179*t604);
    t616 = pow(t470,2);
    t617 = complex<T>(1,0)/t616;
    t625 = complex<T>(1,0)/t604;
    t626 = t567*t625;
    t643 = pow(t93,2);
    t647 = pow(t371,2);
    t657 = pow(t80,2);
    t659 = pow(t587,2);
    t662 = pow(t117,2);
    t665 = complex<T>(5,0);
    t681 = -t211*t587+t218*t470;
    t682 = pow(t681,2);
    t696 = -t597-t117*t470;
    t700 = SS(3,4,5);
    t701 = complex<T>(1,0)/t700;
    t732 = -t93*t219*t61;
    t734 = -t553*t587;
    t752 = pow(t90*t595+t86*t195,2);
    t755 = pow(t498,2);
    t765 = pow(t93,4);
    return(t1*(t4*t6/(t8*t9+t7*t2*t12)+t2*t22*(t24*t30*t32/t22+t37*t3/(-t39+t40
))*t51*t52)*t58/t24*t62*t65+t68*(t70*t49*t30*t73*t64*t83/(t49*t82+(t86*t87*t61-
t76*t90*t61+t93*(t87*t94+t90*t96)+t100*(t87*t101+t90*t72))*t64)*t52+(t113*(t49*
t121*t73/t90*t131*t83+t135*t64)+t148*t73*t64*t52*t155-t113*t135*t64*t40/(t40-
t153))*t165)*t58*t46*t32*t171/t172-t1*((t177*t125+t180*t61)*t6*t184/t177*t62/
t130*t65-(t125*t196*t199/t200-t180*(t204*t214*t165+t204*t220*t73*t155)/t204*t50
*t65/t229+t37*(-t233*t125*t199/(-t152+t200)+t180*t241*t50*t244/(-t153+t229)))*
t184*t131)*t257+recursive*(-t68*t3*t30*t21*(-t261*t263*t64+t274*t280)*t58*t51*
t286*t62/(-t289*t29+t4*(t7*t25-t27*t19))*t65/t280*t52+t303*(-t180*t214*t307*
t257*t50*t346+t303*t350*(-t303*t180*t214*t58*t256*t50*t356*t346+t303*(t80*(-
t363*t46+t209*t373*t307*t50*t379)*t58*t32*t356+(-t388*t214-t390*t113*t3*t213*
t373*t307*t379-(-(t397*(t399*t25+t78*t12*t5+(t39+t370)*t19+t78*t3*t101)+t76*t27
*t213)*t363*t345*t46*t356-t417*t3*t419*t307*(t421*t422+t233*t424-t57*t93*t100*
t25*t341*t307+(-t431*t24+t433*t434)*t350-t57*t86*t322*(-t342+t344))/t444)*t32)*
t452/t454)*t460)*t460)*t460*t184*t471*t244+(-t303*(t475*t80*t477*t486-t451*t9*
t49*(t478*t209+t397*t483)*(t478*t125+t499))*t506*t508+t233*t12*t61*t521*(t1*
t522*t80*t274*t485*t452*t528*t508*t531+t303*t2*t538*(t90*t538*t286-t80*t485*
t536*t531)*t58*t171+t7*t27*t470*(-t303*(-t390*t180*t219*t554*t12*t49-t558*t261*
t3*t560*t399*t470)*t567*t256*t570/t571*t508*t471*t576-(-t1*t196*(-t584+t113*
t587*(t390*t9*t62*t531-t387*t3*t595*(t451+t597*(t600*t595+(t475*t3+t179)*t604)*
t531*t611)*t570*t617)*t611)*t626-t303*t113*t483*(-t558*t3*t600*t61*t483+t451*
t371*(-t598*t477*t195-t9*t498*t575))*t506/t643*t32/t647*t62*t471*t625*t576-t303
*t196*(-t584+(t57*(t80*t595*(-t657*t117*t24*t659-t80*(t451*t662*t24+t665*t117*
t3*t219+t57*t220)*t587*t470-t57*t117*t70*t204*t616)/t682-t117*(-t80*(t209+t117*
t125)*(t93*t325+t100*t19)+(-t211+t69*t125)*t696)*t701)*t570*t471+t390*t80*(-
t117*t69*t9*t61*t470-t597*(-t217*t61+t9*t700))*t62/t681*t701)/t696)*t626+t1*(
t233*(t388*t3*t180*(t732+t734)*(t732+t451*t49*t195+t734)*t570*t508+t451*t554*
t125*t470*t62)+t581*t752*t470*t755*t62*t625/(t200-t229))*t506/t765*t256*t617*
t576)*t184))*t46)*t32*t198*t184*t62*t471*t531));
  }
}

#else

return
(

(complex<T>(0,-1)*((SPA(1,3)*SPB(3,2)*pow(SPB(7,2),2))/(-(-pow(SPA(1,2),2)*SPB(4,2))+SPA(1,2)*SPA(1,3)*SPB(4,3))+
(SPA(1,3)*pow(-SPA(1,2)*SPB(7,2)-SPA(1,3)*SPB(7,3),2)*(-((-pow(SPB(3,2),2)*pow(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2),2))/(SPA(2,3)*pow(-SPA(1,2)*SPB(7,2)-SPA(1,3)*SPB(7,3),2)))-
(-pow(SPB(2,1),2)*SPB(3,2))/(-SPA(2,3)*SPB(3,2)+
SS(1,2,3))))/(SPA(1,2)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SS(1,2,3))))/(complex<T>(2,0)*pow(SPB(3,2),2)*SPB(5,4)*SPB(7,6))+
(complex<T>(0,-1)*SPA(1,3)*(-((-pow(SPA(4,5),2)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*pow(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2),2)*pow(SPB(7,5),2)*SPB(7,6))/((SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))*((-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))+(-(-SPA(1,5)*SPA(4,6)*SPB(5,4))-
SPA(1,4)*SPA(5,6)*SPB(5,4)+SPA(1,6)*(SPA(4,6)*SPB(6,4)+
SPA(5,6)*SPB(6,5))+SPA(1,7)*(SPA(4,6)*SPB(7,4)+
SPA(5,6)*SPB(7,5)))*SPB(7,6))*SS(1,2,3)))+
(pow(SPA(3,4),2)*(-((-(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*pow(SPA(1,5)*SPB(7,1)+
SPA(2,5)*SPB(7,2)+
SPA(3,5)*SPB(7,3)+
SPA(4,5)*SPB(7,4),2)*pow(SPB(7,5),2))/(SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))))+pow(SPB(7,4),2)*SPB(7,6))-
(-pow(-(-SPA(5,6)*(SPA(1,3)*SPB(6,1)+
SPA(2,3)*SPB(6,2)))+SPA(5,7)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)),2)*pow(SPB(7,5),2)*SPB(7,6))/(SS(1,2,3)*(SPA(6,7)*SPB(7,6)-SSS(1,2,3,4)))-
(pow(SPA(3,4),2)*pow(SPB(7,4),2)*SPB(7,6)*SS(1,2,3))/(SS(1,2,3)-SSS(1,2,3,4)))/SSS(1,2,3,4)))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(5,4)*pow(SPB(7,6),2))-
(complex<T>(0,-1)*(-(((-(-(SPA(1,2)*SPB(4,2)+
SPA(1,3)*SPB(4,3))*SPB(5,2))+(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(5,4))*pow(SPB(7,2),2))/(SPB(4,3)*(SPA(1,2)*SPB(4,2)+
SPA(1,3)*SPB(4,3))*SPB(5,4)*(SPA(1,2)*SPB(5,2)+
SPA(1,3)*SPB(5,3)+
SPA(1,4)*SPB(5,4))*SPB(7,6)))-
(-((-SPB(5,2)*pow(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2),2))/(SPA(6,7)*SPB(5,4)*SS(1,6,7)))-
((SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*((pow(SPB(4,2),2)*pow(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))+(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)+SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2))/SSS(1,2,3,4)+
(pow(SPB(4,2),2)*pow(-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2),2)*pow(SPB(7,5),2))/(SPA(6,7)*SPB(7,6)-SSS(1,2,3,4))))/(pow(SPB(4,2),2)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(7,6)*SSS(1,5,6,7))+
pow(SPB(2,1),2)*(-((pow(SPA(1,6),2)*SPB(5,2))/(SPA(6,7)*SPB(5,4)*(-SPA(6,7)*SPB(7,6)+
SS(1,6,7))))-
(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*pow(-SPA(1,2)*SPB(7,2)-SPA(1,3)*SPB(7,3)-SPA(1,4)*SPB(7,4),2))/((-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(7,6)*SSS(1,2,3,4)*(-SSS(1,2,3,4)+
SSS(1,5,6,7)))))/(SPB(4,3)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))))/(complex<T>(2,0)*SPB(3,2))+
recursive*(-((complex<T>(0,-1)*SPA(1,3)*SPB(3,2)*pow(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2),2)*(-SPA(1,2)*SPB(7,2)-SPA(1,3)*SPB(7,3))*(-(pow(SPA(5,6),2)*pow(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3),2)*SPB(7,6))+
(-(-SPA(1,2)*(-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2)))+SPA(1,3)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3)))*(-(-(-(-SPA(1,2)*(-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2)))+SPA(1,3)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3)))*SPB(7,6))+
(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)))))/(complex<T>(2,0)*SPA(1,2)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(-(-SPA(1,2)*(-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2)))+SPA(1,3)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3)))*SPB(5,4)*(-SPA(1,2)*SPB(3,2)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPA(1,3)*SPB(3,2)*(SPA(1,2)*SPB(7,1)-SPA(2,3)*SPB(7,3)))*SPB(7,6)*(-(-(-(-SPA(1,2)*(-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2)))+SPA(1,3)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3)))*SPB(7,6))+
(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)))*SS(1,2,3)))+
(complex<T>(0,1)*(-(((SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*pow(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))+(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)+SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))/(complex<T>(2,0)*SPB(3,2)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(-(SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
SPA(1,7)*SPB(7,1)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))))+
(complex<T>(0,1)*pow(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6),2)*(-((complex<T>(0,1)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*pow(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))+(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)+SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2))/(complex<T>(2,0)*SPB(3,2)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))*(-(SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
SPA(1,7)*SPB(7,1)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))))+
(complex<T>(0,1)*((SPA(3,4)*(-((-(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6)))/SPA(1,2))+
(SPA(3,4)*SPB(4,2)*(-(-(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(7,1))-
SPA(2,3)*SPB(4,3)*SPB(7,2)+SPA(2,3)*SPB(4,2)*SPB(7,3)+(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(7,4))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))/((-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(-(SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(-SPA(2,3)*SPB(3,2)+
SPA(1,7)*SPB(7,1))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))))))/(complex<T>(2,0)*SPA(2,3)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))+
(-(complex<T>(-2,0)*SPA(3,4)*pow(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))+(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)+SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2))-
(complex<T>(-3,0)*pow(SPA(3,4),2)*SPB(3,2)*(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))+(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)+SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4))*(-(-(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(7,1))-
SPA(2,3)*SPB(4,3)*SPB(7,2)+SPA(2,3)*SPB(4,2)*SPB(7,3)+(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(7,4))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))/(-(SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(-SPA(2,3)*SPB(3,2)+
SPA(1,7)*SPB(7,1))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))-
(-((-(SPA(1,3)*SPA(3,4)*(-(-(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*SPB(7,1))+SPA(2,4)*SPB(4,3)*SPB(7,2)+(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))*SPB(7,3)+SPA(2,4)*SPB(3,2)*SPB(7,4))+
SPA(1,4)*SPA(2,3)*(-(-(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))+(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)+SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4)))*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))*(-(SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
SPA(1,7)*SPB(7,1)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))))/(SPA(1,2)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))))+
(-pow(SPA(3,4),3)*SPB(3,2)*pow(-(-(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(7,1))-
SPA(2,3)*SPB(4,3)*SPB(7,2)+SPA(2,3)*SPB(4,2)*SPB(7,3)+(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(7,4),2)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))*(pow(SPA(1,5),2)*pow(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4),2)+
pow(SPA(1,6),2)*pow(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4),2)-
complex<T>(2,0)*SPA(1,6)*SPA(1,7)*SPB(7,1)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))+
(-pow(SPA(2,3),2)*pow(SPB(3,2),2)+
pow(SPA(1,7),2)*pow(SPB(7,1),2))*pow(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6),2)-
complex<T>(2,0)*SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4))*(-(SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4)))+
SPA(1,7)*SPB(7,1)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))))/pow(-(SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(-SPA(2,3)*SPB(3,2)+
SPA(1,7)*SPB(7,1))*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)),3))/SPA(2,3))/(complex<T>(3,0)*pow(-(SPA(1,5)*(-(-(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1))+(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,2)+(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*SPB(7,3)+(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(-(-SPA(1,2)*SPB(6,2)-SPA(1,3)*SPB(6,3)-SPA(1,4)*SPB(6,4))*SPB(7,1))+(-SPA(2,3)*SPB(6,3)-SPA(2,4)*SPB(6,4))*SPB(7,2)+(SPA(2,3)*SPB(6,2)-SPA(3,4)*SPB(6,4))*SPB(7,3)+(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
SPA(1,7)*SPB(7,1)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)),2))))/complex<T>(1,0)))/complex<T>(1,0)))/(complex<T>(1,0)*SPB(4,3)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*SPB(7,6)*SSS(1,2,3,4))+
(-((complex<T>(0,1)*(complex<T>(2,0)*SPA(1,3)*SPA(3,4)*SPB(3,2)*SPB(5,4)*pow(-(-SPA(1,6)*SPA(2,3)*SPB(4,2))+SPA(1,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4)),2)-
complex<T>(3,0)*SPB(4,2)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(SPA(1,6)*SPA(2,3)*SPA(3,4)*SPB(4,2)+
SPA(1,3)*SPA(3,4)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4)))*(-(-SPA(1,6)*SPA(2,3)*SPB(5,2))+SPA(1,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))))/(complex<T>(6,0)*pow(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3),2))-
(pow(SPA(1,6),2)*SPB(4,3)*SPB(5,4)*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))*(-((complex<T>(0,-1)*pow(SPA(1,3),complex<T>(2,0))*SPA(3,4)*(-(-SPA(1,2)*(-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2)))+SPA(1,3)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3)))*(-(SPA(1,6)*SPA(2,3)*SPB(4,2))-
SPA(1,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))))/(complex<T>(3,0)*pow(SPA(1,6),2)*pow(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3),2)*(-(-SPA(1,6)*SPA(2,3)*SPB(6,2))+SPA(1,3)*SPA(1,7)*SPB(7,1)+SPA(1,7)*SPA(2,3)*SPB(7,2)+SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))))+
(complex<T>(0,1)*SPA(1,3)*(SPA(2,3)*SPB(5,2)+
(SPA(1,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/SPA(1,6))*(-((-SPA(5,6)*(SPA(2,3)*SPB(5,2)+
(SPA(1,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/SPA(1,6)))/(-(-SPA(1,2)*(-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2)))+SPA(1,3)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))))+
(-SPA(3,4)*(-(-SPA(1,6)*SPA(2,3)*SPB(4,2))+SPA(1,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))))/(SPA(1,6)*(-(-SPA(1,6)*SPA(2,3)*SPB(6,2))+SPA(1,3)*SPA(1,7)*SPB(7,1)+SPA(1,7)*SPA(2,3)*SPB(7,2)+SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))))/(complex<T>(2,0)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*SPB(5,4))+
SPA(1,2)*SPA(2,3)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(-((complex<T>(0,1)*(-(complex<T>(-3,0)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*(-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2))*pow(-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2),2)*SPB(4,3)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3)))-
complex<T>(2,0)*SPA(3,4)*pow(SPA(5,6),2)*SPB(3,2)*pow(SPB(4,2),3)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))))/(complex<T>(6,0)*pow(SPA(1,6),2)*SPB(3,2)*(-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2))*pow(SPB(4,3),2)*pow(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*SS(2,3,4)))-
(-((complex<T>(0,-1)*pow(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2),2)*(-((complex<T>(-3,0)*SPB(5,2))/(SPB(3,2)*SPB(5,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))))+
(pow(SPA(3,4),2)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*((complex<T>(-3,0)*SPB(4,2))/(SPB(5,4)*(-(-SPA(1,6)*SPA(2,3)*SPB(6,2))+SPA(1,3)*SPA(1,7)*SPB(7,1)+SPA(1,7)*SPA(2,3)*SPB(7,2)+SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))-
(complex<T>(-2,0)*SPB(3,2)*(SPA(1,6)*SPB(6,2)+
SPA(1,7)*SPB(7,2))*(complex<T>(3,0)-
(-SPA(3,4)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*((complex<T>(2,0)*SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*(SPA(1,6)*SPB(6,2)+
SPA(1,7)*SPB(7,2))+
(complex<T>(2,0)*SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))))/((-(-SPA(1,6)*SPA(2,3)*SPB(6,2))+SPA(1,3)*SPA(1,7)*SPB(7,1)+SPA(1,7)*SPA(2,3)*SPB(7,2)+SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))*((SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*(SPA(1,6)*SPB(6,2)+
SPA(1,7)*SPB(7,2))+
SPA(1,4)*SPB(4,2)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))))))/((-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2))*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2))))/((SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*(SPA(1,6)*SPB(6,2)+
SPA(1,7)*SPB(7,2))+
SPA(1,4)*SPB(4,2)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6)))))/(complex<T>(6,0)*pow(SPA(1,6),2)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6)))+
(complex<T>(0,1)*pow(SPA(3,4),2)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))*(-(complex<T>(2,0)*SPA(3,4)*SPB(3,2)*(complex<T>(2,0)*SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(5,4)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4)))+
complex<T>(3,0)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*(-(complex<T>(2,0)*SPA(2,3)*SPB(3,2)*SPB(5,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2)))-
SPB(4,2)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*SS(2,3,4))))/(complex<T>(6,0)*pow(SPA(1,6),complex<T>(2,0))*SPA(2,3)*pow(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3),2)*SPB(5,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))*SS(2,3,4))+
(complex<T>(0,1)*pow(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2),2)*(-((complex<T>(-3,0)*SPB(5,2))/(SPB(3,2)*SPB(5,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))))+
((complex<T>(2,0)*(-((-SPA(3,4)*(SPA(1,6)*SPB(6,2)+
SPA(1,7)*SPB(7,2))*(-(pow(SPA(3,4),complex<T>(2,0))*SPA(3,5)*pow(SPB(3,2),2)*pow(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4),2))-
SPA(3,4)*(complex<T>(3,0)*pow(SPA(3,5),2)*pow(SPB(3,2),2)+
complex<T>(5,0)*SPA(3,5)*SPB(3,2)*(-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2))+
complex<T>(2,0)*pow(-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2),2))*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-
complex<T>(2,0)*SPA(3,5)*pow(SPA(4,5),2)*pow(SPB(4,2),2)*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)))/pow(-SPA(3,4)*SPB(3,2)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))+
SPA(4,5)*SPB(4,2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)),2))-
(SPA(3,5)*(-(SPA(3,4)*(SPA(3,4)*SPB(4,2)+
SPA(3,5)*SPB(5,2))*(SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3)))+
(-SPA(3,4)*SPB(3,2)+
SPA(4,5)*SPB(5,2))*(-(SPA(3,4)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4)))-
SPA(3,5)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))))/SS(3,4,5)))/((-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))+
(complex<T>(-3,0)*SPA(3,4)*(-SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-
SPA(3,4)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*(-(SPA(3,5)*SPB(3,2)*SPB(5,4))+SPB(4,2)*SS(3,4,5))))/(SPB(5,4)*(-SPA(3,4)*SPB(3,2)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))+
SPA(4,5)*SPB(4,2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))*SS(3,4,5)))/(-(SPA(3,4)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4)))-
SPA(3,5)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))))/(complex<T>(6,0)*pow(SPA(1,6),2)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6)))-
(complex<T>(0,-1)*(pow(SPA(1,6),2)*((complex<T>(-2,0)*SPA(3,4)*SPB(3,2)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*(SPA(1,6)*(SPA(3,5)*SPB(3,2)+
SPA(4,5)*SPB(4,2))*SPB(5,4)+
(SPA(3,6)*SPB(3,2)+
SPA(4,6)*SPB(4,2))*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4)))*(SPA(1,6)*(SPA(3,5)*SPB(3,2)+
SPA(4,5)*SPB(4,2))*SPB(5,4)+
complex<T>(3,0)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
(SPA(3,6)*SPB(3,2)+
SPA(4,6)*SPB(4,2))*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))))/((-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2))*pow(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3),2))+
(complex<T>(3,0)*pow(-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2),2)*SPB(5,2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))/SPB(5,4))+
(complex<T>(-3,0)*SPB(5,2)*pow(-(-SPA(5,6)*(SPA(1,6)*SPB(6,2)+
SPA(1,7)*SPB(7,2)))+SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2)),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*pow(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5),2))/(SPB(5,4)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))*(SS(1,6,7)-SSS(1,5,6,7)))))/(complex<T>(6,0)*pow(SPA(1,6),4)*SPB(3,2)*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*SS(2,3,4))))/SPB(4,3))))/SPA(1,2)))/(SPA(2,3)*SPA(6,7)*SPB(4,3)*SPB(5,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(-(-SPA(1,6)*SPA(2,3)*SPB(6,2))+SPA(1,3)*SPA(1,7)*SPB(7,1)+SPA(1,7)*SPA(2,3)*SPB(7,2)+SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))))
)
;
#endif

}



template <class T> complex<T> R2q3g2l_qpppmqbmlmlbp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
//{{qp, p, p, m, qbm, lm, lbp}, nf}

#if _VERBOSE
  _MESSAGE("R2q3g2l_eval :  qpppmqbmlmlbp lc");
#endif

complex<T> recursive;
recursive=complex<T>(1,0);
//new

#if _USE_OPTIMIZED

{
  complex<T> t1;
  complex<T> t109;
  complex<T> t11;
  complex<T> t116;
  complex<T> t12;
  complex<T> t125;
  complex<T> t126;
  complex<T> t128;
  complex<T> t130;
  complex<T> t132;
  complex<T> t133;
  complex<T> t135;
  complex<T> t137;
  complex<T> t138;
  complex<T> t141;
  complex<T> t142;
  complex<T> t143;
  complex<T> t145;
  complex<T> t148;
  complex<T> t149;
  complex<T> t152;
  complex<T> t161;
  complex<T> t162;
  complex<T> t163;
  complex<T> t164;
  complex<T> t168;
  complex<T> t170;
  complex<T> t173;
  complex<T> t174;
  complex<T> t176;
  complex<T> t177;
  complex<T> t178;
  complex<T> t18;
  complex<T> t181;
  complex<T> t183;
  complex<T> t187;
  complex<T> t188;
  complex<T> t19;
  complex<T> t190;
  complex<T> t193;
  complex<T> t198;
  complex<T> t2;
  complex<T> t20;
  complex<T> t200;
  complex<T> t21;
  complex<T> t214;
  complex<T> t216;
  complex<T> t218;
  complex<T> t219;
  complex<T> t22;
  complex<T> t221;
  complex<T> t223;
  complex<T> t224;
  complex<T> t225;
  complex<T> t227;
  complex<T> t229;
  complex<T> t230;
  complex<T> t233;
  complex<T> t234;
  complex<T> t239;
  complex<T> t24;
  complex<T> t246;
  complex<T> t249;
  complex<T> t25;
  complex<T> t252;
  complex<T> t253;
  complex<T> t254;
  complex<T> t256;
  complex<T> t257;
  complex<T> t26;
  complex<T> t263;
  complex<T> t264;
  complex<T> t266;
  complex<T> t269;
  complex<T> t27;
  complex<T> t270;
  complex<T> t272;
  complex<T> t274;
  complex<T> t278;
  complex<T> t28;
  complex<T> t281;
  complex<T> t283;
  complex<T> t284;
  complex<T> t285;
  complex<T> t286;
  complex<T> t288;
  complex<T> t29;
  complex<T> t292;
  complex<T> t293;
  complex<T> t297;
  complex<T> t3;
  complex<T> t30;
  complex<T> t301;
  complex<T> t306;
  complex<T> t307;
  complex<T> t316;
  complex<T> t317;
  complex<T> t32;
  complex<T> t321;
  complex<T> t323;
  complex<T> t328;
  complex<T> t330;
  complex<T> t334;
  complex<T> t341;
  complex<T> t346;
  complex<T> t348;
  complex<T> t349;
  complex<T> t357;
  complex<T> t36;
  complex<T> t364;
  complex<T> t368;
  complex<T> t369;
  complex<T> t37;
  complex<T> t39;
  complex<T> t399;
  complex<T> t4;
  complex<T> t40;
  complex<T> t415;
  complex<T> t431;
  complex<T> t433;
  complex<T> t444;
  complex<T> t448;
  complex<T> t451;
  complex<T> t453;
  complex<T> t454;
  complex<T> t46;
  complex<T> t463;
  complex<T> t470;
  complex<T> t48;
  complex<T> t483;
  complex<T> t487;
  complex<T> t489;
  complex<T> t49;
  complex<T> t494;
  complex<T> t495;
  complex<T> t499;
  complex<T> t50;
  complex<T> t501;
  complex<T> t51;
  complex<T> t513;
  complex<T> t514;
  complex<T> t516;
  complex<T> t517;
  complex<T> t534;
  complex<T> t535;
  complex<T> t536;
  complex<T> t537;
  complex<T> t538;
  complex<T> t542;
  complex<T> t548;
  complex<T> t554;
  complex<T> t558;
  complex<T> t56;
  complex<T> t561;
  complex<T> t566;
  complex<T> t567;
  complex<T> t57;
  complex<T> t59;
  complex<T> t597;
  complex<T> t599;
  complex<T> t6;
  complex<T> t60;
  complex<T> t601;
  complex<T> t602;
  complex<T> t604;
  complex<T> t61;
  complex<T> t611;
  complex<T> t613;
  complex<T> t614;
  complex<T> t624;
  complex<T> t632;
  complex<T> t637;
  complex<T> t638;
  complex<T> t64;
  complex<T> t646;
  complex<T> t647;
  complex<T> t648;
  complex<T> t65;
  complex<T> t655;
  complex<T> t656;
  complex<T> t662;
  complex<T> t663;
  complex<T> t68;
  complex<T> t684;
  complex<T> t685;
  complex<T> t686;
  complex<T> t689;
  complex<T> t693;
  complex<T> t696;
  complex<T> t7;
  complex<T> t703;
  complex<T> t709;
  complex<T> t712;
  complex<T> t72;
  complex<T> t725;
  complex<T> t73;
  complex<T> t733;
  complex<T> t74;
  complex<T> t745;
  complex<T> t75;
  complex<T> t756;
  complex<T> t77;
  complex<T> t770;
  complex<T> t79;
  complex<T> t8;
  complex<T> t80;
  complex<T> t83;
  complex<T> t88;
  complex<T> t9;
  complex<T> t91;
  complex<T> t92;
  complex<T> t93;
  complex<T> t99;
  {
    t1 = complex<T>(0,-1);
    t2 = SPA(3,4);
    t3 = SPA(4,6);
    t4 = pow(t3,2);
    t6 = SPB(5,3);
    t7 = SPA(2,3);
    t8 = t7*t6;
    t9 = SPB(5,4);
    t11 = SPA(2,4);
    t12 = pow(t9,2);
    t18 = SPA(3,6);
    t19 = t18*t6;
    t20 = t3*t9;
    t21 = t19+t20;
    t22 = pow(t21,2);
    t24 = pow(t2,2);
    t25 = SPB(4,3);
    t26 = t3*t25;
    t27 = SPA(5,6);
    t28 = t27*t6;
    t29 = -t26-t28;
    t30 = pow(t29,2);
    t32 = complex<T>(1,0)/t25;
    t36 = SPA(4,5);
    t37 = pow(t36,2);
    t39 = t2*t25;
    t40 = SS(3,4,5);
    t46 = complex<T>(1,0)/t9;
    t48 = -t8-t11*t9;
    t49 = complex<T>(1,0)/t48;
    t50 = t46*t49;
    t51 = complex<T>(1,0)/t40;
    t56 = complex<T>(2,0);
    t57 = complex<T>(1,0)/t56;
    t59 = SPA(1,2);
    t60 = complex<T>(1,0)/t59;
    t61 = pow(t2,2);
    t64 = SPA(6,7);
    t65 = complex<T>(1,0)/t64;
    t68 = complex<T>(0,1);
    t72 = SPA(1,3);
    t73 = t72*t6;
    t74 = SPA(1,4);
    t75 = t74*t9;
    t77 = pow(-t73-t75,2);
    t79 = SPB(7,1);
    t80 = pow(t79,2);
    t83 = SPB(7,2);
    t88 = t11*t83;
    t91 = t6*(t72*t79+t7*t83)+t9*(t74*t79+t88);
    t92 = SPA(1,6);
    t93 = SPB(2,1);
    t99 = t48*(t92*t93+t64*t83)+t64*t91;
    t109 = t2*t6;
    t116 = complex<T>(1,0)/t91;
    t125 = SPB(3,2);
    t126 = t125*t6;
    t128 = SPB(7,5);
    t130 = SPA(2,5);
    t132 = SPA(2,6);
    t133 = SPB(7,6);
    t135 = -t130*t128-t132*t133;
    t137 = -t11*t25*t128+t6*t135;
    t138 = pow(t137,2);
    t141 = complex<T>(-3,0);
    t142 = t141*t11;
    t143 = t74*t25;
    t145 = SPA(1,5);
    t148 = -t145*t128-t92*t133;
    t149 = t6*t148;
    t152 = t11*t125;
    t161 = complex<T>(6,0);
    t162 = complex<T>(1,0)/t161;
    t163 = pow(t48,2);
    t164 = complex<T>(1,0)/t163;
    t168 = pow(t128,2);
    t170 = SPB(6,5);
    t173 = SPA(4,7);
    t174 = t173*t25;
    t176 = SPA(5,7);
    t177 = t176*t128;
    t178 = t64*t133;
    t181 = -t26*t170-t28*t170-t174*t128-t6*(t177+t178);
    t183 = pow(t6,2);
    t187 = complex<T>(3,0);
    t188 = complex<T>(1,0)/t187;
    t190 = complex<T>(1,0)/t168;
    t193 = complex<T>(-1,0)/t181;
    t198 = complex<T>(1,0)/t128;
    t200 = t143-t149*t198;
    t214 = SPA(1,7);
    t216 = t92*t170+t214*t128;
    t218 = t11*t93;
    t219 = SPB(3,1);
    t221 = -t218-t2*t219;
    t223 = SPB(5,2);
    t224 = t11*t223;
    t225 = t224+t109;
    t227 = SPB(7,3);
    t229 = t88+t2*t227;
    t230 = pow(t229,2);
    t233 = complex<T>(-2,0);
    t234 = pow(t11,3);
    t239 = t7*t223-t2*t9;
    t246 = pow(t7,2);
    t249 = complex<T>(1,0)/t2;
    t252 = complex<T>(1,0)/t221;
    t253 = t252*t164;
    t254 = complex<T>(1,0)/t216;
    t256 = SS(2,3,4);
    t257 = complex<T>(1,0)/t256;
    t263 = -t36*t128-t3*t133;
    t264 = pow(t263,2);
    t266 = t187*t74;
    t269 = t266*t60*t249*t254;
    t270 = pow(t125,2);
    t272 = SPA(2,7);
    t274 = t132*t170+t272*t128;
    t278 = t56*t2;
    t281 = t3*t170+t173*t128;
    t283 = t125*t274;
    t284 = t7*t125;
    t285 = SPB(4,2);
    t286 = t11*t285;
    t288 = t284+t286+t278*t25;
    t292 = t27*t170;
    t293 = t292+t177+t178;
    t297 = t284+t286;
    t301 = complex<T>(1,0)/(t297*t281+t224*t293);
    t306 = pow(t216,2);
    t307 = complex<T>(1,0)/t306;
    t316 = complex<T>(1,0)/t293;
    t317 = t162*t190*t316;
    t321 = pow(t11,2);
    t323 = pow(t93,2);
    t328 = pow(t219,2);
    t330 = complex<T>(5,0);
    t334 = pow(t221,2);
    t341 = pow(t274,2);
    t346 = t2*t125;
    t348 = t218*t216-t346*t274;
    t349 = pow(t348,2);
    t357 = SPA(3,7);
    t364 = -t219*t216-t283;
    t368 = SS(1,2,3);
    t369 = complex<T>(1,0)/t368;
    t399 = t233*t59;
    t415 = pow(t297,2);
    t431 = -t59*t221*t128;
    t433 = -t229*t274;
    t444 = pow(t148,2);
    t448 = SPB(5,1);
    t451 = pow(-t79*t281+t448*t263,2);
    t453 = SS(5,6,7);
    t454 = SSS(1,5,6,7);
    t463 = pow(t128,4);
    t470 = complex<T>(1,0)/t7;
    t483 = complex<T>(1,0)/t133;
    t487 = t132*t2;
    t489 = t11*t18;
    t494 = t487*t125-t489*t125-t3*(t286+t39)-t27*t225;
    t495 = pow(t494,2);
    t499 = t92*t448+t64*t128;
    t501 = t57*t249;
    t513 = t59*t223;
    t514 = -t513-t73-t75;
    t516 = -t132*(t72*t125+t74*t285)-t3*(-t59*t285-t72*t25)-t18*(-t59*t125+t143
)-t27*t514;
    t517 = t448*t516;
    t534 = -t132*(-t357*t125-t173*t285)-t3*(t272*t285+t357*t25)-t18*(t272*t125-
t174)-t27*(t272*t223+t357*t6+t173*t9);
    t535 = t534*t128;
    t536 = t292*t499;
    t537 = -t517+t535+t536;
    t538 = complex<T>(1,0)/t537;
    t542 = pow(t499,2);
    t548 = complex<T>(1,0)/t499;
    t554 = t92*t219+t64*t227;
    t558 = t7*t3;
    t561 = -t132*t297-t489*t25+t558*t25-t27*t48;
    t566 = -t517+t535+(-t39+t292)*t499;
    t567 = complex<T>(1,0)/t566;
    t597 = pow(t125,3);
    t599 = pow(t561,2);
    t601 = pow(t448,2);
    t602 = pow(t516,2);
    t604 = pow(t534,2);
    t611 = pow(t25,2);
    t613 = pow(t27,2);
    t614 = pow(t170,2);
    t624 = pow(t566,3);
    t632 = pow(t537,2);
    t637 = complex<T>(1,0);
    t638 = complex<T>(1,0)/t637;
    t646 = t470*t65;
    t647 = SSS(2,3,4,5);
    t648 = complex<T>(1,0)/t647;
    t655 = pow(t92,2);
    t656 = t655*t64;
    t662 = -t18*t125-t3*t285-t27*t223;
    t663 = complex<T>(1,0)/t662;
    t684 = pow(t132,2);
    t685 = t684*t64;
    t686 = pow(t92,2);
    t689 = SPB(4,1);
    t693 = pow(-t132*t93-t18*t219-t3*t689-t27*t448,2);
    t696 = complex<T>(1,0)/t514;
    t703 = SPB(6,1);
    t709 = pow(-t29*t703-(-t174-t176*t6)*t79,2);
    t712 = complex<T>(1,0)/(t178-t647);
    t725 = pow(t64,2);
    t733 = -t48;
    t745 = t60*t483;
    t756 = pow(t11,2);
    t770 = pow(t132*t223+t19+t20,2);
    return(t1*(t2*t4*t6/(t8*t9+t11*t12)+t6*t22*(t24*t30*t32/t22+t2*t37/(-t39+
t40))*t50*t51)*t57*t60/t61*t65+recursive*(t68*t2*t6*t30*t21*(-t64*t77*t80+t91*
t99)*t57*t60*t65*t50/(-t2*t29*t9+t109*(t18*t25-t27*t9))*t116/t99*t51+(-t68*(t56
*t59*t2*t126*t138+t142*t48*(-t143*t128+t149)*(-t152*t25*t128+t126*t135))*t162*
t164+t59*t7*t168*t181*(-t68*t125*t183*t91*t137*t188*t164*t190*t193+t68*t6*t200*
(-t79*t200*t116+t125*t137*t198*t193)*t57*t60*t49+t25*t9*t216*(t1*(t187*t7*t221*
t225*t48*t230+t233*t234*t2*t125*t239*t80*t216)*t162/t246*t249*t253*t190*t254*
t257+(t1*t264*(t269+t270*t274*(t142*t60*t193+t278*t281*(t187+t283*(t288*t281+(
t224+t278*t6)*t293)*t193*t301)*t252*t307)*t301)*t317+t68*t264*(t269+(t56*(t125*
t281*(t233*t321*t323*t219*t306-(t187*t24*t328+t330*t2*t219*t221+t56*t334)*t125*
t216*t274-t24*t219*t270*t341)/t349-t219*(-t125*(t74*t219+t152)*(t18*t170+t357*
t128)+(t74*t93-t346)*t364)*t369)*t252*t254+t141*t125*(-t59*t11*t93*t219*t216-
t283*(-t59*t2*t219+t11*t368))*t60/t348*t369)/t364)*t317+t68*t270*t135*(t399*t2*
t125*t288*t135+t187*t297*(t399*t39*t263-t11*t148*t256))*t162*t60/t415*t32*t190*
t254*t316*t257+t68*(t168*(t266*t230*t216*t60+t233*t2*t125*t225*(-t431+t433)*(-
t431+t433+t187*t48*t263)*t253)+t141*t74*t216*t444*t451*t60*t316/(t453-t454))*
t162*t249/t463*t307*t257)*t470))*t46)*t60*t470*t32*t254*t483*t193+t68*(-t225*
t495*t499*t501*t49*t538+t68*t542*(t1*t225*t495*t57*t249*t49*t548*t538+t68*(t125
*(t554*t46+t152*t561*t499*t49*t567)*t57*t32*t548+(t56*t125*t495+t187*t2*t270*
t494*t561*t499*t567+(-(t25*t223*t494+t126*(-t487*t285-t558*t285-t18*(t284+t39)-
t27*t239))*t554*t537*t46*t548+t2*t597*t599*t499*(t601*t602+t604*t168+t56*t27*
t534*t170*t128*t499+(-t24*t611+t613*t614)*t542+t233*t448*t516*(t535+t536))/t624
)*t32)*t188/t632)*t638)*t638)*t638*t646*t254*t648)+t1*t6*(t656*t323*t30*t48*
t663/(t662*t48+t64*(t513*t79-t59*t448*t83-t170*(t92*t79+t132*t83)-(t214*t79+
t272*t83)*t128))*t51+(t270*(t685-t686*t693*t48*t663*t696/t79)+t656*t709*t51*
t712-t685*t270*t40/(t40-t647))*t648)*t57*t60/t725*t32*t46*t49+t68*(t4*(t59*t225
+t74*t733)*t60*t646/t514/t733+(-t74*t264*t745/t453+t225*(t321*t495*t648+t655*
t321*t334*t712)/t756*t65*t49/t454-t37*(-t74*t168*t745/(-t178+t453)+t225*t770*
t65*t49*t648/(t454-t647)))*t470*t696)*t501);
  }
}

#else

return
 (

(complex<T>(0,-1)*((SPA(3,4)*pow(SPA(4,6),2)*SPB(5,3))/(SPA(2,3)*SPB(5,3)*SPB(5,4)+
SPA(2,4)*pow(SPB(5,4),2))+
(SPB(5,3)*pow(SPA(3,6)*SPB(5,3)+
SPA(4,6)*SPB(5,4),2)*((pow(SPA(3,4),2)*pow(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3),2))/(SPB(4,3)*pow(SPA(3,6)*SPB(5,3)+
SPA(4,6)*SPB(5,4),2))+
(SPA(3,4)*pow(SPA(4,5),2))/(-SPA(3,4)*SPB(4,3)+
SS(3,4,5))))/(SPB(5,4)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SS(3,4,5))))/(complex<T>(2,0)*SPA(1,2)*pow(SPA(3,4),complex<T>(2,0))*SPA(6,7))+
recursive*((complex<T>(0,1)*SPA(3,4)*SPB(5,3)*pow(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3),2)*(SPA(3,6)*SPB(5,3)+
SPA(4,6)*SPB(5,4))*(-SPA(6,7)*pow(-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4),2)*pow(SPB(7,1),2)+
(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))*((-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
SPA(6,7)*(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2))))))/(complex<T>(2,0)*SPA(1,2)*SPA(6,7)*SPB(5,4)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*(-SPA(3,4)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))*SPB(5,4)+
SPA(3,4)*SPB(5,3)*(SPA(3,6)*SPB(4,3)-SPA(5,6)*SPB(5,4)))*(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))*((-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
SPA(6,7)*(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2))))*SS(3,4,5))+
(-((complex<T>(0,1)*(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*SPB(3,2)*SPB(5,3)*pow(-SPA(2,4)*SPB(4,3)*SPB(7,5)+
SPB(5,3)*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6)),2)+
complex<T>(-3,0)*SPA(2,4)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*(-SPA(1,4)*SPB(4,3)*SPB(7,5)+
SPB(5,3)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))*(-SPA(2,4)*SPB(3,2)*SPB(4,3)*SPB(7,5)+
SPB(3,2)*SPB(5,3)*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6)))))/(complex<T>(6,0)*pow(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4),2))+
(-SPA(1,2)*SPA(2,3)*pow(SPB(7,5),2)*(-SPA(4,6)*SPB(4,3)*SPB(6,5)-SPA(5,6)*SPB(5,3)*SPB(6,5)-SPA(4,7)*SPB(4,3)*SPB(7,5)-SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*((complex<T>(0,1)*SPB(3,2)*pow(SPB(5,3),2)*(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))*(SPA(2,4)*SPB(4,3)*SPB(7,5)-SPB(5,3)*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))))/(complex<T>(3,0)*pow(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4),2)*pow(SPB(7,5),2)*(SPA(4,6)*SPB(4,3)*SPB(6,5)+
SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(4,7)*SPB(4,3)*SPB(7,5)+
SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))+
(complex<T>(0,1)*SPB(5,3)*(SPA(1,4)*SPB(4,3)+
(-SPB(5,3)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))/SPB(7,5))*((-SPB(7,1)*(SPA(1,4)*SPB(4,3)+
(-SPB(5,3)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6)))/SPB(7,5)))/(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))+
(SPB(3,2)*(-SPA(2,4)*SPB(4,3)*SPB(7,5)+
SPB(5,3)*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))))/(SPB(7,5)*(SPA(4,6)*SPB(4,3)*SPB(6,5)+
SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(4,7)*SPB(4,3)*SPB(7,5)+
SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))))/(complex<T>(2,0)*SPA(1,2)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)))+
SPB(4,3)*SPB(5,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*((complex<T>(0,-1)*(complex<T>(3,0)*SPA(2,3)*(-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1))*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*pow(SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3),2)+
complex<T>(-2,0)*pow(SPA(2,4),3)*SPA(3,4)*SPB(3,2)*(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))*pow(SPB(7,1),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))))/(complex<T>(6,0)*pow(SPA(2,3),complex<T>(2,0))*SPA(3,4)*(-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1))*pow(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4),2)*pow(SPB(7,5),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*SS(2,3,4))+
((complex<T>(0,-1)*pow(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6),2)*((complex<T>(3,0)*SPA(1,4))/(SPA(1,2)*SPA(3,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))+
(pow(SPB(3,2),2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))*((complex<T>(-3,0)*SPA(2,4))/(SPA(1,2)*(SPA(4,6)*SPB(4,3)*SPB(6,5)+
SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(4,7)*SPB(4,3)*SPB(7,5)+
SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))+
(complex<T>(2,0)*SPA(3,4)*(SPA(4,6)*SPB(6,5)+
SPA(4,7)*SPB(7,5))*(complex<T>(3,0)+
(SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))*((SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2)+
complex<T>(2,0)*SPA(3,4)*SPB(4,3))*(SPA(4,6)*SPB(6,5)+
SPA(4,7)*SPB(7,5))+
(SPA(2,4)*SPB(5,2)+
complex<T>(2,0)*SPA(3,4)*SPB(5,3))*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))/((SPA(4,6)*SPB(4,3)*SPB(6,5)+
SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(4,7)*SPB(4,3)*SPB(7,5)+
SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*((SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*(SPA(4,6)*SPB(6,5)+
SPA(4,7)*SPB(7,5))+
SPA(2,4)*SPB(5,2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))))/((-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1))*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2))))/((SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*(SPA(4,6)*SPB(6,5)+
SPA(4,7)*SPB(7,5))+
SPA(2,4)*SPB(5,2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))))/(complex<T>(6,0)*pow(SPB(7,5),2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
(complex<T>(0,1)*pow(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6),2)*((complex<T>(3,0)*SPA(1,4))/(SPA(1,2)*SPA(3,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))+
((complex<T>(2,0)*((SPB(3,2)*(SPA(4,6)*SPB(6,5)+
SPA(4,7)*SPB(7,5))*(complex<T>(-2,0)*pow(SPA(2,4),2)*pow(SPB(2,1),2)*SPB(3,1)*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)-(complex<T>(3,0)*pow(SPA(3,4),2)*pow(SPB(3,1),2)+
complex<T>(5,0)*SPA(3,4)*SPB(3,1)*(-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1))+
complex<T>(2,0)*pow(-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1),2))*SPB(3,2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))-pow(SPA(3,4),2)*SPB(3,1)*pow(SPB(3,2),2)*pow(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5),2)))/pow(SPA(2,4)*SPB(2,1)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-SPA(3,4)*SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5)),2)+
(-SPB(3,1)*(-SPB(3,2)*(SPA(1,4)*SPB(3,1)+
SPA(2,4)*SPB(3,2))*(SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5))+
(SPA(1,4)*SPB(2,1)-SPA(3,4)*SPB(3,2))*(-SPB(3,1)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5)))))/SS(1,2,3)))/((-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))+
(complex<T>(-3,0)*SPB(3,2)*(-SPA(1,2)*SPA(2,4)*SPB(2,1)*SPB(3,1)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))*(-SPA(1,2)*SPA(3,4)*SPB(3,1)+
SPA(2,4)*SS(1,2,3))))/(SPA(1,2)*(SPA(2,4)*SPB(2,1)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-SPA(3,4)*SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5)))*SS(1,2,3)))/(-SPB(3,1)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5)))))/(complex<T>(6,0)*pow(SPB(7,5),2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
(complex<T>(0,1)*pow(SPB(3,2),2)*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))*(complex<T>(-2,0)*SPA(1,2)*SPA(3,4)*SPB(3,2)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2)+
complex<T>(2,0)*SPA(3,4)*SPB(4,3))*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))+
complex<T>(3,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*(complex<T>(-2,0)*SPA(1,2)*SPA(3,4)*SPB(4,3)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))-SPA(2,4)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))*SS(2,3,4))))/(complex<T>(6,0)*SPA(1,2)*pow(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2),2)*SPB(4,3)*pow(SPB(7,5),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))*SS(2,3,4))+
(complex<T>(0,1)*(pow(SPB(7,5),2)*((complex<T>(3,0)*SPA(1,4)*pow(SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))/SPA(1,2)+
(complex<T>(-2,0)*SPA(3,4)*SPB(3,2)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*(-SPA(1,2)*(SPA(2,4)*SPB(2,1)+
SPA(3,4)*SPB(3,1))*SPB(7,5)+
(-SPA(2,4)*SPB(7,2)-SPA(3,4)*SPB(7,3))*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5)))*(-SPA(1,2)*(SPA(2,4)*SPB(2,1)+
SPA(3,4)*SPB(3,1))*SPB(7,5)+
(-SPA(2,4)*SPB(7,2)-SPA(3,4)*SPB(7,3))*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))+
complex<T>(3,0)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))))/((-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1))*pow(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4),2)))+
(complex<T>(-3,0)*SPA(1,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*pow(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6),2)*pow(-SPB(7,1)*(SPA(4,6)*SPB(6,5)+
SPA(4,7)*SPB(7,5))+
SPB(5,1)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)),2))/(SPA(1,2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))*(SS(5,6,7)-SSS(1,5,6,7)))))/(complex<T>(6,0)*SPA(3,4)*pow(SPB(7,5),4)*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*SS(2,3,4)))/SPA(2,3))))/SPB(5,4)))/(SPA(1,2)*SPA(2,3)*SPB(4,3)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*SPB(7,6)*(SPA(4,6)*SPB(4,3)*SPB(6,5)+
SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(4,7)*SPB(4,3)*SPB(7,5)+
SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))+
(complex<T>(0,1)*(-(((SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/(complex<T>(2,0)*SPA(3,4)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*(-SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))+
(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))))+
(complex<T>(0,1)*pow(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5),2)*((complex<T>(0,-1)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2))/(complex<T>(2,0)*SPA(3,4)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(-SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))+
(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))+
(complex<T>(0,1)*((SPB(3,2)*((SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3))/SPB(5,4)+
(SPA(2,4)*SPB(3,2)*(-SPA(2,6)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))-SPA(2,4)*SPA(3,6)*SPB(4,3)+
SPA(2,3)*SPA(4,6)*SPB(4,3)-SPA(5,6)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/((-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*(-SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))+
(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(-SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))))/(complex<T>(2,0)*SPB(4,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))+
(complex<T>(2,0)*SPB(3,2)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)+
(complex<T>(3,0)*SPA(3,4)*pow(SPB(3,2),2)*(SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)))*(-SPA(2,6)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))-SPA(2,4)*SPA(3,6)*SPB(4,3)+
SPA(2,3)*SPA(4,6)*SPB(4,3)-SPA(5,6)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/(-SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))+
(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(-SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))+
(-(((SPB(4,3)*SPB(5,2)*(SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)))+
SPB(3,2)*SPB(5,3)*(-SPA(2,6)*SPA(3,4)*SPB(4,2)-SPA(2,3)*SPA(4,6)*SPB(4,2)-SPA(3,6)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4))))*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3))*(-SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))+
(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/(SPB(5,4)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))+
(-SPA(3,4)*pow(SPB(3,2),3)*pow(-SPA(2,6)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))-SPA(2,4)*SPA(3,6)*SPB(4,3)+
SPA(2,3)*SPA(4,6)*SPB(4,3)-SPA(5,6)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)),2)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(pow(SPB(5,1),2)*pow(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)),2)+
pow(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)),2)*pow(SPB(7,5),2)+
complex<T>(2,0)*SPA(5,6)*(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(6,5)*SPB(7,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))+
(-pow(SPA(3,4),2)*pow(SPB(4,3),2)+
pow(SPA(5,6),2)*pow(SPB(6,5),2))*pow(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5),2)+
complex<T>(-2,0)*SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))*((-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))))/pow(-SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))+
(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(-SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)),3)))/SPB(4,3))/(complex<T>(3,0)*pow(-SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))+
(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)),2))))/complex<T>(1,0)))/complex<T>(1,0)))/(complex<T>(1,0)*SPA(2,3)*SPA(6,7)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*SSS(2,3,4,5)))+
(complex<T>(0,-1)*SPB(5,3)*((pow(SPA(1,6),complex<T>(2,0))*SPA(6,7)*pow(SPB(2,1),2)*pow(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3),2)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)))/((-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2))*((-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2))*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))+
SPA(6,7)*(SPA(1,2)*SPB(5,2)*SPB(7,1)-SPA(1,2)*SPB(5,1)*SPB(7,2)-SPB(6,5)*(SPA(1,6)*SPB(7,1)+
SPA(2,6)*SPB(7,2))-(SPA(1,7)*SPB(7,1)+
SPA(2,7)*SPB(7,2))*SPB(7,5)))*SS(3,4,5))+
(pow(SPB(3,2),2)*(pow(SPA(2,6),complex<T>(2,0))*SPA(6,7)+
(-pow(SPA(1,6),2)*pow(-SPA(2,6)*SPB(2,1)-SPA(3,6)*SPB(3,1)-SPA(4,6)*SPB(4,1)-SPA(5,6)*SPB(5,1),2)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)))/((-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2))*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*SPB(7,1)))+
(pow(SPA(1,6),complex<T>(2,0))*SPA(6,7)*pow(-(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))*SPB(6,1)-(-SPA(4,7)*SPB(4,3)-SPA(5,7)*SPB(5,3))*SPB(7,1),2))/(SS(3,4,5)*(SPA(6,7)*SPB(7,6)-SSS(2,3,4,5)))+
(-pow(SPA(2,6),complex<T>(2,0))*SPA(6,7)*pow(SPB(3,2),2)*SS(3,4,5))/(SS(3,4,5)-SSS(2,3,4,5)))/SSS(2,3,4,5)))/(complex<T>(2,0)*SPA(1,2)*pow(SPA(6,7),2)*SPB(4,3)*SPB(5,4)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)))+
(complex<T>(0,1)*((-pow(SPA(4,6),2)*(SPA(1,2)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))+
SPA(1,4)*(SPA(2,3)*SPB(5,3)+
SPA(2,4)*SPB(5,4))))/(SPA(1,2)*SPA(2,3)*SPA(6,7)*(SPA(1,2)*SPB(5,2)+
SPA(1,3)*SPB(5,3)+
SPA(1,4)*SPB(5,4))*(SPA(2,3)*SPB(5,3)+
SPA(2,4)*SPB(5,4)))+
(-((SPA(1,4)*pow(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6),2))/(SPA(1,2)*SPB(7,6)*SS(5,6,7))+
(-(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*((pow(SPA(2,4),2)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2))/SSS(2,3,4,5)+
(pow(SPA(1,6),complex<T>(2,0))*pow(SPA(2,4),2)*pow(-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1),2))/(SPA(6,7)*SPB(7,6)-SSS(2,3,4,5))))/(pow(SPA(2,4),complex<T>(2,0))*SPA(6,7)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SSS(1,5,6,7))+
pow(SPA(4,5),2)*((-SPA(1,4)*pow(SPB(7,5),2))/(SPA(1,2)*SPB(7,6)*(-SPA(6,7)*SPB(7,6)+
SS(5,6,7)))+
((SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*pow(SPA(2,6)*SPB(5,2)+
SPA(3,6)*SPB(5,3)+
SPA(4,6)*SPB(5,4),2))/(SPA(6,7)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SSS(2,3,4,5)*(SSS(1,5,6,7)-SSS(2,3,4,5))))))/(SPA(2,3)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))))/(complex<T>(2,0)*SPA(3,4))

	)

	;

#endif
}



template <class T> complex<T> R2q3g2l_qppppqbmlmlbp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
//{{qp, p, p, p, qbm, lm, lbp}, lc}

#if _VERBOSE
  _MESSAGE("R2q3g2l_eval :  qppppqbmlmlbp lc");
#endif

complex<T> recursive;
recursive=complex<T>(1,0);
//new
//here
#if _USE_OPTIMIZED

{
  complex<T> t1;
  complex<T> t108;
  complex<T> t109;
  complex<T> t11;
  complex<T> t111;
  complex<T> t116;
  complex<T> t12;
  complex<T> t121;
  complex<T> t123;
  complex<T> t124;
  complex<T> t125;
  complex<T> t129;
  complex<T> t13;
  complex<T> t136;
  complex<T> t137;
  complex<T> t139;
  complex<T> t140;
  complex<T> t141;
  complex<T> t143;
  complex<T> t144;
  complex<T> t145;
  complex<T> t147;
  complex<T> t149;
  complex<T> t15;
  complex<T> t150;
  complex<T> t151;
  complex<T> t152;
  complex<T> t156;
  complex<T> t161;
  complex<T> t164;
  complex<T> t165;
  complex<T> t166;
  complex<T> t167;
  complex<T> t17;
  complex<T> t172;
  complex<T> t176;
  complex<T> t177;
  complex<T> t19;
  complex<T> t194;
  complex<T> t2;
  complex<T> t20;
  complex<T> t201;
  complex<T> t203;
  complex<T> t21;
  complex<T> t22;
  complex<T> t221;
  complex<T> t222;
  complex<T> t223;
  complex<T> t225;
  complex<T> t226;
  complex<T> t23;
  complex<T> t232;
  complex<T> t24;
  complex<T> t242;
  complex<T> t244;
  complex<T> t247;
  complex<T> t248;
  complex<T> t25;
  complex<T> t250;
  complex<T> t251;
  complex<T> t255;
  complex<T> t259;
  complex<T> t26;
  complex<T> t261;
  complex<T> t263;
  complex<T> t264;
  complex<T> t267;
  complex<T> t272;
  complex<T> t275;
  complex<T> t277;
  complex<T> t28;
  complex<T> t284;
  complex<T> t285;
  complex<T> t29;
  complex<T> t292;
  complex<T> t295;
  complex<T> t3;
  complex<T> t31;
  complex<T> t310;
  complex<T> t314;
  complex<T> t316;
  complex<T> t320;
  complex<T> t33;
  complex<T> t337;
  complex<T> t341;
  complex<T> t348;
  complex<T> t349;
  complex<T> t35;
  complex<T> t353;
  complex<T> t355;
  complex<T> t365;
  complex<T> t37;
  complex<T> t388;
  complex<T> t389;
  complex<T> t395;
  complex<T> t396;
  complex<T> t398;
  complex<T> t4;
  complex<T> t41;
  complex<T> t410;
  complex<T> t411;
  complex<T> t412;
  complex<T> t413;
  complex<T> t414;
  complex<T> t416;
  complex<T> t42;
  complex<T> t43;
  complex<T> t440;
  complex<T> t443;
  complex<T> t448;
  complex<T> t449;
  complex<T> t460;
  complex<T> t468;
  complex<T> t5;
  complex<T> t50;
  complex<T> t51;
  complex<T> t52;
  complex<T> t53;
  complex<T> t55;
  complex<T> t57;
  complex<T> t59;
  complex<T> t6;
  complex<T> t60;
  complex<T> t61;
  complex<T> t65;
  complex<T> t66;
  complex<T> t68;
  complex<T> t71;
  complex<T> t76;
  complex<T> t78;
  complex<T> t8;
  complex<T> t80;
  complex<T> t86;
  complex<T> t88;
  complex<T> t9;
  {
    t1 = complex<T>(0,-1);
    t2 = complex<T>(3,0);
    t3 = SPA(2,3);
    t4 = t2*t3;
    t5 = SPA(1,5);
    t6 = SPB(7,5);
    t8 = SPA(1,6);
    t9 = SPB(7,6);
    t11 = -t5*t6-t8*t9;
    t12 = SPB(7,3);
    t13 = SPA(3,5);
    t15 = SPA(3,6);
    t17 = -t13*t6-t15*t9;
    t19 = SPA(5,6);
    t20 = SPB(6,5);
    t21 = t19*t20;
    t22 = SPA(5,7);
    t23 = t22*t6;
    t24 = SPA(6,7);
    t25 = t24*t9;
    t26 = t21+t23+t25;
    t28 = SPB(2,1);
    t29 = SPA(2,5);
    t31 = SPA(2,6);
    t33 = -t29*t6-t31*t9;
    t35 = SPB(3,1);
    t37 = SPB(7,1);
    t41 = SPB(3,2);
    t42 = t41*t17;
    t43 = SPB(7,2);
    t50 = SPA(1,2);
    t51 = complex<T>(1,0)/t50;
    t52 = SPA(3,4);
    t53 = complex<T>(1,0)/t52;
    t55 = SPA(4,5);
    t57 = SPA(4,6);
    t59 = -t55*t6-t57*t9;
    t60 = complex<T>(1,0)/t59;
    t61 = complex<T>(1,0)/t26;
    t65 = complex<T>(-2,0);
    t66 = SPA(1,3);
    t68 = SPB(4,3);
    t71 = pow(t52,2);
    t76 = SPB(4,2);
    t78 = t68*t17;
    t80 = pow(-t76*t33-t78,2);
    t86 = SPB(5,1);
    t88 = SPA(1,7);
    t108 = complex<T>(6,0);
    t109 = complex<T>(1,0)/t108;
    t111 = pow(t3,2);
    t116 = complex<T>(1,0)/(t8*t20+t88*t6);
    t121 = pow(t19,2);
    t123 = complex<T>(-3,0);
    t124 = t123*t3;
    t125 = SPB(5,4);
    t129 = complex<T>(2,0);
    t136 = SPA(2,4);
    t137 = t8*t28;
    t139 = t137+t24*t43;
    t140 = t136*t139;
    t141 = t8*t35;
    t143 = t141+t24*t12;
    t144 = t52*t143;
    t145 = -t140-t144;
    t147 = t31*t41;
    t149 = t57*t68;
    t150 = SPB(5,3);
    t151 = t19*t150;
    t152 = -t149-t151;
    t156 = SS(1,2,3);
    t161 = complex<T>(1,0)/t156;
    t164 = pow(t31,2);
    t165 = t164*t71;
    t166 = pow(t41,2);
    t167 = t52*t152;
    t172 = pow(t8,2);
    t176 = pow(-t136*t28-t52*t35,2);
    t177 = t172*t176;
    t194 = t3*t52;
    t201 = pow(t136*t43+t52*t12,2);
    t203 = t31*t52;
    t221 = pow(t3,2);
    t222 = complex<T>(1,0)/t221;
    t223 = complex<T>(1,0)/t55;
    t225 = complex<T>(1,0)/t24;
    t226 = complex<T>(1,0)/t125;
    t232 = t15*t68-t19*t125;
    t242 = pow(t52,2);
    t244 = complex<T>(1,0)/t242*t225;
    t247 = -t3*t150-t136*t125;
    t248 = complex<T>(1,0)/t247;
    t250 = SS(3,4,5);
    t251 = complex<T>(1,0)/t250;
    t255 = complex<T>(0,1);
    t259 = t8*t86+t24*t6;
    t261 = t136*t15;
    t263 = t136*t76;
    t264 = t52*t68;
    t267 = SPB(5,2);
    t272 = t203*t41-t261*t41-t57*(t263+t264)-t19*(t136*t267+t52*t150);
    t275 = t3*t57;
    t277 = t3*t41;
    t284 = -t203*t76-t275*t76-t15*(t277+t264)-t19*(t3*t267-t52*t125);
    t285 = t284*t259;
    t292 = -t31*(t277+t263)-t261*t68+t275*t68-t19*t247;
    t295 = SPA(1,4);
    t310 = -t50*t267-t66*t150-t295*t125;
    t314 = SPA(3,7);
    t316 = SPA(4,7);
    t320 = SPA(2,7);
    t337 = -t86*(-t31*(t66*t41+t295*t76)-t57*(-t50*t76-t66*t68)-t15*(-t50*t41+
t295*t68)-t19*t310)+(-t31*(-t314*t41-t316*t76)-t57*(t320*t76+t314*t68)-t15*(
t320*t41-t316*t68)-t19*(t320*t267+t314*t150+t316*t125))*t6+t21*t259;
    t341 = complex<T>(1,0)/t272;
    t348 = t31*t20*t259+(t50*t86+t320*t6)*t259;
    t349 = complex<T>(1,0)/t348;
    t353 = SPB(4,1);
    t355 = SPB(7,4);
    t365 = pow(t272,2);
    t388 = SSS(2,3,4,5);
    t389 = complex<T>(1,0)/t388;
    t395 = pow(t15,2);
    t396 = pow(t68,2);
    t398 = complex<T>(-1,0);
    t410 = t15*t41;
    t411 = t57*t76;
    t412 = t19*t267;
    t413 = -t410-t411-t412;
    t414 = t136*t413;
    t416 = t52*(t147-t149-t151);
    t440 = pow(t31,2);
    t443 = pow(t8,2);
    t448 = SS(2,3,4);
    t449 = pow(t448,2);
    t460 = pow(-t29*t28-t13*t35-t55*t353,2);
    t468 = pow(-t13*t41-t55*t76,2);
    return(recursive*(t1*(t4*t11*(-t12*t17*t26+t11*(t28*t33+t35*t17-t37*t26)+
t33*(t42-t43*t26))*t51*t53*t60*t61+t65*(-t66*t3*t68*t17*t59*t51/t71-t6*t80/(-
t20*(t8*t37+t19*t6)+t86*t11-t6*(t88*t37+t23+t25))+t17*(-t78*t59+t33*(-t42-t76*
t59))*t53*t60)*t61)*t109/t111*t116/t9+t1*t121*(t124*t125*t51*t53+(t129*(t125*(-
t66*t35*t37-t3*t35*t43)*t145+t52*(t147*t125+t152*t125)*t37*t156)*t60*t161+(t124
*t125*(t165*t166/(-t167-t140-t144)+t177/(t140+t144+t24*t59))+t4*t125*(t136*t52*
t24*t41*t37*(-t50*t43-t66*t12)-t194*t41*t143*t59+(t24*t201+t203*t41*t59)*t156)*
t60*t161)*t51*t53)/t145)*t109*t222*t223*t225*t226+t1*t232*(t167*t125+t123*t52*
t150*t232)*t109*t51*t244*t226*t248*t251+t255*(t4*t52*t259*(-t68*t272*t285+t292*
t139*t337)*t341*t349+t129*(-t194*t139*(t8*t353+t24*t355)-t285*(-t52*t41*t292*
t284*t348-t3*t68*t365*(t15*t20*t259+(t66*t86+t314*t6)*t259))*t341/t337*t349))*
t109*t222*t244*t116*t389)+t255*(-t395*t396*(t398*t51*t226*t251-t55*t51*t251/(-
t55*t125+t250))+(t121*(t165*t166/(t414-t167+t416)+t177/(t295*(-t31*t28-t15*t35-
t57*t353-t19*t86)-t414+t136*(t137-t410-t411-t412)-t416+t52*(t141+t147-t149-t151
)))*t223/(-t414-t416)+(-t440*t413*t310+t443*t24*t247*t37)*t449*t341/t310*t248*
t389+(-t172*t460/(-t25+t388)-t440*t468/(t250-t388))*t223*t389)*t51)/t129/t3*t53
*t225);
  }
}



#else


return

(
recursive*((complex<T>(0,-1)*((complex<T>(3,0)*SPA(2,3)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))*(-SPB(7,3)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))+
(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))*(SPB(2,1)*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))+
SPB(3,1)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))-SPB(7,1)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))*(SPB(3,2)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))-SPB(7,2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))))/(SPA(1,2)*SPA(3,4)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
(complex<T>(-2,0)*((-SPA(1,3)*SPA(2,3)*SPB(4,3)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)))/(SPA(1,2)*pow(SPA(3,4),2))+
(-SPB(7,5)*pow(-SPB(4,2)*(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))-SPB(4,3)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6)),2))/(-SPB(6,5)*(SPA(1,6)*SPB(7,1)+
SPA(5,6)*SPB(7,5))+
SPB(5,1)*(-SPA(1,5)*SPB(7,5)-SPA(1,6)*SPB(7,6))-SPB(7,5)*(SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
((-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))*(-SPB(4,3)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))+
(-SPA(2,5)*SPB(7,5)-SPA(2,6)*SPB(7,6))*(-SPB(3,2)*(-SPA(3,5)*SPB(7,5)-SPA(3,6)*SPB(7,6))-SPB(4,2)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)))))/(SPA(3,4)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)))))/(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))/(complex<T>(6,0)*pow(SPA(2,3),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*SPB(7,6))+
(complex<T>(0,-1)*pow(SPA(5,6),2)*((complex<T>(-3,0)*SPA(2,3)*SPB(5,4))/(SPA(1,2)*SPA(3,4))+
((complex<T>(2,0)*(SPB(5,4)*(-SPA(1,3)*SPB(3,1)*SPB(7,1)-SPA(2,3)*SPB(3,1)*SPB(7,2))*(-SPA(2,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-SPA(3,4)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))+
SPA(3,4)*(SPA(2,6)*SPB(3,2)*SPB(5,4)+
(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))*SPB(5,4))*SPB(7,1)*SS(1,2,3)))/((-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))*SS(1,2,3))+
(complex<T>(-3,0)*SPA(2,3)*SPB(5,4)*((pow(SPA(2,6),complex<T>(2,0))*pow(SPA(3,4),2)*pow(SPB(3,2),2))/(-SPA(3,4)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))-SPA(2,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-SPA(3,4)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))+
(pow(SPA(1,6),2)*pow(-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1),2))/(SPA(2,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
SPA(3,4)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3))+
SPA(6,7)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))))+
(complex<T>(3,0)*SPA(2,3)*SPB(5,4)*(SPA(2,4)*SPA(3,4)*SPA(6,7)*SPB(3,2)*SPB(7,1)*(-SPA(1,2)*SPB(7,2)-SPA(1,3)*SPB(7,3))-SPA(2,3)*SPA(3,4)*SPB(3,2)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3))*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))+
(SPA(6,7)*pow(SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3),2)+
SPA(2,6)*SPA(3,4)*SPB(3,2)*(-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6)))*SS(1,2,3)))/((-SPA(4,5)*SPB(7,5)-SPA(4,6)*SPB(7,6))*SS(1,2,3)))/(SPA(1,2)*SPA(3,4)))/(-SPA(2,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-SPA(3,4)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))))/(complex<T>(6,0)*pow(SPA(2,3),complex<T>(2,0))*SPA(4,5)*SPA(6,7)*SPB(5,4))+
(complex<T>(0,-1)*(SPA(3,6)*SPB(4,3)-SPA(5,6)*SPB(5,4))*(SPA(3,4)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))*SPB(5,4)+
complex<T>(-3,0)*SPA(3,4)*SPB(5,3)*(SPA(3,6)*SPB(4,3)-SPA(5,6)*SPB(5,4))))/(complex<T>(6,0)*SPA(1,2)*pow(SPA(3,4),complex<T>(2,0))*SPA(6,7)*SPB(5,4)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SS(3,4,5))+
(complex<T>(0,1)*((complex<T>(3,0)*SPA(2,3)*SPA(3,4)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(-SPB(4,3)*(SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)))*(-SPA(2,6)*SPA(3,4)*SPB(4,2)-SPA(2,3)*SPA(4,6)*SPB(4,2)-SPA(3,6)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4)))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))+
(-SPA(2,6)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))-SPA(2,4)*SPA(3,6)*SPB(4,3)+
SPA(2,3)*SPA(4,6)*SPB(4,3)-SPA(5,6)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)))*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(-SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))+
(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))))/((SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)))*(SPA(2,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))+
(SPA(1,2)*SPB(5,1)+
SPA(2,7)*SPB(7,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))+
complex<T>(2,0)*(-SPA(2,3)*SPA(3,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))+
((-SPA(2,6)*SPA(3,4)*SPB(4,2)-SPA(2,3)*SPA(4,6)*SPB(4,2)-SPA(3,6)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4)))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(-SPA(3,4)*SPB(3,2)*(-SPA(2,6)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))-SPA(2,4)*SPA(3,6)*SPB(4,3)+
SPA(2,3)*SPA(4,6)*SPB(4,3)-SPA(5,6)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4)))*(-SPA(2,6)*SPA(3,4)*SPB(4,2)-SPA(2,3)*SPA(4,6)*SPB(4,2)-SPA(3,6)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,3)*SPB(5,2)-SPA(3,4)*SPB(5,4)))*(SPA(2,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))+
(SPA(1,2)*SPB(5,1)+
SPA(2,7)*SPB(7,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))-SPA(2,3)*SPB(4,3)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)*(SPA(3,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))+
(SPA(1,3)*SPB(5,1)+
SPA(3,7)*SPB(7,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))))/((SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)))*(SPB(5,1)*(-SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))-SPA(4,6)*(-SPA(1,2)*SPB(4,2)-SPA(1,3)*SPB(4,3))-SPA(3,6)*(-SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))-SPA(5,6)*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4)))-(-SPA(2,6)*(-SPA(3,7)*SPB(3,2)-SPA(4,7)*SPB(4,2))-SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))-SPA(3,6)*(SPA(2,7)*SPB(3,2)-SPA(4,7)*SPB(4,3))-SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)-SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))*(SPA(2,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))+
(SPA(1,2)*SPB(5,1)+
SPA(2,7)*SPB(7,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))))))/(complex<T>(6,0)*pow(SPA(2,3),complex<T>(2,0))*pow(SPA(3,4),complex<T>(2,0))*SPA(6,7)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*SSS(2,3,4,5)))+
(complex<T>(0,1)*(-pow(SPA(3,6),2)*pow(SPB(4,3),2)*(complex<T>(-1,0)/(SPA(1,2)*SPB(5,4)*SS(3,4,5))+
(-SPA(4,5))/(SPA(1,2)*SS(3,4,5)*(-SPA(4,5)*SPB(5,4)+
SS(3,4,5))))+
((pow(SPA(5,6),2)*((pow(SPA(2,6),complex<T>(2,0))*pow(SPA(3,4),2)*pow(SPB(3,2),2))/(SPA(2,4)*(-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2))-SPA(3,4)*(-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))+
SPA(3,4)*(SPA(2,6)*SPB(3,2)-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3)))+
(pow(SPA(1,6),2)*pow(-SPA(2,4)*SPB(2,1)-SPA(3,4)*SPB(3,1),2))/(SPA(1,4)*(-SPA(2,6)*SPB(2,1)-SPA(3,6)*SPB(3,1)-SPA(4,6)*SPB(4,1)-SPA(5,6)*SPB(5,1))-SPA(2,4)*(-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2))+
SPA(2,4)*(SPA(1,6)*SPB(2,1)-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2))-SPA(3,4)*(SPA(2,6)*SPB(3,2)-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))+
SPA(3,4)*(SPA(1,6)*SPB(3,1)+
SPA(2,6)*SPB(3,2)-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3)))))/(SPA(4,5)*(-SPA(2,4)*(-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2))-SPA(3,4)*(SPA(2,6)*SPB(3,2)-SPA(4,6)*SPB(4,3)-SPA(5,6)*SPB(5,3))))+
((-pow(SPA(2,6),2)*(-SPA(3,6)*SPB(3,2)-SPA(4,6)*SPB(4,2)-SPA(5,6)*SPB(5,2))*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))+
pow(SPA(1,6),complex<T>(2,0))*SPA(6,7)*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SPB(7,1))*pow(SS(2,3,4),2))/((SPA(2,6)*SPA(3,4)*SPB(3,2)-SPA(2,4)*SPA(3,6)*SPB(3,2)-SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))-SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)))*(-SPA(1,2)*SPB(5,2)-SPA(1,3)*SPB(5,3)-SPA(1,4)*SPB(5,4))*(-SPA(2,3)*SPB(5,3)-SPA(2,4)*SPB(5,4))*SSS(2,3,4,5))+
(-((pow(SPA(1,6),2)*pow(-SPA(2,5)*SPB(2,1)-SPA(3,5)*SPB(3,1)-SPA(4,5)*SPB(4,1),2))/(-SPA(6,7)*SPB(7,6)+
SSS(2,3,4,5))+
(pow(SPA(2,6),2)*pow(-SPA(3,5)*SPB(3,2)-SPA(4,5)*SPB(4,2),2))/(SS(3,4,5)-SSS(2,3,4,5))))/(SPA(4,5)*SSS(2,3,4,5)))/SPA(1,2)))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*SPA(6,7))
)

;
#endif
}




 // *************** table of switch values *************

#define _R_qpmmmqbmlmlbp_lc R2q3g2l_265682_lc
#define _R_qpmmmqbmlmlbp_nf R2q3g2l_265682_nf
#define _R_qpmmpqbmlmlbp_lc R2q3g2l_266330_lc
#define _R_qpmmpqbmlmlbp_nf R2q3g2l_266330_nf
#define _R_qpmmqbmmlmlbp_nf R2q3g2l_264602_nf
#define _R_qpmmqbmmlmlbp_sl R2q3g2l_264602_sl
#define _R_qpmmqbmplmlbp_nf R2q3g2l_268490_nf
#define _R_qpmmqbmplmlbp_sl R2q3g2l_268490_sl
#define _R_qpmpmqbmlmlbp_lc R2q3g2l_265790_lc
#define _R_qpmpmqbmlmlbp_nf R2q3g2l_265790_nf
#define _R_qpmppqbmlmlbp_lc R2q3g2l_266438_lc
#define _R_qpmppqbmlmlbp_nf R2q3g2l_266438_nf
#define _R_qpmpqbmmlmlbp_nf R2q3g2l_264710_nf
#define _R_qpmpqbmmlmlbp_sl R2q3g2l_264710_sl
#define _R_qpmpqbmplmlbp_nf R2q3g2l_268598_nf
#define _R_qpmpqbmplmlbp_sl R2q3g2l_268598_sl
#define _R_qpmqbmmmlmlbp_nf R2q3g2l_264422_nf
#define _R_qpmqbmmmlmlbp_sl R2q3g2l_264422_sl
#define _R_qpmqbmmplmlbp_nf R2q3g2l_268310_nf
#define _R_qpmqbmmplmlbp_sl R2q3g2l_268310_sl
#define _R_qpmqbmpmlmlbp_nf R2q3g2l_265070_nf
#define _R_qpmqbmpmlmlbp_sl R2q3g2l_265070_sl
#define _R_qpmqbmpplmlbp_nf R2q3g2l_268958_nf
#define _R_qpmqbmpplmlbp_sl R2q3g2l_268958_sl
#define _R_qppmmqbmlmlbp_lc R2q3g2l_265700_lc
#define _R_qppmmqbmlmlbp_nf R2q3g2l_265700_nf
#define _R_qppmpqbmlmlbp_lc R2q3g2l_266348_lc
#define _R_qppmpqbmlmlbp_nf R2q3g2l_266348_nf
#define _R_qppmqbmmlmlbp_nf R2q3g2l_264620_nf
#define _R_qppmqbmmlmlbp_sl R2q3g2l_264620_sl
#define _R_qppmqbmplmlbp_nf R2q3g2l_268508_nf
#define _R_qppmqbmplmlbp_sl R2q3g2l_268508_sl
#define _R_qpppmqbmlmlbp_lc R2q3g2l_265808_lc
#define _R_qpppmqbmlmlbp_nf R2q3g2l_265808_nf
#define _R_qppppqbmlmlbp_lc R2q3g2l_266456_lc
#define _R_qppppqbmlmlbp_nf R2q3g2l_266456_nf
#define _R_qpppqbmmlmlbp_nf R2q3g2l_264728_nf
#define _R_qpppqbmmlmlbp_sl R2q3g2l_264728_sl
#define _R_qpppqbmplmlbp_nf R2q3g2l_268616_nf
#define _R_qpppqbmplmlbp_sl R2q3g2l_268616_sl
#define _R_qppqbmmmlmlbp_nf R2q3g2l_264440_nf
#define _R_qppqbmmmlmlbp_sl R2q3g2l_264440_sl
#define _R_qppqbmmplmlbp_sl R2q3g2l_268328_sl
#define _R_qppqbmpmlmlbp_nf R2q3g2l_265088_nf
#define _R_qppqbmpmlmlbp_sl R2q3g2l_265088_sl
#define _R_qppqbmpplmlbp_nf R2q3g2l_268976_nf
#define _R_qppqbmpplmlbp_sl R2q3g2l_268976_sl
#define _R_qpqbmmmmlmlbp_nf R2q3g2l_264392_nf
#define _R_qpqbmmmmlmlbp_sl R2q3g2l_264392_sl
#define _R_qpqbmmmplmlbp_nf R2q3g2l_268280_nf
#define _R_qpqbmmmplmlbp_sl R2q3g2l_268280_sl
#define _R_qpqbmmpmlmlbp_nf R2q3g2l_265040_nf
#define _R_qpqbmmpmlmlbp_sl R2q3g2l_265040_sl
#define _R_qpqbmmpplmlbp_nf R2q3g2l_268928_nf
#define _R_qpqbmmpplmlbp_sl R2q3g2l_268928_sl
#define _R_qpqbmpmmlmlbp_nf R2q3g2l_264500_nf
#define _R_qpqbmpmmlmlbp_sl R2q3g2l_264500_sl
#define _R_qpqbmpmplmlbp_nf R2q3g2l_268388_nf
#define _R_qpqbmpmplmlbp_sl R2q3g2l_268388_sl
#define _R_qpqbmppmlmlbp_nf R2q3g2l_265148_nf
#define _R_qpqbmppmlmlbp_sl R2q3g2l_265148_sl
#define _R_qpqbmppplmlbp_nf R2q3g2l_269036_nf
#define _R_qpqbmppplmlbp_sl R2q3g2l_269036_sl


 // *************** more macro definitions *************

#define _CASE_qpmmmqbmlmlbp_lc case 265682 : \
          return &R2q3g2l_265682_lc
#define _CASE_qpmmmqbmlmlbp_nf case 265682 : \
          return &R2q3g2l_265682_nf
#define _CASE_qpmmpqbmlmlbp_lc case 266330 : \
          return &R2q3g2l_266330_lc
#define _CASE_qpmmpqbmlmlbp_nf case 266330 : \
          return &R2q3g2l_266330_nf
#define _CASE_qpmmqbmmlmlbp_nf case 264602 : \
          return &R2q3g2l_264602_nf
#define _CASE_qpmmqbmmlmlbp_sl case 264602 : \
          return &R2q3g2l_264602_sl
#define _CASE_qpmmqbmplmlbp_nf case 268490 : \
          return &R2q3g2l_268490_nf
#define _CASE_qpmmqbmplmlbp_sl case 268490 : \
          return &R2q3g2l_268490_sl
#define _CASE_qpmpmqbmlmlbp_lc case 265790 : \
          return &R2q3g2l_265790_lc
#define _CASE_qpmpmqbmlmlbp_nf case 265790 : \
          return &R2q3g2l_265790_nf
#define _CASE_qpmppqbmlmlbp_lc case 266438 : \
          return &R2q3g2l_266438_lc
#define _CASE_qpmppqbmlmlbp_nf case 266438 : \
          return &R2q3g2l_266438_nf
#define _CASE_qpmpqbmmlmlbp_nf case 264710 : \
          return &R2q3g2l_264710_nf
#define _CASE_qpmpqbmmlmlbp_sl case 264710 : \
          return &R2q3g2l_264710_sl
#define _CASE_qpmpqbmplmlbp_nf case 268598 : \
          return &R2q3g2l_268598_nf
#define _CASE_qpmpqbmplmlbp_sl case 268598 : \
          return &R2q3g2l_268598_sl
#define _CASE_qpmqbmmmlmlbp_nf case 264422 : \
          return &R2q3g2l_264422_nf
#define _CASE_qpmqbmmmlmlbp_sl case 264422 : \
          return &R2q3g2l_264422_sl
#define _CASE_qpmqbmmplmlbp_nf case 268310 : \
          return &R2q3g2l_268310_nf
#define _CASE_qpmqbmmplmlbp_sl case 268310 : \
          return &R2q3g2l_268310_sl
#define _CASE_qpmqbmpmlmlbp_nf case 265070 : \
          return &R2q3g2l_265070_nf
#define _CASE_qpmqbmpmlmlbp_sl case 265070 : \
          return &R2q3g2l_265070_sl
#define _CASE_qpmqbmpplmlbp_nf case 268958 : \
          return &R2q3g2l_268958_nf
#define _CASE_qpmqbmpplmlbp_sl case 268958 : \
          return &R2q3g2l_268958_sl
#define _CASE_qppmmqbmlmlbp_lc case 265700 : \
          return &R2q3g2l_265700_lc
#define _CASE_qppmmqbmlmlbp_nf case 265700 : \
          return &R2q3g2l_265700_nf
#define _CASE_qppmpqbmlmlbp_lc case 266348 : \
          return &R2q3g2l_266348_lc
#define _CASE_qppmpqbmlmlbp_nf case 266348 : \
          return &R2q3g2l_266348_nf
#define _CASE_qppmqbmmlmlbp_nf case 264620 : \
          return &R2q3g2l_264620_nf
#define _CASE_qppmqbmmlmlbp_sl case 264620 : \
          return &R2q3g2l_264620_sl
#define _CASE_qppmqbmplmlbp_nf case 268508 : \
          return &R2q3g2l_268508_nf
#define _CASE_qppmqbmplmlbp_sl case 268508 : \
          return &R2q3g2l_268508_sl
#define _CASE_qpppmqbmlmlbp_lc case 265808 : \
          return &R2q3g2l_265808_lc
#define _CASE_qpppmqbmlmlbp_nf case 265808 : \
          return &R2q3g2l_265808_nf
#define _CASE_qppppqbmlmlbp_lc case 266456 : \
          return &R2q3g2l_266456_lc
#define _CASE_qppppqbmlmlbp_nf case 266456 : \
          return &R2q3g2l_266456_nf
#define _CASE_qpppqbmmlmlbp_nf case 264728 : \
          return &R2q3g2l_264728_nf
#define _CASE_qpppqbmmlmlbp_sl case 264728 : \
          return &R2q3g2l_264728_sl
#define _CASE_qpppqbmplmlbp_nf case 268616 : \
          return &R2q3g2l_268616_nf
#define _CASE_qpppqbmplmlbp_sl case 268616 : \
          return &R2q3g2l_268616_sl
#define _CASE_qppqbmmmlmlbp_nf case 264440 : \
          return &R2q3g2l_264440_nf
#define _CASE_qppqbmmmlmlbp_sl case 264440 : \
          return &R2q3g2l_264440_sl
#define _CASE_qppqbmmplmlbp_sl case 268328 : \
          return &R2q3g2l_268328_sl
#define _CASE_qppqbmpmlmlbp_nf case 265088 : \
          return &R2q3g2l_265088_nf
#define _CASE_qppqbmpmlmlbp_sl case 265088 : \
          return &R2q3g2l_265088_sl
#define _CASE_qppqbmpplmlbp_nf case 268976 : \
          return &R2q3g2l_268976_nf
#define _CASE_qppqbmpplmlbp_sl case 268976 : \
          return &R2q3g2l_268976_sl
#define _CASE_qpqbmmmmlmlbp_nf case 264392 : \
          return &R2q3g2l_264392_nf
#define _CASE_qpqbmmmmlmlbp_sl case 264392 : \
          return &R2q3g2l_264392_sl
#define _CASE_qpqbmmmplmlbp_nf case 268280 : \
          return &R2q3g2l_268280_nf
#define _CASE_qpqbmmmplmlbp_sl case 268280 : \
          return &R2q3g2l_268280_sl
#define _CASE_qpqbmmpmlmlbp_nf case 265040 : \
          return &R2q3g2l_265040_nf
#define _CASE_qpqbmmpmlmlbp_sl case 265040 : \
          return &R2q3g2l_265040_sl
#define _CASE_qpqbmmpplmlbp_nf case 268928 : \
          return &R2q3g2l_268928_nf
#define _CASE_qpqbmmpplmlbp_sl case 268928 : \
          return &R2q3g2l_268928_sl
#define _CASE_qpqbmpmmlmlbp_nf case 264500 : \
          return &R2q3g2l_264500_nf
#define _CASE_qpqbmpmmlmlbp_sl case 264500 : \
          return &R2q3g2l_264500_sl
#define _CASE_qpqbmpmplmlbp_nf case 268388 : \
          return &R2q3g2l_268388_nf
#define _CASE_qpqbmpmplmlbp_sl case 268388 : \
          return &R2q3g2l_268388_sl
#define _CASE_qpqbmppmlmlbp_nf case 265148 : \
          return &R2q3g2l_265148_nf
#define _CASE_qpqbmppmlmlbp_sl case 265148 : \
          return &R2q3g2l_265148_sl
#define _CASE_qpqbmppplmlbp_nf case 269036 : \
          return &R2q3g2l_269036_nf
#define _CASE_qpqbmppplmlbp_sl case 269036 : \
          return &R2q3g2l_269036_sl


 // *************** function definitions using macros *************

template <class T> complex<T> _R_qpmmmqbmlmlbp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g2l_qpmmmqbmlmlbp_lc(ep,mpc);}


template <class T> complex<T> _R_qppmmqbmlmlbp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g2l_qppmmqbmlmlbp_lc(ep,mpc);}

template <class T> complex<T> _R_qpppmqbmlmlbp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g2l_qpppmqbmlmlbp_lc(ep,mpc);}

template <class T> complex<T> _R_qppppqbmlmlbp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q3g2l_qppppqbmlmlbp_lc(ep,mpc);}


 // *************** define pointers *************

template <class T> complex<T> ( *R2q3g2l_L_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qpmmmqbmlmlbp_lc;
       _CASE_qppmmqbmlmlbp_lc;
       _CASE_qpppmqbmlmlbp_lc;
       _CASE_qppppqbmlmlbp_lc;

       default: return 0;
        }
 }


 // *************** definitions for template *************

template complex<R> ( *R2q3g2l_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q3g2l_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q3g2l_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q3g2l_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif



}

