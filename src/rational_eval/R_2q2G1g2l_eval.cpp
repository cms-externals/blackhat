/*
*R_2q3g2l_num.cpp
*
* Created on 11/13, 2008
*/

#include "R_2q2G1g2l_eval.h"
#include "eval_param.h"

using namespace std;

namespace BH  {


#define _VERBOSE 0

#define _USE_OPTIMIZED 1

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)
#define SSS(i,j,k,l) ep.s(i-1,j-1,k-1,l-1)
#define SSSS(i,j,k,l,m) ep.s(i-1,j-1,k-1,l-1,m-1)


template <class T> complex<T> R2q2Q1g2l_qppQbpQmqbmlmlbp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
//{{qp, p, Qbp, Qm, qbm, lm, lbp}, lc}

#if _VERBOSE
  _MESSAGE("R2q2Q1g2l :  qppQbpQmqbmlmlbp lc");
#endif

complex<T> recursive;
recursive=complex<T>(1,0);

#if _USE_OPTIMIZED
  complex<T> t1;
  complex<T> t10;
  complex<T> t100;
  complex<T> t1007;
  complex<T> t101;
  complex<T> t1011;
  complex<T> t1019;
  complex<T> t102;
  complex<T> t1024;
  complex<T> t1025;
  complex<T> t1028;
  complex<T> t1029;
  complex<T> t1033;
  complex<T> t1036;
  complex<T> t1038;
  complex<T> t1050;
  complex<T> t1057;
  complex<T> t1059;
  complex<T> t106;
  complex<T> t1062;
  complex<T> t1067;
  complex<T> t107;
  complex<T> t109;
  complex<T> t1090;
  complex<T> t1101;
  complex<T> t1108;
  complex<T> t1112;
  complex<T> t1128;
  complex<T> t1149;
  complex<T> t115;
  complex<T> t1150;
  complex<T> t1156;
  complex<T> t1157;
  complex<T> t1158;
  complex<T> t116;
  complex<T> t1164;
  complex<T> t1166;
  complex<T> t1168;
  complex<T> t117;
  complex<T> t118;
  complex<T> t1187;
  complex<T> t1188;
  complex<T> t1189;
  complex<T> t1193;
  complex<T> t1197;
  complex<T> t1199;
  complex<T> t12;
  complex<T> t1213;
  complex<T> t1220;
  complex<T> t1222;
  complex<T> t1232;
  complex<T> t1238;
  complex<T> t1259;
  complex<T> t1263;
  complex<T> t1264;
  complex<T> t1267;
  complex<T> t128;
  complex<T> t1281;
  complex<T> t1284;
  complex<T> t129;
  complex<T> t13;
  complex<T> t1314;
  complex<T> t132;
  complex<T> t133;
  complex<T> t136;
  complex<T> t137;
  complex<T> t138;
  complex<T> t141;
  complex<T> t142;
  complex<T> t146;
  complex<T> t148;
  complex<T> t149;
  complex<T> t15;
  complex<T> t151;
  complex<T> t152;
  complex<T> t159;
  complex<T> t16;
  complex<T> t160;
  complex<T> t163;
  complex<T> t164;
  complex<T> t167;
  complex<T> t17;
  complex<T> t170;
  complex<T> t173;
  complex<T> t182;
  complex<T> t185;
  complex<T> t187;
  complex<T> t188;
  complex<T> t189;
  complex<T> t19;
  complex<T> t190;
  complex<T> t196;
  complex<T> t197;
  complex<T> t198;
  complex<T> t2;
  complex<T> t20;
  complex<T> t204;
  complex<T> t209;
  complex<T> t21;
  complex<T> t211;
  complex<T> t217;
  complex<T> t22;
  complex<T> t221;
  complex<T> t223;
  complex<T> t228;
  complex<T> t23;
  complex<T> t247;
  complex<T> t248;
  complex<T> t249;
  complex<T> t25;
  complex<T> t250;
  complex<T> t251;
  complex<T> t252;
  complex<T> t254;
  complex<T> t258;
  complex<T> t26;
  complex<T> t262;
  complex<T> t263;
  complex<T> t264;
  complex<T> t267;
  complex<T> t272;
  complex<T> t275;
  complex<T> t278;
  complex<T> t28;
  complex<T> t282;
  complex<T> t285;
  complex<T> t289;
  complex<T> t29;
  complex<T> t290;
  complex<T> t291;
  complex<T> t292;
  complex<T> t294;
  complex<T> t296;
  complex<T> t297;
  complex<T> t299;
  complex<T> t3;
  complex<T> t301;
  complex<T> t302;
  complex<T> t306;
  complex<T> t307;
  complex<T> t309;
  complex<T> t31;
  complex<T> t311;
  complex<T> t312;
  complex<T> t313;
  complex<T> t320;
  complex<T> t322;
  complex<T> t327;
  complex<T> t331;
  complex<T> t332;
  complex<T> t335;
  complex<T> t347;
  complex<T> t35;
  complex<T> t36;
  complex<T> t360;
  complex<T> t362;
  complex<T> t365;
  complex<T> t366;
  complex<T> t373;
  complex<T> t38;
  complex<T> t39;
  complex<T> t392;
  complex<T> t393;
  complex<T> t394;
  complex<T> t395;
  complex<T> t397;
  complex<T> t4;
  complex<T> t40;
  complex<T> t405;
  complex<T> t406;
  complex<T> t412;
  complex<T> t413;
  complex<T> t414;
  complex<T> t421;
  complex<T> t424;
  complex<T> t428;
  complex<T> t433;
  complex<T> t435;
  complex<T> t437;
  complex<T> t439;
  complex<T> t443;
  complex<T> t46;
  complex<T> t460;
  complex<T> t462;
  complex<T> t463;
  complex<T> t466;
  complex<T> t467;
  complex<T> t471;
  complex<T> t472;
  complex<T> t474;
  complex<T> t475;
  complex<T> t476;
  complex<T> t477;
  complex<T> t48;
  complex<T> t480;
  complex<T> t487;
  complex<T> t488;
  complex<T> t496;
  complex<T> t498;
  complex<T> t50;
  complex<T> t503;
  complex<T> t51;
  complex<T> t52;
  complex<T> t520;
  complex<T> t526;
  complex<T> t528;
  complex<T> t541;
  complex<T> t543;
  complex<T> t549;
  complex<T> t554;
  complex<T> t556;
  complex<T> t566;
  complex<T> t57;
  complex<T> t58;
  complex<T> t584;
  complex<T> t585;
  complex<T> t6;
  complex<T> t60;
  complex<T> t600;
  complex<T> t602;
  complex<T> t61;
  complex<T> t617;
  complex<T> t618;
  complex<T> t62;
  complex<T> t620;
  complex<T> t621;
  complex<T> t624;
  complex<T> t625;
  complex<T> t632;
  complex<T> t633;
  complex<T> t65;
  complex<T> t652;
  complex<T> t66;
  complex<T> t671;
  complex<T> t682;
  complex<T> t683;
  complex<T> t685;
  complex<T> t69;
  complex<T> t694;
  complex<T> t697;
  complex<T> t699;
  complex<T> t7;
  complex<T> t70;
  complex<T> t701;
  complex<T> t711;
  complex<T> t712;
  complex<T> t715;
  complex<T> t716;
  complex<T> t717;
  complex<T> t72;
  complex<T> t724;
  complex<T> t727;
  complex<T> t731;
  complex<T> t733;
  complex<T> t741;
  complex<T> t75;
  complex<T> t76;
  complex<T> t77;
  complex<T> t770;
  complex<T> t773;
  complex<T> t78;
  complex<T> t782;
  complex<T> t785;
  complex<T> t791;
  complex<T> t797;
  complex<T> t80;
  complex<T> t811;
  complex<T> t82;
  complex<T> t820;
  complex<T> t83;
  complex<T> t835;
  complex<T> t836;
  complex<T> t84;
  complex<T> t843;
  complex<T> t847;
  complex<T> t848;
  complex<T> t851;
  complex<T> t853;
  complex<T> t858;
  complex<T> t862;
  complex<T> t863;
  complex<T> t87;
  complex<T> t872;
  complex<T> t873;
  complex<T> t879;
  complex<T> t885;
  complex<T> t889;
  complex<T> t89;
  complex<T> t890;
  complex<T> t891;
  complex<T> t895;
  complex<T> t896;
  complex<T> t9;
  complex<T> t910;
  complex<T> t912;
  complex<T> t913;
  complex<T> t92;
  complex<T> t923;
  complex<T> t93;
  complex<T> t931;
  complex<T> t933;
  complex<T> t936;
  complex<T> t937;
  complex<T> t94;
  complex<T> t945;
  complex<T> t946;
  complex<T> t947;
  complex<T> t948;
  complex<T> t95;
  complex<T> t96;
  complex<T> t960;
  complex<T> t965;
  complex<T> t97;
  complex<T> t970;
  complex<T> t98;
  complex<T> t981;
  complex<T> t983;
  complex<T> t99;
  complex<T> t993;
  complex<T> t998;
  {
    t1 = complex<T>(0,1);
    t2 = SPA(3,4);
    t3 = SPA(4,6);
    t4 = pow(t3,2);
    t6 = SPA(2,3);
    t7 = SPB(5,3);
    t9 = SPA(2,4);
    t10 = SPB(5,4);
    t12 = t6*t7+t9*t10;
    t13 = complex<T>(1,0)/t12;
    t15 = SPA(3,6);
    t16 = t15*t7;
    t17 = t3*t10;
    t19 = pow(t16+t17,2);
    t20 = pow(t2,2);
    t21 = complex<T>(-1,0);
    t22 = t21*t3;
    t23 = SPB(4,3);
    t25 = SPA(5,6);
    t26 = t21*t25;
    t28 = t22*t23+t26*t7;
    t29 = pow(t28,2);
    t31 = complex<T>(1,0)/t23;
    t35 = SPA(4,5);
    t36 = pow(t35,2);
    t38 = t21*t2;
    t39 = t38*t23;
    t40 = SS(3,4,5);
    t46 = t21*t6;
    t48 = t21*t9;
    t50 = t46*t7+t48*t10;
    t51 = complex<T>(1,0)/t50;
    t52 = complex<T>(1,0)/t40;
    t57 = complex<T>(2,0);
    t58 = complex<T>(1,0)/t57;
    t60 = SPA(1,2);
    t61 = complex<T>(1,0)/t60;
    t62 = pow(t2,2);
    t65 = SPA(6,7);
    t66 = complex<T>(1,0)/t65;
    t69 = SPB(7,5);
    t70 = pow(t69,2);
    t72 = SPB(5,2);
    t75 = t6*t72+t38*t10;
    t76 = SPB(7,2);
    t77 = t9*t76;
    t78 = SPB(7,3);
    t80 = t77+t2*t78;
    t82 = SPA(1,4);
    t83 = t21*t82;
    t84 = SPB(7,1);
    t87 = SPB(7,6);
    t89 = t83*t84+t35*t69+t3*t87;
    t92 = complex<T>(1,0)/t6;
    t93 = complex<T>(1,0)/t2;
    t94 = t92*t93;
    t95 = complex<T>(1,0)/t70;
    t96 = SPA(1,6);
    t97 = SPB(6,5);
    t98 = t96*t97;
    t99 = SPA(1,7);
    t100 = t99*t69;
    t101 = t98+t100;
    t102 = complex<T>(1,0)/t101;
    t106 = complex<T>(0,-1);
    t107 = t106*t9;
    t109 = SPB(3,2);
    t115 = t2*t7;
    t116 = t9*t72+t115;
    t117 = complex<T>(3,0);
    t118 = complex<T>(1,0)/t117;
    t128 = t21*t35*t69+t22*t87;
    t129 = complex<T>(1,0)/t69;
    t132 = t21*t72;
    t133 = SPA(2,7);
    t136 = t21*t60;
    t137 = SPA(4,7);
    t138 = t137*t84;
    t141 = SPA(2,6);
    t142 = t141*t97;
    t146 = t97*t84;
    t148 = SPB(5,1);
    t149 = t148*t128;
    t151 = t21*t141;
    t152 = t97*t128;
    t159 = t21*t7;
    t160 = SPA(3,7);
    t163 = SPA(1,3);
    t164 = t21*t163;
    t167 = t15*t97;
    t170 = t163*t3;
    t173 = t21*t15;
    t182 = complex<T>(1,0)/t116;
    t185 = pow(t21*t10*t128*t129+t21*(t132*(t82*t133*t84+t136*t138+t133*t128+
t21*(t83*t142*t84+t60*t3*t146+t136*t149+t151*t152)*t129)+t159*(t82*t160*t84+
t164*t138+t160*t128+t21*(t83*t167*t84+t170*t146+t164*t149+t173*t152)*t129))*
t182,3);
    t187 = complex<T>(1,0);
    t188 = complex<T>(1,0)/t187;
    t189 = pow(t101,2);
    t190 = complex<T>(1,0)/t189;
    t196 = t138+t21*(t22*t146+t149)*t129;
    t197 = complex<T>(1,0)/t196;
    t198 = t10*t128;
    t204 = t21*t133;
    t209 = t3*t97*t84;
    t211 = t60*t148;
    t217 = t83*t133*t84+t60*t137*t84+t204*t128+t21*(t82*t141*t146+t136*t209+
t211*t128+t142*t128)*t129;
    t221 = t163*t137;
    t223 = t21*t160;
    t228 = t163*t148;
    t247 = t25*t97;
    t248 = t99*t84;
    t249 = SPA(5,7);
    t250 = t249*t69;
    t251 = t65*t87;
    t252 = t21*t96;
    t254 = SPA(1,5);
    t258 = t21*t254*t69+t252*t87;
    t262 = t21*(t252*t146+t148*t258)*t129;
    t263 = t247+t248+t250+t251+t262;
    t264 = complex<T>(1,0)/t263;
    t267 = t1*t109;
    t272 = t82*t84;
    t275 = t7*(t163*t84+t6*t76)+t10*(t272+t77);
    t278 = SPA(2,5);
    t282 = t21*t278*t69+t151*t87;
    t285 = t9*t23+t159*t282*t129;
    t289 = pow(t50,2);
    t290 = complex<T>(1,0)/t289;
    t291 = t290*t129;
    t292 = t23*t97;
    t294 = t7*t97;
    t296 = t21*t137;
    t297 = t23*t69;
    t299 = t250+t251;
    t301 = t22*t292+t26*t294+t296*t297+t159*t299;
    t302 = complex<T>(1,0)/t301;
    t306 = t1*t82;
    t307 = SPB(2,1);
    t309 = SPB(3,1);
    t311 = t48*t307+t38*t309;
    t312 = pow(t311,2);
    t313 = pow(t258,2);
    t320 = complex<T>(1,0)/(t247+t250+t251);
    t322 = complex<T>(1,0)/(t248+t262);
    t327 = t21*t311;
    t331 = t26*t97;
    t332 = t21*t249;
    t335 = t21*t65*t87;
    t347 = complex<T>(1,0)/t9;
    t360 = complex<T>(1,0)/t109;
    t362 = t264*t197;
    t365 = t106*t84;
    t366 = t9*t109;
    t373 = t21*t99;
    t392 = complex<T>(-3,0);
    t393 = t392*t82;
    t394 = t393*t6;
    t395 = t163*t9;
    t397 = (t394+t395)*t109;
    t405 = complex<T>(6,0);
    t406 = complex<T>(1,0)/t405;
    t412 = complex<T>(-2,0);
    t413 = t412*t60;
    t414 = t9*t2;
    t421 = t82*t23;
    t424 = t421+t159*t258*t129;
    t428 = pow(t9,2);
    t433 = t7*t282;
    t435 = t60*t2;
    t437 = t392*t60;
    t439 = SPA(3,5);
    t443 = t21*t439*t69+t173*t87;
    t460 = pow(t128,2);
    t462 = t21*t109;
    t463 = t9*t307;
    t466 = t6*t109;
    t467 = t136*t307+t466;
    t471 = t133*t69;
    t472 = t142+t471;
    t474 = t160*t69;
    t475 = t167+t474;
    t476 = t463*t475;
    t477 = t437*t476;
    t480 = pow(t196,2);
    t487 = t60*t307+t46*t109;
    t488 = t487*t101;
    t496 = SPB(7,4);
    t498 = t97*t496;
    t503 = t247+t137*t496+t250+t251+t21*(t22*t498+t198)*t129;
    t520 = t117*t60;
    t526 = t412*t9;
    t528 = t109*t101;
    t541 = t2*t3*t292;
    t543 = t2*t25*t294;
    t549 = t2*t137*t297;
    t554 = t115*t299;
    t556 = pow(t97,2);
    t566 = pow(t541+t543+t373*t209+t252*t137*t97*t84+t549+t373*t138*t69+t99*
t148*t128+t554+t21*(t96*t3*t556*t84+t252*t148*t97*t128)*t129,2);
    t584 = pow(t82,2);
    t585 = pow(t84,2);
    t600 = pow(t424,2);
    t602 = complex<T>(1,0)/t275;
    t617 = t26*t97*t76+t307*t258+t21*t76*t299;
    t618 = pow(t617,2);
    t620 = t39+t247+t248+t250+t251+t262;
    t621 = t620*t95;
    t624 = t6*t2;
    t625 = pow(t109,2);
    t632 = t2*t23;
    t633 = t632+t247+t248+t250+t251+t262;
    t652 = pow(t633,2);
    t671 = pow(t620,2);
    t682 = complex<T>(9,0);
    t683 = t682*t82;
    t685 = t109*t472;
    t694 = complex<T>(1,0)/t311;
    t697 = complex<T>(23,0);
    t699 = pow(t163,2);
    t701 = complex<T>(1,0)/t503;
    t711 = t117*t163;
    t712 = t60*t6;
    t715 = t82*t307+t38*t109;
    t716 = SPB(4,2);
    t717 = t3*t716;
    t724 = t717*t97+t25*t72*t97+t137*t716*t69+t72*t299;
    t727 = t699*t311;
    t731 = complex<T>(1,0)/t20;
    t733 = pow(t503,2);
    t741 = pow(t301,2);
    t770 = t96*t9*t307*t97;
    t773 = t96*t2*t309*t97;
    t782 = t99*t9*t307*t69;
    t785 = t99*t2*t309*t69;
    t791 = t770+t773+t38*t3*t23*t97+t38*t25*t7*t97+t782+t785+t38*t137*t23*t69+
t38*t7*t299;
    t797 = pow(t309,2);
    t811 = t163*t7;
    t820 = pow(t791,2);
    t835 = complex<T>(18,0);
    t836 = complex<T>(1,0)/t835;
    t843 = complex<T>(1,0)/t87;
    t847 = complex<T>(14,0);
    t848 = complex<T>(-9,0);
    t851 = t136*t72;
    t853 = pow(t851+t98+t100,2);
    t858 = t96*t307+t65*t76;
    t862 = t136*t148+t151*t97+t204*t69;
    t863 = t858*t862;
    t872 = complex<T>(1,0)/t862;
    t873 = t31*t872;
    t879 = t96*t148+t65*t69;
    t885 = t9*t716;
    t889 = t141*t2*t109+t48*t15*t109+t22*(t885+t632)+t26*t116;
    t890 = pow(t889,2);
    t891 = t21*t890;
    t895 = t164*t148+t173*t97+t223*t69;
    t896 = t891*t895;
    t910 = t851+t164*t7+t83*t10;
    t912 = t151*(t163*t109+t82*t716)+t22*(t136*t716+t164*t23)+t173*(t136*t109+
t421)+t26*t910;
    t913 = t148*t912;
    t923 = t296*t23;
    t931 = t151*(t223*t109+t296*t716)+t22*(t133*t716+t160*t23)+t173*(t133*t109+
t923)+t26*(t133*t72+t160*t7+t137*t10);
    t933 = t21*t931*t69;
    t936 = t913+t933+t26*t97*t879;
    t937 = complex<T>(1,0)/t936;
    t945 = pow(t10,2);
    t946 = complex<T>(1,0)/t945;
    t947 = pow(t936,2);
    t948 = complex<T>(1,0)/t947;
    t960 = t151*(t466+t885)+t48*t15*t23+t6*t3*t23+t26*t50;
    t965 = t83*t148+t22*t97+t296*t69;
    t970 = complex<T>(1,0)/(t913+t933+(t632+t331)*t879);
    t981 = t151*t2*t716+t46*t717+t173*(t466+t632)+t26*t75;
    t983 = t211+t142+t471;
    t993 = t96*t309+t65*t78;
    t998 = pow(t15,2);
    t1007 = t228+t474;
    t1011 = pow(t1007,2);
    t1019 = complex<T>(1,0)/t10;
    t1024 = complex<T>(1,0)/t879;
    t1025 = t872*t1024;
    t1028 = t117*t23;
    t1029 = pow(t895,3);
    t1033 = t872*t948;
    t1036 = pow(t993,2);
    t1038 = t228+t167+t474;
    t1050 = pow(t879,2);
    t1057 = t21*t148;
    t1059 = t931*t69;
    t1062 = t1057*t912+t1059+(t39+t247)*t879;
    t1067 = complex<T>(-23,0);
    t1090 = pow(t1038,2);
    t1101 = t39+t331;
    t1108 = t228+t412*t2*t10+t167+t474;
    t1112 = t21*t97;
    t1128 = pow(t1062,2);
    t1149 = SSS(2,3,4,5);
    t1150 = complex<T>(1,0)/t1149;
    t1156 = pow(t96,2);
    t1157 = t1156*t65;
    t1158 = pow(t307,2);
    t1164 = t173*t109+t22*t716+t26*t72;
    t1166 = t50/t1164;
    t1168 = t60*t72;
    t1187 = pow(t141,2);
    t1188 = t1187*t65;
    t1189 = pow(t96,2);
    t1193 = SPB(4,1);
    t1197 = pow(t151*t307+t173*t309+t22*t1193+t26*t148,2);
    t1199 = complex<T>(1,0)/t910;
    t1213 = SPB(6,1);
    t1220 = pow(t21*t28*t1213+t21*(t923+t332*t7)*t84,2);
    t1222 = t21*t1149;
    t1232 = pow(t65,2);
    t1238 = pow(t9,2);
    t1259 = t136*t624*t109+t6*(t164*t9+t435)*t109;
    t1263 = pow(t6,2);
    t1264 = complex<T>(1,0)/t1263;
    t1267 = SS(5,6,7);
    t1281 = SSS(1,5,6,7);
    t1284 = t51/t1281*t1150;
    t1314 = pow(t141*t72+t16+t17,2);
    return(t1*(t2*t4*t13+t19*(t20*t29*t31/t19+t2*t36/(t39+t40))*t51*t52)*t58*
t61/t62*t66+recursive*(t21*t70*(t1*(-t75*t80*t89*t58*t94*t95*t102+t1*(t107*t2*
t109*t58*t92+t107*t109*t116*t118*t51)*t185*t188*t190*t197/(t198*t129+t21*(t132*
t217+t159*(t83*t160*t84+t221*t84+t223*t128+t21*(t82*t15*t146+t164*t209+t228*
t128+t167*t128)*t129))*t182))*t188*t51*t264+(t267*t7*t275*t285*t118*t31*t291*
t302+(t306*t312*t313*t58*t61*t93*t95*t320*t322*t264+t1*(t327*t258*t128*t95+t272
*t80*(t331+t332*t69+t335)*t95)*t58*t92*t322*t264)*t347+(t1*(t83*t109*t84*t89*
t196+t109*t89*t128*t196)*t58*t347*t360*t95*t362+(t365*(t366*t50*(t83*t99*t6*t84
+t60*t99*t2*t84+t373*t6*t128+t21*(t82*t96*t6*t97*t84+t136*t96*t2*t97*t84+t96*t6
*t152)*t129)+(t397*t50+t57*t6*t366*t101)*t217)*t406*t291*t362+(t106*t109*(t413*
t414*t10*t275*t285*t129+t21*t50*(t46*t2*t424*t285+t275*(t163*t428*t23+t413*t414
*t23+t21*(t395*t433+t435*t433+t437*t9*t7*t443)*t129)*t129))*t406*t31*t290*t302+
t1*t460*(t462*t189*(t21*(t164*t463+t2*t467)*t472+t477)*t480+t21*t20*t307*t301*(
(t48*t488+t397*t472)*t301+t57*t9*t109*t101*t472*t503)+t38*t101*t196*((t463*t467
*t101+t462*(t21*(t393*t6*t307+t412*t163*t463+t2*t487)*t472+t520*t476))*t301+
t526*t307*t528*t472*t503))*t406*t95*t320*t302*t197/t566)*t93)*t61)*t92)*t102+((
(t106*t82*t460*t58*t95*t320+t306*(t584*t585*t95+t57*t82*t84*t128*t95+t460*t95)*
t58*t264)*t347*t93+t21*((t365*t600*t58*t602+t267*t424*t285*t58*t302)*t51+t106*(
t136*t414*t23*t618*t621+t624*t625*t258*t282*t263*t621+t462*t617*(t526*t2*t282*
t633*(t99*t496+t21*(t252*t498+t10*t258)*t129)+t48*t620*(t624*t23*t258+t520*t443
*t263)+t21*t282*(t164*t9*t652+t394*t632*t620+t136*t2*t263*t620))*t95)*t406*t92*
t93*t320*t264/t671)*t31)*t102+t106*t460*(t21*(t683*t307*t685/(t21*t312*t101+t2*
t311*t301)+(t683*t307*t694+(t697*t163+t392*t699*t309*t701)*t92)*t93)*t102+(t711
*(t712*t715*t724+t727*t301)*t731*t190/t733+((t117*(t712*t309*t715*t101*t724+
t727*t741)*t190*t302*t701+t392*(t437*t312*t528*t475+t414*t307*t488*t301+t327*
t109*(t21*t101*(t38*t467*t472+t477)+t437*t2*t475*t301))*t102*t302/t791)*t93+
t117*(t711*t797*t302+(t699*t312*t731*t701+t463*t685*(t770+t773+t541+t543+t782+
t785+t549+t554)*(t164*t311*t101+t2*(t170*t292+t163*t25*t294+t221*t297+t811*t299
+t412*t101*t503))*t93*t302/t820)*t102))*t694)*t92)*t836*t95*t320)*t61)*t843+
t106*t29*(t847+t848*t65*t585*t853*t602/(t65*t275+t863))*t836*t61*t66*t873*t52+
t1*t879*(t896*t94*t872*t937+t879*(t117*t716*t72*t890*t895*t731*t946*t948+t117*(
t392*t109*t960*t858*t965*t879*t970+t889*(t117*t109*t10*t981*t983+t117*t6*t23*
t10*t858*t965+t993*(t6*t141*t716*t72*t97+t21*t998*t23*t556+t6*t716*t72*(t211+
t471)+t412*t15*t292*t1007+t21*t23*t1011))*t92*t93*t31*t1019)*t1025*t937+(t1028*
t890*t1029*t731*t946*t1033+t392*(t117*t1036*t1038+t109*t960*t879*(t858*t983+t46
*t109*t879)*t970)*t873/t1050+(t682*t9*t109*t981*t858*t31/t1062+t896*(t1067*t148
*t10*t912+t697*t10*t1059+t879*(t21*(t392*t15*t23+t1067*t25*t10)*t97+t1028*t1007
))*t1019*t1033)*t93*t1024+t117*(t891*t1090*t731*t1019*t937+t48*t858*(t863+t466*
t879)*t970+t366*t960*t858*(t913+t933+t1101*t879)*(t1057*t912*t1108+t1059*t1108+
t879*(t1112*(t2*(t173*t23+t57*t25*t10)+t173*t247)+t21*t1101*t1007))*t93*t31*
t937/t1128)*t1025)*t92)/t682)*t58*t66/(t252*t97+t373*t69)*t1150)+t1*(t1157*
t1158*t29*t1166/(t1164*t50+t65*(t1168*t84+t136*t148*t76+t1112*(t96*t84+t141*t76
)+t21*(t248+t133*t76)*t69))*t52+(t625*(t1188+t21*t1189*t1197*t1166*t1199/t84)+
t1188*t625*t40/(t21*t40+t1149)+t1157*t1220*t52/(t251+t1222))*t1150)*t58*t61/
t1232*t31*t51+t1*(t21*t1238*t4*(t163*t12+t60*t75)*t61*t92*t66/(t1168+t811+t82*
t10)*t13+t9*t1259*t460*t61*t1264*t360*t1199*t843/t1267+t9*t890*(t466*t116+t2*
t109*t50)*t92*t66*t360*t1199*t1284+t21*(t48*t36*t1259*t70*t61*t1264*t360*t843/(
t251+t21*t1267)+t462*t116*(t21*t360+t38*t50*t92*t360*t182)*(t1156*t428*t312*
t1149/(t335+t1149)+t1238*t36*t1314*t1281/(t1281+t1222))*t347*t66*t1284)*t1199)*
t58/t1238*t93);
  }
#else

return
(

(complex<T>(0,1)*((SPA(3,4)*pow(SPA(4,6),2))/(SPA(2,3)*SPB(5,3)+
SPA(2,4)*SPB(5,4))+
(pow(SPA(3,6)*SPB(5,3)+
SPA(4,6)*SPB(5,4),2)*((pow(SPA(3,4),2)*pow(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3),2))/(SPB(4,3)*pow(SPA(3,6)*SPB(5,3)+
SPA(4,6)*SPB(5,4),2))+
(SPA(3,4)*pow(SPA(4,5),2))/(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SS(3,4,5))))/((complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SS(3,4,5))))/(complex<T>(2,0)*SPA(1,2)*pow(SPA(3,4),complex<T>(2,0))*SPA(6,7))+
recursive*((complex<T>(-1,0)*pow(SPB(7,5),2)*((complex<T>(0,1)*(-(((SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*(SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))*(complex<T>(-1,0)*SPA(1,4)*SPB(7,1)+
SPA(4,5)*SPB(7,5)+
SPA(4,6)*SPB(7,6)))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*pow(SPB(7,5),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))))+
(complex<T>(0,1)*((complex<T>(0,-1)*SPA(2,4)*SPA(3,4)*SPB(3,2))/(complex<T>(2,0)*SPA(2,3))+
(complex<T>(0,-1)*SPA(2,4)*SPB(3,2)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)))/(complex<T>(3,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))))*pow((complex<T>(-1,0)*SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6)))/SPB(7,5)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPB(5,2)*(SPA(1,4)*SPA(2,7)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,2)*SPA(4,7)*SPB(7,1)+
SPA(2,7)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,4)*SPA(2,6)*SPB(6,5)*SPB(7,1)+
SPA(1,2)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,2)*SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))+
complex<T>(-1,0)*SPB(5,3)*(SPA(1,4)*SPA(3,7)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,3)*SPA(4,7)*SPB(7,1)+
SPA(3,7)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,4)*SPA(3,6)*SPB(6,5)*SPB(7,1)+
SPA(1,3)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
complex<T>(-1,0)*SPA(3,6)*SPB(6,5)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))))/(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),3))/(complex<T>(1,0)*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*(SPA(4,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))*((SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6)))/SPB(7,5)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPB(5,2)*(complex<T>(-1,0)*SPA(1,4)*SPA(2,7)*SPB(7,1)+
SPA(1,2)*SPA(4,7)*SPB(7,1)+
complex<T>(-1,0)*SPA(2,7)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
(complex<T>(-1,0)*(SPA(1,4)*SPA(2,6)*SPB(6,5)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,2)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPA(1,2)*SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
SPA(2,6)*SPB(6,5)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))+
complex<T>(-1,0)*SPB(5,3)*(complex<T>(-1,0)*SPA(1,4)*SPA(3,7)*SPB(7,1)+
SPA(1,3)*SPA(4,7)*SPB(7,1)+
complex<T>(-1,0)*SPA(3,7)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
(complex<T>(-1,0)*(SPA(1,4)*SPA(3,6)*SPB(6,5)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,3)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPA(1,3)*SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
SPA(3,6)*SPB(6,5)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))))/(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))))))/(complex<T>(1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5)))+
((complex<T>(0,1)*SPB(3,2)*SPB(5,3)*(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))*(SPA(2,4)*SPB(4,3)+
(complex<T>(-1,0)*SPB(5,3)*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6)))/SPB(7,5)))/(complex<T>(3,0)*SPB(4,3)*pow(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4),2)*SPB(7,5)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))+
((complex<T>(0,1)*SPA(1,4)*pow(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1),2)*pow(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6),2))/(complex<T>(2,0)*SPA(1,2)*SPA(3,4)*pow(SPB(7,5),2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))*(SPA(1,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5)))+
(complex<T>(0,1)*((complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6)))/pow(SPB(7,5),2)+
(SPA(1,4)*SPB(7,1)*(SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))*(complex<T>(-1,0)*SPA(5,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,7)*SPB(7,5)+
complex<T>(-1,0)*SPA(6,7)*SPB(7,6)))/pow(SPB(7,5),2)))/(complex<T>(2,0)*SPA(2,3)*(SPA(1,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))))/SPA(2,4)+
((complex<T>(0,1)*(complex<T>(-1,0)*SPA(1,4)*SPB(3,2)*SPB(7,1)*(complex<T>(-1,0)*SPA(1,4)*SPB(7,1)+
SPA(4,5)*SPB(7,5)+
SPA(4,6)*SPB(7,6))*(SPA(4,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))+
SPB(3,2)*(complex<T>(-1,0)*SPA(1,4)*SPB(7,1)+
SPA(4,5)*SPB(7,5)+
SPA(4,6)*SPB(7,6))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))*(SPA(4,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))))/(complex<T>(2,0)*SPA(2,4)*SPB(3,2)*pow(SPB(7,5),2)*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*(SPA(4,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5)))+
((complex<T>(0,-1)*SPB(7,1)*(SPA(2,4)*SPB(3,2)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*(complex<T>(-1,0)*SPA(1,4)*SPA(1,7)*SPA(2,3)*SPB(7,1)+
SPA(1,2)*SPA(1,7)*SPA(3,4)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,7)*SPA(2,3)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
(complex<T>(-1,0)*(SPA(1,4)*SPA(1,6)*SPA(2,3)*SPB(6,5)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,2)*SPA(1,6)*SPA(3,4)*SPB(6,5)*SPB(7,1)+
SPA(1,6)*SPA(2,3)*SPB(6,5)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))+
((complex<T>(-3,0)*SPA(1,4)*SPA(2,3)+
SPA(1,3)*SPA(2,4))*SPB(3,2)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))+
complex<T>(2,0)*SPA(2,3)*SPA(2,4)*SPB(3,2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))*(complex<T>(-1,0)*SPA(1,4)*SPA(2,7)*SPB(7,1)+
SPA(1,2)*SPA(4,7)*SPB(7,1)+
complex<T>(-1,0)*SPA(2,7)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
(complex<T>(-1,0)*(SPA(1,4)*SPA(2,6)*SPB(6,5)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,2)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPA(1,2)*SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
SPA(2,6)*SPB(6,5)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))))/(complex<T>(6,0)*pow(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4),2)*SPB(7,5)*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*(SPA(4,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5)))+
((complex<T>(0,-1)*SPB(3,2)*((complex<T>(-2,0)*SPA(1,2)*SPA(2,4)*SPA(3,4)*SPB(5,4)*(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))*(SPA(2,4)*SPB(4,3)+
(complex<T>(-1,0)*SPB(5,3)*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6)))/SPB(7,5)))/SPB(7,5)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*(complex<T>(-1,0)*SPA(2,3)*SPA(3,4)*(SPA(1,4)*SPB(4,3)+
(complex<T>(-1,0)*SPB(5,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))/SPB(7,5))*(SPA(2,4)*SPB(4,3)+
(complex<T>(-1,0)*SPB(5,3)*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6)))/SPB(7,5))+
((SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))*(SPA(1,3)*pow(SPA(2,4),2)*SPB(4,3)+
complex<T>(-2,0)*SPA(1,2)*SPA(2,4)*SPA(3,4)*SPB(4,3)+
(complex<T>(-1,0)*(SPA(1,3)*SPA(2,4)*SPB(5,3)*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6))+
SPA(1,2)*SPA(3,4)*SPB(5,3)*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6))+
complex<T>(-3,0)*SPA(1,2)*SPA(2,4)*SPB(5,3)*(complex<T>(-1,0)*SPA(3,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(3,6)*SPB(7,6))))/SPB(7,5)))/SPB(7,5))))/(complex<T>(6,0)*SPB(4,3)*pow(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4),2)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))+
(complex<T>(0,1)*pow(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6),2)*(complex<T>(-1,0)*SPB(3,2)*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,3)*SPA(2,4)*SPB(2,1)+
SPA(3,4)*(complex<T>(-1,0)*SPA(1,2)*SPB(2,1)+
SPA(2,3)*SPB(3,2)))*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))+
complex<T>(-3,0)*SPA(1,2)*SPA(2,4)*SPB(2,1)*(SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5)))*pow(SPA(4,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5),2)+
complex<T>(-1,0)*pow(SPA(3,4),2)*SPB(2,1)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*((complex<T>(-1,0)*SPA(2,4)*(SPA(1,2)*SPB(2,1)+
complex<T>(-1,0)*SPA(2,3)*SPB(3,2))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
(complex<T>(-3,0)*SPA(1,4)*SPA(2,3)+
SPA(1,3)*SPA(2,4))*SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5)))*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
complex<T>(2,0)*SPA(2,4)*SPB(3,2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))*(SPA(5,6)*SPB(6,5)+
SPA(4,7)*SPB(7,4)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,4)+
SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5)))+
complex<T>(-1,0)*SPA(3,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(4,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))*((SPA(2,4)*SPB(2,1)*(complex<T>(-1,0)*SPA(1,2)*SPB(2,1)+
SPA(2,3)*SPB(3,2))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
complex<T>(-1,0)*SPB(3,2)*(complex<T>(-1,0)*(complex<T>(-3,0)*SPA(1,4)*SPA(2,3)*SPB(2,1)+
complex<T>(-2,0)*SPA(1,3)*SPA(2,4)*SPB(2,1)+
SPA(3,4)*(SPA(1,2)*SPB(2,1)+
complex<T>(-1,0)*SPA(2,3)*SPB(3,2)))*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))+
complex<T>(3,0)*SPA(1,2)*SPA(2,4)*SPB(2,1)*(SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5))))*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
complex<T>(-2,0)*SPA(2,4)*SPB(2,1)*SPB(3,2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))*(SPA(5,6)*SPB(6,5)+
SPA(4,7)*SPB(7,4)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,4)+
SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5)))))/(complex<T>(6,0)*pow(SPB(7,5),2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*(SPA(4,7)*SPB(7,1)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5))*pow(SPA(3,4)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
SPA(3,4)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(1,7)*SPA(4,6)*SPB(6,5)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,6)*SPA(4,7)*SPB(6,5)*SPB(7,1)+
SPA(3,4)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,7)*SPA(4,7)*SPB(7,1)*SPB(7,5)+
SPA(1,7)*SPB(5,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))+
SPA(3,4)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))+
(complex<T>(-1,0)*(SPA(1,6)*SPA(4,6)*pow(SPB(6,5),2)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,6)*SPB(5,1)*SPB(6,5)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5),2)))/SPA(3,4))/SPA(1,2))/SPA(2,3))/(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
((((complex<T>(0,-1)*SPA(1,4)*pow(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6),2))/(complex<T>(2,0)*pow(SPB(7,5),2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
(complex<T>(0,1)*SPA(1,4)*((pow(SPA(1,4),2)*pow(SPB(7,1),2))/pow(SPB(7,5),2)+
(complex<T>(2,0)*SPA(1,4)*SPB(7,1)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6)))/pow(SPB(7,5),2)+
pow(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6),2)/pow(SPB(7,5),2)))/(complex<T>(2,0)*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))))/(SPA(2,4)*SPA(3,4))+
(complex<T>(-1,0)*(((complex<T>(0,-1)*SPB(7,1)*pow(SPA(1,4)*SPB(4,3)+
(complex<T>(-1,0)*SPB(5,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))/SPB(7,5),2))/(complex<T>(2,0)*(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2))))+
(complex<T>(0,1)*SPB(3,2)*(SPA(1,4)*SPB(4,3)+
(complex<T>(-1,0)*SPB(5,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))/SPB(7,5))*(SPA(2,4)*SPB(4,3)+
(complex<T>(-1,0)*SPB(5,3)*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6)))/SPB(7,5)))/(complex<T>(2,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))))/(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))+
(complex<T>(0,-1)*((complex<T>(-1,0)*SPA(1,2)*SPA(2,4)*SPA(3,4)*SPB(4,3)*pow(complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*SPB(7,2)+
SPB(2,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))+
complex<T>(-1,0)*SPB(7,2)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)),2)*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5)))/pow(SPB(7,5),2)+
(SPA(2,3)*SPA(3,4)*pow(SPB(3,2),2)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6))*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5)))/pow(SPB(7,5),2)+
(complex<T>(-1,0)*SPB(3,2)*(complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*SPB(7,2)+
SPB(2,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))+
complex<T>(-1,0)*SPB(7,2)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*(complex<T>(-2,0)*SPA(2,4)*SPA(3,4)*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6))*(SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*(SPA(1,7)*SPB(7,4)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,4)+
SPB(5,4)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))+
complex<T>(-1,0)*SPA(2,4)*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*(SPA(2,3)*SPA(3,4)*SPB(4,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))+
complex<T>(3,0)*SPA(1,2)*(complex<T>(-1,0)*SPA(3,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(3,6)*SPB(7,6))*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(2,6)*SPB(7,6))*(complex<T>(-1,0)*SPA(1,3)*SPA(2,4)*pow(SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5),2)+
complex<T>(-3,0)*SPA(1,4)*SPA(2,3)*SPA(3,4)*SPB(4,3)*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))+
complex<T>(-1,0)*SPA(1,2)*SPA(3,4)*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5)))))/pow(SPB(7,5),2)))/(complex<T>(6,0)*SPA(2,3)*SPA(3,4)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))*(SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5))*pow(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)*SPB(7,1)+
SPB(5,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/SPB(7,5),2))))/SPB(4,3))/(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
(complex<T>(0,-1)*pow(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6),2)*((complex<T>(-1,0)*((complex<T>(9,0)*SPA(1,4)*SPB(2,1)*SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5)))/(complex<T>(-1,0)*pow(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
SPA(3,4)*(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1))*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))+
((complex<T>(9,0)*SPA(1,4)*SPB(2,1))/(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1))+
(complex<T>(23,0)*SPA(1,3)+
(complex<T>(-3,0)*pow(SPA(1,3),2)*SPB(3,1))/(SPA(5,6)*SPB(6,5)+
SPA(4,7)*SPB(7,4)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,4)+
SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5)))/SPA(2,3))/SPA(3,4)))/(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
((complex<T>(3,0)*SPA(1,3)*(SPA(1,2)*SPA(2,3)*(SPA(1,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,2))*(SPA(4,6)*SPB(4,2)*SPB(6,5)+
SPA(5,6)*SPB(5,2)*SPB(6,5)+
SPA(4,7)*SPB(4,2)*SPB(7,5)+
SPB(5,2)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
pow(SPA(1,3),2)*(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1))*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))))/(pow(SPA(3,4),2)*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*pow(SPA(5,6)*SPB(6,5)+
SPA(4,7)*SPB(7,4)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,4)+
SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5),2))+
(((complex<T>(3,0)*(SPA(1,2)*SPA(2,3)*SPB(3,1)*(SPA(1,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,2))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(4,6)*SPB(4,2)*SPB(6,5)+
SPA(5,6)*SPB(5,2)*SPB(6,5)+
SPA(4,7)*SPB(4,2)*SPB(7,5)+
SPB(5,2)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
pow(SPA(1,3),2)*(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1))*pow(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)),2)))/(pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*(SPA(5,6)*SPB(6,5)+
SPA(4,7)*SPB(7,4)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,4)+
SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5)))+
(complex<T>(-3,0)*(complex<T>(-3,0)*SPA(1,2)*pow(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1),2)*SPB(3,2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5))+
SPA(2,4)*SPA(3,4)*SPB(2,1)*(SPA(1,2)*SPB(2,1)+
complex<T>(-1,0)*SPA(2,3)*SPB(3,2))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1))*SPB(3,2)*(complex<T>(-1,0)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(complex<T>(-1,0)*SPA(3,4)*(complex<T>(-1,0)*SPA(1,2)*SPB(2,1)+
SPA(2,3)*SPB(3,2))*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))+
complex<T>(-3,0)*SPA(1,2)*SPA(2,4)*SPB(2,1)*(SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5)))+
complex<T>(-3,0)*SPA(1,2)*SPA(3,4)*(SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5))*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))))/((SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*(SPA(1,6)*SPA(2,4)*SPB(2,1)*SPB(6,5)+
SPA(1,6)*SPA(3,4)*SPB(3,1)*SPB(6,5)+
complex<T>(-1,0)*SPA(3,4)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(3,4)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(1,7)*SPA(2,4)*SPB(2,1)*SPB(7,5)+
SPA(1,7)*SPA(3,4)*SPB(3,1)*SPB(7,5)+
complex<T>(-1,0)*SPA(3,4)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))))/SPA(3,4)+
complex<T>(3,0)*((complex<T>(3,0)*SPA(1,3)*pow(SPB(3,1),2))/(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))+
((pow(SPA(1,3),2)*pow(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1),2))/(pow(SPA(3,4),2)*(SPA(5,6)*SPB(6,5)+
SPA(4,7)*SPB(7,4)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,4)+
SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5)))+
(SPA(2,4)*SPB(2,1)*SPB(3,2)*(SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))*(SPA(1,6)*SPA(2,4)*SPB(2,1)*SPB(6,5)+
SPA(1,6)*SPA(3,4)*SPB(3,1)*SPB(6,5)+
SPA(3,4)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
SPA(3,4)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(1,7)*SPA(2,4)*SPB(2,1)*SPB(7,5)+
SPA(1,7)*SPA(3,4)*SPB(3,1)*SPB(7,5)+
SPA(3,4)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
SPA(3,4)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*(complex<T>(-1,0)*SPA(1,3)*(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
SPA(3,4)*(SPA(1,3)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
SPA(1,3)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(1,3)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
SPA(1,3)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))+
complex<T>(-2,0)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(5,6)*SPB(6,5)+
SPA(4,7)*SPB(7,4)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(6,5)*SPB(7,4)+
SPB(5,4)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/SPB(7,5)))))/(SPA(3,4)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)))*pow(SPA(1,6)*SPA(2,4)*SPB(2,1)*SPB(6,5)+
SPA(1,6)*SPA(3,4)*SPB(3,1)*SPB(6,5)+
complex<T>(-1,0)*SPA(3,4)*SPA(4,6)*SPB(4,3)*SPB(6,5)+
complex<T>(-1,0)*SPA(3,4)*SPA(5,6)*SPB(5,3)*SPB(6,5)+
SPA(1,7)*SPA(2,4)*SPB(2,1)*SPB(7,5)+
SPA(1,7)*SPA(3,4)*SPB(3,1)*SPB(7,5)+
complex<T>(-1,0)*SPA(3,4)*SPA(4,7)*SPB(4,3)*SPB(7,5)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,3)*(SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6)),2)))/(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))))/(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1)))/SPA(2,3)))/(complex<T>(18,0)*pow(SPB(7,5),2)*(SPA(5,6)*SPB(6,5)+
SPA(5,7)*SPB(7,5)+
SPA(6,7)*SPB(7,6))))/SPA(1,2)))/SPB(7,6)+
(complex<T>(0,-1)*pow(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3),2)*(complex<T>(14,0)+
(complex<T>(-9,0)*SPA(6,7)*pow(SPB(7,1),2)*pow(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2))/((SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))*(SPA(6,7)*(SPB(5,3)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2))+
SPB(5,4)*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)))+
(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))))))/(complex<T>(18,0)*SPA(1,2)*SPA(6,7)*SPB(4,3)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))*SS(3,4,5))+
(complex<T>(0,1)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*((complex<T>(-1,0)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)*(complex<T>(-1,0)*SPA(1,3)*SPB(5,1)+
complex<T>(-1,0)*SPA(3,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(3,7)*SPB(7,5)))/(SPA(2,3)*SPA(3,4)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))*(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))+
((SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*((complex<T>(3,0)*SPB(4,2)*SPB(5,2)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)*(complex<T>(-1,0)*SPA(1,3)*SPB(5,1)+
complex<T>(-1,0)*SPA(3,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(3,7)*SPB(7,5)))/(pow(SPA(3,4),2)*pow(SPB(5,4),2)*pow(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)),2))+
(complex<T>(3,0)*((complex<T>(-3,0)*SPB(3,2)*(complex<T>(-1,0)*SPA(2,6)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(4,3)+
SPA(2,3)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4)))*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(complex<T>(-1,0)*SPA(1,4)*SPB(5,1)+
complex<T>(-1,0)*SPA(4,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(7,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(SPA(3,4)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))+
((SPA(2,6)*SPA(3,4)*SPB(3,2)+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)))*(complex<T>(3,0)*SPB(3,2)*SPB(5,4)*(complex<T>(-1,0)*SPA(2,6)*SPA(3,4)*SPB(4,2)+
complex<T>(-1,0)*SPA(2,3)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4)))*(SPA(1,2)*SPB(5,1)+
SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))+
complex<T>(3,0)*SPA(2,3)*SPB(4,3)*SPB(5,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(complex<T>(-1,0)*SPA(1,4)*SPB(5,1)+
complex<T>(-1,0)*SPA(4,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(4,7)*SPB(7,5))+
(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3))*(SPA(2,3)*SPA(2,6)*SPB(4,2)*SPB(5,2)*SPB(6,5)+
complex<T>(-1,0)*pow(SPA(3,6),2)*SPB(4,3)*pow(SPB(6,5),2)+
SPA(2,3)*SPB(4,2)*SPB(5,2)*(SPA(1,2)*SPB(5,1)+
SPA(2,7)*SPB(7,5))+
complex<T>(-2,0)*SPA(3,6)*SPB(4,3)*SPB(6,5)*(SPA(1,3)*SPB(5,1)+
SPA(3,7)*SPB(7,5))+
complex<T>(-1,0)*SPB(4,3)*pow(SPA(1,3)*SPB(5,1)+
SPA(3,7)*SPB(7,5),2))))/(SPA(2,3)*SPA(3,4)*SPB(4,3)*SPB(5,4))))/((complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))+
((complex<T>(3,0)*SPB(4,3)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)*pow(complex<T>(-1,0)*SPA(1,3)*SPB(5,1)+
complex<T>(-1,0)*SPA(3,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(3,7)*SPB(7,5),3))/(pow(SPA(3,4),2)*pow(SPB(5,4),2)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))*pow(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)),2))+
(complex<T>(-3,0)*(complex<T>(3,0)*pow(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3),2)*(SPA(1,3)*SPB(5,1)+
SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5))+
(SPB(3,2)*(complex<T>(-1,0)*SPA(2,6)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(4,3)+
SPA(2,3)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4)))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*((SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(SPA(1,2)*SPB(5,1)+
SPA(2,6)*SPB(6,5)+
SPA(2,7)*SPB(7,5))+
complex<T>(-1,0)*SPA(2,3)*SPB(3,2)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(SPA(3,4)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))))/(SPB(4,3)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))*pow(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5),2))+
((complex<T>(9,0)*SPA(2,4)*SPB(3,2)*(complex<T>(-1,0)*SPA(2,6)*SPA(3,4)*SPB(4,2)+
complex<T>(-1,0)*SPA(2,3)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4)))*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2)))/(SPB(4,3)*(complex<T>(-1,0)*SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))+
(complex<T>(-1,0)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)*(complex<T>(-1,0)*SPA(1,3)*SPB(5,1)+
complex<T>(-1,0)*SPA(3,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(3,7)*SPB(7,5))*(complex<T>(-23,0)*SPB(5,1)*SPB(5,4)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(23,0)*SPB(5,4)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(complex<T>(-1,0)*(complex<T>(-3,0)*SPA(3,6)*SPB(4,3)+
complex<T>(-23,0)*SPA(5,6)*SPB(5,4))*SPB(6,5)+
complex<T>(3,0)*SPB(4,3)*(SPA(1,3)*SPB(5,1)+
SPA(3,7)*SPB(7,5)))))/(SPB(5,4)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))*pow(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)),2)))/(SPA(3,4)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))+
(complex<T>(3,0)*((complex<T>(-1,0)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)*pow(SPA(1,3)*SPB(5,1)+
SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5),2))/(pow(SPA(3,4),2)*SPB(5,4)*(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))+
(complex<T>(-1,0)*SPA(2,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*((SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))+
SPA(2,3)*SPB(3,2)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(SPA(3,4)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))+
(SPA(2,4)*SPB(3,2)*(complex<T>(-1,0)*SPA(2,6)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(4,3)+
SPA(2,3)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4)))*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))*(complex<T>(-1,0)*SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))*(SPA(1,3)*SPB(5,1)+
complex<T>(-2,0)*SPA(3,4)*SPB(5,4)+
SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5))+
(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)*(SPA(1,3)*SPB(5,1)+
complex<T>(-2,0)*SPA(3,4)*SPB(5,4)+
SPA(3,6)*SPB(6,5)+
SPA(3,7)*SPB(7,5))+
(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(complex<T>(-1,0)*SPB(6,5)*(SPA(3,4)*(complex<T>(-1,0)*SPA(3,6)*SPB(4,3)+
complex<T>(2,0)*SPA(5,6)*SPB(5,4))+
complex<T>(-1,0)*SPA(3,6)*SPA(5,6)*SPB(6,5))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5))*(SPA(1,3)*SPB(5,1)+
SPA(3,7)*SPB(7,5)))))/(SPA(3,4)*SPB(4,3)*(SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
complex<T>(-1,0)*SPA(5,6)*SPB(6,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))*pow(complex<T>(-1,0)*SPB(5,1)*(complex<T>(-1,0)*SPA(2,6)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4)))+
(complex<T>(-1,0)*SPA(2,6)*(complex<T>(-1,0)*SPA(3,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,2))+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,7)*SPB(4,2)+
SPA(3,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(3,6)*(SPA(2,7)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,7)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,7)*SPB(5,2)+
SPA(3,7)*SPB(5,3)+
SPA(4,7)*SPB(5,4)))*SPB(7,5)+
(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(5,6)*SPB(6,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)),2))))/((complex<T>(-1,0)*SPA(1,2)*SPB(5,1)+
complex<T>(-1,0)*SPA(2,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(2,7)*SPB(7,5))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(2,3)))/complex<T>(9,0)))/(complex<T>(2,0)*SPA(6,7)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,5))*SSS(2,3,4,5)))+
(complex<T>(0,1)*((pow(SPA(1,6),complex<T>(2,0))*SPA(6,7)*pow(SPB(2,1),2)*pow(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3),2)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4)))/((complex<T>(-1,0)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2))*((complex<T>(-1,0)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2))*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))+
SPA(6,7)*(SPA(1,2)*SPB(5,2)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,2)*SPB(5,1)*SPB(7,2)+
complex<T>(-1,0)*SPB(6,5)*(SPA(1,6)*SPB(7,1)+
SPA(2,6)*SPB(7,2))+
complex<T>(-1,0)*(SPA(1,7)*SPB(7,1)+
SPA(2,7)*SPB(7,2))*SPB(7,5)))*SS(3,4,5))+
(pow(SPB(3,2),2)*(pow(SPA(2,6),complex<T>(2,0))*SPA(6,7)+
(complex<T>(-1,0)*pow(SPA(1,6),2)*pow(complex<T>(-1,0)*SPA(2,6)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,6)*SPB(3,1)+
complex<T>(-1,0)*SPA(4,6)*SPB(4,1)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,1),2)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4)))/((complex<T>(-1,0)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2))*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1)))+
(pow(SPA(2,6),complex<T>(2,0))*SPA(6,7)*pow(SPB(3,2),2)*SS(3,4,5))/(complex<T>(-1,0)*SS(3,4,5)+
SSS(2,3,4,5))+
(pow(SPA(1,6),complex<T>(2,0))*SPA(6,7)*pow(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3))*SPB(6,1)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,7)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,7)*SPB(5,3))*SPB(7,1),2))/(SS(3,4,5)*(SPA(6,7)*SPB(7,6)+
complex<T>(-1,0)*SSS(2,3,4,5))))/SSS(2,3,4,5)))/(complex<T>(2,0)*SPA(1,2)*pow(SPA(6,7),2)*SPB(4,3)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4)))+
(complex<T>(0,1)*((complex<T>(-1,0)*pow(SPA(2,4),complex<T>(2,0))*pow(SPA(4,6),2)*(SPA(1,3)*(SPA(2,3)*SPB(5,3)+
SPA(2,4)*SPB(5,4))+
SPA(1,2)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))))/(SPA(1,2)*SPA(2,3)*SPA(6,7)*(SPA(1,2)*SPB(5,2)+
SPA(1,3)*SPB(5,3)+
SPA(1,4)*SPB(5,4))*(SPA(2,3)*SPB(5,3)+
SPA(2,4)*SPB(5,4)))+
(SPA(2,4)*(complex<T>(-1,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(3,2)+
SPA(2,3)*(complex<T>(-1,0)*SPA(1,3)*SPA(2,4)+
SPA(1,2)*SPA(3,4))*SPB(3,2))*pow(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6),2))/(SPA(1,2)*pow(SPA(2,3),2)*SPB(3,2)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,6)*SS(5,6,7))+
(SPA(2,4)*pow(SPA(2,6)*SPA(3,4)*SPB(3,2)+
complex<T>(-1,0)*SPA(2,4)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))+
complex<T>(-1,0)*SPA(5,6)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3)),2)*(SPA(2,3)*SPB(3,2)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))+
SPA(3,4)*SPB(3,2)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))))/(SPA(2,3)*SPA(6,7)*SPB(3,2)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SSS(1,5,6,7)*SSS(2,3,4,5))+
(complex<T>(-1,0)*((complex<T>(-1,0)*SPA(2,4)*pow(SPA(4,5),2)*(complex<T>(-1,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPB(3,2)+
SPA(2,3)*(complex<T>(-1,0)*SPA(1,3)*SPA(2,4)+
SPA(1,2)*SPA(3,4))*SPB(3,2))*pow(SPB(7,5),2))/(SPA(1,2)*pow(SPA(2,3),2)*SPB(3,2)*SPB(7,6)*(SPA(6,7)*SPB(7,6)+
complex<T>(-1,0)*SS(5,6,7)))+
(complex<T>(-1,0)*SPB(3,2)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*(complex<T>(-1,0)/SPB(3,2)+
(complex<T>(-1,0)*SPA(3,4)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4)))/(SPA(2,3)*SPB(3,2)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))))*((pow(SPA(1,6),complex<T>(2,0))*pow(SPA(2,4),2)*pow(complex<T>(-1,0)*SPA(2,4)*SPB(2,1)+
complex<T>(-1,0)*SPA(3,4)*SPB(3,1),2)*SSS(2,3,4,5))/(complex<T>(-1,0)*SPA(6,7)*SPB(7,6)+
SSS(2,3,4,5))+
(pow(SPA(2,4),complex<T>(2,0))*pow(SPA(4,5),2)*pow(SPA(2,6)*SPB(5,2)+
SPA(3,6)*SPB(5,3)+
SPA(4,6)*SPB(5,4),2)*SSS(1,5,6,7))/(SSS(1,5,6,7)+
complex<T>(-1,0)*SSS(2,3,4,5))))/(SPA(2,4)*SPA(6,7)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SSS(1,5,6,7)*SSS(2,3,4,5))))/(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))))/(complex<T>(2,0)*pow(SPA(2,4),complex<T>(2,0))*SPA(3,4))
)
;

#endif

      }



template <class T> complex<T> R2q2Q1g2l_qpQbpQmmqbmlmlbp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
//{{qp, Qbp, Qm, m, qbm, lm, lbp}, lc}

#if _VERBOSE
  _MESSAGE("R2q2Q1g2l :  qpQbpQmmqbmlmlbp lc");
#endif

complex<T> recursive;
recursive=complex<T>(1,0);


#if _USE_OPTIMIZED
complex<T> t1;
  complex<T> t10;
  complex<T> t100;
  complex<T> t1002;
  complex<T> t1007;
  complex<T> t1009;
  complex<T> t101;
  complex<T> t1018;
  complex<T> t102;
  complex<T> t1021;
  complex<T> t1022;
  complex<T> t1024;
  complex<T> t1033;
  complex<T> t104;
  complex<T> t1041;
  complex<T> t105;
  complex<T> t1050;
  complex<T> t1057;
  complex<T> t1058;
  complex<T> t106;
  complex<T> t1061;
  complex<T> t1070;
  complex<T> t1071;
  complex<T> t1074;
  complex<T> t1080;
  complex<T> t1084;
  complex<T> t11;
  complex<T> t110;
  complex<T> t1102;
  complex<T> t1104;
  complex<T> t111;
  complex<T> t1117;
  complex<T> t1121;
  complex<T> t1127;
  complex<T> t1133;
  complex<T> t114;
  complex<T> t115;
  complex<T> t1153;
  complex<T> t118;
  complex<T> t1180;
  complex<T> t1181;
  complex<T> t1182;
  complex<T> t1187;
  complex<T> t1191;
  complex<T> t1195;
  complex<T> t1196;
  complex<T> t1218;
  complex<T> t1223;
  complex<T> t1227;
  complex<T> t1231;
  complex<T> t1247;
  complex<T> t1251;
  complex<T> t1255;
  complex<T> t1266;
  complex<T> t128;
  complex<T> t1291;
  complex<T> t1295;
  complex<T> t1296;
  complex<T> t1299;
  complex<T> t13;
  complex<T> t130;
  complex<T> t1304;
  complex<T> t1313;
  complex<T> t1315;
  complex<T> t1345;
  complex<T> t138;
  complex<T> t14;
  complex<T> t150;
  complex<T> t153;
  complex<T> t155;
  complex<T> t156;
  complex<T> t157;
  complex<T> t164;
  complex<T> t165;
  complex<T> t169;
  complex<T> t17;
  complex<T> t171;
  complex<T> t175;
  complex<T> t177;
  complex<T> t18;
  complex<T> t181;
  complex<T> t19;
  complex<T> t195;
  complex<T> t2;
  complex<T> t20;
  complex<T> t203;
  complex<T> t21;
  complex<T> t210;
  complex<T> t211;
  complex<T> t212;
  complex<T> t213;
  complex<T> t214;
  complex<T> t215;
  complex<T> t216;
  complex<T> t220;
  complex<T> t221;
  complex<T> t223;
  complex<T> t224;
  complex<T> t225;
  complex<T> t226;
  complex<T> t229;
  complex<T> t23;
  complex<T> t231;
  complex<T> t233;
  complex<T> t235;
  complex<T> t238;
  complex<T> t239;
  complex<T> t24;
  complex<T> t240;
  complex<T> t241;
  complex<T> t248;
  complex<T> t249;
  complex<T> t250;
  complex<T> t251;
  complex<T> t252;
  complex<T> t254;
  complex<T> t257;
  complex<T> t26;
  complex<T> t261;
  complex<T> t262;
  complex<T> t263;
  complex<T> t264;
  complex<T> t275;
  complex<T> t278;
  complex<T> t28;
  complex<T> t281;
  complex<T> t284;
  complex<T> t288;
  complex<T> t289;
  complex<T> t291;
  complex<T> t293;
  complex<T> t295;
  complex<T> t297;
  complex<T> t299;
  complex<T> t3;
  complex<T> t30;
  complex<T> t300;
  complex<T> t304;
  complex<T> t305;
  complex<T> t307;
  complex<T> t308;
  complex<T> t31;
  complex<T> t310;
  complex<T> t311;
  complex<T> t313;
  complex<T> t319;
  complex<T> t32;
  complex<T> t321;
  complex<T> t332;
  complex<T> t333;
  complex<T> t346;
  complex<T> t348;
  complex<T> t358;
  complex<T> t362;
  complex<T> t364;
  complex<T> t367;
  complex<T> t368;
  complex<T> t37;
  complex<T> t371;
  complex<T> t372;
  complex<T> t375;
  complex<T> t38;
  complex<T> t382;
  complex<T> t383;
  complex<T> t394;
  complex<T> t4;
  complex<T> t40;
  complex<T> t407;
  complex<T> t408;
  complex<T> t41;
  complex<T> t414;
  complex<T> t415;
  complex<T> t42;
  complex<T> t421;
  complex<T> t422;
  complex<T> t429;
  complex<T> t430;
  complex<T> t433;
  complex<T> t447;
  complex<T> t450;
  complex<T> t461;
  complex<T> t463;
  complex<T> t465;
  complex<T> t466;
  complex<T> t467;
  complex<T> t468;
  complex<T> t469;
  complex<T> t470;
  complex<T> t471;
  complex<T> t473;
  complex<T> t475;
  complex<T> t479;
  complex<T> t480;
  complex<T> t481;
  complex<T> t487;
  complex<T> t490;
  complex<T> t498;
  complex<T> t5;
  complex<T> t505;
  complex<T> t506;
  complex<T> t51;
  complex<T> t515;
  complex<T> t52;
  complex<T> t53;
  complex<T> t543;
  complex<T> t547;
  complex<T> t549;
  complex<T> t556;
  complex<T> t566;
  complex<T> t568;
  complex<T> t58;
  complex<T> t59;
  complex<T> t590;
  complex<T> t591;
  complex<T> t6;
  complex<T> t606;
  complex<T> t608;
  complex<T> t61;
  complex<T> t617;
  complex<T> t62;
  complex<T> t622;
  complex<T> t63;
  complex<T> t632;
  complex<T> t633;
  complex<T> t649;
  complex<T> t65;
  complex<T> t650;
  complex<T> t654;
  complex<T> t66;
  complex<T> t663;
  complex<T> t664;
  complex<T> t667;
  complex<T> t684;
  complex<T> t694;
  complex<T> t695;
  complex<T> t7;
  complex<T> t70;
  complex<T> t706;
  complex<T> t708;
  complex<T> t71;
  complex<T> t711;
  complex<T> t714;
  complex<T> t715;
  complex<T> t72;
  complex<T> t728;
  complex<T> t73;
  complex<T> t732;
  complex<T> t734;
  complex<T> t735;
  complex<T> t737;
  complex<T> t739;
  complex<T> t74;
  complex<T> t742;
  complex<T> t748;
  complex<T> t752;
  complex<T> t766;
  complex<T> t77;
  complex<T> t791;
  complex<T> t794;
  complex<T> t8;
  complex<T> t802;
  complex<T> t805;
  complex<T> t808;
  complex<T> t81;
  complex<T> t815;
  complex<T> t82;
  complex<T> t831;
  complex<T> t84;
  complex<T> t840;
  complex<T> t85;
  complex<T> t856;
  complex<T> t857;
  complex<T> t86;
  complex<T> t866;
  complex<T> t869;
  complex<T> t871;
  complex<T> t873;
  complex<T> t881;
  complex<T> t884;
  complex<T> t885;
  complex<T> t894;
  complex<T> t901;
  complex<T> t906;
  complex<T> t910;
  complex<T> t914;
  complex<T> t916;
  complex<T> t918;
  complex<T> t919;
  complex<T> t924;
  complex<T> t928;
  complex<T> t93;
  complex<T> t942;
  complex<T> t943;
  complex<T> t95;
  complex<T> t955;
  complex<T> t964;
  complex<T> t965;
  complex<T> t968;
  complex<T> t969;
  complex<T> t97;
  complex<T> t976;
  complex<T> t977;
  complex<T> t979;
  complex<T> t980;
  complex<T> t988;
  complex<T> t99;
  complex<T> t997;
  {
    t1 = complex<T>(0,-1);
    t2 = complex<T>(-1,0);
    t3 = SPB(3,2);
    t4 = t2*t3;
    t5 = SPB(7,2);
    t6 = pow(t5,2);
    t7 = SPA(1,2);
    t8 = SPB(4,2);
    t10 = SPA(1,3);
    t11 = SPB(4,3);
    t13 = t7*t8+t10*t11;
    t14 = complex<T>(1,0)/t13;
    t17 = t2*t7;
    t18 = t17*t5;
    t19 = t2*t10;
    t20 = SPB(7,3);
    t21 = t19*t20;
    t23 = pow(t18+t21,2);
    t24 = pow(t3,2);
    t26 = SPB(7,1);
    t28 = SPA(2,3);
    t30 = t10*t26+t28*t5;
    t31 = pow(t30,2);
    t32 = complex<T>(1,0)/t28;
    t37 = SPB(2,1);
    t38 = pow(t37,2);
    t40 = t2*t28;
    t41 = t40*t3;
    t42 = SS(1,2,3);
    t51 = t17*t8+t19*t11;
    t52 = complex<T>(1,0)/t51;
    t53 = complex<T>(1,0)/t42;
    t58 = complex<T>(2,0);
    t59 = complex<T>(1,0)/t58;
    t61 = complex<T>(1,0)/t24;
    t62 = SPB(5,4);
    t63 = complex<T>(1,0)/t62;
    t65 = SPB(7,6);
    t66 = complex<T>(1,0)/t65;
    t70 = SPA(1,6);
    t71 = pow(t70,2);
    t72 = complex<T>(0,1);
    t73 = SPA(3,4);
    t74 = t1*t73;
    t77 = complex<T>(1,0)/t11;
    t81 = t10*t3;
    t82 = SPA(1,4);
    t84 = t81+t82*t8;
    t85 = complex<T>(3,0);
    t86 = complex<T>(1,0)/t85;
    t93 = SPA(6,7);
    t95 = t70*t37+t93*t5;
    t97 = complex<T>(1,0)/t70;
    t99 = SPA(5,6);
    t100 = SPB(5,3);
    t101 = t99*t100;
    t102 = SPB(6,2);
    t104 = t2*t99;
    t105 = SPB(5,2);
    t106 = SPB(6,3);
    t110 = SPA(1,7);
    t111 = t2*t110;
    t114 = SPA(1,5);
    t115 = t2*t114;
    t118 = t110*t99;
    t128 = t99*t62;
    t130 = SPB(6,4);
    t138 = SPB(7,4);
    t150 = complex<T>(1,0)/t84;
    t153 = pow(t7*t95*t97+t2*(-t10*(-t101*t102-t104*t105*t106-t106*t95-t2*(t111
*t101*t5+t115*t100*t95+t118*t105*t20+t111*t95*t20)*t97)-t82*(-t128*t102-t104*
t105*t130-t130*t95-t2*(t111*t128*t5+t115*t62*t95+t118*t105*t138+t111*t95*t138)*
t97))*t150,3);
    t155 = complex<T>(1,0);
    t156 = complex<T>(1,0)/t155;
    t157 = t99*t102;
    t164 = t157-t2*(-t111*t99*t5-t115*t95)*t97;
    t165 = complex<T>(1,0)/t164;
    t169 = t100*t102;
    t171 = t99*t105;
    t175 = t100*t5;
    t177 = t114*t100;
    t181 = t110*t95;
    t195 = t114*t62;
    t203 = -t104*t62*t102-t171*t130-t2*t130*t95-t2*(t118*t62*t5+t195*t95+t111*
t171*t138+t181*t138)*t97;
    t210 = SPB(6,5);
    t211 = t70*t210;
    t212 = SPB(7,5);
    t213 = t110*t212;
    t214 = t211+t213;
    t215 = pow(t214,2);
    t216 = complex<T>(1,0)/t215;
    t220 = SPA(3,6);
    t221 = t2*t220;
    t223 = SPA(4,6);
    t224 = t2*t223;
    t225 = t224*t8;
    t226 = t221*t3+t225;
    t229 = t17*t3+t82*t11;
    t231 = t2*t70;
    t233 = t2*t93;
    t235 = t231*t37+t171+t233*t5;
    t238 = complex<T>(1,0)/t71;
    t239 = complex<T>(1,0)/t3;
    t240 = t238*t239;
    t241 = complex<T>(1,0)/t214;
    t248 = SPB(6,1);
    t249 = t70*t248;
    t250 = t99*t210;
    t251 = t110*t26;
    t252 = t99*t212;
    t254 = SPB(5,1);
    t257 = t70*t254+t93*t212;
    t261 = t2*(-t111*t252-t115*t257)*t97;
    t262 = t93*t65;
    t263 = t249+t250+t251-t261+t262;
    t264 = complex<T>(1,0)/t263;
    t275 = -t17*(t225+t104*t105)-t19*(t224*t11+t104*t100);
    t278 = SPB(4,1);
    t281 = t70*t278+t93*t138;
    t284 = t28*t8+t10*t281*t97;
    t288 = pow(t51,2);
    t289 = complex<T>(1,0)/t288;
    t291 = t70*t28;
    t293 = t10*t110;
    t295 = t110*t28;
    t297 = t249+t262;
    t299 = -t291*t102-t293*t26-t295*t5-t10*t297;
    t300 = complex<T>(1,0)/t299;
    t304 = SPA(3,5);
    t305 = t2*t304;
    t307 = SPA(4,5);
    t308 = t2*t307;
    t310 = t305*t3+t308*t8;
    t311 = pow(t310,2);
    t313 = pow(t257,2);
    t319 = complex<T>(1,0)/(t250-t261);
    t321 = complex<T>(1,0)/(t249+t251+t262);
    t332 = t111*t26;
    t333 = t233*t65;
    t346 = complex<T>(1,0)/t8;
    t348 = t2*t73;
    t358 = pow(t70,2);
    t362 = complex<T>(1,0)/t73;
    t364 = t165*t264;
    t367 = t1*t99;
    t368 = complex<T>(-3,0);
    t371 = t8*t100;
    t372 = t368*t11*t105+t371;
    t375 = complex<T>(-2,0);
    t382 = t73*t8;
    t383 = t11*t105;
    t394 = t3*t62;
    t407 = complex<T>(6,0);
    t408 = complex<T>(1,0)/t407;
    t414 = t375*t7;
    t415 = t3*t8;
    t421 = t2*t51;
    t422 = pow(t8,2);
    t429 = t8*t62;
    t430 = SPB(3,1);
    t433 = t70*t430+t93*t20;
    t447 = t28*t105;
    t450 = t447+t10*t257*t97;
    t461 = pow(t95,2);
    t463 = pow(t164,2);
    t465 = t85*t307;
    t466 = t70*t106;
    t467 = t110*t20;
    t468 = t466+t467;
    t469 = t429*t468;
    t470 = t465*t469;
    t471 = t307*t8;
    t473 = t73*t11;
    t475 = t473+t308*t62;
    t479 = t70*t130;
    t480 = t110*t138;
    t481 = t479+t480;
    t487 = t58*t73;
    t490 = SPA(2,6);
    t498 = t249+t490*t102+t251-t2*(-t111*t490*t5-t17*t95)*t97+t262;
    t505 = t348*t11+t307*t62;
    t506 = t8*t505;
    t515 = t8*t481;
    t543 = t291*t3*t102;
    t547 = t293*t3*t26;
    t549 = t295*t3*t5;
    t556 = pow(t110,2);
    t566 = t81*t297;
    t568 = pow(t543+t231*t157*t210+t547+t549+t111*t250*t5+t115*t210*t95+t111*
t157*t212-t2*(-t556*t99*t5*t212-t114*t110*t95*t212)*t97+t566,2);
    t590 = pow(t99,2);
    t591 = pow(t105,2);
    t606 = pow(t450,2);
    t608 = complex<T>(1,0)/t275;
    t617 = pow(t73,2);
    t622 = t41+t249+t250+t251-t261+t262;
    t632 = -t111*t223*t26-t308*t257-t224*t297;
    t633 = pow(t632,2);
    t649 = t28*t3;
    t650 = t649+t249+t250+t251-t261+t262;
    t654 = t3*t11;
    t663 = t2*t8;
    t664 = pow(t650,2);
    t667 = t368*t28;
    t684 = pow(t622,2);
    t694 = complex<T>(-9,0);
    t695 = t694*t73;
    t706 = complex<T>(9,0);
    t708 = complex<T>(1,0)/t310;
    t711 = complex<T>(-23,0);
    t714 = pow(t100,2);
    t715 = complex<T>(1,0)/t498;
    t728 = t310*t714;
    t732 = t348*t3+t307*t105;
    t734 = SPA(2,4);
    t735 = t734*t102;
    t737 = t2*t82;
    t739 = t734*t5;
    t742 = -t231*t735-t737*t251-t111*t739-t737*t297;
    t748 = pow(t498,2);
    t752 = pow(t299,2);
    t766 = t62*t468;
    t791 = t70*t304*t3*t210;
    t794 = t70*t307*t8*t210;
    t802 = t110*t304*t3*t212;
    t805 = t110*t307*t8*t212;
    t808 = t231*t649*t102+t791+t794+t19*t110*t3*t26+t111*t649*t5+t802+t805+t19*
t3*t297;
    t815 = pow(t304,2);
    t831 = t10*t100;
    t840 = pow(t808,2);
    t856 = complex<T>(18,0);
    t857 = complex<T>(1,0)/t856;
    t866 = complex<T>(1,0)/t93;
    t869 = complex<T>(14,0);
    t871 = t737*t62;
    t873 = pow(t871+t211+t213,2);
    t881 = t115*t62+t231*t130+t111*t138;
    t884 = t308*t212+t224*t65;
    t885 = t881*t884;
    t894 = complex<T>(1,0)/t881;
    t901 = t115*t212+t231*t65;
    t906 = t115*t100+t231*t106+t111*t20;
    t910 = t734*t8;
    t914 = t8*t20;
    t916 = t73*t3;
    t918 = -t2*t84*t26-t2*(t649+t910)*t5-t348*t914-t916*t138;
    t919 = pow(t918,2);
    t924 = t17*t105+t19*t100+t871;
    t928 = t2*t734;
    t942 = -t2*t924*t26-t2*(t40*t100+t928*t62)*t5-t2*(t447+t348*t62)*t20-t2*(
t734*t105+t73*t100)*t138;
    t943 = t115*t942;
    t955 = t28*t102;
    t964 = -t2*(t17*t102+t19*t106+t737*t130)*t26-t2*(t40*t106+t928*t130)*t5-t2*
(t955+t348*t130)*t20-t2*(t735+t73*t106)*t138;
    t965 = t231*t964;
    t968 = -t943-t965+t111*t26*t901;
    t969 = complex<T>(1,0)/t968;
    t976 = pow(t7,2);
    t977 = complex<T>(1,0)/t976;
    t979 = pow(t968,2);
    t980 = complex<T>(1,0)/t979;
    t988 = t115*t105+t231*t102+t111*t5;
    t997 = -t421*t26-t28*t11*t5-t40*t914-t2*(t910+t473)*t138;
    t1002 = complex<T>(1,0)/(-t943-t965+(t649+t332)*t901);
    t1007 = t195+t479+t480;
    t1009 = t2*t229;
    t1018 = -t1009*t26-t928*t11*t5-t2*(t649+t473)*t20-t928*t3*t138;
    t1021 = t177+t466;
    t1022 = pow(t1021,2);
    t1024 = t734*t11;
    t1033 = pow(t20,2);
    t1041 = t305*t212+t221*t65;
    t1050 = complex<T>(1,0)/t7;
    t1057 = complex<T>(1,0)/t901;
    t1058 = t894*t1057;
    t1061 = pow(t906,3);
    t1070 = t114*t942;
    t1071 = t70*t964;
    t1074 = -t1070-t1071+(t41+t251)*t901;
    t1080 = t711*t7;
    t1084 = complex<T>(23,0);
    t1102 = t177+t466+t467;
    t1104 = pow(t1041,2);
    t1117 = pow(t901,2);
    t1121 = pow(t1102,2);
    t1127 = t41+t332;
    t1133 = t414*t3+t177+t466+t467;
    t1153 = pow(t1074,2);
    t1180 = SSS(1,2,3,4);
    t1181 = complex<T>(1,0)/t1180;
    t1182 = t66*t1181;
    t1187 = pow(t307,2);
    t1191 = pow(t212,2);
    t1195 = t82*t26+t739+t73*t20;
    t1196 = complex<T>(1,0)/t1195;
    t1218 = SPA(2,5);
    t1223 = pow(t114*t26+t1218*t5+t304*t20+t307*t138,2);
    t1227 = complex<T>(1,0)/t924;
    t1231 = pow(t138,2);
    t1247 = SPA(5,7);
    t1251 = pow(-t104*(t10*t248+t955)-t2*t1247*t30,2);
    t1255 = t2*t1180;
    t1266 = pow(t65,2);
    t1291 = t348*t654*t62+t473*(t663*t100+t394);
    t1295 = pow(t11,2);
    t1296 = complex<T>(1,0)/t1295;
    t1299 = SS(1,6,7);
    t1304 = t73*t84;
    t1313 = SSS(1,5,6,7);
    t1315 = t1182/t1313;
    t1345 = pow(t18+t21+t737*t138,2);
    return(t1*(-t4*t6*t14+t23*(-t2*t24*t31*t32/t23-t2*t38*t3/(t41+t42))*t52*t53
)*t59*t61*t63*t66+recursive*t2*(-t71*(t72*(t72*(t74*t3*t8*t59*t77+t74*t8*t84*
t86*t52)*t153*t156*t165/(t17*t95*t97+t2*(-t10*(-t104*t169-t171*t106-t2*t106*t95
-t2*(t118*t175+t177*t95+t111*t171*t20+t181*t20)*t97)-t82*t203)*t150)*t216-t226*
t229*t235*t59*t240*t77*t241)*t156*t52*t264+(t1*t10*t73*t275*t284*t86*t97*t32*
t289*t300-t2*(-t1*t311*t105*t313*t59*t240*t63*t319*t321*t264-t1*(t2*t310*t95*
t257*t238+t104*t226*t105*(t231*t248+t332+t333)*t238)*t59*t77*t319*t264)*t346-t2
*(t72*(-t348*t99*t105*t235*t164-t348*t95*t235*t164)*t59/t358*t362*t346*t364-t2*
(t367*(t203*(-t348*t51*t372-t375*t73*t8*t11*t214)+t382*t51*(t104*t383*t210+t99*
t3*t62*t210+t11*t210*t95-t2*(-t118*t383*t212-t111*t99*t394*t212-t111*t11*t95*
t212)*t97))*t408*t97*t289*t364-t2*(t74*(-t414*t415*t275*t62*t284*t97+t421*(-
t275*(t28*t422*t100+t375*t28*t415*t62-t2*(-t85*t10*t429*t433-t19*t371*t281-t19*
t394*t281)*t97)*t97+t4*t11*t284*t450))*t408*t32*t289*t300+t72*t461*(-t73*t463*(
-t470+t2*(-t471*t100-t4*t475)*t481)*t215-t307*t24*t299*(t487*t8*t481*t214*t498+
(-t348*t372*t481-t506*t214)*t299)-t3*t164*t214*(-t487*t307*t515*t214*t498+(-t73
*(-t368*t307*t469+t2*(-t465*t383-t58*t307*t371-t4*t505)*t481)+t471*t475*t214)*
t299))*t408*t238*t165*t321*t300/t568)*t239)*t63)*t77)*t241-t2*(((-t72*t105*t461
*t59*t238*t321-t1*t105*(t590*t591*t238+t375*t99*t105*t95*t238+t461*t238)*t59*
t264)*t239*t346-((-t367*t606*t59*t608-t74*t284*t450*t59*t300)*t52+t1*(t617*t3*
t11*t281*t257*t263*t622*t238+t40*t415*t62*t622*t633*t238-t73*t632*(t375*t3*t8*
t281*(t490*t210-t2*(-t111*t490*t212-t17*t257)*t97)*t650-t8*t622*(-t40*t654*t257
-t368*t62*t433*t263)+t2*t281*(t663*t100*t664+t667*t3*t383*t622+t4*t62*t263*t622
))*t238)*t408*t239*t77*t321*t264/t684)*t32)*t241+t1*t461*(t2*(-t695*t307*t105*
t481/(t2*t311*t214-t4*t310*t299)-t2*(t706*t307*t105*t708-t2*(-t711*t100-t85*
t304*t714*t715)*t77)*t239)*t241-t2*(-t368*t100*(t728*t299+t11*t732*t62*t742)*
t61*t216/t748+(-t2*(t85*(t728*t752-t305*t11*t732*t62*t214*t742)*t216*t715*t300+
t368*(t368*t73*t311*t766*t214-t308*t3*t506*t214*t299-t73*t310*(t2*(-t470-t3*
t475*t481)*t214+t368*t3*t766*t299))*t241*t300/t808)*t239+t85*(-t368*t815*t100*
t300+(t311*t714*t61*t715+t73*t307*t515*(t543+t791+t794+t547+t549+t802+t805+t566
)*(-t310*t100*t214-t4*(t291*t169+t293*t100*t26+t295*t175+t831*t297+t375*t214*
t498))*t239*t300/t840)*t241))*t708)*t77)*t857*t238*t321)*t63)*t866-t72*t31*(
t869-t706*t590*t873*t65*t608/(-t2*t275*t65+t885))*t857*t32*t63*t894*t66*t53-t1*
t901*(t2*t906*t919*t239*t77*t894*t969+t901*(t85*t82*t734*t906*t919*t977*t61*
t980+t85*(-t85*t73*t988*t997*t901*t884*t1002+t918*(t85*t7*t73*t1007*t1018+(-t28
*t1022-t737*t1024*(t195+t479)-t58*t110*t28*t1021*t20-t556*t28*t1033-t737*t110*
t1024*t138)*t1041-t368*t7*t28*t11*t988*t884)*t1050*t32*t239*t77)*t1058*t969-t2*
(-t667*t1061*t919*t977*t61*t894*t980-t2*(-t695*t8*t1018*t884*t32/t1074-t906*
t919*(t1080*t1070+t1080*t1071+(-t667*t1021-t111*(t1084*t7*t26+t85*t28*t20))*
t901)*t1050*t894*t980)*t239*t1057-t85*(t85*t1102*t1104-t348*t997*t901*(t348*t11
*t901+t1007*t884)*t1002)*t32*t894/t1117+t85*(-t1121*t919*t1050*t61*t969+t382*
t997*t884*(-t943-t965+t1127*t901)*(-t114*t1133*t942-t70*t1133*t964+(t2*t1021*
t1127-t111*(-t111*t26*t20-t4*(t414*t26+t28*t20)))*t901)*t32*t239*t969/t1153-t8*
t884*(t473*t901+t885)*t1002)*t1058)*t77)/t706)*t59/(t231*t210+t111*t212)*t1182)
+t1*(-t2*t1187*t51*t31*t1191*t65*t1196/(t51*t1195-t2*(-t115*t223*t62-t82*t99*
t62-t231*(t223*t130+t250)-t111*(t223*t138+t252))*t65)*t53+(t617*(-t421*t1223*
t1191/t99*t1227*t1196-t2*t1231*t65)-t2*t617*t1231*t65*t42/(t2*t42+t1180)-t2*
t1251*t1191*t65*t53/(t262+t1255))*t1181)*t59*t32*t52*t63/t1266-t72*(-t422*(-t2*
t13*t100-t1009*t62)*t6*t77*t14*t63/(t7*t105+t831+t82*t62)*t66+t8*t1291*t461*
t362*t866*t1296*t63*t1227/t1299+t8*(t1304*t11+t916*t51)*t919*t362*t77*t52*t1227
*t1315+t2*(t2*t71*t38*t8*t1291*t362*t866*t1296*t63/(t262+t2*t1299)-t1304*(-t155
*t362-t3*t51*t362*t150*t77)*(t422*t311*t1191*t1180/(t333+t1180)+t38*t422*t1345*
t1313/(t1255+t1313))*t346*t52*t1315)*t1227)*t59*t239/t422);
  }
#else

return
(


(complex<T>(0,-1)*(-((complex<T>(-1,0)*SPB(3,2)*pow(SPB(7,2),2))/(SPA(1,2)*SPB(4,2)+
SPA(1,3)*SPB(4,3)))+
(pow(complex<T>(-1,0)*SPA(1,2)*SPB(7,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(7,3),2)*(-((complex<T>(-1,0)*pow(SPB(3,2),2)*pow(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2),2))/(SPA(2,3)*pow(complex<T>(-1,0)*SPA(1,2)*SPB(7,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(7,3),2)))-
(complex<T>(-1,0)*pow(SPB(2,1),2)*SPB(3,2))/(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SS(1,2,3))))/((complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*SS(1,2,3))))/(complex<T>(2,0)*pow(SPB(3,2),2)*SPB(5,4)*SPB(7,6))+
recursive*complex<T>(-1,0)*(-((pow(SPA(1,6),2)*((complex<T>(0,1)*((complex<T>(0,1)*((complex<T>(0,-1)*SPA(3,4)*SPB(3,2)*SPB(4,2))/(complex<T>(2,0)*SPB(4,3))+
(complex<T>(0,-1)*SPA(3,4)*SPB(4,2)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2)))/(complex<T>(3,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))))*pow((SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2)))/SPA(1,6)+
(complex<T>(-1,0)*(-(SPA(1,3)*(-(SPA(5,6)*SPB(5,3)*SPB(6,2))-
complex<T>(-1,0)*SPA(5,6)*SPB(5,2)*SPB(6,3)-
SPB(6,3)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(5,3)*SPB(7,2)+
complex<T>(-1,0)*SPA(1,5)*SPB(5,3)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
SPA(1,7)*SPA(5,6)*SPB(5,2)*SPB(7,3)+
complex<T>(-1,0)*SPA(1,7)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,3)))/SPA(1,6)))-
SPA(1,4)*(-(SPA(5,6)*SPB(5,4)*SPB(6,2))-
complex<T>(-1,0)*SPA(5,6)*SPB(5,2)*SPB(6,4)-
SPB(6,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(5,4)*SPB(7,2)+
complex<T>(-1,0)*SPA(1,5)*SPB(5,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
SPA(1,7)*SPA(5,6)*SPB(5,2)*SPB(7,4)+
complex<T>(-1,0)*SPA(1,7)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,4)))/SPA(1,6))))/(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2)),3))/(complex<T>(1,0)*(SPA(5,6)*SPB(6,2)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6))*((complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2)))/SPA(1,6)+
(complex<T>(-1,0)*(-(SPA(1,3)*(-(complex<T>(-1,0)*SPA(5,6)*SPB(5,3)*SPB(6,2))-
SPA(5,6)*SPB(5,2)*SPB(6,3)-
complex<T>(-1,0)*SPB(6,3)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-
(complex<T>(-1,0)*(SPA(1,7)*SPA(5,6)*SPB(5,3)*SPB(7,2)+
SPA(1,5)*SPB(5,3)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(5,2)*SPB(7,3)+
SPA(1,7)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,3)))/SPA(1,6)))-
SPA(1,4)*(-(complex<T>(-1,0)*SPA(5,6)*SPB(5,4)*SPB(6,2))-
SPA(5,6)*SPB(5,2)*SPB(6,4)-
complex<T>(-1,0)*SPB(6,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-
(complex<T>(-1,0)*(SPA(1,7)*SPA(5,6)*SPB(5,4)*SPB(7,2)+
SPA(1,5)*SPB(5,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(5,2)*SPB(7,4)+
SPA(1,7)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,4)))/SPA(1,6))))/(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2)))*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2))-
((complex<T>(-1,0)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*SPB(4,2))*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*(complex<T>(-1,0)*SPA(1,6)*SPB(2,1)+
SPA(5,6)*SPB(5,2)+
complex<T>(-1,0)*SPA(6,7)*SPB(7,2)))/(complex<T>(2,0)*pow(SPA(1,6),2)*SPB(3,2)*SPB(4,3)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))))/(complex<T>(1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))+
((complex<T>(0,-1)*SPA(1,3)*SPA(3,4)*(-(complex<T>(-1,0)*SPA(1,2)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2)))-
complex<T>(-1,0)*SPA(1,3)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)))*(SPA(2,3)*SPB(4,2)+
(SPA(1,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4)))/SPA(1,6)))/(complex<T>(3,0)*SPA(1,6)*SPA(2,3)*pow(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3),2)*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))-
(complex<T>(-1,0)*(-((complex<T>(0,-1)*pow(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2),2)*SPB(5,2)*pow(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5),2))/(complex<T>(2,0)*pow(SPA(1,6),2)*SPB(3,2)*SPB(5,4)*(SPA(5,6)*SPB(6,5)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6))*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))))-
(complex<T>(0,-1)*((complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2))*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/pow(SPA(1,6),2)+
(complex<T>(-1,0)*SPA(5,6)*(complex<T>(-1,0)*SPA(3,6)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,6)*SPB(4,2))*SPB(5,2)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,1)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1)+
complex<T>(-1,0)*SPA(6,7)*SPB(7,6)))/pow(SPA(1,6),2)))/(complex<T>(2,0)*SPB(4,3)*(SPA(5,6)*SPB(6,5)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6))*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))))/SPB(4,2)-
(complex<T>(-1,0)*((complex<T>(0,1)*(-(complex<T>(-1,0)*SPA(3,4)*SPA(5,6)*SPB(5,2)*(complex<T>(-1,0)*SPA(1,6)*SPB(2,1)+
SPA(5,6)*SPB(5,2)+
complex<T>(-1,0)*SPA(6,7)*SPB(7,2))*(SPA(5,6)*SPB(6,2)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6)))-
complex<T>(-1,0)*SPA(3,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*(complex<T>(-1,0)*SPA(1,6)*SPB(2,1)+
SPA(5,6)*SPB(5,2)+
complex<T>(-1,0)*SPA(6,7)*SPB(7,2))*(SPA(5,6)*SPB(6,2)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6))))/(complex<T>(2,0)*pow(SPA(1,6),complex<T>(2,0))*SPA(3,4)*SPB(4,2)*(SPA(5,6)*SPB(6,2)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6))*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))-
(complex<T>(-1,0)*((complex<T>(0,-1)*SPA(5,6)*((-(complex<T>(-1,0)*SPA(5,6)*SPB(5,4)*SPB(6,2))-
SPA(5,6)*SPB(5,2)*SPB(6,4)-
complex<T>(-1,0)*SPB(6,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-
(complex<T>(-1,0)*(SPA(1,7)*SPA(5,6)*SPB(5,4)*SPB(7,2)+
SPA(1,5)*SPB(5,4)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(5,2)*SPB(7,4)+
SPA(1,7)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,4)))/SPA(1,6))*(-(complex<T>(-1,0)*SPA(3,4)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*(complex<T>(-3,0)*SPB(4,3)*SPB(5,2)+
SPB(4,2)*SPB(5,3)))-
complex<T>(-2,0)*SPA(3,4)*SPB(4,2)*SPB(4,3)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))+
SPA(3,4)*SPB(4,2)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*(complex<T>(-1,0)*SPA(5,6)*SPB(4,3)*SPB(5,2)*SPB(6,5)+
SPA(5,6)*SPB(3,2)*SPB(5,4)*SPB(6,5)+
SPB(4,3)*SPB(6,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))-
(complex<T>(-1,0)*(-(SPA(1,7)*SPA(5,6)*SPB(4,3)*SPB(5,2)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(3,2)*SPB(5,4)*SPB(7,5)-
complex<T>(-1,0)*SPA(1,7)*SPB(4,3)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,5)))/SPA(1,6))))/(complex<T>(6,0)*SPA(1,6)*pow(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3),2)*(SPA(5,6)*SPB(6,2)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6))*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))-
(complex<T>(-1,0)*((complex<T>(0,-1)*SPA(3,4)*(-((complex<T>(-2,0)*SPA(1,2)*SPB(3,2)*SPB(4,2)*(-(complex<T>(-1,0)*SPA(1,2)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2)))-
complex<T>(-1,0)*SPA(1,3)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)))*SPB(5,4)*(SPA(2,3)*SPB(4,2)+
(SPA(1,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4)))/SPA(1,6)))/SPA(1,6))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*(-(((-(complex<T>(-1,0)*SPA(1,2)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2)))-
complex<T>(-1,0)*SPA(1,3)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)))*(SPA(2,3)*pow(SPB(4,2),2)*SPB(5,3)+
complex<T>(-2,0)*SPA(2,3)*SPB(3,2)*SPB(4,2)*SPB(5,4)-
(complex<T>(-1,0)*(-(complex<T>(3,0)*SPA(1,3)*SPB(4,2)*SPB(5,4)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3)))-
complex<T>(-1,0)*SPA(1,3)*SPB(4,2)*SPB(5,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))-
complex<T>(-1,0)*SPA(1,3)*SPB(3,2)*SPB(5,4)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))))/SPA(1,6)))/SPA(1,6))+
complex<T>(-1,0)*SPB(3,2)*SPB(4,3)*(SPA(2,3)*SPB(4,2)+
(SPA(1,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4)))/SPA(1,6))*(SPA(2,3)*SPB(5,2)+
(SPA(1,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/SPA(1,6)))))/(complex<T>(6,0)*SPA(2,3)*pow(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3),2)*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))+
(complex<T>(0,1)*pow(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2),2)*(-(SPA(3,4)*pow(SPA(5,6)*SPB(6,2)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6),2)*(-(complex<T>(3,0)*SPA(4,5)*SPB(4,2)*SPB(5,4)*(SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3)))+
complex<T>(-1,0)*(-(SPA(4,5)*SPB(4,2)*SPB(5,3))-
complex<T>(-1,0)*SPB(3,2)*(SPA(3,4)*SPB(4,3)+
complex<T>(-1,0)*SPA(4,5)*SPB(5,4)))*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4)))*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2))-
SPA(4,5)*pow(SPB(3,2),2)*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))*(complex<T>(2,0)*SPA(3,4)*SPB(4,2)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(1,6)*SPB(6,1)+
SPA(2,6)*SPB(6,2)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(2,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))+
(-(complex<T>(-1,0)*SPA(3,4)*(complex<T>(-3,0)*SPB(4,3)*SPB(5,2)+
SPB(4,2)*SPB(5,3))*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4)))-
SPB(4,2)*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(4,5)*SPB(5,4))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))-
SPB(3,2)*(SPA(5,6)*SPB(6,2)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(-(complex<T>(2,0)*SPA(3,4)*SPA(4,5)*SPB(4,2)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(1,6)*SPB(6,1)+
SPA(2,6)*SPB(6,2)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(2,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))+
(-(SPA(3,4)*(-(complex<T>(-3,0)*SPA(4,5)*SPB(4,2)*SPB(5,4)*(SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3)))+
complex<T>(-1,0)*(-(complex<T>(3,0)*SPA(4,5)*SPB(4,3)*SPB(5,2))-
complex<T>(2,0)*SPA(4,5)*SPB(4,2)*SPB(5,3)-
complex<T>(-1,0)*SPB(3,2)*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(4,5)*SPB(5,4)))*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))))+
SPA(4,5)*SPB(4,2)*(SPA(3,4)*SPB(4,3)+
complex<T>(-1,0)*SPA(4,5)*SPB(5,4))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))))/(complex<T>(6,0)*pow(SPA(1,6),2)*(SPA(5,6)*SPB(6,2)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6))*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))*pow(SPA(1,6)*SPA(2,3)*SPB(3,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,6)*SPA(5,6)*SPB(6,2)*SPB(6,5)+
SPA(1,3)*SPA(1,7)*SPB(3,2)*SPB(7,1)+
SPA(1,7)*SPA(2,3)*SPB(3,2)*SPB(7,2)+
complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(6,5)*SPB(7,2)+
complex<T>(-1,0)*SPA(1,5)*SPB(6,5)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))+
complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(6,2)*SPB(7,5)-
(complex<T>(-1,0)*(-(pow(SPA(1,7),complex<T>(2,0))*SPA(5,6)*SPB(7,2)*SPB(7,5))-
SPA(1,5)*SPA(1,7)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))*SPB(7,5)))/SPA(1,6)+
SPA(1,3)*SPB(3,2)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)),2))))/SPB(3,2)))/SPB(5,4)))/SPB(4,3))/(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-
(complex<T>(-1,0)*(((-((complex<T>(0,1)*SPB(5,2)*pow(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2),2))/(complex<T>(2,0)*pow(SPA(1,6),2)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))))-
(complex<T>(0,-1)*SPB(5,2)*((pow(SPA(5,6),2)*pow(SPB(5,2),2))/pow(SPA(1,6),2)+
(complex<T>(-2,0)*SPA(5,6)*SPB(5,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2)))/pow(SPA(1,6),2)+
pow(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2),2)/pow(SPA(1,6),2)))/(complex<T>(2,0)*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))))/(SPB(3,2)*SPB(4,2))-
((-((complex<T>(0,-1)*SPA(5,6)*pow(SPA(2,3)*SPB(5,2)+
(SPA(1,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/SPA(1,6),2))/(complex<T>(2,0)*(-(complex<T>(-1,0)*SPA(1,2)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2)))-
complex<T>(-1,0)*SPA(1,3)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)))))-
(complex<T>(0,-1)*SPA(3,4)*(SPA(2,3)*SPB(4,2)+
(SPA(1,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4)))/SPA(1,6))*(SPA(2,3)*SPB(5,2)+
(SPA(1,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))/SPA(1,6)))/(complex<T>(2,0)*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))))/(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))+
(complex<T>(0,-1)*((pow(SPA(3,4),2)*SPB(3,2)*SPB(4,3)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))*(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))/pow(SPA(1,6),2)+
(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)*SPB(4,2)*SPB(5,4)*(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))*pow(-(complex<T>(-1,0)*SPA(1,7)*SPA(4,6)*SPB(7,1))-
complex<T>(-1,0)*SPA(4,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))-
complex<T>(-1,0)*SPA(4,6)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)),2))/pow(SPA(1,6),2)-
(SPA(3,4)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(4,6)*SPB(7,1))-
complex<T>(-1,0)*SPA(4,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))-
complex<T>(-1,0)*SPA(4,6)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))*(complex<T>(-2,0)*SPB(3,2)*SPB(4,2)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))*(SPA(2,6)*SPB(6,5)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(2,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6))*(SPA(2,3)*SPB(3,2)+
SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))-
SPB(4,2)*(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))*(-(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)*SPB(4,3)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5)))-
complex<T>(-3,0)*SPB(5,4)*(SPA(1,6)*SPB(3,1)+
SPA(6,7)*SPB(7,3))*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))+
complex<T>(-1,0)*(SPA(1,6)*SPB(4,1)+
SPA(6,7)*SPB(7,4))*(complex<T>(-1,0)*SPB(4,2)*SPB(5,3)*pow(SPA(2,3)*SPB(3,2)+
SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6),2)+
complex<T>(-3,0)*SPA(2,3)*SPB(3,2)*SPB(4,3)*SPB(5,2)*(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))+
complex<T>(-1,0)*SPB(3,2)*SPB(5,4)*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))*(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))))/pow(SPA(1,6),2)))/(complex<T>(6,0)*SPB(3,2)*SPB(4,3)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6))*(SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))*pow(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SPA(1,6)*SPB(6,1)+
SPA(5,6)*SPB(6,5)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(5,6)*SPB(7,5))-
complex<T>(-1,0)*SPA(1,5)*(SPA(1,6)*SPB(5,1)+
SPA(6,7)*SPB(7,5))))/SPA(1,6)+
SPA(6,7)*SPB(7,6),2)))/SPA(2,3))/(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
(complex<T>(0,-1)*pow(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2),2)*((complex<T>(-1,0)*(-((complex<T>(-9,0)*SPA(3,4)*SPA(4,5)*SPB(5,2)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4)))/(complex<T>(-1,0)*pow(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2),2)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-
complex<T>(-1,0)*SPB(3,2)*(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2))*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))))-
(complex<T>(-1,0)*((complex<T>(9,0)*SPA(4,5)*SPB(5,2))/(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2))-
(complex<T>(-1,0)*(-(complex<T>(-23,0)*SPB(5,3))-
(complex<T>(3,0)*SPA(3,5)*pow(SPB(5,3),2))/(SPA(1,6)*SPB(6,1)+
SPA(2,6)*SPB(6,2)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(2,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))))/SPB(4,3)))/SPB(3,2)))/(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-
(complex<T>(-1,0)*(-((complex<T>(-3,0)*SPB(5,3)*((complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2))*pow(SPB(5,3),2)*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))+
SPB(4,3)*(complex<T>(-1,0)*SPA(3,4)*SPB(3,2)+
SPA(4,5)*SPB(5,2))*SPB(5,4)*(-(complex<T>(-1,0)*SPA(1,6)*SPA(2,4)*SPB(6,2))-
complex<T>(-1,0)*SPA(1,4)*SPA(1,7)*SPB(7,1)-
complex<T>(-1,0)*SPA(1,7)*SPA(2,4)*SPB(7,2)-
complex<T>(-1,0)*SPA(1,4)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))))/(pow(SPB(3,2),2)*pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*pow(SPA(1,6)*SPB(6,1)+
SPA(2,6)*SPB(6,2)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(2,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6)+
SPA(6,7)*SPB(7,6),2)))+
(-((complex<T>(-1,0)*((complex<T>(3,0)*((complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2))*pow(SPB(5,3),2)*pow(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)),2)-
complex<T>(-1,0)*SPA(3,5)*SPB(4,3)*(complex<T>(-1,0)*SPA(3,4)*SPB(3,2)+
SPA(4,5)*SPB(5,2))*SPB(5,4)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(-(complex<T>(-1,0)*SPA(1,6)*SPA(2,4)*SPB(6,2))-
complex<T>(-1,0)*SPA(1,4)*SPA(1,7)*SPB(7,1)-
complex<T>(-1,0)*SPA(1,7)*SPA(2,4)*SPB(7,2)-
complex<T>(-1,0)*SPA(1,4)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))))/(pow(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*(SPA(1,6)*SPB(6,1)+
SPA(2,6)*SPB(6,2)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(2,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6)+
SPA(6,7)*SPB(7,6))*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))+
(complex<T>(-3,0)*(complex<T>(-3,0)*SPA(3,4)*pow(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2),2)*SPB(5,4)*(SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))-
complex<T>(-1,0)*SPA(4,5)*SPB(3,2)*SPB(4,2)*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)+
SPA(4,5)*SPB(5,4))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))-
SPA(3,4)*(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2))*(complex<T>(-1,0)*(-(complex<T>(3,0)*SPA(4,5)*SPB(4,2)*SPB(5,4)*(SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3)))-
SPB(3,2)*(SPA(3,4)*SPB(4,3)+
complex<T>(-1,0)*SPA(4,5)*SPB(5,4))*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4)))*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))+
complex<T>(-3,0)*SPB(3,2)*SPB(5,4)*(SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3))*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))))/((SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))*(complex<T>(-1,0)*SPA(1,6)*SPA(2,3)*SPB(3,2)*SPB(6,2)+
SPA(1,6)*SPA(3,5)*SPB(3,2)*SPB(6,5)+
SPA(1,6)*SPA(4,5)*SPB(4,2)*SPB(6,5)+
complex<T>(-1,0)*SPA(1,3)*SPA(1,7)*SPB(3,2)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,7)*SPA(2,3)*SPB(3,2)*SPB(7,2)+
SPA(1,7)*SPA(3,5)*SPB(3,2)*SPB(7,5)+
SPA(1,7)*SPA(4,5)*SPB(4,2)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,3)*SPB(3,2)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))))/SPB(3,2))+
complex<T>(3,0)*(-((complex<T>(-3,0)*pow(SPA(3,5),2)*SPB(5,3))/(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))))+
((pow(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2),2)*pow(SPB(5,3),2))/(pow(SPB(3,2),2)*(SPA(1,6)*SPB(6,1)+
SPA(2,6)*SPB(6,2)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(2,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))+
(SPA(3,4)*SPA(4,5)*SPB(4,2)*(SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*(SPA(1,6)*SPA(2,3)*SPB(3,2)*SPB(6,2)+
SPA(1,6)*SPA(3,5)*SPB(3,2)*SPB(6,5)+
SPA(1,6)*SPA(4,5)*SPB(4,2)*SPB(6,5)+
SPA(1,3)*SPA(1,7)*SPB(3,2)*SPB(7,1)+
SPA(1,7)*SPA(2,3)*SPB(3,2)*SPB(7,2)+
SPA(1,7)*SPA(3,5)*SPB(3,2)*SPB(7,5)+
SPA(1,7)*SPA(4,5)*SPB(4,2)*SPB(7,5)+
SPA(1,3)*SPB(3,2)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))*(-((complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2))*SPB(5,3)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5)))-
complex<T>(-1,0)*SPB(3,2)*(SPA(1,6)*SPA(2,3)*SPB(5,3)*SPB(6,2)+
SPA(1,3)*SPA(1,7)*SPB(5,3)*SPB(7,1)+
SPA(1,7)*SPA(2,3)*SPB(5,3)*SPB(7,2)+
SPA(1,3)*SPB(5,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6))+
complex<T>(-2,0)*(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))*(SPA(1,6)*SPB(6,1)+
SPA(2,6)*SPB(6,2)+
SPA(1,7)*SPB(7,1)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,7)*SPA(2,6)*SPB(7,2))-
complex<T>(-1,0)*SPA(1,2)*(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2))))/SPA(1,6)+
SPA(6,7)*SPB(7,6)))))/(SPB(3,2)*(-(SPA(1,6)*SPA(2,3)*SPB(6,2))-
SPA(1,3)*SPA(1,7)*SPB(7,1)-
SPA(1,7)*SPA(2,3)*SPB(7,2)-
SPA(1,3)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)))*pow(complex<T>(-1,0)*SPA(1,6)*SPA(2,3)*SPB(3,2)*SPB(6,2)+
SPA(1,6)*SPA(3,5)*SPB(3,2)*SPB(6,5)+
SPA(1,6)*SPA(4,5)*SPB(4,2)*SPB(6,5)+
complex<T>(-1,0)*SPA(1,3)*SPA(1,7)*SPB(3,2)*SPB(7,1)+
complex<T>(-1,0)*SPA(1,7)*SPA(2,3)*SPB(3,2)*SPB(7,2)+
SPA(1,7)*SPA(3,5)*SPB(3,2)*SPB(7,5)+
SPA(1,7)*SPA(4,5)*SPB(4,2)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,3)*SPB(3,2)*(SPA(1,6)*SPB(6,1)+
SPA(6,7)*SPB(7,6)),2)))/(SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5))))/(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2))))/SPB(4,3)))/(complex<T>(18,0)*pow(SPA(1,6),2)*(SPA(1,6)*SPB(6,1)+
SPA(1,7)*SPB(7,1)+
SPA(6,7)*SPB(7,6)))))/SPB(5,4)))/SPA(6,7))-
(complex<T>(0,1)*pow(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2),2)*(complex<T>(14,0)-
(complex<T>(9,0)*pow(SPA(5,6),2)*pow(complex<T>(-1,0)*SPA(1,4)*SPB(5,4)+
SPA(1,6)*SPB(6,5)+
SPA(1,7)*SPB(7,5),2)*SPB(7,6))/((-(complex<T>(-1,0)*SPA(1,2)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2)))-
complex<T>(-1,0)*SPA(1,3)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)))*(-(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,2)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,2)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,2)))-
complex<T>(-1,0)*SPA(1,3)*(complex<T>(-1,0)*SPA(4,6)*SPB(4,3)+
complex<T>(-1,0)*SPA(5,6)*SPB(5,3)))*SPB(7,6))+
(complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))))/(complex<T>(18,0)*SPA(2,3)*SPB(5,4)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*SPB(7,6)*SS(1,2,3))-
(complex<T>(0,-1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))*((complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,3))*pow(-(complex<T>(-1,0)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2))/(SPB(3,2)*SPB(4,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))+
((complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))*((complex<T>(3,0)*SPA(1,4)*SPA(2,4)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,3))*pow(-(complex<T>(-1,0)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2))/(pow(SPA(1,2),2)*pow(SPB(3,2),2)*pow(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)),2))+
(complex<T>(3,0)*(-((complex<T>(3,0)*SPA(3,4)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,2))*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*SPB(7,1))-
SPA(2,3)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*SPA(2,3)*SPB(4,2)*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(7,4))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6)))/(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(SPA(2,3)*SPB(3,2)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))+
((-(complex<T>(-1,0)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4))*(complex<T>(3,0)*SPA(1,2)*SPA(3,4)*(SPA(1,5)*SPB(5,4)+
SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*SPB(7,1))-
complex<T>(-1,0)*SPA(2,4)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))*SPB(7,3)-
complex<T>(-1,0)*SPA(2,4)*SPB(3,2)*SPB(7,4))+
(-(SPA(2,3)*pow(SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3),2))-
complex<T>(-1,0)*SPA(1,4)*SPA(2,4)*SPB(4,3)*(SPA(1,5)*SPB(5,4)+
SPA(1,6)*SPB(6,4))-
complex<T>(2,0)*SPA(1,7)*SPA(2,3)*(SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3))*SPB(7,3)-
pow(SPA(1,7),complex<T>(2,0))*SPA(2,3)*pow(SPB(7,3),2)-
complex<T>(-1,0)*SPA(1,4)*SPA(1,7)*SPA(2,4)*SPB(4,3)*SPB(7,4))*(complex<T>(-1,0)*SPA(3,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(3,6)*SPB(7,6))-
complex<T>(-3,0)*SPA(1,2)*SPA(2,3)*SPB(4,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,2))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/(SPA(1,2)*SPA(2,3)*SPB(3,2)*SPB(4,3))))/((complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))*(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))-
(complex<T>(-1,0)*(-((complex<T>(-3,0)*SPA(2,3)*pow(complex<T>(-1,0)*SPA(1,5)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,3),3)*pow(-(complex<T>(-1,0)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2))/(pow(SPA(1,2),2)*pow(SPB(3,2),2)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*pow(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)),2)))-
(complex<T>(-1,0)*(-((complex<T>(-9,0)*SPA(3,4)*SPB(4,2)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*SPB(7,1))-
complex<T>(-1,0)*SPA(2,4)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(3,4)*SPB(4,3))*SPB(7,3)-
complex<T>(-1,0)*SPA(2,4)*SPB(3,2)*SPB(7,4))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6)))/(SPA(2,3)*(-(SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SPA(1,7)*SPB(7,1))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))))-
((complex<T>(-1,0)*SPA(1,5)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,3))*pow(-(complex<T>(-1,0)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2)*(complex<T>(-23,0)*SPA(1,2)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4))+
complex<T>(-23,0)*SPA(1,2)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(-(complex<T>(-3,0)*SPA(2,3)*(SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3)))-
complex<T>(-1,0)*SPA(1,7)*(complex<T>(23,0)*SPA(1,2)*SPB(7,1)+
complex<T>(3,0)*SPA(2,3)*SPB(7,3)))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/(SPA(1,2)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*pow(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)),2))))/(SPB(3,2)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))-
(complex<T>(3,0)*(complex<T>(3,0)*(SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3))*pow(complex<T>(-1,0)*SPA(3,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(3,6)*SPB(7,6),2)-
(complex<T>(-1,0)*SPA(3,4)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*SPB(7,1))-
SPA(2,3)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*SPA(2,3)*SPB(4,2)*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(7,4))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))*(complex<T>(-1,0)*SPA(3,4)*SPB(4,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))+
(SPA(1,5)*SPB(5,4)+
SPA(1,6)*SPB(6,4)+
SPA(1,7)*SPB(7,4))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(SPA(2,3)*SPB(3,2)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))))/(SPA(2,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*pow(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6),2))+
(complex<T>(3,0)*(-((pow(SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3),2)*pow(-(complex<T>(-1,0)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2))/(SPA(1,2)*pow(SPB(3,2),2)*(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))))+
(SPA(3,4)*SPB(4,2)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*SPB(7,1))-
SPA(2,3)*SPB(4,3)*SPB(7,2)-
complex<T>(-1,0)*SPA(2,3)*SPB(4,2)*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(4,2)+
SPA(3,4)*SPB(4,3))*SPB(7,4))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))*(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))*(-(SPA(1,5)*(complex<T>(-2,0)*SPA(1,2)*SPB(3,2)+
SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3))*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(complex<T>(-2,0)*SPA(1,2)*SPB(3,2)+
SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3)+
SPA(1,7)*SPB(7,3))*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(complex<T>(-1,0)*(SPA(1,5)*SPB(5,3)+
SPA(1,6)*SPB(6,3))*(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1))-
complex<T>(-1,0)*SPA(1,7)*(-(complex<T>(-1,0)*SPA(1,7)*SPB(7,1)*SPB(7,3))-
complex<T>(-1,0)*SPB(3,2)*(complex<T>(-2,0)*SPA(1,2)*SPB(7,1)+
SPA(2,3)*SPB(7,3))))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))))/(SPA(2,3)*SPB(3,2)*(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))*pow(-(SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(complex<T>(-1,0)*SPA(2,3)*SPB(3,2)+
SPA(1,7)*SPB(7,1))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)),2))-
(SPB(4,2)*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))*(SPA(3,4)*SPB(4,3)*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6))+
(complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*(complex<T>(-1,0)*SPA(4,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(4,6)*SPB(7,6))))/(-(complex<T>(-1,0)*SPA(1,5)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(5,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(5,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(5,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(5,2)+
SPA(3,4)*SPB(5,3))*SPB(7,4)))-
complex<T>(-1,0)*SPA(1,6)*(-(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(6,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(6,4))*SPB(7,1))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,3)*SPB(6,3)+
complex<T>(-1,0)*SPA(2,4)*SPB(6,4))*SPB(7,2)-
complex<T>(-1,0)*(SPA(2,3)*SPB(6,2)+
complex<T>(-1,0)*SPA(3,4)*SPB(6,4))*SPB(7,3)-
complex<T>(-1,0)*(SPA(2,4)*SPB(6,2)+
SPA(3,4)*SPB(6,3))*SPB(7,4))+
(SPA(2,3)*SPB(3,2)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,1))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))))/((complex<T>(-1,0)*SPA(1,5)*SPB(5,4)+
complex<T>(-1,0)*SPA(1,6)*SPB(6,4)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,4))*(complex<T>(-1,0)*SPA(1,5)*SPB(7,5)+
complex<T>(-1,0)*SPA(1,6)*SPB(7,6)))))/SPB(4,3)))/complex<T>(9,0)))/(complex<T>(2,0)*(complex<T>(-1,0)*SPA(1,6)*SPB(6,5)+
complex<T>(-1,0)*SPA(1,7)*SPB(7,5))*SPB(7,6)*SSS(1,2,3,4)))+
(complex<T>(0,-1)*(-((complex<T>(-1,0)*pow(SPA(4,5),2)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*pow(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2),2)*pow(SPB(7,5),2)*SPB(7,6))/((SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))*((complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))-
complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(1,5)*SPA(4,6)*SPB(5,4))-
SPA(1,4)*SPA(5,6)*SPB(5,4)-
complex<T>(-1,0)*SPA(1,6)*(SPA(4,6)*SPB(6,4)+
SPA(5,6)*SPB(6,5))-
complex<T>(-1,0)*SPA(1,7)*(SPA(4,6)*SPB(7,4)+
SPA(5,6)*SPB(7,5)))*SPB(7,6))*SS(1,2,3)))+
(pow(SPA(3,4),2)*(-((complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*pow(SPA(1,5)*SPB(7,1)+
SPA(2,5)*SPB(7,2)+
SPA(3,5)*SPB(7,3)+
SPA(4,5)*SPB(7,4),2)*pow(SPB(7,5),2))/(SPA(5,6)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*(SPA(1,4)*SPB(7,1)+
SPA(2,4)*SPB(7,2)+
SPA(3,4)*SPB(7,3))))-
complex<T>(-1,0)*pow(SPB(7,4),2)*SPB(7,6))-
(complex<T>(-1,0)*pow(SPA(3,4),2)*pow(SPB(7,4),2)*SPB(7,6)*SS(1,2,3))/(complex<T>(-1,0)*SS(1,2,3)+
SSS(1,2,3,4))-
(complex<T>(-1,0)*pow(-(complex<T>(-1,0)*SPA(5,6)*(SPA(1,3)*SPB(6,1)+
SPA(2,3)*SPB(6,2)))-
complex<T>(-1,0)*SPA(5,7)*(SPA(1,3)*SPB(7,1)+
SPA(2,3)*SPB(7,2)),2)*pow(SPB(7,5),2)*SPB(7,6))/(SS(1,2,3)*(SPA(6,7)*SPB(7,6)+
complex<T>(-1,0)*SSS(1,2,3,4))))/SSS(1,2,3,4)))/(complex<T>(2,0)*SPA(2,3)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*SPB(5,4)*pow(SPB(7,6),2))-
(complex<T>(0,1)*(-((pow(SPB(4,2),2)*(-(complex<T>(-1,0)*(SPA(1,2)*SPB(4,2)+
SPA(1,3)*SPB(4,3))*SPB(5,3))-
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(1,2)*SPB(3,2)+
SPA(1,4)*SPB(4,3))*SPB(5,4))*pow(SPB(7,2),2))/(SPB(4,3)*(SPA(1,2)*SPB(4,2)+
SPA(1,3)*SPB(4,3))*SPB(5,4)*(SPA(1,2)*SPB(5,2)+
SPA(1,3)*SPB(5,3)+
SPA(1,4)*SPB(5,4))*SPB(7,6)))+
(SPB(4,2)*(complex<T>(-1,0)*SPA(3,4)*SPB(3,2)*SPB(4,3)*SPB(5,4)+
SPA(3,4)*SPB(4,3)*(complex<T>(-1,0)*SPB(4,2)*SPB(5,3)+
SPB(3,2)*SPB(5,4)))*pow(SPA(1,6)*SPB(2,1)+
SPA(6,7)*SPB(7,2),2))/(SPA(3,4)*SPA(6,7)*pow(SPB(4,3),2)*SPB(5,4)*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SS(1,6,7))+
(SPB(4,2)*(SPA(3,4)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(4,3)+
SPA(3,4)*SPB(3,2)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3)))*pow(-(complex<T>(-1,0)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(7,1))-
complex<T>(-1,0)*(SPA(2,3)*SPB(3,2)+
SPA(2,4)*SPB(4,2))*SPB(7,2)-
complex<T>(-1,0)*SPA(3,4)*SPB(4,2)*SPB(7,3)-
SPA(3,4)*SPB(3,2)*SPB(7,4),2))/(SPA(3,4)*SPB(4,3)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))*SPB(7,6)*SSS(1,2,3,4)*SSS(1,5,6,7))+
(complex<T>(-1,0)*((complex<T>(-1,0)*pow(SPA(1,6),2)*pow(SPB(2,1),2)*SPB(4,2)*(complex<T>(-1,0)*SPA(3,4)*SPB(3,2)*SPB(4,3)*SPB(5,4)+
SPA(3,4)*SPB(4,3)*(complex<T>(-1,0)*SPB(4,2)*SPB(5,3)+
SPB(3,2)*SPB(5,4))))/(SPA(3,4)*SPA(6,7)*pow(SPB(4,3),2)*SPB(5,4)*(SPA(6,7)*SPB(7,6)+
complex<T>(-1,0)*SS(1,6,7)))-
(SPA(3,4)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*(-(complex<T>(1,0)/SPA(3,4))-
(SPB(3,2)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3)))/(SPA(3,4)*(SPA(1,3)*SPB(3,2)+
SPA(1,4)*SPB(4,2))*SPB(4,3)))*((pow(SPB(4,2),2)*pow(complex<T>(-1,0)*SPA(3,5)*SPB(3,2)+
complex<T>(-1,0)*SPA(4,5)*SPB(4,2),2)*pow(SPB(7,5),2)*SSS(1,2,3,4))/(complex<T>(-1,0)*SPA(6,7)*SPB(7,6)+
SSS(1,2,3,4))+
(pow(SPB(2,1),2)*pow(SPB(4,2),2)*pow(complex<T>(-1,0)*SPA(1,2)*SPB(7,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(7,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(7,4),2)*SSS(1,5,6,7))/(complex<T>(-1,0)*SSS(1,2,3,4)+
SSS(1,5,6,7))))/(SPB(4,2)*(complex<T>(-1,0)*SPA(1,2)*SPB(4,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(4,3))*SPB(7,6)*SSS(1,2,3,4)*SSS(1,5,6,7))))/(complex<T>(-1,0)*SPA(1,2)*SPB(5,2)+
complex<T>(-1,0)*SPA(1,3)*SPB(5,3)+
complex<T>(-1,0)*SPA(1,4)*SPB(5,4))))/(complex<T>(2,0)*SPB(3,2)*pow(SPB(4,2),2))

 )
;

#endif

}


template <class T> complex<T> R2q2Q1g2l_qmmQbmQpqbplmlbp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
//{{qm, m, Qbm, Qp, qbp, lm, lbp}, lc}

#if _VERBOSE
  _MESSAGE("R2q2Q1g2l :  qmmQbmQpqbplmlbp lc");
#endif

complex<T> recursive;
recursive=complex<T>(1,0);

#if _USE_OPTIMIZED

  complex<T> t1;
  complex<T> t10;
  complex<T> t1001;
  complex<T> t1007;
  complex<T> t1012;
  complex<T> t1024;
  complex<T> t1025;
  complex<T> t1029;
  complex<T> t103;
  complex<T> t1035;
  complex<T> t1036;
  complex<T> t1039;
  complex<T> t104;
  complex<T> t1045;
  complex<T> t1047;
  complex<T> t105;
  complex<T> t106;
  complex<T> t1060;
  complex<T> t1068;
  complex<T> t1071;
  complex<T> t1076;
  complex<T> t1098;
  complex<T> t1109;
  complex<T> t111;
  complex<T> t1115;
  complex<T> t112;
  complex<T> t113;
  complex<T> t1135;
  complex<T> t114;
  complex<T> t1156;
  complex<T> t1157;
  complex<T> t1163;
  complex<T> t1164;
  complex<T> t1171;
  complex<T> t1172;
  complex<T> t1197;
  complex<T> t12;
  complex<T> t120;
  complex<T> t1202;
  complex<T> t1205;
  complex<T> t1209;
  complex<T> t122;
  complex<T> t1223;
  complex<T> t1227;
  complex<T> t1230;
  complex<T> t124;
  complex<T> t1240;
  complex<T> t125;
  complex<T> t1267;
  complex<T> t1271;
  complex<T> t1272;
  complex<T> t1274;
  complex<T> t128;
  complex<T> t1287;
  complex<T> t129;
  complex<T> t1290;
  complex<T> t13;
  complex<T> t130;
  complex<T> t1319;
  complex<T> t132;
  complex<T> t136;
  complex<T> t137;
  complex<T> t140;
  complex<T> t143;
  complex<T> t144;
  complex<T> t147;
  complex<T> t15;
  complex<T> t155;
  complex<T> t157;
  complex<T> t158;
  complex<T> t16;
  complex<T> t163;
  complex<T> t17;
  complex<T> t176;
  complex<T> t179;
  complex<T> t181;
  complex<T> t182;
  complex<T> t183;
  complex<T> t184;
  complex<T> t187;
  complex<T> t19;
  complex<T> t193;
  complex<T> t194;
  complex<T> t195;
  complex<T> t199;
  complex<T> t2;
  complex<T> t20;
  complex<T> t208;
  complex<T> t21;
  complex<T> t210;
  complex<T> t215;
  complex<T> t219;
  complex<T> t22;
  complex<T> t228;
  complex<T> t23;
  complex<T> t230;
  complex<T> t248;
  complex<T> t25;
  complex<T> t250;
  complex<T> t253;
  complex<T> t257;
  complex<T> t258;
  complex<T> t259;
  complex<T> t26;
  complex<T> t260;
  complex<T> t261;
  complex<T> t262;
  complex<T> t263;
  complex<T> t266;
  complex<T> t272;
  complex<T> t275;
  complex<T> t278;
  complex<T> t28;
  complex<T> t281;
  complex<T> t285;
  complex<T> t286;
  complex<T> t287;
  complex<T> t29;
  complex<T> t292;
  complex<T> t294;
  complex<T> t295;
  complex<T> t299;
  complex<T> t3;
  complex<T> t300;
  complex<T> t301;
  complex<T> t302;
  complex<T> t304;
  complex<T> t305;
  complex<T> t307;
  complex<T> t308;
  complex<T> t31;
  complex<T> t314;
  complex<T> t317;
  complex<T> t318;
  complex<T> t327;
  complex<T> t328;
  complex<T> t340;
  complex<T> t342;
  complex<T> t35;
  complex<T> t352;
  complex<T> t356;
  complex<T> t359;
  complex<T> t36;
  complex<T> t360;
  complex<T> t364;
  complex<T> t367;
  complex<T> t38;
  complex<T> t383;
  complex<T> t387;
  complex<T> t388;
  complex<T> t39;
  complex<T> t391;
  complex<T> t392;
  complex<T> t397;
  complex<T> t398;
  complex<T> t4;
  complex<T> t403;
  complex<T> t404;
  complex<T> t405;
  complex<T> t411;
  complex<T> t412;
  complex<T> t415;
  complex<T> t421;
  complex<T> t423;
  complex<T> t426;
  complex<T> t430;
  complex<T> t433;
  complex<T> t436;
  complex<T> t453;
  complex<T> t455;
  complex<T> t456;
  complex<T> t457;
  complex<T> t46;
  complex<T> t460;
  complex<T> t461;
  complex<T> t465;
  complex<T> t466;
  complex<T> t467;
  complex<T> t468;
  complex<T> t469;
  complex<T> t472;
  complex<T> t477;
  complex<T> t478;
  complex<T> t48;
  complex<T> t480;
  complex<T> t485;
  complex<T> t488;
  complex<T> t49;
  complex<T> t491;
  complex<T> t50;
  complex<T> t502;
  complex<T> t508;
  complex<T> t535;
  complex<T> t537;
  complex<T> t541;
  complex<T> t543;
  complex<T> t546;
  complex<T> t55;
  complex<T> t555;
  complex<T> t557;
  complex<T> t558;
  complex<T> t56;
  complex<T> t560;
  complex<T> t578;
  complex<T> t579;
  complex<T> t58;
  complex<T> t59;
  complex<T> t594;
  complex<T> t596;
  complex<T> t6;
  complex<T> t60;
  complex<T> t605;
  complex<T> t609;
  complex<T> t610;
  complex<T> t617;
  complex<T> t62;
  complex<T> t620;
  complex<T> t621;
  complex<T> t63;
  complex<T> t634;
  complex<T> t635;
  complex<T> t638;
  complex<T> t640;
  complex<T> t649;
  complex<T> t650;
  complex<T> t653;
  complex<T> t66;
  complex<T> t67;
  complex<T> t670;
  complex<T> t68;
  complex<T> t681;
  complex<T> t682;
  complex<T> t69;
  complex<T> t692;
  complex<T> t695;
  complex<T> t698;
  complex<T> t699;
  complex<T> t7;
  complex<T> t710;
  complex<T> t713;
  complex<T> t714;
  complex<T> t72;
  complex<T> t722;
  complex<T> t725;
  complex<T> t73;
  complex<T> t730;
  complex<T> t739;
  complex<T> t74;
  complex<T> t746;
  complex<T> t75;
  complex<T> t77;
  complex<T> t772;
  complex<T> t775;
  complex<T> t777;
  complex<T> t780;
  complex<T> t784;
  complex<T> t79;
  complex<T> t791;
  complex<T> t797;
  complex<T> t80;
  complex<T> t81;
  complex<T> t821;
  complex<T> t836;
  complex<T> t837;
  complex<T> t84;
  complex<T> t844;
  complex<T> t848;
  complex<T> t849;
  complex<T> t851;
  complex<T> t853;
  complex<T> t86;
  complex<T> t860;
  complex<T> t863;
  complex<T> t864;
  complex<T> t872;
  complex<T> t873;
  complex<T> t880;
  complex<T> t885;
  complex<T> t889;
  complex<T> t89;
  complex<T> t891;
  complex<T> t897;
  complex<T> t898;
  complex<T> t899;
  complex<T> t9;
  complex<T> t90;
  complex<T> t91;
  complex<T> t910;
  complex<T> t918;
  complex<T> t92;
  complex<T> t921;
  complex<T> t922;
  complex<T> t93;
  complex<T> t930;
  complex<T> t94;
  complex<T> t944;
  complex<T> t945;
  complex<T> t948;
  complex<T> t949;
  complex<T> t95;
  complex<T> t950;
  complex<T> t956;
  complex<T> t957;
  complex<T> t959;
  complex<T> t96;
  complex<T> t960;
  complex<T> t961;
  complex<T> t965;
  complex<T> t976;
  complex<T> t979;
  complex<T> t98;
  complex<T> t984;
  complex<T> t988;
  complex<T> t99;
  complex<T> t995;
  complex<T> t996;
  {
    t1 = complex<T>(0,1);
    t2 = SPB(3,4);
    t3 = SPB(4,7);
    t4 = pow(t3,2);
    t6 = SPA(5,3);
    t7 = SPB(2,3);
    t9 = SPA(5,4);
    t10 = SPB(2,4);
    t12 = t6*t7+t9*t10;
    t13 = complex<T>(1,0)/t12;
    t15 = SPB(3,7);
    t16 = t6*t15;
    t17 = t9*t3;
    t19 = pow(t16+t17,2);
    t20 = pow(t2,2);
    t21 = complex<T>(-1,0);
    t22 = SPA(4,3);
    t23 = t21*t22;
    t25 = t21*t6;
    t26 = SPB(5,7);
    t28 = t23*t3+t25*t26;
    t29 = pow(t28,2);
    t31 = complex<T>(1,0)/t22;
    t35 = SPB(4,5);
    t36 = pow(t35,2);
    t38 = t23*t2;
    t39 = SS(3,4,5);
    t46 = t21*t9;
    t48 = t25*t7+t46*t10;
    t49 = complex<T>(1,0)/t48;
    t50 = complex<T>(1,0)/t39;
    t55 = complex<T>(2,0);
    t56 = complex<T>(1,0)/t55;
    t58 = SPB(1,2);
    t59 = complex<T>(1,0)/t58;
    t60 = complex<T>(1,0)/t20;
    t62 = SPB(6,7);
    t63 = complex<T>(1,0)/t62;
    t66 = SPA(6,5);
    t67 = pow(t66,2);
    t68 = t21*t67;
    t69 = SPA(5,2);
    t72 = t69*t7+t46*t2;
    t73 = SPA(6,2);
    t74 = t73*t10;
    t75 = SPA(6,3);
    t77 = t74+t75*t2;
    t79 = SPA(6,1);
    t80 = t21*t79;
    t81 = SPB(1,4);
    t84 = SPA(7,6);
    t86 = t80*t81+t66*t35-t84*t3;
    t89 = complex<T>(1,0)/t67;
    t90 = SPB(1,6);
    t91 = t66*t90;
    t92 = SPA(7,5);
    t93 = SPB(1,7);
    t94 = t92*t93;
    t95 = t91+t94;
    t96 = complex<T>(1,0)/t95;
    t98 = complex<T>(1,0)/t7;
    t99 = complex<T>(1,0)/t2;
    t103 = complex<T>(0,-1);
    t104 = SPA(3,2);
    t105 = t103*t104;
    t106 = t105*t10;
    t111 = t6*t2;
    t112 = t69*t10+t111;
    t113 = complex<T>(3,0);
    t114 = complex<T>(1,0)/t113;
    t120 = t21*t66;
    t122 = t21*t84;
    t124 = t120*t35-t122*t3;
    t125 = complex<T>(1,0)/t66;
    t128 = t21*t69;
    t129 = t79*t81;
    t130 = SPB(2,6);
    t132 = SPB(4,6);
    t136 = t92*t81;
    t137 = SPB(2,7);
    t140 = t79*t92;
    t143 = SPA(5,1);
    t144 = t21*t143;
    t147 = t21*t92;
    t155 = SPB(3,6);
    t157 = SPB(1,3);
    t158 = t157*t132;
    t163 = t157*t3;
    t176 = complex<T>(1,0)/t112;
    t179 = pow(t46*t124*t125+t21*(t128*(t129*t130+t80*t58*t132+t130*t124+t21*(
t80*t136*t137+t140*t58*t3+t144*t58*t124+t147*t137*t124)*t125)+t25*(t129*t155+
t80*t158+t155*t124+t21*(t80*t136*t15+t140*t163+t144*t157*t124+t147*t15*t124)*
t125))*t176,3);
    t181 = complex<T>(1,0);
    t182 = complex<T>(1,0)/t181;
    t183 = pow(t95,2);
    t184 = complex<T>(1,0)/t183;
    t187 = t92*t3;
    t193 = t79*t132+t21*(t80*t187+t143*t124)*t125;
    t194 = complex<T>(1,0)/t193;
    t195 = t9*t124;
    t199 = t79*t58;
    t208 = t143*t58;
    t210 = t92*t137;
    t215 = t80*t81*t130+t199*t132+t21*t130*t124+t21*(t140*t81*t137+t80*t92*t58*
t3+t208*t124+t210*t124)*t125;
    t219 = t79*t157;
    t228 = t143*t157;
    t230 = t92*t15;
    t248 = t79*t90;
    t250 = SPB(1,5);
    t253 = t120*t250-t122*t93;
    t257 = t21*(t80*t94+t143*t253)*t125;
    t258 = SPB(5,6);
    t259 = t66*t258;
    t260 = t92*t26;
    t261 = t84*t62;
    t262 = t248+t257+t259+t260+t261;
    t263 = complex<T>(1,0)/t262;
    t266 = t1*t104;
    t272 = t6*(t219+t73*t7)+t9*(t129+t74);
    t275 = SPB(2,5);
    t278 = t120*t275-t122*t137;
    t281 = t22*t10+t25*t278*t125;
    t285 = pow(t48,2);
    t286 = complex<T>(1,0)/t285;
    t287 = t125*t286;
    t292 = t259+t261;
    t294 = t23*t66*t132+t23*t187+t25*t260+t25*t292;
    t295 = complex<T>(1,0)/t294;
    t299 = t1*t81;
    t300 = pow(t253,2);
    t301 = SPA(2,1);
    t302 = t21*t301;
    t304 = SPA(3,1);
    t305 = t21*t304;
    t307 = t302*t10+t305*t2;
    t308 = pow(t307,2);
    t314 = complex<T>(1,0)/(t248+t257);
    t317 = complex<T>(1,0)/(t259+t260+t261);
    t318 = t99*t317;
    t327 = t147*t26;
    t328 = t122*t62;
    t340 = complex<T>(1,0)/t10;
    t342 = t21*t104;
    t352 = complex<T>(1,0)/t104;
    t356 = t194*t263;
    t359 = t103*t79;
    t360 = t104*t10;
    t364 = t90*t2;
    t367 = t7*t124;
    t383 = t55*t104;
    t387 = t104*t48;
    t388 = complex<T>(-3,0);
    t391 = t157*t10;
    t392 = t388*t81*t7+t391;
    t397 = complex<T>(6,0);
    t398 = complex<T>(1,0)/t397;
    t403 = complex<T>(-2,0);
    t404 = t403*t9;
    t405 = t58*t10;
    t411 = t21*t48;
    t412 = t22*t81;
    t415 = t412+t25*t253*t125;
    t421 = pow(t10,2);
    t423 = t403*t22;
    t426 = t6*t157;
    t430 = t278*t2;
    t433 = SPB(3,5);
    t436 = t120*t433-t122*t15;
    t453 = pow(t124,2);
    t455 = t66*t130;
    t456 = t455+t210;
    t457 = t21*t456;
    t460 = t104*t7;
    t461 = t302*t58+t460;
    t465 = t388*t301;
    t466 = t66*t155;
    t467 = t466+t230;
    t468 = t405*t467;
    t469 = t465*t468;
    t472 = pow(t193,2);
    t477 = t10*t456;
    t478 = SPA(6,4);
    t480 = t21*t478;
    t485 = t478*t132+t21*(t480*t187+t195)*t125+t259+t260+t261;
    t488 = t21*t95;
    t491 = t301*t58+t342*t7;
    t502 = t403*t301;
    t508 = t301*t95;
    t535 = t22*t66;
    t537 = t535*t2*t132;
    t541 = t22*t92;
    t543 = t541*t2*t3;
    t546 = pow(t92,2);
    t555 = t6*t92;
    t557 = t555*t2*t26;
    t558 = t111*t292;
    t560 = pow(t80*t91*t132+t80*t94*t132+t537+t80*t92*t90*t3+t543+t143*t90*t124
+t21*(t79*t546*t93*t3+t144*t94*t124)*t125+t557+t558,2);
    t578 = pow(t79,2);
    t579 = pow(t81,2);
    t594 = pow(t415,2);
    t596 = complex<T>(1,0)/t272;
    t605 = pow(t104,2);
    t609 = t2*t262;
    t610 = t248+t257+t38+t259+t260+t261;
    t617 = t21*t73;
    t620 = t301*t253+t617*t260+t617*t292;
    t621 = pow(t620,2);
    t634 = t22*t2;
    t635 = t248+t257+t634+t259+t260+t261;
    t638 = t21*t10;
    t640 = t7*t2;
    t649 = t21*t157;
    t650 = pow(t635,2);
    t653 = t388*t22;
    t670 = pow(t610,2);
    t681 = complex<T>(9,0);
    t682 = t681*t301;
    t692 = complex<T>(1,0)/t307;
    t695 = complex<T>(23,0);
    t698 = pow(t157,2);
    t699 = complex<T>(1,0)/t485;
    t710 = t58*t7;
    t713 = t301*t81+t342*t2;
    t714 = SPA(4,2);
    t722 = t714*t66*t132+t714*t92*t3+t69*t92*t26+t69*t292;
    t725 = t698*t307;
    t730 = pow(t485,2);
    t739 = pow(t294,2);
    t746 = t388*t104;
    t772 = t301*t66*t90*t10;
    t775 = t301*t92*t93*t10;
    t777 = t304*t66*t364;
    t780 = t304*t92*t93*t2;
    t784 = t92*t2;
    t791 = t772+t775+t777+t780+t23*t66*t2*t132+t23*t784*t3+t25*t784*t26+t25*t2*
t292;
    t797 = pow(t304,2);
    t821 = pow(t791,2);
    t836 = complex<T>(18,0);
    t837 = complex<T>(1,0)/t836;
    t844 = complex<T>(1,0)/t84;
    t848 = complex<T>(14,0);
    t849 = complex<T>(-9,0);
    t851 = t128*t58;
    t853 = pow(t851+t91+t94,2);
    t860 = t144*t58+t120*t130+t147*t137;
    t863 = t301*t93-t73*t62;
    t864 = t860*t863;
    t872 = t31*t59;
    t873 = complex<T>(1,0)/t860;
    t880 = t143*t93-t66*t62;
    t885 = t144*t157+t120*t155+t147*t15;
    t889 = t10*t15;
    t891 = t714*t10;
    t897 = t104*t137*t2+t342*t889+t21*(t891+t634)*t3+t21*t112*t26;
    t898 = pow(t897,2);
    t899 = t21*t885*t898;
    t910 = t21*t714;
    t918 = t851+t25*t157+t46*t81;
    t921 = t21*(t104*t157+t714*t81)*t137+t21*(t342*t58+t412)*t15+t21*(t910*t58+
t23*t157)*t3+t21*t918*t26;
    t922 = t143*t921;
    t930 = t23*t132;
    t944 = t21*t137*(t342*t155+t910*t132)+t21*t15*(t104*t130+t930)+t21*(t714*
t130+t22*t155)*t3+t21*(t69*t130+t6*t155+t9*t132)*t26;
    t945 = t120*t944;
    t948 = t922+t945+t147*t26*t880;
    t949 = complex<T>(1,0)/t948;
    t950 = t99*t949;
    t956 = pow(t9,2);
    t957 = complex<T>(1,0)/t956;
    t959 = pow(t948,2);
    t960 = complex<T>(1,0)/t959;
    t961 = t60*t960;
    t965 = t208+t455+t210;
    t976 = t910*t137*t2+t21*(t460+t634)*t15+t910*t7*t3+t21*t72*t26;
    t979 = t113*t22;
    t984 = t144*t81+t120*t132+t147*t3;
    t988 = t714*t69;
    t995 = t228+t466;
    t996 = pow(t995,2);
    t1001 = pow(t15,2);
    t1007 = t304*t93-t75*t62;
    t1012 = complex<T>(1,0)/t9;
    t1024 = t21*(t460+t891)*t137+t23*t889+t22*t7*t3+t411*t26;
    t1025 = t1024*t863;
    t1029 = complex<T>(1,0)/(t922+t945+(t634+t327)*t880);
    t1035 = complex<T>(1,0)/t880;
    t1036 = t873*t1035;
    t1039 = pow(t885,3);
    t1045 = t228+t466+t230;
    t1047 = pow(t1007,2);
    t1060 = pow(t880,2);
    t1068 = t66*t944;
    t1071 = t144*t921+t1068+(t38+t260)*t880;
    t1076 = complex<T>(-23,0);
    t1098 = pow(t1045,2);
    t1109 = t38+t327;
    t1115 = t228+t404*t2+t466+t230;
    t1135 = pow(t1071,2);
    t1156 = SSS(2,3,4,5);
    t1157 = complex<T>(1,0)/t1156;
    t1163 = pow(t301,2);
    t1164 = pow(t93,2);
    t1171 = t342*t15+t910*t3+t128*t26;
    t1172 = complex<T>(1,0)/t1171;
    t1197 = SPA(4,1);
    t1202 = pow(t302*t137+t305*t15+t21*t1197*t3+t144*t26,2);
    t1205 = complex<T>(1,0)/t918;
    t1209 = pow(t137,2);
    t1223 = SPA(7,1);
    t1227 = pow(t80*(t930+t25*t258)+t21*t1223*t28,2);
    t1230 = t21*t1156;
    t1240 = pow(t62,2);
    t1267 = t10*(t342*t710*t2+t460*(t649*t10+t58*t2));
    t1271 = pow(t7,2);
    t1272 = complex<T>(1,0)/t1271;
    t1274 = SS(5,6,7);
    t1287 = SSS(1,5,6,7);
    t1290 = t63/t1287*t1157;
    t1319 = pow(t69*t137+t16+t17,2);
    return(t1*(t2*t4*t13+t19*(t20*t29*t31/t19+t2*t36/(t38+t39))*t49*t50)*t56*
t59*t60*t63-recursive*(-t68*(t1*(-t72*t77*t86*t56*t89*t96*t98*t99+t1*(t106*t2*
t56*t98+t106*t112*t114*t49)*t179*t182*t184*t194/(t195*t125+t21*(t128*t215+t25*(
t80*t81*t155+t219*t132+t21*t155*t124+t21*(t140*t81*t15+t80*t92*t157*t3+t228*
t124+t230*t124)*t125))*t176))*t182*t49*t263+(t266*t6*t272*t281*t114*t31*t287*
t295+(t299*t300*t308*t56*t89*t59*t314*t318*t263+t1*(t21*t253*t307*t124*t89+t129
*t77*(t120*t258+t327+t328)*t89)*t56*t314*t98*t263)*t340+(t1*(t342*t79*t81*t86*
t193+t104*t86*t124*t193)*t56*t352*t89*t340*t356+(t359*(t360*t48*(t80*t81*t90*t7
+t199*t364+t21*t90*t367+t21*(t140*t81*t93*t7+t80*t92*t58*t93*t2+t94*t367)*t125)
+(t383*t95*t7*t10+t387*t392)*t215)*t398*t287*t356+(t105*(t404*t405*t272*t281*t2
*t125+t411*(t21*t415*t7*t281*t2+t272*(t22*t157*t421+t423*t405*t2+t21*(t426*t10*
t278+t6*t58*t430+t388*t6*t405*t436)*t125)*t125))*t398*t31*t286*t295+t1*t453*(
t342*t183*(t457*(t302*t391+t461*t2)+t469)*t472+t302*t20*t294*(t383*t95*t477*
t485+(t488*t491*t10+t104*t392*t456)*t294)+t488*t2*t193*(t502*t104*t95*t10*t456*
t485+(t508*t461*t10+t342*(t457*(t465*t81*t7+t502*t391+t491*t2)+t113*t301*t468))
*t294))*t398*t89*t194*t317*t295/t560)*t99)*t59)*t98)*t96+(((t103*t81*t453*t56*
t89*t317+t299*(t578*t579*t89+t55*t79*t81*t124*t89+t453*t89)*t56*t263)*t340*t99+
t21*((t359*t594*t56*t596+t266*t415*t281*t56*t295)*t49+t103*(t605*t253*t7*t278*
t609*t610*t89+t23*t405*t2*t610*t621*t89+t342*t620*(t403*(t478*t90+t21*(t480*t94
+t9*t253)*t125)*t10*t430*t635+t638*t610*(t22*t253*t640+t113*t58*t436*t262)+t21*
t278*(t649*t10*t650+t653*t81*t640*t610+t21*t58*t609*t610))*t89)*t398*t98*t318*
t263/t670)*t31)*t96+t103*t453*(t21*(t682*t104*t81*t456/(t488*t308+t2*t307*t294)
+(t682*t81*t692+(t695*t157+t388*t304*t698*t699)*t98)*t99)*t96+(t113*t157*(t710*
t713*t722+t725*t294)*t184*t60/t730+((t113*(t304*t58*t95*t7*t713*t722+t725*t739)
*t184*t699*t295+t388*(t746*t58*t95*t308*t467+t508*t491*t10*t2*t294+t342*t307*(
t488*(t21*t461*t456*t2+t469)+t388*t58*t2*t467*t294))*t96*t295/t791)*t99+t113*(
t113*t797*t157*t295+(t698*t308*t60*t699+t301*t104*t477*(t772+t775+t777+t780+
t537+t543+t557+t558)*(t649*t95*t307+t2*(t535*t158+t541*t163+t555*t157*t26+t426*
t292+t403*t95*t485))*t99*t295/t821)*t96))*t692)*t98)*t837*t89*t317)*t59)*t844-
t103*t29*(t848-t849*t578*t853*t62*t596/(-t272*t62+t864))*t837*t872*t873*t63*t50
-t1*t880*(t899*t98*t873*t950+t880*(t113*t714*t69*t885*t898*t957*t961+t113*(t897
*(t113*t104*t9*t965*t976+t979*t9*t7*t984*t863+(t988*t7*(t208+t455)+t988*t92*t7*
t137+t23*t996+t423*t92*t995*t15+t23*t546*t1001)*t1007)*t31*t1012*t98*t99+t746*
t984*t1025*t880*t1029)*t1036*t949+(t979*t1039*t898*t957*t873*t961+t388*(t113*
t1045*t1047+t104*t1024*t880*(t965*t863+t342*t7*t880)*t1029)*t31*t873/t1060+(
t681*t104*t10*t976*t863*t31/t1071+t899*(t1076*t143*t9*t921+t695*t9*t1068+(t979*
t995+t147*(t653*t15+t1076*t9*t26))*t880)*t1012*t873*t960)*t99*t1035+t113*(t21*
t1098*t898*t1012*t60*t949+t638*t863*(t864+t460*t880)*t1029+t360*t1025*(t922+
t945+t1109*t880)*(t144*t1115*t921+t66*t1115*t944+(t21*t995*t1109+t147*(t147*t15
*t26+t2*(t23*t15+t55*t9*t26)))*t880)*t31*t950/t1135)*t1036)*t98)/t681)*t56/(
t120*t90+t147*t93)*t63*t1157)-t1*(-t1163*t1164*t48*t29*t62*t1172/(t48*t1171-(
t69*t79*t58+t144*t73*t58+t120*(t248+t73*t130)+t147*(t79*t93+t73*t137))*t62)*t50
+(t605*(t21*t1164*t48*t1202/t79*t1205*t1172-t1209*t62)-t605*t1209*t62*t39/(t21*
t39+t1156)-t1164*t1227*t62*t50/(t261+t1230))*t1157)*t56*t872*t49/t1240-t1*(-t21
*t421*(t157*t12+t58*t72)*t4*t59/(t69*t58+t426+t9*t81)*t98*t13*t63-t1267*t453*
t352*t844*t59*t1205*t1272/t1274-t10*(t387*t2+t460*t112)*t898*t352*t1205*t98*t49
*t1290+t21*(-t68*t1267*t36*t352*t844*t59*t1272/(t261+t21*t1274)-t342*t112*(t21*
t352+t411*t2*t352*t98*t176)*(t1164*t421*t308*t1156/(t328+t1156)+t421*t36*t1319*
t1287/(t1287+t1230))*t340*t49*t1290)*t1205)*t56/t421*t99);
  }
#else

return
(

(complex<T>(0,1)*((SPB(3,4)*pow(SPB(4,7),2))/(SPA(5,3)*SPB(2,3)+
SPA(5,4)*SPB(2,4))+
(pow(SPA(5,3)*SPB(3,7)+
SPA(5,4)*SPB(4,7),2)*((pow(SPB(3,4),2)*pow(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7),2))/(SPA(4,3)*pow(SPA(5,3)*SPB(3,7)+
SPA(5,4)*SPB(4,7),2))+
(SPB(3,4)*pow(SPB(4,5),2))/(complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SS(3,4,5))))/((complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*SS(3,4,5))))/(complex<T>(2,0)*SPB(1,2)*pow(SPB(3,4),2)*SPB(6,7))-
recursive*(-((complex<T>(-1,0)*pow(SPA(6,5),2)*((complex<T>(0,1)*(-(((SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4))*(SPA(6,2)*SPB(2,4)+
SPA(6,3)*SPB(3,4))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)+
SPA(6,5)*SPB(4,5)-
SPA(7,6)*SPB(4,7)))/(complex<T>(2,0)*pow(SPA(6,5),2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(2,3)*SPB(3,4)))+
(complex<T>(0,1)*((complex<T>(0,-1)*SPA(3,2)*SPB(2,4)*SPB(3,4))/(complex<T>(2,0)*SPB(2,3))+
(complex<T>(0,-1)*SPA(3,2)*SPB(2,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4)))/(complex<T>(3,0)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))))*pow((complex<T>(-1,0)*SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7)))/SPA(6,5)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*(SPA(6,1)*SPB(1,4)*SPB(2,6)+
complex<T>(-1,0)*SPA(6,1)*SPB(1,2)*SPB(4,6)+
SPB(2,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,4)*SPB(2,7)+
SPA(6,1)*SPA(7,5)*SPB(1,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,1)*SPB(1,2)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,1)*SPB(1,4)*SPB(3,6)+
complex<T>(-1,0)*SPA(6,1)*SPB(1,3)*SPB(4,6)+
SPB(3,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,4)*SPB(3,7)+
SPA(6,1)*SPA(7,5)*SPB(1,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,1)*SPB(1,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(3,7)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))))/(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4)),3))/(complex<T>(1,0)*pow(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*(SPA(6,1)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(4,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))*((SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7)))/SPA(6,5)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)*SPB(2,6)+
SPA(6,1)*SPB(1,2)*SPB(4,6)+
complex<T>(-1,0)*SPB(2,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
(complex<T>(-1,0)*(SPA(6,1)*SPA(7,5)*SPB(1,4)*SPB(2,7)+
complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,2)*SPB(4,7)+
SPA(5,1)*SPB(1,2)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
SPA(7,5)*SPB(2,7)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))+
complex<T>(-1,0)*SPA(5,3)*(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)*SPB(3,6)+
SPA(6,1)*SPB(1,3)*SPB(4,6)+
complex<T>(-1,0)*SPB(3,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
(complex<T>(-1,0)*(SPA(6,1)*SPA(7,5)*SPB(1,4)*SPB(3,7)+
complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,3)*SPB(4,7)+
SPA(5,1)*SPB(1,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
SPA(7,5)*SPB(3,7)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))))/(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))))))/(complex<T>(1,0)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))+
((complex<T>(0,1)*SPA(3,2)*SPA(5,3)*(SPA(5,3)*(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3))+
SPA(5,4)*(SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4)))*(SPA(4,3)*SPB(2,4)+
(complex<T>(-1,0)*SPA(5,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7)))/SPA(6,5)))/(complex<T>(3,0)*SPA(4,3)*SPA(6,5)*pow(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4),2)*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7))))+
((complex<T>(0,1)*SPB(1,4)*pow(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7),2)*pow(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4),2))/(complex<T>(2,0)*pow(SPA(6,5),2)*SPB(1,2)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5))*SPB(3,4)*(SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))+
(complex<T>(0,1)*((complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4))*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7)))/pow(SPA(6,5),2)+
(SPA(6,1)*SPB(1,4)*(SPA(6,2)*SPB(2,4)+
SPA(6,3)*SPB(3,4))*(complex<T>(-1,0)*SPA(6,5)*SPB(5,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(7,6)*SPB(6,7)))/pow(SPA(6,5),2)))/(complex<T>(2,0)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5))*SPB(2,3)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))))/SPB(2,4)+
((complex<T>(0,1)*(complex<T>(-1,0)*SPA(3,2)*SPA(6,1)*SPB(1,4)*(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)+
SPA(6,5)*SPB(4,5)-
SPA(7,6)*SPB(4,7))*(SPA(6,1)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(4,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))+
SPA(3,2)*(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)+
SPA(6,5)*SPB(4,5)-
SPA(7,6)*SPB(4,7))*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))*(SPA(6,1)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(4,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))))/(complex<T>(2,0)*SPA(3,2)*pow(SPA(6,5),2)*SPB(2,4)*(SPA(6,1)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(4,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))+
((complex<T>(0,-1)*SPA(6,1)*(SPA(3,2)*SPB(2,4)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)*SPB(1,6)*SPB(2,3)+
SPA(6,1)*SPB(1,2)*SPB(1,6)*SPB(3,4)+
complex<T>(-1,0)*SPB(1,6)*SPB(2,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
(complex<T>(-1,0)*(SPA(6,1)*SPA(7,5)*SPB(1,4)*SPB(1,7)*SPB(2,3)+
complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,2)*SPB(1,7)*SPB(3,4)+
SPA(7,5)*SPB(1,7)*SPB(2,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))+
(complex<T>(2,0)*SPA(3,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(2,3)*SPB(2,4)+
SPA(3,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*(complex<T>(-3,0)*SPB(1,4)*SPB(2,3)+
SPB(1,3)*SPB(2,4)))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)*SPB(2,6)+
SPA(6,1)*SPB(1,2)*SPB(4,6)+
complex<T>(-1,0)*SPB(2,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
(complex<T>(-1,0)*(SPA(6,1)*SPA(7,5)*SPB(1,4)*SPB(2,7)+
complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,2)*SPB(4,7)+
SPA(5,1)*SPB(1,2)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
SPA(7,5)*SPB(2,7)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))))/(complex<T>(6,0)*SPA(6,5)*pow(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4),2)*(SPA(6,1)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(4,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))+
((complex<T>(0,-1)*SPA(3,2)*((complex<T>(-2,0)*SPA(5,4)*SPB(1,2)*SPB(2,4)*(SPA(5,3)*(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3))+
SPA(5,4)*(SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4)))*(SPA(4,3)*SPB(2,4)+
(complex<T>(-1,0)*SPA(5,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7)))/SPA(6,5))*SPB(3,4))/SPA(6,5)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*(complex<T>(-1,0)*(SPA(4,3)*SPB(1,4)+
(complex<T>(-1,0)*SPA(5,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7)))/SPA(6,5))*SPB(2,3)*(SPA(4,3)*SPB(2,4)+
(complex<T>(-1,0)*SPA(5,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7)))/SPA(6,5))*SPB(3,4)+
((SPA(5,3)*(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3))+
SPA(5,4)*(SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4)))*(SPA(4,3)*SPB(1,3)*pow(SPB(2,4),2)+
complex<T>(-2,0)*SPA(4,3)*SPB(1,2)*SPB(2,4)*SPB(3,4)+
(complex<T>(-1,0)*(SPA(5,3)*SPB(1,3)*SPB(2,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7))+
SPA(5,3)*SPB(1,2)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7))*SPB(3,4)+
complex<T>(-3,0)*SPA(5,3)*SPB(1,2)*SPB(2,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(3,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(3,7))))/SPA(6,5)))/SPA(6,5))))/(complex<T>(6,0)*SPA(4,3)*pow(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4),2)*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7))))+
(complex<T>(0,1)*pow(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7),2)*(complex<T>(-1,0)*SPA(3,2)*pow(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*(complex<T>(-1,0)*(SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7))*(complex<T>(-1,0)*SPA(2,1)*SPB(1,3)*SPB(2,4)+
(complex<T>(-1,0)*SPA(2,1)*SPB(1,2)+
SPA(3,2)*SPB(2,3))*SPB(3,4))+
complex<T>(-3,0)*SPA(2,1)*SPB(1,2)*SPB(2,4)*(SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7)))*pow(SPA(6,1)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(4,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5),2)+
complex<T>(-1,0)*SPA(2,1)*pow(SPB(3,4),2)*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))*(complex<T>(2,0)*SPA(3,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(2,4)*(SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7))*(SPA(6,4)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,4)*SPA(7,5)*SPB(4,7)+
SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))+
(complex<T>(-1,0)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(SPA(2,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3))*SPB(2,4)+
SPA(3,2)*(complex<T>(-3,0)*SPB(1,4)*SPB(2,3)+
SPB(1,3)*SPB(2,4))*(SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7)))*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7))))+
complex<T>(-1,0)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(3,4)*(SPA(6,1)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(4,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))*(complex<T>(-2,0)*SPA(2,1)*SPA(3,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(2,4)*(SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7))*(SPA(6,4)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,4)*SPA(7,5)*SPB(4,7)+
SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))+
(SPA(2,1)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(complex<T>(-1,0)*SPA(2,1)*SPB(1,2)+
SPA(3,2)*SPB(2,3))*SPB(2,4)+
complex<T>(-1,0)*SPA(3,2)*(complex<T>(-1,0)*(SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7))*(complex<T>(-3,0)*SPA(2,1)*SPB(1,4)*SPB(2,3)+
complex<T>(-2,0)*SPA(2,1)*SPB(1,3)*SPB(2,4)+
(SPA(2,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3))*SPB(3,4))+
complex<T>(3,0)*SPA(2,1)*SPB(1,2)*SPB(2,4)*(SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7))))*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7))))))/(complex<T>(6,0)*pow(SPA(6,5),2)*(SPA(6,1)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(4,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5))*(SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))*pow(complex<T>(-1,0)*SPA(6,1)*SPA(6,5)*SPB(1,6)*SPB(4,6)+
complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)*SPB(4,6)+
SPA(4,3)*SPA(6,5)*SPB(3,4)*SPB(4,6)+
complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,6)*SPB(4,7)+
SPA(4,3)*SPA(7,5)*SPB(3,4)*SPB(4,7)+
SPA(5,1)*SPB(1,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
(complex<T>(-1,0)*(SPA(6,1)*pow(SPA(7,5),2)*SPB(1,7)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,1)*SPA(7,5)*SPB(1,7)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5)+
SPA(5,3)*SPA(7,5)*SPB(3,4)*SPB(5,7)+
SPA(5,3)*SPB(3,4)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)),2)))/SPB(3,4))/SPB(1,2))/SPB(2,3))/(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))+
((((complex<T>(0,-1)*SPB(1,4)*pow(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7),2))/(complex<T>(2,0)*pow(SPA(6,5),2)*(SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))+
(complex<T>(0,1)*SPB(1,4)*((pow(SPA(6,1),2)*pow(SPB(1,4),2))/pow(SPA(6,5),2)+
(complex<T>(2,0)*SPA(6,1)*SPB(1,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7)))/pow(SPA(6,5),2)+
pow(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7),2)/pow(SPA(6,5),2)))/(complex<T>(2,0)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))))/(SPB(2,4)*SPB(3,4))+
(complex<T>(-1,0)*(((complex<T>(0,-1)*SPA(6,1)*pow(SPA(4,3)*SPB(1,4)+
(complex<T>(-1,0)*SPA(5,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7)))/SPA(6,5),2))/(complex<T>(2,0)*(SPA(5,3)*(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3))+
SPA(5,4)*(SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4))))+
(complex<T>(0,1)*SPA(3,2)*(SPA(4,3)*SPB(1,4)+
(complex<T>(-1,0)*SPA(5,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7)))/SPA(6,5))*(SPA(4,3)*SPB(2,4)+
(complex<T>(-1,0)*SPA(5,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7)))/SPA(6,5)))/(complex<T>(2,0)*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))))/(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))+
(complex<T>(0,-1)*((pow(SPA(3,2),2)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*SPB(2,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7))*SPB(3,4)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))/pow(SPA(6,5),2)+
(complex<T>(-1,0)*SPA(4,3)*SPB(1,2)*SPB(2,4)*SPB(3,4)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*pow(SPA(2,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))+
complex<T>(-1,0)*SPA(6,2)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(6,2)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)),2))/pow(SPA(6,5),2)+
(complex<T>(-1,0)*SPA(3,2)*(SPA(2,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))+
complex<T>(-1,0)*SPA(6,2)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(6,2)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))*(complex<T>(-2,0)*(SPA(6,4)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,4)*SPA(7,5)*SPB(1,7)+
SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5))*SPB(2,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7))*SPB(3,4)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(4,3)*SPB(3,4)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))+
complex<T>(-1,0)*SPB(2,4)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*(SPA(4,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*SPB(2,3)*SPB(3,4)+
complex<T>(3,0)*SPB(1,2)*(complex<T>(-1,0)*SPA(6,5)*SPB(3,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(3,7))*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,5)*SPB(2,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(2,7))*(complex<T>(-1,0)*SPB(1,3)*SPB(2,4)*pow(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(4,3)*SPB(3,4)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7),2)+
complex<T>(-3,0)*SPA(4,3)*SPB(1,4)*SPB(2,3)*SPB(3,4)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))+
complex<T>(-1,0)*SPB(1,2)*SPB(3,4)*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))))/pow(SPA(6,5),2)))/(complex<T>(6,0)*SPB(2,3)*SPB(3,4)*(SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*pow(SPA(6,1)*SPB(1,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,1)*SPA(7,5)*SPB(1,7)+
SPA(5,1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))))/SPA(6,5)+
complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7),2))))/SPA(4,3))/(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))+
(complex<T>(0,-1)*pow(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7),2)*((complex<T>(-1,0)*((complex<T>(9,0)*SPA(2,1)*SPA(3,2)*SPB(1,4)*(SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7)))/(complex<T>(-1,0)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*pow(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4),2)+
SPB(3,4)*(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4))*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7))))+
((complex<T>(9,0)*SPA(2,1)*SPB(1,4))/(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4))+
(complex<T>(23,0)*SPB(1,3)+
(complex<T>(-3,0)*SPA(3,1)*pow(SPB(1,3),2))/(SPA(6,4)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,4)*SPA(7,5)*SPB(4,7)+
SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))/SPB(2,3))/SPB(3,4)))/(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))+
((complex<T>(3,0)*SPB(1,3)*(SPB(1,2)*SPB(2,3)*(SPA(2,1)*SPB(1,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(3,4))*(SPA(4,2)*SPA(6,5)*SPB(4,6)+
SPA(4,2)*SPA(7,5)*SPB(4,7)+
SPA(5,2)*SPA(7,5)*SPB(5,7)+
SPA(5,2)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))+
pow(SPB(1,3),2)*(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4))*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))))/(pow(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*pow(SPB(3,4),2)*pow(SPA(6,4)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,4)*SPA(7,5)*SPB(4,7)+
SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7),2))+
(((complex<T>(3,0)*(SPA(3,1)*SPB(1,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(2,3)*(SPA(2,1)*SPB(1,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(3,4))*(SPA(4,2)*SPA(6,5)*SPB(4,6)+
SPA(4,2)*SPA(7,5)*SPB(4,7)+
SPA(5,2)*SPA(7,5)*SPB(5,7)+
SPA(5,2)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))+
pow(SPB(1,3),2)*(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4))*pow(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)),2)))/(pow(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*(SPA(6,4)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,4)*SPA(7,5)*SPB(4,7)+
SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7))))+
(complex<T>(-3,0)*(complex<T>(-3,0)*SPA(3,2)*SPB(1,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*pow(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4),2)*(SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7))+
SPA(2,1)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(SPA(2,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3))*SPB(2,4)*SPB(3,4)*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))+
complex<T>(-1,0)*SPA(3,2)*(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4))*(complex<T>(-1,0)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(2,1)*SPB(1,2)+
SPA(3,2)*SPB(2,3))*(SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7))*SPB(3,4)+
complex<T>(-3,0)*SPA(2,1)*SPB(1,2)*SPB(2,4)*(SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7)))+
complex<T>(-3,0)*SPB(1,2)*SPB(3,4)*(SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7))*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7))))))/((SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))*(SPA(2,1)*SPA(6,5)*SPB(1,6)*SPB(2,4)+
SPA(2,1)*SPA(7,5)*SPB(1,7)*SPB(2,4)+
SPA(3,1)*SPA(6,5)*SPB(1,6)*SPB(3,4)+
SPA(3,1)*SPA(7,5)*SPB(1,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(3,4)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(3,4)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(3,4)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(3,4)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))))/SPB(3,4)+
complex<T>(3,0)*((complex<T>(3,0)*pow(SPA(3,1),2)*SPB(1,3))/(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))+
((pow(SPB(1,3),2)*pow(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4),2))/(pow(SPB(3,4),2)*(SPA(6,4)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,4)*SPA(7,5)*SPB(4,7)+
SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))+
(SPA(2,1)*SPA(3,2)*SPB(2,4)*(SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7))*(SPA(2,1)*SPA(6,5)*SPB(1,6)*SPB(2,4)+
SPA(2,1)*SPA(7,5)*SPB(1,7)*SPB(2,4)+
SPA(3,1)*SPA(6,5)*SPB(1,6)*SPB(3,4)+
SPA(3,1)*SPA(7,5)*SPB(1,7)*SPB(3,4)+
SPA(4,3)*SPA(6,5)*SPB(3,4)*SPB(4,6)+
SPA(4,3)*SPA(7,5)*SPB(3,4)*SPB(4,7)+
SPA(5,3)*SPA(7,5)*SPB(3,4)*SPB(5,7)+
SPA(5,3)*SPB(3,4)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))*(complex<T>(-1,0)*SPB(1,3)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4))+
SPB(3,4)*(SPA(4,3)*SPA(6,5)*SPB(1,3)*SPB(4,6)+
SPA(4,3)*SPA(7,5)*SPB(1,3)*SPB(4,7)+
SPA(5,3)*SPA(7,5)*SPB(1,3)*SPB(5,7)+
SPA(5,3)*SPB(1,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7))+
complex<T>(-2,0)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(SPA(6,4)*SPB(4,6)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,4)*SPA(7,5)*SPB(4,7)+
SPA(5,4)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/SPA(6,5)+
SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)))))/(SPB(3,4)*(complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)))*pow(SPA(2,1)*SPA(6,5)*SPB(1,6)*SPB(2,4)+
SPA(2,1)*SPA(7,5)*SPB(1,7)*SPB(2,4)+
SPA(3,1)*SPA(6,5)*SPB(1,6)*SPB(3,4)+
SPA(3,1)*SPA(7,5)*SPB(1,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(3,4)*SPB(4,6)+
complex<T>(-1,0)*SPA(4,3)*SPA(7,5)*SPB(3,4)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPA(7,5)*SPB(3,4)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(3,4)*(SPA(6,5)*SPB(5,6)+
SPA(7,6)*SPB(6,7)),2)))/(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))))/(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4)))/SPB(2,3)))/(complex<T>(18,0)*pow(SPA(6,5),2)*(SPA(6,5)*SPB(5,6)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7))))/SPB(1,2)))/SPA(7,6))-
(complex<T>(0,-1)*pow(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7),2)*(complex<T>(14,0)-
(complex<T>(-9,0)*pow(SPA(6,1),2)*pow(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*SPB(6,7))/((SPA(5,3)*(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3))+
SPA(5,4)*(SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4)))*(-((SPA(5,3)*(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3))+
SPA(5,4)*(SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4)))*SPB(6,7))+
(complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))))/(complex<T>(18,0)*SPA(4,3)*SPB(1,2)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*SPB(6,7)*SS(3,4,5))-
(complex<T>(0,1)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))*((complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,3)+
complex<T>(-1,0)*SPA(6,5)*SPB(3,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(3,7))*pow(SPA(3,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,4)*SPB(3,7)+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))*SPB(5,7),2))/(SPB(2,3)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*SPB(3,4)*(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))+
((SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))*((complex<T>(3,0)*SPA(4,2)*SPA(5,2)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,3)+
complex<T>(-1,0)*SPA(6,5)*SPB(3,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(3,7))*pow(SPA(3,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,4)*SPB(3,7)+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))*SPB(5,7),2))/(pow(SPA(5,4),2)*pow(SPB(3,4),2)*pow(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)),2))+
(complex<T>(3,0)*(((SPA(3,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,4)*SPB(3,7)+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))*SPB(5,7))*(complex<T>(3,0)*SPA(3,2)*SPA(5,4)*(SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7))*(complex<T>(-1,0)*SPA(4,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*(SPA(3,2)*SPB(2,3)+
SPA(4,3)*SPB(3,4))*SPB(3,7)+
complex<T>(-1,0)*SPA(4,2)*SPB(2,3)*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4))*SPB(5,7))+
complex<T>(3,0)*SPA(4,3)*SPA(5,4)*SPB(2,3)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,4)+
complex<T>(-1,0)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(4,7))*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
(SPA(4,2)*SPA(5,2)*SPB(2,3)*(SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6))+
SPA(4,2)*SPA(5,2)*SPA(7,5)*SPB(2,3)*SPB(2,7)+
complex<T>(-1,0)*SPA(4,3)*pow(SPA(5,1)*SPB(1,3)+
SPA(6,5)*SPB(3,6),2)+
complex<T>(-2,0)*SPA(4,3)*SPA(7,5)*(SPA(5,1)*SPB(1,3)+
SPA(6,5)*SPB(3,6))*SPB(3,7)+
complex<T>(-1,0)*SPA(4,3)*pow(SPA(7,5),2)*pow(SPB(3,7),2))*(SPA(3,1)*SPB(1,7)-
SPA(6,3)*SPB(6,7))))/(SPA(4,3)*SPA(5,4)*SPB(2,3)*SPB(3,4))+
(complex<T>(-3,0)*SPA(3,2)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,4)+
complex<T>(-1,0)*SPA(6,5)*SPB(4,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(4,7))*(complex<T>(-1,0)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))*SPB(2,7)+
complex<T>(-1,0)*SPA(4,3)*SPB(2,4)*SPB(3,7)+
SPA(4,3)*SPB(2,3)*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*SPB(5,7))*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))/(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
(SPA(4,3)*SPB(3,4)+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))))/((complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))*(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))+
((complex<T>(3,0)*SPA(4,3)*pow(complex<T>(-1,0)*SPA(5,1)*SPB(1,3)+
complex<T>(-1,0)*SPA(6,5)*SPB(3,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(3,7),3)*pow(SPA(3,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,4)*SPB(3,7)+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))*SPB(5,7),2))/(pow(SPA(5,4),2)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*pow(SPB(3,4),2)*pow(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)),2))+
(complex<T>(-3,0)*(complex<T>(3,0)*(SPA(5,1)*SPB(1,3)+
SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7))*pow(SPA(3,1)*SPB(1,7)-
SPA(6,3)*SPB(6,7),2)+
(SPA(3,2)*(complex<T>(-1,0)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))*SPB(2,7)+
complex<T>(-1,0)*SPA(4,3)*SPB(2,4)*SPB(3,7)+
SPA(4,3)*SPB(2,3)*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*SPB(5,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))*((SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6)+
SPA(7,5)*SPB(2,7))*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
(SPA(4,3)*SPB(3,4)+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))))/(SPA(4,3)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*pow(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7),2))+
((complex<T>(9,0)*SPA(3,2)*SPB(2,4)*(complex<T>(-1,0)*SPA(4,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*(SPA(3,2)*SPB(2,3)+
SPA(4,3)*SPB(3,4))*SPB(3,7)+
complex<T>(-1,0)*SPA(4,2)*SPB(2,3)*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4))*SPB(5,7))*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7)))/(SPA(4,3)*(complex<T>(-1,0)*SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
(complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(7,5)*SPB(5,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,3)+
complex<T>(-1,0)*SPA(6,5)*SPB(3,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(3,7))*pow(SPA(3,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,4)*SPB(3,7)+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))*SPB(5,7),2)*(complex<T>(-23,0)*SPA(5,1)*SPA(5,4)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(23,0)*SPA(5,4)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
(complex<T>(3,0)*SPA(4,3)*(SPA(5,1)*SPB(1,3)+
SPA(6,5)*SPB(3,6))+
complex<T>(-1,0)*SPA(7,5)*(complex<T>(-3,0)*SPA(4,3)*SPB(3,7)+
complex<T>(-23,0)*SPA(5,4)*SPB(5,7)))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/(SPA(5,4)*(complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*pow(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)),2)))/(SPB(3,4)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))+
(complex<T>(3,0)*((complex<T>(-1,0)*pow(SPA(5,1)*SPB(1,3)+
SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7),2)*pow(SPA(3,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,4)*SPB(3,7)+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))*SPB(5,7),2))/(SPA(5,4)*pow(SPB(3,4),2)*(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))+
(complex<T>(-1,0)*SPB(2,4)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))*((complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
SPA(3,2)*SPB(2,3)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
(SPA(4,3)*SPB(3,4)+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))+
(SPA(3,2)*SPB(2,4)*(complex<T>(-1,0)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))*SPB(2,7)+
complex<T>(-1,0)*SPA(4,3)*SPB(2,4)*SPB(3,7)+
SPA(4,3)*SPB(2,3)*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*SPB(5,7))*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))*(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
(complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))*(complex<T>(-1,0)*SPA(5,1)*(SPA(5,1)*SPB(1,3)+
complex<T>(-2,0)*SPA(5,4)*SPB(3,4)+
SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7))*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
SPA(6,5)*(SPA(5,1)*SPB(1,3)+
complex<T>(-2,0)*SPA(5,4)*SPB(3,4)+
SPA(6,5)*SPB(3,6)+
SPA(7,5)*SPB(3,7))*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
(complex<T>(-1,0)*(SPA(5,1)*SPB(1,3)+
SPA(6,5)*SPB(3,6))*(complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7))+
complex<T>(-1,0)*SPA(7,5)*(complex<T>(-1,0)*SPA(7,5)*SPB(3,7)*SPB(5,7)+
SPB(3,4)*(complex<T>(-1,0)*SPA(4,3)*SPB(3,7)+
complex<T>(2,0)*SPA(5,4)*SPB(5,7))))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/(SPA(4,3)*SPB(3,4)*(SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
complex<T>(-1,0)*SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
complex<T>(-1,0)*SPA(7,5)*SPB(5,7)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))*pow(complex<T>(-1,0)*SPA(5,1)*(complex<T>(-1,0)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(2,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*SPB(3,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(4,7)+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(5,7))+
SPA(6,5)*(complex<T>(-1,0)*SPB(2,7)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,6)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,6))+
complex<T>(-1,0)*SPB(3,7)*(SPA(3,2)*SPB(2,6)+
complex<T>(-1,0)*SPA(4,3)*SPB(4,6))+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,6)+
SPA(4,3)*SPB(3,6))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,6)+
SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SPB(5,7))+
(complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(7,5)*SPB(5,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)),2))))/((complex<T>(-1,0)*SPA(5,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*SPB(2,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(2,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(2,3)))/complex<T>(9,0)))/(complex<T>(2,0)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(1,7))*SPB(6,7)*SSS(2,3,4,5)))-
(complex<T>(0,1)*(-((pow(SPA(2,1),2)*pow(SPB(1,7),2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*pow(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7),2)*SPB(6,7))/((complex<T>(-1,0)*SPA(3,2)*SPB(3,7)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7))*((complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*(complex<T>(-1,0)*SPA(3,2)*SPB(3,7)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7))-
(SPA(5,2)*SPA(6,1)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,1)*SPA(6,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,5)*(SPA(6,1)*SPB(1,6)+
SPA(6,2)*SPB(2,6))+
complex<T>(-1,0)*SPA(7,5)*(SPA(6,1)*SPB(1,7)+
SPA(6,2)*SPB(2,7)))*SPB(6,7))*SS(3,4,5)))+
(pow(SPA(3,2),2)*((complex<T>(-1,0)*pow(SPB(1,7),2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*pow(complex<T>(-1,0)*SPA(2,1)*SPB(2,7)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,7)+
complex<T>(-1,0)*SPA(4,1)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,1)*SPB(5,7),2))/(SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*(complex<T>(-1,0)*SPA(3,2)*SPB(3,7)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7)))-
pow(SPB(2,7),2)*SPB(6,7))-
(pow(SPA(3,2),2)*pow(SPB(2,7),2)*SPB(6,7)*SS(3,4,5))/(complex<T>(-1,0)*SS(3,4,5)+
SSS(2,3,4,5))-
(pow(SPB(1,7),2)*pow(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(4,3)*SPB(4,6)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,6))+
complex<T>(-1,0)*SPA(7,1)*(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7)),2)*SPB(6,7))/(SS(3,4,5)*(SPA(7,6)*SPB(6,7)+
complex<T>(-1,0)*SSS(2,3,4,5))))/SSS(2,3,4,5)))/(complex<T>(2,0)*SPA(4,3)*SPB(1,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*pow(SPB(6,7),2))-
(complex<T>(0,1)*(-((complex<T>(-1,0)*pow(SPB(2,4),2)*(SPB(1,3)*(SPA(5,3)*SPB(2,3)+
SPA(5,4)*SPB(2,4))+
SPB(1,2)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))*pow(SPB(4,7),2))/(SPB(1,2)*(SPA(5,2)*SPB(1,2)+
SPA(5,3)*SPB(1,3)+
SPA(5,4)*SPB(1,4))*SPB(2,3)*(SPA(5,3)*SPB(2,3)+
SPA(5,4)*SPB(2,4))*SPB(6,7)))-
(SPB(2,4)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)*SPB(2,3)*SPB(3,4)+
SPA(3,2)*SPB(2,3)*(complex<T>(-1,0)*SPB(1,3)*SPB(2,4)+
SPB(1,2)*SPB(3,4)))*pow(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7),2))/(SPA(3,2)*SPA(7,6)*SPB(1,2)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*pow(SPB(2,3),2)*SS(5,6,7))-
(SPB(2,4)*(SPA(3,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*SPB(3,4)+
SPA(3,2)*SPB(2,3)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4)))*pow(SPA(3,2)*SPB(2,7)*SPB(3,4)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,4)*SPB(3,7)+
complex<T>(-1,0)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4))*SPB(4,7)+
complex<T>(-1,0)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))*SPB(5,7),2))/(SPA(3,2)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(2,3)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*SPB(6,7)*SSS(1,5,6,7)*SSS(2,3,4,5))+
(complex<T>(-1,0)*(-((complex<T>(-1,0)*pow(SPA(6,5),2)*SPB(2,4)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)*SPB(2,3)*SPB(3,4)+
SPA(3,2)*SPB(2,3)*(complex<T>(-1,0)*SPB(1,3)*SPB(2,4)+
SPB(1,2)*SPB(3,4)))*pow(SPB(4,5),2))/(SPA(3,2)*SPA(7,6)*SPB(1,2)*pow(SPB(2,3),2)*(SPA(7,6)*SPB(6,7)+
complex<T>(-1,0)*SS(5,6,7))))-
(complex<T>(-1,0)*SPA(3,2)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))*(complex<T>(-1,0)/SPA(3,2)+
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*SPB(3,4))/(SPA(3,2)*SPB(2,3)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))))*((pow(SPB(1,7),2)*pow(SPB(2,4),2)*pow(complex<T>(-1,0)*SPA(2,1)*SPB(2,4)+
complex<T>(-1,0)*SPA(3,1)*SPB(3,4),2)*SSS(2,3,4,5))/(complex<T>(-1,0)*SPA(7,6)*SPB(6,7)+
SSS(2,3,4,5))+
(pow(SPB(2,4),2)*pow(SPB(4,5),2)*pow(SPA(5,2)*SPB(2,7)+
SPA(5,3)*SPB(3,7)+
SPA(5,4)*SPB(4,7),2)*SSS(1,5,6,7))/(SSS(1,5,6,7)+
complex<T>(-1,0)*SSS(2,3,4,5))))/(SPB(2,4)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))*SPB(6,7)*SSS(1,5,6,7)*SSS(2,3,4,5))))/(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))))/(complex<T>(2,0)*pow(SPB(2,4),2)*SPB(3,4))
 )
;


#endif

      }



template <class T> complex<T> R2q2Q1g2l_qmQbmQppqbplmlbp_lc
      (const eval_param<T>& ep,
                 const mass_param_coll& mpc){
//{{qm, Qbm, Qp,p, qbp, lm, lbp}, lc}

#if _VERBOSE
  _MESSAGE("R2q2Q1g2l :  qmQbmQppqbplmlbp lc");
#endif

complex<T> recursive;
recursive=complex<T>(1,0);

#if _USE_OPTIMIZED

  complex<T> t1;
  complex<T> t10;
  complex<T> t100;
  complex<T> t1003;
  complex<T> t101;
  complex<T> t1015;
  complex<T> t102;
  complex<T> t1024;
  complex<T> t103;
  complex<T> t1030;
  complex<T> t1031;
  complex<T> t1034;
  complex<T> t104;
  complex<T> t1045;
  complex<T> t1046;
  complex<T> t105;
  complex<T> t1055;
  complex<T> t1056;
  complex<T> t1058;
  complex<T> t1070;
  complex<T> t1079;
  complex<T> t1081;
  complex<T> t1094;
  complex<T> t1099;
  complex<T> t11;
  complex<T> t1105;
  complex<T> t111;
  complex<T> t1121;
  complex<T> t1130;
  complex<T> t1159;
  complex<T> t116;
  complex<T> t1160;
  complex<T> t1166;
  complex<T> t1167;
  complex<T> t1170;
  complex<T> t1173;
  complex<T> t1174;
  complex<T> t1178;
  complex<T> t118;
  complex<T> t1195;
  complex<T> t1197;
  complex<T> t1198;
  complex<T> t120;
  complex<T> t1202;
  complex<T> t1207;
  complex<T> t1208;
  complex<T> t121;
  complex<T> t1224;
  complex<T> t1231;
  complex<T> t1233;
  complex<T> t1243;
  complex<T> t1267;
  complex<T> t1269;
  complex<T> t1270;
  complex<T> t1274;
  complex<T> t1275;
  complex<T> t1281;
  complex<T> t129;
  complex<T> t1291;
  complex<T> t1292;
  complex<T> t13;
  complex<T> t131;
  complex<T> t1318;
  complex<T> t1322;
  complex<T> t133;
  complex<T> t134;
  complex<T> t135;
  complex<T> t137;
  complex<T> t138;
  complex<T> t14;
  complex<T> t142;
  complex<T> t143;
  complex<T> t144;
  complex<T> t147;
  complex<T> t149;
  complex<T> t150;
  complex<T> t152;
  complex<T> t161;
  complex<T> t165;
  complex<T> t167;
  complex<T> t17;
  complex<T> t171;
  complex<T> t18;
  complex<T> t180;
  complex<T> t183;
  complex<T> t185;
  complex<T> t186;
  complex<T> t187;
  complex<T> t188;
  complex<T> t19;
  complex<T> t190;
  complex<T> t192;
  complex<T> t197;
  complex<T> t198;
  complex<T> t199;
  complex<T> t2;
  complex<T> t20;
  complex<T> t205;
  complex<T> t207;
  complex<T> t209;
  complex<T> t21;
  complex<T> t212;
  complex<T> t223;
  complex<T> t227;
  complex<T> t23;
  complex<T> t230;
  complex<T> t236;
  complex<T> t24;
  complex<T> t249;
  complex<T> t250;
  complex<T> t251;
  complex<T> t252;
  complex<T> t253;
  complex<T> t254;
  complex<T> t256;
  complex<T> t259;
  complex<T> t26;
  complex<T> t263;
  complex<T> t264;
  complex<T> t265;
  complex<T> t27;
  complex<T> t277;
  complex<T> t28;
  complex<T> t280;
  complex<T> t283;
  complex<T> t286;
  complex<T> t289;
  complex<T> t290;
  complex<T> t293;
  complex<T> t294;
  complex<T> t296;
  complex<T> t297;
  complex<T> t298;
  complex<T> t299;
  complex<T> t3;
  complex<T> t30;
  complex<T> t300;
  complex<T> t304;
  complex<T> t305;
  complex<T> t307;
  complex<T> t309;
  complex<T> t31;
  complex<T> t310;
  complex<T> t311;
  complex<T> t318;
  complex<T> t320;
  complex<T> t330;
  complex<T> t331;
  complex<T> t334;
  complex<T> t335;
  complex<T> t34;
  complex<T> t348;
  complex<T> t362;
  complex<T> t364;
  complex<T> t367;
  complex<T> t368;
  complex<T> t369;
  complex<T> t37;
  complex<T> t371;
  complex<T> t373;
  complex<T> t376;
  complex<T> t38;
  complex<T> t387;
  complex<T> t39;
  complex<T> t392;
  complex<T> t4;
  complex<T> t40;
  complex<T> t405;
  complex<T> t406;
  complex<T> t41;
  complex<T> t412;
  complex<T> t415;
  complex<T> t416;
  complex<T> t417;
  complex<T> t418;
  complex<T> t419;
  complex<T> t42;
  complex<T> t420;
  complex<T> t421;
  complex<T> t422;
  complex<T> t423;
  complex<T> t425;
  complex<T> t427;
  complex<T> t432;
  complex<T> t438;
  complex<T> t441;
  complex<T> t445;
  complex<T> t446;
  complex<T> t447;
  complex<T> t448;
  complex<T> t450;
  complex<T> t456;
  complex<T> t468;
  complex<T> t49;
  complex<T> t491;
  complex<T> t492;
  complex<T> t494;
  complex<T> t495;
  complex<T> t497;
  complex<T> t498;
  complex<T> t5;
  complex<T> t505;
  complex<T> t508;
  complex<T> t51;
  complex<T> t510;
  complex<T> t512;
  complex<T> t521;
  complex<T> t526;
  complex<T> t527;
  complex<T> t528;
  complex<T> t53;
  complex<T> t536;
  complex<T> t539;
  complex<T> t54;
  complex<T> t542;
  complex<T> t547;
  complex<T> t548;
  complex<T> t55;
  complex<T> t551;
  complex<T> t554;
  complex<T> t589;
  complex<T> t590;
  complex<T> t6;
  complex<T> t60;
  complex<T> t605;
  complex<T> t607;
  complex<T> t62;
  complex<T> t623;
  complex<T> t624;
  complex<T> t626;
  complex<T> t627;
  complex<T> t631;
  complex<T> t64;
  complex<T> t647;
  complex<T> t648;
  complex<T> t65;
  complex<T> t655;
  complex<T> t661;
  complex<T> t664;
  complex<T> t67;
  complex<T> t68;
  complex<T> t681;
  complex<T> t691;
  complex<T> t694;
  complex<T> t695;
  complex<T> t7;
  complex<T> t703;
  complex<T> t705;
  complex<T> t708;
  complex<T> t710;
  complex<T> t712;
  complex<T> t72;
  complex<T> t724;
  complex<T> t725;
  complex<T> t727;
  complex<T> t73;
  complex<T> t730;
  complex<T> t732;
  complex<T> t74;
  complex<T> t740;
  complex<T> t744;
  complex<T> t746;
  complex<T> t750;
  complex<T> t76;
  complex<T> t78;
  complex<T> t786;
  complex<T> t789;
  complex<T> t79;
  complex<T> t792;
  complex<T> t795;
  complex<T> t797;
  complex<T> t8;
  complex<T> t804;
  complex<T> t81;
  complex<T> t82;
  complex<T> t827;
  complex<T> t83;
  complex<T> t843;
  complex<T> t844;
  complex<T> t853;
  complex<T> t856;
  complex<T> t858;
  complex<T> t860;
  complex<T> t866;
  complex<T> t869;
  complex<T> t87;
  complex<T> t870;
  complex<T> t879;
  complex<T> t88;
  complex<T> t880;
  complex<T> t887;
  complex<T> t892;
  complex<T> t895;
  complex<T> t898;
  complex<T> t90;
  complex<T> t902;
  complex<T> t903;
  complex<T> t909;
  complex<T> t92;
  complex<T> t922;
  complex<T> t923;
  complex<T> t924;
  complex<T> t938;
  complex<T> t942;
  complex<T> t943;
  complex<T> t944;
  complex<T> t945;
  complex<T> t95;
  complex<T> t953;
  complex<T> t954;
  complex<T> t955;
  complex<T> t956;
  complex<T> t96;
  complex<T> t963;
  complex<T> t969;
  complex<T> t97;
  complex<T> t973;
  complex<T> t978;
  complex<T> t98;
  complex<T> t983;
  complex<T> t986;
  complex<T> t99;
  complex<T> t992;
  complex<T> t995;
  complex<T> t999;
  {
    t1 = complex<T>(0,-1);
    t2 = complex<T>(-1,0);
    t3 = SPA(3,2);
    t4 = t2*t3;
    t5 = SPA(6,2);
    t6 = pow(t5,2);
    t7 = SPA(4,2);
    t8 = SPB(1,2);
    t10 = SPA(4,3);
    t11 = SPB(1,3);
    t13 = t7*t8+t10*t11;
    t14 = complex<T>(1,0)/t13;
    t17 = t2*t5;
    t18 = t17*t8;
    t19 = SPA(6,3);
    t20 = t2*t19;
    t21 = t20*t11;
    t23 = pow(t18+t21,2);
    t24 = pow(t3,2);
    t26 = SPA(6,1);
    t27 = t26*t11;
    t28 = SPB(2,3);
    t30 = t27+t5*t28;
    t31 = pow(t30,2);
    t34 = complex<T>(1,0)/t28;
    t37 = SPA(2,1);
    t38 = complex<T>(2,0);
    t39 = pow(t37,2);
    t40 = t2*t39;
    t41 = t4*t28;
    t42 = SS(1,2,3);
    t49 = t2*t7;
    t51 = t2*t10;
    t53 = t49*t8+t51*t11;
    t54 = complex<T>(1,0)/t53;
    t55 = complex<T>(1,0)/t42;
    t60 = complex<T>(1,0)/t38;
    t62 = pow(t3,2);
    t64 = SPA(5,4);
    t65 = complex<T>(1,0)/t64;
    t67 = SPA(7,6);
    t68 = complex<T>(1,0)/t67;
    t72 = SPB(1,7);
    t73 = pow(t72,2);
    t74 = complex<T>(0,1);
    t76 = SPB(1,4);
    t78 = t4*t8+t10*t76;
    t79 = SPB(3,7);
    t81 = SPB(4,7);
    t82 = t49*t81;
    t83 = t4*t79+t82;
    t87 = SPA(5,2);
    t88 = SPB(5,7);
    t90 = SPB(6,7);
    t92 = t2*t37*t72+t87*t88-t17*t90;
    t95 = complex<T>(1,0)/t3;
    t96 = complex<T>(1,0)/t10;
    t97 = t95*t96;
    t98 = complex<T>(1,0)/t73;
    t99 = SPA(6,5);
    t100 = SPB(1,6);
    t101 = t99*t100;
    t102 = SPA(7,5);
    t103 = t102*t72;
    t104 = t101+t103;
    t105 = complex<T>(1,0)/t104;
    t111 = SPB(3,4);
    t116 = t3*t11;
    t118 = t116+t7*t76;
    t120 = complex<T>(3,0);
    t121 = complex<T>(1,0)/t120;
    t129 = t37*t72-t5*t90;
    t131 = complex<T>(1,0)/t72;
    t133 = SPA(5,3);
    t134 = SPA(7,2);
    t135 = t133*t134;
    t137 = t2*t87;
    t138 = SPA(7,3);
    t142 = t2*t133;
    t143 = t5*t100;
    t144 = t143*t88;
    t147 = t100*t88;
    t149 = SPB(1,5);
    t150 = t149*t129;
    t152 = t100*t129;
    t161 = SPA(7,4);
    t165 = t2*t64;
    t167 = SPA(6,4);
    t171 = t2*t167;
    t180 = complex<T>(1,0)/t118;
    t183 = pow(t8*t129*t131+t2*(-t11*(-t135*t88-t137*t138*t88-t138*t129-t2*(
t142*t144+t87*t19*t147+t142*t150+t20*t152)*t131)-t76*(-t64*t134*t88-t137*t161*
t88-t161*t129-t2*(t165*t144+t87*t167*t147+t165*t150+t171*t152)*t131))*t180,3)
;
    t185 = complex<T>(1,0);
    t186 = complex<T>(1,0)/t185;
    t187 = pow(t104,2);
    t188 = complex<T>(1,0)/t187;
    t190 = t134*t88;
    t192 = t2*t149;
    t197 = t190-t2*(-t17*t147-t192*t129)*t131;
    t198 = complex<T>(1,0)/t197;
    t199 = t2*t8;
    t205 = t2*t138;
    t207 = t133*t5;
    t209 = t19*t100;
    t212 = t133*t149;
    t223 = t2*t161;
    t227 = t167*t100;
    t230 = t64*t149;
    t236 = -t165*t190-t87*t161*t88-t223*t129-t2*(t64*t5*t147+t137*t227*t88+t230
*t129+t227*t129)*t131;
    t249 = t26*t100;
    t250 = SPA(7,1);
    t251 = t250*t72;
    t252 = t102*t88;
    t253 = t67*t90;
    t254 = t2*t99;
    t256 = SPA(5,1);
    t259 = t256*t72-t99*t90;
    t263 = t2*(-t254*t147-t192*t259)*t131;
    t264 = t249+t251+t252+t253-t263;
    t265 = complex<T>(1,0)/t264;
    t277 = -t199*(t82+t137*t88)-t2*t11*(t51*t81+t142*t88);
    t280 = SPA(4,1);
    t283 = t280*t72-t167*t90;
    t286 = t7*t28+t11*t283*t131;
    t289 = pow(t53,2);
    t290 = complex<T>(1,0)/t289;
    t293 = t27*t100;
    t294 = t143*t28;
    t296 = t134*t72*t28;
    t297 = t251+t253;
    t298 = t11*t297;
    t299 = -t293-t294-t296-t298;
    t300 = complex<T>(1,0)/t299;
    t304 = t1*t87;
    t305 = SPB(3,5);
    t307 = SPB(4,5);
    t309 = t4*t305+t49*t307;
    t310 = pow(t309,2);
    t311 = pow(t259,2);
    t318 = complex<T>(1,0)/(t249+t251+t253);
    t320 = complex<T>(1,0)/(t252-t263);
    t330 = t2*t26;
    t331 = t330*t100;
    t334 = t2*t67;
    t335 = t334*t90;
    t348 = complex<T>(1,0)/t7;
    t362 = complex<T>(1,0)/t111;
    t364 = t198*t265;
    t367 = t1*t88;
    t368 = complex<T>(-3,0);
    t369 = t368*t10;
    t371 = t7*t133;
    t373 = t2*(t369*t87+t371);
    t376 = complex<T>(-2,0);
    t387 = t3*t64;
    t392 = t101*t88;
    t405 = complex<T>(6,0);
    t406 = complex<T>(1,0)/t405;
    t412 = pow(t129,2);
    t415 = t120*t7;
    t416 = t138*t72;
    t417 = t209+t416;
    t418 = t64*t417;
    t419 = t418*t307;
    t420 = t415*t419;
    t421 = t161*t72;
    t422 = t227+t421;
    t423 = t2*t422;
    t425 = t10*t111;
    t427 = t425+t165*t307;
    t432 = pow(t197,2);
    t438 = t7*t104;
    t441 = t51*t111+t64*t307;
    t445 = t38*t7;
    t446 = t445*t422;
    t447 = t104*t111;
    t448 = SPB(2,7);
    t450 = t100*t448;
    t456 = t249+t251+t134*t448+t253-t2*(-t17*t450-t199*t129)*t131;
    t468 = t87*t307;
    t491 = t11*t100;
    t492 = t3*t26*t491;
    t494 = t100*t28;
    t495 = t3*t5*t494;
    t497 = t72*t28;
    t498 = t3*t134*t497;
    t505 = t2*t134;
    t508 = t2*t102;
    t510 = t116*t297;
    t512 = pow(t100,2);
    t521 = pow(t492+t495+t498+t254*t134*t100*t88+t17*t102*t100*t88+t505*t103*
t88+t508*t150+t510-t2*(-t5*t99*t512*t88-t99*t149*t152)*t131,2);
    t526 = t1*t111;
    t527 = t376*t3;
    t528 = t7*t64;
    t536 = t87*t28;
    t539 = t536+t11*t259*t131;
    t542 = pow(t7,2);
    t547 = t64*t11;
    t548 = SPA(3,1);
    t551 = t548*t72-t19*t90;
    t554 = t133*t11;
    t589 = pow(t87,2);
    t590 = pow(t88,2);
    t605 = pow(t539,2);
    t607 = complex<T>(1,0)/t277;
    t623 = -t330*t100*t81-t2*t307*t259-t2*t81*t297;
    t624 = pow(t623,2);
    t626 = t249+t251+t41+t252+t253-t263;
    t627 = t626*t98;
    t631 = pow(t111,2);
    t647 = t3*t28;
    t648 = t249+t251+t647+t252+t253-t263;
    t655 = t368*t64;
    t661 = pow(t648,2);
    t664 = t368*t3;
    t681 = pow(t626,2);
    t691 = complex<T>(-9,0);
    t694 = t111*t307;
    t695 = t2*t104;
    t703 = complex<T>(9,0);
    t705 = complex<T>(1,0)/t309;
    t708 = complex<T>(-23,0);
    t710 = pow(t133,2);
    t712 = complex<T>(1,0)/t456;
    t724 = t368*t133;
    t725 = t710*t309;
    t727 = t10*t64;
    t730 = t76*t100;
    t732 = SPB(2,4);
    t740 = (t4*t111+t468)*(-t330*t730-t17*t100*t732-t505*t72*t732-t2*t76*t297);
    t744 = complex<T>(1,0)/t24;
    t746 = pow(t456,2);
    t750 = pow(t299,2);
    t786 = t3*t99*t100*t305;
    t789 = t3*t102*t72*t305;
    t792 = t7*t99*t100*t307;
    t795 = t7*t102*t72*t307;
    t797 = t4*t293+t4*t294+t4*t296+t786+t789+t792+t795+t4*t298;
    t804 = pow(t305,2);
    t827 = pow(t797,2);
    t843 = complex<T>(18,0);
    t844 = complex<T>(1,0)/t843;
    t853 = complex<T>(1,0)/t90;
    t856 = complex<T>(14,0);
    t858 = t165*t76;
    t860 = pow(t858+t101+t103,2);
    t866 = t165*t149+t171*t100+t223*t72;
    t869 = t254*t307-t334*t81;
    t870 = t866*t869;
    t879 = t65*t68;
    t880 = complex<T>(1,0)/t866;
    t887 = t254*t149-t334*t72;
    t892 = t142*t149+t20*t100+t205*t72;
    t895 = t7*t732;
    t898 = t19*t111;
    t902 = -t330*t118-t17*(t647+t895)-t49*t898-t3*t167*t111;
    t903 = pow(t902,2);
    t909 = t137*t8+t142*t11+t858;
    t922 = -t330*t909-t17*(t142*t28+t165*t732)-t171*(t87*t732+t133*t111)-t20*(
t536+t165*t111);
    t923 = t192*t922;
    t924 = t2*t72;
    t938 = t134*t28;
    t942 = -t330*(t505*t8+t205*t11+t223*t76)-t17*(t205*t28+t223*t732)-t171*(
t134*t732+t138*t111)-t20*(t938+t223*t111);
    t943 = t924*t942;
    t944 = t330*t100*t887-t923-t943;
    t945 = complex<T>(1,0)/t944;
    t953 = pow(t8,2);
    t954 = complex<T>(1,0)/t953;
    t955 = pow(t944,2);
    t956 = complex<T>(1,0)/t955;
    t963 = t137*t149+t17*t100+t505*t72;
    t969 = t19*t28;
    t973 = -t330*t53-t10*t5*t28-t49*t969-t171*(t895+t425);
    t978 = complex<T>(1,0)/(t887*(t331+t647)-t923-t943);
    t983 = t230+t227+t421;
    t986 = t5*t732;
    t992 = -t330*t78-t51*t986-t4*t167*t732-t20*(t647+t425);
    t995 = pow(t19,2);
    t999 = t212+t416;
    t1003 = pow(t999,2);
    t1015 = t254*t305-t334*t79;
    t1024 = complex<T>(1,0)/t8;
    t1030 = complex<T>(1,0)/t887;
    t1031 = t880*t1030;
    t1034 = pow(t892,3);
    t1045 = t2*t100;
    t1046 = complex<T>(23,0);
    t1055 = t708*t8;
    t1056 = t149*t922;
    t1058 = t72*t942;
    t1070 = t887*(t249+t41)-t1056-t1058;
    t1079 = t212+t209+t416;
    t1081 = pow(t1015,2);
    t1094 = pow(t887,2);
    t1099 = pow(t1079,2);
    t1105 = t331+t41;
    t1121 = t527*t8+t212+t209+t416;
    t1130 = pow(t1070,2);
    t1159 = SSS(1,2,3,4);
    t1160 = complex<T>(1,0)/t1159;
    t1166 = pow(t99,2);
    t1167 = t2*t1166;
    t1170 = pow(t307,2);
    t1173 = t26*t76+t986+t898;
    t1174 = complex<T>(1,0)/t1173;
    t1178 = t64*t76;
    t1195 = pow(t167,2);
    t1197 = t2*t1195*t67;
    t1198 = pow(t99,2);
    t1202 = SPB(2,5);
    t1207 = pow(t26*t149+t5*t1202+t19*t305+t167*t307,2);
    t1208 = complex<T>(1,0)/t909;
    t1224 = SPB(5,6);
    t1231 = pow(-t2*t30*t1224-t2*(t250*t11+t938)*t88,2);
    t1233 = t2*t1159;
    t1243 = pow(t67,2);
    t1267 = t4*t727*t111+t10*(t49*t133+t387)*t111;
    t1269 = pow(t10,2);
    t1270 = complex<T>(1,0)/t1269;
    t1274 = t362*t853;
    t1275 = SS(1,6,7);
    t1281 = t3*t53;
    t1291 = SSS(1,5,6,7);
    t1292 = complex<T>(1,0)/t1291;
    t1318 = pow(t7,2);
    t1322 = pow(t18+t21+t171*t76,2);
    return(t1*(-t4*t6*t14+t23*(-t2*t24*t31/t23*t34-t40*t3/(t41+t42))*t54*t55)*
t60/t62*t65*t68-recursive*t2*(t73*(t74*(-t78*t83*t92*t60*t97*t98*t105+t74*(t1*
t3*t7*t111*t60*t96+t1*t7*t118*t111*t121*t54)*t183*t186*t188*t198/(t199*t129*
t131+t2*(-t11*(-t142*t190-t87*t138*t88-t205*t129-t2*(t207*t147+t137*t209*t88+
t212*t129+t209*t129)*t131)-t76*t236)*t180))*t186*t54*t265+(t1*t11*t111*t277*
t286*t121*t290*t131*t34*t300-t2*(-t304*t310*t311*t60*t95*t65*t98*t318*t320*t265
-t1*(t2*t309*t129*t259*t98+t137*t83*t88*(t331+t2*t250*t72+t335)*t98)*t60*t96*
t320*t265)*t348-t2*(t74*(-t137*t111*t88*t92*t197-t2*t111*t129*t92*t197)*t60*
t348*t98*t362*t364-t2*(t367*((-t373*t53*t111-t376*t7*t10*t104*t111)*t236+t7*t53
*t111*(t51*t87*t102*t88+t387*t252+t10*t102*t129-t2*(-t10*t87*t392-t4*t64*t392-
t51*t101*t129)*t131))*t406*t290*t131*t364-t2*(t74*t412*(-t187*t111*(-t420+t423*
(-t371*t307-t4*t427))*t432-t24*t307*t299*((-t373*t422*t111-t438*t441)*t299+t446
*t447*t456)-t3*t104*t197*((t438*t307*t427-t111*(-t368*t7*t419+t423*(-t120*t10*
t468-t445*t133*t307-t4*t441)))*t299-t446*t447*t307*t456))*t406*t98*t318*t300*
t198/t521+t526*(-t527*t528*t8*t277*t286*t131+t2*t53*(t4*t10*t286*t539-t277*(
t542*t133*t28+t527*t528*t28-t2*(-t415*t547*t551-t49*t554*t283-t4*t547*t283)*
t131)*t131))*t406*t290*t34*t300)*t95)*t65)*t96)*t105-t2*(((-t74*t87*t412*t60*
t98*t318-t304*(t589*t590*t98+t376*t87*t88*t129*t98+t412*t98)*t60*t265)*t95*t348
-((-t367*t605*t60*t607-t526*t286*t539*t60*t300)*t54+t1*(t4*t528*t28*t624*t627+
t3*t10*t631*t283*t259*t264*t627-t111*t623*(t527*t7*t283*(t102*t448-t2*(-t254*
t450-t199*t259)*t131)*t648-t7*t626*(-t4*t10*t28*t259-t655*t551*t264)+t2*t283*(
t49*t133*t661+t664*t10*t536*t626+t4*t64*t264*t626))*t98)*t406*t95*t96*t318*t265
/t681)*t34)*t105+t1*t412*(t2*(-t691*t87*t422*t694/(t695*t310-t4*t309*t299)-t2*(
t703*t87*t307*t705-t2*(-t708*t133-t120*t710*t305*t712)*t96)*t95)*t105-t2*(-t724
*(t725*t299+t727*t740)*t744*t188/t746+(-t2*(t120*(t725*t750-t51*t64*t104*t305*
t740)*t188*t300*t712+t368*(t655*t417*t447*t310-t4*t7*t104*t307*t441*t299-t111*
t309*(t695*(-t420-t3*t422*t427)+t664*t418*t299))*t105*t300/t797)*t95+t120*(-
t724*t804*t300+(t710*t310*t744*t712+t7*t422*t694*(t492+t495+t498+t786+t789+t792
+t795+t510)*(-t133*t104*t309-t4*(t133*t26*t491+t207*t494+t135*t497+t554*t297+
t376*t104*t456))*t95*t300/t827)*t105))*t705)*t96)*t844*t98*t318)*t65)*t853+t74*
t31*(t856+t703*t67*t860*t590*t607/(t870+t334*t277))*t844*t879*t880*t34*t55+t1*
t887*(t2*t892*t903*t97*t880*t945+t887*(t120*t76*t892*t732*t903*t744*t954*t956+
t120*(-t120*t963*t887*t111*t973*t869*t978+t902*(t120*t8*t983*t111*t992+(-t995*
t512*t28-t38*t19*t100*t999*t28-t1003*t28-t51*t167*t730*t732-t51*t76*(t230+t421)
*t732)*t1015-t369*t8*t963*t28*t869)*t95*t96*t1024*t34)*t1031*t945-t2*(-t368*
t1034*t28*t903*t744*t954*t880*t956-t2*(-t892*t903*(t887*(-t368*t999*t28-t1045*(
t1046*t26*t8+t120*t19*t28))+t1055*t1056+t1055*t1058)*t1024*t880*t956-t691*t7*
t111*t992*t869*t34/t1070)*t95*t1030-t120*(t120*t1079*t1081-t2*t887*t111*t973*(
t51*t887*t111+t983*t869)*t978)*t880/t1094*t34+t120*(-t1099*t903*t744*t1024*t945
+t7*t111*t973*(t887*t1105-t923-t943)*(t887*(t2*t999*t1105-t1045*(-t330*t209-t4*
(t376*t26*t8+t969)))-t149*t1121*t922-t72*t1121*t942)*t869*t95*t34/t1130*t945-t7
*t869*(t10*t887*t111+t870)*t978)*t1031)*t96)/t703)*t60*t68/(t254*t100+t508*t72)
*t1160)-t1*(t1167*t67*t53*t31*t1170*t1174/(t53*t1173+t334*(-t165*t149*t81-t1178
*t88-t1045*(t167*t81+t99*t88)-t924*(t161*t81+t252)))*t55+(t631*(t1197-t2*t1198*
t53*t1207*t1208*t1174/t88)+t1197*t631*t42/(t2*t42+t1159)+t1167*t67*t1231*t55/(
t253+t1233))*t1160)*t60*t65/t1243*t54*t34+t74*(t542*t6*(-t142*t13-t165*t78)*t96
*t879*t14/(t87*t8+t554+t1178)-t7*t1267*t412*t1270*t65*t1208*t1274/t1275-t7*t903
*(t1281*t111+t10*t118*t111)*t96*t68*t54*t1208*t362*t1160*t1292+t2*(-t40*t7*t73*
t1267*t1270*t65*t1274/(t253+t2*t1275)+t118*(-t185*t362-t1281*t96*t180*t362)*
t111*(t542*t1198*t310*t1159/(t335+t1159)+t39*t1318*t1322*t1291/(t1233+t1291))*
t348*t68*t54*t1160*t1292)*t1208)*t60*t95/t1318);
  }

#else

return
(

-(-((complex<T>(0,-1)*(-((complex<T>(-1,0)*SPA(3,2)*pow(SPA(6,2),2))/(SPA(4,2)*SPB(1,2)+
SPA(4,3)*SPB(1,3)))+
(pow(complex<T>(-1,0)*SPA(6,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,3)*SPB(1,3),2)*(-((complex<T>(-1,0)*pow(SPA(3,2),2)*pow(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3),2))/(pow(complex<T>(-1,0)*SPA(6,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,3)*SPB(1,3),2)*SPB(2,3)))-
(complex<T>(-1,0)*pow(SPA(2,1),complex<T>(2,0))*SPA(3,2))/(complex<T>(-1,0)*SPA(3,2)*SPB(2,3)+
SS(1,2,3))))/((complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SS(1,2,3))))/(complex<T>(2,0)*pow(SPA(3,2),complex<T>(2,0))*SPA(5,4)*SPA(7,6)))+
recursive*complex<T>(-1,0)*((pow(SPB(1,7),2)*((complex<T>(0,1)*(-(((complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))*(complex<T>(-1,0)*SPA(3,2)*SPB(3,7)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,7))*(complex<T>(-1,0)*SPA(2,1)*SPB(1,7)+
SPA(5,2)*SPB(5,7)-
complex<T>(-1,0)*SPA(6,2)*SPB(6,7)))/(complex<T>(2,0)*SPA(3,2)*SPA(4,3)*pow(SPB(1,7),2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))))+
(complex<T>(0,1)*((complex<T>(0,-1)*SPA(3,2)*SPA(4,2)*SPB(3,4))/(complex<T>(2,0)*SPA(4,3))+
(complex<T>(0,-1)*SPA(4,2)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(3,4))/(complex<T>(3,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))))*pow((SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7)))/SPB(1,7)+
(complex<T>(-1,0)*(-(SPB(1,3)*(-(SPA(5,3)*SPA(7,2)*SPB(5,7))-
complex<T>(-1,0)*SPA(5,2)*SPA(7,3)*SPB(5,7)-
SPA(7,3)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))-
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,3)*SPA(6,2)*SPB(1,6)*SPB(5,7)+
SPA(5,2)*SPA(6,3)*SPB(1,6)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
complex<T>(-1,0)*SPA(6,3)*SPB(1,6)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7)))-
SPB(1,4)*(-(SPA(5,4)*SPA(7,2)*SPB(5,7))-
complex<T>(-1,0)*SPA(5,2)*SPA(7,4)*SPB(5,7)-
SPA(7,4)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))-
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,4)*SPA(6,2)*SPB(1,6)*SPB(5,7)+
SPA(5,2)*SPA(6,4)*SPB(1,6)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))))/(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4)),3))/(complex<T>(1,0)*pow(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*(SPA(7,2)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))*((complex<T>(-1,0)*SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7)))/SPB(1,7)+
(complex<T>(-1,0)*(-(SPB(1,3)*(-(complex<T>(-1,0)*SPA(5,3)*SPA(7,2)*SPB(5,7))-
SPA(5,2)*SPA(7,3)*SPB(5,7)-
complex<T>(-1,0)*SPA(7,3)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))-
(complex<T>(-1,0)*(SPA(5,3)*SPA(6,2)*SPB(1,6)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,2)*SPA(6,3)*SPB(1,6)*SPB(5,7)+
SPA(5,3)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
SPA(6,3)*SPB(1,6)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7)))-
SPB(1,4)*(-(complex<T>(-1,0)*SPA(5,4)*SPA(7,2)*SPB(5,7))-
SPA(5,2)*SPA(7,4)*SPB(5,7)-
complex<T>(-1,0)*SPA(7,4)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))-
(complex<T>(-1,0)*(SPA(5,4)*SPA(6,2)*SPB(1,6)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,2)*SPA(6,4)*SPB(1,6)*SPB(5,7)+
SPA(5,4)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
SPA(6,4)*SPB(1,6)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))))/(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))))))/(complex<T>(1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7)))+
((complex<T>(0,-1)*SPB(1,3)*SPB(3,4)*(-(complex<T>(-1,0)*SPB(1,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7)))-
complex<T>(-1,0)*SPB(1,3)*(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7)))*(SPA(4,2)*SPB(2,3)+
(SPB(1,3)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7)))/SPB(1,7)))/(complex<T>(3,0)*pow(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3),2)*SPB(1,7)*SPB(2,3)*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))))-
(complex<T>(-1,0)*(-((complex<T>(0,-1)*SPA(5,2)*pow(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5),2)*pow(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7),2))/(complex<T>(2,0)*SPA(3,2)*SPA(5,4)*pow(SPB(1,7),2)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))*(SPA(7,5)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))))-
(complex<T>(0,-1)*((complex<T>(-1,0)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5))*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))/pow(SPB(1,7),2)+
(complex<T>(-1,0)*SPA(5,2)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,7)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,7))*SPB(5,7)*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,1)*SPB(1,7)+
complex<T>(-1,0)*SPA(7,6)*SPB(6,7)))/pow(SPB(1,7),2)))/(complex<T>(2,0)*SPA(4,3)*(SPA(7,5)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7)))))/SPA(4,2)-
(complex<T>(-1,0)*((complex<T>(0,1)*(-(complex<T>(-1,0)*SPA(5,2)*SPB(3,4)*SPB(5,7)*(complex<T>(-1,0)*SPA(2,1)*SPB(1,7)+
SPA(5,2)*SPB(5,7)-
complex<T>(-1,0)*SPA(6,2)*SPB(6,7))*(SPA(7,2)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7)))-
complex<T>(-1,0)*SPB(3,4)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))*(complex<T>(-1,0)*SPA(2,1)*SPB(1,7)+
SPA(5,2)*SPB(5,7)-
complex<T>(-1,0)*SPA(6,2)*SPB(6,7))*(SPA(7,2)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))))/(complex<T>(2,0)*SPA(4,2)*pow(SPB(1,7),2)*SPB(3,4)*(SPA(7,2)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7)))-
(complex<T>(-1,0)*((complex<T>(0,-1)*SPB(5,7)*((-(complex<T>(-1,0)*(complex<T>(-3,0)*SPA(4,3)*SPA(5,2)+
SPA(4,2)*SPA(5,3))*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(3,4))-
complex<T>(-2,0)*SPA(4,2)*SPA(4,3)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(3,4))*(-(complex<T>(-1,0)*SPA(5,4)*SPA(7,2)*SPB(5,7))-
SPA(5,2)*SPA(7,4)*SPB(5,7)-
complex<T>(-1,0)*SPA(7,4)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))-
(complex<T>(-1,0)*(SPA(5,4)*SPA(6,2)*SPB(1,6)*SPB(5,7)+
complex<T>(-1,0)*SPA(5,2)*SPA(6,4)*SPB(1,6)*SPB(5,7)+
SPA(5,4)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
SPA(6,4)*SPB(1,6)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))+
SPA(4,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(3,4)*(complex<T>(-1,0)*SPA(4,3)*SPA(5,2)*SPA(7,5)*SPB(5,7)+
SPA(3,2)*SPA(5,4)*SPA(7,5)*SPB(5,7)+
SPA(4,3)*SPA(7,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))-
(complex<T>(-1,0)*(-(SPA(4,3)*SPA(5,2)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPA(3,2)*SPA(5,4)*SPA(6,5)*SPB(1,6)*SPB(5,7)-
complex<T>(-1,0)*SPA(4,3)*SPA(6,5)*SPB(1,6)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))))/(complex<T>(6,0)*pow(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3),2)*SPB(1,7)*(SPA(7,2)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7)))-
(complex<T>(-1,0)*((complex<T>(0,1)*pow(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7),2)*(-(pow(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*SPB(3,4)*(-(complex<T>(3,0)*SPA(4,2)*SPA(5,4)*(SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7))*SPB(4,5))+
complex<T>(-1,0)*(SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*(-(SPA(4,2)*SPA(5,3)*SPB(4,5))-
complex<T>(-1,0)*SPA(3,2)*(SPA(4,3)*SPB(3,4)+
complex<T>(-1,0)*SPA(5,4)*SPB(4,5))))*pow(SPA(7,2)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7),2))-
pow(SPA(3,2),2)*SPB(4,5)*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))*((-(complex<T>(-1,0)*(complex<T>(-3,0)*SPA(4,3)*SPA(5,2)+
SPA(4,2)*SPA(5,3))*(SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*SPB(3,4))-
SPA(4,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(5,4)*SPB(4,5)))*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))+
complex<T>(2,0)*SPA(4,2)*(SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(3,4)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,2)*SPB(2,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,7))-
complex<T>(-1,0)*SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7)))-
SPA(3,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(SPA(7,2)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))*((SPA(4,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(4,5)*(SPA(4,3)*SPB(3,4)+
complex<T>(-1,0)*SPA(5,4)*SPB(4,5))-
SPB(3,4)*(-(complex<T>(-3,0)*SPA(4,2)*SPA(5,4)*(SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7))*SPB(4,5))+
complex<T>(-1,0)*(SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*(-(complex<T>(3,0)*SPA(4,3)*SPA(5,2)*SPB(4,5))-
complex<T>(2,0)*SPA(4,2)*SPA(5,3)*SPB(4,5)-
complex<T>(-1,0)*SPA(3,2)*(complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(5,4)*SPB(4,5)))))*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))-
complex<T>(2,0)*SPA(4,2)*(SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(3,4)*SPB(4,5)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,2)*SPB(2,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,7))-
complex<T>(-1,0)*SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7)))))/(complex<T>(6,0)*pow(SPB(1,7),2)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))*(SPA(7,2)*SPB(5,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))*pow(SPA(3,2)*SPA(6,1)*SPB(1,3)*SPB(1,6)+
SPA(3,2)*SPA(6,2)*SPB(1,6)*SPB(2,3)+
SPA(3,2)*SPA(7,2)*SPB(1,7)*SPB(2,3)+
complex<T>(-1,0)*SPA(6,5)*SPA(7,2)*SPB(1,6)*SPB(5,7)+
complex<T>(-1,0)*SPA(6,2)*SPA(7,5)*SPB(1,6)*SPB(5,7)+
complex<T>(-1,0)*SPA(7,2)*SPA(7,5)*SPB(1,7)*SPB(5,7)+
complex<T>(-1,0)*SPA(7,5)*SPB(1,5)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))+
SPA(3,2)*SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))-
(complex<T>(-1,0)*(-(SPA(6,2)*SPA(6,5)*pow(SPB(1,6),2)*SPB(5,7))-
SPA(6,5)*SPB(1,5)*SPB(1,6)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7),2))+
(complex<T>(0,-1)*SPB(3,4)*(-((complex<T>(-2,0)*SPA(3,2)*SPA(4,2)*SPA(5,4)*SPB(1,2)*(-(complex<T>(-1,0)*SPB(1,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7)))-
complex<T>(-1,0)*SPB(1,3)*(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7)))*(SPA(4,2)*SPB(2,3)+
(SPB(1,3)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7)))/SPB(1,7)))/SPB(1,7))+
complex<T>(-1,0)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*(complex<T>(-1,0)*SPA(3,2)*SPA(4,3)*(SPA(4,2)*SPB(2,3)+
(SPB(1,3)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7)))/SPB(1,7))*(SPA(5,2)*SPB(2,3)+
(SPB(1,3)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))/SPB(1,7))-
((-(complex<T>(-1,0)*SPB(1,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7)))-
complex<T>(-1,0)*SPB(1,3)*(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7)))*(pow(SPA(4,2),complex<T>(2,0))*SPA(5,3)*SPB(2,3)+
complex<T>(-2,0)*SPA(3,2)*SPA(4,2)*SPA(5,4)*SPB(2,3)-
(complex<T>(-1,0)*(-(complex<T>(3,0)*SPA(4,2)*SPA(5,4)*SPB(1,3)*(SPA(3,1)*SPB(1,7)-
SPA(6,3)*SPB(6,7)))-
complex<T>(-1,0)*SPA(4,2)*SPA(5,3)*SPB(1,3)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7))-
complex<T>(-1,0)*SPA(3,2)*SPA(5,4)*SPB(1,3)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7))))/SPB(1,7)))/SPB(1,7))))/(complex<T>(6,0)*pow(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3),2)*SPB(2,3)*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))))))/SPA(3,2)))/SPA(5,4)))/SPA(4,3))/(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))-
(complex<T>(-1,0)*(((-((complex<T>(0,1)*SPA(5,2)*pow(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7),2))/(complex<T>(2,0)*pow(SPB(1,7),2)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))))-
(complex<T>(0,-1)*SPA(5,2)*((pow(SPA(5,2),2)*pow(SPB(5,7),2))/pow(SPB(1,7),2)+
(complex<T>(-2,0)*SPA(5,2)*SPB(5,7)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7)))/pow(SPB(1,7),2)+
pow(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7),2)/pow(SPB(1,7),2)))/(complex<T>(2,0)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))))/(SPA(3,2)*SPA(4,2))-
((-((complex<T>(0,-1)*SPB(5,7)*pow(SPA(5,2)*SPB(2,3)+
(SPB(1,3)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))/SPB(1,7),2))/(complex<T>(2,0)*(-(complex<T>(-1,0)*SPB(1,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7)))-
complex<T>(-1,0)*SPB(1,3)*(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7)))))-
(complex<T>(0,-1)*SPB(3,4)*(SPA(4,2)*SPB(2,3)+
(SPB(1,3)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7)))/SPB(1,7))*(SPA(5,2)*SPB(2,3)+
(SPB(1,3)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))/SPB(1,7)))/(complex<T>(2,0)*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))))/(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))+
(complex<T>(0,-1)*((complex<T>(-1,0)*SPA(3,2)*SPA(4,2)*SPA(5,4)*SPB(2,3)*pow(-(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*SPB(4,7))-
complex<T>(-1,0)*SPB(4,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))-
complex<T>(-1,0)*SPB(4,7)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)),2)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7)))/pow(SPB(1,7),2)+
(SPA(3,2)*SPA(4,3)*pow(SPB(3,4),2)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7))*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7)))/pow(SPB(1,7),2)-
(SPB(3,4)*(-(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*SPB(4,7))-
complex<T>(-1,0)*SPB(4,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))-
complex<T>(-1,0)*SPB(4,7)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))*(complex<T>(-2,0)*SPA(3,2)*SPA(4,2)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7))*(SPA(7,5)*SPB(2,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(2,7))-
complex<T>(-1,0)*SPB(1,2)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(3,2)*SPB(2,3)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))-
SPA(4,2)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))*(-(complex<T>(-1,0)*SPA(3,2)*SPA(4,3)*SPB(2,3)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7)))-
complex<T>(-3,0)*SPA(5,4)*(SPA(3,1)*SPB(1,7)-
SPA(6,3)*SPB(6,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7)))+
complex<T>(-1,0)*(SPA(4,1)*SPB(1,7)-
SPA(6,4)*SPB(6,7))*(complex<T>(-1,0)*SPA(4,2)*SPA(5,3)*pow(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(3,2)*SPB(2,3)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7),2)+
complex<T>(-3,0)*SPA(3,2)*SPA(4,3)*SPA(5,2)*SPB(2,3)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))+
complex<T>(-1,0)*SPA(3,2)*SPA(5,4)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7)))))/pow(SPB(1,7),2)))/(complex<T>(6,0)*SPA(3,2)*SPA(4,3)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7))*pow(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3)+
SPA(7,5)*SPB(5,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,5)*(SPA(5,1)*SPB(1,7)-
SPA(6,5)*SPB(6,7))))/SPB(1,7),2)))/SPB(2,3))/(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))+
(complex<T>(0,-1)*pow(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7),2)*((complex<T>(-1,0)*(-((complex<T>(-9,0)*SPA(5,2)*(SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*SPB(3,4)*SPB(4,5))/(complex<T>(-1,0)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*pow(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5),2)-
complex<T>(-1,0)*SPA(3,2)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5))*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))))-
(complex<T>(-1,0)*((complex<T>(9,0)*SPA(5,2)*SPB(4,5))/(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5))-
(complex<T>(-1,0)*(-(complex<T>(-23,0)*SPA(5,3))-
(complex<T>(3,0)*pow(SPA(5,3),2)*SPB(3,5))/(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,2)*SPB(2,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,7))-
complex<T>(-1,0)*SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7))))/SPA(4,3)))/SPA(3,2)))/(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))-
(complex<T>(-1,0)*(-((complex<T>(-3,0)*SPA(5,3)*(pow(SPA(5,3),2)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5))*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))+
SPA(4,3)*SPA(5,4)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,4)+
SPA(5,2)*SPB(4,5))*(-(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)*SPB(1,6))-
complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,4)-
complex<T>(-1,0)*SPA(7,2)*SPB(1,7)*SPB(2,4)-
complex<T>(-1,0)*SPB(1,4)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))))/(pow(SPA(3,2),2)*pow(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*pow(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,2)*SPB(2,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,7))-
complex<T>(-1,0)*SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7),2)))+
(-((complex<T>(-1,0)*((complex<T>(3,0)*(pow(SPA(5,3),2)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5))*pow(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)),2)-
complex<T>(-1,0)*SPA(4,3)*SPA(5,4)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(3,5)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,4)+
SPA(5,2)*SPB(4,5))*(-(complex<T>(-1,0)*SPA(6,1)*SPB(1,4)*SPB(1,6))-
complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,4)-
complex<T>(-1,0)*SPA(7,2)*SPB(1,7)*SPB(2,4)-
complex<T>(-1,0)*SPB(1,4)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))))/(pow(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,2)*SPB(2,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,7))-
complex<T>(-1,0)*SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7)))+
(complex<T>(-3,0)*(complex<T>(-3,0)*SPA(5,4)*(SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7))*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(3,4)*pow(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5),2)-
complex<T>(-1,0)*SPA(3,2)*SPA(4,2)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*SPB(4,5)*(complex<T>(-1,0)*SPA(4,3)*SPB(3,4)+
SPA(5,4)*SPB(4,5))*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))-
SPB(3,4)*(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5))*(complex<T>(-1,0)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(-(complex<T>(3,0)*SPA(4,2)*SPA(5,4)*(SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7))*SPB(4,5))-
SPA(3,2)*(SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*(SPA(4,3)*SPB(3,4)+
complex<T>(-1,0)*SPA(5,4)*SPB(4,5)))+
complex<T>(-3,0)*SPA(3,2)*SPA(5,4)*(SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7))*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))))))/((SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))*(complex<T>(-1,0)*SPA(3,2)*SPA(6,1)*SPB(1,3)*SPB(1,6)+
complex<T>(-1,0)*SPA(3,2)*SPA(6,2)*SPB(1,6)*SPB(2,3)+
complex<T>(-1,0)*SPA(3,2)*SPA(7,2)*SPB(1,7)*SPB(2,3)+
SPA(3,2)*SPA(6,5)*SPB(1,6)*SPB(3,5)+
SPA(3,2)*SPA(7,5)*SPB(1,7)*SPB(3,5)+
SPA(4,2)*SPA(6,5)*SPB(1,6)*SPB(4,5)+
SPA(4,2)*SPA(7,5)*SPB(1,7)*SPB(4,5)+
complex<T>(-1,0)*SPA(3,2)*SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))))))/SPA(3,2))+
complex<T>(3,0)*(-((complex<T>(-3,0)*SPA(5,3)*pow(SPB(3,5),2))/(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))))+
((pow(SPA(5,3),2)*pow(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5),2))/(pow(SPA(3,2),2)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,2)*SPB(2,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,7))-
complex<T>(-1,0)*SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7)))+
(SPA(4,2)*(SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*SPB(3,4)*SPB(4,5)*(SPA(3,2)*SPA(6,1)*SPB(1,3)*SPB(1,6)+
SPA(3,2)*SPA(6,2)*SPB(1,6)*SPB(2,3)+
SPA(3,2)*SPA(7,2)*SPB(1,7)*SPB(2,3)+
SPA(3,2)*SPA(6,5)*SPB(1,6)*SPB(3,5)+
SPA(3,2)*SPA(7,5)*SPB(1,7)*SPB(3,5)+
SPA(4,2)*SPA(6,5)*SPB(1,6)*SPB(4,5)+
SPA(4,2)*SPA(7,5)*SPB(1,7)*SPB(4,5)+
SPA(3,2)*SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))*(-(SPA(5,3)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5)))-
complex<T>(-1,0)*SPA(3,2)*(SPA(5,3)*SPA(6,1)*SPB(1,3)*SPB(1,6)+
SPA(5,3)*SPA(6,2)*SPB(1,6)*SPB(2,3)+
SPA(5,3)*SPA(7,2)*SPB(1,7)*SPB(2,3)+
SPA(5,3)*SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7))+
complex<T>(-2,0)*(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,2)*SPB(2,7)+
SPA(7,6)*SPB(6,7)-
(complex<T>(-1,0)*(-(complex<T>(-1,0)*SPA(6,2)*SPB(1,6)*SPB(2,7))-
complex<T>(-1,0)*SPB(1,2)*(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7))))/SPB(1,7)))))/(SPA(3,2)*(-(SPA(6,1)*SPB(1,3)*SPB(1,6))-
SPA(6,2)*SPB(1,6)*SPB(2,3)-
SPA(7,2)*SPB(1,7)*SPB(2,3)-
SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))*pow(complex<T>(-1,0)*SPA(3,2)*SPA(6,1)*SPB(1,3)*SPB(1,6)+
complex<T>(-1,0)*SPA(3,2)*SPA(6,2)*SPB(1,6)*SPB(2,3)+
complex<T>(-1,0)*SPA(3,2)*SPA(7,2)*SPB(1,7)*SPB(2,3)+
SPA(3,2)*SPA(6,5)*SPB(1,6)*SPB(3,5)+
SPA(3,2)*SPA(7,5)*SPB(1,7)*SPB(3,5)+
SPA(4,2)*SPA(6,5)*SPB(1,6)*SPB(4,5)+
SPA(4,2)*SPA(7,5)*SPB(1,7)*SPB(4,5)+
complex<T>(-1,0)*SPA(3,2)*SPB(1,3)*(SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)),2)))/(SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7))))/(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5))))/SPA(4,3)))/(complex<T>(18,0)*pow(SPB(1,7),2)*(SPA(6,1)*SPB(1,6)+
SPA(7,1)*SPB(1,7)+
SPA(7,6)*SPB(6,7)))))/SPA(5,4)))/SPB(6,7)+
(complex<T>(0,1)*pow(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3),2)*(complex<T>(14,0)+
(complex<T>(9,0)*SPA(7,6)*pow(complex<T>(-1,0)*SPA(5,4)*SPB(1,4)+
SPA(6,5)*SPB(1,6)+
SPA(7,5)*SPB(1,7),2)*pow(SPB(5,7),2))/((-(complex<T>(-1,0)*SPB(1,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7)))-
complex<T>(-1,0)*SPB(1,3)*(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7)))*((complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))+
complex<T>(-1,0)*SPA(7,6)*(-(complex<T>(-1,0)*SPB(1,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,2)*SPB(5,7)))-
complex<T>(-1,0)*SPB(1,3)*(complex<T>(-1,0)*SPA(4,3)*SPB(4,7)+
complex<T>(-1,0)*SPA(5,3)*SPB(5,7)))))))/(complex<T>(18,0)*SPA(5,4)*SPA(7,6)*(complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*SPB(2,3)*SS(1,2,3))+
(complex<T>(0,-1)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*((complex<T>(-1,0)*(complex<T>(-1,0)*SPA(5,3)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,3)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,7))*pow(-(complex<T>(-1,0)*SPA(6,1)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(3,4)-
SPA(3,2)*SPA(6,4)*SPB(3,4),2))/(SPA(3,2)*SPA(4,3)*(complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4)))))+
((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*((complex<T>(3,0)*SPB(1,4)*(complex<T>(-1,0)*SPA(5,3)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,3)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,7))*SPB(2,4)*pow(-(complex<T>(-1,0)*SPA(6,1)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(3,4)-
SPA(3,2)*SPA(6,4)*SPB(3,4),2))/(pow(SPA(3,2),2)*pow(SPB(1,2),2)*pow(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))),2))+
(complex<T>(3,0)*(-((complex<T>(3,0)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,2)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,2)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*SPB(3,4)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3)))-
SPA(4,3)*SPA(6,2)*SPB(2,3)-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(2,3)-
complex<T>(-1,0)*SPA(6,4)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4)))*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7)))/((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)+
SPA(3,2)*SPB(2,3))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4)))))+
((-(complex<T>(-1,0)*SPA(6,1)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(3,4)-
SPA(3,2)*SPA(6,4)*SPB(3,4))*(complex<T>(3,0)*SPB(1,2)*(SPA(5,4)*SPB(1,5)+
SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*SPB(3,4)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4)))-
complex<T>(-1,0)*SPA(4,3)*SPA(6,2)*SPB(2,4)-
complex<T>(-1,0)*SPA(3,2)*SPA(6,4)*SPB(2,4)-
complex<T>(-1,0)*SPA(6,3)*(SPA(3,2)*SPB(2,3)+
SPA(4,3)*SPB(3,4)))+
(-(pow(SPA(6,3),2)*pow(SPB(1,6),2)*SPB(2,3))-
complex<T>(2,0)*SPA(6,3)*SPB(1,6)*(SPA(5,3)*SPB(1,5)+
SPA(7,3)*SPB(1,7))*SPB(2,3)-
pow(SPA(5,3)*SPB(1,5)+
SPA(7,3)*SPB(1,7),2)*SPB(2,3)-
complex<T>(-1,0)*SPA(4,3)*SPA(6,4)*SPB(1,4)*SPB(1,6)*SPB(2,4)-
complex<T>(-1,0)*SPA(4,3)*SPB(1,4)*(SPA(5,4)*SPB(1,5)+
SPA(7,4)*SPB(1,7))*SPB(2,4))*(complex<T>(-1,0)*SPA(6,5)*SPB(3,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(3,7))-
complex<T>(-3,0)*SPA(4,3)*SPB(1,2)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,2)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,2)*SPB(1,7))*SPB(2,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/(SPA(3,2)*SPA(4,3)*SPB(1,2)*SPB(2,3))))/((complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4)))))-
(complex<T>(-1,0)*(-((complex<T>(-3,0)*pow(complex<T>(-1,0)*SPA(5,3)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,3)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,7),3)*SPB(2,3)*pow(-(complex<T>(-1,0)*SPA(6,1)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(3,4)-
SPA(3,2)*SPA(6,4)*SPB(3,4),2))/(pow(SPA(3,2),2)*pow(SPB(1,2),2)*(complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*pow(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))),2)))-
(complex<T>(-1,0)*(-(((complex<T>(-1,0)*SPA(5,3)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,3)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,7))*pow(-(complex<T>(-1,0)*SPA(6,1)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(3,4)-
SPA(3,2)*SPA(6,4)*SPB(3,4),2)*((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(-(complex<T>(-3,0)*(SPA(5,3)*SPB(1,5)+
SPA(7,3)*SPB(1,7))*SPB(2,3))-
complex<T>(-1,0)*SPB(1,6)*(complex<T>(23,0)*SPA(6,1)*SPB(1,2)+
complex<T>(3,0)*SPA(6,3)*SPB(2,3)))+
complex<T>(-23,0)*SPB(1,2)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))+
complex<T>(-23,0)*SPB(1,2)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4)))))/(SPB(1,2)*(complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*pow(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))),2)))-
(complex<T>(-9,0)*SPA(4,2)*SPB(3,4)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4)))-
complex<T>(-1,0)*SPA(4,3)*SPA(6,2)*SPB(2,4)-
complex<T>(-1,0)*SPA(3,2)*SPA(6,4)*SPB(2,4)-
complex<T>(-1,0)*SPA(6,3)*(SPA(3,2)*SPB(2,3)+
SPA(4,3)*SPB(3,4)))*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7)))/(SPB(2,3)*((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(SPA(6,1)*SPB(1,6)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3))-
SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4)))))))/(SPA(3,2)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7)))-
(complex<T>(3,0)*(complex<T>(3,0)*(SPA(5,3)*SPB(1,5)+
SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7))*pow(complex<T>(-1,0)*SPA(6,5)*SPB(3,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(3,7),2)-
(complex<T>(-1,0)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*SPB(3,4)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3)))-
SPA(4,3)*SPA(6,2)*SPB(2,3)-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(2,3)-
complex<T>(-1,0)*SPA(6,4)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4)))*(complex<T>(-1,0)*SPA(4,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*SPB(3,4)+
(SPA(5,4)*SPB(1,5)+
SPA(6,4)*SPB(1,6)+
SPA(7,4)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)+
SPA(3,2)*SPB(2,3))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))))))/((complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*pow(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7),2)*SPB(2,3))+
(complex<T>(3,0)*(-((pow(SPA(5,3)*SPB(1,5)+
SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7),2)*pow(-(complex<T>(-1,0)*SPA(6,1)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(3,4)-
SPA(3,2)*SPA(6,4)*SPB(3,4),2))/(pow(SPA(3,2),2)*SPB(1,2)*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))))))+
(SPA(4,2)*SPB(3,4)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3)))-
SPA(4,3)*SPA(6,2)*SPB(2,3)-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(2,3)-
complex<T>(-1,0)*SPA(6,4)*(SPA(4,2)*SPB(2,4)+
SPA(4,3)*SPB(3,4)))*((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))))*((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(complex<T>(-1,0)*(SPA(5,3)*SPB(1,5)+
SPA(7,3)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3))-
complex<T>(-1,0)*SPB(1,6)*(-(complex<T>(-1,0)*SPA(6,1)*SPA(6,3)*SPB(1,6))-
complex<T>(-1,0)*SPA(3,2)*(complex<T>(-2,0)*SPA(6,1)*SPB(1,2)+
SPA(6,3)*SPB(2,3))))-
SPB(1,5)*(complex<T>(-2,0)*SPA(3,2)*SPB(1,2)+
SPA(5,3)*SPB(1,5)+
SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7))*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
SPB(1,7)*(complex<T>(-2,0)*SPA(3,2)*SPB(1,2)+
SPA(5,3)*SPB(1,5)+
SPA(6,3)*SPB(1,6)+
SPA(7,3)*SPB(1,7))*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))))*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7)))/(SPA(3,2)*SPB(2,3)*pow((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(SPA(6,1)*SPB(1,6)+
complex<T>(-1,0)*SPA(3,2)*SPB(2,3))-
SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))),2)*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4)))))-
(SPA(4,2)*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))*(SPA(4,3)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*SPB(3,4)+
(complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,5)*SPB(4,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(4,7))))/((complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,1)*SPB(1,6)+
SPA(3,2)*SPB(2,3))-
complex<T>(-1,0)*SPB(1,5)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(5,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(5,2)*SPB(2,4)+
SPA(5,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(5,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(3,4)))-
complex<T>(-1,0)*SPB(1,7)*(-(complex<T>(-1,0)*SPA(6,1)*(complex<T>(-1,0)*SPA(7,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(7,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(complex<T>(-1,0)*SPA(7,3)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(2,4))-
complex<T>(-1,0)*SPA(6,4)*(SPA(7,2)*SPB(2,4)+
SPA(7,3)*SPB(3,4))-
complex<T>(-1,0)*SPA(6,3)*(SPA(7,2)*SPB(2,3)+
complex<T>(-1,0)*SPA(7,4)*SPB(3,4))))))/((complex<T>(-1,0)*SPA(5,4)*SPB(1,5)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,4)*SPB(1,7))*(complex<T>(-1,0)*SPA(6,5)*SPB(1,5)-
complex<T>(-1,0)*SPA(7,6)*SPB(1,7)))))/SPA(4,3)))/complex<T>(9,0)))/(complex<T>(2,0)*SPA(7,6)*(complex<T>(-1,0)*SPA(6,5)*SPB(1,6)+
complex<T>(-1,0)*SPA(7,5)*SPB(1,7))*SSS(1,2,3,4)))+
(complex<T>(0,-1)*((complex<T>(-1,0)*pow(SPA(6,5),complex<T>(2,0))*SPA(7,6)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*pow(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3),2)*pow(SPB(4,5),2))/((SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4)+
SPA(6,3)*SPB(3,4))*((complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*(SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4)+
SPA(6,3)*SPB(3,4))+
complex<T>(-1,0)*SPA(7,6)*(-(complex<T>(-1,0)*SPA(5,4)*SPB(1,5)*SPB(4,7))-
SPA(5,4)*SPB(1,4)*SPB(5,7)-
complex<T>(-1,0)*SPB(1,6)*(SPA(6,4)*SPB(4,7)+
SPA(6,5)*SPB(5,7))-
complex<T>(-1,0)*SPB(1,7)*(SPA(7,4)*SPB(4,7)+
SPA(7,5)*SPB(5,7))))*SS(1,2,3))+
(pow(SPB(3,4),2)*(complex<T>(-1,0)*pow(SPA(6,4),complex<T>(2,0))*SPA(7,6)-
(complex<T>(-1,0)*pow(SPA(6,5),2)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*pow(SPA(6,1)*SPB(1,5)+
SPA(6,2)*SPB(2,5)+
SPA(6,3)*SPB(3,5)+
SPA(6,4)*SPB(4,5),2))/((complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*(SPA(6,1)*SPB(1,4)+
SPA(6,2)*SPB(2,4)+
SPA(6,3)*SPB(3,4))*SPB(5,7)))+
(complex<T>(-1,0)*pow(SPA(6,4),complex<T>(2,0))*SPA(7,6)*pow(SPB(3,4),2)*SS(1,2,3))/(complex<T>(-1,0)*SS(1,2,3)+
SSS(1,2,3,4))+
(complex<T>(-1,0)*pow(SPA(6,5),complex<T>(2,0))*SPA(7,6)*pow(-(complex<T>(-1,0)*(SPA(6,1)*SPB(1,3)+
SPA(6,2)*SPB(2,3))*SPB(5,6))-
complex<T>(-1,0)*(SPA(7,1)*SPB(1,3)+
SPA(7,2)*SPB(2,3))*SPB(5,7),2))/(SS(1,2,3)*(SPA(7,6)*SPB(6,7)+
complex<T>(-1,0)*SSS(1,2,3,4))))/SSS(1,2,3,4)))/(complex<T>(2,0)*SPA(5,4)*pow(SPA(7,6),2)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(2,3))-
(complex<T>(0,1)*((pow(SPA(4,2),complex<T>(2,0))*pow(SPA(6,2),2)*(-(complex<T>(-1,0)*SPA(5,3)*(SPA(4,2)*SPB(1,2)+
SPA(4,3)*SPB(1,3)))-
complex<T>(-1,0)*SPA(5,4)*(complex<T>(-1,0)*SPA(3,2)*SPB(1,2)+
SPA(4,3)*SPB(1,4))))/(SPA(4,3)*SPA(5,4)*SPA(7,6)*(SPA(4,2)*SPB(1,2)+
SPA(4,3)*SPB(1,3))*(SPA(5,2)*SPB(1,2)+
SPA(5,3)*SPB(1,3)+
SPA(5,4)*SPB(1,4)))-
(SPA(4,2)*(complex<T>(-1,0)*SPA(3,2)*SPA(4,3)*SPA(5,4)*SPB(3,4)+
SPA(4,3)*(complex<T>(-1,0)*SPA(4,2)*SPA(5,3)+
SPA(3,2)*SPA(5,4))*SPB(3,4))*pow(SPA(2,1)*SPB(1,7)-
SPA(6,2)*SPB(6,7),2))/(pow(SPA(4,3),complex<T>(2,0))*SPA(5,4)*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(3,4)*SPB(6,7)*SS(1,6,7))-
(SPA(4,2)*pow(-(complex<T>(-1,0)*SPA(6,1)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4)))-
complex<T>(-1,0)*SPA(6,2)*(SPA(3,2)*SPB(2,3)+
SPA(4,2)*SPB(2,4))-
complex<T>(-1,0)*SPA(4,2)*SPA(6,3)*SPB(3,4)-
SPA(3,2)*SPA(6,4)*SPB(3,4),2)*(SPA(3,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SPB(3,4)+
SPA(4,3)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(3,4)))/(SPA(4,3)*SPA(7,6)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))*SPB(3,4)*SSS(1,2,3,4)*SSS(1,5,6,7))+
(complex<T>(-1,0)*(-((complex<T>(-1,0)*pow(SPA(2,1),complex<T>(2,0))*SPA(4,2)*pow(SPB(1,7),2)*(complex<T>(-1,0)*SPA(3,2)*SPA(4,3)*SPA(5,4)*SPB(3,4)+
SPA(4,3)*(complex<T>(-1,0)*SPA(4,2)*SPA(5,3)+
SPA(3,2)*SPA(5,4))*SPB(3,4)))/(pow(SPA(4,3),complex<T>(2,0))*SPA(5,4)*SPB(3,4)*SPB(6,7)*(SPA(7,6)*SPB(6,7)+
complex<T>(-1,0)*SS(1,6,7))))+
((SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*(-(complex<T>(1,0)/SPB(3,4))-
(SPA(3,2)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3)))/(SPA(4,3)*(SPA(3,2)*SPB(1,3)+
SPA(4,2)*SPB(1,4))*SPB(3,4)))*SPB(3,4)*((pow(SPA(4,2),complex<T>(2,0))*pow(SPA(6,5),2)*pow(complex<T>(-1,0)*SPA(3,2)*SPB(3,5)+
complex<T>(-1,0)*SPA(4,2)*SPB(4,5),2)*SSS(1,2,3,4))/(complex<T>(-1,0)*SPA(7,6)*SPB(6,7)+
SSS(1,2,3,4))+
(pow(SPA(2,1),complex<T>(2,0))*pow(SPA(4,2),2)*pow(complex<T>(-1,0)*SPA(6,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(6,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(6,4)*SPB(1,4),2)*SSS(1,5,6,7))/(complex<T>(-1,0)*SSS(1,2,3,4)+
SSS(1,5,6,7))))/(SPA(4,2)*SPA(7,6)*(complex<T>(-1,0)*SPA(4,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(4,3)*SPB(1,3))*SSS(1,2,3,4)*SSS(1,5,6,7))))/(complex<T>(-1,0)*SPA(5,2)*SPB(1,2)+
complex<T>(-1,0)*SPA(5,3)*SPB(1,3)+
complex<T>(-1,0)*SPA(5,4)*SPB(1,4))))/(complex<T>(2,0)*SPA(3,2)*pow(SPA(4,2),2))
)
 )
;


#endif

    }






 // *************** table of switch values *************

#define _R_qppQbpQmqbmlmlbp_lc R2q2Q1g2l_2032873_lc
#define _R_qpQbpQmmqbmlmlbp_lc R2q2Q1g2l_2033817_lc
#define _R_qmmQbmQpqbplmlbp_lc R2q2Q1g2l_2037408_lc
#define _R_qmQbmQppqbplmlbp_lc R2q2Q1g2l_2038480_lc



 // *************** more macro definitions *************

#define _CASE_qppQbpQmqbmlmlbp_lc case 2032873 : \
          return &R2q2Q1g2l_2032873_lc
#define _CASE_qpQbpQmmqbmlmlbp_lc case 2033817 : \
          return &R2q2Q1g2l_2033817_lc
#define _CASE_qmmQbmQpqbplmlbp_lc case 2037408 : \
          return &R2q2Q1g2l_2037408_lc
#define _CASE_qmQbmQppqbplmlbp_lc case 2038480 : \
          return &R2q2Q1g2l_2038480_lc


 // *************** function definitions using macros *************

template <class T> complex<T> _R_qppQbpQmqbmlmlbp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g2l_qppQbpQmqbmlmlbp_lc(ep,mpc);}

template <class T> complex<T> _R_qpQbpQmmqbmlmlbp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g2l_qpQbpQmmqbmlmlbp_lc(ep,mpc);}

template <class T> complex<T> _R_qmmQbmQpqbplmlbp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g2l_qmmQbmQpqbplmlbp_lc(ep,mpc);}

template <class T> complex<T> _R_qmQbmQppqbplmlbp_lc(
        const eval_param<T>& ep,  const mass_param_coll& mpc){
          return R2q2Q1g2l_qmQbmQppqbplmlbp_lc(ep,mpc);}


 // *************** define pointers *************

template <class T> complex<T> ( *R2q2Q1g2l_L_Ptr_eval( int hc))
     (const eval_param<T>& ,const mass_param_coll&) {
       switch (hc) {
       _CASE_qppQbpQmqbmlmlbp_lc;
       _CASE_qpQbpQmmqbmlmlbp_lc;
       _CASE_qmmQbmQpqbplmlbp_lc;
       _CASE_qmQbmQppqbplmlbp_lc;

       default: return 0;
        }
 }


 // *************** definitions for template *************

template complex<R> ( *R2q2Q1g2l_L_Ptr_eval(int hc))
             (const eval_param<R>&,const mass_param_coll&);
template complex<RHP> ( *R2q2Q1g2l_L_Ptr_eval(int hc))
             (const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> ( *R2q2Q1g2l_L_Ptr_eval(int hc))
             (const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> ( *R2q2Q1g2l_L_Ptr_eval(int hc))
             (const eval_param<RGMP>&,const mass_param_coll&);
#endif



}

