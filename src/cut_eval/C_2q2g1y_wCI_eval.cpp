/*
*C_2q2g1y.cpp
*
* Created on 11/17, 2008
*      Author: Zvi's script
*/
//
//#include <vector>
//#include "integrals_ep.h"
//#include "cached_integral.h"
//
//using namespace std;
//using BH::CachedIntegral::Cached_Bubble_Integral_User;
//using BH::CachedIntegral::Cached_Triangle_Integral_User;
//using BH::CachedIntegral::Cached_Box_Integral_User;
//
//namespace BH  {
//
//class Index_Vector;
//
//namespace CachedIntegral {
//
//#define _VERBOSE 0
//
//#define SPA(i,j) ep.spa(i-1,j-1)
//#define SPB(i,j) ep.spb(i-1,j-1)
//#define S(i,j) ep.s(i-1,j-1)
//#define SS(i,j,k) ep.s(i-1,j-1,k-1)
//#define Complex complex<R>
//template<class T> static inline complex<T> square(complex<T> x)
//{return(x*x);}
//template<class T> static inline complex<T> cube(complex<T> x)
//{return(x*x*x);}
//
//template static inline complex<R> square(complex<R> x);
//template static inline complex<RHP> square(complex<RHP> x);
//template static inline complex<RVHP> square(complex<RVHP> x);
//
//template static inline complex<R> cube(complex<R> x);
//template static inline complex<RHP> cube(complex<RHP> x);
//template static inline complex<RVHP> cube(complex<RVHP> x);
//
//
//
//
//class C2q2g1y_qpppqmgap_L_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpppqmgap_L_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//
//
//
//class C2q2g1y_qppmqmgap_L_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qppmqmgap_L_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpmpqmgap_L_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpmpqmgap_L_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpmmqmgap_L_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpmmqmgap_L_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpgapqmpp_L_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpgapqmpp_L_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpgapqmpm_L_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpgapqmpm_L_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpgapqmmp_L_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpgapqmmp_L_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpgapqmmm_L_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpgapqmmm_L_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//
//
//
//
//
//
//class C2q2g1y_qppqmpgap_SLC_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qppqmpgap_SLC_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qppqmmgap_SLC_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qppqmmgap_SLC_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpmqmpgap_SLC_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpmqmpgap_SLC_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpmqmmgap_SLC_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpmqmmgap_SLC_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpqmppgap_SLC_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpqmppgap_SLC_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpqmpmgap_SLC_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpqmpmgap_SLC_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpqmmpgap_SLC_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpqmmpgap_SLC_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpqmmmgap_SLC_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpqmmmgap_SLC_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//
//
//
//
//
//class C2q2g1y_qpppqmgap_nf_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpppqmgap_nf_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qppmqmgap_nf_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qppmqmgap_nf_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpmpqmgap_nf_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpmpqmgap_nf_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//
//
//class C2q2g1y_qpmmqmgap_nf_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpmmqmgap_nf_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpgapqmpp_nf_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpgapqmpp_nf_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpgapqmpm_nf_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpgapqmpm_nf_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//class C2q2g1y_qpgapqmmp_nf_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpgapqmmp_nf_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//
//class C2q2g1y_qpgapqmmm_nf_wCI : public Cut_Part_wCI {
//public:
//	C2q2g1y_qpgapqmmm_nf_wCI(const std::vector<int>&);
//	SeriesC<R> eval(const eval_param<R>& ep,const T& mu);
//};
//
//
//C2q2g1y_qpppqmgap_L_wCI::C2q2g1y_qpppqmgap_L_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpppqmgap_L_wCI::eval(const eval_param<R>& ep,
//                 const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//C2q2g1y_qppmqmgap_L_wCI::C2q2g1y_qppmqmgap_L_wCI(const std::vector<int>& ind){
////	cout << "C2q2g1y_qppmqmgap constructor called with ";
////	copy(ind.begin(),ind.end(),ostream_iterator<int>(cout," "));
////	cout << endl;
//	vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//        vector<int> c15 = c51;
//        vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//      	 vector<int> c512;    c512.push_back(5-1);
//                         for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//        vector<int> c125 = c512;
//      	 vector<int> c412;    c412.push_back(4-1);
//                         for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//        vector<int> c124 = c412;
//        vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                           c134.push_back(4-1);
//        vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                           c235.push_back(5-1);
//
//     CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//     CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123));
//
//     CI_users.push_back(new Cached_Triangle_Integral_User(c1,c2,c345));
//     CI_users.push_back(new Cached_Triangle_Integral_User(c1,c5,c234));
//     CI_users.push_back(new Cached_Triangle_Integral_User(c2,c3,c145));
//     CI_users.push_back(new Cached_Triangle_Integral_User(c3,c4,c125));
//
//     CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c3,c45));
//     CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c4,c15));
//     CI_users.push_back(new Cached_Box_Integral_User(c3,c4,c5,c12));
//     CI_users.push_back(new Cached_Box_Integral_User(c5,c1,c2,c34));
//}
//
//
//
//SeriesC<R> C2q2g1y_qppmqmgap_L_wCI::eval(const eval_param<R>& ep,
//                 const T& mu){
//
////	cout << "C2q2g1y_qppmqmgap eval called with ";
////	copy(ind.begin(),ind.end(),ostream_iterator<int>(cout," "));
////	cout << endl;
//
// // #define TimeStamp "Mon 17 Nov 2008 23:49:47 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//complex<R> spa12 = SPA(1,2);
//complex<R> spa13 = SPA(1,3);
//complex<R> spa15 = SPA(1,5);
//complex<R> spa23 = SPA(2,3);
//complex<R> spa34 = SPA(3,4);
//complex<R> spa45 = SPA(4,5);
//complex<R> spa14 = SPA(1,4);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb23 = SPB(2,3);
//complex<R> spa24 = SPA(2,4);
//complex<R> spb25 = SPB(2,5);
//complex<R> spb15 = SPB(1,5);
//complex<R> spb45 = SPB(4,5);
//complex<R> s23 = -(spa23*spb23);
//complex<R> s45 = -(spa45*spb45);
//complex<R> s34 = S(3,4);
//complex<R> s12 = -(spa12*spb12);
//complex<R> s15 = -(spa15*spb15);
//complex<R> t1 = spa12*spa15;
//complex<R> t6 = square(spa34);
//complex<R> t10 = square(spa14);
//complex<R> t13 = spa13*spb12;
//complex<R> t21 = square(spb12);
//complex<R> t25 = spa14*spa34;
//complex<R> t27 = -(spb12*spb25);
//complex<R> d4 = (s12 + s15 - s34)*spa45; d4 = R(1)/d4;
//complex<R> d5 = spa12*spa24*spa45; d5 = R(1)/d5;
//complex<R> d6 = (s12 + s15 - s34)*spa12*spa45; d6 = R(1)/d6;
//complex<R> d9 = spa15*spa45*R(2); d9 = R(1)/d9;
//complex<R> d13 = spa23*spa45*R(2); d13 = R(1)/d13;
//complex<R> d14 = (s12 + s15 - s34)*spa45*R(2); d14 = R(1)/d14;
//complex<R> t8 = -t13;
//complex<R> t14 = s34*t6;
//complex<R> t19 = d6*spb25;
//complex<R> t22 = spa13*t10;
//complex<R> t32 = s15*t6;
//complex<R> d1 = (s23 - s45)*spa45*t1; d1 = R(1)/d1;
//complex<R> d2 = spa45*t1*square(s23 - s45)*R(2); d2 = R(1)/d2;
//complex<R> d3 = spa23*spa45*t1*R(2); d3 = R(1)/d3;
//complex<R> d7 = spa24*spa45*t1; d7 = R(1)/d7;
//complex<R> d8 = spa45*t1; d8 = R(1)/d8;
//complex<R> d10 = spa24*spa45*t1*R(2); d10 = R(1)/d10;
//complex<R> d11 = spa45*t1*R(2); d11 = R(1)/d11;
//complex<R> d12 = spa23*t1*R(2); d12 = R(1)/d12;
//complex<R> t12 = -t19;
//complex<R> t16 = d7*spa14;
//complex<R> t17 = d2*spa23;
//complex<R> t23 = (d10*s23*spa14 + d11*spa13*spb23)*t14;
//complex<R> t33 = t19*t32 + d5*spa14*spb15*t6;
//complex<R> t35 = d14*t27*t32 + d13*spb15*t6*t8;
//complex<R> t4 = -(t17*t21*t22) - d1*t13*t25*R(2) - d3*spa13*t6*R(3);
//complex<R> t5 = (d8*spa13*spb23 + s23*t16)*t6;
//complex<R> t20 = t17*t21*t22 + d1*t13*t25*R(2);
//complex<R> t29 = t14*(t12 + t16);
//complex<R> co1 = d4*t27*t6;
//complex<R> co2 = d9*spb23*t6*t8;
//complex<R> co3 = d12*spa13*spb45*t14;
//complex<R> co4 = complex<R>(0,1);
//SeriesC<R> result =
//	co4*
//	(
//		t20*(*CI_users[0]->get_value(ep,mu))
//		+t4*(*CI_users[1]->get_value(ep,mu))
//		+co1*(*CI_users[2]->get_value(ep,mu))
//		+t33*(*CI_users[3]->get_value(ep,mu))
//		+t5*(*CI_users[4]->get_value(ep,mu))
//		+t29*(*CI_users[5]->get_value(ep,mu))
//		+co2*(*CI_users[6]->get_value(ep,mu))
//		+t23*(*CI_users[7]->get_value(ep,mu))
//		+co3*(*CI_users[8]->get_value(ep,mu))
//		+t35*(*CI_users[9]->get_value(ep,mu))
//	)
//		;
////_PRINT(ind.get_permutation_code());
////_PRINT(result);
////CI_users[0]->print_info();
////cout << "index vector " ;
////copy(ind.begin(),ind.end(),ostream_iterator<int>(cout," "));
////cout << endl;
//return(result);
//}
//
//
//
//
//C2q2g1y_qpmpqmgap_L_wCI::C2q2g1y_qpmpqmgap_L_wCI(const vector<int>& ind){
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c12,c345));
//CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//CI_users.push_back(new Cached_Bubble_Integral_User(c34,c125));
//CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123) );
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c2,c345) );
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c3,c145));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c4,c125) );
//CI_users.push_back(new Cached_Triangle_Integral_User(c4,c5,c123) );
//CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c3,c45) );
//CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c4,c15) );
//CI_users.push_back(new Cached_Box_Integral_User(c3,c2,c1,c45) );
//CI_users.push_back(new Cached_Box_Integral_User(c3,c4,c5,c12) );
//CI_users.push_back(new Cached_Box_Integral_User(c5,c1,c2,c34));
//
//}
//
//
//
//
//
//SeriesC<R> C2q2g1y_qpmpqmgap_L_wCI::eval(const eval_param<R>& ep,const T& mu){
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpmpqmgap L");
//#endif
//
//
//
// // #define TimeStamp "Mon 17 Nov 2008 23:50:30 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//complex<R> spa12 = SPA(1,2);
//complex<R> spa13 = SPA(1,3);
//complex<R> spa14 = SPA(1,4);
//complex<R> spa15 = SPA(1,5);
//complex<R> spa23 = SPA(2,3);
//complex<R> spa45 = SPA(4,5);
//complex<R> spa24 = SPA(2,4);
//complex<R> spa34 = SPA(3,4);
//complex<R> spa25 = SPA(2,5);
//complex<R> spa35 = SPA(3,5);
//complex<R> spb13 = SPB(1,3);
//complex<R> spb15 = SPB(1,5);
//complex<R> spb23 = SPB(2,3);
//complex<R> spb34 = SPB(3,4);
//complex<R> spb35 = SPB(3,5);
//complex<R> spb45 = SPB(4,5);
//complex<R> s12 = S(1,2);
//complex<R> s23 = -(spa23*spb23);
//complex<R> s13 = -(spa13*spb13);
//complex<R> s45 = -(spa45*spb45);
//complex<R> s34 = -(spa34*spb34);
//complex<R> s35 = -(spa35*spb35);
//complex<R> t3 = spa15*spa45;
//complex<R> t4 = square(s12 - s45);
//complex<R> t18 = square(spa24);
//complex<R> t19 = square(spb13);
//complex<R> t20 = square(spb35);
//complex<R> t21 = square(spa14);
//complex<R> t24 = -(spa12*R(2));
//complex<R> t31 = cube(spb13);
//complex<R> t32 = cube(spb35);
//complex<R> t33 = cube(spa24);
//complex<R> t35 = -(spa25*R(2));
//complex<R> t37 = s12*spa14;
//complex<R> t47 = spa24*R(2);
//complex<R> t48 = spa23*spa45;
//complex<R> t54 = spa12*spa23;
//complex<R> t57 = spa24*spb35;
//complex<R> d7 = spa15*square(s12 - s34)*R(2); d7 = R(1)/d7;
//complex<R> d8 = (s12 - s34)*s35*spa15; d8 = R(1)/d8;
//complex<R> d9 = s35*(s12 - s45)*spa15; d9 = R(1)/d9;
//complex<R> d10 = (s12 - s34)*spa15*square(s35); d10 = R(1)/d10;
//complex<R> d11 = s35*spa15*square(s12 - s34)*R(2); d11 = R(1)/d11;
//complex<R> d13 = (s12 - s45)*spa15*square(s35); d13 = R(1)/d13;
//complex<R> d15 = spa15*spa34*spa35*R(2); d15 = R(1)/d15;
//complex<R> d16 = s12 - s45; d16 = R(1)/d16;
//complex<R> d21 = spa15*square(spa35); d21 = R(1)/d21;
//complex<R> d22 = spa15*spa34*spa35; d22 = R(1)/d22;
//complex<R> d26 = spa15*cube(spa35); d26 = R(1)/d26;
//complex<R> d28 = spa15*spa35; d28 = R(1)/d28;
//complex<R> d29 = spa15*cube(spa13); d29 = R(1)/d29;
//complex<R> d30 = spa15*square(spa13); d30 = R(1)/d30;
//complex<R> d31 = spa13*spa15*spa34; d31 = R(1)/d31;
//complex<R> d36 = spa15*cube(spa35)*R(2); d36 = R(1)/d36;
//complex<R> d37 = spa15*spa35*R(2); d37 = R(1)/d37;
//complex<R> d38 = spa15*spa23*R(2); d38 = R(1)/d38;
//complex<R> t22 = -t48;
//complex<R> t34 = -t54;
//complex<R> t49 = t21*t31;
//complex<R> t50 = spa25*t32;
//complex<R> t55 = spa14*t19;
//complex<R> t59 = spa12*t47;
//complex<R> t64 = spa25*t20;
//complex<R> t66 = spb23*t33;
//complex<R> t77 = spa24*t35;
//complex<R> d1 = spa23*spa34*t3; d1 = R(1)/d1;
//complex<R> d2 = t3*t4*R(2); d2 = R(1)/d2;
//complex<R> d3 = s13*(s12 - s45)*t3; d3 = R(1)/d3;
//complex<R> d4 = s13*t3*t4*R(2); d4 = R(1)/d4;
//complex<R> d5 = (s12 - s45)*t3*square(s13); d5 = R(1)/d5;
//complex<R> d6 = spa15*spa35*t4*R(2); d6 = R(1)/d6;
//complex<R> d12 = s35*spa15*t4*R(2); d12 = R(1)/d12;
//complex<R> d14 = spa23*t3; d14 = R(1)/d14;
//complex<R> d17 = s13*(s23 - s45)*t3; d17 = R(1)/d17;
//complex<R> d18 = s13*t3*square(s23 - s45)*R(2); d18 = R(1)/d18;
//complex<R> d19 = (s23 - s45)*t3*square(s13); d19 = R(1)/d19;
//complex<R> d20 = spa23*spa34*t3*R(2); d20 = R(1)/d20;
//complex<R> d23 = t3*cube(spa13); d23 = R(1)/d23;
//complex<R> d24 = t3*square(spa13); d24 = R(1)/d24;
//complex<R> d25 = spa13*spa34*t3; d25 = R(1)/d25;
//complex<R> d27 = spa34*t3; d27 = R(1)/d27;
//complex<R> d32 = spa34*t3*R(2); d32 = R(1)/d32;
//complex<R> d33 = t3*R(2); d33 = R(1)/d33;
//complex<R> d34 = t3*cube(spa13)*R(2); d34 = R(1)/d34;
//complex<R> d35 = spa13*spa34*t3*R(2); d35 = R(1)/d35;
//complex<R> d39 = spa34*t48*R(2); d39 = R(1)/d39;
//complex<R> t11 = s23*(d34*s12*t21*t34 + d24*spa12*spa24*t37 + d35*t18*t37);
//complex<R> t13 = -(d28*spb34*t18) + s34*spa25*(d26*t22 + d21*t47);
//complex<R> t15 = d10*t22*t50 + d11*t22*t50 - d7*spa24*t64 + d8*t20*t77;
//complex<R> t16 = d22*s45*t18 + d31*spa14*spb45*t18 + d26*s45*spa25*t22 + d29*spb45*t21*t34 + d21*s45*spa25*t47 + d30*spa14*spb45*t59;
//complex<R> t17 = d21*s34*s45*spa24*spa25 - d37*s45*spb34*t18 + d36*s34*s45*spa25*t22 - d38*spb34*spb45*t33;
//complex<R> t42 = d6*spb34;
//complex<R> t45 = d18*t34*t49 + d19*t49*t54 + d17*t55*t59;
//complex<R> t62 = d14*spb13;
//complex<R> t69 = spa24*t55;
//complex<R> t75 = d23*t21;
//complex<R> t1 = d4*t34*t49 + d10*t48*t50 + d11*t48*t50 + d12*t48*t50 + d13*t48*t50 + d5*t49*t54 + d15*d16*t22*t57 + t22*t42*t57 + d3*t55*t59 + d7*spa24*t64 + d8*t47*t64 + d9*t47*t64 - d2*spa12*t69 - d1*t33*R(2) + d16*spa12*t18*t62*R(2);
//complex<R> t2 = d20*t33 + d19*t34*t49 + d5*t34*t49 + d12*t22*t50 + d13*t22*t50 + d18*t49*t54 + d4*t49*t54 + d15*d16*t48*t57 + t42*t48*t57 + d16*t18*t24*t62 + d2*spa12*t69 + d17*t24*t69 + d3*t24*t69 + d9*t20*t77;
//complex<R> t12 = d25*s23*spa14*t18 + d24*s23*spa14*t59 + d27*t66 + s23*t34*t75;
//complex<R> t46 = -(d22*s12*t18) + d25*t18*t37 + d26*s12*spa25*t48 + d24*t37*t59 + s12*t34*t75 + d21*s12*t77;
//complex<R> co1 = d32*s12*t66;
//complex<R> co2 = -(d33*spb34*t66);
//complex<R> co3 = d39*s12*spb15*t33;
//complex<R> co4 = Complex(0,1);
//SeriesC<R> result = co4*(
//		t1*(*CI_users[0]->get_value(ep,mu))
//		+ t45*(*CI_users[1]->get_value(ep,mu))
//		+ t15*(*CI_users[2]->get_value(ep,mu))
//		+ t2*(*CI_users[3]->get_value(ep,mu))
//		+ t46*(*CI_users[4]->get_value(ep,mu))
//		+ t12*(*CI_users[5]->get_value(ep,mu))
//		+ t13*(*CI_users[6]->get_value(ep,mu))
//		+ t16*(*CI_users[7]->get_value(ep,mu))
//		+ co1*(*CI_users[8]->get_value(ep,mu))
//		+ co2*(*CI_users[9]->get_value(ep,mu))
//		+ t11*(*CI_users[10]->get_value(ep,mu))
//		+ t17*(*CI_users[11]->get_value(ep,mu))
//		+ co3*(*CI_users[12]->get_value(ep,mu))
//		);
// return(result);
//}
//
//
//C2q2g1y_qpmmqmgap_L_wCI::C2q2g1y_qpmmqmgap_L_wCI(const vector<int>& ind){
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpmmqmgap L");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Mon 17 Nov 2008 23:50:31 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//CI_users.push_back(new Cached_Bubble_Integral_User(c12,c345));
//CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c3,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c1,c5,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c4,c15));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c4,c5,c12));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c5,c1,c23));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c1,c2,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c4,c3,c12));
//}
//
//
//
//SeriesC<R> C2q2g1y_qpmmqmgap_L_wCI::eval(const eval_param<R>& ep,const T& mu){
////{{qp, m, m, qm, gap}, L}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpmmqmgap L");
//#endif
//
//
//
// // #define TimeStamp "Mon 17 Nov 2008 23:50:31 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa23 = SPA(2,3);
//complex<R> spa34 = SPA(3,4);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb15 = SPB(1,5);
//complex<R> spb23 = SPB(2,3);
//complex<R> spa12 = SPA(1,2);
//complex<R> spb34 = SPB(3,4);
//complex<R> spb35 = SPB(3,5);
//complex<R> s45 = S(4,5);
//complex<R> s15 = S(1,5);
//complex<R> s12 = -(spa12*spb12);
//complex<R> t5 = square(spb15);
//complex<R> t7 = square(spa23);
//complex<R> t8 = square(spb35);
//complex<R> t12 = spa23*spb15;
//complex<R> d1 = (s12 - s45)*spb23*spb34; d1 = R(1)/d1;
//complex<R> d2 = spb23*spb34*square(s12 - s45)*R(2); d2 = R(1)/d2;
//complex<R> d3 = spb12*spb23*spb34*R(2); d3 = R(1)/d3;
//complex<R> d4 = spb34*R(2); d4 = R(1)/d4;
//complex<R> d5 = spb23*spb34*R(2); d5 = R(1)/d5;
//complex<R> d6 = spb12*R(2); d6 = R(1)/d6;
//complex<R> d7 = spb12*spb23*R(2); d7 = R(1)/d7;
//complex<R> t3 = -(d2*spb12*t7*t8) - d1*spb35*t12*R(2);
//complex<R> t15 = d3*t5;
//complex<R> t4 = d2*spb12*t7*t8 + d1*spb35*t12*R(2) + t15*R(3);
//complex<R> co1 = d4*spa12*spa23*t5;
//complex<R> co2 = d5*s15*spa12*t5;
//complex<R> co3 = d6*spa23*spa34*t5;
//complex<R> co4 = -(d7*s45*spa34*t5);
//complex<R> co5 = s15*s45*t15;
//complex<R> co6 = -(d5*s15*spa12*t5);
//complex<R> co7 = d7*s45*spa34*t5;
//complex<R> co8 = -co2;
//complex<R> co9 = Complex(0,1);
//SeriesC<R> result = co9*(
//		t3*(*CI_users[0]->get_value(ep,mu))
//		+ t4*(*CI_users[1]->get_value(ep,mu))
//		+ co1*(*CI_users[2]->get_value(ep,mu))
//		+ co2*(*CI_users[3]->get_value(ep,mu))
//		+ co3*(*CI_users[4]->get_value(ep,mu))
//		+ co4*(*CI_users[5]->get_value(ep,mu))
//		+ co5*(*CI_users[6]->get_value(ep,mu))
//		+ co8*(*CI_users[7]->get_value(ep,mu))
//		+ co7*(*CI_users[8]->get_value(ep,mu))
//		);
// return(result);
//}
//
//C2q2g1y_qppqmpgap_SLC_wCI::C2q2g1y_qppqmpgap_SLC_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qppqmpgap_SLC_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//
//C2q2g1y_qppqmmgap_SLC_wCI::C2q2g1y_qppqmmgap_SLC_wCI(const vector<int>& ind){
//
////{{qp, p, qm, m, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qppqmmgap SLC");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Mon 17 Nov 2008 23:51:53 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c14,c235));
//CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//CI_users.push_back(new Cached_Bubble_Integral_User(c35,c124));
//CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c2,c345));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c4,c235));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c5,c234));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c3,c145));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c4,c125));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c5,c124));
//CI_users.push_back(new Cached_Triangle_Integral_User(c4,c5,c123));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c4,c5,c23));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c1,c5,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c5,c14));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c2,c1,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c3,c2,c14));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c4,c1,c23));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c4,c3,c12));
//
//
//}
//
//
//
//
//
//SeriesC<R> C2q2g1y_qppqmmgap_SLC_wCI::eval(const eval_param<R>& ep,const T& mu){
////{{qp, p, qm, m, gap}, SLC}
//
//
// // #define TimeStamp "Mon 17 Nov 2008 23:51:53 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa12 = SPA(1,2);
//complex<R> spa13 = SPA(1,3);
//complex<R> spa15 = SPA(1,5);
//complex<R> spa23 = SPA(2,3);
//complex<R> spa34 = SPA(3,4);
//complex<R> spa35 = SPA(3,5);
//complex<R> spa14 = SPA(1,4);
//complex<R> spb15 = SPB(1,5);
//complex<R> spa45 = SPA(4,5);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb23 = SPB(2,3);
//complex<R> spb25 = SPB(2,5);
//complex<R> spa24 = SPA(2,4);
//complex<R> spb35 = SPB(3,5);
//complex<R> spa25 = SPA(2,5);
//complex<R> spb45 = SPB(4,5);
//complex<R> s34 = S(3,4);
//complex<R> s45 = -(spa45*spb45);
//complex<R> s14 = S(1,4);
//complex<R> s23 = -(spa23*spb23);
//complex<R> s12 = -(spa12*spb12);
//complex<R> s15 = -(spa15*spb15);
//complex<R> s35 = -(spa35*spb35);
//complex<R> t5 = (s14 - s23)*spa12;
//complex<R> t6 = (s23 - s45)*spa23;
//complex<R> t8 = spa15*R(2);
//complex<R> t10 = square(spa34);
//complex<R> t11 = spa12*spa23;
//complex<R> t16 = -(spa14*R(3));
//complex<R> t23 = square(spb15);
//complex<R> t24 = square(spb25);
//complex<R> t25 = square(spa13);
//complex<R> t26 = s14 - s23 + s45;
//complex<R> t29 = -(spa14*spa34);
//complex<R> t31 = s14 - s35;
//complex<R> t34 = s14*s45;
//complex<R> t39 = spa12*R(2);
//complex<R> t40 = cube(spb15);
//complex<R> t41 = cube(spb25);
//complex<R> t43 = spa13*spa34;
//complex<R> t45 = s23*spa24;
//complex<R> t49 = spb12*spb23;
//complex<R> t54 = spa13*spb15;
//complex<R> t59 = spa14*spa45;
//complex<R> t60 = spa23*spa24;
//complex<R> t66 = s45*spa13;
//complex<R> t76 = spa23*spa45;
//complex<R> t105 = spa34*spb25;
//complex<R> d16 = s14 - s23; d16 = R(1)/d16;
//complex<R> d26 = spa15*spa23*spa35; d26 = R(1)/d26;
//complex<R> d27 = (s12 + s15 - s34)*spa15*spa23; d27 = R(1)/d27;
//complex<R> d34 = (s12 + s15 - s34)*spa23; d34 = R(1)/d34;
//complex<R> d43 = (s12 + s15 - s34)*spa23*R(2); d43 = R(1)/d43;
//complex<R> d46 = spa35*spa45*R(2); d46 = R(1)/d46;
//complex<R> t27 = -t76;
//complex<R> t28 = -s23 + t31;
//complex<R> t37 = -t49;
//complex<R> t42 = -t59;
//complex<R> t44 = t11*R(2);
//complex<R> t57 = t25*t40;
//complex<R> t68 = spa34*t24;
//complex<R> t83 = t23*t29;
//complex<R> t88 = spa14*t43;
//complex<R> t89 = spa45*t41;
//complex<R> t91 = t23*t34;
//complex<R> t92 = s35*t45;
//complex<R> t94 = spb25*t10;
//complex<R> t100 = spa14*t10;
//complex<R> t112 = s34*t10;
//complex<R> d3 = spa23*t26*t5; d3 = R(1)/d3;
//complex<R> d4 = spa23*t5*square(t26); d4 = R(1)/d4;
//complex<R> d6 = t39*square(t31); d6 = R(1)/d6;
//complex<R> d13 = spa35*t5*t76*R(2); d13 = R(1)/d13;
//complex<R> d14 = spa45*t11; d14 = R(1)/d14;
//complex<R> d15 = spa25*spa35*t39; d15 = R(1)/d15;
//complex<R> d17 = spa25*t39*square(s14 - s23); d17 = R(1)/d17;
//complex<R> d18 = spa35*spa45*t11; d18 = R(1)/d18;
//complex<R> d19 = spa12*spa15*t6; d19 = R(1)/d19;
//complex<R> d20 = spa12*t26*t6; d20 = R(1)/d20;
//complex<R> d21 = t11*t8*square(s23 - s45); d21 = R(1)/d21;
//complex<R> d22 = spa12*t6*square(t26); d22 = R(1)/d22;
//complex<R> d24 = spa12*spa45*t6*t8; d24 = R(1)/d24;
//complex<R> d25 = spa45*t11*t8; d25 = R(1)/d25;
//complex<R> d28 = spa35*t11*t26; d28 = R(1)/d28;
//complex<R> d29 = t11*square(t26); d29 = R(1)/d29;
//complex<R> d30 = t11*cube(t26); d30 = R(1)/d30;
//complex<R> d35 = spa12*spa35*t26; d35 = R(1)/d35;
//complex<R> d36 = spa12*square(t26); d36 = R(1)/d36;
//complex<R> d37 = spa12*cube(t26); d37 = R(1)/d37;
//complex<R> d38 = spa15*spa35*t11; d38 = R(1)/d38;
//complex<R> d40 = spa15*t11; d40 = R(1)/d40;
//complex<R> d41 = spa35*t11; d41 = R(1)/d41;
//complex<R> d45 = spa45*t8; d45 = R(1)/d45;
//complex<R> d51 = spa35*t11*t8; d51 = R(1)/d51;
//complex<R> t17 = -(d24*(spa14*spa23*spb12 + spa45*t54));
//complex<R> t19 = d13*(-(spa14*spa35*spb15) + spb25*t76);
//complex<R> t50 = d17*spb35;
//complex<R> t52 = d45*t100*t49 + d46*t37*cube(spa34);
//complex<R> t61 = -(d14*R(2));
//complex<R> t69 = t57*t59;
//complex<R> t71 = d14*R(2);
//complex<R> t75 = -t94;
//complex<R> t81 = d26*spa13*spb12*t10 + d27*s12*t94;
//complex<R> t82 = t27*t41;
//complex<R> t86 = d21*t25;
//complex<R> d1 = spa35*spa45*t44; d1 = R(1)/d1;
//complex<R> d2 = t44*square(s14 - s23); d2 = R(1)/d2;
//complex<R> d5 = t26*t44*square(s14 - s23); d5 = R(1)/d5;
//complex<R> d7 = t28*t5; d7 = R(1)/d7;
//complex<R> d8 = spa12*t28*t31; d8 = R(1)/d8;
//complex<R> d9 = t5*square(t28); d9 = R(1)/d9;
//complex<R> d10 = spa12*t31*square(t28); d10 = R(1)/d10;
//complex<R> d11 = t28*t39*square(s14 - s23); d11 = R(1)/d11;
//complex<R> d12 = t28*t39*square(t31); d12 = R(1)/d12;
//complex<R> d23 = t26*t44*square(s23 - s45); d23 = R(1)/d23;
//complex<R> d31 = spa12*spa35*t28; d31 = R(1)/d31;
//complex<R> d32 = spa12*square(t28); d32 = R(1)/d32;
//complex<R> d33 = spa12*cube(t28); d33 = R(1)/d33;
//complex<R> d39 = spa12*t28; d39 = R(1)/d39;
//complex<R> d42 = t44*square(t26); d42 = R(1)/d42;
//complex<R> d44 = t39*square(t28); d44 = R(1)/d44;
//complex<R> d47 = t39*cube(t28); d47 = R(1)/d47;
//complex<R> d48 = t28*t39; d48 = R(1)/d48;
//complex<R> d49 = spa35*t26*t44; d49 = R(1)/d49;
//complex<R> d50 = t44*cube(t26); d50 = R(1)/d50;
//complex<R> t1 = d15*d16*t105*t27 + t105*t27*t50 + d22*t42*t57 + d7*spa24*t68 + d23*t69 + d4*t69 + d5*t69 + d16*spb15*t100*t71 + d3*spa13*t83 + t23*t42*t86 + d2*t23*t88 + d20*t23*t88 + d11*t60*t89 + d9*t60*t89 + d19*spb15*t88*R(2) - d18*cube(spa34)*R(2) + spa14*spa34*t17*R(3) + t10*t19*R(3);
//complex<R> t2 = d4*t42*t57 + d5*t42*t57 + d16*spb15*t100*t61 + d6*spa24*t68 - d7*spa24*t68 - d8*spa24*t68 + d15*d16*t105*t76 + t105*t50*t76 + d10*spa24*t82 + d11*spa24*t82 + d12*spa24*t82 + d9*spa24*t82 + d2*spa13*t83 + d3*t23*t88 + d1*cube(spa34) - t10*t19*R(3);
//complex<R> t9 = spa34*t16*t17 + d23*t42*t57 + d22*t69 + d20*spa13*t83 + t23*t59*t86 - d19*spb15*t88*R(2) + d25*t100*R(3);
//complex<R> t18 = d50*t34*t69 + d29*t88*t91 + d49*t34*t54*square(spa34);
//complex<R> t20 = -(d6*spa24*t68) + d8*spa24*t68 + (d10 + d12)*t60*t89;
//complex<R> t21 = -(d32*t45*t68) + d37*spb23*t69 + d33*t45*t82 + d36*spa13*spb23*t83 - d31*s23*spb25*square(spa34) + d35*spb23*t54*square(spa34);
//complex<R> t22 = d30*s45*t69 + d29*t66*t83 - d41*spb45*cube(spa34) + d40*spa14*spb45*square(spa34) + d28*s45*t54*square(spa34) + d38*t66*square(spa34);
//complex<R> t48 = d32*s35;
//complex<R> t56 = s14*(d28*spa13*spb15*t10 + d32*spa24*t68 + d30*t69 + d29*spa13*t83 + d33*t60*t89 + d31*t94);
//complex<R> t93 = d38*spa13*t112 + d27*s34*t75;
//complex<R> t35 = -t48;
//complex<R> t74 = t45*t48*t68 + d47*t82*t92 + d48*s23*spb35*t94;
//complex<R> t67 = spa24*t35*t68 + d33*s35*spa24*t82 + d39*spb35*t94;
//complex<R> co1 = d34*spb15*t75;
//complex<R> co2 = d42*t16*t43*t91;
//complex<R> co3 = d43*s12*spb15*t75;
//complex<R> co4 = -(d44*t68*t92*R(3));
//complex<R> co5 = d51*t112*t66;
//complex<R> co6 = Complex(0,1);
//SeriesC<R> result = co6*(
//		t2*(*CI_users[0]->get_value(ep,mu))
//		+ t1*(*CI_users[1]->get_value(ep,mu))
//		+ t20*(*CI_users[2]->get_value(ep,mu))
//		+ t9*(*CI_users[3]->get_value(ep,mu))
//		+ t81*(*CI_users[4]->get_value(ep,mu))
//		+ t56*(*CI_users[5]->get_value(ep,mu))
//		+ co1*(*CI_users[6]->get_value(ep,mu))
//		+ t21*(*CI_users[7]->get_value(ep,mu))
//		+ t93*(*CI_users[8]->get_value(ep,mu))
//		+ t67*(*CI_users[9]->get_value(ep,mu))
//		+ t22*(*CI_users[10]->get_value(ep,mu))
//		+ co2*(*CI_users[11]->get_value(ep,mu))
//		+ co3*(*CI_users[12]->get_value(ep,mu))
//		+ co4*(*CI_users[13]->get_value(ep,mu))
//		+ t52*(*CI_users[14]->get_value(ep,mu))
//		+ t74*(*CI_users[15]->get_value(ep,mu))
//		+ t18*(*CI_users[16]->get_value(ep,mu))
//		+ co5*(*CI_users[17]->get_value(ep,mu))
//		);
// return(result);
//}
//
//C2q2g1y_qpmqmpgap_SLC_wCI::C2q2g1y_qpmqmpgap_SLC_wCI(const vector<int>& ind){
//
////{{qp, m, qm, p, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpmqmpgap SLC");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Mon 17 Nov 2008 23:51:55 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c12,c345) );
//CI_users.push_back(new Cached_Bubble_Integral_User(c34,c125) );
//CI_users.push_back(new Cached_Bubble_Integral_User(c35,c124));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c1,c4,c35) );
//CI_users.push_back(new Cached_Box_Integral_User(c2,c1,c5,c34) );
//CI_users.push_back(new Cached_Box_Integral_User(c4,c3,c2,c15) );
//CI_users.push_back(new Cached_Box_Integral_User(c5,c3,c2,c14));
//
//
//}
//
//
//
//SeriesC<R> C2q2g1y_qpmqmpgap_SLC_wCI::eval(const eval_param<R>& ep,const T& mu){
//
//
////{{qp, m, qm, p, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpmqmpgap SLC");
//#endif
//
//
//
// // #define TimeStamp "Mon 17 Nov 2008 23:51:55 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa23 = SPA(2,3);
//complex<R> spa35 = SPA(3,5);
//complex<R> spa45 = SPA(4,5);
//complex<R> spb14 = SPB(1,4);
//complex<R> spa34 = SPA(3,4);
//complex<R> spb15 = SPB(1,5);
//complex<R> spa15 = SPA(1,5);
//complex<R> spb34 = SPB(3,4);
//complex<R> spa14 = SPA(1,4);
//complex<R> spb35 = SPB(3,5);
//complex<R> spa25 = SPA(2,5);
//complex<R> spb45 = SPB(4,5);
//complex<R> spa12 = SPA(1,2);
//complex<R> spa24 = SPA(2,4);
//complex<R> s12 = S(1,2);
//complex<R> s23 = S(2,3);
//complex<R> s34 = -(spa34*spb34);
//complex<R> s35 = -(spa35*spb35);
//complex<R> t12 = square(spa23);
//complex<R> t13 = square(spb45);
//complex<R> t16 = square(spa24);
//complex<R> t17 = square(spa25);
//complex<R> t28 = spa24*spb45;
//complex<R> d1 = (s12 - s35)*spa14*spa45; d1 = R(1)/d1;
//complex<R> d2 = (s12 - s34)*spa15*spa45; d2 = R(1)/d2;
//complex<R> d3 = spa15*spa45*square(s12 - s34)*R(2); d3 = R(1)/d3;
//complex<R> d4 = spa14*spa45*square(s12 - s35)*R(2); d4 = R(1)/d4;
//complex<R> d5 = (s12 - s34)*spa15*spa34*spa45*R(2); d5 = R(1)/d5;
//complex<R> d6 = (s12 - s35)*spa14*spa35*spa45*R(2); d6 = R(1)/d6;
//complex<R> d7 = spa15*spa34*spa45*R(2); d7 = R(1)/d7;
//complex<R> d8 = spa14*spa35*spa45*R(2); d8 = R(1)/d8;
//complex<R> d9 = spa35*spa45*R(2); d9 = R(1)/d9;
//complex<R> d10 = spa34*spa45*R(2); d10 = R(1)/d10;
//complex<R> d11 = spa15*spa45*R(2); d11 = R(1)/d11;
//complex<R> d12 = spa14*spa45*R(2); d12 = R(1)/d12;
//complex<R> t6 = d4*spa35*t13*t16 + d1*spa23*t28*R(2) + d8*t12*R(3) + d6*spa23*(spa12*spa34*spb14 - spa35*t28)*R(3);
//complex<R> t24 = d2*spa25;
//complex<R> t4 = -(d3*spa34*t13*t17) + spa23*spb45*t24*R(2) - d5*spa23*(spa12*spa35*spb15 + spa25*spa34*spb45)*R(3) - d7*t12*R(3);
//complex<R> t5 = -(d4*spa35*t13*t16) + d3*spa34*t13*t17 - spa23*spb45*t24*R(2) - d1*spa23*t28*R(2) + d5*spa23*(spa12*spa35*spb15 + spa25*spa34*spb45)*R(3) + d6*(-(spa12*spa23*spa34*spb14*R(3)) + spa23*spa35*t28*R(3));
//complex<R> co1 = -(d9*s12*spb14*t12);
//complex<R> co2 = d10*s12*spb15*t12;
//complex<R> co3 = d11*s23*spb34*t12;
//complex<R> co4 = -(d12*s23*spb35*t12);
//complex<R> co5 = Complex(0,1);
//SeriesC<R> result = co5*(
//		t5*(*CI_users[0]->get_value(ep,mu))
//		+ t4*(*CI_users[1]->get_value(ep,mu))
//		+ t6*(*CI_users[2]->get_value(ep,mu))
//		+ co1*(*CI_users[3]->get_value(ep,mu))
//		+ co2*(*CI_users[4]->get_value(ep,mu))
//		+ co3*(*CI_users[5]->get_value(ep,mu))
//		+ co4*(*CI_users[6]->get_value(ep,mu))
//		);
// return(result);
//}
//
//
//
//C2q2g1y_qpmqmmgap_SLC_wCI::C2q2g1y_qpmqmmgap_SLC_wCI(const vector<int>& ind){
//
//
////{{qp, m, qm, m, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpmqmmgap SLC");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Mon 17 Nov 2008 23:53:08 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c12,c345));
//CI_users.push_back(new Cached_Bubble_Integral_User(c14,c235));
//CI_users.push_back(new Cached_Bubble_Integral_User(c35,c124));
//CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c2,c345));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c4,c235));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c5,c234));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c3,c145));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c4,c125));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c5,c124));
//CI_users.push_back(new Cached_Triangle_Integral_User(c4,c5,c123));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c3,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c1,c4,c35));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c4,c15));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c5,c4,c12));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c1,c2,c35));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c5,c1,c23));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c5,c3,c12));
//
//
//}
//
//
//SeriesC<R> C2q2g1y_qpmqmmgap_SLC_wCI::eval(const eval_param<R>& ep,const T& mu){
//
//
////{{qp, m, qm, m, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpmqmmgap SLC");
//#endif
//
//
//
// // #define TimeStamp "Mon 17 Nov 2008 23:53:08 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa34 = SPA(3,4);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb15 = SPB(1,5);
//complex<R> spb24 = SPB(2,4);
//complex<R> spa24 = SPA(2,4);
//complex<R> spb23 = SPB(2,3);
//complex<R> spb25 = SPB(2,5);
//complex<R> spa14 = SPA(1,4);
//complex<R> spb13 = SPB(1,3);
//complex<R> spb34 = SPB(3,4);
//complex<R> spa23 = SPA(2,3);
//complex<R> spb35 = SPB(3,5);
//complex<R> spa12 = SPA(1,2);
//complex<R> spb14 = SPB(1,4);
//complex<R> spb45 = SPB(4,5);
//complex<R> spa45 = SPA(4,5);
//complex<R> s23 = -(spa23*spb23);
//complex<R> s12 = -(spa12*spb12);
//complex<R> s14 = -(spa14*spb14);
//complex<R> s35 = S(3,5);
//complex<R> s15 = S(1,5);
//complex<R> s45 = -(spa45*spb45);
//complex<R> s34 = -(spa34*spb34);
//complex<R> t4 = spb12*spb23;
//complex<R> t6 = spb34*R(2);
//complex<R> t10 = square(spb15);
//complex<R> t16 = -(spb35*R(3));
//complex<R> t23 = square(spa24);
//complex<R> t24 = square(spa34);
//complex<R> t25 = square(spb13);
//complex<R> t28 = -(spb15*spb35);
//complex<R> t30 = s12 - s45;
//complex<R> t31 = s14 - s35;
//complex<R> t33 = s35*s45;
//complex<R> t40 = cube(spa24);
//complex<R> t41 = cube(spa34);
//complex<R> t43 = spb13*spb15;
//complex<R> t44 = spb23*R(2);
//complex<R> t46 = s12*spb25;
//complex<R> t49 = spa12*spa23;
//complex<R> t55 = s45*spb13;
//complex<R> t58 = spb35*spb45;
//complex<R> t68 = spb12*spb45;
//complex<R> t70 = -(spa34*R(2));
//complex<R> t77 = spa34*R(2);
//complex<R> t89 = spa24*spb15;
//complex<R> d19 = s12 - s35; d19 = R(1)/d19;
//complex<R> d28 = spb14*spb23*spb34; d28 = R(1)/d28;
//complex<R> d29 = spb23*square(spb34); d29 = R(1)/d29;
//complex<R> d31 = spb23*cube(spb34); d31 = R(1)/d31;
//complex<R> d34 = spb12*spb24*spb34; d34 = R(1)/d34;
//complex<R> d35 = (s15 - s23 + s45)*spb12*spb34; d35 = R(1)/d35;
//complex<R> d36 = spb12*spb24; d36 = R(1)/d36;
//complex<R> d42 = spb14*spb45*R(2); d42 = R(1)/d42;
//complex<R> d45 = spb12*spb24*R(2); d45 = R(1)/d45;
//complex<R> t11 = t4*R(2);
//complex<R> t27 = -t68;
//complex<R> t37 = -t49;
//complex<R> t42 = -t58;
//complex<R> t59 = t25*t41;
//complex<R> t67 = spb25*t40;
//complex<R> t69 = spb35*t43;
//complex<R> t75 = spb15*t23;
//complex<R> t76 = t25*t58;
//complex<R> t84 = -(spa14*t10);
//complex<R> t86 = spb13*t24;
//complex<R> t100 = s15*t10;
//complex<R> t104 = s14*t46;
//complex<R> t108 = s23*t10;
//complex<R> d1 = (s12 - s35)*spb23*(s12 + t31); d1 = R(1)/d1;
//complex<R> d3 = s34*(s12 - s35)*t4; d3 = R(1)/d3;
//complex<R> d4 = s34*t30*t4; d4 = R(1)/d4;
//complex<R> d5 = spb34*t30*t4; d5 = R(1)/d5;
//complex<R> d6 = spb14*spb45*t4; d6 = R(1)/d6;
//complex<R> d7 = spb24*t44*square(s12 - s35); d7 = R(1)/d7;
//complex<R> d8 = (s12 - s35)*spb23*square(s12 + t31); d8 = R(1)/d8;
//complex<R> d9 = (s12 + t31)*t44*square(s12 - s35); d9 = R(1)/d9;
//complex<R> d11 = (s12 - s35)*t4*square(s34); d11 = R(1)/d11;
//complex<R> d13 = t30*t4*square(s34); d13 = R(1)/d13;
//complex<R> d14 = t4*t6*square(t30); d14 = R(1)/d14;
//complex<R> d16 = spb45*t30*t4*t6; d16 = R(1)/d16;
//complex<R> d17 = spb45*t4; d17 = R(1)/d17;
//complex<R> d18 = spb14*spb24*t44; d18 = R(1)/d18;
//complex<R> d20 = t44*square(t31); d20 = R(1)/d20;
//complex<R> d21 = spb23*t31*(s12 + t31); d21 = R(1)/d21;
//complex<R> d22 = spb23*t31*square(s12 + t31); d22 = R(1)/d22;
//complex<R> d23 = (s12 + t31)*t44*square(t31); d23 = R(1)/d23;
//complex<R> d25 = spb45*t4*t6; d25 = R(1)/d25;
//complex<R> d26 = spb14*spb23*(s12 + t31); d26 = R(1)/d26;
//complex<R> d27 = spb23*square(s12 + t31); d27 = R(1)/d27;
//complex<R> d30 = spb23*cube(s12 + t31); d30 = R(1)/d30;
//complex<R> d32 = spb23*(s12 + t31); d32 = R(1)/d32;
//complex<R> d33 = (s15 - s23 + s45)*spb34*t4; d33 = R(1)/d33;
//complex<R> d37 = spb14*spb34*t4; d37 = R(1)/d37;
//complex<R> d38 = t4*square(spb34); d38 = R(1)/d38;
//complex<R> d39 = t4*cube(spb34); d39 = R(1)/d39;
//complex<R> d40 = spb14*t4; d40 = R(1)/d40;
//complex<R> d41 = spb34*t4; d41 = R(1)/d41;
//complex<R> d43 = spb45*t6; d43 = R(1)/d43;
//complex<R> d44 = t44*square(s12 + t31); d44 = R(1)/d44;
//complex<R> d47 = (s12 + t31)*t44; d47 = R(1)/d47;
//complex<R> d48 = t44*cube(s12 + t31); d48 = R(1)/d48;
//complex<R> d49 = (s15 - s23 + s45)*t4*t6; d49 = R(1)/d49;
//complex<R> d50 = spb14*t4*t6; d50 = R(1)/d50;
//complex<R> t18 = d16*(spa23*spb12*spb35 + spa34*spb13*spb45);
//complex<R> t19 = d22*t27*t67 + d23*t27*t67 + (-d20 + d21)*spb25*t75;
//complex<R> t20 = s35*(d38*spb13*t28 + d30*t27*t67 + d27*spb25*t75 + d39*t76 - d26*spa24*square(spb15) + d37*spb13*square(spb15));
//complex<R> t36 = d33*spa14;
//complex<R> t47 = d27*s14;
//complex<R> t48 = d7*spa14;
//complex<R> t51 = d30*spb12;
//complex<R> t53 = d43*spb35*t10*t49 + d42*t37*cube(spb15);
//complex<R> t80 = d35*spa14*spa23*spb13*t10 + d34*t108;
//complex<R> t101 = d17*spb35;
//complex<R> d2 = t11*square(s12 - s35); d2 = R(1)/d2;
//complex<R> d10 = s34*t11*square(s12 - s35); d10 = R(1)/d10;
//complex<R> d12 = s34*t11*square(t30); d12 = R(1)/d12;
//complex<R> d15 = (s12 - s35)*spb14*spb45*t11; d15 = R(1)/d15;
//complex<R> d24 = spb14*spb45*t11; d24 = R(1)/d24;
//complex<R> d46 = t11*square(spb34); d46 = R(1)/d46;
//complex<R> d51 = t11*cube(spb34); d51 = R(1)/d51;
//complex<R> t9 = d12*t58*t59 + d13*t58*t59 + d4*t24*t69 + d5*t69*t70 + d14*t24*t76 + d25*spb35*t10*R(3) + spb15*spb35*t18*R(3);
//complex<R> t17 = d15*(spa34*spb14*spb35 + spa24*t27);
//complex<R> t21 = d38*t28*t55 + d39*s45*t76 - d40*spa45*cube(spb15) + d41*spa45*spb35*square(spb15) + d37*t55*square(spb15) + t36*t55*square(spb15);
//complex<R> t22 = d29*spa12*spb13*t28 + spb45*t40*t46*t51 - d27*t46*t75 + d31*spa12*t76 + d26*s12*spa24*square(spb15) + d28*spa12*spb13*square(spb15);
//complex<R> t35 = -t47;
//complex<R> t81 = d48*t104*t40*t68 + t46*t47*t75 + d47*s12*spa24*t84;
//complex<R> t82 = t33*(d50*spb13*t10 + d38*t69 + d51*t76);
//complex<R> t88 = -(d34*t100) + spb13*t100*t36;
//complex<R> t1 = d10*t58*t59 + d11*t58*t59 + d22*t67*t68 + d23*t67*t68 + d8*t67*t68 + d9*t67*t68 + d3*t24*t69 - d1*spb25*t75 + d20*spb25*t75 - d21*spb25*t75 + d19*t10*t101*t77 + d2*t28*t86 + d18*d19*t27*t89 + t48*t68*t89 + d24*cube(spb15) - t10*t17*R(3);
//complex<R> t2 = spb15*t16*t18 + d14*t24*t25*t42 + d10*t42*t59 + d11*t42*t59 + d12*t42*t59 + d13*t42*t59 + d8*t27*t67 + d9*t27*t67 + d2*t24*t69 + d19*t10*t101*t70 + d1*spb25*t75 + d5*t69*t77 + d3*t28*t86 + d4*t28*t86 + t27*t48*t89 + d18*d19*t68*t89 - d6*cube(spb15)*R(2) + t10*t17*R(3);
//complex<R> t74 = s14*spb45*t51*t67 + spb25*t35*t75 + d32*spa24*t84;
//complex<R> co1 = -(d36*spa34*t10);
//complex<R> co2 = -(d44*t104*t75*R(3));
//complex<R> co3 = -(d45*spa34*t108);
//complex<R> co4 = d46*t16*t33*t43;
//complex<R> co5 = d49*spa14*t100*t55;
//complex<R> co6 = Complex(0,1);
//SeriesC<R> result = co6*(
//		t2*(*CI_users[0]->get_value(ep,mu))
//		+ t19*(*CI_users[1]->get_value(ep,mu))
//		+ t1*(*CI_users[2]->get_value(ep,mu))
//		+ t9*(*CI_users[3]->get_value(ep,mu))
//		+ t22*(*CI_users[4]->get_value(ep,mu))
//		+ t74*(*CI_users[5]->get_value(ep,mu))
//		+ t88*(*CI_users[6]->get_value(ep,mu))
//		+ t80*(*CI_users[7]->get_value(ep,mu))
//		+ co1*(*CI_users[8]->get_value(ep,mu))
//		+ t20*(*CI_users[9]->get_value(ep,mu))
//		+ t21*(*CI_users[10]->get_value(ep,mu))
//		+ t53*(*CI_users[11]->get_value(ep,mu))
//		+ co2*(*CI_users[12]->get_value(ep,mu))
//		+ co3*(*CI_users[13]->get_value(ep,mu))
//		+ co4*(*CI_users[14]->get_value(ep,mu))
//		+ t81*(*CI_users[15]->get_value(ep,mu))
//		+ co5*(*CI_users[16]->get_value(ep,mu))
//		+ t82*(*CI_users[17]->get_value(ep,mu)));
// return(result);
//}
//
//
//
//C2q2g1y_qpqmppgap_SLC_wCI::C2q2g1y_qpqmppgap_SLC_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpqmppgap_SLC_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//
//C2q2g1y_qpqmpmgap_SLC_wCI::C2q2g1y_qpqmpmgap_SLC_wCI(const vector<int>& ind){
////{{qp, qm, p, m, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpqmpmgap SLC");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Tue 18 Nov 2008 00:08:49 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c12,c345));
//CI_users.push_back(new Cached_Bubble_Integral_User(c14,c235));
//CI_users.push_back(new Cached_Bubble_Integral_User(c15,c234));
//CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//CI_users.push_back(new Cached_Bubble_Integral_User(c25,c134));
//CI_users.push_back(new Cached_Bubble_Integral_User(c34,c125));
//CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c2,c345));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c4,c235));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c5,c234));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c3,c145));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c5,c134));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c4,c125));
//CI_users.push_back(new Cached_Triangle_Integral_User(c4,c5,c123));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c3,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c2,c1,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c4,c1,c25));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c4,c5,c12));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c1,c2,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c2,c1,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c4,c1,c23));
//
//}
//
//SeriesC<R> C2q2g1y_qpqmpmgap_SLC_wCI::eval(const eval_param<R>& ep,const T& mu){
//
////{{qp, qm, p, m, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpqmpmgap SLC");
//#endif
//
//
// // #define TimeStamp "Tue 18 Nov 2008 00:08:49 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa12 = SPA(1,2);
//complex<R> spa24 = SPA(2,4);
//complex<R> spa35 = SPA(3,5);
//complex<R> spa14 = SPA(1,4);
//complex<R> spa15 = SPA(1,5);
//complex<R> spa25 = SPA(2,5);
//complex<R> spa23 = SPA(2,3);
//complex<R> spa13 = SPA(1,3);
//complex<R> spa34 = SPA(3,4);
//complex<R> spa45 = SPA(4,5);
//complex<R> spb13 = SPB(1,3);
//complex<R> spb15 = SPB(1,5);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb23 = SPB(2,3);
//complex<R> spb25 = SPB(2,5);
//complex<R> spb34 = SPB(3,4);
//complex<R> spb35 = SPB(3,5);
//complex<R> spb45 = SPB(4,5);
//complex<R> s34 = -(spa34*spb34);
//complex<R> s45 = -(spa45*spb45);
//complex<R> s14 = S(1,4);
//complex<R> s25 = -(spa25*spb25);
//complex<R> s23 = -(spa23*spb23);
//complex<R> s13 = -(spa13*spb13);
//complex<R> s12 = -(spa12*spb12);
//complex<R> s15 = -(spa15*spb15);
//complex<R> s35 = -(spa35*spb35);
//complex<R> t6 = spa23*R(2);
//complex<R> t10 = square(s12 - s45);
//complex<R> t12 = square(s12 + s15 - s34);
//complex<R> t18 = spa15*spa34;
//complex<R> t25 = square(spa24);
//complex<R> t26 = square(spb13);
//complex<R> t27 = square(spb35);
//complex<R> t28 = spa25*spa34;
//complex<R> t29 = spa12*spa23;
//complex<R> t30 = -(s12*spa14);
//complex<R> t31 = spa35*R(2);
//complex<R> t32 = spa15*spa45;
//complex<R> t45 = square(spb15);
//complex<R> t47 = cube(spb13);
//complex<R> t50 = s23 - s45;
//complex<R> t51 = square(spa45);
//complex<R> t55 = cube(spb35);
//complex<R> t56 = s25 - s34;
//complex<R> t61 = s14 - s23;
//complex<R> t62 = s14 - s25;
//complex<R> t69 = s14*s34;
//complex<R> t72 = square(s13);
//complex<R> t73 = square(spa14);
//complex<R> t74 = cube(spb15);
//complex<R> t75 = cube(spa24);
//complex<R> t76 = spa12*spa34;
//complex<R> t79 = -(spa25*R(3));
//complex<R> t81 = s12*s23;
//complex<R> t82 = spa15*R(2);
//complex<R> t83 = s25*spb25;
//complex<R> t85 = s14*s45;
//complex<R> t92 = spa14*spa24;
//complex<R> t93 = spa34*spa35;
//complex<R> t99 = square(spa34);
//complex<R> t100 = -(spa45*spb15);
//complex<R> t107 = spa12*spa14;
//complex<R> t108 = spa25*spa45;
//complex<R> t111 = spa23*spa35;
//complex<R> t166 = spa25*R(3);
//complex<R> t186 = spa13*spa25;
//complex<R> t246 = spa14*spa23;
//complex<R> d7 = s13*(s12 - s45)*spa35*spa45; d7 = R(1)/d7;
//complex<R> d17 = (s12 - s45)*spa15*spa23; d17 = R(1)/d17;
//complex<R> d18 = (s12 - s34)*spa12*spa35; d18 = R(1)/d18;
//complex<R> d19 = (s12 - s45)*spa12*spa35; d19 = R(1)/d19;
//complex<R> d20 = (s12 - s34)*s35*spa12; d20 = R(1)/d20;
//complex<R> d21 = s35*(s12 - s45)*spa12; d21 = R(1)/d21;
//complex<R> d22 = (s12 - s34)*s35*spa12*spa15; d22 = R(1)/d22;
//complex<R> d23 = s35*(s12 - s45)*spa12*spa15; d23 = R(1)/d23;
//complex<R> d29 = (s12 - s34)*spa15*square(s35); d29 = R(1)/d29;
//complex<R> d32 = (s12 - s45)*spa15*square(s35); d32 = R(1)/d32;
//complex<R> d34 = spa12*spa35; d34 = R(1)/d34;
//complex<R> d35 = s12 - s45; d35 = R(1)/d35;
//complex<R> d37 = s12 - s34; d37 = R(1)/d37;
//complex<R> d59 = spa12*spa35*spa45; d59 = R(1)/d59;
//complex<R> d81 = spa35*spa45*square(spa13); d81 = R(1)/d81;
//complex<R> d83 = spa35*spa45*cube(spa13); d83 = R(1)/d83;
//complex<R> d84 = spa15*cube(spa35); d84 = R(1)/d84;
//complex<R> d85 = square(spa35); d85 = R(1)/d85;
//complex<R> d86 = spa15*square(spa35); d86 = R(1)/d86;
//complex<R> d87 = spa15*spa23*square(spa35); d87 = R(1)/d87;
//complex<R> d100 = spa13*spa35*spa45; d100 = R(1)/d100;
//complex<R> d108 = spa12*square(spa35); d108 = R(1)/d108;
//complex<R> d109 = spa12*spa15*square(spa35); d109 = R(1)/d109;
//complex<R> d115 = spa15*square(spa13); d115 = R(1)/d115;
//complex<R> d116 = spa15*spa23*square(spa13); d116 = R(1)/d116;
//complex<R> d117 = spa15*cube(spa13); d117 = R(1)/d117;
//complex<R> d118 = spa35*square(spa13); d118 = R(1)/d118;
//complex<R> d120 = spa35*cube(spa13); d120 = R(1)/d120;
//complex<R> d129 = spa12*square(spa35)*R(2); d129 = R(1)/d129;
//complex<R> t14 = spa25*t56;
//complex<R> t16 = square(t62);
//complex<R> t46 = -t92;
//complex<R> t48 = -t76;
//complex<R> t52 = s12 + t56;
//complex<R> t53 = s34 + t62;
//complex<R> t58 = spa13*spa24 + t246*R(2);
//complex<R> t59 = -t83;
//complex<R> t60 = spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6;
//complex<R> t68 = -(spa14*spa25*spb15) - spb13*t76;
//complex<R> t78 = spa23*t26;
//complex<R> t84 = spa35*t50;
//complex<R> t94 = t26*(spa13*spa24 + spa14*t6);
//complex<R> t98 = -t107;
//complex<R> t112 = -(d34*R(2));
//complex<R> t127 = t107*t74;
//complex<R> t129 = -(spa45*t28);
//complex<R> t130 = d34*R(2);
//complex<R> t146 = t47*t99;
//complex<R> t147 = spa45*t55;
//complex<R> t154 = spa25*t45;
//complex<R> t155 = spa24*t73;
//complex<R> t156 = spa34*t47;
//complex<R> t162 = spa14*t29;
//complex<R> t163 = spb15*t25;
//complex<R> t164 = spa34*t26;
//complex<R> t165 = spa45*t27;
//complex<R> t176 = spa45*t45;
//complex<R> t177 = t47*t73;
//complex<R> t185 = spb13*t25;
//complex<R> t210 = t79*t92;
//complex<R> t224 = t75*t83;
//complex<R> t239 = t107*t99;
//complex<R> t261 = t47*t76;
//complex<R> d1 = spa12*spa35*t28; d1 = R(1)/d1;
//complex<R> d2 = spa12*spa45*t18*t6; d2 = R(1)/d2;
//complex<R> d3 = spa35*spa45*t29; d3 = R(1)/d3;
//complex<R> d4 = s13*(s12 - s45)*t32; d4 = R(1)/d4;
//complex<R> d5 = s13*(s12 - s45)*spa23*t32; d5 = R(1)/d5;
//complex<R> d6 = t10*t32*t6; d6 = R(1)/d6;
//complex<R> d8 = s13*t10*t32*R(2); d8 = R(1)/d8;
//complex<R> d9 = (s12 - s45)*t32*t72; d9 = R(1)/d9;
//complex<R> d10 = s13*spa45*t10*t31; d10 = R(1)/d10;
//complex<R> d11 = (s12 - s45)*spa35*spa45*t72; d11 = R(1)/d11;
//complex<R> d12 = (s12 - s34)*spa23*t18; d12 = R(1)/d12;
//complex<R> d16 = (s12 - s34)*(s12 + s15 - s34)*spa23*t18; d16 = R(1)/d16;
//complex<R> d24 = (s12 - s34)*s35*spa15*t29; d24 = R(1)/d24;
//complex<R> d25 = s35*(s12 - s45)*spa15*t29; d25 = R(1)/d25;
//complex<R> d26 = spa12*t31*square(s12 - s34); d26 = R(1)/d26;
//complex<R> d27 = t10*t82; d27 = R(1)/d27;
//complex<R> d28 = spa12*t10*t31; d28 = R(1)/d28;
//complex<R> d30 = s35*t82*square(s12 - s34); d30 = R(1)/d30;
//complex<R> d31 = s35*t10*t82; d31 = R(1)/d31;
//complex<R> d33 = spa45*t111; d33 = R(1)/d33;
//complex<R> d36 = spa35*t28; d36 = R(1)/d36;
//complex<R> d42 = spa35*t6*square(t61); d42 = R(1)/d42;
//complex<R> d43 = t111*t61*(s45 + t61); d43 = R(1)/d43;
//complex<R> d44 = t111*t61*square(s45 + t61); d44 = R(1)/d44;
//complex<R> d45 = spa35*t6*(s45 + t61)*square(t61); d45 = R(1)/d45;
//complex<R> d46 = (s15 - s34)*spa23*t18; d46 = R(1)/d46;
//complex<R> d47 = (s15 - s34)*(s12 + s15 - s34)*spa23*t18; d47 = R(1)/d47;
//complex<R> d48 = spa12*spa35*spa45*t6; d48 = R(1)/d48;
//complex<R> d49 = spa23*t32*t50; d49 = R(1)/d49;
//complex<R> d50 = t32*square(t50)*R(2); d50 = R(1)/d50;
//complex<R> d51 = s13*t32*t50; d51 = R(1)/d51;
//complex<R> d52 = s13*spa23*t32*t50; d52 = R(1)/d52;
//complex<R> d54 = s13*t32*square(t50)*R(2); d54 = R(1)/d54;
//complex<R> d55 = t32*t50*t72; d55 = R(1)/d55;
//complex<R> d56 = s13*spa45*t31*square(t50); d56 = R(1)/d56;
//complex<R> d58 = spa45*t31; d58 = R(1)/d58;
//complex<R> d60 = square(t50); d60 = R(1)/d60;
//complex<R> d62 = t111*t50*(s45 + t61); d62 = R(1)/d62;
//complex<R> d63 = t111*t50*square(s45 + t61); d63 = R(1)/d63;
//complex<R> d64 = spa35*t6*(s45 + t61)*square(t50); d64 = R(1)/d64;
//complex<R> d66 = spa12*t28*t31; d66 = R(1)/d66;
//complex<R> d75 = spa34*t31; d75 = R(1)/d75;
//complex<R> d76 = spa35*t76; d76 = R(1)/d76;
//complex<R> d77 = square(t56); d77 = R(1)/d77;
//complex<R> d78 = t32*square(spa13); d78 = R(1)/d78;
//complex<R> d79 = spa23*t32*square(spa13); d79 = R(1)/d79;
//complex<R> d80 = t32*cube(spa13); d80 = R(1)/d80;
//complex<R> d82 = spa13*spa45*t111; d82 = R(1)/d82;
//complex<R> d88 = spa23*spa45*t18; d88 = R(1)/d88;
//complex<R> d92 = spa23*t12*t18; d92 = R(1)/d92;
//complex<R> d96 = t111*(s45 + t61); d96 = R(1)/d96;
//complex<R> d97 = t111*square(s45 + t61); d97 = R(1)/d97;
//complex<R> d98 = t111*cube(s45 + t61); d98 = R(1)/d98;
//complex<R> d99 = spa23*spa34*t12; d99 = R(1)/d99;
//complex<R> d101 = spa35*(s45 + t61); d101 = R(1)/d101;
//complex<R> d102 = spa35*square(s45 + t61); d102 = R(1)/d102;
//complex<R> d103 = spa35*cube(s45 + t61); d103 = R(1)/d103;
//complex<R> d110 = spa15*t29*square(spa35); d110 = R(1)/d110;
//complex<R> d114 = spa15*spa23*t12; d114 = R(1)/d114;
//complex<R> d119 = spa13*t111; d119 = R(1)/d119;
//complex<R> d121 = t32*square(spa13)*R(2); d121 = R(1)/d121;
//complex<R> d122 = t32*cube(spa13)*R(2); d122 = R(1)/d122;
//complex<R> d123 = spa45*t31*square(spa13); d123 = R(1)/d123;
//complex<R> d124 = spa45*t31*cube(spa13); d124 = R(1)/d124;
//complex<R> d125 = spa13*spa45*t31; d125 = R(1)/d125;
//complex<R> d130 = spa12*t82*square(spa35); d130 = R(1)/d130;
//complex<R> d131 = spa12*spa15*t6*square(spa35); d131 = R(1)/d131;
//complex<R> d132 = t82*cube(spa35); d132 = R(1)/d132;
//complex<R> d133 = spa34*t12*t6; d133 = R(1)/d133;
//complex<R> d137 = spa35*t6*(s45 + t61); d137 = R(1)/d137;
//complex<R> d138 = spa35*t6*square(s45 + t61); d138 = R(1)/d138;
//complex<R> d139 = spa35*t6*cube(s45 + t61); d139 = R(1)/d139;
//complex<R> t8 = square(t53);
//complex<R> t19 = s34*s45*(d132*t129 + d130*t210 - d129*t25 + d131*spa24*(spa15*spa23*spa24 + spa24*t186 + spa25*t246*R(2)));
//complex<R> t34 = cube(t53);
//complex<R> t39 = d110*(spa15*spa23*spa24 + spa24*t186 + spa25*t246*R(2));
//complex<R> t80 = -t147;
//complex<R> t109 = -t163;
//complex<R> t110 = -t156;
//complex<R> t115 = d58*spb12;
//complex<R> t133 = d75*spb12;
//complex<R> t134 = -(d33*spb13);
//complex<R> t135 = d46*spb25;
//complex<R> t144 = -t155;
//complex<R> t153 = d125*s12*spa34*spb23*t25 + d124*t239*t81 + d123*spa34*t81*t92;
//complex<R> t159 = -(d78*R(3));
//complex<R> t169 = t108*t127;
//complex<R> t172 = spb35*t112;
//complex<R> t180 = spb35*t130;
//complex<R> t184 = t127*t51;
//complex<R> t193 = t74*t98;
//complex<R> t200 = t154*t46;
//complex<R> t231 = t48*t73;
//complex<R> t233 = d6*t26;
//complex<R> d13 = (s12 - s34)*t52*t93; d13 = R(1)/d13;
//complex<R> d14 = (s12 - s34)*t93*square(t52); d14 = R(1)/d14;
//complex<R> d15 = spa34*t31*t52*square(s12 - s34); d15 = R(1)/d15;
//complex<R> d38 = spa25*t16*t31; d38 = R(1)/d38;
//complex<R> d39 = spa25*spa35*t53*t62; d39 = R(1)/d39;
//complex<R> d41 = spa25*t16*t31*t53; d41 = R(1)/d41;
//complex<R> d53 = s13*spa45*t84; d53 = R(1)/d53;
//complex<R> d57 = spa45*t72*t84; d57 = R(1)/d57;
//complex<R> d61 = t6*t84; d61 = R(1)/d61;
//complex<R> d65 = spa12*spa45*t6*t84; d65 = R(1)/d65;
//complex<R> d67 = t14*t31; d67 = R(1)/d67;
//complex<R> d68 = spa35*t14*t53; d68 = R(1)/d68;
//complex<R> d70 = spa25*t31*t53*square(t56); d70 = R(1)/d70;
//complex<R> d71 = t52*t56*t93; d71 = R(1)/d71;
//complex<R> d72 = t56*t93*square(t52); d72 = R(1)/d72;
//complex<R> d73 = spa34*t31*t52*square(t56); d73 = R(1)/d73;
//complex<R> d74 = t14*t31*t76; d74 = R(1)/d74;
//complex<R> d89 = spa35*t28*t52; d89 = R(1)/d89;
//complex<R> d90 = t93*square(t52); d90 = R(1)/d90;
//complex<R> d91 = t93*cube(t52); d91 = R(1)/d91;
//complex<R> d93 = spa25*spa35*t53; d93 = R(1)/d93;
//complex<R> d104 = spa35*t53; d104 = R(1)/d104;
//complex<R> d107 = t52*t93; d107 = R(1)/d107;
//complex<R> d111 = spa25*spa35*t52; d111 = R(1)/d111;
//complex<R> d112 = spa35*square(t52); d112 = R(1)/d112;
//complex<R> d113 = spa35*cube(t52); d113 = R(1)/d113;
//complex<R> d126 = spa25*t31*t53; d126 = R(1)/d126;
//complex<R> d134 = spa34*t31*square(t52); d134 = R(1)/d134;
//complex<R> d135 = spa34*t31*cube(t52); d135 = R(1)/d135;
//complex<R> d136 = spa34*t31*t52; d136 = R(1)/d136;
//complex<R> t5 = d12*spa14*t109 + d10*t107*t146 + d27*spa24*t165 + d28*spa23*spa24*t165 + d16*t224 + d17*spb35*t25 - d18*spb35*t25 - d19*spb35*t25 + d37*(d36*t100 + t180)*t25 + d35*(d33*spa34*t185 + t180*t25) + d22*t210*t27 + d23*t210*t27 - d20*t25*t27 - d21*t25*t27 + d29*t147*t28 + d30*t147*t28 + d31*t147*t28 + d32*t147*t28 + d26*spa24*t27*t28 + d8*t177*t48 + d14*t193*t51 + d15*t193*t51 + d24*spa24*t27*(spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6) + d25*spa24*t27*(spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6) - d1*t75 - d3*t75 + d9*t177*t76 + t233*t46*t76 + d7*t164*t92 + d13*t176*t92 + d5*t26*t58*t92 + d11*t146*t98 - d4*t155*t26*R(3) + d2*spa14*t75*R(3);
//complex<R> t20 = d121*spa24*spb23*t30*t58 + d122*t231*t81 - d121*t155*t81*R(3);
//complex<R> t21 = spb23*(d101*t163 + d103*t169 + d102*t200 + d100*spa34*t25 + d78*t46*t58) + s23*(t155*t159 + d80*t231 + d83*t239 + d81*spa34*t92);
//complex<R> t22 = d86*spb12*t210 - d85*spb12*t25 + d90*spa24*t176*t30 - d33*spb12*t75 - d36*spb12*t75 + d88*spa14*spb12*t75 + s12*(t155*t159 + d89*spa45*t163 + d91*t184 + d80*t231 + d83*t239 - d82*spa34*t25 + d84*spa45*t28 + d92*t59*t75 + d81*spa34*t92 + d79*t58*t92) + d87*spa24*spb12*(spa15*spa23*spa24 + spa24*t186 + spa25*t246*R(2));
//complex<R> t24 = d84*s45*t129 + d96*s45*t163 + d98*s45*t169 + d97*s45*t200 + d109*s45*t210 + d117*spb45*t231 + d120*spb45*t239 - d108*s45*t25 - d119*spa34*spb45*t25 + s45*spa24*t39 + d118*spa34*spb45*t92 + d116*spb45*t58*t92 - d115*spb45*t155*R(3);
//complex<R> t41 = d65*(spa12*spa45*spb15 - spb13*t246);
//complex<R> t87 = -t115;
//complex<R> t141 = d135*s12*s25*t184 + d136*s12*spb25*t100*t25 + d134*s25*spa24*t176*t30;
//complex<R> t175 = (d137*t163 + d139*t169 + d138*t200)*t85;
//complex<R> t192 = d47*t224 + t135*t75;
//complex<R> d40 = spa25*spa35*t62*t8; d40 = R(1)/d40;
//complex<R> d69 = spa35*t14*t8; d69 = R(1)/d69;
//complex<R> d94 = spa25*spa35*t8; d94 = R(1)/d94;
//complex<R> d95 = spa25*spa35*t34; d95 = R(1)/d95;
//complex<R> d105 = spa35*t8; d105 = R(1)/d105;
//complex<R> d106 = spa35*t34; d106 = R(1)/d106;
//complex<R> d127 = spa25*t31*t8; d127 = R(1)/d127;
//complex<R> d128 = spa25*t31*t34; d128 = R(1)/d128;
//complex<R> t1 = d56*t107*t146 + d61*t163 + d44*t169 + d45*t169 + d64*t169 + d49*spa14*t185 + d63*t108*t193 + d43*t200 + d50*t144*t26 + d54*t177*t48 - d48*t75 + d55*t177*t76 + d59*d60*t144*t78 + d60*spa14*t185*t87 + d42*t154*t92 + d62*t154*t92 + d53*t164*t92 + d52*t26*t58*t92 + d57*t146*t98 - d51*t155*t26*R(3) - t25*t41*R(3);
//complex<R> t2 = d77*spa14*t109*t133 + d76*d77*t144*t154 + d69*t110*t162 + d40*t156*t162 + d41*t156*t162 + d70*t156*t162 - d67*t185 + d72*t193*t51 + d73*t193*t51 - d66*t75 + d38*t46*t78 + d68*t46*t78 + d71*t176*t92 + d39*t78*t92 - d74*t25*t68*R(3);
//complex<R> t3 = d61*t109 + d11*t107*t146 + d57*t107*t146 - d27*spa24*t165 - d28*spa23*spa24*t165 + d63*t169 - d49*spa14*t185 + d64*t108*t193 + d62*t200 - d17*spb35*t25 + d19*spb35*t25 + d35*(spa34*t134 + t172)*t25 + d50*t155*t26 + d21*t25*t27 + d53*t164*t46 + d7*t164*t46 + d55*t177*t48 + d9*t177*t48 + d31*t129*t55 + d32*t129*t55 + d5*t26*t46*t58 + d52*t26*t46*t58 - d25*spa24*t27*(spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6) + d54*t177*t76 + d8*t177*t76 + d60*(spa14*t115*t185 + d59*t155*t78) + d23*t166*t27*t92 + t233*t76*t92 + d10*t146*t98 + d56*t146*t98 + d4*t155*t26*R(3) + d51*t155*t26*R(3) + t25*t41*R(3);
//complex<R> t4 = d70*t110*t162 + d69*t156*t162 + d12*spa14*t163 + d77*(d76*t154*t155 + spa14*t133*t163) + d14*t184 + d15*t184 + d72*t184 + d73*t184 + d67*t185 + d18*spb35*t25 + d37*(d36*spa45*t163 + t172*t25) + d20*t25*t27 - d26*spa24*t27*t28 + d13*t176*t46 + d71*t176*t46 + d29*t129*t55 + d30*t129*t55 - d24*spa24*t27*(spa15*spa23*spa24 + spa24*t186 + spa14*spa25*t6) - t135*t75 + d16*t59*t75 + d47*t59*t75 + d22*t166*t27*t92 + d68*t78*t92 + d74*t25*t68*R(3);
//complex<R> t23 = d84*s34*t129 + d95*s34*t156*t162 + d111*spa45*spb34*t163 + d113*spb34*t184 - d93*s34*t185 + d109*s34*t210 - d108*s34*t25 + s34*spa24*t39 + d112*spb34*t176*t46 + d114*spb34*t59*t75 + d94*s34*t78*t92;
//complex<R> t42 = s14*(d98*t169 + d97*t200 + d95*t246*t261 - d93*spb13*square(spa24) + d96*spb15*square(spa24) + d94*spa23*t92*square(spb13));
//complex<R> t44 = d91*s25*t184 + d106*spb25*t246*t261 + d90*s25*t176*t46 - d104*spb13*spb25*square(spa24) + d107*spb25*t100*square(spa24) + d105*spa23*spb25*t92*square(spb13);
//complex<R> t126 = d40*t110*t162 + d41*t110*t162 + d44*t108*t193 + d45*t108*t193 + d42*t200 + d39*t46*t78 + d43*t154*t92 + d38*t78*t92;
//complex<R> t143 = t69*(d128*t156*t162 - d126*t185 + d127*t78*t92);
//complex<R> co1 = d99*spb15*t224;
//complex<R> co2 = d133*s12*spb15*t224;
//complex<R> co3 = Complex(0,1);
//SeriesC<R> result = co3*(
//		t5*(*CI_users[0]->get_value(ep,mu))
//		+ t126*(*CI_users[1]->get_value(ep,mu))
//		+ t192*(*CI_users[2]->get_value(ep,mu))
//		+ t1*(*CI_users[3]->get_value(ep,mu))
//		+ t2*(*CI_users[4]->get_value(ep,mu))
//		+ t4*(*CI_users[5]->get_value(ep,mu))
//		+ t3*(*CI_users[6]->get_value(ep,mu))
//		+ t22*(*CI_users[7]->get_value(ep,mu))
//		+ t42*(*CI_users[8]->get_value(ep,mu))
//		+ co1*(*CI_users[9]->get_value(ep,mu))
//		+ t21*(*CI_users[10]->get_value(ep,mu))
//		+ t44*(*CI_users[11]->get_value(ep,mu))
//		+ t23*(*CI_users[12]->get_value(ep,mu))
//		+ t24*(*CI_users[13]->get_value(ep,mu))
//		+ t20*(*CI_users[14]->get_value(ep,mu))
//		+ t153*(*CI_users[15]->get_value(ep,mu))
//		+ t143*(*CI_users[16]->get_value(ep,mu))
//		+ t19*(*CI_users[17]->get_value(ep,mu))
//		+ co2*(*CI_users[18]->get_value(ep,mu))
//		+ t141*(*CI_users[19]->get_value(ep,mu))
//		+ t175*(*CI_users[20]->get_value(ep,mu))
//
//);
// return(result);
//}
//
//
//C2q2g1y_qpqmmpgap_SLC_wCI::C2q2g1y_qpqmmpgap_SLC_wCI(const vector<int>& ind){
//
////{{qp, qm, m, p, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpqmmpgap SLC");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Tue 18 Nov 2008 00:14:55 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c12,c345));
//CI_users.push_back(new Cached_Bubble_Integral_User(c14,c235));
//CI_users.push_back(new Cached_Bubble_Integral_User(c15,c234));
//CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//CI_users.push_back(new Cached_Bubble_Integral_User(c25,c134));
//CI_users.push_back(new Cached_Bubble_Integral_User(c34,c125));
//CI_users.push_back(new Cached_Bubble_Integral_User(c35,c124));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c2,c345));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c4,c235));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c5,c234));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c3,c145));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c5,c134));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c4,c125));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c5,c124));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c5,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c4,c15));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c5,c14) );
//CI_users.push_back(new Cached_Box_Integral_User(c4,c1,c2,c35));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c1,c2,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c3,c4,c12));
//
//
//
//}
//
//
//SeriesC<R> C2q2g1y_qpqmmpgap_SLC_wCI::eval(const eval_param<R>& ep,const T& mu){
////{{qp, qm, m, p, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpqmmpgap SLC");
//#endif
//
//	//define corner vectors
//
// // #define TimeStamp "Tue 18 Nov 2008 00:14:55 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa15 = SPA(1,5);
//complex<R> spa23 = SPA(2,3);
//complex<R> spa24 = SPA(2,4);
//complex<R> spa45 = SPA(4,5);
//complex<R> spa13 = SPA(1,3);
//complex<R> spa14 = SPA(1,4);
//complex<R> spa25 = SPA(2,5);
//complex<R> spa34 = SPA(3,4);
//complex<R> spb15 = SPB(1,5);
//complex<R> spa12 = SPA(1,2);
//complex<R> spa35 = SPA(3,5);
//complex<R> spb14 = SPB(1,4);
//complex<R> spb24 = SPB(2,4);
//complex<R> spb25 = SPB(2,5);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb45 = SPB(4,5);
//complex<R> spb34 = SPB(3,4);
//complex<R> spb35 = SPB(3,5);
//complex<R> s23 = S(2,3);
//complex<R> s34 = -(spa34*spb34);
//complex<R> s25 = -(spa25*spb25);
//complex<R> s12 = -(spa12*spb12);
//complex<R> s24 = -(spa24*spb24);
//complex<R> s14 = -(spa14*spb14);
//complex<R> s35 = -(spa35*spb35);
//complex<R> s15 = -(spa15*spb15);
//complex<R> t5 = spa25*R(2);
//complex<R> t8 = square(s12 - s34);
//complex<R> t13 = spa15*spa45;
//complex<R> t20 = square(spb45);
//complex<R> t21 = square(spb15);
//complex<R> t23 = spa12*spa35;
//complex<R> t24 = spa14*spa25;
//complex<R> t27 = spa45*R(2);
//complex<R> t30 = s34*s35;
//complex<R> t39 = square(spa23);
//complex<R> t40 = square(spa13);
//complex<R> t41 = s12 - s34 - s35;
//complex<R> t42 = cube(spb15);
//complex<R> t43 = cube(spb45);
//complex<R> t45 = s12 + s25 - s34;
//complex<R> t49 = s14 - s35;
//complex<R> t50 = s12 + s14 - s35;
//complex<R> t52 = s15 - s34;
//complex<R> t53 = s12 + s15 - s34;
//complex<R> t64 = spa14*R(2);
//complex<R> t65 = cube(spa23);
//complex<R> t66 = -(spa13*R(3));
//complex<R> t69 = s25*spb25;
//complex<R> t70 = s24*spb24;
//complex<R> t71 = spa14*spa45;
//complex<R> t74 = spa24*R(2);
//complex<R> t90 = spa12*R(2);
//complex<R> t97 = spa24*spa35;
//complex<R> t115 = spa15*spa24;
//complex<R> t116 = -(spb45*R(3));
//complex<R> t128 = -(spa12*spb14);
//complex<R> t155 = spa25*spa34;
//complex<R> t159 = spa13*spa23;
//complex<R> t170 = spa35*spb15;
//complex<R> d14 = (s12 - s34)*spa12*spa45; d14 = R(1)/d14;
//complex<R> d15 = (s12 - s35)*spa12*spa45; d15 = R(1)/d15;
//complex<R> d31 = spa12*spa45; d31 = R(1)/d31;
//complex<R> d32 = s12 - s34; d32 = R(1)/d32;
//complex<R> d33 = s12 - s35; d33 = R(1)/d33;
//complex<R> d66 = spa24*spa45; d66 = R(1)/d66;
//complex<R> t6 = square(s12 + t49);
//complex<R> t7 = spa14*t41;
//complex<R> t9 = square(s12 + t52);
//complex<R> t22 = -t159;
//complex<R> t44 = -t97;
//complex<R> t46 = -t70;
//complex<R> t47 = -t69;
//complex<R> t48 = spa15*spa23 + spa13*t5;
//complex<R> t51 = spa23*(t115 + t24) + spa13*spa24*t5;
//complex<R> t56 = spa34*t128 + spb45*t97*R(3);
//complex<R> t59 = t116*t155 - spb15*t23;
//complex<R> t67 = -(t23*t42);
//complex<R> t68 = spa34*t43;
//complex<R> t83 = spa23*t20;
//complex<R> t87 = t21*(spa15*spa23 + spa13*spa25*R(2));
//complex<R> t98 = -(t40*R(3));
//complex<R> t99 = -t155;
//complex<R> t109 = spa14*t5;
//complex<R> t112 = spa23*t21;
//complex<R> t113 = t23*t40;
//complex<R> t114 = spa34*t64;
//complex<R> t137 = t39*t69;
//complex<R> t139 = spa34*t20;
//complex<R> t144 = spa13*t39;
//complex<R> t146 = spa34*t70;
//complex<R> t148 = d31*spb45;
//complex<R> t152 = -(t20*t39);
//complex<R> t153 = -t170;
//complex<R> t154 = spa24*t66;
//complex<R> t167 = spa25*t74;
//complex<R> t178 = spa34*t39;
//complex<R> d2 = spa34*t13*t90; d2 = R(1)/d2;
//complex<R> d3 = spa14*t23*t27; d3 = R(1)/d3;
//complex<R> d4 = (s12 - s35)*spa35*t24; d4 = R(1)/d4;
//complex<R> d5 = (s12 - s34)*spa14*spa34*t45; d5 = R(1)/d5;
//complex<R> d6 = (s12 - s34)*spa34*t24*t45; d6 = R(1)/d6;
//complex<R> d8 = (s12 - s34)*spa14*spa34*square(t45); d8 = R(1)/d8;
//complex<R> d10 = (s12 - s35)*spa35*t24*(s12 + t49); d10 = R(1)/d10;
//complex<R> d11 = (s12 - s35)*spa35*(s12 + t49)*t71; d11 = R(1)/d11;
//complex<R> d12 = (s12 - s34)*spa34*t13*(s12 + t52); d12 = R(1)/d12;
//complex<R> d13 = (s12 - s34)*t24; d13 = R(1)/d13;
//complex<R> d16 = (s12 - s34)*spa12*t41; d16 = R(1)/d16;
//complex<R> d17 = (s12 - s35)*spa12*t41; d17 = R(1)/d17;
//complex<R> d22 = t64*t8; d22 = R(1)/d22;
//complex<R> d23 = spa14*t27*square(s12 - s35); d23 = R(1)/d23;
//complex<R> d24 = t13*t8*R(2); d24 = R(1)/d24;
//complex<R> d25 = (s12 - s34)*spa14*square(t41); d25 = R(1)/d25;
//complex<R> d26 = (s12 - s35)*spa14*square(t41); d26 = R(1)/d26;
//complex<R> d29 = (s12 - s34)*spa34*t13*t90; d29 = R(1)/d29;
//complex<R> d30 = (s12 - s35)*spa14*t23*t27; d30 = R(1)/d30;
//complex<R> d34 = spa35*t49*t71; d34 = R(1)/d34;
//complex<R> d35 = spa35*t24*t49; d35 = R(1)/d35;
//complex<R> d36 = spa35*t24*t49*(s12 + t49); d36 = R(1)/d36;
//complex<R> d37 = spa35*t49*(s12 + t49)*t71; d37 = R(1)/d37;
//complex<R> d38 = (s14 - s23)*t71; d38 = R(1)/d38;
//complex<R> d39 = (s14 - s23)*(-s23 + t49)*t71; d39 = R(1)/d39;
//complex<R> d40 = t49*(-s23 + t49)*t71; d40 = R(1)/d40;
//complex<R> d41 = spa34*t13*t52; d41 = R(1)/d41;
//complex<R> d42 = t13*t52; d42 = R(1)/d42;
//complex<R> d43 = spa34*t13*t52*(s12 + t52); d43 = R(1)/d43;
//complex<R> d44 = (s25 - s34)*spa34*t24; d44 = R(1)/d44;
//complex<R> d46 = (s25 - s34)*spa14*spa34*t45; d46 = R(1)/d46;
//complex<R> d47 = (s25 - s34)*spa34*t24*t45; d47 = R(1)/d47;
//complex<R> d48 = (s25 - s34)*spa14*spa34*square(t45); d48 = R(1)/d48;
//complex<R> d50 = spa34*spa35*t24; d50 = R(1)/d50;
//complex<R> d51 = spa34*t13; d51 = R(1)/d51;
//complex<R> d52 = spa35*t71; d52 = R(1)/d52;
//complex<R> d53 = spa14*spa34*square(t45); d53 = R(1)/d53;
//complex<R> d54 = spa34*t24*square(t45); d54 = R(1)/d54;
//complex<R> d55 = spa14*spa34*cube(t45); d55 = R(1)/d55;
//complex<R> d59 = square(t41); d59 = R(1)/d59;
//complex<R> d60 = spa14*square(t41); d60 = R(1)/d60;
//complex<R> d61 = t24*square(t41); d61 = R(1)/d61;
//complex<R> d62 = spa14*cube(t41); d62 = R(1)/d62;
//complex<R> d65 = spa45*square(s23 - t49); d65 = R(1)/d65;
//complex<R> d68 = spa24*t13; d68 = R(1)/d68;
//complex<R> d69 = t71*square(s23 - t49); d69 = R(1)/d69;
//complex<R> d70 = spa14*square(t45); d70 = R(1)/d70;
//complex<R> d71 = t24*square(t45); d71 = R(1)/d71;
//complex<R> d72 = spa14*cube(t45); d72 = R(1)/d72;
//complex<R> d74 = spa12*square(t41); d74 = R(1)/d74;
//complex<R> d75 = spa12*spa14*square(t41); d75 = R(1)/d75;
//complex<R> d76 = spa12*t24*square(t41); d76 = R(1)/d76;
//complex<R> d81 = t13*t74; d81 = R(1)/d81;
//complex<R> d82 = spa14*t27*square(s23 - t49); d82 = R(1)/d82;
//complex<R> d86 = t90*square(t41); d86 = R(1)/d86;
//complex<R> d87 = spa12*t64*square(t41); d87 = R(1)/d87;
//complex<R> d89 = t64*cube(t41); d89 = R(1)/d89;
//complex<R> t36 = d30*(spa34*t128 + spb45*t97);
//complex<R> t76 = d3*t56;
//complex<R> t79 = d34*spb12;
//complex<R> t82 = d2*t59;
//complex<R> t94 = t39*(d38*spb25 + d39*t47);
//complex<R> t100 = -t148;
//complex<R> t101 = spb14*t46;
//complex<R> t103 = -(d41*spb12);
//complex<R> t104 = -(d38*spb25);
//complex<R> t120 = -(d35*spb24);
//complex<R> t123 = s23*(d69*t137 + d68*t39);
//complex<R> t126 = t68*t97;
//complex<R> t132 = spb25*t48;
//complex<R> t138 = t113*t42;
//complex<R> t175 = t112*t98;
//complex<R> t206 = t137*t170;
//complex<R> d1 = spa34*t109*t23; d1 = R(1)/d1;
//complex<R> d7 = spa34*t109*t8; d7 = R(1)/d7;
//complex<R> d9 = t114*t45*t8; d9 = R(1)/d9;
//complex<R> d18 = (s12 - s34)*spa12*t7; d18 = R(1)/d18;
//complex<R> d19 = (s12 - s35)*spa12*t7; d19 = R(1)/d19;
//complex<R> d20 = (s12 - s34)*spa12*spa25*t7; d20 = R(1)/d20;
//complex<R> d21 = (s12 - s35)*spa12*spa25*t7; d21 = R(1)/d21;
//complex<R> d27 = t7*t8*R(2); d27 = R(1)/d27;
//complex<R> d28 = t7*square(s12 - s35)*R(2); d28 = R(1)/d28;
//complex<R> d45 = t114*square(s25 - s34); d45 = R(1)/d45;
//complex<R> d49 = t114*t45*square(s25 - s34); d49 = R(1)/d49;
//complex<R> d56 = spa35*t24*t6; d56 = R(1)/d56;
//complex<R> d57 = spa35*t6*t71; d57 = R(1)/d57;
//complex<R> d58 = spa34*t13*t9; d58 = R(1)/d58;
//complex<R> d63 = spa25*spa35*t6; d63 = R(1)/d63;
//complex<R> d64 = spa35*spa45*t6; d64 = R(1)/d64;
//complex<R> d67 = spa34*spa45*t9; d67 = R(1)/d67;
//complex<R> d73 = t13*t9; d73 = R(1)/d73;
//complex<R> d77 = t24*t6; d77 = R(1)/d77;
//complex<R> d78 = t6*t71; d78 = R(1)/d78;
//complex<R> d79 = t114*square(t45); d79 = R(1)/d79;
//complex<R> d80 = t114*cube(t45); d80 = R(1)/d80;
//complex<R> d83 = spa35*t5*t6; d83 = R(1)/d83;
//complex<R> d84 = spa35*t27*t6; d84 = R(1)/d84;
//complex<R> d85 = spa34*t27*t9; d85 = R(1)/d85;
//complex<R> d88 = spa12*t109*square(t41); d88 = R(1)/d88;
//complex<R> t14 = d44*spb15*t144 - d45*t112*t40 + d48*t40*t67 + d49*t40*t67 + d47*t22*t87 + d46*t112*t40*R(3);
//complex<R> t15 = d55*s12*t138 + d51*spb12*t144 - d52*spb12*t144 + d59*spb12*t152 + d53*s12*t175 + d57*s12*t146*t39 + d58*s12*spa35*t39*t47 - d50*spa13*spb12*t65 + d62*s12*t44*t68 + d56*s12*t65*t70 + d60*spb12*t154*t83 + d61*spb12*(spa13*t167 + spa23*(t115 + t24))*t83 + d54*s12*t159*t87;
//complex<R> t16 = d62*s34*t126 + d72*spb34*t138 + d74*s34*t152 + d70*spb34*t175 + d68*s34*t39 + d73*spa35*spb34*t39*t47 + d75*s34*t154*t83 + d76*s34*(spa13*t167 + spa23*(t115 + t24))*t83 + d71*spb34*t159*t87;
//complex<R> t17 = d62*s35*t126 + d69*s35*t137 + d74*s35*t152 + d78*spb35*t146*t39 + d77*spb35*t65*t70 + d75*s35*t154*t83 + d76*s35*(spa13*t167 + spa23*(t115 + t24))*t83;
//complex<R> t18 = t30*(d89*t126 + d86*t152 + d87*t154*t83 + d88*(spa13*t167 + spa23*(t115 + t24))*t83);
//complex<R> t77 = d79*s12;
//complex<R> t110 = s12*t101*(d84*t178 + d83*t65);
//complex<R> t136 = d43*spa35*t137 + t103*t144 - d42*spb24*t39;
//complex<R> t140 = d18*spa24;
//complex<R> t143 = d67*t206 + d66*spb15*t39;
//complex<R> t149 = d19*spa24;
//complex<R> t151 = d65*spb14*t137 + d64*t101*t178 + d63*t101*t65;
//complex<R> t162 = t132*t21;
//complex<R> t168 = d39*t137 + d40*t137 + t104*t39 + d37*t178*t46 + t120*t65 + d36*t46*t65 + t144*t79;
//complex<R> t1 = d48*t138 + d49*t138 + d8*t138 + d9*t138 + d41*spb12*t144 - d44*spb15*t144 + d46*t175 + d5*t175 + d24*spa35*t155*t20 + d7*t21*t22*t23 + d42*spb24*t39 + d13*spb45*t39 + d14*spb45*t39 + d16*t20*t39 + d32*(d2*t159*(t116*t155 - spb15*t23) + t100*t39) + d45*t112*t40 + d12*spa35*t39*t47 + d43*spa35*t39*t47 + d47*t159*t21*t48 + d6*t159*t21*t48 + d25*t44*t68 + d27*t44*t68 - d22*spa34*t83 - d20*(spa23*(t115 + t24) + spa13*spa24*t5)*t83 + d29*spa23*t66*(-(spb15*t23) + spb45*t99) + t140*t159*t20*R(3);
//complex<R> t2 = -(d4*spb14*t144) + d15*spb45*t39 + d11*t146*t39 + d37*t146*t39 + d17*t20*t39 + d23*t139*t44 + d40*t39*t47 + d35*spb24*t65 + d26*t44*t68 + d28*t44*t68 + d10*t65*t70 + d36*t65*t70 + d33*(t100*t39 + t22*t76) - t144*t79 - d21*(spa23*(t115 + t24) + spa13*spa24*t5)*t83 + t149*t159*t20*R(3) + t159*t36*R(3);
//complex<R> t3 = d25*t126 + d26*t126 + d27*t126 + d28*t126 + d12*spa35*t137 + d4*spb14*t144 + d16*t152 + d17*t152 + d7*t159*t21*t23 - d13*spb45*t39 - d14*spb45*t39 - d15*spb45*t39 + d32*(d2*t22*(t116*t155 - spb15*t23) + t148*t39) + d11*t178*t46 + d6*t21*t22*t48 + d10*t46*t65 + spa23*t36*t66 + d3*t39*t66 + d1*t65*t66 + d8*t40*t67 + d9*t40*t67 + d33*(t148*t39 + t159*t76) + d22*spa34*t83 + d20*(spa23*(t115 + t24) + spa13*spa24*t5)*t83 + d21*(spa23*(t115 + t24) + spa13*spa24*t5)*t83 + t140*t66*t83 + t149*t66*t83 + d23*t139*t97 + d24*spa35*t20*t99 + d2*t144*R(3) + d5*t112*t40*R(3) + d29*t159*(-(spb15*t23) + spb45*t99)*R(3);
//complex<R> t95 = d55*s25*t138 + d53*s25*t175 + d53*t162*t22;
//complex<R> t124 = d80*s12*s25*t138 + s25*t175*t77 + t162*t22*t77;
//complex<R> co1 = d81*s23*s34*t39;
//complex<R> co2 = d82*s23*s35*t137;
//complex<R> co3 = d85*s12*t206;
//complex<R> co4 = Complex(0,1);
//SeriesC<R> result = co4*(
//		t3*(*CI_users[0]->get_value(ep,mu))
//		+ t168*(*CI_users[1]->get_value(ep,mu))
//		+ t136*(*CI_users[2]->get_value(ep,mu))
//		+ t94*(*CI_users[3]->get_value(ep,mu))
//		+ t14*(*CI_users[4]->get_value(ep,mu))
//		+ t1*(*CI_users[5]->get_value(ep,mu))
//		+ t2*(*CI_users[6]->get_value(ep,mu))
//		+ t15*(*CI_users[7]->get_value(ep,mu))
//		+ t151*(*CI_users[8]->get_value(ep,mu))
//		+ t143*(*CI_users[9]->get_value(ep,mu))
//		+ t123*(*CI_users[10]->get_value(ep,mu))
//		+ t95*(*CI_users[11]->get_value(ep,mu))
//		+ t16*(*CI_users[12]->get_value(ep,mu))
//		+ t17*(*CI_users[13]->get_value(ep,mu))
//		+ t124*(*CI_users[14]->get_value(ep,mu))
//		+ co1*(*CI_users[15]->get_value(ep,mu))
//		+ co2*(*CI_users[16]->get_value(ep,mu))
//		+ t110*(*CI_users[17]->get_value(ep,mu))
//		+ co3*(*CI_users[18]->get_value(ep,mu))
//		+ t18*(*CI_users[19]->get_value(ep,mu))
//		);
// return(result);
//}
//
//
//C2q2g1y_qpqmmmgap_SLC_wCI::C2q2g1y_qpqmmmgap_SLC_wCI(const vector<int>& ind){
//
////{{qp, qm, m, m, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpqmmmgap SLC");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Tue 18 Nov 2008 00:25:28 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c12,c345));
//CI_users.push_back(new Cached_Bubble_Integral_User(c14,c235));
//CI_users.push_back(new Cached_Bubble_Integral_User(c15,c234));
//CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//CI_users.push_back(new Cached_Bubble_Integral_User(c25,c134));
//CI_users.push_back(new Cached_Bubble_Integral_User(c35,c124));
//CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c2,c345));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c4,c235));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c5,c234));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c3,c145));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c5,c134));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c5,c124));
//CI_users.push_back(new Cached_Triangle_Integral_User(c4,c5,c123));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c5,c4,c23));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c1,c4,c35));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c2,c1,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c5,c2,c14));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c1,c2,c35));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c5,c3,c12));
//
//}
//
//
//SeriesC<R> C2q2g1y_qpqmmmgap_SLC_wCI::eval(const eval_param<R>& ep,const T& mu){
//
//
////{{qp, qm, m, m, gap}, SLC}
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpqmmmgap SLC");
//#endif
//
//	//define corner vectors
//
// // #define TimeStamp "Tue 18 Nov 2008 00:25:28 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa14 = SPA(1,4);
//complex<R> spb15 = SPB(1,5);
//complex<R> spb23 = SPB(2,3);
//complex<R> spb34 = SPB(3,4);
//complex<R> spa23 = SPA(2,3);
//complex<R> spb14 = SPB(1,4);
//complex<R> spb13 = SPB(1,3);
//complex<R> spb25 = SPB(2,5);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb35 = SPB(3,5);
//complex<R> spb45 = SPB(4,5);
//complex<R> spa12 = SPA(1,2);
//complex<R> spa13 = SPA(1,3);
//complex<R> spa24 = SPA(2,4);
//complex<R> spb24 = SPB(2,4);
//complex<R> spa45 = SPA(4,5);
//complex<R> spa35 = SPA(3,5);
//complex<R> spa34 = SPA(3,4);
//complex<R> s15 = S(1,5);
//complex<R> s23 = -(spa23*spb23);
//complex<R> s14 = -(spa14*spb14);
//complex<R> s45 = -(spa45*spb45);
//complex<R> s25 = S(2,5);
//complex<R> s35 = -(spa35*spb35);
//complex<R> s12 = -(spa12*spb12);
//complex<R> s34 = -(spa34*spb34);
//complex<R> t6 = spb14*R(2);
//complex<R> t8 = square(s12 - s35);
//complex<R> t9 = spb23*spb34;
//complex<R> t21 = square(spb15);
//complex<R> t22 = square(spa24);
//complex<R> t23 = square(spa34);
//complex<R> t24 = spb12*spb45;
//complex<R> t25 = -(spb25*R(3));
//complex<R> t26 = spb14*spb23;
//complex<R> t28 = spb34*R(2);
//complex<R> t38 = cube(spa24);
//complex<R> t39 = square(spb25);
//complex<R> t40 = square(spa23);
//complex<R> t41 = -(spb13*spb15);
//complex<R> t43 = square(spb45);
//complex<R> t45 = s14 - s35;
//complex<R> t47 = cube(spa34);
//complex<R> t49 = spb15*spb24 + spb14*spb25*R(2);
//complex<R> t50 = spb23*R(2);
//complex<R> t51 = s14 - s25;
//complex<R> t52 = s15 - s23 + s45;
//complex<R> t53 = s12*s14;
//complex<R> t65 = cube(spa23);
//complex<R> t67 = -(spb12*spb25);
//complex<R> t70 = spb14*spb34;
//complex<R> t72 = s25*s35;
//complex<R> t88 = cube(spb15);
//complex<R> t90 = -(spa12*spb25);
//complex<R> t97 = spb13*spb25;
//complex<R> t100 = -(spa23*spb35);
//complex<R> t101 = -(spb14*R(2));
//complex<R> t103 = s14*spa14;
//complex<R> t114 = spb13*spb45;
//complex<R> t116 = spb12*spb35;
//complex<R> t152 = spb23*spb35;
//complex<R> t157 = spb25*spb45;
//complex<R> d1 = s34*(s12 - s35)*spb12; d1 = R(1)/d1;
//complex<R> d2 = s34*(s12 - s45)*spb12; d2 = R(1)/d2;
//complex<R> d4 = s34*(s12 - s35)*spb12*spb23; d4 = R(1)/d4;
//complex<R> d5 = s34*(s12 - s45)*spb12*spb23; d5 = R(1)/d5;
//complex<R> d8 = (s12 - s35)*spb12*spb34; d8 = R(1)/d8;
//complex<R> d9 = (s12 - s45)*spb12*spb34; d9 = R(1)/d9;
//complex<R> d24 = (s12 - s35)*spb23*square(s34); d24 = R(1)/d24;
//complex<R> d26 = (s12 - s45)*spb23*square(s34); d26 = R(1)/d26;
//complex<R> d31 = spb12*spb34; d31 = R(1)/d31;
//complex<R> d33 = s12 - s35; d33 = R(1)/d33;
//complex<R> d34 = s12 - s45; d34 = R(1)/d34;
//complex<R> d61 = square(spb34); d61 = R(1)/d61;
//complex<R> d62 = spb23*square(spb34); d62 = R(1)/d62;
//complex<R> d73 = spb23*cube(spb34); d73 = R(1)/d73;
//complex<R> d88 = spb12*square(spb34); d88 = R(1)/d88;
//complex<R> d89 = spb12*spb23*square(spb34); d89 = R(1)/d89;
//complex<R> d108 = spb12*square(spb34)*R(2); d108 = R(1)/d108;
//complex<R> t10 = spb35*t45;
//complex<R> t42 = s12 + t45;
//complex<R> t44 = -t116;
//complex<R> t48 = -(spb13*spb15*spb24) - spb15*t26 + t101*t97;
//complex<R> t57 = spa23*t116 + spa34*t114*R(3);
//complex<R> t58 = -t72;
//complex<R> t66 = -t114;
//complex<R> t69 = -t103;
//complex<R> t74 = -(d31*R(2));
//complex<R> t81 = spb15*t22;
//complex<R> t85 = spb25*(spb15*spb24 + spb25*t6);
//complex<R> t115 = spb15*t23;
//complex<R> t117 = spb25*t38;
//complex<R> t128 = t65*t97;
//complex<R> t130 = spb35*t47;
//complex<R> t131 = spb12*t43;
//complex<R> t140 = spb35*t114;
//complex<R> t150 = t24*t38;
//complex<R> t151 = spb25*t40;
//complex<R> t153 = d89*spb13;
//complex<R> t158 = spa24*t21;
//complex<R> t170 = spb35*t21;
//complex<R> t174 = t38*t43;
//complex<R> t185 = spb25*t21;
//complex<R> t187 = spb15*t25;
//complex<R> t200 = spb15*t40;
//complex<R> d3 = (s12 - s35)*t26; d3 = R(1)/d3;
//complex<R> d6 = s34*(s12 - s35)*spb12*t26; d6 = R(1)/d6;
//complex<R> d7 = s34*(s12 - s45)*spb12*t26; d7 = R(1)/d7;
//complex<R> d12 = t116*t70; d12 = R(1)/d12;
//complex<R> d13 = t50*t8; d13 = R(1)/d13;
//complex<R> d14 = spb12*t28*t8; d14 = R(1)/d14;
//complex<R> d15 = (s12 - s45)*spb45*t26; d15 = R(1)/d15;
//complex<R> d16 = t24*t9*R(2); d16 = R(1)/d16;
//complex<R> d17 = t152*t24*t6; d17 = R(1)/d17;
//complex<R> d18 = (s12 - s45)*spb45*t9; d18 = R(1)/d18;
//complex<R> d19 = t152*t6*t8; d19 = R(1)/d19;
//complex<R> d23 = s34*t50*t8; d23 = R(1)/d23;
//complex<R> d25 = s34*t50*square(s12 - s45); d25 = R(1)/d25;
//complex<R> d27 = t9*square(s12 - s45)*R(2); d27 = R(1)/d27;
//complex<R> d30 = (s12 - s45)*t24*t9*R(2); d30 = R(1)/d30;
//complex<R> d32 = spb35*t70; d32 = R(1)/d32;
//complex<R> d35 = spb34*t45*t6; d35 = R(1)/d35;
//complex<R> d36 = spb34*t6*square(t51); d36 = R(1)/d36;
//complex<R> d37 = (-s25 + t45)*t51*t70; d37 = R(1)/d37;
//complex<R> d38 = t45*(-s25 + t45)*t70; d38 = R(1)/d38;
//complex<R> d39 = spb35*t28; d39 = R(1)/d39;
//complex<R> d40 = spb34*t116; d40 = R(1)/d40;
//complex<R> d41 = square(t45); d41 = R(1)/d41;
//complex<R> d43 = spb35*t50*square(t45); d43 = R(1)/d43;
//complex<R> d46 = spb34*t116*t6; d46 = R(1)/d46;
//complex<R> d47 = t51*t70*square(s25 - t45); d47 = R(1)/d47;
//complex<R> d48 = t45*t70*square(s25 - t45); d48 = R(1)/d48;
//complex<R> d49 = spb34*(-s25 + t45)*t6*square(t51); d49 = R(1)/d49;
//complex<R> d50 = spb34*(-s25 + t45)*t6*square(t45); d50 = R(1)/d50;
//complex<R> d57 = (s15 - s23)*t9; d57 = R(1)/d57;
//complex<R> d58 = (s15 - s23)*t52*t9; d58 = R(1)/d58;
//complex<R> d59 = (s23 - s45)*t52*t9; d59 = R(1)/d59;
//complex<R> d60 = (s23 - s45)*spb45*t9; d60 = R(1)/d60;
//complex<R> d63 = t26*square(spb34); d63 = R(1)/d63;
//complex<R> d66 = t114*t26; d66 = R(1)/d66;
//complex<R> d67 = spb45*t9; d67 = R(1)/d67;
//complex<R> d68 = spb35*spb45*t26; d68 = R(1)/d68;
//complex<R> d69 = t114*t9; d69 = R(1)/d69;
//complex<R> d75 = spb34*(-s25 + t45); d75 = R(1)/d75;
//complex<R> d76 = spb34*square(s25 - t45); d76 = R(1)/d76;
//complex<R> d77 = spb34*cube(-s25 + t45); d77 = R(1)/d77;
//complex<R> d79 = t9*square(t52); d79 = R(1)/d79;
//complex<R> d80 = spb34*square(t52); d80 = R(1)/d80;
//complex<R> d81 = spb14*t114; d81 = R(1)/d81;
//complex<R> d82 = spb34*t114; d82 = R(1)/d82;
//complex<R> d83 = (-s25 + t45)*t70; d83 = R(1)/d83;
//complex<R> d84 = t70*square(s25 - t45); d84 = R(1)/d84;
//complex<R> d85 = t70*cube(-s25 + t45); d85 = R(1)/d85;
//complex<R> d90 = spb12*t26*square(spb34); d90 = R(1)/d90;
//complex<R> d95 = spb13*t26; d95 = R(1)/d95;
//complex<R> d96 = spb13*t9; d96 = R(1)/d96;
//complex<R> d97 = t9*square(t52)*R(2); d97 = R(1)/d97;
//complex<R> d100 = t114*t6; d100 = R(1)/d100;
//complex<R> d101 = t114*t28; d101 = R(1)/d101;
//complex<R> d102 = spb34*(-s25 + t45)*t6; d102 = R(1)/d102;
//complex<R> d103 = spb34*t6*square(s25 - t45); d103 = R(1)/d103;
//complex<R> d104 = spb34*t6*cube(-s25 + t45); d104 = R(1)/d104;
//complex<R> d109 = spb12*t50*square(spb34); d109 = R(1)/d109;
//complex<R> d110 = spb12*spb23*t6*square(spb34); d110 = R(1)/d110;
//complex<R> d111 = t50*cube(spb34); d111 = R(1)/d111;
//complex<R> t5 = square(t42);
//complex<R> t16 = s35*s45*(d111*t140 + d109*spb13*t187 - d108*t21 + d110*spb15*(spb15*t26 - spb24*t41 - t101*t97));
//complex<R> t27 = spb35*t42;
//complex<R> t64 = -t81;
//complex<R> t68 = cube(t42);
//complex<R> t80 = d16*t57;
//complex<R> t99 = -t150;
//complex<R> t105 = d39*spa12;
//complex<R> t121 = -(d57*spa14);
//complex<R> t123 = d14*spb14;
//complex<R> t125 = (d57*spa14 + d58*t103)*t21;
//complex<R> t129 = -t158;
//complex<R> t134 = -(d60*spa13);
//complex<R> t136 = d40*spb14;
//complex<R> t137 = -(d90*t48);
//complex<R> t139 = t39*t81;
//complex<R> t146 = d79*t69;
//complex<R> t164 = t117*t131;
//complex<R> t169 = t150*t39;
//complex<R> t172 = s12*(d101*t100*t21 + d100*spa23*t88);
//complex<R> t173 = t128*t44;
//complex<R> t175 = t21*t69;
//complex<R> t178 = t151*t41;
//complex<R> t208 = t115*R(3);
//complex<R> t210 = spb13*t115;
//complex<R> d42 = t10*t26; d42 = R(1)/d42;
//complex<R> d44 = spb23*t10*t42; d44 = R(1)/d44;
//complex<R> d45 = t10*t26*t42; d45 = R(1)/d45;
//complex<R> d51 = spb12*spb34*t10*t6; d51 = R(1)/d51;
//complex<R> d54 = spb34*t10*t42; d54 = R(1)/d54;
//complex<R> d92 = t42*t70; d92 = R(1)/d92;
//complex<R> t1 = d18*spa13*t170 + d60*spa13*t170 + d59*t175 + d60*spa12*t185 - d15*spa23*t185 + d9*spa34*t21 + d27*t140*t23 + d2*t21*t23 + d25*t140*t47 + d26*t140*t47 - d34*(d31*spa34*t21 + spb15*spb25*t80) - d15*spa13*t88 + d5*t208*t97 + d7*t115*(-(spb15*t26) + spb24*t41 + t101*t97) + d30*spb15*spb25*(spa34*t114 + spa23*t116)*R(3);
//complex<R> t34 = d51*(spa24*spb14*spb25 + spa23*t44);
//complex<R> t35 = s25*(d85*t173 + d84*t178 - d83*spa23*square(spb15));
//complex<R> t61 = -t136;
//complex<R> t76 = -t105;
//complex<R> t113 = d47*t116*t128 + d49*t116*t128 + d36*t178 + d37*t200*t97;
//complex<R> t127 = d102*spa23*t21*t58 + d104*t173*t72 + d103*t178*t72;
//complex<R> t168 = d80*spa23*t175 + d82*t100*t21 + d81*spa23*t88;
//complex<R> t177 = t134*t170 + d58*t175 + t21*(d59*t103 + t121 + d60*t90);
//complex<R> t182 = t157*t64;
//complex<R> t231 = t146*t21;
//complex<R> d10 = (s12 - s35)*spb23*t27; d10 = R(1)/d10;
//complex<R> d11 = (s12 - s35)*t26*t27; d11 = R(1)/d11;
//complex<R> d20 = (s12 - s35)*t152*t5; d20 = R(1)/d20;
//complex<R> d21 = t27*t50*t8; d21 = R(1)/d21;
//complex<R> d22 = (s12 - s35)*spb34*t27; d22 = R(1)/d22;
//complex<R> d28 = (s12 - s35)*spb34*spb35*t5; d28 = R(1)/d28;
//complex<R> d29 = t27*t28*t8; d29 = R(1)/d29;
//complex<R> d52 = spb23*t10*t5; d52 = R(1)/d52;
//complex<R> d53 = t27*t50*square(t45); d53 = R(1)/d53;
//complex<R> d55 = spb34*t10*t5; d55 = R(1)/d55;
//complex<R> d56 = t27*t28*square(t45); d56 = R(1)/d56;
//complex<R> d64 = t152*t5; d64 = R(1)/d64;
//complex<R> d65 = spb35*t26*t5; d65 = R(1)/d65;
//complex<R> d70 = t152*t68; d70 = R(1)/d70;
//complex<R> d71 = t27*t70; d71 = R(1)/d71;
//complex<R> d72 = spb34*spb35*t5; d72 = R(1)/d72;
//complex<R> d74 = spb34*spb35*t68; d74 = R(1)/d74;
//complex<R> d78 = spb34*t27; d78 = R(1)/d78;
//complex<R> d86 = spb23*t5; d86 = R(1)/d86;
//complex<R> d87 = t26*t5; d87 = R(1)/d87;
//complex<R> d91 = spb23*t68; d91 = R(1)/d91;
//complex<R> d93 = spb34*t5; d93 = R(1)/d93;
//complex<R> d94 = spb34*t68; d94 = R(1)/d94;
//complex<R> d98 = spb35*t5*t50; d98 = R(1)/d98;
//complex<R> d99 = spb35*t50*t68; d99 = R(1)/d99;
//complex<R> d105 = t27*t28; d105 = R(1)/d105;
//complex<R> d106 = spb35*t28*t5; d106 = R(1)/d106;
//complex<R> d107 = spb35*t28*t68; d107 = R(1)/d107;
//complex<R> t2 = d13*spb35*t115 + spb35*t115*t123 - d18*spa13*t170 + d15*spa23*t185 + d30*(spa34*t114 + spa23*t116)*t187 + d3*spa34*t21 - d8*spa34*t21 - d9*spa34*t21 + d4*spb13*t187*t23 + d5*spb13*t187*t23 - d1*t21*t23 - d2*t21*t23 + d11*spb25*t49*t64 + d23*t130*t66 + d24*t130*t66 + d25*t130*t66 + d26*t130*t66 + d27*spb35*t23*t66 + d28*t174*t67 + d29*t174*t67 + d34*(d31*spa34*t21 + spb15*spb25*t80) + d22*t157*t81 + d19*spb25*t24*t81 - d12*t88 + d15*spa13*t88 + d17*t25*t88 + d6*t115*(spb15*t26 - spb24*t41 - t101*t97) + d7*t115*(spb15*t26 - spb24*t41 - t101*t97) + d20*t39*t99 + d21*t39*t99 + d33*(d32*spb45*t129 + d31*spa34*t21*R(2)) + d10*t139*R(3) + d16*t185*R(3);
//complex<R> t3 = d42*spb25*t158 + d47*t173 + d48*t173 + d49*t173 + d50*t173 + d37*t178 + d38*t178 + d35*spa23*t21 + d41*t139*t61 + d43*t39*t64 + d45*spb25*t49*t64 + d55*t174*t67 + d56*t174*t67 + d41*spb25*t158*t76 + d54*t157*t81 - d46*t88 + d36*t200*t97 + d52*t39*t99 + d53*t39*t99 + d44*t139*R(3) + t21*t34*R(3);
//complex<R> t4 = -(d13*spb35*t115) - spb35*t115*t123 + d48*t116*t128 + d50*t116*t128 + d42*spb25*t129 + d43*t139 + d41*t136*t139 + d32*d33*spb45*t158 + d41*spb25*t105*t158 + d28*t164 + d29*t164 + d55*t164 + d56*t164 + d20*t169 + d21*t169 + d52*t169 + d53*t169 + d22*t182 + d54*t182 - d35*spa23*t21 - d3*spa34*t21 + d8*spa34*t21 + d1*t21*t23 - d6*spb15*t115*t26 + d6*spb24*t115*t41 + d23*t140*t47 + d24*t140*t47 + d19*spb25*t24*t64 + d33*spa34*t21*t74 + d11*spb25*t49*t81 + d45*spb25*t49*t81 + d6*t101*t115*t97 + d38*t200*t97 + d4*t208*t97 - d10*t139*R(3) - d44*t139*R(3) - t21*t34*R(3);
//complex<R> t14 = d99*t169*t53 + d98*s12*spa14*spb25*t49*t64 - d98*t139*t53*R(3);
//complex<R> t15 = d73*s45*t140 + d96*spa45*t170 + s45*t153*t187 - d88*s45*t21 + s45*t231 - d95*spa45*t88 + d90*s45*spb15*(spb15*t26 - spb24*t41 - t101*t97);
//complex<R> t19 = d73*s35*t140 + d92*spa35*spb45*t158 + d94*spa35*t164 + d91*spa35*t169 + d85*s35*t173 + d84*s35*t178 + d93*spa35*t182 + s35*t153*t187 - d88*s35*t21 - d83*s35*spa23*t21 + d87*spa35*spb25*t49*t81 + d90*s35*spb15*(spb15*t26 - spb24*t41 - t101*t97) - d86*spa35*t139*R(3);
//complex<R> t96 = d105*s12*spa14*spb45*t129 + d107*t164*t53 + d106*t182*t53;
//complex<R> t133 = -(d64*R(3));
//complex<R> t17 = d67*spa12*t185 + d62*spa12*spb13*t187 - d61*spa12*t21 - d63*spa12*spb15*spb24*t41 - d32*spa12*t88 + s12*(t133*t139 + d71*spb45*t158 + d74*t164 + d70*t169 + d69*t170 + d72*t182 + d73*spb35*t66 + d65*spb25*t49*t81 - d66*t88) + d68*t88*t90 - d63*spa12*spb15*t101*t97 + d63*spa12*t26*square(spb15);
//complex<R> t18 = d78*spa14*spb45*t129 + s14*(t133*t139 + d74*t164 + d70*t169 + d72*t182) + spa14*(d77*t173 + d76*t178 - d75*spa23*t21 + d64*spb25*t49*t64);
//complex<R> co1 = s15*t231;
//complex<R> co2 = d97*s15*s45*t175;
//complex<R> co3 = Complex(0,1);
//SeriesC<R> result = co3*(
//		t2*(*CI_users[0]->get_value(ep,mu))
//		+ t3*(*CI_users[1]->get_value(ep,mu))
//		+ t125*(*CI_users[2]->get_value(ep,mu))
//		+ t177*(*CI_users[3]->get_value(ep,mu))
//		+ t113*(*CI_users[4]->get_value(ep,mu))
//		+ t4*(*CI_users[5]->get_value(ep,mu))
//		+ t1*(*CI_users[6]->get_value(ep,mu))
//		+ t17*(*CI_users[7]->get_value(ep,mu))
//		+ t18*(*CI_users[8]->get_value(ep,mu))
//		+ co1*(*CI_users[9]->get_value(ep,mu))
//		+ t168*(*CI_users[10]->get_value(ep,mu))
//		+ t35*(*CI_users[11]->get_value(ep,mu))
//		+ t19*(*CI_users[12]->get_value(ep,mu))
//		+ t15*(*CI_users[13]->get_value(ep,mu))
//		+ co2*(*CI_users[14]->get_value(ep,mu))
//		+ t14*(*CI_users[15]->get_value(ep,mu))
//		+ t172*(*CI_users[16]->get_value(ep,mu))
//		+ t127*(*CI_users[17]->get_value(ep,mu))
//		+ t96*(*CI_users[18]->get_value(ep,mu))
//		+ t16*(*CI_users[19]->get_value(ep,mu))
//		);
// return(result);
//}
//
//
//C2q2g1y_qpppqmgap_nf_wCI::C2q2g1y_qpppqmgap_nf_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpppqmgap_nf_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//C2q2g1y_qppmqmgap_nf_wCI::C2q2g1y_qppmqmgap_nf_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qppmqmgap_nf_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//
//C2q2g1y_qpmmqmgap_nf_wCI::C2q2g1y_qpmmqmgap_nf_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpmmqmgap_nf_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//C2q2g1y_qpmpqmgap_nf_wCI::C2q2g1y_qpmpqmgap_nf_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpmpqmgap_nf_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//
//C2q2g1y_qpgapqmpp_L_wCI::C2q2g1y_qpgapqmpp_L_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpgapqmpp_L_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//
//C2q2g1y_qpgapqmpm_L_wCI::C2q2g1y_qpgapqmpm_L_wCI(const vector<int>& ind){
//
//
//
//
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpgapqmpm L");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Tue 18 Nov 2008 00:26:30 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c15,c234));
//CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//CI_users.push_back(new Cached_Bubble_Integral_User(c34,c125));
//CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c5,c234));
//CI_users.push_back(new Cached_Triangle_Integral_User(c2,c3,c145));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c4,c125));
//CI_users.push_back(new Cached_Triangle_Integral_User(c4,c5,c123));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c3,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c4,c15));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c2,c1,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c4,c5,c12));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c3,c2,c15));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c5,c1,c23));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c1,c2,c34));
//
//}
//
//
//SeriesC<R> C2q2g1y_qpgapqmpm_L_wCI::eval(const eval_param<R>& ep,const T& mu){
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpgapqmpm L");
//#endif
//
//
// // #define TimeStamp "Tue 18 Nov 2008 00:26:30 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa23 = SPA(2,3);
//complex<R> spa34 = SPA(3,4);
//complex<R> spa35 = SPA(3,5);
//complex<R> spa45 = SPA(4,5);
//complex<R> spb12 = SPB(1,2);
//complex<R> spa12 = SPA(1,2);
//complex<R> spa24 = SPA(2,4);
//complex<R> spa25 = SPA(2,5);
//complex<R> spa13 = SPA(1,3);
//complex<R> spb14 = SPB(1,4);
//complex<R> spa15 = SPA(1,5);
//complex<R> spb23 = SPB(2,3);
//complex<R> spb24 = SPB(2,4);
//complex<R> spb34 = SPB(3,4);
//complex<R> spb45 = SPB(4,5);
//complex<R> s15 = S(1,5);
//complex<R> s23 = -(spa23*spb23);
//complex<R> s45 = -(spa45*spb45);
//complex<R> s34 = -(spa34*spb34);
//complex<R> s24 = -(spa24*spb24);
//complex<R> t4 = square(s15 - s23);
//complex<R> t5 = spa12*R(2);
//complex<R> t7 = (s23 - s45)*spa23;
//complex<R> t20 = square(spb14);
//complex<R> t21 = square(spa35);
//complex<R> t22 = square(spb24);
//complex<R> t23 = square(spa13);
//complex<R> t25 = s15 - s23 + s45;
//complex<R> t26 = -(spa15*R(2));
//complex<R> t27 = s15 - s34;
//complex<R> t32 = cube(spb14);
//complex<R> t33 = cube(spa35);
//complex<R> t34 = cube(spb24);
//complex<R> t35 = spa23*spa45;
//complex<R> t36 = -(spa25*R(2));
//complex<R> t42 = s15*s45;
//complex<R> t49 = spa35*R(2);
//complex<R> t53 = spa15*spa45;
//complex<R> t54 = spa35*spb24;
//complex<R> d6 = (s15 - s23)*s24*spa12; d6 = R(1)/d6;
//complex<R> d9 = (s15 - s23)*spa12*square(s24); d9 = R(1)/d9;
//complex<R> d15 = s15 - s23; d15 = R(1)/d15;
//complex<R> d21 = spa12*square(spa24); d21 = R(1)/d21;
//complex<R> d22 = spa12*spa24*spa34; d22 = R(1)/d22;
//complex<R> d23 = spa12*cube(spa24); d23 = R(1)/d23;
//complex<R> d30 = spa12*spa24; d30 = R(1)/d30;
//complex<R> d31 = spa12*spa23*spa34; d31 = R(1)/d31;
//complex<R> d32 = spa34*spa45*R(2); d32 = R(1)/d32;
//complex<R> t24 = -t35;
//complex<R> t30 = -t42;
//complex<R> t46 = -t53;
//complex<R> t48 = spa25*t34;
//complex<R> t52 = t23*t32;
//complex<R> t57 = spa13*t20;
//complex<R> t58 = spa35*t26;
//complex<R> t63 = spa25*t22;
//complex<R> t72 = spb14*t21;
//complex<R> t77 = spa35*t36;
//complex<R> t86 = spb45*t33;
//complex<R> d1 = spa12*spa34*t35; d1 = R(1)/d1;
//complex<R> d2 = spa23*t4*t5; d2 = R(1)/d2;
//complex<R> d3 = (s15 - s23)*spa12*spa23*t25; d3 = R(1)/d3;
//complex<R> d4 = (s15 - s23)*spa12*spa23*square(t25); d4 = R(1)/d4;
//complex<R> d5 = spa23*t25*t4*t5; d5 = R(1)/d5;
//complex<R> d7 = t5*square(t27); d7 = R(1)/d7;
//complex<R> d8 = s24*spa12*t27; d8 = R(1)/d8;
//complex<R> d10 = s24*t4*t5; d10 = R(1)/d10;
//complex<R> d11 = s24*t5*square(t27); d11 = R(1)/d11;
//complex<R> d12 = spa12*t27*square(s24); d12 = R(1)/d12;
//complex<R> d13 = spa12*t35; d13 = R(1)/d13;
//complex<R> d14 = spa24*spa34*t5; d14 = R(1)/d14;
//complex<R> d16 = spa24*t4*t5; d16 = R(1)/d16;
//complex<R> d17 = spa34*t35*t5; d17 = R(1)/d17;
//complex<R> d18 = spa12*t25*t7; d18 = R(1)/d18;
//complex<R> d19 = spa12*t7*square(t25); d19 = R(1)/d19;
//complex<R> d20 = spa23*t25*t5*square(s23 - s45); d20 = R(1)/d20;
//complex<R> d24 = spa12*spa23*spa34*t25; d24 = R(1)/d24;
//complex<R> d25 = spa12*spa23*square(t25); d25 = R(1)/d25;
//complex<R> d26 = spa12*spa23*cube(t25); d26 = R(1)/d26;
//complex<R> d27 = spa12*spa34*t25; d27 = R(1)/d27;
//complex<R> d28 = spa12*square(t25); d28 = R(1)/d28;
//complex<R> d29 = spa12*cube(t25); d29 = R(1)/d29;
//complex<R> d33 = spa45*t5; d33 = R(1)/d33;
//complex<R> d34 = spa23*t5; d34 = R(1)/d34;
//complex<R> d35 = t5*cube(spa24); d35 = R(1)/d35;
//complex<R> d36 = spa24*t5; d36 = R(1)/d36;
//complex<R> d37 = spa23*spa34*t25*t5; d37 = R(1)/d37;
//complex<R> d38 = spa23*t5*cube(t25); d38 = R(1)/d38;
//complex<R> d39 = spa23*spa34*t5; d39 = R(1)/d39;
//complex<R> d40 = spa34*t35*R(2); d40 = R(1)/d40;
//complex<R> t16 = d11*t35*t48 + d12*t35*t48 + d7*spa35*t63 + d8*t49*t63;
//complex<R> t17 = d30*spb34*t21 + d23*s34*spa25*t35 + d21*s34*t77;
//complex<R> t18 = s23*(-(d21*s34*spa25*spa35) + d36*spb34*t21 + d35*s34*spa25*t35);
//complex<R> t43 = d16*spb34;
//complex<R> t62 = t52*t53;
//complex<R> t66 = -t72;
//complex<R> t70 = t57*t58;
//complex<R> t73 = d24*spa13;
//complex<R> t79 = spa15*t57;
//complex<R> t1 = -(d17*t33) + d10*t35*t48 + d9*t35*t48 + d19*t46*t52 + d14*d15*t24*t54 + t24*t43*t54 + d20*t62 + d4*t62 + d5*t62 + d6*t49*t63 + d3*t70 - d2*spa35*t79 + d18*t49*t79 + d13*d15*spa15*t72*R(2);
//complex<R> t2 = d10*t24*t48 + d11*t24*t48 + d12*t24*t48 + d9*t24*t48 + d4*t46*t52 + d5*t46*t52 + d14*d15*t35*t54 + t35*t43*t54 - d7*spa35*t63 + d13*d15*t26*t72 + d6*t22*t77 + d8*t22*t77 + d2*spa35*t79 + d3*t49*t79 + d1*t33*R(2);
//complex<R> t12 = d20*t46*t52 + d19*t62 + d18*t70;
//complex<R> t13 = s15*(d22*t21 + d23*spa25*t24 + d21*spa25*t49 + d26*t62 + d25*t70 + t66*t73);
//complex<R> t14 = -(d22*s23*t21) + d23*s23*spa25*t35 + d29*spb23*t62 + d27*spa13*spb23*t66 + d28*spb23*t70 + d21*s23*t77;
//complex<R> t19 = d26*s45*t62 + d25*s45*t70 + s45*t66*t73 - d31*t86;
//complex<R> t45 = d38*t42*t62 + d37*spa13*t30*t72 + d25*spa35*t30*t79 - d39*s15*t86;
//complex<R> co1 = d32*spb12*spb23*t33;
//complex<R> co2 = d33*spb23*spb34*t33;
//complex<R> co3 = -(d32*spb12*spb23*t33);
//complex<R> co4 = d34*spb34*t86;
//complex<R> co5 = -(d40*s15*spb12*t33);
//complex<R> co6 = -co1;
//complex<R> co7 = Complex(0,1);
//
//SeriesC<R> result = co7*(
//t2*(*CI_users[0]->get_value(ep,mu))
//+ t1*(*CI_users[1]->get_value(ep,mu))
//+ t16*(*CI_users[2]->get_value(ep,mu))
//+ t12*(*CI_users[3]->get_value(ep,mu))
//+ t13*(*CI_users[4]->get_value(ep,mu))
//+ t14*(*CI_users[5]->get_value(ep,mu))
//+ t17*(*CI_users[6]->get_value(ep,mu))
//+ t19*(*CI_users[7]->get_value(ep,mu))
//+ co1*(*CI_users[8]->get_value(ep,mu))
//+ co2*(*CI_users[9]->get_value(ep,mu))
//+ co6*(*CI_users[10]->get_value(ep,mu))
//+ co4*(*CI_users[11]->get_value(ep,mu))
//+ t18*(*CI_users[12]->get_value(ep,mu))
//+ t45*(*CI_users[13]->get_value(ep,mu))
//+ co5*(*CI_users[14]->get_value(ep,mu))
//);
// return(result);
//}
//
//
//
//C2q2g1y_qpgapqmmp_L_wCI::C2q2g1y_qpgapqmmp_L_wCI(const vector<int>& ind){
//
//
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpgapqmmp L");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Tue 18 Nov 2008 00:26:33 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//CI_users.push_back(new Cached_Bubble_Integral_User(c45,c123));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c2,c345));
//CI_users.push_back(new Cached_Triangle_Integral_User(c1,c5,c234));
//CI_users.push_back(new Cached_Triangle_Integral_User(c3,c4,c125));
//CI_users.push_back(new Cached_Triangle_Integral_User(c4,c5,c123));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c3,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c1,c5,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c2,c3,c4,c15));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c2,c1,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c4,c5,c12));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c5,c1,c23));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c1,c2,c34));
//CI_users.push_back(new Cached_Box_Integral_User(c5,c4,c3,c12));
//
//}
//
//SeriesC<R> C2q2g1y_qpgapqmmp_L_wCI::eval(const eval_param<R>& ep,const T& mu){
//
//
// // #define TimeStamp "Tue 18 Nov 2008 00:26:33 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa12 = SPA(1,2);
//complex<R> spa13 = SPA(1,3);
//complex<R> spa15 = SPA(1,5);
//complex<R> spa23 = SPA(2,3);
//complex<R> spa34 = SPA(3,4);
//complex<R> spa35 = SPA(3,5);
//complex<R> spa14 = SPA(1,4);
//complex<R> spa45 = SPA(4,5);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb15 = SPB(1,5);
//complex<R> spb23 = SPB(2,3);
//complex<R> spb25 = SPB(2,5);
//complex<R> spb45 = SPB(4,5);
//complex<R> s34 = S(3,4);
//complex<R> s45 = -(spa45*spb45);
//complex<R> s23 = -(spa23*spb23);
//complex<R> s12 = -(spa12*spb12);
//complex<R> s15 = -(spa15*spb15);
//complex<R> t1 = spa12*spa15;
//complex<R> t6 = square(spa34);
//complex<R> t9 = square(spa13);
//complex<R> t10 = square(spb15);
//complex<R> t11 = -(spa14*spb45);
//complex<R> t24 = spa13*spa14;
//complex<R> t28 = spa34*spb15;
//complex<R> d4 = spa15*spa23*spa35; d4 = R(1)/d4;
//complex<R> d5 = (s12 + s15 - s34)*spa15*spa23; d5 = R(1)/d5;
//complex<R> d6 = (s12 + s15 - s34)*spa23; d6 = R(1)/d6;
//complex<R> d9 = spa15*spa45*R(2); d9 = R(1)/d9;
//complex<R> d10 = (s12 + s15 - s34)*spa23*R(2); d10 = R(1)/d10;
//complex<R> d13 = spa12*spa23*R(2); d13 = R(1)/d13;
//complex<R> d14 = spa23*spa45*R(2); d14 = R(1)/d14;
//complex<R> t14 = spa14*t6;
//complex<R> t19 = d5*spb25;
//complex<R> t22 = s34*t6;
//complex<R> t29 = spa14*t9;
//complex<R> t32 = s12*t6;
//complex<R> d1 = spa23*spa45*t1*R(2); d1 = R(1)/d1;
//complex<R> d2 = (s23 - s45)*spa23*t1; d2 = R(1)/d2;
//complex<R> d3 = spa23*t1*square(s23 - s45)*R(2); d3 = R(1)/d3;
//complex<R> d7 = spa23*spa35*t1; d7 = R(1)/d7;
//complex<R> d8 = spa23*t1; d8 = R(1)/d8;
//complex<R> d11 = spa45*t1*R(2); d11 = R(1)/d11;
//complex<R> d12 = spa23*t1*R(2); d12 = R(1)/d12;
//complex<R> d15 = spa23*spa35*t1*R(2); d15 = R(1)/d15;
//complex<R> t5 = -(t19*t32) - d4*spa13*spb12*t6;
//complex<R> t15 = -(d2*R(2));
//complex<R> t16 = -(d7*spa13);
//complex<R> t18 = d3*spa45;
//complex<R> t4 = -(t10*t18*t29) + d2*t24*t28*R(2);
//complex<R> t20 = (d8*t11 + s45*t16)*t6;
//complex<R> t21 = t15*t24*t28 + t10*t18*t29 + d1*t14*R(3);
//complex<R> t26 = (t16 + t19)*t22;
//complex<R> co1 = d6*spb15*spb25*t6;
//complex<R> co2 = d9*spb12*spb23*t14;
//complex<R> co3 = d10*spb15*spb25*t32;
//complex<R> co4 = -(d11*s34*spb23*t14);
//complex<R> co5 = -(d9*spb12*spb23*t14);
//complex<R> co6 = d12*t11*t22;
//complex<R> co7 = d13*spb15*spb45*t14;
//complex<R> co8 = d14*spb12*spb15*t14;
//complex<R> co9 = -(d15*s45*spa13*t22);
//complex<R> co10 = -co2;
//complex<R> co11 = Complex(0,1);
//
//
//SeriesC<R> result = co11*(
//t21*(*CI_users[0]->get_value(ep,mu))
//+ t4*(*CI_users[1]->get_value(ep,mu))
//+ t5*(*CI_users[2]->get_value(ep,mu))
//+ co1*(*CI_users[3]->get_value(ep,mu))
//+ t26*(*CI_users[4]->get_value(ep,mu))
//+ t20*(*CI_users[5]->get_value(ep,mu))
//+ co2*(*CI_users[6]->get_value(ep,mu))
//+ co3*(*CI_users[7]->get_value(ep,mu))
//+ co4*(*CI_users[8]->get_value(ep,mu))
//+ co10*(*CI_users[9]->get_value(ep,mu))
//+ co6*(*CI_users[10]->get_value(ep,mu))
//+ co7*(*CI_users[11]->get_value(ep,mu))
//+ co8*(*CI_users[12]->get_value(ep,mu))
//+ co9*(*CI_users[13]->get_value(ep,mu))
//);
// return(result);
//}
//
//
//
//C2q2g1y_qpgapqmmm_L_wCI::C2q2g1y_qpgapqmmm_L_wCI(const vector<int>& ind){
//
//
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpgapqmmm L");
//#endif
//
//	//define corner vectors
//	 vector<int> c1;  c1.push_back(1-1);
//	 vector<int> c2;  c2.push_back(2-1);
//	 vector<int> c3;  c3.push_back(3-1);
//	 vector<int> c4;  c4.push_back(4-1);
//	 vector<int> c5;  c5.push_back(5-1);
//
//	 vector<int> c12;  c12.push_back(1-1); c12.push_back(2-1);
//	 vector<int> c23;  c23.push_back(2-1); c23.push_back(3-1);
//	 vector<int> c34;  c34.push_back(3-1); c34.push_back(4-1);
//	 vector<int> c45;  c45.push_back(4-1); c45.push_back(5-1);
//	 vector<int> c51;  c51.push_back(5-1); c51.push_back(1-1);
//         vector<int> c15 = c51;
//         vector<int> c14;  c14.push_back(1-1); c14.push_back(4-1);
//	 vector<int> c25;  c25.push_back(2-1); c25.push_back(5-1);
//	 vector<int> c35;  c35.push_back(3-1); c35.push_back(5-1);
//
//	 vector<int> c123; for(int i = 1; i<=3; i++){c123.push_back(i-1);}
//	 vector<int> c234; for(int i = 2; i<=4; i++){c234.push_back(i-1);}
//	 vector<int> c345; for(int i = 3; i<=5; i++){c345.push_back(i-1);}
//	 vector<int> c451; for(int i = 4; i<=5; i++){c451.push_back(i-1);}
//	                      c451.push_back(1-1);
//	 vector<int> c145 = c451;
//       	 vector<int> c512;    c512.push_back(5-1);
//                          for(int i = 1; i<=2; i++) {c512.push_back(i-1);}
//         vector<int> c125 = c512;
//       	 vector<int> c412;    c412.push_back(4-1);
//                          for(int i = 1; i<=2; i++) {c412.push_back(i-1);}
//         vector<int> c124 = c412;
//         vector<int> c134; c134.push_back(1-1); c134.push_back(3-1);
//                            c134.push_back(4-1);
//         vector<int> c235; c235.push_back(2-1); c235.push_back(3-1);
//                            c235.push_back(5-1);
// // #define TimeStamp "Tue 18 Nov 2008 00:26:34 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//
//
//CI_users.push_back(new Cached_Bubble_Integral_User(c15,c234));
//CI_users.push_back(new Cached_Bubble_Integral_User(c23,c145));
//CI_users.push_back(new Cached_Box_Integral_User(c1,c2,c3,c45));
//CI_users.push_back(new Cached_Box_Integral_User(c3,c4,c5,c12));
//CI_users.push_back(new Cached_Box_Integral_User(c4,c5,c1,c23));
//
//
//
//
//}
//
//SeriesC<R> C2q2g1y_qpgapqmmm_L_wCI::eval(const eval_param<R>& ep,const T& mu){
//
//
//
//
//#if _VERBOSE
//  _MESSAGE("C2q2g1y :  qpgapqmmm L");
//#endif
//
//
// // #define TimeStamp "Tue 18 Nov 2008 00:26:34 on n303"
// // #define Options "extract CSEs: True, shake: True {-2, -5}, simplifyFirst: False"
//#define Complex complex<R>
//complex<R> spa34 = SPA(3,4);
//complex<R> spa45 = SPA(4,5);
//complex<R> spb12 = SPB(1,2);
//complex<R> spb15 = SPB(1,5);
//complex<R> spa15 = SPA(1,5);
//complex<R> spb34 = SPB(3,4);
//complex<R> spb45 = SPB(4,5);
//complex<R> spb24 = SPB(2,4);
//complex<R> s15 = -(spa15*spb15);
//complex<R> s23 = S(2,3);
//complex<R> s12 = S(1,2);
//complex<R> t4 = spb34*spb45;
//complex<R> t5 = square(spb12);
//complex<R> t7 = square(spa45);
//complex<R> t8 = square(spb24);
//complex<R> t11 = spa45*spb12;
//complex<R> d4 = spb15*R(2); d4 = R(1)/d4;
//complex<R> d5 = spb34*R(2); d5 = R(1)/d5;
//complex<R> d1 = (s15 - s23)*t4; d1 = R(1)/d1;
//complex<R> d2 = t4*square(s15 - s23)*R(2); d2 = R(1)/d2;
//complex<R> d3 = spb15*t4*R(2); d3 = R(1)/d3;
//complex<R> t3 = d2*spb15*t7*t8 + d1*spb24*t11*R(2);
//complex<R> t16 = d3*t5;
//complex<R> t2 = -(d2*spb15*t7*t8) - d1*spb24*t11*R(2) - t16*R(3);
//complex<R> co1 = -(s12*s23*t16);
//complex<R> co2 = -(d4*spa34*spa45*t5);
//complex<R> co3 = -(d5*spa15*spa45*t5);
//complex<R> co4 = Complex(0,1);
//
//SeriesC<R> result = co4*(
//t3*(*CI_users[0]->get_value(ep,mu))
//+ t2*(*CI_users[1]->get_value(ep,mu))
//+ co1*(*CI_users[2]->get_value(ep,mu))
//+ co2*(*CI_users[3]->get_value(ep,mu))
//+ co3*(*CI_users[4]->get_value(ep,mu))
//);
// return(result);
//}
//
//
//C2q2g1y_qpgapqmpp_nf_wCI::C2q2g1y_qpgapqmpp_nf_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpgapqmpp_nf_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//
//C2q2g1y_qpgapqmpm_nf_wCI::C2q2g1y_qpgapqmpm_nf_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpgapqmpm_nf_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//
//C2q2g1y_qpgapqmmp_nf_wCI::C2q2g1y_qpgapqmmp_nf_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpgapqmmp_nf_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//C2q2g1y_qpgapqmmm_nf_wCI::C2q2g1y_qpgapqmmm_nf_wCI(const std::vector<int>& ind){
//
//}
//
//
//SeriesC<R> C2q2g1y_qpgapqmmm_nf_wCI::eval(const eval_param<R>& ep,const T& mu){
//    SeriesC<R> res(-2,0);
//    return res;
//
//}
//
//
//
//
// // *************** table of switch values *************
//
//#define _C_qpppqmgap_L C2q2g1y_6824_L
//#define _C_qppmqmgap_L C2q2g1y_6716_L
//#define _C_qpmpqmgap_L C2q2g1y_6806_L
//#define _C_qpmmqmgap_L C2q2g1y_6698_L
//#define _C_qppqmpgap_SLC C2q2g1y_7184_SLC
//#define _C_qppqmmgap_SLC C2q2g1y_6536_SLC
//#define _C_qpmqmpgap_SLC C2q2g1y_7166_SLC
//#define _C_qpmqmmgap_SLC C2q2g1y_6518_SLC
//#define _C_qpqmppgap_SLC C2q2g1y_7244_SLC
//#define _C_qpqmpmgap_SLC C2q2g1y_6596_SLC
//#define _C_qpqmmpgap_SLC C2q2g1y_7136_SLC
//#define _C_qpqmmmgap_SLC C2q2g1y_6488_SLC
//#define _C_qpppqmgap_nf C2q2g1y_6824_nf
//#define _C_qppmqmgap_nf C2q2g1y_6716_nf
//#define _C_qpmpqmgap_nf C2q2g1y_6806_nf
//#define _C_qpmmqmgap_nf C2q2g1y_6698_nf
//#define _C_qpgapqmpp_L C2q2g1y_4604_L
//#define _C_qpgapqmpm_L C2q2g1y_716_L
//#define _C_qpgapqmmp_L C2q2g1y_3956_L
//#define _C_qpgapqmmm_L C2q2g1y_68_L
//#define _C_qpgapqmpp_nf C2q2g1y_4604_nf
//#define _C_qpgapqmpm_nf C2q2g1y_716_nf
//#define _C_qpgapqmmp_nf C2q2g1y_3956_nf
//#define _C_qpgapqmmm_nf C2q2g1y_68_nf
//
//
// // *************** more macro definitions *************
//
//#define _CASE_qpppqmgap_L case 6824
//
//#define _CASE_qppmqmgap_L case 6716
//
//#define _CASE_qpmpqmgap_L case 6806
//
//#define _CASE_qpmmqmgap_L case 6698
//
//#define _CASE_qppqmpgap_SLC case 7184
//
//#define _CASE_qppqmmgap_SLC case 6536
//
//#define _CASE_qpmqmpgap_SLC case 7166
//
//#define _CASE_qpmqmmgap_SLC case 6518
//
//#define _CASE_qpqmppgap_SLC case 7244
//
//#define _CASE_qpqmpmgap_SLC case 6596
//
//#define _CASE_qpqmmpgap_SLC case 7136
//
//#define _CASE_qpqmmmgap_SLC case 6488
//
//#define _CASE_qpppqmgap_nf case 6824
//
//#define _CASE_qppmqmgap_nf case 6716
//
//#define _CASE_qpmpqmgap_nf case 6806
//
//#define _CASE_qpmmqmgap_nf case 6698
//
//#define _CASE_qpgapqmpp_L case 4604
//
//#define _CASE_qpgapqmpm_L case 716
//
//#define _CASE_qpgapqmmp_L case 3956
//
//#define _CASE_qpgapqmmm_L case 68
//
//#define _CASE_qpgapqmpp_nf case 4604
//
//#define _CASE_qpgapqmpm_nf case 716
//
//#define _CASE_qpgapqmmp_nf case 3956
//
//#define _CASE_qpgapqmmm_nf case 68
//
//
//
//
// // *************** define pointers *************
//
//
//Cut_Part_wCI*  CwCI_2q2g1y_L( int hc,const std::vector<int>& ind) {
//       switch (hc) {
//
//       _CASE_qpppqmgap_L: return new C2q2g1y_qpppqmgap_L_wCI(ind);
//       _CASE_qppmqmgap_L: return new C2q2g1y_qppmqmgap_L_wCI(ind);
//       _CASE_qpmpqmgap_L: return new C2q2g1y_qpmpqmgap_L_wCI(ind);
//       _CASE_qpmmqmgap_L: return new C2q2g1y_qpmmqmgap_L_wCI(ind);
//       _CASE_qpgapqmpp_L: return new C2q2g1y_qpgapqmpp_L_wCI(ind);
//       _CASE_qpgapqmpm_L: return new C2q2g1y_qpgapqmpm_L_wCI(ind);
//       _CASE_qpgapqmmp_L: return new C2q2g1y_qpgapqmmp_L_wCI(ind);
//       _CASE_qpgapqmmm_L: return new C2q2g1y_qpgapqmmm_L_wCI(ind);
//
//       default: return 0;
//        }
//}
//
//Cut_Part_wCI*  CwCI_2q2g1y_nf( int hc,const std::vector<int>& ind)  {
//       switch (hc) {
//       _CASE_qpppqmgap_nf: return new C2q2g1y_qpppqmgap_nf_wCI(ind);
//       _CASE_qppmqmgap_nf: return new C2q2g1y_qppmqmgap_nf_wCI(ind);
//       _CASE_qpmpqmgap_nf: return new C2q2g1y_qpmpqmgap_nf_wCI(ind);
//       _CASE_qpmmqmgap_nf: return new C2q2g1y_qpmmqmgap_nf_wCI(ind);
//       _CASE_qpgapqmpp_nf: return new C2q2g1y_qpgapqmpp_nf_wCI(ind);
//       _CASE_qpgapqmpm_nf: return new C2q2g1y_qpgapqmpm_nf_wCI(ind);
//       _CASE_qpgapqmmp_nf: return new C2q2g1y_qpgapqmmp_nf_wCI(ind);
//       _CASE_qpgapqmmm_nf: return new C2q2g1y_qpgapqmmm_nf_wCI(ind);
//
//       default: return 0;
//        }
// }
//
//Cut_Part_wCI*  CwCI_2q2g1y_SLC( int hc,const std::vector<int>& ind)  {
//       switch (hc) {
//       _CASE_qppqmpgap_SLC: return new C2q2g1y_qppqmpgap_SLC_wCI(ind);
//       _CASE_qppqmmgap_SLC: return new C2q2g1y_qppqmmgap_SLC_wCI(ind);
//       _CASE_qpmqmpgap_SLC: return new C2q2g1y_qpmqmpgap_SLC_wCI(ind);
//       _CASE_qpmqmmgap_SLC: return new C2q2g1y_qpmqmmgap_SLC_wCI(ind);
//       _CASE_qpqmppgap_SLC: return new C2q2g1y_qpqmppgap_SLC_wCI(ind);
//       _CASE_qpqmpmgap_SLC: return new C2q2g1y_qpqmpmgap_SLC_wCI(ind);
//       _CASE_qpqmmpgap_SLC: return new C2q2g1y_qpqmmpgap_SLC_wCI(ind);
//       _CASE_qpqmmmgap_SLC: return new C2q2g1y_qpqmmmgap_SLC_wCI(ind);
//
//       default: return 0;
//        }
// }
//
//
//// // *************** definitions for template *************
////
////template SeriesC<R> ( *C2q2g1y_L_Ptr_eval(int hc))
////             (const eval_param<R>&, const R&);
////template SeriesC<RHP> ( *C2q2g1y_L_Ptr_eval(int hc))
////             (const eval_param<RHP>&, const RHP&);
////template SeriesC<RVHP> ( *C2q2g1y_L_Ptr_eval(int hc))
////             (const eval_param<RVHP>&, const RVHP&);
////
////
////template SeriesC<R> ( *C2q2g1y_nf_Ptr_eval(int hc))
////             (const eval_param<R>&, const R&);
////template SeriesC<RHP> ( *C2q2g1y_nf_Ptr_eval(int hc))
////             (const eval_param<RHP>&, const RHP&);
////template SeriesC<RVHP> ( *C2q2g1y_nf_Ptr_eval(int hc))
////             (const eval_param<RVHP>&, const RVHP&);
////
////
////template SeriesC<R> ( *C2q2g1y_SLC_Ptr_eval(int hc))
////             (const eval_param<R>&, const R&);
////template SeriesC<RHP> ( *C2q2g1y_SLC_Ptr_eval(int hc))
////             (const eval_param<RHP>&, const RHP&);
////template SeriesC<RVHP> ( *C2q2g1y_SLC_Ptr_eval(int hc))
////             (const eval_param<RVHP>&, const RVHP&);
////
////
//
//}
//}
