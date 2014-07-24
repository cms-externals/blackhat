/* The 2 scalar amplitudes */

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;

namespace BH {

/*
 *
 *
 * A_eval function that returns zero
 *
 *
 */

template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses){return complex<T>(0,0);}


/*
 *
 *
 * The 2 scalar and a single gluon amplitudes
 *
 *
 */

template <int i0, int i1, int i2, class T> complex<T>  A2s1g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*(ep.ref()->L()*ep.p(i0)->Sm()*ep.p(i1)->Lt())/(ep.ref()->L()*ep.p(i1)->L());
}
template <int i0, int i1, int i2, class T> complex<T>  A2s1g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*(ep.p(i1)->L()*ep.p(i0)->Sm()*ep.ref()->Lt())/(ep.p(i1)->Lt()*ep.ref()->Lt());
}

template <class T> complex<T> (*A2s1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 46://s(0)+s(0)
		return &A2s1g1_eval<0,1,2>;
		break;
	case 34://s(0)-s(0)
		return &A2s1g2_eval<0,1,2>;
		break;
	case 58://s(0)s(0)+
		return &A2s1g1_eval<1,2,0>;
		break;
	case 10://s(0)s(0)-
		return &A2s1g2_eval<1,2,0>;
		break;
	case 43://+s(0)s(0)
		return &A2s1g1_eval<2,0,1>;
		break;
	case 40://-s(0)s(0)
		return &A2s1g2_eval<2,0,1>;
		break;
	default:// We return zero for all other helicity combinations
		return &ZeroF;
	}
}



/*
 *
 *
 * The 2 scalar and two gluon amplitudes
 *
 *
 */

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2s2g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,1)*mass*ep.spb(i1,i2)/(ep.spa(i1,i2)*(-T(2.)*ep.sp(i0,i1)));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2s2g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,1)*mass*ep.spa(i1,i2)/(ep.spb(i1,i2)*(-T(2.)*ep.sp(i0,i1)));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2s2g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,1)*pow(ep.spab(i2,i0,i1),2)/(ep.s(i1,i2)*(-T(2.)*ep.sp(i0,i1)));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2s2g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,1)*pow(ep.spba(i2,i0,i1),2)/(ep.s(i1,i2)*(-T(2.)*ep.sp(i0,i1)));
}

template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s2g5_eval(const eval_param<T>& ep, const mass_param_coll& masses, const complex<T>& mass){
	return complex<T>(0,1)*mass*pow(ep.spb(i1,i3),2)/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}
template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s2g6_eval(const eval_param<T>& ep, const mass_param_coll& masses, const complex<T>& mass){
	return complex<T>(0,1)*mass*pow(ep.spa(i1,i3),2)/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}


template <class T> complex<T> (*A2s2g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 190://s(0)++s(0)
		return &A2s2g1_eval<0,1,2,3>;
		break;
	case 250://s(0)s(0)++
		return &A2s2g1_eval<1,2,3,0>;
		break;
	case 235://+s(0)s(0)+
		return &A2s2g1_eval<2,3,0,1>;
		break;
	case 175://++s(0)s(0)
		return &A2s2g1_eval<3,0,1,2>;
		break;

	case 130://s(0)--s(0)
		return &A2s2g2_eval<0,1,2,3>;
		break;
	case 10://s(0)s(0)--
		return &A2s2g2_eval<1,2,3,0>;
		break;
	case 160://--s(0)s(0)
		return &A2s2g2_eval<3,0,1,2>;
		break;
	case 40://-s(0)s(0)-
		return &A2s2g2_eval<2,3,0,1>;
		break;

	case 142://s(0)+-s(0)
		return &A2s2g3_eval<0,1,2,3>;
		break;
	case 58://s(0)s(0)+-
		return &A2s2g3_eval<1,2,3,0>;
		break;
	case 163://+-s(0)s(0)
		return &A2s2g3_eval<3,0,1,2>;
		break;
	case 232://-s(0)s(0)+
		return &A2s2g3_eval<2,3,0,1>;
		break;

	case 178://s(0)-+s(0)
		return &A2s2g4_eval<0,1,2,3>;
		break;
	case 43://+s(0)s(0)-
		return &A2s2g4_eval<2,3,0,1>;
		break;
	case 202://s(0)s(0)-+
		return &A2s2g4_eval<1,2,3,0>;
		break;
	case 172://-+s(0)s(0)
		return &A2s2g4_eval<3,0,1,2>;
		break;
	default:// We return zero for all other helicity combinations
		return 0;
	}
}


/*
 *
 *
 * The 2 scalar and three gluon amplitudes
 *
 *
 */

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2s3g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,-1)*(pow(ep.spba(i1,i0,i3)*ep.spba(i3,i4,i3)+ep.spba(i1,i0,i2)*ep.spba(i2,i4,i3),2)
			/((-T(2)*ep.sp(i4,i3))*ep.spb(i3,i2)*ep.spb(i2,i1)*(-T(2)*ep.sp(i0,i1))*(ep.spa(i1,i3)*ep.spba(i3,i4,i3)+ep.spa(i1,i2)*ep.spba(i2,i4,i3)))
			-mass*pow(ep.spa(i3,i2),3)/(ep.s(i3,i2,i1)*ep.spa(i2,i1)*(ep.spa(i1,i3)*ep.spba(i3,i4,i3)+ep.spa(i1,i2)*ep.spba(i2,i4,i3))));}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2s3g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,1)*(pow(ep.spba(i2,i0,i1)*ep.spba(i2,i4,i3),2)
			/((-T(2)*ep.sp(i0,i1))*ep.spb(i1,i2)*ep.spb(i2,i3)*(-T(2)*ep.sp(i4,i3))*(ep.spa(i3,i1)*ep.spba(i1,i0,i1)+ep.spa(i3,i2)*ep.spba(i2,i0,i1)))
			-mass*pow(ep.spa(i1,i3),4)/(ep.s(i1,i2,i3)*ep.spa(i1,i2)*ep.spa(i2,i3)*(ep.spa(i3,i1)*ep.spba(i1,i0,i1)+ep.spa(i3,i2)*ep.spba(i2,i0,i1))));}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2s3g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,-1)*(pow(ep.spab(i3,i4,i1)*ep.spab(i1,i0,i1)+ep.spab(i3,i4,i2)*ep.spab(i2,i0,i1),2)
			/((-T(2)*ep.sp(i0,i1))*ep.spa(i1,i2)*ep.spa(i2,i3)*(-T(2)*ep.sp(i4,i3))*(ep.spb(i3,i1)*ep.spab(i1,i0,i1)+ep.spb(i3,i2)*ep.spab(i2,i0,i1)))
			-mass*pow(ep.spb(i1,i2),3)/(ep.s(i1,i2,i3)*ep.spb(i2,i3)*(ep.spb(i3,i1)*ep.spab(i1,i0,i1)+ep.spb(i3,i2)*ep.spab(i2,i0,i1))));}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2s3g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,1)*(pow(ep.spba(i3,i4,i1)*ep.spba(i1,i0,i1)+ep.spba(i3,i4,i2)*ep.spba(i2,i0,i1),2)
			/((-T(2)*ep.sp(i0,i1))*ep.spb(i1,i2)*ep.spb(i2,i3)*(-T(2)*ep.sp(i4,i3))*(ep.spa(i3,i1)*ep.spba(i1,i0,i1)+ep.spa(i3,i2)*ep.spba(i2,i0,i1)))
			-mass*pow(ep.spa(i1,i2),3)/(ep.s(i1,i2,i3)*ep.spa(i2,i3)*(ep.spa(i3,i1)*ep.spba(i1,i0,i1)+ep.spa(i3,i2)*ep.spba(i2,i0,i1))));}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2s3g5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,-1)*(pow(ep.spab(i2,i0,i1)*ep.spab(i2,i4,i3),2)
			/((-T(2)*ep.sp(i0,i1))*ep.spa(i1,i2)*ep.spa(i2,i3)*(-T(2)*ep.sp(i4,i3))*(ep.spb(i3,i1)*ep.spab(i1,i0,i1)+ep.spb(i3,i2)*ep.spab(i2,i0,i1)))
			-mass*pow(ep.spb(i1,i3),4)/(ep.s(i1,i2,i3)*ep.spb(i1,i2)*ep.spb(i2,i3)*(ep.spb(i3,i1)*ep.spab(i1,i0,i1)+ep.spb(i3,i2)*ep.spab(i2,i0,i1))));}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2s3g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,1)*(pow(ep.spab(i1,i0,i3)*ep.spab(i3,i4,i3)+ep.spab(i1,i0,i2)*ep.spab(i2,i4,i3),2)
			/((-T(2)*ep.sp(i4,i3))*ep.spa(i3,i2)*ep.spa(i2,i1)*(-T(2)*ep.sp(i0,i1))*(ep.spb(i1,i3)*ep.spab(i3,i4,i3)+ep.spb(i1,i2)*ep.spab(i2,i4,i3)))
			-mass*pow(ep.spb(i3,i2),3)/(ep.s(i3,i2,i1)*ep.spb(i2,i1)*(ep.spb(i1,i3)*ep.spab(i3,i4,i3)+ep.spb(i1,i2)*ep.spab(i2,i4,i3))));}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2s3g7_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,1)*mass*(ep.spb(i3,i1)*ep.spab(i1,i0,i1)+ep.spb(i3,i2)*ep.spab(i2,i0,i1))
			/((-T(2)*ep.sp(i0,i1))*ep.spa(i1,i2)*ep.spa(i2,i3)*(-T(2)*ep.sp(i4,i3)));}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2s3g8_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,-1)*mass*(ep.spa(i3,i1)*ep.spba(i1,i0,i1)+ep.spa(i3,i2)*ep.spba(i2,i0,i1))
			/((-T(2)*ep.sp(i0,i1))*ep.spb(i1,i2)*ep.spb(i2,i3)*(-T(2)*ep.sp(i4,i3)));}


template <class T> complex<T> (*A2s3g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 526://s(0)+--s(0)
		return &A2s3g1_eval<0,1,2,3,4>;
		break;
	case 232://-s(0)s(0)+-
		return &A2s3g1_eval<2,3,4,0,1>;
		break;
	case 928://--s(0)s(0)+
		return &A2s3g1_eval<3,4,0,1,2>;
		break;
	case 643://+--s(0)s(0)
		return &A2s3g1_eval<4,0,1,2,3>;
		break;
	case 58://s(0)s(0)+--
		return &A2s3g1_eval<1,2,3,4,0>;
		break;

	case 562://s(0)-+-s(0)
		return &A2s3g2_eval<0,1,2,3,4>;
		break;
	case 202://s(0)s(0)-+-
		return &A2s3g2_eval<1,2,3,4,0>;
		break;
	case 808://-s(0)s(0)-+
		return &A2s3g2_eval<2,3,4,0,1>;
		break;
	case 163://+-s(0)s(0)-
		return &A2s3g2_eval<3,4,0,1,2>;
		break;
	case 652://-+-s(0)s(0)
		return &A2s3g2_eval<4,0,1,2,3>;
		break;

	case 574://s(0)++-s(0)
		return &A2s3g3_eval<0,1,2,3,4>;
		break;
	case 931://+-s(0)s(0)+
		return &A2s3g3_eval<3,4,0,1,2>;
		break;
	case 1000://-s(0)s(0)++
		return &A2s3g3_eval<2,3,4,0,1>;
		break;
	case 250://s(0)s(0)++-
		return &A2s3g3_eval<1,2,3,4,0>;
		break;
	case 655://++-s(0)s(0)
		return &A2s3g3_eval<4,0,1,2,3>;
		break;

	case 706://s(0)--+s(0)
		return &A2s3g4_eval<0,1,2,3,4>;
		break;
	case 172://-+s(0)s(0)-
		return &A2s3g4_eval<3,4,0,1,2>;
		break;
	case 688://--+s(0)s(0)
		return &A2s3g4_eval<4,0,1,2,3>;
		break;
	case 43://+s(0)s(0)--
		return &A2s3g4_eval<2,3,4,0,1>;
		break;
	case 778://s(0)s(0)--+
		return &A2s3g4_eval<1,2,3,4,0>;
		break;

	case 718://s(0)+-+s(0)
		return &A2s3g5_eval<0,1,2,3,4>;
		break;
	case 235://+s(0)s(0)+-
		return &A2s3g5_eval<2,3,4,0,1>;
		break;
	case 691://+-+s(0)s(0)
		return &A2s3g5_eval<4,0,1,2,3>;
		break;
	case 826://s(0)s(0)+-+
		return &A2s3g5_eval<1,2,3,4,0>;
		break;
	case 940://-+s(0)s(0)+
		return &A2s3g5_eval<3,4,0,1,2>;
		break;

	case 754://s(0)-++s(0)
		return &A2s3g6_eval<0,1,2,3,4>;
		break;
	case 175://++s(0)s(0)-
		return &A2s3g6_eval<3,4,0,1,2>;
		break;
	case 970://s(0)s(0)-++
		return &A2s3g6_eval<1,2,3,4,0>;
		break;
	case 811://+s(0)s(0)-+
		return &A2s3g6_eval<2,3,4,0,1>;
		break;
	case 700://-++s(0)s(0)
		return &A2s3g6_eval<4,0,1,2,3>;
		break;

	case 766://s(0)+++s(0)
		return &A2s3g7_eval<0,1,2,3,4>;
		break;
	case 1003://+s(0)s(0)++
		return &A2s3g7_eval<2,3,4,0,1>;
		break;
	case 703://+++s(0)s(0)
		return &A2s3g7_eval<4,0,1,2,3>;
		break;
	case 1018://s(0)s(0)+++
		return &A2s3g7_eval<1,2,3,4,0>;
		break;
	case 943://++s(0)s(0)+
		return &A2s3g7_eval<3,4,0,1,2>;
		break;

	case 514://s(0)---s(0)
		return &A2s3g8_eval<0,1,2,3,4>;
		break;
	case 10://s(0)s(0)---
		return &A2s3g8_eval<1,2,3,4,0>;
		break;
	case 160://--s(0)s(0)-
		return &A2s3g8_eval<3,4,0,1,2>;
		break;
	case 640://---s(0)s(0)
		return &A2s3g8_eval<4,0,1,2,3>;
		break;
	case 40://-s(0)s(0)--
		return &A2s3g8_eval<2,3,4,0,1>;
		break;
	default:// We return zero for all other helicity combinations
		return 0;
	}
}


/*
 *
 *
 * The 2 scalar and four gluon amplitudes
 *
 *
 */

template <class T> complex<T>  A2s4g1278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses[0]);

	const complex<T> Q1=(-T(2)*ep.sp(0,1));
	const complex<T> Q2=(-T(2)*(ep.sp(0,1)+ep.sp(0,2)+ep.sp(1,2)));
	const complex<T> Q3=(-T(2)*ep.sp(5,4));
	const complex<T> denfac1=ep.spba(4,5,3)*ep.spb(3,1)*ep.spab(1,0,1)
									+ep.spba(4,5,3)*ep.spb(3,2)*ep.spab(2,0,1)
									+ep.spba(4,5,4)*ep.spb(4,1)*ep.spab(1,0,1)
									+ep.spba(4,5,4)*ep.spb(4,2)*ep.spab(2,0,1);
	const complex<T> denfac2=ep.spa(2,3)*ep.spbb(3,0,0,1)+ep.spa(2,3)*ep.spbb(3,5,0,1)
							+ep.spa(2,4)*ep.spbb(4,0,0,1)+ep.spa(2,4)*ep.spbb(4,5,0,1);

	const complex<T> Q1f=(-T(2)*ep.sp(5,4));
	const complex<T> Q2f=(-T(2)*(ep.sp(5,4)+ep.sp(5,3)+ep.sp(4,3)));
	const complex<T> Q3f=(-T(2)*ep.sp(0,1));
	const complex<T> denfac1f=ep.spba(1,0,2)*ep.spb(2,4)*ep.spab(4,5,4)
												+ep.spba(1,0,2)*ep.spb(2,3)*ep.spab(3,5,4)
												+ep.spba(1,0,1)*ep.spb(1,4)*ep.spab(4,5,4)
												+ep.spba(1,0,1)*ep.spb(1,3)*ep.spab(3,5,4);
	const complex<T> denfac2f=ep.spa(3,2)*ep.spbb(2,5,5,4)+ep.spa(3,2)*ep.spbb(2,0,5,4)
										+ep.spa(3,1)*ep.spbb(1,5,5,4)+ep.spa(3,1)*ep.spbb(1,0,5,4);


		return complex<T>(0,1)*(pow(Q2*(ep.spab(4,5,1)*ep.spab(1,0,1)+ep.spab(4,5,2)*ep.spab(2,0,1)+ep.spab(4,5,3)*ep.spab(3,0,1))
										-mass*ep.spab(4,5,3)*ep.spa(3,2)*ep.spb(2,1),2)
											/(Q1*Q2*Q3*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*denfac1)
								+mass*pow(ep.spb(3,1)*ep.spab(1,0,1)+ep.spb(3,2)*ep.spab(2,0,1),3)
									/(Q1*ep.spa(1,2)*ep.spb(3,4)*denfac1*denfac2)
								-mass*pow(ep.spab(4,2,1)+ep.spab(4,3,1),3)
									/(ep.s(1,2,3,4)*ep.s(2,3,4)*ep.spa(2,3)*ep.spa(3,4)*denfac2));
}
template <class T> complex<T>  A2s4g1854_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses[0]);

	const complex<T> Q1=(-T(2)*ep.sp(0,1));
	const complex<T> Q2=(-T(2)*(ep.sp(0,1)+ep.sp(0,2)+ep.sp(1,2)));
	const complex<T> Q3=(-T(2)*ep.sp(5,4));
	const complex<T> denfac1=ep.spba(4,5,3)*ep.spb(3,1)*ep.spab(1,0,1)
										+ep.spba(4,5,3)*ep.spb(3,2)*ep.spab(2,0,1)
										+ep.spba(4,5,4)*ep.spb(4,1)*ep.spab(1,0,1)
										+ep.spba(4,5,4)*ep.spb(4,2)*ep.spab(2,0,1);
	const complex<T> denfac2=ep.spa(2,3)*ep.spbb(3,0,0,1)+ep.spa(2,3)*ep.spbb(3,5,0,1)
								+ep.spa(2,4)*ep.spbb(4,0,0,1)+ep.spa(2,4)*ep.spbb(4,5,0,1);
	const complex<T> denfac3=ep.spb(3,1)*ep.spab(1,0,1)+ep.spb(3,2)*ep.spab(2,0,1);
	const complex<T> denfac4=ep.spab(4,2,1)+ep.spab(4,3,1);


	const complex<T> Q1f=(-T(2)*ep.sp(5,4));
	const complex<T> Q2f=(-T(2)*(ep.sp(5,4)+ep.sp(5,3)+ep.sp(4,3)));
	const complex<T> Q3f=(-T(2)*ep.sp(0,1));
	const complex<T> denfac1f=ep.spba(1,0,2)*ep.spb(2,4)*ep.spab(4,5,4)
													+ep.spba(1,0,2)*ep.spb(2,3)*ep.spab(3,5,4)
													+ep.spba(1,0,1)*ep.spb(1,4)*ep.spab(4,5,4)
													+ep.spba(1,0,1)*ep.spb(1,3)*ep.spab(3,5,4);
	const complex<T> denfac2f=ep.spa(3,2)*ep.spbb(2,5,5,4)+ep.spa(3,2)*ep.spbb(2,0,5,4)
											+ep.spa(3,1)*ep.spbb(1,5,5,4)+ep.spa(3,1)*ep.spbb(1,0,5,4);
	const complex<T> denfac3f=ep.spb(2,4)*ep.spab(4,5,4)+ep.spb(2,3)*ep.spab(3,5,4);
	const complex<T> denfac4f=ep.spab(1,3,4)+ep.spab(1,2,4);



	return complex<T>(0,1)*(pow((Q2*ep.spab(3,0,1)-mass*ep.spb(3,2)*ep.spb(2,1))*ep.spab(3,5,4),2)
										/(Q1*Q2*Q3*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*denfac1)
									+mass*pow(ep.spb(4,1)*ep.spab(1,0,1)+ep.spb(4,2)*ep.spab(2,0,1),4)
										/(Q1*ep.spa(1,2)*ep.spb(3,4)*denfac1*denfac2*denfac3)
									-mass*pow(ep.spab(3,2,1)+ep.spab(3,4,1),4)
										/(ep.s(1,2,3,4)*ep.s(2,3,4)*ep.spa(2,3)*ep.spa(3,4)*denfac2*denfac4)
									+mass*pow(ep.spb(1,2),3)*(ep.spb(4,1)*ep.spab(1,0,1)+ep.spb(4,2)*ep.spab(2,0,1)+ep.spb(4,3)*ep.spab(3,0,1))
										/(ep.s(1,2,3)*ep.spb(2,3)*Q3*denfac3*denfac4));
}
template <class T> complex<T>  A2s4g1998_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses[0]);

	const complex<T> Q1=(-T(2)*ep.sp(5,4));
	const complex<T> Q2=(-T(2)*(ep.sp(5,4)+ep.sp(5,3)+ep.sp(4,3)));
	const complex<T> Q3=(-T(2)*ep.sp(0,1));
	const complex<T> denfac1=ep.spba(1,0,2)*ep.spb(2,4)*ep.spab(4,5,4)
											+ep.spba(1,0,2)*ep.spb(2,3)*ep.spab(3,5,4)
											+ep.spba(1,0,1)*ep.spb(1,4)*ep.spab(4,5,4)
											+ep.spba(1,0,1)*ep.spb(1,3)*ep.spab(3,5,4);
	const complex<T> denfac2=ep.spa(3,2)*ep.spbb(2,5,5,4)+ep.spa(3,2)*ep.spbb(2,0,5,4)
									+ep.spa(3,1)*ep.spbb(1,5,5,4)+ep.spa(3,1)*ep.spbb(1,0,5,4);
	const complex<T> denfac3=ep.spb(2,4)*ep.spab(4,5,4)+ep.spb(2,3)*ep.spab(3,5,4);
	const complex<T> denfac4=ep.spab(1,3,4)+ep.spab(1,2,4);
	return complex<T>(0,1)*(pow((Q2*ep.spab(2,5,4)-mass*ep.spb(2,3)*ep.spb(3,4))*ep.spab(2,0,1),2)
											/(Q1*Q2*Q3*ep.spa(4,3)*ep.spa(3,2)*ep.spa(2,1)*denfac1)
										+mass*pow(ep.spb(1,4)*ep.spab(4,5,4)+ep.spb(1,3)*ep.spab(3,5,4),4)
											/(Q1*ep.spa(4,3)*ep.spb(2,1)*denfac1*denfac2*denfac3)
										-mass*pow(ep.spab(2,3,4)+ep.spab(2,1,4),4)
											/(ep.s(4,3,2,1)*ep.s(3,2,1)*ep.spa(3,2)*ep.spa(2,1)*denfac2*denfac4)
										+mass*pow(ep.spb(4,3),3)*(ep.spb(1,4)*ep.spab(4,5,4)+ep.spb(1,3)*ep.spab(3,5,4)+ep.spb(1,2)*ep.spab(2,5,4))
											/(ep.s(4,3,2)*ep.spb(3,2)*Q3*denfac3*denfac4));
}
template <class T> complex<T>  A2s4g2034_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses[0]);

	const complex<T> Q1=(-T(2)*ep.sp(5,4));
	const complex<T> Q2=(-T(2)*(ep.sp(5,4)+ep.sp(5,3)+ep.sp(4,3)));
	const complex<T> Q3=(-T(2)*ep.sp(0,1));
	const complex<T> denfac1=ep.spba(1,0,2)*ep.spb(2,4)*ep.spab(4,5,4)
										+ep.spba(1,0,2)*ep.spb(2,3)*ep.spab(3,5,4)
										+ep.spba(1,0,1)*ep.spb(1,4)*ep.spab(4,5,4)
										+ep.spba(1,0,1)*ep.spb(1,3)*ep.spab(3,5,4);
	const complex<T> denfac2=ep.spa(3,2)*ep.spbb(2,5,5,4)+ep.spa(3,2)*ep.spbb(2,0,5,4)
								+ep.spa(3,1)*ep.spbb(1,5,5,4)+ep.spa(3,1)*ep.spbb(1,0,5,4);
	return complex<T>(0,1)*(pow(Q2*(ep.spab(1,0,4)*ep.spab(4,5,4)+ep.spab(1,0,3)*ep.spab(3,5,4)+ep.spab(1,0,2)*ep.spab(2,5,4))
											-mass*ep.spab(1,0,2)*ep.spa(2,3)*ep.spb(3,4),2)
												/(Q1*Q2*Q3*ep.spa(4,3)*ep.spa(3,2)*ep.spa(2,1)*denfac1)
									+mass*pow(ep.spb(2,4)*ep.spab(4,5,4)+ep.spb(2,3)*ep.spab(3,5,4),3)
										/(Q1*ep.spa(4,3)*ep.spb(2,1)*denfac1*denfac2)
									-mass*pow(ep.spab(1,3,4)+ep.spab(1,2,4),3)
										/(ep.s(4,3,2,1)*ep.s(3,2,1)*ep.spa(3,2)*ep.spa(2,1)*denfac2));
}
template <int i0, int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2s4g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	return complex<T>(0,1)*mass*(ep.spb(i4,i1)*ep.spab(i1,i0,i1)+ep.spb(i4,i2)*ep.spab(i2,i0,i1)+ep.spb(i4,i3)*ep.spab(i3,i0,i1)
						-mass*ep.spb(i1,i2)*ep.spa(i2,i3)*ep.spb(i3,i4)/(-T(2)*(ep.sp(i0,i1)+ep.sp(i0,i2)+ep.sp(i1,i2))))
						/((-T(2)*ep.sp(i0,i1))*(-T(2)*ep.sp(i5,i4))*ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4));
}

template <class T> complex<T> (*A2s4g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 3070://s(0)++++s(0)
		return &A2s4g1_eval<0,1,2,3,4,5>;
	break;
	case 2815://++++s(0)s(0)
		return &A2s4g1_eval<5,0,1,2,3,4>;
		break;
	case 3775://+++s(0)s(0)+
		return &A2s4g1_eval<4,5,0,1,2,3>;
		break;
	case 4015://++s(0)s(0)++
		return &A2s4g1_eval<3,4,5,0,1,2>;
		break;
	case 4075://+s(0)s(0)+++
		return &A2s4g1_eval<2,3,4,5,0,1>;
		break;
	case 4090://s(0)s(0)++++
		return &A2s4g1_eval<1,2,3,4,5,0>;
		break;
	default:// We return zero for all other helicity combinations
		return 0;//We do not know the amplitude and so we return null
	}
}


template complex<R> (*A2s1g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2s1g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2s1g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2s2g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2s2g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2s2g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2s3g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2s3g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2s3g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2s4g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2s4g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2s4g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> (*A2s1g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A2s2g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A2s3g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A2s4g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif


}
