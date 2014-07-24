/*
 * A0_2q_massive_eval.cpp
 *
 *  Created on: 14-Jul-2008
 *      Author: Darren
 */

/* The 2 Massive Quarks and gluons amplitudes */

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;

namespace BH {

template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses){return complex<T>(0,0);}

/*
 *
 *
 * The 2 massive quarks and a single gluon amplitudes
 *
 *
 */

template <int i0, int i1, int i2, class T> complex<T>  A2QMg1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);//const complex<T> mass=ep.ref_mass()/T(2);

	const lambdat<T> p1Lt=lat(ep.mom(i0)-(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> pnLt=lat(ep.mom(i2)-(mass/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> i1Lat=ep.p(i1)->Lt();

	return complex<T>(0,1)*pow((i1Lat*p1Lt),2)/(p1Lt*pnLt);
}

template <int i0, int i1, int i2, class T> complex<T>  A2QMg2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);//const complex<T> mass=ep.ref_mass()/T(2);

	const lambda<T> p1L=la(ep.mom(i0)-(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> pnL=la(ep.mom(i2)-(mass/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> i1L=ep.p(i1)->L();

	return complex<T>(0,1)*pow((i1L*p1L),2)/(p1L*pnL);
}

template <int i0, int i1, int i2, class T> complex<T>  A2QMg3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);//const complex<T> mass=ep.ref_mass()/T(2);

	const lambdat<T> p1Lt=lat(ep.mom(i0)-(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> pnLt=lat(ep.mom(i2)-(mass/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return complex<T>(0,1)*pow(i1Lt*pnLt,2)/(p1Lt*pnLt);
}

template <int i0, int i1, int i2, class T> complex<T>  A2QMg4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);//const complex<T> mass=ep.ref_mass()/T(2);

	const lambda<T> p1L=la(ep.mom(i0)-(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> pnL=la(ep.mom(i2)-(mass/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> i1L=ep.p(i1)->L();

	return complex<T>(0,1)*pow(i1L*pnL,2)/(p1L*pnL);
}

template <int i0, int i1, int i2, class T> complex<T>  A2QMg5m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);//const complex<T> mass=ep.ref_mass()/T(2);

	const lambdat<T> p1Lt=lat(ep.mom(i0)-(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> pnLt=lat(ep.mom(i2)-(mass/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> i1Lt=ep.p(i1)->Lt();
	const lambdat<T> qLt=ep.ref()->Lt();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*pow(i1Lt*qLt,2)/((p1Lt*qLt)*(pnLt*qLt));
}
template <int i0, int i1, int i2, class T> complex<T>  A2QMg5p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);//const complex<T> mass=ep.ref_mass()/T(2);

	const lambdat<T> p1Lt=lat(ep.mom(i0)-(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> pnLt=lat(ep.mom(i2)-(mass/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> i1Lt=ep.p(i1)->Lt();
	const lambdat<T> qLt=ep.ref()->Lt();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i0))*pow(i1Lt*qLt,2)/((p1Lt*qLt)*(pnLt*qLt));
}

template <int i0, int i1, int i2, class T> complex<T>  A2QMg6m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);//const complex<T> mass=ep.ref_mass()/T(2);

	const lambda<T> p1L=la(ep.mom(i0)-(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> pnL=la(ep.mom(i2)-(mass/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> i1L=ep.p(i1)->L();
	const lambda<T> qL=ep.ref()->L();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*pow(i1L*qL,2)/((p1L*qL)*(pnL*qL));
}
template <int i0, int i1, int i2, class T> complex<T>  A2QMg6p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);//const complex<T> mass=ep.ref_mass()/T(2);

	const lambda<T> p1L=la(ep.mom(i0)-(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> pnL=la(ep.mom(i2)-(mass/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> i1L=ep.p(i1)->L();
	const lambda<T> qL=ep.ref()->L();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i0))*pow(i1L*qL,2)/((p1L*qL)*(pnL*qL));
}

template <class T> complex<T> (*A2QM1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	// case 2q_massive
	switch (hc) {
	case 0x617: //167: //(Qb+,g+,Q-)
	case 0x815:
		return &A2QMg1_eval<0,1,2>;
		break;
	case 0x761: //142: //(Q-,Qb+,g+)
	case 0x581:
		return &A2QMg1_eval<1,2,0>;
		break;
	case 0x176: //207: //(g+,Q-,Qb+)
	case 0x158:
		return &A2QMg1_eval<2,0,1>;
		break;

	case 0x508: //184: //(Qb-,g-,Q+)
	case 0x706:
		return &A2QMg2_eval<0,1,2>;
		break;
	case 0x850: //29: //(Q+,Qb-,g-)
	case 0x670:
		return &A2QMg2_eval<1,2,0>;
		break;
	case 0x085: //174: //(g-,Q+,Qb-)
	case 0x067:
		return &A2QMg2_eval<2,0,1>;
		break;

	case 0x518: //202: //(Qb-,g+,Q+)
	case 0x716:
		return &A2QMg3_eval<0,1,2>;
		break;
	case 0x851: //137: //(Q+,Qb-,g+)
	case 0x671:
		return &A2QMg3_eval<1,2,0>;
		break;
	case 0x185: //177: //(g+,Q+,Qb-)
	case 0x167:
		return &A2QMg3_eval<2,0,1>;
		break;

	case 0x805: //149: //(Qb+,g-,Q-)
	case 0x607:
		return &A2QMg4_eval<0,1,2>;
		break;
	case 0x580: //34: //(Qb+,g-,Q-)
	case 0x760:
		return &A2QMg4_eval<1,2,0>;
		break;
	case 0x058: //204: //(g-,Q-,Qb+)
	case 0x076:
		return &A2QMg4_eval<2,0,1>;
		break;

	case 0x715: //166: //(Qb-,g+,Q-)
		return &A2QMg5m_eval<0,1,2>;
		break;
	case 0x517:
		return &A2QMg5p_eval<0,1,2>;
		break;
	case 0x571: //136: //(Q-,Qb-,g+)
		return &A2QMg5m_eval<1,2,0>;
		break;
	case 0x751:
		return &A2QMg5p_eval<1,2,0>;
		break;
	case 0x157: //171: //(g+,Q-,Qb-)
		return &A2QMg5m_eval<2,0,1>;
		break;
	case 0x175:
		return &A2QMg5p_eval<2,0,1>;
		break;

	case 0x608: //185: //(Qb+,g-,Q+)
		return &A2QMg6p_eval<0,1,2>;
		break;
	case 0x806:
		return &A2QMg6m_eval<0,1,2>;
		break;
	case 0x860: //35: //(Q+,Qb+,g-)
		return &A2QMg6p_eval<1,2,0>;
		break;
	case 0x680:
		return &A2QMg6m_eval<1,2,0>;
		break;
	case 0x086: //210: //(g-,Q+,Qb+)
		return &A2QMg6p_eval<2,0,1>;
		break;
	case 0x068:
		return &A2QMg6m_eval<2,0,1>;
		break;

	case 0x618: //203://(Qb+,g+,Q+)
	case 0x816:
	case 0x861: //143://(Q+,Qb+,g+)
	case 0x681:
	case 0x186: //213://(g+,Q+,Qb+)
	case 0x168:
	case 0x507: //148://(Qb-,g-,Q-)
	case 0x705:
	case 0x750: //28://(Q-,Qb-,g-)
	case 0x570:
	case 0x075: //168://(g-,Q-,Qb-)
	case 0x057:
		return &ZeroF;
		break;

	default:// We return zero for all other helicity combinations
		cout << "3 pt A2QM1g_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return &ZeroF;
	}
}


/*
 *
 *
 * The 2 massive quarks and two gluon amplitudes
 *
 *
 */

template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s2g_help_1_eval(const eval_param<T>& ep, const complex<T>& mass){
	return complex<T>(0,1)*mass*ep.spb(i1,i2)/(ep.spa(i1,i2)*(-T(2)*ep.sp(i0,i1)));
}
template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s2g_help_2_eval(const eval_param<T>& ep, const complex<T>& mass){
	return complex<T>(0,1)*mass*ep.spa(i1,i2)/(ep.spb(i1,i2)*(-T(2)*ep.sp(i0,i1)));
}
template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s2g_help_3_eval(const eval_param<T>& ep, const complex<T>& mass){
	return complex<T>(0,1)*mass*pow(ep.spb(i1,i3),2)/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}
template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s2g_help_4_eval(const eval_param<T>& ep, const complex<T>& mass){
	return complex<T>(0,1)*mass*pow(ep.spa(i1,i3),2)/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g1_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> qL=ep.ref()->L();

	return ((pnL*qL)/(p1L*qL))*A2s2g_help_1_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g2_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> qL=ep.ref()->L();

	return ((p1L*qL)/(pnL*qL))*A2s2g_help_1_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g3m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -((p1L*pnL)/eval_param<T>::mass(masses.p(i0)))*A2s2g_help_1_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g3p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((p1L*pnL)/eval_param<T>::mass(masses.p(i0)))*A2s2g_help_1_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g4_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return complex<T>(0,-1)*ep.spba(i1,i3,i2)*
					(ep.spba(i1,i3,i2)*(qLt*p1Lt)/(T(2)*ep.sp(i1,i2))
							+(ep.ref()->Lt()*ep.p(i1)->Lt())*(i1Lt*p1Lt)/ep.spb(i1,i2))
					/((qLt*pnLt)*(-T(2)*ep.sp(i0,i1)));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g5_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return complex<T>(0.,-1)*ep.spba(i1,i3,i2)*
			(ep.spba(i1,i3,i2)*(qLt*pnLt)/(T(2)*ep.sp(i1,i2))
					-(ep.ref()->Lt()*ep.p(i1)->Lt())*(i1Lt*pnLt)/ep.spb(i1,i2))
			/((qLt*p1Lt)*(-T(2)*ep.sp(i0,i1)));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g6_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> i1L=ep.p(i1)->L();

	return complex<T>(0,1)*(ep.spab(i1,i3,i2)*
					(ep.spab(i1,i3,i2)*(qL*pnL)/(T(2)*ep.sp(i1,i2))
							-(ep.ref()->L()*ep.p(i1)->L())*(i1L*pnL)/ep.spa(i1,i2))
					/((qL*p1L)*(-T(2)*ep.sp(i0,i1))));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g7_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> i1L=ep.p(i1)->L();

	return complex<T>(0,1)*(ep.spab(i1,i3,i2)*
			(ep.spab(i1,i3,i2)*(qL*p1L)/(T(2)*ep.sp(i1,i2))
					+(ep.ref()->L()*ep.p(i1)->L())*(i1L*p1L)/ep.spa(i1,i2))
			/((qL*pnL)*(-T(2)*ep.sp(i0,i1))));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g8m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return eval_param<T>::mass(masses.p(i0))*ep.spba(i1,i0,i2)*pow(qLt*i1Lt,2)/(ep.sp(i0,i1)*ep.spb(i1,i2)*(qLt*p1Lt)*(pnLt*qLt)*complex<T>(0,-2));

}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g8p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return eval_param<T>::mass(masses.p(i0))*ep.spba(i1,i0,i2)*pow(qLt*i1Lt,2)/(ep.sp(i0,i1)*ep.spb(i1,i2)*(qLt*p1Lt)*(pnLt*qLt)*complex<T>(0,2));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g9m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();


	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i2Lt=ep.p(i2)->Lt();

	return eval_param<T>::mass(masses.p(i0))*ep.spba(i2,i3,i1)*pow(qLt*i2Lt,2)/(ep.sp(i3,i2)*ep.spb(i2,i1)*(qLt*pnLt)*(p1Lt*qLt)*complex<T>(0,2));

}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g9p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();


	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i2Lt=ep.p(i2)->Lt();

	return eval_param<T>::mass(masses.p(i0))*ep.spba(i2,i3,i1)*pow(qLt*i2Lt,2)/(ep.sp(i3,i2)*ep.spb(i2,i1)*(qLt*pnLt)*(p1Lt*qLt)*complex<T>(0,-2));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g10_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return -((p1Lt*qLt)/(pnLt*qLt))*A2s2g_help_2_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g11_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return -((pnLt*qLt)/(p1Lt*qLt))*A2s2g_help_2_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g12m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((p1Lt*pnLt)/eval_param<T>::mass(masses.p(i0)))*A2s2g_help_2_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g12p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -((p1Lt*pnLt)/eval_param<T>::mass(masses.p(i0)))*A2s2g_help_2_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g13m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();


	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> i2L=ep.p(i2)->L();

	return eval_param<T>::mass(masses.p(i0))*ep.spab(i2,i3,i1)*pow(qL*i2L,2)/(ep.sp(i3,i2)*ep.spa(i2,i1)*(qL*pnL)*(p1L*qL)*complex<T>(0,-2));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g13p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();


	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> i2L=ep.p(i2)->L();

	return eval_param<T>::mass(masses.p(i0))*ep.spab(i2,i3,i1)*pow(qL*i2L,2)/(ep.sp(i3,i2)*ep.spa(i2,i1)*(qL*pnL)*(p1L*qL)*complex<T>(0,2));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g14m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> i1L=ep.p(i1)->L();

	return eval_param<T>::mass(masses.p(i0))*ep.spab(i1,i0,i2)*pow(qL*i1L,2)/(ep.sp(i0,i1)*ep.spa(i1,i2)*(qL*p1L)*(pnL*qL)*complex<T>(0,2));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g14p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> i1L=ep.p(i1)->L();

	return eval_param<T>::mass(masses.p(i0))*ep.spab(i1,i0,i2)*pow(qL*i1L,2)/(ep.sp(i0,i1)*ep.spa(i1,i2)*(qL*p1L)*(pnL*qL)*complex<T>(0,-2));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g15_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();

	return ((p1L*qL)/(pnL*qL))*A2s2g_help_3_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g16_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -((p1L*pnL)/eval_param<T>::mass(masses.p(i0)))*A2s2g_help_3_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g17_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return -((p1Lt*qLt)/(pnLt*qLt))*A2s2g_help_4_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g18_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((p1Lt*pnLt)/eval_param<T>::mass(masses.p(i0)))*A2s2g_help_4_eval<i0,i1,i2,i3>(ep,mass);

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g19_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambda<T> qL=ep.ref()->L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();
	const lambda<T> i3L=ep.p(i3)->L();

	return complex<T>(0,-1)*ep.spab(i3,i0,i1)*((p1Lt*i1Lt)*(i3L*pnL)-mass*(ep.ref()->L()*ep.p(i3)->L())*(ep.p(i1)->Lt()*ep.ref()->Lt())/((qL*p1L)*(pnLt*qLt)))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g20_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambda<T> qL=ep.ref()->L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();
	const lambda<T> i3L=ep.p(i3)->L();

	return complex<T>(0,1)*ep.spab(i3,i0,i1)*((p1L*i3L)*(i1Lt*pnLt)-mass*(ep.ref()->L()*ep.p(i3)->L())*(ep.p(i1)->Lt()*ep.ref()->Lt())/((qLt*p1Lt)*pnL*qL))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g23m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambda<T> qL=ep.ref()->L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*ep.spab(i3,i0,i1)*((p1Lt*i1Lt)*(ep.p(i3)->L()*ep.ref()->L())/(pnL*qL)-(ep.ref()->L()*ep.p(i3)->L())*(i1Lt*pnLt)/(qL*p1L))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g23p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambda<T> qL=ep.ref()->L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i0))*ep.spab(i3,i0,i1)*((p1Lt*i1Lt)*(ep.p(i3)->L()*ep.ref()->L())/(pnL*qL)-(ep.ref()->L()*ep.p(i3)->L())*(i1Lt*pnLt)/(qL*p1L))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g24m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambda<T> i1L=ep.p(i1)->L();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i0))*ep.spba(i3,i0,i1)*((p1L*i1L)*(ep.p(i3)->Lt()*ep.ref()->Lt())/(pnLt*qLt)-(ep.ref()->Lt()*ep.p(i3)->Lt())*(i1L*pnL)/(qLt*p1Lt))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));

}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2g24p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambda<T> i1L=ep.p(i1)->L();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*ep.spba(i3,i0,i1)*((p1L*i1L)*(ep.p(i3)->Lt()*ep.ref()->Lt())/(pnLt*qLt)-(ep.ref()->Lt()*ep.p(i3)->Lt())*(i1L*pnL)/(qLt*p1Lt))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));

}


template <class T> complex<T> (*A2QM2g_Tree_Ptr_eval(int hc))(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// case 2q_massive
	switch (hc) {
	case 0x8115://995: // Qb+,g+,g+,Q-
	case 0x6117:
		return &A2QM2g1_eval<0,1,2,3>;
		break;
	case 0x5811://790: //(Q3-,Qb3+,g+,g+)
	case 0x7611:
		return &A2QM2g1_eval<1,2,3,0>;
		break;
	case 0x1581://855: // (g+,Qb3-,Q3+,g+)
	case 0x1761:
		return &A2QM2g1_eval<2,3,0,1>;
		break;
	case 0x1158://1245: // (g+,g+,Qb3-,Q3+)
	case 0x1176:
		return &A2QM2g1_eval<3,0,1,2>;
		break;

	case 0x7116://1210: // Qb-,g+,g+,Q+
	case 0x5118:
		return &A2QM2g2_eval<0,1,2,3>;
		break;
	case 0x6711://785: //(Q3+,Qb3-,g+,g+)
	case 0x8511:
		return &A2QM2g2_eval<1,2,3,0>;
		break;
	case 0x1671://825: // (g+,Qb3+,Q3-,g+)
	case 0x1851:
		return &A2QM2g2_eval<2,3,0,1>;
		break;
	case 0x1167://1065: // (g+,g+,Qb3+,Q3-)
	case 0x1185:
		return &A2QM2g2_eval<3,0,1,2>;
		break;

	case 0x5181: // (Q-,g+,Qb+,g+) NEW
	case 0x7161: // (Qb-,g+,Q+,g+) NEW
		return &A2QM2g15_eval<0,1,2,3>;
		break;
	case 0x1518: // (g+,Q-,g+,Qb+) NEW
	case 0x1716: // (g+,Qb-,g+,Q+) NEW
		return &A2QM2g15_eval<1,2,3,0>;
		break;
	case 0x8151: // (Qb+,g+,Q-,g+) NEW
	case 0x6171: // (Q+,g+,Qb-,g+) NEW
		return &A2QM2g15_eval<2,3,0,1>;
		break;
	case 0x1815: // (g+,Qb+,g+,Q-) NEW
	case 0x1617: // (g+,Q+,g+,Qb-) NEW
		return &A2QM2g15_eval<3,0,1,2>;
		break;

	case 0x5117: //994: // Qb-g+g+Q-
		return &A2QM2g3p_eval<0,1,2,3>;
		break;
	case 0x7115:
		return &A2QM2g3m_eval<0,1,2,3>;
		break;
	case 0x7511: //784: //(Q3-,Qb3-,g+,g+)
		return &A2QM2g3p_eval<1,2,3,0>;
		break;
	case 0x5711:
		return &A2QM2g3m_eval<1,2,3,0>;
		break;
	case 0x1751: //819: // (g+,Qb3-,Q3-,g+)
		return &A2QM2g3p_eval<2,3,0,1>;
		break;
	case 0x1571:
		return &A2QM2g3m_eval<2,3,0,1>;
		break;
	case 0x1175: //1029: // (g+,g+,Qb3-,Q3-)
		return &A2QM2g3p_eval<3,0,1,2>;
		break;
	case 0x1157:
		return &A2QM2g3m_eval<3,0,1,2>;
		break;

	case 0x7151: // (Qb-,g+,Q-,g+) NEW
		return &A2QM2g16_eval<0,1,2,3>;
		break;
	case 0x1715: // (g+,Qb-,g+,Q-) NEW
		return &A2QM2g16_eval<1,2,3,0>;
		break;
	case 0x5171: // (Q-,g+,Qb-,g+) NEW
		return &A2QM2g16_eval<2,3,0,1>;
		break;
	case 0x1517: // (g+,Q-,g+,Qb-) NEW
		return &A2QM2g16_eval<3,0,1,2>;
		break;

	case 0x8105: //887: // Qb+g+g-Q-
	case 0x6107:
		return &A2QM2g4_eval<0,1,2,3>;
		break;
	case 0x5810: //142: //(Q3-,Qb3+,g+,g-)
	case 0x7610:
		return &A2QM2g4_eval<1,2,3,0>;
		break;
	case 0x0581: //852: // (g-,Qb3-,Q3+,g+)
	case 0x0761:
		return &A2QM2g4_eval<2,3,0,1>;
		break;
	case 0x1058: //1227: //(g+,g-,Qb3-,Q3+)
	case 0x1076:
		return &A2QM2g4_eval<3,0,1,2>;
		break;

	case 0x8150: //887: // Qb+g+Q-g-
	case 0x6170:
		return &A2QM2g19_eval<0,1,2,3>;
		break;
	case 0x0815: //142: //(g-,Q3-,Qb3+,g+)
	case 0x0617:
		return &A2QM2g19_eval<1,2,3,0>;
		break;
	case 0x5081: //852: // (Qb3-,g-,Q3+,g+)
	case 0x7061:
		return &A2QM2g19_eval<2,3,0,1>;
		break;
	case 0x1508: //1227: //(g+,Qb3-,g-,Q3+)
	case 0x1706:
		return &A2QM2g19_eval<3,0,1,2>;
		break;

	case 0x7106: //1102: // Qb-g+g-Q+
	case 0x5108:
		return &A2QM2g5_eval<0,1,2,3>;
		break;
	case 0x6710: //137: // (Q3+,Qb3-,g+,g-)
	case 0x8510:
		return &A2QM2g5_eval<1,2,3,0>;
		break;
	case 0x0671: //
	case 0x0851:
		return &A2QM2g5_eval<2,3,0,1>;
		break;
	case 0x1067: //1047: // (g+,g-,Q3+,Qb3-)
	case 0x1085:
		return &A2QM2g5_eval<3,0,1,2>;
		break;

	case 0x7160: //1102: // Qb-g+Q+g-
	case 0x5180:
		return &A2QM2g20_eval<0,1,2,3>;
		break;
	case 0x0716: //137: // (Q3+,Qb3-,g+,g-)
	case 0x0518:
		return &A2QM2g20_eval<1,2,3,0>;
		break;
	case 0x6071: //
	case 0x8051:
		return &A2QM2g20_eval<2,3,0,1>;
		break;
	case 0x1607: //1047: // (g+,g-,Q3+,Qb3-)
	case 0x1805:
		return &A2QM2g20_eval<3,0,1,2>;
		break;

	case 0x8015: //977: // Qb+g-g+Q-
	case 0x6017:
		return &A2QM2g6_eval<0,1,2,3>;
		break;
	case 0x5801:
	case 0x7601:
		return &A2QM2g6_eval<1,2,3,0>;
		break;
	case 0x1580: //207: // (g+,Qb3-,Q3+,g-)
	case 0x1760:
		return &A2QM2g6_eval<2,3,0,1>;
		break;
	case 0x0158: //1242: // (g-,g+,Qb3-,Q3+)
	case 0x0176:
		return &A2QM2g6_eval<3,0,1,2>;
		break;

	case 0x7016: //1192: // Qb-g-Q+g+
	case 0x5018:
		return &A2QM2g7_eval<0,1,2,3>;
		break;
	case 0x6701: //
	case 0x8501:
		return &A2QM2g7_eval<1,2,3,0>;
		break;
	case 0x1670: //177: // (g+,Qb3+,Q3-,g-)
	case 0x1850:
		return &A2QM2g7_eval<2,3,0,1>;
		break;
	case 0x0167: //
	case 0x0185:
		return &A2QM2g7_eval<3,0,1,2>;
		break;

	case 0x7105: //886: //Qb-,g+,g-,Q-
		return &A2QM2g8m_eval<0,1,2,3>;
		break;
	case 0x5107:
		return &A2QM2g8p_eval<0,1,2,3>;
		break;
	case 0x5710: //136: //(Q3-,Qb3-,g+,g-)
		return &A2QM2g8m_eval<1,2,3,0>;
		break;
	case 0x7510:
		return &A2QM2g8p_eval<1,2,3,0>;
		break;
	case 0x0571: //816: //(g-,Qb3-,Q3-,g+)
		return &A2QM2g8m_eval<2,3,0,1>;
		break;
	case 0x0751:
		return &A2QM2g8p_eval<2,3,0,1>;
		break;
	case 0x1057: //1011: //(g+,g-,Q3-,Qb3-)
		return &A2QM2g8m_eval<3,0,1,2>;
		break;
	case 0x1075:
		return &A2QM2g8p_eval<3,0,1,2>;
		break;

	case 0x7015: //976: //Qb-,g-,g+,Q-
		return &A2QM2g9m_eval<0,1,2,3>;
		break;
	case 0x5017:
		return &A2QM2g9p_eval<0,1,2,3>;
		break;
	case 0x5701: //
		return &A2QM2g9m_eval<1,2,3,0>;
		break;
	case 0x7501:
		return &A2QM2g9p_eval<1,2,3,0>;
		break;
	case 0x1570: //171: //(g+,Q3-,Qb3-,g-)
		return &A2QM2g9m_eval<2,3,0,1>;
		break;
	case 0x1750:
		return &A2QM2g9p_eval<2,3,0,1>;
		break;
	case 0x0157: //1026: //(g-,g+,Qb3-,Q3-)
		return &A2QM2g9m_eval<3,0,1,2>;
		break;
	case 0x0175:
		return &A2QM2g9p_eval<3,0,1,2>;
		break;

	case 0x7150: //Qb-,g+,Q-,g-
		return &A2QM2g24m_eval<2,3,0,1>;
		break;
	case 0x0715: //g-,Qb-,g+,Q-
		return &A2QM2g24m_eval<3,0,1,2>;
		break;
	case 0x5071: //Q-,g-,Qb-,g+
		return &A2QM2g24m_eval<0,1,2,3>;
		break;
	case 0x1507: //g+,Q-,g-,Qb-
		return &A2QM2g24m_eval<1,2,3,0>;
		break;

	case 0x5170: //Q-,g+,Qb-,g-
		return &A2QM2g24p_eval<2,3,0,1>;
		break;
	case 0x0517: //g-,Q-,g+,Qb-
		return &A2QM2g24p_eval<3,0,1,2>;
		break;
	case 0x7051: //Qb-,g-,Q-,g+
		return &A2QM2g24p_eval<0,1,2,3>;
		break;
	case 0x1705: //g+,Qb-,g-,Q-
		return &A2QM2g24p_eval<1,2,3,0>;
		break;

	case 0x8005: //869: // Qb+g-g-Q-
	case 0x6007:
		return &A2QM2g10_eval<0,1,2,3>;
		break;
	case 0x5800: //34: // (Q3-,Qb3+,g-,g-)
	case 0x7600:
		return &A2QM2g10_eval<1,2,3,0>;
		break;
	case 0x0580: //204: // (g-,Qb3-,Q3+,g-)
	case 0x0760:
		return &A2QM2g10_eval<2,3,0,1>;
		break;
	case 0x0058: //1224: // (g-,g-,Qb3-,Q3+)
	case 0x0076:
		return &A2QM2g10_eval<3,0,1,2>;
		break;

	case 0x7006: //1084: // Qb-g-g-Q+
	case 0x5008:
		return &A2QM2g11_eval<0,1,2,3>;
		break;
	case 0x6700: //29: // (Q3+,Qb3-,g-,g-)
	case 0x8500:
		return &A2QM2g11_eval<1,2,3,0>;
		break;
	case 0x0670: //174: // (g-,Q3+,Qb3-,g-)
	case 0x0850:
		return &A2QM2g11_eval<2,3,0,1>;
		break;
	case 0x0067: //1044: //(g-,g-,Q3+,Qb3-)
	case 0x0085:
		return &A2QM2g11_eval<3,0,1,2>;
		break;

	case 0x5080: // (Q-,g-,Qb+,g-) NEW
	case 0x7060: // (Qb-,g-,Q+,g-) NEW
		return &A2QM2g17_eval<2,3,0,1>;
		break;
	case 0x0508: // (g-,Q-,g-,Qb+) NEW
	case 0x0706: // (g-,Qb-,g-,Q+) NEW
		return &A2QM2g17_eval<3,0,1,2>;
		break;
	case 0x8050: // (Qb+,g-,Q-,g-) NEW
	case 0x6070: // (Q+,g-,Qb-,g-) NEW
		return &A2QM2g17_eval<0,1,2,3>;
		break;
	case 0x0805: // (g-,Qb+,g-,Q-) NEW
	case 0x0607: // (g-,Q+,g-,Qb-) NEW
		return &A2QM2g17_eval<1,2,3,0>;
		break;

	case 0x8006: //1085: // Qb+g-g-Q+
		return &A2QM2g12m_eval<0,1,2,3>;
		break;
	case 0x6008:
		return &A2QM2g12p_eval<0,1,2,3>;
		break;
	case 0x6800: //35: // (Q3+,Qb3+,g-,g-)
		return &A2QM2g12m_eval<1,2,3,0>;
		break;
	case 0x8600:
		return &A2QM2g12p_eval<1,2,3,0>;
		break;
	case 0x0680: //210: // (g-,Q3+,Qb3+,g-)
		return &A2QM2g12m_eval<2,3,0,1>;
		break;
	case 0x0860:
		return &A2QM2g12p_eval<2,3,0,1>;
		break;
	case 0x0068: //1260: // (g-,g-,Q3+,Qb3+)
		return &A2QM2g12m_eval<3,0,1,2>;
		break;
	case 0x0086:
		return &A2QM2g12p_eval<3,0,1,2>;
		break;

	case 0x8060: // (Qb+,g-,Q+,g-) NEW
		return &A2QM2g18_eval<0,1,2,3>;
		break;
	case 0x0806: // (g-,Q+-,g-,Q+) NEW
		return &A2QM2g18_eval<1,2,3,0>;
		break;
	case 0x6080: // (Q+,g-,Qb+,g-) NEW
		return &A2QM2g18_eval<2,3,0,1>;
		break;
	case 0x0608: // (g-,Q+,g-,Qb+) NEW
		return &A2QM2g18_eval<3,0,1,2>;
		break;

	case 0x8106: //1103: //Qb+,g+,g-,Q+
		return &A2QM2g13m_eval<0,1,2,3>;
		break;
	case 0x6108:
		return &A2QM2g13p_eval<0,1,2,3>;
		break;
	case 0x6810: //143: //(Q3+,Qb3+,g+,g-)
		return &A2QM2g13m_eval<1,2,3,0>;
		break;
	case 0x8610:
		return &A2QM2g13p_eval<1,2,3,0>;
		break;
	case 0x0681: //858: //(g-,Qb3+,Q3+,g+)
		return &A2QM2g13m_eval<2,3,0,1>;
		break;
	case 0x0861:
		return &A2QM2g13p_eval<2,3,0,1>;
		break;
	case 0x1068: //1263: //(g+,g-,Q3+,Qb3+)
		return &A2QM2g13m_eval<3,0,1,2>;
		break;
	case 0x1086:
		return &A2QM2g13p_eval<3,0,1,2>;
		break;

	case 0x8160: //1103: //Qb+,g+,Q+,g-
		return &A2QM2g23p_eval<0,1,2,3>;
		break;
	case 0x0816: //
		return &A2QM2g23p_eval<1,2,3,0>;
		break;
	case 0x6081: //
		return &A2QM2g23p_eval<2,3,0,1>;
		break;
	case 0x1608: //
		return &A2QM2g23p_eval<3,0,1,2>;
		break;

	case 0x6180:
		return &A2QM2g23m_eval<0,1,2,3>;
		break;
	case 0x0618:
		return &A2QM2g23m_eval<1,2,3,0>;
		break;
	case 0x8061:
		return &A2QM2g23m_eval<2,3,0,1>;
		break;
	case 0x1806:
		return &A2QM2g23m_eval<3,0,1,2>;
		break;

	case 0x8016: //1193: //Qb+,g-,g+,Q+
		return &A2QM2g14m_eval<0,1,2,3>;
		break;
	case 0x6018:
		return &A2QM2g14p_eval<0,1,2,3>;
		break;
	case 0x6801: //
		return &A2QM2g14m_eval<1,2,3,0>;
		break;
	case 0x8601:
		return &A2QM2g14p_eval<1,2,3,0>;
		break;
	case 0x1680: //213: //(g+,Q3+,Qb3+,g-)
		return &A2QM2g14m_eval<2,3,0,1>;
		break;
	case 0x1860:
		return &A2QM2g14p_eval<2,3,0,1>;
		break;
	case 0x0168: //1278: //(g-,g+,Qb3+,Q3+)
		return &A2QM2g14m_eval<3,0,1,2>;
		break;
	case 0x0186:
		return &A2QM2g14p_eval<3,0,1,2>;
		break;

	case 0x8116: //1211: // Qb+g+g+Q+
	case 0x6811:
	case 0x1681:
	case 0x1168:
	case 0x6118: //
	case 0x8611:
	case 0x1861:
	case 0x1186:
	case 0x7005: //868: // Qb-g-g-Q-
	case 0x5700:
	case 0x0570:
	case 0x0057:
	case 0x5007: //
	case 0x7500:
	case 0x0750:
	case 0x0075:
		return &ZeroF;
		break;

	default:// We return zero for all other helicity combinations
		return 0;
	}
}


/*
 *
 *
 * The 2 massive quarks and two massless quark amplitudes
 *
 *
 */

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->L()*ep.p(i2)->L())*(ep.ref()->L()*ep.p(i2)->L())*ep.spb(i1,i2)/(complex<T>(0,-2)*(ep.ref()->L()*f1.L())*(ep.ref()->L()*fn.L())*ep.sp(i1,i2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->L()*ep.p(i2)->L())*(ep.ref()->L()*ep.p(i2)->L())*ep.spb(i1,i2)/(complex<T>(0,2)*(ep.ref()->L()*f1.L())*(ep.ref()->L()*fn.L())*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*fn.L())
								-mass*(ep.p(i2)->L()*ep.ref()->L())*(ep.ref()->Lt()*ep.p(i1)->Lt())/((ep.ref()->L()*f1.L())*(fn.Lt()*ep.ref()->Lt())))/(complex<T>(0,2)*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((fn.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*f1.L())
								-mass*(ep.p(i2)->L()*ep.ref()->L())*(ep.ref()->Lt()*ep.p(i1)->Lt())/((ep.ref()->L()*fn.L())*(f1.Lt()*ep.ref()->Lt())))/(complex<T>(0,-2)*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i2)->Lt())*(ep.ref()->Lt()*ep.p(i2)->Lt())/(complex<T>(0,-1)*(ep.ref()->Lt()*f1.Lt())*(ep.ref()->Lt()*fn.Lt())*ep.spb(i2,i1));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i2)->Lt())*(ep.ref()->Lt()*ep.p(i2)->Lt())/(complex<T>(0,1)*(ep.ref()->Lt()*f1.Lt())*(ep.ref()->Lt()*fn.Lt())*ep.spb(i2,i1));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((f1.L()*ep.p(i1)->L())*(ep.p(i2)->Lt()*fn.Lt())
								-mass*(ep.p(i2)->Lt()*ep.ref()->Lt())*(ep.ref()->L()*ep.p(i1)->L())/((ep.ref()->Lt()*f1.Lt())*(fn.L()*ep.ref()->L())))/(complex<T>(0,2)*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((fn.L()*ep.p(i1)->L())*(ep.p(i2)->Lt()*f1.Lt())
								-mass*(ep.p(i2)->Lt()*ep.ref()->Lt())*(ep.ref()->L()*ep.p(i1)->L())/((ep.ref()->Lt()*fn.Lt())*(f1.L()*ep.ref()->L())))/(complex<T>(0,-2)*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q7m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->L()*ep.p(i1)->L())*(ep.ref()->L()*ep.p(i1)->L())/(complex<T>(0,-1)*(ep.ref()->L()*f1.L())*(ep.ref()->L()*fn.L())*ep.spa(i2,i1));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q7p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->L()*ep.p(i1)->L())*(ep.ref()->L()*ep.p(i1)->L())/(complex<T>(0,1)*(ep.ref()->L()*f1.L())*(ep.ref()->L()*fn.L())*ep.spa(i2,i1));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.ref()->Lt()*ep.p(i1)->Lt())/(complex<T>(0,-1)*(ep.ref()->Lt()*f1.Lt())*(ep.ref()->Lt()*fn.Lt())*ep.spb(i2,i1));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.ref()->Lt()*ep.p(i1)->Lt())/(complex<T>(0,1)*(ep.ref()->Lt()*f1.Lt())*(ep.ref()->Lt()*fn.Lt())*ep.spb(i2,i1));
}



template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.p(i2)->Lt()*fn.Lt())*(ep.ref()->L()*ep.p(i1)->L())/(complex<T>(0,4)*(ep.ref()->L()*f1.L())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.p(i2)->Lt()*fn.Lt())*(ep.ref()->L()*ep.p(i1)->L())/(complex<T>(0,-4)*(ep.ref()->L()*f1.L())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_2p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return mass*(ep.p(i2)->Lt()*ep.ref()->Lt())*(ep.ref()->L()*ep.p(i1)->L())/(complex<T>(0,4)*(fn.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f1.L())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_2m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return mass*(ep.p(i2)->Lt()*ep.ref()->Lt())*(ep.ref()->L()*ep.p(i1)->L())/(complex<T>(0,-4)*(fn.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f1.L())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.p(i2)->Lt()*fn.Lt())*(f1.L()*ep.p(i1)->L())/(complex<T>(0,-4)*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,4)*(fn.L()*ep.ref()->L())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,-4)*(fn.L()*ep.ref()->L())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.p(i2)->L()*fn.L())*(f1.Lt()*ep.p(i1)->Lt())/(complex<T>(0,-4)*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_6p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return mass*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,4)*(ep.ref()->Lt()*f1.Lt())*(fn.L()*ep.ref()->L())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_6m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return mass*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,-4)*(ep.ref()->Lt()*f1.Lt())*(fn.L()*ep.ref()->L())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_7m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*fn.L())/(complex<T>(0,4)*(ep.ref()->Lt()*f1.Lt())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_7p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*fn.L())/(complex<T>(0,-4)*(ep.ref()->Lt()*f1.Lt())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i2)->Lt())*(f1.L()*ep.p(i1)->L())/(complex<T>(0,4)*(ep.ref()->Lt()*fn.Lt())*ep.sp(i1,i0));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i2)->Lt())*(f1.L()*ep.p(i1)->L())/(complex<T>(0,-4)*(ep.ref()->Lt()*fn.Lt())*ep.sp(i1,i0));
}


template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_9m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*ep.ref()->Lt())/(complex<T>(0,4)*(fn.Lt()*ep.ref()->Lt())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_9p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*ep.ref()->Lt())/(complex<T>(0,-4)*(fn.Lt()*ep.ref()->Lt())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_10m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*fn.Lt())/(complex<T>(0,4)*(ep.ref()->Lt()*f1.Lt())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_10p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*fn.Lt())/(complex<T>(0,-4)*(ep.ref()->Lt()*f1.Lt())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_11m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(f1.L()*ep.p(i1)->L())*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,-4)*(fn.L()*ep.ref()->L())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_11p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(f1.L()*ep.p(i1)->L())*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,4)*(fn.L()*ep.ref()->L())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_12m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->L()*ep.p(i1)->L())*(ep.p(i2)->L()*fn.L())/(complex<T>(0,-4)*(ep.ref()->L()*f1.L())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_12p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return eval_param<T>::mass(masses.p(i0))*(ep.ref()->L()*ep.p(i1)->L())*(ep.p(i2)->L()*fn.L())/(complex<T>(0,4)*(ep.ref()->L()*f1.L())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_13p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return mass*(ep.ref()->L()*ep.p(i1)->L())*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,-4)*(ep.ref()->L()*f1.L())*(fn.L()*ep.ref()->L())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_13m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return mass*(ep.ref()->L()*ep.p(i1)->L())*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,4)*(ep.ref()->L()*f1.L())*(fn.L()*ep.ref()->L())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_14p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return mass*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*ep.ref()->Lt())/(complex<T>(0,4)*(ep.ref()->Lt()*f1.Lt())*(fn.Lt()*ep.ref()->Lt())*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2f_14m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return mass*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*ep.ref()->Lt())/(complex<T>(0,-4)*(ep.ref()->Lt()*f1.Lt())*(fn.Lt()*ep.ref()->Lt())*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2fp_15_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (f1.L()*ep.p(i1)->L())*(ep.p(i2)->L()*fn.L())/(complex<T>(0,-4)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2fm_15_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (f1.L()*ep.p(i1)->L())*(ep.p(i2)->L()*fn.L())/(complex<T>(0,-4)*ep.sp(i2,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2fp_16_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*fn.Lt())/(complex<T>(0,-4)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2q2fm_16_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	// We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*fn.Lt())/(complex<T>(0,-4)*ep.sp(i2,i3));
}

template <class T> complex<T> (*A2QM2q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	// case 2q2Qs_massive
	switch(hc){
	case 0x8A96: //9731://(Q+,q3+,qb3-,Qb+)
	case 0xE32C: //12131://Qb2+,q+,qb-,Q2+
		return &A2QM2q1m_eval<0,1,2,3>;
		break;
	case 0x6A98:
	case 0xC32E:
		return &A2QM2q1p_eval<0,1,2,3>;
		break;
	case 0x68A9: //4561://(Qb+,Q+,q3+,qb3-)
	case 0xCE32: //1681://(Q2+,Qb2+,q+,qb-)
		return &A2QM2q1m_eval<1,2,3,0>;
		break;
	case 0x86A9:
	case 0xEC32:
		return &A2QM2q1p_eval<1,2,3,0>;
		break;
	case 0x968A: //
	case 0x2CE3: //
		return &A2QM2q1m_eval<2,3,0,1>;
		break;
	case 0x986A:
	case 0x2EC3:
		return &A2QM2q1p_eval<2,3,0,1>;
		break;
	case 0xA968: //10201://(q3+,qb3-,Q+,Qb+)
	case 0x32CE: //13081://(qb+,q-,Q2+,Qb2+)
		return &A2QM2q1m_eval<3,0,1,2>;
		break;
	case 0xA986:
	case 0x32EC:
		return &A2QM2q1p_eval<3,0,1,2>;
		break;

	case 0x8A95: //8400: //(Q+,q3+,qb3-,Qb-)
	case 0x6A97:
	case 0xE32B: //10800://Qb2+,q+,qb-,Q2-
	case 0xC32D:
		return &A2QM2q2_eval<0,1,2,3>;
		break;
	case 0x58A9: //4560://(Qb-,Q+,q3+,qb3-)
	case 0x76A9:
	case 0xBE32: //1680://(Q2-,Qb2+,q+,qb-)
	case 0xDC32:
		return &A2QM2q2_eval<1,2,3,0>;
		break;
	case 0x958A: //
	case 0x976A:
	case 0x2BE3: //
	case 0x2DC3:
		return &A2QM2q2_eval<2,3,0,1>;
		break;
	case 0xA958: //10080://(q3+,qb3-,Q-,Qb+)
	case 0xA976:
	case 0x32BE: //12960://qb+,q-,Q2-,Qb2+
	case 0x32DC:
		return &A2QM2q2_eval<3,0,1,2>;
		break;

	case 0x7A96: //9730://(Qb-,q3+,qb3-,Q+)
	case 0x5A98:
	case 0xD32C: //12130://Qb2-,q+,qb-,Q2+
	case 0xB32E:
		return &A2QM2q3_eval<0,1,2,3>;
		break;
	case 0x67A9: //4550://(Qb+,Q-,q3+,qb3-)
	case 0x85A9:
	case 0xCD32: //1670://(Q2+,Qb2-,q+,qb-)
	case 0xEB32:
		return &A2QM2q3_eval<1,2,3,0>;
		break;
	case 0x967A: //
	case 0x985A:
	case 0x2CD3: //
	case 0x2EB3:
		return &A2QM2q3_eval<2,3,0,1>;
		break;
	case 0xA967: //8870://(q3+,qb3-,Q+,Qb-)
	case 0xA985:
	case 0x32CD: //11750://(qb+,q-,Q2+,Qb2-)
	case 0x32EB:
		return &A2QM2q3_eval<3,0,1,2>;
		break;

	case 0x79A5: //8509://(Qb-,qb3-,q3+,Q-)
	case 0xD23B: //10909://Qb2-,q-,qb+,Q2-
		return &A2QM2q4m_eval<0,1,2,3>;
		break;
	case 0x59A7:
	case 0xB23D:
		return &A2QM2q4p_eval<0,1,2,3>;
		break;
	case 0x579A: //2879://(Qb3-,Q3-,qb-,q+)
	case 0xBD23: //5759://(Qb-,Q-,qb3-,q3+)
		return &A2QM2q4m_eval<1,2,3,0>;
		break;
	case 0x759A:
	case 0xDB23:
		return &A2QM2q4p_eval<1,2,3,0>;
		break;
	case 0xA579: //2389://(q+,Q3-,Qb3-,qb-)
	case 0x3BD2: //4789://(q3+,Q-,Qb-,qb3-)
		return &A2QM2q4m_eval<2,3,0,1>;
		break;
	case 0xA759:
	case 0x3DB2:
		return &A2QM2q4p_eval<2,3,0,1>;
		break;
	case 0x9A57: //8759://(qb3-,q3+,Q-,Qb-)
	case 0x23BD: //11639://(qb-,q+,Q2-,Qb2-)
		return &A2QM2q4m_eval<3,0,1,2>;
		break;
	case 0x9A75:
	case 0x23DB:
		return &A2QM2q4p_eval<3,0,1,2>;
		break;

	case 0x79A6: //9840://(Qb-,qb3-,q3+,Q+)
	case 0x59A8:
	case 0xD23C: //12240://Qb2-,q-,qb+,Q2+
	case 0xB23E:
		return &A2QM2q5_eval<0,1,2,3>;
		break;
	case 0x679A: //2880://(Qb3+,Q3-,qb-,q+)
	case 0x859A:
	case 0xCD23: //5760://(Qb+,Q-,qb3-,q3+)
	case 0xEB23:
		return &A2QM2q5_eval<1,2,3,0>;
		break;
	case 0xA679: //2400://(q+,Q3+,Qb3-,qb-)
	case 0xA859:
	case 0x3CD2: //4800://(q3+,Q+,Qb-,qb3-)
	case 0x3EB2:
		return &A2QM2q5_eval<2,3,0,1>;
		break;
	case 0x9A67: //8880://(qb3-,q3+,Q+,Qb-)
	case 0x9A85:
	case 0x23CD: //11760://(qb-,q+,Q2+,Qb2-)
	case 0x23EB:
		return &A2QM2q5_eval<3,0,1,2>;
		break;

	case 0x89A5: //8510://(Q+,qb3-,q3+,Qb-)
	case 0x69A7:
	case 0xE23B: //10910://Qb2+,q-,qb+,Q2-
	case 0xC23D:
		return &A2QM2q6_eval<0,1,2,3>;
		break;
	case 0x589A: //2890://(Qb3-,Q3+,qb-,q+)
	case 0x769A:
	case 0xBE23: //5770://(Qb-,Q+,qb3-,q3+)
	case 0xDC23:
		return &A2QM2q6_eval<1,2,3,0>;
		break;
	case 0xA589: //2510://(q+,Qb3-,Q3+,qb-)
	case 0xA769:
	case 0x3BE2: //4910://(q3+,Qb-,Q+,qb3-)
	case 0x3DC2:
		return &A2QM2q6_eval<2,3,0,1>;
		break;
	case 0x9A58: //10090://(qb3-,q3+,Q-,Qb+)
	case 0x9A76:
	case 0x23BE: //12970://(qb-,q+,Q2-,Qb2+)
	case 0x23DC:
		return &A2QM2q6_eval<3,0,1,2>;
		break;

	case 0x89A6: //9841://(Qb+,qb3-,q3+,Q+)
	case 0xE23C: //12241://Qb2+,qb-,q+,Q2+
		return &A2QM2q7m_eval<0,1,2,3>;
		break;
	case 0x69A8:
	case 0xC23E:
		return &A2QM2q7p_eval<0,1,2,3>;
		break;
	case 0x689A: //2891://(Qb3+,Q3+,qb-,q+)
	case 0xCE23: //5771://(Qb+,Q+,qb3-,q3+)
		return &A2QM2q7m_eval<1,2,3,0>;
		break;
	case 0x869A:
	case 0xEC23:
		return &A2QM2q7p_eval<1,2,3,0>;
		break;
	case 0xA689: //2521://(q+,Q3+,Qb3+,qb-)
	case 0x3CE2: //4921://(q3+,Q+,Qb+,qb3-)
		return &A2QM2q7m_eval<2,3,0,1>;
		break;
	case 0xA869:
	case 0x3EC2:
		return &A2QM2q7p_eval<2,3,0,1>;
		break;
	case 0x9A68: //10211://(qb3-,q3+,Q+,Qb+)
	case 0x23CE: //13091://(qb-,q+,Q2+,Qb2+)
		return &A2QM2q7m_eval<3,0,1,2>;
		break;
	case 0x9A86:
	case 0x23EC:
		return &A2QM2q7p_eval<3,0,1,2>;
		break;

	case 0x7A95: //8399://(Q-,q3+,qb3-,Qb-)
	case 0xD32B: //10799://Qb2-,qb+,q-,Q2-
		return &A2QM2q8m_eval<0,1,2,3>;
		break;
	case 0x5A97:
	case 0xB32D:
		return &A2QM2q8p_eval<0,1,2,3>;
		break;
	case 0x57A9: //4549://(Qb-,Q-,q3+,qb3-)
	case 0xBD32: //1669://(Q2-,Qb2-,q+,qb-)
		return &A2QM2q8m_eval<1,2,3,0>;
		break;
	case 0x75A9:
	case 0xDB32:
		return &A2QM2q8p_eval<1,2,3,0>;
		break;
	case 0x957A: //
	case 0x2BD3: //
		return &A2QM2q8m_eval<2,3,0,1>;
		break;
	case 0x975A:
	case 0x2DB3:
		return &A2QM2q8p_eval<2,3,0,1>;
		break;
	case 0xA957: //8749://(q3+,qb3-,Q-,Qb-)
	case 0x32BD: //11629://(qb+,q-,Q2-,Qb2-)
		return &A2QM2q8m_eval<3,0,1,2>;
		break;
	case 0xA975:
	case 0x32DB:
		return &A2QM2q8p_eval<3,0,1,2>;
		break;


	case 0x82AC: //(Qb+,q-,qb3+,Q3+)
	case 0xE936: //(Qb3+,q3-,qb+,Q+)
	case 0x82AE: //(Qb+,q-,q3+,Qb3+)
	case 0xE938: //(Qb+,q-,q3+,Qb3+)
		return &A2QM2q2f_1m_eval<0,1,2,3>;
		break;
	case 0x62AE: //(Q+,qb-,q3+,Qb3+)
	case 0xC938: //(Q+,qb-,q3+,Qb3+)
	case 0x62AC: //(Q+,qb-,qb3+,Q3+)
	case 0xC936: //(Q3+,qb3-,qb+,Q+)
		return &A2QM2q2f_1p_eval<0,1,2,3>;
		break;
	case 0xC82A: //12481://(Q+,qb-,q3+,Qb3+)
	case 0x6E93: //9601://(Qb3+,q3-,qb+,Q+)
	case 0xE82A:
	case 0x8E93:
		return &A2QM2q2f_1m_eval<1,2,3,0>;
		break;
	case 0xE62A:
	case 0x8C93:
	case 0xC62A:
	case 0x6C93:
		return &A2QM2q2f_1p_eval<1,2,3,0>;
		break;
	case 0xAC82: //12481://(Q+,qb-,q3+,Qb3+)
	case 0x36E9: //9601://(Qb3+,q3-,qb+,Q+)
	case 0xAE82:
	case 0x38E9:
		return &A2QM2q2f_1m_eval<2,3,0,1>;
		break;
	case 0xAE62:
	case 0x38C9:
	case 0xAC62:
	case 0x36C9:
		return &A2QM2q2f_1p_eval<2,3,0,1>;
		break;
	case 0x2AC8: //12481://(Q+,qb-,q3+,Qb3+)
	case 0x936E: //9601://(Qb3+,q3-,qb+,Q+)
	case 0x2AE8:
	case 0x938E:
		return &A2QM2q2f_1m_eval<3,0,1,2>;
		break;
	case 0x2AE6:
	case 0x938C:
	case 0x2AC6:
	case 0x936C:
		return &A2QM2q2f_1p_eval<3,0,1,2>;
		break;


	case 0x62AD: //(Q+,qb-,q3+,Qb3-)
	case 0xC937: //(Q3+,qb3-,q+,Qb-)
	case 0x82AB: //(Qb+,q-,qb3+,Q3-)
	case 0xE935: //(Qb3+,q3-,qb+,Q-)
		return &A2QM2q2f_2m_eval<0,1,2,3>;
		break;
	case 0x62AB: //(Q+,qb-,qb3+,Q3-)
	case 0xC935: //(Q3+,qb3-,qb+,Q-)
	case 0x82AD: //(Qb+,q-,q3+,Qb3-)
	case 0xE937: //(Qb3+,q3-,q+,Qb-)
		return &A2QM2q2f_2p_eval<0,1,2,3>;
		break;
	case 0xD62A: //
	case 0x7C93:
	case 0xB82A:
	case 0x5E93:
		return &A2QM2q2f_2m_eval<1,2,3,0>;
		break;
	case 0xB62A: //
	case 0x5C93:
	case 0xD82A:
	case 0x7E93:
		return &A2QM2q2f_2p_eval<1,2,3,0>;
		break;
	case 0xAD62: //
	case 0x37C9:
	case 0xAB82:
	case 0x35E9:
		return &A2QM2q2f_2m_eval<2,3,0,1>;
		break;
	case 0xAB62: //
	case 0x35C9:
	case 0xAD82:
	case 0x37E9:
		return &A2QM2q2f_2p_eval<2,3,0,1>;
		break;
	case 0x2AD6: //
	case 0x937C:
	case 0x2AB8:
	case 0x935E:
		return &A2QM2q2f_2m_eval<3,0,1,2>;
		break;
	case 0x2AB6: //
	case 0x935C:
	case 0x2AD8:
	case 0x937E:
		return &A2QM2q2f_2p_eval<3,0,1,2>;
		break;

	case 0x72AC: //(Qb-,q-,qb3+,Q3+)
	case 0x72AE: //(Qb-,q-,q3+,Qb3+)
	case 0x52AC: //(Q-,qb-,qb3+,Q3+)
	case 0x52AE: //(Q-,qb-,q3+,Qb3+)
	case 0xD936: //(Qb3-,q3-,qb+,Q+)
	case 0xD938: //(Qb3-,q3-,q+,Qb+)
	case 0xB936: //(Q3-,qb3-,qb+,Q+)
	case 0xB938: //(Q3-,qb3-,q+,Qb+)
		return &A2QM2q2f_3_eval<0,1,2,3>;
		break;
	case 0xC72A:
	case 0xE72A:
	case 0xC52A:
	case 0xE52A:
	case 0x6D93:
	case 0x8D93:
	case 0x6B93:
	case 0x8B93:
		return &A2QM2q2f_3_eval<1,2,3,0>;
		break;
	case 0xAC72:
	case 0xAE72:
	case 0xAC52:
	case 0xAE52:
	case 0x36D9:
	case 0x38D9:
	case 0x36B9:
	case 0x38B9:
		return &A2QM2q2f_3_eval<2,3,0,1>;
		break;
	case 0x2AC7:
	case 0x2AE7:
	case 0x2AC5:
	case 0x2AE5:
	case 0x936D:
	case 0x938D:
	case 0x936B:
	case 0x938B:
		return &A2QM2q2f_3_eval<3,0,1,2>;
		break;

	case 0x839C: //(Qb+,q+,qb3-,Q3+)
	case 0xEA26: //(Qb3+,q3+,qb-,Q+)
	case 0x639C: //(Q+,qb+,qb3-,Q3+)
	case 0xCA26: //(Q3+,qb3+,qb-,Q+)
		return &A2QM2q2f_4m_eval<0,1,2,3>;
		break;
	case 0x639E: //(Q+,qb+,q3-,Qb3+)
	case 0xCA28: //(Q3+,qb3+,q-,Qb+)
	case 0x839E: //(Qb+,q+,q3-,Qb3+)
	case 0xEA28: //(Qb3+,q3+,q-,Qb+)
		return &A2QM2q2f_4p_eval<0,1,2,3>;
		break;
	case 0xC839: //12371://(Qb+,q+,qb3-,Q3+)
	case 0x6EA2: //9491://(Qb3+,q3+,qb-,Q+)
	case 0xC639: //(Q+,qb+,qb3-,Q3+)
	case 0x6CA2:
		return &A2QM2q2f_4m_eval<1,2,3,0>;
		break;
	case 0xE639:
	case 0x8CA2:
	case 0xE839:
	case 0x8EA2:
		return &A2QM2q2f_4p_eval<1,2,3,0>;
		break;
	case 0x9C83: //12371://(Qb+,q+,qb3-,Q3+)
	case 0x26EA: //9491://(Qb3+,q3+,qb-,Q+)
	case 0x9C63: //(Q+,qb+,qb3-,Q3+)
	case 0x26CA:
		return &A2QM2q2f_4m_eval<2,3,0,1>;
		break;
	case 0x9E63:
	case 0x28CA:
	case 0x9E83:
	case 0x28EA:
		return &A2QM2q2f_4p_eval<2,3,0,1>;
		break;
	case 0x39C8: //12371://(Qb+,q+,qb3-,Q3+)
	case 0xA26E: //9491://(Qb3+,q3+,qb-,Q+)
	case 0x39C6: //(Q+,qb+,qb3-,Q3+)
	case 0xA26C:
		return &A2QM2q2f_4m_eval<3,0,1,2>;
		break;
	case 0x39E6:
	case 0xA28C:
	case 0x39E8:
	case 0xA28E:
		return &A2QM2q2f_4p_eval<3,0,1,2>;
		break;

	case 0x639B: //(Q+,qb+,q3-,Qb3-)
	case 0x839B: //(Qb+,q+,q3-,Qb3-)
	case 0xCA25: //(Q3+,qb3+,qb-,Q-)
	case 0xEA25: //(Qb3+,q3+,qb-,Q-)
	case 0x639D: //(Q+,qb+,qb3-,Q3-)
	case 0x839D: //(Qb+,q+,qb3-,Q3-)
	case 0xCA27: //(Q3+,qb3+,q-,Qb-)
	case 0xEA27: //(Qb3+,q3+,q-,Qb-)
		return &A2QM2q2f_5_eval<0,1,2,3>;
		break;
	case 0xB639:
	case 0xB839:
	case 0x5CA2:
	case 0x5EA2:
	case 0xD639:
	case 0xD839:
	case 0x7CA2:
	case 0x7EA2:
		return &A2QM2q2f_5_eval<1,2,3,0>;
		break;
	case 0x9B63:
	case 0x9B83:
	case 0x25CA:
	case 0x25EA:
	case 0x9D63:
	case 0x9D83:
	case 0x27CA:
	case 0x27EA:
		return &A2QM2q2f_5_eval<2,3,0,1>;
		break;
	case 0x39B6:
	case 0x39B8:
	case 0xA25C:
	case 0xA25E:
	case 0x39D6:
	case 0x39D8:
	case 0xA27C:
	case 0xA27E:
		return &A2QM2q2f_5_eval<3,0,1,2>;
		break;

	case 0x739C: //(Qb-,q+,qb3-,Q3+)
	case 0xDA26: //(Qb3-,q3+,qb-,Q+)
	case 0x539E: //(Q-,qb+,q3-,Qb3+)
	case 0xBA28: //(Q3-,qb3+,q-,Qb+)
		return &A2QM2q2f_6m_eval<0,1,2,3>;
		break;
	case 0x739E: //(Qb-,q+,q3-,Qb3+)
	case 0xDA28: //(Qb3-,q3+,q-,Qb+)
	case 0x539C: //(Q-,qb+,qb3-,Q3+)
	case 0xBA26: //(Q3-,qb3+,qb-,Q+)
		return &A2QM2q2f_6p_eval<0,1,2,3>;
		break;
	case 0xC739: //
	case 0x6DA2:
	case 0xE539:
	case 0x8BA2:
		return &A2QM2q2f_6m_eval<1,2,3,0>;
		break;
	case 0xE739:
	case 0x8DA2:
	case 0xC539:
	case 0x6BA2:
		return &A2QM2q2f_6p_eval<1,2,3,0>;
		break;
	case 0x9C73: //
	case 0x26DA:
	case 0x9E53:
	case 0x28BA:
		return &A2QM2q2f_6m_eval<2,3,0,1>;
		break;
	case 0x9E73:
	case 0x28DA:
	case 0x9C53:
	case 0x26BA:
		return &A2QM2q2f_6p_eval<2,3,0,1>;
		break;
	case 0x39C7: //
	case 0xA26D:
	case 0x39E5:
	case 0xA28B:
		return &A2QM2q2f_6m_eval<3,0,1,2>;
		break;
	case 0x39E7:
	case 0xA28D:
	case 0x39C5:
	case 0xA26B:
		return &A2QM2q2f_6p_eval<3,0,1,2>;
		break;

	case 0x739B: //(Qb-,q+,qb3-,Q3-)
	case 0xDA25: //(Qb3-,q3+,qb-,Q-)
	case 0x739D: //(Qb-,q+,q3-,Qb3-)
	case 0xDA27: //(Qb3-,q3+,q-,Qb-)
		return &A2QM2q2f_7m_eval<0,1,2,3>;
		break;
	case 0x539D:
	case 0xBA27:
	case 0x539B:
	case 0xBA25:
		return &A2QM2q2f_7p_eval<0,1,2,3>;
		break;
	case 0xB739: //(Qb-,q+,qb3-,Q3-)
	case 0x5DA2:
	case 0xD739:
	case 0x7DA2:
		return &A2QM2q2f_7m_eval<1,2,3,0>;
		break;
	case 0xD539:
	case 0x7BA2:
	case 0xB539:
	case 0x5BA2:
		return &A2QM2q2f_7p_eval<1,2,3,0>;
		break;
	case 0x9B73: //(Qb-,q+,qb3-,Q3-)
	case 0x25DA:
	case 0x9D73:
	case 0x27DA:
		return &A2QM2q2f_7m_eval<2,3,0,1>;
		break;
	case 0x9D53:
	case 0x27BA:
	case 0x9B53:
	case 0x25BA:
		return &A2QM2q2f_7p_eval<2,3,0,1>;
		break;
	case 0x39B7: //(Qb-,q+,qb3-,Q3-)
	case 0xA25D:
	case 0x39D7:
	case 0xA27D:
		return &A2QM2q2f_7m_eval<3,0,1,2>;
		break;
	case 0x39D5:
	case 0xA27B:
	case 0x39B5:
	case 0xA25B:
		return &A2QM2q2f_7p_eval<3,0,1,2>;
		break;

	case 0x72AB: //(Qb-,q-,qb3+,Q3-)
	case 0xD935: //(Q3-,qb3-,q+,Qb-)
	case 0x52AB: //(Q-,qb-,qb3+,Q3-)
	case 0xB935: //(Q3-,qb3-,q+,Qb-)
		return &A2QM2q2f_8m_eval<0,1,2,3>;
		break;
	case 0x52AD: //(Q-,qb-,q3+,Qb3-)
	case 0xB937: //(Q3-,qb3-,q+,Qb-)
	case 0x72AD: //(Qb-,q-,qb3+,Qb3-)
	case 0xD937: //(Qb3-,q3-,q+,Qb-)
		return &A2QM2q2f_8p_eval<0,1,2,3>;
		break;
	case 0xB72A: //(Qb-,q-,qb3+,Q3-)
	case 0x5D93:
	case 0xB52A:
	case 0x5B93:
		return &A2QM2q2f_8m_eval<1,2,3,0>;
		break;
	case 0xD52A:
	case 0x7B93:
	case 0xD72A:
	case 0x7D93:
		return &A2QM2q2f_8p_eval<1,2,3,0>;
		break;
	case 0xAB72: //(Qb-,q-,qb3+,Q3-)
	case 0x35D9:
	case 0xAB52:
	case 0x35B9:
		return &A2QM2q2f_8m_eval<2,3,0,1>;
		break;
	case 0xAD52:
	case 0x37B9:
	case 0xAD72:
	case 0x37D9:
		return &A2QM2q2f_8p_eval<2,3,0,1>;
		break;
	case 0x2AB7: //(Qb-,q-,qb3+,Q3-)
	case 0x935D:
	case 0x2AB5:
	case 0x935B:
		return &A2QM2q2f_8m_eval<3,0,1,2>;
		break;
	case 0x2AD5:
	case 0x937B:
	case 0x2AD7:
	case 0x937D:
		return &A2QM2q2f_8p_eval<3,0,1,2>;
		break;

		//NEWEST

	case 0x63AD: //(Q+,qb+,q2+,Qb2-)
	case 0xCA37: //(Q2+,qb2+,q+,Qb-)
	case 0x83AD: //(Qb+,q+,q2+,Qb2-)
	case 0xEA37: //(Qb2+,q2+,q+,Qb-)
		return &A2QM2q2f_9p_eval<0,1,2,3>;
		break;
	case 0xEA35: //(Qb2+,q2+,qb+,Q-)
	case 0x83AB: //(Qb+,q+,qb2+,Q2-)
	case 0xCA35: //(Q2+,qb2+,qb+,Q-)
	case 0x63AB: //(Q+,qb+,qb2+,Q2-)
		return &A2QM2q2f_9m_eval<0,1,2,3>;
		break;
	case 0xD63A: //
	case 0x7CA3:
	case 0xD83A:
	case 0x7EA3:
		return &A2QM2q2f_9p_eval<1,2,3,0>;
		break;
	case 0x5EA3:
	case 0xB83A:
	case 0x5CA3:
	case 0xB63A:
		return &A2QM2q2f_9m_eval<1,2,3,0>;
		break;
	case 0xAD63: //
	case 0x37CA:
	case 0xAD83:
	case 0x37EA:
		return &A2QM2q2f_9p_eval<2,3,0,1>;
		break;
	case 0x35EA:
	case 0xAB83:
	case 0x35CA:
	case 0xAB63:
		return &A2QM2q2f_9m_eval<2,3,0,1>;
		break;
	case 0x3AD6: //
	case 0xA37C:
	case 0x3AD8:
	case 0xA37E:
		return &A2QM2q2f_9p_eval<3,0,1,2>;
		break;
	case 0xA35E:
	case 0x3AB8:
	case 0xA35C:
	case 0x3AB6:
		return &A2QM2q2f_9m_eval<3,0,1,2>;
		break;

	case 0x73AC: //(Qb-,q+,qb2+,Q2+)
	case 0xDA36: //(Qb2-,q2+,qb+,Q+)
	case 0x73AE: //(Qb-,q+,q2+,Qb2+)
	case 0xDA38: //(Qb2-,q2+,q+,Qb+)
		return &A2QM2q2f_10m_eval<0,1,2,3>;
		break;
	case 0x53AE: //(Q-,qb+,q2+,Qb2+)
	case 0xBA38: //(Q2-,qb2+,q+,Qb+)
	case 0x53AC: //(Q-,qb+,qb2+,Q2+)
	case 0xBA36: //(Q2-,qb2+,qb+,Q+)
		return &A2QM2q2f_10p_eval<0,1,2,3>;
		break;
	case 0xC73A: //
	case 0x6DA3:
	case 0xE73A:
	case 0x8DA3:
		return &A2QM2q2f_10m_eval<1,2,3,0>;
		break;
	case 0xE53A:
	case 0x8BA3:
	case 0xC53A:
	case 0x6BA3:
		return &A2QM2q2f_10p_eval<1,2,3,0>;
		break;
	case 0xAC73: //
	case 0x36DA:
	case 0xAE73:
	case 0x38DA:
		return &A2QM2q2f_10m_eval<2,3,0,1>;
		break;
	case 0xAE53:
	case 0x38BA:
	case 0xAC53:
	case 0x36BA:
		return &A2QM2q2f_10p_eval<2,3,0,1>;
		break;
	case 0x3AC7: //
	case 0xA36D:
	case 0x3AE7:
	case 0xA38D:
		return &A2QM2q2f_10m_eval<3,0,1,2>;
		break;
	case 0x3AE5:
	case 0xA38B:
	case 0x3AC5:
	case 0xA36B:
		return &A2QM2q2f_10p_eval<3,0,1,2>;
		break;


	case 0x529E: //(Q-,qb-,q2-,Qb2+)
	case 0xB928: //(Q2-,qb2-,q-,Qb+)
	case 0x729E: //(Qb-,q-,q2-,Qb2+)
	case 0xD928: //(Qb2-,q-,q-,Qb+)
		return &A2QM2q2f_11m_eval<0,1,2,3>;
		break;
	case 0x529C: //(Q-,qb-,qb2-,Q+)
	case 0xB926: //(Q2-,qb2-,qb-,Q+)
	case 0x729C: //(Qb-,q-,qb2-,Q2+)
	case 0xD926: //(Qb2-,q2-,qb-,Q+)
		return &A2QM2q2f_11p_eval<0,1,2,3>;
		break;
	case 0xE529: //
	case 0x8B92:
	case 0xE729:
	case 0x8D92:
		return &A2QM2q2f_11m_eval<1,2,3,0>;
		break;
	case 0xC529:
	case 0x6B92:
	case 0xC729:
	case 0x6D92:
		return &A2QM2q2f_11p_eval<1,2,3,0>;
		break;
	case 0x9E52: //
	case 0x28B9:
	case 0x9E72:
	case 0x28D9:
		return &A2QM2q2f_11m_eval<2,3,0,1>;
		break;
	case 0x9C52:
	case 0x26B9:
	case 0x9C72:
	case 0x26D9:
		return &A2QM2q2f_11p_eval<2,3,0,1>;
		break;
	case 0x29E5://
	case 0x928B:
	case 0x29E7:
	case 0x928D:
		return &A2QM2q2f_11m_eval<3,0,1,2>;
		break;
	case 0x29C5:
	case 0x926B:
	case 0x29C7:
	case 0x926D:
		return &A2QM2q2f_11p_eval<3,0,1,2>;
		break;


	case 0x629D: //(Q+,qb-,q2-,Qb2-)
	case 0xC927: //(Q2+,qb2-,q-,Qb-)
	case 0x629B: //(Q+,qb-,qb2-,Q2-)
	case 0xC925: //(Q2+,qb2-,qb-,Q-)
		return &A2QM2q2f_12m_eval<0,1,2,3>;
		break;
	case 0x829B: //(Qb+,q-,qb2-,Q2-)
	case 0xE925: //(Qb2+,q2-,qb-,Q-)
	case 0x829D: //(Qb+,q-,q2-,Qb2-)
	case 0xE927: //(Qb2+,q2-,q-,Qb-)
		return &A2QM2q2f_12p_eval<0,1,2,3>;
		break;
	case 0xD629://
	case 0x7C92:
	case 0xB629:
	case 0x5C92:
		return &A2QM2q2f_12m_eval<1,2,3,0>;
		break;
	case 0xB829:
	case 0x5E92:
	case 0xD829:
	case 0x7E92:
		return &A2QM2q2f_12p_eval<1,2,3,0>;
		break;
	case 0x9D62://
	case 0x27C9:
	case 0x9B62:
	case 0x25C9:
		return &A2QM2q2f_12m_eval<2,3,0,1>;
		break;
	case 0x9B82:
	case 0x25E9:
	case 0x9D82:
	case 0x27E9:
		return &A2QM2q2f_12p_eval<2,3,0,1>;
		break;
	case 0x29D6://
	case 0x927C:
	case 0x29B6:
	case 0x925C:
		return &A2QM2q2f_12m_eval<3,0,1,2>;
		break;
	case 0x29B8:
	case 0x925E:
	case 0x29D8:
	case 0x927E:
		return &A2QM2q2f_12p_eval<3,0,1,2>;
		break;



	case 0x829C: //(Qb+,q-,qb2-,Q2+)
	case 0xE926: //(Qb2+,q2-,qb-,Q+)
	case 0x629E: //(Q+,qb-,q2-,Qb2+)
	case 0xC928: //(Q2+,qb2-,q-,Qb+)
		return &A2QM2q2f_13p_eval<0,1,2,3>;
		break;
	case 0x829E: //(Qb+,q-,q2-,Qb2+)
	case 0xE928: //(Qb2+,q2-,q-,Qb+)
	case 0x629C: //(Q+,qb-,qb-,Q2+)
	case 0xC926: //(Q2+,qb2-,qb-,Q+)
		return &A2QM2q2f_13m_eval<0,1,2,3>;
		break;
	case 0xC829: //
	case 0x6E92:
	case 0xE629:
	case 0x8C92:
		return &A2QM2q2f_13p_eval<1,2,3,0>;
		break;
	case 0xE829:
	case 0x8E92:
	case 0xC629:
	case 0x6C92:
		return &A2QM2q2f_13m_eval<1,2,3,0>;
		break;
	case 0x9C82: //
	case 0x26E9:
	case 0x9E62:
	case 0x28C9:
		return &A2QM2q2f_13p_eval<2,3,0,1>;
		break;
	case 0x9E82:
	case 0x28E9:
	case 0x9C62:
	case 0x26C9:
		return &A2QM2q2f_13m_eval<2,3,0,1>;
		break;
	case 0x29C8: //
	case 0x926E:
	case 0x29E6:
	case 0x928C:
		return &A2QM2q2f_13p_eval<3,0,1,2>;
		break;
	case 0x29E8:
	case 0x928E:
	case 0x29C6:
	case 0x926C:
		return &A2QM2q2f_13m_eval<3,0,1,2>;
		break;


	case 0x53AD: //(Q-,qb+,q+,Qb2-)
	case 0xBA37: //(Q2-,qb2+,q+,Qb-)
	case 0x73AB: //(Qb-,q+,qb2+,Q2-)
	case 0xDA35: //(Qb2-,q2+,qb+,Q-)
		return &A2QM2q2f_14m_eval<0,1,2,3>;
		break;
	case 0x73AD: //(Qb-,q+,q2+,Qb2-)
	case 0xDA37: //(Qb2-,q2+,q+,Qb-)
	case 0x53AB: //(Q-,qb+,qb2+,Q2-)
	case 0xBA35: //(Q2-,qb2+,qb+,Q-)
		return &A2QM2q2f_14p_eval<0,1,2,3>;
		break;
	case 0xD53A: //
	case 0x7BA3:
	case 0xB73A:
	case 0x5DA3:
		return &A2QM2q2f_14m_eval<1,2,3,0>;
		break;
	case 0xD73A:
	case 0x7DA3:
	case 0xB53A:
	case 0x5BA3:
		return &A2QM2q2f_14p_eval<1,2,3,0>;
		break;
	case 0xAD53: //
	case 0x37BA:
	case 0xAB73:
	case 0x35DA:
		return &A2QM2q2f_14m_eval<2,3,0,1>;
		break;
	case 0xAD73:
	case 0x37DA:
	case 0xAB53:
	case 0x35BA:
		return &A2QM2q2f_14p_eval<2,3,0,1>;
		break;
	case 0x3AD5: //
	case 0xA37B:
	case 0x3AB7:
	case 0xA35D:
		return &A2QM2q2f_14m_eval<3,0,1,2>;
		break;
	case 0x3AD7:
	case 0xA37D:
	case 0x3AB5:
	case 0xA35B:
		return &A2QM2q2f_14p_eval<3,0,1,2>;
		break;

	case 0x729B: //(Q-,qb-,q2-,Qb2-)
	case 0xD925: //(Qb2-,q2-,qb-,Q-)
	case 0x529D:
	case 0xB927:
		return &A2QM2q2fp_15_eval<0,1,2,3>;
		break;
	case 0x729D:
	case 0xD927:
	case 0x529B:
	case 0xB925:
		return &A2QM2q2fm_15_eval<0,1,2,3>;
		break;
	case 0xB729:
	case 0x5D92:
	case 0xD529:
	case 0x7B92:
		return &A2QM2q2fp_15_eval<1,2,3,0>;
		break;
	case 0xD729:
	case 0x7D92:
	case 0xB529:
	case 0x5B92:
		return &A2QM2q2fm_15_eval<1,2,3,0>;
		break;
	case 0x9B72:
	case 0x25D9:
	case 0x9D52:
	case 0x27B9:
		return &A2QM2q2fp_15_eval<2,3,0,1>;
		break;
	case 0x9D72:
	case 0x27D9:
	case 0x9B52:
	case 0x25B9:
		return &A2QM2q2fm_15_eval<2,3,0,1>;
		break;
	case 0x29B7:
	case 0x925D:
	case 0x29D5:
	case 0x927B:
		return &A2QM2q2fp_15_eval<3,0,1,2>;
		break;
	case 0x29D7:
	case 0x927D:
	case 0x29B5:
	case 0x925B:
		return &A2QM2q2fm_15_eval<3,0,1,2>;
		break;

	case 0x83AC: //(Qb+,q+,qb2+,Q2+)
	case 0xEA36: //(Qb2+,q2+,qb+,Q+)
	case 0x63AE: //(Q+,qb+,q2+,Qb2+)
	case 0xCA38: //(Q2+,qb2+,q+,Qb+)
		return &A2QM2q2fp_16_eval<0,1,2,3>;
		break;
	case 0x83AE: //(Qb+,q+,q2+,Qb2+)
	case 0xEA38: //(Qb2+,q2+,q+,Qb+)
	case 0x63AC: //(Q+,qb+,qb2+,Q2+)
	case 0xCA36: //(Q2+,qb2+,qb+,Q+)
		return &A2QM2q2fm_16_eval<0,1,2,3>;
		break;
	case 0xC83A:
	case 0x6EA3:
	case 0xE63A:
	case 0x8CA3:
		return &A2QM2q2fp_16_eval<1,2,3,0>;
		break;
	case 0xE83A:
	case 0x8EA3:
	case 0xC63A:
	case 0x6CA3:
		return &A2QM2q2fm_16_eval<1,2,3,0>;
		break;
	case 0xAC83:
	case 0x36EA:
	case 0xAE63:
	case 0x38CA:
		return &A2QM2q2fp_16_eval<2,3,0,1>;
		break;
	case 0xAE83:
	case 0x38EA:
	case 0xAC63:
	case 0x36CA:
		return &A2QM2q2fm_16_eval<2,3,0,1>;
		break;
	case 0x3AC8:
	case 0xA36E:
	case 0x3AE6:
	case 0xA38C:
		return &A2QM2q2fp_16_eval<3,0,1,2>;
		break;
	case 0x3AE8:
	case 0xA38E:
	case 0x3AC6:
	case 0xA36C:
		return &A2QM2q2fm_16_eval<3,0,1,2>;
		break;

	default:// We return zero for all other helicity combinations
		cout << "4 pt A2QM2q_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return &ZeroF;
	}
}


/*
 *
 *
 * The 1 massive quark, 1 quark and one scalar amplitudes
 *
 *
 */


template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(eval_param<T>::mass2(masses.p(i2))/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i2))*(ep.p(i1)->L()*qL)/((QL*qL)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(eval_param<T>::mass2(masses.p(i2))/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i2))*(ep.p(i1)->L()*qL)/((QL*qL)*sqrt(T(2)));
}

template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(eval_param<T>::mass2(masses.p(i2))/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());

	return complex<T>(0,1)*(ep.p(i1)->L()*QL)/sqrt(T(2));
}

template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(eval_param<T>::mass2(masses.p(i0))/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());

	return complex<T>(0,1)*(QL*ep.p(i1)->L())/sqrt(T(2));
}

template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(eval_param<T>::mass2(masses.p(i0))/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*(qL*ep.p(i1)->L())/((qL*QL)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(eval_param<T>::mass2(masses.p(i0))/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i0))*(qL*ep.p(i1)->L())/((qL*QL)*sqrt(T(2)));
}

template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs5m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(eval_param<T>::mass2(masses.p(i0))/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*(qLt*ep.p(i1)->Lt())/((qLt*QLt)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs5p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(eval_param<T>::mass2(masses.p(i0))/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i0))*(qLt*ep.p(i1)->Lt())/((qLt*QLt)*sqrt(T(2)));
}

template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(eval_param<T>::mass2(masses.p(i0))/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());

	return complex<T>(0,1)*(QLt*ep.p(i1)->Lt())/sqrt(T(2));
}

template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs7_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(eval_param<T>::mass2(masses.p(i2))/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());

	return complex<T>(0,1)*(ep.p(i1)->Lt()*QLt)/sqrt(T(2));
}

template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(eval_param<T>::mass2(masses.p(i2))/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i2))*(ep.p(i1)->Lt()*qLt)/((QLt*qLt)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1QM1qs8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(eval_param<T>::mass2(masses.p(i2))/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i2))*(ep.p(i1)->Lt()*qLt)/((QLt*qLt)*sqrt(T(2)));
}

template <class T> complex<T> (*A1QM1qs_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	// case 0x: //2qs_massive
	switch (hc) {
	case 0x428: //258: //(s,q3-,Qb+)
		return &A1QM1qs1m_eval<0,1,2>;
		break;
	case 0x426:
		return &A1QM1qs1p_eval<0,1,2>;
		break;
	case 0x284: //330: //(qb-,Q+,s)
		return &A1QM1qs1m_eval<2,0,1>;
		break;
	case 0x264:
		return &A1QM1qs1p_eval<2,0,1>;
		break;
	case 0x842: //96: //(Q+,s,qb-)
		return &A1QM1qs1m_eval<1,2,0>;
		break;
	case 0x642:
		return &A1QM1qs1p_eval<1,2,0>;
		break;

	case 0x427: //209: //(s,q3-,Qb-)
	case 0x425:
		return &A1QM1qs2_eval<0,1,2>;
		break;
	case 0x274: //323: //(qb-,Q-,s)
	case 0x254:
		return &A1QM1qs2_eval<2,0,1>;
		break;
	case 0x742: //95: //(qb-,Q-,s)
	case 0x542:
		return &A1QM1qs2_eval<1,2,0>;
		break;

	case 0x724: //305: //(Q-,q3-,s)
	case 0x524:
		return &A1QM1qs3_eval<0,1,2>;
		break;
	case 0x247: //239: //(qb-,s,Q-)
	case 0x245:
		return &A1QM1qs3_eval<2,0,1>;
		break;
	case 0x472: //83: //(s,Q3-,qb3-)
	case 0x452:
		return &A1QM1qs3_eval<1,2,0>;
		break;

	case 0x824: //306: //(Q+,q3-,s)
		return &A1QM1qs4m_eval<0,1,2>;
		break;
	case 0x624:
		return &A1QM1qs4p_eval<0,1,2>;
		break;
	case 0x248: //288: //(qb-,s,Q+)
		return &A1QM1qs4m_eval<2,0,1>;
		break;
	case 0x246:
		return &A1QM1qs4p_eval<2,0,1>;
		break;
	case 0x482: //90: //(s,Q3+,qb3-)
		return &A1QM1qs4m_eval<1,2,0>;
		break;
	case 0x462:
		return &A1QM1qs4p_eval<1,2,0>;
		break;

	case 0x734: //312: //(Q-,qb+,s)
		return &A1QM1qs5m_eval<0,1,2>;
		break;
	case 0x534:
		return &A1QM1qs5p_eval<0,1,2>;
		break;
	case 0x347: //240: //(q+,s,Qb-)
		return &A1QM1qs5m_eval<2,0,1>;
		break;
	case 0x345:
		return &A1QM1qs5p_eval<2,0,1>;
		break;
	case 0x473: //132: //(s,Qb-,q+)
		return &A1QM1qs5m_eval<1,2,0>;
		break;
	case 0x453:
		return &A1QM1qs5p_eval<1,2,0>;
		break;

	case 0x834: //313: //(Q+,qb+,s)
	case 0x634:
		return &A1QM1qs6_eval<0,1,2>;
		break;
	case 0x348: //289: //(q+,s,Qb+)
	case 0x346:
		return &A1QM1qs6_eval<2,0,1>;
		break;
	case 0x483: //139: //(s,Qb+,q+)
	case 0x463:
		return &A1QM1qs6_eval<1,2,0>;
		break;

	case 0x438: //265: //(s,qb+,Qb+)
	case 0x436:
		return &A1QM1qs7_eval<0,1,2>;
		break;
	case 0x384: //331: //(q+,Qb+,s)
	case 0x364:
		return &A1QM1qs7_eval<2,0,1>;
		break;
	case 0x843: //145: //(Qb+,s,q+)
	case 0x643:
		return &A1QM1qs7_eval<1,2,0>;
		break;

	case 0x437: //216: //(s,qb+,Qb-)
		return &A1QM1qs8m_eval<0,1,2>;
		break;
	case 0x435:
		return &A1QM1qs8p_eval<0,1,2>;
		break;
	case 0x374: //324: //(q+,Qb-,s)
		return &A1QM1qs8m_eval<2,0,1>;
		break;
	case 0x354:
		return &A1QM1qs8p_eval<2,0,1>;
		break;
	case 0x743: //144: //(Qb-,s,q+)
		return &A1QM1qs8m_eval<1,2,0>;
		break;
	case 0x543:
		return &A1QM1qs8p_eval<1,2,0>;
		break;

	default:// We return zero for all other helicity combinations
		cout << "3 pt A1QM1qs_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return &ZeroF;
	}
}



/*
 *
 *
 * The 1 massive quark, 1 quark, 1 gluon and 1 scalar amplitudes
 *
 *
 */


template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	Cmom<T> fQ(ep.mom(i3)-(eval_param<T>::mass2(masses.p(i3))/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spab(i2,i0,i1)*(ep.p(i2)->L()*fQ.L())/(complex<T>(0,2)*ep.spa(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (mass*ep.spb(i1,i2)*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,2)*ep.spa(i2,i1)*(fQ.L()*ep.ref()->L())*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs3m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i3))*ep.spab(i2,i0,i1)*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,-2)*(fQ.L()*ep.ref()->L())*ep.spa(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs3p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i3))*ep.spab(i2,i0,i1)*(ep.p(i2)->L()*ep.ref()->L())/(complex<T>(0,-2)*(fQ.L()*ep.ref()->L())*ep.spa(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i3))*(mass*(ep.p(i1)->Lt()*ep.ref()->Lt())/(fQ.Lt()*ep.ref()->Lt())
							+(ep.p(i1)->Lt()*ep.p(i0)->Sm()*fQ.L()))/(complex<T>(0,2)*ep.spa(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i3))*(mass*(ep.p(i1)->Lt()*ep.ref()->Lt())/(fQ.Lt()*ep.ref()->Lt())
							+(ep.p(i1)->Lt()*ep.p(i0)->Sm()*fQ.L()))/(complex<T>(0,2)*ep.spa(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (mass*ep.spa(i1,i2)*(ep.p(i2)->Lt()*ep.ref()->Lt())/(complex<T>(0,-2)*ep.spb(i2,i1)*(fQ.Lt()*ep.ref()->Lt())*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	Cmom<T> fQ(ep.mom(i3)-(eval_param<T>::mass2(masses.p(i3))/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spba(i2,i0,i1)*(ep.p(i2)->Lt()*fQ.Lt())/(complex<T>(0,-2)*ep.spb(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs7m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i3))*(mass*(ep.p(i1)->L()*ep.ref()->L())/(fQ.L()*ep.ref()->L())
					+(ep.p(i1)->L()*ep.p(i0)->Sm()*fQ.Lt()))/(complex<T>(0,-2)*ep.spb(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs7p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i3))*(mass*(ep.p(i1)->L()*ep.ref()->L())/(fQ.L()*ep.ref()->L())
					+(ep.p(i1)->L()*ep.p(i0)->Sm()*fQ.Lt()))/(complex<T>(0,-2)*ep.spb(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i3))*ep.spba(i2,i0,i1)*(ep.p(i2)->Lt()*ep.ref()->Lt())/(complex<T>(0,2)*(fQ.Lt()*ep.ref()->Lt())*ep.spb(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i3))*ep.spba(i2,i0,i1)*(ep.p(i2)->Lt()*ep.ref()->Lt())/(complex<T>(0,2)*(fQ.Lt()*ep.ref()->Lt())*ep.spb(i2,i1)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs9_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	Cmom<T> fQ(ep.mom(i0)-T(0.5)*(eval_param<T>::mass2(masses.p(i0))/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());

	return (ep.spab(i1,i3,i2)*(ep.p(i1)->L()*fQ.L())/(complex<T>(0,2)*ep.spa(i1,i2)*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs10_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (mass*ep.spb(i2,i1)*(ep.p(i1)->L()*ep.ref()->L())/(complex<T>(0,2)*ep.spa(i1,i2)*(fQ.L()*ep.ref()->L())*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs11m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i0))*ep.spab(i1,i3,i2)*(ep.p(i1)->L()*ep.ref()->L())/(complex<T>(0,2)*(fQ.L()*ep.ref()->L())*ep.spa(i1,i2)*ep.sp(i3,i2)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs11p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i0))*ep.spab(i1,i3,i2)*(ep.p(i1)->L()*ep.ref()->L())/(complex<T>(0,2)*(fQ.L()*ep.ref()->L())*ep.spa(i1,i2)*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs12m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i0))*(mass*(ep.p(i2)->Lt()*ep.ref()->Lt())/(fQ.Lt()*ep.ref()->Lt())
					+(ep.p(i2)->Lt()*ep.p(i3)->Sm()*fQ.L()))/(complex<T>(0,-2)*ep.spa(i1,i2)*ep.sp(i0,i1)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs12p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i0))*(mass*(ep.p(i2)->Lt()*ep.ref()->Lt())/(fQ.Lt()*ep.ref()->Lt())
					+(ep.p(i2)->Lt()*ep.p(i3)->Sm()*fQ.L()))/(complex<T>(0,-2)*ep.spa(i1,i2)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs13_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (mass*ep.spa(i2,i1)*(ep.p(i1)->Lt()*ep.ref()->Lt())/(complex<T>(0,-2)*ep.spb(i1,i2)*(fQ.Lt()*ep.ref()->Lt())*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs14_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	Cmom<T> fQ(ep.mom(i0)-T(0.5)*(eval_param<T>::mass2(masses.p(i0))/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());

	return (ep.spba(i1,i3,i2)*(ep.p(i1)->Lt()*fQ.Lt())/(complex<T>(0,-2)*ep.spb(i1,i2)*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs15m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i0))*(mass*(ep.p(i2)->L()*ep.ref()->L())/(fQ.L()*ep.ref()->L())
						+(ep.p(i2)->L()*ep.p(i3)->Sm()*fQ.Lt()))/(complex<T>(0,2)*ep.spb(i1,i2)*ep.sp(i0,i1)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs15p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i0))*(mass*(ep.p(i2)->L()*ep.ref()->L())/(fQ.L()*ep.ref()->L())
						+(ep.p(i2)->L()*ep.p(i3)->Sm()*fQ.Lt()))/(complex<T>(0,2)*ep.spb(i1,i2)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs16m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i0))*ep.spba(i1,i3,i2)*(ep.p(i1)->Lt()*ep.ref()->Lt())/(complex<T>(0,-2)*(fQ.Lt()*ep.ref()->Lt())*ep.spb(i1,i2)*ep.sp(i3,i2)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs16p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i0))*ep.spba(i1,i3,i2)*(ep.p(i1)->Lt()*ep.ref()->Lt())/(complex<T>(0,-2)*(fQ.Lt()*ep.ref()->Lt())*ep.spb(i1,i2)*ep.sp(i3,i2)))/sqrt(T(2));
}



template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs1S_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	Cmom<T> fQ(ep.mom(i3)-(eval_param<T>::mass2(masses.p(i3))/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spab(i1,i0,i2)*(ep.p(i1)->L()*fQ.L())
				/(complex<T>(0,2)*ep.spa(i1,i2)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs2Sm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i3))*ep.spab(i1,i0,i2)*(ep.p(i1)->L()*ep.ref()->L())
				/(complex<T>(0,2)*ep.spa(i1,i2)*(ep.ref()->L()*fQ.L())*ep.sp(i0,i1)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs2Sp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i3))*ep.spab(i1,i0,i2)*(ep.p(i1)->L()*ep.ref()->L())
				/(complex<T>(0,2)*ep.spa(i1,i2)*(ep.ref()->L()*fQ.L())*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs3S_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	// This amplitude was generated from the Feynman diagrams and includes both orderings of the position of the scalar
	return -(ep.spab(i1,i0,i2)*(ep.p(i1)->Lt()*ep.p(i0)->Sm()*ep.ref()->L())/(complex<T>(0,2)*ep.sp(i0,i1)*ep.spa(i1,i2)*(fQ.L()*ep.ref()->L()))
			+complex<T>(0,-1)*(ep.p(i2)->Lt()*ep.p(i0)->Sm()*ep.ref()->L())/(ep.spa(i2,i1)*(fQ.L()*ep.ref()->L())))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs4Sm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(((T(2)*ep.sp(i0,i1)+mass)*ep.spb(i2,i1)*(fQ.L()*ep.p(i1)->L())
									-mass*T(2)*ep.sp(i0,i1)*(ep.ref()->Lt()*ep.p(i2)->Lt())/(ep.ref()->Lt()*fQ.Lt()))
									/(complex<T>(0,2)*ep.sp(i0,i1)*ep.spa(i1,i2))
				+complex<T>(0,-1)*(ep.p(i2)->Lt()*ep.p(i0)->Sm()*fQ.L())/(ep.spa(i2,i1)))/eval_param<T>::mass(masses.p(i3))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs4Sp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (((T(2)*ep.sp(i0,i1)+mass)*ep.spb(i2,i1)*(fQ.L()*ep.p(i1)->L())
									-mass*T(2)*ep.sp(i0,i1)*(ep.ref()->Lt()*ep.p(i2)->Lt())/(ep.ref()->Lt()*fQ.Lt()))
									/(complex<T>(0,2)*ep.sp(i0,i1)*ep.spa(i1,i2))
				+complex<T>(0,-1)*(ep.p(i2)->Lt()*ep.p(i0)->Sm()*fQ.L())/(ep.spa(i2,i1)))/eval_param<T>::mass(masses.p(i3))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs5S_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spba(i1,i0,i2)*(ep.p(i1)->L()*ep.p(i0)->Sm()*ep.ref()->Lt())/(complex<T>(0,2)*ep.sp(i0,i1)*ep.spb(i1,i2)*(fQ.Lt()*ep.ref()->Lt()))
			+complex<T>(0,-1)*(ep.p(i2)->L()*ep.p(i0)->Sm()*ep.ref()->Lt())/(ep.spb(i2,i1)*(fQ.Lt()*ep.ref()->Lt())))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs6Sm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(((T(2)*ep.sp(i0,i1)+mass)*ep.spa(i2,i1)*(fQ.Lt()*ep.p(i1)->Lt())
			-mass*T(2)*ep.sp(i0,i1)*(ep.ref()->L()*ep.p(i2)->L())/(ep.ref()->L()*fQ.L()))
			/(complex<T>(0,-2)*ep.sp(i0,i1)*ep.spb(i1,i2))
			+complex<T>(0,1)*(ep.p(i2)->L()*ep.p(i0)->Sm()*fQ.Lt())/(ep.spb(i2,i1)))/eval_param<T>::mass(masses.p(i3))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs6Sp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (((T(2)*ep.sp(i0,i1)+mass)*ep.spa(i2,i1)*(fQ.Lt()*ep.p(i1)->Lt())
			-mass*T(2)*ep.sp(i0,i1)*(ep.ref()->L()*ep.p(i2)->L())/(ep.ref()->L()*fQ.L()))
			/(complex<T>(0,-2)*ep.sp(i0,i1)*ep.spb(i1,i2))
			+complex<T>(0,1)*(ep.p(i2)->L()*ep.p(i0)->Sm()*fQ.Lt())/(ep.spb(i2,i1)))/eval_param<T>::mass(masses.p(i3))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs7S_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spba(i1,i0,i2)*(ep.p(i1)->Lt()*fQ.Lt())
					/(complex<T>(0,-2)*ep.spb(i1,i2)*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs8Sm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i3))*ep.spba(i1,i0,i2)*(ep.p(i1)->Lt()*ep.ref()->Lt())
				/(complex<T>(0,-2)*ep.spb(i1,i2)*(ep.ref()->Lt()*fQ.Lt())*ep.sp(i0,i1)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs8Sp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i3))*ep.spba(i1,i0,i2)*(ep.p(i1)->Lt()*ep.ref()->Lt())
				/(complex<T>(0,-2)*ep.spb(i1,i2)*(ep.ref()->Lt()*fQ.Lt())*ep.sp(i0,i1)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs9Sm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (((T(2)*ep.sp(i3,i2)+mass)*ep.spb(i1,i2)*(fQ.L()*ep.p(i2)->L())
			-mass*T(2)*ep.sp(i3,i2)*(ep.ref()->Lt()*ep.p(i1)->Lt())/(ep.ref()->Lt()*fQ.Lt()))
			/(complex<T>(0,-2)*ep.sp(i3,i2)*ep.spa(i2,i1))
			+complex<T>(0,1)*(ep.p(i1)->Lt()*ep.p(i3)->Sm()*fQ.L())/(ep.spa(i1,i2)))/eval_param<T>::mass(masses.p(i0))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs9Sp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(((T(2)*ep.sp(i3,i2)+mass)*ep.spb(i1,i2)*(fQ.L()*ep.p(i2)->L())
			-mass*T(2)*ep.sp(i3,i2)*(ep.ref()->Lt()*ep.p(i1)->Lt())/(ep.ref()->Lt()*fQ.Lt()))
			/(complex<T>(0,-2)*ep.sp(i3,i2)*ep.spa(i2,i1))
			+complex<T>(0,1)*(ep.p(i1)->Lt()*ep.p(i3)->Sm()*fQ.L())/(ep.spa(i1,i2)))/eval_param<T>::mass(masses.p(i0))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs10S_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spab(i2,i3,i1)*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(complex<T>(0,-2)*ep.sp(i3,i2)*ep.spa(i2,i1)*(fQ.L()*ep.ref()->L()))
				+complex<T>(0,1)*(ep.p(i1)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(ep.spa(i1,i2)*(fQ.L()*ep.ref()->L())))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs11S_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spba(i2,i3,i1)*(ep.p(i2)->Lt()*fQ.Lt())
				/(complex<T>(0,-2)*ep.spb(i2,i1)*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs12Sm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i0))*ep.spba(i2,i3,i1)*(ep.p(i2)->Lt()*ep.ref()->Lt())
				/(complex<T>(0,2)*ep.spb(i2,i1)*(ep.ref()->Lt()*fQ.Lt())*ep.sp(i3,i2)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs12Sp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i0))*ep.spba(i2,i3,i1)*(ep.p(i2)->Lt()*ep.ref()->Lt())
				/(complex<T>(0,2)*ep.spb(i2,i1)*(ep.ref()->Lt()*fQ.Lt())*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs13S_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spab(i2,i3,i1)*(ep.p(i2)->L()*fQ.L())
				/(complex<T>(0,2)*ep.spa(i2,i1)*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs14Sm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (eval_param<T>::mass(masses.p(i0))*ep.spab(i2,i3,i1)*(ep.p(i2)->L()*ep.ref()->L())
			/(complex<T>(0,-2)*ep.spa(i2,i1)*(ep.ref()->L()*fQ.L())*ep.sp(i3,i2)))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs14Sp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(eval_param<T>::mass(masses.p(i0))*ep.spab(i2,i3,i1)*(ep.p(i2)->L()*ep.ref()->L())
			/(complex<T>(0,-2)*ep.spa(i2,i1)*(ep.ref()->L()*fQ.L())*ep.sp(i3,i2)))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs15Sm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (((T(2)*ep.sp(i3,i2)+mass)*ep.spa(i1,i2)*(fQ.Lt()*ep.p(i2)->Lt())
			-mass*T(2)*ep.sp(i3,i2)*(ep.ref()->L()*ep.p(i1)->L())/(ep.ref()->L()*fQ.L()))
			/(complex<T>(0,2)*ep.sp(i3,i2)*ep.spb(i2,i1))
			+complex<T>(0,-1)*(ep.p(i1)->L()*ep.p(i3)->Sm()*fQ.Lt())/(ep.spb(i1,i2)))/eval_param<T>::mass(masses.p(i0))/sqrt(T(2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs15Sp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return -(((T(2)*ep.sp(i3,i2)+mass)*ep.spa(i1,i2)*(fQ.Lt()*ep.p(i2)->Lt())
			-mass*T(2)*ep.sp(i3,i2)*(ep.ref()->L()*ep.p(i1)->L())/(ep.ref()->L()*fQ.L()))
			/(complex<T>(0,2)*ep.sp(i3,i2)*ep.spb(i2,i1))
			+complex<T>(0,-1)*(ep.p(i1)->L()*ep.p(i3)->Sm()*fQ.Lt())/(ep.spb(i1,i2)))/eval_param<T>::mass(masses.p(i0))/sqrt(T(2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A1QM1q1gs16S_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
	Cmom<T> fQ(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());

	return (ep.spba(i2,i3,i1)*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(complex<T>(0,2)*ep.sp(i3,i2)*ep.spb(i2,i1)*(fQ.Lt()*ep.ref()->Lt()))
				+complex<T>(0.,-1)*(ep.p(i1)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(ep.spb(i1,i2)*(fQ.Lt()*ep.ref()->Lt())))/sqrt(T(2));
}


template <class T> complex<T> (*A1QM1q1gs_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	// case 0x: //2qs_massive
	switch (hc) {
	case 0x4127: //1448://(s,g+,qb-,Q-)
	case 0x4125:
		return &A1QM1q1gs1_eval<0,1,2,3>;
		break;
	case 0x7412: //536://(Q-,s,g+,qb-)
	case 0x5412:
		return &A1QM1q1gs1_eval<1,2,3,0>;
		break;
	case 0x2741: //
	case 0x2541:
		return &A1QM1q1gs1_eval<2,3,0,1>;
		break;
	case 0x1274: //2264://(g+,qb-,Q-,s)
	case 0x1254:
		return &A1QM1q1gs1_eval<3,0,1,2>;
		break;

	case 0x4138: //1840://(s,g+,q+,Qb+)
	case 0x4136:
		return &A1QM1q1gs2_eval<0,1,2,3>;
		break;
	case 0x8413: //880://(Qb+,s,g+,q+)
	case 0x6413:
		return &A1QM1q1gs2_eval<1,2,3,0>;
		break;
	case 0x3841: //
	case 0x3641:
		return &A1QM1q1gs2_eval<2,3,0,1>;
		break;
	case 0x1384: //2320://(g+,q+,Qb+,s)
	case 0x1364:
		return &A1QM1q1gs2_eval<3,0,1,2>;
		break;

	case 0x4128: //1791://(s,g+,qb-,Q+)
		return &A1QM1q1gs3m_eval<0,1,2,3>;
		break;
	case 0x4126:
		return &A1QM1q1gs3p_eval<0,1,2,3>;
		break;
	case 0x8412: //537://(Q+,s,g+,qb-)
		return &A1QM1q1gs3m_eval<1,2,3,0>;
		break;
	case 0x6412:
		return &A1QM1q1gs3p_eval<1,2,3,0>;
		break;
	case 0x2841: //
		return &A1QM1q1gs3m_eval<2,3,0,1>;
		break;
	case 0x2641:
		return &A1QM1q1gs3p_eval<2,3,0,1>;
		break;
	case 0x1284: //2313://(g+,qb-,Q+,s)
		return &A1QM1q1gs3m_eval<3,0,1,2>;
		break;
	case 0x1264:
		return &A1QM1q1gs3p_eval<3,0,1,2>;
		break;

	case 0x4137: //1497://(s,g+,q+,Qb-)
		return &A1QM1q1gs4m_eval<0,1,2,3>;
		break;
	case 0x4135:
		return &A1QM1q1gs4p_eval<0,1,2,3>;
		break;
	case 0x7413: //879://(Qb-,s,g+,q+)
		return &A1QM1q1gs4m_eval<1,2,3,0>;
		break;
	case 0x5413:
		return &A1QM1q1gs4p_eval<1,2,3,0>;
		break;
	case 0x3741: //
		return &A1QM1q1gs4m_eval<2,3,0,1>;
		break;
	case 0x3541:
		return &A1QM1q1gs4p_eval<2,3,0,1>;
		break;
	case 0x1374: //2271://(g+,q+,Qb-,s)
		return &A1QM1q1gs4m_eval<3,0,1,2>;
		break;
	case 0x1354:
		return &A1QM1q1gs4p_eval<3,0,1,2>;
		break;

	case 0x4027: //1427://(s,g-,qb-,Q-)
	case 0x4025:
		return &A1QM1q1gs5_eval<0,1,2,3>;
		break;
	case 0x7402: //389://(Qb-,s,g-,q-)
	case 0x5402:
		return &A1QM1q1gs5_eval<1,2,3,0>;
		break;
	case 0x2740: //323://(qb-,Q-,s,g-)
	case 0x2540:
		return &A1QM1q1gs5_eval<2,3,0,1>;
		break;
	case 0x0274: //2261://(g-,q-,Qb-,s)
	case 0x0254:
		return &A1QM1q1gs5_eval<3,0,1,2>;
		break;

	case 0x4038: //1819://(s,g-,q+,Qb+)
	case 0x4036:
		return &A1QM1q1gs6_eval<0,1,2,3>;
		break;
	case 0x8403: //
	case 0x6403:
		return &A1QM1q1gs6_eval<1,2,3,0>;
		break;
	case 0x3840: //331://(q+,Qb+,s,g-)
	case 0x3640:
		return &A1QM1q1gs6_eval<2,3,0,1>;
		break;
	case 0x0384: //2317://(g-,q+,Qb+,s)
	case 0x0364:
		return &A1QM1q1gs6_eval<3,0,1,2>;
		break;

	case 0x4028: //1770://(s,g-,qb-,Q+)
		return &A1QM1q1gs7m_eval<0,1,2,3>;
		break;
	case 0x4026:
		return &A1QM1q1gs7p_eval<0,1,2,3>;
		break;
	case 0x8402: //390://(Qb+,s,g-,q-)
		return &A1QM1q1gs7m_eval<1,2,3,0>;
		break;
	case 0x6402:
		return &A1QM1q1gs7p_eval<1,2,3,0>;
		break;
	case 0x2840: //330://(qb-,Q+,s,g-)
		return &A1QM1q1gs7m_eval<2,3,0,1>;
		break;
	case 0x2640:
		return &A1QM1q1gs7p_eval<2,3,0,1>;
		break;
	case 0x0284: //2310://(g-,q-,Qb+,s)
		return &A1QM1q1gs7m_eval<3,0,1,2>;
		break;
	case 0x0264:
		return &A1QM1q1gs7p_eval<3,0,1,2>;
		break;

	case 0x4037: //1476://(s,g-,q+,Qb-)
		return &A1QM1q1gs8m_eval<0,1,2,3>;
		break;
	case 0x4035:
		return &A1QM1q1gs8p_eval<0,1,2,3>;
		break;
	case 0x7403: //
		return &A1QM1q1gs8m_eval<1,2,3,0>;
		break;
	case 0x5403:
		return &A1QM1q1gs8p_eval<1,2,3,0>;
		break;
	case 0x3740: //324://(q+,Qb-,s,g-)
		return &A1QM1q1gs8m_eval<2,3,0,1>;
		break;
	case 0x3540:
		return &A1QM1q1gs8p_eval<2,3,0,1>;
		break;
	case 0x0374: //2268://(g-,q+,Qb-,s)
		return &A1QM1q1gs8m_eval<3,0,1,2>;
		break;
	case 0x0354:
		return &A1QM1q1gs8p_eval<3,0,1,2>;
		break;

	case 0x7214: //2216://(Q-,qb-,g+,s)
	case 0x5214:
		return &A1QM1q1gs9_eval<0,1,2,3>;
		break;
	case 0x4721: //
	case 0x4521:
		return &A1QM1q1gs9_eval<1,2,3,0>;
		break;
	case 0x1472: //584://(g+,s,Qb-,q-)
	case 0x1452:
		return &A1QM1q1gs9_eval<2,3,0,1>;
		break;
	case 0x2147: //1688://(q-,g+,s,Qb-)
	case 0x2145:
		return &A1QM1q1gs9_eval<3,0,1,2>;
		break;

	case 0x8314: //2224://(Qb+,q+,g+,s)
	case 0x6314:
		return &A1QM1q1gs10_eval<0,1,2,3>;
		break;
	case 0x4831: //
	case 0x4631:
		return &A1QM1q1gs10_eval<1,2,3,0>;
		break;
	case 0x1483: //976://(g+,s,Qb+,q+)
	case 0x1463:
		return &A1QM1q1gs10_eval<2,3,0,1>;
		break;
	case 0x3148: //2032://(q+,g+,s,Qb+)
	case 0x3146:
		return &A1QM1q1gs10_eval<3,0,1,2>;
		break;

	case 0x8214: //2217://(Q+,qb-,g+,s)
		return &A1QM1q1gs11m_eval<0,1,2,3>;
		break;
	case 0x6214:
		return &A1QM1q1gs11p_eval<0,1,2,3>;
		break;
	case 0x4821: //
		return &A1QM1q1gs11m_eval<1,2,3,0>;
		break;
	case 0x4621:
		return &A1QM1q1gs11p_eval<1,2,3,0>;
		break;
	case 0x1482: //633://(g+,s,Qb+,q-)
		return &A1QM1q1gs11m_eval<2,3,0,1>;
		break;
	case 0x1462:
		return &A1QM1q1gs11p_eval<2,3,0,1>;
		break;
	case 0x2148: //2031://(q-,g+,s,Qb+)
		return &A1QM1q1gs11m_eval<3,0,1,2>;
		break;
	case 0x2146:
		return &A1QM1q1gs11p_eval<3,0,1,2>;
		break;

	case 0x7314: //2223://(Qb-,q+,g+,s)
		return &A1QM1q1gs12m_eval<0,1,2,3>;
		break;
	case 0x5314:
		return &A1QM1q1gs12p_eval<0,1,2,3>;
		break;
	case 0x4731: //
		return &A1QM1q1gs12m_eval<1,2,3,0>;
		break;
	case 0x4531:
		return &A1QM1q1gs12p_eval<1,2,3,0>;
		break;
	case 0x1473: //927://(g+,s,Qb-,q+)
		return &A1QM1q1gs12m_eval<2,3,0,1>;
		break;
	case 0x1453:
		return &A1QM1q1gs12p_eval<2,3,0,1>;
		break;
	case 0x3147: //1689://(q+,g+,s,Qb-)
		return &A1QM1q1gs12m_eval<3,0,1,2>;
		break;
	case 0x3145:
		return &A1QM1q1gs12p_eval<3,0,1,2>;
		break;

	case 0x7204: //2069://(Q-,qb-,g-,s)
	case 0x5204:
		return &A1QM1q1gs13_eval<0,1,2,3>;
		break;
	case 0x4720: //83://(s,Qb-,q-,g-)
	case 0x4520:
		return &A1QM1q1gs13_eval<1,2,3,0>;
		break;
	case 0x0472: //581://(g-,s,Qb-,q-)
	case 0x0452:
		return &A1QM1q1gs13_eval<2,3,0,1>;
		break;
	case 0x2047: //1667://(q-,g-,s,Qb-)
	case 0x2045:
		return &A1QM1q1gs13_eval<3,0,1,2>;
		break;

	case 0x8304: //2077://(Qb+,q+,g-,s)
	case 0x6304:
		return &A1QM1q1gs14_eval<0,1,2,3>;
		break;
	case 0x4830: //139://(s,Qb+,q+,g-)
	case 0x4630:
		return &A1QM1q1gs14_eval<1,2,3,0>;
		break;
	case 0x0483: //973://(g-,s,Qb+,q+)
	case 0x0463:
		return &A1QM1q1gs14_eval<2,3,0,1>;
		break;
	case 0x3048: //2011://(q+,g-,s,Qb+)
	case 0x3046:
		return &A1QM1q1gs14_eval<3,0,1,2>;
		break;

	case 0x8204: //2070://(Q+,qb-,g-,s)
		return &A1QM1q1gs15m_eval<0,1,2,3>;
		break;
	case 0x6204:
		return &A1QM1q1gs15p_eval<0,1,2,3>;
		break;
	case 0x4820: //90://(s,Qb+,q-,g-)
		return &A1QM1q1gs15m_eval<1,2,3,0>;
		break;
	case 0x4620:
		return &A1QM1q1gs15p_eval<1,2,3,0>;
		break;
	case 0x0482: //630://(g-,s,Qb+,q-)
		return &A1QM1q1gs15m_eval<2,3,0,1>;
		break;
	case 0x0462:
		return &A1QM1q1gs15p_eval<2,3,0,1>;
		break;
	case 0x2048: //2010://(q-,g-,s,Qb+)
		return &A1QM1q1gs15m_eval<3,0,1,2>;
		break;
	case 0x2046:
		return &A1QM1q1gs15p_eval<3,0,1,2>;
		break;

	case 0x7304: //2076://(Qb-,q+,g-,s)
		return &A1QM1q1gs16m_eval<0,1,2,3>;
		break;
	case 0x5304:
		return &A1QM1q1gs16p_eval<0,1,2,3>;
		break;
	case 0x4730: //132://(s,Qb-,q+,g-)
		return &A1QM1q1gs16m_eval<1,2,3,0>;
		break;
	case 0x4530:
		return &A1QM1q1gs16p_eval<1,2,3,0>;
		break;
	case 0x0473: //924://(g-,s,Qb-,q+)
		return &A1QM1q1gs16m_eval<2,3,0,1>;
		break;
	case 0x0453:
		return &A1QM1q1gs16p_eval<2,3,0,1>;
		break;
	case 0x3047: //1668://(q+,g-,s,Qb-)
		return &A1QM1q1gs16m_eval<3,0,1,2>;
		break;
	case 0x3045:
		return &A1QM1q1gs16p_eval<3,0,1,2>;
		break;

	case 0x4217: //1532://(s,qb-,g+,Q-)
	case 0x4215:
		return &A1QM1q1gs1S_eval<0,1,2,3>;
		break;
	case 0x7421: //
	case 0x5421:
		return &A1QM1q1gs1S_eval<1,2,3,0>;
		break;
	case 0x1742: //668://(g+,Qb-,s,q-)
	case 0x1542:
		return &A1QM1q1gs1S_eval<2,3,0,1>;
		break;
	case 0x2174: //2276://(q-,g+,Qb-,s)
	case 0x2154:
		return &A1QM1q1gs1S_eval<3,0,1,2>;
		break;

	case 0x4218: //1875://(s,qb-,g+,Q+)
		return &A1QM1q1gs2Sm_eval<0,1,2,3>;
		break;
	case 0x4216:
		return &A1QM1q1gs2Sp_eval<0,1,2,3>;
		break;
	case 0x8421: //
		return &A1QM1q1gs2Sm_eval<1,2,3,0>;
		break;
	case 0x6421:
		return &A1QM1q1gs2Sp_eval<1,2,3,0>;
		break;
	case 0x1842: //675://(g+,Qb+,s,q-)
		return &A1QM1q1gs2Sm_eval<2,3,0,1>;
		break;
	case 0x1642:
		return &A1QM1q1gs2Sp_eval<2,3,0,1>;
		break;
	case 0x2184: //2325://(q-,g+,Qb+,s)
		return &A1QM1q1gs2Sm_eval<3,0,1,2>;
		break;
	case 0x2164:
		return &A1QM1q1gs2Sp_eval<3,0,1,2>;
		break;

	case 0x4318: //1882://(s,q+,g+,Qb+)
	case 0x4316:
		return &A1QM1q1gs3S_eval<0,1,2,3>;
		break;
	case 0x8431: //
	case 0x6431:
		return &A1QM1q1gs3S_eval<1,2,3,0>;
		break;
	case 0x1843: //1018://(g+,Qb+,s,q+)
	case 0x1643:
		return &A1QM1q1gs3S_eval<2,3,0,1>;
		break;
	case 0x3184: //2326://(q+,g+,Qb+,s)
	case 0x3164:
		return &A1QM1q1gs3S_eval<3,0,1,2>;
		break;

	case 0x4317: //1539://(s,q+,g+,Qb-)
		return &A1QM1q1gs4Sm_eval<0,1,2,3>;
		break;
	case 0x4315:
		return &A1QM1q1gs4Sp_eval<0,1,2,3>;
		break;
	case 0x7431: //
		return &A1QM1q1gs4Sm_eval<1,2,3,0>;
		break;
	case 0x5431:
		return &A1QM1q1gs4Sp_eval<1,2,3,0>;
		break;
	case 0x1743: //1011://(g+,Qb-,s,q+)
		return &A1QM1q1gs4Sm_eval<2,3,0,1>;
		break;
	case 0x1543:
		return &A1QM1q1gs4Sp_eval<2,3,0,1>;
		break;
	case 0x3174: //2277://(q+,g+,Qb-,s)
		return &A1QM1q1gs4Sm_eval<3,0,1,2>;
		break;
	case 0x3154:
		return &A1QM1q1gs4Sp_eval<3,0,1,2>;
		break;

	case 0x4207: //1385://(s,qb-,g-,Q-)
	case 0x4205:
		return &A1QM1q1gs5S_eval<0,1,2,3>;
		break;
	case 0x7420: //95://(Qb-,s,q-,g-)
	case 0x5420:
		return &A1QM1q1gs5S_eval<1,2,3,0>;
		break;
	case 0x0742: //665://(g-,Qb-,s,q-)
	case 0x0542:
		return &A1QM1q1gs5S_eval<2,3,0,1>;
		break;
	case 0x2074: //2255://(q-,g-,Qb-,s)
	case 0x2054:
		return &A1QM1q1gs5S_eval<3,0,1,2>;
		break;

	case 0x4208: //1728://(s,qb-,g-,Q+)
		return &A1QM1q1gs6Sm_eval<0,1,2,3>;
		break;
	case 0x4206:
		return &A1QM1q1gs6Sp_eval<0,1,2,3>;
		break;
	case 0x8420: //96://(Qb+,s,q-,g-)
		return &A1QM1q1gs6Sm_eval<1,2,3,0>;
		break;
	case 0x6420:
		return &A1QM1q1gs6Sp_eval<1,2,3,0>;
		break;
	case 0x0842: //672://(g-,Qb+,s,q-)
		return &A1QM1q1gs6Sm_eval<2,3,0,1>;
		break;
	case 0x0642:
		return &A1QM1q1gs6Sp_eval<2,3,0,1>;
		break;
	case 0x2084: //2304://(q-,g-,Qb+,s)
		return &A1QM1q1gs6Sm_eval<3,0,1,2>;
		break;
	case 0x2064:
		return &A1QM1q1gs6Sp_eval<3,0,1,2>;
		break;

	case 0x4308: //1735://(s,q+,g-,Qb+)
	case 0x4306:
		return &A1QM1q1gs7S_eval<0,1,2,3>;
		break;
	case 0x8430: //145://(Qb+,s,q+,g-)
	case 0x6430:
		return &A1QM1q1gs7S_eval<1,2,3,0>;
		break;
	case 0x0843: //1015://(g-,Qb+,s,q+)
	case 0x0643:
		return &A1QM1q1gs7S_eval<2,3,0,1>;
		break;
	case 0x3084: //2305://(q+,g-,Qb+,s)
	case 0x3064:
		return &A1QM1q1gs7S_eval<3,0,1,2>;
		break;

	case 0x4307: //1392://(s,q+,g-,Qb-)
		return &A1QM1q1gs8Sm_eval<0,1,2,3>;
		break;
	case 0x4305:
		return &A1QM1q1gs8Sp_eval<0,1,2,3>;
		break;
	case 0x7430: //144://(Qb-,s,q+,g-)
		return &A1QM1q1gs8Sm_eval<1,2,3,0>;
		break;
	case 0x5430:
		return &A1QM1q1gs8Sp_eval<1,2,3,0>;
		break;
	case 0x0743: //1008://(g-,Qb-,s,q+)
		return &A1QM1q1gs8Sm_eval<2,3,0,1>;
		break;
	case 0x0543:
		return &A1QM1q1gs8Sp_eval<2,3,0,1>;
		break;
	case 0x3074: //2256://(q+,g-,Qb-,s)
		return &A1QM1q1gs8Sm_eval<3,0,1,2>;
		break;
	case 0x3054:
		return &A1QM1q1gs8Sp_eval<3,0,1,2>;
		break;

	case 0x7134: //2181://(Qb-,g+,q+,s)
		return &A1QM1q1gs9Sm_eval<0,1,2,3>;
		break;
	case 0x5134:
		return &A1QM1q1gs9Sp_eval<0,1,2,3>;
		break;
	case 0x4713: //867://(s,Qb-,g+,q+)
		return &A1QM1q1gs9Sm_eval<1,2,3,0>;
		break;
	case 0x4513:
		return &A1QM1q1gs9Sp_eval<1,2,3,0>;
		break;
	case 0x3471: //
		return &A1QM1q1gs9Sm_eval<2,3,0,1>;
		break;
	case 0x3451:
		return &A1QM1q1gs9Sp_eval<2,3,0,1>;
		break;
	case 0x1347: //1683://(g+,q+,s,Qb-)
		return &A1QM1q1gs9Sm_eval<3,0,1,2>;
		break;
	case 0x1345:
		return &A1QM1q1gs9Sp_eval<3,0,1,2>;
		break;

	case 0x8134: //2182://(Qb+,g+,q+,s)
	case 0x6134:
		return &A1QM1q1gs10S_eval<0,1,2,3>;
		break;
	case 0x4813: //874://(s,Qb+,g+,q+)
	case 0x4613:
		return &A1QM1q1gs10S_eval<1,2,3,0>;
		break;
	case 0x3481: //
	case 0x3461:
		return &A1QM1q1gs10S_eval<2,3,0,1>;
		break;
	case 0x1348: //2026://(g+,q+,s,Qb+)
	case 0x1346:
		return &A1QM1q1gs10S_eval<3,0,1,2>;
		break;

	case 0x8034: //2161://(Qb+,g-,q+,s)
	case 0x6034:
		return &A1QM1q1gs11S_eval<0,1,2,3>;
		break;
	case 0x4803: //
	case 0x4603:
		return &A1QM1q1gs11S_eval<1,2,3,0>;
		break;
	case 0x3480: //289://(q+,s,Qb+,g-)
	case 0x3460:
		return &A1QM1q1gs11S_eval<2,3,0,1>;
		break;
	case 0x0348: //2023://(g-,q+,s,Qb+)
	case 0x0346:
		return &A1QM1q1gs11S_eval<3,0,1,2>;
		break;

	case 0x7034: //2160://(Qb-,g-,q+,s)
		return &A1QM1q1gs12Sm_eval<0,1,2,3>;
		break;
	case 0x5034:
		return &A1QM1q1gs12Sp_eval<0,1,2,3>;
		break;
	case 0x4703: //
		return &A1QM1q1gs12Sm_eval<1,2,3,0>;
		break;
	case 0x4503:
		return &A1QM1q1gs12Sp_eval<1,2,3,0>;
		break;
	case 0x3470: //240://(q+,s,Qb-,g-)
		return &A1QM1q1gs12Sm_eval<2,3,0,1>;
		break;
	case 0x3450:
		return &A1QM1q1gs12Sp_eval<2,3,0,1>;
		break;
	case 0x0347: //1680://(g-,q+,s,Qb-)
		return &A1QM1q1gs12Sm_eval<3,0,1,2>;
		break;
	case 0x0345:
		return &A1QM1q1gs12Sp_eval<3,0,1,2>;
		break;

	case 0x7124: //2132://(Q-,g+,qb-,s)
	case 0x5124:
		return &A1QM1q1gs13S_eval<0,1,2,3>;
		break;
	case 0x4712: //524://(s,Q-,g+,qb-)
	case 0x4512:
		return &A1QM1q1gs13S_eval<1,2,3,0>;
		break;
	case 0x2471: //
	case 0x2451:
		return &A1QM1q1gs13S_eval<2,3,0,1>;
		break;
	case 0x1247: //1676://(g+,qb-,s,Q-)
	case 0x1245:
		return &A1QM1q1gs13S_eval<3,0,1,2>;
		break;

	case 0x8124: //2133://(Q+,g+,qb-,s)
		return &A1QM1q1gs14Sm_eval<0,1,2,3>;
		break;
	case 0x6124:
		return &A1QM1q1gs14Sp_eval<0,1,2,3>;
		break;
	case 0x4812: //531://(s,Q+,g+,qb-)
		return &A1QM1q1gs14Sm_eval<1,2,3,0>;
		break;
	case 0x4612:
		return &A1QM1q1gs14Sp_eval<1,2,3,0>;
		break;
	case 0x2481: //
		return &A1QM1q1gs14Sm_eval<2,3,0,1>;
		break;
	case 0x2461:
		return &A1QM1q1gs14Sp_eval<2,3,0,1>;
		break;
	case 0x1248: //2019://(g+,qb-,s,Q+)
		return &A1QM1q1gs14Sm_eval<3,0,1,2>;
		break;
	case 0x1246:
		return &A1QM1q1gs14Sp_eval<3,0,1,2>;
		break;

	case 0x8024: //2112://(Q+,g-,qb-,s)
		return &A1QM1q1gs15Sm_eval<0,1,2,3>;
		break;
	case 0x6024:
		return &A1QM1q1gs15Sp_eval<0,1,2,3>;
		break;
	case 0x4802: //384://(s,Qb+,g-,q-)
		return &A1QM1q1gs15Sm_eval<1,2,3,0>;
		break;
	case 0x4602:
		return &A1QM1q1gs15Sp_eval<1,2,3,0>;
		break;
	case 0x2480: //288://(qb-,s,Q+,g-)
		return &A1QM1q1gs15Sm_eval<2,3,0,1>;
		break;
	case 0x2460:
		return &A1QM1q1gs15Sp_eval<2,3,0,1>;
		break;
	case 0x0248: //2016://(g-,q-,s,Qb+)
		return &A1QM1q1gs15Sm_eval<3,0,1,2>;
		break;
	case 0x0246:
		return &A1QM1q1gs15Sp_eval<3,0,1,2>;
		break;

	case 0x7024: //2111://(Q-,g-,qb-,s)
	case 0x5024:
		return &A1QM1q1gs16S_eval<0,1,2,3>;
		break;
	case 0x4702: //377://(s,Qb-,g-,q-)
	case 0x4502:
		return &A1QM1q1gs16S_eval<1,2,3,0>;
		break;
	case 0x2470: //239://(qb-,s,Q-,g-)
	case 0x2450:
		return &A1QM1q1gs16S_eval<2,3,0,1>;
		break;
	case 0x0247: //1673://(g-,q-,s,Qb-)
	case 0x0245:
		return &A1QM1q1gs16S_eval<3,0,1,2>;
		break;

	default:// We return zero for all other helicity combinations
		return 0;
	}
}



/*
 *
 *
 *
 * Two scalar two quark amplitudes
 *
 *
 *
 */

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2s2q1_eval(const eval_param<T>& ep, const mass_param_coll& masses){

	// The first term here corresponds to the S channel gluon exchange bewteen a fermion and a scalar line
	//  the second line is a T channel exchange of a massive quark between scalar massless quark vertices
	//	return complex<T>(0.,1)*ep.spab(i1,i0,i2)/(ep.sp(i1,i2))
	//					+complex<T>(0.,0.5)*ep.spab(i1,i0,i2)/ep.sp(i0,i1);

	return ep.spab(i1,i0,i2)*(T(1)/ep.sp(i1,i2)+T(1)/(T(2)*ep.sp(i0,i1)))/complex<T>(0,-2);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2s2q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return ep.spab(i2,i0,i1)*(T(1)/ep.sp(i1,i2)+T(0.5)/ep.sp(i0,i1))/complex<T>(0,2);
}

template <class T> complex<T> (*A2s2q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	// case 0x: //2qs_massive
	switch (hc) {
	case 0x4234: //2169: //(s,qb-,q+,s)
	case 0x49A4:
		return &A2s2q1_eval<0,1,2,3>;
		break;
	case 0x4423: //783: //(s,s,qb-,q+)
	case 0x449A:
		return &A2s2q1_eval<1,2,3,0>;
		break;
	case 0x3442: //681: //(q+,s,s,qb-)
	case 0xA449:
		return &A2s2q1_eval<2,3,0,1>;
		break;
	case 0x2344: //2367: //(qb-,q+,s,s)
	case 0x9A44:
		return &A2s2q1_eval<3,0,1,2>;
		break;

	case 0x4324: //2127: //(s,q+,qb-,s)
	case 0x4A94:
		return &A2s2q2_eval<0,1,2,3>;
		break;
	case 0x4432: //489: //(s,s,q+,qb-)
	case 0x44A9:
		return &A2s2q2_eval<1,2,3,0>;
		break;
	case 0x2443: //1023: //(qb-,s,s,q+)
	case 0x944A:
		return &A2s2q2_eval<2,3,0,1>;
		break;
	case 0x3244: //2361: //(q+,qb-,s,s)
	case 0xA944:
		return &A2s2q2_eval<3,0,1,2>;
		break;

	case 0x42A4: //2169: //(s,qb-,q+,s)
	case 0x4934:
	case 0x442A: //783: //(s,s,qb-,q+)
	case 0x4493:
	case 0xA442: //681: //(q+,s,s,qb-)
	case 0x3449:
	case 0x2A44: //2367: //(qb-,q+,s,s)
	case 0x9344:
	case 0x4A24: //2127: //(s,q+,qb-,s)
	case 0x4394:
	case 0x44A2: //489: //(s,s,q+,qb-)
	case 0x4439:
	case 0x244A: //1023: //(qb-,s,s,q+)
	case 0x9443:
	case 0xA244: //2361: //(q+,qb-,s,s)
	case 0x3944:
		return &ZeroF;
		break;

	default:// We return zero for all other helicity combinations
		return 0;
	}
}




template complex<R> (*A2QM1g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2QM1g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2QM1g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> (*A2QM1g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A2QM2g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A2QM2q_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1QM1qs_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1QM1q1gs_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A2s2q_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

template complex<R> (*A2QM2g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2QM2g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2QM2g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2QM2q_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2QM2q_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2QM2q_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1QM1qs_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1QM1qs_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1QM1qs_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1QM1q1gs_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1QM1q1gs_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1QM1q1gs_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2s2q_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2s2q_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2s2q_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

}
