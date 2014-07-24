/*
*  A0_1phi_Q_massive.cpp
*  BlackHat
*
*  Created by Darren Forde on 30/04/2010.
*  Copyright 2010 BlackHat Collaboration. All rights reserved.
*
*/
#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

#define _HIGGS_SCALAR_SCALAR_PARAM T(1) // T(12)ph,Q,q,q,Q higgs attached to scalar prop
#define _OTHER_FAC T(1)  // ph,g,Q,quarks higgs attached to gluon prop

using namespace std;

namespace BH {
	
// Defined in A0_1phi_eval.cpp
template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses);

#define S(i,j) (ep.p(i)->P()*ep.p(j)->P())
	
/*
 *
 *
 * The 1 complex higgs and 2 quarks and 2 massive quarks amplitudes
 *
 *
 */

//(Q+,q3+,qb3-,Qb+)
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i1)->Lt(),ep.p(i2)->L());
	momentum<complex<T> > mp03Q1=(mass/(f1.L()*ep.ref()->L()))*PfLLt(fn.Lt(),ep.ref()->L());
	momentum<complex<T> > mp03Q2=(mass/(fn.L()*ep.ref()->L()))*PfLLt(f1.Lt(),ep.ref()->L());
	complex<T> den1=m12.square()*m03.square();
	
	
//	_MESSAGE2("Case 1m:",complex<T>(0,-2)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/(den1));
	
	return complex<T>(0,-4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i1)->Lt(),ep.p(i2)->L());
	momentum<complex<T> > mp03Q1=(mass/(f1.L()*ep.ref()->L()))*PfLLt(fn.Lt(),ep.ref()->L());
	momentum<complex<T> > mp03Q2=(mass/(fn.L()*ep.ref()->L()))*PfLLt(f1.Lt(),ep.ref()->L());
	complex<T> den1=m12.square()*m03.square();
	
//		_MESSAGE2("Case 1p:",complex<T>(0,2)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}

	
//(Q+,q3+,qb3-,Qb-)	
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i1)->Lt(),ep.p(i2)->L());
	momentum<complex<T> > mp03Q1=PfLLt(f1.Lt(),fn.L());
	momentum<complex<T> > mp03Q2=(-mass2/((f1.L()*ep.ref()->L())*(ep.ref()->Lt()*fn.Lt())))*ep.ref()->P();
	complex<T> den1=m12.square()*m03.square();
	
//	_MESSAGE2("Case 2:",complex<T>(0,1)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}

	
//(Qb-,q3+,qb3-,Q+)
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i1)->Lt(),ep.p(i2)->L());
	momentum<complex<T> > mp03Q1=PfLLt(fn.Lt(),f1.L());
	momentum<complex<T> > mp03Q2=(-mass2/((fn.L()*ep.ref()->L())*(ep.ref()->Lt()*f1.Lt())))*ep.ref()->P();
	complex<T> den1=m12.square()*m03.square();
	
//		_MESSAGE2("Case 3:",complex<T>(0,1)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,-4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}

	
//(Qb-,qb3-,q3+,Q-)	
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));

	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i2)->Lt(),ep.p(i1)->L());
	momentum<complex<T> > mp03Q1=(mass/(f1.Lt()*ep.ref()->Lt()))*PfLLt(ep.ref()->Lt(),fn.L());
	momentum<complex<T> > mp03Q2=(mass/(fn.Lt()*ep.ref()->Lt()))*PfLLt(ep.ref()->Lt(),f1.L());
	complex<T> den1=m12.square()*m03.square();
	
//		_MESSAGE2("Case 4m:",complex<T>(0,-2)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,-4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));

	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i2)->Lt(),ep.p(i1)->L());
	momentum<complex<T> > mp03Q1=(mass/(f1.Lt()*ep.ref()->Lt()))*PfLLt(ep.ref()->Lt(),fn.L());
	momentum<complex<T> > mp03Q2=(mass/(fn.Lt()*ep.ref()->Lt()))*PfLLt(ep.ref()->Lt(),f1.L());
	complex<T> den1=m12.square()*m03.square();
	
//		_MESSAGE2("Case 4p:",complex<T>(0,2)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}

	
//(Qb-,qb3-,q3+,Q+)
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i2)->Lt(),ep.p(i1)->L());
	momentum<complex<T> > mp03Q1=PfLLt(fn.Lt(),f1.L());
	momentum<complex<T> > mp03Q2=(-mass2/((fn.L()*ep.ref()->L())*(ep.ref()->Lt()*f1.Lt())))*ep.ref()->P();
	complex<T> den1=m12.square()*m03.square();
	
//	_MESSAGE2("Case 5:",complex<T>(0,1)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}

	
//(Q+,qb3-,q3+,Qb-)	
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i2)->Lt(),ep.p(i1)->L());
	momentum<complex<T> > mp03Q1=PfLLt(f1.Lt(),fn.L());
	momentum<complex<T> > mp03Q2=(-mass2/((f1.L()*ep.ref()->L())*(ep.ref()->Lt()*fn.Lt())))*ep.ref()->P();
	complex<T> den1=m12.square()*m03.square();
	
//	_MESSAGE2("Case 6:",complex<T>(0,1)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,-4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}

	
//(Qb+,qb3-,q3+,Q+)	
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM7m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i2)->Lt(),ep.p(i1)->L());
	momentum<complex<T> > mp03Q1=(mass/(f1.L()*ep.ref()->L()))*PfLLt(fn.Lt(),ep.ref()->L());
	momentum<complex<T> > mp03Q2=(mass/(fn.L()*ep.ref()->L()))*PfLLt(f1.Lt(),ep.ref()->L());
	complex<T> den1=m12.square()*m03.square();
	
//		_MESSAGE2("Case 7m:",complex<T>(0,2)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM7p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i2)->Lt(),ep.p(i1)->L());
	momentum<complex<T> > mp03Q1=(mass/(f1.L()*ep.ref()->L()))*PfLLt(fn.Lt(),ep.ref()->L());
	momentum<complex<T> > mp03Q2=(mass/(fn.L()*ep.ref()->L()))*PfLLt(f1.Lt(),ep.ref()->L());
	complex<T> den1=m12.square()*m03.square();
	
//		_MESSAGE2("Case 7p:",complex<T>(0,-2)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,-4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}

	
//(Q-,q3+,qb3-,Qb-)
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i1)->Lt(),ep.p(i2)->L());
	momentum<complex<T> > mp03Q1=(mass/(f1.Lt()*ep.ref()->Lt()))*PfLLt(ep.ref()->Lt(),fn.L());
	momentum<complex<T> > mp03Q2=(mass/(fn.Lt()*ep.ref()->Lt()))*PfLLt(ep.ref()->Lt(),f1.L());
	complex<T> den1=m12.square()*m03.square();
	
//		_MESSAGE2("Case 8m:",complex<T>(0,2)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2q2QM8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	momentum<complex<T> > m12(ep.p(i1)->P()+ep.p(i2)->P());
	momentum<complex<T> > m03(ep.p(i3)->P()+ep.p(i0)->P());
	momentum<complex<T> > mp12=PfLLt(ep.p(i1)->Lt(),ep.p(i2)->L());
	momentum<complex<T> > mp03Q1=(mass/(f1.Lt()*ep.ref()->Lt()))*PfLLt(ep.ref()->Lt(),fn.L());
	momentum<complex<T> > mp03Q2=(mass/(fn.Lt()*ep.ref()->Lt()))*PfLLt(ep.ref()->Lt(),f1.L());
	complex<T> den1=m12.square()*m03.square();
	
//		_MESSAGE2("Case 8p:",complex<T>(0,-2)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1);
	
	return complex<T>(0,-4)*(((mp03Q1+mp03Q2)*mp12)*(m12*m03)-((mp03Q1+mp03Q2)*m12)*(mp12*m03))/den1;
}

	
	// Original
	

	//(Q+,q-,q2+,Q2+)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	//_MESSAGE8("1m: Q+,q-,q2+,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,1)*(ep.p(i2)->Lt()*fn.Lt())*(ep.ref()->L()*ep.p(i1)->L())/(T(8)*(ep.ref()->L()*f1.L())*ep.sp(i1,i0)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	//_MESSAGE8("1p: Q+,q-,q2+,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,-1)*(ep.p(i2)->Lt()*fn.Lt())*(ep.ref()->L()*ep.p(i1)->L())/(T(8)*(ep.ref()->L()*f1.L())*ep.sp(i1,i0)*ep.sp(i2,i3));
}

	 //(Q+,q-,q2+,Q2-)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_2p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("2p: Q+,q-,q2+,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*pow(mass,2)*complex<T>(0,1)*(ep.p(i2)->Lt()*ep.ref()->Lt())*(ep.ref()->L()*ep.p(i1)->L())/(T(8)*(fn.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f1.L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_2m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("2m: Q+,q-,q2+,Q2- : ",i0,", ",i1,", ",i2,", ",i3);	
	
	return -_HIGGS_SCALAR_SCALAR_PARAM*pow(mass,2)*complex<T>(0,-1)*(ep.p(i2)->Lt()*ep.ref()->Lt())*(ep.ref()->L()*ep.p(i1)->L())/(T(8)*(fn.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f1.L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}

	//(Q-,q-,q2+,Q2+)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("3: Q-,q-,q2+,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*complex<T>(0,-1)*mass*(ep.p(i2)->Lt()*fn.Lt())*(f1.L()*ep.p(i1)->L())/(T(8)*ep.sp(i1,i0)*ep.sp(i2,i3));
}

	//(Q+,q+,q2-,Q2+)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("4m: Q+,q+,q2-,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,1)*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/(T(8)*(fn.L()*ep.ref()->L())*ep.sp(i1,i0)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("4p: Q+,q+,q2-,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,-1)*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/(T(8)*(fn.L()*ep.ref()->L())*ep.sp(i1,i0)*ep.sp(i2,i3));
}

	//(Q+,q+,q2-,Q2-)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("5: Q+,q+,q2-,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*complex<T>(0,-1)*mass*(ep.p(i2)->L()*fn.L())*(f1.Lt()*ep.p(i1)->Lt())/(T(8)*ep.sp(i1,i0)*ep.sp(i2,i3));
}

	//(Q-,q+,q2-,Q2+)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_6p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("6p: Q-,q+,q2-,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*pow(mass,2)*complex<T>(0,1)*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/(T(8)*(ep.ref()->Lt()*f1.Lt())*(fn.L()*ep.ref()->L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_6m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("6m: Q-,q+,q2-,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*pow(mass,2)*complex<T>(0,-1)*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/(T(8)*(ep.ref()->Lt()*f1.Lt())*(fn.L()*ep.ref()->L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}

	//(Q-,q+,q2-,Q2-)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_7m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("7m: Q-,q+,q2-,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,1)*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*fn.L())/(T(8)*(ep.ref()->Lt()*f1.Lt())*ep.sp(i1,i0)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_7p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("7p: Q-,q+,q2-,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,-1)*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->L()*fn.L())/(T(8)*(ep.ref()->Lt()*f1.Lt())*ep.sp(i1,i0)*ep.sp(i2,i3));
}

	//(Q-,q-,q2+,Q2-)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("8m: Q-,q+,q2+,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,1)*(ep.ref()->Lt()*ep.p(i2)->Lt())*(f1.L()*ep.p(i1)->L())/(T(8)*(ep.ref()->Lt()*fn.Lt())*ep.sp(i1,i0)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("8p: Q-,q+,q2+,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,-1)*(ep.ref()->Lt()*ep.p(i2)->Lt())*(f1.L()*ep.p(i1)->L())/(T(8)*(ep.ref()->Lt()*fn.Lt())*ep.sp(i1,i0)*ep.sp(i2,i3));
}

	//(Q+,q+,q2+,Q2-)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_9m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("9m: Q+,q+,q2+,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,1)*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*ep.ref()->Lt())/(T(8)*(fn.Lt()*ep.ref()->Lt())*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_9p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("9p: Q+,q+,q2+,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,-1)*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*ep.ref()->Lt())/(T(8)*(fn.Lt()*ep.ref()->Lt())*ep.sp(i0,i1)*ep.sp(i2,i3));
}

	//(Q-,q+,q2+,Q2+)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_10m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("10m: Q-,q+,q2+,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,1)*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*fn.Lt())/(T(8)*(ep.ref()->Lt()*f1.Lt())*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_10p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("10p: Q-,q+,q2+,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,-1)*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*fn.Lt())/(T(8)*(ep.ref()->Lt()*f1.Lt())*ep.sp(i0,i1)*ep.sp(i2,i3));
}

	//(Q-,q-,q2-,Q2+)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_11m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("11m: Q-,q-,q2-,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,-1)*(f1.L()*ep.p(i1)->L())*(ep.p(i2)->L()*ep.ref()->L())/(T(8)*(fn.L()*ep.ref()->L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_11p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("11p: Q-,q-,q2-,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,1)*(f1.L()*ep.p(i1)->L())*(ep.p(i2)->L()*ep.ref()->L())/(T(8)*(fn.L()*ep.ref()->L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}

	//(Q+,q-,q2-,Q2-)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_12m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("12m: Q+,q-,q2-,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,-1)*(ep.ref()->L()*ep.p(i1)->L())*(ep.p(i2)->L()*fn.L())/(T(8)*(ep.ref()->L()*f1.L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_12p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("12p: Q+,q-,q2-,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*eval_param<T>::mass(masses.p(i0))*complex<T>(0,1)*(ep.ref()->L()*ep.p(i1)->L())*(ep.p(i2)->L()*fn.L())/(T(8)*(ep.ref()->L()*f1.L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}

	//(Q+,q-,q2-,Q2+)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_13p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("13p: Q+,q-,q2-,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*pow(mass,2)*complex<T>(0,-1)*(ep.ref()->L()*ep.p(i1)->L())*(ep.p(i2)->L()*ep.ref()->L())/(T(8)*(ep.ref()->L()*f1.L())*(fn.L()*ep.ref()->L())*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_13m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("13m: Q+,q-,q2-,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*pow(mass,2)*complex<T>(0,1)*(ep.ref()->L()*ep.p(i1)->L())*(ep.p(i2)->L()*ep.ref()->L())/(T(8)*(ep.ref()->L()*f1.L())*(fn.L()*ep.ref()->L())*ep.sp(i0,i1)*ep.sp(i2,i3));

}
	
	//(Q-,q+,q+,Q2-)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_14p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE10("14p: Q-,q+,q2+,Q2- : ",i0,", ",i1,", ",i2,", ",i3,", mass^2=",mass);
	return -_HIGGS_SCALAR_SCALAR_PARAM*pow(mass,2)*complex<T>(0,1)*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*ep.ref()->Lt())/(T(8)*(ep.ref()->Lt()*f1.Lt())*(fn.Lt()*ep.ref()->Lt())*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2f_14m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE10("14m: Q-,q+,q2+,Q2- : ",i0,", ",i1,", ",i2,", ",i3,", mass^2=",mass);
	return -_HIGGS_SCALAR_SCALAR_PARAM*pow(mass,2)*complex<T>(0,-1)*(ep.ref()->Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*ep.ref()->Lt())/(T(8)*(ep.ref()->Lt()*f1.Lt())*(fn.Lt()*ep.ref()->Lt())*ep.sp(i0,i1)*ep.sp(i2,i3));
}

	//(Q-,q-,q2-,Q2-)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2fp_15_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("15p: Q-,q-,q2-,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*complex<T>(0,-1)*(f1.L()*ep.p(i1)->L())*(ep.p(i2)->L()*fn.L())/(T(8)*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2fm_15_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("15m: Q-,q-,q2-,Q2- : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*complex<T>(0,-1)*(f1.L()*ep.p(i1)->L())*(ep.p(i2)->L()*fn.L())/(T(8)*ep.sp(i0,i1)*ep.sp(i2,i3));}

	//(Q+,q+,q2+,Q2+)
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2fp_16_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("16p: Q+,q+,q2+,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*complex<T>(0,-1)*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*fn.Lt())/(T(8)*ep.sp(i0,i1)*ep.sp(i2,i3));
}
template <int i0, int i1, int i2, int i3,  class T> complex<T>  A1ph2q2QM2fm_16_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));//const complex<T> mass=ep.ref_mass();
															 // We use the negative helicity leg as the reference vector
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> fn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	//_MESSAGE8("16m: Q+,q+,q2+,Q2+ : ",i0,", ",i1,", ",i2,", ",i3);
	return -_HIGGS_SCALAR_SCALAR_PARAM*mass*complex<T>(0,-1)*(f1.Lt()*ep.p(i1)->Lt())*(ep.p(i2)->Lt()*fn.Lt())/(T(8)*ep.sp(i0,i1)*ep.sp(i2,i3));}

template <class T, int QMPOS, int QPOS, int QPOS2, int QMPOS2> complex<T> (*A1ph2q2QM_Tree_Ptr_non_higgs_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
//	cout << hex << " hqQq2Q2:" << hc << dec << endl;
	switch(hc){
		case 0x07090805: //8631://(Q+,q3+,qb3-,Qb+)
		case 0x0D03020B: //12131://Qb2+,q+,qb-,Q2+
			return &A1ph2q2QM1m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x05090807:
		case 0x0B03020D:
			return &A1ph2q2QM1p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x05070908: //4451://(Qb+,Q+,q3+,qb3-)
		case 0x0B0D0302: //1571://(Q2+,Qb2+,q+,qb-)
			return &A1ph2q2QM1m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x07050908:
		case 0x0D0B0302:
			return &A1ph2q2QM1p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x08050709: //
		case 0x020B0D03: //
			return &A1ph2q2QM1m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x08070509:
		case 0x020D0B03:
			return &A1ph2q2QM1p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x09080507: //10201://(q3+,qb3-,Q+,Qb+)
		case 0x03020B0D: //13071://(qb+,q-,Q2+,Qb2+)
			return &A1ph2q2QM1m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
		case 0x09080705:
		case 0x03020D0B:
			return &A1ph2q2QM1p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
			
		case 0x07090804: //7400: //(Q+,q3+,qb3-,Qb-)
		case 0x05090806:
		case 0x0D03020A: //10700://Qb2+,q+,qb-,Q2-
		case 0x0B03020C:
			return &A1ph2q2QM2_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x04070908: //4450://(Qb-,Q+,q3+,qb3-)
		case 0x06050908:
		case 0x0A0D0302: //1570://(Q2-,Qb2+,q+,qb-)
		case 0x0C0B0302:
			return &A1ph2q2QM2_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x08040709: //
		case 0x08060509:
		case 0x020A0D03: //
		case 0x020C0B03:
			return &A1ph2q2QM2_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x09080407: //10070://(q3+,qb3-,Q-,Qb+)
		case 0x09080605:
		case 0x03020A0D: //12850://qb+,q-,Q2-,Qb2+
		case 0x03020C0B:
			return &A1ph2q2QM2_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
			
		case 0x06090805: //8630://(Qb-,q3+,qb3-,Q+)
		case 0x04090807:
		case 0x0C03020B: //12130://Qb2-,q+,qb-,Q2+
		case 0x0A03020D:
			return &A1ph2q2QM3_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x05060908: //4440://(Qb+,Q-,q3+,qb3-)
		case 0x07040908:
		case 0x0B0C0302: //1560://(Q2+,Qb2-,q+,qb-)
		case 0x0D0A0302:
			return &A1ph2q2QM3_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x08050609: //
		case 0x08070409:
		case 0x020B0C03: //
		case 0x020D0A03:
			return &A1ph2q2QM3_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x09080506: //7760://(q3+,qb3-,Q+,Qb-)
		case 0x09080704:
		case 0x03020B0C: //11640://(qb+,q-,Q2+,Qb2-)
		case 0x03020D0A:
			return &A1ph2q2QM3_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
			
		case 0x06080904: //7408://(Qb-,qb3-,q3+,Q-)
		case 0x0C02030A: //10808://Qb2-,q-,qb+,Q2-
			return &A1ph2q2QM4m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x04080906:
		case 0x0A02030C:
			return &A1ph2q2QM4p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x04060809: //2768://(Qb3-,Q3-,qb-,q+)
		case 0x0A0C0203: //4648://(Qb-,Q-,qb3-,q3+)
			return &A1ph2q2QM4m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x06040809:
		case 0x0C0A0203:
			return &A1ph2q2QM4p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x09040608: //2378://(q+,Q3-,Qb3-,qb-)
		case 0x030A0C02: //4678://(q3+,Q-,Qb-,qb3-)
			return &A1ph2q2QM4m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x09060408:
		case 0x030C0A02:
			return &A1ph2q2QM4p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x08090406: //7648://(qb3-,q3+,Q-,Qb-)
		case 0x02030A0C: //11538://(qb-,q+,Q2-,Qb2-)
			return &A1ph2q2QM4m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
		case 0x08090604:
		case 0x02030C0A:
			return &A1ph2q2QM4p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
			
		case 0x06080905: //8740://(Qb-,qb3-,q3+,Q+)
		case 0x04080907:
		case 0x0C02030B: //12240://Qb2-,q-,qb+,Q2+
		case 0x0A02030D:
			return &A1ph2q2QM5_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x05060809: //2770://(Qb3+,Q3-,qb-,q+)
		case 0x07040809:
		case 0x0B0C0203: //4650://(Qb+,Q-,qb3-,q3+)
		case 0x0D0A0203:
			return &A1ph2q2QM5_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x09050608: //2400://(q+,Q3+,Qb3-,qb-)
		case 0x09070408:
		case 0x030B0C02: //4700://(q3+,Q+,Qb-,qb3-)
		case 0x030D0A02:
			return &A1ph2q2QM5_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x08090506: //7770://(qb3-,q3+,Q+,Qb-)
		case 0x08090704:
		case 0x02030B0C: //11650://(qb-,q+,Q2+,Qb2-)
		case 0x02030D0A:
			return &A1ph2q2QM5_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
			
		case 0x07080904: //7410://(Q+,qb3-,q3+,Qb-)
		case 0x05080906:
		case 0x0D02030A: //10810://Qb2+,q-,qb+,Q2-
		case 0x0B02030C:
			return &A1ph2q2QM6_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x04070809: //2780://(Qb3-,Q3+,qb-,q+)
		case 0x06050809:
		case 0x0A0D0203: //4660://(Qb-,Q+,qb3-,q3+)
		case 0x0C0B0203:
			return &A1ph2q2QM6_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x09040708: //2410://(q+,Qb3-,Q3+,qb-)
		case 0x09060508:
		case 0x030A0D02: //4810://(q3+,Qb-,Q+,qb3-)
		case 0x030C0B02:
			return &A1ph2q2QM6_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x08090407: //10080://(qb3-,q3+,Q-,Qb+)
		case 0x08090605:
		case 0x02030A0D: //12860://(qb-,q+,Q2-,Qb2+)
		case 0x02030C0B:
			return &A1ph2q2QM6_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
			
		case 0x07080905: //8741://(Qb+,qb3-,q3+,Q+)
		case 0x0D02030B: //12241://Qb2+,qb-,q+,Q2+
			return &A1ph2q2QM7m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x05080907:
		case 0x0B02030D:
			return &A1ph2q2QM7p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x05070809: //2781://(Qb3+,Q3+,qb-,q+)
		case 0x0B0D0203: //4661://(Qb+,Q+,qb3-,q3+)
			return &A1ph2q2QM7m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x07050809:
		case 0x0D0B0203:
			return &A1ph2q2QM7p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x09050708: //2421://(q+,Q3+,Qb3+,qb-)
		case 0x030B0D02: //4821://(q3+,Q+,Qb+,qb3-)
			return &A1ph2q2QM7m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x09070508:
		case 0x030D0B02:
			return &A1ph2q2QM7p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x08090507: //10211://(qb3-,q3+,Q+,Qb+)
		case 0x02030B0D: //13081://(qb-,q+,Q2+,Qb2+)
			return &A1ph2q2QM7m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
		case 0x08090705:
		case 0x02030D0B:
			return &A1ph2q2QM7p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
			
		case 0x06090804: //7388://(Q-,q3+,qb3-,Qb-)
		case 0x0C03020A: //10688://Q2-,qb+,q-,Qb2-
			return &A1ph2q2QM8m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x04090806:
		case 0x0A03020C:
			return &A1ph2q2QM8p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
			break;
		case 0x04060908: //4448://(Qb-,Q-,q3+,qb3-)
		case 0x0A0C0302: //1558://(Q2-,Qb2-,q+,qb-)
			return &A1ph2q2QM8m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x06040908:
		case 0x0C0A0302:
			return &A1ph2q2QM8p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
			break;
		case 0x08040609: //
		case 0x020A0C03: //
			return &A1ph2q2QM8m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x08060409:
		case 0x020C0A03:
			return &A1ph2q2QM8p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
			break;
		case 0x09080406: //7648://(q3+,qb3-,Q-,Qb-)
		case 0x03020A0C: //11528://(qb+,q-,Q2-,Qb2-)
			return &A1ph2q2QM8m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
		case 0x09080604:
		case 0x03020C0A:
			return &A1ph2q2QM8p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
			break;
			
	// COMMENTED OUT FOR HIGGS DEBUGGING		
//			// Original
//			
//			//1
//			
//		case 0x0702090B: //(Qb+,q-,qb3+,Q3+)
//		case 0x0D080305: //(Qb3+,q3-,qb+,Q+)
//		case 0x0702090D: //(Qb+,q-,q3+,Qb3+)
//		case 0x0D080307: //(Qb+,q-,q3+,Qb3+)
//			return &A1ph2q2QM2f_1m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0502090D: //(Q+,qb-,q3+,Qb3+)
//		case 0x0B080307: //(Q+,qb-,q3+,Qb3+)
//		case 0x0502090B: //(Q+,qb-,qb3+,Q3+)
//		case 0x0B080305: //(Q3+,qb3-,qb+,Q+)
//			return &A1ph2q2QM2f_1p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0B070209: //12471://(Q+,qb-,q3+,Qb3+)
//		case 0x050D0803: //8501://(Qb3+,q3-,qb+,Q+)
//		case 0x0D070209:
//		case 0x070D0803:
//			return &A1ph2q2QM2f_1m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0D050209:
//		case 0x070B0803:
//		case 0x0B050209:
//		case 0x050B0803:
//			return &A1ph2q2QM2f_1p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x090B0702: //12471://(Q+,qb-,q3+,Qb3+)
//		case 0x03050D08: //8501://(Qb3+,q3-,qb+,Q+)
//		case 0x090D0702:
//		case 0x03070D08:
//			return &A1ph2q2QM2f_1m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x090D0502:
//		case 0x03070B08:
//		case 0x090B0502:
//		case 0x03050B08:
//			return &A1ph2q2QM2f_1p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x02090B07: //12471://(Q+,qb-,q3+,Qb3+)
//		case 0x0803050D: //8501://(Qb3+,q3-,qb+,Q+)
//		case 0x02090D07:
//		case 0x0803070D:
//			return &A1ph2q2QM2f_1m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x02090D05:
//		case 0x0803070B:
//		case 0x02090B05:
//		case 0x0803050B:
//			return &A1ph2q2QM2f_1p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//						
//			//2
//			
//		case 0x0502090C: //(Q+,qb-,q3+,Qb3-)
//		case 0x0B080306: //(Q3+,qb3-,q+,Qb-)
//		case 0x0702090A: //(Qb+,q-,qb3+,Q3-)
//		case 0x0D080304: //(Qb3+,q3-,qb+,Q-)
//			return &A1ph2q2QM2f_2m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0502090A: //(Q+,qb-,qb3+,Q3-)
//		case 0x0B080304: //(Q3+,qb3-,qb+,Q-)
//		case 0x0702090C: //(Qb+,q-,q3+,Qb3-)
//		case 0x0D080306: //(Qb3+,q3-,q+,Qb-)
//			return &A1ph2q2QM2f_2p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0C050209: //
//		case 0x060B0803:
//		case 0x0A070209:
//		case 0x040D0803:
//			return &A1ph2q2QM2f_2m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0A050209: //
//		case 0x040B0803:
//		case 0x0C070209:
//		case 0x060D0803:
//			return &A1ph2q2QM2f_2p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x090C0502: //
//		case 0x03060B08:
//		case 0x090A0702:
//		case 0x03040D08:
//			return &A1ph2q2QM2f_2m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x090A0502: //
//		case 0x03040B08:
//		case 0x090C0702:
//		case 0x03060D08:
//			return &A1ph2q2QM2f_2p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x02090C05: //
//		case 0x0803060B:
//		case 0x02090A07:
//		case 0x0803040D:
//			return &A1ph2q2QM2f_2m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x02090A05: //
//		case 0x0803040B:
//		case 0x02090C07:
//		case 0x0803060D:
//			return &A1ph2q2QM2f_2p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//						
//			//3
//			
//		case 0x0602090B: //(Qb-,q-,qb3+,Q3+)
//		case 0x0602090D: //(Qb-,q-,q3+,Qb3+)
//		case 0x0402090B: //(Q-,qb-,qb3+,Q3+)
//		case 0x0402090D: //(Q-,qb-,q3+,Qb3+)
//		case 0x0C080305: //(Qb3-,q3-,qb+,Q+)
//		case 0x0C080307: //(Qb3-,q3-,q+,Qb+)
//		case 0x0A080305: //(Q3-,qb3-,qb+,Q+)
//		case 0x0A080307: //(Q3-,qb3-,q+,Qb+)
//			return &A1ph2q2QM2f_3_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0B060209:
//		case 0x0D060209:
//		case 0x0B040209:
//		case 0x0D040209:
//		case 0x050C0803:
//		case 0x070C0803:
//		case 0x050A0803:
//		case 0x070A0803:
//			return &A1ph2q2QM2f_3_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x090B0602:
//		case 0x090D0602:
//		case 0x090B0402:
//		case 0x090D0402:
//		case 0x03050C08:
//		case 0x03070C08:
//		case 0x03050A08:
//		case 0x03070A08:
//			return &A1ph2q2QM2f_3_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x02090B06:
//		case 0x02090D06:
//		case 0x02090B04:
//		case 0x02090D04:
//		case 0x0803050C:
//		case 0x0803070C:
//		case 0x0803050A:
//		case 0x0803070A:
//			return &A1ph2q2QM2f_3_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//						
//			//4
//			
//			
//		case 0x0703080B: //(Qb+,q+,qb3-,Q3+)
//		case 0x0D090205: //(Qb3+,q3+,qb-,Q+)
//		case 0x0503080B: //(Q+,qb+,qb3-,Q3+)
//		case 0x0B090205: //(Q3+,qb3+,qb-,Q+)
//			return &A1ph2q2QM2f_4m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0503080D: //(Q+,qb+,q3-,Qb3+)
//		case 0x0B090207: //(Q3+,qb3+,q-,Qb+)
//		case 0x0703080D: //(Qb+,q+,q3-,Qb3+)
//		case 0x0D090207: //(Qb3+,q3+,q-,Qb+)
//			return &A1ph2q2QM2f_4p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0B070308: //12361://(Qb+,q+,qb3-,Q3+)
//		case 0x050D0902: //8481://(Qb3+,q3+,qb-,Q+)
//		case 0x0B050308: //(Q+,qb+,qb3-,Q3+)
//		case 0x050B0902:
//			return &A1ph2q2QM2f_4m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0D050308:
//		case 0x070B0902:
//		case 0x0D070308:
//		case 0x070D0902:
//			return &A1ph2q2QM2f_4p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x080B0703: //12361://(Qb+,q+,qb3-,Q3+)
//		case 0x02050D09: //8481://(Qb3+,q3+,qb-,Q+)
//		case 0x080B0503: //(Q+,qb+,qb3-,Q3+)
//		case 0x02050B09:
//			return &A1ph2q2QM2f_4m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x080D0503:
//		case 0x02070B09:
//		case 0x080D0703:
//		case 0x02070D09:
//			return &A1ph2q2QM2f_4p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03080B07: //12361://(Qb+,q+,qb3-,Q3+)
//		case 0x0902050D: //8481://(Qb3+,q3+,qb-,Q+)
//		case 0x03080B05: //(Q+,qb+,qb3-,Q3+)
//		case 0x0902050B:
//			return &A1ph2q2QM2f_4m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x03080D05:
//		case 0x0902070B:
//		case 0x03080D07:
//		case 0x0902070D:
//			return &A1ph2q2QM2f_4p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//5
//			
//		case 0x0503080A: //(Q+,qb+,q3-,Qb3-)
//		case 0x0703080A: //(Qb+,q+,q3-,Qb3-)
//		case 0x0B090204: //(Q3+,qb3+,qb-,Q-)
//		case 0x0D090204: //(Qb3+,q3+,qb-,Q-)
//		case 0x0503080C: //(Q+,qb+,qb3-,Q3-)
//		case 0x0703080C: //(Qb+,q+,qb3-,Q3-)
//		case 0x0B090206: //(Q3+,qb3+,q-,Qb-)
//		case 0x0D090206: //(Qb3+,q3+,q-,Qb-)
//			return &A1ph2q2QM2f_5_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0A050308:
//		case 0x0A070308:
//		case 0x040B0902:
//		case 0x040D0902:
//		case 0x0C050308:
//		case 0x0C070308:
//		case 0x060B0902:
//		case 0x060D0902:
//			return &A1ph2q2QM2f_5_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x080A0503:
//		case 0x080A0703:
//		case 0x02040B09:
//		case 0x02040D09:
//		case 0x080C0503:
//		case 0x080C0703:
//		case 0x02060B09:
//		case 0x02060D09:
//			return &A1ph2q2QM2f_5_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03080A05:
//		case 0x03080A07:
//		case 0x0902040B:
//		case 0x0902040D:
//		case 0x03080C05:
//		case 0x03080C07:
//		case 0x0902060B:
//		case 0x0902060D:
//			return &A1ph2q2QM2f_5_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//6
//			
//		case 0x0603080B: //(Qb-,q+,qb3-,Q3+)
//		case 0x0C090205: //(Qb3-,q3+,qb-,Q+)
//		case 0x0403080D: //(Q-,qb+,q3-,Qb3+)
//		case 0x0A090207: //(Q3-,qb3+,q-,Qb+)
//			return &A1ph2q2QM2f_6m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0603080D: //(Qb-,q+,q3-,Qb3+)
//		case 0x0C090207: //(Qb3-,q3+,q-,Qb+)
//		case 0x0403080B: //(Q-,qb+,qb3-,Q3+)
//		case 0x0A090205: //(Q3-,qb3+,qb-,Q+)
//			return &A1ph2q2QM2f_6p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0B060308: //
//		case 0x050C0902:
//		case 0x0D040308:
//		case 0x070A0902:
//			return &A1ph2q2QM2f_6m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0D060308:
//		case 0x070C0902:
//		case 0x0B040308:
//		case 0x050A0902:
//			return &A1ph2q2QM2f_6p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x080B0603: //
//		case 0x02050C09:
//		case 0x080D0403:
//		case 0x02070A09:
//			return &A1ph2q2QM2f_6m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x080D0603:
//		case 0x02070C09:
//		case 0x080B0403:
//		case 0x02050A09:
//			return &A1ph2q2QM2f_6p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03080B06: //
//		case 0x0902050C:
//		case 0x03080D04:
//		case 0x0902070A:
//			return &A1ph2q2QM2f_6m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x03080D06:
//		case 0x0902070C:
//		case 0x03080B04:
//		case 0x0902050A:
//			return &A1ph2q2QM2f_6p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//7
//			
//		case 0x0603080A: //(Qb-,q+,qb3-,Q3-)
//		case 0x0C090204: //(Qb3-,q3+,qb-,Q-)
//		case 0x0603080C: //(Qb-,q+,q3-,Qb3-)
//		case 0x0C090206: //(Qb3-,q3+,q-,Qb-)
//			return &A1ph2q2QM2f_7m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0403080C:
//		case 0x0A090206:
//		case 0x0403080A:
//		case 0x0A090204:
//			return &A1ph2q2QM2f_7p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0A060308: //(Qb-,q+,qb3-,Q3-)
//		case 0x040C0902:
//		case 0x0C060308:
//		case 0x060C0902:
//			return &A1ph2q2QM2f_7m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0C040308:
//		case 0x060A0902:
//		case 0x0A040308:
//		case 0x040A0902:
//			return &A1ph2q2QM2f_7p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x080A0603: //(Qb-,q+,qb3-,Q3-)
//		case 0x02040C09:
//		case 0x080C0603:
//		case 0x02060C09:
//			return &A1ph2q2QM2f_7m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x080C0403:
//		case 0x02060A09:
//		case 0x080A0403:
//		case 0x02040A09:
//			return &A1ph2q2QM2f_7p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03080A06: //(Qb-,q+,qb3-,Q3-)
//		case 0x0902040C:
//		case 0x03080C06:
//		case 0x0902060C:
//			return &A1ph2q2QM2f_7m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x03080C04:
//		case 0x0902060A:
//		case 0x03080A04:
//		case 0x0902040A:
//			return &A1ph2q2QM2f_7p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//8
//			
//		case 0x0602090A: //(Qb-,q-,qb3+,Q3-)
//		case 0x0C080304: //(Q3-,qb3-,q+,Qb-)
//		case 0x0402090A: //(Q-,qb-,qb3+,Q3-)
//		case 0x0A080304: //(Q3-,qb3-,q+,Qb-)
//			return &A1ph2q2QM2f_8m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0402090C: //(Q-,qb-,q3+,Qb3-)
//		case 0x0A080306: //(Q3-,qb3-,q+,Qb-)
//		case 0x0602090C: //(Qb-,q-,qb3+,Qb3-)
//		case 0x0C080306: //(Qb3-,q3-,q+,Qb-)
//			return &A1ph2q2QM2f_8p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0A060209: //(Qb-,q-,qb3+,Q3-)
//		case 0x040C0803:
//		case 0x0A040209:
//		case 0x040A0803:
//			return &A1ph2q2QM2f_8m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0C040209:
//		case 0x060A0803:
//		case 0x0C060209:
//		case 0x060C0803:
//			return &A1ph2q2QM2f_8p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x090A0602: //(Qb-,q-,qb3+,Q3-)
//		case 0x03040C08:
//		case 0x090A0402:
//		case 0x03040A08:
//			return &A1ph2q2QM2f_8m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x090C0402:
//		case 0x03060A08:
//		case 0x090C0602:
//		case 0x03060C08:
//			return &A1ph2q2QM2f_8p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x02090A06: //(Qb-,q-,qb3+,Q3-)
//		case 0x0803040C:
//		case 0x02090A04:
//		case 0x0803040A:
//			return &A1ph2q2QM2f_8m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x02090C04:
//		case 0x0803060A:
//		case 0x02090C06:
//		case 0x0803060C:
//			return &A1ph2q2QM2f_8p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//9
//			
//		case 0x0503090C: //(Q+,qb+,q2+,Qb2-)
//		case 0x0B090306: //(Q2+,qb2+,q+,Qb-)
//		case 0x0703090C: //(Qb+,q+,q2+,Qb2-)
//		case 0x0D090306: //(Qb2+,q2+,q+,Qb-)
//			return &A1ph2q2QM2f_9p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0D090304: //(Qb2+,q2+,qb+,Q-)
//		case 0x0703090A: //(Qb+,q+,qb2+,Q2-)
//		case 0x0B090304: //(Q2+,qb2+,qb+,Q-)
//		case 0x0503090A: //(Q+,qb+,qb2+,Q2-)
//			return &A1ph2q2QM2f_9m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0C050309: //
//		case 0x060B0903:
//		case 0x0C070309:
//		case 0x060D0903:
//			return &A1ph2q2QM2f_9p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x040D0903:
//		case 0x0A070309:
//		case 0x040B0903:
//		case 0x0A050309:
//			return &A1ph2q2QM2f_9m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x090C0503: //
//		case 0x03060B09:
//		case 0x090C0703:
//		case 0x03060D09:
//			return &A1ph2q2QM2f_9p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03040D09:
//		case 0x090A0703:
//		case 0x03040B09:
//		case 0x090A0503:
//			return &A1ph2q2QM2f_9m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03090C05: //
//		case 0x0903060B:
//		case 0x03090C07:
//		case 0x0903060D:
//			return &A1ph2q2QM2f_9p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x0903040D:
//		case 0x03090A07:
//		case 0x0903040B:
//		case 0x03090A05:
//			return &A1ph2q2QM2f_9m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//10
//			
//		case 0x0603090B: //(Qb-,q+,qb2+,Q2+)
//		case 0x0C090305: //(Qb2-,q2+,qb+,Q+)
//		case 0x0603090D: //(Qb-,q+,q2+,Qb2+)
//		case 0x0C090307: //(Qb2-,q2+,q+,Qb+)
//			return &A1ph2q2QM2f_10m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0403090D: //(Q-,qb+,q2+,Qb2+)
//		case 0x0A090307: //(Q2-,qb2+,q+,Qb+)
//		case 0x0403090B: //(Q-,qb+,qb2+,Q2+)
//		case 0x0A090305: //(Q2-,qb2+,qb+,Q+)
//			return &A1ph2q2QM2f_10p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0B060309: //
//		case 0x050C0903:
//		case 0x0D060309:
//		case 0x070C0903:
//			return &A1ph2q2QM2f_10m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0D040309:
//		case 0x070A0903:
//		case 0x0B040309:
//		case 0x050A0903:
//			return &A1ph2q2QM2f_10p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x090B0603: //
//		case 0x03050C09:
//		case 0x090D0603:
//		case 0x03070C09:
//			return &A1ph2q2QM2f_10m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x090D0403:
//		case 0x03070A09:
//		case 0x090B0403:
//		case 0x03050A09:
//			return &A1ph2q2QM2f_10p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03090B06: //
//		case 0x0903050C:
//		case 0x03090D06:
//		case 0x0903070C:
//			return &A1ph2q2QM2f_10m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x03090D04:
//		case 0x0903070A:
//		case 0x03090B04:
//		case 0x0903050A:
//			return &A1ph2q2QM2f_10p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//11
//			
//		case 0x0402080D: //(Q-,qb-,q2-,Qb2+)
//		case 0x0A080207: //(Q2-,qb2-,q-,Qb+)
//		case 0x0602080D: //(Qb-,q-,q2-,Qb2+)
//		case 0x0C080207: //(Qb2-,q-,q-,Qb+)
//			return &A1ph2q2QM2f_11m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0402080B: //(Q-,qb-,qb2-,Q+)
//		case 0x0A080205: //(Q2-,qb2-,qb-,Q+)
//		case 0x0602080B: //(Qb-,q-,qb2-,Q2+)
//		case 0x0C080205: //(Qb2-,q2-,qb-,Q+)
//			return &A1ph2q2QM2f_11p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0D040208: //
//		case 0x070A0802:
//		case 0x0D060208:
//		case 0x070C0802:
//			return &A1ph2q2QM2f_11m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0B040208:
//		case 0x050A0802:
//		case 0x0B060208:
//		case 0x050C0802:
//			return &A1ph2q2QM2f_11p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x080D0402: //
//		case 0x02070A08:
//		case 0x080D0602:
//		case 0x02070C08:
//			return &A1ph2q2QM2f_11m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x080B0402:
//		case 0x02050A08:
//		case 0x080B0602:
//		case 0x0205C008:
//			return &A1ph2q2QM2f_11p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x02080D04://
//		case 0x0802070A:
//		case 0x02080D06:
//		case 0x0802070C:
//			return &A1ph2q2QM2f_11m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x02080B04:
//		case 0x0802050A:
//		case 0x02080B06:
//		case 0x0802050C:
//			return &A1ph2q2QM2f_11p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			
//		case 0x0502080C: //(Q+,qb-,q2-,Qb2-)
//		case 0x0B080206: //(Q2+,qb2-,q-,Qb-)
//		case 0x0502080A: //(Q+,qb-,qb2-,Q2-)
//		case 0x0B080204: //(Q2+,qb2-,qb-,Q-)
//			return &A1ph2q2QM2f_12m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0702080A: //(Qb+,q-,qb2-,Q2-)
//		case 0x0D080204: //(Qb2+,q2-,qb-,Q-)
//		case 0x0702080C: //(Qb+,q-,q2-,Qb2-)
//		case 0x0D080206: //(Qb2+,q2-,q-,Qb-)
//			return &A1ph2q2QM2f_12p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0C050208://
//		case 0x060B0802:
//		case 0x0A050208:
//		case 0x040B0802:
//			return &A1ph2q2QM2f_12m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0A070208:
//		case 0x040D0802:
//		case 0x0C070208:
//		case 0x060D0802:
//			return &A1ph2q2QM2f_12p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x080C0502://
//		case 0x02060B08:
//		case 0x080A0502:
//		case 0x02040B08:
//			return &A1ph2q2QM2f_12m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x080A0702:
//		case 0x02040D08:
//		case 0x080C0702:
//		case 0x02060D08:
//			return &A1ph2q2QM2f_12p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x02080C05://
//		case 0x0802060B:
//		case 0x02080A05:
//		case 0x0802040B:
//			return &A1ph2q2QM2f_12m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x02080A07:
//		case 0x0802040D:
//		case 0x02080C07:
//		case 0x0802060D:
//			return &A1ph2q2QM2f_12p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//13
//			
//		case 0x0702080B: //(Qb+,q-,qb2-,Q2+)
//		case 0x0D080205: //(Qb2+,q2-,qb-,Q+)
//		case 0x0502080D: //(Q+,qb-,q2-,Qb2+)
//		case 0x0B080207: //(Q2+,qb2-,q-,Qb+)
//			return &A1ph2q2QM2f_13p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0702080D: //(Qb+,q-,q2-,Qb2+)
//		case 0x0D080207: //(Qb2+,q2-,q-,Qb+)
//		case 0x0502080B: //(Q+,qb-,qb-,Q2+)
//		case 0x0B080205: //(Q2+,qb2-,qb-,Q+)
//			return &A1ph2q2QM2f_13m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0B070208: //
//		case 0x050D0802:
//		case 0x0D050208:
//		case 0x070B0802:
//			return &A1ph2q2QM2f_13p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0D070208:
//		case 0x070D0802:
//		case 0x0B050208:
//		case 0x050B0802:
//			return &A1ph2q2QM2f_13m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x080B0702: //
//		case 0x02050D08:
//		case 0x080D0502:
//		case 0x02070B08:
//			return &A1ph2q2QM2f_13p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x080D0702:
//		case 0x02070D08:
//		case 0x080B0502:
//		case 0x02050B08:
//			return &A1ph2q2QM2f_13m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x02080B07: //
//		case 0x0802050D:
//		case 0x02080D05:
//		case 0x0802070B:
//			return &A1ph2q2QM2f_13p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x02080D07:
//		case 0x0802070D:
//		case 0x0208B005:
//		case 0x0802050B:
//			return &A1ph2q2QM2f_13m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//14
//			
//		case 0x0403090C: //(Q-,qb+,q+,Qb2-)
//		case 0x0A090306: //(Q2-,qb2+,q+,Qb-)
//		case 0x0603090A: //(Qb-,q+,qb2+,Q2-)
//		case 0x0C090304: //(Qb2-,q2+,qb+,Q-)
//			return &A1ph2q2QM2f_14m_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0603090C: //(Qb-,q+,q2+,Qb2-)
//		case 0x0C090306: //(Qb2-,q2+,q+,Qb-)
//		case 0x0403090A: //(Q-,qb+,qb2+,Q2-)
//		case 0x0A090304: //(Q2-,qb2+,qb+,Q-)
//			return &A1ph2q2QM2f_14p_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0C040309: //
//		case 0x060A0903:
//		case 0x0A060309:
//		case 0x040C0903:
//			return &A1ph2q2QM2f_14m_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0C060309:
//		case 0x060C0903:
//		case 0x0A040309:
//		case 0x040A0903:
//			return &A1ph2q2QM2f_14p_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x090C0403: //
//		case 0x03060A09:
//		case 0x090A0603:
//		case 0x03040C09:
//			return &A1ph2q2QM2f_14m_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x090C0603:
//		case 0x03060C09:
//		case 0x090A4003:
//		case 0x03040A09:
//			return &A1ph2q2QM2f_14p_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03090C04: //
//		case 0x0903060A:
//		case 0x03090A06:
//		case 0x0903040C:
//			return &A1ph2q2QM2f_14m_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x03090C06:
//		case 0x0903060C:
//		case 0x03090A04:
//		case 0x0903040A:
//			return &A1ph2q2QM2f_14p_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//15
//			
//		case 0x0602080A: //(Q-,qb-,q2-,Qb2-)
//		case 0x0C080204: //(Qb2-,q2-,qb-,Q-)
//		case 0x0402080C:
//		case 0x0A080206:
//			return &A1ph2q2QM2fp_15_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0602080C:
//		case 0x0C080206:
//		case 0x0402080A:
//		case 0x0A080204:
//			return &A1ph2q2QM2fm_15_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0A060208:
//		case 0x040C0802:
//		case 0x0C040208:
//		case 0x060A0802:
//			return &A1ph2q2QM2fp_15_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0C060208:
//		case 0x060C0802:
//		case 0x0A040208:
//		case 0x040A0802:
//			return &A1ph2q2QM2fm_15_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x080A0602:
//		case 0x02040C08:
//		case 0x080C0402:
//		case 0x02060A08:
//			return &A1ph2q2QM2fp_15_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x080C0602:
//		case 0x02060C08:
//		case 0x080A0402:
//		case 0x02040A08:
//			return &A1ph2q2QM2fm_15_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x02080A06:
//		case 0x0802040C:
//		case 0x02080C04:
//		case 0x0802060A:
//			return &A1ph2q2QM2fp_15_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x02080C06:
//		case 0x0802060C:
//		case 0x02080A04:
//		case 0x0802040A:
//			return &A1ph2q2QM2fm_15_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//			
//			//16
//			
//		case 0x0703090B: //(Qb+,q+,qb2+,Q2+)
//		case 0x0D090305: //(Qb2+,q2+,qb+,Q+)
//		case 0x0503090D: //(Q+,qb+,q2+,Qb2+)
//		case 0x0B093007: //(Q2+,qb2+,q+,Qb+)
//			return &A1ph2q2QM2fp_16_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0703090D: //(Qb+,q+,q2+,Qb2+)
//		case 0x0D090307: //(Qb2+,q2+,q+,Qb+)
//		case 0x0503090B: //(Q+,qb+,qb2+,Q2+)
//		case 0x0B090305: //(Q2+,qb2+,qb+,Q+)
//			return &A1ph2q2QM2fm_16_eval<QMPOS,QPOS,QPOS2,QMPOS2>;
//			break;
//		case 0x0B070309:
//		case 0x050D0903:
//		case 0x0D050309:
//		case 0x070B0903:
//			return &A1ph2q2QM2fp_16_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x0D070309:
//		case 0x070D0903:
//		case 0x0B050309:
//		case 0x050B9003:
//			return &A1ph2q2QM2fm_16_eval<QPOS,QPOS2,QMPOS2,QMPOS>;
//			break;
//		case 0x090B0703:
//		case 0x03050D09:
//		case 0x090D0503:
//		case 0x03070B09:
//			return &A1ph2q2QM2fp_16_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x090D0703:
//		case 0x03070D09:
//		case 0x090B0503:
//		case 0x03050B09:
//			return &A1ph2q2QM2fm_16_eval<QPOS2,QMPOS2,QMPOS,QPOS>;
//			break;
//		case 0x03090B07:
//		case 0x0903050D:
//		case 0x03090D05:
//		case 0x0903070B:
//			return &A1ph2q2QM2fp_16_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
//		case 0x03090D07:
//		case 0x0903070D:
//		case 0x03090B05:
//		case 0x0903050B:
//			return &A1ph2q2QM2fm_16_eval<QMPOS2,QMPOS,QPOS,QPOS2>;
//			break;
			
		default:// We return zero for all other helicity combinations
//				cout << "5 pt A1ph2q2QM_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
			return 0; // Until everything is implemented return 0 for unknow results.
	}
}

template <class T> complex<T> (*A1ph2q2QM_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	// We first extract the position of the higgs and then analyse the colour ordered particles in the next level down
	long rem=hc;
	long newhc=0;
	int epos;
	long fac=0x1;
	for(int i=0;i<5;i++){
		int efac=rem%0x100;
		rem=rem/0x100;
		if(efac==0x0E||efac==0x0F||efac==0x10){ // As the higgs couples to a scalar only we treat all phis as the same
			epos=i;
		}
		else{
			newhc+=efac*fac;
			fac*=0x100;
		}
	}
		
	switch (epos) {
		case 0:
			return A1ph2q2QM_Tree_Ptr_non_higgs_eval<T,0,1,2,3>(newhc);
			break;
		case 1:
			return A1ph2q2QM_Tree_Ptr_non_higgs_eval<T,0,1,2,4>(newhc);
			break;
		case 2:
			return A1ph2q2QM_Tree_Ptr_non_higgs_eval<T,0,1,3,4>(newhc);
			break;
		case 3:
			return A1ph2q2QM_Tree_Ptr_non_higgs_eval<T,0,2,3,4>(newhc);
			break;
		case 4:
			return A1ph2q2QM_Tree_Ptr_non_higgs_eval<T,1,2,3,4>(newhc);
			break;
			
		default:// Should never get here as there is always a higgs. 
			return &ZeroF;
			break;
	}
}




/*
 *
 *
 * The 1 complex higgs, two massive quarks and a gluon amplitude
 *
 *
 */

//Q+phmQ-
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	_MESSAGE("Calling 1");
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> f4(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	Cmom<T> rf=*ep.ref();
	//Cmom<T> rf=f1;
	
	complex<T> fac1=((rf.Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())+(rf.Lt()*ep.p(i0)->Sm()*ep.p(i2)->L()))/T(2);
	complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=(ep.mom(i3)+ep.mom(i0)).square();
	
	return _OTHER_FAC*complex<T>(0,2)*(
							(ep.p(i2)->L()*f4.L())/(ep.p(i2)->Lt()*rf.Lt())
							*((f1.Lt()*rf.Lt())*(fac2)-(f1.Lt()*ep.p(i2)->Lt())*(fac1))
							-mass*((ep.ref()->Lt()*ep.p(i2)->Lt())*(ep.ref()->L()*ep.p(i2)->L())*(fac1)/((ep.p(i2)->Lt()*rf.Lt())*(f1.L()*ep.ref()->L())*(ep.ref()->Lt()*f4.Lt()))
								   // The commented out lines below are only needed when rf!=ep.ref
								   //-(ep.ref()->Lt()*rf.Lt())*(ep.ref()->L()*ep.p(i2)->L())*(fac2)/((ep.p(i2)->Lt()*rf.Lt())*(f1.L()*ep.ref()->L())*(ep.ref()->Lt()*f4.Lt()))
								)
							)/den;
	
}
//Q-phmQ+
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	_MESSAGE("Calling 2");
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> f4(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	Cmom<T> rf=*ep.ref();
	//		Cmom<T> rf=f1;
	
	complex<T> fac1=((rf.Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())+(rf.Lt()*ep.p(i0)->Sm()*ep.p(i2)->L()))/T(2);
	complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=(ep.mom(i3)+ep.mom(i0)).square();
	
	return _OTHER_FAC*complex<T>(0,-2)*(
							 (ep.p(i2)->L()*f1.L())/(ep.p(i2)->Lt()*rf.Lt())
							 *((f4.Lt()*rf.Lt())*(fac2)-(f4.Lt()*ep.p(i2)->Lt())*(fac1))
							 -mass*((ep.ref()->Lt()*ep.p(i2)->Lt())*(ep.ref()->L()*ep.p(i2)->L())*(fac1)/((ep.p(i2)->Lt()*rf.Lt())*(f4.L()*ep.ref()->L())*(ep.ref()->Lt()*f1.Lt()))
									// The commented out lines below are only needed when rf!=ep.ref
									//-(ep.ref()->Lt()*rf.Lt())*(ep.ref()->L()*ep.p(i2)->L())*(fac2)/((ep.p(i2)->Lt()*rf.Lt())*(f4.L()*ep.ref()->L())*(ep.ref()->Lt()*f1.Lt()))
								)
							 )/den;
	
}

//Q-phmQ-
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g7m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	_MESSAGE("Calling 3m");
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const Cmom<T> f1(ep.p(i0)->P()-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> f4(ep.p(i3)->P()-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> rf=*ep.ref();
	//	Cmom<T> rf=f4;
	
	momentum<complex<T> > prop=ep.mom(i3)+ep.mom(i0);
    complex<T> fac1=((rf.Lt()*f4.Lt())*(f4.L()*ep.p(i2)->L())+(rf.Lt()*f1.Lt())*(f1.L()*ep.p(i2)->L()))/T(2);
	//complex<T> fac2=S(i2,i3)+S(i2,i0);
 	
	complex<T> den=prop.square();
    
	return _OTHER_FAC*complex<T>(0,2)*(mass/den)*((f1.L()*ep.p(i2)->L())*(ep.p(i2)->Lt()*ep.ref()->Lt())/((ep.p(i2)->Lt()*rf.Lt())*(ep.ref()->Lt()*f4.Lt()))*fac1
									   +(ep.ref()->Lt()*ep.p(i2)->Lt())*(ep.p(i2)->L()*f4.L())*fac1/((rf.Lt()*ep.p(i2)->Lt())*(f1.Lt()*ep.ref()->Lt()))
									   // The commented out lines below are only needed when rf!=ep.ref
									   //-(f1.L()*ep.p(i2)->L())*(rf.Lt()*ep.ref()->Lt())/((ep.p(i2)->Lt()*rf.Lt())*(ep.ref()->Lt()*f4.Lt()))*fac2
									   //-(ep.ref()->Lt()*rf.Lt())*(ep.p(i2)->L()*f4.L())*fac2/((rf.Lt()*ep.p(i2)->Lt())*(f1.Lt()*ep.ref()->Lt()))
									   );
}
//Q-phmQ-
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g7p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	_MESSAGE("Calling 3p");
    
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const Cmom<T> f1(ep.p(i0)->P()-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> f4(ep.p(i3)->P()-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> rf=*ep.ref();
	//		Cmom<T> rf=f4;
	
	momentum<complex<T> > prop=ep.mom(i3)+ep.mom(i0);
	//complex<T> fac1=(rf.Lt()*smatrix<T>(prop)*ep.p(i2)->L())/T(2);
    complex<T> fac1=((rf.Lt()*f4.Lt())*(f4.L()*ep.p(i2)->L())+(rf.Lt()*f1.Lt())*(f1.L()*ep.p(i2)->L()))/T(2);
	//complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=prop.square();	
    
	return _OTHER_FAC*complex<T>(0,-2)*(mass/den)*((f1.L()*ep.p(i2)->L())*(ep.p(i2)->Lt()*ep.ref()->Lt())/((ep.p(i2)->Lt()*rf.Lt())*(ep.ref()->Lt()*f4.Lt()))*fac1
										+(ep.ref()->Lt()*ep.p(i2)->Lt())*(ep.p(i2)->L()*f4.L())*fac1/((rf.Lt()*ep.p(i2)->Lt())*(f1.Lt()*ep.ref()->Lt()))
										// The commented out lines below are only needed when rf!=ep.ref
									    //-(f1.L()*ep.p(i2)->L())*(rf.Lt()*ep.ref()->Lt())/((ep.p(i2)->Lt()*rf.Lt())*(ep.ref()->Lt()*f4.Lt()))*fac2
										//-(ep.ref()->Lt()*rf.Lt())*(ep.p(i2)->L()*f4.L())*fac2/((rf.Lt()*ep.p(i2)->Lt())*(f1.Lt()*ep.ref()->Lt()))
										);
}

//Q+phmQ+
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	_MESSAGE("Calling 4m");
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const Cmom<T> f1(ep.p(i0)->P()-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> f4(ep.p(i3)->P()-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> rf=*ep.ref();
	//		Cmom<T> rf=f4;
	
	complex<T> fac1=((rf.Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())+(rf.Lt()*ep.p(i0)->Sm()*ep.p(i2)->L()))/T(2);
	complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=(ep.mom(i3)+ep.mom(i0)).square();
	
	return _OTHER_FAC*complex<T>(0,2)*(mass/den)*((f1.Lt()*rf.Lt())*(ep.p(i2)->L()*ep.ref()->L())/((ep.p(i2)->Lt()*rf.Lt())*(ep.ref()->L()*f4.L()))*fac2
									   -(f1.Lt()*ep.p(i2)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/((ep.p(i2)->Lt()*rf.Lt())*(ep.ref()->L()*f4.L()))*fac1
									   +(ep.ref()->L()*ep.p(i2)->L())*(rf.Lt()*f4.Lt())*fac2/((rf.Lt()*ep.p(i2)->Lt())*(f1.L()*ep.ref()->L()))
									   -(ep.ref()->L()*ep.p(i2)->L())*(ep.p(i2)->Lt()*f4.Lt())*fac1/((rf.Lt()*ep.p(i2)->Lt())*(f1.L()*ep.ref()->L())));
}
//Q+phmQ+
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	_MESSAGE("Calling 4p");
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const Cmom<T> f1(ep.p(i0)->P()-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> f4(ep.p(i3)->P()-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> rf=*ep.ref();
	//		Cmom<T> rf=f4;
	
	complex<T> fac1=((rf.Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())+(rf.Lt()*ep.p(i0)->Sm()*ep.p(i2)->L()))/T(2);
	complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=(ep.mom(i3)+ep.mom(i0)).square();
	
	return _OTHER_FAC*complex<T>(0,-2)*(mass/den)*((f1.Lt()*rf.Lt())*(ep.p(i2)->L()*ep.ref()->L())/((ep.p(i2)->Lt()*rf.Lt())*(ep.ref()->L()*f4.L()))*fac2
										-(f1.Lt()*ep.p(i2)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/((ep.p(i2)->Lt()*rf.Lt())*(ep.ref()->L()*f4.L()))*fac1
										+(ep.ref()->L()*ep.p(i2)->L())*(rf.Lt()*f4.Lt())*fac2/((rf.Lt()*ep.p(i2)->Lt())*(f1.L()*ep.ref()->L()))
										-(ep.ref()->L()*ep.p(i2)->L())*(ep.p(i2)->Lt()*f4.Lt())*fac1/((rf.Lt()*ep.p(i2)->Lt())*(f1.L()*ep.ref()->L())));
}
	
	//phd
	
//Q-phdpQ+
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g14_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	_MESSAGE("Calling 5");
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
	
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> f4(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	Cmom<T> rf=*ep.ref();
	//		Cmom<T> rf=f1;
	
	complex<T> fac1=((rf.L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())+(rf.L()*ep.p(i0)->Sm()*ep.p(i2)->Lt()))/T(2);
	complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=(ep.mom(i3)+ep.mom(i0)).square();
	
	return _OTHER_FAC*complex<T>(0,2)*(
							(ep.p(i2)->Lt()*f4.Lt())/(ep.p(i2)->L()*rf.L())
							*((f1.L()*rf.L())*(fac2)-(f1.L()*ep.p(i2)->L())*(fac1))
							-mass*((ep.ref()->L()*ep.p(i2)->L())*(ep.ref()->Lt()*ep.p(i2)->Lt())*(fac1)/((ep.p(i2)->L()*rf.L())*(f1.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f4.L()))
								   // The commented out lines below are only needed when rf!=ep.ref
                                   //-(ep.ref()->L()*rf.L())*(ep.ref()->Lt()*ep.p(i2)->Lt())*(fac2)/((ep.p(i2)->L()*rf.L())*(f1.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f4.L()))
								   //
								   )
							)/den;
}
//Q+phdpQ-
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g13_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));
//	_MESSAGE("Calling 6");
	Cmom<T> f1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	Cmom<T> f4(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	
	Cmom<T> rf=*ep.ref();
	//		Cmom<T> rf=f1;
	
	complex<T> fac1=((rf.L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())+(rf.L()*ep.p(i0)->Sm()*ep.p(i2)->Lt()))/T(2);
	complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=(ep.mom(i3)+ep.mom(i0)).square();
	
	return _OTHER_FAC*complex<T>(0,-2)*(
							 (ep.p(i2)->Lt()*f1.Lt())/(ep.p(i2)->L()*rf.L())
							 *((f4.L()*rf.L())*(fac2)-(f4.L()*ep.p(i2)->L())*(fac1))
							 -mass*((ep.ref()->L()*ep.p(i2)->L())*(ep.ref()->Lt()*ep.p(i2)->Lt())*(fac1)/((ep.p(i2)->L()*rf.L())*(f4.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f1.L()))
									// The commented out lines below are only needed when rf!=ep.ref
                                    //-(ep.ref()->L()*rf.L())*(ep.ref()->Lt()*ep.p(i2)->Lt())*(fac2)/((ep.p(i2)->L()*rf.L())*(f4.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f1.L()))
									//
									)
							 )/den;
	
}

//Q+phdpQ+
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g16m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    _MESSAGE("Calling 7m");
    
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const Cmom<T> f1(ep.p(i0)->P()-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> f4(ep.p(i3)->P()-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> rf=*ep.ref();
	//	Cmom<T> rf=f4;

	momentum<complex<T> > prop=ep.mom(i3)+ep.mom(i0);
	//complex<T> fac2=S(i2,i3)+S(i2,i0);
    complex<T> fac1=((rf.L()*f4.L())*(f4.Lt()*ep.p(i2)->Lt())+(rf.L()*f1.L())*(f1.Lt()*ep.p(i2)->Lt()))/T(2);
 	
	complex<T> den=prop.square();
	
	return _OTHER_FAC*complex<T>(0,2)*(mass/den)*((f1.Lt()*ep.p(i2)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/((ep.p(i2)->L()*rf.L())*(ep.ref()->L()*f4.L()))*fac1
									   +(ep.ref()->L()*ep.p(i2)->L())*(ep.p(i2)->Lt()*f4.Lt())*fac1/((rf.L()*ep.p(i2)->L())*(f1.L()*ep.ref()->L()))
                                        // The commented out lines below are only needed when rf!=ep.ref
                                        //-(f1.Lt()*ep.p(i2)->Lt())*(rf.L()*ep.ref()->L())/((ep.p(i2)->L()*rf.L())*(ep.ref()->L()*f4.L()))*fac2
									   //-(ep.ref()->L()*rf.L())*(ep.p(i2)->Lt()*f4.Lt())*fac2/((rf.L()*ep.p(i2)->L())*(f1.L()*ep.ref()->L()))
									   //
									   );
}
//Q+phdpQ+
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g16p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    _MESSAGE("Calling 7p");
     
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const Cmom<T> f1(ep.p(i0)->P()-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> f4(ep.p(i3)->P()-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> rf=*ep.ref();
	//	Cmom<T> rf=f4;
	
	momentum<complex<T> > prop=ep.mom(i3)+ep.mom(i0);
    complex<T> fac1=((rf.L()*f4.L())*(f4.Lt()*ep.p(i2)->Lt())+(rf.L()*f1.L())*(f1.Lt()*ep.p(i2)->Lt()))/T(2);
 	
	complex<T> den=prop.square();
	
	return _OTHER_FAC*complex<T>(0,-2)*(mass/den)*((f1.Lt()*ep.p(i2)->Lt())*(ep.p(i2)->L()*ep.ref()->L())/((ep.p(i2)->L()*rf.L())*(ep.ref()->L()*f4.L()))*fac1
									   +(ep.ref()->L()*ep.p(i2)->L())*(ep.p(i2)->Lt()*f4.Lt())*fac1/((rf.L()*ep.p(i2)->L())*(f1.L()*ep.ref()->L()))
                                        // The commented out lines below are only needed when rf!=ep.ref
                                        //-(f1.Lt()*ep.p(i2)->Lt())*(rf.L()*ep.ref()->L())/((ep.p(i2)->L()*rf.L())*(ep.ref()->L()*f4.L()))*fac2
									   //-(ep.ref()->L()*rf.L())*(ep.p(i2)->Lt()*f4.Lt())*fac2/((rf.L()*ep.p(i2)->L())*(f1.L()*ep.ref()->L()))
										//
										);
}

//Q-phdpQ-
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g15m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    _MESSAGE("Calling 8m");
    
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const Cmom<T> f1(ep.p(i0)->P()-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> f4(ep.p(i3)->P()-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> rf=*ep.ref();
	//		Cmom<T> rf=f4;
	
	complex<T> fac1=((rf.L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())+(rf.L()*ep.p(i0)->Sm()*ep.p(i2)->Lt()))/T(2);
	complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=(ep.mom(i3)+ep.mom(i0)).square();
	
	return _OTHER_FAC*complex<T>(0,2)*(mass/den)*((f1.L()*rf.L())*(ep.p(i2)->Lt()*ep.ref()->Lt())/((ep.p(i2)->L()*rf.L())*(ep.ref()->Lt()*f4.Lt()))*fac2
									   -(f1.L()*ep.p(i2)->L())*(ep.p(i2)->Lt()*ep.ref()->Lt())/((ep.p(i2)->L()*rf.L())*(ep.ref()->Lt()*f4.Lt()))*fac1
									   +(ep.ref()->Lt()*ep.p(i2)->Lt())*(rf.L()*f4.L())*fac2/((rf.L()*ep.p(i2)->L())*(f1.Lt()*ep.ref()->Lt()))
									   -(ep.ref()->Lt()*ep.p(i2)->Lt())*(ep.p(i2)->L()*f4.L())*fac1/((rf.L()*ep.p(i2)->L())*(f1.Lt()*ep.ref()->Lt())));
}
//Q-phdpQ-
template <int i0, int i2, int i3, class T> complex<T>  A1ph2QM1g15p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    _MESSAGE("Calling 8p");
    
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const Cmom<T> f1(ep.p(i0)->P()-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> f4(ep.p(i3)->P()-(mass2/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
	Cmom<T> rf=*ep.ref();
	//		Cmom<T> rf=f4;
	
	complex<T> fac1=((rf.L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())+(rf.L()*ep.p(i0)->Sm()*ep.p(i2)->Lt()))/T(2);
	complex<T> fac2=S(i2,i3)+S(i2,i0);
	
	complex<T> den=(ep.mom(i3)+ep.mom(i0)).square();
	
	return _OTHER_FAC*complex<T>(0,-2)*(mass/den)*((f1.L()*rf.L())*(ep.p(i2)->Lt()*ep.ref()->Lt())/((ep.p(i2)->L()*rf.L())*(ep.ref()->Lt()*f4.Lt()))*fac2
									   -(f1.L()*ep.p(i2)->L())*(ep.p(i2)->Lt()*ep.ref()->Lt())/((ep.p(i2)->L()*rf.L())*(ep.ref()->Lt()*f4.Lt()))*fac1
									   +(ep.ref()->Lt()*ep.p(i2)->Lt())*(rf.L()*f4.L())*fac2/((rf.L()*ep.p(i2)->L())*(f1.Lt()*ep.ref()->Lt()))
									   -(ep.ref()->Lt()*ep.p(i2)->Lt())*(ep.p(i2)->L()*f4.L())*fac1/((rf.L()*ep.p(i2)->L())*(f1.Lt()*ep.ref()->Lt())));
}

template <class T, int QPOS, int GLUE, int QPOS2> complex<T> (*A1ph2QM1g_phi_simp_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	//	cout << hex << hc << dec << endl;
	
	switch (hc) {
			// The ph amplitudes
		case 0x607: //167: //(phd,Qb+,g-,Q-)
		case 0x805:
			return &A1ph2QM1g5_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x760: //142: //(Q-,Qb+,g-)
		case 0x580:
			return &A1ph2QM1g5_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x076: //207: //(g-,Q-,Qb+)
		case 0x058:
			return &A1ph2QM1g5_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x508: //202: //(phd,Qb-,g-,Q+)
		case 0x706:
			return &A1ph2QM1g6_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x850: //137: //(Q+,Qb-,g-)
		case 0x670:
			return &A1ph2QM1g6_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x085: //177: //(g-,Q+,Qb-)
		case 0x067:
			return &A1ph2QM1g6_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x705: //166: //(phd,Qb-,g-,Q-)
			return &A1ph2QM1g7m_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x507:
			return &A1ph2QM1g7p_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x570: //136: //(Q-,Qb-,g-)
			return &A1ph2QM1g7m_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x750:
			return &A1ph2QM1g7p_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x057: //171: //(g-,Q-,Qb-)
			return &A1ph2QM1g7m_eval<QPOS2,QPOS,GLUE>;
			break;
		case 0x075:
			return &A1ph2QM1g7p_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x806: //166: //(phd,Qb+,g-,Q+)
			return &A1ph2QM1g8m_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x608:
			return &A1ph2QM1g8p_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x680:
			return &A1ph2QM1g8m_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x860:
			return &A1ph2QM1g8p_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x068:
			return &A1ph2QM1g8m_eval<QPOS2,QPOS,GLUE>;
			break;
		case 0x086:
			return &A1ph2QM1g8p_eval<QPOS2,QPOS,GLUE>;
			break;
			
		default:// We return zero for all other helicity combinations
//			cout << "Missing helcode ph: " << hex << hc << dec << endl;
			return &ZeroF;
	}
}
	
template <class T, int QPOS, int GLUE, int QPOS2> complex<T> (*A1ph2QM1g_phi_simp_dag_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	//	cout << hex << hc << dec << endl;
	
	switch (hc) {
			// The ph amplitudes
		case 0x617: //167: //(phd,Qb+,g+,Q-)
		case 0x815:
			return &A1ph2QM1g13_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x761: //142: //(Q-,Qb+,g+)
		case 0x581:
			return &A1ph2QM1g13_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x176: //207: //(g+,Q-,Qb+)
		case 0x158:
			return &A1ph2QM1g13_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x518: //202: //(phd,Qb-,g+,Q+)
		case 0x716:
			return &A1ph2QM1g14_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x851: //137: //(Q+,Qb-,g+)
		case 0x671:
			return &A1ph2QM1g14_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x185: //177: //(g+,Q+,Qb-)
		case 0x167:
			return &A1ph2QM1g14_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x715: //166: //(phd,Qb-,g+,Q-)
			return &A1ph2QM1g15m_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x517:
			return &A1ph2QM1g15p_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x571: //136: //(Q-,Qb-,g+)
			return &A1ph2QM1g15m_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x751:
			return &A1ph2QM1g15p_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x157: //171: //(g+,Q-,Qb-)
			return &A1ph2QM1g15m_eval<QPOS2,QPOS,GLUE>;
			break;
		case 0x175:
			return &A1ph2QM1g15p_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x816: //166: //(phd,Qb+,g+,Q+)
			return &A1ph2QM1g16m_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x618:
			return &A1ph2QM1g16p_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x681:
			return &A1ph2QM1g16m_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x861:
			return &A1ph2QM1g16p_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x168:
			return &A1ph2QM1g16m_eval<QPOS2,QPOS,GLUE>;
			break;
		case 0x186:
			return &A1ph2QM1g16p_eval<QPOS2,QPOS,GLUE>;
			break;
		default:// We return zero for all other helicity combinations
//			cout << "Missing helcode phd: " << hex << hc << dec << endl;
			return &ZeroF;
	}
}
			
			
template <class T, int QPOS, int GLUE, int QPOS2> complex<T> (*A1ph2QM1g_phi_simp_full_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	//	cout << hex << hc << dec << endl;
	
	switch (hc) {
			// The H amplitudes
		case 0x607: //167: //(phd,Qb+,g-,Q-)
		case 0x805:
			return &A1ph2QM1g5_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x760: //142: //(Q-,Qb+,g-)
		case 0x580:
			return &A1ph2QM1g5_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x076: //207: //(g-,Q-,Qb+)
		case 0x058:
			return &A1ph2QM1g5_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x508: //202: //(phd,Qb-,g-,Q+)
		case 0x706:
			return &A1ph2QM1g6_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x850: //137: //(Q+,Qb-,g-)
		case 0x670:
			return &A1ph2QM1g6_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x085: //177: //(g-,Q+,Qb-)
		case 0x067:
			return &A1ph2QM1g6_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x705: //166: //(phd,Qb-,g-,Q-)
			return &A1ph2QM1g7m_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x507:
			return &A1ph2QM1g7p_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x570: //136: //(Q-,Qb-,g-)
			return &A1ph2QM1g7m_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x750:
			return &A1ph2QM1g7p_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x057: //171: //(g-,Q-,Qb-)
			return &A1ph2QM1g7m_eval<QPOS2,QPOS,GLUE>;
			break;
		case 0x075:
			return &A1ph2QM1g7p_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x806: //166: //(phd,Qb+,g-,Q+)
			return &A1ph2QM1g8m_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x608:
			return &A1ph2QM1g8p_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x680:
			return &A1ph2QM1g8m_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x860:
			return &A1ph2QM1g8p_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x068:
			return &A1ph2QM1g8m_eval<QPOS2,QPOS,GLUE>;
			break;
		case 0x086:
			return &A1ph2QM1g8p_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x617: //167: //(phd,Qb+,g+,Q-)
		case 0x815:
			return &A1ph2QM1g13_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x761: //142: //(Q-,Qb+,g+)
		case 0x581:
			return &A1ph2QM1g13_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x176: //207: //(g+,Q-,Qb+)
		case 0x158:
			return &A1ph2QM1g13_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x518: //202: //(phd,Qb-,g+,Q+)
		case 0x716:
			return &A1ph2QM1g14_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x851: //137: //(Q+,Qb-,g+)
		case 0x671:
			return &A1ph2QM1g14_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x185: //177: //(g+,Q+,Qb-)
		case 0x167:
			return &A1ph2QM1g14_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x715: //166: //(phd,Qb-,g+,Q-)
			return &A1ph2QM1g15m_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x517:
			return &A1ph2QM1g15p_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x571: //136: //(Q-,Qb-,g+)
			return &A1ph2QM1g15m_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x751:
			return &A1ph2QM1g15p_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x157: //171: //(g+,Q-,Qb-)
			return &A1ph2QM1g15m_eval<QPOS2,QPOS,GLUE>;
			break;
		case 0x175:
			return &A1ph2QM1g15p_eval<QPOS2,QPOS,GLUE>;
			break;
			
			
		case 0x816: //166: //(phd,Qb+,g+,Q+)
			return &A1ph2QM1g16m_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x618:
			return &A1ph2QM1g16p_eval<QPOS,GLUE,QPOS2>;
			break;
		case 0x681:
			return &A1ph2QM1g16m_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x861:
			return &A1ph2QM1g16p_eval<GLUE,QPOS2,QPOS>;
			break;
		case 0x168:
			return &A1ph2QM1g16m_eval<QPOS2,QPOS,GLUE>;
			break;
		case 0x186:
			return &A1ph2QM1g16p_eval<QPOS2,QPOS,GLUE>;
			break;
		default:// We return zero for all other helicity combinations
//			cout << "Missing helcode H: " << hex << hc << dec << endl;
			return &ZeroF;
	}
}
	
template <class T> complex<T> (*A1ph2QM1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	// We first extract the position of the higgs and then analyse the colour ordered particles in the next level down
	long rem=hc;
	long newhc=0;
	int epos;
	long fac=0x1;
	bool dagger=false;
	bool full_higgs=false;
	for(int i=0;i<4;i++){
		int efac=rem%0x10;
		rem=rem/0x10;
		if(efac==0x9||efac==0xA||efac==0xE){
			epos=i;
			if(efac==0xA){// If this is a phi dagger then record this
				dagger=true;
			}
			if(efac==0xE){// If this is the full higgs then record this
				full_higgs=true;
			}
		}
		else{
			newhc+=efac*fac;
			fac*=0x10;
		}
	}
	
	// We now have the position of the higgs (from the right hand side, i.e. 9XXX=>3) in epos and a new helcode minus the higgs in newhc
	
	if(!full_higgs){
		if(!dagger){
			switch (epos) {
				case 0:
					return A1ph2QM1g_phi_simp_Tree_Ptr_eval<T,0,1,2>(newhc);
					break;
				case 1:
					return A1ph2QM1g_phi_simp_Tree_Ptr_eval<T,0,1,3>(newhc);
					break;
				case 2:
					return A1ph2QM1g_phi_simp_Tree_Ptr_eval<T,0,2,3>(newhc);
					break;
				case 3:
					return A1ph2QM1g_phi_simp_Tree_Ptr_eval<T,1,2,3>(newhc);
					break;
					
				default:// Should never get here as there is always a higgs. 
					return &ZeroF;
					break;
			}
		}
		else{
			switch (epos) {
				case 0:
					return A1ph2QM1g_phi_simp_dag_Tree_Ptr_eval<T,0,1,2>(newhc);
					break;
				case 1:
					return A1ph2QM1g_phi_simp_dag_Tree_Ptr_eval<T,0,1,3>(newhc);
					break;
				case 2:
					return A1ph2QM1g_phi_simp_dag_Tree_Ptr_eval<T,0,2,3>(newhc);
					break;
				case 3:
					return A1ph2QM1g_phi_simp_dag_Tree_Ptr_eval<T,1,2,3>(newhc);
					break;
					
				default:// Should never get here as there is always a higgs. 
					return &ZeroF;
					break;
			}
		}
	}
	else{
		switch (epos) {
			case 0:
				return A1ph2QM1g_phi_simp_full_Tree_Ptr_eval<T,0,1,2>(newhc);
				break;
			case 1:
				return A1ph2QM1g_phi_simp_full_Tree_Ptr_eval<T,0,1,3>(newhc);
				break;
			case 2:
				return A1ph2QM1g_phi_simp_full_Tree_Ptr_eval<T,0,2,3>(newhc);
				break;
			case 3:
				return A1ph2QM1g_phi_simp_full_Tree_Ptr_eval<T,1,2,3>(newhc);
				break;
				
			default:// Should never get here as there is always a higgs. 
				return &ZeroF;
				break;
		}
	}

}


template complex<R> (*A1ph2q2QM_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2q2QM_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2q2QM_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1ph2QM1g_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2QM1g_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2QM1g_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> (*A1ph2q2QM_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2QM1g_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);

#endif

}


