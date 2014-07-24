/* The 2 gluon massive +/- amplitudes */

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


    

// Defined in A0_1phi_scalar.cpp
momentum<complex<R> > ddim2_R(complex<R>(-1,0),complex<R>(0,0),complex<R>(0,0),complex<R>(0,0));
momentum<complex<RHP> > ddim2_RHP(complex<RHP>(1,0),complex<RHP>(0,0),complex<RHP>(0,0),complex<RHP>(0,0));
momentum<complex<RVHP> > ddim2_RVHP(complex<RVHP>(1,0),complex<RVHP>(0,0),complex<RVHP>(0,0),complex<RVHP>(0,0));
template <typename T> momentum<complex<T> >& get_ddim2(){};
template<> momentum<complex<R> >& get_ddim2(){return ddim2_R;};
template<> momentum<complex<RHP> >& get_ddim2(){return ddim2_RHP;};
template<> momentum<complex<RVHP> >& get_ddim2(){return ddim2_RVHP;};
#ifdef BH_USE_GMP
momentum<complex<RGMP> > ddim2_RGMP(complex<RGMP>(1,0),complex<RGMP>(0,0),complex<RGMP>(0,0),complex<RGMP>(0,0));
template<> momentum<complex<RGMP> >& get_ddim2(){return ddim2_RGMP;};
#endif


/*
 *
 *
 * The 2 gluon massive +/- and a single gluon amplitudes
 *
 *
 */
    
// G+ g+ G-
template <int i0, int i1, int i2, class T> complex<T>  A2Gsc1g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
        
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
    const Cmom<T> p3L(ep.mom(i2)-(mass2/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());

	momentum<complex<T> > eps1=-PfLLt(p1L.Lt(),ep.ref()->L())/(ep.ref()->L()*p1L.L());
	momentum<complex<T> > eps2=-PfLLt(ep.p(i1)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i1)->L());
	momentum<complex<T> > eps3=PfLLt(ep.ref()->Lt(),p3L.L())/(ep.ref()->Lt()*p3L.Lt());
	
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=zero;
	
	
	//SgS
	complex<T> res1(complex<T>(0,-2)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));

	
//	_MESSAGE4(" glue A2Gsc1g1 : ",res1,"==",complex<T>(0,-1)*pow(p1L.Lt()*ep.p(i1)->Lt(),3)/((ep.p(i1)->Lt()*p3L.Lt())*(p3L.Lt()*p1L.Lt())));
    
//     return complex<T>(0,-1)*pow(p1L.Lt()*ep.p(i1)->Lt(),3)/((ep.p(i1)->Lt()*p3L.Lt())*(p3L.Lt()*p1L.Lt()));
	return res1;
}
// G+ g- G-
template <int i0, int i1, int i2, class T> complex<T>  A2Gsc1g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){    
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);
//    
//    const Cmom<T> nref=*ep.ref();
//	const Cmom<T> p1L(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
//    const Cmom<T> p3L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
//
//	_MESSAGE2(" glue A2Gsc1g2 : ",complex<T>(0,1)*pow(ep.p(i1)->L()*p3L.L(),3)/((p1L.L()*ep.p(i1)->L())*(p3L.L()*p1L.L())));
//	return complex<T>(0,1)*pow(ep.p(i1)->L()*p3L.L(),3)/((p1L.L()*ep.p(i1)->L())*(p3L.L()*p1L.L()));
	
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
    const Cmom<T> p3L(ep.mom(i2)-(mass2/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=-PfLLt(p1L.Lt(),ep.ref()->L())/(ep.ref()->L()*p1L.L());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i1)->L())/(ep.ref()->Lt()*ep.p(i1)->Lt());
	momentum<complex<T> > eps3=PfLLt(ep.ref()->Lt(),p3L.L())/(ep.ref()->Lt()*p3L.Lt());
	
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=zero;
	
	
	//SgS
	complex<T> res1(complex<T>(0,-2)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A2Gsc1g2 : ",res1,"==",complex<T>(0,1)*pow(ep.p(i1)->L()*p3L.L(),3)/((p1L.L()*ep.p(i1)->L())*(p3L.L()*p1L.L())));
    
	return res1;

}
    
// G- g+ G+
template <int i0, int i1, int i2, class T> complex<T>  A2Gsc1g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);
// 
//    const Cmom<T> nref=*ep.ref();
//	const Cmom<T> p1L(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
//    const Cmom<T> p3L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
//    
//		_MESSAGE2(" glue A2Gsc1g3 : ",complex<T>(0,-1)*pow(ep.p(i1)->Lt()*p3L.Lt(),3)/((p1L.Lt()*ep.p(i1)->Lt())*(p3L.Lt()*p1L.Lt())));
//    return complex<T>(0,-1)*pow(ep.p(i1)->Lt()*p3L.Lt(),3)/((p1L.Lt()*ep.p(i1)->Lt())*(p3L.Lt()*p1L.Lt()));
	
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
    const Cmom<T> p3L(ep.mom(i2)-(mass2/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=PfLLt(ep.ref()->Lt(),p1L.L())/(ep.ref()->Lt()*p1L.Lt());
	momentum<complex<T> > eps2=-PfLLt(ep.p(i1)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i1)->L());
	momentum<complex<T> > eps3=-PfLLt(p3L.Lt(),ep.ref()->L())/(ep.ref()->L()*p3L.L());
	
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=zero;
	
	
	//SgS
	complex<T> res1(complex<T>(0,-2)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A2Gsc1g3 : ",res1,"==",complex<T>(0,-1)*pow(ep.p(i1)->Lt()*p3L.Lt(),3)/((p1L.Lt()*ep.p(i1)->Lt())*(p3L.Lt()*p1L.Lt())));
    
	return res1;

}
// G- g- G+
template <int i0, int i1, int i2, class T> complex<T>  A2Gsc1g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);
//
//    const Cmom<T> nref=*ep.ref();
//	const Cmom<T> p1L(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
//    const Cmom<T> p3L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
//    
//		_MESSAGE2(" glue A2Gsc1g4 : ",complex<T>(0,1)*pow(p1L.L()*ep.p(i1)->L(),3)/((ep.p(i1)->L()*p3L.L())*(p3L.L()*p1L.L())));
//    return complex<T>(0,1)*pow(p1L.L()*ep.p(i1)->L(),3)/((ep.p(i1)->L()*p3L.L())*(p3L.L()*p1L.L()));
	
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
    const Cmom<T> p3L(ep.mom(i2)-(mass2/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=PfLLt(ep.ref()->Lt(),p1L.L())/(ep.ref()->Lt()*p1L.Lt());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i1)->L())/(ep.ref()->Lt()*ep.p(i1)->Lt());
	momentum<complex<T> > eps3=-PfLLt(p3L.Lt(),ep.ref()->L())/(ep.ref()->L()*p3L.L());
	
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=zero;
	
	
	//SgS
	complex<T> res1(complex<T>(0,-2)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A2Gsc1g4 : ",res1,"==",complex<T>(0,1)*pow(p1L.L()*ep.p(i1)->L(),3)/((ep.p(i1)->L()*p3L.L())*(p3L.L()*p1L.L())));
    
	return res1;
}
    
// G- g+ G-
template <int i0, int i1, int i2, class T> complex<T>  A2Gsc1g5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);
//
//    const Cmom<T> nref=*ep.ref();
//	const Cmom<T> p1L(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
//    const Cmom<T> p3L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
//    
//		_MESSAGE2(" glue A2Gsc1g5 : ",complex<T>(0,-1)*pow(p1L.L()*p3L.L(),3)/((p1L.L()*ep.p(i1)->L())*(ep.p(i1)->L()*p3L.L())));
//    return complex<T>(0,-1)*pow(p1L.L()*p3L.L(),3)/((p1L.L()*ep.p(i1)->L())*(ep.p(i1)->L()*p3L.L()));
	
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
    const Cmom<T> p3L(ep.mom(i2)-(mass2/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=PfLLt(ep.ref()->Lt(),p1L.L())/(ep.ref()->Lt()*p1L.Lt());
	momentum<complex<T> > eps2=-PfLLt(ep.p(i1)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i1)->L());
	momentum<complex<T> > eps3=PfLLt(ep.ref()->Lt(),p3L.L())/(ep.ref()->Lt()*p3L.Lt());
	
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=zero;
	
	
	//SgS
	complex<T> res1(complex<T>(0,-2)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A2Gsc1g5 : ",res1,"==",complex<T>(0,-1)*pow(p1L.L()*p3L.L(),3)/((p1L.L()*ep.p(i1)->L())*(ep.p(i1)->L()*p3L.L())));
    
	return res1;
}
// G+ g- G+
template <int i0, int i1, int i2, class T> complex<T>  A2Gsc1g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i0))/T(2);
//    
//	const Cmom<T> nref=*ep.ref();
//	const Cmom<T> p1L(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
//    const Cmom<T> p3L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());	
//    
//		_MESSAGE2(" glue A2Gsc1g6 : ",complex<T>(0,1)*pow(p1L.Lt()*p3L.Lt(),3)/((p1L.Lt()*ep.p(i1)->Lt())*(ep.p(i1)->Lt()*p3L.Lt())));
//    return complex<T>(0,1)*pow(p1L.Lt()*p3L.Lt(),3)/((p1L.Lt()*ep.p(i1)->Lt())*(ep.p(i1)->Lt()*p3L.Lt()));
	
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i0));
	const complex<T> mass=eval_param<T>::mass(masses.p(i0));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
    const Cmom<T> p3L(ep.mom(i2)-(mass2/(T(2)*(ep.p(i2)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=-PfLLt(p1L.Lt(),ep.ref()->L())/(ep.ref()->L()*p1L.L());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i1)->L())/(ep.ref()->Lt()*ep.p(i1)->Lt());
	momentum<complex<T> > eps3=-PfLLt(p3L.Lt(),ep.ref()->L())/(ep.ref()->L()*p3L.L());
	
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=zero;
	
	
	//SgS
	complex<T> res1(complex<T>(0,-2)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A2Gsc1g6 : ",res1,"==",complex<T>(0,1)*pow(p1L.Lt()*p3L.Lt(),3)/((p1L.Lt()*ep.p(i1)->Lt())*(ep.p(i1)->Lt()*p3L.Lt())));
    
	return res1;

}
    
    
/*
 *
 *
 * The 1 gluon massive +/- 1 gluon scalar and a single gluon amplitudes
 *
 *
 */

//G-g+Gsc
template <int i0, int i1, int i2, class T> complex<T>  A1Gsc1gM1g1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
    
    const Cmom<T> pGm(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
//    const Cmom<T> pGs(ep.mom(i2)-(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
//    
//	
//	_MESSAGE2(" glue A1Gsc1gM1g1p : ",complex<T>(0,1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt())/((pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L()))
//			  *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P())));
//    return complex<T>(0,1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt())/((pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L()))
//    *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P()));
//    
//    //return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*pow(pGs.Lt()*ep.p(i1)->Lt(),2)/((pGm.Lt()*ep.p(i1)->Lt())*(pGs.Lt()*pGm.Lt()))*((pGm.L()*pGs.L())/(pGm.L()*ep.p(i1)->L()));
	
	
//	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=PfLLt(ep.ref()->Lt(),p1L.L())/(ep.ref()->Lt()*p1L.Lt());
	momentum<complex<T> > eps2=-PfLLt(ep.p(i1)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i1)->L());
	momentum<complex<T> > eps3=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
//	momentum<complex<T> > eps3=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=-ddim;
	
	
	//SgS
	complex<T> res1(sqrt(T(2))*complex<T>(0,1)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A1Gsc1gM1g1p : ",res1,"==",complex<T>(0,1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*
//		(pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt())/((pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L()))
//			*((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P())));
    
	return res1;

}
template <int i0, int i1, int i2, class T> complex<T>  A1Gsc1gM1g1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
    
    const Cmom<T> pGm(ep.mom(i1)-(mass2/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
//    const Cmom<T> pGs(ep.mom(i2)-(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
//    
//	_MESSAGE2(" glue A1Gsc1gM1g1m : ",complex<T>(0,-1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt())/((pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L()))
//			  *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P())));
//    return complex<T>(0,-1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt())/((pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L()))
//    *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P()));
//    
//    //return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*pow(pGs.Lt()*ep.p(i1)->Lt(),2)/((pGm.Lt()*ep.p(i1)->Lt())*(pGs.Lt()*pGm.Lt()))*((pGm.L()*pGs.L())/(pGm.L()*ep.p(i1)->L()));
	
//	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i1)-(mass2/(T(2)*(ep.p(i1)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=-PfLLt(ep.p(i0)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i0)->L());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),p1L.L())/(ep.ref()->Lt()*p1L.Lt());
	momentum<complex<T> > eps3=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
//	momentum<complex<T> > eps3=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > mom1D(zero);
	momentum<complex<T> > mom2D(mass*ddim);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=-ddim;
	
	//SgS
	complex<T> res1(sqrt(T(2))*complex<T>(0,1)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A1Gsc1gM1g1m : ",res1,"==",complex<T>(0,1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i1))*(pGm.L()*ep.ref()->L())*(ep.p(i0)->Lt()*ep.ref()->Lt())/((pGm.Lt()*ep.ref()->Lt())*(ep.p(i0)->L()*ep.ref()->L()))
//			  *((ep.ref()->P()*ep.p(i0)->P())/(ep.ref()->P()*ep.p(i2)->P())));
    
	return res1;
}
    
//G+g-Gsc
template <int i0, int i1, int i2, class T> complex<T>  A1Gsc1gM1g2p_eval(const eval_param<T>& ep, const mass_param_coll& masses){    
    const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
    
    const Cmom<T> pGm(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
//    const Cmom<T> pGs(ep.mom(i2)-(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
//        
//	_MESSAGE2(" glue A1Gsc1gM1g2p : ",complex<T>(0,1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L())/((pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt()))
//			  *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P())));
//    return complex<T>(0,1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L())/((pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt()))
//    *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P()));
//    
//    //return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*pow(pGs.Lt()*ep.p(i1)->Lt(),2)/((pGm.Lt()*ep.p(i1)->Lt())*(pGs.Lt()*pGm.Lt()))*((pGm.L()*pGs.L())/(pGm.L()*ep.p(i1)->L()));
	
//	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i0)-(mass2/(T(2)*(ep.p(i0)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=-PfLLt(p1L.Lt(),ep.ref()->L())/(ep.ref()->L()*p1L.L());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i1)->L())/(ep.ref()->Lt()*ep.p(i1)->Lt());
	momentum<complex<T> > eps3=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
//	momentum<complex<T> > eps3=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=-ddim;
	
	
	//SgS
	complex<T> res1(sqrt(T(2))*complex<T>(0,1)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A1Gsc1gM1g2p : ",res1,"==",complex<T>(0,1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L())/((pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt()))
//			  *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P())));
    
	return res1;

}
template <int i0, int i1, int i2, class T> complex<T>  A1Gsc1gM1g2m_eval(const eval_param<T>& ep, const mass_param_coll& masses){    
    const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
    
    const Cmom<T> pGm(ep.mom(i1)-(mass2/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
//    const Cmom<T> pGs(ep.mom(i2)-(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
//    
//	_MESSAGE2(" glue A1Gsc1gM1g2m : ",complex<T>(0,-1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L())/((pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt()))
//			  *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P())));
//    return complex<T>(0,-1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i0))*(pGm.Lt()*ep.ref()->Lt())*(ep.p(i1)->L()*ep.ref()->L())/((pGm.L()*ep.ref()->L())*(ep.p(i1)->Lt()*ep.ref()->Lt()))
//    *((ep.ref()->P()*ep.p(i1)->P())/(ep.ref()->P()*ep.p(i2)->P()));
//    
//    //return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*pow(pGs.Lt()*ep.p(i1)->Lt(),2)/((pGm.Lt()*ep.p(i1)->Lt())*(pGs.Lt()*pGm.Lt()))*((pGm.L()*pGs.L())/(pGm.L()*ep.p(i1)->L()));
	
//	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	
    const Cmom<T> nref=*ep.ref();
	const Cmom<T> p1L(ep.mom(i1)-(mass2/(T(2)*(ep.p(i1)->P()*nref.P())))*nref.P());
	
	momentum<complex<T> > eps1=PfLLt(ep.ref()->Lt(),ep.p(i0)->L())/(ep.ref()->Lt()*ep.p(i0)->Lt());
	momentum<complex<T> > eps2=-PfLLt(p1L.Lt(),ep.ref()->L())/(ep.ref()->L()*p1L.L());
	momentum<complex<T> > eps3=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
//	momentum<complex<T> > eps3=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim2<T>());
//	momentum<complex<T> > ddim(zero);
	momentum<complex<T> > mom1D(zero);
	momentum<complex<T> > mom2D(mass*ddim);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=-ddim;	
	
	//SgS
	complex<T> res1(sqrt(T(2))*complex<T>(0,1)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
									 +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
									 +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
	
//	_MESSAGE4(" glue A1Gsc1gM1g2m : ",res1,"==",complex<T>(0,1)*sqrt(T(2))*eval_param<T>::mass(masses.p(i1))*(pGm.Lt()*ep.ref()->Lt())*(ep.p(i0)->L()*ep.ref()->L())/((pGm.L()*ep.ref()->L())*(ep.p(i0)->Lt()*ep.ref()->Lt()))
//			*((ep.ref()->P()*ep.p(i0)->P())/(ep.ref()->P()*ep.p(i2)->P())));
    
	return res1;

}

    
/*
 *
 *
 * The 2 scalar and a single gluon amplitudes
 *
 *
 */

template <int i0, int i1, int i2, class T> complex<T>  A2Gsc1g7_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
//	momentum<complex<T> > eps1=(ep.mom(i0)/mass)-(mass*ep.ref()->P())/(ep.mom(i0)*ep.ref()->P());
	momentum<complex<T> > eps1=(mass*ep.ref()->P())/(ep.mom(i0)*ep.ref()->P());
	momentum<complex<T> > eps2=-PfLLt(ep.p(i1)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i1)->L());
//	momentum<complex<T> > eps3=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	momentum<complex<T> > eps3=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=ddim;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=-ddim;

	
	//SgS
	complex<T> res1(complex<T>(0,-1)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
					+((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
					+((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
//    _MESSAGE2(" glue A1Gsc1gM1g7 : ",res1);

	
	return res1;
//    return complex<T>(0,-1)*(ep.ref()->L()*ep.p(i0)->Sm()*ep.p(i1)->Lt())/(ep.ref()->L()*ep.p(i1)->L());
}
template <int i0, int i1, int i2, class T> complex<T>  A2Gsc1g8_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
//	momentum<complex<T> > eps1=(ep.mom(i0)/mass)-(mass*ep.ref()->P())/(ep.mom(i0)*ep.ref()->P());
	momentum<complex<T> > eps1=(mass*ep.ref()->P())/(ep.mom(i0)*ep.ref()->P());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i1)->L())/(ep.ref()->Lt()*ep.p(i1)->Lt());
//	momentum<complex<T> > eps3=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	momentum<complex<T> > eps3=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());

	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim2<T>());
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(zero);
	momentum<complex<T> > mom3D(-mass*ddim);
	momentum<complex<T> > eps1D=ddim;
	momentum<complex<T> > eps2D=zero;
	momentum<complex<T> > eps3D=-ddim;
	
	
	//SgS
	complex<T> res1(complex<T>(0,-1)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
					+((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
					+((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
	
//    _MESSAGE2(" glue A1Gsc1gM1g8 : ",res1);

	return res1;
//    return complex<T>(0,-1)*(ep.p(i1)->L()*ep.p(i0)->Sm()*ep.ref()->Lt())/(ep.p(i1)->Lt()*ep.ref()->Lt());
}



template <class T> complex<T> (*A2Gsc1g_Tree_Ptr_eval(long long hc))(const eval_param<T>&, const mass_param_coll&) {
//    cout << "Found massive three-point vertex:" << hex << hc << dec << endl;
	switch (hc) {// Note that the labels are in the reverse order of the respective process
	case 0x413://G+g+G-
		return &A2Gsc1g1_eval<0,1,2>;
		break;
	case 0x341://G-G+g+
		return &A2Gsc1g1_eval<1,2,0>;
		break;
	case 0x134://g+G-G+
		return &A2Gsc1g1_eval<2,0,1>;
		break;
            
	case 0x403://G+g-G-
		return &A2Gsc1g2_eval<0,1,2>;
		break;
	case 0x340://G-G+g-
		return &A2Gsc1g2_eval<1,2,0>;
		break;
	case 0x034://g-G-G+
		return &A2Gsc1g2_eval<2,0,1>;
		break;
            
	case 0x314://G+g+G-
		return &A2Gsc1g3_eval<0,1,2>;
		break;
	case 0x431://G-G+g+
		return &A2Gsc1g3_eval<1,2,0>;
		break;
	case 0x143://g+G-G+
		return &A2Gsc1g3_eval<2,0,1>;
		break;
		
	case 0x304://G+g-G-
		return &A2Gsc1g4_eval<0,1,2>;
		break;
	case 0x430://G-G+g-
		return &A2Gsc1g4_eval<1,2,0>;
		break;
	case 0x043://g-G-G+
		return &A2Gsc1g4_eval<2,0,1>;
		break;
            
    case 0x313://G-g+G-
        return &A2Gsc1g5_eval<0,1,2>;
        break;
    case 0x331://G-G-g+
        return &A2Gsc1g5_eval<1,2,0>;
        break;
    case 0x133://g+G-G-
        return &A2Gsc1g5_eval<2,0,1>;
        break;
        
    case 0x404://G+g-G+
        return &A2Gsc1g6_eval<0,1,2>;
        break;
    case 0x440://G+G+g-
        return &A2Gsc1g6_eval<1,2,0>;
        break;
    case 0x044://g-G+G+
        return &A2Gsc1g6_eval<2,0,1>;
        break;
            
            // 1 scalar, 1 helicity
            
    case 0x312://G-g+Gsc
        return &A1Gsc1gM1g1p_eval<0,1,2>;
        break;
    case 0x123://g+GscG-
        return &A1Gsc1gM1g1p_eval<2,0,1>;
        break;
    case 0x231://GscG-g+
        return &A1Gsc1gM1g1p_eval<1,2,0>;
        break;
        
	case 0x402://G+g-Gsc
		return &A1Gsc1gM1g2p_eval<0,1,2>;
		break;
	case 0x024://g-GscG+
		return &A1Gsc1gM1g2p_eval<2,0,1>;
		break;
	case 0x240://GscG+g-
		return &A1Gsc1gM1g2p_eval<1,2,0>;
		break;
			
	case 0x132://g+G-Gsc
		return &A1Gsc1gM1g1m_eval<0,1,2>;
		break;
	case 0x321://G-Gscg+
		return &A1Gsc1gM1g1m_eval<2,0,1>;
		break;
	case 0x213://Gscg+G-
		return &A1Gsc1gM1g1m_eval<1,2,0>;
		break;
		
	case 0x042://g-G+Gsc
		return &A1Gsc1gM1g2m_eval<0,1,2>;
		break;
	case 0x420://G+Gscg-
		return &A1Gsc1gM1g2m_eval<2,0,1>;
		break;
	case 0x204://Gscg-G+
		return &A1Gsc1gM1g2m_eval<1,2,0>;
		break;
            
            // 2 scalars

    case 0x212://s(0)+s(0)
        return &A2Gsc1g7_eval<0,1,2>;
        break;
    case 0x202://s(0)-s(0)
        return &A2Gsc1g8_eval<0,1,2>;
        break;
    case 0x221://s(0)s(0)+
        return &A2Gsc1g7_eval<1,2,0>;
        break;
    case 0x220://s(0)s(0)-
        return &A2Gsc1g8_eval<1,2,0>;
        break;
    case 0x122://+s(0)s(0)
        return &A2Gsc1g7_eval<2,0,1>;
        break;
    case 0x022://-s(0)s(0)
        return &A2Gsc1g8_eval<2,0,1>;
        break;

            
	default:// We return zero for all other helicity combinations
//        cout << "NOT:" << hex << hc << dec << endl;
        return &ZeroF;
	}
}





template complex<R> (*A2Gsc1g_Tree_Ptr_eval(long long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2Gsc1g_Tree_Ptr_eval(long long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2Gsc1g_Tree_Ptr_eval(long long hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> (*A2Gsc1g_Tree_Ptr_eval(long long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif


}
