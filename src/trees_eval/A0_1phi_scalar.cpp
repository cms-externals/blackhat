/*
*  A0_1phi_scalar.cpp
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

using namespace std;

#define _HIGGS_SCALAR_SCALAR_PARAM_3 T(1)//T(-12) // ph,sc,sc
#define _HIGGS_SCALAR_SCALAR_PARAM_X (_HIGGS_SCALAR_SCALAR_PARAM_3)//T(12) // ph,sc,sc
#define _HIGGS_SCALAR_SCALAR_PARAM T(_HIGGS_SCALAR_SCALAR_PARAM_3)//T(12) // ph,sc,sc
#define _HIGGS_SCALAR_SCALAR_PARAM_D T(_HIGGS_SCALAR_SCALAR_PARAM_3)//T(12) // ph,sc,sc

#define _HIGGS_GLUON_GLUON_PARAM T(1) // Usually to include the gluon-gluon-higgs vertex

//#define BIT2 T(0) //From D-Dim result
//#define BIT1 T(0) //From D-Dim result
//#define BIT0 T(1)

//#define _HIGGS_GLUON_FAC T(1) // T(1)ph,s,g,s higgs coupling to the gluon
//#define _MASSIVE_GLUON_PIECE T(0)

namespace BH {

// Defined in A0_1phi_eval.cpp
template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses);

/*
 *
 *
 * The 1 complex higgs and two massive scalars amplitudes
 *
 *
 */
    
    
momentum<complex<R> > ddim_R(complex<R>(-1,0),complex<R>(0,0),complex<R>(0,0),complex<R>(0,0));
momentum<complex<RHP> > ddim_RHP(complex<RHP>(1,0),complex<RHP>(0,0),complex<RHP>(0,0),complex<RHP>(0,0));
momentum<complex<RVHP> > ddim_RVHP(complex<RVHP>(1,0),complex<RVHP>(0,0),complex<RVHP>(0,0),complex<RVHP>(0,0));

template <typename T> momentum<complex<T> >& get_ddim(){};
template<> momentum<complex<R> >& get_ddim(){return ddim_R;};
template<> momentum<complex<RHP> >& get_ddim(){return ddim_RHP;};
template<> momentum<complex<RVHP> >& get_ddim(){return ddim_RVHP;};

	
// Function that returns the scalar epsilon^{\mu1\mu2\mu3\mu4}p1_\mu1 p2_\mu2 p3_\mu3 p4_\mu4
template <class T> void epsS(complex<T>& pR, const momentum<complex<T> >& p0, const momentum<complex<T> >& p1,const momentum<complex<T> >& p2, const momentum<complex<T> >& p3){
	complex<T> v0,v1,v2,v3;
	v0=-p1.Z()*p2.Y()*p3.X() + p1.Y()*p2.Z()*p3.X()
	+ p1.Z()*p2.X()*p3.Y() -p1.X()*p2.Z()*p3.Y()
	-p1.Y()*p2.X()*p3.Z() +p1.X()*p2.Y()*p3.Z();
	v1=p1.Z()*p2.Y()*p3.E() - p1.Y()*p2.Z()*p3.E()
	- p1.Z()*p2.E()*p3.Y() +p1.E()*p2.Z()*p3.Y()
	+p1.Y()*p2.E()*p3.Z() -p1.E()*p2.Y()*p3.Z();
	v2=-p1.Z()*p2.X()*p3.E() + p1.X()*p2.Z()*p3.E()
	+ p1.Z()*p2.E()*p3.X() -p1.E()*p2.Z()*p3.X()
	-p1.X()*p2.E()*p3.Z() +p1.E()*p2.X()*p3.Z();
	v3=p1.Y()*p2.X()*p3.E() - p1.X()*p2.Y()*p3.E()
	- p1.Y()*p2.E()*p3.X() +p1.E()*p2.Y()*p3.X()
	+p1.X()*p2.E()*p3.Y() -p1.E()*p2.X()*p3.Y();
	// signs arranged to get covariant vector
	pR=(momentum<complex<T> >(v0,-v1,-v2,-v3)*p0);
}


//ph(d)ss
template <int i0, int i1, int i2, class T> complex<T>  A1ph2sc1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
    
    
//	return complex<T>(0,1)*T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM_3*(
//        BIT0*mass
//        -(BIT1*(ep.p(i1)->P()*ep.p(i2)->P())
//        -BIT2*mass*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//            /((ep.p(i2)->P()*ep.ref()->P())*(ep.p(i1)->P()*ep.ref()->P())))
//        );
	
//	momentum<complex<T> > eps1=(ep.mom(i1)/mass)-(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
//	momentum<complex<T> > eps2=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	
	momentum<complex<T> > ddim(get_ddim<T>());
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(-mass*ddim);
	
	
	momentum<complex<T> > eps1=(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
	momentum<complex<T> > eps2=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	momentum<complex<T> > eps1D=ddim;
	momentum<complex<T> > eps2D=-ddim;
    	
//	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM*((ep.mom(i1)*ep.mom(i2))*(eps1*eps2)-(ep.mom(i1)*eps2)*(eps1*ep.mom(i2)));
	
//	_PRINT((ep.mom(i1)*ep.mom(i2))*(eps1*eps2)-(ep.mom(i1)*eps2)*(eps1*ep.mom(i2)));
//	_PRINT((mom1D*mom2D)*(eps1D*eps2D)-(mom1D*eps2D)*(eps1D*mom2D));
	
	return complex<T>(0,-2)*_HIGGS_SCALAR_SCALAR_PARAM*(
		((ep.mom(i1)*ep.mom(i2))-(mom1D*mom1D))*((eps1*eps2)-(eps1D*eps2D))
		-((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)));

}
	
//Hss
template <int i0, int i1, int i2, class T> complex<T>  A1ph2sc1H_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	
	momentum<complex<T> > eps1=(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
	momentum<complex<T> > eps2=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());

	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim<T>());
	momentum<complex<T> > mom1D(-mass*ddim);
	momentum<complex<T> > mom2D(mass*ddim);
	momentum<complex<T> > eps1D=-ddim;
	momentum<complex<T> > eps2D=ddim;
	
//	_MESSAGE4(" A1ph2sc1_eval:",complex<T>(0,-2)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																			 ((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
//																			 -((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)))
//			  ,"=",complex<T>(0,-2)*_HIGGS_SCALAR_SCALAR_PARAM*(-mass2
//																-mass2*((ep.mom(i1)*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P())+(ep.mom(i2)*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P()))
//																+(ep.mom(i1)*ep.mom(i2))
//																));    
    
	return complex<T>(0,-2)*_HIGGS_SCALAR_SCALAR_PARAM*(
		((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
		-((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)));
	
//	return complex<T>(0,2)*_HIGGS_SCALAR_SCALAR_PARAM_3*(mass2);
//	return complex<T>(0,-2)*_HIGGS_SCALAR_SCALAR_PARAM*(
//			-mass2
//			-mass2*((ep.mom(i1)*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P())+(ep.mom(i2)*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P()))
//			+(ep.mom(i1)*ep.mom(i2))
//			);

}
	
//phgmMs
template <int i0, int i1, int i2, class T> complex<T>  A1ph2sc2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i2));
//    const complex<T> massSqrt=eval_param<T>::mass(masses.p(i2));
//    
//    const Cmom<T> p1L(ep.mom(i1)-(mass/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
//    const Cmom<T> p2L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*ep.ref()->P())))*ep.ref()->P());
//
//	return complex<T>(0,SIGN)*sqrt(T(2))*massSqrt*(p1L.L()*p2L.L())*(ep.ref()->L()*p1L.L())/(ep.ref()->L()*p2L.L());
	
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	
	const Cmom<T> pL(ep.mom(i1)-(mass2/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
	momentum<complex<T> > eps1=PfLLt(ep.ref()->Lt(),pL.L())/(ep.ref()->Lt()*pL.Lt());
	momentum<complex<T> > eps2=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim<T>());
	momentum<complex<T> > mom1D(-mass*ddim);
	momentum<complex<T> > mom2D(mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=ddim;
	
//	_MESSAGE4(" A1ph2sc2_eval:",sqrt(T(2))*complex<T>(0,-2)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																			 ((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
//																			 -((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)))
//			  ,"=",mass);

	
	return sqrt(T(2))*complex<T>(0,2)*_HIGGS_SCALAR_SCALAR_PARAM*(
														((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
														-((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)));
}

    
//phdgpMs
template <int i0, int i1, int i2, class T> complex<T>  A1ph2sc3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i2));
//    const complex<T> massSqrt=eval_param<T>::mass(masses.p(i2));
//    
//    const Cmom<T> p1L(ep.mom(i1)-(mass/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
//    const Cmom<T> p2L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*ep.ref()->P())))*ep.ref()->P());
//    
//    
//	return complex<T>(0,SIGN)*sqrt(T(2))*massSqrt*(p1L.Lt()*p2L.Lt())*(ep.ref()->Lt()*p1L.Lt())/(ep.ref()->Lt()*p2L.Lt());
	
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	
	const Cmom<T> pL(ep.mom(i1)-(mass2/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
	momentum<complex<T> > eps1=-PfLLt(pL.Lt(),ep.ref()->L())/(ep.ref()->L()*pL.L());
	momentum<complex<T> > eps2=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim<T>());
	momentum<complex<T> > mom1D(-mass*ddim);
	momentum<complex<T> > mom2D(mass*ddim);
	momentum<complex<T> > eps1D=-ddim;
	momentum<complex<T> > eps2D=zero;
	
//	_MESSAGE4(" A1ph2sc3_eval:",sqrt(T(2))*complex<T>(0,-2)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																			 ((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
//																			 -((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)))
//			  ,"=",mass);

	
	return sqrt(T(2))*complex<T>(0,2)*_HIGGS_SCALAR_SCALAR_PARAM*(
														((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
														-((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)));

}

//phgmMgmM
template <int i0, int i1, int i2, class T> complex<T>  A1ph2sc4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i2));
//
//    const Cmom<T> p1L(ep.mom(i1)-(mass/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
//    const Cmom<T> p2L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*ep.ref()->P())))*ep.ref()->P());
//    
//    return complex<T>(0,-1)*pow(p1L.L()*p2L.L(),2);
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	
	const Cmom<T> p1L(ep.mom(i1)-(mass2/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> p2L(ep.mom(i2)-(mass2/(T(2)*(ep.p(i2)->P()*ep.ref()->P())))*ep.ref()->P());
	momentum<complex<T> > eps1=PfLLt(ep.ref()->Lt(),p1L.L())/(ep.ref()->Lt()*p1L.Lt());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),p2L.L())/(ep.ref()->Lt()*p2L.Lt());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim<T>());
	momentum<complex<T> > mom1D(-mass*ddim);
	momentum<complex<T> > mom2D(mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;
    	
//	_MESSAGE4(" A1ph2sc4_eval:",complex<T>(0,-4)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																			 ((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
//																			 -((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)))
//			  ,"=",complex<T>(0,-1)*pow(p1L.L()*p2L.L(),2));

	
	return complex<T>(0,-4)*_HIGGS_SCALAR_SCALAR_PARAM*(
														((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
														-((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)));	
}

//phdgpMgpM
template <int i0, int i1, int i2, class T> complex<T>  A1ph2sc5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i2));
//    
//    const Cmom<T> p1L(ep.mom(i1)-(mass/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
//    const Cmom<T> p2L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*ep.ref()->P())))*ep.ref()->P());
//            
//    return complex<T>(0,-1)*pow(p1L.Lt()*p2L.Lt(),2);
	
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	
	const Cmom<T> p1L(ep.mom(i1)-(mass2/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> p2L(ep.mom(i2)-(mass2/(T(2)*(ep.p(i2)->P()*ep.ref()->P())))*ep.ref()->P());
	momentum<complex<T> > eps1=-PfLLt(p1L.Lt(),ep.ref()->L())/(ep.ref()->L()*p1L.L());
	momentum<complex<T> > eps2=-PfLLt(p2L.Lt(),ep.ref()->L())/(ep.ref()->L()*p2L.L());
	
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(get_ddim<T>());
	momentum<complex<T> > mom1D(-mass*ddim);
	momentum<complex<T> > mom2D(mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;

    	
//	_MESSAGE4(" A1ph2sc5_eval:",complex<T>(0,-4)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																			 ((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
//																			 -((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)))
//			  ,"=",complex<T>(0,-1)*pow(p1L.Lt()*p2L.Lt(),2));
    
	
	return complex<T>(0,-4)*_HIGGS_SCALAR_SCALAR_PARAM*(
														((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
														-((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)));
}
	
//phdgmMgpM
template <int i0, int i1, int i2, class T> complex<T>  A1ph2sc6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	//    const complex<T> mass=eval_param<T>::mass2(masses.p(i2));
	//    
	//    const Cmom<T> p1L(ep.mom(i1)-(mass/(T(2)*(ep.p(i1)->P()*ep.ref()->P())))*ep.ref()->P());
	//    const Cmom<T> p2L(ep.mom(i2)-(mass/(T(2)*(ep.p(i2)->P()*ep.ref()->P())))*ep.ref()->P());
	//            
	//    return complex<T>(0,-1)*pow(p1L.Lt()*p2L.Lt(),2);
	
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	
	const Cmom<T> p1L(ep.mom(i1)-(mass2/(T(2)*(ep.mom(i1)*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> p2L(ep.mom(i2)-(mass2/(T(2)*(ep.mom(i2)*ep.ref()->P())))*ep.ref()->P());
	momentum<complex<T> > eps1=PfLLt(ep.ref()->Lt(),p1L.L())/(ep.ref()->Lt()*p1L.Lt());
	momentum<complex<T> > eps2=-PfLLt(p2L.Lt(),ep.ref()->L())/(ep.ref()->L()*p2L.L());
	
	momentum<complex<T> > ddim(get_ddim<T>());
	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > mom1D(-mass*ddim);
	momentum<complex<T> > mom2D(mass*ddim);
	momentum<complex<T> > eps1D=zero;
	momentum<complex<T> > eps2D=zero;

    
//    _MESSAGE2(" A1ph2sc6_eval:",complex<T>(0,-4)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																			 ((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
//																			 -((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D))));
//
    
	return complex<T>(0,-4)*_HIGGS_SCALAR_SCALAR_PARAM*(
														((ep.mom(i1)*ep.mom(i2))-(mom1D*mom2D))*((eps1*eps2)-(eps1D*eps2D))
														-((ep.mom(i1)*eps2)-(mom1D*eps2D))*((eps1*ep.mom(i2))-(eps1D*mom2D)));
}


    
template <class T> complex<T> (*A1ph2sc_Tree_Ptr_eval(long long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 0x944://phss
        case 0xA44://phdss
			return &A1ph2sc1_eval<0,1,2>;
			break;
	case 0xE44://Hss
			return &A1ph2sc1H_eval<0,1,2>;
			break;
	case 0x494://sphs
        case 0x4A4://sphds
			return &A1ph2sc1_eval<1,0,2>;
			break;
        case 0x4E4://sHs
			return &A1ph2sc1H_eval<1,0,2>;
			break;
	case 0x449://ssph
        case 0x44A://ssphd
			return &A1ph2sc1_eval<2,0,1>;
			break;
        case 0x44E://ssH
			return &A1ph2sc1H_eval<2,0,1>;
			break;            
            
            
//        case 0x9B4://phgmMs
        case 0xEB4://phdgmMs
			return &A1ph2sc2_eval<0,1,2>;
			break;
//		case 0xB94://gmMphs
        case 0xBE4://gmMphds
			return &A1ph2sc2_eval<1,0,2>;
			break;
//		case 0xB49://gmMsph
        case 0xB4E://gmMsphd
			return &A1ph2sc2_eval<2,0,1>;
			break; 
            
//        case 0x94B://phsgmM
        case 0xE4B://phdsgmM
			return &A1ph2sc2_eval<0,2,1>;
			break;
//		case 0x49B://sgmMph
        case 0x4EB://sgmMphd
			return &A1ph2sc2_eval<1,2,0>;
			break;
//		case 0x4B9://sgmMph
        case 0x4BE://sgmMphd
			return &A1ph2sc2_eval<2,1,0>;
			break; 
            
            
        case 0xEC4://HgpMs
//        case 0xAC4://phdgpMs
			return &A1ph2sc3_eval<0,1,2>;
			break;
		case 0xCE4://gpMHs
//        case 0xCA4://gpMphds
			return &A1ph2sc3_eval<1,0,2>;
			break;
		case 0xC4E://gpMsH
//        case 0xC4A://gpMsphd
			return &A1ph2sc3_eval<2,0,1>;
			break; 
            
        case 0xE4C://HsgpM
//        case 0xA4C://phdsgpM
			return &A1ph2sc3_eval<0,2,1>;
			break;
		case 0x4EC://sHgpM
//        case 0x4AC://sgpMphd
			return &A1ph2sc3_eval<1,2,0>;
			break;
		case 0x4CE://sgpMH
//        case 0x4CA://sgpMphd
			return &A1ph2sc3_eval<2,1,0>;
			break; 
            
            
//        case 0x9BB://phgmMgmM
		case 0xEBB://HgmMgmM
			return &A1ph2sc4_eval<0,1,2>;
			break;
//		case 0xB9B://gmMphgmM
		case 0xBEB://gmMGgmM
			return &A1ph2sc4_eval<1,0,2>;
			break;
//		case 0xBB9://gmMgmMph
		case 0xBBE://gmMgmMG
			return &A1ph2sc4_eval<2,0,1>;
			break; 
            
//        case 0xACC://phdgpMgpM
        case 0xECC://HgpMgpM
			return &A1ph2sc5_eval<0,1,2>;
			break;
//        case 0xCAC://gpMphdgpM
        case 0xCEC://gpMHgpM
			return &A1ph2sc5_eval<1,0,2>;
			break;
//        case 0xCCA://gpMgpMphd
        case 0xCCE://gpMgpMH
			return &A1ph2sc5_eval<2,0,1>;
			break; 
			
		case 0xEBC://HgmMgpM
			return &A1ph2sc6_eval<0,1,2>;
			break;
		case 0xBEC://gmMHgpM
			return &A1ph2sc6_eval<1,0,2>;
			break;
		case 0xBCE://gmMgpMH
			return &A1ph2sc6_eval<2,0,1>;
			break; 
			
		case 0xECB://HgpMgmM
			return &A1ph2sc6_eval<0,2,1>;
			break;
		case 0xCEB://gpMHgmM
			return &A1ph2sc6_eval<1,2,0>;
			break;
		case 0xCBE://gpMgmMH
			return &A1ph2sc6_eval<2,1,0>;
			break; 


						
		default:// We return zero for all other helicity combinations
			cout << hex << " missing 3pt:" << hc << dec << endl;
			return &ZeroF;
	}
}


/*
 *
 *
 * The 1 complex higgs, two massive scalars and a gluon amplitude
 *
 *
 */
	
// Function that returns the vector v_ii^\mu1 = epsilon^{\mu1\mu2\mu3\mu4}p1_\mu2 p2_\mu3 p3_\mu4
template <class T> void epsC(momentum<complex<T> >& p0, const momentum<complex<T> >& p1,const momentum<complex<T> >& p2, const momentum<complex<T> >& p3){
	complex<T> v0,v1,v2,v3;
	v0=-p1.Z()*p2.Y()*p3.X() + p1.Y()*p2.Z()*p3.X()
	+ p1.Z()*p2.X()*p3.Y() -p1.X()*p2.Z()*p3.Y()
	-p1.Y()*p2.X()*p3.Z() +p1.X()*p2.Y()*p3.Z();
	v1=p1.Z()*p2.Y()*p3.E() - p1.Y()*p2.Z()*p3.E()
	- p1.Z()*p2.E()*p3.Y() +p1.E()*p2.Z()*p3.Y()
	+p1.Y()*p2.E()*p3.Z() -p1.E()*p2.Y()*p3.Z();
	v2=-p1.Z()*p2.X()*p3.E() + p1.X()*p2.Z()*p3.E()
	+ p1.Z()*p2.E()*p3.X() -p1.E()*p2.Z()*p3.X()
	-p1.X()*p2.E()*p3.Z() +p1.E()*p2.X()*p3.Z();
	v3=p1.Y()*p2.X()*p3.E() - p1.X()*p2.Y()*p3.E()
	- p1.Y()*p2.E()*p3.X() +p1.E()*p2.Y()*p3.X()
	+p1.X()*p2.E()*p3.Y() -p1.E()*p2.X()*p3.Y();
	// signs arranged to get covariant vector
	p0=momentum<complex<T> >(v0,-v1,-v2,-v3);
}

#define S(i,j) (ep.p(i)->P()*ep.p(j)->P())

//phsms
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i1));
	const complex<T> mass=eval_param<T>::mass(masses.p(i1));

	const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();
    const momentum<complex<T> > m12=ep.p(i1)->P()+ep.p(i2)->P();
    const momentum<complex<T> > m23=ep.p(i2)->P()+ep.p(i3)->P();
    const smatrix<T> Sm12=Cmom<T>(m12).Sm();
    const smatrix<T> Sm23=Cmom<T>(m23).Sm();
    
//	return complex<T>(0,1)*(
//     T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM*(
//      (BIT0*mass2+(BIT1*(m23*ep.p(i1)->P())
//       -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//        /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))))
//      *(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(T(2)*S(i3,i2))
//      -(BIT0*mass2+(BIT1*(m12*ep.p(i3)->P())
//       -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//       /((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))))
//      *(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())/(T(2)*S(i1,i2))        
//        )/(ep.p(i2)->Lt()*ep.ref()->Lt())
//        +_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->Lt()*ep.ref()->Lt())))*(
//        (ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())
//        -(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.p(i2)->Lt()))
//        +_MASSIVE_GLUON_PIECE*(mass2*pow(ep.p(i2)->L()*ep.ref()->L(),2)
//            /((ep.ref()->P()*ep.p(i1)->P())*(ep.ref()->P()*ep.p(i2)->P()))
//            *(ep.ref()->L()*ep.p(i3)->Sm()*Sm12*ep.ref()->L())/(T(2)*S(i1,i2))
//            +(ep.ref()->L()*ep.p(i1)->Sm()*Sm23*ep.ref()->L())/(T(2)*S(i3,i2)))
//        );

	
	momentum<complex<T> > eps1=(ep.mom(i1)/mass)-(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i2)->L())/(ep.ref()->Lt()*ep.p(i2)->Lt());
	momentum<complex<T> > eps3=(ep.mom(i3)/mass)-(mass*ep.ref()->P())/(ep.mom(i3)*ep.ref()->P());
	
	//SHS-SSg
	momentum<complex<T> > propleg1=-(ep.mom(i2)+ep.mom(i3));
	momentum<complex<T> > lc_higgsvert1;
	epsC(lc_higgsvert1,eps1,propleg1,ep.mom(i1));
	momentum<complex<T> > propeps1=(propleg1/mass)-(mass*ep.ref()->P())/(propleg1*ep.ref()->P());
	complex<T> higgsvert1=_HIGGS_SCALAR_SCALAR_PARAM*(propeps1*((propleg1*ep.mom(i1))*eps1-(eps1*propleg1)*ep.mom(i1)+complex<T>(0,1)*lc_higgsvert1))/(propleg1.square()-mass2);
	complex<T> res1(higgsvert1*((propeps1*eps2)*((propleg1-ep.mom(i2))*eps3)
					+(eps2*eps3)*((ep.mom(i2)-ep.mom(i3))*propeps1)
					+(eps3*propeps1)*((ep.mom(i3)-propleg1)*eps2)));
	//gHg-gSS
	momentum<complex<T> > propleg2=-(ep.mom(i1)+ep.mom(i3));
	momentum<complex<T> > lc_higgsvert2;
	epsC(lc_higgsvert2,eps2,propleg2,ep.mom(i2));
	momentum<complex<T> > higgsvert2=((propleg2*ep.mom(i2))*eps2-(eps2*propleg2)*ep.mom(i2)+complex<T>(0,1)*lc_higgsvert2)/propleg2.square();
	complex<T> res2new((eps1*higgsvert2)*((ep.mom(i1)-propleg2)*eps3)
					+(higgsvert2*eps3)*((propleg2-ep.mom(i3))*eps1)
					+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*higgsvert2));
//	_MESSAGE4("res2 bits:",complex<T>(0,1)*((eps1*eps2)*((ep.mom(i1)-propleg2)*eps3)
//											+(eps2*eps3)*((propleg2-ep.mom(i3))*eps1)
//											+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2))
//								,"=?=",complex<T>(0,-1)*((ep.mom(i1)-ep.mom(i3))*eps2));
	complex<T> res2(_HIGGS_GLUON_GLUON_PARAM*((ep.mom(i1)-ep.mom(i3))*higgsvert2));
	//SHS-SSg
	momentum<complex<T> > propleg3=-(ep.mom(i1)+ep.mom(i2));
	momentum<complex<T> > lc_higgsvert3;
	epsC(lc_higgsvert3,eps3,propleg3,ep.mom(i3));
	momentum<complex<T> > propeps3=(propleg3/mass)-(mass*ep.ref()->P())/(propleg3*ep.ref()->P());
	complex<T> higgsvert3=_HIGGS_SCALAR_SCALAR_PARAM*(propeps3*((propleg3*ep.mom(i3))*eps3-(eps3*propleg3)*ep.mom(i3)+complex<T>(0,1)*lc_higgsvert3))/(propleg3.square()-mass2);
	complex<T> res3(higgsvert3*((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*propeps3)
					+(eps2*propeps3)*((ep.mom(i2)-propleg3)*eps1)
					+(propeps3*eps1)*((propleg3-ep.mom(i1))*eps2)));
	//SHgS
	complex<T> res4((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*eps3)
					+(eps2*eps3)*((ep.mom(i2)-ep.mom(i2))*eps1)
					+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2));

//	_PRINT(complex<T>(0,1)*(res1+res3));
////	_MESSAGE2(" res1+res3 simplify to:",complex<T>(0,-2)*((ep.mom(i3)*eps2)*(mass2)/(propleg1.square()-mass2)
//										   -(ep.mom(i1)*eps2)*(mass2)/(propleg3.square()-mass2)));
//
//	_PRINT(complex<T>(0,1)*res2);
//	_PRINT(complex<T>(0,-1)*(ep.mom(i1)-ep.mom(i3))*higgsvert2);
//	_PRINT(complex<T>(0,1)*res4);
	
//	_MESSAGE6(" OLD phSMS:",complex<T>(0,1)*(
//											 T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																					 (BIT0*mass2+(BIT1*(m23*ep.p(i1)->P())
//																								  -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																								  /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))))
//																					 *(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(T(2)*S(i3,i2))
//																					 -(BIT0*mass2+(BIT1*(m12*ep.p(i3)->P())
//																								   -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																								   /((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))))
//																					 *(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())/(T(2)*S(i1,i2))        
//																					 )/(ep.p(i2)->Lt()*ep.ref()->Lt())
//											 +_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->Lt()*ep.ref()->Lt())))*(
//																											 (ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())
//																											 -(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.p(i2)->Lt()))
//											 ),
//			  " NEW:",complex<T>(0,1)*(res1+res3),"+",complex<T>(0,1)*res2);
	
	return complex<T>(0,1)*(res1+res2+res3);
}

//phsps
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i1));
	const complex<T> mass=eval_param<T>::mass(masses.p(i1));
    
    const momentum<complex<T> > m12=ep.p(i1)->P()+ep.p(i2)->P();
    const momentum<complex<T> > m23=ep.p(i2)->P()+ep.p(i3)->P();
            
//    return complex<T>(0,-1)*(
//    T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM_D*(
//      (BIT0*mass2
//      +(BIT1*(m23*ep.p(i1)->P())
//       -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//         /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))
//         )
//        )
//       *(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//       -(BIT0*mass2
//        +(BIT1*(m12*ep.p(i3)->P())
//       -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//         /((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))
//          )
//        )
//       *(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))
//      )/(ep.p(i2)->L()*ep.ref()->L())
//    );
	
	momentum<complex<T> > eps1=(ep.mom(i1)/mass)-(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
	momentum<complex<T> > eps2=-PfLLt(ep.p(i2)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i2)->L());
	momentum<complex<T> > eps3=(ep.mom(i3)/mass)-(mass*ep.ref()->P())/(ep.mom(i3)*ep.ref()->P());
	
	//SHS-SSg
	momentum<complex<T> > propleg1=-(ep.mom(i2)+ep.mom(i3));
	momentum<complex<T> > lc_higgsvert1;
	epsC(lc_higgsvert1,eps1,propleg1,ep.mom(i1));
	momentum<complex<T> > propeps1=(propleg1/mass)-(mass*ep.ref()->P())/(propleg1*ep.ref()->P());
	complex<T> higgsvert1=_HIGGS_SCALAR_SCALAR_PARAM*(propeps1*((propleg1*ep.mom(i1))*eps1-(eps1*propleg1)*ep.mom(i1)+complex<T>(0,1)*lc_higgsvert1))/(propleg1.square()-mass2);
	complex<T> res1(higgsvert1*((propeps1*eps2)*((propleg1-ep.mom(i2))*eps3)
								+(eps2*eps3)*((ep.mom(i2)-ep.mom(i3))*propeps1)
								+(eps3*propeps1)*((ep.mom(i3)-propleg1)*eps2)));
	//gHg-gSS
	momentum<complex<T> > propleg2=-(ep.mom(i1)+ep.mom(i3));
	momentum<complex<T> > lc_higgsvert2;
	epsC(lc_higgsvert2,eps2,propleg2,ep.mom(i2));
	momentum<complex<T> > higgsvert2=((propleg2*ep.mom(i2))*eps2-(eps2*propleg2)*ep.mom(i2)+complex<T>(0,1)*lc_higgsvert2)/propleg2.square();
	complex<T> res2new((eps1*higgsvert2)*((ep.mom(i1)-propleg2)*eps3)
					   +(higgsvert2*eps3)*((propleg2-ep.mom(i3))*eps1)
					   +(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*higgsvert2));
	//	_MESSAGE4("res2 bits:",complex<T>(0,1)*((eps1*eps2)*((ep.mom(i1)-propleg2)*eps3)
	//											+(eps2*eps3)*((propleg2-ep.mom(i3))*eps1)
	//											+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2))
	//								,"=?=",complex<T>(0,-1)*((ep.mom(i1)-ep.mom(i3))*eps2));
	complex<T> res2(_HIGGS_GLUON_GLUON_PARAM*((ep.mom(i1)-ep.mom(i3))*higgsvert2));
	//SHS-SSg
	momentum<complex<T> > propleg3=-(ep.mom(i1)+ep.mom(i2));
	momentum<complex<T> > lc_higgsvert3;
	epsC(lc_higgsvert3,eps3,propleg3,ep.mom(i3));
	momentum<complex<T> > propeps3=(propleg3/mass)-(mass*ep.ref()->P())/(propleg3*ep.ref()->P());
	complex<T> higgsvert3=_HIGGS_SCALAR_SCALAR_PARAM*(propeps3*((propleg3*ep.mom(i3))*eps3-(eps3*propleg3)*ep.mom(i3)+complex<T>(0,1)*lc_higgsvert3))/(propleg3.square()-mass2);
	complex<T> res3(higgsvert3*((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*propeps3)
								+(eps2*propeps3)*((ep.mom(i2)-propleg3)*eps1)
								+(propeps3*eps1)*((propleg3-ep.mom(i1))*eps2)));
	//SHgS
	complex<T> res4((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*eps3)
					+(eps2*eps3)*((ep.mom(i2)-ep.mom(i2))*eps1)
					+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2));
	
//	_PRINT(complex<T>(0,1)*(res1+res3));
	//	_MESSAGE2(" res1+res3 simplify to:",complex<T>(0,-2)*((ep.mom(i3)*eps2)*(mass2)/(propleg1.square()-mass2)
	//										   -(ep.mom(i1)*eps2)*(mass2)/(propleg3.square()-mass2)));
	
//	_PRINT(complex<T>(0,1)*res2);
	//	_PRINT(complex<T>(0,-1)*(ep.mom(i1)-ep.mom(i3))*higgsvert2);
//	_PRINT(complex<T>(0,1)*res4);
	
//	_MESSAGE4(" OLD phSPS:",complex<T>(0,-1)*(
//											  T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM_D*(
//																						(BIT0*mass2
//																						 +(BIT1*(m23*ep.p(i1)->P())
//																						   -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																						   /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))
//																						   )
//																						 )
//																						*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//																						-(BIT0*mass2
//																						  +(BIT1*(m12*ep.p(i3)->P())
//																							-BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																							/((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))
//																							)
//																						  )
//																						*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))
//																						)/(ep.p(i2)->L()*ep.ref()->L())
//											  ),
//			  " NEW:",complex<T>(0,1)*(res1+res2+res3));
	
	return complex<T>(0,1)*(res1+res2+res3);
}
    
//phdsms
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i1));
	const complex<T> mass=eval_param<T>::mass(masses.p(i1));
    
    const momentum<complex<T> > m12=ep.p(i1)->P()+ep.p(i2)->P();
    const momentum<complex<T> > m23=ep.p(i2)->P()+ep.p(i3)->P();
    
//    return complex<T>(0,1)*(
//    T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM_D*(
//       (BIT0*mass2
//        +(BIT1*(m23*ep.p(i1)->P())
//         -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//           /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))
//          )
//        )
//         *(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(T(2)*S(i3,i2))
//       -(BIT0*mass2
//         +(BIT1*(m12*ep.p(i3)->P())
//          -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//         /((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))
//           )
//         )
//        *(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())/(T(2)*S(i1,i2))
//    )/(ep.p(i2)->Lt()*ep.ref()->Lt())
//    );
	
	momentum<complex<T> > eps1=(ep.mom(i1)/mass)-(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i2)->L())/(ep.ref()->Lt()*ep.p(i2)->Lt());
	momentum<complex<T> > eps3=(ep.mom(i3)/mass)-(mass*ep.ref()->P())/(ep.mom(i3)*ep.ref()->P());
	
	//SHS-SSg
	momentum<complex<T> > propleg1=-(ep.mom(i2)+ep.mom(i3));
	momentum<complex<T> > lc_higgsvert1;
	epsC(lc_higgsvert1,eps1,propleg1,ep.mom(i1));
	momentum<complex<T> > propeps1=(propleg1/mass)-(mass*ep.ref()->P())/(propleg1*ep.ref()->P());
	complex<T> higgsvert1=_HIGGS_SCALAR_SCALAR_PARAM*(propeps1*((propleg1*ep.mom(i1))*eps1-(eps1*propleg1)*ep.mom(i1)-complex<T>(0,1)*lc_higgsvert1))/(propleg1.square()-mass2);
	complex<T> res1(higgsvert1*((propeps1*eps2)*((propleg1-ep.mom(i2))*eps3)
								+(eps2*eps3)*((ep.mom(i2)-ep.mom(i3))*propeps1)
								+(eps3*propeps1)*((ep.mom(i3)-propleg1)*eps2)));
	//gHg-gSS
	momentum<complex<T> > propleg2=-(ep.mom(i1)+ep.mom(i3));
	momentum<complex<T> > lc_higgsvert2;
	epsC(lc_higgsvert2,eps2,propleg2,ep.mom(i2));
	momentum<complex<T> > higgsvert2=((propleg2*ep.mom(i2))*eps2-(eps2*propleg2)*ep.mom(i2)-complex<T>(0,1)*lc_higgsvert2)/propleg2.square();
	complex<T> res2new((eps1*higgsvert2)*((ep.mom(i1)-propleg2)*eps3)
					   +(higgsvert2*eps3)*((propleg2-ep.mom(i3))*eps1)
					   +(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*higgsvert2));
	//	_MESSAGE4("res2 bits:",complex<T>(0,1)*((eps1*eps2)*((ep.mom(i1)-propleg2)*eps3)
	//											+(eps2*eps3)*((propleg2-ep.mom(i3))*eps1)
	//											+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2))
	//								,"=?=",complex<T>(0,-1)*((ep.mom(i1)-ep.mom(i3))*eps2));
	complex<T> res2(_HIGGS_GLUON_GLUON_PARAM*((ep.mom(i1)-ep.mom(i3))*higgsvert2));
	//SHS-SSg
	momentum<complex<T> > propleg3=-(ep.mom(i1)+ep.mom(i2));
	momentum<complex<T> > lc_higgsvert3;
	epsC(lc_higgsvert3,eps3,propleg3,ep.mom(i3));
	momentum<complex<T> > propeps3=(propleg3/mass)-(mass*ep.ref()->P())/(propleg3*ep.ref()->P());
	complex<T> higgsvert3=_HIGGS_SCALAR_SCALAR_PARAM*(propeps3*((propleg3*ep.mom(i3))*eps3-(eps3*propleg3)*ep.mom(i3)-complex<T>(0,1)*lc_higgsvert3))/(propleg3.square()-mass2);
	complex<T> res3(higgsvert3*((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*propeps3)
								+(eps2*propeps3)*((ep.mom(i2)-propleg3)*eps1)
								+(propeps3*eps1)*((propleg3-ep.mom(i1))*eps2)));
	//SHgS
	complex<T> res4((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*eps3)
					+(eps2*eps3)*((ep.mom(i2)-ep.mom(i2))*eps1)
					+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2));
	
//	_PRINT(complex<T>(0,1)*(res1+res3));
	//	_MESSAGE2(" res1+res3 simplify to:",complex<T>(0,-2)*((ep.mom(i3)*eps2)*(mass2)/(propleg1.square()-mass2)
	//										   -(ep.mom(i1)*eps2)*(mass2)/(propleg3.square()-mass2)));
	
//	_PRINT(complex<T>(0,1)*res2);
	//	_PRINT(complex<T>(0,-1)*(ep.mom(i1)-ep.mom(i3))*higgsvert2);
//	_PRINT(complex<T>(0,1)*res4);
	
//	_MESSAGE4(" OLD phdSMS:",complex<T>(0,1)*(
//											 T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM_D*(
//																					   (BIT0*mass2
//																						+(BIT1*(m23*ep.p(i1)->P())
//																						  -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																						  /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))
//																						  )
//																						)
//																					   *(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(T(2)*S(i3,i2))
//																					   -(BIT0*mass2
//																						 +(BIT1*(m12*ep.p(i3)->P())
//																						   -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																						   /((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))
//																						   )
//																						 )
//																					   *(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())/(T(2)*S(i1,i2))
//																					   )/(ep.p(i2)->Lt()*ep.ref()->Lt())
//											 ),
//			  " NEW:",complex<T>(0,1)*(res1+res2+res3));
	
	return complex<T>(0,1)*(res1+res2+res3);
	
}
	
//phdsps
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i1));
	const complex<T> mass=eval_param<T>::mass(masses.p(i1));

    const momentum<complex<T> > m12=ep.p(i1)->P()+ep.p(i2)->P();
    const momentum<complex<T> > m23=ep.p(i2)->P()+ep.p(i3)->P();
        
	const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();
    
//	return complex<T>(0,-1)*(
//    T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM*(
//      (BIT0*mass2+(BIT1*(m23*ep.p(i1)->P())
//           -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//           /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))))
//        *(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//      -(BIT0*mass2+(BIT1*(m12*ep.p(i3)->P())
//           -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//           /((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))))
//        *(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))      
//       )/(ep.p(i2)->L()*ep.ref()->L())
//     +_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->L()*ep.ref()->L())))*(
//        (ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())
//        -(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.p(i2)->L()))
//     );
	
	momentum<complex<T> > eps1=(ep.mom(i1)/mass)-(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
	momentum<complex<T> > eps2=-PfLLt(ep.p(i2)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i2)->L());
	momentum<complex<T> > eps3=(ep.mom(i3)/mass)-(mass*ep.ref()->P())/(ep.mom(i3)*ep.ref()->P());
	
	//SHS-SSg
	momentum<complex<T> > propleg1=-(ep.mom(i2)+ep.mom(i3));
	momentum<complex<T> > lc_higgsvert1;
	epsC(lc_higgsvert1,eps1,propleg1,ep.mom(i1));
	momentum<complex<T> > propeps1=(propleg1/mass)-(mass*ep.ref()->P())/(propleg1*ep.ref()->P());
	complex<T> higgsvert1=_HIGGS_SCALAR_SCALAR_PARAM*(propeps1*((propleg1*ep.mom(i1))*eps1-(eps1*propleg1)*ep.mom(i1)-complex<T>(0,1)*lc_higgsvert1))/(propleg1.square()-mass2);
	complex<T> res1(higgsvert1*((propeps1*eps2)*((propleg1-ep.mom(i2))*eps3)
								+(eps2*eps3)*((ep.mom(i2)-ep.mom(i3))*propeps1)
								+(eps3*propeps1)*((ep.mom(i3)-propleg1)*eps2)));
	//gHg-gSS
	momentum<complex<T> > propleg2=-(ep.mom(i1)+ep.mom(i3));
	momentum<complex<T> > lc_higgsvert2;
	epsC(lc_higgsvert2,eps2,propleg2,ep.mom(i2));
	momentum<complex<T> > higgsvert2=((propleg2*ep.mom(i2))*eps2-(eps2*propleg2)*ep.mom(i2)-complex<T>(0,1)*lc_higgsvert2)/propleg2.square();
	complex<T> res2new((eps1*higgsvert2)*((ep.mom(i1)-propleg2)*eps3)
					   +(higgsvert2*eps3)*((propleg2-ep.mom(i3))*eps1)
					   +(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*higgsvert2));
	//	_MESSAGE4("res2 bits:",complex<T>(0,1)*((eps1*eps2)*((ep.mom(i1)-propleg2)*eps3)
	//											+(eps2*eps3)*((propleg2-ep.mom(i3))*eps1)
	//											+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2))
	//								,"=?=",complex<T>(0,-1)*((ep.mom(i1)-ep.mom(i3))*eps2));
	complex<T> res2(_HIGGS_GLUON_GLUON_PARAM*((ep.mom(i1)-ep.mom(i3))*higgsvert2));
	//SHS-SSg
	momentum<complex<T> > propleg3=-(ep.mom(i1)+ep.mom(i2));
	momentum<complex<T> > lc_higgsvert3;
	epsC(lc_higgsvert3,eps3,propleg3,ep.mom(i3));
	momentum<complex<T> > propeps3=(propleg3/mass)-(mass*ep.ref()->P())/(propleg3*ep.ref()->P());
	complex<T> higgsvert3=_HIGGS_SCALAR_SCALAR_PARAM*(propeps3*((propleg3*ep.mom(i3))*eps3-(eps3*propleg3)*ep.mom(i3)-complex<T>(0,1)*lc_higgsvert3))/(propleg3.square()-mass2);
	complex<T> res3(higgsvert3*((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*propeps3)
								+(eps2*propeps3)*((ep.mom(i2)-propleg3)*eps1)
								+(propeps3*eps1)*((propleg3-ep.mom(i1))*eps2)));
	//SHgS
	complex<T> res4((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*eps3)
					+(eps2*eps3)*((ep.mom(i2)-ep.mom(i2))*eps1)
					+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2));
	
//	_PRINT(complex<T>(0,1)*(res1+res3));
	//	_MESSAGE2(" res1+res3 simplify to:",complex<T>(0,-2)*((ep.mom(i3)*eps2)*(mass2)/(propleg1.square()-mass2)
	//										   -(ep.mom(i1)*eps2)*(mass2)/(propleg3.square()-mass2)));
	
//	_PRINT(complex<T>(0,1)*res2);
	//	_PRINT(complex<T>(0,-1)*(ep.mom(i1)-ep.mom(i3))*higgsvert2);
//	_PRINT(complex<T>(0,1)*res4);
//	_MESSAGE4(" OLD phdSPS:",complex<T>(0,-1)*(
//											  T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																					  (BIT0*mass2+(BIT1*(m23*ep.p(i1)->P())
//																								   -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																								   /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))))
//																					  *(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//																					  -(BIT0*mass2+(BIT1*(m12*ep.p(i3)->P())
//																									-BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																									/((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))))
//																					  *(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))      
//																					  )/(ep.p(i2)->L()*ep.ref()->L())
//											  +_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->L()*ep.ref()->L())))*(
//																											(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())
//																											-(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.p(i2)->L()))
//											  ),
//			  " NEW:",complex<T>(0,1)*(res1+res2+res3));
	
	
	return complex<T>(0,1)*(res1+res2+res3);
}
	
//Hsms
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g1H_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i1));
	const complex<T> mass=eval_param<T>::mass(masses.p(i1));

	const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();

	//	momentum<complex<T> > ddim(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > ddim(complex<T>(1,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	momentum<complex<T> > mom1D(mass*ddim);
	momentum<complex<T> > mom2D(-mass*ddim);
	momentum<complex<T> > eps1D=ddim;
	momentum<complex<T> > eps2D=-ddim;
	
	momentum<complex<T> > eps1=(ep.mom(i1)/mass)-(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i2)->L())/(ep.ref()->Lt()*ep.p(i2)->Lt());
	momentum<complex<T> > eps3=(ep.mom(i3)/mass)-(mass*ep.ref()->P())/(ep.mom(i3)*ep.ref()->P());
//	momentum<complex<T> > eps1=(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
//	momentum<complex<T> > eps2=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());

	
	//SHS-SSg
	momentum<complex<T> > propleg1=-(ep.mom(i2)+ep.mom(i3));
	momentum<complex<T> > propeps1=(propleg1/mass)-(mass*ep.ref()->P())/(propleg1*ep.ref()->P());
	complex<T> higgsvert1=_HIGGS_SCALAR_SCALAR_PARAM*(propeps1*((propleg1*ep.mom(i1))*eps1-(eps1*propleg1)*ep.mom(i1)))/(propleg1.square()-mass2);
	complex<T> res1(higgsvert1*((propeps1*eps2)*((propleg1-ep.mom(i2))*eps3)
								+(eps2*eps3)*((ep.mom(i2)-ep.mom(i3))*propeps1)
								+(eps3*propeps1)*((ep.mom(i3)-propleg1)*eps2)));
	//gHg-gSS
	momentum<complex<T> > propleg2=-(ep.mom(i1)+ep.mom(i3));
	momentum<complex<T> > higgsvert2=((propleg2*ep.mom(i2))*eps2-(eps2*propleg2)*ep.mom(i2))/propleg2.square();
	complex<T> res2new((eps1*higgsvert2)*((ep.mom(i1)-propleg2)*eps3)
					   +(higgsvert2*eps3)*((propleg2-ep.mom(i3))*eps1)
					   +(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*higgsvert2));
//	_MESSAGE4("res2 bits:",complex<T>(0,1)*((eps1)*((ep.mom(i1)-propleg2)*eps3))*(ep.mom(i1)),
//			  complex<T>(0,1)*((eps3)*((propleg2-ep.mom(i3))*eps1))*(ep.mom(i1)),
//			  complex<T>(0,1)*((eps3*eps1)*((ep.mom(i3)-ep.mom(i1))))*(ep.mom(i1)));
//	_MESSAGE4("=",complex<T>(0,1)*((eps1)*((ep.mom(i1)-propleg2)*eps3)
//							   +(eps3)*((propleg2-ep.mom(i3))*eps1)
//							   +(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))))*(ep.mom(i1))
//								,"=?=",complex<T>(0,-1)*((ep.mom(i1)-ep.mom(i3))*(ep.mom(i1))));
	complex<T> res2(_HIGGS_GLUON_GLUON_PARAM*((ep.mom(i1)-ep.mom(i3))*higgsvert2));
	//SHS-SSg
	momentum<complex<T> > propleg3=-(ep.mom(i1)+ep.mom(i2));
	momentum<complex<T> > propeps3=(propleg3/mass)-(mass*ep.ref()->P())/(propleg3*ep.ref()->P());
	complex<T> higgsvert3=_HIGGS_SCALAR_SCALAR_PARAM*(propeps3*((propleg3*ep.mom(i3))*eps3-(eps3*propleg3)*ep.mom(i3)))/(propleg3.square()-mass2);
	complex<T> res3(higgsvert3*((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*propeps3)
								+(eps2*propeps3)*((ep.mom(i2)-propleg3)*eps1)
								+(propeps3*eps1)*((propleg3-ep.mom(i1))*eps2)));
	//SHgS
	complex<T> res4((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*eps3)
					+(eps2*eps3)*((ep.mom(i2)-ep.mom(i2))*eps1)
					+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2));
	
//	_PRINT(complex<T>(0,1)*(res1+res3));
//	_MESSAGE2(" res1+res3 simplify to:",complex<T>(0,-2)*((ep.mom(i3)*eps2)*(mass2)/(propleg1.square()-mass2)
//										   -(ep.mom(i1)*eps2)*(mass2)/(propleg3.square()-mass2)));

//	_PRINT(complex<T>(0,1)*res2);
//	_PRINT(complex<T>(0,-1)*(ep.mom(i1)-ep.mom(i3))*higgsvert2);
//	_PRINT(complex<T>(0,1)*res4);
//	_MESSAGE4(" OLD phdSPS:",complex<T>(0,-1)*(
//											  T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																					  (BIT0*mass2+(BIT1*(m23*ep.p(i1)->P())
//																								   -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																								   /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))))
//																					  *(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//																					  -(BIT0*mass2+(BIT1*(m12*ep.p(i3)->P())
//																									-BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																									/((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))))
//																					  *(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))      
//																					  )/(ep.p(i2)->L()*ep.ref()->L())
//											  +_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->L()*ep.ref()->L())))*(
//																											(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())
//																											-(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.p(i2)->L()))
//											  ),
//			  " NEW:",complex<T>(0,1)*(res1+res2+res3));
	
	
	return complex<T>(0,2)*(res1+res2+res3);
	
//	return complex<T>(0,-1)*(
//		-_HIGGS_SCALAR_SCALAR_PARAM_X*T(2)*(
//			(mass2)*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(T(2)*S(i3,i2))
//			-(mass2)*(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())/(T(2)*S(i1,i2))        
//			)/(ep.p(i2)->Lt()*ep.ref()->Lt())
//		-_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->Lt()*ep.ref()->Lt())))*(
//			(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())
//			-(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.p(i2)->Lt()))			 
//		);

}

//Hsps
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g2H_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i1));
	const complex<T> mass=eval_param<T>::mass(masses.p(i1));

	const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();	
	const complex<T> mass1=eval_param<T>::mass(masses.p(i1));
	const complex<T> mass3=eval_param<T>::mass(masses.p(i3));
//	const momentum<complex<T> > e01(ep.mom(i1)*(T(1)/(mass1))-mass1/(ep.mom(i1)*ep.ref()->P())*ep.ref()->P());
//	const momentum<complex<T> > e03(ep.mom(i3)*(T(1)/(mass3))-mass3/(ep.mom(i3)*ep.ref()->P())*ep.ref()->P());
	const momentum<complex<T> > e01(ep.mom(i1));
	const momentum<complex<T> > e03(ep.mom(i3));
	const momentum<complex<T> > e02=PfLLt(ep.p(i2)->Lt(),ep.ref()->L())/(ep.p(i2)->L()*ep.ref()->L());
	
//	_MESSAGE6("Bits:",(e01*e02)*((ep.mom(i1)-ep.mom(i2))*e03)
//					  ,", ",(e02*e03)*((ep.mom(i2)-ep.mom(i3))*e01)
//					  ,", ",(e03*e01)*((ep.mom(i3)-ep.mom(i1))*e02));
	
	
//	return complex<T>(0,1)*(
//		-_HIGGS_SCALAR_SCALAR_PARAM_X*T(2)*(
//			 (mass2)*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//			 -(mass2)*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))        
//			 )/(ep.p(i2)->L()*ep.ref()->L())
//		-_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->L()*ep.ref()->L())))*(
//			 (ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())
//			 -(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.p(i2)->L()))
//		);
	
	momentum<complex<T> > eps1=(ep.mom(i1)/mass)-(mass*ep.ref()->P())/(ep.mom(i1)*ep.ref()->P());
	momentum<complex<T> > eps2=-PfLLt(ep.p(i2)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i2)->L());
	momentum<complex<T> > eps3=(ep.mom(i3)/mass)-(mass*ep.ref()->P())/(ep.mom(i3)*ep.ref()->P());
	
	//SHS-SSg
	momentum<complex<T> > propleg1=-(ep.mom(i2)+ep.mom(i3));
	momentum<complex<T> > propeps1=(propleg1/mass)-(mass*ep.ref()->P())/(propleg1*ep.ref()->P());
	complex<T> higgsvert1=_HIGGS_SCALAR_SCALAR_PARAM*(propeps1*((propleg1*ep.mom(i1))*eps1-(eps1*propleg1)*ep.mom(i1)))/(propleg1.square()-mass2);
	complex<T> res1(higgsvert1*((propeps1*eps2)*((propleg1-ep.mom(i2))*eps3)
								+(eps2*eps3)*((ep.mom(i2)-ep.mom(i3))*propeps1)
								+(eps3*propeps1)*((ep.mom(i3)-propleg1)*eps2)));
	//gHg-gSS
	momentum<complex<T> > propleg2=-(ep.mom(i1)+ep.mom(i3));
	momentum<complex<T> > higgsvert2=((propleg2*ep.mom(i2))*eps2-(eps2*propleg2)*ep.mom(i2))/propleg2.square();
	complex<T> res2new((eps1*higgsvert2)*((ep.mom(i1)-propleg2)*eps3)
					   +(higgsvert2*eps3)*((propleg2-ep.mom(i3))*eps1)
					   +(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*higgsvert2));
//	_MESSAGE4("res2 bits:",complex<T>(0,1)*((eps1*eps2)*((ep.mom(i1)-propleg2)*eps3)
//											+(eps2*eps3)*((propleg2-ep.mom(i3))*eps1)
//											+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2))
//								,"=?=",complex<T>(0,-1)*((ep.mom(i1)-ep.mom(i3))*eps2));
	complex<T> res2(_HIGGS_GLUON_GLUON_PARAM*((ep.mom(i1)-ep.mom(i3))*higgsvert2));
	//SHS-SSg
	momentum<complex<T> > propleg3=-(ep.mom(i1)+ep.mom(i2));
	momentum<complex<T> > propeps3=(propleg3/mass)-(mass*ep.ref()->P())/(propleg3*ep.ref()->P());
	complex<T> higgsvert3=_HIGGS_SCALAR_SCALAR_PARAM*(propeps3*((propleg3*ep.mom(i3))*eps3-(eps3*propleg3)*ep.mom(i3)))/(propleg3.square()-mass2);
	complex<T> res3(higgsvert3*((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*propeps3)
								+(eps2*propeps3)*((ep.mom(i2)-propleg3)*eps1)
								+(propeps3*eps1)*((propleg3-ep.mom(i1))*eps2)));
	//SHgS
	complex<T> res4((eps1*eps2)*((ep.mom(i1)-ep.mom(i2))*eps3)
					+(eps2*eps3)*((ep.mom(i2)-ep.mom(i2))*eps1)
					+(eps3*eps1)*((ep.mom(i3)-ep.mom(i1))*eps2));

//	_PRINT(complex<T>(0,1)*(res1+res3));
//	_MESSAGE2(" res1+res3 simplify to:",complex<T>(0,-2)*((ep.mom(i3)*eps2)*(mass2)/(propleg1.square()-mass2)
//										   -(ep.mom(i1)*eps2)*(mass2)/(propleg3.square()-mass2)));

//	_PRINT(complex<T>(0,1)*res2);
//	_PRINT(complex<T>(0,-1)*(ep.mom(i1)-ep.mom(i3))*higgsvert2);
//	_PRINT(complex<T>(0,1)*res4);
//	_MESSAGE4(" OLD phdSPS:",complex<T>(0,-1)*(
//											  T(1/*0.5*/)*_HIGGS_SCALAR_SCALAR_PARAM*(
//																					  (BIT0*mass2+(BIT1*(m23*ep.p(i1)->P())
//																								   -BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																								   /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))))
//																					  *(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//																					  -(BIT0*mass2+(BIT1*(m12*ep.p(i3)->P())
//																									-BIT2*mass2*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//																									/((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))))
//																					  *(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))      
//																					  )/(ep.p(i2)->L()*ep.ref()->L())
//											  +_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->L()*ep.ref()->L())))*(
//																											(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())
//																											-(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.p(i2)->L()))
//											  ),
//			  " NEW:",complex<T>(0,1)*(res1+res2+res3));
	
	
	return complex<T>(0,2)*(res1+res2+res3);
}
    
    
//phsg+G-
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
        
    return complex<T>(0,0);
}
    
//phsg+G-
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
    
    return complex<T>(0,0);
}

//phsg+G+
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g7_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
    
    return complex<T>(0,0);
}

//phsg-G+
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g8_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
    
    return complex<T>(0,0);
}
    
    
    
//phdsg+G-
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g9_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
    
    return complex<T>(0,0);
}

//phdsg+G-
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g10_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
    
    return complex<T>(0,0);
}

//phdsg+G+
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g11_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
    
    return complex<T>(0,0);
}

//phdsg-G+
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1g12_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
    
    return complex<T>(0,0);
}
    


template <class T> complex<T> (*A1ph2sc1g_Tree_Ptr_eval(long long hc))(const eval_param<T>&, const mass_param_coll&) {
//		cout << hex << hc << dec << endl;
	switch (hc) {
			
			//
			//ph
			//
			
		case 0x9404://phssg
			return &A1ph2sc1g1_eval<0,1,2,3>;
			break;
		case 0x4904://sphsg
			return &A1ph2sc1g1_eval<1,0,2,3>;
			break;
		case 0x4094://ssphg
			return &A1ph2sc1g1_eval<2,0,1,3>;
			break;
		case 0x4049://ssgph
			return &A1ph2sc1g1_eval<3,0,1,2>;
			break;
			
		case 0x9414://phssg
			return &A1ph2sc1g2_eval<0,1,2,3>;
			break;
		case 0x4914://sphsg
			return &A1ph2sc1g2_eval<1,0,2,3>;
			break;
		case 0x4194://ssphg
			return &A1ph2sc1g2_eval<2,0,1,3>;
			break;
		case 0x4149://ssgph
			return &A1ph2sc1g2_eval<3,0,1,2>;
			break;
            
        case 0x9440://phssg
			return &A1ph2sc1g1_eval<0,2,3,1>;
			break;
		case 0x4940://sphsg
			return &A1ph2sc1g1_eval<1,2,3,0>;
			break;
		case 0x4490://ssphg
			return &A1ph2sc1g1_eval<2,1,3,0>;
			break;
		case 0x4409://ssgph
			return &A1ph2sc1g1_eval<3,1,2,0>;
			break;
			
		case 0x9441://phssg
			return &A1ph2sc1g2_eval<0,2,3,1>;
			break;
		case 0x4941://sphsg
			return &A1ph2sc1g2_eval<1,2,3,0>;
			break;
		case 0x4491://ssphg
			return &A1ph2sc1g2_eval<2,1,3,0>;
			break;
		case 0x4419://ssgph
			return &A1ph2sc1g2_eval<3,1,2,0>;
			break;
            
			
			
		case 0x9044://phssg
			return &A1ph2sc1g1_eval<0,3,1,2>;
			break;
		case 0x0944://sphsg
			return &A1ph2sc1g1_eval<1,3,0,2>;
			break;
		case 0x0494://ssphg
			return &A1ph2sc1g1_eval<2,3,0,1>;
			break;
		case 0x0449://ssgph
			return &A1ph2sc1g1_eval<3,2,0,1>;
			break;
			
		case 0x9144://phssg
			return &A1ph2sc1g2_eval<0,3,1,2>;
			break;
		case 0x1944://sphsg
			return &A1ph2sc1g2_eval<1,3,0,2>;
			break;
		case 0x1494://ssphg
			return &A1ph2sc1g2_eval<2,3,0,1>;
			break;
		case 0x1449://ssgph
			return &A1ph2sc1g2_eval<3,2,0,1>;
			break;
            
//            //
//            // gmM, gpM and gsc ph amplitudes
//            //
//            
//  		case 0x94C0://phsG+g-
//			return &A1ph2sc1g5_eval<0,1,2,3>;
//			break;
//		case 0x49C0://sphG+g-
//			return &A1ph2sc1g5_eval<1,0,2,3>;
//			break;
//		case 0x4C90://G+sphg-
//			return &A1ph2sc1g5_eval<2,0,1,3>;
//			break;
//		case 0x4C09://G+sg-ph
//			return &A1ph2sc1g5_eval<3,0,1,2>;
//			break;
//            
//        case 0x904C://phg-sG+
//			return &A1ph2sc1g5_eval<0,2,3,1>;
//			break;
//		case 0x094C://g-phsG+
//			return &A1ph2sc1g5_eval<1,2,3,0>;
//			break;
//		case 0x049C://g-sphG+
//			return &A1ph2sc1g5_eval<2,1,3,0>;
//			break;
//		case 0x04C9://g-sG+ph
//			return &A1ph2sc1g5_eval<3,1,2,0>;
//			break;
//            
//        case 0x940C://phG+g-s
//			return &A1ph2sc1g5_eval<0,3,1,2>;
//			break;
//		case 0x490C://G+phg-s
//			return &A1ph2sc1g5_eval<1,3,0,2>;
//			break;
//		case 0x409C://G+g-phs
//			return &A1ph2sc1g5_eval<2,3,0,1>;
//			break;
//		case 0x40C9://G+g-sph
//			return &A1ph2sc1g5_eval<3,2,0,1>;
//			break;
//            
//            
//        case 0x94C1://phsG+g+
//			return &A1ph2sc1g6_eval<0,1,2,3>;
//			break;
//		case 0x49C1://sphG+g+
//			return &A1ph2sc1g6_eval<1,0,2,3>;
//			break;
//		case 0x4C91://G+sphg+
//			return &A1ph2sc1g6_eval<2,0,1,3>;
//			break;
//		case 0x4C19://G+sg+ph
//			return &A1ph2sc1g6_eval<3,0,1,2>;
//			break;
//            
//        case 0x914C://phg+sG+
//			return &A1ph2sc1g6_eval<0,2,3,1>;
//			break;
//		case 0x194C://g+phsG+
//			return &A1ph2sc1g6_eval<1,2,3,0>;
//			break;
//		case 0x149C://g+sphG+
//			return &A1ph2sc1g6_eval<2,1,3,0>;
//			break;
//		case 0x14C9://g+sG+ph
//			return &A1ph2sc1g6_eval<3,1,2,0>;
//			break;
//            
//        case 0x941C://phG+g+s
//			return &A1ph2sc1g6_eval<0,3,1,2>;
//			break;
//		case 0x491C://G+phg+s
//			return &A1ph2sc1g6_eval<1,3,0,2>;
//			break;
//		case 0x419C://G+g+phs
//			return &A1ph2sc1g6_eval<2,3,0,1>;
//			break;
//		case 0x41C9://G+g+sph
//			return &A1ph2sc1g6_eval<3,2,0,1>;
//			break;
//
//			
//          
//        case 0x94B0://phsG-g-
//			return &A1ph2sc1g7_eval<0,1,2,3>;
//			break;
//		case 0x49B0://sphG-g-
//			return &A1ph2sc1g7_eval<1,0,2,3>;
//			break;
//		case 0x4B90://G-sphg-
//			return &A1ph2sc1g7_eval<2,0,1,3>;
//			break;
//		case 0x4B09://G-sg-ph
//			return &A1ph2sc1g7_eval<3,0,1,2>;
//			break;
//            
//        case 0x904B://phg-sG-
//			return &A1ph2sc1g7_eval<0,2,3,1>;
//			break;
//		case 0x094B://g-phsG-
//			return &A1ph2sc1g7_eval<1,2,3,0>;
//			break;
//		case 0x049B://g-sphG-
//			return &A1ph2sc1g7_eval<2,1,3,0>;
//			break;
//		case 0x04B9://g-sG-ph
//			return &A1ph2sc1g7_eval<3,1,2,0>;
//			break;
//            
//        case 0x940B://phG-g-s
//			return &A1ph2sc1g7_eval<0,3,1,2>;
//			break;
//		case 0x490B://G-phg-s
//			return &A1ph2sc1g7_eval<1,3,0,2>;
//			break;
//		case 0x409B://G-g-phs
//			return &A1ph2sc1g7_eval<2,3,0,1>;
//			break;
//		case 0x40B9://G-g-sph
//			return &A1ph2sc1g7_eval<3,2,0,1>;
//			break;
//            
//        
//        case 0x94B1://phsG-g+
//			return &A1ph2sc1g8_eval<0,1,2,3>;
//			break;
//		case 0x49B1://sphG-g+
//			return &A1ph2sc1g8_eval<1,0,2,3>;
//			break;
//		case 0x4B91://G-sphg+
//			return &A1ph2sc1g8_eval<2,0,1,3>;
//			break;
//		case 0x4B19://G-sg+ph
//			return &A1ph2sc1g8_eval<3,0,1,2>;
//			break;
//            
//        case 0x914B://phg+sG-
//			return &A1ph2sc1g8_eval<0,2,3,1>;
//			break;
//		case 0x194B://g+phsG-
//			return &A1ph2sc1g8_eval<1,2,3,0>;
//			break;
//		case 0x149B://g+sphG-
//			return &A1ph2sc1g8_eval<2,1,3,0>;
//			break;
//		case 0x14B9://g+sG-ph
//			return &A1ph2sc1g8_eval<3,1,2,0>;
//			break;
//            
//        case 0x941B://phG-g+s
//			return &A1ph2sc1g8_eval<0,3,1,2>;
//			break;
//		case 0x491B://G-phg+s
//			return &A1ph2sc1g8_eval<1,3,0,2>;
//			break;
//		case 0x419B://G-g+phs
//			return &A1ph2sc1g8_eval<2,3,0,1>;
//			break;
//		case 0x41B9://G-g+sph
//			return &A1ph2sc1g8_eval<3,2,0,1>;
//			break;
//            
            
            
            
//        case 0x94B1://phG-G-g+
//			return &A1ph2sc1g8_eval<0,1,2,3>;
//			break;
//		case 0x49B1://sphG-g+
//			return &A1ph2sc1g8_eval<1,0,2,3>;
//			break;
//		case 0x4B91://G-sphg+
//			return &A1ph2sc1g8_eval<2,0,1,3>;
//			break;
//		case 0x4B19://G-sg+ph
//			return &A1ph2sc1g8_eval<3,0,1,2>;
//			break;

			
			
			//
			// phd
			//
			
		case 0xA440://phdssg
			return &A1ph2sc1g4_eval<0,2,3,1>;
			break;
		case 0x4A40://spdhsg
			return &A1ph2sc1g4_eval<1,2,3,0>;
			break;
		case 0x44A0://ssphdg
			return &A1ph2sc1g4_eval<2,1,3,0>;
			break;
		case 0x440A://ssgphd
			return &A1ph2sc1g4_eval<3,1,2,0>;
			break;
			
		case 0xA441://phdssg
			return &A1ph2sc1g3_eval<0,2,3,1>;
			break;
		case 0x4A41://sphdsg
			return &A1ph2sc1g3_eval<1,2,3,0>;
			break;
		case 0x44A1://ssphdg
			return &A1ph2sc1g3_eval<2,1,3,0>;
			break;
		case 0x441A://ssgphd
			return &A1ph2sc1g3_eval<3,1,2,0>;
			break;
			
			
		case 0xA404://phdssg
			return &A1ph2sc1g4_eval<0,1,2,3>;
			break;
		case 0x4A04://sphdsg
			return &A1ph2sc1g4_eval<1,0,2,3>;
			break;
		case 0x40A4://ssphdg
			return &A1ph2sc1g4_eval<2,0,1,3>;
			break;
		case 0x404A://ssgphd
			return &A1ph2sc1g4_eval<3,0,1,2>;
			break;
			
		case 0xA414://phdssg
			return &A1ph2sc1g3_eval<0,1,2,3>;
			break;
		case 0x4A14://sphdsg
			return &A1ph2sc1g3_eval<1,0,2,3>;
			break;
		case 0x41A4://ssphdg
			return &A1ph2sc1g3_eval<2,0,1,3>;
			break;
		case 0x414A://ssgphd
			return &A1ph2sc1g3_eval<3,0,1,2>;
			break;
			
			
		case 0xA044://phdssg
			return &A1ph2sc1g4_eval<0,3,1,2>;
			break;
		case 0x0A44://sphdsg
			return &A1ph2sc1g4_eval<1,3,0,2>;
			break;
		case 0x04A4://ssphdg
			return &A1ph2sc1g4_eval<2,3,0,1>;
			break;
		case 0x044A://ssgphd
			return &A1ph2sc1g4_eval<3,2,0,1>;
			break;
			
		case 0xA144://phdssg
			return &A1ph2sc1g3_eval<0,3,1,2>;
			break;
		case 0x1A44://sphdsg
			return &A1ph2sc1g3_eval<1,3,0,2>;
			break;
		case 0x14A4://ssphdg
			return &A1ph2sc1g3_eval<2,3,0,1>;
			break;
		case 0x144A://ssgpdh
			return &A1ph2sc1g3_eval<3,2,0,1>;
			break;
            
            
            
//            //
//            // gmM, gpM and gsc phd amplitudes
//            //
//            
//  		case 0xA4C0://phdsG+g-
//			return &A1ph2sc1g9_eval<0,1,2,3>;
//			break;
//		case 0x4AC0://sphdG+g-
//			return &A1ph2sc1g9_eval<1,0,2,3>;
//			break;
//		case 0x4CA0://G+sphdg-
//			return &A1ph2sc1g9_eval<2,0,1,3>;
//			break;
//		case 0x4C0A://G+sg-phd
//			return &A1ph2sc1g9_eval<3,0,1,2>;
//			break;
//            
//        case 0xA04C://phdg-sG+
//			return &A1ph2sc1g9_eval<0,2,3,1>;
//			break;
//		case 0x0A4C://g-phdsG+
//			return &A1ph2sc1g9_eval<1,2,3,0>;
//			break;
//		case 0x04AC://g-sphdG+
//			return &A1ph2sc1g9_eval<2,1,3,0>;
//			break;
//		case 0x04CA://g-sG+phd
//			return &A1ph2sc1g9_eval<3,1,2,0>;
//			break;
//            
//        case 0xA40C://phdG+g-s
//			return &A1ph2sc1g9_eval<0,3,1,2>;
//			break;
//		case 0x4A0C://G+phdg-s
//			return &A1ph2sc1g9_eval<1,3,0,2>;
//			break;
//		case 0x40AC://G+g-phds
//			return &A1ph2sc1g9_eval<2,3,0,1>;
//			break;
//		case 0x40CA://G+g-sphd
//			return &A1ph2sc1g9_eval<3,2,0,1>;
//			break;
//            
//            
//        case 0xA4C1://phdsG+g+
//			return &A1ph2sc1g10_eval<0,1,2,3>;
//			break;
//		case 0x4AC1://sphdG+g+
//			return &A1ph2sc1g10_eval<1,0,2,3>;
//			break;
//		case 0x4CA1://G+sphdg+
//			return &A1ph2sc1g10_eval<2,0,1,3>;
//			break;
//		case 0x4C1A://G+sg+phd
//			return &A1ph2sc1g10_eval<3,0,1,2>;
//			break;
//            
//        case 0xA14C://phdg+sG+
//			return &A1ph2sc1g10_eval<0,2,3,1>;
//			break;
//		case 0x1A4C://g+phdsG+
//			return &A1ph2sc1g10_eval<1,2,3,0>;
//			break;
//		case 0x14AC://g+sphdG+
//			return &A1ph2sc1g10_eval<2,1,3,0>;
//			break;
//		case 0x14CA://g+sG+phd
//			return &A1ph2sc1g10_eval<3,1,2,0>;
//			break;
//            
//        case 0xA41C://phdG+g+s
//			return &A1ph2sc1g10_eval<0,3,1,2>;
//			break;
//		case 0x4A1C://G+phdg+s
//			return &A1ph2sc1g10_eval<1,3,0,2>;
//			break;
//		case 0x41AC://G+g+phds
//			return &A1ph2sc1g10_eval<2,3,0,1>;
//			break;
//		case 0x41CA://G+g+sphd
//			return &A1ph2sc1g10_eval<3,2,0,1>;
//			break;
//            
//			
//            
//        case 0xA4B0://phdsG-g-
//			return &A1ph2sc1g11_eval<0,1,2,3>;
//			break;
//		case 0x4AB0://sphdG-g-
//			return &A1ph2sc1g11_eval<1,0,2,3>;
//			break;
//		case 0x4BA0://G-sphdg-
//			return &A1ph2sc1g11_eval<2,0,1,3>;
//			break;
//		case 0x4B0A://G-sg-phd
//			return &A1ph2sc1g11_eval<3,0,1,2>;
//			break;
//            
//        case 0xA04B://phdg-sG-
//			return &A1ph2sc1g11_eval<0,2,3,1>;
//			break;
//		case 0x0A4B://g-phdsG-
//			return &A1ph2sc1g11_eval<1,2,3,0>;
//			break;
//		case 0x04AB://g-sphdG-
//			return &A1ph2sc1g11_eval<2,1,3,0>;
//			break;
//		case 0x04BA://g-sG-phd
//			return &A1ph2sc1g11_eval<3,1,2,0>;
//			break;
//            
//        case 0xA40B://phdG-g-s
//			return &A1ph2sc1g11_eval<0,3,1,2>;
//			break;
//		case 0x4A0B://G-phdg-s
//			return &A1ph2sc1g11_eval<1,3,0,2>;
//			break;
//		case 0x40AB://G-g-phds
//			return &A1ph2sc1g11_eval<2,3,0,1>;
//			break;
//		case 0x40BA://G-g-sphd
//			return &A1ph2sc1g11_eval<3,2,0,1>;
//			break;
//            
//            
//        case 0xA4B1://phdsG-g+
//			return &A1ph2sc1g12_eval<0,1,2,3>;
//			break;
//		case 0x4AB1://sphdG-g+
//			return &A1ph2sc1g12_eval<1,0,2,3>;
//			break;
//		case 0x4BA1://G-sphdg+
//			return &A1ph2sc1g12_eval<2,0,1,3>;
//			break;
//		case 0x4B1A://G-sg+phd
//			return &A1ph2sc1g12_eval<3,0,1,2>;
//			break;
//            
//        case 0xA14B://phdg+sG-
//			return &A1ph2sc1g12_eval<0,2,3,1>;
//			break;
//		case 0x1A4B://g+phdsG-
//			return &A1ph2sc1g12_eval<1,2,3,0>;
//			break;
//		case 0x14AB://g+sphdG-
//			return &A1ph2sc1g12_eval<2,1,3,0>;
//			break;
//		case 0x14BA://g+sG-phd
//			return &A1ph2sc1g12_eval<3,1,2,0>;
//			break;
//            
//        case 0xA41B://phdG-g+s
//			return &A1ph2sc1g12_eval<0,3,1,2>;
//			break;
//		case 0x4A1B://G-phdg+s
//			return &A1ph2sc1g12_eval<1,3,0,2>;
//			break;
//		case 0x41AB://G-g+phds
//			return &A1ph2sc1g12_eval<2,3,0,1>;
//			break;
//		case 0x41BA://G-g+sphd
//			return &A1ph2sc1g12_eval<3,2,0,1>;
//			break;

			
			
			
		//
		//H
		//
		
	case 0xE404://phssg
		return &A1ph2sc1g1H_eval<0,1,2,3>;
		break;
	case 0x4E04://sphsg
		return &A1ph2sc1g1H_eval<1,0,2,3>;
		break;
	case 0x40E4://ssphg
		return &A1ph2sc1g1H_eval<2,0,1,3>;
		break;
	case 0x404E://ssgph
		return &A1ph2sc1g1H_eval<3,0,1,2>;
		break;
		
	case 0xE414://phssg
		return &A1ph2sc1g2H_eval<0,1,2,3>;
		break;
	case 0x4E14://sphsg
		return &A1ph2sc1g2H_eval<1,0,2,3>;
		break;
	case 0x41E4://ssphg
		return &A1ph2sc1g2H_eval<2,0,1,3>;
		break;
	case 0x414E://ssgph
		return &A1ph2sc1g2H_eval<3,0,1,2>;
		break;
		
	case 0xE440://phssg
		return &A1ph2sc1g1H_eval<0,2,3,1>;
		break;
	case 0x4E40://sphsg
		return &A1ph2sc1g1H_eval<1,2,3,0>;
		break;
	case 0x44E0://ssphg
		return &A1ph2sc1g1H_eval<2,1,3,0>;
		break;
	case 0x440E://ssgph
		return &A1ph2sc1g1H_eval<3,1,2,0>;
		break;
		
	case 0xE441://phssg
		return &A1ph2sc1g2H_eval<0,2,3,1>;
		break;
	case 0x4E41://sphsg
		return &A1ph2sc1g2H_eval<1,2,3,0>;
		break;
	case 0x44E1://ssphg
		return &A1ph2sc1g2H_eval<2,1,3,0>;
		break;
	case 0x441E://ssgph
		return &A1ph2sc1g2H_eval<3,1,2,0>;
		break;
		
		
		
	case 0xE044://phssg
		return &A1ph2sc1g1H_eval<0,3,1,2>;
		break;
	case 0x0E44://sphsg
		return &A1ph2sc1g1H_eval<1,3,0,2>;
		break;
	case 0x04E4://ssphg
		return &A1ph2sc1g1H_eval<2,3,0,1>;
		break;
	case 0x044E://ssgph
		return &A1ph2sc1g1H_eval<3,2,0,1>;
		break;
		
	case 0xE144://phssg
		return &A1ph2sc1g2H_eval<0,3,1,2>;
		break;
	case 0x1E44://sphsg
		return &A1ph2sc1g2H_eval<1,3,0,2>;
		break;
	case 0x14E4://ssphg
		return &A1ph2sc1g2H_eval<2,3,0,1>;
		break;
	case 0x144E://ssgph
		return &A1ph2sc1g2H_eval<3,2,0,1>;
		break;
		
		//            //
		//            // gmM, gpM and gsc ph amplitudes
		//            //
		//            
		//  		case 0xE4C0://phsG+g-
		//			return &A1ph2sc1g5_eval<0,1,2,3>;
		//			break;
		//		case 0x4EC0://sphG+g-
		//			return &A1ph2sc1g5_eval<1,0,2,3>;
		//			break;
		//		case 0x4CE0://G+sphg-
		//			return &A1ph2sc1g5_eval<2,0,1,3>;
		//			break;
		//		case 0x4C0E://G+sg-ph
		//			return &A1ph2sc1g5_eval<3,0,1,2>;
		//			break;
		//            
		//        case 0xE04C://phg-sG+
		//			return &A1ph2sc1g5_eval<0,2,3,1>;
		//			break;
		//		case 0x0E4C://g-phsG+
		//			return &A1ph2sc1g5_eval<1,2,3,0>;
		//			break;
		//		case 0x04EC://g-sphG+
		//			return &A1ph2sc1g5_eval<2,1,3,0>;
		//			break;
		//		case 0x04CE://g-sG+ph
		//			return &A1ph2sc1g5_eval<3,1,2,0>;
		//			break;
		//            
		//        case 0xE40C://phG+g-s
		//			return &A1ph2sc1g5_eval<0,3,1,2>;
		//			break;
		//		case 0x4E0C://G+phg-s
		//			return &A1ph2sc1g5_eval<1,3,0,2>;
		//			break;
		//		case 0x40EC://G+g-phs
		//			return &A1ph2sc1g5_eval<2,3,0,1>;
		//			break;
		//		case 0x40CE://G+g-sph
		//			return &A1ph2sc1g5_eval<3,2,0,1>;
		//			break;
		//            
		//            
		//        case 0xE4C1://phsG+g+
		//			return &A1ph2sc1g6_eval<0,1,2,3>;
		//			break;
		//		case 0x4EC1://sphG+g+
		//			return &A1ph2sc1g6_eval<1,0,2,3>;
		//			break;
		//		case 0x4CE1://G+sphg+
		//			return &A1ph2sc1g6_eval<2,0,1,3>;
		//			break;
		//		case 0x4C1E://G+sg+ph
		//			return &A1ph2sc1g6_eval<3,0,1,2>;
		//			break;
		//            
		//        case 0xE14C://phg+sG+
		//			return &A1ph2sc1g6_eval<0,2,3,1>;
		//			break;
		//		case 0x1E4C://g+phsG+
		//			return &A1ph2sc1g6_eval<1,2,3,0>;
		//			break;
		//		case 0x14EC://g+sphG+
		//			return &A1ph2sc1g6_eval<2,1,3,0>;
		//			break;
		//		case 0x14CE://g+sG+ph
		//			return &A1ph2sc1g6_eval<3,1,2,0>;
		//			break;
		//            
		//        case 0xE41C://phG+g+s
		//			return &A1ph2sc1g6_eval<0,3,1,2>;
		//			break;
		//		case 0x4E1C://G+phg+s
		//			return &A1ph2sc1g6_eval<1,3,0,2>;
		//			break;
		//		case 0x41EC://G+g+phs
		//			return &A1ph2sc1g6_eval<2,3,0,1>;
		//			break;
		//		case 0x41CE://G+g+sph
		//			return &A1ph2sc1g6_eval<3,2,0,1>;
		//			break;
		//
		//			
		//          
		//        case 0xE4B0://phsG-g-
		//			return &A1ph2sc1g7_eval<0,1,2,3>;
		//			break;
		//		case 0x4EB0://sphG-g-
		//			return &A1ph2sc1g7_eval<1,0,2,3>;
		//			break;
		//		case 0x4BE0://G-sphg-
		//			return &A1ph2sc1g7_eval<2,0,1,3>;
		//			break;
		//		case 0x4B0E://G-sg-ph
		//			return &A1ph2sc1g7_eval<3,0,1,2>;
		//			break;
		//            
		//        case 0xE04B://phg-sG-
		//			return &A1ph2sc1g7_eval<0,2,3,1>;
		//			break;
		//		case 0x0E4B://g-phsG-
		//			return &A1ph2sc1g7_eval<1,2,3,0>;
		//			break;
		//		case 0x04EB://g-sphG-
		//			return &A1ph2sc1g7_eval<2,1,3,0>;
		//			break;
		//		case 0x04BE://g-sG-ph
		//			return &A1ph2sc1g7_eval<3,1,2,0>;
		//			break;
		//            
		//        case 0xE40B://phG-g-s
		//			return &A1ph2sc1g7_eval<0,3,1,2>;
		//			break;
		//		case 0x4E0B://G-phg-s
		//			return &A1ph2sc1g7_eval<1,3,0,2>;
		//			break;
		//		case 0x40EB://G-g-phs
		//			return &A1ph2sc1g7_eval<2,3,0,1>;
		//			break;
		//		case 0x40BE://G-g-sph
		//			return &A1ph2sc1g7_eval<3,2,0,1>;
		//			break;
		//            
		//        
		//        case 0xE4B1://phsG-g+
		//			return &A1ph2sc1g8_eval<0,1,2,3>;
		//			break;
		//		case 0x4EB1://sphG-g+
		//			return &A1ph2sc1g8_eval<1,0,2,3>;
		//			break;
		//		case 0x4BE1://G-sphg+
		//			return &A1ph2sc1g8_eval<2,0,1,3>;
		//			break;
		//		case 0x4B1E://G-sg+ph
		//			return &A1ph2sc1g8_eval<3,0,1,2>;
		//			break;
		//            
		//        case 0xE14B://phg+sG-
		//			return &A1ph2sc1g8_eval<0,2,3,1>;
		//			break;
		//		case 0x1E4B://g+phsG-
		//			return &A1ph2sc1g8_eval<1,2,3,0>;
		//			break;
		//		case 0x14EB://g+sphG-
		//			return &A1ph2sc1g8_eval<2,1,3,0>;
		//			break;
		//		case 0x14BE://g+sG-ph
		//			return &A1ph2sc1g8_eval<3,1,2,0>;
		//			break;
		//            
		//        case 0xE41B://phG-g+s
		//			return &A1ph2sc1g8_eval<0,3,1,2>;
		//			break;
		//		case 0x4E1B://G-phg+s
		//			return &A1ph2sc1g8_eval<1,3,0,2>;
		//			break;
		//		case 0x41EB://G-g+phs
		//			return &A1ph2sc1g8_eval<2,3,0,1>;
		//			break;
		//		case 0x41BE://G-g+sph
		//			return &A1ph2sc1g8_eval<3,2,0,1>;
		//			break;
		//            
		
		
		
		//        case 0xE4B1://phG-G-g+
		//			return &A1ph2sc1g8_eval<0,1,2,3>;
		//			break;
		//		case 0x4EB1://sphG-g+
		//			return &A1ph2sc1g8_eval<1,0,2,3>;
		//			break;
		//		case 0x4BE1://G-sphg+
		//			return &A1ph2sc1g8_eval<2,0,1,3>;
		//			break;
		//		case 0x4B1E://G-sg+ph
		//			return &A1ph2sc1g8_eval<3,0,1,2>;
		//			break;
		
		
//	case 0xE440://phdssg
//		return &A1ph2sc1g4_eval<0,2,3,1>;
//		break;
//	case 0x4E40://spdhsg
//		return &A1ph2sc1g4_eval<1,2,3,0>;
//		break;
//	case 0x44E0://ssphdg
//		return &A1ph2sc1g4_eval<2,1,3,0>;
//		break;
//	case 0x440E://ssgphd
//		return &A1ph2sc1g4_eval<3,1,2,0>;
//		break;
//		
//	case 0xE441://phdssg
//		return &A1ph2sc1g3_eval<0,2,3,1>;
//		break;
//	case 0x4E41://sphdsg
//		return &A1ph2sc1g3_eval<1,2,3,0>;
//		break;
//	case 0x44E1://ssphdg
//		return &A1ph2sc1g3_eval<2,1,3,0>;
//		break;
//	case 0x441E://ssgphd
//		return &A1ph2sc1g3_eval<3,1,2,0>;
//		break;
//		
//		
//	case 0xE404://phdssg
//		return &A1ph2sc1g4_eval<0,1,2,3>;
//		break;
//	case 0x4E04://sphdsg
//		return &A1ph2sc1g4_eval<1,0,2,3>;
//		break;
//	case 0x40E4://ssphdg
//		return &A1ph2sc1g4_eval<2,0,1,3>;
//		break;
//	case 0x404E://ssgphd
//		return &A1ph2sc1g4_eval<3,0,1,2>;
//		break;
//		
//	case 0xE414://phdssg
//		return &A1ph2sc1g3_eval<0,1,2,3>;
//		break;
//	case 0x4E14://sphdsg
//		return &A1ph2sc1g3_eval<1,0,2,3>;
//		break;
//	case 0x41E4://ssphdg
//		return &A1ph2sc1g3_eval<2,0,1,3>;
//		break;
//	case 0x414E://ssgphd
//		return &A1ph2sc1g3_eval<3,0,1,2>;
//		break;
//		
//		
//	case 0xE044://phdssg
//		return &A1ph2sc1g4_eval<0,3,1,2>;
//		break;
//	case 0x0E44://sphdsg
//		return &A1ph2sc1g4_eval<1,3,0,2>;
//		break;
//	case 0x04E4://ssphdg
//		return &A1ph2sc1g4_eval<2,3,0,1>;
//		break;
//	case 0x044E://ssgphd
//		return &A1ph2sc1g4_eval<3,2,0,1>;
//		break;
//		
//	case 0xE144://phdssg
//		return &A1ph2sc1g3_eval<0,3,1,2>;
//		break;
//	case 0x1E44://sphdsg
//		return &A1ph2sc1g3_eval<1,3,0,2>;
//		break;
//	case 0x14E4://ssphdg
//		return &A1ph2sc1g3_eval<2,3,0,1>;
//		break;
//	case 0x144E://ssgpdh
//		return &A1ph2sc1g3_eval<3,2,0,1>;
//		break;
		
		
		
		//            //
		//            // gmM, gpM and gsc phd amplitudes
		//            //
		//            
		//  		case 0xA4C0://phdsG+g-
		//			return &A1ph2sc1gE_eval<0,1,2,3>;
		//			break;
		//		case 0x4AC0://sphdG+g-
		//			return &A1ph2sc1gE_eval<1,0,2,3>;
		//			break;
		//		case 0x4CA0://G+sphdg-
		//			return &A1ph2sc1gE_eval<2,0,1,3>;
		//			break;
		//		case 0x4C0A://G+sg-phd
		//			return &A1ph2sc1gE_eval<3,0,1,2>;
		//			break;
		//            
		//        case 0xA04C://phdg-sG+
		//			return &A1ph2sc1g9_eval<0,2,3,1>;
		//			break;
		//		case 0x0A4C://g-phdsG+
		//			return &A1ph2sc1g9_eval<1,2,3,0>;
		//			break;
		//		case 0x04AC://g-sphdG+
		//			return &A1ph2sc1g9_eval<2,1,3,0>;
		//			break;
		//		case 0x04CA://g-sG+phd
		//			return &A1ph2sc1g9_eval<3,1,2,0>;
		//			break;
		//            
		//        case 0xA40C://phdG+g-s
		//			return &A1ph2sc1g9_eval<0,3,1,2>;
		//			break;
		//		case 0x4A0C://G+phdg-s
		//			return &A1ph2sc1g9_eval<1,3,0,2>;
		//			break;
		//		case 0x40AC://G+g-phds
		//			return &A1ph2sc1g9_eval<2,3,0,1>;
		//			break;
		//		case 0x40CA://G+g-sphd
		//			return &A1ph2sc1g9_eval<3,2,0,1>;
		//			break;
		//            
		//            
		//        case 0xA4C1://phdsG+g+
		//			return &A1ph2sc1g10_eval<0,1,2,3>;
		//			break;
		//		case 0x4AC1://sphdG+g+
		//			return &A1ph2sc1g10_eval<1,0,2,3>;
		//			break;
		//		case 0x4CA1://G+sphdg+
		//			return &A1ph2sc1g10_eval<2,0,1,3>;
		//			break;
		//		case 0x4C1A://G+sg+phd
		//			return &A1ph2sc1g10_eval<3,0,1,2>;
		//			break;
		//            
		//        case 0xA14C://phdg+sG+
		//			return &A1ph2sc1g10_eval<0,2,3,1>;
		//			break;
		//		case 0x1A4C://g+phdsG+
		//			return &A1ph2sc1g10_eval<1,2,3,0>;
		//			break;
		//		case 0x14AC://g+sphdG+
		//			return &A1ph2sc1g10_eval<2,1,3,0>;
		//			break;
		//		case 0x14CA://g+sG+phd
		//			return &A1ph2sc1g10_eval<3,1,2,0>;
		//			break;
		//            
		//        case 0xA41C://phdG+g+s
		//			return &A1ph2sc1g10_eval<0,3,1,2>;
		//			break;
		//		case 0x4A1C://G+phdg+s
		//			return &A1ph2sc1g10_eval<1,3,0,2>;
		//			break;
		//		case 0x41AC://G+g+phds
		//			return &A1ph2sc1g10_eval<2,3,0,1>;
		//			break;
		//		case 0x41CA://G+g+sphd
		//			return &A1ph2sc1g10_eval<3,2,0,1>;
		//			break;
		//            
		//			
		//            
		//        case 0xA4B0://phdsG-g-
		//			return &A1ph2sc1g11_eval<0,1,2,3>;
		//			break;
		//		case 0x4AB0://sphdG-g-
		//			return &A1ph2sc1g11_eval<1,0,2,3>;
		//			break;
		//		case 0x4BA0://G-sphdg-
		//			return &A1ph2sc1g11_eval<2,0,1,3>;
		//			break;
		//		case 0x4B0A://G-sg-phd
		//			return &A1ph2sc1g11_eval<3,0,1,2>;
		//			break;
		//            
		//        case 0xA04B://phdg-sG-
		//			return &A1ph2sc1g11_eval<0,2,3,1>;
		//			break;
		//		case 0x0A4B://g-phdsG-
		//			return &A1ph2sc1g11_eval<1,2,3,0>;
		//			break;
		//		case 0x04AB://g-sphdG-
		//			return &A1ph2sc1g11_eval<2,1,3,0>;
		//			break;
		//		case 0x04BA://g-sG-phd
		//			return &A1ph2sc1g11_eval<3,1,2,0>;
		//			break;
		//            
		//        case 0xA40B://phdG-g-s
		//			return &A1ph2sc1g11_eval<0,3,1,2>;
		//			break;
		//		case 0x4A0B://G-phdg-s
		//			return &A1ph2sc1g11_eval<1,3,0,2>;
		//			break;
		//		case 0x40AB://G-g-phds
		//			return &A1ph2sc1g11_eval<2,3,0,1>;
		//			break;
		//		case 0x40BA://G-g-sphd
		//			return &A1ph2sc1g11_eval<3,2,0,1>;
		//			break;
		//            
		//            
		//        case 0xA4B1://phdsG-g+
		//			return &A1ph2sc1g12_eval<0,1,2,3>;
		//			break;
		//		case 0x4AB1://sphdG-g+
		//			return &A1ph2sc1g12_eval<1,0,2,3>;
		//			break;
		//		case 0x4BA1://G-sphdg+
		//			return &A1ph2sc1g12_eval<2,0,1,3>;
		//			break;
		//		case 0x4B1A://G-sg+phd
		//			return &A1ph2sc1g12_eval<3,0,1,2>;
		//			break;
		//            
		//        case 0xA14B://phdg+sG-
		//			return &A1ph2sc1g12_eval<0,2,3,1>;
		//			break;
		//		case 0x1A4B://g+phdsG-
		//			return &A1ph2sc1g12_eval<1,2,3,0>;
		//			break;
		//		case 0x14AB://g+sphdG-
		//			return &A1ph2sc1g12_eval<2,1,3,0>;
		//			break;
		//		case 0x14BA://g+sG-phd
		//			return &A1ph2sc1g12_eval<3,1,2,0>;
		//			break;
		//            
		//        case 0xA41B://phdG-g+s
		//			return &A1ph2sc1g12_eval<0,3,1,2>;
		//			break;
		//		case 0x4A1B://G-phdg+s
		//			return &A1ph2sc1g12_eval<1,3,0,2>;
		//			break;
		//		case 0x41AB://G-g+phds
		//			return &A1ph2sc1g12_eval<2,3,0,1>;
		//			break;
		//		case 0x41BA://G-g+sphd
		//			return &A1ph2sc1g12_eval<3,2,0,1>;
		//			break;
			

			
		default:// We return zero for all other helicity combinations
			cout << hex << "missing helcode A1ph2sc1g_Tree_Ptr_eval:" << hc << dec << endl; 	
			return 0;
	}
}

	
template complex<R> (*A1ph2sc_Tree_Ptr_eval(long long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2sc_Tree_Ptr_eval(long long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2sc_Tree_Ptr_eval(long long hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1ph2sc1g_Tree_Ptr_eval(long long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2sc1g_Tree_Ptr_eval(long long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2sc1g_Tree_Ptr_eval(long long hc))(const eval_param<RVHP>&, const mass_param_coll&);
			
#if BH_USE_GMP
template complex<RGMP> (*A1ph2sc_Tree_Ptr_eval(long long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2sc1g_Tree_Ptr_eval(long long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}	


