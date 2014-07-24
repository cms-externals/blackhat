/*
 *  A0_sc_massless.cpp
 *  BlackHat
 *
 *  Created by Darren Forde on 12/07/2010.
 *  Copyright 2010 BlackHat Collaboration. All rights reserved.
 *
 */

/* The 2 scalar amplitudes */

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;


//#define _HIGGS_SCALAR_SCALAR_PARAM_3 T(1)//T(12) // ph,sc,sc
//#define _HIGGS_SCALAR_SCALAR_PARAM T(1)//T(12) // ph,sc,sc
//#define _HIGGS_SCALAR_SCALAR_PARAM_D T(1)//T(12) // ph,sc,sc
//#define _HIGGS_GLUON_FAC T(1) // ph,s,g,s higgs coupling to the gluon
//#define BIT2 T(0)
//#define BIT1 T(0)
//#define BIT0 T(1)

#define _HIGGS_SCALAR_SCALAR_PARAM_1 T(12) // ph,sc,q,Q
#define _HIGGS_PARAM_2 T(12) // ph,sc,sc,q,q gluon prop
#define _HIGGS_PARAM_3 T(12) // ph,sc,sc,q,q massive quark prop
#define _OTHER_PARAM_2 T(1) // ph,sc,sc,q,q higgs attached to gluon prop


namespace BH {

/*
 *
 *
 * eval function that returns zero
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

template <int i0, int i1, int i2, class T> complex<T>  A2sc1g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
//	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
////	momentum<complex<T> > eps1=(ep.mom(i0)/mass)-(mass*ep.ref()->P())/(ep.mom(i0)*ep.ref()->P());
//	momentum<complex<T> > eps2=-PfLLt(ep.p(i1)->Lt(),ep.ref()->L())/(ep.ref()->L()*ep.p(i1)->L());
////	momentum<complex<T> > eps3=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
//	
//	momentum<complex<T> > ddim(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(1,0));
//	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
//	momentum<complex<T> > mom1D(mass*ddim);
//	momentum<complex<T> > mom2D(zero);
//	momentum<complex<T> > mom3D(mass*ddim);
//	
//	
//	momentum<complex<T> > eps1=(mass*ep.ref()->P())/(ep.mom(i0)*ep.ref()->P());
//	momentum<complex<T> > eps3=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
//	momentum<complex<T> > eps1D=ddim;
//	momentum<complex<T> > eps2D=zero;
//	momentum<complex<T> > eps3D=ddim;
//
//	
//	//SSg
//	complex<T> res1(complex<T>(0,-1)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
//		+((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
//		+((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
//	
//	_MESSAGE4("compare orig a:",complex<T>(0,-1)*(ep.ref()->L()*ep.p(i0)->Sm()*ep.p(i1)->Lt())/(ep.ref()->L()*ep.p(i1)->L()),"new:",res1);
	

	
    return complex<T>(0,-1)*(ep.ref()->L()*ep.p(i0)->Sm()*ep.p(i1)->Lt())/(ep.ref()->L()*ep.p(i1)->L());
}
	
template <int i0, int i1, int i2, class T> complex<T>  A2sc1g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
//	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
////	momentum<complex<T> > eps1=(ep.mom(i0)/mass)-(mass*ep.ref()->P())/(ep.mom(i0)*ep.ref()->P());
//	momentum<complex<T> > eps2=PfLLt(ep.ref()->Lt(),ep.p(i1)->L())/(ep.ref()->Lt()*ep.p(i1)->Lt());
////	momentum<complex<T> > eps3=(ep.mom(i2)/mass)-(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
//	
//	momentum<complex<T> > ddim(complex<T>(1,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
//	momentum<complex<T> > zero(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
//	momentum<complex<T> > mom1D(mass*ddim);
//	momentum<complex<T> > mom2D(zero);
//	momentum<complex<T> > mom3D(-mass*ddim);
//	
//	
//	momentum<complex<T> > eps1=(mass*ep.ref()->P())/(ep.mom(i0)*ep.ref()->P());
//	momentum<complex<T> > eps3=(mass*ep.ref()->P())/(ep.mom(i2)*ep.ref()->P());
//	momentum<complex<T> > eps1D=ddim;
//	momentum<complex<T> > eps2D=zero;
//	momentum<complex<T> > eps3D=-ddim;
//	
//	
//	//SSg
//	complex<T> res1(complex<T>(0,-1)*(((eps1*eps2)-(eps1D*eps2D))*(((ep.mom(i0)-ep.mom(i1))*eps3)-((mom1D-mom2D)*eps3D))
//									  +((eps2*eps3)-(eps2D*eps3D))*(((ep.mom(i1)-ep.mom(i2))*eps1)-((mom2D-mom3D)*eps1D))
//									  +((eps3*eps1)-(eps3D*eps1D))*(((ep.mom(i2)-ep.mom(i0))*eps2)-((mom3D-mom1D)*eps2D))));
//	
////	_MESSAGE4("compare orig b:",complex<T>(0,-1)*(ep.p(i1)->L()*ep.p(i0)->Sm()*ep.ref()->Lt())/(ep.p(i1)->Lt()*ep.ref()->Lt()),"new:",res1);

	return complex<T>(0,-1)*(ep.p(i1)->L()*ep.p(i0)->Sm()*ep.ref()->Lt())/(ep.p(i1)->Lt()*ep.ref()->Lt());
}

template <class T> complex<T> (*A2sc1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
		case 0x414://s(0)+s(0)
			return &A2sc1g1_eval<0,1,2>;
			break;
		case 0x404://s(0)-s(0)
			return &A2sc1g2_eval<0,1,2>;
			break;
		case 0x441://s(0)s(0)+
			return &A2sc1g1_eval<1,2,0>;
			break;
		case 0x440://s(0)s(0)-
			return &A2sc1g2_eval<1,2,0>;
			break;
		case 0x144://+s(0)s(0)
			return &A2sc1g1_eval<2,0,1>;
			break;
		case 0x044://-s(0)s(0)
			return &A2sc1g2_eval<2,0,1>;
			break;
		default:// We return zero for all other helicity combinations
			//cout << "missing:" << hex << hc << dec << endl;
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

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2sc2g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	//return complex<T>(0,1)*pow(ep.spab(i2,i0,i1),2)/(ep.s(i1,i2)*(-T(2.)*ep.sp(i0,i1)));
	return complex<T>(0,1)*ep.spa(i0,i2)*ep.spa(i0,i2)*ep.spa(i2,i3)/(ep.spa(i0,i1)*ep.spa(i1,i2)*ep.spa(i3,i0));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2sc2g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	//return complex<T>(0,1)*pow(ep.spba(i2,i0,i1),2)/(ep.s(i1,i2)*(-T(2.)*ep.sp(i0,i1)));
	return complex<T>(0,1)*ep.spa(i0,i1)*ep.spa(i1,i3)*ep.spa(i1,i3)/(ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i0));
}


template <class T> complex<T> (*A2sc2g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {			
		case 0x4104://s(0)+-s(0)
			return &A2sc2g3_eval<0,1,2,3>;
			break;
		case 0x4410://s(0)s(0)+-
			return &A2sc2g3_eval<1,2,3,0>;
			break;
		case 0x1044://+-s(0)s(0)
			return &A2sc2g3_eval<3,0,1,2>;
			break;
		case 0x0441://-s(0)s(0)+
			return &A2sc2g3_eval<2,3,0,1>;
			break;
			
		case 0x4014://s(0)-+s(0)
			return &A2sc2g4_eval<0,1,2,3>;
			break;
		case 0x1440://+s(0)s(0)-
			return &A2sc2g4_eval<2,3,0,1>;
			break;
		case 0x4401://s(0)s(0)-+
			return &A2sc2g4_eval<1,2,3,0>;
			break;
		case 0x0144://-+s(0)s(0)
			return &A2sc2g4_eval<3,0,1,2>;
			break;
            
		default:// We return zero for all other helicity combinations
			//cout << "missing:" << hex << hc << dec << endl;
			return &ZeroF;
	}
}


	
/*
 *
 *
 * The 1 complex higgs, two massless scalars and a gluon amplitude
 *
 *
 */
    
//phsms
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2scm1g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    
    return complex<T>(0,1)*ep.spa(i2,i1)*ep.spa(i2,i3)/ep.spa(i1,i3);
}
    
    
//phdsps
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2scm1g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    
    return complex<T>(0,1)*ep.spb(i2,i1)*ep.spb(i2,i3)/ep.spb(i1,i3);
}
    

template <class T> complex<T> (*A1ph2scm1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	//cout << hex << hc << dec << endl;
	switch (hc) {
			
			//
			//ph
			//
			
		case 0x9440://phssg
			return &A1ph2scm1g1_eval<0,2,3,1>;
			break;
		case 0x4940://sphsg
			return &A1ph2scm1g1_eval<1,2,3,0>;
			break;
		case 0x4490://ssphg
			return &A1ph2scm1g1_eval<2,1,3,0>;
			break;
		case 0x4409://ssgph
			return &A1ph2scm1g1_eval<3,1,2,0>;
			break;
						
			
		case 0x9404://phssg
			return &A1ph2scm1g1_eval<0,1,2,3>;
			break;
		case 0x4904://sphsg
			return &A1ph2scm1g1_eval<1,0,2,3>;
			break;
		case 0x4094://ssphg
			return &A1ph2scm1g1_eval<2,0,1,3>;
			break;
		case 0x4049://ssgph
			return &A1ph2scm1g1_eval<3,0,1,2>;
			break;
			
			
		case 0x9044://phssg
			return &A1ph2scm1g1_eval<0,3,1,2>;
			break;
		case 0x0944://sphsg
			return &A1ph2scm1g1_eval<1,3,0,2>;
			break;
		case 0x0494://ssphg
			return &A1ph2scm1g1_eval<2,3,0,1>;
			break;
		case 0x0449://ssgph
			return &A1ph2scm1g1_eval<3,2,0,1>;
			break;
				
            //
			// phd
			//
						
		case 0xA441://phdssg
			return &A1ph2scm1g3_eval<0,1,3,2>;
			break;
		case 0x4A41://sphdsg
			return &A1ph2scm1g3_eval<1,0,3,2>;
			break;
		case 0x44A1://ssphdg
			return &A1ph2scm1g3_eval<2,0,3,1>;
			break;
		case 0x441A://ssgphd
			return &A1ph2scm1g3_eval<3,0,2,1>;
			break;
			
						
		case 0xA414://phdssg
			return &A1ph2scm1g3_eval<0,1,2,3>;
			break;
		case 0x4A14://sphdsg
			return &A1ph2scm1g3_eval<1,0,2,3>;
			break;
		case 0x41A4://ssphdg
			return &A1ph2scm1g3_eval<2,0,1,3>;
			break;
		case 0x414A://ssgphd
			return &A1ph2scm1g3_eval<3,0,1,2>;
			break;
			
						
		case 0xA144://phdssg
			return &A1ph2scm1g3_eval<0,3,1,2>;
			break;
		case 0x1A44://sphdsg
			return &A1ph2scm1g3_eval<1,3,0,2>;
			break;
		case 0x14A4://ssphdg
			return &A1ph2scm1g3_eval<2,3,0,1>;
			break;
		case 0x144A://ssgpdh
			return &A1ph2scm1g3_eval<3,2,0,1>;
			break;
			
			
			
			
			//
			//H
			//
			
		case 0xB440://phssg
			return &A1ph2scm1g1_eval<0,2,3,1>;
			break;
		case 0x4B40://sphsg
			return &A1ph2scm1g1_eval<1,2,3,0>;
			break;
		case 0x44B0://ssphg
			return &A1ph2scm1g1_eval<2,1,3,0>;
			break;
		case 0x440B://ssgph
			return &A1ph2scm1g1_eval<3,1,2,0>;
			break;
			
			
		case 0xB404://phssg
			return &A1ph2scm1g1_eval<0,1,2,3>;
			break;
		case 0x4B04://sphsg
			return &A1ph2scm1g1_eval<1,0,2,3>;
			break;
		case 0x40B4://ssphg
			return &A1ph2scm1g1_eval<2,0,1,3>;
			break;
		case 0x404B://ssgph
			return &A1ph2scm1g1_eval<3,0,1,2>;
			break;
			
			
		case 0xB044://phssg
			return &A1ph2scm1g1_eval<0,3,1,2>;
			break;
		case 0x0B44://sphsg
			return &A1ph2scm1g1_eval<1,3,0,2>;
			break;
		case 0x04B4://ssphg
			return &A1ph2scm1g1_eval<2,3,0,1>;
			break;
		case 0x044B://ssgph
			return &A1ph2scm1g1_eval<3,2,0,1>;
			break;
			
 			
		case 0xB441://phdssg
			return &A1ph2scm1g3_eval<0,1,3,2>;
			break;
		case 0x4B41://sphdsg
			return &A1ph2scm1g3_eval<1,0,3,2>;
			break;
		case 0x44B1://ssphdg
			return &A1ph2scm1g3_eval<2,0,3,1>;
			break;
		case 0x441B://ssgphd
			return &A1ph2scm1g3_eval<3,0,2,1>;
			break;
			
			
		case 0xB414://phdssg
			return &A1ph2scm1g3_eval<0,1,2,3>;
			break;
		case 0x4B14://sphdsg
			return &A1ph2scm1g3_eval<1,0,2,3>;
			break;
		case 0x41B4://ssphdg
			return &A1ph2scm1g3_eval<2,0,1,3>;
			break;
		case 0x414B://ssgphd
			return &A1ph2scm1g3_eval<3,0,1,2>;
			break;
			
			
		case 0xB144://phdssg
			return &A1ph2scm1g3_eval<0,3,1,2>;
			break;
		case 0x1B44://sphdsg
			return &A1ph2scm1g3_eval<1,3,0,2>;
			break;
		case 0x14B4://ssphdg
			return &A1ph2scm1g3_eval<2,3,0,1>;
			break;
		case 0x144B://ssgpdh
			return &A1ph2scm1g3_eval<3,2,0,1>;
			break;
            
            			
		default:// We return zero for all other helicity combinations
			//_MESSAGE2("Missed",hc);
			return &ZeroF;
	}
}
    
    
    
    
/*
 *
 *
 * The 1 complex higgs and two massive scalars amplitudes
 *
 *
 */

//ph(d)ss
template <int i0, int i1, int i2, class T> complex<T>  A1ph2SM1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
    const complex<T> mass=eval_param<T>::mass2(masses.p(i2));
    
    return complex<T>(0,1)*mass;
}

template <class T> complex<T> (*A1ph2SM_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
    switch (hc) {
        case 0x944://phss
        case 0xA44://phdss
            return &A1ph2SM1_eval<0,1,2>;
            break;
        case 0x494://sphs
        case 0x4A4://sphds
            return &A1ph2SM1_eval<1,2,0>;
            break;
        case 0x449://ssph
        case 0x44A://ssphd
            return &A1ph2SM1_eval<2,0,1>;
            break;            
            
        default:// We return zero for all other helicity combinations
            return &ZeroF;
    }
}


///*
// *
// *
// * The 1 complex higgs, two massive scalars and a gluon amplitude
// *
// *
// */
//
//#define S(i,j) (ep.p(i)->P()*ep.p(j)->P())
//
////phsms
//template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2SM1g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));	
//    const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();
//    
//    const momentum<complex<T> > m12=ep.p(i1)->P()+ep.p(i2)->P();
//    const momentum<complex<T> > m23=ep.p(i2)->P()+ep.p(i3)->P();
//    const smatrix<T> Sm12=Cmom<T>(m12).Sm();
//    const smatrix<T> Sm23=Cmom<T>(m23).Sm();
//    
//    return complex<T>(0,-1)*(
//                             T(0.5)*_HIGGS_SCALAR_SCALAR_PARAM*(
//                                   (BIT0*mass)*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(T(2)*S(i3,i2))
//                                -(BIT0*mass)*(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())/(T(2)*S(i1,i2))        
//                                                                )/(ep.p(i2)->Lt()*ep.ref()->Lt())
//                             +_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->Lt()*ep.ref()->Lt())))*(
//                                       (ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())
//                                     -(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.p(i2)->Lt())
//                             ));
//}
//
////phsps
//template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2SM1g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
//    
//    const momentum<complex<T> > m12=ep.p(i1)->P()+ep.p(i2)->P();
//    const momentum<complex<T> > m23=ep.p(i2)->P()+ep.p(i3)->P();
//    
//    return complex<T>(0,1)*(
//                            T(0.5)*_HIGGS_SCALAR_SCALAR_PARAM_D*(
//                                (BIT0*mass)*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//                            -(BIT0*mass)*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))
//                                )/(ep.p(i2)->L()*ep.ref()->L())
//                            );
//}
//
////phdsms
//template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2SM1g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
//    
//    const momentum<complex<T> > m12=ep.p(i1)->P()+ep.p(i2)->P();
//    const momentum<complex<T> > m23=ep.p(i2)->P()+ep.p(i3)->P();
//    
//    return complex<T>(0,-1)*(
//                             T(0.5)*_HIGGS_SCALAR_SCALAR_PARAM_D*(
//                                                                  (BIT0*mass
//                                                                   +(BIT1*(m23*ep.p(i1)->P())
//                                                                     -BIT2*mass*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//                                                                     /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))
//                                                                     )
//                                                                   )
//                                                                  *(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())/(T(2)*S(i3,i2))
//                                                                  -(BIT0*mass
//                                                                    +(BIT1*(m12*ep.p(i3)->P())
//                                                                      -BIT2*mass*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//                                                                      /((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))
//                                                                      )
//                                                                    )
//                                                                  *(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())/(T(2)*S(i1,i2))
//                                                                  )/(ep.p(i2)->Lt()*ep.ref()->Lt())
//                             );
//}
//
////phdsps
//template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2SM1g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//    const complex<T> mass=eval_param<T>::mass2(masses.p(i1));
//    
//    const momentum<complex<T> > m12=ep.p(i1)->P()+ep.p(i2)->P();
//    const momentum<complex<T> > m23=ep.p(i2)->P()+ep.p(i3)->P();
//    
//    const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();
//    
//    return complex<T>(0,1)*(
//                            T(0.5)*_HIGGS_SCALAR_SCALAR_PARAM*(
//                                                               (BIT0*mass+(BIT1*(m23*ep.p(i1)->P())
//                                                                           -BIT2*mass*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//                                                                           /((ep.p(i1)->P()*ep.ref()->P())*(m23*ep.ref()->P()))))
//                                                               *(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())/(T(2)*S(i3,i2))
//                                                               -(BIT0*mass+(BIT1*(m12*ep.p(i3)->P())
//                                                                            -BIT2*mass*pow(ep.p(i0)->P()*ep.ref()->P(),2)
//                                                                            /((ep.p(i3)->P()*ep.ref()->P())*(m12*ep.ref()->P()))))
//                                                               *(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())/(T(2)*S(i1,i2))      
//                                                               )/(ep.p(i2)->L()*ep.ref()->L())
//                            );
//}
//
//template <class T> complex<T> (*A1ph2SM1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
//    //		cout << hex << hc << dec << endl;
//    switch (hc) {
//            
//            //
//            //ph
//            //
//            
//        case 0x9404://phssg
//            return &A1ph2SM1g1_eval<0,1,2,3>;
//            break;
//        case 0x4904://sphsg
//            return &A1ph2SM1g1_eval<1,0,2,3>;
//            break;
//        case 0x4094://ssphg
//            return &A1ph2SM1g1_eval<2,0,1,3>;
//            break;
//        case 0x4049://ssgph
//            return &A1ph2SM1g1_eval<3,0,1,2>;
//            break;
//            
//        case 0x9414://phssg
//            return &A1ph2SM1g2_eval<0,1,2,3>;
//            break;
//        case 0x4914://sphsg
//            return &A1ph2SM1g2_eval<1,0,2,3>;
//            break;
//        case 0x4194://ssphg
//            return &A1ph2SM1g2_eval<2,0,1,3>;
//            break;
//        case 0x4149://ssgph
//            return &A1ph2SM1g2_eval<3,0,1,2>;
//            break;
//            
//        case 0x9440://phssg
//            return &A1ph2SM1g1_eval<0,2,3,1>;
//            break;
//        case 0x4940://sphsg
//            return &A1ph2SM1g1_eval<1,2,3,0>;
//            break;
//        case 0x4490://ssphg
//            return &A1ph2SM1g1_eval<2,1,3,0>;
//            break;
//        case 0x4409://ssgph
//            return &A1ph2SM1g1_eval<3,1,2,0>;
//            break;
//            
//        case 0x9441://phssg
//            return &A1ph2SM1g2_eval<0,2,3,1>;
//            break;
//        case 0x4941://sphsg
//            return &A1ph2SM1g2_eval<1,2,3,0>;
//            break;
//        case 0x4491://ssphg
//            return &A1ph2SM1g2_eval<2,1,3,0>;
//            break;
//        case 0x4419://ssgph
//            return &A1ph2SM1g2_eval<3,1,2,0>;
//            break;
//            
//            
//            
//        case 0x9044://phssg
//            return &A1ph2SM1g1_eval<0,3,1,2>;
//            break;
//        case 0x0944://sphsg
//            return &A1ph2SM1g1_eval<1,3,0,2>;
//            break;
//        case 0x0494://ssphg
//            return &A1ph2SM1g1_eval<2,3,0,1>;
//            break;
//        case 0x0449://ssgph
//            return &A1ph2SM1g1_eval<3,2,0,1>;
//            break;
//            
//        case 0x9144://phssg
//            return &A1ph2SM1g2_eval<0,3,1,2>;
//            break;
//        case 0x1944://sphsg
//            return &A1ph2SM1g2_eval<1,3,0,2>;
//            break;
//        case 0x1494://ssphg
//            return &A1ph2SM1g2_eval<2,3,0,1>;
//            break;
//        case 0x1449://ssgph
//            return &A1ph2SM1g2_eval<3,2,0,1>;
//            break;
//
//            
//            
//            //
//            // phd
//            //
//            
//        case 0xA440://phdssg
//            return &A1ph2SM1g4_eval<0,2,3,1>;
//            break;
//        case 0x4A40://spdhsg
//            return &A1ph2SM1g4_eval<1,2,3,0>;
//            break;
//        case 0x44A0://ssphdg
//            return &A1ph2SM1g4_eval<2,1,3,0>;
//            break;
//        case 0x440A://ssgphd
//            return &A1ph2SM1g4_eval<3,1,2,0>;
//            break;
//            
//        case 0xA441://phdssg
//            return &A1ph2SM1g3_eval<0,2,3,1>;
//            break;
//        case 0x4A41://sphdsg
//            return &A1ph2SM1g3_eval<1,2,3,0>;
//            break;
//        case 0x44A1://ssphdg
//            return &A1ph2SM1g3_eval<2,1,3,0>;
//            break;
//        case 0x441A://ssgphd
//            return &A1ph2SM1g3_eval<3,1,2,0>;
//            break;
//            
//            
//        case 0xA404://phdssg
//            return &A1ph2SM1g4_eval<0,1,2,3>;
//            break;
//        case 0x4A04://sphdsg
//            return &A1ph2SM1g4_eval<1,0,2,3>;
//            break;
//        case 0x40A4://ssphdg
//            return &A1ph2SM1g4_eval<2,0,1,3>;
//            break;
//        case 0x404A://ssgphd
//            return &A1ph2SM1g4_eval<3,0,1,2>;
//            break;
//            
//        case 0xA414://phdssg
//            return &A1ph2SM1g3_eval<0,1,2,3>;
//            break;
//        case 0x4A14://sphdsg
//            return &A1ph2SM1g3_eval<1,0,2,3>;
//            break;
//        case 0x41A4://ssphdg
//            return &A1ph2SM1g3_eval<2,0,1,3>;
//            break;
//        case 0x414A://ssgphd
//            return &A1ph2SM1g3_eval<3,0,1,2>;
//            break;
//            
//            
//        case 0xA044://phdssg
//            return &A1ph2SM1g4_eval<0,3,1,2>;
//            break;
//        case 0x0A44://sphdsg
//            return &A1ph2SM1g4_eval<1,3,0,2>;
//            break;
//        case 0x04A4://ssphdg
//            return &A1ph2SM1g4_eval<2,3,0,1>;
//            break;
//        case 0x044A://ssgphd
//            return &A1ph2SM1g4_eval<3,2,0,1>;
//            break;
//            
//        case 0xA144://phdssg
//            return &A1ph2SM1g3_eval<0,3,1,2>;
//            break;
//        case 0x1A44://sphdsg
//            return &A1ph2SM1g3_eval<1,3,0,2>;
//            break;
//        case 0x14A4://ssphdg
//            return &A1ph2SM1g3_eval<2,3,0,1>;
//            break;
//        case 0x144A://ssgpdh
//            return &A1ph2SM1g3_eval<3,2,0,1>;
//            break;
//            
//            
//
//            
//            
//        default:// We return zero for all other helicity combinations
//            //			_MESSAGE2("Missed",hc);
//            return 0;
//    }
//}
	


/*
 *
 *
 * The 1 complex higgs, two massive scalars and a gluon amplitude
 *
 *
 */

#define S(i,j) (ep.p(i)->P()*ep.p(j)->P())

//phsms
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2SM1g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();
	
	return complex<T>(0,1)*((T(1)/(den*(ep.p(i2)->Lt()*ep.ref()->Lt())))*((ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())
										-(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.p(i2)->Lt())));
}

//phdsps
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2SM1g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();
	
	return complex<T>(0,-1)*((T(1)/(den*(ep.p(i2)->L()*ep.ref()->L())))*((ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())
					  -(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.p(i2)->L())));
}


template <class T> complex<T> (*A1ph2SM1g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	//	cout << hex << hc << dec << endl;
	switch (hc) {
			
			//
			//ph
			//
			
		case 0x9550://phssg
			return &A1ph2SM1g1_eval<0,2,3,1>;
			break;
		case 0x5950://sphsg
			return &A1ph2SM1g1_eval<1,2,3,0>;
			break;
		case 0x5590://ssphg
			return &A1ph2SM1g1_eval<2,1,3,0>;
			break;
		case 0x5509://ssgph
			return &A1ph2SM1g1_eval<3,1,2,0>;
			break;	
			
		case 0x9505://phssg
			return &A1ph2SM1g1_eval<0,1,2,3>;
			break;
		case 0x5905://sphsg
			return &A1ph2SM1g1_eval<1,0,2,3>;
			break;
		case 0x5095://ssphg
			return &A1ph2SM1g1_eval<2,0,1,3>;
			break;
		case 0x5059://ssgph
			return &A1ph2SM1g1_eval<3,0,1,2>;
			break;
						
		case 0x9055://phssg
			return &A1ph2SM1g1_eval<0,3,1,2>;
			break;
		case 0x0955://sphsg
			return &A1ph2SM1g1_eval<1,3,0,2>;
			break;
		case 0x0595://ssphg
			return &A1ph2SM1g1_eval<2,3,0,1>;
			break;
		case 0x0559://ssgph
			return &A1ph2SM1g1_eval<3,2,0,1>;
			break;
			

			//
			// phd
			//
						
		case 0xA551://phdssg
			return &A1ph2SM1g3_eval<0,2,3,1>;
			break;
		case 0x5A51://sphdsg
			return &A1ph2SM1g3_eval<1,2,3,0>;
			break;
		case 0x55A1://ssphdg
			return &A1ph2SM1g3_eval<2,1,3,0>;
			break;
		case 0x551A://ssgphd
			return &A1ph2SM1g3_eval<3,1,2,0>;
			break;
			
			
		case 0xA515://phdssg
			return &A1ph2SM1g3_eval<0,1,2,3>;
			break;
		case 0x5A15://sphdsg
			return &A1ph2SM1g3_eval<1,0,2,3>;
			break;
		case 0x51A5://ssphdg
			return &A1ph2SM1g3_eval<2,0,1,3>;
			break;
		case 0x515A://ssgphd
			return &A1ph2SM1g3_eval<3,0,1,2>;
			break;
			
						
		case 0xA155://phdssg
			return &A1ph2SM1g3_eval<0,3,1,2>;
			break;
		case 0x1A55://sphdsg
			return &A1ph2SM1g3_eval<1,3,0,2>;
			break;
		case 0x15A5://ssphdg
			return &A1ph2SM1g3_eval<2,3,0,1>;
			break;
		case 0x155A://ssgpdh
			return &A1ph2SM1g3_eval<3,2,0,1>;
			break;
			
			
			
			//
			//H
			//
			
		case 0xB550://phssg
			return &A1ph2SM1g1_eval<0,2,3,1>;
			break;
		case 0x5B50://sphsg
			return &A1ph2SM1g1_eval<1,2,3,0>;
			break;
		case 0x55B0://ssphg
			return &A1ph2SM1g1_eval<2,1,3,0>;
			break;
		case 0x550B://ssgph
			return &A1ph2SM1g1_eval<3,1,2,0>;
			break;	
			
		case 0xB505://phssg
			return &A1ph2SM1g1_eval<0,1,2,3>;
			break;
		case 0x5B05://sphsg
			return &A1ph2SM1g1_eval<1,0,2,3>;
			break;
		case 0x50B5://ssphg
			return &A1ph2SM1g1_eval<2,0,1,3>;
			break;
		case 0x505B://ssgph
			return &A1ph2SM1g1_eval<3,0,1,2>;
			break;
			
		case 0xB055://phssg
			return &A1ph2SM1g1_eval<0,3,1,2>;
			break;
		case 0x0B55://sphsg
			return &A1ph2SM1g1_eval<1,3,0,2>;
			break;
		case 0x05B5://ssphg
			return &A1ph2SM1g1_eval<2,3,0,1>;
			break;
		case 0x055B://ssgph
			return &A1ph2SM1g1_eval<3,2,0,1>;
			break;
			
		case 0xB551://phdssg
			return &A1ph2SM1g3_eval<0,2,3,1>;
			break;
		case 0x5B51://sphdsg
			return &A1ph2SM1g3_eval<1,2,3,0>;
			break;
		case 0x55B1://ssphdg
			return &A1ph2SM1g3_eval<2,1,3,0>;
			break;
		case 0x551B://ssgphd
			return &A1ph2SM1g3_eval<3,1,2,0>;
			break;
			
			
		case 0xB515://phdssg
			return &A1ph2SM1g3_eval<0,1,2,3>;
			break;
		case 0x5B15://sphdsg
			return &A1ph2SM1g3_eval<1,0,2,3>;
			break;
		case 0x51B5://ssphdg
			return &A1ph2SM1g3_eval<2,0,1,3>;
			break;
		case 0x515B://ssgphd
			return &A1ph2SM1g3_eval<3,0,1,2>;
			break;
			
			
		case 0xB155://phdssg
			return &A1ph2SM1g3_eval<0,3,1,2>;
			break;
		case 0x1B55://sphdsg
			return &A1ph2SM1g3_eval<1,3,0,2>;
			break;
		case 0x15B5://ssphdg
			return &A1ph2SM1g3_eval<2,3,0,1>;
			break;
		case 0x155B://ssgpdh
			return &A1ph2SM1g3_eval<3,2,0,1>;
			break;

			
			
		default:// We return zero for all other helicity combinations
//			cout << hex<< "Missed: " <<hc<<dec<<endl;
			return &ZeroF;
	}
}
	
	
/*
 *
 *
 * The 1 complex higgs, 1 massive quark a massive scalar and a massless quark amplitude
 *
 *
 */
                         
//(s,q-,Q+,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
 const lambda<T> qL=ep.ref()->L();
 
 const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
 
 return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->L()*qL)/(den*(QL*qL)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
 const lambda<T> qL=ep.ref()->L();
 
 const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
 
 return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->L()*qL)/(den*(QL*qL)*sqrt(T(2)));
}

//(s,q-,Q-,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
 
 const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
 
 return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(ep.p(i1)->L()*QL)/(den*sqrt(T(2)));
}

//(Q-,q-,s,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
 
 const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
 
 return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(QL*ep.p(i1)->L())/(den*sqrt(T(2)));
}

//(Q+,q-,s,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
 const lambda<T> qL=ep.ref()->L();
 
 const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
 
 return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qL*ep.p(i1)->L())/(den*(qL*QL)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
 const lambda<T> qL=ep.ref()->L();
 
 const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
 
 return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qL*ep.p(i1)->L())/(den*(qL*QL)*sqrt(T(2)));
}

//(Q-,q+,s,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q5m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
 const lambdat<T> qLt=ep.ref()->Lt();
 
 const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
 
 return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qLt*ep.p(i1)->Lt())/(den*(qLt*QLt)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q5p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
 const lambdat<T> qLt=ep.ref()->Lt();
 
 const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
 
 return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qLt*ep.p(i1)->Lt())/(den*(qLt*QLt)*sqrt(T(2)));
}

//(Q+,q+,s,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
 
 const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
 
 return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(QLt*ep.p(i1)->Lt())/(den*sqrt(T(2)));
}

//(s,q+,Q+,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q7_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
 
 const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
 
 return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(ep.p(i1)->Lt()*QLt)/(den*sqrt(T(2)));
}

//(s,q+,Q-,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
 const lambdat<T> qLt=ep.ref()->Lt();
 
 const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
 
 return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->Lt()*qLt)/(den*(QLt*qLt)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1SM1q8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
 const lambdat<T> qLt=ep.ref()->Lt();
 
 const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
 
 return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->Lt()*qLt)/(den*(QLt*qLt)*sqrt(T(2)));
}
 
 
template <class T, int SPOS, int QPOS, int QMPOS> complex<T> (*A1ph1QM1SM1q_phi_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
     //	cout << hex << hc << dec << endl;
     switch (hc) {
         case 0x428: //258: //(s,q3-,Qb+)
             return &A1ph1QM1SM1q1m_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x426:
             return &A1ph1QM1SM1q1p_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x284: //330: //(qb-,Q+,s)
             return &A1ph1QM1SM1q1m_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x264:
             return &A1ph1QM1SM1q1p_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x842: //96: //(Q+,s,qb-)
             return &A1ph1QM1SM1q1m_eval<QPOS,QMPOS,SPOS>;
             break;
         case 0x642:
             return &A1ph1QM1SM1q1p_eval<QPOS,QMPOS,SPOS>;
             break;
             
             
         case 0x427: //209: //(s,q3-,Qb-)
         case 0x425:
             return &A1ph1QM1SM1q2_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x274: //323: //(qb-,Q-,s)
         case 0x254:
             return &A1ph1QM1SM1q2_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x742: //95: //(qb-,Q-,s)
         case 0x542:
             return &A1ph1QM1SM1q2_eval<QPOS,QMPOS,SPOS>;
             break;
             
             
             
         case 0x724: //305: //(Q-,q3-,s)
         case 0x524:
             return &A1ph1QM1SM1q3_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x247: //239: //(qb-,s,Q-)
         case 0x245:
             return &A1ph1QM1SM1q3_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x472: //83: //(s,Q3-,qb3-)
         case 0x452:
             return &A1ph1QM1SM1q3_eval<QPOS,QMPOS,SPOS>;
             break;
             
             
             
         case 0x824: //306: //(Q+,q3-,s)
             return &A1ph1QM1SM1q4m_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x624:
             return &A1ph1QM1SM1q4p_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x248: //288: //(qb-,s,Q+)
             return &A1ph1QM1SM1q4m_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x246:
             return &A1ph1QM1SM1q4p_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x482: //90: //(s,Q3+,qb3-)
             return &A1ph1QM1SM1q4m_eval<QPOS,QMPOS,SPOS>;
             break;
         case 0x462:
             return &A1ph1QM1SM1q4p_eval<QPOS,QMPOS,SPOS>;
             break;
             
             
         case 0x734: //312: //(Q-,qb+,s)
             return &A1ph1QM1SM1q5m_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x534:
             return &A1ph1QM1SM1q5p_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x347: //240: //(q+,s,Qb-)
             return &A1ph1QM1SM1q5m_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x345:
             return &A1ph1QM1SM1q5p_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x473: //132: //(s,Qb-,q+)
             return &A1ph1QM1SM1q5m_eval<QPOS,QMPOS,SPOS>;
             break;
         case 0x453:
             return &A1ph1QM1SM1q5p_eval<QPOS,QMPOS,SPOS>;
             break;
             
             
         case 0x834: //313: //(Q+,qb+,s)
         case 0x634:
             return &A1ph1QM1SM1q6_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x348: //289: //(q+,s,Qb+)
         case 0x346:
             return &A1ph1QM1SM1q6_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x483: //139: //(s,Qb+,q+)
         case 0x463:
             return &A1ph1QM1SM1q6_eval<QPOS,QMPOS,SPOS>;
             break;
             
             
         case 0x438: //265: //(s,qb+,Qb+)
         case 0x436:
             return &A1ph1QM1SM1q7_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x384: //331: //(q+,Qb+,s)
         case 0x364:
             return &A1ph1QM1SM1q7_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x843: //145: //(Qb+,s,q+)
         case 0x643:
             return &A1ph1QM1SM1q7_eval<QPOS,QMPOS,SPOS>;
             break;
             
             
         case 0x437: //216: //(s,qb+,Qb-)
             return &A1ph1QM1SM1q8m_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x435:
             return &A1ph1QM1SM1q8p_eval<SPOS,QPOS,QMPOS>;
             break;
         case 0x374: //324: //(q+,Qb-,s)
             return &A1ph1QM1SM1q8m_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x354:
             return &A1ph1QM1SM1q8p_eval<QMPOS,SPOS,QPOS>;
             break;
         case 0x743: //144: //(Qb-,s,q+)
             return &A1ph1QM1SM1q8m_eval<QPOS,QMPOS,SPOS>;
             break;
         case 0x543:
             return &A1ph1QM1SM1q8p_eval<QPOS,QMPOS,SPOS>;
             break;
             
             
         default:// We return zero for all other helicity combinations
             //cout << "Missed vertex:" << hex << hc << dec << endl;
             return &ZeroF;
     }
}
 
template <class T> complex<T> (*A1ph1QM1SM1q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
     // We first extract the position of the higgs and then analyse the colour ordered particles in the next level down
     long rem=hc;
     long newhc=0;
     int epos;
     long fac=0x1;
     bool dagger=false;
	 bool full_higgs=true;
     for(int i=0;i<4;i++){
         int efac=rem%0x10;
         rem=rem/0x10;
         if(efac==0x9||efac==0xA||efac==0xB){
             epos=i;
             if(efac==0xA){// If this is a phi dagger then set record this
                 dagger=true;
             }
			 if(efac==0xB){// If this is a phi dagger then set record this
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
					 return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,1,2>(newhc);
					 break;
				 case 1:
					 return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,1,3>(newhc);
					 break;
				 case 2:
					 return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,2,3>(newhc);
					 break;
				 case 3:
					 return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,1,2,3>(newhc);
					 break;
					 
				 default:// Should never get here as there is always a higgs. 
					 return &ZeroF;
					 break;
			 }
		 }
		 else{
			 switch (epos) {
				 case 0:
					 return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,1,2>(newhc);
					 break;
				 case 1:
					 return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,1,3>(newhc);
					 break;
				 case 2:
					 return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,2,3>(newhc);
					 break;
				 case 3:
					 return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,1,2,3>(newhc);
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
				return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,1,2>(newhc);
				break;
			case 1:
				return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,1,3>(newhc);
				break;
			case 2:
				return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,0,2,3>(newhc);
				break;
			case 3:
				return A1ph1QM1SM1q_phi_Tree_Ptr_eval<T,1,2,3>(newhc);
				break;
				
			default:// Should never get here as there is always a higgs. 
				return &ZeroF;
				break;
		}
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
                         
//phsqmqps
template <int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A1ph2SM2q1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 
 const momentum<complex<T> > sum34(ep.p(i3)->P()+ep.p(i4)->P());
 const momentum<complex<T> > sum25(ep.p(i2)->P()+ep.p(i5)->P());
 const momentum<complex<T> > sum15(ep.p(i1)->P()+ep.p(i5)->P());
 const momentum<complex<T> > sum12(ep.p(i1)->P()+ep.p(i2)->P());
 const momentum<complex<T> > sum45(ep.p(i4)->P()+ep.p(i5)->P());
 const momentum<complex<T> > sum23(ep.p(i2)->P()+ep.p(i3)->P());
 const complex<T> Isqrt12=T(1)/sqrt(sum12.square());
 const complex<T> Isqrt15=T(1)/sqrt(sum15.square());
 const complex<T> sI34=T(2)/sum34.square();
 const complex<T> sI25=T(1)/sum25.square();
 const complex<T> sI15=T(1)/(sum15.square()-mass2);
 const complex<T> sI12=T(1)/(sum12.square()-mass2);
 const complex<T> sI23=T(1)/(sum23.square()-mass2);
 const complex<T> sI45=T(1)/(sum45.square()-mass2);
 
// _MESSAGE("1:");
// _PRINT(ep.spab(i3,i2,i4)*(sI15*(_HIGGS_PARAM_2)));
// _PRINT(ep.spab(i3,i2,i4)*(sI23*(_HIGGS_PARAM_3)));
// _PRINT(ep.spab(i3,i2,i4)*(_OTHER_PARAM_2*(sum34*ep.mom(i5))*sI34*sI25));
// 
// _PRINT(-ep.spab(i3,i5,i4)*(sI12*(_HIGGS_PARAM_2)));
// _PRINT(-ep.spab(i3,i5,i4)*(sI45*(_HIGGS_PARAM_3)));
// _PRINT(-ep.spab(i3,i5,i4)*(_OTHER_PARAM_2*(sum34*ep.mom(i2))*sI34*sI25));    
 
	return complex<T>(0,0);
	
 return complex<T>(0,-1)*(ep.spab(i3,i2,i4)*(
                                             _HIGGS_PARAM_2*sI15/**sI23*sum23.square()*/
                                             +_HIGGS_PARAM_3*sI23
                                             +_OTHER_PARAM_2*(sum34*ep.mom(i5))*sI34*sI25)
                          -ep.spab(i3,i5,i4)*(
                                              _HIGGS_PARAM_2*sI12/**sI45*sum45.square()*/
                                              +_HIGGS_PARAM_3*sI45
                                              +_OTHER_PARAM_2*(sum34*ep.mom(i2))*sI34*sI25));
}


//phsqpqms
template <int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A1ph2SM2q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
 const complex<T> mass=eval_param<T>::mass(masses.p(i2));
 
 const momentum<complex<T> > sum34(ep.p(i3)->P()+ep.p(i4)->P());
 const momentum<complex<T> > sum25(ep.p(i2)->P()+ep.p(i5)->P());
 const momentum<complex<T> > sum15(ep.p(i1)->P()+ep.p(i5)->P());
 const momentum<complex<T> > sum12(ep.p(i1)->P()+ep.p(i2)->P());
 const momentum<complex<T> > sum45(ep.p(i4)->P()+ep.p(i5)->P());
 const momentum<complex<T> > sum23(ep.p(i2)->P()+ep.p(i3)->P());
 const complex<T> Isqrt12=T(1)/sqrt(sum12.square());
 const complex<T> Isqrt15=T(1)/sqrt(sum15.square());
 const complex<T> sI34=T(2)/sum34.square();
 const complex<T> sI25=T(1)/sum25.square();
 const complex<T> sI15=T(1)/(sum15.square()-mass2);
 const complex<T> sI12=T(1)/(sum12.square()-mass2);
 const complex<T> sI23=T(1)/(sum23.square()-mass2);
 const complex<T> sI45=T(1)/(sum45.square()-mass2);
 
// _MESSAGE("2:");
// _PRINT(ep.spba(i3,i2,i4)*(sI15*(_HIGGS_PARAM_2)));
// _PRINT(ep.spba(i3,i2,i4)*(sI23*(_HIGGS_PARAM_3)));
// _PRINT(ep.spba(i3,i2,i4)*(_OTHER_PARAM_2*(sum34*ep.mom(i5))*sI34*sI25));
// 
// _PRINT(-ep.spba(i3,i5,i4)*(sI12*(_HIGGS_PARAM_2)));
// _PRINT(-ep.spba(i3,i5,i4)*(sI45*(_HIGGS_PARAM_3)));
// _PRINT(-ep.spba(i3,i5,i4)*(_OTHER_PARAM_2*(sum34*ep.mom(i2))*sI34*sI25));
 
	return complex<T>(0,0);

	
 return complex<T>(0,1)*(
                         ep.spba(i3,i2,i4)*(
                                            _HIGGS_PARAM_2*sI15/**sI23*sum23.square()*/
                                            +_HIGGS_PARAM_3*sI23
                                            +_OTHER_PARAM_2*(sum34*ep.mom(i5))*sI34*sI25)
                         -ep.spba(i3,i5,i4)*(
                                             _HIGGS_PARAM_2*sI12/**sI45*sum45.square()*/
                                             +_HIGGS_PARAM_3*sI45
                                             +_OTHER_PARAM_2*(sum34*ep.mom(i2))*sI34*sI25));
 
}

template <class T> complex<T> (*A1ph2SM2q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
    switch (hc) {
        case 0x94234: //2169: //(ph,s,qb-,q+,s)
            return &A1ph2SM2q1_eval<0,1,2,3,4>;
            break;
        case 0x49234: //2169: //(s,ph,qb-,q+,s)
            return &A1ph2SM2q1_eval<1,0,2,3,4>;
            break;
        case 0x42934: //2169: //(s,qb-,ph,q+,s)
            return &A1ph2SM2q1_eval<2,0,1,3,4>;
            break;
        case 0x42394: //2169: //(s,qb-,q+,ph,s)
            return &A1ph2SM2q1_eval<3,0,1,2,4>;
            break;
        case 0x42349: //2169: //(s,qb-,q+,s,ph)
            return &A1ph2SM2q1_eval<4,0,1,2,3>;
            break;
            
        case 0x94423: //783: //(ph,s,s,qb-,q+)
            return &A1ph2SM2q1_eval<0,2,3,4,1>;
            break;
        case 0x49423: //783: //(s,ph,s,qb-,q+)
            return &A1ph2SM2q1_eval<1,2,3,4,0>;
            break;
        case 0x44923: //783: //(s,s,ph,qb-,q+)
            return &A1ph2SM2q1_eval<2,1,3,4,0>;
            break;
        case 0x44293: //783: //(s,s,qb-,ph,q+)
            return &A1ph2SM2q1_eval<3,1,2,4,0>;
            break;
        case 0x44239: //783: //(ph,s,s,qb-,q+,ph)
            return &A1ph2SM2q1_eval<4,1,2,3,0>;
            break;
            
        case 0x93442: //681: //(ph,q+,s,s,qb-)
            return &A1ph2SM2q1_eval<0,3,4,1,2>;
            break;
        case 0x39442: //681: //(q+,ph,s,s,qb-)
            return &A1ph2SM2q1_eval<1,3,4,0,2>;
            break;
        case 0x34942: //681: //(q+,s,ph,s,qb-)
            return &A1ph2SM2q1_eval<2,3,4,0,1>;
            break;
        case 0x34492: //681: //(q+,s,s,ph,qb-)
            return &A1ph2SM2q1_eval<3,2,4,0,1>;
            break;
        case 0x34429: //681: //(q+,s,s,qb-,ph)
            return &A1ph2SM2q1_eval<4,2,3,0,1>;
            break;
            
        case 0x92344: //2367: //(ph,qb-,q+,s,s)
            return &A1ph2SM2q1_eval<0,4,1,2,3>;
            break;
        case 0x29344: //2367: //(qb-,ph,q+,s,s)
            return &A1ph2SM2q1_eval<1,4,0,2,3>;
            break;
        case 0x23944: //2367: //(qb-,q+,ph,s,s)
            return &A1ph2SM2q1_eval<2,4,0,1,3>;
            break;
        case 0x23494: //2367: //(qb-,q+,s,ph,s)
            return &A1ph2SM2q1_eval<3,4,0,1,2>;
            break;
        case 0x23449: //2367: //(qb-,q+,s,s,ph)
            return &A1ph2SM2q1_eval<4,3,0,1,2>;
            break;
            
            
        case 0x94324: //2127: //(ph,s,q+,qb-,s)
            return &A1ph2SM2q2_eval<0,1,2,3,4>;
            break;
        case 0x49324: //2127: //(s,ph,q+,qb-,s)
            return &A1ph2SM2q2_eval<1,0,2,3,4>;
            break;
        case 0x43924: //2127: //(s,q+,ph,qb-,s)
            return &A1ph2SM2q2_eval<2,0,1,3,4>;
            break;
        case 0x43294: //2127: //(s,q+,qb-,ph,s)
            return &A1ph2SM2q2_eval<3,0,1,2,4>;
            break;
        case 0x43249: //2127: //(s,q+,qb-,s,ph)
            return &A1ph2SM2q2_eval<4,0,1,2,3>;
            break;
            
        case 0x94432: //489: //(ph,s,s,q+,qb-)
            return &A1ph2SM2q2_eval<0,2,3,4,1>;
            break;
        case 0x49432: //489: //(s,ph,s,q+,qb-)
            return &A1ph2SM2q2_eval<1,2,3,4,0>;
            break;
        case 0x44932: //489: //(s,s,ph,q+,qb-)
            return &A1ph2SM2q2_eval<2,1,3,4,0>;
            break;
        case 0x44392: //489: //(s,s,q+,ph,qb-)
            return &A1ph2SM2q2_eval<3,1,2,4,0>;
            break;
        case 0x44329: //489: //(s,s,q+,qb-,ph)
            return &A1ph2SM2q2_eval<4,1,2,3,0>;
            break;
            
        case 0x92443: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<0,3,4,1,2>;
            break;
        case 0x29443: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<1,3,4,0,2>;
            break;
        case 0x24943: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<2,3,4,0,1>;
            break;
        case 0x24493: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<3,2,4,0,1>;
            break;
        case 0x24439: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<4,2,3,0,1>;
            break;
            
        case 0x93244: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<0,4,1,2,3>;
            break;
        case 0x39244: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<1,4,0,2,3>;
            break;
        case 0x32944: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<2,4,0,1,3>;
            break;
        case 0x32494: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<3,4,0,1,2>;
            break;
        case 0x32449: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<4,3,0,1,2>;
            break;
            
            
            //
            //phd
            //
            
        case 0xA4234: //216A: //(ph,s,qb-,q+,s)
            return &A1ph2SM2q1_eval<0,1,2,3,4>;
            break;
        case 0x4A234: //216A: //(s,ph,qb-,q+,s)
            return &A1ph2SM2q1_eval<1,0,2,3,4>;
            break;
        case 0x42A34: //216A: //(s,qb-,ph,q+,s)
            return &A1ph2SM2q1_eval<2,0,1,3,4>;
            break;
        case 0x423A4: //216A: //(s,qb-,q+,ph,s)
            return &A1ph2SM2q1_eval<3,0,1,2,4>;
            break;
        case 0x4234A: //216A: //(s,qb-,q+,s,ph)
            return &A1ph2SM2q1_eval<4,0,1,2,3>;
            break;
            
        case 0xA4423: //783: //(ph,s,s,qb-,q+)
            return &A1ph2SM2q1_eval<0,2,3,4,1>;
            break;
        case 0x4A423: //783: //(s,ph,s,qb-,q+)
            return &A1ph2SM2q1_eval<1,2,3,4,0>;
            break;
        case 0x44A23: //783: //(s,s,ph,qb-,q+)
            return &A1ph2SM2q1_eval<2,1,3,4,0>;
            break;
        case 0x442A3: //783: //(s,s,qb-,ph,q+)
            return &A1ph2SM2q1_eval<3,1,2,4,0>;
            break;
        case 0x4423A: //783: //(ph,s,s,qb-,q+,ph)
            return &A1ph2SM2q1_eval<4,1,2,3,0>;
            break;
            
        case 0xA3442: //681: //(ph,q+,s,s,qb-)
            return &A1ph2SM2q1_eval<0,3,4,1,2>;
            break;
        case 0x3A442: //681: //(q+,ph,s,s,qb-)
            return &A1ph2SM2q1_eval<1,3,4,0,2>;
            break;
        case 0x34A42: //681: //(q+,s,ph,s,qb-)
            return &A1ph2SM2q1_eval<2,3,4,0,1>;
            break;
        case 0x344A2: //681: //(q+,s,s,ph,qb-)
            return &A1ph2SM2q1_eval<3,2,4,0,1>;
            break;
        case 0x3442A: //681: //(q+,s,s,qb-,ph)
            return &A1ph2SM2q1_eval<4,2,3,0,1>;
            break;
            
        case 0xA2344: //2367: //(ph,qb-,q+,s,s)
            return &A1ph2SM2q1_eval<0,4,1,2,3>;
            break;
        case 0x2A344: //2367: //(qb-,ph,q+,s,s)
            return &A1ph2SM2q1_eval<1,4,0,2,3>;
            break;
        case 0x23A44: //2367: //(qb-,q+,ph,s,s)
            return &A1ph2SM2q1_eval<2,4,0,1,3>;
            break;
        case 0x234A4: //2367: //(qb-,q+,s,ph,s)
            return &A1ph2SM2q1_eval<3,4,0,1,2>;
            break;
        case 0x2344A: //2367: //(qb-,q+,s,s,ph)
            return &A1ph2SM2q1_eval<4,3,0,1,2>;
            break;
            
            
        case 0xA4324: //2127: //(ph,s,q+,qb-,s)
            return &A1ph2SM2q2_eval<0,1,2,3,4>;
            break;
        case 0x4A324: //2127: //(s,ph,q+,qb-,s)
            return &A1ph2SM2q2_eval<1,0,2,3,4>;
            break;
        case 0x43A24: //2127: //(s,q+,ph,qb-,s)
            return &A1ph2SM2q2_eval<2,0,1,3,4>;
            break;
        case 0x432A4: //2127: //(s,q+,qb-,ph,s)
            return &A1ph2SM2q2_eval<3,0,1,2,4>;
            break;
        case 0x4324A: //2127: //(s,q+,qb-,s,ph)
            return &A1ph2SM2q2_eval<4,0,1,2,3>;
            break;
            
        case 0xA4432: //48A: //(ph,s,s,q+,qb-)
            return &A1ph2SM2q2_eval<0,2,3,4,1>;
            break;
        case 0x4A432: //48A: //(s,ph,s,q+,qb-)
            return &A1ph2SM2q2_eval<1,2,3,4,0>;
            break;
        case 0x44A32: //48A: //(s,s,ph,q+,qb-)
            return &A1ph2SM2q2_eval<2,1,3,4,0>;
            break;
        case 0x443A2: //48A: //(s,s,q+,ph,qb-)
            return &A1ph2SM2q2_eval<3,1,2,4,0>;
            break;
        case 0x4432A: //48A: //(s,s,q+,qb-,ph)
            return &A1ph2SM2q2_eval<4,1,2,3,0>;
            break;
            
        case 0xA2443: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<0,3,4,1,2>;
            break;
        case 0x2A443: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<1,3,4,0,2>;
            break;
        case 0x24A43: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<2,3,4,0,1>;
            break;
        case 0x244A3: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<3,2,4,0,1>;
            break;
        case 0x2443A: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<4,2,3,0,1>;
            break;
            
        case 0xA3244: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<0,4,1,2,3>;
            break;
        case 0x3A244: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<1,4,0,2,3>;
            break;
        case 0x32A44: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<2,4,0,1,3>;
            break;
        case 0x324A4: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<3,4,0,1,2>;
            break;
        case 0x3244A: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<4,3,0,1,2>;
            break;
            
			
			//
			//H
			//
			
			
		case 0xB4234: //2169: //(ph,s,qb-,q+,s)
            return &A1ph2SM2q1_eval<0,1,2,3,4>;
            break;
        case 0x4B234: //2169: //(s,ph,qb-,q+,s)
            return &A1ph2SM2q1_eval<1,0,2,3,4>;
            break;
        case 0x42B34: //2169: //(s,qb-,ph,q+,s)
            return &A1ph2SM2q1_eval<2,0,1,3,4>;
            break;
        case 0x423B4: //2169: //(s,qb-,q+,ph,s)
            return &A1ph2SM2q1_eval<3,0,1,2,4>;
            break;
        case 0x4234B: //2169: //(s,qb-,q+,s,ph)
            return &A1ph2SM2q1_eval<4,0,1,2,3>;
            break;
            
        case 0xB4423: //783: //(ph,s,s,qb-,q+)
            return &A1ph2SM2q1_eval<0,2,3,4,1>;
            break;
        case 0x4B423: //783: //(s,ph,s,qb-,q+)
            return &A1ph2SM2q1_eval<1,2,3,4,0>;
            break;
        case 0x44B23: //783: //(s,s,ph,qb-,q+)
            return &A1ph2SM2q1_eval<2,1,3,4,0>;
            break;
        case 0x442B3: //783: //(s,s,qb-,ph,q+)
            return &A1ph2SM2q1_eval<3,1,2,4,0>;
            break;
        case 0x4423B: //783: //(ph,s,s,qb-,q+,ph)
            return &A1ph2SM2q1_eval<4,1,2,3,0>;
            break;
            
        case 0xB3442: //681: //(ph,q+,s,s,qb-)
            return &A1ph2SM2q1_eval<0,3,4,1,2>;
            break;
        case 0x3B442: //681: //(q+,ph,s,s,qb-)
            return &A1ph2SM2q1_eval<1,3,4,0,2>;
            break;
        case 0x34B42: //681: //(q+,s,ph,s,qb-)
            return &A1ph2SM2q1_eval<2,3,4,0,1>;
            break;
        case 0x344B2: //681: //(q+,s,s,ph,qb-)
            return &A1ph2SM2q1_eval<3,2,4,0,1>;
            break;
        case 0x3442B: //681: //(q+,s,s,qb-,ph)
            return &A1ph2SM2q1_eval<4,2,3,0,1>;
            break;
            
        case 0xB2344: //2367: //(ph,qb-,q+,s,s)
            return &A1ph2SM2q1_eval<0,4,1,2,3>;
            break;
        case 0x2B344: //2367: //(qb-,ph,q+,s,s)
            return &A1ph2SM2q1_eval<1,4,0,2,3>;
            break;
        case 0x23B44: //2367: //(qb-,q+,ph,s,s)
            return &A1ph2SM2q1_eval<2,4,0,1,3>;
            break;
        case 0x234B4: //2367: //(qb-,q+,s,ph,s)
            return &A1ph2SM2q1_eval<3,4,0,1,2>;
            break;
        case 0x2344B: //2367: //(qb-,q+,s,s,ph)
            return &A1ph2SM2q1_eval<4,3,0,1,2>;
            break;
            
            
        case 0xB4324: //2127: //(ph,s,q+,qb-,s)
            return &A1ph2SM2q2_eval<0,1,2,3,4>;
            break;
        case 0x4B324: //2127: //(s,ph,q+,qb-,s)
            return &A1ph2SM2q2_eval<1,0,2,3,4>;
            break;
        case 0x43B24: //2127: //(s,q+,ph,qb-,s)
            return &A1ph2SM2q2_eval<2,0,1,3,4>;
            break;
        case 0x432B4: //2127: //(s,q+,qb-,ph,s)
            return &A1ph2SM2q2_eval<3,0,1,2,4>;
            break;
        case 0x4324B: //2127: //(s,q+,qb-,s,ph)
            return &A1ph2SM2q2_eval<4,0,1,2,3>;
            break;
            
        case 0xB4432: //489: //(ph,s,s,q+,qb-)
            return &A1ph2SM2q2_eval<0,2,3,4,1>;
            break;
        case 0x4B432: //489: //(s,ph,s,q+,qb-)
            return &A1ph2SM2q2_eval<1,2,3,4,0>;
            break;
        case 0x44B32: //489: //(s,s,ph,q+,qb-)
            return &A1ph2SM2q2_eval<2,1,3,4,0>;
            break;
        case 0x443B2: //489: //(s,s,q+,ph,qb-)
            return &A1ph2SM2q2_eval<3,1,2,4,0>;
            break;
        case 0x4432B: //489: //(s,s,q+,qb-,ph)
            return &A1ph2SM2q2_eval<4,1,2,3,0>;
            break;
            
        case 0xB2443: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<0,3,4,1,2>;
            break;
        case 0x2B443: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<1,3,4,0,2>;
            break;
        case 0x24B43: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<2,3,4,0,1>;
            break;
        case 0x244B3: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<3,2,4,0,1>;
            break;
        case 0x2443B: //1023: //(qb-,s,s,q+)
            return &A1ph2SM2q2_eval<4,2,3,0,1>;
            break;
            
        case 0xB3244: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<0,4,1,2,3>;
            break;
        case 0x3B244: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<1,4,0,2,3>;
            break;
        case 0x32B44: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<2,4,0,1,3>;
            break;
        case 0x324B4: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<3,4,0,1,2>;
            break;
        case 0x3244B: //2361: //(q+,qb-,s,s)
            return &A1ph2SM2q2_eval<4,3,0,1,2>;
            break;
                        
//        case 0xB4234: //216A: //(ph,s,qb-,q+,s)
//            return &A1ph2SM2q1_eval<0,1,2,3,4>;
//            break;
//        case 0x4B234: //216A: //(s,ph,qb-,q+,s)
//            return &A1ph2SM2q1_eval<1,0,2,3,4>;
//            break;
//        case 0x42B34: //216A: //(s,qb-,ph,q+,s)
//            return &A1ph2SM2q1_eval<2,0,1,3,4>;
//            break;
//        case 0x423B4: //216A: //(s,qb-,q+,ph,s)
//            return &A1ph2SM2q1_eval<3,0,1,2,4>;
//            break;
//        case 0x4234B: //216A: //(s,qb-,q+,s,ph)
//            return &A1ph2SM2q1_eval<4,0,1,2,3>;
//            break;
//            
//        case 0xB4423: //783: //(ph,s,s,qb-,q+)
//            return &A1ph2SM2q1_eval<0,2,3,4,1>;
//            break;
//        case 0x4B423: //783: //(s,ph,s,qb-,q+)
//            return &A1ph2SM2q1_eval<1,2,3,4,0>;
//            break;
//        case 0x44B23: //783: //(s,s,ph,qb-,q+)
//            return &A1ph2SM2q1_eval<2,1,3,4,0>;
//            break;
//        case 0x442B3: //783: //(s,s,qb-,ph,q+)
//            return &A1ph2SM2q1_eval<3,1,2,4,0>;
//            break;
//        case 0x4423B: //783: //(ph,s,s,qb-,q+,ph)
//            return &A1ph2SM2q1_eval<4,1,2,3,0>;
//            break;
//            
//        case 0xB3442: //681: //(ph,q+,s,s,qb-)
//            return &A1ph2SM2q1_eval<0,3,4,1,2>;
//            break;
//        case 0x3B442: //681: //(q+,ph,s,s,qb-)
//            return &A1ph2SM2q1_eval<1,3,4,0,2>;
//            break;
//        case 0x34B42: //681: //(q+,s,ph,s,qb-)
//            return &A1ph2SM2q1_eval<2,3,4,0,1>;
//            break;
//        case 0x344B2: //681: //(q+,s,s,ph,qb-)
//            return &A1ph2SM2q1_eval<3,2,4,0,1>;
//            break;
//        case 0x3442B: //681: //(q+,s,s,qb-,ph)
//            return &A1ph2SM2q1_eval<4,2,3,0,1>;
//            break;
//            
//        case 0xB2344: //2367: //(ph,qb-,q+,s,s)
//            return &A1ph2SM2q1_eval<0,4,1,2,3>;
//            break;
//        case 0x2B344: //2367: //(qb-,ph,q+,s,s)
//            return &A1ph2SM2q1_eval<1,4,0,2,3>;
//            break;
//        case 0x23B44: //2367: //(qb-,q+,ph,s,s)
//            return &A1ph2SM2q1_eval<2,4,0,1,3>;
//            break;
//        case 0x234B4: //2367: //(qb-,q+,s,ph,s)
//            return &A1ph2SM2q1_eval<3,4,0,1,2>;
//            break;
//        case 0x2344B: //2367: //(qb-,q+,s,s,ph)
//            return &A1ph2SM2q1_eval<4,3,0,1,2>;
//            break;
//            
//            
//        case 0xB4324: //2127: //(ph,s,q+,qb-,s)
//            return &A1ph2SM2q2_eval<0,1,2,3,4>;
//            break;
//        case 0x4B324: //2127: //(s,ph,q+,qb-,s)
//            return &A1ph2SM2q2_eval<1,0,2,3,4>;
//            break;
//        case 0x43B24: //2127: //(s,q+,ph,qb-,s)
//            return &A1ph2SM2q2_eval<2,0,1,3,4>;
//            break;
//        case 0x432B4: //2127: //(s,q+,qb-,ph,s)
//            return &A1ph2SM2q2_eval<3,0,1,2,4>;
//            break;
//        case 0x4324B: //2127: //(s,q+,qb-,s,ph)
//            return &A1ph2SM2q2_eval<4,0,1,2,3>;
//            break;
//            
//        case 0xB4432: //48A: //(ph,s,s,q+,qb-)
//            return &A1ph2SM2q2_eval<0,2,3,4,1>;
//            break;
//        case 0x4B432: //48A: //(s,ph,s,q+,qb-)
//            return &A1ph2SM2q2_eval<1,2,3,4,0>;
//            break;
//        case 0x44B32: //48A: //(s,s,ph,q+,qb-)
//            return &A1ph2SM2q2_eval<2,1,3,4,0>;
//            break;
//        case 0x443B2: //48A: //(s,s,q+,ph,qb-)
//            return &A1ph2SM2q2_eval<3,1,2,4,0>;
//            break;
//        case 0x4432B: //48A: //(s,s,q+,qb-,ph)
//            return &A1ph2SM2q2_eval<4,1,2,3,0>;
//            break;
//            
//        case 0xB2443: //1023: //(qb-,s,s,q+)
//            return &A1ph2SM2q2_eval<0,3,4,1,2>;
//            break;
//        case 0x2B443: //1023: //(qb-,s,s,q+)
//            return &A1ph2SM2q2_eval<1,3,4,0,2>;
//            break;
//        case 0x24B43: //1023: //(qb-,s,s,q+)
//            return &A1ph2SM2q2_eval<2,3,4,0,1>;
//            break;
//        case 0x244B3: //1023: //(qb-,s,s,q+)
//            return &A1ph2SM2q2_eval<3,2,4,0,1>;
//            break;
//        case 0x2443B: //1023: //(qb-,s,s,q+)
//            return &A1ph2SM2q2_eval<4,2,3,0,1>;
//            break;
//            
//        case 0xB3244: //2361: //(q+,qb-,s,s)
//            return &A1ph2SM2q2_eval<0,4,1,2,3>;
//            break;
//        case 0x3B244: //2361: //(q+,qb-,s,s)
//            return &A1ph2SM2q2_eval<1,4,0,2,3>;
//            break;
//        case 0x32B44: //2361: //(q+,qb-,s,s)
//            return &A1ph2SM2q2_eval<2,4,0,1,3>;
//            break;
//        case 0x324B4: //2361: //(q+,qb-,s,s)
//            return &A1ph2SM2q2_eval<3,4,0,1,2>;
//            break;
//        case 0x3244B: //2361: //(q+,qb-,s,s)
//            return &A1ph2SM2q2_eval<4,3,0,1,2>;
//            break;
//            
            
        default:// We return zero for all other helicity combinations
            cout << "Missing helcode:" << hex << hc << dec << endl;
            return &ZeroF;
    }
}

                             
                             
template complex<R> (*A2sc1g_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2sc1g_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2sc1g_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2sc2g_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2sc2g_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2sc2g_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

  //template complex<R> (*A2sc3g_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
  //template complex<RHP> (*A2sc3g_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
  //template complex<RVHP> (*A2sc3g_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);
	
template complex<R> (*A1ph2scm1g_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2scm1g_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2scm1g_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1ph2SM_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2SM_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2SM_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);
                             
template complex<R> (*A1ph2SM1g_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2SM1g_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2SM1g_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);
                             
template complex<R> (*A1ph1QM1SM1q_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph1QM1SM1q_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph1QM1SM1q_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1ph2SM2q_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2SM2q_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2SM2q_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

                             
#if BH_USE_GMP
    
template complex<RGMP> (*A2sc1g_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A2sc2g_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2scm1g_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2SM_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2SM1g_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph1QM1SM1q_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2SM2q_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}
