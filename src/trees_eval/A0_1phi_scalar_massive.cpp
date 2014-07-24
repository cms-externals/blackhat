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

#define _HIGGS_GLUON_FAC T(1)//T(1) // ph,s,g,s higgs coupling to the gluon

#define _HIGGS_PARAM_2 T(0)//T(12) // ph,sc,sc,q,q gluon prop
#define _HIGGS_PARAM_3 T(0)//T(12) // ph,sc,sc,q,q massive quark prop
#define _OTHER_PARAM_2 T(1)//T(1) // ph,sc,sc,q,q higgs attached to gluon prop

namespace BH {

// Defined in A0_1phi_eval.cpp
template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses);



/*
 *
 *
 * The 1 complex higgs, two massive scalars and a gluon amplitude
 *
 *
 */

#define S(i,j) (ep.p(i)->P()*ep.p(j)->P())

//phsms
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1gM1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i1));	
	const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();
	
	return complex<T>(0,-1)*(_HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->Lt()*ep.ref()->Lt())))*(
        (ep.p(i2)->L()*ep.p(i1)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.p(i2)->Lt())
        -(ep.p(i2)->L()*ep.p(i3)->Sm()*ep.ref()->Lt())*(ep.p(i2)->L()*ep.p(i1)->Sm()*ep.p(i2)->Lt()))
        );
}
    
	
//phdsps
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph2sc1gM3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass=eval_param<T>::mass2(masses.p(i1));        
	const complex<T> den=(ep.p(i1)->P()+ep.p(i3)->P()).square();
    
	return complex<T>(0,1)*(
            _HIGGS_GLUON_FAC*(T(1)/(den*(ep.p(i2)->L()*ep.ref()->L())))*(
        (ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.p(i2)->L())
        -(ep.p(i2)->Lt()*ep.p(i3)->Sm()*ep.ref()->L())*(ep.p(i2)->Lt()*ep.p(i1)->Sm()*ep.p(i2)->L()))
     );
}


	


template <class T> complex<T> (*A1ph2sc1gM_Tree_Ptr_eval(long long hc))(const eval_param<T>&, const mass_param_coll&) {
	//	cout << hex << hc << dec << endl;
	switch (hc) {
			
			//
			//ph
			//
			
		case 0x9440://phssg
			return &A1ph2sc1gM1_eval<0,2,3,1>;
			break;
		case 0x4940://sphsg
			return &A1ph2sc1gM1_eval<1,2,3,0>;
			break;
		case 0x4490://ssphg
			return &A1ph2sc1gM1_eval<2,1,3,0>;
			break;
		case 0x4409://ssgph
			return &A1ph2sc1gM1_eval<3,1,2,0>;
			break;
			
//		case 0x9441://phssg
//			return &A1ph2sc1gM2_eval<0,2,3,1>;
//			break;
//		case 0x4941://sphsg
//			return &A1ph2sc1gM2_eval<1,2,3,0>;
//			break;
//		case 0x4491://ssphg
//			return &A1ph2sc1gM2_eval<2,1,3,0>;
//			break;
//		case 0x4419://ssgph
//			return &A1ph2sc1gM2_eval<3,1,2,0>;
//			break;
			
			
		case 0x9404://phssg
			return &A1ph2sc1gM1_eval<0,1,2,3>;
			break;
		case 0x4904://sphsg
			return &A1ph2sc1gM1_eval<1,0,2,3>;
			break;
		case 0x4094://ssphg
			return &A1ph2sc1gM1_eval<2,0,1,3>;
			break;
		case 0x4049://ssgph
			return &A1ph2sc1gM1_eval<3,0,1,2>;
			break;
			
//		case 0x9414://phssg
//			return &A1ph2sc1gM2_eval<0,1,2,3>;
//			break;
//		case 0x4914://sphsg
//			return &A1ph2sc1gM2_eval<1,0,2,3>;
//			break;
//		case 0x4194://ssphg
//			return &A1ph2sc1gM2_eval<2,0,1,3>;
//			break;
//		case 0x4149://ssgph
//			return &A1ph2sc1gM2_eval<3,0,1,2>;
//			break;
			
			
		case 0x9044://phssg
			return &A1ph2sc1gM1_eval<0,3,1,2>;
			break;
		case 0x0944://sphsg
			return &A1ph2sc1gM1_eval<1,3,0,2>;
			break;
		case 0x0494://ssphg
			return &A1ph2sc1gM1_eval<2,3,0,1>;
			break;
		case 0x0449://ssgph
			return &A1ph2sc1gM1_eval<3,2,0,1>;
			break;
			
//		case 0x9144://phssg
//			return &A1ph2sc1gM2_eval<0,3,1,2>;
//			break;
//		case 0x1944://sphsg
//			return &A1ph2sc1gM2_eval<1,3,0,2>;
//			break;
//		case 0x1494://ssphg
//			return &A1ph2sc1gM2_eval<2,3,0,1>;
//			break;
//		case 0x1449://ssgph
//			return &A1ph2sc1gM2_eval<3,2,0,1>;
//			break;
                      
			
			
			//
			// phd
			//
			
//		case 0xA440://phdssg
//			return &A1ph2sc1gM4_eval<0,2,3,1>;
//			break;
//		case 0x4A40://spdhsg
//			return &A1ph2sc1gM4_eval<1,2,3,0>;
//			break;
//		case 0x44A0://ssphdg
//			return &A1ph2sc1gM4_eval<2,1,3,0>;
//			break;
//		case 0x440A://ssgphd
//			return &A1ph2sc1gM4_eval<3,1,2,0>;
//			break;
			
		case 0xA441://phdssg
			return &A1ph2sc1gM3_eval<0,2,3,1>;
			break;
		case 0x4A41://sphdsg
			return &A1ph2sc1gM3_eval<1,2,3,0>;
			break;
		case 0x44A1://ssphdg
			return &A1ph2sc1gM3_eval<2,1,3,0>;
			break;
		case 0x441A://ssgphd
			return &A1ph2sc1gM3_eval<3,1,2,0>;
			break;
			
			
//		case 0xA404://phdssg
//			return &A1ph2sc1gM4_eval<0,1,2,3>;
//			break;
//		case 0x4A04://sphdsg
//			return &A1ph2sc1gM4_eval<1,0,2,3>;
//			break;
//		case 0x40A4://ssphdg
//			return &A1ph2sc1gM4_eval<2,0,1,3>;
//			break;
//		case 0x404A://ssgphd
//			return &A1ph2sc1gM4_eval<3,0,1,2>;
//			break;
			
		case 0xA414://phdssg
			return &A1ph2sc1gM3_eval<0,1,2,3>;
			break;
		case 0x4A14://sphdsg
			return &A1ph2sc1gM3_eval<1,0,2,3>;
			break;
		case 0x41A4://ssphdg
			return &A1ph2sc1gM3_eval<2,0,1,3>;
			break;
		case 0x414A://ssgphd
			return &A1ph2sc1gM3_eval<3,0,1,2>;
			break;
			
			
//		case 0xA044://phdssg
//			return &A1ph2sc1gM4_eval<0,3,1,2>;
//			break;
//		case 0x0A44://sphdsg
//			return &A1ph2sc1gM4_eval<1,3,0,2>;
//			break;
//		case 0x04A4://ssphdg
//			return &A1ph2sc1gM4_eval<2,3,0,1>;
//			break;
//		case 0x044A://ssgphd
//			return &A1ph2sc1gM4_eval<3,2,0,1>;
//			break;
			
		case 0xA144://phdssg
			return &A1ph2sc1gM3_eval<0,3,1,2>;
			break;
		case 0x1A44://sphdsg
			return &A1ph2sc1gM3_eval<1,3,0,2>;
			break;
		case 0x14A4://ssphdg
			return &A1ph2sc1gM3_eval<2,3,0,1>;
			break;
		case 0x144A://ssgpdh
			return &A1ph2sc1gM3_eval<3,2,0,1>;
			break;
			
			
		default:// We return zero for all other helicity combinations
			_MESSAGE2("Missed",hc);
			return &ZeroF;
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
template <int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A1ph2sc2qM1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

//    _MESSAGE("1:");
//    _PRINT(ep.spab(i3,i2,i4)*(sI15*(_HIGGS_PARAM_2)));
//    _PRINT(ep.spab(i3,i2,i4)*(sI23*(_HIGGS_PARAM_3)));
//    _PRINT(ep.spab(i3,i2,i4)*(_OTHER_PARAM_2*(sum34*ep.mom(i5))*sI34*sI25));
//           
//    _PRINT(-ep.spab(i3,i5,i4)*(sI12*(_HIGGS_PARAM_2)));
//    _PRINT(-ep.spab(i3,i5,i4)*(sI45*(_HIGGS_PARAM_3)));
//    _PRINT(-ep.spab(i3,i5,i4)*(_OTHER_PARAM_2*(sum34*ep.mom(i2))*sI34*sI25));    
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
template <int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A1ph2sc2qM2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
    
//    _MESSAGE("2:");
//    _PRINT(ep.spba(i3,i2,i4)*(sI15*(_HIGGS_PARAM_2)));
//    _PRINT(ep.spba(i3,i2,i4)*(sI23*(_HIGGS_PARAM_3)));
//    _PRINT(ep.spba(i3,i2,i4)*(_OTHER_PARAM_2*(sum34*ep.mom(i5))*sI34*sI25));
//    
//    _PRINT(-ep.spba(i3,i5,i4)*(sI12*(_HIGGS_PARAM_2)));
//    _PRINT(-ep.spba(i3,i5,i4)*(sI45*(_HIGGS_PARAM_3)));
//    _PRINT(-ep.spba(i3,i5,i4)*(_OTHER_PARAM_2*(sum34*ep.mom(i2))*sI34*sI25));
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

template <class T> complex<T> (*A1ph2sc2qM_Tree_Ptr_eval(long long hc))(const eval_param<T>&, const mass_param_coll&)
{
    switch (hc) {
		case 0x94234: //2169: //(ph,s,qb-,q+,s)
			return &A1ph2sc2qM1_eval<0,1,2,3,4>;
			break;
		case 0x49234: //2169: //(s,ph,qb-,q+,s)
			return &A1ph2sc2qM1_eval<1,0,2,3,4>;
			break;
		case 0x42934: //2169: //(s,qb-,ph,q+,s)
			return &A1ph2sc2qM1_eval<2,0,1,3,4>;
			break;
		case 0x42394: //2169: //(s,qb-,q+,ph,s)
			return &A1ph2sc2qM1_eval<3,0,1,2,4>;
			break;
		case 0x42349: //2169: //(s,qb-,q+,s,ph)
			return &A1ph2sc2qM1_eval<4,0,1,2,3>;
			break;
			
		case 0x94423: //783: //(ph,s,s,qb-,q+)
			return &A1ph2sc2qM1_eval<0,2,3,4,1>;
			break;
		case 0x49423: //783: //(s,ph,s,qb-,q+)
			return &A1ph2sc2qM1_eval<1,2,3,4,0>;
			break;
		case 0x44923: //783: //(s,s,ph,qb-,q+)
			return &A1ph2sc2qM1_eval<2,1,3,4,0>;
			break;
		case 0x44293: //783: //(s,s,qb-,ph,q+)
			return &A1ph2sc2qM1_eval<3,1,2,4,0>;
			break;
		case 0x44239: //783: //(ph,s,s,qb-,q+,ph)
			return &A1ph2sc2qM1_eval<4,1,2,3,0>;
			break;
			
		case 0x93442: //681: //(ph,q+,s,s,qb-)
			return &A1ph2sc2qM1_eval<0,3,4,1,2>;
			break;
		case 0x39442: //681: //(q+,ph,s,s,qb-)
			return &A1ph2sc2qM1_eval<1,3,4,0,2>;
			break;
		case 0x34942: //681: //(q+,s,ph,s,qb-)
			return &A1ph2sc2qM1_eval<2,3,4,0,1>;
			break;
		case 0x34492: //681: //(q+,s,s,ph,qb-)
			return &A1ph2sc2qM1_eval<3,2,4,0,1>;
			break;
		case 0x34429: //681: //(q+,s,s,qb-,ph)
			return &A1ph2sc2qM1_eval<4,2,3,0,1>;
			break;
			
		case 0x92344: //2367: //(ph,qb-,q+,s,s)
			return &A1ph2sc2qM1_eval<0,4,1,2,3>;
			break;
		case 0x29344: //2367: //(qb-,ph,q+,s,s)
			return &A1ph2sc2qM1_eval<1,4,0,2,3>;
			break;
		case 0x23944: //2367: //(qb-,q+,ph,s,s)
			return &A1ph2sc2qM1_eval<2,4,0,1,3>;
			break;
		case 0x23494: //2367: //(qb-,q+,s,ph,s)
			return &A1ph2sc2qM1_eval<3,4,0,1,2>;
			break;
		case 0x23449: //2367: //(qb-,q+,s,s,ph)
			return &A1ph2sc2qM1_eval<4,3,0,1,2>;
			break;
			
			
		case 0x94324: //2127: //(ph,s,q+,qb-,s)
			return &A1ph2sc2qM2_eval<0,1,2,3,4>;
			break;
		case 0x49324: //2127: //(s,ph,q+,qb-,s)
			return &A1ph2sc2qM2_eval<1,0,2,3,4>;
			break;
		case 0x43924: //2127: //(s,q+,ph,qb-,s)
			return &A1ph2sc2qM2_eval<2,0,1,3,4>;
			break;
		case 0x43294: //2127: //(s,q+,qb-,ph,s)
			return &A1ph2sc2qM2_eval<3,0,1,2,4>;
			break;
		case 0x43249: //2127: //(s,q+,qb-,s,ph)
			return &A1ph2sc2qM2_eval<4,0,1,2,3>;
			break;
			
		case 0x94432: //489: //(ph,s,s,q+,qb-)
			return &A1ph2sc2qM2_eval<0,2,3,4,1>;
			break;
		case 0x49432: //489: //(s,ph,s,q+,qb-)
			return &A1ph2sc2qM2_eval<1,2,3,4,0>;
			break;
		case 0x44932: //489: //(s,s,ph,q+,qb-)
			return &A1ph2sc2qM2_eval<2,1,3,4,0>;
			break;
		case 0x44392: //489: //(s,s,q+,ph,qb-)
			return &A1ph2sc2qM2_eval<3,1,2,4,0>;
			break;
		case 0x44329: //489: //(s,s,q+,qb-,ph)
			return &A1ph2sc2qM2_eval<4,1,2,3,0>;
			break;
			
		case 0x92443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<0,3,4,1,2>;
			break;
		case 0x29443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<1,3,4,0,2>;
			break;
		case 0x24943: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<2,3,4,0,1>;
			break;
		case 0x24493: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<3,2,4,0,1>;
			break;
		case 0x24439: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<4,2,3,0,1>;
			break;
			
		case 0x93244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<0,4,1,2,3>;
			break;
		case 0x39244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<1,4,0,2,3>;
			break;
		case 0x32944: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<2,4,0,1,3>;
			break;
		case 0x32494: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<3,4,0,1,2>;
			break;
		case 0x32449: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<4,3,0,1,2>;
			break;
			
			
			//
			//phd
			//
			
		case 0xA4234: //216A: //(ph,s,qb-,q+,s)
			return &A1ph2sc2qM1_eval<0,1,2,3,4>;
			break;
		case 0x4A234: //216A: //(s,ph,qb-,q+,s)
			return &A1ph2sc2qM1_eval<1,0,2,3,4>;
			break;
		case 0x42A34: //216A: //(s,qb-,ph,q+,s)
			return &A1ph2sc2qM1_eval<2,0,1,3,4>;
			break;
		case 0x423A4: //216A: //(s,qb-,q+,ph,s)
			return &A1ph2sc2qM1_eval<3,0,1,2,4>;
			break;
		case 0x4234A: //216A: //(s,qb-,q+,s,ph)
			return &A1ph2sc2qM1_eval<4,0,1,2,3>;
			break;
			
		case 0xA4423: //783: //(ph,s,s,qb-,q+)
			return &A1ph2sc2qM1_eval<0,2,3,4,1>;
			break;
		case 0x4A423: //783: //(s,ph,s,qb-,q+)
			return &A1ph2sc2qM1_eval<1,2,3,4,0>;
			break;
		case 0x44A23: //783: //(s,s,ph,qb-,q+)
			return &A1ph2sc2qM1_eval<2,1,3,4,0>;
			break;
		case 0x442A3: //783: //(s,s,qb-,ph,q+)
			return &A1ph2sc2qM1_eval<3,1,2,4,0>;
			break;
		case 0x4423A: //783: //(ph,s,s,qb-,q+,ph)
			return &A1ph2sc2qM1_eval<4,1,2,3,0>;
			break;
			
		case 0xA3442: //681: //(ph,q+,s,s,qb-)
			return &A1ph2sc2qM1_eval<0,3,4,1,2>;
			break;
		case 0x3A442: //681: //(q+,ph,s,s,qb-)
			return &A1ph2sc2qM1_eval<1,3,4,0,2>;
			break;
		case 0x34A42: //681: //(q+,s,ph,s,qb-)
			return &A1ph2sc2qM1_eval<2,3,4,0,1>;
			break;
		case 0x344A2: //681: //(q+,s,s,ph,qb-)
			return &A1ph2sc2qM1_eval<3,2,4,0,1>;
			break;
		case 0x3442A: //681: //(q+,s,s,qb-,ph)
			return &A1ph2sc2qM1_eval<4,2,3,0,1>;
			break;
			
		case 0xA2344: //2367: //(ph,qb-,q+,s,s)
			return &A1ph2sc2qM1_eval<0,4,1,2,3>;
			break;
		case 0x2A344: //2367: //(qb-,ph,q+,s,s)
			return &A1ph2sc2qM1_eval<1,4,0,2,3>;
			break;
		case 0x23A44: //2367: //(qb-,q+,ph,s,s)
			return &A1ph2sc2qM1_eval<2,4,0,1,3>;
			break;
		case 0x234A4: //2367: //(qb-,q+,s,ph,s)
			return &A1ph2sc2qM1_eval<3,4,0,1,2>;
			break;
		case 0x2344A: //2367: //(qb-,q+,s,s,ph)
			return &A1ph2sc2qM1_eval<4,3,0,1,2>;
			break;
			
			
		case 0xA4324: //2127: //(ph,s,q+,qb-,s)
			return &A1ph2sc2qM2_eval<0,1,2,3,4>;
			break;
		case 0x4A324: //2127: //(s,ph,q+,qb-,s)
			return &A1ph2sc2qM2_eval<1,0,2,3,4>;
			break;
		case 0x43A24: //2127: //(s,q+,ph,qb-,s)
			return &A1ph2sc2qM2_eval<2,0,1,3,4>;
			break;
		case 0x432A4: //2127: //(s,q+,qb-,ph,s)
			return &A1ph2sc2qM2_eval<3,0,1,2,4>;
			break;
		case 0x4324A: //2127: //(s,q+,qb-,s,ph)
			return &A1ph2sc2qM2_eval<4,0,1,2,3>;
			break;
			
		case 0xA4432: //48A: //(ph,s,s,q+,qb-)
			return &A1ph2sc2qM2_eval<0,2,3,4,1>;
			break;
		case 0x4A432: //48A: //(s,ph,s,q+,qb-)
			return &A1ph2sc2qM2_eval<1,2,3,4,0>;
			break;
		case 0x44A32: //48A: //(s,s,ph,q+,qb-)
			return &A1ph2sc2qM2_eval<2,1,3,4,0>;
			break;
		case 0x443A2: //48A: //(s,s,q+,ph,qb-)
			return &A1ph2sc2qM2_eval<3,1,2,4,0>;
			break;
		case 0x4432A: //48A: //(s,s,q+,qb-,ph)
			return &A1ph2sc2qM2_eval<4,1,2,3,0>;
			break;
			
		case 0xA2443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<0,3,4,1,2>;
			break;
		case 0x2A443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<1,3,4,0,2>;
			break;
		case 0x24A43: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<2,3,4,0,1>;
			break;
		case 0x244A3: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<3,2,4,0,1>;
			break;
		case 0x2443A: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2qM2_eval<4,2,3,0,1>;
			break;
			
		case 0xA3244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<0,4,1,2,3>;
			break;
		case 0x3A244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<1,4,0,2,3>;
			break;
		case 0x32A44: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<2,4,0,1,3>;
			break;
		case 0x324A4: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<3,4,0,1,2>;
			break;
		case 0x3244A: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2qM2_eval<4,3,0,1,2>;
			break;
			
			
			
		default:// We return zero for all other helicity combinations
            cout << "Missing helcode:" << hex << hc << dec << endl;
			return &ZeroF;
	}
}
	
template complex<R> (*A1ph2sc1gM_Tree_Ptr_eval(long long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2sc1gM_Tree_Ptr_eval(long long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2sc1gM_Tree_Ptr_eval(long long hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1ph2sc2qM_Tree_Ptr_eval(long long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2sc2qM_Tree_Ptr_eval(long long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2sc2qM_Tree_Ptr_eval(long long hc))(const eval_param<RVHP>&, const mass_param_coll&);
			
#if BH_USE_GMP
template complex<RGMP> (*A1ph2sc1gM_Tree_Ptr_eval(long long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2sc2qM_Tree_Ptr_eval(long long hc))(const eval_param<RGMP>&, const mass_param_coll&);

#endif

}	


