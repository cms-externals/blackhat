//
//  A0_1phi_scalar_2.cpp
//  BlackHat
//
//  Created by Darren Forde on 20/01/2011.
//  Copyright 2011 BlackHat Collaboration. All rights reserved.
//

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;

#define _HIGGS_SCALAR_SCALAR_PARAM_3 T(1)//T(-12) // ph,sc,sc
#define _HIGGS_SCALAR_SCALAR_PARAM_1 T(_HIGGS_SCALAR_SCALAR_PARAM_3) // T(12) ph,sc,q,Q

#define _HIGGS_PARAM_2 T(0) // ph,sc,sc,q,q gluon prop
#define _HIGGS_PARAM_3 T(0) // ph,sc,sc,q,q massive quark prop
#define _OTHER_PARAM_2 T(1) // ph,sc,sc,q,q higgs attached to gluon prop

namespace BH {

// Defined in A0_1phi_eval.cpp
template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses);


			
			
/*
 *
 *
 * The 1 complex higgs, 1 massive quark a massive scalar and a massless quark amplitude
 *
 *
 */

//(s,q-,Q+,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	
	const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 1m:",complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->L()*qL)/(den*(QL*qL)*sqrt(T(2))));
	
	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->L()*qL)/(den*(QL*qL)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	
	const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 1p:",complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->L()*qL)/(den*(QL*qL)*sqrt(T(2))));
	
	return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->L()*qL)/(den*(QL*qL)*sqrt(T(2)));
}

//(s,q-,Q-,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const lambda<T> QL=la(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	
	const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 2:",complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(ep.p(i1)->L()*QL)/(den*sqrt(T(2))));
	
	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(ep.p(i1)->L()*QL)/(den*sqrt(T(2)));
}

//(Q-,q-,s,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	
	const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 3:",complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(QL*ep.p(i1)->L())/(den*sqrt(T(2))));
	
	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(QL*ep.p(i1)->L())/(den*sqrt(T(2)));
}

//(Q+,q-,s,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	
	const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 4m:",complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qL*ep.p(i1)->L())/(den*(qL*QL)*sqrt(T(2))));
	
	return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qL*ep.p(i1)->L())/(den*(qL*QL)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const lambda<T> QL=la(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();
	
	const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 4p:",complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qL*ep.p(i1)->L())/(den*(qL*QL)*sqrt(T(2))));
	
	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qL*ep.p(i1)->L())/(den*(qL*QL)*sqrt(T(2)));
}

//(Q-,q+,s,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q5m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	
	const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 5m:",complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qLt*ep.p(i1)->Lt())/(den*(qLt*QLt)*sqrt(T(2))));
	
	return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qLt*ep.p(i1)->Lt())/(den*(qLt*QLt)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q5p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	
	const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 5p:",complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qLt*ep.p(i1)->Lt())/(den*(qLt*QLt)*sqrt(T(2))));
	
	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(qLt*ep.p(i1)->Lt())/(den*(qLt*QLt)*sqrt(T(2)));
}

//(Q+,q+,s,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const lambdat<T> QLt=lat(ep.mom(i0)-T(0.5)*(mass2/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	
	const complex<T> den=(ep.p(i0)->P()+ep.p(i1)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 6:",complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(QLt*ep.p(i1)->Lt())/(den*sqrt(T(2))));
	
	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(QLt*ep.p(i1)->Lt())/(den*sqrt(T(2)));
}

//(s,q+,Q+,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q7_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	
	const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 7:",complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(ep.p(i1)->Lt()*QLt)/(den*sqrt(T(2))));
	
	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*(ep.p(i1)->Lt()*QLt)/(den*sqrt(T(2)));
}

//(s,q+,Q-,ph)
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	
	const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 8m:",complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->Lt()*qLt)/(den*(QLt*qLt)*sqrt(T(2))));
	
	return complex<T>(0,-1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->Lt()*qLt)/(den*(QLt*qLt)*sqrt(T(2)));
}
template <int i0, int i1, int i2, class T> complex<T>  A1ph1QM1sc1q8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	const complex<T> mass2=eval_param<T>::mass2(masses.p(i2));
	const complex<T> mass=eval_param<T>::mass(masses.p(i2));
	const lambdat<T> QLt=lat(ep.mom(i2)-T(0.5)*(mass2/(ep.p(i2)->P()*ep.ref()->P()))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();
	
	const complex<T> den=(ep.p(i1)->P()+ep.p(i2)->P()).square()-mass2;
	
	//	_MESSAGE2("Case 8p:",complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->Lt()*qLt)/(den*(QLt*qLt)*sqrt(T(2))));
	
	return complex<T>(0,1)*_HIGGS_SCALAR_SCALAR_PARAM_1*mass2*mass*(ep.p(i1)->Lt()*qLt)/(den*(QLt*qLt)*sqrt(T(2)));
}


template <class T, int SPOS, int QPOS, int QMPOS> complex<T> (*A1ph1QM1sc1q_phi_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	//	cout << hex << hc << dec << endl;
	switch (hc) {
		case 0xd28: //258: //(s,q3-,Qb+)
			return &A1ph1QM1sc1q1m_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0xd26:
			return &A1ph1QM1sc1q1p_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x28d: //330: //(qb-,Q+,s)
			return &A1ph1QM1sc1q1m_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0x26d:
			return &A1ph1QM1sc1q1p_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0x8d2: //96: //(Q+,s,qb-)
			return &A1ph1QM1sc1q1m_eval<QPOS,QMPOS,SPOS>;
			break;
		case 0x6d2:
			return &A1ph1QM1sc1q1p_eval<QPOS,QMPOS,SPOS>;
			break;
			
			
		case 0xd27: //209: //(s,q3-,Qb-)
		case 0xd25:
			return &A1ph1QM1sc1q2_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x27d: //323: //(qb-,Q-,s)
		case 0x25d:
			return &A1ph1QM1sc1q2_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0x7d2: //95: //(qb-,Q-,s)
		case 0x5d2:
			return &A1ph1QM1sc1q2_eval<QPOS,QMPOS,SPOS>;
			break;
			
			
			
		case 0x72d: //305: //(Q-,q3-,s)
		case 0x52d:
			return &A1ph1QM1sc1q3_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x2d7: //239: //(qb-,s,Q-)
		case 0x2d5:
			return &A1ph1QM1sc1q3_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0xd72: //83: //(s,Q3-,qb3-)
		case 0xd52:
			return &A1ph1QM1sc1q3_eval<QPOS,QMPOS,SPOS>;
			break;
			
			
			
		case 0x82d: //306: //(Q+,q3-,s)
			return &A1ph1QM1sc1q4m_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x62d:
			return &A1ph1QM1sc1q4p_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x2d8: //288: //(qb-,s,Q+)
			return &A1ph1QM1sc1q4m_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0x2d6:
			return &A1ph1QM1sc1q4p_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0xd82: //90: //(s,Q3+,qb3-)
			return &A1ph1QM1sc1q4m_eval<QPOS,QMPOS,SPOS>;
			break;
		case 0xd62:
			return &A1ph1QM1sc1q4p_eval<QPOS,QMPOS,SPOS>;
			break;
			
			
		case 0x73d: //312: //(Q-,qb+,s)
			return &A1ph1QM1sc1q5m_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x53d:
			return &A1ph1QM1sc1q5p_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x3d7: //240: //(q+,s,Qb-)
			return &A1ph1QM1sc1q5m_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0x3d5:
			return &A1ph1QM1sc1q5p_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0xd73: //132: //(s,Qb-,q+)
			return &A1ph1QM1sc1q5m_eval<QPOS,QMPOS,SPOS>;
			break;
		case 0xd53:
			return &A1ph1QM1sc1q5p_eval<QPOS,QMPOS,SPOS>;
			break;
			
			
		case 0x83d: //313: //(Q+,qb+,s)
		case 0x63d:
			return &A1ph1QM1sc1q6_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x3d8: //289: //(q+,s,Qb+)
		case 0x3d6:
			return &A1ph1QM1sc1q6_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0xd83: //139: //(s,Qb+,q+)
		case 0xd63:
			return &A1ph1QM1sc1q6_eval<QPOS,QMPOS,SPOS>;
			break;
			
			
		case 0xd38: //265: //(s,qb+,Qb+)
		case 0xd36:
			return &A1ph1QM1sc1q7_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x38d: //331: //(q+,Qb+,s)
		case 0x36d:
			return &A1ph1QM1sc1q7_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0x8d3: //145: //(Qb+,s,q+)
		case 0x6d3:
			return &A1ph1QM1sc1q7_eval<QPOS,QMPOS,SPOS>;
			break;
			
			
		case 0xd37: //216: //(s,qb+,Qb-)
			return &A1ph1QM1sc1q8m_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0xd35:
			return &A1ph1QM1sc1q8p_eval<SPOS,QPOS,QMPOS>;
			break;
		case 0x37d: //324: //(q+,Qb-,s)
			return &A1ph1QM1sc1q8m_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0x35d:
			return &A1ph1QM1sc1q8p_eval<QMPOS,SPOS,QPOS>;
			break;
		case 0x7d3: //144: //(Qb-,s,q+)
			return &A1ph1QM1sc1q8m_eval<QPOS,QMPOS,SPOS>;
			break;
		case 0x5d3:
			return &A1ph1QM1sc1q8p_eval<QPOS,QMPOS,SPOS>;
			break;
			
			
		default:// We return zero for all other helicity combinations
			//cout << "Missed vertex:" << hex << hc << dec << endl;
			return &ZeroF;
	}
}

template <class T> complex<T> (*A1ph1QM1sc1q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
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
			if(efac==0xA){// If this is a phi dagger then set record this
				dagger=true;
			}
			if(efac==0xE){// If this is a phi dagger then set record this
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
					return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,1,2>(newhc);
					break;
				case 1:
					return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,1,3>(newhc);
					break;
				case 2:
					return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,2,3>(newhc);
					break;
				case 3:
					return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,1,2,3>(newhc);
					break;
					
				default:// Should never get here as there is always a higgs. 
					return &ZeroF;
					break;
			}
		}
		else{
			switch (epos) {
				case 0:
					return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,1,2>(newhc);
					break;
				case 1:
					return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,1,3>(newhc);
					break;
				case 2:
					return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,2,3>(newhc);
					break;
				case 3:
					return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,1,2,3>(newhc);
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
				return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,1,2>(newhc);
				break;
			case 1:
				return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,1,3>(newhc);
				break;
			case 2:
				return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,0,2,3>(newhc);
				break;
			case 3:
				return A1ph1QM1sc1q_phi_Tree_Ptr_eval<T,1,2,3>(newhc);
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
template <int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A1ph2sc2q1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
	
	//	return complex<T>(0,0);
	
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
template <int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A1ph2sc2q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
	
	//	return complex<T>(0,0);
	
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

template <class T> complex<T> (*A1ph2sc2q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {
		case 0x94234: //2169: //(ph,s,qb-,q+,s)
			return &A1ph2sc2q1_eval<0,1,2,3,4>;
			break;
		case 0x49234: //2169: //(s,ph,qb-,q+,s)
			return &A1ph2sc2q1_eval<1,0,2,3,4>;
			break;
		case 0x42934: //2169: //(s,qb-,ph,q+,s)
			return &A1ph2sc2q1_eval<2,0,1,3,4>;
			break;
		case 0x42394: //2169: //(s,qb-,q+,ph,s)
			return &A1ph2sc2q1_eval<3,0,1,2,4>;
			break;
		case 0x42349: //2169: //(s,qb-,q+,s,ph)
			return &A1ph2sc2q1_eval<4,0,1,2,3>;
			break;
			
		case 0x94423: //783: //(ph,s,s,qb-,q+)
			return &A1ph2sc2q1_eval<0,2,3,4,1>;
			break;
		case 0x49423: //783: //(s,ph,s,qb-,q+)
			return &A1ph2sc2q1_eval<1,2,3,4,0>;
			break;
		case 0x44923: //783: //(s,s,ph,qb-,q+)
			return &A1ph2sc2q1_eval<2,1,3,4,0>;
			break;
		case 0x44293: //783: //(s,s,qb-,ph,q+)
			return &A1ph2sc2q1_eval<3,1,2,4,0>;
			break;
		case 0x44239: //783: //(ph,s,s,qb-,q+,ph)
			return &A1ph2sc2q1_eval<4,1,2,3,0>;
			break;
			
		case 0x93442: //681: //(ph,q+,s,s,qb-)
			return &A1ph2sc2q1_eval<0,3,4,1,2>;
			break;
		case 0x39442: //681: //(q+,ph,s,s,qb-)
			return &A1ph2sc2q1_eval<1,3,4,0,2>;
			break;
		case 0x34942: //681: //(q+,s,ph,s,qb-)
			return &A1ph2sc2q1_eval<2,3,4,0,1>;
			break;
		case 0x34492: //681: //(q+,s,s,ph,qb-)
			return &A1ph2sc2q1_eval<3,2,4,0,1>;
			break;
		case 0x34429: //681: //(q+,s,s,qb-,ph)
			return &A1ph2sc2q1_eval<4,2,3,0,1>;
			break;
			
		case 0x92344: //2367: //(ph,qb-,q+,s,s)
			return &A1ph2sc2q1_eval<0,4,1,2,3>;
			break;
		case 0x29344: //2367: //(qb-,ph,q+,s,s)
			return &A1ph2sc2q1_eval<1,4,0,2,3>;
			break;
		case 0x23944: //2367: //(qb-,q+,ph,s,s)
			return &A1ph2sc2q1_eval<2,4,0,1,3>;
			break;
		case 0x23494: //2367: //(qb-,q+,s,ph,s)
			return &A1ph2sc2q1_eval<3,4,0,1,2>;
			break;
		case 0x23449: //2367: //(qb-,q+,s,s,ph)
			return &A1ph2sc2q1_eval<4,3,0,1,2>;
			break;
			
			
		case 0x94324: //2127: //(ph,s,q+,qb-,s)
			return &A1ph2sc2q2_eval<0,1,2,3,4>;
			break;
		case 0x49324: //2127: //(s,ph,q+,qb-,s)
			return &A1ph2sc2q2_eval<1,0,2,3,4>;
			break;
		case 0x43924: //2127: //(s,q+,ph,qb-,s)
			return &A1ph2sc2q2_eval<2,0,1,3,4>;
			break;
		case 0x43294: //2127: //(s,q+,qb-,ph,s)
			return &A1ph2sc2q2_eval<3,0,1,2,4>;
			break;
		case 0x43249: //2127: //(s,q+,qb-,s,ph)
			return &A1ph2sc2q2_eval<4,0,1,2,3>;
			break;
			
		case 0x94432: //489: //(ph,s,s,q+,qb-)
			return &A1ph2sc2q2_eval<0,2,3,4,1>;
			break;
		case 0x49432: //489: //(s,ph,s,q+,qb-)
			return &A1ph2sc2q2_eval<1,2,3,4,0>;
			break;
		case 0x44932: //489: //(s,s,ph,q+,qb-)
			return &A1ph2sc2q2_eval<2,1,3,4,0>;
			break;
		case 0x44392: //489: //(s,s,q+,ph,qb-)
			return &A1ph2sc2q2_eval<3,1,2,4,0>;
			break;
		case 0x44329: //489: //(s,s,q+,qb-,ph)
			return &A1ph2sc2q2_eval<4,1,2,3,0>;
			break;
			
		case 0x92443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<0,3,4,1,2>;
			break;
		case 0x29443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<1,3,4,0,2>;
			break;
		case 0x24943: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<2,3,4,0,1>;
			break;
		case 0x24493: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<3,2,4,0,1>;
			break;
		case 0x24439: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<4,2,3,0,1>;
			break;
			
		case 0x93244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<0,4,1,2,3>;
			break;
		case 0x39244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<1,4,0,2,3>;
			break;
		case 0x32944: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<2,4,0,1,3>;
			break;
		case 0x32494: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<3,4,0,1,2>;
			break;
		case 0x32449: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<4,3,0,1,2>;
			break;
			
			
			//
			//phd
			//
			
		case 0xA4234: //216A: //(ph,s,qb-,q+,s)
			return &A1ph2sc2q1_eval<0,1,2,3,4>;
			break;
		case 0x4A234: //216A: //(s,ph,qb-,q+,s)
			return &A1ph2sc2q1_eval<1,0,2,3,4>;
			break;
		case 0x42A34: //216A: //(s,qb-,ph,q+,s)
			return &A1ph2sc2q1_eval<2,0,1,3,4>;
			break;
		case 0x423A4: //216A: //(s,qb-,q+,ph,s)
			return &A1ph2sc2q1_eval<3,0,1,2,4>;
			break;
		case 0x4234A: //216A: //(s,qb-,q+,s,ph)
			return &A1ph2sc2q1_eval<4,0,1,2,3>;
			break;
			
		case 0xA4423: //783: //(ph,s,s,qb-,q+)
			return &A1ph2sc2q1_eval<0,2,3,4,1>;
			break;
		case 0x4A423: //783: //(s,ph,s,qb-,q+)
			return &A1ph2sc2q1_eval<1,2,3,4,0>;
			break;
		case 0x44A23: //783: //(s,s,ph,qb-,q+)
			return &A1ph2sc2q1_eval<2,1,3,4,0>;
			break;
		case 0x442A3: //783: //(s,s,qb-,ph,q+)
			return &A1ph2sc2q1_eval<3,1,2,4,0>;
			break;
		case 0x4423A: //783: //(ph,s,s,qb-,q+,ph)
			return &A1ph2sc2q1_eval<4,1,2,3,0>;
			break;
			
		case 0xA3442: //681: //(ph,q+,s,s,qb-)
			return &A1ph2sc2q1_eval<0,3,4,1,2>;
			break;
		case 0x3A442: //681: //(q+,ph,s,s,qb-)
			return &A1ph2sc2q1_eval<1,3,4,0,2>;
			break;
		case 0x34A42: //681: //(q+,s,ph,s,qb-)
			return &A1ph2sc2q1_eval<2,3,4,0,1>;
			break;
		case 0x344A2: //681: //(q+,s,s,ph,qb-)
			return &A1ph2sc2q1_eval<3,2,4,0,1>;
			break;
		case 0x3442A: //681: //(q+,s,s,qb-,ph)
			return &A1ph2sc2q1_eval<4,2,3,0,1>;
			break;
			
		case 0xA2344: //2367: //(ph,qb-,q+,s,s)
			return &A1ph2sc2q1_eval<0,4,1,2,3>;
			break;
		case 0x2A344: //2367: //(qb-,ph,q+,s,s)
			return &A1ph2sc2q1_eval<1,4,0,2,3>;
			break;
		case 0x23A44: //2367: //(qb-,q+,ph,s,s)
			return &A1ph2sc2q1_eval<2,4,0,1,3>;
			break;
		case 0x234A4: //2367: //(qb-,q+,s,ph,s)
			return &A1ph2sc2q1_eval<3,4,0,1,2>;
			break;
		case 0x2344A: //2367: //(qb-,q+,s,s,ph)
			return &A1ph2sc2q1_eval<4,3,0,1,2>;
			break;
			
			
		case 0xA4324: //2127: //(ph,s,q+,qb-,s)
			return &A1ph2sc2q2_eval<0,1,2,3,4>;
			break;
		case 0x4A324: //2127: //(s,ph,q+,qb-,s)
			return &A1ph2sc2q2_eval<1,0,2,3,4>;
			break;
		case 0x43A24: //2127: //(s,q+,ph,qb-,s)
			return &A1ph2sc2q2_eval<2,0,1,3,4>;
			break;
		case 0x432A4: //2127: //(s,q+,qb-,ph,s)
			return &A1ph2sc2q2_eval<3,0,1,2,4>;
			break;
		case 0x4324A: //2127: //(s,q+,qb-,s,ph)
			return &A1ph2sc2q2_eval<4,0,1,2,3>;
			break;
			
		case 0xA4432: //48A: //(ph,s,s,q+,qb-)
			return &A1ph2sc2q2_eval<0,2,3,4,1>;
			break;
		case 0x4A432: //48A: //(s,ph,s,q+,qb-)
			return &A1ph2sc2q2_eval<1,2,3,4,0>;
			break;
		case 0x44A32: //48A: //(s,s,ph,q+,qb-)
			return &A1ph2sc2q2_eval<2,1,3,4,0>;
			break;
		case 0x443A2: //48A: //(s,s,q+,ph,qb-)
			return &A1ph2sc2q2_eval<3,1,2,4,0>;
			break;
		case 0x4432A: //48A: //(s,s,q+,qb-,ph)
			return &A1ph2sc2q2_eval<4,1,2,3,0>;
			break;
			
		case 0xA2443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<0,3,4,1,2>;
			break;
		case 0x2A443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<1,3,4,0,2>;
			break;
		case 0x24A43: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<2,3,4,0,1>;
			break;
		case 0x244A3: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<3,2,4,0,1>;
			break;
		case 0x2443A: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<4,2,3,0,1>;
			break;
			
		case 0xA3244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<0,4,1,2,3>;
			break;
		case 0x3A244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<1,4,0,2,3>;
			break;
		case 0x32A44: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<2,4,0,1,3>;
			break;
		case 0x324A4: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<3,4,0,1,2>;
			break;
		case 0x3244A: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<4,3,0,1,2>;
			break;
			
			
			//
			// H
			//
			
			
		case 0xB4234: //216E: //(ph,s,qb-,q+,s)
			return &A1ph2sc2q1_eval<0,1,2,3,4>;
			break;
		case 0x4B234: //216E: //(s,ph,qb-,q+,s)
			return &A1ph2sc2q1_eval<1,0,2,3,4>;
			break;
		case 0x42B34: //216E: //(s,qb-,ph,q+,s)
			return &A1ph2sc2q1_eval<2,0,1,3,4>;
			break;
		case 0x423B4: //216E: //(s,qb-,q+,ph,s)
			return &A1ph2sc2q1_eval<3,0,1,2,4>;
			break;
		case 0x4234B: //216E: //(s,qb-,q+,s,ph)
			return &A1ph2sc2q1_eval<4,0,1,2,3>;
			break;
			
		case 0xB4423: //783: //(ph,s,s,qb-,q+)
			return &A1ph2sc2q1_eval<0,2,3,4,1>;
			break;
		case 0x4B423: //783: //(s,ph,s,qb-,q+)
			return &A1ph2sc2q1_eval<1,2,3,4,0>;
			break;
		case 0x44B23: //783: //(s,s,ph,qb-,q+)
			return &A1ph2sc2q1_eval<2,1,3,4,0>;
			break;
		case 0x442B3: //783: //(s,s,qb-,ph,q+)
			return &A1ph2sc2q1_eval<3,1,2,4,0>;
			break;
		case 0x4423B: //783: //(ph,s,s,qb-,q+,ph)
			return &A1ph2sc2q1_eval<4,1,2,3,0>;
			break;
			
		case 0xB3442: //681: //(ph,q+,s,s,qb-)
			return &A1ph2sc2q1_eval<0,3,4,1,2>;
			break;
		case 0x3B442: //681: //(q+,ph,s,s,qb-)
			return &A1ph2sc2q1_eval<1,3,4,0,2>;
			break;
		case 0x34B42: //681: //(q+,s,ph,s,qb-)
			return &A1ph2sc2q1_eval<2,3,4,0,1>;
			break;
		case 0x344B2: //681: //(q+,s,s,ph,qb-)
			return &A1ph2sc2q1_eval<3,2,4,0,1>;
			break;
		case 0x3442B: //681: //(q+,s,s,qb-,ph)
			return &A1ph2sc2q1_eval<4,2,3,0,1>;
			break;
			
		case 0xB2344: //2367: //(ph,qb-,q+,s,s)
			return &A1ph2sc2q1_eval<0,4,1,2,3>;
			break;
		case 0x2B344: //2367: //(qb-,ph,q+,s,s)
			return &A1ph2sc2q1_eval<1,4,0,2,3>;
			break;
		case 0x23B44: //2367: //(qb-,q+,ph,s,s)
			return &A1ph2sc2q1_eval<2,4,0,1,3>;
			break;
		case 0x234B4: //2367: //(qb-,q+,s,ph,s)
			return &A1ph2sc2q1_eval<3,4,0,1,2>;
			break;
		case 0x2344B: //2367: //(qb-,q+,s,s,ph)
			return &A1ph2sc2q1_eval<4,3,0,1,2>;
			break;
			
			
		case 0xB4324: //2127: //(ph,s,q+,qb-,s)
			return &A1ph2sc2q2_eval<0,1,2,3,4>;
			break;
		case 0x4B324: //2127: //(s,ph,q+,qb-,s)
			return &A1ph2sc2q2_eval<1,0,2,3,4>;
			break;
		case 0x43B24: //2127: //(s,q+,ph,qb-,s)
			return &A1ph2sc2q2_eval<2,0,1,3,4>;
			break;
		case 0x432B4: //2127: //(s,q+,qb-,ph,s)
			return &A1ph2sc2q2_eval<3,0,1,2,4>;
			break;
		case 0x4324B: //2127: //(s,q+,qb-,s,ph)
			return &A1ph2sc2q2_eval<4,0,1,2,3>;
			break;
			
		case 0xB4432: //48E: //(ph,s,s,q+,qb-)
			return &A1ph2sc2q2_eval<0,2,3,4,1>;
			break;
		case 0x4B432: //48E: //(s,ph,s,q+,qb-)
			return &A1ph2sc2q2_eval<1,2,3,4,0>;
			break;
		case 0x44B32: //48E: //(s,s,ph,q+,qb-)
			return &A1ph2sc2q2_eval<2,1,3,4,0>;
			break;
		case 0x443B2: //48E: //(s,s,q+,ph,qb-)
			return &A1ph2sc2q2_eval<3,1,2,4,0>;
			break;
		case 0x4432B: //48E: //(s,s,q+,qb-,ph)
			return &A1ph2sc2q2_eval<4,1,2,3,0>;
			break;
			
		case 0xB2443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<0,3,4,1,2>;
			break;
		case 0x2B443: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<1,3,4,0,2>;
			break;
		case 0x24B43: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<2,3,4,0,1>;
			break;
		case 0x244B3: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<3,2,4,0,1>;
			break;
		case 0x2443B: //1023: //(qb-,s,s,q+)
			return &A1ph2sc2q2_eval<4,2,3,0,1>;
			break;
			
		case 0xB3244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<0,4,1,2,3>;
			break;
		case 0x3B244: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<1,4,0,2,3>;
			break;
		case 0x32B44: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<2,4,0,1,3>;
			break;
		case 0x324B4: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<3,4,0,1,2>;
			break;
		case 0x3244B: //2361: //(q+,qb-,s,s)
			return &A1ph2sc2q2_eval<4,3,0,1,2>;
			break;
			
			
			//			//
			//			//phd
			//			//
			//			
			//		case 0xE4234: //216E: //(ph,s,qb-,q+,s)
			//			return &A1ph2sc2q1_eval<0,1,2,3,4>;
			//			break;
			//		case 0x4E234: //216A: //(s,ph,qb-,q+,s)
			//			return &A1ph2sc2q1_eval<1,0,2,3,4>;
			//			break;
			//		case 0x42E34: //216A: //(s,qb-,ph,q+,s)
			//			return &A1ph2sc2q1_eval<2,0,1,3,4>;
			//			break;
			//		case 0x423E4: //216A: //(s,qb-,q+,ph,s)
			//			return &A1ph2sc2q1_eval<3,0,1,2,4>;
			//			break;
			//		case 0x4234E: //216A: //(s,qb-,q+,s,ph)
			//			return &A1ph2sc2q1_eval<4,0,1,2,3>;
			//			break;
			//			
			//		case 0xE4423: //783: //(ph,s,s,qb-,q+)
			//			return &A1ph2sc2q1_eval<0,2,3,4,1>;
			//			break;
			//		case 0x4E423: //783: //(s,ph,s,qb-,q+)
			//			return &A1ph2sc2q1_eval<1,2,3,4,0>;
			//			break;
			//		case 0x44E23: //783: //(s,s,ph,qb-,q+)
			//			return &A1ph2sc2q1_eval<2,1,3,4,0>;
			//			break;
			//		case 0x442E3: //783: //(s,s,qb-,ph,q+)
			//			return &A1ph2sc2q1_eval<3,1,2,4,0>;
			//			break;
			//		case 0x4423E: //783: //(ph,s,s,qb-,q+,ph)
			//			return &A1ph2sc2q1_eval<4,1,2,3,0>;
			//			break;
			//			
			//		case 0xE3442: //681: //(ph,q+,s,s,qb-)
			//			return &A1ph2sc2q1_eval<0,3,4,1,2>;
			//			break;
			//		case 0x3E442: //681: //(q+,ph,s,s,qb-)
			//			return &A1ph2sc2q1_eval<1,3,4,0,2>;
			//			break;
			//		case 0x34E42: //681: //(q+,s,ph,s,qb-)
			//			return &A1ph2sc2q1_eval<2,3,4,0,1>;
			//			break;
			//		case 0x344E2: //681: //(q+,s,s,ph,qb-)
			//			return &A1ph2sc2q1_eval<3,2,4,0,1>;
			//			break;
			//		case 0x3442E: //681: //(q+,s,s,qb-,ph)
			//			return &A1ph2sc2q1_eval<4,2,3,0,1>;
			//			break;
			//			
			//		case 0xE2344: //2367: //(ph,qb-,q+,s,s)
			//			return &A1ph2sc2q1_eval<0,4,1,2,3>;
			//			break;
			//		case 0x2E344: //2367: //(qb-,ph,q+,s,s)
			//			return &A1ph2sc2q1_eval<1,4,0,2,3>;
			//			break;
			//		case 0x23E44: //2367: //(qb-,q+,ph,s,s)
			//			return &A1ph2sc2q1_eval<2,4,0,1,3>;
			//			break;
			//		case 0x234E4: //2367: //(qb-,q+,s,ph,s)
			//			return &A1ph2sc2q1_eval<3,4,0,1,2>;
			//			break;
			//		case 0x2344E: //2367: //(qb-,q+,s,s,ph)
			//			return &A1ph2sc2q1_eval<4,3,0,1,2>;
			//			break;
			//			
			//			
			//		case 0xE4324: //2127: //(ph,s,q+,qb-,s)
			//			return &A1ph2sc2q2_eval<0,1,2,3,4>;
			//			break;
			//		case 0x4E324: //2127: //(s,ph,q+,qb-,s)
			//			return &A1ph2sc2q2_eval<1,0,2,3,4>;
			//			break;
			//		case 0x43E24: //2127: //(s,q+,ph,qb-,s)
			//			return &A1ph2sc2q2_eval<2,0,1,3,4>;
			//			break;
			//		case 0x432E4: //2127: //(s,q+,qb-,ph,s)
			//			return &A1ph2sc2q2_eval<3,0,1,2,4>;
			//			break;
			//		case 0x4324E: //2127: //(s,q+,qb-,s,ph)
			//			return &A1ph2sc2q2_eval<4,0,1,2,3>;
			//			break;
			//			
			//		case 0xE4432: //48A: //(ph,s,s,q+,qb-)
			//			return &A1ph2sc2q2_eval<0,2,3,4,1>;
			//			break;
			//		case 0x4E432: //48A: //(s,ph,s,q+,qb-)
			//			return &A1ph2sc2q2_eval<1,2,3,4,0>;
			//			break;
			//		case 0x44E32: //48A: //(s,s,ph,q+,qb-)
			//			return &A1ph2sc2q2_eval<2,1,3,4,0>;
			//			break;
			//		case 0x443E2: //48A: //(s,s,q+,ph,qb-)
			//			return &A1ph2sc2q2_eval<3,1,2,4,0>;
			//			break;
			//		case 0x4432E: //48A: //(s,s,q+,qb-,ph)
			//			return &A1ph2sc2q2_eval<4,1,2,3,0>;
			//			break;
			//			
			//		case 0xE2443: //1023: //(qb-,s,s,q+)
			//			return &A1ph2sc2q2_eval<0,3,4,1,2>;
			//			break;
			//		case 0x2E443: //1023: //(qb-,s,s,q+)
			//			return &A1ph2sc2q2_eval<1,3,4,0,2>;
			//			break;
			//		case 0x24E43: //1023: //(qb-,s,s,q+)
			//			return &A1ph2sc2q2_eval<2,3,4,0,1>;
			//			break;
			//		case 0x244E3: //1023: //(qb-,s,s,q+)
			//			return &A1ph2sc2q2_eval<3,2,4,0,1>;
			//			break;
			//		case 0x2443E: //1023: //(qb-,s,s,q+)
			//			return &A1ph2sc2q2_eval<4,2,3,0,1>;
			//			break;
			//			
			//		case 0xE3244: //2361: //(q+,qb-,s,s)
			//			return &A1ph2sc2q2_eval<0,4,1,2,3>;
			//			break;
			//		case 0x3E244: //2361: //(q+,qb-,s,s)
			//			return &A1ph2sc2q2_eval<1,4,0,2,3>;
			//			break;
			//		case 0x32E44: //2361: //(q+,qb-,s,s)
			//			return &A1ph2sc2q2_eval<2,4,0,1,3>;
			//			break;
			//		case 0x324E4: //2361: //(q+,qb-,s,s)
			//			return &A1ph2sc2q2_eval<3,4,0,1,2>;
			//			break;
			//		case 0x3244E: //2361: //(q+,qb-,s,s)
			//			return &A1ph2sc2q2_eval<4,3,0,1,2>;
			//			break;
			
			
			
			
		default:// We return zero for all other helicity combinations
			cout << "Missing helcode:" << hex << hc << dec << endl;
			return &ZeroF;
	}
}

template complex<R> (*A1ph1QM1sc1q_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph1QM1sc1q_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph1QM1sc1q_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1ph2sc2q_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2sc2q_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2sc2q_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> (*A1ph1QM1sc1q_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2sc2q_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}	
			
			
