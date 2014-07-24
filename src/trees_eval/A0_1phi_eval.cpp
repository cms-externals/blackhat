/* The 2 scalar amplitudes */

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;

#define _HIGGS_GLUON_GLUON_PARAM T(1)
#define _HIGGS_GLUON_GLUON_PARAM_D _HIGGS_GLUON_GLUON_PARAM

namespace BH {

/*
 *
 *
 * A function that returns zero
 *
 *
 */

template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses){return complex<T>(0,0);}


/*
 *
 *
 * The 1 complex higgs and two gluons amplitudes
 *
 *
 */

//ph--
template <int i0, int i1, int i2, class T> complex<T>  A1ph2g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){	
	return complex<T>(0,-1)*_HIGGS_GLUON_GLUON_PARAM*pow(ep.spa(i1,i2),2);
}
//phd++
template <int i0, int i1, int i2, class T> complex<T>  A1ph2g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*_HIGGS_GLUON_GLUON_PARAM_D*pow(ep.spb(i1,i2),2);
}

template <class T> complex<T> (*A1ph2g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 0x900://ph--
		return &A1ph2g1_eval<0,1,2>;
		break;
	case 0x090://-ph-
		return &A1ph2g1_eval<1,2,0>;
		break;
	case 0x009://--ph
		return &A1ph2g1_eval<2,0,1>;
		break;

	case 0xA11://phd++
		return &A1ph2g2_eval<0,1,2>;
		break;
	case 0x1A1://+phd+
		return &A1ph2g2_eval<1,0,2>;
		break;
	case 0x11A://++phd
		return &A1ph2g2_eval<2,0,1>;
		break;
			
			
	case 0xE00://H--
		return &A1ph2g1_eval<0,1,2>;
		break;
	case 0x0E0://-H-
		return &A1ph2g1_eval<1,2,0>;
		break;
	case 0x00E://--H
		return &A1ph2g1_eval<2,0,1>;
		break;			
	case 0xE11://H++
		return &A1ph2g2_eval<0,1,2>;
		break;
	case 0x1E1://+H+
		return &A1ph2g2_eval<1,0,2>;
		break;
	case 0x11E://++H
		return &A1ph2g2_eval<2,0,1>;
		break;

	default:// We return zero for all other helicity combinations
		return &ZeroF;
	}
}



/*
 *
 *
 * The 1 complex higgs and three gluons amplitudes
 *
 *
 */

//ph--+
template <int i1, int i2, int i3, class T> complex<T>  A1ph3g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,1)*pow(ep.spa(i1,i2),3)/(ep.spa(i2,i3)*ep.spa(i3,i1));
}
//ph---
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph3g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*pow(eval_param<T>::mass2(masses.p(i0)),2)/(ep.spb(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i1));
}
//phd++-
template <int i1, int i2, int i3, class T> complex<T>  A1ph3g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*pow(ep.spb(i1,i2),3)/(ep.spb(i2,i3)*ep.spb(i3,i1));
}
//phd+++
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph3g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,1)*pow(eval_param<T>::mass2(masses.p(i0)),2)/(ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i1));
}


template <class T, int HIGGS, int G1, int G2, int G3> complex<T> (*A1ph3g_Tree_Ptr_higgs_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
		case 0x001://ph--+
			return &A1ph3g1_eval<G1,G2,G3>;
			break;
			
		case 0x010://ph-+-
			return &A1ph3g1_eval<G3,G1,G2>;//5
			break;
			
		case 0x100://ph+--
			return &A1ph3g1_eval<G2,G3,G1>;//6
			break;
			
		case 0x000://ph---
			return &A1ph3g2_eval<HIGGS,G1,G2,G3>;
			break;
						
		default:// We return zero for all other helicity combinations
			return &ZeroF;
	}
}

template <class T, int HIGGS, int G1, int G2, int G3> complex<T> (*A1ph3g_Tree_Ptr_higgs_dag_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {			
		case 0x110://phd++-
			return &A1ph3g3_eval<G1,G2,G3>;
			break;
			
		case 0x101://phd+-+
			return &A1ph3g3_eval<G3,G1,G2>;//7
			break;
			
		case 0x011://phd-++
			return &A1ph3g3_eval<G2,G3,G1>;//8
			break;
			
		case 0x111://phd+++
			return &A1ph3g4_eval<HIGGS,G1,G2,G3>;
			break;
			
		default:// We return zero for all other helicity combinations
			return &ZeroF;
	}
}
	
	
template <class T, int HIGGS, int G1, int G2, int G3> complex<T> (*A1ph3g_Tree_Ptr_higgs_full_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
		case 0x001://H--+
			return &A1ph3g1_eval<G1,G2,G3>;
			break;
			
		case 0x010://H-+-
			return &A1ph3g1_eval<G3,G1,G2>;//5
			break;
			
		case 0x100://H+--
			return &A1ph3g1_eval<G2,G3,G1>;//6
			break;
			
		case 0x000://H---
			return &A1ph3g2_eval<HIGGS,G1,G2,G3>;
			break;
			
		case 0x110://H++-
			return &A1ph3g3_eval<G1,G2,G3>;
			break;
			
		case 0x101://H+-+
			return &A1ph3g3_eval<G3,G1,G2>;//7
			break;
			
		case 0x011://H-++
			return &A1ph3g3_eval<G2,G3,G1>;//8
			break;
			
		case 0x111://H+++
			return &A1ph3g4_eval<HIGGS,G1,G2,G3>;
			break;
			
		default:// We return zero for all other helicity combinations
			return &ZeroF;
	}
}

	
	
template <class T> complex<T> (*A1ph3g_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	long rem=hc;
	long newhc=0;
	int epos;
	long fac=0x1;
	bool dagger=false;
	bool full_higgs=false;
	for(int i=0;i<4;i++){
		int efac=rem%0x10;
		rem=rem/0x10;
		if(efac==0x9||efac==0xA||0xE){
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
		if(!dagger&&!full_higgs){
			switch (epos) {
				case 0:
					return A1ph3g_Tree_Ptr_higgs_eval<T,3,0,1,2>(newhc);
					break;
				case 1:
					return A1ph3g_Tree_Ptr_higgs_eval<T,2,0,1,3>(newhc);
					break;
				case 2:
					return A1ph3g_Tree_Ptr_higgs_eval<T,1,0,2,3>(newhc);
					break;
				case 3:
					return A1ph3g_Tree_Ptr_higgs_eval<T,0,1,2,3>(newhc);
					break;
					
				default:// Should never get here as there is always a higgs. 
					return &ZeroF;
					break;
			}
		}
		else{
			switch (epos) {
				case 0:
					return A1ph3g_Tree_Ptr_higgs_dag_eval<T,3,0,1,2>(newhc);
					break;
				case 1:
					return A1ph3g_Tree_Ptr_higgs_dag_eval<T,2,0,1,3>(newhc);
					break;
				case 2:
					return A1ph3g_Tree_Ptr_higgs_dag_eval<T,1,0,2,3>(newhc);
					break;
				case 3:
					return A1ph3g_Tree_Ptr_higgs_dag_eval<T,0,1,2,3>(newhc);
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
				return A1ph3g_Tree_Ptr_higgs_full_eval<T,3,0,1,2>(newhc);
				break;
			case 1:
				return A1ph3g_Tree_Ptr_higgs_full_eval<T,2,0,1,3>(newhc);
				break;
			case 2:
				return A1ph3g_Tree_Ptr_higgs_full_eval<T,1,0,2,3>(newhc);
				break;
			case 3:
				return A1ph3g_Tree_Ptr_higgs_full_eval<T,0,1,2,3>(newhc);
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
 * The 1 complex higgs, one gluon and 2 quarks amplitudes
 *
 *
 */

//ph,q-,q+,-
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph1g2q1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spa(i1,i3)*ep.spa(i1,i3)/ep.spa(i1,i2);
}
//ph,q-,-,q+
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph1g2q3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spa(i1,i2)*ep.spa(i1,i2)/ep.spa(i3,i1);
}
//ph,q+,q-,-
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph1g2q4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spa(i2,i3)*ep.spa(i2,i3)/ep.spa(i1,i2);
}
//ph,q+,-,q-
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph1g2q5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spa(i3,i2)*ep.spa(i3,i2)/ep.spa(i3,i1);
}

//phd,q-,q+,+
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph1g2q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spb(i2,i3)*ep.spb(i2,i3)/ep.spb(i1,i2);
}
//phd,q-,+,q+
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph1g2q6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spb(i2,i3)*ep.spb(i2,i3)/ep.spb(i3,i1);
}
//phd,q+,q-,+
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph1g2q7_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spb(i1,i3)*ep.spb(i1,i3)/ep.spb(i1,i2);
}
//phd,q+,+,q-
template <int i0, int i1, int i2, int i3, class T> complex<T>  A1ph1g2q8_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spb(i2,i1)*ep.spb(i2,i1)/ep.spb(i3,i1);
}

template <class T> complex<T> (*A1ph1g2q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 0x9230://phq-q+-
		return &A1ph1g2q1_eval<0,1,2,3>;
		break;
	case 0x0923://-phq-q+
		return &A1ph1g2q1_eval<1,2,3,0>;
		break;
	case 0x3092://q+-phq-
		return &A1ph1g2q1_eval<2,3,0,1>;
		break;
	case 0x2309://q-q+-ph
		return &A1ph1g2q1_eval<3,0,1,2>;
		break;
	case 0x9320://phq+q--
		return &A1ph1g2q4_eval<0,1,2,3>;
		break;
	case 0x0932://-phq+q-
		return &A1ph1g2q4_eval<1,2,3,0>;
		break;
	case 0x2093://q--phq+
		return &A1ph1g2q4_eval<2,3,0,2>;
		break;
	case 0x3209://q+q--ph
		return &A1ph1g2q4_eval<3,0,1,2>;
		break;

	case 0x9203://phq--q+
		return &A1ph1g2q3_eval<0,1,2,3>;
		break;
	case 0x3920://q+phq--
		return &A1ph1g2q3_eval<1,2,3,0>;
		break;
	case 0x0392://-q+phq-
		return &A1ph1g2q3_eval<2,3,0,1>;
		break;
	case 0x2039://q--q+ph
		return &A1ph1g2q3_eval<3,0,1,2>;
		break;
	case 0x9302://phq+-q-
		return &A1ph1g2q5_eval<0,1,2,3>;
		break;
	case 0x2930://q-phq+-
		return &A1ph1g2q5_eval<1,2,3,0>;
		break;
	case 0x0293://-q-phq+
		return &A1ph1g2q5_eval<2,3,0,1>;
		break;
	case 0x3029://q+-q-ph
		return &A1ph1g2q5_eval<3,0,1,2>;
		break;

	case 0x9023://ph-q-q+
		return &A1ph1g2q1_eval<0,2,3,1>;
		break;
	case 0x3902://q+ph-q-
		return &A1ph1g2q1_eval<1,3,0,2>;
		break;
	case 0x2390://q-q+ph-
		return &A1ph1g2q1_eval<2,0,1,3>;
		break;
	case 0x0239://-q-q+ph
		return &A1ph1g2q1_eval<3,1,2,0>;
		break;
	case 0x9032://ph-q+q-
		return &A1ph1g2q4_eval<0,2,3,1>;
		break;
	case 0x2903://q-ph-q+
		return &A1ph1g2q4_eval<1,3,0,2>;
		break;
	case 0x3290://q+q-ph-
		return &A1ph1g2q4_eval<2,0,1,3>;
		break;
	case 0x0329://-q+q-ph
		return &A1ph1g2q4_eval<3,1,2,0>;
		break;

	case 0xA231://phdq-q++
		return &A1ph1g2q2_eval<0,1,2,3>;
		break;
	case 0x1A23://+phdq-q+
		return &A1ph1g2q2_eval<1,2,3,0>;
		break;
	case 0x31A2://q++phdq-
		return &A1ph1g2q2_eval<2,3,0,1>;
		break;
	case 0x231A://q-q++phd
		return &A1ph1g2q2_eval<3,0,1,2>;
		break;
	case 0xA321://phdq+q-+
		return &A1ph1g2q7_eval<0,1,2,3>;
		break;
	case 0x1A32://+phdq+q-
		return &A1ph1g2q7_eval<1,2,3,0>;
		break;
	case 0x21A3://q-+phdq+
		return &A1ph1g2q7_eval<2,3,0,1>;
		break;
	case 0x321A://q+q-+phd
		return &A1ph1g2q7_eval<0,1,2,3>;
		break;

	case 0xA213://phq--q+
		return &A1ph1g2q6_eval<0,1,2,3>;
		break;
	case 0x3A21://q+phq--
		return &A1ph1g2q6_eval<1,2,3,0>;
		break;
	case 0x13A2://-q+phq-
		return &A1ph1g2q6_eval<2,3,0,1>;
		break;
	case 0x213A://q--q+ph
		return &A1ph1g2q6_eval<3,0,1,2>;
		break;
	case 0xA312://phq+-q-
		return &A1ph1g2q8_eval<0,1,2,3>;
		break;
	case 0x2A31://q-phq+-
		return &A1ph1g2q8_eval<1,2,3,0>;
		break;
	case 0x12A3://-q-phq+
		return &A1ph1g2q8_eval<2,3,0,1>;
		break;
	case 0x312A://q+-q-ph
		return &A1ph1g2q8_eval<3,0,1,2>;
		break;

	case 0xA123://ph-q-q+
		return &A1ph1g2q2_eval<0,2,3,1>;
		break;
	case 0x3A12://q+ph-q-
		return &A1ph1g2q2_eval<1,3,0,2>;
		break;
	case 0x23A1://q-q+ph-
		return &A1ph1g2q2_eval<2,0,1,3>;
		break;
	case 0x123A://-q-q+ph
		return &A1ph1g2q2_eval<3,1,2,0>;
		break;
	case 0xA132://ph-q+q-
		return &A1ph1g2q7_eval<0,2,3,1>;
		break;
	case 0x2A13://q-ph-q+
		return &A1ph1g2q7_eval<1,3,0,2>;
		break;
	case 0x32A1://q+q-ph-
		return &A1ph1g2q7_eval<2,0,1,3>;
		break;
	case 0x132A://-q+q-ph
		return &A1ph1g2q7_eval<3,1,2,0>;
		break;

	default:// We return zero for all other helicity combinations
//		return &ZeroF;
		return 0;
	}
}




/*
 *
 *
 * The 1 complex higgs, two gluons and 2 quarks amplitudes
 *
 *
 */

////ph--
//template <int i0, int i1, int i2, class T> complex<T>  A1ph2g2q1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	return complex<T>(0,-1)*ep.spa(i1,i2)*ep.spa(i1,i2);
//}
////phd++
//template <int i0, int i1, int i2, class T> complex<T>  A1ph2g2q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
//	return complex<T>(0,1)*ep.spb(i1,i2)*ep.spb(i1,i2);
//}

template <class T> complex<T> (*A1ph2g2q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {

	default:// We return zero for all other helicity combinations
//		return &ZeroF;
		return 0; // We don't do anything so return telling the calling function that we did nothing
	}
}



/*
 *
 *
 * The 1 complex higgs and 4 quarks amplitudes
 *
 *
 */

//ph,qm,qp,q2p,q2m
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*pow(ep.spa(i1,i4),2)/(ep.spa(i1,i2)*ep.spa(i3,i4));
}
//ph,qm,qp,q2m,q2p
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*pow(ep.spa(i1,i3),2)/(ep.spa(i1,i2)*ep.spa(i3,i4));
}
//ph,qp,qm,q2m,q2p
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*pow(ep.spa(i2,i3),2)/(ep.spa(i1,i2)*ep.spa(i3,i4));
}
//ph,qp,qm,q2p,q2m
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*pow(ep.spa(i2,i4),2)/(ep.spa(i1,i2)*ep.spa(i3,i4));
}

//phd,qm,qp,q2p,q2m
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spb(i2,i3)*ep.spb(i2,i3)/(ep.spb(i1,i2)*ep.spb(i3,i4));
}
//phd,qm,qp,q2m,q2p
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q4_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spb(i2,i4)*ep.spb(i2,i4)/(ep.spb(i1,i2)*ep.spb(i3,i4));
}
//phd,qp,qm,q2m,q2p
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q7_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spb(i1,i4)*ep.spb(i1,i4)/(ep.spb(i1,i2)*ep.spb(i3,i4));
}
//phd,qp,qm,q2p,q2m
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q8_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spb(i1,i3)*ep.spb(i1,i3)/(ep.spb(i1,i2)*ep.spb(i3,i4));
}

	
//H,qm,qp,q2p,q2m
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q9_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
	return complex<T>(0,-1)*(pow(ep.spa(i1,i4),2)/(ep.spa(i1,i2)*ep.spa(i3,i4))+pow(ep.spb(i2,i3),2)/(ep.spb(i1,i2)*ep.spb(i3,i4)));
}
//H,qm,qp,q2m,q2p
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q10_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
	return complex<T>(0,-1)*(pow(ep.spa(i1,i3),2)/(ep.spa(i1,i2)*ep.spa(i3,i4))+pow(ep.spb(i2,i4),2)/(ep.spb(i1,i2)*ep.spb(i3,i4)));
}
//H,qp,qm,q2m,q2p
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q11_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
	return complex<T>(0,-1)*(pow(ep.spa(i2,i3),2)/(ep.spa(i1,i2)*ep.spa(i3,i4))+pow(ep.spb(i1,i4),2)/(ep.spb(i1,i2)*ep.spb(i3,i4)));
}
//H,qp,qm,q2p,q2m
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A1ph2q2Q12_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	
	return complex<T>(0,-1)*(pow(ep.spa(i2,i4),2)/(ep.spa(i1,i2)*ep.spa(i3,i4))+pow(ep.spb(i1,i3),2)/(ep.spb(i1,i2)*ep.spb(i3,i4)));
}


template <class T> complex<T> (*A1ph2q2Q_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {
//	cout << hex << hc << dec << endl;
	switch (hc) {
	case 0x0E02030908://ph,qm,qp,q2p,q2m
		return &A1ph2q2Q1_eval<0,1,2,3,4>;
		break;
	case 0x080E020309:
		return &A1ph2q2Q1_eval<1,2,3,4,0>;
		break;
	case 0x09080E0203:
		return &A1ph2q2Q1_eval<2,3,4,0,1>;
		break;
	case 0x0309080E02:
		return &A1ph2q2Q1_eval<3,4,0,1,2>;
		break;
	case 0x020309080E:
		return &A1ph2q2Q1_eval<4,0,1,2,3>;
		break;
	case 0x0E08090302://ph,qm,qp,q2p,q2m
		return &A1ph2q2Q1_eval<0,1,2,3,4>;
		break;
	case 0x020E080903:
		return &A1ph2q2Q1_eval<1,2,3,4,0>;
		break;
	case 0x03020E0809:
		return &A1ph2q2Q1_eval<2,3,4,0,1>;
		break;
	case 0x0903020E08:
		return &A1ph2q2Q1_eval<3,4,0,1,2>;
		break;
	case 0x080903020E:
		return &A1ph2q2Q1_eval<4,0,1,2,3>;
		break;


	case 0x0E02030809://ph,qm,qp,q2m,q2p
		return &A1ph2q2Q3_eval<0,1,2,3,4>;
		break;
	case 0x090E020308:
		return &A1ph2q2Q3_eval<1,2,3,4,0>;
		break;
	case 0x08090E0203:
		return &A1ph2q2Q3_eval<2,3,4,0,1>;
		break;
	case 0x0308090E02:
		return &A1ph2q2Q3_eval<3,4,0,1,2>;
		break;
	case 0x020308090E:
		return &A1ph2q2Q3_eval<4,0,1,2,3>;
		break;
	case 0x0E08090203://ph,qm,qp,q2m,q2p
		return &A1ph2q2Q3_eval<0,1,2,3,4>;
		break;
	case 0x030E080902:
		return &A1ph2q2Q3_eval<1,2,3,4,0>;
		break;
	case 0x02030E0809:
		return &A1ph2q2Q3_eval<2,3,4,0,1>;
		break;
	case 0x0902030E08:
		return &A1ph2q2Q3_eval<3,4,0,1,2>;
		break;
	case 0x080902030E:
		return &A1ph2q2Q3_eval<4,0,1,2,3>;
		break;


	case 0x0E03020809://ph,qp,qm,q2m,q2p
		return &A1ph2q2Q5_eval<0,1,2,3,4>;
		break;
	case 0x090E030208:
		return &A1ph2q2Q5_eval<1,2,3,4,0>;
		break;
	case 0x08090E0302:
		return &A1ph2q2Q5_eval<2,3,4,0,1>;
		break;
	case 0x0208090E03:
		return &A1ph2q2Q5_eval<3,4,0,1,2>;
		break;
	case 0x030208090E:
		return &A1ph2q2Q5_eval<4,0,1,2,3>;
		break;
	case 0x0E09080203://ph,qp,qm,q2m,q2p
		return &A1ph2q2Q5_eval<0,1,2,3,4>;
		break;
	case 0x030E090802:
		return &A1ph2q2Q5_eval<1,2,3,4,0>;
		break;
	case 0x02030E0908:
		return &A1ph2q2Q5_eval<2,3,4,0,1>;
		break;
	case 0x0802030E09:
		return &A1ph2q2Q5_eval<3,4,0,1,2>;
		break;
	case 0x090802030E:
		return &A1ph2q2Q5_eval<4,0,1,2,3>;
		break;


	case 0x0E03020908://ph,qp,qm,q2p,q2m
		return &A1ph2q2Q6_eval<0,1,2,3,4>;
		break;
	case 0x080E030209:
		return &A1ph2q2Q6_eval<1,2,3,4,0>;
		break;
	case 0x09080E0302:
		return &A1ph2q2Q6_eval<2,3,4,0,1>;
		break;
	case 0x0209080E03:
		return &A1ph2q2Q6_eval<3,4,0,1,2>;
		break;
	case 0x030209080E:
		return &A1ph2q2Q6_eval<4,0,1,2,3>;
		break;
	case 0x0E09080302://ph,qp,qm,q2p,q2m
		return &A1ph2q2Q6_eval<0,1,2,3,4>;
		break;
	case 0x020E090803:
		return &A1ph2q2Q6_eval<1,2,3,4,0>;
		break;
	case 0x03020E0908:
		return &A1ph2q2Q6_eval<2,3,4,0,1>;
		break;
	case 0x0803020E09:
		return &A1ph2q2Q6_eval<3,4,0,1,2>;
		break;
	case 0x090803020E:
		return &A1ph2q2Q6_eval<4,0,1,2,3>;
		break;


	case 0x0F02030908://phd,qm,qp,q2p,q2m
		return &A1ph2q2Q2_eval<0,1,2,3,4>;
		break;
	case 0x080F020309:
		return &A1ph2q2Q2_eval<1,2,3,4,0>;
		break;
	case 0x09080F0203:
		return &A1ph2q2Q2_eval<2,3,4,0,1>;
		break;
	case 0x0309080F02:
		return &A1ph2q2Q2_eval<3,4,0,1,2>;
		break;
	case 0x020309080F:
		return &A1ph2q2Q2_eval<4,0,1,2,3>;
		break;
	case 0x0F08090302://phd,qm,qp,q2p,q2m
		return &A1ph2q2Q2_eval<0,1,2,3,4>;
		break;
	case 0x020F080903:
		return &A1ph2q2Q2_eval<1,2,3,4,0>;
		break;
	case 0x03020F0809:
		return &A1ph2q2Q2_eval<2,3,4,0,1>;
		break;
	case 0x0903020F08:
		return &A1ph2q2Q2_eval<3,4,0,1,2>;
		break;
	case 0x080903020F:
		return &A1ph2q2Q2_eval<4,0,1,2,3>;
		break;


	case 0x0F02030809://phd,qm,qp,q2m,q2p
		return &A1ph2q2Q4_eval<0,1,2,3,4>;
		break;
	case 0x090F020308:
		return &A1ph2q2Q4_eval<1,2,3,4,0>;
		break;
	case 0x08090F0203:
		return &A1ph2q2Q4_eval<2,3,4,0,1>;
		break;
	case 0x0308090F02:
		return &A1ph2q2Q4_eval<3,4,0,1,2>;
		break;
	case 0x020308090F:
		return &A1ph2q2Q4_eval<4,0,1,2,3>;
		break;
	case 0x0F08090203://phd,qm,qp,q2m,q2p
		return &A1ph2q2Q4_eval<0,1,2,3,4>;
		break;
	case 0x030F080902:
		return &A1ph2q2Q4_eval<1,2,3,4,0>;
		break;
	case 0x02030F0809:
		return &A1ph2q2Q4_eval<2,3,4,0,1>;
		break;
	case 0x0902030F08:
		return &A1ph2q2Q4_eval<3,4,0,1,2>;
		break;
	case 0x080902030F:
		return &A1ph2q2Q4_eval<4,0,1,2,3>;
		break;


	case 0x0F03020809://phd,qp,qm,q2m,q2p
		return &A1ph2q2Q7_eval<0,1,2,3,4>;
		break;
	case 0x090F030208:
		return &A1ph2q2Q7_eval<1,2,3,4,0>;
		break;
	case 0x08090F0302:
		return &A1ph2q2Q7_eval<2,3,4,0,1>;
		break;
	case 0x0208090F03:
		return &A1ph2q2Q7_eval<3,4,0,1,2>;
		break;
	case 0x030208090F:
		return &A1ph2q2Q7_eval<4,0,1,2,3>;
		break;
	case 0x0F09080203://phd,qp,qm,q2m,q2p
		return &A1ph2q2Q7_eval<0,1,2,3,4>;
		break;
	case 0x030F090802:
		return &A1ph2q2Q7_eval<1,2,3,4,0>;
		break;
	case 0x02030F0908:
		return &A1ph2q2Q7_eval<2,3,4,0,1>;
		break;
	case 0x0802030F09:
		return &A1ph2q2Q7_eval<3,4,0,1,2>;
		break;
	case 0x090802030F:
		return &A1ph2q2Q7_eval<4,0,1,2,3>;
		break;


	case 0x0F03020908://phd,qp,qm,q2p,q2m
		return &A1ph2q2Q8_eval<0,1,2,3,4>;
		break;
	case 0x080F030209:
		return &A1ph2q2Q8_eval<1,2,3,4,0>;
		break;
	case 0x09080F0302:
		return &A1ph2q2Q8_eval<2,3,4,0,1>;
		break;
	case 0x0209080F03:
		return &A1ph2q2Q8_eval<3,4,0,1,2>;
		break;
	case 0x030209080F:
		return &A1ph2q2Q8_eval<4,0,1,2,3>;
		break;
	case 0x0F09080302://phd,qp,qm,q2p,q2m
		return &A1ph2q2Q8_eval<0,1,2,3,4>;
		break;
	case 0x020F090803:
		return &A1ph2q2Q8_eval<1,2,3,4,0>;
		break;
	case 0x03020F0908:
		return &A1ph2q2Q8_eval<2,3,4,0,1>;
		break;
	case 0x0803020F09:
		return &A1ph2q2Q8_eval<3,4,0,1,2>;
		break;
	case 0x090803020F:
		return &A1ph2q2Q8_eval<4,0,1,2,3>;
		break;
			
			
			
		//H
			
	case 0x1002030908://ph,qm,qp,q2p,q2m
		return &A1ph2q2Q9_eval<0,1,2,3,4>;
		break;
	case 0x0810020309:
		return &A1ph2q2Q9_eval<1,2,3,4,0>;
		break;
	case 0x0908100203:
		return &A1ph2q2Q9_eval<2,3,4,0,1>;
		break;
	case 0x0309081002:
		return &A1ph2q2Q9_eval<3,4,0,1,2>;
		break;
	case 0x0203090810:
		return &A1ph2q2Q9_eval<4,0,1,2,3>;
		break;
	case 0x1008090302://ph,qm,qp,q2p,q2m
		return &A1ph2q2Q9_eval<0,1,2,3,4>;
		break;
	case 0x0210080903:
		return &A1ph2q2Q9_eval<1,2,3,4,0>;
		break;
	case 0x0302100809:
		return &A1ph2q2Q9_eval<2,3,4,0,1>;
		break;
	case 0x0903021008:
		return &A1ph2q2Q9_eval<3,4,0,1,2>;
		break;
	case 0x0809030210:
		return &A1ph2q2Q9_eval<4,0,1,2,3>;
		break;
		
		
	case 0x1002030809://ph,qm,qp,q2m,q2p
		return &A1ph2q2Q10_eval<0,1,2,3,4>;
		break;
	case 0x0910020308:
		return &A1ph2q2Q10_eval<1,2,3,4,0>;
		break;
	case 0x0809100203:
		return &A1ph2q2Q10_eval<2,3,4,0,1>;
		break;
	case 0x0308091002:
		return &A1ph2q2Q10_eval<3,4,0,1,2>;
		break;
	case 0x0203080910:
		return &A1ph2q2Q10_eval<4,0,1,2,3>;
		break;
	case 0x1008090203://ph,qm,qp,q2m,q2p
		return &A1ph2q2Q10_eval<0,1,2,3,4>;
		break;
	case 0x0310080902:
		return &A1ph2q2Q10_eval<1,2,3,4,0>;
		break;
	case 0x0203100809:
		return &A1ph2q2Q10_eval<2,3,4,0,1>;
		break;
	case 0x0902031008:
		return &A1ph2q2Q10_eval<3,4,0,1,2>;
		break;
	case 0x0809020310:
		return &A1ph2q2Q10_eval<4,0,1,2,3>;
		break;
		
		
	case 0x1003020809://ph,qp,qm,q2m,q2p
		return &A1ph2q2Q11_eval<0,1,2,3,4>;
		break;
	case 0x0910030208:
		return &A1ph2q2Q11_eval<1,2,3,4,0>;
		break;
	case 0x0809100302:
		return &A1ph2q2Q11_eval<2,3,4,0,1>;
		break;
	case 0x0208091003:
		return &A1ph2q2Q11_eval<3,4,0,1,2>;
		break;
	case 0x0302080910:
		return &A1ph2q2Q11_eval<4,0,1,2,3>;
		break;
	case 0x1009080203://ph,qp,qm,q2m,q2p
		return &A1ph2q2Q11_eval<0,1,2,3,4>;
		break;
	case 0x0310090802:
		return &A1ph2q2Q11_eval<1,2,3,4,0>;
		break;
	case 0x0203100908:
		return &A1ph2q2Q11_eval<2,3,4,0,1>;
		break;
	case 0x0802031009:
		return &A1ph2q2Q11_eval<3,4,0,1,2>;
		break;
	case 0x0908020310:
		return &A1ph2q2Q11_eval<4,0,1,2,3>;
		break;
		
		
	case 0x1003020908://ph,qp,qm,q2p,q2m
		return &A1ph2q2Q12_eval<0,1,2,3,4>;
		break;
	case 0x0810030209:
		return &A1ph2q2Q12_eval<1,2,3,4,0>;
		break;
	case 0x0908100302:
		return &A1ph2q2Q12_eval<2,3,4,0,1>;
		break;
	case 0x0209081003:
		return &A1ph2q2Q12_eval<3,4,0,1,2>;
		break;
	case 0x0302090810:
		return &A1ph2q2Q12_eval<4,0,1,2,3>;
		break;
	case 0x1009080302://ph,qp,qm,q2p,q2m
		return &A1ph2q2Q12_eval<0,1,2,3,4>;
		break;
	case 0x0210090803:
		return &A1ph2q2Q12_eval<1,2,3,4,0>;
		break;
	case 0x0302100908:
		return &A1ph2q2Q12_eval<2,3,4,0,1>;
		break;
	case 0x0803021009:
		return &A1ph2q2Q12_eval<3,4,0,1,2>;
		break;
	case 0x0908030210:
		return &A1ph2q2Q12_eval<4,0,1,2,3>;
		break;
		
		
//	case 0x1002030908://phd,qm,qp,q2p,q2m
//		return &A1ph2q2Q2_eval<0,1,2,3,4>;
//		break;
//	case 0x0810020309:
//		return &A1ph2q2Q2_eval<1,2,3,4,0>;
//		break;
//	case 0x0908100203:
//		return &A1ph2q2Q2_eval<2,3,4,0,1>;
//		break;
//	case 0x0309081002:
//		return &A1ph2q2Q2_eval<3,4,0,1,2>;
//		break;
//	case 0x0203090810:
//		return &A1ph2q2Q2_eval<4,0,1,2,3>;
//		break;
//	case 0x1008090302://phd,qm,qp,q2p,q2m
//		return &A1ph2q2Q2_eval<0,1,2,3,4>;
//		break;
//	case 0x0210080903:
//		return &A1ph2q2Q2_eval<1,2,3,4,0>;
//		break;
//	case 0x0302100809:
//		return &A1ph2q2Q2_eval<2,3,4,0,1>;
//		break;
//	case 0x0903021008:
//		return &A1ph2q2Q2_eval<3,4,0,1,2>;
//		break;
//	case 0x0809030210:
//		return &A1ph2q2Q2_eval<4,0,1,2,3>;
//		break;
//		
//		
//	case 0x1002030809://phd,qm,qp,q2m,q2p
//		return &A1ph2q2Q4_eval<0,1,2,3,4>;
//		break;
//	case 0x0910020308:
//		return &A1ph2q2Q4_eval<1,2,3,4,0>;
//		break;
//	case 0x0809100203:
//		return &A1ph2q2Q4_eval<2,3,4,0,1>;
//		break;
//	case 0x0308091002:
//		return &A1ph2q2Q4_eval<3,4,0,1,2>;
//		break;
//	case 0x0203080910:
//		return &A1ph2q2Q4_eval<4,0,1,2,3>;
//		break;
//	case 0x1008090203://phd,qm,qp,q2m,q2p
//		return &A1ph2q2Q4_eval<0,1,2,3,4>;
//		break;
//	case 0x0310080902:
//		return &A1ph2q2Q4_eval<1,2,3,4,0>;
//		break;
//	case 0x0203100809:
//		return &A1ph2q2Q4_eval<2,3,4,0,1>;
//		break;
//	case 0x0902031008:
//		return &A1ph2q2Q4_eval<3,4,0,1,2>;
//		break;
//	case 0x0809020310:
//		return &A1ph2q2Q4_eval<4,0,1,2,3>;
//		break;
//		
//		
//	case 0x1003020809://phd,qp,qm,q2m,q2p
//		return &A1ph2q2Q7_eval<0,1,2,3,4>;
//		break;
//	case 0x0910030208:
//		return &A1ph2q2Q7_eval<1,2,3,4,0>;
//		break;
//	case 0x0809100302:
//		return &A1ph2q2Q7_eval<2,3,4,0,1>;
//		break;
//	case 0x0208091003:
//		return &A1ph2q2Q7_eval<3,4,0,1,2>;
//		break;
//	case 0x0302080910:
//		return &A1ph2q2Q7_eval<4,0,1,2,3>;
//		break;
//	case 0x1009080203://phd,qp,qm,q2m,q2p
//		return &A1ph2q2Q7_eval<0,1,2,3,4>;
//		break;
//	case 0x0310090802:
//		return &A1ph2q2Q7_eval<1,2,3,4,0>;
//		break;
//	case 0x0203100908:
//		return &A1ph2q2Q7_eval<2,3,4,0,1>;
//		break;
//	case 0x0802031009:
//		return &A1ph2q2Q7_eval<3,4,0,1,2>;
//		break;
//	case 0x0908020310:
//		return &A1ph2q2Q7_eval<4,0,1,2,3>;
//		break;
//		
//		
//	case 0x1003020908://phd,qp,qm,q2p,q2m
//		return &A1ph2q2Q8_eval<0,1,2,3,4>;
//		break;
//	case 0x0810030209:
//		return &A1ph2q2Q8_eval<1,2,3,4,0>;
//		break;
//	case 0x0908100302:
//		return &A1ph2q2Q8_eval<2,3,4,0,1>;
//		break;
//	case 0x0209081003:
//		return &A1ph2q2Q8_eval<3,4,0,1,2>;
//		break;
//	case 0x0302090810:
//		return &A1ph2q2Q8_eval<4,0,1,2,3>;
//		break;
//	case 0x1009080302://phd,qp,qm,q2p,q2m
//		return &A1ph2q2Q8_eval<0,1,2,3,4>;
//		break;
//	case 0x0210090803:
//		return &A1ph2q2Q8_eval<1,2,3,4,0>;
//		break;
//	case 0x0302100908:
//		return &A1ph2q2Q8_eval<2,3,4,0,1>;
//		break;
//	case 0x0803021009:
//		return &A1ph2q2Q8_eval<3,4,0,1,2>;
//		break;
//	case 0x0908030210:
//		return &A1ph2q2Q8_eval<4,0,1,2,3>;
//		break;

	default:// We return zero for all other helicity combinations
//		return &ZeroF;
		return 0;
	}
}
	


template complex<R> (*A1ph2g_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2g_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2g_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);



template complex<R> (*A1ph3g_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph3g_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph3g_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);



template complex<R> (*A1ph1g2q_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph1g2q_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph1g2q_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);



template complex<R> (*A1ph2g2q_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2g2q_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2g2q_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A1ph2q2Q_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A1ph2q2Q_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A1ph2q2Q_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> (*A1ph2g_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph3g_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph1g2q_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2g2q_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
template complex<RGMP> (*A1ph2q2Q_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);

#endif
		
}
