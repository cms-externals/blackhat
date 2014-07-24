/*
 * A0_2q_massive_lepton_eval.cpp
 *
 *  Created on: Aug_eval 28, 2008
 *      Author: Darren
 *
 */
#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"
#include <iostream>


using namespace std;

namespace BH {

template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses){return complex<T>(0,0);}

/*
 *
 *
 *
 * Two scalar two lepton amplitudes
 *
 *
 *
 */

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2s2l1_eval(const eval_param<T>& ep, const mass_param_coll& masses){

	// The first term here corresponds to the S channel gluon exchange bewteen a fermion and a scalar line
	//  the second line is a T channel exchange of a massive quark between scalar massless quark vertices
	//	return complex<T>(0.,1.)*ep.spab(i1,i0,i2)/(ep.sp(i1,i2))
	//					+complex<T>(0.,0.5)*ep.spab(i1,i0,i2)/ep.sp(i0,i1);

	return complex<T>(0,1)*ep.spab(i1,i0,i2)*(T(1.)/ep.sp(i1,i2)/*+T(1)/(T(2)*ep.sp(i0,i1))*/);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2s2l2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return complex<T>(0,-1)*ep.spab(i2,i0,i1)*(T(1.)/ep.sp(i1,i2)/*+T(1)/(T(2)*ep.sp(i0,i1))*/);
}

template <class T> complex<T> (*A2s2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	// case 0x: //2qs_massive
	switch (hc) {
	case 38324://0x040F1004: //2169: //(s,lb-,l+,s)
		return &A2s2l1_eval<0,1,2,3>;
		break;
	case 33916://0x04040F10: //783: //(s,s,lb-,l+)
		return &A2s2l1_eval<1,2,3,0>;
		break;
	case 129695://0x1004040F: //681: //(l+,s,s,lb-)
		return &A2s2l1_eval<2,3,0,1>;
		break;
	case 126484://0x0F100404: //2367: //(lb-,l+,s,s)
		return &A2s2l1_eval<3,0,1,2>;
		break;

	case 38704://0x04100F04: //2127: //(s,l+,lb-,s)
		return &A2s2l2_eval<0,1,2,3>;
		break;
	case 33935://0x0404100F: //489: //(s,s,l+,lb-)
		return &A2s2l2_eval<1,2,3,0>;
		break;
	case 121696://0x0F040410: //1023: //(lb-,s,s,l+)
		return &A2s2l2_eval<2,3,0,1>;
		break;
	case 134084://0x100F0404: //2361: //(l+,lb-,s,s)
		return &A2s2l2_eval<3,0,1,2>;
		break;

	default:// We return zero for all other helicity combinations
		cout << "4 pt A2s2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return &ZeroF;
	}
}


/*
 *
 *
 *
 * Two scalar two lepton amplitudes
 *
 *
 *
 */

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_1m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambdat<T> i1Lat=ep.p(i1)->Lt();
	lambda<T> i2La=ep.p(i2)->L();
	lambdat<T> i2Lat=ep.p(i2)->Lt();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*pow(qLa*i2La,2)*(i1Lat*i2Lat)/(complex<T>(0,-2)*(qLa*p1flatc.L())*(qLa*pnflatc.L())*ep.sp(i1,i2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_1p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambdat<T> i1Lat=ep.p(i1)->Lt();
	lambda<T> i2La=ep.p(i2)->L();
	lambdat<T> i2Lat=ep.p(i2)->Lt();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*pow(qLa*i2La,2)*(i1Lat*i2Lat)/(complex<T>(0,2)*(qLa*p1flatc.L())*(qLa*pnflatc.L())*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambdat<T> i1Lat=ep.p(i1)->Lt();
	lambda<T> i2La=ep.p(i2)->L();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return ((p1flatc.Lt()*i1Lat)*(i2La*pnflatc.L())
									-mass*(i2La*qLa)*(qLat*i1Lat)/((qLa*p1flatc.L())*(pnflatc.Lt()*qLat)))/(complex<T>(0,2)*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_3_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambdat<T> i1Lat=ep.p(i1)->Lt();
	lambda<T> i2La=ep.p(i2)->L();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return ((pnflatc.Lt()*i1Lat)*(i2La*p1flatc.L())-mass*(i2La*qLa)*(qLat*i1Lat)/((qLa*pnflatc.L())*(p1flatc.Lt()*qLat)))/(complex<T>(0,-2)*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambda<T> i1La=ep.p(i1)->L();
	lambda<T> i2La=ep.p(i2)->L();
	lambdat<T> i2Lat=ep.p(i2)->Lt();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*pow(qLat*i2Lat,2)*(i1La*i2La)/(complex<T>(0,-2)*(qLat*p1flatc.Lt())*(qLat*pnflatc.Lt())*ep.sp(i1,i2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));



	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambda<T> i1La=ep.p(i1)->L();
	lambda<T> i2La=ep.p(i2)->L();
	lambdat<T> i2Lat=ep.p(i2)->Lt();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*pow(qLat*i2Lat,2)*(i1La*i2La)/(complex<T>(0,2)*(qLat*p1flatc.Lt())*(qLat*pnflatc.Lt())*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambda<T> i1La=ep.p(i1)->L();
	lambdat<T> i2Lat=ep.p(i2)->Lt();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return ((p1flatc.L()*i1La)*(i2Lat*pnflatc.Lt())
									-mass*(i2Lat*qLat)*(qLa*i1La)/((qLat*p1flatc.Lt())*(pnflatc.L()*qLa)))/(complex<T>(0,2)*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	//Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambda<T> i1La=ep.p(i1)->L();
	lambdat<T> i2Lat=ep.p(i2)->Lt();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return ((pnflatc.L()*i1La)*(i2Lat*p1flatc.Lt())-mass*(i2Lat*qLat)*(qLa*i1La)/((qLat*pnflatc.Lt())*(p1flatc.L()*qLa)))/(complex<T>(0,-2)*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_7m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambdat<T> i1Lat=ep.p(i1)->Lt();
	lambda<T> i1La=ep.p(i1)->L();
	lambdat<T> i2Lat=ep.p(i2)->Lt();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*pow(qLa*i1La,2)*(i1Lat*i2Lat)/(complex<T>(0,-2)*(qLa*p1flatc.L())*(qLa*pnflatc.L())*ep.sp(i1,i2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_7p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambdat<T> i1Lat=ep.p(i1)->Lt();
	lambda<T> i1La=ep.p(i1)->L();
	lambdat<T> i2Lat=ep.p(i2)->Lt();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*pow(qLa*i1La,2)*(i1Lat*i2Lat)/(complex<T>(0,2)*(qLa*p1flatc.L())*(qLa*pnflatc.L())*ep.sp(i1,i2));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambda<T> i1La=ep.p(i1)->L();
	lambdat<T> i1Lat=ep.p(i1)->Lt();
	lambda<T> i2La=ep.p(i2)->L();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*pow(qLat*i1Lat,2)*(i1La*i2La)/(complex<T>(0,-2)*(qLat*p1flatc.Lt())*(qLat*pnflatc.Lt())*ep.sp(i1,i2));
}
template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM2l_8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	Cmom<T> p1flatc(ep.mom(i0)-T(0.5)*(mass/(ep.p(i0)->P()*ep.ref()->P()))*ep.ref()->P());
	Cmom<T> pnflatc(ep.mom(i3)-T(0.5)*(mass/(ep.p(i3)->P()*ep.ref()->P()))*ep.ref()->P());
	lambda<T> i1La=ep.p(i1)->L();
	lambdat<T> i1Lat=ep.p(i1)->Lt();
	lambda<T> i2La=ep.p(i2)->L();
	lambdat<T> qLat=ep.ref()->Lt();
	lambda<T> qLa=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*pow(qLat*i1Lat,2)*(i1La*i2La)/(complex<T>(0,2)*(qLat*p1flatc.Lt())*(qLat*pnflatc.Lt())*ep.sp(i1,i2));
}

template <class T> complex<T> (*A2QM2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch(hc){
	case 70326: //0x080F1006
		return &A2QM2l_7m_eval<0,1,2,3>;
		break;
	case 54328: //0x060F1008
		return &A2QM2l_7p_eval<0,1,2,3>;
		break;
	case 51516: //0x06080F10
		return &A2QM2l_7m_eval<1,2,3,0>;
		break;
	case 66716: //0x08060F10
		return &A2QM2l_7p_eval<1,2,3,0>;
		break;
	case 130575: //0x1006080F
		return &A2QM2l_7m_eval<2,3,0,1>;
		break;
	case 131335: //0x1008060F
		return &A2QM2l_7p_eval<2,3,0,1>;
		break;
	case 126528: //0x0F100608
		return &A2QM2l_7m_eval<3,0,1,2>;
		break;
	case 126566: //0x0F100806
		return &A2QM2l_7p_eval<3,0,1,2>;
		break;

	case 70325: //0x080F1005
	case 54327: //0x060F1007
		return &A2QM2l_6_eval<0,1,2,3>;
		break;
	case 43516: //0x05080F10
	case 58716: //0x07060F10
		return &A2QM2l_6_eval<1,2,3,0>;
		break;
	case 130175: //0x1005080F
	case 130935: //0x1007060F
		return &A2QM2l_6_eval<2,3,0,1>;
		break;
	case 126508: //0x0F100508
	case 126546: //0x0F100706
		return &A2QM2l_6_eval<3,0,1,2>;
		break;

	case 62326: //0x070F1006
	case 46328: //0x050F1008
		return &A2QM2l_5_eval<0,1,2,3>;
		break;
	case 51116: //0x06070F10
	case 66316: //0x08050F10
		return &A2QM2l_5_eval<1,2,3,0>;
		break;
	case 130555: //0x1006070F
	case 131315: //0x1008050F
		return &A2QM2l_5_eval<2,3,0,1>;
		break;
	case 126527: //0x0F100607
	case 126565: //0x0F100805
		return &A2QM2l_5_eval<3,0,1,2>;
		break;

	case 70706: //0x08100F06
		return &A2QM2l_1m_eval<0,1,2,3>;
		break;
	case 54708: //0x06100F08
		return &A2QM2l_1p_eval<0,1,2,3>;
		break;
	case 51535: //0x0608100F
		return &A2QM2l_1m_eval<1,2,3,0>;
		break;
	case 66735: //0x0806100F
		return &A2QM2l_1p_eval<1,2,3,0>;
		break;
	case 122576: //0x0F060810
		return &A2QM2l_1m_eval<2,3,0,1>;
		break;
	case 123336: //0x0F080610
		return &A2QM2l_1p_eval<2,3,0,1>;
		break;
	case 134128: //0x100F0608
		return &A2QM2l_1m_eval<3,0,1,2>;
		break;
	case 134166: //0x100F0806
		return &A2QM2l_1p_eval<3,0,1,2>;
		break;

	case 70705: //0x08100F05
	case 54707: //0x06100F07
		return &A2QM2l_2_eval<0,1,2,3>;
		break;
	case 43535: //0x0508100F
	case 58735: //0x0706100F
		return &A2QM2l_2_eval<1,2,3,0>;
		break;
	case 122176: //0x0F050810
	case 122936: //0x0F070610
		return &A2QM2l_2_eval<2,3,0,1>;
		break;
	case 134108: //0x100F0508
	case 134146: //0x100F0706
		return &A2QM2l_2_eval<3,0,1,2>;
		break;

	case 62706: //0x07100F06
	case 46708: //0x05100F08
		return &A2QM2l_3_eval<0,1,2,3>;
		break;
	case 51135: //0x0607100F
	case 66335: //0x0805100F
		return &A2QM2l_3_eval<1,2,3,0>;
		break;
	case 122556: //0x0F060710
	case 123316: //0x0F080510
		return &A2QM2l_3_eval<2,3,1,0>;
		break;
	case 134127: //0x100F0607
	case 134165: //0x100F0805
		return &A2QM2l_3_eval<3,0,1,2>;
		break;

	case 62705: //0x07100F05
		return &A2QM2l_8m_eval<0,1,2,3>;
		break;
	case 46707: //0x05100F07
		return &A2QM2l_8p_eval<0,1,2,3>;
		break;
	case 43135: //0x0507100F
		return &A2QM2l_8m_eval<1,2,3,0>;
		break;
	case 58335: //0x0705100F
		return &A2QM2l_8p_eval<1,2,3,0>;
		break;
	case 122156: //0x0F050710
		return &A2QM2l_8m_eval<2,3,0,1>;
		break;
	case 122916: //0x0F070510
		return &A2QM2l_8p_eval<2,3,0,1>;
		break;
	case 134107: //0x100F0507
		return &A2QM2l_8m_eval<3,0,1,2>;
		break;
	case 134145: //0x100F0705
		return &A2QM2l_8p_eval<3,0,1,2>;
		break;

	case 62325: //0x070F1005
		return &A2QM2l_4m_eval<0,1,2,3>;
		break;
	case 46327: //0x050F1007
		return &A2QM2l_4p_eval<0,1,2,3>;
		break;
	case 43116: //0x05070F10
		return &A2QM2l_4m_eval<1,2,3,0>;
		break;
	case 58316: //0x07050F10
		return &A2QM2l_4p_eval<1,2,3,0>;
		break;
	case 130155: //0x1005070F
		return &A2QM2l_4m_eval<2,3,0,1>;
		break;
	case 130915: //0x1007050F
		return &A2QM2l_4p_eval<2,3,0,1>;
		break;
	case 126507: //0x0F100507
		return &A2QM2l_4m_eval<3,0,1,2>;
		break;
	case 126545: //0x0F100705
		return &A2QM2l_4p_eval<3,0,1,2>;
		break;

	default:// We return zero for all other helicity combinations
		return &ZeroF;
	}
}


/*
 *
 *
 * The 2 massive quarks and one gluon and one photon amplitudes
 *
 *
 */

template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s1g1y_help_1(const eval_param<T>& ep, const complex<T>& mass){
	return complex<T>(0,1)*mass*ep.spb(i1,i2)/(ep.spa(i1,i2)*(-T(2)*ep.sp(i0,i1)));
}
template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s1g1y_help_2(const eval_param<T>& ep, const complex<T>& mass){
	return complex<T>(0,1)*mass*ep.spa(i1,i2)/(ep.spb(i1,i2)*(-T(2)*ep.sp(i0,i1)));
}
template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s1g1y_help_3(const eval_param<T>& ep, const complex<T>& mass){
	return complex<T>(0,1)*mass*pow(ep.spb(i1,i3),2)/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}
template <int i0, int i1, int i2, int i3, class T> inline complex<T>  A2s1g1y_help_4(const eval_param<T>& ep, const complex<T>& mass){
	return complex<T>(0,1)*mass*pow(ep.spa(i1,i3),2)/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y1_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> qL=ep.ref()->L();

	return ((pnL*qL)/(p1L*qL))*(A2s1g1y_help_1<i0,i1,i2,i3>(ep,mass)+A2s1g1y_help_1<i0,i2,i1,i3>(ep,mass));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y2_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> qL=ep.ref()->L();

	return ((p1L*qL)/(pnL*qL))*(A2s1g1y_help_1<i0,i1,i2,i3>(ep,mass)+A2s1g1y_help_1<i0,i2,i1,i3>(ep,mass));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y3m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -((p1L*pnL)/eval_param<T>::mass(masses.p(i0)))*(A2s1g1y_help_1<i0,i1,i2,i3>(ep,mass)+A2s1g1y_help_1<i0,i2,i1,i3>(ep,mass));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y3p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	lambda<T> p1L=la(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((p1L*pnL)/eval_param<T>::mass(masses.p(i0)))*(A2s1g1y_help_1<i0,i1,i2,i3>(ep,mass)+A2s1g1y_help_1<i0,i2,i1,i3>(ep,mass));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y46_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambda<T> qL=ep.ref()->L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambda<T> i1L=ep.p(i1)->L();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();
	const lambda<T> i2L=ep.p(i2)->L();
	const lambdat<T> i2Lt=ep.p(i2)->Lt();

	return complex<T>(0,-1)*(ep.spba(i1,i3,i2)*
					(ep.spba(i1,i3,i2)*(qLt*p1Lt)/(T(2)*ep.sp(i1,i2))
							+(qLt*i1Lt)*(i1Lt*p1Lt)/ep.spb(i1,i2))
					/((qLt*pnLt)*(-T(2)*ep.sp(i0,i1)))
					-ep.spab(i2,i3,i1)*
					(ep.spab(i2,i3,i1)*(qL*pnL)/(T(2)*ep.sp(i2,i1))
							-(qL*i2L)*(i2L*pnL)/ep.spa(i2,i1))
					/((qL*p1L)*(-T(2)*ep.sp(i0,i2))));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y57_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambda<T> qL=ep.ref()->L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambda<T> i1L=ep.p(i1)->L();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();
	const lambda<T> i2L=ep.p(i2)->L();
	const lambdat<T> i2Lt=ep.p(i2)->Lt();

	return complex<T>(0,-1)*(ep.spba(i1,i3,i2)*
			(ep.spba(i1,i3,i2)*(qLt*pnLt)/(T(2)*ep.sp(i1,i2))
					-(qLt*i1Lt)*(i1Lt*pnLt)/ep.spb(i1,i2))
			/((qLt*p1Lt)*(-T(2)*ep.sp(i0,i1)))
			-ep.spab(i2,i3,i1)*
			(ep.spab(i2,i3,i1)*(qL*p1L)/(T(2)*ep.sp(i2,i1))
					+(qL*i2L)*(i2L*p1L)/ep.spa(i2,i1))
			/((qL*pnL)*(-T(2)*ep.sp(i0,i2))));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y89m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();


	return eval_param<T>::mass(masses.p(i0))*(ep.spba(i1,i0,i2)*pow(qLt*ep.p(i1)->Lt(),2)/(ep.sp(i0,i1)*ep.spb(i1,i2)*(qLt*p1Lt)*(pnLt*qLt))
									-ep.spba(i1,i3,i2)*pow(qLt*ep.p(i1)->Lt(),2)/(ep.sp(i3,i1)*ep.spb(i1,i2)*(qLt*pnLt)*(p1Lt*qLt)))/complex<T>(0,-2);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y89p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return eval_param<T>::mass(masses.p(i0))*(ep.spba(i1,i0,i2)*pow((qLt*ep.p(i1)->Lt()),2)/(ep.sp(i0,i1)*ep.spb(i1,i2)*(qLt*p1Lt)*(pnLt*qLt))
					-ep.spba(i1,i3,i2)*pow((qLt*ep.p(i1)->Lt()),2)/(ep.sp(i3,i1)*ep.spb(i1,i2)*(qLt*pnLt)*(p1Lt*qLt)))/complex<T>(0,2);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y10_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return -((p1Lt*qLt)/(pnLt*qLt))*(A2s1g1y_help_2<i0,i1,i2,i3>(ep,mass)+A2s1g1y_help_2<i0,i2,i1,i3>(ep,mass));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y11_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return -((pnLt*qLt)/(p1Lt*qLt))*(A2s1g1y_help_2<i0,i1,i2,i3>(ep,mass)+A2s1g1y_help_2<i0,i2,i1,i3>(ep,mass));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y12m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((p1Lt*pnLt)/eval_param<T>::mass(masses.p(i0)))*(A2s1g1y_help_2<i0,i1,i2,i3>(ep,mass)+A2s1g1y_help_2<i0,i2,i1,i3>(ep,mass));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y12p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -((p1Lt*pnLt)/eval_param<T>::mass(masses.p(i0)))*(A2s1g1y_help_2<i0,i1,i2,i3>(ep,mass)+A2s1g1y_help_2<i0,i2,i1,i3>(ep,mass));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y1314m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*(ep.spab(i2,i3,i1)*pow((qL*ep.p(i2)->L()),2)/(ep.sp(i3,i2)*ep.spa(i2,i1)*(qL*pnL)*(p1L*qL))
			-ep.spab(i2,i0,i1)*pow((qL*ep.p(i2)->L()),2)/(ep.sp(i0,i2)*ep.spa(i2,i1)*(qL*p1L)*(pnL*qL)))/complex<T>(0,-2);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y1314p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));


	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i0))*(ep.spab(i2,i3,i1)*pow((qL*ep.p(i2)->L()),2)/(ep.sp(i3,i2)*ep.spa(i2,i1)*(qL*pnL)*(p1L*qL))
						-ep.spab(i2,i0,i1)*pow((qL*ep.p(i2)->L()),2)/(ep.sp(i0,i2)*ep.spa(i2,i1)*(qL*p1L)*(pnL*qL)))/complex<T>(0,2);
}


template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y15_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambda<T> qL=ep.ref()->L();

	return ((p1L*qL)/(pnL*qL))*A2s1g1y_help_3<i0,i1,i2,i3>(ep,mass);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y16_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambda<T> p1L=la((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambda<T> pnL=la(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return -((p1L*pnL)/eval_param<T>::mass(masses.p(i0)))*A2s1g1y_help_3<i0,i1,i2,i3>(ep,mass);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y17_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> qLt=ep.ref()->Lt();

	return -((p1Lt*qLt)/(pnLt*qLt))*A2s1g1y_help_4<i0,i1,i2,i3>(ep,mass);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y18_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const lambdat<T> p1Lt=lat((ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P()));
	const lambdat<T> pnLt=lat(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());

	return ((p1Lt*pnLt)/eval_param<T>::mass(masses.p(i0)))*A2s1g1y_help_4<i0,i1,i2,i3>(ep,mass);
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y19_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

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

	return complex<T>(0,-1)*ep.spab(i3,i0,i1)*((p1Lt*i1Lt)*(i3L*pnL)-mass*(qL*ep.p(i3)->L())*(ep.p(i1)->Lt()*qLt)/((qL*p1L)*(pnLt*qLt)))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y20_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

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

	return complex<T>(0,1)*ep.spab(i3,i0,i1)*((p1L*i3L)*(i1Lt*pnLt)-mass*(qL*ep.p(i3)->L())*(ep.p(i1)->Lt()*qLt)/((qLt*p1Lt)*pnL*qL))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y23m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambda<T> qL=ep.ref()->L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*ep.spab(i3,i0,i1)*((p1Lt*i1Lt)*(ep.p(i3)->L()*qL)/(pnL*qL)-(qL*ep.p(i3)->L())*(i1Lt*pnLt)/(qL*p1L))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y23p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambda<T> qL=ep.ref()->L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambdat<T> i1Lt=ep.p(i1)->Lt();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i0))*ep.spab(i3,i0,i1)*((p1Lt*i1Lt)*(ep.p(i3)->L()*qL)/(pnL*qL)-(qL*ep.p(i3)->L())*(i1Lt*pnLt)/(qL*p1L))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y24m_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambda<T> i1L=ep.p(i1)->L();

	return complex<T>(0,1)*eval_param<T>::mass(masses.p(i0))*ep.spba(i3,i0,i1)*((p1L*i1L)*(ep.p(i3)->Lt()*qLt)/(pnLt*qLt)-(qLt*ep.p(i3)->Lt())*(i1L*pnL)/(qLt*p1Lt))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}

template <int i0, int i1, int i2, int i3, class T> complex<T>  A2QM1g1y24p_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i0));

	const Cmom<T> p1(ep.mom(i0)-(mass/(T(2)*(ep.p(i0)->P()*ep.ref()->P())))*ep.ref()->P());
	const Cmom<T> pn(ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P());
	const lambdat<T> p1Lt=p1.Lt();
	const lambda<T> p1L=p1.L();
	const lambdat<T> pnLt=pn.Lt();
	const lambda<T> pnL=pn.L();
	const lambdat<T> qLt=ep.ref()->Lt();
	const lambda<T> i1L=ep.p(i1)->L();

	return complex<T>(0,-1)*eval_param<T>::mass(masses.p(i0))*ep.spba(i3,i0,i1)*((p1L*i1L)*(ep.p(i3)->Lt()*qLt)/(pnLt*qLt)-(qLt*ep.p(i3)->Lt())*(i1L*pnL)/(qLt*p1Lt))
					/(T(4)*ep.sp(i0,i1)*ep.sp(i0,i3));
}


template <class T> complex<T> (*A2QM1g1y_Tree_Ptr_eval(int hc))(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// case 2q_massive
	switch (hc) {
	case 0x8315://995: // Qb+,g+,g+,Q-
	case 0x6317:
	case 0x8135:
	case 0x6137:
		return &A2QM1g1y1_eval<0,1,2,3>;
		break;
	case 0x5831://790: //(Q3-,Qb3+,g+,g+)
	case 0x7631:
	case 0x5813:
	case 0x7613:
		return &A2QM1g1y1_eval<1,2,3,0>;
		break;
	case 0x1583://855: // (g+,Qb3-,Q3+,g+)
	case 0x1763:
	case 0x3581:
	case 0x3761:
		return &A2QM1g1y1_eval<2,3,0,1>;
		break;
	case 0x3158://1245: // (g+,g+,Qb3-,Q3+)
	case 0x3176:
	case 0x1358:
	case 0x1376:
		return &A2QM1g1y1_eval<3,0,1,2>;
		break;


	case 0x7316://1210: // Qb-,g+,g+,Q+
	case 0x5318:
	case 0x7136:
	case 0x5138:
		return &A2QM1g1y2_eval<0,1,2,3>;
		break;
	case 0x6731://785: //(Q3+,Qb3-,g+,g+)
	case 0x8531:
	case 0x6713:
	case 0x8513:
		return &A2QM1g1y2_eval<1,2,3,0>;
		break;
	case 0x1673://825: // (g+,Qb3+,Q3-,g+)
	case 0x1853:
	case 0x3671:
	case 0x3851:
		return &A2QM1g1y2_eval<2,3,0,1>;
		break;
	case 0x1367://1065: // (g+,g+,Qb3+,Q3-)
	case 0x1385:
	case 0x3167:
	case 0x3185:
		return &A2QM1g1y2_eval<3,0,1,2>;
		break;

	case 0x5381: // (Q-,g+,Qb+,g+) NEW
	case 0x7361: // (Qb-,g+,Q+,g+) NEW
	case 0x5183: // (Q-,g+,Qb+,g+) NEW
	case 0x7163: // (Qb-,g+,Q+,g+) NEW
		return &A2QM1g1y15_eval<0,1,2,3>;
		break;
	case 0x3518: // (g+,Q-,g+,Qb+) NEW
	case 0x3716: // (g+,Qb-,g+,Q+) NEW
	case 0x1538: // (g+,Q-,g+,Qb+) NEW
	case 0x1736: // (g+,Qb-,g+,Q+) NEW
		return &A2QM1g1y15_eval<1,2,3,0>;
		break;
	case 0x8351: // (Qb+,g+,Q-,g+) NEW
	case 0x6371: // (Q+,g+,Qb-,g+) NEW
	case 0x8153: // (Qb+,g+,Q-,g+) NEW
	case 0x6173: // (Q+,g+,Qb-,g+) NEW
		return &A2QM1g1y15_eval<2,3,0,1>;
		break;
	case 0x3815: // (g+,Qb+,g+,Q-) NEW
	case 0x3617: // (g+,Q+,g+,Qb-) NEW
	case 0x1835: // (g+,Qb+,g+,Q-) NEW
	case 0x1637: // (g+,Q+,g+,Qb-) NEW
		return &A2QM1g1y15_eval<3,0,1,2>;
		break;


	case 0x5317: //994: // Qb-g+g+Q-
	case 0x5137:
		return &A2QM1g1y3p_eval<0,1,2,3>;
		break;
	case 0x7315:
	case 0x7135:
		return &A2QM1g1y3m_eval<0,1,2,3>;
		break;
	case 0x7531: //784: //(Q3-,Qb3-,g+,g+)
	case 0x7513:
		return &A2QM1g1y3p_eval<1,2,3,0>;
		break;
	case 0x5731:
	case 0x5713:
		return &A2QM1g1y3m_eval<1,2,3,0>;
		break;
	case 0x3751: //819: // (g+,Qb3-,Q3-,g+)
	case 0x1753:
		return &A2QM1g1y3p_eval<2,3,0,1>;
		break;
	case 0x3571:
	case 0x1573:
		return &A2QM1g1y3m_eval<2,3,0,1>;
		break;
	case 0x1375: //1029: // (g+,g+,Qb3-,Q3-)
	case 0x3175:
		return &A2QM1g1y3p_eval<3,0,1,2>;
		break;
	case 0x3157:
	case 0x1357:
		return &A2QM1g1y3m_eval<3,0,1,2>;
		break;

	case 0x7351: // (Qb-,g+,Q-,g+) NEW
	case 0x7153: // (Qb-,g+,Q-,g+) NEW
		return &A2QM1g1y16_eval<0,1,2,3>;
		break;
	case 0x3715: // (g+,Qb-,g+,Q-) NEW
	case 0x1735: // (g+,Qb-,g+,Q-) NEW
		return &A2QM1g1y16_eval<1,2,3,0>;
		break;
	case 0x5371: // (Q-,g+,Qb-,g+) NEW
	case 0x5173: // (Q-,g+,Qb-,g+) NEW
		return &A2QM1g1y16_eval<2,3,0,1>;
		break;
	case 0x3517: // (g+,Q-,g+,Qb-) NEW
	case 0x1537: // (g+,Q-,g+,Qb-) NEW
		return &A2QM1g1y16_eval<3,0,1,2>;
		break;

	case 0x8305: //887: // Qb+g+g-Q-
	case 0x6307:
	case 0x8125:
	case 0x6127:
		return &A2QM1g1y46_eval<0,1,2,3>;
		break;
	case 0x8035: // Qb+g-g+Q-
	case 0x6037:
	case 0x8215:
	case 0x6217:
		return &A2QM1g1y46_eval<0,2,1,3>;
		break;
	case 0x5830: //142: //(Q3-,Qb3+,g+,g-)
	case 0x7630:
	case 0x5812:
	case 0x7612:
		return &A2QM1g1y46_eval<1,2,3,0>;
		break;
	case 0x5803: //(Q3-,Qb3+,g-,g+)
	case 0x7603:
	case 0x5821:
	case 0x7621:
		return &A2QM1g1y46_eval<1,3,2,0>;
		break;
	case 0x0583: //852: // (g-,Qb3-,Q3+,g+)
	case 0x0763:
	case 0x2581:
	case 0x2761:
		return &A2QM1g1y46_eval<2,3,0,1>;
		break;
	case 0x3580: // (g+,Qb3-,Q3+,g-)
	case 0x3760:
	case 0x1582:
	case 0x1762:
		return &A2QM1g1y46_eval<2,0,3,1>;
		break;
	case 0x3058: //1227: //(g+,g-,Qb3-,Q3+)
	case 0x3076:
	case 0x1258:
	case 0x1276:
		return &A2QM1g1y46_eval<3,0,1,2>;
		break;
	case 0x0358: //(g-,g+,Qb3-,Q3+)
	case 0x0376:
	case 0x2158:
	case 0x2176:
		return &A2QM1g1y46_eval<3,1,0,2>;
		break;

	case 0x8152: //887: // Qb+g+Q-g-
	case 0x6172:
	case 0x8350: //887: // Qb+g+Q-g-
	case 0x6370:
		return &A2QM1g1y19_eval<0,1,2,3>;
		break;
	case 0x2815: //142: //(g-,Q3-,Qb3+,g+)
	case 0x2617:
	case 0x0835: //142: //(g-,Q3-,Qb3+,g+)
	case 0x0637:
		return &A2QM1g1y19_eval<1,2,3,0>;
		break;
	case 0x5281: //852: // (Qb3-,g-,Q3+,g+)
	case 0x7261:
	case 0x5083: //852: // (Qb3-,g-,Q3+,g+)
	case 0x7063:
		return &A2QM1g1y19_eval<2,3,0,1>;
		break;
	case 0x1528: //1227: //(g+,Qb3-,g-,Q3+)
	case 0x1726:
	case 0x3508: //1227: //(g+,Qb3-,g-,Q3+)
	case 0x3706:
		return &A2QM1g1y19_eval<3,0,1,2>;
		break;

	case 0x7306: //1102: // Qb-g+g-Q+
	case 0x5308:
	case 0x7126:
	case 0x5128:
		return &A2QM1g1y57_eval<0,1,2,3>;
		break;
	case 0x7036: //1192: // Qb-g-g+Q+
	case 0x5038:
	case 0x7216:
	case 0x5218:
		return &A2QM1g1y57_eval<0,2,1,3>;
		break;
	case 0x6730: //137: // (Q3+,Qb3-,g+,g-)
	case 0x8530:
	case 0x6712:
	case 0x8512:
		return &A2QM1g1y57_eval<1,2,3,0>;
		break;
	case 0x6703: //
	case 0x8503:
	case 0x6721: //
	case 0x8521:
		return &A2QM1g1y57_eval<1,3,2,0>;
		break;
	case 0x0673: //
	case 0x0853:
	case 0x2671:
	case 0x2851:
		return &A2QM1g1y57_eval<2,3,0,1>;
		break;
	case 0x3670: //177: // (g+,Qb3+,Q3-,g-)
	case 0x3850:
	case 0x1672:
	case 0x1852:
		return &A2QM1g1y57_eval<2,0,3,1>;
		break;
	case 0x3067: //1047: // (g+,g-,Q3+,Qb3-)
	case 0x3085:
	case 0x1267:
	case 0x1285:
		return &A2QM1g1y57_eval<3,0,1,2>;
		break;
	case 0x0367: //
	case 0x0385:
	case 0x2167:
	case 0x2185:
		return &A2QM1g1y57_eval<3,1,0,2>;
		break;

	case 0x7162: //1102: // Qb-g+Q+g-
	case 0x5182:
	case 0x7360: //1102: // Qb-g+Q+g-
	case 0x5380:
		return &A2QM1g1y20_eval<0,1,2,3>;
		break;
	case 0x2716: //137: // (Q3+,Qb3-,g+,g-)
	case 0x2518:
	case 0x0736: //137: // (Q3+,Qb3-,g+,g-)
	case 0x0538:
		return &A2QM1g1y20_eval<1,2,3,0>;
		break;
	case 0x6271: //
	case 0x8251:
	case 0x6073: //
	case 0x8053:
		return &A2QM1g1y20_eval<2,3,0,1>;
		break;
	case 0x1627: //1047: // (g+,g-,Q3+,Qb3-)
	case 0x1825:
	case 0x3607: //1047: // (g+,g-,Q3+,Qb3-)
	case 0x3805:
		return &A2QM1g1y20_eval<3,0,1,2>;
		break;


	case 0x7305: //886: //Qb-,g+,g-,Q-
	case 0x7125:
		return &A2QM1g1y89m_eval<0,1,2,3>;
		break;
	case 0x7035: //976: //Qb-,g-,g+,Q-
	case 0x7215:
		return &A2QM1g1y89m_eval<0,2,1,3>;
		break;
	case 0x5307:
	case 0x5127:
		return &A2QM1g1y89p_eval<0,1,2,3>;
		break;
	case 0x5037:
	case 0x5217:
		return &A2QM1g1y89p_eval<0,2,1,3>;
		break;
	case 0x5730: //136: //(Q3-,Qb3-,g+,g-)
	case 0x5712:
		return &A2QM1g1y89m_eval<1,2,3,0>;
		break;
	case 0x5703:
	case 0x5721:
		return &A2QM1g1y89m_eval<1,3,2,0>;
		break;
	case 0x7530:
	case 0x7512:
		return &A2QM1g1y89p_eval<1,2,3,0>;
		break;
	case 0x7503:
	case 0x7521:
		return &A2QM1g1y89p_eval<1,3,2,0>;
		break;
	case 0x0573: //816: //(g-,Qb3-,Q3-,g+)
	case 0x2571:
		return &A2QM1g1y89m_eval<2,3,0,1>;
		break;
	case 0x3570: //171: //(g+,Q3-,Qb3-,g-)
	case 0x1572:
		return &A2QM1g1y89m_eval<2,0,3,1>;
		break;
	case 0x0753:
	case 0x2751:
		return &A2QM1g1y89p_eval<2,3,0,1>;
		break;
	case 0x3750:
	case 0x1752:
		return &A2QM1g1y89p_eval<2,0,3,1>;
		break;
	case 0x3057: //1011: //(g+,g-,Q3-,Qb3-)
	case 0x1257:
		return &A2QM1g1y89m_eval<3,0,1,2>;
		break;
	case 0x0357: //1026: //(g-,g+,Qb3-,Q3-)
	case 0x2157:
		return &A2QM1g1y89m_eval<3,1,0,2>;
		break;
	case 0x3075:
	case 0x1275:
		return &A2QM1g1y89p_eval<3,0,1,2>;
		break;
	case 0x0375:
	case 0x2175:
		return &A2QM1g1y89p_eval<3,1,0,2>;
		break;

	case 0x7152: //886: //Qb-,g+,Q-,g-
	case 0x7350: //886: //Qb-,g+,Q-,g-
		return &A2QM1g1y24m_eval<2,3,0,1>;
		break;
	case 0x2715: //
	case 0x0735: //
		return &A2QM1g1y24m_eval<3,0,1,2>;
		break;
	case 0x5271: //
	case 0x5073: //
		return &A2QM1g1y24m_eval<0,1,2,3>;
		break;
	case 0x1527: //
	case 0x3507: //
		return &A2QM1g1y24m_eval<1,2,3,0>;
		break;

	case 0x5172:
	case 0x5370:
		return &A2QM1g1y24p_eval<2,3,0,1>;
		break;
	case 0x2517:
	case 0x0537:
		return &A2QM1g1y24p_eval<3,0,1,2>;
		break;
	case 0x7251:
	case 0x7053:
		return &A2QM1g1y24p_eval<0,1,2,3>;
		break;
	case 0x1725:
	case 0x3705:
		return &A2QM1g1y24p_eval<1,2,3,0>;
		break;

	case 0x8306: //1103: //Qb+,g+,g-,Q+
	case 0x8126:
		return &A2QM1g1y1314m_eval<0,1,2,3>;
		break;
	case 0x8036: //1193: //Qb+,g-,g+,Q+
	case 0x8216:
		return &A2QM1g1y1314m_eval<0,2,1,3>;
		break;
	case 0x6308:
	case 0x6128:
		return &A2QM1g1y1314p_eval<0,1,2,3>;
		break;
	case 0x6038:
	case 0x6218:
		return &A2QM1g1y1314p_eval<0,2,1,3>;
		break;
	case 0x6830: //143: //(Q3+,Qb3+,g+,g-)
	case 0x6812:
		return &A2QM1g1y1314m_eval<1,2,3,0>;
		break;
	case 0x6803: //
	case 0x6821: //
		return &A2QM1g1y1314m_eval<1,3,2,0>;
		break;
	case 0x8630:
	case 0x8612:
		return &A2QM1g1y1314p_eval<1,2,3,0>;
		break;
	case 0x8603:
	case 0x8621:
		return &A2QM1g1y1314p_eval<1,3,2,0>;
		break;
	case 0x0683: //858: //(g-,Qb3+,Q3+,g+)
	case 0x2681:
		return &A2QM1g1y1314m_eval<2,3,0,1>;
		break;
	case 0x3680: //213: //(g+,Q3+,Qb3+,g-)
	case 0x1682:
		return &A2QM1g1y1314m_eval<2,0,3,1>;
		break;
	case 0x0863:
	case 0x2861:
		return &A2QM1g1y1314p_eval<2,3,0,1>;
		break;
	case 0x3860:
	case 0x1862:
		return &A2QM1g1y1314p_eval<2,0,3,1>;
		break;
	case 0x3068: //1263: //(g+,g-,Q3+,Qb3+)
	case 0x1268:
		return &A2QM1g1y1314m_eval<3,0,1,2>;
		break;
	case 0x0368: //1278: //(g-,g+,Qb3+,Q3+)
	case 0x2168:
		return &A2QM1g1y1314m_eval<3,1,0,2>;
		break;
	case 0x3086:
	case 0x1286:
		return &A2QM1g1y1314p_eval<3,0,1,2>;
		break;
	case 0x0386:
	case 0x2186:
		return &A2QM1g1y1314p_eval<3,1,0,2>;
		break;

	case 0x8360: //1103: //Qb+,g+,Q+,g-
	case 0x8162: //1103: //Qb+,g+,Q+,g-
		return &A2QM1g1y23p_eval<0,1,2,3>;
		break;
	case 0x0836: //
	case 0x2816: //
		return &A2QM1g1y23p_eval<1,2,3,0>;
		break;
	case 0x6083: //
	case 0x6281: //
		return &A2QM1g1y23p_eval<2,3,0,1>;
		break;
	case 0x3608: //
	case 0x1628: //
		return &A2QM1g1y23p_eval<3,0,1,2>;
		break;

	case 0x6380:
	case 0x6182:
		return &A2QM1g1y23m_eval<0,1,2,3>;
		break;
	case 0x0638:
	case 0x2618:
		return &A2QM1g1y23m_eval<1,2,3,0>;
		break;
	case 0x8063:
	case 0x8261:
		return &A2QM1g1y23m_eval<2,3,0,1>;
		break;
	case 0x3806:
	case 0x1826:
		return &A2QM1g1y23m_eval<3,0,1,2>;
		break;

	case 0x8205: //869: // Qb+g-g-Q-
	case 0x6207:
	case 0x8025:
	case 0x6027:
		return &A2QM1g1y10_eval<0,1,2,3>;
		break;
	case 0x5820: //34: // (Q3-,Qb3+,g-,g-)
	case 0x7620:
	case 0x5802:
	case 0x7602:
		return &A2QM1g1y10_eval<1,2,3,0>;
		break;
	case 0x0582: //204: // (g-,Qb3-,Q3+,g-)
	case 0x0762:
	case 0x2580:
	case 0x2760:
		return &A2QM1g1y10_eval<2,3,0,1>;
		break;
	case 0x2058: //1224: // (g-,g-,Qb3-,Q3+)
	case 0x2076:
	case 0x0258:
	case 0x0276:
		return &A2QM1g1y10_eval<3,0,1,2>;
		break;

	case 0x7206: //1084: // Qb-g-g-Q+
	case 0x5208:
	case 0x7026:
	case 0x5028:
		return &A2QM1g1y11_eval<0,1,2,3>;
		break;
	case 0x6720: //29: // (Q3+,Qb3-,g-,g-)
	case 0x8520:
	case 0x6702:
	case 0x8502:
		return &A2QM1g1y11_eval<1,2,3,0>;
		break;
	case 0x2670: //174: // (g-,Q3+,Qb3-,g-)
	case 0x2850:
	case 0x0672:
	case 0x0852:
		return &A2QM1g1y11_eval<2,3,0,1>;
		break;
	case 0x2067: //1044: //(g-,g-,Q3+,Qb3-)
	case 0x2085:
	case 0x0267:
	case 0x0285:
		return &A2QM1g1y11_eval<3,0,1,2>;
		break;

	case 0x5280: // (Q-,g-,Qb+,g-) NEW
	case 0x7260: // (Qb-,g-,Q+,g-) NEW
	case 0x5082: // (Q-,g-,Qb+,g-) NEW
	case 0x7062: // (Qb-,g-,Q+,g-) NEW
		return &A2QM1g1y17_eval<2,3,0,1>;
		break;
	case 0x2508: // (g-,Q-,g-,Qb+) NEW
	case 0x2706: // (g-,Qb-,g-,Q+) NEW
	case 0x0528: // (g-,Q-,g-,Qb+) NEW
	case 0x0726: // (g-,Qb-,g-,Q+) NEW
		return &A2QM1g1y17_eval<3,0,1,2>;
		break;
	case 0x8250: // (Qb+,g-,Q-,g-) NEW
	case 0x6270: // (Q+,g-,Qb-,g-) NEW
	case 0x8052: // (Qb+,g-,Q-,g-) NEW
	case 0x6072: // (Q+,g-,Qb-,g-) NEW
		return &A2QM1g1y17_eval<0,1,2,3>;
		break;
	case 0x2805: // (g-,Qb+,g-,Q-) NEW
	case 0x2607: // (g-,Q+,g-,Qb-) NEW
	case 0x0825: // (g-,Qb+,g-,Q-) NEW
	case 0x0627: // (g-,Q+,g-,Qb-) NEW
		return &A2QM1g1y17_eval<1,2,3,0>;
		break;

	case 0x8206: //1085: // Qb+g-g-Q+
	case 0x8026:
		return &A2QM1g1y12m_eval<0,1,2,3>;
		break;
	case 0x6208:
	case 0x6028:
		return &A2QM1g1y12p_eval<0,1,2,3>;
		break;
	case 0x6820: //35: // (Q3+,Qb3+,g-,g-)
	case 0x6802:
		return &A2QM1g1y12m_eval<1,2,3,0>;
		break;
	case 0x8620:
	case 0x8602:
		return &A2QM1g1y12p_eval<1,2,3,0>;
		break;
	case 0x2680: //210: // (g-,Q3+,Qb3+,g-)
	case 0x0682:
		return &A2QM1g1y12m_eval<2,3,0,1>;
		break;
	case 0x2860:
	case 0x0862:
		return &A2QM1g1y12p_eval<2,3,0,1>;
		break;
	case 0x2068: //1260: // (g-,g-,Q3+,Qb3+)
	case 0x0268:
		return &A2QM1g1y12m_eval<3,0,1,2>;
		break;
	case 0x2086:
	case 0x0286:
		return &A2QM1g1y12p_eval<3,0,1,2>;
		break;

	case 0x8260: // (Qb+,g-,Q+,g-) NEW
	case 0x8062: // (Qb+,g-,Q+,g-) NEW
		return &A2QM1g1y18_eval<0,1,2,3>;
		break;
	case 0x2806: // (g-,Q+-,g-,Q+) NEW
	case 0x0826: // (g-,Q+-,g-,Q+) NEW
		return &A2QM1g1y18_eval<1,2,3,0>;
		break;
	case 0x6280: // (Q+,g-,Qb+,g-) NEW
	case 0x6082: // (Q+,g-,Qb+,g-) NEW
		return &A2QM1g1y18_eval<2,3,0,1>;
		break;
	case 0x2608: // (g-,Q+,g-,Qb+) NEW
	case 0x0628: // (g-,Q+,g-,Qb+) NEW
		return &A2QM1g1y18_eval<3,0,1,2>;
		break;


	case 0x8136:// Qb+g+g+Q+
	case 0x6813:
	case 0x3681:
	case 0x1368:
	case 0x6138:
	case 0x8613:
	case 0x3861:
	case 0x1386:
	case 0x7025:// Qb-g-g-Q-
	case 0x5702:
	case 0x2570:
	case 0x0257:
	case 0x5027:
	case 0x7502:
	case 0x2750:
	case 0x0275:
	case 0x8316:// Qb+g+g+Q+
	case 0x6831:
	case 0x1683:
	case 0x3168:
	case 0x6318:
	case 0x8631:
	case 0x1863:
	case 0x3186:
	case 0x7205:// Qb-g-g-Q-
	case 0x5720:
	case 0x0572:
	case 0x2057:
	case 0x5207:
	case 0x7520:
	case 0x0752:
	case 0x2075:
		return &ZeroF;
		break;

	default:// We return zero for all other helicity combinations
		return 0;
	}
}




template <class T> complex<T> (*A2QM1g2l_Tree_Ptr_eval(int hc))(const eval_param<T>& ep, const mass_param_coll& masses)
{
	// case 2q_massive
	switch (hc) {
	default:// We return zero for all other helicity combinations
		return 0;
	}
}


/*
 *
 *
 *
 * Two scalar two lepton two quark amplitudes (DOES NOT WORK
 *
 *
 *
 */



template <class T> complex<T> (*A2s2q2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch(hc){

	default:// We return zero for all other helicity combinations
		return &ZeroF;
	}
}


/*
 *
 *
 * The 2 massive quarks and two massless quarks and a photon amplitudes (with the photon only coupling to the massless quarks)
 *
 *     q\  /y
 *       \/       /Q
 *        \___g__/
 *        /      \
 *     qb/        \Qb
 *
 *
 *
 */

#define SIGN1 T(-1)
#define SIGN2 1

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_1_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> qLt=ep.ref()->Lt();

	return ep.spb(i2,i0)*((ep.spa(i1,i0)*(ep.p(i0)->Lt()*f4.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f4.Lt()))*(f5.L()*ep.p(i1)->L())
				-SIGN1*mass*(ep.spa(i1,i0)*(ep.p(i0)->Lt()*qLt)+ep.spa(i1,i2)*ep.spb(i2,5))*(qL*ep.p(i1)->L())/((f5.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f4.L())))
				/(complex<T>(0,-4)*ep.spa(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_2_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> qLt=ep.ref()->Lt();

	return ep.spb(i2,i0)*((ep.spa(i1,i0)*(ep.p(i0)->Lt()*f5.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f5.Lt()))*(f4.L()*ep.p(i1)->L())
						-SIGN1*mass*(ep.spa(i1,i0)*(ep.p(i0)->Lt()*qLt)+ep.spa(i1,i2)*ep.spb(i2,5))*(qL*ep.p(i1)->L())/((f5.L()*ep.ref()->L())*(ep.ref()->Lt()*f4.Lt())))
						/(complex<T>(0,4)*ep.spa(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_3m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*ep.spb(i2,i0)*(qL*ep.p(i1)->L())
			*((ep.spa(i1,i0)*(ep.p(i0)->Lt()*f4.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f4.Lt()))/(f5.L()*ep.ref()->L())
			+(ep.spa(i1,i0)*(ep.p(i0)->Lt()*f5.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f5.Lt()))/(ep.ref()->L()*f4.L()))
			/(complex<T>(0,-4)*ep.spa(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_3p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*ep.spb(i2,i0)*(qL*ep.p(i1)->L())
			*((ep.spa(i1,i0)*(ep.p(i0)->Lt()*f4.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f4.Lt()))/(f5.L()*ep.ref()->L())
			+(ep.spa(i1,i0)*(ep.p(i0)->Lt()*f5.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f5.Lt()))/(ep.ref()->L()*f4.L()))
			/(complex<T>(0,4)*ep.spa(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}


template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_4m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*ep.spb(i2,i0)*(qL*ep.p(i1)->L())
			*((ep.spa(i1,i0)*(ep.p(i0)->Lt()*f5.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f5.Lt()))/(f4.L()*ep.ref()->L())
			+(ep.spa(i1,i0)*(ep.p(i0)->Lt()*f4.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f4.Lt()))/(ep.ref()->L()*f5.L()))
			/(complex<T>(0,-4)*ep.spa(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_4p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*ep.spb(i2,i0)*(qL*ep.p(i1)->L())
			*((ep.spa(i1,i0)*(ep.p(i0)->Lt()*f5.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f5.Lt()))/(f4.L()*ep.ref()->L())
			+(ep.spa(i1,i0)*(ep.p(i0)->Lt()*f4.Lt())+ep.spa(i1,i2)*(ep.p(i2)->Lt()*f4.Lt()))/(ep.ref()->L()*f5.L()))
			/(complex<T>(0,4)*ep.spa(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_5_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> qLt=ep.ref()->Lt();

	return ((ep.p(i0)->Lt()*f4.Lt())*(f5.L()*ep.p(i2)->L())
			+ep.spa(i2,i1)*(ep.p(i2)->L()*f5.L())*(f4.Lt()*ep.p(i1)->Lt())/ep.spa(i2,i0)
			+SIGN1*mass*((ep.p(i0)->Lt()*qLt)*((qL*ep.p(i2)->L()))
				+ep.spa(i2,i1)*(ep.p(i2)->L()*qL)*ep.spb(5,i1)/ep.spa(i2,i0))
				/((ep.ref()->L()*f4.L())*(ep.ref()->Lt()*f5.Lt()))
			)/(complex<T>(0,2)*ep.spa(i1,i0)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();
	const lambda<T> qLt=ep.ref()->Lt();

	return ((ep.p(i0)->Lt()*f5.Lt())*(f4.L()*ep.p(i2)->L())
			+ep.spa(i2,i1)*(ep.p(i2)->L()*f4.L())*(f5.Lt()*ep.p(i1)->Lt())/ep.spa(i2,i0)
			+SIGN1*mass*((ep.p(i0)->Lt()*qLt)*((qL*ep.p(i2)->L()))
				+ep.spa(i2,i1)*(ep.p(i2)->L()*qL)*ep.spb(5,i1)/ep.spa(i2,i0))
				/((ep.ref()->L()*f5.L())*(ep.ref()->Lt()*f4.Lt()))
			)/(complex<T>(0,-2)*ep.spa(i1,i0)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_7m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*((ep.p(i0)->Lt()*f4.Lt())*((qL*ep.p(i2)->L()))/(f5.L()*ep.ref()->L())
			+ep.spa(i2,i1)*(ep.p(i1)->Lt()*f4.Lt())*((qL*ep.p(i2)->L()))/(ep.spa(i2,i0)*(f5.L()*ep.ref()->L()))
			+(ep.p(i0)->Lt()*f5.Lt())*((qL*ep.p(i2)->L()))/(ep.ref()->L()*f4.L())
			+ep.spa(i2,i1)*(ep.p(i1)->Lt()*f5.Lt())*((qL*ep.p(i2)->L()))/(ep.spa(i2,i0)*(ep.ref()->L()*f4.L()))
			)/(complex<T>(0,2)*ep.spa(i1,i0)*ep.sp(i3,i4));
}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_7p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*((ep.p(i0)->Lt()*f4.Lt())*((qL*ep.p(i2)->L()))/(f5.L()*ep.ref()->L())
			+ep.spa(i2,i1)*(ep.p(i1)->Lt()*f4.Lt())*((qL*ep.p(i2)->L()))/(ep.spa(i2,i0)*(f5.L()*ep.ref()->L()))
			+(ep.p(i0)->Lt()*f5.Lt())*((qL*ep.p(i2)->L()))/(ep.ref()->L()*f4.L())
			+ep.spa(i2,i1)*(ep.p(i1)->Lt()*f5.Lt())*((qL*ep.p(i2)->L()))/(ep.spa(i2,i0)*(ep.ref()->L()*f4.L()))
			)/(complex<T>(0,-2)*ep.spa(i1,i0)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_8m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*((ep.p(i0)->Lt()*f5.Lt())*((qL*ep.p(i2)->L()))/(f4.L()*ep.ref()->L())
			+ep.spa(i2,i1)*(ep.p(i1)->Lt()*f5.Lt())*((qL*ep.p(i2)->L()))/(ep.spa(i2,i0)*(f4.L()*ep.ref()->L()))
			+(ep.p(i0)->Lt()*f4.Lt())*((qL*ep.p(i2)->L()))/(ep.ref()->L()*f5.L())
			+ep.spa(i2,i1)*(ep.p(i1)->Lt()*f4.Lt())*((qL*ep.p(i2)->L()))/(ep.spa(i2,i0)*(ep.ref()->L()*f5.L()))
			)/(complex<T>(0,2)*ep.spa(i1,i0)*ep.sp(i3,i4));
}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_8p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*((ep.p(i0)->Lt()*f5.Lt())*((qL*ep.p(i2)->L()))/(f4.L()*ep.ref()->L())
			+ep.spa(i2,i1)*(ep.p(i1)->Lt()*f5.Lt())*((qL*ep.p(i2)->L()))/(ep.spa(i2,i0)*(f4.L()*ep.ref()->L()))
			+(ep.p(i0)->Lt()*f4.Lt())*((qL*ep.p(i2)->L()))/(ep.ref()->L()*f5.L())
			+ep.spa(i2,i1)*(ep.p(i1)->Lt()*f4.Lt())*((qL*ep.p(i2)->L()))/(ep.spa(i2,i0)*(ep.ref()->L()*f5.L()))
			)/(complex<T>(0,-2)*ep.spa(i1,i0)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_9_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return ep.spa(i2,i0)*((ep.spb(i1,i0)*(ep.p(i0)->L()*f5.L())+ep.spb(i1,i2)*(ep.p(i2)->L()*f5.L()))*(f4.Lt()*ep.p(i1)->Lt())
						-mass*(ep.spb(i1,i0)*(ep.p(i0)->L()*qL)+ep.spb(i1,i2)*(ep.p(i2)->L()*qL))*ep.spb(5,i1)/((f5.Lt()*ep.ref()->Lt())*(ep.ref()->L()*f4.L())))
						/(complex<T>(0,-4)*ep.spb(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_10_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return ep.spa(i2,i0)*((ep.spb(i1,i0)*(ep.p(i0)->L()*f4.L())+ep.spb(i1,i2)*(ep.p(i2)->L()*f4.L()))*(f5.Lt()*ep.p(i1)->Lt())
				-SIGN1*mass*(ep.spb(i1,i0)*(ep.p(i0)->L()*qL)+ep.spb(i1,i2)*(ep.p(i2)->L()*qL))*ep.spb(5,i1)/((f5.L()*ep.ref()->L())*(ep.ref()->Lt()*f4.Lt())))
				/(complex<T>(0,4)*ep.spb(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_11m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*ep.spa(i2,i0)*(ep.spb(i2,i0)*(ep.p(i0)->L()*qL)+ep.spb(i1,i2)*(ep.p(i2)->L()*qL))
		*((f4.Lt()*ep.p(i1)->Lt())/(f5.L()*ep.ref()->L())+(f5.Lt()*ep.p(i1)->Lt())/(ep.ref()->L()*f4.L()))
		/(complex<T>(0,-4)*ep.spa(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_11p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*ep.spa(i2,i0)*(ep.spb(i2,i0)*(ep.p(i0)->L()*qL)+ep.spb(i1,i2)*(ep.p(i2)->L()*qL))
		*((f4.Lt()*ep.p(i1)->Lt())/(f5.L()*ep.ref()->L())+(f5.Lt()*ep.p(i1)->Lt())/(ep.ref()->L()*f4.L()))
		/(complex<T>(0,4)*ep.spa(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_12m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	return eval_param<T>::mass(masses.p(i3))*ep.spa(i2,i0)*ep.spb(5,i1)
			*((ep.spb(i1,i0)*(ep.p(i0)->L()*f4.L())+ep.spb(i1,i2)*(ep.p(i2)->L()*f4.L()))/(f5.Lt()*ep.ref()->Lt())
			+(ep.spb(i1,i0)*(ep.p(i0)->L()*f5.L())+ep.spb(i1,i2)*(ep.p(i2)->L()*f5.L()))/(ep.ref()->Lt()*f4.Lt()))
			/(complex<T>(0,-4)*ep.spb(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_12p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	return eval_param<T>::mass(masses.p(i3))*ep.spa(i2,i0)*ep.spb(5,i1)
			*((ep.spb(i1,i0)*(ep.p(i0)->L()*f4.L())+ep.spb(i1,i2)*(ep.p(i2)->L()*f4.L()))/(f5.Lt()*ep.ref()->Lt())
			+(ep.spb(i1,i0)*(ep.p(i0)->L()*f5.L())+ep.spb(i1,i2)*(ep.p(i2)->L()*f5.L()))/(ep.ref()->Lt()*f4.Lt()))
			/(complex<T>(0,4)*ep.spb(i1,i0)*ep.sp(i0,i2)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_13_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return ((ep.p(i0)->L()*f5.L())*(f4.Lt()*ep.p(i2)->Lt())
			+ep.spb(i2,i1)*(ep.p(i2)->Lt()*f4.Lt())*(f5.L()*ep.p(i1)->L())/ep.spb(i2,i0)
			+SIGN1*mass*((ep.p(i0)->L()*qL)*ep.spb(5,i2)
				+ep.spb(i2,i1)*ep.spb(i2,5)*(qL*ep.p(i1)->L())/ep.spb(i2,i0))
				/((ep.ref()->Lt()*f5.Lt())*(ep.ref()->L()*f4.L()))
			)/(complex<T>(0,2)*ep.spb(i1,i0)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_14_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return ((ep.p(i0)->L()*f4.L())*(f5.Lt()*ep.p(i2)->Lt())
			+ep.spb(i2,i1)*(ep.p(i2)->Lt()*f5.Lt())*(f4.L()*ep.p(i1)->L())/ep.spb(i2,i0)
			+SIGN1*mass*((ep.p(i0)->L()*qL)*ep.spb(5,i2)
				+ep.spb(i2,i1)*ep.spb(i2,5)*(qL*ep.p(i1)->L())/ep.spb(i2,i0))
				/((ep.ref()->Lt()*f4.Lt())*(ep.ref()->L()*f5.L()))
			)/(complex<T>(0,-2)*ep.spb(i1,i0)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_15m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*((ep.p(i0)->L()*qL)*(f4.Lt()*ep.p(i2)->Lt())/(f5.L()*ep.ref()->L())
			+ep.spb(i2,i1)*(ep.p(i1)->L()*qL)*(f4.Lt()*ep.p(i2)->Lt())/(ep.spb(i2,i0)*(f5.L()*ep.ref()->L()))
			+((ep.p(i0)->L()*qL)*(f5.Lt()*ep.p(i2)->Lt())/(ep.ref()->L()*f4.L())
				+ep.spb(i2,i1)*(ep.p(i1)->L()*qL)*(f5.Lt()*ep.p(i2)->Lt())/(ep.spb(i2,i0)*(ep.ref()->L()*f4.L())))
				/((ep.ref()->Lt()*f5.Lt())*(ep.ref()->L()*f4.L()))
			)/(complex<T>(0,-2)*ep.spb(i1,i0)*ep.sp(i3,i4));
}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_15p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	const lambda<T> qL=ep.ref()->L();

	return eval_param<T>::mass(masses.p(i3))*((ep.p(i0)->L()*qL)*(f4.Lt()*ep.p(i2)->Lt())/(f5.L()*ep.ref()->L())
			+ep.spb(i2,i1)*(ep.p(i1)->L()*qL)*(f4.Lt()*ep.p(i2)->Lt())/(ep.spb(i2,i0)*(f5.L()*ep.ref()->L()))
			+((ep.p(i0)->L()*qL)*(f5.Lt()*ep.p(i2)->Lt())/(ep.ref()->L()*f4.L())
				+ep.spb(i2,i1)*(ep.p(i1)->L()*qL)*(f5.Lt()*ep.p(i2)->Lt())/(ep.spb(i2,i0)*(ep.ref()->L()*f4.L())))
				/((ep.ref()->Lt()*f5.Lt())*(ep.ref()->L()*f4.L()))
			)/(complex<T>(0,2)*ep.spb(i1,i0)*ep.sp(i3,i4));
}

template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_16m_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	return eval_param<T>::mass(masses.p(i3))*((ep.p(i0)->L()*f4.L())*ep.spb(5,i2)/(f5.Lt()*ep.ref()->Lt())
			+ep.spb(i2,i1)*(ep.p(i1)->L()*f4.L())*ep.spb(5,i2)/(ep.spb(i2,i0)*(f5.Lt()*ep.ref()->Lt()))
			+((ep.p(i0)->L()*f5.L())*ep.spb(5,i2)/(ep.ref()->Lt()*f4.Lt())
				+ep.spb(i2,i1)*(ep.p(i1)->L()*f5.L())*ep.spb(5,i2)/(ep.spa(i2,i0)*(ep.ref()->Lt()*f4.Lt())))
				/((ep.ref()->Lt()*f5.Lt())*(ep.ref()->L()*f4.L()))
			)/(complex<T>(0,-2)*ep.spb(i1,i0)*ep.sp(i3,i4));
}
template <int i0, int i1, int i2, int i3, int i4, class T> complex<T>  A2QM2q1y_16p_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	// Construct the flatted momenta for the massive 0 and 2 legs
	const complex<T> mass=eval_param<T>::mass2(masses.p(i3));
	// We use the negative helicity leg as the reference vector
	Cmom<T> f4=ep.mom(i3)-(mass/(T(2)*(ep.p(i3)->P()*ep.ref()->P())))*ep.ref()->P();
	Cmom<T> f5=ep.mom(i4)-(mass/(T(2)*(ep.p(i4)->P()*ep.ref()->P())))*ep.ref()->P();
	return eval_param<T>::mass(masses.p(i3))*((ep.p(i0)->L()*f4.L())*ep.spb(5,i2)/(f5.Lt()*ep.ref()->Lt())
			+ep.spb(i2,i1)*(ep.p(i1)->L()*f4.L())*ep.spb(5,i2)/(ep.spb(i2,i0)*(f5.Lt()*ep.ref()->Lt()))
			+((ep.p(i0)->L()*f5.L())*ep.spb(5,i2)/(ep.ref()->Lt()*f4.Lt())
				+ep.spb(i2,i1)*(ep.p(i1)->L()*f5.L())*ep.spb(5,i2)/(ep.spa(i2,i0)*(ep.ref()->Lt()*f4.Lt())))
				/((ep.ref()->Lt()*f5.Lt())*(ep.ref()->L()*f4.L()))
			)/(complex<T>(0,2)*ep.spb(i1,i0)*ep.sp(i3,i4));
}



template <class T> complex<T> (*A2QM2q1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch(hc){
	case 0x12385: //y+qb-q+Qb+Q-
	case 0x12367: //y+qb-q+Q+Qb-
		return &A2QM2q1y_1_eval<0,1,2,3,4>;
		break;
	case 0x51238:
	case 0x71236:
		return &A2QM2q1y_1_eval<1,2,3,4,0>;
		break;
	case 0x85123:
	case 0x67123:
		return &A2QM2q1y_1_eval<2,3,4,0,1>;
		break;
	case 0x38512:
	case 0x36712:
		return &A2QM2q1y_1_eval<3,4,0,1,2>;
		break;
	case 0x23851:
	case 0x23671:
		return &A2QM2q1y_1_eval<4,0,1,2,3>;
		break;

	case 0x12376: //y+qb-q+Qb-Q+
	case 0x12358: //y+qb-q+Qb-Q+
		return &A2QM2q1y_2_eval<0,1,2,3,4>;
		break;
	case 0x61237:
	case 0x81235:
		return &A2QM2q1y_2_eval<1,2,3,4,0>;
		break;
	case 0x76123:
	case 0x58123:
		return &A2QM2q1y_2_eval<2,3,4,0,1>;
		break;
	case 0x37612:
	case 0x35812:
		return &A2QM2q1y_2_eval<3,4,0,1,2>;
		break;
	case 0x23761:
	case 0x23581:
		return &A2QM2q1y_2_eval<4,0,1,2,3>;
		break;

	case 0x12386: //y+qb-q+Qb+Q+
		return &A2QM2q1y_3m_eval<0,1,2,3,4>;
		break;
	case 0x12368: //y+qb-q+Q+Qb+
		return &A2QM2q1y_3p_eval<0,1,2,3,4>;
		break;
	case 0x61238:
		return &A2QM2q1y_3m_eval<1,2,3,4,0>;
		break;
	case 0x81236:
		return &A2QM2q1y_3p_eval<1,2,3,4,0>;
		break;
	case 0x86123:
		return &A2QM2q1y_3m_eval<2,3,4,0,1>;
		break;
	case 0x68123:
		return &A2QM2q1y_3p_eval<2,3,4,0,1>;
		break;
	case 0x38612:
		return &A2QM2q1y_3m_eval<3,4,0,1,2>;
		break;
	case 0x36812:
		return &A2QM2q1y_3p_eval<3,4,0,1,2>;
		break;
	case 0x23861:
		return &A2QM2q1y_3m_eval<4,0,1,2,3>;
		break;
	case 0x23681:
		return &A2QM2q1y_3p_eval<4,0,1,2,3>;
		break;

	case 0x12375: //y+qb-q+Qb-Q-
		return &A2QM2q1y_4m_eval<0,1,2,3,4>;
		break;
	case 0x12357: //y+qb-q+Q-Qb-
		return &A2QM2q1y_4p_eval<0,1,2,3,4>;
		break;
	case 0x51237:
		return &A2QM2q1y_4m_eval<1,2,3,4,0>;
		break;
	case 0x71235:
		return &A2QM2q1y_4p_eval<1,2,3,4,0>;
		break;
	case 0x75123:
		return &A2QM2q1y_4m_eval<2,3,4,0,1>;
		break;
	case 0x57123:
		return &A2QM2q1y_4p_eval<2,3,4,0,1>;
		break;
	case 0x37512:
		return &A2QM2q1y_4m_eval<3,4,0,1,2>;
		break;
	case 0x35712:
		return &A2QM2q1y_4p_eval<3,4,0,1,2>;
		break;
	case 0x23751:
		return &A2QM2q1y_4m_eval<4,0,1,2,3>;
		break;
	case 0x23571:
		return &A2QM2q1y_4p_eval<4,0,1,2,3>;
		break;

	case 0x13285: //y+qb+q-Qb+Q-
	case 0x13267: //y+qb+q-Q+Qb-
		return &A2QM2q1y_5_eval<0,1,2,3,4>;
		break;
	case 0x51328:
	case 0x71326:
		return &A2QM2q1y_5_eval<1,2,3,4,0>;
		break;
	case 0x85132:
	case 0x67132:
		return &A2QM2q1y_5_eval<2,3,4,0,1>;
		break;
	case 0x28513:
	case 0x26713:
		return &A2QM2q1y_5_eval<3,4,0,1,2>;
		break;
	case 0x32851:
	case 0x32671:
		return &A2QM2q1y_5_eval<4,0,1,2,3>;
		break;

	case 0x13276: //y+qb+q-Qb-Q+
	case 0x13258: //y+qb+q-Qb-Q+
		return &A2QM2q1y_6_eval<0,1,2,3,4>;
		break;
	case 0x61327:
	case 0x81325:
		return &A2QM2q1y_6_eval<1,2,3,4,0>;
		break;
	case 0x76132:
	case 0x58132:
		return &A2QM2q1y_6_eval<2,3,4,0,1>;
		break;
	case 0x27613:
	case 0x25813:
		return &A2QM2q1y_6_eval<3,4,0,1,2>;
		break;
	case 0x32761:
	case 0x32581:
		return &A2QM2q1y_6_eval<4,0,1,2,3>;
		break;

	case 0x13286: //y+qb+q-Qb+Q+
		return &A2QM2q1y_7m_eval<0,1,2,3,4>;
		break;
	case 0x13268: //y+qb+q-Q+Qb+
		return &A2QM2q1y_7p_eval<0,1,2,3,4>;
		break;
	case 0x61328:
		return &A2QM2q1y_7m_eval<1,2,3,4,0>;
		break;
	case 0x81326:
		return &A2QM2q1y_7p_eval<1,2,3,4,0>;
		break;
	case 0x86132:
		return &A2QM2q1y_7m_eval<2,3,4,0,1>;
		break;
	case 0x68132:
		return &A2QM2q1y_7p_eval<2,3,4,0,1>;
		break;
	case 0x28613:
		return &A2QM2q1y_7m_eval<2,3,4,0,1>;
		break;
	case 0x26813:
		return &A2QM2q1y_7p_eval<3,4,0,1,2>;
		break;
	case 0x32861:
		return &A2QM2q1y_7m_eval<3,4,0,1,2>;
		break;
	case 0x32681:
		return &A2QM2q1y_7p_eval<4,0,1,2,3>;
		break;

	case 0x13275: //y+qb+q-Qb-Q-
		return &A2QM2q1y_8m_eval<0,1,2,3,4>;
		break;
	case 0x13257: //y+qb+q-Q-Qb-
		return &A2QM2q1y_8p_eval<0,1,2,3,4>;
		break;
	case 0x51327:
		return &A2QM2q1y_8m_eval<1,2,3,4,0>;
		break;
	case 0x71325:
		return &A2QM2q1y_8p_eval<1,2,3,4,0>;
		break;
	case 0x75132:
		return &A2QM2q1y_8m_eval<2,3,4,0,1>;
		break;
	case 0x57132:
		return &A2QM2q1y_8p_eval<2,3,4,0,1>;
		break;
	case 0x27513:
		return &A2QM2q1y_8m_eval<3,4,0,1,2>;
		break;
	case 0x25713:
		return &A2QM2q1y_8p_eval<3,4,0,1,2>;
		break;
	case 0x32751:
		return &A2QM2q1y_8m_eval<4,0,1,2,3>;
		break;
	case 0x32571:
		return &A2QM2q1y_8p_eval<4,0,1,2,3>;
		break;

	case 0x03285: //y-qb+q-Qb+Q-
	case 0x03267: //y-qb+q-Q+Qb-
		return &A2QM2q1y_9_eval<0,1,2,3,4>;
		break;
	case 0x50328:
	case 0x70326:
		return &A2QM2q1y_9_eval<1,2,3,4,0>;
		break;
	case 0x85032:
	case 0x67032:
		return &A2QM2q1y_9_eval<2,3,4,0,1>;
		break;
	case 0x28503:
	case 0x26703:
		return &A2QM2q1y_9_eval<3,4,0,1,2>;
		break;
	case 0x32850:
	case 0x32670:
		return &A2QM2q1y_9_eval<4,0,1,2,3>;
		break;

	case 0x03276: //y-qb+q-Qb-Q+
	case 0x03258: //y-qb+q-Qb-Q+
		return &A2QM2q1y_10_eval<0,1,2,3,4>;
		break;
	case 0x60327:
	case 0x80325:
		return &A2QM2q1y_10_eval<1,2,3,4,0>;
		break;
	case 0x76032:
	case 0x58032:
		return &A2QM2q1y_10_eval<2,3,4,0,1>;
		break;
	case 0x27603:
	case 0x25803:
		return &A2QM2q1y_10_eval<3,4,0,1,2>;
		break;
	case 0x32760:
	case 0x32580:
		return &A2QM2q1y_10_eval<4,0,1,2,3>;
		break;

	case 0x03286: //y-qb+q-Qb+Q+
		return &A2QM2q1y_11m_eval<0,1,2,3,4>;
		break;
	case 0x03268: //y-qb+q-Q+Qb+
		return &A2QM2q1y_11p_eval<0,1,2,3,4>;
		break;
	case 0x60328:
		return &A2QM2q1y_11m_eval<1,2,3,4,0>;
		break;
	case 0x80326:
		return &A2QM2q1y_11p_eval<1,2,3,4,0>;
		break;
	case 0x86032:
		return &A2QM2q1y_11m_eval<2,3,4,0,1>;
		break;
	case 0x68032:
		return &A2QM2q1y_11p_eval<2,3,4,0,1>;
		break;
	case 0x28603:
		return &A2QM2q1y_11m_eval<3,4,0,1,2>;
		break;
	case 0x26803:
		return &A2QM2q1y_11p_eval<3,4,0,1,2>;
		break;
	case 0x32860:
		return &A2QM2q1y_11m_eval<4,0,1,2,3>;
		break;
	case 0x32680:
		return &A2QM2q1y_11p_eval<4,0,1,2,3>;
		break;

	case 0x03275: //y-qb+q-Qb-Q-
		return &A2QM2q1y_12m_eval<0,1,2,3,4>;
		break;
	case 0x03257: //y-qb+q-Q-Qb-
		return &A2QM2q1y_12p_eval<0,1,2,3,4>;
		break;
	case 0x50327:
		return &A2QM2q1y_12m_eval<1,2,3,4,0>;
		break;
	case 0x70325:
		return &A2QM2q1y_12p_eval<1,2,3,4,0>;
		break;
	case 0x75032:
		return &A2QM2q1y_12m_eval<2,3,4,0,1>;
		break;
	case 0x57032:
		return &A2QM2q1y_12p_eval<2,3,4,0,1>;
		break;
	case 0x27503:
		return &A2QM2q1y_12m_eval<3,4,0,1,2>;
		break;
	case 0x25703:
		return &A2QM2q1y_12p_eval<3,4,0,1,2>;
		break;
	case 0x32750:
		return &A2QM2q1y_12m_eval<4,0,1,2,3>;
		break;
	case 0x32570:
		return &A2QM2q1y_12p_eval<4,0,1,2,3>;
		break;

	case 0x02385: //y-qb-q+Qb+Q-
	case 0x02367: //y-qb-q+Q+Qb-
		return &A2QM2q1y_13_eval<0,1,2,3,4>;
		break;
	case 0x50238:
	case 0x70236:
		return &A2QM2q1y_13_eval<1,2,3,4,0>;
		break;
	case 0x85023:
	case 0x67023:
		return &A2QM2q1y_13_eval<2,3,4,0,1>;
		break;
	case 0x38502:
	case 0x36702:
		return &A2QM2q1y_13_eval<3,4,0,1,2>;
		break;
	case 0x23850:
	case 0x23670:
		return &A2QM2q1y_13_eval<4,0,1,2,3>;
		break;

	case 0x02376: //y-qb-q+Qb-Q+
	case 0x02358: //y-qb-q+Qb-Q+
		return &A2QM2q1y_14_eval<0,1,2,3,4>;
		break;
	case 0x60237:
	case 0x80235:
		return &A2QM2q1y_14_eval<1,2,3,4,0>;
		break;
	case 0x76023:
	case 0x58023:
		return &A2QM2q1y_14_eval<2,3,4,0,1>;
		break;
	case 0x37602:
	case 0x35802:
		return &A2QM2q1y_14_eval<3,4,0,1,2>;
		break;
	case 0x23760:
	case 0x23580:
		return &A2QM2q1y_14_eval<4,0,1,2,3>;
		break;

	case 0x02386: //y-qb-q+Qb+Q+
		return &A2QM2q1y_15m_eval<0,1,2,3,4>;
		break;
	case 0x02368: //y-qb-q+Q+Qb+
		return &A2QM2q1y_15p_eval<0,1,2,3,4>;
		break;
	case 0x60238:
		return &A2QM2q1y_15m_eval<1,2,3,4,0>;
		break;
	case 0x80236:
		return &A2QM2q1y_15p_eval<1,2,3,4,0>;
		break;
	case 0x86023:
		return &A2QM2q1y_15m_eval<2,3,4,0,1>;
		break;
	case 0x68023:
		return &A2QM2q1y_15p_eval<2,3,4,0,1>;
		break;
	case 0x38602:
		return &A2QM2q1y_15m_eval<3,4,0,1,2>;
		break;
	case 0x36802:
		return &A2QM2q1y_15p_eval<3,4,0,1,2>;
		break;
	case 0x23860:
		return &A2QM2q1y_15m_eval<4,0,1,2,3>;
		break;
	case 0x23680:
		return &A2QM2q1y_15p_eval<4,0,1,2,3>;
		break;

	case 0x02375: //y-qb-q+Qb-Q-
		return &A2QM2q1y_16m_eval<0,1,2,3,4>;
		break;
	case 0x02357: //y-qb-q+Q-Qb-
		return &A2QM2q1y_16p_eval<0,1,2,3,4>;
		break;
	case 0x50237:
		return &A2QM2q1y_16m_eval<1,2,3,4,0>;
		break;
	case 0x70235:
		return &A2QM2q1y_16p_eval<1,2,3,4,0>;
		break;
	case 0x75023:
		return &A2QM2q1y_16m_eval<2,3,4,0,1>;
		break;
	case 0x57023:
		return &A2QM2q1y_16p_eval<2,3,4,0,1>;
		break;
	case 0x37502:
		return &A2QM2q1y_16m_eval<3,4,0,1,2>;
		break;
	case 0x35702:
		return &A2QM2q1y_16p_eval<3,4,0,1,2>;
		break;
	case 0x23750:
		return &A2QM2q1y_16m_eval<4,0,1,2,3>;
		break;
	case 0x23570:
		return &A2QM2q1y_16p_eval<4,0,1,2,3>;
		break;

		//*********************************************************

	case 0x18532: //y+Qb+Q-qb+q-
	case 0x16732: //y+Q+Qb-qb+q-
		return &A2QM2q1y_5_eval<0,3,4,1,2>;
		break;
	case 0x21853:
	case 0x21673:
		return &A2QM2q1y_5_eval<1,4,0,2,3>;
		break;
	case 0x32185:
	case 0x32167:
		return &A2QM2q1y_5_eval<2,0,1,3,4>;
		break;
	case 0x53218:
	case 0x73216:
		return &A2QM2q1y_5_eval<3,1,2,4,0>;
		break;
	case 0x85321:
	case 0x67321:
		return &A2QM2q1y_5_eval<4,2,3,0,1>;
		break;

	case 0x17632: //y+Qb-Q+qb+q-
	case 0x15832: //y+Qb-Q+qb+q-
		return &A2QM2q1y_6_eval<0,3,4,1,2>;
		break;
	case 0x21763:
	case 0x21583:
		return &A2QM2q1y_6_eval<1,4,0,2,3>;
		break;
	case 0x32176:
	case 0x32158:
		return &A2QM2q1y_6_eval<2,0,1,3,4>;
		break;
	case 0x63217:
	case 0x83215:
		return &A2QM2q1y_6_eval<3,1,2,4,0>;
		break;
	case 0x76321:
	case 0x58321:
		return &A2QM2q1y_6_eval<4,2,3,0,1>;
		break;

	case 0x18632: //y+Qb+Q+qb+q-
		return &A2QM2q1y_7m_eval<0,3,4,1,2>;
		break;
	case 0x16832: //y+Q+Qb+qb+q-
		return &A2QM2q1y_7p_eval<0,3,4,1,2>;
		break;
	case 0x21863:
		return &A2QM2q1y_7m_eval<1,4,0,2,3>;
		break;
	case 0x21683:
		return &A2QM2q1y_7p_eval<1,4,0,2,3>;
		break;
	case 0x32186:
		return &A2QM2q1y_7m_eval<2,0,1,3,4>;
		break;
	case 0x32168:
		return &A2QM2q1y_7p_eval<2,0,1,3,4>;
		break;
	case 0x63218:
		return &A2QM2q1y_7m_eval<3,1,2,4,0>;
		break;
	case 0x83216:
		return &A2QM2q1y_7p_eval<3,1,2,4,0>;
		break;
	case 0x86321:
		return &A2QM2q1y_7m_eval<4,2,3,0,1>;
		break;
	case 0x68321:
		return &A2QM2q1y_7p_eval<4,2,3,0,1>;
		break;

	case 0x17532: //y+Qb-Q-qb+q-
		return &A2QM2q1y_8m_eval<0,3,4,1,2>;
		break;
	case 0x15732: //y+Q-Qb-qb+q-
		return &A2QM2q1y_8p_eval<0,3,4,1,2>;
		break;
	case 0x21753:
		return &A2QM2q1y_8m_eval<1,4,0,2,3>;
		break;
	case 0x21573:
		return &A2QM2q1y_8p_eval<1,4,0,2,3>;
		break;
	case 0x32175:
		return &A2QM2q1y_8m_eval<2,0,1,3,4>;
		break;
	case 0x32157:
		return &A2QM2q1y_8p_eval<2,0,1,3,4>;
		break;
	case 0x53217:
		return &A2QM2q1y_8m_eval<3,1,2,4,0>;
		break;
	case 0x73215:
		return &A2QM2q1y_8p_eval<3,1,2,4,0>;
		break;
	case 0x75321:
		return &A2QM2q1y_8m_eval<4,2,3,0,1>;
		break;
	case 0x57321:
		return &A2QM2q1y_8p_eval<4,2,3,0,1>;
		break;

	case 0x18523: //y+Qb+Q-qb-q+
	case 0x16723: //y+Q+Qb-qb-q+
		return &A2QM2q1y_1_eval<0,3,4,1,2>;
		break;
	case 0x31852:
	case 0x31672:
		return &A2QM2q1y_1_eval<1,4,0,2,3>;
		break;
	case 0x23185:
	case 0x23167:
		return &A2QM2q1y_1_eval<2,0,1,3,4>;
		break;
	case 0x52318:
	case 0x72316:
		return &A2QM2q1y_1_eval<3,1,2,4,0>;
		break;
	case 0x85231:
	case 0x67231:
		return &A2QM2q1y_1_eval<4,2,3,0,1>;
		break;

	case 0x17623: //y+Qb-Q+qb-q+
	case 0x15823: //y+Qb-Q+qb-q+
		return &A2QM2q1y_2_eval<0,3,4,1,2>;
		break;
	case 0x31762:
	case 0x31582:
		return &A2QM2q1y_2_eval<1,4,0,2,3>;
		break;
	case 0x23176:
	case 0x23158:
		return &A2QM2q1y_2_eval<2,0,1,3,4>;
		break;
	case 0x62317:
	case 0x82315:
		return &A2QM2q1y_2_eval<3,1,2,4,0>;
		break;
	case 0x76231:
	case 0x58231:
		return &A2QM2q1y_2_eval<4,2,3,0,1>;
		break;

	case 0x18623: //y+Qb+Q+qb-q+
		return &A2QM2q1y_3m_eval<0,3,4,1,2>;
		break;
	case 0x16823: //y+Q+Qb+qb-q+
		return &A2QM2q1y_3p_eval<0,3,4,1,2>;
		break;
	case 0x31862:
		return &A2QM2q1y_3m_eval<1,4,0,2,3>;
		break;
	case 0x31682:
		return &A2QM2q1y_3p_eval<1,4,0,2,3>;
		break;
	case 0x23186:
		return &A2QM2q1y_3m_eval<2,0,1,3,4>;
		break;
	case 0x23168:
		return &A2QM2q1y_3p_eval<2,0,1,3,4>;
		break;
	case 0x62318:
		return &A2QM2q1y_3m_eval<3,1,2,4,0>;
		break;
	case 0x82316:
		return &A2QM2q1y_3p_eval<3,1,2,4,0>;
		break;
	case 0x86231:
		return &A2QM2q1y_3m_eval<4,2,3,0,1>;
		break;
	case 0x68231:
		return &A2QM2q1y_3p_eval<4,2,3,0,1>;
		break;

	case 0x17523: //y+Qb-Q-qb-q+
		return &A2QM2q1y_4m_eval<0,3,4,1,2>;
		break;
	case 0x15723: //y+Q-Qb-qb-q+
		return &A2QM2q1y_4p_eval<0,3,4,1,2>;
		break;
	case 0x31752:
		return &A2QM2q1y_4m_eval<1,4,0,2,3>;
		break;
	case 0x31572:
		return &A2QM2q1y_4p_eval<1,4,0,2,3>;
		break;
	case 0x23175:
		return &A2QM2q1y_4m_eval<2,0,1,3,4>;
		break;
	case 0x23157:
		return &A2QM2q1y_4p_eval<2,0,1,3,4>;
		break;
	case 0x52317:
		return &A2QM2q1y_4m_eval<3,1,2,4,0>;
		break;
	case 0x72315:
		return &A2QM2q1y_4p_eval<3,1,2,4,0>;
		break;
	case 0x75231:
		return &A2QM2q1y_4m_eval<4,2,3,0,1>;
		break;
	case 0x57231:
		return &A2QM2q1y_4p_eval<4,2,3,0,1>;
		break;

	case 0x08532: //y-Qb+Q-qb+q-
	case 0x06732: //y-Q+Qb-qb+q-
		return &A2QM2q1y_9_eval<0,3,4,1,2>;
		break;
	case 0x20853:
	case 0x20673:
		return &A2QM2q1y_9_eval<1,4,0,2,3>;
		break;
	case 0x32085:
	case 0x32067:
		return &A2QM2q1y_9_eval<2,0,1,3,4>;
		break;
	case 0x53208:
	case 0x73206:
		return &A2QM2q1y_9_eval<3,1,2,4,0>;
		break;
	case 0x85320:
	case 0x67320:
		return &A2QM2q1y_9_eval<4,2,3,0,1>;
		break;

	case 0x07632: //y-Qb-Q+qb+q-
	case 0x05832: //y-Qb-Q+qb+q-
		return &A2QM2q1y_10_eval<0,3,4,1,2>;
		break;
	case 0x20763:
	case 0x20583:
		return &A2QM2q1y_10_eval<1,4,0,2,3>;
		break;
	case 0x32076:
	case 0x32058:
		return &A2QM2q1y_10_eval<2,0,1,3,4>;
		break;
	case 0x63207:
	case 0x83205:
		return &A2QM2q1y_10_eval<3,1,2,4,0>;
		break;
	case 0x76320:
	case 0x58320:
		return &A2QM2q1y_10_eval<4,2,3,0,1>;
		break;

	case 0x08632: //y-Qb+Q+qb+q-
		return &A2QM2q1y_11m_eval<0,3,4,1,2>;
		break;
	case 0x06832: //y-Q+Qb+qb+q-
		return &A2QM2q1y_11p_eval<0,3,4,1,2>;
		break;
	case 0x20863:
		return &A2QM2q1y_11m_eval<1,4,0,2,3>;
		break;
	case 0x20683:
		return &A2QM2q1y_11p_eval<1,4,0,2,3>;
		break;
	case 0x32086:
		return &A2QM2q1y_11m_eval<2,0,1,3,4>;
		break;
	case 0x32068:
		return &A2QM2q1y_11p_eval<2,0,1,3,4>;
		break;
	case 0x63208:
		return &A2QM2q1y_11m_eval<3,1,2,4,0>;
		break;
	case 0x83206:
		return &A2QM2q1y_11p_eval<3,1,2,4,0>;
		break;
	case 0x86320:
		return &A2QM2q1y_11m_eval<4,2,3,0,1>;
		break;
	case 0x68320:
		return &A2QM2q1y_11p_eval<4,2,3,0,1>;
		break;

	case 0x07532: //y-Qb-Q-qb+q-
		return &A2QM2q1y_12m_eval<0,3,4,1,2>;
		break;
	case 0x05732: //y-Q-Qb-qb+q-
		return &A2QM2q1y_12p_eval<0,3,4,1,2>;
		break;
	case 0x20753:
		return &A2QM2q1y_12m_eval<1,4,0,2,3>;
		break;
	case 0x20573:
		return &A2QM2q1y_12p_eval<1,4,0,2,3>;
		break;
	case 0x32075:
		return &A2QM2q1y_12m_eval<2,0,1,3,4>;
		break;
	case 0x32057:
		return &A2QM2q1y_12p_eval<2,0,1,3,4>;
		break;
	case 0x53207:
		return &A2QM2q1y_12m_eval<3,1,2,4,0>;
		break;
	case 0x73205:
		return &A2QM2q1y_12p_eval<3,1,2,4,0>;
		break;
	case 0x75320:
		return &A2QM2q1y_12m_eval<4,2,3,0,1>;
		break;
	case 0x57320:
		return &A2QM2q1y_12p_eval<4,2,3,0,1>;
		break;

	case 0x08523: //y-Qb+Q-qb-q+
	case 0x06723: //y-Q+Qb-qb-q+
		return &A2QM2q1y_13_eval<0,3,4,1,2>;
		break;
	case 0x30852:
	case 0x30672:
		return &A2QM2q1y_13_eval<1,4,0,2,3>;
		break;
	case 0x23085:
	case 0x23067:
		return &A2QM2q1y_13_eval<2,0,1,3,4>;
		break;
	case 0x52308:
	case 0x72306:
		return &A2QM2q1y_13_eval<3,1,2,4,0>;
		break;
	case 0x85230:
	case 0x67230:
		return &A2QM2q1y_13_eval<4,2,3,0,1>;
		break;

	case 0x07623: //y-Qb-Q+qb-q+
	case 0x05823: //y-Qb-Q+qb-q+
		return &A2QM2q1y_14_eval<0,3,4,1,2>;
		break;
	case 0x30762:
	case 0x30582:
		return &A2QM2q1y_14_eval<1,4,0,2,3>;
		break;
	case 0x23076:
	case 0x23058:
		return &A2QM2q1y_14_eval<2,0,1,3,4>;
		break;
	case 0x62307:
	case 0x82305:
		return &A2QM2q1y_14_eval<3,1,2,4,0>;
		break;
	case 0x76230:
	case 0x58230:
		return &A2QM2q1y_14_eval<4,2,3,0,1>;
		break;

	case 0x08623: //y-Qb+Q+qb-q+
		return &A2QM2q1y_15m_eval<0,3,4,1,2>;
		break;
	case 0x06823: //y-Q+Qb+qb-q+
		return &A2QM2q1y_15p_eval<0,3,4,1,2>;
		break;
	case 0x30862:
		return &A2QM2q1y_15m_eval<1,4,0,2,3>;
		break;
	case 0x30682:
		return &A2QM2q1y_15p_eval<1,4,0,2,3>;
		break;
	case 0x23086:
		return &A2QM2q1y_15m_eval<2,0,1,3,4>;
		break;
	case 0x23068:
		return &A2QM2q1y_15p_eval<2,0,1,3,4>;
		break;
	case 0x62308:
		return &A2QM2q1y_15m_eval<3,1,2,4,0>;
		break;
	case 0x82306:
		return &A2QM2q1y_15p_eval<3,1,2,4,0>;
		break;
	case 0x86230:
		return &A2QM2q1y_15m_eval<4,2,3,0,1>;
		break;
	case 0x68230:
		return &A2QM2q1y_15p_eval<4,2,3,0,1>;
		break;

	case 0x07523: //y-Qb-Q-qb-q+
		return &A2QM2q1y_16m_eval<0,3,4,1,2>;
		break;
	case 0x05723: //y-Q-Qb-qb-q+
		return &A2QM2q1y_16p_eval<0,3,4,1,2>;
		break;
	case 0x30752:
		return &A2QM2q1y_16m_eval<1,4,0,2,3>;
		break;
	case 0x30572:
		return &A2QM2q1y_16p_eval<1,4,0,2,3>;
		break;
	case 0x23075:
		return &A2QM2q1y_16m_eval<2,0,1,3,4>;
		break;
	case 0x23057:
		return &A2QM2q1y_16p_eval<2,0,1,3,4>;
		break;
	case 0x52307:
		return &A2QM2q1y_16m_eval<3,1,2,4,0>;
		break;
	case 0x72305:
		return &A2QM2q1y_16p_eval<3,1,2,4,0>;
		break;
	case 0x75230:
		return &A2QM2q1y_16m_eval<4,2,3,0,1>;
		break;
	case 0x57230:
		return &A2QM2q1y_16p_eval<4,2,3,0,1>;
		break;

	default:// We return zero for all other helicity combinations
		cout << "5 pt A2QM2q1y_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}





template complex<R> (*A2s2l_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2s2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2s2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2QM2l_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2QM2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2QM2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2QM1g1y_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2QM1g1y_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2QM1g1y_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2QM1g2l_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2QM1g2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2QM1g2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2s2q2l_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2s2q2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2s2q2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

template complex<R> (*A2QM2q1y_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2QM2q1y_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2QM2q1y_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> (*A2s2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);

template complex<RGMP> (*A2QM2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);

template complex<RGMP> (*A2QM1g1y_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);

template complex<RGMP> (*A2QM1g2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);

template complex<RGMP> (*A2s2q2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);

template complex<RGMP> (*A2QM2q1y_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);

#endif



}

