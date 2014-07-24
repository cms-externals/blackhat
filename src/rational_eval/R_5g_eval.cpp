

#include "R_5g_eval.h"
#include "eval_param.h"
#include "tree_amp.h"

using namespace std;
namespace BH  {


#define _ONLY_X_PART 1
#define X 2
//#define _VERBOSE 1

int helcode_g(const process& pro);

template <class T> complex<T> R5g0(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R5 : -----");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
	// DF : 07/03/2008 changed the overall sign to minus to match my extraction code, does not seem to effect the 6 point result
	- complex<T>(-1,0)/complex<T>(0,3)*(
			ep.spa(0,1)*ep.spa(1,2)/(ep.spb(2,3)*ep.spb(3,4)*ep.spb(4,0))
			+ep.spa(3,4)*ep.spa(4,0)/(ep.spb(0,1)*ep.spb(1,2)*ep.spb(2,3))
			+ep.spa(1,4)*ep.spa(2,3)/(ep.spb(0,1)*ep.spb(2,3)*ep.spb(4,0)))
#endif
;
}

template <class T> complex<T> R5g31(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R5 : +++++");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
	+ complex<T>(1,0)/complex<T>(0,-3)*(
			ep.spb(0,1)*ep.spb(1,2)/(ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,0))
			+ep.spb(3,4)*ep.spb(4,0)/(ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3))
			+ep.spb(1,4)*ep.spb(2,3)/(ep.spa(0,1)*ep.spa(2,3)*ep.spa(4,0)))
#endif
	;
}

template <int i1, int i2, int i3, int i4, int i5, class T> complex<T> R5g1p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R5 : +----");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
	// DF : 07/03/2008 changed the overall sign to minus to match my extraction code, does not seem to effect the 6 point result
	- complex<T>(0,1)/(complex<T>(3,0)*pow(ep.spb(i3,i4),2))*(
					-pow(ep.spa(i2,i5),3)/(ep.spa(i1,i2)*ep.spa(i5,i1))
					+pow(ep.spb(i1,i4),3)*ep.spa(i4,i5)*ep.spb(i3,i5)/(ep.spb(i1,i2)*ep.spb(i2,i3)*pow(ep.spb(i4,i5),2))
					-pow(ep.spb(i1,i3),3)*ep.spa(i3,i2)*ep.spb(i4,i2)/(ep.spb(i1,i5)*ep.spb(i5,i4)*pow(ep.spb(i3,i2),2)))
#endif
	;
}

template <int i1, int i2, int i3, int i4, int i5, class T> complex<T> R5g1m(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE6("R5 : -++++",i1,i2,i3,i4,i5);
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
	+ complex<T>(0,1)/(complex<T>(3,0)*pow(ep.spa(i3,i4),2))*(
					-pow(ep.spb(i2,i5),3)/(ep.spb(i1,i2)*ep.spb(i5,i1))
					+pow(ep.spa(i1,i4),3)*ep.spb(i4,i5)*ep.spa(i3,i5)/(ep.spa(i1,i2)*ep.spa(i2,i3)*pow(ep.spa(i4,i5),2))
					-pow(ep.spa(i1,i3),3)*ep.spb(i3,i2)*ep.spa(i4,i2)/(ep.spa(i1,i5)*ep.spa(i5,i4)*pow(ep.spa(i3,i2),2)))

#endif
	;
}

template <int i1, int i2, int i3, int i4, int i5, class T> complex<T> R5g2m(const eval_param<T>& ep,const mass_param_coll& mpc)
{
complex<T> result(0,0);
complex<T> Xpart=complex<T>(X,0)/complex<T>(9,0)*A05g2m_eval<i1,i2>(ep,mpc);
	result+=Xpart;
#if _ONLY_X_PART

	complex<T> Rpart=-(complex<T>(0,1)/complex<T>(3,0)*(
				ep.spa(i3,i5)*pow(ep.spb(i3,i5),3)/(ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spb(i5,i1)*ep.spb(i1,i2)*ep.spb(i2,i3))
				-pow(ep.spb(i3,i5),2)*ep.spa(i1,i2)/(ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spb(i5,i1)*ep.spb(i2,i3))
				+complex<T>(1,0)/complex<T>(2,0)*ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spa(i4,i1)*ep.spa(i4,i2)*ep.spa(i1,i2)/(ep.spa(i3,i4)*ep.spa(i4,i5)*ep.s(i5,i1)*ep.s(i2,i3))));
				/*CRhat*/
complex<T> Cpart=-(-complex<T>(0,1)/complex<T>(6,0)*(ep.spb(i3,i4)*ep.spa(i4,i1)*ep.spa(i2,i4)*ep.spb(i4,i5)*(ep.spa(i2,i3)*ep.spb(i3,i4)*ep.spa(i4,i1)+ep.spa(i2,i4)*ep.spb(i4,i5)*ep.spa(i5,i1))
				/(ep.spa(i3,i4)*ep.spa(i4,i5))*(ep.s(i2,i3)/ep.s(i5,i1)-ep.s(i5,i1)/ep.s(i2,i3))/pow(ep.s(i5,i1)-ep.s(i2,i3),3)));
#if _VERBOSE
	_MESSAGE6("R5 : --+++",i1,i2,i3,i4,i5);
	_PRINT(Xpart);
	_PRINT(Rpart);
	_PRINT(Cpart);
#endif

	result+=Rpart;
	result+=Cpart;

	#endif
	;

	return result;

}

template <int i1, int i2, int i3, int i4, int i5, class T> complex<T> R5g2p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE6("R5 : ++---",i1,i2,i3,i4,i5);
#endif
	return 			+complex<T>(X,0)/complex<T>(9,0)*A05g2p_eval<i1,i2>(ep,mpc)
#if _ONLY_X_PART


	+(complex<T>(0,1)/complex<T>(3,0)*(
				ep.spb(i3,i5)*pow(ep.spa(i3,i5),3)/(ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spa(i5,i1)*ep.spa(i1,i2)*ep.spa(i2,i3))
				-pow(ep.spa(i3,i5),2)*ep.spb(i1,i2)/(ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spa(i5,i1)*ep.spa(i2,i3))
				+complex<T>(1,0)/complex<T>(2,0)*ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spb(i4,i1)*ep.spb(i4,i2)*ep.spb(i1,i2)/(ep.spb(i3,i4)*ep.spb(i4,i5)*ep.s(i5,i1)*ep.s(i2,i3)))
					/*CRhat*/
				-complex<T>(0,1)/complex<T>(6,0)*(ep.spa(i3,i4)*ep.spb(i4,i1)*ep.spb(i2,i4)*ep.spa(i4,i5)*(ep.spb(i2,i3)*ep.spa(i3,i4)*ep.spb(i4,i1)+ep.spb(i2,i4)*ep.spa(i4,i5)*ep.spb(i5,i1))
				/(ep.spb(i3,i4)*ep.spb(i4,i5))*(ep.s(i2,i3)/ep.s(i5,i1)-ep.s(i5,i1)/ep.s(i2,i3))/pow(ep.s(i5,i1)-ep.s(i2,i3),3)));
#endif
	;

}

template <int i1, int i2, int i3, int i4, int i5, class T> complex<T> R5gmpm(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE6("R5 : -+-++",i1,i2,i3,i4,i5);
#endif
	return (complex<T>(0,1)/complex<T>(3,0)*(
			complex<T>(X,0)/complex<T>(3,0)*pow(ep.spa(i1,i3),4)/(ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i5,i1))
#if _ONLY_X_PART
			+pow(ep.spb(i2,i4),2)*pow(ep.spb(i2,i5),2)/(ep.spb(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i4)*ep.spa(i4,i5)*ep.spb(i5,i1))
			-ep.spa(i1,i2)*pow(ep.spa(i4,i1),2)*pow(ep.spb(i2,i4),3)/(ep.spa(i4,i5)*ep.spa(i5,i1)*ep.spa(i2,i4)*ep.spb(i2,i3)*ep.spb(i3,i4)*ep.s(i5,i1))
			+ep.spa(i3,i2)*pow(ep.spa(i5,i3),2)*pow(ep.spb(i2,i5),3)/(ep.spa(i5,i4)*ep.spa(i4,i3)*ep.spa(i2,i5)*ep.spb(i2,i1)*ep.spb(i1,i5)*ep.s(i3,i4))
			+complex<T>(1,0)/complex<T>(2,0)*pow(ep.spa(i1,i3),2)*ep.spb(i2,i4)*ep.spb(i2,i5)/(ep.s(i3,i4)*ep.spa(i4,i5)*ep.s(i5,i1))
#endif

	)
#if _ONLY_X_PART
	/*CRhat*/
			+complex<T>(0,1)*(-(((T(1)/(T(1) - ep.s(i2,i3)/ep.s(i5,i1)) + T(1)/(T(1) - ep.s(i3,i4)/ep.s(i5,i1)))*ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*
	        pow(ep.spa(i4,i1),2)*pow(ep.spb(i2,i4),2))/
	      (pow(ep.s(i5,i1),2)*pow(ep.spa(i2,i4),2)*ep.spa(i4,i5)*ep.spa(i5,i1))) -
	   ((ep.s(i2,i3)/ep.s(i5,i1) - ep.s(i5,i1)/ep.s(i2,i3))*pow(ep.spa(i2,i3),2)*pow(ep.spa(i4,i1),3)*
	      pow(ep.spb(i2,i4),3))/(T(3)*pow(ep.s(i5,i1) - ep.s(i2,i3),3)*ep.spa(i2,i4)*ep.spa(i4,i5)*ep.spa(i5,i1)) +
	   ((T(1)/(T(1) - ep.s(i1,i2)/ep.s(i3,i4)) + T(1)/(T(1) - ep.s(i5,i1)/ep.s(i3,i4)))*ep.spa(i1,i5)*ep.spa(i2,i1)*ep.spa(i3,i2)*
	      pow(ep.spa(i5,i3),2)*pow(ep.spb(i2,i5),2))/(pow(ep.s(i3,i4),2)*pow(ep.spa(i2,i5),2)*ep.spa(i4,i3)*ep.spa(i5,i4))
	     + ((ep.s(i1,i2)/ep.s(i3,i4) - ep.s(i3,i4)/ep.s(i1,i2))*pow(ep.spa(i2,i1),2)*pow(ep.spa(i5,i3),3)*
	      pow(ep.spb(i2,i5),3))/(T(3)*pow(ep.s(i3,i4) - ep.s(i1,i2),3)*ep.spa(i2,i5)*ep.spa(i4,i3)*ep.spa(i5,i4)) -
	   ((ep.s(i3,i4)/ep.s(i5,i1) - ep.s(i5,i1)/ep.s(i3,i4))*
	      ((T(2)*pow(ep.spa(i1,i2),2)*pow(ep.spa(i3,i4),2)*ep.spa(i4,i1)*pow(ep.spb(i2,i4),3))/
	         (T(3)*ep.spa(i2,i4)*ep.spa(i4,i5)*ep.spa(i5,i1)) -
	        (T(2)*pow(ep.spa(i1,i5),2)*pow(ep.spa(i3,i2),2)*ep.spa(i5,i3)*pow(ep.spb(i2,i5),3))/
	         (T(3)*ep.spa(i2,i5)*ep.spa(i4,i3)*ep.spa(i5,i4)) +
	        (ep.spa(i1,i3)*ep.spb(i2,i4)*ep.spb(i2,i5)*
	           (-(ep.spa(i2,i1)*ep.spa(i3,i4)*ep.spb(i4,i2)) + ep.spa(i1,i5)*ep.spa(i2,i3)*ep.spb(i5,i2)))/(T(3)*ep.spa(i4,i5))))/
	    (T(2)*pow(ep.s(i5,i1) - ep.s(i3,i4),3)))
#endif
	   ) ;
}

template <int i1, int i2, int i3, int i4, int i5, class T> complex<T> R5gpmp(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE6("R5 : +-+--",i1,i2,i3,i4,i5);
#endif
	return -(complex<T>(0,1)/complex<T>(3,0)*(
			complex<T>(X,0)/complex<T>(3,0)*pow(ep.spb(i1,i3),4)/(ep.spb(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i5,i1))
#if _ONLY_X_PART
			+pow(ep.spa(i2,i4),2)*pow(ep.spa(i2,i5),2)/(ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*ep.spb(i4,i5)*ep.spa(i5,i1))
			-ep.spb(i1,i2)*pow(ep.spb(i4,i1),2)*pow(ep.spa(i2,i4),3)/(ep.spb(i4,i5)*ep.spb(i5,i1)*ep.spb(i2,i4)*ep.spa(i2,i3)*ep.spa(i3,i4)*ep.s(i5,i1))
			+ep.spb(i3,i2)*pow(ep.spb(i5,i3),2)*pow(ep.spa(i2,i5),3)/(ep.spb(i5,i4)*ep.spb(i4,i3)*ep.spb(i2,i5)*ep.spa(i2,i1)*ep.spa(i1,i5)*ep.s(i3,i4))
			+complex<T>(1,0)/complex<T>(2,0)*pow(ep.spb(i1,i3),2)*ep.spa(i2,i4)*ep.spa(i2,i5)/(ep.s(i3,i4)*ep.spb(i4,i5)*ep.s(i5,i1))
#endif
			)
#if _ONLY_X_PART
		/*CRhat*/
			+complex<T>(0,1)*(-(((T(1)/(T(1) - ep.s(i2,i3)/ep.s(i5,i1)) + T(1)/(T(1) - ep.s(i3,i4)/ep.s(i5,i1)))*ep.spb(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i4)*
		        pow(ep.spb(i4,i1),2)*pow(ep.spa(i2,i4),2))/
		      (pow(ep.s(i5,i1),2)*pow(ep.spb(i2,i4),2)*ep.spb(i4,i5)*ep.spb(i5,i1))) -
		   ((ep.s(i2,i3)/ep.s(i5,i1) - ep.s(i5,i1)/ep.s(i2,i3))*pow(ep.spb(i2,i3),2)*pow(ep.spb(i4,i1),3)*
		      pow(ep.spa(i2,i4),3))/(T(3)*pow(ep.s(i5,i1) - ep.s(i2,i3),3)*ep.spb(i2,i4)*ep.spb(i4,i5)*ep.spb(i5,i1)) +
		   ((T(1)/(T(1) - ep.s(i1,i2)/ep.s(i3,i4)) + T(1)/(T(1) - ep.s(i5,i1)/ep.s(i3,i4)))*ep.spb(i1,i5)*ep.spb(i2,i1)*ep.spb(i3,i2)*
		      pow(ep.spb(i5,i3),2)*pow(ep.spa(i2,i5),2))/(pow(ep.s(i3,i4),2)*pow(ep.spb(i2,i5),2)*ep.spb(i4,i3)*ep.spb(i5,i4))
		     + ((ep.s(i1,i2)/ep.s(i3,i4) - ep.s(i3,i4)/ep.s(i1,i2))*pow(ep.spb(i2,i1),2)*pow(ep.spb(i5,i3),3)*
		      pow(ep.spa(i2,i5),3))/(T(3)*pow(ep.s(i3,i4) - ep.s(i1,i2),3)*ep.spb(i2,i5)*ep.spb(i4,i3)*ep.spb(i5,i4)) -
		   ((ep.s(i3,i4)/ep.s(i5,i1) - ep.s(i5,i1)/ep.s(i3,i4))*
		      ((T(2)*pow(ep.spb(i1,i2),2)*pow(ep.spb(i3,i4),2)*ep.spb(i4,i1)*pow(ep.spa(i2,i4),3))/
		         (T(3)*ep.spb(i2,i4)*ep.spb(i4,i5)*ep.spb(i5,i1)) -
		        (T(2)*pow(ep.spb(i1,i5),2)*pow(ep.spb(i3,i2),2)*ep.spb(i5,i3)*pow(ep.spa(i2,i5),3))/
		         (T(3)*ep.spb(i2,i5)*ep.spb(i4,i3)*ep.spb(i5,i4)) +
		        (ep.spb(i1,i3)*ep.spa(i2,i4)*ep.spa(i2,i5)*
		           (-(ep.spb(i2,i1)*ep.spb(i3,i4)*ep.spa(i4,i2)) + ep.spb(i1,i5)*ep.spb(i2,i3)*ep.spa(i5,i2)))/(T(3)*ep.spb(i4,i5))))/
		    (T(2)*pow(ep.s(i5,i1) - ep.s(i3,i4),3)))
#endif
    );
}




//R5g0 and R5g15 already in rat_ampl;
			template <class T> complex<T> R5g1(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1p<0,1,2,3,4>(ep,mpc);}
			template <class T> complex<T> R5g2(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1p<1,2,3,4,0>(ep,mpc);}
			template <class T> complex<T> R5g3(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2p<0,1,2,3,4>(ep,mpc);}
			template <class T> complex<T> R5g4(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1p<2,3,4,0,1>(ep,mpc);}
			template <class T> complex<T> R5g5(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gpmp<0,1,2,3,4>(ep,mpc);}
			template <class T> complex<T> R5g6(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2p<1,2,3,4,0>(ep,mpc);}
			template <class T> complex<T> R5g7(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2m<3,4,0,1,2>(ep,mpc);}
			template <class T> complex<T> R5g8(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1p<3,4,0,1,2>(ep,mpc);}
			template <class T> complex<T> R5g9(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gpmp<3,4,0,1,2>(ep,mpc);}
			template <class T> complex<T> R5g10(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gpmp<1,2,3,4,0>(ep,mpc);}
			template <class T> complex<T> R5g11(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gmpm<2,3,4,0,1>(ep,mpc);}
			template <class T> complex<T> R5g12(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2p<2,3,4,0,1>(ep,mpc);}
			template <class T> complex<T> R5g13(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gmpm<4,0,1,2,3>(ep,mpc);}
			template <class T> complex<T> R5g14(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2m<4,0,1,2,3>(ep,mpc);}
			template <class T> complex<T> R5g15(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1m<4,0,1,2,3>(ep,mpc);}
			template <class T> complex<T> R5g16(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1p<4,0,1,2,3>(ep,mpc);}
			template <class T> complex<T> R5g17(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2p<4,0,1,2,3>(ep,mpc);}
			template <class T> complex<T> R5g18(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gpmp<4,0,1,2,3>(ep,mpc);}
			template <class T> complex<T> R5g19(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2m<2,3,4,0,1>(ep,mpc);}
			template <class T> complex<T> R5g20(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gpmp<2,3,4,0,1>(ep,mpc);}
			template <class T> complex<T> R5g21(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gmpm<1,2,3,4,0>(ep,mpc);}
			template <class T> complex<T> R5g22(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gmpm<3,4,0,1,2>(ep,mpc);}
			template <class T> complex<T> R5g23(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1m<3,4,0,1,2>(ep,mpc);}
			template <class T> complex<T> R5g24(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2p<3,4,0,1,2>(ep,mpc);}
			template <class T> complex<T> R5g25(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2m<1,2,3,4,0>(ep,mpc);}
			template <class T> complex<T> R5g26(const eval_param<T>& ep,const mass_param_coll& mpc){return R5gmpm<0,1,2,3,4>(ep,mpc);}
			template <class T> complex<T> R5g27(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1m<2,3,4,0,1>(ep,mpc);}
			template <class T> complex<T> R5g28(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g2m<0,1,2,3,4>(ep,mpc);}
			template <class T> complex<T> R5g29(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1m<1,2,3,4,0>(ep,mpc);}
			template <class T> complex<T> R5g30(const eval_param<T>& ep,const mass_param_coll& mpc){return R5g1m<0,1,2,3,4>(ep,mpc);}



			//------------------  NF   -----------------------------




template <class T> complex<T> (*R5g_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	switch (hc) {

	case 0 : return &R5g0;
	case 31 : return &R5g31;

	case 1 : return &R5g1;
	case 2 : return &R5g2;
	case 4 : return &R5g4;
	case 8 : return &R5g8;
	case 16 : return &R5g16;

	case 15 : return &R5g15;
	case 23 : return &R5g23;
	case 27 : return &R5g27;
	case 29 : return &R5g29;
	case 30 : return &R5g30;

#if _RECURSION_TESTING_MODE == 0

	case 3 : return &R5g3;
	case 5 : return &R5g5;
	case 6 : return &R5g6;
	case 7 : return &R5g7;
	case 9 : return &R5g9;
	case 10 : return &R5g10;
	case 11 : return &R5g11;
	case 12 : return &R5g12;
	case 13 : return &R5g13;
	case 14 : return &R5g14;
	case 17 : return &R5g17;
	case 18 : return &R5g18;
	case 19 : return &R5g19;
	case 20 : return &R5g20;
	case 21 : return &R5g21;
	case 22 : return &R5g22;
	case 24 : return &R5g24;
	case 25 : return &R5g25;
	case 26 : return &R5g26;
	case 28 : return &R5g28;

#endif

	}
}


template <class T> complex<T> R5g0_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g0(ep,mpc);}
template <class T> complex<T> R5g1_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g1(ep,mpc);}
template <class T> complex<T> R5g2_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g2(ep,mpc);}
template <class T> complex<T> R5g3_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g3(ep,mpc);}
template <class T> complex<T> R5g4_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g4(ep,mpc);}
template <class T> complex<T> R5g5_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g5(ep,mpc);}
template <class T> complex<T> R5g6_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g6(ep,mpc);}
template <class T> complex<T> R5g7_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g7(ep,mpc);}
template <class T> complex<T> R5g8_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g8(ep,mpc);}
template <class T> complex<T> R5g9_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g9(ep,mpc);}
template <class T> complex<T> R5g10_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g10(ep,mpc);}
template <class T> complex<T> R5g11_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g11(ep,mpc);}
template <class T> complex<T> R5g12_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g12(ep,mpc);}
template <class T> complex<T> R5g13_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g13(ep,mpc);}
template <class T> complex<T> R5g14_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g14(ep,mpc);}
template <class T> complex<T> R5g15_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g15(ep,mpc);}
template <class T> complex<T> R5g16_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g16(ep,mpc);}
template <class T> complex<T> R5g17_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g17(ep,mpc);}
template <class T> complex<T> R5g18_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g18(ep,mpc);}
template <class T> complex<T> R5g19_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g19(ep,mpc);}
template <class T> complex<T> R5g20_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g20(ep,mpc);}
template <class T> complex<T> R5g21_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g21(ep,mpc);}
template <class T> complex<T> R5g22_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g22(ep,mpc);}
template <class T> complex<T> R5g23_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g23(ep,mpc);}
template <class T> complex<T> R5g24_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g24(ep,mpc);}
template <class T> complex<T> R5g25_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g25(ep,mpc);}
template <class T> complex<T> R5g26_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g26(ep,mpc);}
template <class T> complex<T> R5g27_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g27(ep,mpc);}
template <class T> complex<T> R5g28_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g28(ep,mpc);}
template <class T> complex<T> R5g29_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g29(ep,mpc);}
template <class T> complex<T> R5g30_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g30(ep,mpc);}
template <class T> complex<T> R5g31_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R5g31(ep,mpc);}




template <class T> complex<T> (*R5g_nf_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	switch (hc) {
	case 0 : return &R5g0_nf;
	case 31 : return &R5g31_nf;

	case 1 : return &R5g1_nf;
	case 2 : return &R5g2_nf;
	case 4 : return &R5g4_nf;
	case 8 : return &R5g8_nf;
	case 16 : return &R5g16_nf;

	case 15 : return &R5g15_nf;
	case 23 : return &R5g23_nf;
	case 27 : return &R5g27_nf;
	case 29 : return &R5g29_nf;
	case 30 : return &R5g30_nf;

#if _RECURSION_TESTING_MODE == 0

	case 3 : return &R5g3_nf;
	case 5 : return &R5g5_nf;
	case 6 : return &R5g6_nf;
	case 7 : return &R5g7_nf;
	case 9 : return &R5g9_nf;
	case 10 : return &R5g10_nf;
	case 11 : return &R5g11_nf;
	case 12 : return &R5g12_nf;
	case 13 : return &R5g13_nf;
	case 14 : return &R5g14_nf;
	case 17 : return &R5g17_nf;
	case 18 : return &R5g18_nf;
	case 19 : return &R5g19_nf;
	case 20 : return &R5g20_nf;
	case 21 : return &R5g21_nf;
	case 22 : return &R5g22_nf;
	case 24 : return &R5g24_nf;
	case 25 : return &R5g25_nf;
	case 26 : return &R5g26_nf;
	case 28 : return &R5g28_nf;
#endif
	default:  return 0;
	}
}




template complex<R> (*R5g_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R5g_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R5g_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R5g_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif



template complex<R> (*R5g_nf_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R5g_nf_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R5g_nf_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R5g_nf_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

}


