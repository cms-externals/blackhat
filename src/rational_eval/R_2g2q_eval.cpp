

#include "R_2g2q_eval.h"
#include "eval_param.h"
#include "helcode.h"

using namespace std;
namespace BH  {


#define _ONLY_X_PART 1
#define X 2
#define _VERBOSE 0


template <class T> complex<T> R2g2q0(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q2g : ----");
#endif
	return complex<T>(0,0);
}
template <class T> complex<T> R2g2q15(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q2g : ++++");
#endif
	return complex<T>(0,0);
}


template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mppp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -+++ L");
#endif
	return complex<T>(0,1)/complex<T>(2,0)*ep.spa(i1,i2)*ep.spb(i4,i2)/(ep.spa(i2,i3)*ep.spa(i3,i4)) + complex<T>(0,1)/complex<T>(3,0)*ep.spb(i3,i2)*ep.spb(i4,i2)/(ep.spa(i3,i4)*ep.spb(i2,i1)) ;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +--- L");
#endif
	return -(
			complex<T>(0,1)/complex<T>(2,0)*ep.spb(i2,i1)*ep.spa(i2,i4)/(ep.spb(i3,i2)*ep.spb(i4,i3)) + complex<T>(0,1)/complex<T>(3,0)*ep.spa(i2,i3)*ep.spa(i2,i4)/(ep.spb(i4,i3)*ep.spa(i1,i2))
			);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mpmp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -+-+ L");
#endif
	return
(complex<T>(0,-1)/complex<T>(2,0)*pow(ep.spa(i1,i3),3))/
    (ep.spa(i1,i2)*ep.spa(i1,i4)*ep.spa(i3,i4))
     - (complex<T>(0,1)/complex<T>(2,0)*ep.s(i1,i2)*
      pow(ep.spa(i1,i3),3))/
    (ep.s(i1,i3)*ep.spa(i1,i2)*ep.spa(i1,i4)*
      ep.spa(i3,i4));
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmpm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +-+- L");
#endif
	return -(
		(complex<T>(0,-1)/complex<T>(2,0)*pow(ep.spb(i3,i1),3))/
			(ep.spb(i2,i1)*ep.spb(i4,i1)*ep.spb(i4,i3))
			 - (complex<T>(0,1)/complex<T>(2,0)*ep.s(i1,i2)*
			  pow(ep.spb(i3,i1),3))/
			(ep.s(i1,i3)*ep.spb(i2,i1)*ep.spb(i4,i1)*
			  ep.spb(i4,i3))
	)
      ;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mppm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -++- L");
#endif
	return (complex<T>(0,-1)/complex<T>(2,0)*pow(ep.spa(i1,i4),2)*ep.spa(i2,i4))/
	   (ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4));
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmmp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +--+ L");
#endif
	return (complex<T>(0,1)/complex<T>(2,0)*pow(ep.spb(i4,i1),2)*ep.spb(i4,i2))/
	   (ep.spb(i2,i1)*ep.spb(i3,i2)*ep.spb(i4,i3));
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mpmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -+-- L");
#endif
	return (complex<T>(0,-1)/complex<T>(2,0)*ep.spa(i1,i3)*ep.spb(i2,i1))/
    (ep.spb(i4,i1)*ep.spb(i4,i3)) -
   (complex<T>(0,1)/complex<T>(3,0)*ep.spa(i1,i3)*
      ep.spb(i2,i1)*ep.s(i1,i4))/
    (ep.spb(i4,i1)*ep.spb(i4,i3)*ep.s(i1,i2));
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmpp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +-++ L");
#endif
	return -(
			(complex<T>(0,-1)/complex<T>(2,0)*ep.spb(i3,i1)*ep.spa(i1,i2))/
			    (ep.spa(i1,i4)*ep.spa(i3,i4)) -
			   (complex<T>(0,1)/complex<T>(3,0)*ep.spb(i3,i1)*
			      ep.spa(i1,i2)*ep.s(i1,i4))/
			    (ep.spa(i1,i4)*ep.spa(i3,i4)*ep.s(i1,i2))
    );
	}




template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_sl_mppp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q sl : -+++ L");
#endif
	return complex<T>(0,-1)/complex<T>(2,0)*ep.spa(i1,i3)*ep.spb(i2,i4)/(ep.spa(i2,i3)*ep.spa(i3,i4));
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_sl_pmmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q sl : +--- L");
#endif
	return - (
			complex<T>(0,-1)/complex<T>(2,0)*ep.spb(i1,i3)*ep.spa(i2,i4)/(ep.spb(i2,i3)*ep.spb(i3,i4))
			);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_sl_ppmp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q sl : ++-+ L");
#endif
	return -R2g2q_sl_mppp_L<i3,i4,i1,i2>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_sl_mmpm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q sl : --+- L");
#endif
			return -R2g2q_sl_pmmm_L<i3,i4,i1,i2>(ep,mpc);
}




//===================  RT  =================================

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mppp_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -+++ SLC");
#endif
	return complex<T>(0,-1)/complex<T>(2,0)*ep.spa(i1,i2)*ep.spb(i4,i2)/(ep.spa(i2,i3)*ep.spa(i3,i4)) ;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmmm_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +--- SLC");
#endif
	return -complex<T>(0,-1)/complex<T>(2,0)*ep.spb(i1,i2)*ep.spa(i4,i2)/(ep.spb(i2,i3)*ep.spb(i3,i4));
}


template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mpmp_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -+-+ SLC ");
#endif
	return	-R2g2q_mpmp_L<i1,i2,i3,i4>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmpm_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +-+- SLC ");
#endif
	return	-R2g2q_pmpm_L<i1,i2,i3,i4>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mppm_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -++- SLC ");
#endif
	return	-R2g2q_mppm_L<i1,i2,i3,i4>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmmp_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +--+ SLC ");
#endif
	return	-R2g2q_pmmp_L<i1,i2,i3,i4>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mpmm_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -+-- SLC");
#endif
	return (complex<T>(0,1)/complex<T>(2,0)*ep.spa(i1,i3)*ep.spb(i2,i1))/
    (ep.spb(i4,i1)*ep.spb(i4,i3)) ;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmpp_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +-++ SLC");
#endif
	return -(
			(complex<T>(0,1)/complex<T>(2,0)*ep.spa(i1,i2)*ep.spb(i3,i1))/
    (ep.spa(i1,i4)*ep.spa(i3,i4))
    );
	}


template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_sl_mppp_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q sl : -+-+ L");
#endif
	return -R2g2q_sl_mppp_L<i1,i2,i3,i4,T>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_sl_pmmm_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q sl : -+-+ L");
#endif

			return -R2g2q_sl_pmmm_L<i1,i2,i3,i4,T>(ep,mpc);
			;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_sl_ppmp_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q sl : ++-+ L");
#endif
	return -R2g2q_sl_mppp_SLC<i3,i4,i1,i2>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_sl_mmpm_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q sl : --+- L");
#endif
			return -R2g2q_sl_pmmm_SLC<i3,i4,i1,i2>(ep,mpc);
}

//===================  nf  =================================

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mppp_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -+++ nf");
#endif
	return - complex<T>(0,1)/complex<T>(3,0)*ep.spb(i3,i2)*ep.spb(i4,i2)/(ep.spa(i3,i4)*ep.spb(i2,i1))  ;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmmm_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +--- nf");
#endif
	return + complex<T>(0,1)/complex<T>(3,0)*ep.spa(i3,i2)*ep.spa(i4,i2)/(ep.spb(i3,i4)*ep.spa(i2,i1));
}


template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_mpmm_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : -+-- nf");
#endif
	return
	   (complex<T>(0,1)/complex<T>(3,0)*ep.spa(i1,i3)*
	      ep.spb(i2,i1)*ep.s(i1,i4))/
	    (ep.spb(i4,i1)*ep.spb(i4,i3)*ep.s(i1,i2)) ;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R2g2q_pmpp_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2g2q : +-++ nf");
#endif
	return -(complex<T>(0,1)/complex<T>(3,0)*ep.spb(i1,i3)*
		      ep.spa(i2,i1)*ep.s(i1,i4))/
		    (ep.spa(i4,i1)*ep.spa(i4,i3)*ep.s(i1,i2)) ;
	}


// case with qqb

template <class T> complex<T> Rqqbgg_mppp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R qm qbp p p  sl : -+++ L");
#endif
	return R2g2q_mppp_SLC(ep,mpc) ;
}

template <class T> complex<T> Rqqbgg_pmmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("qp qbm m m : +--- L");
#endif
	return R2g2q_pmmm_SLC(ep,mpc) ;
}

template <class T> complex<T> Rqqbgg_pmpp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("qp qbm p p : +-++ L");
#endif
	return R2g2q_pmpp_SLC(ep,mpc) ;
}

template <class T> complex<T> Rqqbgg_mpmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("qm qbp m m : -+-- SLC");
#endif
	return R2g2q_mpmm_SLC(ep,mpc) ;
}


template <class T> complex<T> Rqqbgg_mppp_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R qm qbp p p  sl : -+++ SLC");
#endif
	return R2g2q_mppp_L(ep,mpc) ;
}

template <class T> complex<T> Rqqbgg_pmmm_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("qp qbm m m : +--- SLC");
#endif
	return R2g2q_pmmm_L(ep,mpc) ;
}

template <class T> complex<T> Rqqbgg_pmpp_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("qp qbm p p : +-++ SLC");
#endif
	return R2g2q_pmpp_L(ep,mpc) ;
}

template <class T> complex<T> Rqqbgg_mpmm_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("qm qbp m m : -+-- SLC");
#endif
	return R2g2q_mpmm_L(ep,mpc) ;
}




/* *************** table of switch values ************* */

/*
144: mmqmqp
96: mmqpqm
132: mqmmqp
36: mqmqpm
228: mqmqpp
180: mqmpqp
72: mqpmqm
24: mqpqmm
216: mqpqmp
120: mqppqm
156: mpqmqp
108: mpqpqm
129: qmmmqp
33: qmmqpm
225: qmmqpp
177: qmmpqp
9: qmqpmm
201: qmqpmp
57: qmqppm
249: qmqppp
                         141: qmpmqp
45: qmpqpm
237: qmpqpp
189: qmppqp
66: qpmmqm
18: qpmqmm
210: qpmqmp
114: qpmpqm
6: qpqmmm
198: qpqmmp
54: qpqmpm
246: qpqmpp
78: qppmqm
30: qppqmm
222: qppqmp
126: qpppqm
147: pmqmqp
99: pmqpqm
135: pqmmqp
39: pqmqpm
231: pqmqpp
183: pqmpqp
75: pqpmqm
27: pqpqmm
219: pqpqmp
123: pqppqm
159: ppqmqp
111: ppqpqm
*/

template <class T> complex<T> R2g2q_L_249(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_L<0,1,2,3>(ep,mpc);}   // qm qp p p
template <class T> complex<T> R2g2q_L_231(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_L<1,2,3,0>(ep,mpc);}   // p qm qp p
template <class T> complex<T> R2g2q_L_159(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_L<2,3,0,1>(ep,mpc);}   // p p qm qp
template <class T> complex<T> R2g2q_L_126(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_L<3,0,1,2>(ep,mpc);}   // qp p p qm

template <class T> complex<T> R2g2q_L_6  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_L<0,1,2,3>(ep,mpc);}   // qp qm m m
template <class T> complex<T> R2g2q_L_24 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_L<1,2,3,0>(ep,mpc);}   // m qp qm m
template <class T> complex<T> R2g2q_L_96 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_L<2,3,0,1>(ep,mpc);}   // m m qp qm
template <class T> complex<T> R2g2q_L_129(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_L<3,0,1,2>(ep,mpc);}   // qm m m qp

template <class T> complex<T> R2g2q_L_201(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmp_L<0,1,2,3>(ep,mpc);}   // qm qp m p
template <class T> complex<T> R2g2q_L_39 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmp_L<1,2,3,0>(ep,mpc);}   // p qm qp m
template <class T> complex<T> R2g2q_L_156(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmp_L<2,3,0,1>(ep,mpc);}   // m p qm qp
template <class T> complex<T> R2g2q_L_114(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmp_L<3,0,1,2>(ep,mpc);}   // qp m p qm

template <class T> complex<T> R2g2q_L_54 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpm_L<0,1,2,3>(ep,mpc);}   // qp qm p m
template <class T> complex<T> R2g2q_L_216(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpm_L<1,2,3,0>(ep,mpc);}   // m qp qm p
template <class T> complex<T> R2g2q_L_99 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpm_L<2,3,0,1>(ep,mpc);}   // p m qp qm
template <class T> complex<T> R2g2q_L_141(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpm_L<3,0,1,2>(ep,mpc);}   // qm p m qp

template <class T> complex<T> R2g2q_L_57 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppm_L<0,1,2,3>(ep,mpc);}   // qm qp p m
template <class T> complex<T> R2g2q_L_228(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppm_L<1,2,3,0>(ep,mpc);}   // m qm qp p
template <class T> complex<T> R2g2q_L_147(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppm_L<2,3,0,1>(ep,mpc);}   // p m qm qp
template <class T> complex<T> R2g2q_L_78 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppm_L<3,0,1,2>(ep,mpc);}   // qp p m qm

template <class T> complex<T> R2g2q_L_198(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmp_L<0,1,2,3>(ep,mpc);}   // qp qm m p
template <class T> complex<T> R2g2q_L_27 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmp_L<1,2,3,0>(ep,mpc);}   // p qp qm m
template <class T> complex<T> R2g2q_L_108(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmp_L<2,3,0,1>(ep,mpc);}   // m p qp qm
template <class T> complex<T> R2g2q_L_177 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmp_L<3,0,1,2>(ep,mpc);}   // qm m p qp

template <class T> complex<T> R2g2q_L_9   (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_L<0,1,2,3>(ep,mpc);}   // qm qp m m
template <class T> complex<T> R2g2q_L_36  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_L<1,2,3,0>(ep,mpc);}   // m qm qp m
template <class T> complex<T> R2g2q_L_144 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_L<2,3,0,1>(ep,mpc);}   // m m qm qp
template <class T> complex<T> R2g2q_L_66  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_L<3,0,1,2>(ep,mpc);}   // qp m m qm

template <class T> complex<T> R2g2q_L_246 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_L<0,1,2,3>(ep,mpc);}   // qp qm p p
template <class T> complex<T> R2g2q_L_219 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_L<1,2,3,0>(ep,mpc);}   // p qp qm p
template <class T> complex<T> R2g2q_L_111 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_L<2,3,0,1>(ep,mpc);}   // p p qp qm
template <class T> complex<T> R2g2q_L_189 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_L<3,0,1,2>(ep,mpc);}   // qm p p qp

template <class T> complex<T> R2g2q_L_237 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_sl_mppp_L<0,1,2,3>(ep,mpc);}   // qm p qp p
template <class T> complex<T> R2g2q_L_183 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_sl_mppp_L<1,2,3,0>(ep,mpc);}   // p qm p qp
template <class T> complex<T> R2g2q_L_222 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_sl_mppp_L<2,3,0,1>(ep,mpc);}   // qp p qm p
template <class T> complex<T> R2g2q_L_123 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_sl_mppp_L<3,0,1,2>(ep,mpc);}   // p qp p qm

template <class T> complex<T> R2g2q_L_18 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_sl_pmmm_L<0,1,2,3>(ep,mpc);}   // qp m qm m
template <class T> complex<T> R2g2q_L_72 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_sl_pmmm_L<1,2,3,0>(ep,mpc);}   // m qp m qm
template <class T> complex<T> R2g2q_L_33 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_sl_pmmm_L<2,3,0,1>(ep,mpc);}   // qm m qp m
template <class T> complex<T> R2g2q_L_132(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_sl_pmmm_L<3,0,1,2>(ep,mpc);}   // m qm m qp



template <class T> complex<T> R2g2q_SLC_249(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_SLC<0,1,2,3>(ep,mpc);}   // qm qp p p
template <class T> complex<T> R2g2q_SLC_231(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_SLC<1,2,3,0>(ep,mpc);}   // p qm qp p
template <class T> complex<T> R2g2q_SLC_159(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_SLC<2,3,0,1>(ep,mpc);}   // p p qm qp
template <class T> complex<T> R2g2q_SLC_126(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_SLC<3,0,1,2>(ep,mpc);}   // qp p p qm

template <class T> complex<T> R2g2q_SLC_6  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_SLC<0,1,2,3>(ep,mpc);}   // qp qm m m
template <class T> complex<T> R2g2q_SLC_24 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_SLC<1,2,3,0>(ep,mpc);}   // m qp qm m
template <class T> complex<T> R2g2q_SLC_96 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_SLC<2,3,0,1>(ep,mpc);}   // m m qp qm
template <class T> complex<T> R2g2q_SLC_129(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_SLC<3,0,1,2>(ep,mpc);}   // qm m m qp

template <class T> complex<T> R2g2q_SLC_201(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmp_SLC<0,1,2,3>(ep,mpc);}   // qm qp m p
template <class T> complex<T> R2g2q_SLC_39 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmp_SLC<1,2,3,0>(ep,mpc);}   // p qm qp m
template <class T> complex<T> R2g2q_SLC_156(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmp_SLC<2,3,0,1>(ep,mpc);}   // m p qm qp
template <class T> complex<T> R2g2q_SLC_114(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmp_SLC<3,0,1,2>(ep,mpc);}   // qp m p qm

template <class T> complex<T> R2g2q_SLC_54 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpm_SLC<0,1,2,3>(ep,mpc);}   // qp qm p m
template <class T> complex<T> R2g2q_SLC_216(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpm_SLC<1,2,3,0>(ep,mpc);}   // m qp qm p
template <class T> complex<T> R2g2q_SLC_99 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpm_SLC<2,3,0,1>(ep,mpc);}   // p m qp qm
template <class T> complex<T> R2g2q_SLC_141(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpm_SLC<3,0,1,2>(ep,mpc);}   // qm p m qp

template <class T> complex<T> R2g2q_SLC_57 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppm_SLC<0,1,2,3>(ep,mpc);}   // qm qp p m
template <class T> complex<T> R2g2q_SLC_228(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppm_SLC<1,2,3,0>(ep,mpc);}   // m qm qp p
template <class T> complex<T> R2g2q_SLC_147(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppm_SLC<2,3,0,1>(ep,mpc);}   // p m qm qp
template <class T> complex<T> R2g2q_SLC_78 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppm_SLC<3,0,1,2>(ep,mpc);}   // qp p m qm

template <class T> complex<T> R2g2q_SLC_198(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmp_SLC<0,1,2,3>(ep,mpc);}   // qp qm m p
template <class T> complex<T> R2g2q_SLC_27 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmp_SLC<1,2,3,0>(ep,mpc);}   // p qp qm m
template <class T> complex<T> R2g2q_SLC_108(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmp_SLC<2,3,0,1>(ep,mpc);}   // m p qp qm
template <class T> complex<T> R2g2q_SLC_177 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmp_SLC<3,0,1,2>(ep,mpc);}   // qm m p qp

template <class T> complex<T> R2g2q_SLC_9   (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_SLC<0,1,2,3>(ep,mpc);}   // qm qp m m
template <class T> complex<T> R2g2q_SLC_36  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_SLC<1,2,3,0>(ep,mpc);}   // m qm qp m
template <class T> complex<T> R2g2q_SLC_144 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_SLC<2,3,0,1>(ep,mpc);}   // m m qm qp
template <class T> complex<T> R2g2q_SLC_66  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_SLC<3,0,1,2>(ep,mpc);}   // qp m m qm

template <class T> complex<T> R2g2q_SLC_246 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_SLC<0,1,2,3>(ep,mpc);}   // qp qm p p
template <class T> complex<T> R2g2q_SLC_219 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_SLC<1,2,3,0>(ep,mpc);}   // p qp qm p
template <class T> complex<T> R2g2q_SLC_111 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_SLC<2,3,0,1>(ep,mpc);}   // p p qp qm
template <class T> complex<T> R2g2q_SLC_189 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_SLC<3,0,1,2>(ep,mpc);}   // qm p p qp

template <class T> complex<T> R2g2q_SLC_237 (const eval_param<T>& ep,const mass_param_coll& mpc){return -R2g2q_sl_mppp_L<0,1,2,3>(ep,mpc);}   // qm p qp p
template <class T> complex<T> R2g2q_SLC_183 (const eval_param<T>& ep,const mass_param_coll& mpc){return -R2g2q_sl_mppp_L<1,2,3,0>(ep,mpc);}   // p qm p qp
template <class T> complex<T> R2g2q_SLC_222 (const eval_param<T>& ep,const mass_param_coll& mpc){return -R2g2q_sl_mppp_L<2,3,0,1>(ep,mpc);}   // qp p qm p
template <class T> complex<T> R2g2q_SLC_123 (const eval_param<T>& ep,const mass_param_coll& mpc){return -R2g2q_sl_mppp_L<3,0,1,2>(ep,mpc);}   // p qp p qm

template <class T> complex<T> R2g2q_SLC_18 (const eval_param<T>& ep,const mass_param_coll& mpc){return -R2g2q_sl_pmmm_L<0,1,2,3>(ep,mpc);}   // qp m qm m
template <class T> complex<T> R2g2q_SLC_72 (const eval_param<T>& ep,const mass_param_coll& mpc){return -R2g2q_sl_pmmm_L<1,2,3,0>(ep,mpc);}   // m qp m qm
template <class T> complex<T> R2g2q_SLC_33 (const eval_param<T>& ep,const mass_param_coll& mpc){return -R2g2q_sl_pmmm_L<2,3,0,1>(ep,mpc);}   // qm m qp m
template <class T> complex<T> R2g2q_SLC_132(const eval_param<T>& ep,const mass_param_coll& mpc){return -R2g2q_sl_pmmm_L<3,0,1,2>(ep,mpc);}   // m qm m qp


template <class T> complex<T> R2g2q_nf_249(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_nf<0,1,2,3>(ep,mpc);}   // qm qp p p
template <class T> complex<T> R2g2q_nf_231(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_nf<1,2,3,0>(ep,mpc);}   // p qm qp p
template <class T> complex<T> R2g2q_nf_159(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_nf<2,3,0,1>(ep,mpc);}   // p p qm qp
template <class T> complex<T> R2g2q_nf_126(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mppp_nf<3,0,1,2>(ep,mpc);}   // qp p p qm

template <class T> complex<T> R2g2q_nf_6  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_nf<0,1,2,3>(ep,mpc);}   // qp qm m m
template <class T> complex<T> R2g2q_nf_24 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_nf<1,2,3,0>(ep,mpc);}   // m qp qm m
template <class T> complex<T> R2g2q_nf_96 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_nf<2,3,0,1>(ep,mpc);}   // m m qp qm
template <class T> complex<T> R2g2q_nf_129(const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmmm_nf<3,0,1,2>(ep,mpc);}   // qm m m qp

template <class T> complex<T> R2g2q_nf_201(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // qm qp m p
template <class T> complex<T> R2g2q_nf_39 (const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // p qm qp m
template <class T> complex<T> R2g2q_nf_156(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // m p qm qp
template <class T> complex<T> R2g2q_nf_114(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // qp m p qm

template <class T> complex<T> R2g2q_nf_54 (const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // qp qm p m
template <class T> complex<T> R2g2q_nf_216(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // m qp qm p
template <class T> complex<T> R2g2q_nf_99 (const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // p m qp qm
template <class T> complex<T> R2g2q_nf_141(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // qm p m qp

template <class T> complex<T> R2g2q_nf_57 (const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // qm qp p m
template <class T> complex<T> R2g2q_nf_228(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // m qm qp p
template <class T> complex<T> R2g2q_nf_147(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // p m qm qp
template <class T> complex<T> R2g2q_nf_78 (const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // qp p m qm

template <class T> complex<T> R2g2q_nf_198(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // qp qm m p
template <class T> complex<T> R2g2q_nf_27 (const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // p qp qm m
template <class T> complex<T> R2g2q_nf_108(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // m p qp qm
template <class T> complex<T> R2g2q_nf_177(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);}   // qm m p qp

template <class T> complex<T> R2g2q_nf_9   (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_nf<0,1,2,3>(ep,mpc);}   // qm qp m m
template <class T> complex<T> R2g2q_nf_36  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_nf<1,2,3,0>(ep,mpc);}   // m qm qp m
template <class T> complex<T> R2g2q_nf_144 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_nf<2,3,0,1>(ep,mpc);}   // m m qm qp
template <class T> complex<T> R2g2q_nf_66  (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_mpmm_nf<3,0,1,2>(ep,mpc);}   // qp m m qm

template <class T> complex<T> R2g2q_nf_246 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_nf<0,1,2,3>(ep,mpc);}   // qp qm p p
template <class T> complex<T> R2g2q_nf_219 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_nf<1,2,3,0>(ep,mpc);}   // p qp qm p
template <class T> complex<T> R2g2q_nf_111 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_nf<2,3,0,1>(ep,mpc);}   // p p qp qm
template <class T> complex<T> R2g2q_nf_189 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2g2q_pmpp_nf<3,0,1,2>(ep,mpc);}   // qm p p qp



#define _CASE_R2g2q_L_Ptr_eval(K) case K : return &R2g2q_L_ ## K
#define _CASE_R2g2q_SLC_Ptr_eval(K) case K : return &R2g2q_SLC_ ## K

template <class T> complex<T> (*R2g2q_L_Ptr_eval(const process&  pro))(const eval_param<T>& ,const mass_param_coll&) {
	switch (helcode_2q(pro)) {
	case 249 : return &R2g2q_L_249;
//	case 231 : return &R2g2q_L_231;
//	case 159 : return &R2g2q_L_159;
	case 126 : return &R2g2q_SLC_126;

	case   6 : return &R2g2q_L_6;
//	case  24 : return &R2g2q_L_24;
//	case  96 : return &R2g2q_L_96;
	case 129 : return &R2g2q_SLC_129;

	case 201 : return &R2g2q_L_201;
//	case  39 : return &R2g2q_L_39;
//	case 156 : return &R2g2q_L_156;
	case 114 : return &R2g2q_SLC_114;

	case 54  : return &R2g2q_L_54;
//	case 216 : return &R2g2q_L_216;
//	case 99  : return &R2g2q_L_99;
	case 141 : return &R2g2q_SLC_141;

	case 57  : return &R2g2q_L_57;
//	case 228 : return &R2g2q_L_228;
//	case 147 : return &R2g2q_L_147;
	case 78  : return &R2g2q_SLC_78;

	case 198 : return &R2g2q_L_198;
//	case 27  : return &R2g2q_L_27;
//	case 108 : return &R2g2q_L_108;
	case 177 : return &R2g2q_SLC_177;
	_CASE_R2g2q_L_Ptr_eval(9);
//	_CASE_R2g2q_L_Ptr_eval(36);
//	_CASE_R2g2q_L_Ptr_eval(144);
	_CASE_R2g2q_SLC_Ptr_eval(66);

	_CASE_R2g2q_L_Ptr_eval(246);
//	_CASE_R2g2q_L_Ptr_eval(219);
//	_CASE_R2g2q_L_Ptr_eval(111);
	_CASE_R2g2q_SLC_Ptr_eval(189);

	//finite amplitudes

//	 _CASE_R2g2q_L_Ptr_eval(237); // qm p qp p
//	 _CASE_R2g2q_L_Ptr_eval(183); // p qm p qp
//	 _CASE_R2g2q_L_Ptr_eval(222); // qp p qm p
//	 _CASE_R2g2q_L_Ptr_eval(123); // p qp p qm
//
//	 _CASE_R2g2q_L_Ptr_eval(18); // qp m qm m
//	 _CASE_R2g2q_L_Ptr_eval(72); // m qp m qm
//	 _CASE_R2g2q_L_Ptr_eval(33); // qm m qp m
//	 _CASE_R2g2q_L_Ptr_eval(132); // m qm m qp
	// vanishing R pieces:
	case 225: // qmmqpp
	case 135: // pqmmqp
	case 30:  // qppqmm
	case 120: // mqppqm
	case 210:// qpmqmp
	case 180:// mqmpqp
	case 75: //pqpmqm
	case 45: // qmpqpm

	          return &R2g2q0;

	}

	if ( pro == process(qbm,p,qp,p)) return &R2g2q_sl_mppp_L<0,1,2,3,T>;
//	if ( pro == process(p,qbm,p,qp)) return &R2g2q_sl_mppp_L<1,2,3,0,T>;
//	if ( pro == process(qp,p,qbm,p)) return &R2g2q_sl_mppp_L<2,3,0,1,T>;
//	if ( pro == process(p,qp,p,qbm)) return &R2g2q_sl_mppp_L<3,0,1,2,T>;

	if ( pro == process(qbp,m,qm,m)) return &R2g2q_sl_pmmm_L<0,1,2,3,T>;
//	if ( pro == process(m,qbp,m,qm)) return &R2g2q_sl_pmmm_L<1,2,3,0,T>;
//	if ( pro == process(qm,m,qbp,m)) return &R2g2q_sl_pmmm_L<2,3,0,1,T>;
//	if ( pro == process(m,qm,m,qbp)) return &R2g2q_sl_pmmm_L<3,0,1,2,T>;

	if ( pro == process(qbm,m,qp,m)) return &R2g2q_sl_mmpm_L<0,1,2,3,T>;
//	if ( pro == process(m,qbm,m,qp)) return &R2g2q_sl_mmpm_L<1,2,3,0,T>;
//	if ( pro == process(qp,m,qbm,m)) return &R2g2q_sl_mmpm_L<2,3,0,1,T>;
//	if ( pro == process(m,qp,m,qbm)) return &R2g2q_sl_mmpm_L<3,0,1,2,T>;

	if ( pro == process(qbp,p,qm,p)) return &R2g2q_sl_ppmp_L<0,1,2,3,T>;
//	if ( pro == process(p,qbp,p,qm)) return &R2g2q_sl_ppmp_L<1,2,3,0,T>;
//	if ( pro == process(qm,p,qbp,p)) return &R2g2q_sl_ppmp_L<2,3,0,1,T>;
//	if ( pro == process(p,qm,p,qbp)) return &R2g2q_sl_ppmp_L<3,0,1,2,T>;


	return 0;
}





template <class T> complex<T> (*R2g2q_SLC_Ptr_eval(const process& pro))(const eval_param<T>& ,const mass_param_coll&) {
	switch (helcode_2q(pro)) {
	case 249 : return &R2g2q_SLC_249;
//	case 231 : return &R2g2q_SLC_231;
//	case 159 : return &R2g2q_SLC_159;
	case 126 : return &R2g2q_L_126;

	case   6 : return &R2g2q_SLC_6;
//	case  24 : return &R2g2q_SLC_24;
//	case  96 : return &R2g2q_SLC_96;
	case 129 : return &R2g2q_L_129;

	case 201 : return &R2g2q_SLC_201;
//	case  39 : return &R2g2q_SLC_39;
//	case 156 : return &R2g2q_SLC_156;
	case 114 : return &R2g2q_L_114;

	case 54  : return &R2g2q_SLC_54;
//	case 216 : return &R2g2q_SLC_216;
//	case 99  : return &R2g2q_SLC_99;
	case 141 : return &R2g2q_L_141;

	case 57  : return &R2g2q_SLC_57;
//	case 228 : return &R2g2q_SLC_228;
//	case 147 : return &R2g2q_SLC_147;
	case 78  : return &R2g2q_L_78;

	case 198 : return &R2g2q_SLC_198;
//	case 27  : return &R2g2q_SLC_27;
//	case 108 : return &R2g2q_SLC_108;
	case 177 : return &R2g2q_L_177;

	_CASE_R2g2q_SLC_Ptr_eval(9);
//	_CASE_R2g2q_SLC_Ptr_eval(36);
//	_CASE_R2g2q_SLC_Ptr_eval(144);
	_CASE_R2g2q_L_Ptr_eval(66);
	_CASE_R2g2q_SLC_Ptr_eval(246);
//	_CASE_R2g2q_SLC_Ptr_eval(219);
//	_CASE_R2g2q_SLC_Ptr_eval(111);
	_CASE_R2g2q_L_Ptr_eval(189);

	//finite amplitudes
//		 _CASE_R2g2q_SLC_Ptr_eval(237); // qm p qp p
//		 _CASE_R2g2q_SLC_Ptr_eval(183); // p qm p qp
//		 _CASE_R2g2q_SLC_Ptr_eval(222); // qp p qm p
//		 _CASE_R2g2q_SLC_Ptr_eval(123); // p qp p qm
//
//		 _CASE_R2g2q_SLC_Ptr_eval(18); // qp m qm m
//		 _CASE_R2g2q_SLC_Ptr_eval(72); // m qp m qm
//		 _CASE_R2g2q_SLC_Ptr_eval(33); // qm m qp m
//		 _CASE_R2g2q_SLC_Ptr_eval(132); // m qm m qp
	// vanishing R pieces:
	case 225: // qmmqpp
	case 135: // pqmmqp
	case 30:  // qppqmm
	case 120: // mqppqm
	case 210:// qpmqmp
	case 180:// mqmpqp
	case 75: //pqpmqm
	case 45: // qmpqpm
	          return &R2g2q0;

	}

	if ( pro == process(qbm,p,qp,p)) return &R2g2q_sl_mppp_SLC<0,1,2,3,T>;
//	if ( pro == process(p,qbm,p,qp)) return &R2g2q_sl_mppp_SLC<1,2,3,0,T>;
//	if ( pro == process(qp,p,qbm,p)) return &R2g2q_sl_mppp_SLC<2,3,0,1,T>;
//	if ( pro == process(p,qp,p,qbm)) return &R2g2q_sl_mppp_SLC<3,0,1,2,T>;

	if ( pro == process(qbp,m,qm,m)) return &R2g2q_sl_pmmm_SLC<0,1,2,3,T>;
//	if ( pro == process(m,qbp,m,qm)) return &R2g2q_sl_pmmm_SLC<1,2,3,0,T>;
//	if ( pro == process(qm,m,qbp,m)) return &R2g2q_sl_pmmm_SLC<2,3,0,1,T>;
//	if ( pro == process(m,qm,m,qbp)) return &R2g2q_sl_pmmm_SLC<3,0,1,2,T>;

	if ( pro == process(qbm,m,qp,m)) return &R2g2q_sl_mmpm_SLC<0,1,2,3,T>;
//	if ( pro == process(m,qbm,m,qp)) return &R2g2q_sl_mmpm_SLC<1,2,3,0,T>;
//	if ( pro == process(qp,m,qbm,m)) return &R2g2q_sl_mmpm_SLC<2,3,0,1,T>;
//	if ( pro == process(m,qp,m,qbm)) return &R2g2q_sl_mmpm_SLC<3,0,1,2,T>;

	if ( pro == process(qbp,p,qm,p)) return &R2g2q_sl_ppmp_SLC<0,1,2,3,T>;
//	if ( pro == process(p,qbp,p,qm)) return &R2g2q_sl_ppmp_SLC<1,2,3,0,T>;
//	if ( pro == process(qm,p,qbp,p)) return &R2g2q_sl_ppmp_SLC<2,3,0,1,T>;
//	if ( pro == process(p,qm,p,qbp)) return &R2g2q_sl_ppmp_SLC<3,0,1,2,T>;

return 0;
}


#define _CASE_R2g2q_nf_Ptr_eval(K) case K : return &R2g2q_nf_ ## K

template <class T> complex<T> (*R2g2q_nf_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	switch (hc) {
	case 249 : return &R2g2q_nf_249;
	case 231 : return &R2g2q_nf_231;
	case 159 : return &R2g2q_nf_159;
	case 126 : return &R2g2q_nf_126;
	case   6 : return &R2g2q_nf_6;
	case  24 : return &R2g2q_nf_24;
	case  96 : return &R2g2q_nf_96;
	case 129 : return &R2g2q_nf_129;
	case 201 : return &R2g2q_nf_201;
	case  39 : return &R2g2q_nf_39;
	case 156 : return &R2g2q_nf_156;
	case 114 : return &R2g2q_nf_114;
	case 54  : return &R2g2q_nf_54;
	case 216 : return &R2g2q_nf_216;
	case 99  : return &R2g2q_nf_99;
	case 141 : return &R2g2q_nf_141;
	case 57  : return &R2g2q_nf_57;
	case 228 : return &R2g2q_nf_228;
	case 147 : return &R2g2q_nf_147;
	case 78  : return &R2g2q_nf_78;
	case 198 : return &R2g2q_nf_198;
	case 27  : return &R2g2q_nf_27;
	case 108 : return &R2g2q_nf_108;
	case 177 : return &R2g2q_nf_177;
	_CASE_R2g2q_nf_Ptr_eval(9);
	_CASE_R2g2q_nf_Ptr_eval(36);
	_CASE_R2g2q_nf_Ptr_eval(144);
	_CASE_R2g2q_nf_Ptr_eval(66);
	_CASE_R2g2q_nf_Ptr_eval(246);
	_CASE_R2g2q_nf_Ptr_eval(219);
	_CASE_R2g2q_nf_Ptr_eval(111);
	_CASE_R2g2q_nf_Ptr_eval(189);
	// vanishing R pieces:
	case 225: // qmmqpp
	case 135: // pqmmqp
	case 30:  // qppqmm
	case 120: // mqppqm
	case 210:// qpmqmp
	case 180:// mqmpqp
	case 75: //pqpmqm
	case 45: // qmpqpm

		//finite amplitudes
			 case 237 :// qm p qp p
			 case 183 : // p qm p qp
			 case 222 : // qp p qm p
			 case 123 : // p qp p qm

			 case 18: // qp m qm m
			 case 72: // m qp m qm
			 case 33: // qm m qp m
			 case 132: // m qm m qp

	          return &R2g2q0;
default: return &R2g2q0;
	}
}





template complex<R> (*R2g2q_L_Ptr_eval(const process& pro))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R2g2q_L_Ptr_eval(const process& pro))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R2g2q_L_Ptr_eval(const process& pro))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R2g2q_L_Ptr_eval(const process& pro))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

 template complex<R> (*R2g2q_SLC_Ptr_eval(const process& pro))(const eval_param<R>&,const mass_param_coll&) ;
  template complex<RHP> (*R2g2q_SLC_Ptr_eval(const process& pro))(const eval_param<RHP>&,const mass_param_coll&) ;
  template complex<RVHP> (*R2g2q_SLC_Ptr_eval(const process& pro))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

  template complex<RGMP> (*R2g2q_SLC_Ptr_eval(const process& pro))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

  template complex<R> (*R2g2q_nf_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
   template complex<RHP> (*R2g2q_nf_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
   template complex<RVHP> (*R2g2q_nf_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

   template complex<RGMP> (*R2g2q_nf_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

}

