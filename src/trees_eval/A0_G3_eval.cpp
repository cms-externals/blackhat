
/*
 * A0_G3_eval.cpp
 *
 * created with Harald's scripts from trees computed by Lance and Yael
 * notes:
 * this file is used only for BHsetting: USE_G3_COUPLING yes
 * *) include Tr(G^3) coupling at tree level as in hep-ph/9312363 
 *
 */

#include <complex>
#include <vector>
#include "constants.h"
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"
#include "BH_debug.h"

using namespace std;

#define SPA(i,j) ep.spa(i,j)
#define SPB(i,j) ep.spb(i,j)
#define S(i,j) ep.s(i,j)

namespace BH {
template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_mmmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmmm");
	return((complex<T>(0,6)*S(i1,i2)*S(i1,i3)*S(i2,i3))/(Lam2*SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_mmmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmmp");
	return((complex<T>(0,-3)*SPA(i1,i2)*SPA(i2,i3)*pow(SPA(i3,i1),2))/(Lam2*SPA(i3,i4)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_mmpm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmpm");
	return((complex<T>(0,-3)*SPA(i1,i2)*pow(SPA(i2,i4),2)*SPA(i4,i1))/(Lam2*SPA(i2,i3)*SPA(i3,i4)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_mmpp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmpp");
	return((complex<T>(0,1)*pow(SPA(i1,i2),3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_mpmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mpmm");
	return((complex<T>(0,-3)*pow(SPA(i1,i3),2)*SPA(i3,i4)*SPA(i4,i1))/(Lam2*SPA(i1,i2)*SPA(i2,i3)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_mpmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mpmp");
	return((complex<T>(0,1)*pow(SPA(i1,i3),4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_mppm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mppm");
	return((complex<T>(0,1)*pow(SPB(i2,i3),3))/(SPB(i1,i2)*SPB(i3,i4)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_mppp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mppp");
	return((complex<T>(0,-3)*SPB(i2,i3)*SPB(i3,i4)*pow(SPB(i4,i2),2))/(Lam2*SPB(i1,i2)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_pmmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmmm");
	return((complex<T>(0,-3)*SPA(i2,i3)*SPA(i3,i4)*pow(SPA(i4,i2),2))/(Lam2*SPA(i1,i2)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_pmmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmmp");
	return((complex<T>(0,1)*pow(SPA(i2,i3),3))/(SPA(i1,i2)*SPA(i3,i4)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_pmpm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmpm");
	return((complex<T>(0,1)*pow(SPB(i1,i3),4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_pmpp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmpp");
	return((complex<T>(0,-3)*pow(SPB(i1,i3),2)*SPB(i3,i4)*SPB(i4,i1))/(Lam2*SPB(i1,i2)*SPB(i2,i3)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_ppmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","ppmm");
	return((complex<T>(0,1)*pow(SPB(i1,i2),3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_ppmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","ppmp");
	return((complex<T>(0,-3)*SPB(i1,i2)*pow(SPB(i2,i4),2)*SPB(i4,i1))/(Lam2*SPB(i2,i3)*SPB(i3,i4)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_pppm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pppm");
	return((complex<T>(0,-3)*SPB(i1,i2)*SPB(i2,i3)*pow(SPB(i3,i1),2))/(Lam2*SPB(i3,i4)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A4g_pppp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pppp");
	return((complex<T>(0,6)*S(i1,i2)*S(i1,i3)*S(i2,i3))/(Lam2*SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i1)));
}


template <class T> complex<T> (*A4g_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
        switch (hc) {
    case	0:	return &A4g_mmmm_G3_eval<0,1,2,3>;
    case	8:	return &A4g_mmmp_G3_eval<0,1,2,3>;
    case	4:	return &A4g_mmpm_G3_eval<0,1,2,3>;
    case	12:	return &A4g_mmpp_G3_eval<0,1,2,3>;
    case	2:	return &A4g_mpmm_G3_eval<0,1,2,3>;
    case	10:	return &A4g_mpmp_G3_eval<0,1,2,3>;
    case	6:	return &A4g_mppm_G3_eval<0,1,2,3>;
    case	14:	return &A4g_mppp_G3_eval<0,1,2,3>;
    case	1:	return &A4g_pmmm_G3_eval<0,1,2,3>;
    case	9:	return &A4g_pmmp_G3_eval<0,1,2,3>;
    case	5:	return &A4g_pmpm_G3_eval<0,1,2,3>;
    case	13:	return &A4g_pmpp_G3_eval<0,1,2,3>;
    case	3:	return &A4g_ppmm_G3_eval<0,1,2,3>;
    case	11:	return &A4g_ppmp_G3_eval<0,1,2,3>;
    case	7:	return &A4g_pppm_G3_eval<0,1,2,3>;
    case	15:	return &A4g_pppp_G3_eval<0,1,2,3>;
  default: // We return the zero pointer for all other helicity combinations
                //cout << " A4g_G3_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
                return 0;
        }
}

template complex<R> (*A4g_G3_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A4g_G3_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A4g_G3_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);
#if BH_USE_GMP
template complex<RGMP> (*A4g_G3_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif


/*
*
*
* 2q2g amplitudes
*
*
*/


template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2g_qmmmqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmmmqbp");
	return((complex<T>(0,3)*SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i1)*SPA(i4,i3))/(Lam2*SPA(i3,i4)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2g_qmmpqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmmpqbp");
	return((complex<T>(0,-1)*SPB(i1,i3)*pow(SPB(i4,i3),3))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2g_qmpmqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmpmqbp");
	return((complex<T>(0,-1)*pow(SPB(i4,i2),3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2g_qmppqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmppqbp");
	return((complex<T>(0,-3)*SPB(i2,i3)*SPB(i2,i4)*pow(SPB(i4,i3),2))/(Lam2*SPB(i3,i4)*SPB(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2g_qpmmqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qpmmqbm");
	return((complex<T>(0,3)*SPA(i2,i3)*SPA(i2,i4)*pow(SPA(i4,i3),2))/(Lam2*SPA(i3,i4)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2g_qpmpqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qpmpqbm");
	return((complex<T>(0,1)*pow(SPA(i4,i2),3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2g_qppmqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qppmqbm");
	return((complex<T>(0,1)*SPA(i1,i3)*pow(SPA(i4,i3),3))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i1)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2g_qpppqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qpppqbm");
	return((complex<T>(0,-3)*SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i1)*SPB(i4,i3))/(Lam2*SPB(i3,i4)*SPB(i4,i1)));
}


template <class T> complex<T> (*A2q2g_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
        switch (hc) {
    case	129:	return &A2q2g_qmmmqbp_G3_eval<0,1,2,3>;
    case	177:	return &A2q2g_qmmpqbp_G3_eval<0,1,2,3>;
    case	141:	return &A2q2g_qmpmqbp_G3_eval<0,1,2,3>;
    case	189:	return &A2q2g_qmppqbp_G3_eval<0,1,2,3>;
    case	66:	return &A2q2g_qpmmqbm_G3_eval<0,1,2,3>;
    case	114:	return &A2q2g_qpmpqbm_G3_eval<0,1,2,3>;
    case	78:	return &A2q2g_qppmqbm_G3_eval<0,1,2,3>;
    case	126:	return &A2q2g_qpppqbm_G3_eval<0,1,2,3>;
  default: // We return the zero pointer for all other helicity combinations
                //cout << " A2q2g_G3_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
                return 0;
        }
}

template complex<R> (*A2q2g_G3_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q2g_G3_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q2g_G3_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);
#if BH_USE_GMP
template complex<RGMP> (*A2q2g_G3_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

/*
*
*
* 4q amplitudes
*
*
*/


template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2Q_q1mq2bmq2pq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bmq2pq1bp");
	return((complex<T>(0,1)*pow(SPA(i1,i2),2))/(SPA(i1,i4)*SPA(i2,i3)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2Q_q1mq2bpq2mq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bpq2mq1bp");
	return((complex<T>(0,-1)*pow(SPA(i3,i1),2))/(SPA(i1,i4)*SPA(i3,i2)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2Q_q1pq2bmq2pq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bmq2pq1bm");
	return((complex<T>(0,1)*pow(SPA(i4,i2),2))/(SPA(i1,i4)*SPA(i2,i3)));
}

template < int i1, int i2, int i3, int i4, class T> complex<T>  A2q2Q_q1pq2bpq2mq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bpq2mq1bm");
	return((complex<T>(0,-1)*pow(SPA(i3,i4),2))/(SPA(i1,i4)*SPA(i3,i2)));
}


template <class T> complex<T> (*A2q2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
        switch (hc) {
    case	637:	return &A2q2Q_q1mq2bmq2pq1bp_G3_eval<0,1,2,3>;
    case	607:	return &A2q2Q_q1mq2bpq2mq1bp_G3_eval<0,1,2,3>;
    case	422:	return &A2q2Q_q1pq2bmq2pq1bm_G3_eval<0,1,2,3>;
    case	392:	return &A2q2Q_q1pq2bpq2mq1bm_G3_eval<0,1,2,3>;
  default: // We return the zero pointer for all other helicity combinations
                //cout << " A2q2Q_G3_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
                return 0;
        }
}

template complex<R> (*A2q2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> (*A2q2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif
/*
*
*
* 5g amplitudes
*
*
*/


template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mmmmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmmmm");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mmmmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmmmp");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mmmpm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmmpm");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mmmpp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmmpp");
	return((complex<T>(0,-3)*SPA(i1,i2)*SPA(i2,i3)*pow(SPA(i3,i1),2))/(Lam2*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i4,i5),3))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mmpmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmpmm");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mmpmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmpmp");
	return((complex<T>(0,-3)*SPA(i1,i2)*pow(SPA(i2,i4),2)*pow(SPA(i4,i1),2))/(Lam2*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i3,i5),4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mmppm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmppm");
	return((complex<T>(0,-3)*SPA(i1,i2)*pow(SPA(i2,i5),2)*SPA(i5,i1))/(Lam2*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)) - (complex<T>(0,1)*pow(SPB(i3,i4),3))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mmppp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mmppp");
	return((complex<T>(0,1)*pow(SPA(i1,i2),3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i3,i4)*SPB(i4,i5)*pow(SPB(i5,i3),2))/(Lam2*SPB(i1,i2)*SPB(i2,i3)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mpmmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mpmmm");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mpmmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mpmmp");
	return((complex<T>(0,-3)*pow(SPA(i1,i3),2)*SPA(i3,i4)*pow(SPA(i4,i1),2))/(Lam2*SPA(i1,i2)*SPA(i2,i3)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i2,i5),4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mpmpm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mpmpm");
	return((complex<T>(0,-3)*pow(SPA(i1,i3),2)*pow(SPA(i3,i5),2)*SPA(i5,i1))/(Lam2*SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)) - (complex<T>(0,1)*pow(SPB(i2,i4),4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mpmpp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mpmpp");
	return((complex<T>(0,1)*pow(SPA(i1,i3),4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*pow(SPB(i2,i4),2)*SPB(i4,i5)*pow(SPB(i5,i2),2))/(Lam2*SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mppmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mppmm");
	return((complex<T>(0,-3)*pow(SPA(i1,i4),2)*SPA(i4,i5)*SPA(i5,i1))/(Lam2*SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)) - (complex<T>(0,1)*pow(SPB(i2,i3),3))/(SPB(i1,i2)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mppmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mppmp");
	return((complex<T>(0,1)*pow(SPA(i1,i4),4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i2,i3)*pow(SPB(i3,i5),2)*pow(SPB(i5,i2),2))/(Lam2*SPB(i1,i2)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mpppm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mpppm");
	return((complex<T>(0,1)*pow(SPA(i1,i5),4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i2,i3)*SPB(i3,i4)*pow(SPB(i4,i2),2))/(Lam2*SPB(i1,i2)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_mpppp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","mpppp");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pmmmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmmmm");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pmmmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmmmp");
	return((complex<T>(0,-3)*SPA(i2,i3)*SPA(i3,i4)*pow(SPA(i4,i2),2))/(Lam2*SPA(i1,i2)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i1,i5),4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pmmpm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmmpm");
	return((complex<T>(0,-3)*SPA(i2,i3)*pow(SPA(i3,i5),2)*pow(SPA(i5,i2),2))/(Lam2*SPA(i1,i2)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i1,i4),4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pmmpp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmmpp");
	return((complex<T>(0,1)*pow(SPA(i2,i3),3))/(SPA(i1,i2)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*pow(SPB(i1,i4),2)*SPB(i4,i5)*SPB(i5,i1))/(Lam2*SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pmpmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmpmm");
	return((complex<T>(0,-3)*pow(SPA(i2,i4),2)*SPA(i4,i5)*pow(SPA(i5,i2),2))/(Lam2*SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i1,i3),4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pmpmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmpmp");
	return((complex<T>(0,1)*pow(SPA(i2,i4),4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*pow(SPB(i1,i3),2)*pow(SPB(i3,i5),2)*SPB(i5,i1))/(Lam2*SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pmppm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmppm");
	return((complex<T>(0,1)*pow(SPA(i2,i5),4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*pow(SPB(i1,i3),2)*SPB(i3,i4)*pow(SPB(i4,i1),2))/(Lam2*SPB(i1,i2)*SPB(i2,i3)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pmppp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pmppp");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_ppmmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","ppmmm");
	return((complex<T>(0,-3)*SPA(i3,i4)*SPA(i4,i5)*pow(SPA(i5,i3),2))/(Lam2*SPA(i1,i2)*SPA(i2,i3)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i1,i2),3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_ppmmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","ppmmp");
	return((complex<T>(0,1)*pow(SPA(i3,i4),3))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i1,i2)*pow(SPB(i2,i5),2)*SPB(i5,i1))/(Lam2*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_ppmpm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","ppmpm");
	return((complex<T>(0,1)*pow(SPA(i3,i5),4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i1,i2)*pow(SPB(i2,i4),2)*pow(SPB(i4,i1),2))/(Lam2*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_ppmpp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","ppmpp");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pppmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pppmm");
	return((complex<T>(0,1)*pow(SPA(i4,i5),3))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i1,i2)*SPB(i2,i3)*pow(SPB(i3,i1),2))/(Lam2*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_pppmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","pppmp");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_ppppm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","ppppm");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A5g_ppppp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","ppppp");
	return(complex<T>(0,0));
}


template <class T> complex<T> (*A5g_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
        switch (hc) {
    case	0:	return &A5g_mmmmm_G3_eval<0,1,2,3,4>;
    case	16:	return &A5g_mmmmp_G3_eval<0,1,2,3,4>;
    case	8:	return &A5g_mmmpm_G3_eval<0,1,2,3,4>;
    case	24:	return &A5g_mmmpp_G3_eval<0,1,2,3,4>;
    case	4:	return &A5g_mmpmm_G3_eval<0,1,2,3,4>;
    case	20:	return &A5g_mmpmp_G3_eval<0,1,2,3,4>;
    case	12:	return &A5g_mmppm_G3_eval<0,1,2,3,4>;
    case	28:	return &A5g_mmppp_G3_eval<0,1,2,3,4>;
    case	2:	return &A5g_mpmmm_G3_eval<0,1,2,3,4>;
    case	18:	return &A5g_mpmmp_G3_eval<0,1,2,3,4>;
    case	10:	return &A5g_mpmpm_G3_eval<0,1,2,3,4>;
    case	26:	return &A5g_mpmpp_G3_eval<0,1,2,3,4>;
    case	6:	return &A5g_mppmm_G3_eval<0,1,2,3,4>;
    case	22:	return &A5g_mppmp_G3_eval<0,1,2,3,4>;
    case	14:	return &A5g_mpppm_G3_eval<0,1,2,3,4>;
    case	30:	return &A5g_mpppp_G3_eval<0,1,2,3,4>;
    case	1:	return &A5g_pmmmm_G3_eval<0,1,2,3,4>;
    case	17:	return &A5g_pmmmp_G3_eval<0,1,2,3,4>;
    case	9:	return &A5g_pmmpm_G3_eval<0,1,2,3,4>;
    case	25:	return &A5g_pmmpp_G3_eval<0,1,2,3,4>;
    case	5:	return &A5g_pmpmm_G3_eval<0,1,2,3,4>;
    case	21:	return &A5g_pmpmp_G3_eval<0,1,2,3,4>;
    case	13:	return &A5g_pmppm_G3_eval<0,1,2,3,4>;
    case	29:	return &A5g_pmppp_G3_eval<0,1,2,3,4>;
    case	3:	return &A5g_ppmmm_G3_eval<0,1,2,3,4>;
    case	19:	return &A5g_ppmmp_G3_eval<0,1,2,3,4>;
    case	11:	return &A5g_ppmpm_G3_eval<0,1,2,3,4>;
    case	27:	return &A5g_ppmpp_G3_eval<0,1,2,3,4>;
    case	7:	return &A5g_pppmm_G3_eval<0,1,2,3,4>;
    case	23:	return &A5g_pppmp_G3_eval<0,1,2,3,4>;
    case	15:	return &A5g_ppppm_G3_eval<0,1,2,3,4>;
    case	31:	return &A5g_ppppp_G3_eval<0,1,2,3,4>;
  default: // We return the zero pointer for all other helicity combinations
                //cout << " A5g_G3_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
                return 0;
        }
}

template complex<R> (*A5g_G3_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A5g_G3_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A5g_G3_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);
#if BH_USE_GMP
template complex<RGMP> (*A5g_G3_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

/*
*
*
* 2q3g amplitudes
*
*
*/


template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qmmmmqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmmmmqbp");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qmmmpqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmmmpqbp");
	return((complex<T>(0,3)*SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i1)*SPA(i5,i3))/(Lam2*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,1)*SPB(i1,i4)*pow(SPB(i5,i4),3))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qmmpmqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmmpmqbp");
	return((complex<T>(0,3)*SPA(i1,i2)*pow(SPA(i2,i4),2)*SPA(i4,i1)*SPA(i5,i4))/(Lam2*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,1)*SPB(i1,i3)*pow(SPB(i5,i3),3))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qmmppqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmmppqbp");
	return((complex<T>(0,1)*pow(SPA(i1,i2),2)*SPA(i5,i2))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i1,i3)*SPB(i3,i5)*pow(SPB(i4,i3),2)*pow(SPB(i5,i4),2))/(Lam2*SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qmpmmqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmpmmqbp");
	return((complex<T>(0,3)*pow(SPA(i1,i3),2)*SPA(i3,i4)*SPA(i4,i1)*SPA(i5,i4))/(Lam2*SPA(i1,i2)*SPA(i2,i3)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,1)*pow(SPB(i5,i2),3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qmpmpqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmpmpqbp");
	return((complex<T>(0,1)*pow(SPA(i1,i3),3)*SPA(i5,i3))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i2,i5)*pow(SPB(i4,i2),2)*pow(SPB(i5,i4),2))/(Lam2*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qmppmqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmppmqbp");
	return((complex<T>(0,1)*pow(SPA(i1,i4),3)*SPA(i5,i4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i2,i5)*pow(SPB(i3,i2),2)*pow(SPB(i5,i3),2))/(Lam2*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qmpppqbp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qmpppqbp");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qpmmmqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qpmmmqbm");
	return(complex<T>(0,0));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qpmmpqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qpmmpqbm");
	return((complex<T>(0,3)*SPA(i2,i5)*pow(SPA(i3,i2),2)*pow(SPA(i5,i3),2))/(Lam2*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,1)*pow(SPB(i1,i4),3)*SPB(i5,i4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qpmpmqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qpmpmqbm");
	return((complex<T>(0,3)*SPA(i2,i5)*pow(SPA(i4,i2),2)*pow(SPA(i5,i4),2))/(Lam2*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,1)*pow(SPB(i1,i3),3)*SPB(i5,i3))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qpmppqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qpmppqbm");
	return((complex<T>(0,1)*pow(SPA(i5,i2),3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*pow(SPB(i1,i3),2)*SPB(i3,i4)*SPB(i4,i1)*SPB(i5,i4))/(Lam2*SPB(i1,i2)*SPB(i2,i3)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qppmmqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qppmmqbm");
	return((complex<T>(0,3)*SPA(i1,i3)*SPA(i3,i5)*pow(SPA(i4,i3),2)*pow(SPA(i5,i4),2))/(Lam2*SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,1)*pow(SPB(i1,i2),2)*SPB(i5,i2))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qppmpqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qppmpqbm");
	return((complex<T>(0,1)*SPA(i1,i3)*pow(SPA(i5,i3),3))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i1,i2)*pow(SPB(i2,i4),2)*SPB(i4,i1)*SPB(i5,i4))/(Lam2*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qpppmqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qpppmqbm");
	return((complex<T>(0,1)*SPA(i1,i4)*pow(SPA(i5,i4),3))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i1)*SPB(i5,i3))/(Lam2*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q3g_qppppqbm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","qppppqbm");
	return(complex<T>(0,0));
}


template <class T> complex<T> (*A2q3g_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
        switch (hc) {
    case	513:	return &A2q3g_qmmmmqbp_G3_eval<0,1,2,3,4>;
    case	705:	return &A2q3g_qmmmpqbp_G3_eval<0,1,2,3,4>;
    case	561:	return &A2q3g_qmmpmqbp_G3_eval<0,1,2,3,4>;
    case	753:	return &A2q3g_qmmppqbp_G3_eval<0,1,2,3,4>;
    case	525:	return &A2q3g_qmpmmqbp_G3_eval<0,1,2,3,4>;
    case	717:	return &A2q3g_qmpmpqbp_G3_eval<0,1,2,3,4>;
    case	573:	return &A2q3g_qmppmqbp_G3_eval<0,1,2,3,4>;
    case	765:	return &A2q3g_qmpppqbp_G3_eval<0,1,2,3,4>;
    case	258:	return &A2q3g_qpmmmqbm_G3_eval<0,1,2,3,4>;
    case	450:	return &A2q3g_qpmmpqbm_G3_eval<0,1,2,3,4>;
    case	306:	return &A2q3g_qpmpmqbm_G3_eval<0,1,2,3,4>;
    case	498:	return &A2q3g_qpmppqbm_G3_eval<0,1,2,3,4>;
    case	270:	return &A2q3g_qppmmqbm_G3_eval<0,1,2,3,4>;
    case	462:	return &A2q3g_qppmpqbm_G3_eval<0,1,2,3,4>;
    case	318:	return &A2q3g_qpppmqbm_G3_eval<0,1,2,3,4>;
    case	510:	return &A2q3g_qppppqbm_G3_eval<0,1,2,3,4>;
  default: // We return the zero pointer for all other helicity combinations
                //cout << " A2q3g_G3_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
                return 0;
        }
}

template complex<R> (*A2q3g_G3_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q3g_G3_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q3g_G3_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);
#if BH_USE_GMP
template complex<RGMP> (*A2q3g_G3_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

/*
*
*
* 4q1g amplitudes
*
*
*/


template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mmq2bmq2pq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mmq2bmq2pq1bp");
	return((complex<T>(0,3)*SPA(i1,i2)*SPA(i1,i3)*SPA(i3,i2))/(Lam2*SPA(i1,i5)*SPA(i3,i4)) - (complex<T>(0,1)*SPB(i1,i3)*pow(SPB(i5,i4),3))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mmq2bpq2mq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mmq2bpq2mq1bp");
	return((complex<T>(0,-3)*SPA(i1,i2)*SPA(i1,i4)*SPA(i4,i2))/(Lam2*SPA(i1,i5)*SPA(i4,i3)) - (complex<T>(0,1)*SPB(i1,i3)*pow(SPB(i5,i3),2)*SPB(i5,i4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mpq2bmq2pq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mpq2bmq2pq1bp");
	return((complex<T>(0,1)*pow(SPA(i1,i3),3)*SPA(i5,i4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,3)*SPB(i4,i2)*SPB(i5,i2)*SPB(i5,i4))/(Lam2*SPB(i1,i5)*SPB(i3,i4)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mpq2bpq2mq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mpq2bpq2mq1bp");
	return((complex<T>(0,1)*SPA(i1,i3)*pow(SPA(i1,i4),2)*SPA(i5,i4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i3,i2)*SPB(i5,i2)*SPB(i5,i3))/(Lam2*SPB(i1,i5)*SPB(i4,i3)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bmmq2pq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bmmq2pq1bp");
	return((complex<T>(0,-1)*pow(SPB(i5,i4),3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bmpq2pq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bmpq2pq1bp");
	return((complex<T>(0,1)*pow(SPA(i1,i2),2)*SPA(i5,i4))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bmq2pmq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bmq2pmq1bp");
	return((complex<T>(0,-3)*SPA(i1,i2)*SPA(i1,i4)*SPA(i2,i4))/(Lam2*SPA(i1,i5)*SPA(i2,i3)) - (complex<T>(0,1)*pow(SPB(i5,i3),3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bmq2ppq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bmq2ppq1bp");
	return((complex<T>(0,1)*pow(SPA(i1,i2),2)*SPA(i5,i3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i3,i4)*SPB(i5,i3)*SPB(i5,i4))/(Lam2*SPB(i1,i5)*SPB(i2,i3)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bmq2pq1bpm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bmq2pq1bpm");
	return((complex<T>(0,-1)*pow(SPB(i4,i3),3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bmq2pq1bpp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bmq2pq1bpp");
	return((complex<T>(0,1)*pow(SPA(i1,i2),2)*SPA(i4,i3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bpmq2mq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bpmq2mq1bp");
	return((complex<T>(0,-1)*pow(SPB(i5,i2),2)*SPB(i5,i4))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bppq2mq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bppq2mq1bp");
	return((complex<T>(0,1)*pow(SPA(i1,i4),2)*SPA(i5,i4))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bpq2mmq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bpq2mmq1bp");
	return((complex<T>(0,3)*SPA(i1,i3)*SPA(i1,i4)*SPA(i3,i4))/(Lam2*SPA(i1,i5)*SPA(i3,i2)) - (complex<T>(0,1)*pow(SPB(i5,i2),2)*SPB(i5,i3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bpq2mpq1bp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bpq2mpq1bp");
	return((complex<T>(0,1)*pow(SPA(i1,i3),2)*SPA(i5,i3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,3)*SPB(i2,i4)*SPB(i5,i2)*SPB(i5,i4))/(Lam2*SPB(i1,i5)*SPB(i3,i2)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bpq2mq1bpm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bpq2mq1bpm");
	return((complex<T>(0,-1)*pow(SPB(i4,i2),2)*SPB(i4,i3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1mq2bpq2mq1bpp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1mq2bpq2mq1bpp");
	return((complex<T>(0,1)*pow(SPA(i1,i3),2)*SPA(i4,i3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pmq2bmq2pq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pmq2bmq2pq1bm");
	return((complex<T>(0,-3)*SPA(i3,i2)*SPA(i5,i2)*SPA(i5,i3))/(Lam2*SPA(i3,i4)*SPA(i5,i1)) - (complex<T>(0,1)*SPB(i1,i3)*pow(SPB(i1,i4),2)*SPB(i5,i4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pmq2bpq2mq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pmq2bpq2mq1bm");
	return((complex<T>(0,3)*SPA(i4,i2)*SPA(i5,i2)*SPA(i5,i4))/(Lam2*SPA(i4,i3)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i1,i3),3)*SPB(i5,i4))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1ppq2bmq2pq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1ppq2bmq2pq1bm");
	return((complex<T>(0,1)*SPA(i1,i3)*pow(SPA(i5,i3),2)*SPA(i5,i4))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i1,i2)*SPB(i1,i4)*SPB(i4,i2))/(Lam2*SPB(i3,i4)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1ppq2bpq2mq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1ppq2bpq2mq1bm");
	return((complex<T>(0,1)*SPA(i1,i3)*pow(SPA(i5,i4),3))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,3)*SPB(i1,i2)*SPB(i1,i3)*SPB(i3,i2))/(Lam2*SPB(i4,i3)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bmmq2pq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bmmq2pq1bm");
	return((complex<T>(0,-1)*pow(SPB(i1,i4),2)*SPB(i5,i4))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bmpq2pq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bmpq2pq1bm");
	return((complex<T>(0,1)*pow(SPA(i5,i2),2)*SPA(i5,i4))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bmq2pmq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bmq2pmq1bm");
	return((complex<T>(0,3)*SPA(i2,i4)*SPA(i5,i2)*SPA(i5,i4))/(Lam2*SPA(i2,i3)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i1,i3),2)*SPB(i5,i3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bmq2ppq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bmq2ppq1bm");
	return((complex<T>(0,1)*pow(SPA(i5,i2),2)*SPA(i5,i3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) - (complex<T>(0,3)*SPB(i1,i3)*SPB(i1,i4)*SPB(i3,i4))/(Lam2*SPB(i2,i3)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bmq2pq1bmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bmq2pq1bmm");
	return((complex<T>(0,-1)*pow(SPB(i1,i3),2)*SPB(i4,i3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bmq2pq1bmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bmq2pq1bmp");
	return((complex<T>(0,1)*pow(SPA(i4,i2),2)*SPA(i4,i3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bpmq2mq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bpmq2mq1bm");
	return((complex<T>(0,-1)*pow(SPB(i1,i2),2)*SPB(i5,i4))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bppq2mq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bppq2mq1bm");
	return((complex<T>(0,1)*pow(SPA(i5,i4),3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bpq2mmq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bpq2mmq1bm");
	return((complex<T>(0,-3)*SPA(i3,i4)*SPA(i5,i3)*SPA(i5,i4))/(Lam2*SPA(i3,i2)*SPA(i5,i1)) - (complex<T>(0,1)*pow(SPB(i1,i2),2)*SPB(i5,i3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bpq2mpq1bm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bpq2mpq1bm");
	return((complex<T>(0,1)*pow(SPA(i5,i3),3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)) + (complex<T>(0,3)*SPB(i1,i2)*SPB(i1,i4)*SPB(i2,i4))/(Lam2*SPB(i3,i2)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bpq2mq1bmm_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bpq2mq1bmm");
	return((complex<T>(0,-1)*pow(SPB(i1,i2),2)*SPB(i4,i3))/(SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i1)));
}

template < int i1, int i2, int i3, int i4, int i5, class T> complex<T>  A2q1g2Q_q1pq2bpq2mq1bmp_G3_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
	complex<T> Lam2 = pow(complex<T>(constants::s_GeV,0),2)*complex<T>(constants::G3_Lambda2,0);
	BH_DEBUG_MESSAGE2("new tree: ","q1pq2bpq2mq1bmp");
	return((complex<T>(0,1)*pow(SPA(i4,i3),3))/(SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i1)));
}


template <class T> complex<T> (*A2q1g2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
        switch (hc) {
    case	3817:	return &A2q1g2Q_q1mmq2bmq2pq1bp_G3_eval<0,1,2,3,4>;
    case	3637:	return &A2q1g2Q_q1mmq2bpq2mq1bp_G3_eval<0,1,2,3,4>;
    case	3835:	return &A2q1g2Q_q1mpq2bmq2pq1bp_G3_eval<0,1,2,3,4>;
    case	3655:	return &A2q1g2Q_q1mpq2bpq2mq1bp_G3_eval<0,1,2,3,4>;
    case	3697:	return &A2q1g2Q_q1mq2bmmq2pq1bp_G3_eval<0,1,2,3,4>;
    case	3805:	return &A2q1g2Q_q1mq2bmpq2pq1bp_G3_eval<0,1,2,3,4>;
    case	2797:	return &A2q1g2Q_q1mq2bmq2pmq1bp_G3_eval<0,1,2,3,4>;
    case	3445:	return &A2q1g2Q_q1mq2bmq2ppq1bp_G3_eval<0,1,2,3,4>;
    case	637:	return &A2q1g2Q_q1mq2bmq2pq1bpm_G3_eval<0,1,2,3,4>;
    case	4525:	return &A2q1g2Q_q1mq2bmq2pq1bpp_G3_eval<0,1,2,3,4>;
    case	3487:	return &A2q1g2Q_q1mq2bpmq2mq1bp_G3_eval<0,1,2,3,4>;
    case	3595:	return &A2q1g2Q_q1mq2bppq2mq1bp_G3_eval<0,1,2,3,4>;
    case	2767:	return &A2q1g2Q_q1mq2bpq2mmq1bp_G3_eval<0,1,2,3,4>;
    case	3415:	return &A2q1g2Q_q1mq2bpq2mpq1bp_G3_eval<0,1,2,3,4>;
    case	607:	return &A2q1g2Q_q1mq2bpq2mq1bpm_G3_eval<0,1,2,3,4>;
    case	4495:	return &A2q1g2Q_q1mq2bpq2mq1bpp_G3_eval<0,1,2,3,4>;
    case	2522:	return &A2q1g2Q_q1pmq2bmq2pq1bm_G3_eval<0,1,2,3,4>;
    case	2342:	return &A2q1g2Q_q1pmq2bpq2mq1bm_G3_eval<0,1,2,3,4>;
    case	2540:	return &A2q1g2Q_q1ppq2bmq2pq1bm_G3_eval<0,1,2,3,4>;
    case	2360:	return &A2q1g2Q_q1ppq2bpq2mq1bm_G3_eval<0,1,2,3,4>;
    case	2402:	return &A2q1g2Q_q1pq2bmmq2pq1bm_G3_eval<0,1,2,3,4>;
    case	2510:	return &A2q1g2Q_q1pq2bmpq2pq1bm_G3_eval<0,1,2,3,4>;
    case	1502:	return &A2q1g2Q_q1pq2bmq2pmq1bm_G3_eval<0,1,2,3,4>;
    case	2150:	return &A2q1g2Q_q1pq2bmq2ppq1bm_G3_eval<0,1,2,3,4>;
    case	422:	return &A2q1g2Q_q1pq2bmq2pq1bmm_G3_eval<0,1,2,3,4>;
    case	4310:	return &A2q1g2Q_q1pq2bmq2pq1bmp_G3_eval<0,1,2,3,4>;
    case	2192:	return &A2q1g2Q_q1pq2bpmq2mq1bm_G3_eval<0,1,2,3,4>;
    case	2300:	return &A2q1g2Q_q1pq2bppq2mq1bm_G3_eval<0,1,2,3,4>;
    case	1472:	return &A2q1g2Q_q1pq2bpq2mmq1bm_G3_eval<0,1,2,3,4>;
    case	2120:	return &A2q1g2Q_q1pq2bpq2mpq1bm_G3_eval<0,1,2,3,4>;
    case	392:	return &A2q1g2Q_q1pq2bpq2mq1bmm_G3_eval<0,1,2,3,4>;
    case	4280:	return &A2q1g2Q_q1pq2bpq2mq1bmp_G3_eval<0,1,2,3,4>;
  default: // We return the zero pointer for all other helicity combinations
                //cout << " A2q1g2Q_G3_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
                return 0;
        }
}

template complex<R> (*A2q1g2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q1g2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q1g2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);
#if BH_USE_GMP
template complex<RGMP> (*A2q1g2Q_G3_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}
