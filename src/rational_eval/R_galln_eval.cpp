	// all-n expressions for rational pieces - glue
//
// last update by CFB - 10/7/08
// all-p and all-m amplitudes, all-but-one p and all-but-one m amplitudes
// from Bern, Dixon, Kosower, PRD 72:125003 (2005) - "the last of the Mohicans"

#include "R_alln_eval.h"

using namespace std;
namespace BH  {

// all-plus amplitude
template <class T> complex<T> Rallp(const eval_param<T>& ep,const mass_param_coll& mpc)
{
    complex<T> temp = complex<T>(0,0);
    size_t n = ep.size();

    for(size_t i1 = 0; i1 < n-3; i1++)
    	for(size_t i2 = i1+1; i2 < n-2; i2++)
    		for(size_t i3 = i2+1; i3 < n-1; i3++)
    			for (size_t i4 = i3+1; i4 < n; i4++)
    				temp += ep.spa(1,2)*ep.spb(2,3)* \
    				        ep.spa(3,4)*ep.spb(4,1);
    return (-temp/complex<T>(3,0)*denomang(ep,mpc)*complex<T>(0,1));

}

// all-minus amplitude
template <class T> complex<T> Rallm(const eval_param<T>& ep,const mass_param_coll& mpc)
{
    complex<T> temp = complex<T>(0,0);
    size_t n = ep.size();

    for(size_t i1 = 0; i1 < n-3; i1++)
    	for(size_t i2 = i1+1; i2 < n-2; i2++)
    		for(size_t i3 = i2+1; i3 < n-1; i3++)
    			for (size_t i4 = i3+1; i4 < n; i4++)
    				temp += ep.spb(1,2)*ep.spa(2,3)* \
    				        ep.spb(3,4)*ep.spa(4,1);
    return (-temp/complex<T>(3,0)*denomsqu(ep,mpc)*complex<T>(0,1));

}

//// mpppppppppppppppp
//// note: 1st = minus
//// BDK, PRD 72:125003 (2005), eqs. (5.44)-(5.46)
//template <class T> complex<T> Rallmp(const eval_param<T>& ep,const mass_param_coll& mpc)
//{
//    complex<T> T1 = complex<T>(0,0);
//    complex<T> T2 = complex<T>(0,0);
//
//    size_t n = ep.size();
//
//    // eq. (5.45) note mpcices -1 in C++
//    for (size_t l = 1; l <= n-2; l++ )
//    {
//    	vector<int> mpclp1n;
//    	for (size_t ll = l+1; ll <= n-1; ll++)
//    		mpclp1n.push_back(ll);
//        int Kllp1 = ep.Sum(l,l+1);
//        int Klp1n = ep.Sum(mpclp1n);
//    	T1 += ep.spa(0,l)*ep.spa(0,l+1)/ \
//    	      ep.spa(l,l+1)*ep.spaa(0,Kllp1,Klp1n,0);
//    }
//
//    // eq. (5.46)
//    for (size_t l2 = 2; l2 <= n-3; l2++)
//    {
//    	vector<int> mpc2lm1;
//    	for (size_t a1 = 1; a1 <= l2-1; a1++)
//    		mpc2lm1.push_back(mpc.at(a1));
//    	int K2lm1 = ep.Sum(mpc2lm1);
//    	for (size_t p = l2+1; p <= n-2; p++)
//    	{
//    		vector<int> mpclp;
//    		vector<int> mpcpp1n;
//    		for (size_t a2 = l2; a2 <= p; a2++)
//    			mpclp.push_back(mpc.at(a2));
//    		int Klp = ep.Sum(mpclp);
//    		for (size_t a3 = p+1; a3 <= n-1; a3++)
//    			mpcpp1n.push_back(mpc.at(a3));
//    		int Kpp1n = ep.Sum(mpcpp1n);
//
//    		complex<T> aux = complex<T>(0,0);
//    		for (size_t i = l2; i <= p-1; i++)
//    			for (size_t m = i+1; m <= p; m++)
//    				for (size_t it = l2; it <= p-1; it++)
//    					for (size_t mt = it+1; mt <= p; mt++)
//    						aux += ep.spaa(0,K2lm1,i,m)* \
//    						       ep.spba(m,t,mpc.at(mt),Kpp1n,0);
//
//
//    		T2 += ep.spa(mpc.at(l2-1),mpc.at(l2))*ep.spa(p,mpc.at(p+1))/ \
//    		      ep.m2(Klp)/ep.spaa(0,Kpp1n,Klp,mpc.at(l2-1))/ \
//    		      ep.spaa(0,Kpp1n,Klp,mpc.at(l2))/ \
//    		      ep.spaa(0,K2lm1,Klp,p)/ \
//    		      ep.spaa(0,K2lm1,Klp,mpc.at(p+1))*aux* \
//    		      pow(ep.spaa(0,Klp,Kpp1n,0),3);
//
//    	}
//    }
//
//    // eq. (5.44)
//    return ( (T1 + T2)/complex<T>(3,0)*denomang(ep,mpc)*complex<T>(0,1) );
//}
//
//// pmmmmmmmmmmmmmmmm
//// note - 1st = plus
//// parity conjugate from above
//template <class T> complex<T> Rallpm(const eval_param<T>& ep,const mass_param_coll& mpc)
//{
//    complex<T> T1 = complex<T>(0,0);
//    complex<T> T2 = complex<T>(0,0);
//
//    size_t n = ep.size();
//
//    // eq. (5.45)
//    for (size_t l = 1; l <= n-2; l++ )
//    {
//    	vector<int> mpclp1n;
//    	for (size_t ll = l+1; ll <= n-1; ll++)
//    		mpclp1n.push_back(mpc.at(ll));
//        int Kllp1 = ep.Sum(l,mpc.at(l+1));
//        int Klp1n = ep.Sum(mpclp1n);
//    	T1 += ep.spb(0,l)*ep.spb(0,mpc.at(l+1))/ \
//    	      ep.spb(l,mpc.at(l+1))*ep.spbb(0,Kllp1,Klp1n,0);
//    }
//
//    // eq. (5.46)
//    for (size_t l2 = 2; l2 <= n-3; l2++)
//    {
//    	vector<int> mpc2lm1;
//    	for (size_t a1 = 1; a1 <= l2-1; a1++)
//    		mpc2lm1.push_back(mpc.at(a1));
//    	int K2lm1 = ep.Sum(mpc2lm1);
//    	for (size_t p = l2+1; p <= n-2; p++)
//    	{
//    		vector<int> mpclp;
//    		vector<int> mpcpp1n;
//    		for (size_t a2 = l2; a2 <= p; a2++)
//    			mpclp.push_back(mpc.at(a2));
//    		int Klp = ep.Sum(mpclp);
//    		for (size_t a3 = p+1; a3 <= n-1; a3++)
//    			mpcpp1n.push_back(mpc.at(a3));
//    		int Kpp1n = ep.Sum(mpcpp1n);
//
//    		complex<T> aux = complex<T>(0,0);
//    		for (size_t i = l2; i <= p-1; i++)
//    			for (size_t m = i+1; m <= p; m++)
//    				for (size_t it = l2; it <= p-1; it++)
//    					for (size_t mt = it+1; mt <= p; mt++)
//    						aux -= ep.spbb(0,K2lm1,i,m)* \
//    						       ep.spab(m,t,mpc.at(mt),Kpp1n,0);
//
//    		T2 -= ep.spb(mpc.at(l2-1),mpc.at(l2))*ep.spb(p,mpc.at(p+1))/ \
//    		      ep.m2(Klp)/ep.spbb(0,Kpp1n,Klp,mpc.at(l2-1))/ \
//    		      ep.spbb(0,Kpp1n,Klp,mpc.at(l2))/ \
//    		      ep.spbb(0,K2lm1,Klp,p)/ \
//    		      ep.spbb(0,K2lm1,Klp,mpc.at(p+1))*aux* \
//    		      pow(ep.spbb(0,Klp,Kpp1n,0),3);
//
//    	}
//    }
//
//    // eq. (5.44)
//    return ( (T1 + T2)/complex<T>(3,0)*denomsqu(ep,mpc)*complex<T>(0,1) );
//
//}


// auxiliary stuff ---------------------

// angle bracket denominator
template <class T> complex<T> denomang(const eval_param<T>& ep,const mass_param_coll& mpc)
{
	complex<T> temp = complex<T>(1,0);
	for (size_t i = 0; i < ep.size()-1; i++)
		temp *= ep.spa(i,i+1);
	temp *= ep.spa(ep.size()-1,0);
	return (T(1)/temp);
}

// square bracket denominator
template <class T> complex<T> denomsqu(const eval_param<T>& ep,const mass_param_coll& mpc)
{
	complex<T> temp = complex<T>(1,0);
	for (size_t i = 0; i < ep.size()-1; i++)
		temp *= ep.spb(i+1,i);
	temp *= ep.spb(0,ep.size()-1);
	return (T(1)/temp);
}

// explicit instantiations
template complex<R> denomang(const eval_param<R>& ep, const mass_param_coll& mpc);
template complex<R> denomsqu(const eval_param<R>& ep, const mass_param_coll& mpc);
template complex<R> Rallp(const eval_param<R>& ep, const mass_param_coll& mpc);
template complex<R> Rallm(const eval_param<R>& ep, const mass_param_coll& mpc);
//template complex<R> Rallmp(const eval_param<R>& ep, const mass_param_coll& mpc);
//template complex<R> Rallpm(const eval_param<R>& ep, const mass_param_coll& mpc);


template complex<RHP> denomang(const eval_param<RHP>& ep, const mass_param_coll& mpc);
template complex<RHP> denomsqu(const eval_param<RHP>& ep, const mass_param_coll& mpc);
template complex<RHP> Rallp(const eval_param<RHP>& ep, const mass_param_coll& mpc);
template complex<RHP> Rallm(const eval_param<RHP>& ep, const mass_param_coll& mpc);
//template complex<RHP> Rallmp(const eval_param<RHP>& ep, const mass_param_coll& mpc);
//template complex<RHP> Rallpm(const eval_param<RHP>& ep, const mass_param_coll& mpc);


template complex<RVHP> denomang(const eval_param<RVHP>& ep, const mass_param_coll& mpc);
template complex<RVHP> denomsqu(const eval_param<RVHP>& ep, const mass_param_coll& mpc);
template complex<RVHP> Rallp(const eval_param<RVHP>& ep, const mass_param_coll& mpc);
template complex<RVHP> Rallm(const eval_param<RVHP>& ep, const mass_param_coll& mpc);
//template complex<RVHP> Rallmp(const eval_param<RVHP>& ep, const mass_param_coll& mpc);
//template complex<RVHP> Rallpm(const eval_param<RVHP>& ep, const mass_param_coll& mpc);

#if BH_USE_GMP
template complex<RGMP> denomang(const eval_param<RGMP>& ep, const mass_param_coll& mpc);
template complex<RGMP> denomsqu(const eval_param<RGMP>& ep, const mass_param_coll& mpc);
template complex<RGMP> Rallp(const eval_param<RGMP>& ep, const mass_param_coll& mpc);
template complex<RGMP> Rallm(const eval_param<RGMP>& ep, const mass_param_coll& mpc);
#endif

}

