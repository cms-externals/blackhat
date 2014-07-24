/*!\file splitamp_tree.cpp
\brief Implementation of tree splitting amplitudes
       \n as listed, e.g. in Bern, Dixon, Dunbar, Kosower,
                             Nucl. Phys. B 425, 217 (1994)
*/

// by CFB, last revision 12/21/2008
// 6/23/08 - update to conform with new particle ID
// 12/21/08 -- added gluinos


#include "factorization.h"
#include "BH_utilities.h"
#include "qd_suppl.h"
#include "process.h"

#if BH_USE_GMP
#include "gmp_r.h"
#endif

namespace BH {
template <class T> complex<T> Split0(const process& p, momentum_configuration<T>& ps,
		int inda, int indb){
    if (p.p(1).is_a(lepton) || p.p(2).is_a(lepton) || p.p(3).is_a(lepton) )
	    return complex<T>(0,0);
	if (p.p(1).is_a(gluon) && p.p(2).is_a(gluon) && p.p(3).is_a(gluon) )
		return Sggg(p, ps, inda, indb);
	if (p.p(1).is_a(quark) && p.p(2).is_a(quark) && p.p(3).is_a(gluon) )
		return Sqqg(p, ps, inda, indb);
	if (p.p(1).is_a(gluino) && p.p(2).is_a(gluino) && p.p(3).is_a(gluon) )
		return Sqqg(p, ps, inda, indb);
	if (p.p(1).is_a(quark) && p.p(2).is_a(gluon) && p.p(3).is_a(quark) )
		return Sqgq(p, ps, inda, indb);
	if (p.p(1).is_a(gluino) && p.p(2).is_a(gluon) && p.p(3).is_a(gluino) )
		return Sqgq(p, ps, inda, indb);
	if (p.p(1).is_a(gluon) && p.p(2).is_a(quark) && p.p(3).is_a(quark) )
		return Sgqq(p, ps, inda, indb);
	if (p.p(1).is_a(gluon) && p.p(2).is_a(gluino) && p.p(3).is_a(gluino) )
		return Sgqq(p, ps, inda, indb);

	// old code commented out below - in case it's needed...
//	long code=scode(p);

//	_PRINT(code);
//	switch ( code ){
//	case 36: return Sggg(p, ps, inda, indb); break;
//	case 20: return Sqqg(p, ps, inda, indb); break;
//	case 12: return Sqgq(p, ps, inda, indb); break;
//	case 11: return Sgqq(p, ps, inda, indb); break;

//	default:
		_WARNING("Unknown tree splitting amplitude for process:");
		cerr << p <<endl ;
		return complex<T>(0,0);
//	}
}

template <class T> complex<T> Sggg(const process& p,momentum_configuration<T>& ps,
		int inda, int indb){
	if (p.p(1).helicity() == p.p(2).helicity() && p.p(1).helicity() == p.p(3).helicity())
	{	return complex<T>(0,0);
	}
	// note this assumes that z = Ka/(Ka+Kb)!!!!
	complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

	if (p.p(1).helicity() == p.p(2).helicity() && p.p(1).helicity() == 1)
		return (T(1)/sqrt(z)/sqrt(T(1)-z)/ps.spa(inda,indb));
	if (p.p(1).helicity() == p.p(2).helicity() && p.p(1).helicity() == -1)
		return (-T(1)/sqrt(z)/sqrt(T(1)-z)/ps.spb(inda,indb));
	if (p.p(1).helicity() == 1 && p.p(2).helicity() == -1)
	{
		if (p.p(3).helicity() == -1)
			return (-z*z/sqrt(z)/sqrt(T(1)-z)/ps.spb(inda,indb));
		if (p.p(3).helicity() == 1)
			return ((T(1)-z)*(T(1)-z)/sqrt(z)/sqrt(T(1)-z)/ps.spa(inda,indb));
	}
	if (p.p(1).helicity() == -1 && p.p(2).helicity() == 1)
	{
		if (p.p(3).helicity() == 1)
			return (z*z/sqrt(z)/sqrt(T(1)-z)/ps.spa(inda,indb));
		if (p.p(3).helicity() == -1)
			return (-(T(1)-z)*(T(1)-z)/sqrt(z)/sqrt(T(1)-z)/ps.spb(inda,indb));
	}

	_WARNING("Unknown tree splitting amplitude for process:");
	cerr << p <<endl ;
	return complex<T>(0,0);
}

template <class T> complex<T> Sqqg(const process& p,momentum_configuration<T>& ps,
		int inda, int indb){
	if (p.p(1).helicity() == p.p(2).helicity())
		// q and qb must have opposite helicity!
		return complex<T>(0,0);

	// note this assumes that z = Ka/(Ka+Kb)!!!!
	complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

	// qbar_1 q_2 -> g
	if ((p.p(1).helicity() == -1) && (p.p(1).is_anti_particle() == true) )
	{
		if (p.p(3).helicity() == -1)
			return ((T(1)-z)/ps.spb(inda,indb));
		if (p.p(3).helicity() == 1)
			return (z/ps.spa(inda,indb));
	}
	if ((p.p(1).helicity() == 1) && (p.p(1).is_anti_particle() == true) )
	{
		if (p.p(3).helicity() == 1)
			return ((T(1)-z)/ps.spa(inda,indb));
		if (p.p(3).helicity() == -1)
			return (z/ps.spb(inda,indb));
	}

	// q_1 qbar_2 -> g
	if ((p.p(1).helicity() == -1) && (p.p(1).is_anti_particle() == false) )
	{
		if (p.p(3).helicity() == 1)
			return ((z)/ps.spa(inda,indb));
		if (p.p(3).helicity() == -1)
			return ((T(1)-z)/ps.spb(inda,indb));
	
	}
	if ((p.p(1).helicity() == 1) && (p.p(1).is_anti_particle() == false) )
	{
		if (p.p(3).helicity() == -1)
			return ((z)/ps.spb(inda,indb));
		if (p.p(3).helicity() == 1)
			return ((T(1)-z)/ps.spa(inda,indb));
	}

	_WARNING("Unknown tree splitting amplitude for process:");
	cerr << p <<endl ;
	return complex<T>(0,0);

}

template <class T> complex<T> Sqgq(const process& p,momentum_configuration<T>& ps,
		int inda, int indb){


	if (p.p(1).helicity() == p.p(3).helicity())
		// q -> q  must have opposite helicity!
		return complex<T>(0,0);

	// note this assumes that z = Ka/(Ka+Kb)!!!!
	complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

	if (p.p(1).helicity() == 1 && p.p(2).helicity() == 1)
		return (T(1)/sqrt(T(1)-z)/ps.spa(inda,indb));
	if (p.p(1).helicity() == 1 && p.p(2).helicity() == -1)
		return (-z/sqrt(T(1)-z)/ps.spb(inda,indb));

	// parity conjugates
	if (p.p(1).helicity() == -1 && p.p(2).helicity() == -1)
		return (-T(1)/sqrt(T(1)-z)/ps.spb(inda,indb));
	if (p.p(1).helicity() == -1 && p.p(2).helicity() == 1)
		return (z/sqrt(T(1)-z)/ps.spa(inda,indb));

	_WARNING("Unknown tree splitting amplitude for process:");
	cerr << p <<endl ;
	return complex<T>(0,0);
}

template <class T> complex<T> Sgqq(const process& p,momentum_configuration<T>& ps,
		int inda, int indb){
	if (p.p(2).helicity() == p.p(3).helicity())
		// qb -> qb  must have opposite helicity!
		return complex<T>(0,0);

	// note this assumes that z = Ka/(Ka+Kb)!!!!
	complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

	if (p.p(2).helicity() == 1 && p.p(1).helicity() == 1)
		return (T(1)/sqrt(z)/ps.spa(inda,indb));
	if (p.p(2).helicity() == 1 && p.p(1).helicity() == -1)
		return (-(T(1)-z)/sqrt(z)/ps.spb(inda,indb));

	// parity conjugates
	if (p.p(2).helicity() == -1 && p.p(1).helicity() == -1)
		return (-T(1)/sqrt(z)/ps.spb(inda,indb));
	if (p.p(2).helicity() == -1 && p.p(1).helicity() == 1)
		return ((T(1)-z)/sqrt(z)/ps.spa(inda,indb));

	_WARNING("Unknown tree splitting amplitude for process:");
	cerr << p <<endl ;
	return complex<T>(0,0);
}

//long scode(const process& p){
//   if (p.n() > 3)
//    	return (-1);
//    long code = 0;
//    for (int i = 1; i <= 3; i++)
//    {
//    	code += p.p(i).type()->particle_nbr()*i;
//    }
//    return code;
//}

// explicit instantiations
template complex<R> Split0(const process& p, momentum_configuration<R>& ps,
		int inda, int indb);
template complex<R> Sggg(const process& p,momentum_configuration<R>& ps,
		int inda, int indb);
template complex<R> Sqqg(const process& p,momentum_configuration<R>& ps,
		int inda, int indb);
template complex<R> Sqgq(const process& p,momentum_configuration<R>& ps,
		int inda, int indb);
template complex<R> Sgqq(const process& p,momentum_configuration<R>& ps,
		int inda, int indb);

template complex<RHP> Split0(const process& p, momentum_configuration<RHP>& ps,
		int inda, int indb);
template complex<RHP> Sggg(const process& p,momentum_configuration<RHP>& ps,
		int inda, int indb);
template complex<RHP> Sqqg(const process& p,momentum_configuration<RHP>& ps,
		int inda, int indb);
template complex<RHP> Sqgq(const process& p,momentum_configuration<RHP>& ps,
		int inda, int indb);
template complex<RHP> Sgqq(const process& p,momentum_configuration<RHP>& ps,
		int inda, int indb);

template complex<RVHP> Split0(const process& p, momentum_configuration<RVHP>& ps,
		int inda, int indb);
template complex<RVHP> Sggg(const process& p,momentum_configuration<RVHP>& ps,
		int inda, int indb);
template complex<RVHP> Sqqg(const process& p,momentum_configuration<RVHP>& ps,
		int inda, int indb);
template complex<RVHP> Sqgq(const process& p,momentum_configuration<RVHP>& ps,
		int inda, int indb);
template complex<RVHP> Sgqq(const process& p,momentum_configuration<RVHP>& ps,
		int inda, int indb);


#if BH_USE_GMP
template complex<RGMP> Split0(const process& p, momentum_configuration<RGMP>& ps,
		int inda, int indb);
template complex<RGMP> Sggg(const process& p,momentum_configuration<RGMP>& ps,
		int inda, int indb);
template complex<RGMP> Sqqg(const process& p,momentum_configuration<RGMP>& ps,
		int inda, int indb);
template complex<RGMP> Sqgq(const process& p,momentum_configuration<RGMP>& ps,
		int inda, int indb);
template complex<RGMP> Sgqq(const process& p,momentum_configuration<RGMP>& ps,
		int inda, int indb);
#endif


}
