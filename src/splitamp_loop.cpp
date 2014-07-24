/*!\file splitamp_loop.cpp
\brief Implementation of loop splitting amplitudes
       \n as listed, e.g. in Bern, Dixon, Dunbar, Kosower,
                             Nucl. Phys. B 425, 217 (1994)
                             and 9409393 (primitive versus color-ordered splitamps)
       \n currently using Lance's dilog function
*/

// by CFB, last revision 1/6/2009
// 8/28/08 added separate functions for cut/rational parts
// 11/18/08 -- leading_c, LT and RT updated, new color structures primitive versus 
//             color-ordered
// 12/21/08 - added gluinos
// note: s_in_loop = [0], f_in_loop = [1/2], v_in_loop = [1] (cf. 9403226)
//       f_left = [L], f_right = [R], [s] = [0] = s_in_loop, [f] = ferm_c (cf. 9409393)
//       primitive versus color-ordered amplitudes
// 1/6/09 - bug fixed

#include "factorization.h"
#include "particles.h"
#include <cmath>
#include "BH_utilities.h"
#include "qd_suppl.h"
#include "process.h"
#include "polylog.h"  //  for pi
#include "BH_debug.h"

#if BH_USE_GMP
#include "gmp_r.h"
#endif


#define DIMREG 0
#define NNc 3   // number of colors
#define Theta(r,v) ( (r) >= T(0) ? (v) : T(0) )

namespace BH {

template <class T> Series<complex<T> > SplitS1(const process& p, int in_loop, momentum_configuration<T>& ps, int inda, int indb, T mu2, int full_cut_rat){
		Series<complex<T> > res=Series<complex<T> >(-2,0,Split1(p,in_loop,ps,inda,indb,-2,mu2,full_cut_rat), \
			Split1(p,in_loop,ps,inda,indb,-1,mu2,full_cut_rat), \
			Split1(p,in_loop,ps,inda,indb,0,mu2,full_cut_rat));
		BH_DEBUG_MESSAGE2("SplitS1",res);

		return res;
};


template <class T> complex<T> Split1(const process& p, int in_loop, momentum_configuration<T>& ps, int inda, int indb, int oeps, T mu2, int full_cut_rat){

	if (oeps > 0)
	{
		_WARNING("Implemented only till order eps^0");
		return complex<T>(0,0);
	}
	if (oeps < -2)
		return complex<T>(0,0);

	if (oeps !=0 && full_cut_rat==2)
		return complex<T>(0,0);

	if (p.p(1).is_a(gluon) && p.p(2).is_a(gluon) && p.p(3).is_a(gluon) ){
		if(full_cut_rat==0) return Sggg1(p, in_loop, ps, inda, indb, oeps, mu2);
		if(full_cut_rat==1) return Sggg1_Cut(p, in_loop, ps, inda, indb, oeps, mu2);
		if(full_cut_rat==2) return Sggg1_Rat(p, in_loop, ps, inda, indb);
	}
	if ((p.p(1).is_a(gluino) && p.p(2).is_a(gluino) && p.p(3).is_a(gluon) && p.p(1).flavor()==p.p(2).flavor())||
		(p.p(1).is_a(quark) && p.p(2).is_a(quark) && p.p(3).is_a(gluon) && p.p(1).flavor()==p.p(2).flavor())){
		if(full_cut_rat==0) return Sqqg1(p, in_loop, ps, inda, indb, oeps, mu2);
		if(full_cut_rat==1) return Sqqg1_Cut(p, in_loop, ps, inda, indb, oeps, mu2);
		if(full_cut_rat==2) return Sqqg1_Rat(p, in_loop, ps, inda, indb);
	}
	if ((p.p(1).is_a(gluino) && p.p(2).is_a(gluon) && p.p(3).is_a(gluino) && p.p(1).flavor()==p.p(3).flavor()) ||
		(p.p(1).is_a(quark) && p.p(2).is_a(gluon) && p.p(3).is_a(quark) && p.p(1).flavor()==p.p(3).flavor())){
		if(full_cut_rat==0) return Sqgq1(p, in_loop, ps, inda, indb, oeps, mu2);
		if(full_cut_rat==1) return Sqgq1_Cut(p, in_loop, ps, inda, indb, oeps, mu2);
		if(full_cut_rat==2) return Sqgq1_Rat(p, in_loop, ps, inda, indb);
	}
	if ((p.p(1).is_a(gluon) && p.p(2).is_a(gluino) && p.p(3).is_a(gluino) && p.p(2).flavor()==p.p(3).flavor())
		|| (p.p(1).is_a(gluon) && p.p(2).is_a(quark) && p.p(3).is_a(quark) && p.p(2).flavor()==p.p(3).flavor())){
		if(full_cut_rat==0) return Sgqq1(p, in_loop, ps, inda, indb, oeps, mu2);
		if(full_cut_rat==1) return Sgqq1_Cut(p, in_loop, ps, inda, indb, oeps, mu2);
		if(full_cut_rat==2) return Sgqq1_Rat(p, in_loop, ps, inda, indb);
	}
		_WARNING("Unknown loop splitting amplitude for process SplitS1:");
		cerr << p <<endl ;
		return complex<T>(0,0);

}



template <class T> complex<T> Sggg1_Cut(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	BH_DEBUG_MESSAGE4("Sggg1_Cut ",in_loop," ",p);
	switch (in_loop){
	case s_in_loop:  return Sggg1s_Cut(p, ps, inda, indb, oeps, mu2); break;
	case f_in_loop: return -Sggg1s_Cut(p, ps, inda, indb, oeps, mu2); break;
	case v_in_loop:  return Sggg1v_Cut(p, ps, inda, indb, oeps, mu2); break;
	case f_left: return Sggg1v_Cut(p, ps, inda, indb, oeps, mu2); break;
	case ferm_c: return complex<T>(0,0); break;
	case leading_c: return Sggg1v_Cut(p, ps, inda, indb, oeps, mu2); break;
	case sub_leading_c: return complex<T>(0,0); break;
	default:
		_WARNING("States not yet implemented");
		return complex<T>(0,0);
	}

	_WARNING("Unknown loop splitting amplitude for process:");
	cerr << p <<endl ;
	return complex<T>(0,0);
}

template <class T> complex<T> Sggg1_Rat(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb)
{
	BH_DEBUG_MESSAGE4("Sggg1_Rat ",in_loop," ",p);
	switch (in_loop){
	case s_in_loop:  return Sggg1s_Rat(p, ps, inda, indb); break;
	case f_in_loop: return -Sggg1s_Rat(p, ps, inda, indb); break;
	case v_in_loop:  return Sggg1v_Rat(p, ps, inda, indb); break;
	case f_left: return Sggg1v_Rat(p, ps, inda, indb); break;
	case ferm_c: return complex<T>(0,0); break;
	case sub_leading_c: return complex<T>(0,0); break;
	case leading_c: return Sggg1v_Rat(p, ps, inda, indb); break;
	default:
		_WARNING("States not yet implemented");
		return complex<T>(0,0);
	}

	_WARNING("Unknown loop splitting amplitude for process:");
	cerr << p <<endl ;
	return complex<T>(0,0);
}


template <class T> complex<T> Sggg1(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	BH_DEBUG_MESSAGE4("Sggg1",in_loop," ",p);
	switch (in_loop){
	case s_in_loop:  return Sggg1s(p, ps, inda, indb, oeps, mu2); break;
	case f_in_loop: return -Sggg1s(p, ps, inda, indb, oeps, mu2); break;
	case v_in_loop:  return Sggg1v(p, ps, inda, indb, oeps, mu2); break;
	case f_left: return Sggg1v(p, ps, inda, indb, oeps, mu2); break;
	case ferm_c: return complex<T>(0,0); break;
	case sub_leading_c: return complex<T>(0,0); break;
	case leading_c: return Sggg1v(p, ps, inda, indb, oeps, mu2); break;
	default:
		_WARNING("States not yet implemented");
		return complex<T>(0,0);
	}

	_WARNING("Unknown loop splitting amplitude for process:");
	cerr << p <<endl ;
	return complex<T>(0,0);
}

template <class T> complex<T> Sggg1s(const process& p, momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
    if (oeps < 0)
    	return complex<T>(0,0);
    // no 1/eps^2 or 1/eps terms here

    if (p.p(1).helicity() != p.p(2).helicity())
    	return complex<T>(0,0);

	// note this assumes that z = Ka/(Ka+Kb)!!!!
	complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (p.p(1).helicity() == p.p(2).helicity() && p.p(1).helicity() == p.p(3).helicity())
    {
    	if (p.p(1).helicity() == 1)
    		return (-T(1)/T(3)*sqrt(z)*sqrt(T(1)-z)*ps.spb(inda,indb)/ \
    				ps.spa(inda,indb)/ps.spa(inda,indb));
    	if (p.p(1).helicity() == -1)
    		return (T(1)/T(3)*sqrt(z)*sqrt(T(1)-z)*ps.spa(inda,indb)/ \
    				ps.spb(inda,indb)/ps.spb(inda,indb));
    }

    if (p.p(1).helicity() == p.p(2).helicity())
    	return (T(1)/T(3)*z*(T(1)-z)*Split0(p, ps, inda, indb));
    // -++ or +-- (which is decided in Split0)

	_WARNING("Unknown loop splitting amplitude for process:");
	cerr << p <<endl ;
	return complex<T>(0,0);
}

template <class T> complex<T> Sggg1s_Cut(const process& p, momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
   	return complex<T>(0,0);
}

template <class T> complex<T> Sggg1s_Rat(const process& p, momentum_configuration<T>& ps,
		int inda, int indb)
{
    return Sggg1s(p,ps,inda,indb,0,T(1));
    // purely rational
}


template <class T> complex<T> Sggg1v(const process& p, momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	BH_DEBUG_MESSAGE3("Sggg1v"," ",p);
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (p.p(1).helicity() == p.p(2).helicity() && p.p(1).helicity() == p.p(3).helicity())
    // our lovely nonfactorizing function...
    {
    	if (oeps < 0)
    		return complex<T>(0,0);
    	if (p.p(1).helicity() == 1 && oeps == 0)
    		return (-T(1)/T(3)*sqrt(z)*sqrt(T(1)-z)*ps.spb(inda,indb)/ \
    				ps.spa(inda,indb)/ps.spa(inda,indb));
    	if (p.p(1).helicity() == -1 && oeps == 0)
    		return (T(1)/T(3)*sqrt(z)*sqrt(T(1)-z)*ps.spa(inda,indb)/ \
    				ps.spb(inda,indb)/ps.spb(inda,indb));
    }
    if (p.p(1).helicity() == p.p(2).helicity())
    // ++ -> - or -- -> + (which is decided in Split0)
    {	switch (oeps) {
            case (-2): return (- Split0(p,ps,inda,indb));
            case (-1):
            	       return ( -(-Clog(ps.s(inda,indb),mu2)-(log(z*(T(1)-z))))* \
            		                Split0(p,ps,inda,indb));
            case (0): return ( -( (log(z*(T(1)-z))+Clog(ps.s(inda,indb),mu2))*   \
            		(log(z*(T(1)-z))+Clog(ps.s(inda,indb),mu2))/T(2) -T(2)*log(z)*log(T(1)-z) - \
            		 T(1)/T(3)*z*(T(1)-z)+T(1)/T(6)*(pi<T>())*(pi<T>()))*Split0(p,ps,inda,indb));
        }
    }

    if (p.p(1).helicity() != p.p(2).helicity())
    // +- -> - or +, -+ -> - or +, again decided in Split0
    // almost identical to the above, except for the 1/3 z (1-z) in eps^0
    {
    	switch (oeps) {
    	    case (-2): return (-Split0(p,ps,inda,indb));
    	    case (-1): return (-(-Clog(ps.s(inda,indb),mu2)-(log(z*(T(1)-z))))* \
    	     		                Split0(p,ps,inda,indb));
    	    case (0): return ( -( (log(z*(T(1)-z))+Clog(ps.s(inda,indb),mu2))*   \
            		(log(z*(T(1)-z))+Clog(ps.s(inda,indb),mu2))/T(2) -T(2)*log(z)*log(T(1)-z) + \
    	          		 T(1)/T(6)*(pi<T>())*(pi<T>()))*Split0(p,ps,inda,indb));
    	}
    }

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    return complex<T>(0,0);
}

template <class T> complex<T> Sggg1v_Cut(const process& p, momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	BH_DEBUG_MESSAGE3("Sggg1v_Cut"," ",p);
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (p.p(1).helicity() == p.p(2).helicity() && p.p(1).helicity() == p.p(3).helicity())
    // our lovely nonfactorizing function...
    {
    	return complex<T>(0,0);
    }
    if (p.p(1).helicity() == p.p(2).helicity())
    // ++ -> - or -- -> + (which is decided in Split0)
    {	switch (oeps) {
            case (-2): return (- Split0(p,ps,inda,indb));
            case (-1):
            	       return ( -(-Clog(ps.s(inda,indb),mu2)-(log(z*(T(1)-z))))* \
            		                Split0(p,ps,inda,indb));
            case (0): return ( -( (log(z*(T(1)-z))+Clog(ps.s(inda,indb),mu2))*   \
            		(log(z*(T(1)-z))+Clog(ps.s(inda,indb),mu2))/T(2) -T(2)*log(z)*log(T(1)-z) + \
            		 T(1)/T(6)*(pi<T>())*(pi<T>()))*Split0(p,ps,inda,indb));
        }
    }

    if (p.p(1).helicity() != p.p(2).helicity())
    // +- -> - or +, -+ -> - or +, again decided in Split0
    // almost identical to the above, except for the 1/3 z (1-z) in eps^0
    {
    	switch (oeps) {
    	    case (-2): return (- Split0(p,ps,inda,indb));
    	    case (-1): return (-(-Clog(ps.s(inda,indb),mu2)-(log(z*(T(1)-z))))* \
    	     		                Split0(p,ps,inda,indb));
    	    case (0): return ( -( (log(z*(T(1)-z))+Clog(ps.s(inda,indb),mu2))*   \
            		(log(z*(T(1)-z))+Clog(ps.s(inda,indb),mu2))/T(2) -T(2)*log(z)*log(T(1)-z) + \
    	          		 T(1)/T(6)*(pi<T>())*(pi<T>()))*Split0(p,ps,inda,indb));
    	}
    }

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    return complex<T>(0,0);
}

template <class T> complex<T> Sggg1v_Rat(const process& p, momentum_configuration<T>& ps,
		int inda, int indb)
{
	BH_DEBUG_MESSAGE3("Sggg1v_Rat"," ",p);
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (p.p(1).helicity() == p.p(2).helicity() && p.p(1).helicity() == p.p(3).helicity())
    // our lovely nonfactorizing function...
    {

    	if (p.p(1).helicity() == 1)
    		return (-T(1)/T(3)*sqrt(z)*sqrt(T(1)-z)*ps.spb(inda,indb)/ \
    				ps.spa(inda,indb)/ps.spa(inda,indb));
    	if (p.p(1).helicity() == -1)
    		return (T(1)/T(3)*sqrt(z)*sqrt(T(1)-z)*ps.spa(inda,indb)/ \
    				ps.spb(inda,indb)/ps.spb(inda,indb));
    }
    if (p.p(1).helicity() == p.p(2).helicity())
    // ++ -> - or -- -> + (which is decided in Split0)
    {	 return ( (  T(1)/T(3)*z*(T(1)-z))*Split0(p,ps,inda,indb));
        
    }

    if (p.p(1).helicity() != p.p(2).helicity())
    // +- -> - or +, -+ -> - or +, again decided in Split0
    {
    	return ( complex<T>(0,0));
    }

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    return complex<T>(0,0);
}

template <class T> complex<T> Sqgq1(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());
	BH_DEBUG_MESSAGE8("Sqgq1 ",in_loop," ",p," z: ",z," mu: ",mu2);
    
    if (in_loop == s_in_loop)
    	return complex<T>(0,0);
    if (in_loop == f_in_loop)
    	return complex<T>(0,0);    
    if (in_loop == ferm_c)
    	return complex<T>(0,0);

	//gluinos included 
   if (p.p(1).helicity() != p.p(2).helicity() && (! p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
        	return (ff(T(1)-z, f_left, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    if (p.p(1).helicity() != p.p(2).helicity() && ( p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
            	return (ff(T(1)-z, f_left, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    if (p.p(1).helicity() == p.p(2).helicity() && (! p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
    	return ( (ff(T(1)-z, f_left, ps.s(inda,indb),mu2,oeps) + (oeps == 0 ? (T(1))/T(2)*(T(1)-z) : T(0)))* \
    			Split0(p,ps,inda,indb));
    if (p.p(1).helicity() == p.p(2).helicity() && ( p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
    	return ( (ff(T(1)-z, f_left, ps.s(inda,indb),mu2,oeps) - (oeps == 0 ? (T(1))/T(2)*(T(1)-z) : T(0)))* \
    			Split0(p,ps,inda,indb));    


    
    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("Specify in loop, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Sgqq1(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	BH_DEBUG_MESSAGE4("Sgqq1 ",in_loop," ",p);
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (in_loop == s_in_loop)
    	return complex<T>(0,0);
    if (in_loop == f_in_loop)
    	return complex<T>(0,0);    
    if (in_loop == ferm_c)
    	return complex<T>(0,0);

	//gluinos included 
    if (p.p(1).helicity() != p.p(2).helicity() && ( p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
        	return (ff(z, f_left, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    if (p.p(1).helicity() != p.p(2).helicity() && ( !p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
            	return (ff(z, f_right, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    if (p.p(1).helicity() == p.p(2).helicity() && ( p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
    	return ( (ff(z, f_left, ps.s(inda,indb),mu2,oeps) + (oeps == 0 ? (T(1))/T(2)*(z) : T(0)))* \
    			Split0(p,ps,inda,indb));
    if (p.p(1).helicity() == p.p(2).helicity() && ( ! p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
    	return ( (ff(z, f_right, ps.s(inda,indb),mu2,oeps) - (oeps == 0 ? (T(1))/T(2)*(z) : T(0)))* \
    			Split0(p,ps,inda,indb));    

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("Specify in loop, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Sqgq1_Cut(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());
	BH_DEBUG_MESSAGE8("Sqgq1_Cut ",in_loop," ",p," z: ",z," mu: ",mu2);
    if (in_loop == s_in_loop)
    	return complex<T>(0,0);
    if (in_loop == f_in_loop)
    	return complex<T>(0,0);    
    if (in_loop == ferm_c)
    	return complex<T>(0,0);
   
 
    if (p.p(1).helicity() != p.p(2).helicity() && (! p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
        	return (ff(T(1)-z, f_left, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    if (p.p(1).helicity() != p.p(2).helicity() && ( p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
            	return (ff(T(1)-z, f_right, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    if (p.p(1).helicity() == p.p(2).helicity() && (! p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
    	return (ff(T(1)-z, f_left, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    if (p.p(1).helicity() == p.p(2).helicity() && ( p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
    	return (ff(T(1)-z, f_right, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));    

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("Specify in loop, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Sqgq1_Rat(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb)
{
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());
        
    if (in_loop == s_in_loop)
    	return complex<T>(0,0);
    if (in_loop == f_in_loop)
    	return complex<T>(0,0);    
    if (in_loop == ferm_c)
    	return complex<T>(0,0);

    if (p.p(1).helicity() != p.p(2).helicity())
        	return (complex<T>(0,0));
    //if (p.p(1).helicity() == p.p(2).helicity() && (p.p(1).is_a(qm) || p.p(1).is_a(qp)))
    if (p.p(1).helicity() == p.p(2).helicity() && (! p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
    	return ((T(1))/T(2)*(T(1)-z)*Split0(p,ps,inda,indb));
    //if (p.p(1).helicity() == p.p(2).helicity() && (p.p(1).is_a(qbm) || p.p(1).is_a(qbp)))
    if (p.p(1).helicity() == p.p(2).helicity() && ( p.p(1).is_anti_particle() && (p.p(1).is_a(quark) || p.p(1).is_a(gluino))))
    	return (-((T(1))/T(2)*(T(1)-z))*Split0(p,ps,inda,indb));    

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("Specify in loop, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Sgqq1_Cut(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	BH_DEBUG_MESSAGE4("Sgqq1_Cut ",in_loop," ",p);
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());
    
    if (in_loop == s_in_loop)
    	return complex<T>(0,0);
    if (in_loop == f_in_loop)
    	return complex<T>(0,0);    
    if (in_loop == ferm_c)
    	return complex<T>(0,0);
    
    //if (p.p(1).helicity() != p.p(2).helicity() && (p.p(1).is_a(qbm) || p.p(1).is_a(qbp)))
    if (p.p(1).helicity() != p.p(2).helicity() && ( p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
        	return (ff(z, f_left, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    //if (p.p(1).helicity() != p.p(2).helicity() && (p.p(1).is_a(qm) || p.p(1).is_a(qp)))
    if (p.p(1).helicity() != p.p(2).helicity() && ( !p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
            	return (ff(z, f_right, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    //if (p.p(1).helicity() == p.p(2).helicity() && (p.p(1).is_a(qbm) || p.p(1).is_a(qbp)))
    if (p.p(1).helicity() == p.p(2).helicity() && ( p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
    	return (ff(z, f_left, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));
    //if (p.p(1).helicity() == p.p(2).helicity() && (p.p(1).is_a(qm) || p.p(1).is_a(qp)))
    if (p.p(1).helicity() == p.p(2).helicity() && ( ! p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
    	return (ff(z, f_right, ps.s(inda,indb),mu2,oeps)*Split0(p,ps,inda,indb));    

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("Specify in loop, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Sgqq1_Rat(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb)
{
	BH_DEBUG_MESSAGE4("Sgqq1_Rat ",in_loop," ",p);
	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (in_loop == s_in_loop)
    	return complex<T>(0,0);
    if (in_loop == f_in_loop)
    	return complex<T>(0,0);    
    if (in_loop == ferm_c)
    	return complex<T>(0,0);

    if (p.p(1).helicity() != p.p(2).helicity())
        	return (complex<T>(0,0));
    //if (p.p(1).helicity() == p.p(2).helicity() && (p.p(1).is_a(qm) || p.p(1).is_a(qp)))
    if (p.p(1).helicity() == p.p(2).helicity() && ( ! p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
    	return (-(T(1))/T(2)*(z)*Split0(p,ps,inda,indb));
    //if (p.p(1).helicity() == p.p(2).helicity() && (p.p(1).is_a(qbm) || p.p(1).is_a(qbp)))
    if (p.p(1).helicity() == p.p(2).helicity() && ( p.p(2).is_anti_particle() && ( p.p(2).is_a(quark) || p.p(2).is_a(gluino))))
    	return (((T(1))/T(2)*(z))*Split0(p,ps,inda,indb));    

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("Specify in loop, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}


template <class T> complex<T> ff(complex<T> z, int in_loop, complex<T> s, T mu2, int oeps)
{
	BH_DEBUG_MESSAGE8("ff: z: ",z," in_loop: ",in_loop," s:",s," mu2: ",mu2);
	BH_DEBUG_MESSAGE4("f_left: ",f_left," f_right: ",f_right);

	T zz = z.real();

    if (in_loop == v_in_loop)
    	return (ff(z,f_left,s,mu2,oeps)-T(1)/T(NNc)/T(NNc)*ff(z,f_right,s,mu2,oeps));
	
	if (in_loop == f_left || in_loop == leading_c)
	{
		switch (oeps){
			case (-2): return (-complex<T>(1,0));
			case (-1): return ( Clog(s,mu2) + log(z) );
			case (0): {
				complex<T> res=(-T(1)/T(2)*(Clog(s,mu2) + log(z))*(Clog(s,mu2)+ log(z)) - li2(T(1)-zz) );
				BH_DEBUG_MESSAGE2("f[z,s,mu2]: ",res);
				return res;
			}
			default: _WARNING("Implemented only till order eps^0");
			         return complex<T>(0,0);
		}
	}

	if (in_loop == f_right || in_loop == sub_leading_c)
	{
		switch (oeps){
			case (-2): return complex<T>(0,0);
			case (-1): return ( (log(T(1)-z)) );
			case (0): return ( -T(1)/T(2)*log(T(1)-z)*log(T(1)-z)-Clog(s,mu2)*log(T(1)-z)- li2(zz));
			default: _WARNING("Implemented only till order eps^0");
			         return complex<T>(0,0);
		}
	}
    _WARNING("Unknown particle in loop in f(z)!");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Sqqg1(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
    if (p.p(1).helicity() == p.p(2).helicity())
    	return complex<T>(0,0);

	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (in_loop == s_in_loop)
    {
    	switch (oeps){
    	case (-2): return complex<T>(0,0);
    	case (-1): return (-T(1)/T(3)*Split0(p,ps,inda,indb));
    	case (0):  return ( (T(1)/T(3)*Clog(ps.s(inda,indb),mu2)-T(8)/T(9) )* \
    			Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}
    }
    if (in_loop == f_in_loop)
    {
    	switch (oeps){
    	case (-2): return complex<T>(0,0);
    	case (-1): return (-T(2)/T(3)*Split0(p,ps,inda,indb));
    	case (0):  return ( (T(2)/T(3)*Clog(ps.s(inda,indb),mu2)-T(10)/T(9) )* \
    			Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}
    }
    if (in_loop == f_left)
    {
    	switch (oeps){
    	case (-2): return complex<T>(0,0);
    	case (-1): return ((log(z) + log(T(1)-z) + T(13)/T(6))*Split0(p,ps,inda,indb));
    	case (0): return ( (T(1)/T(2)*(- (log(z))*(log(z)) - (log(T(1)-z))*(log(T(1)-z))) - \
    			Clog(ps.s(inda,indb),mu2)*(log(z*(T(1)-z))+ T(13)/T(6)) - T(DIMREG)/T(6) +\
    			(log(z))*(log(T(1)-z)) - pi<T>()*pi<T>()/T(6)+T(83)/T(18) )* \
    			 Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}

    }
    if (in_loop == leading_c) // == f_left
    {
    	switch (oeps){
    	case (-2): return complex<T>(0,0);
    	case (-1): return ((log(z) + log(T(1)-z) + T(13)/T(6))*Split0(p,ps,inda,indb));
    	case (0): return ( (T(1)/T(2)*(- (log(z))*(log(z)) - (log(T(1)-z))*(log(T(1)-z))) - \
    			Clog(ps.s(inda,indb),mu2)*(log(z*(T(1)-z))+ T(13)/T(6)) - T(DIMREG)/T(6) +\
    			(log(z))*(log(T(1)-z)) - pi<T>()*pi<T>()/T(6)+T(83)/T(18) )* \
    			 Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}

    }
    if (in_loop == f_right || in_loop == sub_leading_c)
    {
    	switch (oeps){
    	case (-2): return (-Split0(p,ps,inda,indb));
    	case (-1): return (-(-Clog(ps.s(inda,indb),mu2) + T(3)/T(2))* \
    	            Split0(p,ps,inda,indb));
    	case (0): return ( - (T(1)/T(2)*(-Clog(ps.s(inda,indb),mu2))* \
    			(-Clog(ps.s(inda,indb),mu2)) +\
    			 T(3)/T(2)*(Clog(ps.s(inda,indb),mu2))+ T(DIMREG)/T(2) + T(7)/T(2))* \
    			 Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}
    }
    if (in_loop == v_in_loop)
    	return (Sqqg1(p, leading_c, ps, inda, indb, oeps, mu2)-T(1)/T(NNc)/T(NNc)* \
    			Sqqg1(p, f_right, ps, inda, indb, oeps, mu2));
    if (in_loop == ferm_c)
    	return (-Sqqg1(p, s_in_loop, ps, inda, indb, oeps, mu2)- \
    	    			Sqqg1(p, f_in_loop, ps, inda, indb, oeps, mu2));
    
    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("If vector [1], specify f_left or f_right, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Sqqg1_Cut(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2)
{
	BH_DEBUG_MESSAGE4("Sqqg1_Cut ",in_loop," ",p);
	BH_DEBUG_MESSAGE2("leading_color: ",leading_c); 
    if (p.p(1).helicity() == p.p(2).helicity())
    	return complex<T>(0,0);

	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (in_loop == s_in_loop)
    {
    	switch (oeps){
    	case (-2): return complex<T>(0,0);
    	case (-1): return (-T(1)/T(3)*Split0(p,ps,inda,indb));
    	case (0):  return ( (T(1)/T(3)*Clog(ps.s(inda,indb),mu2) )* \
    			Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}
    }
    if (in_loop == f_in_loop)
    {
    	switch (oeps){
    	case (-2): return complex<T>(0,0);
    	case (-1): return (-T(2)/T(3)*Split0(p,ps,inda,indb));
    	case (0):  return ( (T(2)/T(3)*Clog(ps.s(inda,indb),mu2)-T(4)/T(3))* \
    			Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}
    }
    if (in_loop == f_left)
    {
    	switch (oeps){
    	case (-2): return complex<T>(0,0);
    	case (-1): return ((log(z) + log(T(1)-z) + T(13)/T(6))*Split0(p,ps,inda,indb));
    	case (0): return ( (T(1)/T(2)*(- (log(z))*(log(z)) - (log(T(1)-z))*(log(T(1)-z))) - \
    			Clog(ps.s(inda,indb),mu2)*(log(z*(T(1)-z))+ T(13)/T(6))  +\
    			(log(z))*(log(T(1)-z)) - pi<T>()*pi<T>()/T(6) )* \
    			 Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}

    }
    if (in_loop == leading_c)
    {
    	switch (oeps){
    	case (-2): return complex<T>(0,0);
    	case (-1): return ((log(z) + log(T(1)-z) + T(13)/T(6))*Split0(p,ps,inda,indb));
    	case (0):{

		complex<T> rs=T(1)/T(2)*(-log(z)*log(z)-log(T(1)-z)*log(T(1)-z))
    			- Clog(ps.s(inda,indb),mu2)*(log(z*(T(1)-z))+ T(13)/T(6)) 
    			+ (log(z))*(log(T(1)-z)) + T(83)/T(18) - T(5)/T(18) /* def or rational */ - pi<T>()*pi<T>()/T(6) ; 


		BH_DEBUG_MESSAGE8("r_s: ",rs," z=",z," mu2=",mu2," s=",ps.s(inda,indb));

		return (rs*Split0(p,ps,inda,indb));
	}


    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}

    }
    if (in_loop == f_right || in_loop == sub_leading_c)
    {
    	switch (oeps){
    	case (-2): return (-Split0(p,ps,inda,indb));
    	case (-1): return (-(-Clog(ps.s(inda,indb),mu2) + T(3)/T(2))* \
    	            Split0(p,ps,inda,indb));
    	case (0): return ( - (T(1)/T(2)*(-Clog(ps.s(inda,indb),mu2))* \
    			(-Clog(ps.s(inda,indb),mu2)) +\
    			 T(3)/T(2)*(Clog(ps.s(inda,indb),mu2)))* \
    			 Split0(p,ps,inda,indb));
    	default: _WARNING("Implemented only till order eps^0");
    				         return complex<T>(0,0);
    	}
    }
    if (in_loop == v_in_loop)
    	return (Sqqg1_Cut(p, leading_c, ps, inda, indb, oeps, mu2)-T(1)/T(NNc)/T(NNc)* \
    			Sqqg1_Cut(p, f_right, ps, inda, indb, oeps, mu2));
    if (in_loop == ferm_c)
    	return (-Sqqg1_Cut(p, s_in_loop, ps, inda, indb, oeps, mu2)- \
    	    			Sqqg1_Cut(p, f_in_loop, ps, inda, indb, oeps, mu2));

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("If vector [1], specify f_left or f_right, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Sqqg1_Rat(const process& p, int in_loop,
		momentum_configuration<T>& ps,
		int inda, int indb)
{
	BH_DEBUG_MESSAGE4("Sqqg1_Rat ",in_loop," ",p);
    if (p.p(1).helicity() == p.p(2).helicity())
    	return complex<T>(0,0);

	// note this assumes that z = Ka/(Ka+Kb)!!!!
    complex<T> z = ps.p(inda).E()/(ps.p(inda).E()+ps.p(indb).E());

    if (in_loop == s_in_loop)
    {
        return ( (-T(8)/T(9) )*Split0(p,ps,inda,indb));
    }
    if (in_loop == f_in_loop)
    {
    	return ( (-T(10)/T(9) + T(4)/T(3) )*Split0(p,ps,inda,indb));
    }
    if (in_loop == f_left)
    {
        return ( -(  T(DIMREG)/T(6) -T(83)/T(18) )*Split0(p,ps,inda,indb));
    }
    if (in_loop == leading_c)
    {
        return ( -(  T(DIMREG)/T(6) -T(5)/T(18) )*Split0(p,ps,inda,indb));
    }
    if (in_loop == f_right || in_loop == sub_leading_c)
    {
        return ( - ( T(DIMREG)/T(2) + T(7)/T(2))*Split0(p,ps,inda,indb));
    }
    if (in_loop == v_in_loop)
    	return (Sqqg1_Rat(p, leading_c, ps, inda, indb)-T(1)/T(NNc)/T(NNc)* \
    			Sqqg1_Rat(p, f_right, ps, inda, indb));
    if (in_loop == ferm_c)
    	return (-Sqqg1_Rat(p, s_in_loop, ps, inda, indb)- \
    	    			Sqqg1_Rat(p, f_in_loop, ps, inda, indb));

    _WARNING("Unknown loop splitting amplitude for process:");
    cerr << p <<endl ;
    _WARNING("If vector [1], specify f_left or f_right, I don't know ");
    cerr << in_loop << endl;
    return complex<T>(0,0);
}

template <class T> complex<T> Clog(complex<T> s1, T mu2)
{
 T s1r = real(s1);
 complex<T> res=complex<T>(log(abs(s1r/mu2)),Theta(s1r,-pi<T>()));
 return res;
}

template <class T> T li2(T x){ return ReLi2(x); }

// explicit instantiations
template complex<R> Clog(complex<R> s1, R mu2);
template R li2(R x);
template complex<R> ff(complex<R> z, int in_loop, complex<R> s, R mu2, int oeps);
template Series<complex<R> > SplitS1(const process& p, int in_loop, momentum_configuration<R>& ps, int inda, int indb, R mu2,int full_cut_rat);
template complex<R> Split1(const process& p, int in_loop,momentum_configuration<R>& ps,	int inda, int indb, int oeps, R mu2,int full_cut_rat);
template complex<R> Sggg1(const process& p, int in_loop,momentum_configuration<R>& ps,	int inda, int indb, int oeps, R mu2);
template complex<R> Sggg1s(const process& p, momentum_configuration<R>& ps,int inda, int indb, int oeps, R mu2);
template complex<R> Sggg1v(const process& p, momentum_configuration<R>& ps,int inda, int indb, int oeps, R mu2);
template complex<R> Sqgq1(const process& p, int in_loop,momentum_configuration<R>& ps, int inda, int indb, int oeps, R mu2);
template complex<R> Sgqq1(const process& p, int in_loop, momentum_configuration<R>& ps, int inda, int indb, int oeps, R mu2);
template complex<R> Sqqg1(const process& p, int in_loop,momentum_configuration<R>& ps, int inda, int indb, int oeps, R mu2);

template complex<R> Sggg1_Cut(const process& p, int in_loop,momentum_configuration<R>& ps,int inda, int indb, int oeps, R mu2);
template complex<R> Sggg1s_Cut(const process& p, momentum_configuration<R>& ps,	int inda, int indb, int oeps, R mu2);
template complex<R> Sggg1v_Cut(const process& p, momentum_configuration<R>& ps,	int inda, int indb, int oeps, R mu2);
template complex<R> Sqgq1_Cut(const process& p, int in_loop,momentum_configuration<R>& ps,int inda, int indb, int oeps, R mu2);
template complex<R> Sgqq1_Cut(const process& p, int in_loop,momentum_configuration<R>& ps,int inda, int indb, int oeps, R mu2);
template complex<R> Sqqg1_Cut(const process& p, int in_loop,momentum_configuration<R>& ps,int inda, int indb, int oeps, R mu2);

template complex<R> Sggg1_Rat(const process& p, int in_loop,momentum_configuration<R>& ps,int inda, int indb);
template complex<R> Sggg1s_Rat(const process& p, momentum_configuration<R>& ps,	int inda, int indb);
template complex<R> Sggg1v_Rat(const process& p, momentum_configuration<R>& ps,	int inda, int indb);
template complex<R> Sqgq1_Rat(const process& p, int in_loop,momentum_configuration<R>& ps,int inda, int indb);
template complex<R> Sgqq1_Rat(const process& p, int in_loop,momentum_configuration<R>& ps,int inda, int indb);
template complex<R> Sqqg1_Rat(const process& p, int in_loop,momentum_configuration<R>& ps,int inda, int indb);

template RHP li2(RHP x);
template complex<RHP> Clog(complex<RHP> s1, RHP mu2);
template complex<RHP> ff(complex<RHP> z, int in_loop, complex<RHP> s, RHP mu2, int oeps);
template Series<complex<RHP> > SplitS1(const process& p, int in_loop,momentum_configuration<RHP>& ps,int inda, int indb, RHP mu2,int full_cut_rat);
template complex<RHP> Split1(const process& p, int in_loop,momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2,int full_cut_rat);
template complex<RHP> Sggg1(const process& p, int in_loop,momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sggg1s(const process& p, momentum_configuration<RHP>& ps,	int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sggg1v(const process& p, momentum_configuration<RHP>& ps,	int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sqgq1(const process& p, int in_loop,momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sgqq1(const process& p, int in_loop,momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sqqg1(const process& p, int in_loop, momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);

template complex<RHP> Sggg1_Cut(const process& p, int in_loop,	momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sggg1s_Cut(const process& p, momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sggg1v_Cut(const process& p, momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sqgq1_Cut(const process& p, int in_loop,momentum_configuration<RHP>& ps, 	int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sgqq1_Cut(const process& p, int in_loop,	momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);
template complex<RHP> Sqqg1_Cut(const process& p, int in_loop,	momentum_configuration<RHP>& ps,int inda, int indb, int oeps, RHP mu2);

template complex<RHP> Sggg1s_Rat(const process& p, momentum_configuration<RHP>& ps,int inda, int indb);
template complex<RHP> Sggg1v_Rat(const process& p, momentum_configuration<RHP>& ps,int inda, int indb);
template complex<RHP> Sqgq1_Rat(const process& p, int in_loop,	momentum_configuration<RHP>& ps,int inda, int indb);
template complex<RHP> Sgqq1_Rat(const process& p, int in_loop,	momentum_configuration<RHP>& ps,int inda, int indb);
template complex<RHP> Sqqg1_Rat(const process& p, int in_loop,	momentum_configuration<RHP>& ps,int inda, int indb);


template RVHP li2(RVHP x);
template complex<RVHP> Clog(complex<RVHP> s1, RVHP mu2);
template complex<RVHP> ff(complex<RVHP> z, int in_loop, complex<RVHP> s, RVHP mu2, int oeps);
template Series<complex<RVHP> > SplitS1(const process& p, int in_loop,	momentum_configuration<RVHP>& ps,int inda, int indb, RVHP mu2, int full_cut_rat);
template complex<RVHP> Split1(const process& p, int in_loop,momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2, int full_cut_rat);
template complex<RVHP> Sggg1(const process& p, int in_loop,momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sggg1s(const process& p, momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sggg1v(const process& p, momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sqgq1(const process& p, int in_loop,momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sgqq1(const process& p, int in_loop,momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sqqg1(const process& p, int in_loop,momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);

template complex<RVHP> Sggg1_Cut(const process& p, int in_loop,momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sggg1s_Cut(const process& p, momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sggg1v_Cut(const process& p, momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sqgq1_Cut(const process& p, int in_loop,	momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sgqq1_Cut(const process& p, int in_loop,	momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);
template complex<RVHP> Sqqg1_Cut(const process& p, int in_loop,	momentum_configuration<RVHP>& ps,int inda, int indb, int oeps, RVHP mu2);

template complex<RVHP> Sggg1_Rat(const process& p, int in_loop,	momentum_configuration<RVHP>& ps,int inda, int indb);
template complex<RVHP> Sggg1s_Rat(const process& p, momentum_configuration<RVHP>& ps,int inda, int indb);
template complex<RVHP> Sggg1v_Rat(const process& p, momentum_configuration<RVHP>& ps,int inda, int indb);
template complex<RVHP> Sqgq1_Rat(const process& p, int in_loop,	momentum_configuration<RVHP>& ps,int inda, int indb);
template complex<RVHP> Sgqq1_Rat(const process& p, int in_loop,	momentum_configuration<RVHP>& ps,int inda, int indb);
template complex<RVHP> Sqqg1_Rat(const process& p, int in_loop,	momentum_configuration<RVHP>& ps,int inda, int indb);


#if BH_USE_GMP
template RGMP li2(RGMP x);
template complex<RGMP> Clog(complex<RGMP> s1, RGMP mu2);
template complex<RGMP> ff(complex<RGMP> z, int in_loop, complex<RGMP> s, RGMP mu2, int oeps);
template Series<complex<RGMP> > SplitS1(const process& p, int in_loop,	momentum_configuration<RGMP>& ps,int inda, int indb, RGMP mu2, int full_cut_rat);
template complex<RGMP> Split1(const process& p, int in_loop,momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2, int full_cut_rat);
template complex<RGMP> Sggg1(const process& p, int in_loop,momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sggg1s(const process& p, momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sggg1v(const process& p, momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sqgq1(const process& p, int in_loop,momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sgqq1(const process& p, int in_loop,momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sqqg1(const process& p, int in_loop,momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);

template complex<RGMP> Sggg1_Cut(const process& p, int in_loop,momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sggg1s_Cut(const process& p, momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sggg1v_Cut(const process& p, momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sqgq1_Cut(const process& p, int in_loop,	momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sgqq1_Cut(const process& p, int in_loop,	momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);
template complex<RGMP> Sqqg1_Cut(const process& p, int in_loop,	momentum_configuration<RGMP>& ps,int inda, int indb, int oeps, RGMP mu2);

template complex<RGMP> Sggg1_Rat(const process& p, int in_loop,	momentum_configuration<RGMP>& ps,int inda, int indb);
template complex<RGMP> Sggg1s_Rat(const process& p, momentum_configuration<RGMP>& ps,int inda, int indb);
template complex<RGMP> Sggg1v_Rat(const process& p, momentum_configuration<RGMP>& ps,int inda, int indb);
template complex<RGMP> Sqgq1_Rat(const process& p, int in_loop,	momentum_configuration<RGMP>& ps,int inda, int indb);
template complex<RGMP> Sgqq1_Rat(const process& p, int in_loop,	momentum_configuration<RGMP>& ps,int inda, int indb);
template complex<RGMP> Sqqg1_Rat(const process& p, int in_loop,	momentum_configuration<RGMP>& ps,int inda, int indb);
#endif



}
