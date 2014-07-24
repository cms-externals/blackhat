/*!\file fact_momenta.cpp
\brief Implementation of collinear and soft momentum configurations
*/

// by CFB, last revision 10/20/2008
//
// Long explanation for new version - 9/4/08
// Lance's code for generating coll. kinematics does not exactly obey momentum
// conservation -- off by 10^(-6) generally, which sometimes can cause instabilities
// in the numerical code for the cut parts of the collinear target
// A 2nd version of the collinear kinematics was implemented which is based on the 
// dipole formalism -- see Catani, Seymour, PLB 378, 287 (96), eq. (16)-(18)
// added a function which recomputes the SAME coll config with a different
// coll mass from an EXISTING coll config

#include "BH_typedefs.h"
#include "BH_utilities.h"
#include "factorization.h"
#include "spinor.h"
#include "partitions.h"
#include "mom_conf.h"
#include "qd_suppl.h"


#if BH_USE_GMP
#include "gmp_r.h"
#endif


#include <ctime>
using namespace std;
namespace BH {

#define _MAXTRY 100
// to avoid stack overflow

#define  _VERBOSE 0

// this version follows eqs. (16)-(18) in Catani, Seymour, PLB 378, 287 (96)
// The (n+1)th momentum is p_ij, the (n+2)th momentum is \tilde{p}_k,
// and k is taken to simply be b+1 (i.e. adjacent)
// invmass is a NEW invmass, different from the existing invmass
// and simply inverts Catani-Seymour, p_ij, p_j, and \tilde{p}_k
// are kept fixed, but p_i and p_k are varied to produce a momconf
// with momentum conservation
template <class T> momentum_configuration<T>collkinematicsCSredo(
		momentum_configuration<T>& oldCS, 		
		int a, int b, T z, T invmass){
	int nn = (int) oldCS.n();
	
	// compute the emitter and spectator according to inverted Catani-Seymour (16)-(17)
	complex<T> yijk;
    yijk = invmass/T(2)/(dott(oldCS.p(nn),oldCS.p(nn-1)));
    
	Cmom<T> Kk = (T(1) - yijk)*oldCS.p(nn);
	Cmom<T> Kij = oldCS.p(nn-1) + yijk*oldCS.p(nn);
	T abx = real(Kij.X());
	T aby = real(Kij.Y());
	T abz = real(Kij.Z());
	T abE = real(Kij.E());
	
	// a + b = P is always a null vector with
	// P = a/z, where a = z* ab + eps, b = (1-z)*ab - eps,
	// eps : 2 P.eps = (1-2 z)invmass and eps^2 = -z(1-z)invmass
    T csq = (z*(T(1)-z)-(T(1)-T(2)*z)*(T(1)-T(2)*z)*invmass/T(4)/abx/abx)/ \
		        (aby*aby+abz*abz);
	T aE = z*abE;
	T ax = z*abx - (T(1)-T(2)*z)*invmass/T(2)/abx;
	T ay = z*aby + sqrt(csq*invmass)*abz;
	T az = z*abz - sqrt(csq*invmass)*aby;
	Cmom<T> Ka(aE,ax,ay,az);
	Cmom<T> Kb = Kij - Ka;


    // copy over old kin into new one incl. P_ij and Pktilde except a and k
	vector<Cmom<T> > temparray(nn);
	for (int i = 0; i < nn;i++)
	{	
		
	    if (i == a-1)
	    	temparray.at(i) = Ka;
	    else if (i == b-1)
	    	temparray.at(i) = Kb;
	    else if (i == b)
	    	temparray.at(i) = Kk;
	    else
	    	temparray.at(i) = oldCS.p(i+1);
	}
    
// insert everything into a momentum configuration, append KP at end */
    momentum_configuration<T> tempc(temparray);
    temparray.clear();
    return tempc;
}


// this version follows eqs. (16)-(18) in Catani, Seymour, PLB 378, 287 (96)
// The (n+1)th momentum is p_ij, the (n+2)th momentum is \tilde{p}_k,
// and k is taken to simply be b+1 (i.e. adjacent), which shouldn't matter
// note that in this case invmass should not be too small
template <class T> momentum_configuration<T>collkinematicsCS(int n,
		int a, int b, T z, T invmass){
// following in parts, for the generation of p_i and p_j, Lance's maple code
	T abx = randdouble<T>();
	T aby = randdouble<T>();
	T abz = randdouble<T>();
	T abE = sqrt(invmass + abx*abx + aby*aby + abz*abz);
	Cmom<T> Kab(abE,abx,aby,abz);
// a + b = P is always a null vector with
// P = a/z, where a = z* ab + eps, b = (1-z)*ab - eps,
// eps : 2 P.eps = (1-2 z)invmass and eps^2 = -z(1-z)invmass
	T csq = (z*(T(1)-z)-(T(1)-T(2)*z)*(T(1)-T(2)*z)*invmass/T(4)/abx/abx)/ \
	        (aby*aby+abz*abz);
	if (csq < T(0))
	{
#if _VERBOSE
		_WARNING("Imaginary kinematics, we'll do this again...");
#endif
	    return (collkinematicsCS(n,a,b,z,invmass));
	}
    T aE = z*abE;
    T ax = z*abx - (T(1)-T(2)*z)*invmass/T(2)/abx;
    T ay = z*aby + sqrt(csq*invmass)*abz;
    T az = z*abz - sqrt(csq*invmass)*aby;
    Cmom<T> Ka(aE,ax,ay,az);
    Cmom<T> Kb = Kab - Ka;

// now generate the remaining momenta, except for two randomly chosen ones
// which are computed from "decaying" Ktot.
	int n1 = static_cast<int>(static_cast<double>(rand())/RAND_MAX*(n-4));
	int n2 = static_cast<int>(static_cast<double>(rand())/RAND_MAX*(n-4));
	if (n1 == n2)
		n2 = n1+1;
	if (n2 < n1)
	{	int m = n1;
	    n1 = n2;
	    n2 = m;
	}

	vector<Cmom<T> > temparray(n-4);
	Cmom<T> Ktot = Ka + Kb;
	for (int i = 0; i < n-4;i++)
	{	double signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
		Cmom<T> temp = randmom(T(0), signE > 0 ? 1 : (-1));
		temparray.at(i) = temp;
	    Ktot = Ktot +  temp;
	}

    vector<Cmom<T> > tempv(2);
    tempv = decay(Ktot,0);
    if (tempv.at(0).E() == T(0))
    	return (collkinematicsCS(n,a,b,z,invmass));
    // too many iterations, try again

    temparray.insert(temparray.begin()+n1,tempv.at(0));
    temparray.insert(temparray.begin()+n2,tempv.at(1));

    if (a < b)
    {  temparray.insert(temparray.begin()+a-1,Ka);
       temparray.insert(temparray.begin()+b-1,Kb);
    }
    else
    {  temparray.insert(temparray.begin()+b-1,Kb);
       temparray.insert(temparray.begin()+a-1,Ka);
    }
    
// now compute the emitter and spectator according to Catani-Seymour (16)-(17)
    complex<T> yijk;
    int indk = b;
    if (indk == n)
    	indk = 0;
    yijk = dott(Ka,Kb)/(dott(Ka,Kb)+dott(Kb,temparray.at(indk))+dott(Ka,temparray.at(indk)));
    Cmom<T> KP = Ka + Kb - yijk/(T(1)-yijk)*temparray.at(indk);
    Cmom<T> Knew = T(1)/(T(1)-yijk)*temparray.at(indk);

// insert everything into a momentum configuration, append KP at end */
    momentum_configuration<T> tempc(temparray);
    temparray.clear();
    int nw = tempc.insert(KP);
    int nw2 = tempc.insert(Knew);
    if (nw < 1 || nw2 < 1)
    	_WARNING("Error in generating collinear momenta!\n");
    return tempc;
}


// note, old version
template <class T> momentum_configuration<T>collkinematics(int n,
		int a, int b, T z, T invmass){
// following closely Lance's Maple code
	T abx = randdouble<T>();
	T aby = randdouble<T>();
	T abz = randdouble<T>();
	T abE = sqrt(invmass + abx*abx + aby*aby + abz*abz);
	Cmom<T> Kab(abE,abx,aby,abz);
//	_DEBUG(Kab);
// a + b = P is always a null vector with
// P = a/z, where a = z* ab + eps, b = (1-z)*ab - eps,
// eps : 2 P.eps = (1-2 z)invmass and eps^2 = -z(1-z)invmass
	T csq = (z*(T(1)-z)-(T(1)-T(2)*z)*(T(1)-T(2)*z)*invmass/T(4)/abx/abx)/ \
	        (aby*aby+abz*abz);
	if (csq < 0.)
	{
#if _VERBOSE
		_WARNING("Imaginary kinematics, we'll do this again...");
#endif
	    return (collkinematics(n,a,b,z,invmass));
	}
    T aE = z*abE;
    T ax = z*abx - (T(1)-T(2)*z)*invmass/T(2)/abx;
    T ay = z*aby + sqrt(csq*invmass)*abz;
    T az = z*abz - sqrt(csq*invmass)*aby;
    Cmom<T> Ka(aE,ax,ay,az);
//    _DEBUG(Ka);
    Cmom<T> Kb = Kab - Ka;
//    _DEBUG(Kb);
	Cmom<T> KP = Ka/z;

// now generate the remaining momenta, except for two randomly chosen ones
// which are computed from "decaying" Ktot.
	int n1 = static_cast<int>(static_cast<double>(rand())/RAND_MAX*(n-4));
	int n2 = static_cast<int>(static_cast<double>(rand())/RAND_MAX*(n-4));
	if (n1 == n2)
		n2 = n1+1;
	if (n2 < n1)
	{	int m = n1;
	    n1 = n2;
	    n2 = m;
	}

	vector<Cmom<T> > temparray(n-4);
	Cmom<T> Ktot = Ka + Kb;
	for (int i = 0; i < n-4;i++)
	{	double signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
		Cmom<T> temp = randmom(T(0), signE > 0 ? 1 : (-1));
		temparray.at(i) = temp;
	    Ktot = Ktot +  temp;
	}

    vector<Cmom<T> > tempv(2);
    tempv = decay(Ktot,0);
    if (tempv.at(0).E() == T(0))
    	return (collkinematics(n,a,b,z,invmass));
    // too many iterations, try again

    temparray.insert(temparray.begin()+n1,tempv.at(0));
    temparray.insert(temparray.begin()+n2,tempv.at(1));

    if (a < b)
    {  temparray.insert(temparray.begin()+a-1,Ka);
       temparray.insert(temparray.begin()+b-1,Kb);
    }
    else
    {  temparray.insert(temparray.begin()+b-1,Kb);
       temparray.insert(temparray.begin()+a-1,Ka);
    }

// insert everything into a momentum configuration, append KP at end */
    momentum_configuration<T> tempc(temparray);
    temparray.clear();
    int nw = tempc.insert(KP);
    if (nw < 1)
    	_WARNING("Error in generating collinear momenta!\n");
    return tempc;
}


template <class T> complex<T> dott(Cmom<T> Ki, Cmom<T> Kj)
{
	return( Ki.E()*Kj.E()-Ki.X()*Kj.X()- Ki.Y()*Kj.Y() - Ki.Z()*Kj.Z() );
}


template <class T> momentum_configuration<T> softkinematics(int n,
		int a, int b, int c, T invmass){
	// following closely Lance's Maple code
	double signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
	Cmom<T> Ka = randmom(T(0),signE > 0 ? 1 : (-1));
	signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
	Cmom<T> Kb = randmom(T(0),signE > 0 ? 1 : (-1));
	Cmom<T> Kt = Ka + Kb;
	T Ktsquare = real(Kt.E()*Kt.E()-Kt.X()*Kt.X()-Kt.Y()*Kt.Y()-Kt.Z()*Kt.Z());

	T AA = real(Kt.E()*Kt.E())-real(Kt.Y()*Kt.Y());
	T BB = real(Kt.E())*(Ktsquare+invmass);
	T CC = (Ktsquare+invmass)*(Ktsquare+invmass)/T(4) + invmass*real(Kt.Y()*Kt.Y());
	T discrim = BB*BB - T(4)*AA*CC;
	if (discrim < 0.)
	{
#if _VERBOSE
		_WARNING("Imaginary kinematics (loop 1), we'll do this again...");
#endif
	     return(softkinematics(n,a,b,c,invmass));
	}
	if (fabs(AA) < 0.00000000001)
	{
#if _VERBOSE
		_WARNING("Division by almost 0, we'll do this again...");
#endif
	    return(softkinematics(n,a,b,c,invmass));
	}
	T EE = (-BB + sqrt(discrim))/T(2)/AA;
	Cmom<T> Kc(-Kt.E()-EE,-Kt.X(),-Kt.Y()- sqrt(EE*EE-invmass),-Kt.Z());

    // may have to use the other sign of the square root if momentum not
    // massless

	if (fabs(real(Kc.E()*Kc.E()-Kc.X()*Kc.X()-Kc.Y()*Kc.Y()-Kc.Z()*Kc.Z())/T(2)) > 0.00000000001)
	{
		EE = (-BB - sqrt(discrim))/T(2)/AA;
	    Kc = Cmom<T>(-Kt.E()-EE,-Kt.X(),-Kt.Y()- sqrt(EE*EE-invmass),-Kt.Z());
	}

	if (fabs(real(Kc.E()*Kc.E()-Kc.X()*Kc.X()-Kc.Y()*Kc.Y()-Kc.Z()*Kc.Z())/T(2)) > 0.00000000001)
	{
#if _VERBOSE
		_WARNING("Imaginary kinematics (loop 2), we'll do this again...");
#endif
        return(softkinematics(n,a,b,c,invmass));
	}
    Cmom<T> KP(-EE,T(0),-sqrt(EE*EE-invmass),T(0));

	// now generate the remaining momenta, except for two randomly chosen ones
	// which are computed from "decaying" Ktot.
	int n1 = static_cast<int>(static_cast<double>(rand())/RAND_MAX*(n-5));
	int n2 = static_cast<int>(static_cast<double>(rand())/RAND_MAX*(n-5));
	if (n1 == n2)
		n2 = n1+1;
	if (n2 < n1)
	{	int m = n1;
	    n1 = n2;
	    n2 = m;
	}

	vector<Cmom<T> > temparray(n-5);
	Cmom<T> Ktot = Ka + Kb + Kc;
	for (int i = 0; i < n-5;i++)
	{	double signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
		Cmom<T> temp = randmom(T(0), signE > 0 ? 1 : (-1));
		temparray.at(i) = temp;
	    Ktot = Ktot +  temp;
	}

	vector<Cmom<T> > tempv(2);
	tempv = decay(Ktot,0);
    if (tempv.at(0).E() == T(0))
    	return (softkinematics(n,a,b,c,invmass));
    // too many iterations, try again

	temparray.insert(temparray.begin()+n1,tempv.at(0));
	temparray.insert(temparray.begin()+n2,tempv.at(1));

	if (a < b < c)
	{   temparray.insert(temparray.begin()+a-1,Ka);
        temparray.insert(temparray.begin()+b-1,Kb);
        temparray.insert(temparray.begin()+c-1,Kc);
	}
	else if (a < c < b)
	{
		temparray.insert(temparray.begin()+a-1,Ka);
		temparray.insert(temparray.begin()+c-1,Kc);
		temparray.insert(temparray.begin()+b-1,Kb);
	}
	else if (c < a < b)
	{
		temparray.insert(temparray.begin()+c-1,Kc);
		temparray.insert(temparray.begin()+a-1,Ka);
		temparray.insert(temparray.begin()+b-1,Kb);
	}
	else if (c < b < a)
	{
		temparray.insert(temparray.begin()+c-1,Kc);
		temparray.insert(temparray.begin()+b-1,Kb);
		temparray.insert(temparray.begin()+a-1,Ka);
	}
	else if (b < a < c)
	{
		temparray.insert(temparray.begin()+b-1,Kb);
		temparray.insert(temparray.begin()+a-1,Ka);
		temparray.insert(temparray.begin()+c-1,Kc);
	}
	else   // b < c < a
	{
		temparray.insert(temparray.begin()+b-1,Kb);
		temparray.insert(temparray.begin()+c-1,Kc);
		temparray.insert(temparray.begin()+a-1,Ka);
	}

	// insert everything into a momentum configuration, append KP at end
	momentum_configuration<T> tempc(temparray);
	temparray.clear();
	int nw = tempc.insert(KP);
	if (nw < 1)
	   	_WARNING("Error in generating soft momenta!\n");
	return tempc;
}


template <class T> momentum_configuration<T> soft4kinematics(int n,
		int a, int b, int c, int d, T invmass){
	// following closely Lance's Maple code
	double signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
	Cmom<T> Ka = randmom(T(0),signE > 0 ? 1 : (-1));
	signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
	Cmom<T> Kb = randmom(T(0),signE > 0 ? 1 : (-1));
	signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
	Cmom<T> Kc = randmom(T(0),signE > 0 ? 1 : (-1));
	Cmom<T> Kt = Ka + Kb + Kc;

	T Ktsquare = real(Kt.E()*Kt.E()-Kt.X()*Kt.X()-Kt.Y()*Kt.Y()-Kt.Z()*Kt.Z());
	T AA = real(Kt.E()*Kt.E()-Kt.Y()*Kt.Y());
	T BB = real(Kt.E())*(Ktsquare+invmass);
	T CC = (Ktsquare+invmass)*(Ktsquare+invmass)/T(4) + invmass*real(Kt.Y()*Kt.Y());
	T discrim = BB*BB - T(4)*AA*CC;
	if (discrim < 0.)
	{
#if _VERBOSE
		_WARNING("Imaginary kinematics (loop 1), we'll do this again...");
#endif
	     return(soft4kinematics(n,a,b,c,d,invmass));
	}
	if (fabs(AA) < 0.00000000001)
	{
#if _VERBOSE
		_WARNING("Division by almost 0, we'll do this again...");
#endif
	    return(soft4kinematics(n,a,b,c,d,invmass));
	}
	T EE = (-BB + sqrt(discrim))/T(2)/AA;
	Cmom<T> Kd(-Kt.E()-EE,-Kt.X(),-Kt.Y()- sqrt(EE*EE-invmass),-Kt.Z());

    // may have to use the other sign of the square root if momentum not
    // massless
	if (fabs(real(Kd.E()*Kd.E()-Kd.X()*Kd.X()-Kd.Y()*Kd.Y()-Kd.Z()*Kd.Z())/T(2)) > 0.00000000001)
	{
		EE = (-BB - sqrt(discrim))/T(2)/AA;
	    Kd = Cmom<T>(-Kt.E()-EE,-Kt.X(),-Kt.Y()- sqrt(EE*EE-invmass),-Kt.Z());
	}

	if (fabs(real(Kd.E()*Kd.E()-Kd.X()*Kd.X()-Kd.Y()*Kd.Y()-Kd.Z()*Kd.Z())/T(2)) > 0.00000000001)
	{   _WARNING("Imaginary kinematics (loop 2), we'll do this again...");
        return(soft4kinematics(n,a,b,c,d,invmass));
	}
    Cmom<T> KP(-EE,T(0),-sqrt(EE*EE-invmass),T(0));

	// now generate the remaining momenta, except for two randomly chosen ones
	// which are computed from "decaying" Ktot.
	int n1 = static_cast<int>(static_cast<double>(rand())/RAND_MAX*(n-6));
	int n2 = static_cast<int>(static_cast<double>(rand())/RAND_MAX*(n-6));
	if (n1 == n2)
		n2 = n1+1;
	if (n2 < n1)
	{	int m = n1;
	    n1 = n2;
	    n2 = m;
	}

	vector<Cmom<T>  > temparray(n-6);
	Cmom<T> Ktot = Ka + Kb + Kc + Kd;
	for (int i = 0; i < n-6;i++)
	{	double signE = static_cast<double>(rand())/(RAND_MAX)-0.5;
		Cmom<T> temp = randmom(T(0), signE > 0 ? 1 : (-1));
		temparray.at(i) = temp;
	    Ktot = Ktot +  temp;
	}

	vector<Cmom<T> > tempv(2);
	tempv = decay(Ktot,0);
    if (tempv.at(0).E() == T(0))
    	return (soft4kinematics(n,a,b,c,d,invmass));

	temparray.insert(temparray.begin()+n1,tempv.at(0));
	temparray.insert(temparray.begin()+n2,tempv.at(1));

	// order the momenta first
	bool swapped = true;
	vector<Cmom<T> > Ki;
	vector<int> ii;
	Ki.push_back(Ka);  ii.push_back(a);
	Ki.push_back(Kb);  ii.push_back(b);
	Ki.push_back(Kc);  ii.push_back(c);
	Ki.push_back(Kd);  ii.push_back(d);

	while(swapped)
	{
		swapped = false;
		for(int i = 0; i < 3; i++)
		{
			if(ii.at(i) > ii.at(i+1))
			{
				// swap the two elements
				swapped = true;
				int temp = ii.at(i);
				ii.at(i) = ii.at(i+1);
				ii.at(i+1) = temp;
				Cmom<T> tempM = Ki.at(i);
				Ki.at(i) = Ki.at(i+1);
				Ki.at(i+1) = tempM;
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{   temparray.insert(temparray.begin()+ii.at(i)-1,Ki.at(i));
	}

	// insert everything into a momentum configuration, append KP at end
	momentum_configuration<T> tempc(temparray);
	temparray.clear();
	int nw = tempc.insert(KP);
	if (nw < 1)
	   	_WARNING("Error in generating 4 soft momenta!\n");
	return tempc;
}


template <class T> vector<Cmom<T> > decay(Cmom<T> Ktot, int trial){
	vector<Cmom<T> > decayed(2);
	if (trial > _MAXTRY)
	{
		Cmom<T> K0(T(0),T(0),T(0),T(0));
		decayed.at(0) = K0;
		decayed.at(1) = K0;
		return decayed;
	}
	T auxy = T(2)*randdouble<T>()-T(1);
	T auxz = T(2)*randdouble<T>()-T(1);
	T alpha = auxy*auxy + auxz*auxz;
	T beta = auxy*real(Ktot.Y()) + auxz*real(Ktot.Z());
	T square = real(Ktot.E()*Ktot.E()-Ktot.X()*Ktot.X()-Ktot.Y()*Ktot.Y()-Ktot.Z()*Ktot.Z());
	T gamma = beta - 1/T(2)*(square);
	T delta = real(Ktot.E()*Ktot.E()) - real(Ktot.X()*Ktot.X());
	T discrim = -alpha*delta + gamma*gamma;
	if (discrim < 0.)
	{
#if _VERBOSE
		_WARNING("Imaginary decay kinematics -- will try decay again!\n");
#endif
	   return decay(Ktot,trial+1);
	}
	_DEBUG("Real decay achieved!\n");
	T auxx = (real(Ktot.X())*gamma + real(Ktot.E())*sqrt(discrim))/delta;
	T auxE = (auxx*real(Ktot.X()) + auxy*real(Ktot.Y()) + auxz*real(Ktot.Z()) - \
	     1/T(2)*(square))/real(Ktot.E());
    // 7/14/2008 -- check for dangerous configuration (instability in FHZ)
	if (fabs(auxx-auxE) < 0.01)
	{
#if _VERBOSE
		_WARNING("Dangerous decay configuration, we'll decay again...");
#endif
	    return (decay(Ktot,trial+1));
	}
	Cmom<T> Ka(auxE,auxx,auxy,auxz);
	decayed.at(0) = Ka;
	decayed.at(1) = - Ktot - Ka;

	return decayed;
}

template <class T> Cmom<T> randmom(T mass, short signE){
// could be improved...
	_DEBUG(signE);
	T x = T(2)*randdouble<T>()-T(1);
	T y = T(2)*randdouble<T>()-T(1);
	T z = T(2)*randdouble<T>()-T(1);
	if (mass == T(0))
		return Cmom<T>(T(signE)*sqrt(x*x+y*y+z*z),x,y,z);
	_WARNING("Error, massive momenta not yet implemented!\n");
	return Cmom<T>(T(0),T(0),T(0),T(0));
}

template <class T> T randdouble(void){
	return static_cast<T>(static_cast<T>(rand())/(RAND_MAX+T(1)));
}

// explicit instantiations
template momentum_configuration<R> collkinematics(int n, int a, int b,
		R zfrac, R invmass);
template momentum_configuration<R> softkinematics(int n, int a, int b, int c,
		R invmass);
template momentum_configuration<R> soft4kinematics(int n, int a, int b, int c, int d,
		R invmass);
template vector<Cmom<R> > decay(Cmom<R> Ktot, int);
template Cmom<R> randmom(R mass, short signE);
template momentum_configuration<R> collkinematicsCS(int n,
		int a, int b, R z, R invmass);
template complex<R> dott(Cmom<R> Ki, Cmom<R> Kj);
template R randdouble();
template momentum_configuration<R> collkinematicsCSredo(
		momentum_configuration<R>& oldCS,
		int a, int b, R z, R invmass);

template momentum_configuration<RHP> collkinematics(int n, int a, int b,
		RHP zfrac, RHP invmass);
template momentum_configuration<RHP> softkinematics(int n, int a, int b, int c,
		RHP invmass);
template momentum_configuration<RHP> soft4kinematics(int n, int a, int b, int c, int d,
		RHP invmass);
template vector<Cmom<RHP> > decay(Cmom<RHP> Ktot, int);
template Cmom<RHP> randmom(RHP mass, short signE);
template RHP randdouble();
template momentum_configuration<RHP> collkinematicsCS(int n,
		int a, int b, RHP z, RHP invmass);
template complex<RHP> dott(Cmom<RHP> Ki, Cmom<RHP> Kj);
template momentum_configuration<RHP> collkinematicsCSredo(
		momentum_configuration<RHP>& oldCS,
		int a, int b, RHP z, RHP invmass);

template momentum_configuration<RVHP> collkinematics(int n, int a, int b,
		RVHP zfrac, RVHP invmass);
template momentum_configuration<RVHP> softkinematics(int n, int a, int b, int c,
		RVHP invmass);
template momentum_configuration<RVHP> soft4kinematics(int n, int a, int b, int c, int d,
		RVHP invmass);
template vector<Cmom<RVHP> > decay(Cmom<RVHP> Ktot, int);
template Cmom<RVHP> randmom(RVHP mass, short signE);
template RVHP randdouble();
template momentum_configuration<RVHP> collkinematicsCS(int n,
		int a, int b, RVHP z, RVHP invmass);
template complex<RVHP> dott(Cmom<RVHP> Ki, Cmom<RVHP> Kj);
template momentum_configuration<RVHP> collkinematicsCSredo(
		momentum_configuration<RVHP>& oldCS,
		int a, int b, RVHP z, RVHP invmass);
#if BH_USE_GMP
template momentum_configuration<RGMP> collkinematics(int n, int a, int b,
		RGMP zfrac, RGMP invmass);
template momentum_configuration<RGMP> softkinematics(int n, int a, int b, int c,
		RGMP invmass);
template momentum_configuration<RGMP> soft4kinematics(int n, int a, int b, int c, int d,
		RGMP invmass);
template vector<Cmom<RGMP> > decay(Cmom<RGMP> Ktot, int);
template Cmom<RGMP> randmom(RGMP mass, short signE);
template RGMP randdouble();
template momentum_configuration<RGMP> collkinematicsCS(int n,
		int a, int b, RGMP z, RGMP invmass);
template complex<RGMP> dott(Cmom<RGMP> Ki, Cmom<RGMP> Kj);
template momentum_configuration<RGMP> collkinematicsCSredo(
		momentum_configuration<RGMP>& oldCS,
		int a, int b, RGMP z, RGMP invmass);
#endif


}
