/*
 * extend.cpp
 *
 *  Created on: Aug 8, 2008
 *      Author: daniel apdapted from David's extend.cc
 */



/*  extend.cc  */

/* David A. Kosower, June 3, 2008 */

/* Extend precision of momentum configurations.

   To extend the precision of an n-particle configuration, first recompute
   the energy components of those of the 1st n-2 momenta which are supposed
   to be massless (as judged by falling below a threshold).  Then redecay
   the sum of these momenta in the center-of-mass frame of the original
   last 2 momenta, and boost back to the frame of the n-2 momenta.

*/

#include "mom_conf.h"
#include "BH_debug.h"


using namespace std;

namespace BH {

Cmom<RVHP> to_higher_precision(const Cmom<RHP>& P_HP){
	return Cmom<RVHP>(P_HP.E(),P_HP.X(),P_HP.Y(),P_HP.Z());
}
Cmom<RHP> to_higher_precision(const Cmom<R>& P_HP){
	return Cmom<RHP>(P_HP.E(),P_HP.X(),P_HP.Y(),P_HP.Z());
}


// Is it smaller in absolute value than zeroThreshold?
// Be careful that zT^2 is representable!
const double zeroThreshold = 1e-15;
const RHP zeroThresholdHP = RHP("1e-30");  // What should this value be???
const RVHP zeroThresholdVHP = RVHP("1e-60");  // What should this value be???

inline bool IsZero(const complex<double>& value)
{return(real(value)*real(value)+imag(value)*imag(value) <
        zeroThreshold*zeroThreshold);}

inline bool IsZero(const complex<RHP>& value)
{return(real(value)*real(value)+imag(value)*imag(value) <
        zeroThresholdHP*zeroThresholdHP);}


inline bool IsZero(const R& value)
{return(value < zeroThreshold);}

inline bool IsZero(const RHP& value)
{return(value < zeroThresholdHP);}

inline bool IsZero(const RVHP& value)
{return(value < zeroThresholdVHP);}

template<class T> inline T Sign(const complex<T>& value)
{if (IsZero(imag(value)) and real(value) < 0) return T(-1); else return T(1);}

/* Transform the vector "p" from the frame in which "transform" is at rest
   (if it is positive-mass) or has zero energy component (if it is
   negative-mass) to one in which it has the given form. */
template<class Te,class T0> static momentum<complex<Te> >
  Boost(const momentum<complex<T0> >& p,
        const momentum<complex<Te> >& transform)
{complex<Te> massSq = transform*transform;
 if (IsZero(imag(massSq)) and real(massSq) < 0)
    {complex<Te> mass = sqrt(-transform*transform);
    momentum<complex<Te> > base(Te(0),transform.X(),transform.Y(),transform.Z());
    return Boost(p,momentum<complex<Te> >(base*base/mass,Te(0),Te(0),Te(0))
                 -transform.E()/mass*base);}
 Te sign = Sign(transform.E());
 complex<Te> mass = sqrt(transform*transform);
 complex<Te> gamma = sign * transform.E()/mass;
 complex<Te> oneOverOnePlusGamma = Te(1)/(Te(1)+gamma);
 momentum<complex<Te> > boost(Te(0),transform.X(),transform.Y(),transform.Z());
 boost *= sign/mass;
 complex<Te> gammaS = -boost*p;
 complex<Te> boostMix = p.E() + oneOverOnePlusGamma * gammaS;
#if 0
 cout << gamma << endl;
 cout << "p: " << p << endl;
 cout << "transform: " << transform << endl;
 cout << "boost: " << boost << endl;
 cout << gammaS << endl;
 cout << boostMix << endl;
#endif
 momentum<complex<Te> > result(gamma * p.E() + gammaS,
                               p.X()+boostMix*boost.X(),
                               p.Y()+boostMix*boost.Y(),
                               p.Z()+boostMix*boost.Z());
 return result;
}

// this is not the original version, unlike the original version, it returns a mom_conf vith only as many entries as in the index vector. The new vectors are set ind[0] -> mc[1] ...

template<class Te,class T0> momentum_configuration<Te>
  extend(momentum_configuration<T0>& k,
         const vector<int>& indices /* of the momenta within "k" */,
         const Te& dummy /* To signal target type */)
{int n = indices.size();
 // Ensure proper rounding to avoid loss of precision
 unsigned int old_cw;
 fpu_fix_start(&old_cw);
 momentum<complex<Te> > sum;
 vector<Cmom<Te> > result;
 // vector indexing is from 0, mom_conf from 1:
 const int Offset = -1;
 // Although reserve is supposed to leave existing vector elements alone,
 // it seems to zero them out, so reserve space ahead of time...
 // and we need to use resize to reset the size, not just reserve to get space

 for (int i = 0;  i < n-2;  i += 1){
	 int j = indices[i];
    momentum<complex<Te> > extended(complex<Te>(0),complex<Te>(k[j].X()),complex<Te>(k[j].Y()),complex<Te>(k[j].Z()));
    if (IsZero(k.m2(j)))
      extended = momentum<complex<Te> >
        (Te(Sign(k[j].E()))*sqrt(-extended*extended),
         extended.X(),extended.Y(),extended.Z());
    else  {// Just use original mass
      _MESSAGE3(j," is massive:",k.m2(j));extended = momentum<complex<Te> >(complex<Te>(k[j].E()),extended.X(),extended.Y(),
                                        extended.Z());}
    result.push_back( extended);
    sum += extended;
    }

 // Last two momenta
 int j = indices[n-2];
 momentum<complex<Te> > p1(k[j].E(),k[j].X(),k[j].Y(),k[j].Z());
 bool massless1 = IsZero(k[j]*k[j]);
 j = indices[n-1];
 momentum<complex<Te> > p2(k[j].E(),k[j].X(),k[j].Y(),k[j].Z());
 bool massless2 = IsZero(k[j]*k[j]);

 momentum<complex<Te> > boost(sum.E(),-sum.X(),-sum.Y(),-sum.Z());

 p1 = Boost(p1,boost); // to center of mass
 p2 = Boost(p2,boost);

 // Make sure that both are massless
 if (massless1)
    {Te sign = Sign(p1.E());
    p1 = momentum<complex<Te> >(Te(0),p1.X(),p1.Y(),p1.Z());
    p1 = momentum<complex<Te> >(Sign(p1.E())*sqrt(-p1*p1),p1.X(),p1.Y(),p1.Z());}
 if (massless2)
    {Te sign = Sign(p2.E());
    p2 = momentum<complex<Te> >(Te(0),p2.X(),p2.Y(),p2.Z());
    p2 = momentum<complex<Te> >(Sign(p2.E())*sqrt(-p2*p2),p2.X(),p2.Y(),p2.Z());}

 // and that momentum is conserved to new precision -- must be done 2nd
 complex<Te> rescale = sqrt(sum*sum/((p1+p2)*(p1+p2)));
 p1 *= rescale;
 p2 = momentum<complex<Te> >(rescale*p2.E(),-p1.X(),-p1.Y(),-p1.Z());

 // Boost back
 p1 = Boost(p1,sum);
 p2 = Boost(p2,sum);

 j = indices[n-2];
 result.push_back( p1);
 //    cout << "ex " << j << ": " << result[j] << endl;
 j = indices[n-1];
 result.push_back( p2);
 //    cout << "ex " << j << ": " << result[j] << endl;

#if 1
 momentum<complex<Te> > the_sum(complex<Te>(0,0),complex<Te>(0,0),complex<Te>(0,0),complex<Te>(0,0));
 for (int j = 1;  j <= indices.size() ;  j += 1){
    the_sum=the_sum+result[j+Offset].P();


	 std::cout << "before " << j << ": " << k.p(j) << " mass: " << k.m2(j) << std::endl;
	 std::cout << "after " << j << ": " <<  result[j+Offset]  << " mass: " << result[j+Offset]*result[j+Offset] << std::endl;
	 std::cout << "diff " << j << ": " << result[j+Offset] - to_higher_precision(k.p(j)) << std::endl;
 }
std::cout << "sum: " << the_sum << std::endl;
#endif

 return momentum_configuration<Te>(result);
 fpu_fix_end(&old_cw);
}

// this is the original version it returns a mom_conf vith as many entries as the maximum value in the index vector. For example if ind = {1,2,3,4,12}, it will return a mom_conf with 12 vectors
template<class Te,class T0> momentum_configuration<Te>
  extend_and_keep(momentum_configuration<T0>& k,
         const vector<int>& indices /* of the momenta within "k" */,
         const Te& dummy /* To signal target type */)
{int n = indices.size();
 // Ensure proper rounding to avoid loss of precision
 unsigned int old_cw;
 fpu_fix_start(&old_cw);
 momentum<complex<Te> > sum;
 vector<Cmom<Te> > result;
 // vector indexing is from 0, mom_conf from 1:
 const int Offset = -1;
 // Although reserve is supposed to leave existing vector elements alone,
 // it seems to zero them out, so reserve space ahead of time...
 // and we need to use resize to reset the size, not just reserve to get space
 int max = 0;
 for (int i = 0;  i < n;  i += 1){
	 int j = indices[i];
	 if (j > max) max = j; /* Need [j-1] to be valid */
}
 result.resize(max);

 for (int i = 0;  i < n-2;  i += 1)
    {int j = indices[i];
    momentum<complex<Te> > extended(Te(0),k[j].X(),k[j].Y(),k[j].Z());
    if (IsZero(k.m2(j)))
      extended = momentum<complex<Te> >
        (Te(Sign(k[j].E()))*sqrt(-extended*extended),
         extended.X(),extended.Y(),extended.Z());
    else // Just use original mass
      extended = momentum<complex<Te> >(k[j].E(),extended.X(),extended.Y(),
                                        extended.Z());
    result[j+Offset] = extended;
    sum += extended;}

 // Last two momenta
 int j = indices[n-2];
 momentum<complex<Te> > p1(k[j].E(),k[j].X(),k[j].Y(),k[j].Z());
 bool massless1 = IsZero(k[j]*k[j]);
 j = indices[n-1];
 momentum<complex<Te> > p2(k[j].E(),k[j].X(),k[j].Y(),k[j].Z());
 bool massless2 = IsZero(k[j]*k[j]);

 momentum<complex<Te> > boost(sum.E(),-sum.X(),-sum.Y(),-sum.Z());

 p1 = Boost(p1,boost); // to center of mass
 p2 = Boost(p2,boost);

 // Make sure that both are massless
 if (massless1)
    {Te sign = Sign(p1.E());
    p1 = momentum<complex<Te> >(Te(0),p1.X(),p1.Y(),p1.Z());
    p1 = momentum<complex<Te> >(Sign(p1.E())*sqrt(-p1*p1),p1.X(),p1.Y(),p1.Z());}
 if (massless2)
    {Te sign = Sign(p2.E());
    p2 = momentum<complex<Te> >(Te(0),p2.X(),p2.Y(),p2.Z());
    p2 = momentum<complex<Te> >(Sign(p2.E())*sqrt(-p2*p2),p2.X(),p2.Y(),p2.Z());}

 // and that momentum is conserved to new precision -- must be done 2nd
 complex<Te> rescale = sqrt(sum*sum/((p1+p2)*(p1+p2)));
 p1 *= rescale;
 p2 = momentum<complex<Te> >(rescale*p2.E(),-p1.X(),-p1.Y(),-p1.Z());

 // Boost back
 p1 = Boost(p1,sum);
 p2 = Boost(p2,sum);

 j = indices[n-2];
 result[j+Offset] = p1;
 //    cout << "ex " << j << ": " << result[j] << endl;
 j = indices[n-1];
 result[j+Offset] = p2;
 //    cout << "ex " << j << ": " << result[j] << endl;

#if 0
 for (int j = 1;  j <= 6;  j += 1)
    cout << "ex " << j << ": " << result[j+Offset] << endl;
 momentum_configuration<Te> test(result);
 cout << test.n() << "; " << result.size() << endl;
 for (int j = 1;  j <= 6;  j += 1)
   cout << "ex " << j << ": " << test[j] << endl;
#endif
 fpu_fix_end(&old_cw);
 // Need a way to ensure that the final indices are the same as the original
 // ones -- does this do it?
 return momentum_configuration<Te>(result);
}


template <class T_low,class T_high> momentum_configuration<T_high> extend_daniel(const momentum_configuration<T_low>& mc,const vector<int>& ind){
	CPU_FIX CF;
	vector<Cmom<T_high> > momenta;
	momentum<complex<T_high> > sum;
	size_t n=ind.size();
	for (size_t i=0;i<n-2;i++){
		momenta.push_back(Cmom<T_high>(lambdat<T_high>(mc.Lt(ind[i])),lambda<T_high>(mc.L(ind[i]))));
		sum+=momenta.back().P();
	}
//_PRINT(sum);
	smatrix<T_high> P(-sum);

	lambda<T_high> LA1(mc.L(ind[n-2]));
	lambda<T_high> LA2(mc.L(ind[n-1]));

	lambdat<T_high> Z1(complex<T_high>(1,0),complex<T_high>(0,0));
	lambdat<T_high> Z2(complex<T_high>(0,0),complex<T_high>(1,0));
	lambda<T_high> Z3(complex<T_high>(1,0),complex<T_high>(0,0));
	lambda<T_high> Z4(complex<T_high>(0,0),complex<T_high>(1,0));


	complex<T_high> a1=-LA1*Z2;
	complex<T_high> a2=LA1*Z1;
	complex<T_high> A1=-LA2*Z2;
	complex<T_high> A2=LA2*Z1;

	complex<T_high> P1=Z4*P*Z2;
	complex<T_high> P2=-Z4*P*Z1;
	complex<T_high> P3=-Z3*P*Z2;
	complex<T_high> P4=Z3*P*Z1;


	complex<T_high> b1=-((-(A2*P1) + A1*P3)/(-(A1*a2) + a1*A2));
	complex<T_high> b2=-((-(A2*P2) + A1*P4)/(-(A1*a2) + a1*A2));
	complex<T_high> B1=-((-(a2*P1) + a1*P3)/(A1*a2 - a1*A2));
	complex<T_high> B2=-((-(a2*P2) + a1*P4)/(A1*a2 - a1*A2));

//	_PRINT(a1);
//	_PRINT(a2);
//	_PRINT(A1);
//	_PRINT(A2);
//	_PRINT(P1);
//	_PRINT(P2);
//	_PRINT(P3);
//	_PRINT(P4);
//	_PRINT(b1);
//	_PRINT(b2);
//	_PRINT(B1);
//	_PRINT(B2);

	momenta.push_back(Cmom<T_high>(lambdat<T_high>(b1,b2),lambda<T_high>(a1,a2)));
	momenta.push_back(Cmom<T_high>(lambdat<T_high>(B1,B2),lambda<T_high>(A1,A2)));

//	_PRINT(momenta[4]*momenta[4]);
//	_PRINT(momenta[5]*momenta[5]);

//	_PRINT(momenta[0]+momenta[1]+momenta[2]+momenta[3]+momenta[4]+momenta[5]);
	return momentum_configuration<T_high>(momenta);
}
// Fernando and harald's extend
template <class T_low,class T_high> momentum_configuration<T_high> extend_real(const momentum_configuration<T_low>& mc,const vector<int>& ind){
	lambdat<T_high> Z1(complex<T_high>(1,0),complex<T_high>(0,0));
	lambdat<T_high> Z2(complex<T_high>(0,0),complex<T_high>(1,0));
	lambda<T_high> Z3(complex<T_high>(1,0),complex<T_high>(0,0));
	lambda<T_high> Z4(complex<T_high>(0,0),complex<T_high>(1,0));

//	return extend_daniel<T_low,T_high>(mc,ind);
	CPU_FIX CF;
	vector<Cmom<T_high> > momenta;
	momentum<complex<T_high> > sum;
	size_t n=ind.size();
	for (size_t i=0;i<n-2;i++){
		// to avoid issues with new PHASE conventions, it is assumed momenta is real, and built the extended momenta out of lambdat
#if _OLD_PHASE_CONVENTION
		momenta.push_back(Cmom<T_high>(lambdat<T_high>(mc.Lt(ind[i])),lambda<T_high>(mc.L(ind[i]))));
#else
		// check for complex momenta
#define SMALL_COMPLEX_PART 1e-10
		if(abs(imag(mc.p(ind[i]).E()))>SMALL_COMPLEX_PART ||
				abs(imag(mc.p(ind[i]).X()))>SMALL_COMPLEX_PART ||
				abs(imag(mc.p(ind[i]).Y()))>SMALL_COMPLEX_PART ||
				abs(imag(mc.p(ind[i]).Z()))>SMALL_COMPLEX_PART){
			_WARNING("CAREFUL: Using extend_real for a momentum with large complex part");
			momenta.push_back(Cmom<T_high>(lambdat<T_high>(mc.Lt(ind[i])),lambda<T_high>(mc.L(ind[i]))));
		}
		else{
			if(real(mc.p(ind[i]).E())<=0){
				momenta.push_back(Cmom<T_high>(lambdat<T_high>(mc.Lt(ind[i])),
						lambda<T_high>(-conj(Z4*lambdat<T_high>(mc.Lt(ind[i]))),conj(Z3*lambdat<T_high>(mc.Lt(ind[i]))))));
			}
			else{
				momenta.push_back(Cmom<T_high>(lambdat<T_high>(mc.Lt(ind[i])),
						lambda<T_high>(conj(Z4*lambdat<T_high>(mc.Lt(ind[i]))),-conj(Z3*lambdat<T_high>(mc.Lt(ind[i]))))));
			}
		}
#endif
		sum+=momenta.back().P();
	}

	// completed k_{n-1}
#if _OLD_PHASE_CONVENTION
	Cmom<T_high> knm1c(lambdat<T_high>(mc.Lt(ind[n-2])),lambda<T_high>(mc.L(ind[n-2])));
#else
	Cmom<T_high> knm1c;
	if(abs(imag(mc.p(ind[n-2]).E()))>SMALL_COMPLEX_PART ||
			abs(imag(mc.p(ind[n-2]).X()))>SMALL_COMPLEX_PART ||
			abs(imag(mc.p(ind[n-2]).Y()))>SMALL_COMPLEX_PART ||
			abs(imag(mc.p(ind[n-2]).Z()))>SMALL_COMPLEX_PART){
		_WARNING("CAREFUL: Using extend_real for a momentum with large complex part");
		knm1c=Cmom<T_high>(lambdat<T_high>(mc.Lt(ind[n-2])),lambda<T_high>(mc.L(ind[n-2])));
	}
	else{
		if(real(mc.p(ind[n-2]).E())<=0){
			knm1c=Cmom<T_high>(lambdat<T_high>(mc.Lt(ind[n-2])),
					lambda<T_high>(-conj(Z4*lambdat<T_high>(mc.Lt(ind[n-2]))),conj(Z3*lambdat<T_high>(mc.Lt(ind[n-2])))));
		}
		else{
			knm1c=Cmom<T_high>(lambdat<T_high>(mc.Lt(ind[n-2])),
					lambda<T_high>(conj(Z4*lambdat<T_high>(mc.Lt(ind[n-2]))),-conj(Z3*lambdat<T_high>(mc.Lt(ind[n-2])))));
		}
	}
#endif
	// notice that the sqrt is always taken in the vicinity of 1, away of its branch cut
	complex<T_high> sqrtt=sqrt(-sum*sum/(complex<T_high>(2,0)*sum*knm1c.P()));
	// final extended k_{n-1}
#if _OLD_PHASE_CONVENTION
	Cmom<T_high> knm1f(sqrtt*lambdat<T_high>(mc.Lt(ind[n-2])),sqrtt*lambda<T_high>(mc.L(ind[n-2])));
#else
	Cmom<T_high> knm1f;
	if(abs(imag(mc.p(ind[n-2]).E()))>SMALL_COMPLEX_PART ||
			abs(imag(mc.p(ind[n-2]).X()))>SMALL_COMPLEX_PART ||
			abs(imag(mc.p(ind[n-2]).Y()))>SMALL_COMPLEX_PART ||
			abs(imag(mc.p(ind[n-2]).Z()))>SMALL_COMPLEX_PART){
		_WARNING("CAREFUL: Using extend_real for a momentum with large complex part");
		knm1f=Cmom<T_high>(sqrtt*lambdat<T_high>(mc.Lt(ind[n-2])),sqrtt*lambda<T_high>(mc.L(ind[n-2])));
	}
	else{
		if(real(mc.p(ind[n-2]).E())<=0){
			knm1f=Cmom<T_high>(sqrtt*lambdat<T_high>(mc.Lt(ind[n-2])),
					sqrtt*lambda<T_high>(-conj(Z4*lambdat<T_high>(mc.Lt(ind[n-2]))),conj(Z3*lambdat<T_high>(mc.Lt(ind[n-2])))));
		}
		else{
			knm1f=Cmom<T_high>(sqrtt*lambdat<T_high>(mc.Lt(ind[n-2])),
					sqrtt*lambda<T_high>(conj(Z4*lambdat<T_high>(mc.Lt(ind[n-2]))),-conj(Z3*lambdat<T_high>(mc.Lt(ind[n-2])))));
		}
	}
#endif
	sum+=knm1f.P();
	// extended k_{n}, just by momentum conservation
	Cmom<T_high> knf=(-sum);
	// to make sure we didn't hit a branch cut
	// completed k_{n}
	Cmom<T_high> knc(lambdat<T_high>(mc.Lt(ind[n-1])),lambda<T_high>(mc.L(ind[n-1])));

	complex<T_high> phase(1,0);
	complex<T_high> phaseinv(1,0);

	if( (abs((knf.L()-knc.L())*Z1)>T_high(1)/T_high(100000) ||
		abs((knf.L()-knc.L())*Z2)>T_high(1)/T_high(100000) ||
		abs(Z3*(knf.Lt()-knc.Lt()))>T_high(1)/T_high(100000) ||
		abs(Z4*(knf.Lt()-knc.Lt()))>T_high(1)/T_high(100000) )
	){
		BH_DEBUG_MESSAGE("\nfixing phases in extended momentum --- extend of complex momenta\n");
		complex<T_high> knca=(knc.L()*Z1+knc.L()*Z2)/T_high(2);
		complex<T_high> knfa=(knf.L()*Z1+knf.L()*Z2)/T_high(2);
		phase=knca/knfa;
		phaseinv=knfa/knca;

		BH_DEBUG_PRINT(knc.L());
		BH_DEBUG_PRINT(knc.Lt());
		BH_DEBUG_PRINT(knf.L());
		BH_DEBUG_PRINT(knf.Lt());
		BH_DEBUG_PRINT(phase*knf.L());
		BH_DEBUG_PRINT(phaseinv*knf.Lt());
		BH_DEBUG_PRINT(mc.p(ind[n-1]));
		BH_DEBUG_MESSAGE("\n");
	}

	momenta.push_back(knm1f);
	momenta.push_back(Cmom<T_high>(phaseinv*knf.Lt(),phase*knf.L()));


	return momentum_configuration<T_high>(momenta);
}



// Explicit instantiation
template  momentum_configuration<RHP>
  extend(momentum_configuration<R>& k,
         const vector<int>& indices /* of the momenta within "k" */,
         const RHP& dummy /* To signal target type */);

template  momentum_configuration<RVHP>
  extend(momentum_configuration<RHP>& k,
         const vector<int>& indices /* of the momenta within "k" */,
         const RVHP& dummy /* To signal target type */);

template momentum_configuration<RHP> extend_daniel<R,RHP>(const momentum_configuration<R>& mc,const vector<int>& ind);
template momentum_configuration<RVHP> extend_daniel<R,RVHP>(const momentum_configuration<R>& mc,const vector<int>& ind);
template momentum_configuration<RVHP> extend_daniel<RHP,RVHP>(const momentum_configuration<RHP>& mc,const vector<int>& ind);
template momentum_configuration<RHP> extend_real<R,RHP>(const momentum_configuration<R>& mc,const vector<int>& ind);
template momentum_configuration<RVHP> extend_real<R,RVHP>(const momentum_configuration<R>& mc,const vector<int>& ind);
template momentum_configuration<RVHP> extend_real<RHP,RVHP>(const momentum_configuration<RHP>& mc,const vector<int>& ind);

#if BH_USE_GMP
template momentum_configuration<RGMP> extend_real<RGMP,RGMP>(const momentum_configuration<RGMP>& mc,const vector<int>& ind);
#endif


}

