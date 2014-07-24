/* Integrals_ep.cpp */

/*  David A. Kosower,  November 15, 2007  */

/*  Integral functions for one-loop amplitudes

    Based on older code by Lance Dixon; on hep-ph/0502165v2 [Hameren,
    Vollinga and Weinzierl]; on hep-ph/9403226v2 (except for the
    four-mass box, whose basic form is taken from hep-ph/9306240v2
    equiv to hep-ph/9409265v1, fixing sign errors in 9403226v2); and
    on my Mathematica code.  Uses the BlackHat library.

    All externally-visible functions take the following argument:
      a momentum configuration (MomentumConfiguration = mom_conf);
      a scale mu^2;
      either: (a) the appropriate number of vectors of indices [into the
                  momentum configuration], each specifying the
                  set of momenta entering at a corner;
              (b) a partition specifying these vectors.
      a specification of the desired orders in epsilon.
    They return a Series object, giving the complex values of
    the coefficients of different orders in epsilon.

    Vectors specifying massless arguments must contain a single entry, but
    the masslessness of the corresponding entry in the momentum configuration
    will not be checked.  All vectors should be real, but this will not be
    checked either.  The invariants are all taken (implicitly) to have
    an infinitesimal positive imaginary part.

    The integrals have a factor of c_\Gamma removed, so that those
    with 1/e^2 divergences start with precisely -(-1)^n n/e^2/denominator
    as their leading term, where n is the number of massless-massless
    antennae.  As e->0, c_\Gamma -> 1/(4 \pi)^2, which is the only
    additional factor required in the evaluation of a virtual amplitude.
    (Note that the factor of c_\Gamma is included, in expanded form, in
     the formulae in hep-ph/0502165v2; one must add e^2 Zeta[2]/2 * leading
     singular term to remove it).  The integrals are (non-calligraphic) "I"
     functions of our papers, having factors of (-1)^n adjusted from
     the full momentum-space integrals (this changed on 2/6/08 from
     earlier versions).

   1/22/08: Updated with Daniel Maitre's changes for integration into BlackHat
   2/6/08: Updated with more of Daniel Maitre's changes to handle extended-
           precision types: C -> complex<T> and CSeries -> Series<complex<T> >,
           MomentumConfiguration -> momentum_configuration<T>;
           also, use the caches in momentum_configuration now.

   Modifications by Darren
   21/5/09: Altered for use with the eval_param code, only small changes were made and the
		    index based routines removed as the eval param does not use indices.

*/

//#include "header.h"
#include<map>
#include "integrals_ep.h"
#include "integrals.h"   /// for epsilon<T>
#include "Series.h"
#include "polylog.h"
#include "partitions.h"
#include "genkey.h"
#include "BH_debug.h"

#define isnt !=
#define is ==

#define IntCaching 0
#define SaveSymmetry 0
#define BasicFunctionCaching 0
#define GenKey GenKey1
#define Direct 1 // Direct-invariant routines (see Int() for explanation)


#define Re(v) real(v)
#define Im(v) imag(v)
static int version = 2;
using namespace std;

namespace BH {
//const double tolerance = std::numeric_limits<const double>::epsilon()*2;
// The above seems to give 0 (!)
const double tolerance = 1.e-10;

template<class T> inline T epsilon();
template <> inline R epsilon<R>(){return(1.0e-13);};
template <> inline RHP epsilon<RHP>(){return(RHP("1.0e-52"));};
template <> inline RVHP epsilon<RVHP>(){return(RVHP("1.0e-82"));}

#if BH_USE_GMP
template <> inline RGMP epsilon<RGMP>(){ return exp10(-RGMP(RGMP::get_current_nbr_digits()-3));}
#endif


#define ln log

// \theta(r); it is a macro rather than a function because we want to avoid
// evaluation of v unless r >= 0 (the evaluation might be invalid)
#define Theta(r,v) ( (r) >= 0 ? (v) : 0 )

// \theta(r1)-\theta(r2)
#define Theta2(r1,r2,v) \
 ( ((r1) >= 0 and (r2) < 0) ? (v) \
   : ((r1) < 0 and (r2) >= 0) ? (-v) : 0 )

// 3-mass triangle Gram det
template <class T> T Delta(T s1, T s2,T s3)
{return(s1*s1+s2*s2+s3*s3-T(2)*s1*s2-T(2)*s2*s3-T(2)*s3*s1);}
	
	
// Direct (uncacheable) version of above routine, function of invariants
template<class T> complex<T> CLn(const T& s1, const T& s2)
{complex<T> result;
   // Eq. (115a) of hep-ph/0502165v2
   result = complex<T>(ln(fabs(s1/s2)),Theta2(s1,s2,-pi<T>()));

   return result;
}

// Direct (uncacheable) version of above routine, function of invariants
template <class T> complex<T> CLnM(const T& s1, const T& musq /* should be > 0 */)
{
 complex<T> result;

 // Eq. (115a) of hep-ph/0502165v2
 result = complex<T>(ln(abs(s1/musq)),Theta(s1,-pi<T>()));

 return result;
}

// Direct (uncacheable) version of above routine, function of invariants
template <class T> complex<T> CLi2b(const T& s1, const T& s2)
  // Li_2(1-(-s_{i1})/(-s_{i2})), with appropriate branch cut
{
 complex<T> result;

   // Eq. (115b) of hep-ph/0502165v2
   result = complex<T>(ReLi2(T(1)-s1/s2),
                       Theta(-s1/s2,Theta2(s1,s2,pi<T>()*ln(T(1)-s1/s2))));

 return result;
}

template <class T> complex<T> CLi2b(const T& s1, const T& s2,
                                     const T& s3, const T& s4)
  // Li_2(1-(-s_{i1})(-s_{i2})/((-s_{i3})(-s_{i4}))), with appropriate branch cut
{
 complex<T> result;

   // Eq. (115b) of hep-ph/0502165v2
   // Eq. following eq. (115) in hep-ph/0502165v2

   // -Im ln( (-s_1)(-s_2)/((-s_3)(-s_4)) )
   T imln = Theta2(s1,s3,pi<T>())+Theta2(s2,s4,pi<T>());
   complex<T> logb(log(abs(T(1)-s1*s2/(s3*s4))),
                   Theta(s1*s2/(s3*s4)-T(1),-imln*T(0.5)));
   result = complex<T>(ReLi2(T(1)-s1*s2/(s3*s4))-Im(logb)*imln,Re(logb)*imln);
 return result;
}


template<class T> T square(const T& x);//Defined in integrals.cpp

// Square sum of momenta (explicit local computation); should only be
// used for sums of external momenta, which will always be real
template <class T> T SqSum(const eval_param<T> & k, const vector<int>& v1)
{momentum<complex<T> > sum;

  for (int j = 0;  j < v1.size();  j += 1)
    sum += k.mom(v1[j]);

  return Re(sum*sum);
}

template <class T> T SqSum(const eval_param<T> & k, const vector<int>& v1,
      const vector<int>& v2)
{momentum<complex<T> > sum;

  for (int j = 0;  j < v1.size();  j += 1)
    sum += k.mom(v1[j]);

  for (int j = 0;  j < v2.size();  j += 1)
    sum += k.mom(v2[j]);

  return Re(sum*sum);
}

// Box functions

template <class T>  T Sign(T x); //Defined in integrals.cpp

template <class T> complex<T> CLi2(T arg,T branchSide)
  /* Complex dilog; "branchSide" is coefficient of I e in "arg"
     which determines the sign of the imaginary part */
{return(complex<T>(ReLi2(arg),Theta(arg-T(1),pi<T>()*ln(arg)*Sign(branchSide))));}


template <class T> inline T eta(T a, T b, T ab)
// args are Im a, Im b, and Im (a b) compared to eta in Denner, Nierste, & Scharf
// return value doesn't include 2 Pi I factor
{return( a < 0 and b < 0 and ab >= 0 ? 1
         : (a >= 0 and b >= 0 and ab < 0) ? -1 : 0 );}


// 11/20/08
// Direct versions of box functions, arguments are invariants rather
// than indices; will eventually replace above functions
template <class T> complex<T> I4w0m(int order,
                                    const T& musq, const T& s, const T& t)
{
 T result;
	BH_DEBUG_MESSAGE4("I4w0m called for s=",s," t= ",t);

 switch (order) {
 case -2:
   result = (T(4)/(s*t));
   return(result);
 case -1:
   return(T(-2)/(s*t) * (CLnM(s,musq)+CLnM(t,musq)));
 case 0:
   static const T PiSquared = square(pi<T>());
   return( (T(2)*CLnM(s,musq)*CLnM(t,musq)-PiSquared)/(s*t) );
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
 }}

template <class T> complex<T> I4w1m(int order,
                                    const T& musq, const T& s, const T& t,
                                    const T& msq)
{complex<T> result;
BH_DEBUG_MESSAGE6("I4w1m called for s=",s," t= ",t," msq= ",msq);

 switch (order) {
 case -2:
   return(T(2)/(s*t));
 case -1:
   return( T(-2)/(s*t)*(CLnM(s,musq)+CLnM(t,musq)-CLnM(msq,musq)) );
 case 0:
   static const T PiSquaredOverThree = square(pi<T>())/T(3);
   result = (T(1)/(s*t) *(T(2)*CLnM(s,musq)*CLnM(t,musq)-square(CLnM(msq,musq))
            -T(2)*(CLi2b(msq,s)+CLi2b(msq,t)) - PiSquaredOverThree) );
   return result;
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
 }
}

template <class T> complex<T> I4w2me(int order,
				     const T& musq, const T& s, const T& t,
                                     const T& m1sq, const T& m3sq)
{
	BH_DEBUG_MESSAGE8("I4w1m called for s=",s," t= ",t," m1sq= ",m1sq," m3sq= ",m3sq);

	switch (order) {
 case -2:
   return complex<T>(0,0) ;
 case -1:
   return(T(-2)/(s*t-m1sq*m3sq)*(CLnM(s,musq)+CLnM(t,musq)-CLnM(m1sq,musq)-CLnM(m3sq,musq)));
 case 0:
   return (T(1)/(s*t-m1sq*m3sq) *
            (T(2)*CLnM(s,musq)*CLnM(t,musq)
             -square(CLnM(m1sq,musq))-square(CLnM(m3sq,musq))
             -T(2)*(CLi2b(m1sq,s)+CLi2b(m1sq,t)+CLi2b(m3sq,s)+CLi2b(m3sq,t))
             +T(2)*CLi2b(m1sq,m3sq,s,t)) );
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
 }
}

template <class T> complex<T> I4w2mh(int order,
				     const T& musq, const T& s, const T& t,
                                     const T& m1sq, const T& m2sq)
{
	BH_DEBUG_MESSAGE8("I4w1m called for s=",s," t= ",t," m1sq= ",m1sq," m2sq= ",m2sq);

	switch (order) {
 case -2:
   return (T(1)/(s*t));
 case -1:
   return( T(-1)/(s*t) * (CLnM(s,musq)+T(2)*CLnM(t,musq)
                          -CLnM(m1sq,musq)-CLnM(m2sq,musq)) );
 case 0:
   {complex<T> lnS = CLnM(s,musq), ln1 = CLnM(m1sq,musq),
       ln2 = CLnM(m2sq,musq);

   return( T(1)/(s*t) *
            (T(2)*lnS*CLnM(t,musq)
             +square(lnS-ln1-ln2)/T(2)-square(ln1)-square(ln2)
             -T(2)*(CLi2b(m1sq,t)+CLi2b(m2sq,t))) );}
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
 }
}

template <class T> complex<T> I4w3m(int order,
				    const T& musq, const T& s, const T& t,
				    const T& m1sq, const T& m2sq,
				    const T& m3sq)
{switch (order) {
 case -2:
   return complex<T>(0,0);
 case -1:
   return(-T(1)/(s*t-m1sq*m3sq)*(CLnM(s,musq)+CLnM(t,musq)
                                 -CLnM(m1sq,musq)-CLnM(m3sq,musq)));
 case 0:
   {complex<T> lnS = CLnM(s,musq), lnT = CLnM(t,musq),
       ln1 = CLnM(m1sq,musq), ln2 = CLnM(m2sq,musq), ln3 = CLnM(m3sq,musq);

   return(T(1)/(s*t-m1sq*m3sq) *
            (T(2)*lnS*lnT
             +square(lnS-ln1-ln2)/T(2)+square(lnT-ln2-ln3)/T(2)
              -square(ln1)-square(ln2)-square(ln3)
             -T(2)*(CLi2b(m3sq,s)+CLi2b(m1sq,t))
             +T(2)*CLi2b(m1sq,m3sq,s,t)) );}
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
 }
}

// Now correct 10/3/08 (copied from above)
template <class T> complex<T> I4w4m(int order,
				    const T& musq, const T& s, const T& t,
				    const T& m1sq, const T& m2sq,
				    const T& m3sq, const T& m4sq)
{switch (order) {
 case -2:
 case -1:
   return complex<T>(0,0);
 case 0:
   {T lambda1 = m1sq*m3sq/(s*t), lambda2 = m2sq*m4sq/(s*t);
   // rho^2 can be negative, so need to keep everything complex
     T rhoSq = square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2;
   complex<T> result;
   string tag, masstag;

#if 0
   masstag += m1sq > 0 ? "+" : "-";
   masstag += m2sq > 0 ? "+" : "-";
   masstag += m3sq > 0 ? "+" : "-";
   masstag += m4sq > 0 ? "+" : "-";
   masstag += s > 0 ? "+" : "-";
   masstag += t > 0 ? "+" : "-";
   masstag += " ";
#endif

   // But we must treat real & complex case for rho separately
   // because of branch cut handling
   if (rhoSq >= 0)
      {//T rho = sqrt(rhoSq);
        complex<T> rho;

        tag = "rho^2 > 0 "+masstag;
#if DumpArguments
   cout << "rho^2 > 0" << endl;
   cout << "m^2_i,s,t: ";
   cout << m1sq << " ";
   cout << m2sq << " ";
   cout << m3sq << " ";
   cout << m4sq << " ";
   cout << s << " ";
   cout << t;
   cout << endl;
#endif

   complex<T> LiTerms, etaTerms, lnTerms;
   complex<T> ln1, ln2;
   complex<T> Ieps = complex<T>(0,epsilon<T>());
   complex<T> lambda1, lambda2, arg1a, arg2a, arg1b, arg2b;
   complex<T> critA, critB;
   complex<T> x1, x2;
   T sign = 1;
   complex<T> norm;

#define Sqrt sqrt
#define Power pow
   if (version is 1)
      {// 2007 -- with I e in m^2 & s, t, x_i updated
        lambda1 = (m1sq+Ieps)*(m3sq+Ieps)/((s+Ieps)*(t+Ieps));
        lambda2 = (m2sq+Ieps)*(m4sq+Ieps)/((s+Ieps)*(t+Ieps));
        rho = complex<T>(Re(sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2)),0);

        //     complex<T> x1 = (s+Ieps)*(arg1a-1.)/(m3sq+Ieps)+Ieps*m2sq,
        x1 = (s)*(-T(1)-lambda1+lambda2+rho)/(T(2)*m3sq)
          +Ieps*m2sq*s*t/(T(4)*rhoSq);
        x2 = (s)*(-T(1)-lambda1+lambda2-rho)/(T(2)*m3sq)
          -Ieps*m2sq*s*t/(T(4)*rhoSq);
        norm = s*t*rho;}
   else if (version is 2)
      {// 2008 version, only I eps as in Denner et al, none in m_i^2
        lambda1 = complex<T>(m1sq*m3sq/(s*t),0);
        lambda2 = complex<T>(m2sq*m4sq/(s*t),0);
        rho = sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2);
        x1 = -s/(T(2)*m3sq) - m1sq/(T(2)*t) + (m2sq*m4sq)/(T(2)*m3sq*t) +
          Sqrt(square(m1sq*m3sq) - T(2)*m1sq*m2sq*m3sq*m4sq +
               square(m2sq*m4sq) + T(4)*Ieps*m2sq*m3sq*t -
               T(2)*m1sq*m3sq*s*t - T(2)*m2sq*m4sq*s*t
               + square(s*t))/(T(2)*m3sq*t);
        x2 =  -s/(T(2)*m3sq) - m1sq/(T(2)*t) + (m2sq*m4sq)/(T(2)*m3sq*t) -
          Sqrt(square(m1sq*m3sq) - T(2)*m1sq*m2sq*m3sq*m4sq +
               square(m2sq*m4sq) + T(4)*Ieps*m2sq*m3sq*t -
               T(2)*m1sq*m3sq*s*t - T(2)*m2sq*m4sq*s*t
               + square(s*t))/(T(2)*m3sq*t);
        norm = t*m3sq*(x1-x2);
        // Additional sign from s t?  Why?
        sign = Sign(s*t);
#if 0
        x1 = (s)*(-T(1)-lambda1+lambda2+rho)/(T(2)*(m3sq+Ieps))
          +Ieps*m2sq*s*t/(T(4)*rhoSq);
        x2 = (s)*(-T(1)-lambda1+lambda2-rho)/(T(2)*(m3sq+Ieps))
          -Ieps*m2sq*s*t/(T(4)*rhoSq);
#endif
      }
   else if (version is 3)
      {// Use Davydychev--Ussyukina form (PLB 298:262 (1993)), without Ieps
        T s1 = m1sq*m3sq, s2 = m2sq*m4sq, s3 = s*t;
        // Why is there a sign?
        T prefactor = T(-1);
        T Delta3 = -Delta(s1,s2,s3);
        complex<T> result(0,0);
        T delta[3] = {s2+s3-s1,s1+s3-s2,s1+s2-s3};
        complex<T> sqrtDelta;
        int count = 0;
#if DumpTerms
        cout << "Delta3: " << Delta3 << endl;
        cout << "signs: " << Sign(s1) << " " << Sign(s2) << " " << Sign(s3) << endl;
#endif
        if (Delta3 >= 0)
           {T sqrtDelta = sqrt(Delta3);
             /* In this case, the arguments of the dilogs are generic complex
                numbers, and we don't need to do anything special */
             for (int j = 0;  j < 3;  j += 1)
                {complex<T> arg = (sqrtDelta+complex<T>(0,delta[j]))/(-sqrtDelta+complex<T>(0,delta[j]));
                  result += Li2(arg)-Li2(T(1)/arg);}
             result *= complex<T>(0,-1)/sqrtDelta;}
        else { // This branch seems to be problematic: rho^2 > 0 -> Delta3 < 0
          /* Here, the arguments of the dilogs become real, and we have to
             be careful about staying on the right side of the branch cut.
             On the other hand, which sign we choose for sqrt(Delta3) is
             irrelevant. */
          T ImSqrtDelta = sqrt(-Delta3);
#define Crit(s1,s2,s3) (s1*s2 - s2*s2 + s1*s3 - s3*s3)
          T crit[3] = {Crit(s1,s2,s3),Crit(s2,s3,s1),Crit(s3,s1,s2)};
#undef Aux
#define Aux(s1,s2,s3) -(s1+s2+s3)
//          T aux[3] = {Aux(s1,s2,s3),Aux(s2,s3,s1),Aux(s3,s1,s2)};
#if 0
          tag += " ";
          tag += s1 > 0 ? "+" : "-";
          tag += s2 > 0 ? "+" : "-";
          tag += s3 > 0 ? "+" : "-";
          tag += " ";
          tag += aux[0] > 0 ? "+" : "-";
          tag += aux[1] > 0 ? "+" : "-";
          tag += aux[2] > 0 ? "+" : "-";
          tag += " ";
#endif
        T crit1, // coefficient of I eps in x_i
          crit2; // coefficient of I eps' from m^2->m^2+I eps', in x_i;

        crit1 = 1/T(2)* m2sq/(m3sq* t* sqrt(-Delta3));

        crit2 = 1/T(2)*(m1sq*square(m3sq) - m2sq*m3sq*m4sq + m2sq*m3sq*t - square(m3sq)*t - m2sq*m4sq*t +
                        m3sq*m4sq*t - m3sq*square(t) + s*square(t))/(square(m3sq)*square(t))
          +1/T(2)*(square(m1sq*m3sq)*m3sq - T(2)*m1sq*m2sq*square(m3sq)*m4sq + square(m2sq*m4sq)*m3sq +
                   m1sq*m2sq*square(m3sq)*t - m1sq*square(m3sq)*m3sq*t - m1sq*m2sq*m3sq*m4sq*t -
                   square(m2sq)*m3sq*m4sq*t + m1sq*square(m3sq)*m4sq*t + m2sq*square(m3sq)*m4sq*t +
                   square(m2sq*m4sq)*t - m2sq*m3sq*square(m4sq)*t - m1sq*square(m3sq)*s*t -
                   m2sq*m3sq*m4sq*s*t + m1sq*square(m3sq)*square(t) + m2sq*m3sq*m4sq*square(t) -
                   m1sq*m3sq*s*square(t) + m2sq*m3sq*s*square(t) + square(m3sq)*s*square(t) - T(2)*m2sq*m4sq*s*square(t) +
                   m3sq*m4sq*s*square(t) - m3sq*s*t*square(t) + square(s*t)*t)/(square(m3sq)*square(t));

#if DumpTerms
        cout << "crit1: " << crit1 << endl;
        cout << "crit2: " << crit2 << endl;
#endif

        tag += " C";
        tag += Sign(crit1) > 0 ? "+" : "-";
        tag += Sign(crit2) > 0 ? "+" : "-";
        tag += Sign(crit1+crit2)/Sign(crit1) > 0 ? "v" : "x";
        tag += " ";

        if (m1sq > 0) count += 1;
        if (m2sq > 0) count += 1;
        if (m3sq > 0) count += 1;
        if (m4sq > 0) count += 1;

        tag += "S ";
          {
        x1 = -s/(T(2)*m3sq) - m1sq/(T(2)*t) + (m2sq*m4sq)/(T(2)*m3sq*t) +
          Sqrt(square(m1sq*m3sq) - T(2)*m1sq*m2sq*m3sq*m4sq +
               square(m2sq*m4sq) + T(4)*Ieps*m2sq*m3sq*t -
               T(2)*m1sq*m3sq*s*t - T(2)*m2sq*m4sq*s*t
               + square(s*t))/(T(2)*m3sq*t);
        x2 = -s/(T(2)*m3sq) - m1sq/(T(2)*t) + (m2sq*m4sq)/(T(2)*m3sq*t) -
          Sqrt(square(m1sq*m3sq) - T(2)*m1sq*m2sq*m3sq*m4sq +
               square(m2sq*m4sq) + T(4)*Ieps*m2sq*m3sq*t -
               T(2)*m1sq*m3sq*s*t - T(2)*m2sq*m4sq*s*t
               + square(s*t))/(T(2)*m3sq*t);
        norm = t*m3sq*(x1-x2);
#if 0
        tag += Sign(Im(x1)) > 0 ? "+" : "-";
        tag += " ";
#endif
          }
          T sign = count is 1 ? -1 : 1;
          for (int j = 0;  j < 3;  j += 1)
             {T arg = (ImSqrtDelta+delta[j])/(-ImSqrtDelta+delta[j]);
               // It seems that the D-U formula sometimes just needs the opposite imaginary part!
               complex<T> aux;
               result += CLi2(arg,sign*crit[j])-CLi2(T(1)/arg,-sign*crit[j]);}
          result /= -ImSqrtDelta;}
#if DumpTag
        cout << tag << count << "  ";
#endif
        return (prefactor*result);}
   else if (version is 4)
      {// 2008 version, I eps as in Denner et al, but also as in m_i^2
        complex<T> M1sq(m1sq,epsilon<T>()), M2sq(m2sq,epsilon<T>()),
          M3sq(m3sq,epsilon<T>()), M4sq(m4sq,epsilon<T>());
        complex<T> S(s,epsilon<T>()), Ti(t,epsilon<T>());

        lambda1 = M1sq*M3sq/(S*Ti);
        lambda2 = M2sq*M4sq/(S*Ti);
        rho = sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2);
        x1 = -S/(T(2)*M3sq) - M1sq/(T(2)*Ti) + (M2sq*M4sq)/(T(2)*M3sq*Ti) +
          Sqrt(square(M1sq*M3sq) - T(2)*M1sq*M2sq*M3sq*M4sq +
               square(M2sq*M4sq) + T(4)*Ieps*M2sq*M3sq*t -
               T(2)*M1sq*M3sq*S*Ti - T(2)*M2sq*M4sq*S*Ti
               + square(S*Ti))/(T(2)*M3sq*Ti);
        x2 = -S/(T(2)*M3sq) - M1sq/(T(2)*Ti) + (M2sq*M4sq)/(T(2)*M3sq*Ti) -
          Sqrt(square(M1sq*M3sq) - T(2)*M1sq*M2sq*M3sq*M4sq +
               square(M2sq*M4sq) + T(4)*Ieps*M2sq*M3sq*t -
               T(2)*M1sq*M3sq*S*Ti - T(2)*M2sq*M4sq*S*Ti
               + square(S*Ti))/(T(2)*M3sq*Ti);
        norm = t*m3sq*(x1-x2);
        // Additional sign from s t?  Why?
        sign = Sign(s*t);
      }
   else if (version is 5)
      {// Use Davydychev--Ussyukina form (PLB 298:262 (1993)), with Ieps in m^2
        complex<T> s1 = complex<T>(m1sq,epsilon<T>())*complex<T>(m3sq,epsilon<T>()),
          s2 = complex<T>(m2sq,epsilon<T>())*complex<T>(m4sq,epsilon<T>()),
          s3 = complex<T>(s,epsilon<T>())*complex<T>(t,epsilon<T>());
        // Why is there a sign?
        T prefactor = T(-1);
        int count = 0;
        complex<T> Delta3 = -Delta(s1,s2,s3);
        complex<T> result(0,0);
        complex<T> delta[3] = {s2+s3-s1,s1+s3-s2,s1+s2-s3};
#if DumpTerms
        cout << "Delta3: " << Delta3 << endl;
        cout << "signs: " << Sign(Re(s1)) << " " << Sign(Re(s2)) << " " << Sign(Re(s3)) << endl;
#endif
        tag += " ";
        tag += Re(s1) > 0 ? "+" : "-";
        tag += Re(s2) > 0 ? "+" : "-";
        tag += Re(s3) > 0 ? "+" : "-";
        tag += " ";

        T crit1, // coefficient of I eps in x_i
          crit2; // coefficient of I eps' from m^2->m^2+I eps', in x_i;

        crit1 = 1/T(2)* m2sq/(m3sq* t* Re(sqrt(-Delta3)));

        crit2 = 1/T(2)*(m1sq*square(m3sq) - m2sq*m3sq*m4sq + m2sq*m3sq*t - square(m3sq)*t - m2sq*m4sq*t +
                        m3sq*m4sq*t - m3sq*square(t) + s*square(t))/(square(m3sq)*square(t))
          +1/T(2)*(square(m1sq*m3sq)*m3sq - T(2)*m1sq*m2sq*square(m3sq)*m4sq + square(m2sq*m4sq)*m3sq +
                   m1sq*m2sq*square(m3sq)*t - m1sq*square(m3sq)*m3sq*t - m1sq*m2sq*m3sq*m4sq*t -
                   square(m2sq)*m3sq*m4sq*t + m1sq*square(m3sq)*m4sq*t + m2sq*square(m3sq)*m4sq*t +
                   square(m2sq*m4sq)*t - m2sq*m3sq*square(m4sq)*t - m1sq*square(m3sq)*s*t -
                   m2sq*m3sq*m4sq*s*t + m1sq*square(m3sq)*square(t) + m2sq*m3sq*m4sq*square(t) -
                   m1sq*m3sq*s*square(t) + m2sq*m3sq*s*square(t) + square(m3sq)*s*square(t) - T(2)*m2sq*m4sq*s*square(t) +
                   m3sq*m4sq*s*square(t) - m3sq*s*t*square(t) + square(s*t)*t)/(square(m3sq)*square(t));

#if DumpTerms
        cout << "crit1: " << crit1 << endl;
        cout << "crit2: " << crit2 << endl;
#endif

        tag += " C";
        tag += Sign(crit1) > 0 ? "+" : "-";
        tag += Sign(crit2) > 0 ? "+" : "-";
        tag += Sign(crit1+crit2)/Sign(crit1) > 0 ? "v" : "x";
        tag += " ";

        if (m1sq > 0) count += 1;
        if (m2sq > 0) count += 1;
        if (m3sq > 0) count += 1;
        if (m4sq > 0) count += 1;

        tag += "S ";

        complex<T> ImSqrtDelta = sqrt(-Delta3);
        /* In this case, the arguments of the dilogs are generic complex
           numbers, and we don't need to do anything special */
        for (int j = 0;  j < 3;  j += 1)
           {complex<T> arg = (ImSqrtDelta+delta[j])/(-ImSqrtDelta+delta[j]);
            // One term has the wrong sign of imaginary part; but why?
            tag += Im( ((Li2(arg)-Li2(T(1)/arg))/ImSqrtDelta) ) > 0 ? "+" : "-";
#if DumpTerms
            cout << j << ": " << arg << "; "
                  << ((Li2(arg))/ImSqrtDelta)
                  << (-Li2(T(1)/arg)/ImSqrtDelta)
                  << endl;
#endif
            result += Li2(arg)-Li2(T(1)/arg);}
        tag += " ";
        result *= T(-1)/ImSqrtDelta;
#if DumpTag
        cout << tag << count << "  ";;
#endif
        return (prefactor*result);}
   else if (version is 6)
      {// 2008 version, formulae from Denner et al, but I eps only in m_i^2
        complex<T> M1sq(m1sq,epsilon<T>()), M2sq(m2sq,epsilon<T>()),
          M3sq(m3sq,epsilon<T>()), M4sq(m4sq,epsilon<T>());
        complex<T> S(s,epsilon<T>()), Ti(t,epsilon<T>());

        lambda1 = M1sq*M3sq/(S*Ti);
        lambda2 = M2sq*M4sq/(S*Ti);
        rho = sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2);
        x1 = -S/(T(2)*M3sq) - M1sq/(T(2)*Ti) + (M2sq*M4sq)/(T(2)*M3sq*Ti) +
          Sqrt(square(M1sq*M3sq) - T(2)*M1sq*M2sq*M3sq*M4sq +
               square(M2sq*M4sq) -
               T(2)*M1sq*M3sq*S*Ti - T(2)*M2sq*M4sq*S*Ti
               + square(S*Ti))/(T(2)*M3sq*Ti);
        x2 = -S/(T(2)*M3sq) - M1sq/(T(2)*Ti) + (M2sq*M4sq)/(T(2)*M3sq*Ti) -
          Sqrt(square(M1sq*M3sq) - T(2)*M1sq*M2sq*M3sq*M4sq +
               square(M2sq*M4sq) -
               T(2)*M1sq*M3sq*S*Ti - T(2)*M2sq*M4sq*S*Ti
               + square(S*Ti))/(T(2)*M3sq*Ti);
        norm = t*m3sq*(x1-x2);
        // Additional sign from s t?  Why?
        sign = Sign(s*t);
      }
#if DumpArguments
   cout << "-x1: " << (-x1) << ", " << "-x2: " << (-x2) << endl;
#endif
   // 2nd arguments to eta()
   critA = (m3sq+Ieps)/(s+Ieps);
   critB = (t+Ieps)/(m1sq+Ieps);
   // Arguments of 1st polylog (and corresponding ln in eta*ln) for k = 1,2 in eq. (41)
   arg1a = T(1)+critA*x1;
   arg2a = T(1)+critA*x2;
   // Arguments of 2nd polylog (and corresponding ln in eta*ln) for k = 1,2 in eq. (41)
   arg1b = T(1)+critB*x1;
   arg2b = T(1)+critB*x2;

   T crit1, // coefficient of I eps in x_i
     crit2; // coefficient of I eps' from m^2->m^2+I eps', in x_i;

     crit1 = 1/T(2)* m2sq/(m3sq* t* Re(rho));

   crit2 = 1/T(2)*(m1sq*square(m3sq) - m2sq*m3sq*m4sq + m2sq*m3sq*t - square(m3sq)*t - m2sq*m4sq*t +
         m3sq*m4sq*t - m3sq*square(t) + s*square(t))/(square(m3sq)*square(t))
     +1/T(2)*(square(m1sq*m3sq)*m3sq - T(2)*m1sq*m2sq*square(m3sq)*m4sq + square(m2sq*m4sq)*m3sq +
              m1sq*m2sq*square(m3sq)*t - m1sq*square(m3sq)*m3sq*t - m1sq*m2sq*m3sq*m4sq*t -
              square(m2sq)*m3sq*m4sq*t + m1sq*square(m3sq)*m4sq*t + m2sq*square(m3sq)*m4sq*t +
              square(m2sq*m4sq)*t - m2sq*m3sq*square(m4sq)*t - m1sq*square(m3sq)*s*t -
   m2sq*m3sq*m4sq*s*t + m1sq*square(m3sq)*square(t) + m2sq*m3sq*m4sq*square(t) -
   m1sq*m3sq*s*square(t) + m2sq*m3sq*s*square(t) + square(m3sq)*s*square(t) - T(2)*m2sq*m4sq*s*square(t) +
              m3sq*m4sq*s*square(t) - m3sq*s*t*square(t) + square(s*t)*t)/(square(m3sq)*square(t));

#if DumpTerms
   cout << "crit1: " << crit1 << endl;
   cout << "crit2: " << crit2 << endl;
#endif

   tag += " C";
   tag += Sign(crit1) > 0 ? "+" : "-";
   tag += Sign(crit2) > 0 ? "+" : "-";
   tag += Sign(crit1+crit2)/Sign(crit1) > 0 ? "v" : "x";
   tag += " ";

#if DumpArguments
     cout << "arg1a: " << arg1a << "; ";
     cout << "arg2a: " << arg2a << "; ";
     cout << "arg1b: " << arg1b << "; ";
     cout << "arg2b: " << arg2b << "; ";
     cout << endl;

     cout << "critA: " << critA << ", "
          << "critB: " << critB << " " << endl;
     cout << "eta(-x_1,critA): " << eta(Im(-x1),Im(critA),Im(-x1*critA)) << " vs " << (ln(-x1*critA)-ln(-x1)-ln(critA))/(2*pi<T>()*complex<T>(0,1)) << endl;
     cout << "eta(-x_1,critB): " << eta(Im(-x1),Im(critB),Im(-x1*critB)) << " vs " << (ln(-x1*critB)-ln(-x1)-ln(critB))/(2*pi<T>()*complex<T>(0,1)) << endl;
     cout << "eta(-x_2,critA): " << eta(Im(-x2),Im(critA),Im(-x2*critA)) << " vs " << (ln(-x2*critA)-ln(-x2)-ln(critA))/(2*pi<T>()*complex<T>(0,1)) << endl;
     cout << "eta(-x_2,critB): " << eta(Im(-x2),Im(critB),Im(-x2*critB)) << " vs " << (ln(-x2*critB)-ln(-x2)-ln(critB))/(2*pi<T>()*complex<T>(0,1)) << endl;
#endif
     etaTerms =
       eta(Im(-x1),Im(critA),Im(-x1*critA))*ln(arg1a)
       //          *complex<T>(ln(abs(arg1a)),Theta(-arg1a,pi<T>()*Sign(critA)))
       -eta(Im(-x2),Im(critA),Im(-x2*critA))*ln(arg2a)
       //          *complex<T>(ln(abs(arg2a)),Theta(-arg2a,pi<T>()*Sign(critA)))
       +eta(Im(-x1),Im(critB),Im(-x1*critB))*ln(arg1b)
       //          *complex<T>(ln(abs(arg1b)),Theta(-arg1b,pi<T>()*Sign(critB)))
       -eta(Im(-x2),Im(critB),Im(-x2*critB))*ln(arg2b);
       //       *complex<T>(ln(abs(arg2b)),Theta(-arg2b,pi<T>()*Sign(critB)));
     etaTerms *= complex<T>(0,T(2)*pi<T>());
     LiTerms = Li2(arg1a)-Li2(arg2a)+Li2(arg1b)-Li2(arg2b);
#if DumpTerms
     cout << "Li terms: " << LiTerms/(s*t*rho) << endl;
     cout << "eta terms: " << etaTerms/(s*t*rho) << endl;
#endif
     ln1 = ln(-x1);  ln2 = ln(-x2);
#if DumpTerms
     cout << "ln1: " << ln1 << endl;
     cout << "ln2: " << ln2 << endl;
#endif
     lnTerms = (ln1-ln2)*((ln1+ln2)/T(2) - ln(-m1sq-Ieps) - ln(-s-Ieps)
                          + ln(-m4sq-Ieps) + ln(-m2sq-Ieps));
#if DumpTerms
     cout << "ln terms: " << lnTerms/(s*t*rho) << endl;
     cout << "  ln^2 term: " <<
       ((ln1-ln2)*((ln1+ln2)/T(2))/(s*t*rho)) << endl;
     cout << "  ln term1: " << ((ln1-ln2)*(-ln(-m1sq-Ieps))/(s*t*rho)) << endl;
     cout << "  ln term2: " << ((ln1-ln2)*(-ln(-s-Ieps))/(s*t*rho)) << endl;
     cout << "  ln term3: " << ((ln1-ln2)*(ln(-m4sq-Ieps))/(s*t*rho)) << endl;
     cout << "  ln term4: " << ((ln1-ln2)*(ln(-m2sq-Ieps))/(s*t*rho)) << endl;
#endif

#if DumpTag
   cout << tag;
#endif
   return( T(1)/norm * (LiTerms + etaTerms + lnTerms) );}
   else {// rhoSq < 0
     tag = "rho^2 < 0 "+masstag;
     complex<T> rho = sqrt(complex<T>(rhoSq,0));  // Yes, rho can be imaginary...

#if DumpArguments
   cout << "rho^2 < 0" << endl;
   cout << "m^2_i,s,t: ";
   cout << m1sq << " ";
   cout << m2sq << " ";
   cout << m3sq << " ";
   cout << m4sq << " ";
   cout << s << " ";
   cout << t;
   cout << endl;
#endif

   if (false) {// 2007 version
     /* The code below seems to give a real (as expected) answer but
        with totally wrong value */
#if 0
     cout << "imaginary rho" << endl;
#endif

     complex<T> arg1a = (1-lambda1+lambda2+rho)/T(2);
     complex<T> arg2a = (1-lambda1+lambda2-rho)/T(2);
     complex<T> arg1b = -(1-lambda1-lambda2-rho)/lambda1/T(2);
     complex<T> arg2b = -(1-lambda1-lambda2+rho)/lambda1/T(2);
     complex<T> x1 = s*(arg1a-T(1))/m3sq, x2 = s*(arg2a-T(1))/m3sq;

     complex<T> ln1 = ln(x1),
       ln2 = ln(x2);

#if DumpTerms
     cout << "arg1a: " << arg1a << endl;
     cout << "arg2a: " << arg2a << endl;
     cout << "arg1b: " << arg1b << endl;
     cout << "arg2b: " << arg2b << endl;
#endif

     return( T(1)/(s*t*rho) *
           (Li2(arg1a)-Li2(arg2a)+Li2(arg1b)-Li2(arg2b)
            /* eta terms -- specialize to only 2nd arg having infinitesimal
               imaginary part */
            +eta(-Im(x1),s-m3sq,-m3sq*Im(x1)/s)*ln(arg1a)
            -eta(-Im(x2),s-m3sq,-m3sq*Im(x2)/s)*ln(arg2a)
            +eta(-Im(x1),m1sq-t,-t*Im(x1)/m1sq)*ln(arg1b)
            -eta(-Im(x2),m1sq-t,-t*Im(x2)/m1sq)*ln(arg2b)
            // ln terms
            +(ln1-ln2)*((ln1+ln2)/T(2)
                        -complex<T>(ln(abs(m1sq)),Theta(m1sq,-pi<T>()))
                        -complex<T>(ln(abs(s)),Theta(s,-pi<T>()))
                        +complex<T>(ln(abs(m4sq)),Theta(m4sq,-pi<T>()))
                        +complex<T>(ln(abs(m2sq)),Theta(m2sq,-pi<T>())))) );}
   else {// 2008 versions
   complex<T> LiTerms, etaTerms, lnTerms;
   complex<T> Ieps = complex<T>(0,epsilon<T>());
   complex<T> lambda1, lambda2, arg1a, arg2a, arg1b, arg2b;
   complex<T> x1, x2;
   complex<T> ln1, ln2;
   complex<T> critA, critB;
   complex<T> norm;

   if (version is 1)
      {// 2007 -- with I e in m^2 & s, t, x_i updated
        lambda1 = (m1sq+Ieps)*(m3sq+Ieps)/((s+Ieps)*(t+Ieps));
        lambda2 = (m2sq+Ieps)*(m4sq+Ieps)/((s+Ieps)*(t+Ieps));
        rho = sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2);

        //     complex<T> x1 = (s+Ieps)*(arg1a-1.)/(m3sq+Ieps)+Ieps*m2sq,
        x1 = (s)*(-T(1)-lambda1+lambda2+rho)/(T(2)*m3sq)
          +Ieps*m2sq*s*t/(T(4)*rhoSq);
        x2 = (s)*(-T(1)-lambda1+lambda2-rho)/(T(2)*m3sq)
          -Ieps*m2sq*s*t/(T(4)*rhoSq);
        norm = s*t*rho;}
   else if (version is 2)
      {// 2008 version, only I eps as in Denner et al, none in m_i^2
        lambda1 = complex<T>(m1sq*m3sq/(s*t),0);
        lambda2 = complex<T>(m2sq*m4sq/(s*t),0);
        rho = sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2);
        x1 = -s/(T(2)*m3sq) - m1sq/(T(2)*t) + (m2sq*m4sq)/(T(2)*m3sq*t) +
          Sqrt(square(m1sq*m3sq) - T(2)*m1sq*m2sq*m3sq*m4sq +
               square(m2sq*m4sq) + T(4)*Ieps*m2sq*m3sq*t -
               T(2)*m1sq*m3sq*s*t - T(2)*m2sq*m4sq*s*t
               + square(s*t))/(T(2)*m3sq*t);
        x2 = -s/(T(2)*m3sq) - m1sq/(T(2)*t) + (m2sq*m4sq)/(T(2)*m3sq*t) -
          Sqrt(square(m1sq*m3sq) - T(2)*m1sq*m2sq*m3sq*m4sq +
               square(m2sq*m4sq) + T(4)*Ieps*m2sq*m3sq*t -
               T(2)*m1sq*m3sq*s*t - T(2)*m2sq*m4sq*s*t
               + square(s*t))/(T(2)*m3sq*t);
        norm = t*m3sq*(x1-x2);}
   else if (version is /* old 2 */ -2 && false)
      {// 2008 version, only I eps as in Denner et al, none in m_i^2
        lambda1 = complex<T> (m1sq*m3sq/(s*t),0);
        lambda2 = complex<T>(m2sq*m4sq/(s*t),0);
        rho = sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2);
        x1 = (s)*(-T(1)-lambda1+lambda2+rho)/(T(2)*(m3sq+Ieps))
          +Ieps*m2sq*s*t/(T(4)*rhoSq);
        x2 = (s)*(-T(1)-lambda1+lambda2-rho)/(T(2)*(m3sq+Ieps))
          -Ieps*m2sq*s*t/(T(4)*rhoSq);}
   else if (version is 3)
      {// Use Davydychev--Ussyukina form (PLB 298:262 (1993)), without Ieps
        T s1 = m1sq*m3sq, s2 = m2sq*m4sq, s3 = s*t,
          // Why is there a sign?
          prefactor = T(-1);
        T Delta3 = -Delta(s1,s2,s3);
        complex<T> result(0,0);
        T delta[3] = {s2+s3-s1,s1+s3-s2,s1+s2-s3};
        complex<T> sqrtDelta;
#if DumpTerms
        cout << "Delta3: " << Delta3 << endl;
        cout << "signs: " << Sign(s1) << " " << Sign(s2) << " " << Sign(s3) << endl;
#endif
        if (Delta3 >= 0)
           {T sqrtDelta = sqrt(Delta3);
             /* In this case, the arguments of the dilogs are generic complex
                numbers, and we don't need to do anything special */
             for (int j = 0;  j < 3;  j += 1)
                {complex<T> arg = (sqrtDelta+complex<T>(0,delta[j]))/(-sqrtDelta+complex<T>(0,delta[j]));
                  result += Li2(arg)-Li2(T(1)/arg);}
             result *= complex<T>(0,-1)/sqrtDelta;}
        else {// rho^2 < 0 -> Delta3 > 0, so this branch isn't reached
          /* Here, the arguments of the dilogs become real, and we have to
             be careful about staying on the right side of the branch cut.
             On the other hand, which sign we choose for sqrt(Delta3) is
             irrelevant. */
          cout << "Shouldn't be here!!!" << endl;
          T ImSqrtDelta = sqrt(-Delta3);
#define Crit(s1,s2,s3) (s1*s2 - s2*s2 + s1*s3 - s3*s3)
          T crit[3] = {Crit(s1,s2,s3),Crit(s2,s3,s1),Crit(s3,s1,s2)};
          for (int j = 0;  j < 3;  j += 1)
             {T arg = (ImSqrtDelta+delta[j])/(-ImSqrtDelta+delta[j]);
               result += CLi2(arg,crit[j])-CLi2(T(1)/arg,-crit[j]);}
          result /= -ImSqrtDelta;}
#if DumpTag
        cout << tag;
#endif
        return (prefactor*result);}
   else if (version is 4)
      {// 2008 version, I eps as in Denner et al, but also as in m_i^2
        complex<T> M1sq(m1sq,epsilon<T>()), M2sq(m2sq,epsilon<T>()),
          M3sq(m3sq,epsilon<T>()), M4sq(m4sq,epsilon<T>());
        complex<T> S(s,epsilon<T>()), Ti(t,epsilon<T>());

        lambda1 = M1sq*M3sq/(S*Ti);
        lambda2 = M2sq*M4sq/(S*Ti);
        rho = sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2);
        x1 = -S/(T(2)*M3sq) - M1sq/(T(2)*t) + (M2sq*M4sq)/(T(2)*M3sq*t) +
          Sqrt(square(M1sq*M3sq) - T(2)*M1sq*M2sq*M3sq*M4sq +
               square(M2sq*M4sq) + T(4)*Ieps*M2sq*M3sq*t -
               T(2)*M1sq*M3sq*S*Ti - T(2)*M2sq*M4sq*S*Ti
               + square(S*Ti))/(T(2)*M3sq*Ti);
        x2 = -S/(T(2)*M3sq) - M1sq/(T(2)*t) + (M2sq*M4sq)/(T(2)*M3sq*t) -
          Sqrt(square(M1sq*M3sq) - T(2)*M1sq*M2sq*M3sq*M4sq +
               square(M2sq*M4sq) + T(4)*Ieps*M2sq*M3sq*t -
               T(2)*M1sq*M3sq*S*Ti - T(2)*M2sq*M4sq*S*Ti
               + square(S*Ti))/(T(2)*M3sq*Ti);
        norm = t*m3sq*(x1-x2);
      }
   else if (version is 5)
      {// Use Davydychev--Ussyukina form (PLB 298:262 (1993)), with Ieps in m^2
        complex<T> I(0,1);
        complex<T> s1 = complex<T>(m1sq,epsilon<T>())*complex<T>(m3sq,epsilon<T>()),
          s2 = complex<T>(m2sq,epsilon<T>())*complex<T>(m4sq,epsilon<T>()),
          s3 = complex<T>(s,epsilon<T>())*complex<T>(t,epsilon<T>());
        // Why is there a sign?
        T prefactor = T(-1);
        complex<T> Delta3 = -Delta(s1,s2,s3);
        complex<T> result(0,0);
        complex<T> delta[3] = {s2+s3-s1,s1+s3-s2,s1+s2-s3};
#if DumpArguments
        cout << "Delta3: " << Delta3 << endl;
        cout << "signs: " << Sign(Re(s1)) << " " << Sign(Re(s2)) << " " << Sign(Re(s3)) << endl;
#endif
        complex<T> sqrtDelta = sqrt(Delta3);
        /* In this case, the arguments of the dilogs are generic complex
           numbers, and we don't need to do anything special */
        for (int j = 0;  j < 3;  j += 1)
           {complex<T> arg = (sqrtDelta+I*delta[j])/(-sqrtDelta+I*delta[j]);
             result += Li2(arg)-Li2(T(1)/arg);}
        result *= complex<T>(0,-1)/sqrtDelta;
#if DumpTag
        cout << tag;
#endif
        return (prefactor*result);}
   else if (version is 6)
      {// 2008 version, formulae from Denner et al, but I eps only in m_i^2
        complex<T> M1sq(m1sq,epsilon<T>()), M2sq(m2sq,epsilon<T>()),
          M3sq(m3sq,epsilon<T>()), M4sq(m4sq,epsilon<T>());
        complex<T> S(s,epsilon<T>()), Ti(t,epsilon<T>());

        lambda1 = M1sq*M3sq/(S*Ti);
        lambda2 = M2sq*M4sq/(S*Ti);
        rho = sqrt(square(T(1)-lambda1-lambda2)-T(4)*lambda1*lambda2);
        x1 = -S/(T(2)*M3sq) - M1sq/(T(2)*t) + (M2sq*M4sq)/(T(2)*M3sq*t) +
          Sqrt(square(M1sq*M3sq) - T(2)*M1sq*M2sq*M3sq*M4sq +
               square(M2sq*M4sq) -
               T(2)*M1sq*M3sq*S*Ti - T(2)*M2sq*M4sq*S*Ti
               + square(S*Ti))/(T(2)*M3sq*Ti);
        x2 = -S/(T(2)*M3sq) - M1sq/(T(2)*t) + (M2sq*M4sq)/(T(2)*M3sq*t) -
          Sqrt(square(M1sq*M3sq) - T(2)*M1sq*M2sq*M3sq*M4sq +
               square(M2sq*M4sq) -
               T(2)*M1sq*M3sq*S*Ti - T(2)*M2sq*M4sq*S*Ti
               + square(S*Ti))/(T(2)*M3sq*Ti);
        norm = t*m3sq*(x1-x2);
      }
#if DumpArguments
   cout << "-x1: " << (-x1) << ", " << "-x2: " << (-x2) << endl;
#endif
   // 2nd arguments to eta()
   critA = (m3sq+Ieps)/(s+Ieps);
   critB = (t+Ieps)/(m1sq+Ieps);
   // Arguments of 1st polylog (and corresponding ln in eta*ln) for k = 1,2 in eq. (41)
   arg1a = T(1)+critA*x1;
   arg2a = T(1)+critA*x2;
   // Arguments of 2nd polylog (and corresponding ln in eta*ln) for k = 1,2 in eq. (41)
   arg1b = T(1)+critB*x1;
   arg2b = T(1)+critB*x2;

#if DumpArguments
   cout << "arg1a: " << arg1a << "; ";
   cout << "arg2a: " << arg2a << "; ";
   cout << "arg1b: " << arg1b << "; ";
   cout << "arg2b: " << arg2b << "; ";
   cout << endl;

   cout << "critA: " << critA << ", "
        << "critB: " << critB << endl;
   cout << "eta(-x_1,critA): " << eta(Im(-x1),Im(critA),Im(-x1*critA)) << endl;
   cout << "eta(-x_1,critB): " << eta(Im(-x1),Im(critB),Im(-x1*critB)) << endl;
   cout << "eta(-x_2,critA): " << eta(Im(-x2),Im(critA),Im(-x2*critA)) << endl;
   cout << "eta(-x_2,critB): " << eta(Im(-x2),Im(critB),Im(-x2*critB)) << endl;
#endif
     etaTerms =
       eta(Im(-x1),Im(critA),Im(-x1*critA))*ln(arg1a)
       //          *complex<T>(ln(abs(arg1a)),Theta(-arg1a,pi<T>()*Sign(critA)))
       -eta(Im(-x2),Im(critA),Im(-x2*critA))*ln(arg2a)
       //          *complex<T>(ln(abs(arg2a)),Theta(-arg2a,pi<T>()*Sign(critA)))
       +eta(Im(-x1),Im(critB),Im(-x1*critB))*ln(arg1b)
       //          *complex<T>(ln(abs(arg1b)),Theta(-arg1b,pi<T>()*Sign(critB)))
       -eta(Im(-x2),Im(critB),Im(-x2*critB))*ln(arg2b);
       //       *complex<T>(ln(abs(arg2b)),Theta(-arg2b,pi<T>()*Sign(critB)));
     etaTerms *= complex<T>(0,T(2)*pi<T>());
     LiTerms = Li2(arg1a)-Li2(arg2a)+Li2(arg1b)-Li2(arg2b);
#if DumpArguments
     cout << "s t rho: " << (s*t*rho) << endl;
     cout << "a (x1-x2): " << (t*m3sq*(x1-x2)) << endl;
#endif
#if DumpTerms
     cout << "Li terms: " << LiTerms/(s*t*rho) << endl;
     cout << "eta terms: " << etaTerms/(s*t*rho) << endl;
#endif
     ln1 = ln(-x1);  ln2 = ln(-x2);
#if DumpTerms
     cout << "ln1: " << ln1 << endl;
     cout << "ln2: " << ln2 << endl;
#endif
     lnTerms = (ln1-ln2)*((ln1+ln2)/T(2) - ln(-m1sq-Ieps) - ln(-s-Ieps)
                          + ln(-m4sq-Ieps) + ln(-m2sq-Ieps));
#if DumpTerms
     cout << "ln terms: " << lnTerms/(s*t*rho) << endl;
     cout << "  ln^2 term: " <<
       ((ln1-ln2)*((ln1+ln2)/T(2))/(s*t*rho)) << endl;
     cout << "  ln term1: " << ((ln1-ln2)*(-ln(-m1sq-Ieps))/(s*t*rho)) << endl;
     cout << "  ln term2: " << ((ln1-ln2)*(-ln(-s-Ieps))/(s*t*rho)) << endl;
     cout << "  ln term3: " << ((ln1-ln2)*(ln(-m4sq-Ieps))/(s*t*rho)) << endl;
     cout << "  ln term4: " << ((ln1-ln2)*(ln(-m2sq-Ieps))/(s*t*rho)) << endl;
#endif

#if DumpTag
   cout << tag;
#endif
   return( T(1)/norm * (LiTerms + etaTerms + lnTerms) );}
   }
   }
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
 }
}


// Triangle functions
/* Note that to get the physical (full) integral from the I (non-script)
   functions in (e.g.) hep-ph/9409265v1, we must multiply by (-1)^n, which
   gives an extra minus sign in the triangles.  That extra minus sign
   is NOT included here (changed from inclusion on 2/6/08)
*/

// 11/20/08
// Direct versions of triangle functions, arguments are invariants rather
// than indices; will eventually replace above functions
template <class T> complex<T> I3w1m(int order,
				    const T& musq, const T& s)
{switch (order) {
 case -2:
   return(-T(1)/s);
 case -1:
   return(T(1)/s*CLnM(s,musq));
 case 0:
   return(-1/(T(2)*s)*square(CLnM(s,musq)));
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
   }}

template <class T> complex<T> I3w2m(int order,
				    const T& musq, const T& s1, const T& s2)
{switch (order) {
 case -2:
   return complex<T>(0,0);
 case -1:
   return(T(1)/(s2-s1)*(CLnM(s2,musq)-CLnM(s1,musq)));
 case 0:
   return(-T(1)/(T(2)*(s2-s1))*(square(CLnM(s2,musq))-square(CLnM(s1,musq))));
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
   }}

// Requires fully complex dilog
template <class T> complex<T> I3w3m(int order,
				    const T& musq, const T& s1, const T& s2,
				    const T& s3)
{switch (order) {
 case -2:
 case -1:
   return complex<T>(0,0);
 case 0:
   {T Delta3 = -Delta(s1,s2,s3);
   complex<T> result(0,0);
   T delta[3] = {s2+s3-s1,s1+s3-s2,s1+s2-s3};
   complex<T> sqrtDelta;
#if 0
   cout << "Delta3: " << Delta3 << endl;
   cout << "signs: " << Sign(s1) << " " << Sign(s2) << " " << Sign(s3) << endl;
#endif
   if (Delta3 >= 0)
      {T sqrtDelta = sqrt(Delta3);
      /* In this case, the arguments of the dilogs are generic complex
         numbers, and we don't need to do anything special */
      for (int j = 0;  j < 3;  j += 1)
         {complex<T> arg = (sqrtDelta+complex<T>(0,delta[j]))/(-sqrtDelta+complex<T>(0,delta[j]));
         result += Li2(arg)-Li2(T(1)/arg);}
      result *= complex<T>(0,-1)/sqrtDelta;}
   else {
     /* Here, the arguments of the dilogs become real, and we have to
        be careful about staying on the right side of the branch cut.
        On the other hand, which sign we choose for sqrt(Delta3) is
        irrelevant. */
     T ImSqrtDelta = sqrt(-Delta3);
     #define Crit(s1,s2,s3) (s1*s2 - s2*s2 + s1*s3 - s3*s3)
     T crit[3] = {Crit(s1,s2,s3),Crit(s2,s3,s1),Crit(s3,s1,s2)};
     for (int j = 0;  j < 3;  j += 1)
         {T arg = (ImSqrtDelta+delta[j])/(-ImSqrtDelta+delta[j]);
         result += CLi2(arg,crit[j])-CLi2(T(1)/arg,-crit[j]);}
      result /= -ImSqrtDelta;}
#if 0
   cout << "result: " << result << endl;
#endif
   return(result);}
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   else return complex<T>(0,0);
#endif
 }}



// 11/20/08
// Direct version of bubble function, arguments are invariants rather
// than indices; will eventually replace above functions
template <class T> complex<T> I2(int order,
				 const T& musq, const T& s)
{switch (order) {
 case -2:
   return(complex<T>(0,0));
 case -1:
   return complex<T>(1,0);
 case 0:
   return(-CLnM(s,musq)+T(2));
 default:
   if (order < -2) return complex<T>(0,0);
#if _USE_GCC
   else return(std::numeric_limits<C>::signaling_NaN());
#else
   return complex<T>(0,0);
#endif

 }}



// Older spurious-denominator-free functions

template <class T> inline int IsMassless(const eval_param<T> & k, const vector<int>& v, const vector<int>& masses)
  // Returns 1 if "v" represents a massless vector (as an index into "k")
 // Test for masslessness should be replaced by examination of a flag
  // when one becomes available
{	if (v.size() == 1 && masses[v[0]]==0) return 1;
 return 0;}

template <class T> inline int IsMassless(const eval_param<T> & k, const vector<int>& v)
  // Returns 1 if "v" represents a massless vector (as an index into "k")
 // Test for masslessness should be replaced by examination of a flag
  // when one becomes available
{	if (v.size() == 1) return 1;
 return 0;}


#define Massless(k,v) 1

/*
string GenKey(char* tag,int,int,const vector<int>&);
string GenKey(char* tag,int,int,const vector<int>&,const vector<int>&);
string GenKey(char* tag,int,int,const vector<int>&,const vector<int>&,
              const vector<int>&);
string GenKey(char* tag,int,int,const vector<int>&,const vector<int>&,
              const vector<int>&,const vector<int>&);
*/

// Box wrappers
template <class T> complex<T> IntM(int order,
      const eval_param<T> & k, const T& mu,
      const vector<int>& corner1,
      const vector<int>& corner2,
      const vector<int>& corner3,
      const vector<int>& corner4, const vector<int>& masses)
{
	 complex<T> result;
#if IntCaching
	 string key = GenKey(string("Ia"),order,mu,corner1,corner2,corner3,corner4);
 if ( ! k.get_value(key,result) ){
#endif

		 int pattern = (IsMassless(k,corner1,masses) << 3)
		   |(IsMassless(k,corner2,masses) << 2)
		   |(IsMassless(k,corner3,masses) << 1)
		   |IsMassless(k,corner4,masses);


		 // cout << "pattern: " << hex << pattern << dec << endl;

// Direct-invariant routines 11/20/08
   /* Use these to avoid overhead implicit in k.Sum, don't save intermediate
      momenta coming from summing momenta at different corners */
   switch (pattern) {
   case 0x0:
     result = I4w4m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner1),SqSum(k,corner2),SqSum(k,corner3),SqSum(k,corner4));
     break;
   case 0x1:
     result = I4w3m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner1),SqSum(k,corner2),SqSum(k,corner3));
     break;
   case 0x2:
     result = I4w3m(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                    SqSum(k,corner4),SqSum(k,corner1),SqSum(k,corner2));
     break;
   case 0x4:
     result = I4w3m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner3),SqSum(k,corner4),SqSum(k,corner1));
     break;
   case 0x8:
     result = I4w3m(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                    SqSum(k,corner2),SqSum(k,corner3),SqSum(k,corner4));
     break;
   case 0x3:
     result = I4w2mh(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                     SqSum(k,corner1),SqSum(k,corner2));
     break;
   case 0x6:
     result = I4w2mh(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                     SqSum(k,corner4),SqSum(k,corner1));
     break;
   case 0xC:
     result = I4w2mh(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                     SqSum(k,corner3),SqSum(k,corner4));
     break;
   case 0x9:
     result = I4w2mh(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                     SqSum(k,corner2),SqSum(k,corner3));
     break;
   case 0x5:
     result = I4w2me(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                     SqSum(k,corner1),SqSum(k,corner3));
     break;
   case 0xA:
     result = I4w2me(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                     SqSum(k,corner2),SqSum(k,corner4));
     break;
   case 0x7:
     result = I4w1m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner1));
     break;
   case 0xE:
     result = I4w1m(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                    SqSum(k,corner4));
     break;
   case 0xD:
     result = I4w1m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner3));
     break;
   case 0xB:
     result = I4w1m(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                    SqSum(k,corner2));
     break;
   case 0xF:
     result = I4w0m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3));
     break;
   }

#if IntCaching
   k.put_value(key,result);
#if SaveSymmetry
   // By rotation
   key = GenKey("Ia",order,mu,corner2,corner3,corner4,corner1);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner3,corner4,corner1,corner2);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner4,corner1,corner2,corner3);
   k.put_value(key,result);
   // And by reflection
   key = GenKey("Ia",order,mu,corner4,corner3,corner2,corner1);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner3,corner2,corner1,corner4);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner2,corner1,corner4,corner3);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner1,corner4,corner3,corner2);
   k.put_value(key,result);
#endif
   }
#endif
 return result;
}

template <class T> complex<T> Int(int order,
      const eval_param<T> & k, const T& mu,
      const vector<int>& corner1,
      const vector<int>& corner2,
      const vector<int>& corner3,
      const vector<int>& corner4)
{
	 complex<T> result;
#if IntCaching
	 string key = GenKey(string("Ia"),order,mu,corner1,corner2,corner3,corner4);
 if ( ! k.get_value(key,result) ){
#endif

		 int pattern = (IsMassless(k,corner1) << 3)
		   |(IsMassless(k,corner2) << 2)
		   |(IsMassless(k,corner3) << 1)
		   |IsMassless(k,corner4);


		 // cout << "pattern: " << hex << pattern << dec << endl;

// Direct-invariant routines 11/20/08
   /* Use these to avoid overhead implicit in k.Sum, don't save intermediate
      momenta coming from summing momenta at different corners */
   switch (pattern) {
   case 0x0:
     result = I4w4m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner1),SqSum(k,corner2),SqSum(k,corner3),SqSum(k,corner4));
     break;
   case 0x1:
     result = I4w3m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner1),SqSum(k,corner2),SqSum(k,corner3));
     break;
   case 0x2:
     result = I4w3m(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                    SqSum(k,corner4),SqSum(k,corner1),SqSum(k,corner2));
     break;
   case 0x4:
     result = I4w3m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner3),SqSum(k,corner4),SqSum(k,corner1));
     break;
   case 0x8:
     result = I4w3m(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                    SqSum(k,corner2),SqSum(k,corner3),SqSum(k,corner4));
     break;
   case 0x3:
     result = I4w2mh(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                     SqSum(k,corner1),SqSum(k,corner2));
     break;
   case 0x6:
     result = I4w2mh(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                     SqSum(k,corner4),SqSum(k,corner1));
     break;
   case 0xC:
     result = I4w2mh(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                     SqSum(k,corner3),SqSum(k,corner4));
     break;
   case 0x9:
     result = I4w2mh(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                     SqSum(k,corner2),SqSum(k,corner3));
     break;
   case 0x5:
     result = I4w2me(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                     SqSum(k,corner1),SqSum(k,corner3));
     break;
   case 0xA:
     result = I4w2me(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                     SqSum(k,corner2),SqSum(k,corner4));
     break;
   case 0x7:
     result = I4w1m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner1));
     break;
   case 0xE:
     result = I4w1m(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                    SqSum(k,corner4));
     break;
   case 0xD:
     result = I4w1m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3),
                    SqSum(k,corner3));
     break;
   case 0xB:
     result = I4w1m(order,mu,SqSum(k,corner4,corner1),SqSum(k,corner1,corner2),
                    SqSum(k,corner2));
     break;
   case 0xF:
     result = I4w0m(order,mu,SqSum(k,corner1,corner2),SqSum(k,corner2,corner3));
     break;
   }

#if IntCaching
   k.put_value(key,result);
#if SaveSymmetry
   // By rotation
   key = GenKey("Ia",order,mu,corner2,corner3,corner4,corner1);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner3,corner4,corner1,corner2);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner4,corner1,corner2,corner3);
   k.put_value(key,result);
   // And by reflection
   key = GenKey("Ia",order,mu,corner4,corner3,corner2,corner1);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner3,corner2,corner1,corner4);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner2,corner1,corner4,corner3);
   k.put_value(key,result);
   key = GenKey("Ia",order,mu,corner1,corner4,corner3,corner2);
   k.put_value(key,result);
#endif
   }
#endif
 return result;
}


#define HighestPole -2
#define FiniteTerm 0
template <class T> Series<complex<T> >
IntM(const eval_param<T> & k, const T& mu,
            const vector<int>& corner1,
            const vector<int>& corner2,
            const vector<int>& corner3,
            const vector<int>& corner4, const vector<int>& masses)
{vector<complex<T> > coeffs(FiniteTerm-HighestPole+1);

 for (int j = HighestPole;  j <= FiniteTerm;  j += 1)
   coeffs[j-HighestPole] = IntM(j,k,mu,corner1,corner2,corner3,corner4,masses);
 Series<complex<T> > result(HighestPole,FiniteTerm,coeffs);
 // cout << "Int 4: " << result << endl;
 return result;}

template <class T> Series<complex<T> >
Int(const eval_param<T> & k, const T& mu,
            const vector<int>& corner1,
            const vector<int>& corner2,
            const vector<int>& corner3,
            const vector<int>& corner4)
{vector<complex<T> > coeffs(FiniteTerm-HighestPole+1);

 for (int j = HighestPole;  j <= FiniteTerm;  j += 1)
   coeffs[j-HighestPole] = Int(j,k,mu,corner1,corner2,corner3,corner4);
 Series<complex<T> > result(HighestPole,FiniteTerm,coeffs);
 // cout << "Int 4: " << result << endl;
 return result;}



// Triangle wrappers
template <class T> complex<T> IntM(int order,
      const eval_param<T> & k, const T& mu,
      const vector<int>& corner1,
      const vector<int>& corner2,
      const vector<int>& corner3, const vector<int>& masses)
{

	 complex<T> result;
#if IntCaching
	 string key = GenKey(string("Ib"),order,mu,corner1,corner2,corner3);
	 if ( !k.get_value(key,result) ){
#endif
		int pattern = (IsMassless(k,corner1,masses) << 2)
			|(IsMassless(k,corner2,masses) << 1)
			|IsMassless(k,corner3,masses);

 // Direct-invariant routines 11/20/08
   /* Use these to avoid overhead implicit in k.Sum, don't save intermediate
      momenta coming from summing momenta at different corners */
   switch (pattern) {
   case 0x0:
     result = I3w3m(order,mu,SqSum(k,corner1),SqSum(k,corner2),SqSum(k,corner3));
     break;
   case 0x1:
     result = I3w2m(order,mu,SqSum(k,corner1),SqSum(k,corner2));
     break;
   case 0x2:
     result = I3w2m(order,mu,SqSum(k,corner3),SqSum(k,corner1));
     break;
   case 0x4:
     result = I3w2m(order,mu,SqSum(k,corner2),SqSum(k,corner3));
     break;
   case 0x3:
     result = I3w1m(order,mu,SqSum(k,corner1));
     break;
   case 0x5:
     result = I3w1m(order,mu,SqSum(k,corner2));
     break;
   case 0x6:
     result = I3w1m(order,mu,SqSum(k,corner3));
     break;
   }
#if IntCaching
   k.put_value(key,result);
#if SaveSymmetry
   // By rotation
   key = GenKey("Ib",order,mu,corner2,corner3,corner1);
   k.put_value(key,result);
   key = GenKey("Ib",order,mu,corner3,corner1,corner2);
   k.put_value(key,result);
   // And by reflection
   key = GenKey("Ib",order,mu,corner3,corner2,corner1);
   k.put_value(key,result);
   key = GenKey("Ib",order,mu,corner2,corner1,corner3);
   k.put_value(key,result);
   key = GenKey("Ib",order,mu,corner1,corner3,corner2);
   k.put_value(key,result);
#endif
   }
#endif

 return result;
 }

template <class T> complex<T> Int(int order,
      const eval_param<T> & k, const T& mu,
      const vector<int>& corner1,
      const vector<int>& corner2,
      const vector<int>& corner3)
{

	 complex<T> result;
#if IntCaching
	 string key = GenKey(string("Ib"),order,mu,corner1,corner2,corner3);
	 if ( !k.get_value(key,result) ){
#endif
		int pattern = (IsMassless(k,corner1) << 2)
			|(IsMassless(k,corner2) << 1)
			|IsMassless(k,corner3);

 // Direct-invariant routines 11/20/08
   /* Use these to avoid overhead implicit in k.Sum, don't save intermediate
      momenta coming from summing momenta at different corners */
   switch (pattern) {
   case 0x0:
     result = I3w3m(order,mu,SqSum(k,corner1),SqSum(k,corner2),SqSum(k,corner3));
     break;
   case 0x1:
     result = I3w2m(order,mu,SqSum(k,corner1),SqSum(k,corner2));
     break;
   case 0x2:
     result = I3w2m(order,mu,SqSum(k,corner3),SqSum(k,corner1));
     break;
   case 0x4:
     result = I3w2m(order,mu,SqSum(k,corner2),SqSum(k,corner3));
     break;
   case 0x3:
     result = I3w1m(order,mu,SqSum(k,corner1));
     break;
   case 0x5:
     result = I3w1m(order,mu,SqSum(k,corner2));
     break;
   case 0x6:
     result = I3w1m(order,mu,SqSum(k,corner3));
     break;
   }
#if IntCaching
   k.put_value(key,result);
#if SaveSymmetry
   // By rotation
   key = GenKey("Ib",order,mu,corner2,corner3,corner1);
   k.put_value(key,result);
   key = GenKey("Ib",order,mu,corner3,corner1,corner2);
   k.put_value(key,result);
   // And by reflection
   key = GenKey("Ib",order,mu,corner3,corner2,corner1);
   k.put_value(key,result);
   key = GenKey("Ib",order,mu,corner2,corner1,corner3);
   k.put_value(key,result);
   key = GenKey("Ib",order,mu,corner1,corner3,corner2);
   k.put_value(key,result);
#endif
   }
#endif

 return result;
 }



template <class T> Series<complex<T> >
IntM(const eval_param<T> & k, const T& mu,
            const vector<int>& corner1,
            const vector<int>& corner2,
            const vector<int>& corner3, const vector<int>& masses)
{vector<complex<T> > coeffs(FiniteTerm-HighestPole+1);

 for (int j = HighestPole;  j <= FiniteTerm;  j += 1)
   coeffs[j-HighestPole] = IntM(j,k,mu,corner1,corner2,corner3,masses);
 Series<complex<T> > result(HighestPole,FiniteTerm,coeffs);
 return result;}

template <class T> Series<complex<T> >
Int(const eval_param<T> & k, const T& mu,
            const vector<int>& corner1,
            const vector<int>& corner2,
            const vector<int>& corner3)
{vector<complex<T> > coeffs(FiniteTerm-HighestPole+1);

 for (int j = HighestPole;  j <= FiniteTerm;  j += 1)
   coeffs[j-HighestPole] = Int(j,k,mu,corner1,corner2,corner3);
 Series<complex<T> > result(HighestPole,FiniteTerm,coeffs);
 return result;}


// Bubble wrapper
template <class T> complex<T> IntM(int order,
      const eval_param<T> & k, const T& mu,
      const vector<int>& corner1,
      const vector<int>& corner2, const vector<int>& masses)
{
	 complex<T> result;
#if IntCaching
	string key = GenKey(string("Ib"),order,mu,corner1,corner2);
 if ( ! k.get_value(key,result) ){
#endif
#if !Direct // original code 11/20/08
   result = I2(order,k,mu,k.Sum(corner1));;
#else // Direct-invariant routines 11/20/08
   /* Use these to avoid overhead implicit in k.Sum, don't save intermediate
      momenta coming from summing momenta at different corners */
   result = I2(order,mu,SqSum(k,corner1));
#endif
#if IntCaching
   k.put_value(key,result);
#if SaveSymmetry
   // By rotation or reflection
   key = GenKey("Ib",order,mu,corner2,corner1);
   k.put_value(key,result);
#endif
   }
#endif
 return result;
}

template <class T> complex<T> Int(int order,
      const eval_param<T> & k, const T& mu,
      const vector<int>& corner1,
      const vector<int>& corner2)
{
	 complex<T> result;
#if IntCaching
	string key = GenKey(string("Ib"),order,mu,corner1,corner2);
 if ( ! k.get_value(key,result) ){
#endif
#if !Direct // original code 11/20/08
   result = I2(order,k,mu,k.Sum(corner1));;
#else // Direct-invariant routines 11/20/08
   /* Use these to avoid overhead implicit in k.Sum, don't save intermediate
      momenta coming from summing momenta at different corners */
   result = I2(order,mu,SqSum(k,corner1));
#endif
#if IntCaching
   k.put_value(key,result);
#if SaveSymmetry
   // By rotation or reflection
   key = GenKey("Ib",order,mu,corner2,corner1);
   k.put_value(key,result);
#endif
   }
#endif
 return result;
}


template <class T> complex<T> IntM(int order,
      const eval_param<T> & k, const T& mu,
      const vector<int>& corner1, const vector<int>& masses)
{
	 complex<T> result;
#if IntCaching
	 string key = GenKey(string("Ib"),order,mu,corner1);
 if ( ! k.get_value(key,result) ){
#endif

	   result = I2(order,mu,SqSum(k,corner1));
#if IntCaching
   k.put_value(key,result);}
#endif

 return result;
}

template <class T> complex<T> Int(int order,
      const eval_param<T> & k, const T& mu,
      const vector<int>& corner1)
{
	 complex<T> result;
#if IntCaching
	 string key = GenKey(string("Ib"),order,mu,corner1);
 if ( ! k.get_value(key,result) ){
#endif

	   result = I2(order,mu,SqSum(k,corner1));
#if IntCaching
   k.put_value(key,result);}
#endif

 return result;
}


template <class T> Series<complex<T> >
IntM(const eval_param<T> & k, const T& mu,
            const vector<int>& corner1,
            const vector<int>& corner2, const vector<int>& masses)
{vector<complex<T> > coeffs(FiniteTerm-HighestPole+1);

 for (int j = HighestPole;  j <= FiniteTerm;  j += 1)
   coeffs[j-HighestPole] = IntM(j,k,mu,corner1,corner2,masses);
 Series<complex<T> > result(HighestPole,FiniteTerm,coeffs);
 return result;}

template <class T> Series<complex<T> >
Int(const eval_param<T> & k, const T& mu,
            const vector<int>& corner1,
            const vector<int>& corner2)
{vector<complex<T> > coeffs(FiniteTerm-HighestPole+1);

 for (int j = HighestPole;  j <= FiniteTerm;  j += 1)
   coeffs[j-HighestPole] = Int(j,k,mu,corner1,corner2);
 Series<complex<T> > result(HighestPole,FiniteTerm,coeffs);
 return result;}


template <class T> Series<complex<T> >
IntM(const eval_param<T> & k, const T& mu,
            const vector<int>& corner1, const vector<int>& masses)
{vector<complex<T> > coeffs(FiniteTerm-HighestPole+1);

 for (int j = HighestPole;  j <= FiniteTerm;  j += 1)
   coeffs[j-HighestPole] = IntM(j,k,mu,corner1,masses);
 Series<complex<T> > result(HighestPole,FiniteTerm,coeffs);
 return result;}

template <class T> Series<complex<T> >
Int(const eval_param<T> & k, const T& mu,
            const vector<int>& corner1)
{vector<complex<T> > coeffs(FiniteTerm-HighestPole+1);

 for (int j = HighestPole;  j <= FiniteTerm;  j += 1)
   coeffs[j-HighestPole] = Int(j,k,mu,corner1);
 Series<complex<T> > result(HighestPole,FiniteTerm,coeffs);
 return result;}

//EXPLICIT INSTANTITION

template Series<C> IntM(const eval_param<R>& k , const R& mu, const vector<int>& corner1,const vector<int>& corner2, const vector<int>& masses);
template Series<CHP> IntM(const eval_param<RHP>& k , const RHP& mu, const vector<int>& corner1,const vector<int>& corner2, const vector<int>& masses);
template Series<CVHP> IntM(const eval_param<RVHP>& k , const RVHP& mu, const vector<int>& corner1,const vector<int>& corner2, const vector<int>& masses);
template Series<C> IntM(const eval_param<R>& k , const R& mu, const vector<int>& corner1,const vector<int>& corner2,const vector<int>& corner3, const vector<int>& masses);
template Series<CHP> IntM(const eval_param<RHP>& k , const RHP& mu, const vector<int>&,const vector<int>&,const vector<int>&, const vector<int>& masses);
template Series<CVHP> IntM(const eval_param<RVHP>& k , const RVHP& mu, const vector<int>&,const vector<int>&,const vector<int>&, const vector<int>& masses);
template Series<C> IntM(const eval_param<R>& k , const R& mu, const vector<int>& corner1,const vector<int>& corner2,const vector<int>& corner3,const vector<int>&, const vector<int>& masses);
template Series<CHP> IntM(const eval_param<RHP>& k , const RHP& mu, const vector<int>&,const vector<int>&,const vector<int>&,const vector<int>&, const vector<int>& masses);
template Series<CVHP> IntM(const eval_param<RVHP>& k , const RVHP& mu, const vector<int>&,const vector<int>&,const vector<int>&,const vector<int>&, const vector<int>& masses);

template Series<C> Int(const eval_param<R>& k , const R& mu, const vector<int>& corner1,const vector<int>& corner2);
template Series<CHP> Int(const eval_param<RHP>& k , const RHP& mu, const vector<int>& corner1,const vector<int>& corner2);
template Series<CVHP> Int(const eval_param<RVHP>& k , const RVHP& mu, const vector<int>& corner1,const vector<int>& corner2);
template Series<C> Int(const eval_param<R>& k , const R& mu, const vector<int>& corner1,const vector<int>& corner2,const vector<int>& corner3);
template Series<CHP> Int(const eval_param<RHP>& k , const RHP& mu, const vector<int>&,const vector<int>&,const vector<int>&);
template Series<CVHP> Int(const eval_param<RVHP>& k , const RVHP& mu, const vector<int>&,const vector<int>&,const vector<int>&);
template Series<C> Int(const eval_param<R>& k , const R& mu, const vector<int>& corner1,const vector<int>& corner2,const vector<int>& corner3,const vector<int>&);
template Series<CHP> Int(const eval_param<RHP>& k , const RHP& mu, const vector<int>&,const vector<int>&,const vector<int>&,const vector<int>&);
template Series<CVHP> Int(const eval_param<RVHP>& k , const RVHP& mu, const vector<int>&,const vector<int>&,const vector<int>&,const vector<int>&);

#if BH_USE_GMP
template Series<CGMP> IntM(const eval_param<RGMP>& k , const RGMP& mu, const vector<int>& corner1,const vector<int>& corner2, const vector<int>& masses);
template Series<CGMP> IntM(const eval_param<RGMP>& k , const RGMP& mu, const vector<int>&,const vector<int>&,const vector<int>&, const vector<int>& masses);
template Series<CGMP> IntM(const eval_param<RGMP>& k , const RGMP& mu, const vector<int>&,const vector<int>&,const vector<int>&,const vector<int>&, const vector<int>& masses);


template Series<CGMP> Int(const eval_param<RGMP>& k , const RGMP& mu, const vector<int>& corner1,const vector<int>& corner2);
template Series<CGMP> Int(const eval_param<RGMP>& k , const RGMP& mu, const vector<int>&,const vector<int>&,const vector<int>&);
template Series<CGMP> Int(const eval_param<RGMP>& k , const RGMP& mu, const vector<int>&,const vector<int>&,const vector<int>&,const vector<int>&);

#endif

}
