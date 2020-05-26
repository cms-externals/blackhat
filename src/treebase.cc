/* treebase.cc */

/*  David A. Kosower, June 5, 2008

Entry point for tree routines.

*/

#define BRENDSGIELE_H_
#define BERENDSGIELE_IMPL_H_
#include <vector>
#include "spinor.h"
#include "mom_conf.h"
#include "Tree.h"
#include "process.h"

#define isnt !=
#define is ==

using namespace std;

namespace BH {

namespace Tree {
#include "zero.h"

    void PrintVector(vector<particle_ID> v);

template<class T> complex<T> A(momentum_configuration<T>& k,
                               const vector<int>& arg0 /* indices of arguments */,
                               const vector<particle_ID>& leg,
                               int algorithm,
                               // Electroweak vectors
                               const vector<int>& vectorK /* momenta */,
                               const vector<int>& polarization,
                               const vector<int>& coupleTo /* quark flavor */,
                              const vector<int>& massValue)
{vector<int> arg(leg.size()); // "leg" determines the number of arguments
 for (int i = 0;  i < arg.size();  i += 1) arg[i] = arg0[i];
 // cout << "A: " << leg.size() << ": ";
 // PrintVector(leg);  cout << endl;
if (algorithm is 0)
  return (J(k,-1,leg.back(),arg,leg,0,leg.size()-2,massValue.back(),massValue));
 else if (algorithm is 1)
   return (Aosrr(k,arg,leg,vectorK,polarization,coupleTo,massValue));
 else if (algorithm is 2)
    {int mass = (massValue.size() > 0 and massValue.back() >= 0
                 and !IsZero(k.m2(massValue[massValue.back()]))) ?
        massValue.back() : defaultMass;
      return (Jc(k,-1,leg.back(),arg,leg,0,leg.size()-2,0,
                 vectorK,polarization,coupleTo,mass,massValue));}
}

template<class T> complex<T> A(momentum_configuration<T>& k,
                               const vector<int>& arg0 /* indices of arguments */,
                               const process& p,
                               int algorithm,
                               // Electroweak vectors implicit in "p" & "arg0"
                               //                               const vector<int>& vectorK /* momenta */ = empty,
                               //                               const vector<int>& polarization = empty,
                               //                               const vector<int>& coupleTo /* quark flavor */ = empty,
                               const vector<int>& massValue)
{vector<particle_ID> leg;
  vector<int> vectorK(0), polarization(0), coupleTo(0);
  // Assume leptons are paired sequentially, all couple to 1st quark flavor
 int firstQuarkFlavor = -1;
 for (int i = 1;  i <= p.n();  i += 1) // indices run from 1..n
    if (p.p(i).is_a(quark)) {firstQuarkFlavor = p.p(i).flavor(); break;}

 int priorLepton = -1;
 for (int i = 1;  i <= p.n();  i += 1)
    {if (p.p(i).is_a(lepton))
       if (priorLepton is -1) priorLepton = arg0[i-1];
       else {vectorK.push_back(k.Sum(priorLepton,arg0[i-1]));
         if (p.p(i).helicity() > 0)
           polarization.push_back(k.insert(k.L(priorLepton),k.Lt(arg0[i-1])));
         else polarization.push_back(k.insert(k.L(arg0[i-1]),k.Lt(priorLepton)));
         coupleTo.push_back(firstQuarkFlavor);
         priorLepton = -1;}
      else leg.push_back(p.p(i));
     }

 return A(k,arg0,leg,algorithm,vectorK,polarization,coupleTo,massValue);}

// Explicit instantiation
template class complex<R> A(momentum_configuration<R>& k,
                            const vector<int>& arg0 /* indices of arguments */,
                            const vector<particle_ID>& leg,
                            int algorithm,
                            // Electroweak vectors
                            const vector<int>& vectorK /* momenta */ ,
                            const vector<int>& polarization ,
                            const vector<int>& coupleTo /* quark flavor */,
                            const vector<int>& massValue);

template class complex<R> A(momentum_configuration<R>& k,
                            const vector<int>& arg0 /* indices of arguments */,
                            const process& p,
                            int algorithm,
                            // Electroweak vectors
                            //                            const vector<int>& vectorK /* momenta */ = empty,
                            //                            const vector<int>& polarization = empty,
                            //                            const vector<int>& coupleTo /* quark flavor */ = empty,
                            const vector<int>& massValue );

}}


