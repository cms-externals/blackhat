/*!\file factorization.h
\brief Header file for fact_momenta.cpp, splitamp_tree.cpp and _loop.cpp
       and splittingtest.cpp
       implementing factorization and collinear kinematics and
       tree and loop splitting amplitudes. 
*/

// by CFB, last revision 11/18/2008
// 9/5/08 new collinear kinematics following Catani-Seymour
// 10/20/08 coll regeneration with different inv mass from existing coll kinematics
// 10/27/08 different factories now allowed at loop level
// 11/18/08 new color structures for primitive versus color-ordered structures
//          cf. 9403226 and 9409393

#ifndef FACTORIZATION_H_
#define FACTORIZATION_H_

#include <cstdlib> 
#include <math.h>
#include "mom_conf.h"
#include "particles.h"
#include "BH_utilities.h"
#include "qd_suppl.h"
#include "integrals.h"
#include "rational.h"
#include "cut_part_factory.h"
#include "OneLoopHelAmpl.h"
using namespace std;


namespace BH {

enum p_in_loop { s_in_loop, f_in_loop, v_in_loop, f_left, f_right, leading_c, ferm_c, sub_leading_c };
// note: s_in_loop = [0], f_in_loop = [1/2], v_in_loop = [1] (cf. 9403226)
//       f_left = [L], f_right = [R], [s] = [0] = s_in_loop, [f] = ferm_c (cf. 9409393)
//       primitive versus color-ordered amplitudes

//! collkinematicsCS returns a configuration of n particles with a || b
/** \param invmass: P = a+b has invariant mass sqrt(invmass) (not too small!)
 *  \param a and b indices of collinear momenta
 *  \param oldCS: old coll configuration where the collinear momenta are recomputed
 *  \n with a different new invariant mass
 *  \n momentum P will be given as the (n+1)th momentum
 *  \n the shifted new momentum k = b+1 for convenience will be
 *  \n given as the (n+2)nd momentum
 *  \n This version follows Catani-Seymour so that momentum conservation
 *  \n is maintained also for the merged momenta -- which necessitates
 *  \n modification of momentum k as well.
 */
template <class T> momentum_configuration<T>collkinematicsCSredo(
		momentum_configuration<T>& oldCS,
		int a, int b, T z, T invmass);

//! collkinematicsCS returns a configuration of n particles with a || b
/** \param zfrac: fraction of momentum carried by a (b carries (1-z))
 *  \param invmass: P = a+b has invariant mass sqrt(invmass) (not too small!)
 *  \param a and b indices of collinear momenta
 *  \n momentum P will be given as the (n+1)th momentum
 *  \n the shifted new momentum k = b+1 for convenience will be
 *  \n given as the (n+2)nd momentum
 *  \n This version follows Catani-Seymour so that momentum conservation
 *  \n is maintained also for the merged momenta -- which necessitates
 *  \n modification of momentum k as well.
 */
template <class T> momentum_configuration<T>collkinematicsCS(int n,
		int a, int b, T z, T invmass);

//! collkinematics returns a configuration of n particles with a || b
/** \param zfrac: fraction of momentum carried by a (b carries (1-z))
 *  \param invmass: P = a+b has invariant mass sqrt(invmass)
 *  \param a and b indices of collinear momenta
 *  \n momentum P will be given as the (n+1)th momentum
 */
template <class T> momentum_configuration<T> collkinematics(int n, int a, int b, 
		T zfrac, T invmass);


template <class T> complex<T> dott(Cmom<T> Ki, Cmom<T> Kj);		
		
//! softkinematics returns a configuration of n particles with a, b, c soft
/** \param invmass: P = a+b+c has invariant mass sqrt(invmass)
 *  \param a, b, c indices of soft momenta
 *  \n momentum P will be given as the (n+1)th momentum
 */
template <class T> momentum_configuration<T> softkinematics(int n, int a, int b, int c, 
		T invmass);

//! soft4kinematics returns a configuration of n particles with a, b, c, d soft
/** \param invmass: P = a+b+c+d has invariant mass sqrt(invmass)
 *  \param a, b, c, d indices of soft momenta
 *  \n momentum P will be given as the (n+1)th momentum
 */
template <class T> momentum_configuration<T> soft4kinematics(int n, int a, int b, int c, int d,
		T invmass);

//! decay decays a momentum Ktot into two null momenta a and b returned as a vector
template <class T> vector<Cmom<T> > decay(Cmom<T> Ktot, int trial);


//! random momentum generator
/** \param signE: sign of time-like component
 *  \param mass: 0 for massless particles
 */
template <class T> Cmom<T> randmom(T mass, short signE);

//! random double number generator, between 0 and 1
template <class T> T randdouble(void);

template <class T> complex<T> Split0(const process& p, momentum_configuration<T>& ps, 
		int inda, int indb);
//! Split0 returns the tree splitting amplitude 
/** \param p is the process ab -> c, in the order a, b, c
 *  \param ps is the momentum configuration of type mom_conf
 *  \param inda is the index of particle a in momentum configuration ps
 *  \param indb is the index of particle b in momentum configuration ps
 *  Note that this assumes that ps is set up such that a = z P, and 
 *  \n b = (1-z) P, where P = a + b. User responsibility!! 
 *  \n This is automatically the case if ps is generated via
 *  collkinematics()
 */

// internal functions for the tree splitting amplitudes
// ****************************************************
template <class T> complex<T> Sggg(const process& p,momentum_configuration<T>& ps,
		int inda, int indb);
// gg -> g

template <class T> complex<T> Sqqg(const process& p,momentum_configuration<T>& ps,
		int inda, int indb);
// q qbar -> g

template <class T> complex<T> Sqgq(const process& p,momentum_configuration<T>& ps,
		int inda, int indb);
// q g -> q

template <class T> complex<T> Sgqq(const process& p,momentum_configuration<T>& ps,
		int inda, int indb);
// g qbar -> qbar

//long scode(const process& p);
// returns an enumeration/encoding of the process in one long
// not needed any longer due to new particle ID

template <class T> Series<complex<T> > SplitS1(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb, T mu2, int);
//! SplitS1 is a series-wrapper for the loop splitting amplitude Split1

template <class T> complex<T> Split1(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2,int);
//! Split1 returns the unrenormalized(!) loop splitting amplitude 
/** \param p is the process ab -> c, in the order a, b, c
 *  \param in_loop is the particle in the loop: scalar, fermion or vector,
 *  \n or in the case of q splitting amplitudes, whether it's a left- or
 *     right-moving fermion
 *  \n Note that for the right-moving fermions we DO NOT include a factor
 *     of -1/Nc^2 (external quarks), which in some cases (external gluinos)
 *     is => +1  
 *  \param ps is the momentum configuration of type mom_conf
 *  \param inda is the index of particle a in momentum configuration ps
 *  \param indb is the index of particle b in momentum configuration ps
 *  \param oeps is coefficient of which order in eps to be returned (-2,-1,0)
 *  \param mu2 is the scale mu^2
 *  Note that this assumes that ps is set up such that a = z P, and 
 *  \n b = (1-z) P, where P = a + b. User responsibility!! 
 *  \n This is automatically the case if ps is generated via
 *  collkinematics()
 *  \n Note also that the regularization scheme is hardcoded, and
 *  currently set to the four-dimensional helicity scheme. \n
 *  This can be changed by recompiling the code with the constant
 *  DIMREG set to 1.
 *  \n And finally, note that a factor of c_\Gamma = (4 pi)^eps/(16 pi^2)
 *  Gamma(1+eps) Gamma^2(1-eps)/Gamma(1-2eps) has been omitted.
 */

/*template <class T> Series<complex<T> > SplitS1_Cut(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb, T mu2);
//! SplitS1_Cut is a series-wrapper for the loop splitting amplitude Split1_Cut

template <class T> complex<T> Split1_Cut(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
*/
//! Split1_Cut returns the unrenormalized(!) cut part of the loop splitting amplitude 
/** \param p is the process ab -> c, in the order a, b, c
 *  \param in_loop is the particle in the loop: scalar, fermion or vector,
 *  \n or in the case of q splitting amplitudes, whether it's a left- or
 *     right-moving fermion
 *  \n Note that for the right-moving fermions we DO NOT include a factor
 *     of -1/Nc^2 (external quarks), which in some cases (external gluinos)
 *     is => +1  
 *  \param ps is the momentum configuration of type mom_conf
 *  \param inda is the index of particle a in momentum configuration ps
 *  \param indb is the index of particle b in momentum configuration ps
 *  \param oeps is coefficient of which order in eps to be returned (-2,-1,0)
 *  \param mu2 is the scale mu^2
 *  Note that this assumes that ps is set up such that a = z P, and 
 *  \n b = (1-z) P, where P = a + b. User responsibility!! 
 *  \n This is automatically the case if ps is generated via
 *  collkinematics()
 *  \n Note also that the regularization scheme is hardcoded, and
 *  currently set to the four-dimensional helicity scheme. \n
 *  This can be changed by recompiling the code with the constant
 *  DIMREG set to 1.
 *  \n And finally, note that a factor of c_\Gamma = (4 pi)^eps/(16 pi^2)
 *  Gamma(1+eps) Gamma^2(1-eps)/Gamma(1-2eps) has been omitted.
 */

/*template <class T> complex<T> Split1_Rat(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb);
*/
//! Split1_Rat returns the unrenormalized(!) ratl. part of the loop splitting amplitude 
/** \param p is the process ab -> c, in the order a, b, c
 *  \param in_loop is the particle in the loop: scalar, fermion or vector,
 *  \n or in the case of q splitting amplitudes, whether it's a left- or
 *     right-moving fermion
 *  \n Note that for the right-moving fermions we DO NOT include a factor
 *     of -1/Nc^2 (external quarks), which in some cases (external gluinos)
 *     is => +1  
 *  \param ps is the momentum configuration of type mom_conf
 *  \param inda is the index of particle a in momentum configuration ps
 *  \param indb is the index of particle b in momentum configuration ps
 *  Note that this assumes that ps is set up such that a = z P, and 
 *  \n b = (1-z) P, where P = a + b. User responsibility!! 
 *  \n This is automatically the case if ps is generated via
 *  collkinematics()
 *  \n Note also that the regularization scheme is hardcoded, and
 *  currently set to the four-dimensional helicity scheme. \n
 *  This can be changed by recompiling the code with the constant
 *  DIMREG set to 1.
 *  \n And finally, note that a factor of c_\Gamma = (4 pi)^eps/(16 pi^2)
 *  Gamma(1+eps) Gamma^2(1-eps)/Gamma(1-2eps) has been omitted.
 */

// internal functions for the loop splitting amplitudes
// ****************************************************
template <class T> complex<T> Sggg1(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
// gg -> g

template <class T> complex<T> Sggg1s(const process& p, momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
// gg -> g, scalar in the loop, note the fermion one is - the scalar one

template <class T> complex<T> Sggg1v(const process& p, momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
// gg -> g, vector in the loop, adjoint (= gluon), note this is not simply = [0]

template <class T> complex<T> Sqgq1(const process& p, int in_loop, 
		momentum_configuration<T>& ps, 
		int inda, int indb, int oeps, T mu2);
// qg -> q, particles in loop obvious

template <class T> complex<T> Sgqq1(const process& p, int in_loop, 
		momentum_configuration<T>& ps, 
		int inda, int indb, int oeps, T mu2);
// gqbar -> qbar, particles in loop obvious

// internal functions for the loop splitting amplitudes
// ****************************************************
template <class T> complex<T> Sggg1_Cut(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
// gg -> g
template <class T> complex<T> Sggg1_Rat(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb);
// gg -> g

template <class T> complex<T> Sggg1s_Cut(const process& p, momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
// gg -> g, scalar in the loop, note the fermion one is - the scalar one
template <class T> complex<T> Sggg1s_Rat(const process& p, momentum_configuration<T>& ps,
		int inda, int indb);
// gg -> g, scalar in the loop, note the fermion one is - the scalar one

template <class T> complex<T> Sggg1v_Cut(const process& p, momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
// gg -> g, vector in the loop, adjoint (= gluon), note this is not simply = [0]
template <class T> complex<T> Sggg1v_Rat(const process& p, momentum_configuration<T>& ps,
		int inda, int indb);
// gg -> g, vector in the loop, adjoint (= gluon), note this is not simply = [0]


template <class T> complex<T> Sqgq1_Cut(const process& p, int in_loop, 
		momentum_configuration<T>& ps, 
		int inda, int indb, int oeps, T mu2);
// qg -> q, particles in loop obvious
template <class T> complex<T> Sqgq1_Rat(const process& p, int in_loop, 
		momentum_configuration<T>& ps, 
		int inda, int indb);
// qg -> q, particles in loop obvious

template <class T> complex<T> Sgqq1_Cut(const process& p, int in_loop, 
		momentum_configuration<T>& ps, 
		int inda, int indb, int oeps);
// gqbar -> qbar, particles in loop obvious
template <class T> complex<T> Sgqq1_Rat(const process& p, int in_loop, 
		momentum_configuration<T>& ps, 
		int inda, int indb);
// gqbar -> qbar, particles in loop obvious

template <class T> complex<T> ff(complex<T> z, int in_loop, complex<T> s, T mu2, int oeps);
// auxiliary function

template <class T> complex<T> Sqqg1(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
// q qbar -> g

template <class T> complex<T> Sqqg1_Cut(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb, int oeps, T mu2);
// q qbar -> g

template <class T> complex<T> Sqqg1_Rat(const process& p, int in_loop, 
		momentum_configuration<T>& ps,
		int inda, int indb);
// q qbar -> g

template <class T> T li2(T x);
// Lance's dilog function for real arguments

template <class T> complex<T> Clog(complex<T> s1, T mu2);
// log with correct branch (i Pi)

int TestCollinear(const process& pro,color_structure cs, int a, int b, int full_cut_rat=0);
template <class CUT,class RAT> int TestCollinear(const process& pro,color_structure cs, 
		Rational_factory<RAT>* RRF, cut_part_factory<CUT>* CPF, int a, int b);
int TestCollinearTree(const process& pro,color_structure cs, int a, int b,int full_cut_rat=0);

template <class T,class RAT,class CUT> int TestCollinear(const process& pro,color_structure cs, 
		Rational_factory<RAT>* RRF, cut_part_factory<CUT>* CPF,	T mass,int full_cut_rat=0);
template <class T> int TestCollinear(const process& pro,color_structure cs, T mass,double tolerance,int full_cut_rat=0);
//! TestCollinear tests all collinear limits for the given process and prints out a message
/** \param pro is the process
 *  \param cs is the color structure for the process (currently only glue and leading_color)
 *  \param mass is the collinear "mass" to be specified by the user
 */

/*template <class T, class RAT, class CUT> int TestCollFixedKin(const process& pro,color_structure cs,
        Rational_factory<RAT>* RRF, cut_part_factory<CUT>* CPF,
        momentum_configuration<T>& collkin, int a, int b, T z, T mass);
*/
//! TestCollFixedKin tests the collinear limits for the given process and prints out a message
/** \param pro is the process
 *  \param cs is the color structure for the process (currently only glue and leading_color)
 *  \param coll is an original already constructed collinear configuration which is to be modified
 *  \param a,b are the two legs that become collinear
 *  \param z is the collinear momentum fraction of P_a = z(P_a + P_b)
 *  \param mass is the collinear "mass" to be specified by the user
 */


//template <class T> int TestCollinearCut(const process& pro,color_structure cs, T mass);
//! TestCollinearCut tests all collinear limits for the given process and prints out a message
//! only for the Cut part
/** \param pro is the process
 *  \param mass is the collinear "mass" to be specified by the user
 */ 

//template <class T> int TestCollinearRat(const process& pro, color_structure cs,T mass);
//! TestCollinearRat tests all collinear limits for the given process and prints out a message
//! only for the rational part
/** \param pro is the process
 *  \param mass is the collinear "mass" to be specified by the user
 */ 

template <class T> int TestCollinearTree(const process& pro, T mass);
//! TestCollinearTree tests all collinear limits for the given process and prints out a message
//! at tree level
/** \param pro is the process
 *  \param mass is the collinear "mass" to be specified by the user
 */ 

//template <class T> int TestCollFixedTree(const process& pro, 
//		momentum_configuration<T>& collkin, int a, int b, T z, T mass);


}
#endif /*FACTORIZATION_H_*/
