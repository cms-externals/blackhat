/*!\file BH_particle_names.h
\brief Header for the particle names
*/
#ifndef BH_PARTICLE_NAMES_H_
#define BH_PARTICLE_NAMES_H_

namespace BH {
//!   particle types
enum particle_type {
	q, 		/*!< quark */
	qb, 	/*!< anti-quark */
	q2,		/*!< 2nd flavor quark */
	qb2,	/*!< 2nd flavor anti-quark */
	l,		/*!< lepton */
	l2,		/*!< 2nd flavor lepton */
	g,		/*!< gluon */
	sc,		/*!< scalar */
	qM,		/*!< massive quark */
	qbM,		/*!< massive antiquark */
	undef_particle		/*!< undefined particle */
}; 



//!   particle-helicity types
/**
the particle-helicity types are used as short-cuts in constructors. 
 */
enum ph_type {
	qp,  		/*!< positive helicity quark */
	qm,   		/*!< negative helicity quark */
	qbp,   		/*!< positive helicity anti-quark */
	qbm,   		/*!< negative helicity anti-quark */
	q2p,   		/*!< positive helicity quark of the 2nd type */
	q2m,   		/*!< positive helicity quark of the 2nd type */
	qb2p,   	/*!< positive helicity quark of the 2nd type */
	qb2m,   	/*!< positive helicity quark of the 2nd type */
	lp,	    	/*!< positive helicity lepton  */
	lm,	    	/*!< negative helicity lepton  */
	l2p,    	/*!< positive helicity lepton of the 2nd type */
	l2m,    	/*!< negative helicity lepton of the 2nd type */
	p,  	  	/*!< positive helicity gluon */
	m,   		/*!< negative helicity gluon */
	scp,  	  	/*!< positive "helicity" scalar */
	scm,   		/*!< negative "helicity" scalar */
	qMp,  		/*!< positive helicity massive quark */
	qMm,   		/*!< negative helicity massive quark */
	qbMp,   		/*!< positive helicity massive anti-quark */
	qbMm,   		/*!< negative helicity massive anti-quark */
	undef_ph_type   		/*!< undefined particle */
}; 

}


#endif /*BH_PARTICLE_NAMES_H_*/
