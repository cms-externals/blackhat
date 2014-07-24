/*
 * rec_BB.h
 *
 *  Created on: 1 Jul 2009
 *      Author: daniel
 */

#ifndef REC_BB_H_
#define REC_BB_H_

#include "computable.h"

namespace BH {



//! class for building blocks of recursion
class Rec_BB : public computable<std::complex> {
protected:
	std::vector<Rec_BB*> daughters;
public:
#ifdef SWIG
	Rec_BB(){};
#endif
	Rec_BB* get_daughter(size_t i){return daughters[i-1];}
	void set_daughter(size_t i,Rec_BB* new_bb){daughters[i-1]=new_bb;}
	size_t nbr_daughters(){return daughters.size();}
	virtual ~Rec_BB();
};

}

#endif /* REC_BB_H_ */
