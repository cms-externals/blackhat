/*
 * momenta_evaluator_ratext.h
 *
 *  Created on: 21 Aug 2009
 *      Author: darrenforde
 */

#ifndef MOMENTA_EVALUATOR_RATEXT_H_
#define MOMENTA_EVALUATOR_RATEXT_H_

#include "BH_typedefs.h"

namespace BH {

class Rec_Tree;

namespace ratext {
	
#define MUBOXPOINTS_HIGGS 3
#define MUTRIPOINTS_HIGGS 3
#define MUBUBPOINTS_HIGGS 2
#define CPOINTS_HIGGS 9
#define CBUBPOINTS_HIGGS 4
#define YPOINTS_HIGGS 4

class MomentumEvaluator {
public:
	MomentumEvaluator(){};
	~MomentumEvaluator(){};
};

} /* ratext */

} /* BH */


#endif /* MOMENTA_EVALUATOR_RATEXT_H_ */
