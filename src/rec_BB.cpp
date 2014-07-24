/*
 * Rec_BB.cpp
 *
 *  Created on: 1 Jul 2009
 *      Author: daniel
 */

#include "rec_BB.h"

namespace BH {

Rec_BB::~Rec_BB(){
	for (int i=0;i<daughters.size();i++){
		delete daughters[i];
	}
}


}
