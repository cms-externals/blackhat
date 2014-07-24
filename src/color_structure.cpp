/*
 * color_structure.cpp
 *
 *  Created on: 22-Apr-2009
 *      Author: daniel
 */

#include "BH_typedefs.h"
#include <iostream>

namespace BH {

std::ostream& operator<<(std::ostream& s, color_structure cs){
			switch (cs){
			case nf: return s << "nf" ;
			case nf_top: return s << "nf_top" ;
			case LT:return s << "LT" ;
			case RT:return s << "RT" ;
			case leading_color:return s << "leading_color" ;
			case sub_leading_color:return s << "sub_leading_color" ;
			case slc_q:return s << "slc_q" ;
			case slc_G:return s << "slc_G" ;
			case glue: return s << "glue" ;
			case VECT: return s << "VECT" ;
			case AX: return s << "AX" ;
			case AXSL: return s << "AXSL" ;
			case LLT:return s<< "LLT";
			case RRT:return s<< "RRT";
			case LRT:return s<< "LRT";
			case RLT:return s<< "RLT";
			case LLLT:return s<< "LLLT";
			case RLLT:return s<< "RLLT";
			case LRLT:return s<< "LRLT";
			case LLRT:return s<< "LLRT";
			case RRLT:return s<< "RRLT";
			case RLRT:return s<< "RLRT";
			case LRRT:return s<< "LRRT";
			case RRRT:return s<< "RRRT";
			case nfLT:return s<< "nfLT";
			case nfRT:return s<< "nfRT";
			case nfLLT:return s<< "nfLLT";
			case nfRRT:return s<< "nfRRT";
			case nfLRT:return s<< "nfLRT";
			case nfRLT:return s<< "nfRLT";
			case nfLLLT:return s<< "nfLLLT";
			case nfRLLT:return s<< "nfRLLT";
			case nfLRLT:return s<< "nfLRLT";
			case nfLLRT:return s<< "nfLLRT";
			case nfRRLT:return s<< "nfRRLT";
			case nfRLRT:return s<< "nfRLRT";
			case nfLRRT:return s<< "nfLRRT";
			case nfRRRT:return s<< "nfRRRT";
			case unspecified: return s << "unspecified" ;
			}
}

}
