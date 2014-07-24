/*
 * Known_Cut_wCI.cpp
 *
 *  Created on: 23 Sep 2009
 *      Author: daniel
 */

#include "cut/Known_Cut_wCI.h"
#include "process.h"
#include "helcode.h"
#include "process_utils.h"


#include "cut/C_2q1g1y_wCI.h"
#include "cut/C_2q2g1y_wCI.h"
#include "cut/C_2q1g2l_wCI.h"
#include "cut/C_2q2g2l_wCI.h"
#include "cut/C_2q2Gl2l_wCI.h"
#include "cut/C_2q2Q_wCI.h"
#include "cut/C_2q2Q1g_wCI.h"
#include "cut/C_2q3g_wCI.h"
#include "cut/C_5g_wCI.h"
//#include "cut/C_2q2l_wCI.h"
#include "cut/C_2q2Q2l_wCI.h"
#include "cut/C_2q2Q1y_wCI.h"
#include "cut/C_6g_wCI.h"
#include "cut/C_4g_wCI.h"
#include "cut/C_2q2g_wCI.h"

#include "cut/C_3g1ph_wCI.h"
#include "cut/C_4g1ph_wCI.h"
#include "cut/C_2q2g1ph_wCI.h"
#include "cut/C_2q2Q1ph_wCI.h"

namespace BH {

class process;

namespace CachedIntegral {

Cut_Part_wCI* CwCI_Ptr(const process& pro,color_structure cs,const std::vector<int>& ind){
//   std::cout << "in CwCI_Ptr  " << pro << pro.pcode() << "  " << cs << std::endl;

	switch(pro.pcode()){


	case 100000003: // Higgs plus 3 gluons
		switch (cs) {
		case glue: return CwCI_3g1ph_G(helcode_Ng1ph(pro),ind);
		}

	case 100000004: // Higgs plus 4 gluons
		switch (cs) {
		case glue: return CwCI_4g1ph_G(helcode_Ng1ph(pro),ind);
		}

	case 100000022: // Higgs 2 quark 2 gluon
		switch (cs) {
		case LT: return CwCI_2q2g1ph_LT(helcode_phi_1q(pro),ind);
		case RT: return CwCI_2q2g1ph_RT(helcode_phi_1q(pro),ind);
		case nfLT: return CwCI_2q2g1ph_nfLT(helcode_phi_1q(pro),ind);
		}
	case 100000040: // Higgs 4 quark
		switch (cs) {
		case leading_color: return CwCI_2q2Q1ph_lc(helcode_phi_2q2Q(pro),ind);
		case sub_leading_color: return CwCI_2q2Q1ph_slc(helcode_phi_2q2Q(pro),ind);
		case nf: return CwCI_2q2Q1ph_nf(helcode_phi_2q2Q(pro),ind);
		}
	case 4: // added by KJO
		switch (cs) {
		case glue: return CwCI_4g_G(helcode_g(pro),ind);
		case nf: return CwCI_4g_nf(helcode_g(pro),ind);
		}
	case 5:
        switch (cs) {
        case glue: return CwCI_5g_G(helcode_g(pro),ind);
        case nf: return CwCI_5g_nf(helcode_g(pro),ind);
        }
	case 6: // added by KJO
		switch (cs) {
		case glue: return CwCI_6g_G(helcode_g(pro),ind);
		case nf: return CwCI_6g_nf(helcode_g(pro),ind);
		}
	case 22: // added by KJO
		switch (cs) {
		case LT: return CwCI_2q2g_LT(helcode_2q(pro),ind);
		case nfLT: return CwCI_2q2g_nfLT(helcode_2q(pro),ind);
		}
	case 40:{
		switch (cs) {  // added by KJO
		case LLT: return CwCI_2q2Q_LLT(helcode_2q2Q(pro),ind);
		case LRT: return CwCI_2q2Q_LRT(helcode_2q2Q(pro),ind);
		case nfLLT: return CwCI_2q2Q_nfLLT(helcode_2q2Q(pro),ind);
		case nfLRT: return CwCI_2q2Q_nfLRT(helcode_2q2Q(pro),ind);
// 		case leading_color: return CachedIntegral::CwCI_2q2Q_L(helcode_2q2Q(fix_flavors(pro)),ind);
// 		case sub_leading_color: return CachedIntegral::CwCI_2q2Q_SLC(helcode_2q2Q(fix_flavors(pro)),ind);
// 		case nf: return CachedIntegral::CwCI_2q2Q_nf(helcode_2q2Q(fix_flavors(pro)),ind);
		}
		break;
	}
	case 23:
		switch (cs) {
		case LT: return CwCI_2q3g_L(helcode_2q(pro),ind);
//		case sub_leading_color: return C2q3g_SLC_Ptr<T>(helcode_2q(pro));
		case nfLT: return CwCI_2q3g_nf(helcode_2q(pro),ind);
		case nf: return CwCI_2q3g_nf(helcode_2q(pro),ind);
		// case nf: return CwCI_2q3g_L(helcode_2q(pro),ind);
		}
	case 41: {
		switch (cs) {
		case LLT: return CwCI_2q2Q1g_LLT(helcode_2q2Q(pro),ind);
		case nfLLT: return CwCI_2q2Q1g_nfLLT(helcode_2q2Q(pro),ind);
		case LRT: return CwCI_2q2Q1g_LRT(helcode_2q2Q(pro),ind);
// 		case leading_color: return CwCI_2q2Q1g_PentP(helcode_2q2Q(fix_flavors(pro)),ind);
	    } break;
	}
//	case 220:
//		switch (cs) {
//			case leading_color: return C2q2l_L_Ptr<T>(helcode_2q2l(pro));
//	}
//	break;
	case 221: {
		switch (cs) {
		case leading_color: return CachedIntegral::CwCI_2q1g2l_L(helcode_2q2l(pro),ind);
		case sub_leading_color: return CachedIntegral::CwCI_2q1g2l_SLC(helcode_2q2l(pro),ind);
		case nf: return CachedIntegral::CwCI_2q1g2l_nf(helcode_2q2l(pro),ind);
		case AX: return CachedIntegral::CwCI_2q1g2l_nf(helcode_2q2l(pro),ind);
		}
		break;
	}
	case 222: {
		switch (cs) {
		case leading_color: return CachedIntegral::CwCI_2q2g2l_L(helcode_2q2l(pro),ind);
		case sub_leading_color: return CachedIntegral::CwCI_2q2g2l_SLC(helcode_2q2l(pro),ind);
		case nf: return CachedIntegral::CwCI_2q2g2l_nf(helcode_2q2l(pro),ind);
		case nf_top: return CachedIntegral::CwCI_2q2g2l_nf_top(helcode_2q2l(pro),ind);
		case AX: return CachedIntegral::CwCI_2q2g2l_AX(helcode_2q2l(pro),ind);
		case AXSL: return CachedIntegral::CwCI_2q2g2l_AXSL(helcode_2q2l(pro),ind);
		case VECT: return CachedIntegral::CwCI_2q2g2l_VECT(helcode_2q2l(pro),ind);
		}
	break;
	}
	case 240:
		switch (cs) {
		case leading_color: return CwCI_2q2Q2l_L(helcode_2q2l2Q(pro),ind);
		case sub_leading_color: return CwCI_2q2Q2l_sl(helcode_2q2l2Q(pro),ind);
		case nf: return CwCI_2q2Q2l_nf(helcode_2q2l2Q(pro),ind);
		case nf_top: return CwCI_2q2Q2l_nf_top(helcode_2q2l2Q(pro),ind);
		case AX: return CwCI_2q2Q2l_AX(helcode_2q2l2Q(pro),ind);
		}
		break;
	case 100040:
		switch (cs) {
			case leading_color: return CwCI_2q2Q1y_L(helcode_2q2Q(pro),ind);
			case sub_leading_color: return CwCI_2q2Q1y_sl(helcode_2q2Q(pro),ind);
			case nf: return CwCI_2q2Q1y_nf(helcode_2q2Q(pro),ind);
		}
		break;
	case 100021: {
		switch (cs) {
		case leading_color: return CachedIntegral::CwCI_2q1g1y_L(helcode_2q1y(pro),ind);
		case sub_leading_color: return CachedIntegral::CwCI_2q1g1y_SLC(helcode_2q1y(pro),ind);
		case nf: return CachedIntegral::CwCI_2q1g1y_nf(helcode_2q1y(pro),ind);
		}
		break;
	}
	case 100022: {
		switch (cs) {
		case leading_color: return CachedIntegral::CwCI_2q2g1y_L(helcode_2q1y(pro),ind);
		case sub_leading_color: return CachedIntegral::CwCI_2q2g1y_SLC(helcode_2q1y(pro),ind);
		case nf: return CachedIntegral::CwCI_2q2g1y_nf(helcode_2q1y(pro),ind);
		}
		break;
	}
	case 2100020:{
	  return CwCI_Ptr(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1),cs,ind);
	}
	case 2000002: case 2000004: case 2000003:
	case 2000020: case 2000021: case 2000022: 
	case 2000040: case 2000060: {
		return CwCI_Ptr(replace_gluino_with_quark(pro),cs,ind);
    }
	case 2000220:
            switch (cs) {
                    case leading_color: return CachedIntegral::CwCI_2q2G2l_L(helcode_2q2l2G(pro),ind);
                    case sub_leading_color: return CachedIntegral::CwCI_2q2G2l_sl(helcode_2q2l2G(pro),ind);
                    case nf: return CachedIntegral::CwCI_2q2G2l_nf(helcode_2q2l2G(pro),ind);
//                    case nf_top: return CachedIntegral::CwCI_2q2G2l_nf_top(helcode_2q2l2G(pro),ind);
                    case AX: return CachedIntegral::CwCI_2q2G2l_AX(helcode_2q2l2G(pro),ind);
            }

            
            
            
            break;



			default: return 0;
	}
	return 0;
}


} /* CachedIntegral */

} /* BH */

