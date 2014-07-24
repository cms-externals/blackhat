/*
 * C_2q2G2l_wCI.h
 *
 *  Created on: 24 Sep 2009
 *      Author: daniel
 */

#ifndef C_2Q2GL2L_WCI_H_
#define C_2Q2GL2L_WCI_H_


#include <vector>

namespace BH {

namespace CachedIntegral {

class Cut_Part_wCI;

Cut_Part_wCI* CwCI_2q2G2l_AX( int hc,const std::vector<int>& ind);

Cut_Part_wCI* CwCI_2q2G2l_L( int hc,const std::vector<int>& ind) ;

Cut_Part_wCI* CwCI_2q2G2l_nf( int hc,const std::vector<int>& ind);

Cut_Part_wCI* CwCI_2q2G2l_sl( int hc,const std::vector<int>& ind) ;

}


}



#endif /* C_2Q2GL2L_WCI_H_ */
