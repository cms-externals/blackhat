/*
*C_3g1ph_wCI.cpp
*
* Created on 12/9, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C3G1PH_WCI_H_
#define C3G1PH_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_3g1ph_G
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_3g1ph_nf
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C3G1PH_WCI_H_ */
