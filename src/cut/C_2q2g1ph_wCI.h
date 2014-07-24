/*
*C_2q2g1ph_wCI.cpp
*
* Created on 12/11, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q2G1PH_WCI_H_
#define C2Q2G1PH_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q2g1ph_LT
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g1ph_nfLT
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g1ph_RT
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q2G1PH_WCI_H_ */
