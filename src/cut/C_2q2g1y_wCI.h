/*
*C_2q2g1y_wCI.cpp
*
* Created on 2/5, 2011
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q2G1Y_WCI_H_
#define C2Q2G1Y_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q2g1y_L
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g1y_nf
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g1y_SLC
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q2G1Y_WCI_H_ */
