/*
*C_2q2Q1ph_wCI.cpp
*
* Created on 12/7, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q2Q1PH_WCI_H_
#define C2Q2Q1PH_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q2Q1ph_lc
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q1ph_nf
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q1ph_slc
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q2Q1PH_WCI_H_ */
