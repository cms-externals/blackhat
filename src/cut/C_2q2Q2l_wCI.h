/*
*C_2q2Q2l_wCI.cpp
*
* Created on 9/24, 2009
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q2Q2L_WCI_H_
#define C2Q2Q2L_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q2Q2l_AX
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q2l_L
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q2l_nf
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q2l_nf_top
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q2l_sl
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q2Q2L_WCI_H_ */
