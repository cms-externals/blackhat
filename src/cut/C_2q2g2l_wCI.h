/*
*C_2q2g2l_wCI.cpp
*
* Created on 9/25, 2009
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q2G2L_WCI_H_
#define C2Q2G2L_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q2g2l_AX
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g2l_AXSL
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g2l_L
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g2l_nf
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g2l_nf_top
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g2l_SLC
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g2l_VECT
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q2G2L_WCI_H_ */
