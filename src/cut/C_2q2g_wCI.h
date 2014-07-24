/*
*C_2q2g_wCI.cpp
*
* Created on 11/12, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q2G_WCI_H_
#define C2Q2G_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q2g_LT
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2g_nfLT
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q2G_WCI_H_ */
