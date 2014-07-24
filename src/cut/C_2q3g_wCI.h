/*
*C_2q3g_wCI.cpp
*
* Created on 9/24, 2009
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q3G_WCI_H_
#define C2Q3G_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q3g_L
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q3g_nf
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q3G_WCI_H_ */
