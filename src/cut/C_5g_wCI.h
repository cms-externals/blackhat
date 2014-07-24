/*
*C_5g_wCI.cpp
*
* Created on 9/24, 2009
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C5G_WCI_H_
#define C5G_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_5g_G
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_5g_nf
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C5G_WCI_H_ */
