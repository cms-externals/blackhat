/*
*C_4g_wCI.cpp
*
* Created on 11/8, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C4G_WCI_H_
#define C4G_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_4g_G
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_4g_nf
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C4G_WCI_H_ */
