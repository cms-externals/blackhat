/*
*C_6g_wCI.cpp
*
* Created on 11/8, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C6G_WCI_H_
#define C6G_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_6g_G
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_6g_nf
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C6G_WCI_H_ */
