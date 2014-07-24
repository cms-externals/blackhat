/*
*C_2q2Q1g_wCI.cpp
*
* Created on 12/22, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q2Q1G_WCI_H_
#define C2Q2Q1G_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q2Q1g_LLT
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q1g_LRT
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q1g_nfLLT
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q2Q1G_WCI_H_ */
