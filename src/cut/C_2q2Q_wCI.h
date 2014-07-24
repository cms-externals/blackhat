/*
*C_2q2Q_wCI.cpp
*
* Created on 11/15, 2010
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#ifndef C2Q2Q_WCI_H_
#define C2Q2Q_WCI_H_
 
#include <vector> 
 
namespace BH {
 
namespace CachedIntegral {
 
class Cut_Part_wCI;
 
Cut_Part_wCI* CwCI_2q2Q_LLT
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q_LRT
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q_nfLLT
      (int hc,const std::vector<int>& ind); 
Cut_Part_wCI* CwCI_2q2Q_nfLRT
      (int hc,const std::vector<int>& ind); 
 
}
 
 
}
 
 
#endif /* C2Q2Q_WCI_H_ */
