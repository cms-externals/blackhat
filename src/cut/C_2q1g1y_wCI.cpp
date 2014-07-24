/*
*C_2q1g1y_wCI.cpp
*
* Created on 9/24, 2009
*    Author: Zvi's doconvertlooptocpp_wCI.m script 
*/
 
#include <vector>
#include "integrals.h"
#include "cached_integral.h"
using namespace std;

using BH::CachedIntegral::Cached_Bubble_Integral_User;
using BH::CachedIntegral::Cached_Triangle_Integral_User;
using BH::CachedIntegral::Cached_Box_Integral_User;
 
namespace BH  {
 
class Index_Vector;
namespace CachedIntegral {
 
#define _VERBOSE 0
 
#define SPA(i,j) mc.spa(ind[i-1],ind[j-1])
#define SPB(i,j) mc.spb(ind[i-1],ind[j-1])
#define S(i,j) mc.s(ind[i-1],ind[j-1])
#define SS(i,j,k) mc.s(ind[i-1],ind[j-1],ind[k-1])

template<class T> static inline complex<T> square(complex<T> x) 
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x) 
{return(x*x*x);}
 


GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qppqmgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qppqmgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpmqmgap_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpmqmgam_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpgapqmp_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpgapqmm_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpgamqmp_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpgamqmm_L_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qppqmgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qppqmgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpmqmgap_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpmqmgam_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpgapqmp_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpgapqmm_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpgamqmp_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpgamqmm_nf_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpqmpgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpqmpgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpqmmgap_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpqmmgam_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpqmgapp_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpqmgapm_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpqmgamp_SLC_wCI)
GENERATE_CUT_PART_WCI_CLASS_DEFINITION(C2q1g1y_qpqmgamm_SLC_wCI)

 
C2q1g1y_qppqmgap_L_wCI::C2q1g1y_qppqmgap_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qppqmgap_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qppqmgap L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qppqmgam_L_wCI::C2q1g1y_qppqmgam_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qppqmgam_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qppqmgam L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpmqmgap_L_wCI::C2q1g1y_qpmqmgap_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpmqmgap_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpmqmgap L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpmqmgam_L_wCI::C2q1g1y_qpmqmgam_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpmqmgam_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpmqmgam L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpgapqmp_L_wCI::C2q1g1y_qpgapqmp_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpgapqmp_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpgapqmp L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpgapqmm_L_wCI::C2q1g1y_qpgapqmm_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpgapqmm_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpgapqmm L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpgamqmp_L_wCI::C2q1g1y_qpgamqmp_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpgamqmp_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpgamqmp L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpgamqmm_L_wCI::C2q1g1y_qpgamqmm_L_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpgamqmm_L_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpgamqmm L");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qppqmgap_nf_wCI::C2q1g1y_qppqmgap_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qppqmgap_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qppqmgap nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qppqmgam_nf_wCI::C2q1g1y_qppqmgam_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qppqmgam_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qppqmgam nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpmqmgap_nf_wCI::C2q1g1y_qpmqmgap_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpmqmgap_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpmqmgap nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpmqmgam_nf_wCI::C2q1g1y_qpmqmgam_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpmqmgam_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpmqmgam nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpgapqmp_nf_wCI::C2q1g1y_qpgapqmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpgapqmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpgapqmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpgapqmm_nf_wCI::C2q1g1y_qpgapqmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpgapqmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpgapqmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpgamqmp_nf_wCI::C2q1g1y_qpgamqmp_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpgamqmp_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpgamqmp nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpgamqmm_nf_wCI::C2q1g1y_qpgamqmm_nf_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpgamqmm_nf_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpgamqmm nf");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpqmpgap_SLC_wCI::C2q1g1y_qpqmpgap_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpqmpgap_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpqmpgap SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpqmpgam_SLC_wCI::C2q1g1y_qpqmpgam_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpqmpgam_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpqmpgam SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpqmmgap_SLC_wCI::C2q1g1y_qpqmmgap_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpqmmgap_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpqmmgap SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpqmmgam_SLC_wCI::C2q1g1y_qpqmmgam_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpqmmgam_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpqmmgam SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpqmgapp_SLC_wCI::C2q1g1y_qpqmgapp_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpqmgapp_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpqmgapp SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpqmgapm_SLC_wCI::C2q1g1y_qpqmgapm_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpqmgapm_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpqmgapm SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpqmgamp_SLC_wCI::C2q1g1y_qpqmgamp_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpqmgamp_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpqmgamp SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 


C2q1g1y_qpqmgamm_SLC_wCI::C2q1g1y_qpqmgamm_SLC_wCI
        (const std::vector<int>& ind){
}


template <class T> SeriesC<T> 
       C2q1g1y_qpqmgamm_SLC_wCI::eval_fn
     (momentum_configuration<T>& mc,const Index_Vector& ind,int mu){
#if _VERBOSE
  _MESSAGE("C2q1g1y :  qpqmgamm SLC");
#endif
    SeriesC<T> res(-2,0);
    return res;
  
} 
 
 
 // *************** table of switch values ************* 
 
#define _C_qppqmgap_L C2q1g1y_1136_L
#define _C_qppqmgam_L C2q1g1y_920_L
#define _C_qpmqmgap_L C2q1g1y_1118_L
#define _C_qpmqmgam_L C2q1g1y_902_L
#define _C_qpgapqmp_L C2q1g1y_716_L
#define _C_qpgapqmm_L C2q1g1y_68_L
#define _C_qpgamqmp_L C2q1g1y_710_L
#define _C_qpgamqmm_L C2q1g1y_62_L
#define _C_qppqmgap_nf C2q1g1y_1136_nf
#define _C_qppqmgam_nf C2q1g1y_920_nf
#define _C_qpmqmgap_nf C2q1g1y_1118_nf
#define _C_qpmqmgam_nf C2q1g1y_902_nf
#define _C_qpgapqmp_nf C2q1g1y_716_nf
#define _C_qpgapqmm_nf C2q1g1y_68_nf
#define _C_qpgamqmp_nf C2q1g1y_710_nf
#define _C_qpgamqmm_nf C2q1g1y_62_nf
#define _C_qpqmpgap_SLC C2q1g1y_1196_SLC
#define _C_qpqmpgam_SLC C2q1g1y_980_SLC
#define _C_qpqmmgap_SLC C2q1g1y_1088_SLC
#define _C_qpqmmgam_SLC C2q1g1y_872_SLC
#define _C_qpqmgapp_SLC C2q1g1y_836_SLC
#define _C_qpqmgapm_SLC C2q1g1y_188_SLC
#define _C_qpqmgamp_SLC C2q1g1y_800_SLC
#define _C_qpqmgamm_SLC C2q1g1y_152_SLC
 
 
 // *************** more macro definitions ************* 
 
#define _CASE_qppqmgap_L case 1136
 
#define _CASE_qppqmgam_L case 920
 
#define _CASE_qpmqmgap_L case 1118
 
#define _CASE_qpmqmgam_L case 902
 
#define _CASE_qpgapqmp_L case 716
 
#define _CASE_qpgapqmm_L case 68
 
#define _CASE_qpgamqmp_L case 710
 
#define _CASE_qpgamqmm_L case 62
 
#define _CASE_qppqmgap_nf case 1136
 
#define _CASE_qppqmgam_nf case 920
 
#define _CASE_qpmqmgap_nf case 1118
 
#define _CASE_qpmqmgam_nf case 902
 
#define _CASE_qpgapqmp_nf case 716
 
#define _CASE_qpgapqmm_nf case 68
 
#define _CASE_qpgamqmp_nf case 710
 
#define _CASE_qpgamqmm_nf case 62
 
#define _CASE_qpqmpgap_SLC case 1196
 
#define _CASE_qpqmpgam_SLC case 980
 
#define _CASE_qpqmmgap_SLC case 1088
 
#define _CASE_qpqmmgam_SLC case 872
 
#define _CASE_qpqmgapp_SLC case 836
 
#define _CASE_qpqmgapm_SLC case 188
 
#define _CASE_qpqmgamp_SLC case 800
 
#define _CASE_qpqmgamm_SLC case 152
 
 
 // *************** define pointers ************* 
 
Cut_Part_wCI* CwCI_2q1g1y_L( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qppqmgap_L: return new 
                       C2q1g1y_qppqmgap_L_wCI(ind);
    _CASE_qppqmgam_L: return new 
                       C2q1g1y_qppqmgam_L_wCI(ind);
    _CASE_qpmqmgap_L: return new 
                       C2q1g1y_qpmqmgap_L_wCI(ind);
    _CASE_qpmqmgam_L: return new 
                       C2q1g1y_qpmqmgam_L_wCI(ind);
    _CASE_qpgapqmp_L: return new 
                       C2q1g1y_qpgapqmp_L_wCI(ind);
    _CASE_qpgapqmm_L: return new 
                       C2q1g1y_qpgapqmm_L_wCI(ind);
    _CASE_qpgamqmp_L: return new 
                       C2q1g1y_qpgamqmp_L_wCI(ind);
    _CASE_qpgamqmm_L: return new 
                       C2q1g1y_qpgamqmm_L_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q1g1y_nf( int hc,const std::vector<int>& ind) { 
    switch (hc) {
    _CASE_qppqmgap_nf: return new 
                       C2q1g1y_qppqmgap_nf_wCI(ind);
    _CASE_qppqmgam_nf: return new 
                       C2q1g1y_qppqmgam_nf_wCI(ind);
    _CASE_qpmqmgap_nf: return new 
                       C2q1g1y_qpmqmgap_nf_wCI(ind);
    _CASE_qpmqmgam_nf: return new 
                       C2q1g1y_qpmqmgam_nf_wCI(ind);
    _CASE_qpgapqmp_nf: return new 
                       C2q1g1y_qpgapqmp_nf_wCI(ind);
    _CASE_qpgapqmm_nf: return new 
                       C2q1g1y_qpgapqmm_nf_wCI(ind);
    _CASE_qpgamqmp_nf: return new 
                       C2q1g1y_qpgamqmp_nf_wCI(ind);
    _CASE_qpgamqmm_nf: return new 
                       C2q1g1y_qpgamqmm_nf_wCI(ind);
 
       default: return 0;
                   }
      }
 
Cut_Part_wCI* CwCI_2q1g1y_SLC( int hc,const std::vector<int>& ind) { 

    switch (hc) {
    _CASE_qpqmpgap_SLC: return new 
                       C2q1g1y_qpqmpgap_SLC_wCI(ind);
    _CASE_qpqmpgam_SLC: return new 
                       C2q1g1y_qpqmpgam_SLC_wCI(ind);
    _CASE_qpqmmgap_SLC: return new 
                       C2q1g1y_qpqmmgap_SLC_wCI(ind);
    _CASE_qpqmmgam_SLC: return new 
                       C2q1g1y_qpqmmgam_SLC_wCI(ind);
    _CASE_qpqmgapp_SLC: return new 
                       C2q1g1y_qpqmgapp_SLC_wCI(ind);
    _CASE_qpqmgapm_SLC: return new 
                       C2q1g1y_qpqmgapm_SLC_wCI(ind);
    _CASE_qpqmgamp_SLC: return new 
                       C2q1g1y_qpqmgamp_SLC_wCI(ind);
    _CASE_qpqmgamm_SLC: return new 
                       C2q1g1y_qpqmgamm_SLC_wCI(ind);
 
       default: return 0;
                   }
      }
 
 
 }
 }
