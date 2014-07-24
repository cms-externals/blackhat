/*
 * A0_2q3g2l_eval.cpp
 *
 * created with Harald's scripts from trees computed by Lance
 * notes:
 * *) trees are manifest dual super conformal expressions
 * *) speed-up from idendifying commone sumbexpr. still possible
 *  
 *
 */

/* Implementation of super dual conformal trees  */

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;

#define SPA(i,j) ep.spa(i,j)
#define SPB(i,j) ep.spb(i,j)

/* 
 * more spab & spaa  macros 
*/
#include "A0_2q3g2l_eval_fwd.hpp" 

namespace BH {

template<class T> static inline complex<T> square(complex<T> x)
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x)
{return(x*x*x);}

//

/*
 *
 *
 * The 2 quarks 3g 2 lepton amplitudes
 *
 *
 */

template <class T> complex<T> (*A2q3g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {
/*
 * some cases commented in order to match ordinary tree-amplitudes
*/
//case 116705:	 return &A2q3g2l_qmmmpqbplplbm_eval<2,3,4,5,6,0,1>;//lp lbm qm m m p qbp
case 116710:	 return &A2q3g2l_qmmmpqbplmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qm m m p qbp
//case 117353:	 return &A2q3g2l_qmpmpqbplplbm_eval<2,3,4,5,6,0,1>;//lp lbm qm p m p qbp
case 117358:	 return &A2q3g2l_qmpmpqbplmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qm p m p qbp
//case 120593:	 return &A2q3g2l_qmmppqbplplbm_eval<2,3,4,5,6,0,1>;//lp lbm qm m p p qbp
case 120598:	 return &A2q3g2l_qmmppqbplmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qm m p p qbp
//case 121241:	 return &A2q3g2l_qmpppqbplplbm_eval<2,3,4,5,6,0,1>;//lp lbm qm p p p qbp
case 121246:	 return &A2q3g2l_qmpppqbplmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qm p p p qbp
//case 131112:	 return &A2q3g2l_qpmmmqbmlplbm_eval<6,0,1,2,3,4,5>;//m m m qbm lp lbm qp
//case 131115:	 return &A2q3g2l_qppmmqbmlplbm_eval<6,0,1,2,3,4,5>;//p m m qbm lp lbm qp
//case 131130:	 return &A2q3g2l_qpmpmqbmlplbm_eval<6,0,1,2,3,4,5>;//m p m qbm lp lbm qp
//case 131133:	 return &A2q3g2l_qpppmqbmlplbm_eval<6,0,1,2,3,4,5>;//p p m qbm lp lbm qp
//case 131220:	 return &A2q3g2l_qpmmpqbmlplbm_eval<6,0,1,2,3,4,5>;//m m p qbm lp lbm qp
//case 131223:	 return &A2q3g2l_qppmpqbmlplbm_eval<6,0,1,2,3,4,5>;//p m p qbm lp lbm qp
//case 131238:	 return &A2q3g2l_qpmppqbmlplbm_eval<6,0,1,2,3,4,5>;//m p p qbm lp lbm qp
//case 131241:	 return &A2q3g2l_qppppqbmlplbm_eval<6,0,1,2,3,4,5>;//p p p qbm lp lbm qp
case 137592:	 return &A2q3g2l_qpmmmqbmlmlbp_eval<6,0,1,2,3,4,5>;//m m m qbm lm lbp qp
case 137595:	 return &A2q3g2l_qppmmqbmlmlbp_eval<6,0,1,2,3,4,5>;//p m m qbm lm lbp qp
case 137610:	 return &A2q3g2l_qpmpmqbmlmlbp_eval<6,0,1,2,3,4,5>;//m p m qbm lm lbp qp
case 137613:	 return &A2q3g2l_qpppmqbmlmlbp_eval<6,0,1,2,3,4,5>;//p p m qbm lm lbp qp
case 137700:	 return &A2q3g2l_qpmmpqbmlmlbp_eval<6,0,1,2,3,4,5>;//m m p qbm lm lbp qp
case 137703:	 return &A2q3g2l_qppmpqbmlmlbp_eval<6,0,1,2,3,4,5>;//p m p qbm lm lbp qp
case 137718:	 return &A2q3g2l_qpmppqbmlmlbp_eval<6,0,1,2,3,4,5>;//m p p qbm lm lbp qp
case 137721:	 return &A2q3g2l_qppppqbmlmlbp_eval<6,0,1,2,3,4,5>;//p p p qbm lm lbp qp
//case 140360:	 return &A2q3g2l_qmmmpqbplplbm_eval<3,4,5,6,0,1,2>;//qbp lp lbm qm m m p
case 140390:	 return &A2q3g2l_qmmmpqbplmlbp_eval<3,4,5,6,0,1,2>;//qbp lm lbp qm m m p
//case 140575:	 return &A2q3g2l_qpmmpqbmlplbm_eval<3,4,5,6,0,1,2>;//qbm lp lbm qp m m p
case 140605:	 return &A2q3g2l_qpmmpqbmlmlbp_eval<3,4,5,6,0,1,2>;//qbm lm lbp qp m m p
//case 14112:	 return &A2q3g2l_qmmmmqbplplbm_eval<5,6,0,1,2,3,4>;//m m qbp lp lbm qm m
//case 14115:	 return &A2q3g2l_qmmpmqbplplbm_eval<5,6,0,1,2,3,4>;//p m qbp lp lbm qm m
//case 14130:	 return &A2q3g2l_qmmmpqbplplbm_eval<5,6,0,1,2,3,4>;//m p qbp lp lbm qm m
//case 14133:	 return &A2q3g2l_qmmppqbplplbm_eval<5,6,0,1,2,3,4>;//p p qbp lp lbm qm m
//case 142320:	 return &A2q3g2l_qmmpmqbplplbm_eval<4,5,6,0,1,2,3>;//m qbp lp lbm qm m p
//case 142323:	 return &A2q3g2l_qmmppqbplplbm_eval<4,5,6,0,1,2,3>;//p qbp lp lbm qm m p
case 142500:	 return &A2q3g2l_qmmpmqbplmlbp_eval<4,5,6,0,1,2,3>;//m qbp lm lbp qm m p
case 142503:	 return &A2q3g2l_qmmppqbplmlbp_eval<4,5,6,0,1,2,3>;//p qbp lm lbp qm m p
//case 143610:	 return &A2q3g2l_qpmpmqbmlplbm_eval<4,5,6,0,1,2,3>;//m qbm lp lbm qp m p
//case 143613:	 return &A2q3g2l_qpmppqbmlplbm_eval<4,5,6,0,1,2,3>;//p qbm lp lbm qp m p
case 143790:	 return &A2q3g2l_qpmpmqbmlmlbp_eval<4,5,6,0,1,2,3>;//m qbm lm lbp qp m p
case 143793:	 return &A2q3g2l_qpmppqbmlmlbp_eval<4,5,6,0,1,2,3>;//p qbm lm lbp qp m p
//case 144248:	 return &A2q3g2l_qmpmpqbplplbm_eval<3,4,5,6,0,1,2>;//qbp lp lbm qm p m p
case 144278:	 return &A2q3g2l_qmpmpqbplmlbp_eval<3,4,5,6,0,1,2>;//qbp lm lbp qm p m p
//case 144463:	 return &A2q3g2l_qppmpqbmlplbm_eval<3,4,5,6,0,1,2>;//qbm lp lbm qp p m p
case 144493:	 return &A2q3g2l_qppmpqbmlmlbp_eval<3,4,5,6,0,1,2>;//qbm lm lbp qp p m p
case 15192:	 return &A2q3g2l_qmmmmqbplmlbp_eval<5,6,0,1,2,3,4>;//m m qbp lm lbp qm m
case 15195:	 return &A2q3g2l_qmmpmqbplmlbp_eval<5,6,0,1,2,3,4>;//p m qbp lm lbp qm m
case 15210:	 return &A2q3g2l_qmmmpqbplmlbp_eval<5,6,0,1,2,3,4>;//m p qbp lm lbp qm m
case 15213:	 return &A2q3g2l_qmmppqbplmlbp_eval<5,6,0,1,2,3,4>;//p p qbp lm lbp qm m
//case 154080:	 return &A2q3g2l_qmpmmqbplplbm_eval<5,6,0,1,2,3,4>;//m m qbp lp lbm qm p
//case 154083:	 return &A2q3g2l_qmppmqbplplbm_eval<5,6,0,1,2,3,4>;//p m qbp lp lbm qm p
//case 154098:	 return &A2q3g2l_qmpmpqbplplbm_eval<5,6,0,1,2,3,4>;//m p qbp lp lbm qm p
//case 154101:	 return &A2q3g2l_qmpppqbplplbm_eval<5,6,0,1,2,3,4>;//p p qbp lp lbm qm p
case 155160:	 return &A2q3g2l_qmpmmqbplmlbp_eval<5,6,0,1,2,3,4>;//m m qbp lm lbp qm p
case 155163:	 return &A2q3g2l_qmppmqbplmlbp_eval<5,6,0,1,2,3,4>;//p m qbp lm lbp qm p
case 155178:	 return &A2q3g2l_qmpmpqbplmlbp_eval<5,6,0,1,2,3,4>;//m p qbp lm lbp qm p
case 155181:	 return &A2q3g2l_qmpppqbplmlbp_eval<5,6,0,1,2,3,4>;//p p qbp lm lbp qm p
//case 161820:	 return &A2q3g2l_qppmmqbmlplbm_eval<5,6,0,1,2,3,4>;//m m qbm lp lbm qp p
//case 161823:	 return &A2q3g2l_qpppmqbmlplbm_eval<5,6,0,1,2,3,4>;//p m qbm lp lbm qp p
//case 161838:	 return &A2q3g2l_qppmpqbmlplbm_eval<5,6,0,1,2,3,4>;//m p qbm lp lbm qp p
//case 161841:	 return &A2q3g2l_qppppqbmlplbm_eval<5,6,0,1,2,3,4>;//p p qbm lp lbm qp p
case 162900:	 return &A2q3g2l_qppmmqbmlmlbp_eval<5,6,0,1,2,3,4>;//m m qbm lm lbp qp p
case 162903:	 return &A2q3g2l_qpppmqbmlmlbp_eval<5,6,0,1,2,3,4>;//p m qbm lm lbp qp p
case 162918:	 return &A2q3g2l_qppmpqbmlmlbp_eval<5,6,0,1,2,3,4>;//m p qbm lm lbp qp p
case 162921:	 return &A2q3g2l_qppppqbmlmlbp_eval<5,6,0,1,2,3,4>;//p p qbm lm lbp qp p
//case 163688:	 return &A2q3g2l_qmmppqbplplbm_eval<3,4,5,6,0,1,2>;//qbp lp lbm qm m p p
case 163718:	 return &A2q3g2l_qmmppqbplmlbp_eval<3,4,5,6,0,1,2>;//qbp lm lbp qm m p p
//case 163903:	 return &A2q3g2l_qpmppqbmlplbm_eval<3,4,5,6,0,1,2>;//qbm lp lbm qp m p p
case 163933:	 return &A2q3g2l_qpmppqbmlmlbp_eval<3,4,5,6,0,1,2>;//qbm lm lbp qp m p p
//case 165648:	 return &A2q3g2l_qmppmqbplplbm_eval<4,5,6,0,1,2,3>;//m qbp lp lbm qm p p
//case 165651:	 return &A2q3g2l_qmpppqbplplbm_eval<4,5,6,0,1,2,3>;//p qbp lp lbm qm p p
case 165828:	 return &A2q3g2l_qmppmqbplmlbp_eval<4,5,6,0,1,2,3>;//m qbp lm lbp qm p p
case 165831:	 return &A2q3g2l_qmpppqbplmlbp_eval<4,5,6,0,1,2,3>;//p qbp lm lbp qm p p
//case 166938:	 return &A2q3g2l_qpppmqbmlplbm_eval<4,5,6,0,1,2,3>;//m qbm lp lbm qp p p
//case 166941:	 return &A2q3g2l_qppppqbmlplbm_eval<4,5,6,0,1,2,3>;//p qbm lp lbm qp p p
case 167118:	 return &A2q3g2l_qpppmqbmlmlbp_eval<4,5,6,0,1,2,3>;//m qbm lm lbp qp p p
case 167121:	 return &A2q3g2l_qppppqbmlmlbp_eval<4,5,6,0,1,2,3>;//p qbm lm lbp qp p p
//case 167576:	 return &A2q3g2l_qmpppqbplplbm_eval<3,4,5,6,0,1,2>;//qbp lp lbm qm p p p
case 167606:	 return &A2q3g2l_qmpppqbplmlbp_eval<3,4,5,6,0,1,2>;//qbp lm lbp qm p p p
//case 167791:	 return &A2q3g2l_qppppqbmlplbm_eval<3,4,5,6,0,1,2>;//qbm lp lbm qp p p p
case 167821:	 return &A2q3g2l_qppppqbmlmlbp_eval<3,4,5,6,0,1,2>;//qbm lm lbp qp p p p
case 194417:	 return &A2q3g2l_qpmmmqbmlmlbp_eval<1,2,3,4,5,6,0>;//lbp qp m m m qbm lm
case 194525:	 return &A2q3g2l_qppmmqbmlmlbp_eval<1,2,3,4,5,6,0>;//lbp qp p m m qbm lm
case 195065:	 return &A2q3g2l_qpmpmqbmlmlbp_eval<1,2,3,4,5,6,0>;//lbp qp m p m qbm lm
case 195173:	 return &A2q3g2l_qpppmqbmlmlbp_eval<1,2,3,4,5,6,0>;//lbp qp p p m qbm lm
case 198305:	 return &A2q3g2l_qpmmpqbmlmlbp_eval<1,2,3,4,5,6,0>;//lbp qp m m p qbm lm
case 198413:	 return &A2q3g2l_qppmpqbmlmlbp_eval<1,2,3,4,5,6,0>;//lbp qp p m p qbm lm
case 198953:	 return &A2q3g2l_qpmppqbmlmlbp_eval<1,2,3,4,5,6,0>;//lbp qp m p p qbm lm
case 199061:	 return &A2q3g2l_qppppqbmlmlbp_eval<1,2,3,4,5,6,0>;//lbp qp p p p qbm lm
case 202187:	 return &A2q3g2l_qmmmmqbplmlbp_eval<1,2,3,4,5,6,0>;//lbp qm m m m qbp lm
case 202295:	 return &A2q3g2l_qmpmmqbplmlbp_eval<1,2,3,4,5,6,0>;//lbp qm p m m qbp lm
case 202835:	 return &A2q3g2l_qmmpmqbplmlbp_eval<1,2,3,4,5,6,0>;//lbp qm m p m qbp lm
case 202943:	 return &A2q3g2l_qmppmqbplmlbp_eval<1,2,3,4,5,6,0>;//lbp qm p p m qbp lm
case 206075:	 return &A2q3g2l_qmmmpqbplmlbp_eval<1,2,3,4,5,6,0>;//lbp qm m m p qbp lm
case 206183:	 return &A2q3g2l_qmpmpqbplmlbp_eval<1,2,3,4,5,6,0>;//lbp qm p m p qbp lm
case 206723:	 return &A2q3g2l_qmmppqbplmlbp_eval<1,2,3,4,5,6,0>;//lbp qm m p p qbp lm
case 206831:	 return &A2q3g2l_qmpppqbplmlbp_eval<1,2,3,4,5,6,0>;//lbp qm p p p qbp lm
//case 21852:	 return &A2q3g2l_qpmmmqbmlplbm_eval<5,6,0,1,2,3,4>;//m m qbm lp lbm qp m
//case 21855:	 return &A2q3g2l_qpmpmqbmlplbm_eval<5,6,0,1,2,3,4>;//p m qbm lp lbm qp m
//case 21870:	 return &A2q3g2l_qpmmpqbmlplbm_eval<5,6,0,1,2,3,4>;//m p qbm lp lbm qp m
//case 21873:	 return &A2q3g2l_qpmppqbmlplbm_eval<5,6,0,1,2,3,4>;//p p qbm lp lbm qp m
//case 226802:	 return &A2q3g2l_qpmmmqbmlplbm_eval<0,1,2,3,4,5,6>;//qp m m m qbm lp lbm
//case 226820:	 return &A2q3g2l_qppmmqbmlplbm_eval<0,1,2,3,4,5,6>;//qp p m m qbm lp lbm
//case 226910:	 return &A2q3g2l_qpmpmqbmlplbm_eval<0,1,2,3,4,5,6>;//qp m p m qbm lp lbm
//case 226928:	 return &A2q3g2l_qpppmqbmlplbm_eval<0,1,2,3,4,5,6>;//qp p p m qbm lp lbm
//case 227450:	 return &A2q3g2l_qpmmpqbmlplbm_eval<0,1,2,3,4,5,6>;//qp m m p qbm lp lbm
//case 227468:	 return &A2q3g2l_qppmpqbmlplbm_eval<0,1,2,3,4,5,6>;//qp p m p qbm lp lbm
//case 227558:	 return &A2q3g2l_qpmppqbmlplbm_eval<0,1,2,3,4,5,6>;//qp m p p qbm lp lbm
//case 227576:	 return &A2q3g2l_qppppqbmlplbm_eval<0,1,2,3,4,5,6>;//qp p p p qbm lp lbm
//case 228097:	 return &A2q3g2l_qmmmmqbplplbm_eval<0,1,2,3,4,5,6>;//qm m m m qbp lp lbm
//case 228115:	 return &A2q3g2l_qmpmmqbplplbm_eval<0,1,2,3,4,5,6>;//qm p m m qbp lp lbm
//case 228205:	 return &A2q3g2l_qmmpmqbplplbm_eval<0,1,2,3,4,5,6>;//qm m p m qbp lp lbm
//case 228223:	 return &A2q3g2l_qmppmqbplplbm_eval<0,1,2,3,4,5,6>;//qm p p m qbp lp lbm
//case 228745:	 return &A2q3g2l_qmmmpqbplplbm_eval<0,1,2,3,4,5,6>;//qm m m p qbp lp lbm
//case 228763:	 return &A2q3g2l_qmpmpqbplplbm_eval<0,1,2,3,4,5,6>;//qm p m p qbp lp lbm
//case 228853:	 return &A2q3g2l_qmmppqbplplbm_eval<0,1,2,3,4,5,6>;//qm m p p qbp lp lbm
//case 228871:	 return &A2q3g2l_qmpppqbplplbm_eval<0,1,2,3,4,5,6>;//qm p p p qbp lp lbm
case 22932:	 return &A2q3g2l_qpmmmqbmlmlbp_eval<5,6,0,1,2,3,4>;//m m qbm lm lbp qp m
case 22935:	 return &A2q3g2l_qpmpmqbmlmlbp_eval<5,6,0,1,2,3,4>;//p m qbm lm lbp qp m
case 22950:	 return &A2q3g2l_qpmmpqbmlmlbp_eval<5,6,0,1,2,3,4>;//m p qbm lm lbp qp m
case 22953:	 return &A2q3g2l_qpmppqbmlmlbp_eval<5,6,0,1,2,3,4>;//p p qbm lm lbp qp m
//case 2352:	 return &A2q3g2l_qmmmmqbplplbm_eval<4,5,6,0,1,2,3>;//m qbp lp lbm qm m m
//case 2355:	 return &A2q3g2l_qmmmpqbplplbm_eval<4,5,6,0,1,2,3>;//p qbp lp lbm qm m m
//case 23720:	 return &A2q3g2l_qmmpmqbplplbm_eval<3,4,5,6,0,1,2>;//qbp lp lbm qm m p m
case 23750:	 return &A2q3g2l_qmmpmqbplmlbp_eval<3,4,5,6,0,1,2>;//qbp lm lbp qm m p m
//case 23935:	 return &A2q3g2l_qpmpmqbmlplbm_eval<3,4,5,6,0,1,2>;//qbm lp lbm qp m p m
case 23965:	 return &A2q3g2l_qpmpmqbmlmlbp_eval<3,4,5,6,0,1,2>;//qbm lm lbp qp m p m
//case 241072:	 return &A2q3g2l_qpmmmqbmlplbm_eval<1,2,3,4,5,6,0>;//lbm qp m m m qbm lp
//case 241180:	 return &A2q3g2l_qppmmqbmlplbm_eval<1,2,3,4,5,6,0>;//lbm qp p m m qbm lp
//case 241720:	 return &A2q3g2l_qpmpmqbmlplbm_eval<1,2,3,4,5,6,0>;//lbm qp m p m qbm lp
//case 241828:	 return &A2q3g2l_qpppmqbmlplbm_eval<1,2,3,4,5,6,0>;//lbm qp p p m qbm lp
//case 244960:	 return &A2q3g2l_qpmmpqbmlplbm_eval<1,2,3,4,5,6,0>;//lbm qp m m p qbm lp
//case 245068:	 return &A2q3g2l_qppmpqbmlplbm_eval<1,2,3,4,5,6,0>;//lbm qp p m p qbm lp
//case 245608:	 return &A2q3g2l_qpmppqbmlplbm_eval<1,2,3,4,5,6,0>;//lbm qp m p p qbm lp
//case 245716:	 return &A2q3g2l_qppppqbmlplbm_eval<1,2,3,4,5,6,0>;//lbm qp p p p qbm lp
//case 248842:	 return &A2q3g2l_qmmmmqbplplbm_eval<1,2,3,4,5,6,0>;//lbm qm m m m qbp lp
//case 248950:	 return &A2q3g2l_qmpmmqbplplbm_eval<1,2,3,4,5,6,0>;//lbm qm p m m qbp lp
//case 249490:	 return &A2q3g2l_qmmpmqbplplbm_eval<1,2,3,4,5,6,0>;//lbm qm m p m qbp lp
//case 249598:	 return &A2q3g2l_qmppmqbplplbm_eval<1,2,3,4,5,6,0>;//lbm qm p p m qbp lp
//case 252730:	 return &A2q3g2l_qmmmpqbplplbm_eval<1,2,3,4,5,6,0>;//lbm qm m m p qbp lp
//case 252838:	 return &A2q3g2l_qmpmpqbplplbm_eval<1,2,3,4,5,6,0>;//lbm qm p m p qbp lp
case 2532:	 return &A2q3g2l_qmmmmqbplmlbp_eval<4,5,6,0,1,2,3>;//m qbp lm lbp qm m m
//case 253378:	 return &A2q3g2l_qmmppqbplplbm_eval<1,2,3,4,5,6,0>;//lbm qm m p p qbp lp
//case 253486:	 return &A2q3g2l_qmpppqbplplbm_eval<1,2,3,4,5,6,0>;//lbm qm p p p qbp lp
case 2535:	 return &A2q3g2l_qmmmpqbplmlbp_eval<4,5,6,0,1,2,3>;//p qbp lm lbp qm m m
//case 25680:	 return &A2q3g2l_qmpmmqbplplbm_eval<4,5,6,0,1,2,3>;//m qbp lp lbm qm p m
//case 25683:	 return &A2q3g2l_qmpmpqbplplbm_eval<4,5,6,0,1,2,3>;//p qbp lp lbm qm p m
case 25860:	 return &A2q3g2l_qmpmmqbplmlbp_eval<4,5,6,0,1,2,3>;//m qbp lm lbp qm p m
case 25863:	 return &A2q3g2l_qmpmpqbplmlbp_eval<4,5,6,0,1,2,3>;//p qbp lm lbp qm p m
case 265682:	 return &A2q3g2l_qpmmmqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp m m m qbm lm lbp
case 265700:	 return &A2q3g2l_qppmmqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp p m m qbm lm lbp
case 265790:	 return &A2q3g2l_qpmpmqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp m p m qbm lm lbp
case 265808:	 return &A2q3g2l_qpppmqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp p p m qbm lm lbp
case 266330:	 return &A2q3g2l_qpmmpqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp m m p qbm lm lbp
case 266348:	 return &A2q3g2l_qppmpqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp p m p qbm lm lbp
case 266438:	 return &A2q3g2l_qpmppqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp m p p qbm lm lbp
case 266456:	 return &A2q3g2l_qppppqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp p p p qbm lm lbp
case 266977:	 return &A2q3g2l_qmmmmqbplmlbp_eval<0,1,2,3,4,5,6>;//qm m m m qbp lm lbp
case 266995:	 return &A2q3g2l_qmpmmqbplmlbp_eval<0,1,2,3,4,5,6>;//qm p m m qbp lm lbp
case 267085:	 return &A2q3g2l_qmmpmqbplmlbp_eval<0,1,2,3,4,5,6>;//qm m p m qbp lm lbp
case 267103:	 return &A2q3g2l_qmppmqbplmlbp_eval<0,1,2,3,4,5,6>;//qm p p m qbp lm lbp
case 267625:	 return &A2q3g2l_qmmmpqbplmlbp_eval<0,1,2,3,4,5,6>;//qm m m p qbp lm lbp
case 267643:	 return &A2q3g2l_qmpmpqbplmlbp_eval<0,1,2,3,4,5,6>;//qm p m p qbp lm lbp
case 267733:	 return &A2q3g2l_qmmppqbplmlbp_eval<0,1,2,3,4,5,6>;//qm m p p qbp lm lbp
case 267751:	 return &A2q3g2l_qmpppqbplmlbp_eval<0,1,2,3,4,5,6>;//qm p p p qbp lm lbp
//case 26970:	 return &A2q3g2l_qppmmqbmlplbm_eval<4,5,6,0,1,2,3>;//m qbm lp lbm qp p m
//case 26973:	 return &A2q3g2l_qppmpqbmlplbm_eval<4,5,6,0,1,2,3>;//p qbm lp lbm qp p m
case 27150:	 return &A2q3g2l_qppmmqbmlmlbp_eval<4,5,6,0,1,2,3>;//m qbm lm lbp qp p m
case 27153:	 return &A2q3g2l_qppmpqbmlmlbp_eval<4,5,6,0,1,2,3>;//p qbm lm lbp qp p m
//case 27608:	 return &A2q3g2l_qmppmqbplplbm_eval<3,4,5,6,0,1,2>;//qbp lp lbm qm p p m
case 27638:	 return &A2q3g2l_qmppmqbplmlbp_eval<3,4,5,6,0,1,2>;//qbp lm lbp qm p p m
//case 27823:	 return &A2q3g2l_qpppmqbmlplbm_eval<3,4,5,6,0,1,2>;//qbm lp lbm qp p p m
case 27853:	 return &A2q3g2l_qpppmqbmlmlbp_eval<3,4,5,6,0,1,2>;//qbm lm lbp qp p p m
//case 3642:	 return &A2q3g2l_qpmmmqbmlplbm_eval<4,5,6,0,1,2,3>;//m qbm lp lbm qp m m
//case 3645:	 return &A2q3g2l_qpmmpqbmlplbm_eval<4,5,6,0,1,2,3>;//p qbm lp lbm qp m m
case 3822:	 return &A2q3g2l_qpmmmqbmlmlbp_eval<4,5,6,0,1,2,3>;//m qbm lm lbp qp m m
case 3825:	 return &A2q3g2l_qpmmpqbmlmlbp_eval<4,5,6,0,1,2,3>;//p qbm lm lbp qp m m
//case 392:	 return &A2q3g2l_qmmmmqbplplbm_eval<3,4,5,6,0,1,2>;//qbp lp lbm qm m m m
case 422:	 return &A2q3g2l_qmmmmqbplmlbp_eval<3,4,5,6,0,1,2>;//qbp lm lbp qm m m m
//case 4280:	 return &A2q3g2l_qmpmmqbplplbm_eval<3,4,5,6,0,1,2>;//qbp lp lbm qm p m m
case 4310:	 return &A2q3g2l_qmpmmqbplmlbp_eval<3,4,5,6,0,1,2>;//qbp lm lbp qm p m m
//case 4495:	 return &A2q3g2l_qppmmqbmlplbm_eval<3,4,5,6,0,1,2>;//qbm lp lbm qp p m m
case 4525:	 return &A2q3g2l_qppmmqbmlmlbp_eval<3,4,5,6,0,1,2>;//qbm lm lbp qp p m m
//case 46757:	 return &A2q3g2l_qpmmmqbmlplbm_eval<2,3,4,5,6,0,1>;//lp lbm qp m m m qbm
case 46762:	 return &A2q3g2l_qpmmmqbmlmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qp m m m qbm
//case 47405:	 return &A2q3g2l_qppmmqbmlplbm_eval<2,3,4,5,6,0,1>;//lp lbm qp p m m qbm
case 47410:	 return &A2q3g2l_qppmmqbmlmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qp p m m qbm
//case 50645:	 return &A2q3g2l_qpmpmqbmlplbm_eval<2,3,4,5,6,0,1>;//lp lbm qp m p m qbm
case 50650:	 return &A2q3g2l_qpmpmqbmlmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qp m p m qbm
//case 51293:	 return &A2q3g2l_qpppmqbmlplbm_eval<2,3,4,5,6,0,1>;//lp lbm qp p p m qbm
case 51298:	 return &A2q3g2l_qpppmqbmlmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qp p p m qbm
//case 607:	 return &A2q3g2l_qpmmmqbmlplbm_eval<3,4,5,6,0,1,2>;//qbm lp lbm qp m m m
case 637:	 return &A2q3g2l_qpmmmqbmlmlbp_eval<3,4,5,6,0,1,2>;//qbm lm lbp qp m m m
//case 70085:	 return &A2q3g2l_qpmmpqbmlplbm_eval<2,3,4,5,6,0,1>;//lp lbm qp m m p qbm
case 70090:	 return &A2q3g2l_qpmmpqbmlmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qp m m p qbm
//case 70733:	 return &A2q3g2l_qppmpqbmlplbm_eval<2,3,4,5,6,0,1>;//lp lbm qp p m p qbm
case 70738:	 return &A2q3g2l_qppmpqbmlmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qp p m p qbm
//case 73973:	 return &A2q3g2l_qpmppqbmlplbm_eval<2,3,4,5,6,0,1>;//lp lbm qp m p p qbm
case 73978:	 return &A2q3g2l_qpmppqbmlmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qp m p p qbm
//case 74621:	 return &A2q3g2l_qppppqbmlplbm_eval<2,3,4,5,6,0,1>;//lp lbm qp p p p qbm
case 74626:	 return &A2q3g2l_qppppqbmlmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qp p p p qbm
//case 84672:	 return &A2q3g2l_qmmmmqbplplbm_eval<6,0,1,2,3,4,5>;//m m m qbp lp lbm qm
//case 84675:	 return &A2q3g2l_qmpmmqbplplbm_eval<6,0,1,2,3,4,5>;//p m m qbp lp lbm qm
//case 84690:	 return &A2q3g2l_qmmpmqbplplbm_eval<6,0,1,2,3,4,5>;//m p m qbp lp lbm qm
//case 84693:	 return &A2q3g2l_qmppmqbplplbm_eval<6,0,1,2,3,4,5>;//p p m qbp lp lbm qm
//case 84780:	 return &A2q3g2l_qmmmpqbplplbm_eval<6,0,1,2,3,4,5>;//m m p qbp lp lbm qm
//case 84783:	 return &A2q3g2l_qmpmpqbplplbm_eval<6,0,1,2,3,4,5>;//p m p qbp lp lbm qm
//case 84798:	 return &A2q3g2l_qmmppqbplplbm_eval<6,0,1,2,3,4,5>;//m p p qbp lp lbm qm
//case 84801:	 return &A2q3g2l_qmpppqbplplbm_eval<6,0,1,2,3,4,5>;//p p p qbp lp lbm qm
case 91152:	 return &A2q3g2l_qmmmmqbplmlbp_eval<6,0,1,2,3,4,5>;//m m m qbp lm lbp qm
case 91155:	 return &A2q3g2l_qmpmmqbplmlbp_eval<6,0,1,2,3,4,5>;//p m m qbp lm lbp qm
case 91170:	 return &A2q3g2l_qmmpmqbplmlbp_eval<6,0,1,2,3,4,5>;//m p m qbp lm lbp qm
case 91173:	 return &A2q3g2l_qmppmqbplmlbp_eval<6,0,1,2,3,4,5>;//p p m qbp lm lbp qm
case 91260:	 return &A2q3g2l_qmmmpqbplmlbp_eval<6,0,1,2,3,4,5>;//m m p qbp lm lbp qm
case 91263:	 return &A2q3g2l_qmpmpqbplmlbp_eval<6,0,1,2,3,4,5>;//p m p qbp lm lbp qm
case 91278:	 return &A2q3g2l_qmmppqbplmlbp_eval<6,0,1,2,3,4,5>;//m p p qbp lm lbp qm
case 91281:	 return &A2q3g2l_qmpppqbplmlbp_eval<6,0,1,2,3,4,5>;//p p p qbp lm lbp qm
//case 93377:	 return &A2q3g2l_qmmmmqbplplbm_eval<2,3,4,5,6,0,1>;//lp lbm qm m m m qbp
case 93382:	 return &A2q3g2l_qmmmmqbplmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qm m m m qbp
//case 94025:	 return &A2q3g2l_qmpmmqbplplbm_eval<2,3,4,5,6,0,1>;//lp lbm qm p m m qbp
case 94030:	 return &A2q3g2l_qmpmmqbplmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qm p m m qbp
//case 97265:	 return &A2q3g2l_qmmpmqbplplbm_eval<2,3,4,5,6,0,1>;//lp lbm qm m p m qbp
case 97270:	 return &A2q3g2l_qmmpmqbplmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qm m p m qbp
//case 97913:	 return &A2q3g2l_qmppmqbplplbm_eval<2,3,4,5,6,0,1>;//lp lbm qm p p m qbp
case 97918:	 return &A2q3g2l_qmppmqbplmlbp_eval<2,3,4,5,6,0,1>;//lm lbp qm p p m qbp
	
	default: // We return the zero pointer for all other helicity combinations
		//cout << " A2q3g2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}

template complex<R> (*A2q3g2l_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q3g2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q3g2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> (*A2q3g2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}

